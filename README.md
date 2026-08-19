# rcppml

Omnibenchmark **module** wrapping [zdebruine/RcppML](https://github.com/zdebruine/RcppML)
1.0.0. Two entrypoints, both emitting the standard `{dataset}_pcas.tsv` +
`{dataset}_loadings.tsv` pair:

| entrypoint | function | `--solver` / `--loss` |
|---|---|---|
| `pca` | `RcppML::svd(center = TRUE)` | `lanczos`, `krylov` (also `irlba`, `randomized`, `deflation`) |
| `nmf` | `RcppML::nmf` | `mse`, `gp`, `nb`, `gamma`, `inverse_gaussian`, `tweedie` |

No algorithm lives here — this is I/O plus the orientation bookkeeping.

## Orientation

RcppML takes the gene-by-cell matrix as read; `center = TRUE` subtracts **row**
means, i.e. per-gene means, which is what scanpy / scrapper / irlba also do.
Nothing is transposed, so sparsity survives.

    svd:  A (genes x cells) = u d v'   →  scores = v %*% diag(d),  loadings = u
    nmf:  A                 = w d h    →  scores = t(d * h),       loadings = w

Both entrypoints fold `d` into the scores, so the embedding carries the factor
scale. Column names stay `PC1..PCk` in both, including for NMF factors, because
the downstream metrics and validators key off that prefix.

## Dependencies

`pixi.toml` is the source of truth; `envs/rcppml.yml` is generated from it with
`pixi run export-env` and is what omnibenchmark's conda backend actually
consumes (the plan sets `software_backend: conda`). Do not hand-edit the yml.
`pixi run check` verifies the env imports.

### RcppML itself is installed at run time, not by conda

conda-forge stops at `r-rcppml` 0.3.7.1, which predates `svd()`, `pca()` and the
lanczos/krylov backends entirely. `envs/rcppml.yml` therefore ships only the
toolchain, and `src/rcppml.R` installs the pinned commit
(`df69ddd`, v1.0.0) into a module-local `.Rlib/` on first use — the same
build-on-first-use trick `omni-rmt-spca` uses for its Rust crate. Cached per
module clone, so it is paid once, not once per job. First run costs ~10 min of
C++ compilation.

## Thread count

`RcppML`'s own default is `threads = 0` = every core, which would ignore the
stage's `resources: cores:`. Both entrypoints pass `OMP_NUM_THREADS` through
instead, which is what snakemake exports from the rule's thread count.

## NMF refuses negative input

NMF is only defined on non-negative data. `nr-scanpy`'s
`analyticalPearsonResiduals` normalisation produces negatives, and RcppML would
fit it anyway, so `nmf.R` exits non-zero instead of contributing a meaningless
embedding. Pair it with a `log1pCP10k`-style NORM arm.

## Measured on be1-fixture

2000 genes x 1715 cells, `n_components: 50`, `random_seed: 42`,
`OMP_NUM_THREADS=4`, against the sibling runs on the identical input:

| arm | wall | peak RSS | min &#124;cor&#124; vs scanpy arpack, PC1-50 |
|---|---|---|---|
| `pca --solver lanczos` | 6.1 s | 545 MB | 0.999886 (median 1.000000) |
| `pca --solver krylov`  | 5.8 s | 545 MB | 0.999886 (median 1.000000) |
| `nmf --loss mse`       | 7.4 s | 554 MB | n/a |

Same numbers against `pc-scrapper solver-exact`. Only PC50 falls below
1.000000 — the usual tail-component wobble, not a disagreement about the
subspace. lanczos and krylov agree with each other to 1.000000 on all 50, so at
this size krylov's refinement has nothing left to do after its Lanczos seed.

The NMF embedding is 48% exact zeros, both factors non-negative. Mean
silhouette of the true cell lines: 0.265 (PCA) vs 0.275 (NMF) over all 50 dims,
0.392 vs 0.172 over the first 10 — NMF factors are not variance-ordered, so
truncating them is not the same operation as truncating PCs.

## Test

    tests/smoke.sh <dataset>_normalized_selected.h5

Runs both entrypoints and checks the two TSVs against the stage contract.

## Wiring into split-stages-plan

Add `envs/rcppml.yml` to the plan's `software_environments`, then, under the
`PCA` stage:

```yaml
      - id: pc-rcppml
        name: "RcppML truncated SVD"
        software_environment: "rcppml"
        exclude:
          - d-feat-sel
          - d-annotation
        repository:
          url: https://github.com/omni-scrna/rcppml
          commit: <sha>
          entrypoint: pca
        parameters:
          - solver: [lanczos, krylov]
            n_components: 50
            random_seed: 42
```

The `nmf` entrypoint takes the same single input and emits the same two
outputs, so it drops into the `PCA` stage the same way (`entrypoint: nmf`,
`loss: [mse]`). Putting it in `CNTFCT` instead means also accepting that
stage's other two inputs (`--rawdata_h5ad`, `--filtered_cellids`), which the
module does not currently declare.
