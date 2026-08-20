# rcppml

Omnibenchmark **module** wrapping [zdebruine/RcppML](https://github.com/zdebruine/RcppML)
1.0.0. Two entrypoints, both emitting the standard `{dataset}_pcas.tsv` +
`{dataset}_loadings.tsv` pair:

| entrypoint | function | `--solver` / `--loss` |
|---|---|---|
| `pca` | `RcppML::svd(center = TRUE)` | `lanczos`, `krylov` (also `irlba`, `randomized`, `deflation`) |
| `nmf` | `RcppML::nmf` | `mse`, `gp`, `nb`, `gamma`, `inverse_gaussian`, `tweedie` |

Both also take `--backend cpu|gpu` (default `cpu`).

## `--backend`

Maps onto RcppML's `resource`. The module never passes RcppML's own `"auto"`,
which picks a device from whatever the host happens to have. Default is `cpu`.

`--backend gpu` needs `RcppML_gpu.so` plus a visible CUDA device. The
`r-rcppml` conda package ships that library (a plain `R CMD INSTALL` /
`install_github` build does not, which is why this used to be unusable). RcppML
itself only *warns* when a GPU is requested and missing, then computes on the
CPU — so both entrypoints check `gpu_available()` and exit non-zero instead,
before reading the matrix. A CPU run recorded as a GPU arm is worse than a
failed job.

Measured on an RTX 2000 Ada laptop GPU vs 8 CPU threads, k = 50, solver time
only (the h5 read is unchanged and dominates wall clock on small inputs):

| input (genes x cells, nnz) | svd cpu → gpu | nmf cpu → gpu |
|---|---|---|
| be1-fixture 2000 x 1715, 1.1M | 0.17 → 0.47 s (slower) | 1.2 → 0.7 s |
| tenx-0020k 2000 x 19696, 13.4M | 1.34 → 0.44 s (3x) | 8.2 → 1.4 s (6x) |
| pbmc 2000 x 156881, 96.4M | 9.7 → 2.5 s (3.8x) | 121 → 14.4 s (8.4x) |

Below ~10M nonzeros the CUDA context setup (~0.3 s) eats the gain; NMF profits
more than SVD at every size. CPU and GPU embeddings agree to |cor| > 1 - 1e-7
per PC on tenx-0020k.

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
`pixi run check` verifies the env imports and reports `gpu_available`.

RcppML 1.0.0 comes from the [`almost-conductor`][chan] channel, listed ahead of
conda-forge, which stops at 0.3.7.1 — that version predates `svd()`, `pca()` and
the lanczos/krylov backends entirely. The module used to build the pinned
GitHub commit into a local `.Rlib/` on first use; the package replaced that.

[chan]: https://prefix.dev/channels/almost-conductor/packages/r-rcppml

## Thread count

`RcppML`'s own default is `threads = 0` = every core, which would ignore the
stage's `resources: cores:`. Both entrypoints pass `OMP_NUM_THREADS` through
instead, which is what snakemake exports from the rule's thread count.

## NMF refuses negative input

NMF is only defined on non-negative data. `nr-scanpy`'s
`analyticalPearsonResiduals` normalisation produces negatives, and RcppML would
fit it anyway, so `nmf.R` exits non-zero instead of contributing a meaningless
embedding. Pair it with a `log1pCP10k`-style NORM arm.

## Test

    tests/smoke.sh <dataset>_normalized_selected.h5

Runs both entrypoints and checks the two TSVs against the stage contract.

## Usage in omni-scrna

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
            backend: cpu
```

The `nmf` entrypoint takes the same single input and emits the same two
outputs, so it drops into the `PCA` stage the same way (`entrypoint: nmf`,
`loss: [mse]`). 
