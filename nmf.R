#!/usr/bin/env Rscript
# NMF module (RcppML-backed) for omnibenchmark.
#
# RcppML::nmf() factorises A = w d h by alternating least squares with a
# coordinate-descent NNLS solver. A is the gene-by-cell matrix as read, so:
#   loadings = w              (genes x k, non-negative gene programmes)
#   scores   = t(diag(d) h)   (cells x k, per-cell programme usage)
# The d diagonal is folded into the scores so the embedding carries the factor
# scale, matching what the PCA entrypoint emits.
#
# --loss goes straight to RcppML, which validates it: mse (Gaussian) plus the
# count-aware IRLS losses gp / nb / gamma / inverse_gaussian / tweedie.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

source("src/common/cli.R")
source("src/rcppml.R")

p <- arg_parser("NMF module (RcppML)")
p <- add_base_args(p)                 # --output_dir, --name
p <- add_stage_args(p, "PCA")         # --normalized_selected_h5
p <- add_argument(p, "--n_components", type = "integer", help = "rank k")
p <- add_argument(p, "--random_seed", type = "integer", help = "seed")
p <- add_argument(p, "--backend", type = "character", default = "cpu",
                  help = "compute backend: cpu or gpu")
p <- add_argument(p, "--loss", type = "character", default = "mse",
                  help = "mse, gp, nb, gamma, inverse_gaussian or tweedie")
args <- parse_args(p)

run_nmf <- function(X, args) {
  # NMF is only defined on non-negative input. Pearson-residual normalisation
  # (nr-scanpy analyticalPearsonResiduals) produces negatives, and RcppML would
  # otherwise fit it silently, so fail loudly instead of benchmarking nonsense.
  if (min(X) < 0)
    stop("input has negative entries (min ", min(X), "); NMF needs a ",
         "non-negative normalisation such as log1pCP10k", call. = FALSE)

  # verbose = FALSE on purpose: RcppML prints a per-iteration timing table,
  # ~100 lines a job, which buries the stage log.
  fit <- RcppML::nmf(X, k = args$n_components, loss = args$loss,
                     seed = args$random_seed,
                     resource = args$backend,
                     threads = omp_threads(), verbose = FALSE)
  list(embedding = t(fit@d * fit@h),   # d recycles down the k rows of h
       loadings  = fit@w)
}

main <- function() {
  log_args(args)
  load_rcppml()
  args$backend <- resolve_backend(args$backend)   # before the matrix read, not after
  m <- read_tenx(args$normalized_selected_h5)

  res <- run_nmf(m, args)
  rownames(res$embedding) <- colnames(m)
  rownames(res$loadings)  <- rownames(m)
  write_factors(res$embedding, res$loadings, args)
}

if (sys.nframe() == 0L) main()
