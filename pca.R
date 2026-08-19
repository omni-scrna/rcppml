#!/usr/bin/env Rscript
# PCA module (RcppML-backed) for omnibenchmark.
#
# RcppML::svd() takes the gene-by-cell matrix as-is: `center = TRUE` subtracts
# *row* means, i.e. per-gene means, which is the standard scRNA PCA convention
# the sibling modules (scanpy / scrapper / irlba) also use. No transpose, so
# sparsity survives.
#
# With A (genes x cells) = u d v', u is genes x k and v is cells x k:
#   scores   = v %*% diag(d)   (cells x k)
#   loadings = u               (genes x k)
# That is the mirror of the cells-as-rows modules, where u and v swap roles.
#
# --solver is handed straight to RcppML's `method` -- it validates and names the
# accepted set in its own error. The two of interest here:
#   lanczos  unconstrained Lanczos bidiagonalisation
#   krylov   Krylov-Seeded Projected Refinement, block: all k factors at once

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

source("src/common/cli.R")
source("src/rcppml.R")

p <- arg_parser("PCA module (RcppML)")
p <- add_base_args(p)                 # --output_dir, --name
p <- add_stage_args(p, "PCA")         # --normalized_selected_h5
p <- add_argument(p, "--solver", type = "character",
                  help = "RcppML method: lanczos, krylov, irlba, randomized, deflation")
p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
p <- add_argument(p, "--random_seed", type = "integer", help = "seed")
args <- parse_args(p)

run_pca <- function(X, args) {
  fit <- RcppML::svd(X, k = args$n_components, center = TRUE,
                     method = args$solver, seed = args$random_seed,
                     threads = omp_threads(), verbose = TRUE)
  list(embedding = fit@v %*% diag(fit@d, nrow = length(fit@d)),
       loadings  = fit@u)
}

main <- function() {
  log_args(args)
  load_rcppml()
  m <- read_tenx(args$normalized_selected_h5)

  res <- run_pca(m, args)
  rownames(res$embedding) <- colnames(m)
  rownames(res$loadings)  <- rownames(m)
  write_factors(res$embedding, res$loadings, args)
}

if (sys.nframe() == 0L) main()
