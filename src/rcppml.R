# Shared helpers for the RcppML entrypoints (pca.R, nmf.R).

RCPPML_MIN <- "1.0.0"   # svd()/pca() and the lanczos/krylov backends start here

load_rcppml <- function() {
  suppressPackageStartupMessages(library(RcppML, quietly = TRUE))
  if (utils::packageVersion("RcppML") < RCPPML_MIN)
    stop("RcppML ", utils::packageVersion("RcppML"), " < ", RCPPML_MIN,
         "; the env must pull r-rcppml from the almost-conductor channel, not ",
         "conda-forge (which stops at 0.3.7.1).", call. = FALSE)
  cat(sprintf("LOG: RcppML %s\n", utils::packageVersion("RcppML")))
}

# --backend maps onto RcppML's `resource`. We never pass its "auto": that picks a
# device from whatever the host happens to have, so the same invocation would
# mean different things on different machines -- which a benchmark cannot have.
#
# RcppML only *warns* when a GPU is requested and missing, then runs on the CPU
# anyway. Refuse instead: a CPU run recorded as a GPU arm is worse than a failed
# job. Same stance omni-rmt-spca takes on its biwhitening fallback.
resolve_backend <- function(backend) {
  if (!backend %in% c("cpu", "gpu"))
    stop("--backend must be cpu or gpu, got: ", backend, call. = FALSE)
  if (backend == "gpu") {
    if (!RcppML::gpu_available())
      stop("--backend gpu, but RcppML reports no usable GPU. It needs ",
           "RcppML_gpu", .Platform$dynlib.ext, " plus a visible CUDA device. The ",
           "r-rcppml conda package ships that library; a source build via ",
           "`R CMD INSTALL` does not. See the README.", call. = FALSE)
    info <- RcppML::gpu_info()
    for (i in seq_len(nrow(info)))
      cat(sprintf("LOG: gpu %.0f: %s, %.0f MB total, %.0f MB free\n",
                  info$device[i], info$name[i], info$total_mb[i], info$free_mb[i]))
  }
  backend
}

# OpenMP threads. Snakemake exports OMP_NUM_THREADS = <rule threads>; RcppML's
# own default is 0 = grab every core, which would ignore the stage's `cores`.
omp_threads <- function() as.integer(Sys.getenv("OMP_NUM_THREADS", "0"))

read_tenx <- function(path) {
  m <- as(HDF5Array::TENxMatrix(path, group = "matrix"), "dgCMatrix")
  # Both output TSVs are keyed by these; a nameless matrix would silently emit
  # a table with the id column dropped.
  stopifnot(!is.null(rownames(m)), !is.null(colnames(m)))
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(m), ncol(m)))
  m
}

log_args <- function(args) {
  cat(sprintf("Full command: %s\n",
              paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  cat("LOG: command line args\n----------------------------------\n")
  for (i in seq_along(args)) cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
  cat("----------------------------------\n")
}

# Both stages (PCA, CNTFCT) declare the same two outputs. Column names stay
# PC*-prefixed even for NMF factors so the downstream metrics/validators, which
# key off that prefix, read either module unchanged.
write_factors <- function(embedding, loadings, args) {
  stopifnot(nrow(embedding) > 0, ncol(embedding) == ncol(loadings))
  colnames(embedding) <- colnames(loadings) <- paste0("PC", seq_len(ncol(embedding)))
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
              nrow(embedding), ncol(embedding), nrow(loadings), ncol(loadings)))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  out_e <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
  data.table::fwrite(data.frame(cell_id = rownames(embedding), embedding),
                     out_e, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_e))

  out_l <- file.path(args$output_dir, sprintf("%s_loadings.tsv", args$name))
  data.table::fwrite(data.frame(gene_id = rownames(loadings), loadings),
                     out_l, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_l))
}
