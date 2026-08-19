# Shared helpers for the RcppML entrypoints (pca.R, nmf.R).
#
# RcppML 1.0.0 is where svd()/pca() and the lanczos/krylov backends live, and
# it is GitHub-only -- conda-forge still ships 0.3.7.1, which has none of them.
# So the env carries the toolchain and we install the pinned commit into a
# module-local lib on first use, the same build-on-first-use trick
# omni-rmt-spca uses for its Rust crate. Cached in the module clone, so it is
# paid once per clone, not once per job.

RCPPML_REF <- "df69dddbe37962cd0179b6b6c6ace94e400b4f60"  # zdebruine/RcppML @ v1.0.0
RCPPML_MIN <- "1.0.0"

module_dir <- function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(normalizePath(sub("^--file=", "", cargs[[m]]))) else getwd()
}

# ponytail: version check, not ref check. Bumping RCPPML_REF within 1.0.0 will
# not invalidate a cached lib -- wipe .Rlib by hand if that ever matters.
load_rcppml <- function() {
  lib <- file.path(module_dir(), ".Rlib")
  dir.create(lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(lib, .libPaths()))
  have <- requireNamespace("RcppML", quietly = TRUE) &&
    utils::packageVersion("RcppML") >= RCPPML_MIN
  if (!have) {
    cat(sprintf("LOG: installing RcppML@%s into %s (first use, ~minutes)\n",
                substr(RCPPML_REF, 1, 7), lib))
    remotes::install_github(paste0("zdebruine/RcppML@", RCPPML_REF),
                            lib = lib, upgrade = "never")
  }
  # No lib.loc: .Rlib is already first on .libPaths(), and pinning to it would
  # reject an RcppML supplied by conda (r-rcppml on almost-conductor), which is
  # the preferred source once it is published.
  suppressPackageStartupMessages(library(RcppML, quietly = TRUE))
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
           "RcppML_gpu", .Platform$dynlib.ext, ", which comes from a separate ",
           "`make -f src/Makefile.gpu` against a CUDA toolkit -- neither the ",
           "conda r-rcppml package nor this module's .Rlib build produces it. ",
           "See the README.", call. = FALSE)
    info <- RcppML::gpu_info()
    for (i in seq_len(nrow(info)))
      cat(sprintf("LOG: gpu %d: %s, %d MB total, %d MB free\n",
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
