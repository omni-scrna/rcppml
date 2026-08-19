#!/usr/bin/env Rscript
# Env smoke check: the conda side only. RcppML itself is not a conda dep -- it is
# built on first use into .Rlib (see src/rcppml.R), so this just reports it.
suppressPackageStartupMessages({
  library(Matrix); library(HDF5Array); library(argparser)
  library(data.table); library(jsonlite); library(remotes)
})
source("src/rcppml.R")
.libPaths(c(file.path(module_dir(), ".Rlib"), .libPaths()))
cat(sprintf("RcppML: %s\n", tryCatch(
  sprintf("%s (%s)", utils::packageVersion("RcppML"),
          dirname(dirname(system.file(package = "RcppML")))),
  error = function(e) sprintf("not present (would build pinned %s)",
                              substr(RCPPML_REF, 1, 7)))))
cat("OK\n")
