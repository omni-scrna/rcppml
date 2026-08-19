#!/usr/bin/env Rscript
# Env smoke check: the conda side only. RcppML itself is not a conda dep -- it is
# built on first use into .Rlib (see src/rcppml.R), so this just reports it.
suppressPackageStartupMessages({
  library(Matrix); library(HDF5Array); library(argparser)
  library(data.table); library(jsonlite); library(remotes)
})
source("src/rcppml.R")
lib <- file.path(module_dir(), ".Rlib")
cat(sprintf("RcppML: %s\n", tryCatch(
  as.character(utils::packageVersion("RcppML", lib.loc = lib)),
  error = function(e) sprintf("not built yet (pinned %s)", substr(RCPPML_REF, 1, 7)))))
cat("OK\n")
