#!/usr/bin/env Rscript
# Env smoke check: everything the entrypoints load, plus the RcppML version gate.
suppressPackageStartupMessages({
  library(Matrix); library(HDF5Array); library(argparser)
  library(data.table); library(jsonlite)
})
source("src/rcppml.R")
load_rcppml()
cat(sprintf("gpu_available: %s\n", RcppML::gpu_available()))
cat("OK\n")
