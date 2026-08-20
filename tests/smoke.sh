#!/bin/sh
# One runnable check: both entrypoints against a real normalized_selected.h5.
#
#   tests/smoke.sh path/to/<dataset>_normalized_selected.h5
#
# Run from the module root, in the rcppml env.
set -eu
cd "$(dirname "$0")/.."
H5=${1:?usage: tests/smoke.sh <normalized_selected.h5>}
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

Rscript pca.R --output_dir "$OUT" --name smoke --normalized_selected_h5 "$H5" \
  --solver lanczos --n_components 10 --random_seed 42
Rscript nmf.R --output_dir "$OUT" --name smoke --normalized_selected_h5 "$H5" \
  --n_components 10 --random_seed 42 --loss mse

# The stage contract: two TSVs, cell_id/gene_id key plus PC1..PC10, no NaNs.
Rscript -e '
  for (f in c("pcas", "loadings")) {
    d <- read.delim(file.path("'"$OUT"'", paste0("smoke_", f, ".tsv")), check.names = FALSE)
    stopifnot(nrow(d) > 0, ncol(d) == 11, colnames(d)[11] == "PC10",
              !anyNA(d), is.character(d[[1]]))
    cat("OK", f, nrow(d), "x", ncol(d), "\n")
  }'
echo "smoke OK"
