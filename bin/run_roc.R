#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", all_args[grep("--file=", all_args)]))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: run_roc.R <path/to/padj>")

truth_path <- args[1]
negatives_path <- args[2]
path_to_padj <- args[3]

truth <- readLines(truth_path)
negatives <- readLines(negatives_path)


thresholds <- 10^seq(-6, 0, length.out = 40)
files <- list.files(path_to_padj, pattern = "\\.rds$", full.names = TRUE)

roc_rows <- list()
for (f in files) {
  padj  <- readRDS(f)
  parts <- strsplit(sub("\\.rds$", "", basename(f)), "_")[[1]]   # e.g. "draw05" "n3" "10M"
  for (t in thresholds) {
    called <- names(padj)[which(padj < t)]
    roc_rows[[length(roc_rows) + 1]] <- data.frame(
    draw  = as.integer(sub("draw", "", parts[1])),
    n = as.integer(sub("n", "", parts[2])),
    depth = parts[3],
    threshold = t,
    tpr = length(intersect(called, truth))     / length(truth),
    fpr = length(intersect(called, negatives)) / length(negatives))
  }
}
roc <- do.call(rbind, roc_rows)
write.csv(roc, "roc_runs.csv", row.names = FALSE)
