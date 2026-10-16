rm(list=ls())
library(future.apply)
library(progressr)
handlers("progress")
options(progressr.enable = TRUE)

source("R/common.R")

ctrl_ids <- coldata$sample[coldata$condition == "Control"]
e2_ids   <- coldata$sample[coldata$condition == "E2"]

getSamples <- function(n_replicates, seed) {
  set.seed(seed)
  return (c(sample(ctrl_ids, n_replicates), sample(e2_ids, n_replicates)))
}

dds_truth <- getDESeq("25M", coldata$sample)
res_truth <- results(dds_truth, name = "condition_E2_vs_Control", alpha=0.001)
truth <- rownames(res_truth)[which(res_truth$padj < 0.001)]

thresholds <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5)
negatives <- rownames(res_truth)[which(res_truth$padj >= 0.1)]

N_DRAWS = 20
MAX_REPLICATES = 7
dir.create("results/deseq2/padj", showWarnings = FALSE, recursive = TRUE)

run_job <- function(draw, n_rep) {
  picks <- getSamples(n_rep, seed = draw)
  out <- list()
  for (depth in setdiff(names(depths), "full")) {
    res_sim <- results(getDESeq(depth, picks), name = "condition_E2_vs_Control", alpha = 0.05)
    
    # save every gene's padj for this run
    padj <- setNames(res_sim$padj, rownames(res_sim))
    saveRDS(padj, sprintf("results/deseq2/padj/draw%02d_n%d_%s.rds", draw, n_rep, depth))
    
    called  <- rownames(res_sim)[which(res_sim$padj < 0.05)]
    tp <- length(intersect(called, truth))
    fp <- length(intersect(called, negatives))
    
    # populate dataframe
    out[[length(out) + 1]] <- data.frame(
      depth = depth, n = n_rep, draw = draw,
      samples   = paste(picks, collapse = ","),
      n_called  = length(called), 
      tp = tp,
      fp = fp,
      fpr = fp/length(negatives), #false positive rate FP/(FP+TN)
      power = tp / length(truth), # tpr
      precision = if (tp + fp > 0) tp / (tp + fp) else NA)
  }
  do.call(rbind, out)
}

jobs <- expand.grid(draw = 1:N_DRAWS, n_rep = 2:MAX_REPLICATES)
jobs <- jobs[!(jobs$n_rep == 7 & jobs$draw > 1), ]    # n = 7 only once

plan(multisession, workers = 10)

with_progress({
  p <- progressor(steps = nrow(jobs))
  results_list <- future_lapply(seq_len(nrow(jobs)), function(i) {
    out <- run_job(jobs$draw[i], jobs$n_rep[i])
    p()                       # tick the progress bar once per finished job
    out
  }, future.seed = TRUE, future.packages = c("tximport", "DESeq2"))
})

plan(sequential)                               # shut the workers down

runs <- do.call(rbind, results_list)
write.csv(runs, "results/deseq2/power_runs.csv", row.names = FALSE)


#################ROC################################
thresholds <- 10^seq(-6, 0, length.out = 40)
files <- list.files("results/deseq2/padj", full.names = TRUE)

roc_rows <- list()
for (f in files) {
  padj  <- readRDS(f)
  parts <- strsplit(sub("\\.rds$", "", basename(f)), "_")[[1]]   # e.g. "draw05" "n3" "10M"
  for (t in thresholds) {
    called <- names(padj)[which(padj < t)]
    roc_rows[[length(roc_rows) + 1]] <- data.frame(
      draw  = as.integer(sub("draw", "", parts[1])),
      n     = as.integer(sub("n", "", parts[2])),
      depth = parts[3],
      threshold = t,
      tpr = length(intersect(called, truth))     / length(truth),
      fpr = length(intersect(called, negatives)) / length(negatives))
  }
}
roc <- do.call(rbind, roc_rows)
write.csv(roc, "results/deseq2/roc_runs.csv", row.names = FALSE)
