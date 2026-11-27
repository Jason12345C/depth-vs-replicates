#!/usr/bin/env Rscript

#compare total DE genes vs. truth set for varying depths
#Usage: power_analysis.R -mode truth <samplesheet.csv> <tx2gene>
#Usage: power_analysis.R -mode job <samplesheet.csv> <tx2gene> <draw> <n> <truth.txt> <negatives.txt>
#Usage: power_analysis.R -mode combine <truth.txt> <negatives.txt>
suppressMessages({library(rhdf5)})
all_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", all_args[grep("--file=", all_args)]))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("usage: power_analysis.R <truth|job|combine> [args...]")
mode <- args[1]
rest <- args[-1]
depths <- c("2.5M" = "2500000",
            "5M"   = "5000000",
            "10M"  = "10000000",
            "15M"  = "15000000",
            "20M"  = "20000000",
            "25M"  = "25000000",
            "full" = "full")
sim_depths <- setdiff(names(depths), "full")
run_truth   <- function(samplesheet, tx2gene_path) {
  coldata <- load_coldata(samplesheet)
  tx2gene <- read.csv(tx2gene_path)
  dds_truth <- getDESeq("25M", coldata$sample, coldata, tx2gene)
  res_truth <- results(dds_truth,
                       name = "condition_E2_vs_Control",
                       alpha = 0.001)
  truth <- rownames(res_truth)[which(res_truth$padj < 0.001)]
  negatives <- rownames(res_truth)[which(res_truth$padj >= 0.1)]
  writeLines(truth, "truth.txt")
  writeLines(negatives, "negatives.txt")
}

run_job <- function(samplesheet,
                    tx2gene_path,
                    truth_path,
                    negatives_path,
                    draw, n_rep) {

  coldata <- load_coldata(samplesheet)
  tx2gene <- read.csv(tx2gene_path)

  # Split the 14 samples into the two groups
  ctrl_ids <- coldata$sample[coldata$condition == "Control"]
  e2_ids   <- coldata$sample[coldata$condition == "E2"]

  # Randomly pick n samples (this simulates # of replicates) from each group
  # seeding by draw makes each draw repeatable
  set.seed(draw)
  picks <- c(sample(ctrl_ids, n_rep), sample(e2_ids, n_rep))

  # read truth and negatives from file, need both to calculate tp and fp
  truth <- readLines(truth_path)
  negatives <- readLines(negatives_path)

  out <- list()
  for (depth in sim_depths) {
    res_sim <- results(getDESeq(depth, picks, coldata, tx2gene),
                       name = "condition_E2_vs_Control",
                       alpha = 0.05)
    # save every gene's padj for this run
    padj <- setNames(res_sim$padj, rownames(res_sim))
    saveRDS(padj,
            sprintf("draw%02d_n%d_%s.rds", draw, n_rep, depth))
    called  <- rownames(res_sim)[which(res_sim$padj < 0.05)]
    tp <- length(intersect(called, truth))
    fp <- length(intersect(called, negatives))
    # populate dataframe
    out[[length(out) + 1]] <- data.frame(
      depth = depth,
      n = n_rep,
      draw = draw,
      samples   = paste(picks, collapse = ","),
      n_called  = length(called), 
      tp = tp,
      fp = fp,
      fpr = fp/length(negatives), #false positive rate FP/(FP+TN)
      power = tp / length(truth), # tpr
      precision = if (tp + fp > 0) tp / (tp + fp) else NA)
  }
  job_table <- do.call(rbind, out)
  write.csv(job_table,
            sprintf("power_draw%02d_n%d.csv", draw, n_rep),
            row.names = FALSE)

}

switch(mode,
  truth = run_truth(samplesheet = rest[1], tx2gene_path = rest[2]),
  job   = run_job(samplesheet = rest[1], tx2gene_path = rest[2],
                  draw = as.integer(rest[3]), n_rep = as.integer(rest[4]),
                  truth_path = rest[5], negatives_path = rest[6]),
  stop("unknown mode: ", mode)
)

# jobs <- expand.grid(draw = 1:N_DRAWS, n_rep = 2:MAX_REPLICATES)
# jobs <- jobs[!(jobs$n_rep == 7 & jobs$draw > 1), ]    # n = 7 only once

# plan(multisession, workers = 10)

# with_progress({
#   p <- progressor(steps = nrow(jobs))
#   results_list <- future_lapply(seq_len(nrow(jobs)), function(i) {
#     out <- run_job(jobs$draw[i], jobs$n_rep[i])
#     p()                       # tick the progress bar once per finished job
#     out
#   }, future.seed = TRUE, future.packages = c("tximport", "DESeq2"))
# })

# plan(sequential)                               # shut the workers down

# runs <- do.call(rbind, results_list)
# write.csv(runs, "results/deseq2/power_runs.csv", row.names = FALSE)
