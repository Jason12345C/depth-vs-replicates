# common.R: shared helpers for the DESeq2 scripts.
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

# One row per biological sample, with Control as the reference level
load_coldata <- function(samplesheet) {
  sheet   <- read.csv(samplesheet)
  coldata <- unique(sheet[, c("sample", "condition")])
  stopifnot(nrow(coldata) == 14)
  rownames(coldata) <- coldata$sample
  coldata$condition <- factor(coldata$condition, levels = c("Control", "E2"))
  coldata
}

# helper for different depths runs
getDESeq <- function(depth, samples, coldata, tx2gene, quant_dir = ".") {
  subs  <- coldata[samples, , drop = FALSE]
  files <- file.path(quant_dir, paste0(subs$sample, "_", depths[depth]), "abundance.tsv")
  names(files) <- subs$sample
  stopifnot(all(file.exists(files)))
  txi <- suppressMessages(tximport(files, type = "kallisto", tx2gene = tx2gene[, 1:2]))
  dds <- suppressMessages(DESeqDataSetFromTximport(txi, colData = subs, design = ~condition))
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  DESeq(dds, quiet = TRUE)
}
