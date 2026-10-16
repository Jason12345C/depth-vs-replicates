library(tximport)
library(DESeq2)

if (file.exists("ref/tx2gene_r116.csv")) {
  tx2gene <- read.csv("ref/tx2gene_r116.csv")
  
} else {source("R/tx2gene.R")}

sample_sheet <- read.csv("sample_sheets/samplesheet_all.csv")
coldata <- unique(sample_sheet[, c("sample", "condition")]) # get all unique samples (should be 14)
stopifnot(nrow(coldata) == 14)
rownames(coldata) <- coldata$sample # set the index to unique sample names (ex: e2_0, e2_1)
coldata$condition <- factor(coldata$condition, levels = c("Control", "E2"))   # Control is the reference

depths <- c("2.5M" = "2500000",
            "5M"   = "5000000",
            "10M"  = "10000000",
            "15M"  = "15000000",
            "20M"  = "20000000",
            "25M"  = "25000000",
            "full" = "full")

# helper for different depths runs
getDESeq <- function(depth, samples) {
  subs <- coldata[samples, , drop = FALSE]  
  files <- file.path("results/kallisto", paste0(subs$sample, "_", depths[depth]), "abundance.tsv")
  names(files) <- subs$sample
  stopifnot(all(file.exists(files)))  
  txi <- suppressMessages(tximport(files, type = "kallisto", tx2gene = tx2gene[, 1:2])) # txi wants only transcript and gene ID
  dds <- suppressMessages(DESeqDataSetFromTximport(txi, colData = subs, design = ~condition))
  dds <- dds[rowSums(counts(dds)) >= 10, ] # add up gene reads for the samples and filter out with < 10 reads
  return(DESeq(dds, quiet = TRUE)) #suppress outputs
}
