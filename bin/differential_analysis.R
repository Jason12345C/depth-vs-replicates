#!/usr/bin/env Rscript
# Full-depth DE analysis: E2 vs Control, all 14 samples.
# Usage: differential_analysis.R <samplesheet.csv> <tx2gene.csv>
# Expects the <sample>_full kallisto folders in the current folder.

all_args   <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", all_args[grep("--file=", all_args)]))
source(file.path(script_dir, "common.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: differential_analysis.R <samplesheet.csv> <tx2gene.csv>")
coldata <- load_coldata(args[1])
tx2gene <- read.csv(args[2])
symbols <- unique(tx2gene[, c("GENEID", "SYMBOL")])

depths <- c("2.5M" = "2500000",
            "5M"   = "5000000",
            "10M"  = "10000000",
            "15M"  = "15000000",
            "20M"  = "20000000",
            "25M"  = "25000000",
            "full" = "full")
            
# Fit and test
dds_full <- getDESeq("full", coldata$sample, coldata, tx2gene)
res_full <- results(dds_full, name = "condition_E2_vs_Control", alpha = 0.05)

df <- as.data.frame(res_full)
df$GENEID <- rownames(df)
df <- merge(df, symbols, by = "GENEID")
write.csv(df, "res_full.csv", row.names = FALSE)

# Shrunken fold changes, for ranking and the MA plot
resLFC <- as.data.frame(lfcShrink(dds_full, coef = "condition_E2_vs_Control",
                                  res = res_full, type = "apeglm"))
resLFC$GENEID <- rownames(resLFC)
resLFC <- merge(resLFC, symbols, by = "GENEID")

genes_to_highlight <- c("GREB1", "PGR", "TFF1")
resLFC$category <- "Not Significant"
resLFC$category[which(resLFC$padj < 0.05 & resLFC$log2FoldChange > 0)] <- "Up"
resLFC$category[which(resLFC$padj < 0.05 & resLFC$log2FoldChange < 0)] <- "Down"
resLFC$category[resLFC$SYMBOL %in% genes_to_highlight] <- "Known ER Targets"

p <- ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color = category)) +
  geom_point(alpha = 0.6) +
  geom_point(data = subset(resLFC, SYMBOL %in% genes_to_highlight), size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Not Significant" = "grey", "Up" = "red",
                                "Down" = "blue", "Known ER Targets" = "darkred")) +
  theme_minimal() +
  geom_text_repel(data = subset(resLFC, SYMBOL %in% genes_to_highlight),
                  aes(label = SYMBOL), color = "black", fontface = "bold") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        text = element_text(color = "black"), axis.text = element_text(color = "black")) +
  labs(x = "Mean Normalized Counts", y = "log2 fold change (E2 vs Control)", color = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
ggsave("MAplot.png", p, width = 8, height = 6, dpi = 200)

# Significant genes, ranked by shrunken fold change
sig <- subset(resLFC, padj < 0.05)
cols <- c("GENEID", "SYMBOL", "log2FoldChange", "padj", "baseMean")
write.csv(sig[order(-sig$log2FoldChange), cols], "DE_genes.csv", row.names = FALSE)
