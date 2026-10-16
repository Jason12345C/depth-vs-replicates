rm(list=ls())
library(ggplot2)
library(ggrepel) # For non-overlapping labels
source("R/common.R")

dds_full <- getDESeq("full", coldata$sample)
res_full <- results(dds_full, name = "condition_E2_vs_Control", alpha=0.05)

df <- as.data.frame(res_full)
df$GENEID <- rownames(df)
symbols <- unique(tx2gene[, c("GENEID", "SYMBOL")])
df <- merge(df, symbols, by = "GENEID") # add gene symbol so we can rank
write.csv(as.data.frame(df), "results/deseq2/res_full.csv")

#png("results/deseq2/pca.png", width = 1600, height = 1200, res = 200)
#plotPCA(vst(dds_full), intgroup = "condition")
#dev.off()

resLFC <- lfcShrink(dds_full, coef = "condition_E2_vs_Control", res=res_full, type = "apeglm")
resLFC <- as.data.frame(resLFC)
resLFC$GENEID <- rownames(resLFC)
symbols <- unique(tx2gene[, c("GENEID", "SYMBOL")])
resLFC <- merge(resLFC, symbols, by = "GENEID") # add gene symbol so we can rank




png("results/deseq2/MAplot.png", width = 1600, height = 1200, res = 200)

# Highlight specific genes of interest by name
genes_to_highlight <- c("GREB1", "PGR", "TFF1")

resLFC$category[(resLFC$padj < 0.05) & resLFC$log2FoldChange > 0]  <- "Up"
resLFC$category[(resLFC$padj < 0.05) & resLFC$log2FoldChange < 0]  <- "Down"
resLFC$category[(is.na(resLFC$padj)) | resLFC$padj > 0.05] <- "Not Significant"

resLFC$category[resLFC$SYMBOL %in% genes_to_highlight] <- "Known ER Targets"

p <- ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color = category)) +
  geom_point(alpha = 0.6) +
  geom_point(data = subset(resLFC, SYMBOL %in% genes_to_highlight),
             size = 3) +
  scale_x_log10() + # MA plots typically use a log scale for expression
  scale_color_manual(values = c("Not Significant" = "grey",
                                "Up" = "red",
                                "Down" = "blue",
                                "Known ER Targets" = "darkred")) +
  theme_minimal() +
  geom_text_repel(data = subset(resLFC, resLFC$SYMBOL %in% genes_to_highlight),
                  aes(label = SYMBOL), color = "black", fontface = "bold") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"))+
  labs(x = "Mean Normalized Counts",
       y = "log2 fold change (E2 vs Control)",
       color=NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
print(p)
dev.off()


sig <- subset(resLFC, padj < 0.05) # select only genes that are significant (p < 0.001)
cols <- c("GENEID", "SYMBOL", "log2FoldChange", "padj", "baseMean")
DE_genes <- sig[order(-sig$log2FoldChange), cols]
write.csv(DE_genes, "results/deseq2/DE_genes.csv", row.names=FALSE)
