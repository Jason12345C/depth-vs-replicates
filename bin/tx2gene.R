#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: tx2gene.R <fasta.gz> [more fasta.gz ...]")
cmd <- paste("gunzip -c", paste(shQuote(args), collapse = " "), "| grep '^>'")
hdr <- system(cmd, intern = TRUE)

tx2gene <- data.frame(
  TXNAME = sub("^>(\\S+).*", "\\1", hdr), # the first word, which kallisto uses as the ID
  GENEID = sub(".*gene:(\\S+).*", "\\1", hdr),
  SYMBOL = ifelse(grepl("gene_symbol:", hdr),
                  sub(".*gene_symbol:(\\S+).*", "\\1", hdr), NA)
)

write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)
