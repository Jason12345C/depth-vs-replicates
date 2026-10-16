# unzip the reference transcriptome
hdr <- system("gunzip -c ref/Homo_sapiens.GRCh38.r116.transcriptome.fa.gz | grep '^>'", intern = TRUE) #invokes the OS command specified

tx2gene <- data.frame(
  TXNAME = sub("^>(\\S+).*", "\\1", hdr), # the first word, which kallisto uses as the ID
  GENEID = sub(".*gene:(\\S+).*", "\\1", hdr),
  SYMBOL = ifelse(grepl("gene_symbol:", hdr),
                  sub(".*gene_symbol:(\\S+).*", "\\1", hdr), NA)
)

write.csv(tx2gene, "ref/tx2gene_r116.csv", row.names = FALSE)
