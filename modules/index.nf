process KALLISTO_INDEX {
    input:
    path cdna
    path ncrna

    output:
    path "ensembl_r116.idx"

    script:
    """
    cat ${cdna} ${ncrna} > "Homo_sapiens.GRCh38.r116.transcriptome.fa.gz"
    kallisto index -i ensembl_r116.idx -t ${task.cpus} Homo_sapiens.GRCh38.r116.transcriptome.fa.gz
    """
}
