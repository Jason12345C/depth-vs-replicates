process TX2GENE {
    label 'r'
    
    input:
    path cdna
    path ncrna

    output:
    path "tx2gene.csv"

    script:
    """
    tx2gene.R ${cdna} ${ncrna}
    """
}
