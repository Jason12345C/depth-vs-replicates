process RUN_TRUTH {
    label 'r'
    
    input:
    path samplesheet
    path tx2gene
    path quant_dirs

    output:
    path "truth.txt", emit: truth
    path "negatives.txt", emit: negatives

    script:
    """
    power_analysis.R truth ${samplesheet} ${tx2gene}
    """
}
