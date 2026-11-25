process RUN_JOB {
    label 'r'

    input:
    path samplesheet
    path tx2gene
    tuple val(draw), val(n)
    path truth
    path negatives
    path quant_dir

    output:
    path "draw*.rds", emit: padj
    path "power_*.csv", emit: table

    script:
    """
    power_analysis.R job ${samplesheet} \\
                         ${tx2gene} \\
                         ${draw} \\
                         ${n}  \\
                         ${truth} \\
                         ${negatives}
    """
}
