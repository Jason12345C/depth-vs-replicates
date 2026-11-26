process DE_FULL {
    label 'r'

    input:
    path samplesheet
    path tx2gene
    path quant_dirs                // all 14 full-depth kallisto folders

    output:
    path "res_full.csv", emit: results
    path "DE_genes.csv", emit: de_genes
    path "MAplot.png",   emit: ma_plot

    script:
    """
    differential_analysis.R ${samplesheet} ${tx2gene}
    """
}
