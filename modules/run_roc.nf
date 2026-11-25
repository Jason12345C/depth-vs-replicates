process RUN_ROC {
    label 'r'

    input:
    path truth
    path negatives
    path padj_files           // all 606 padj files, from every job

    output:
    path "roc_runs.csv", emit: roc

    script:
    """
    run_roc.R . ${truth} ${negatives}
    """
}
