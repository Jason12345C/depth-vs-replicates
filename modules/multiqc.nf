process MULTIQC {
    input:
    path reports                                 // a list of every report file

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """
}
