process MAKE_SAMPLESHEET {
    input:
    path runtable
    path run_url
    path md5

    output:
    path "samplesheet.csv"

    script:
    """
    make_samplesheet.py --input ${runtable} \\
                        --urls ${run_url} \\
                        --md5 ${md5} \\
                        --out samplesheet.csv

    """
}
