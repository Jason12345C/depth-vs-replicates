process DOWNLOAD_FASTQ {
    maxForks 4
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), val(run), val(url), val(md5)

    output:
    tuple val(sample), path("${run}.fastq.gz")

    script:
    """
    aria2c -x 4 -s 4 --file-allocation=none \\
           --max-tries=5 --retry-wait=10 \\
           -o ${run}.fastq.gz ${url}

    echo "${md5}  ${run}.fastq.gz" | md5sum -c -
    """
}
