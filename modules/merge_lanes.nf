process MERGE_LANES {
    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("${sample}.fastq.gz")

    script:
    """
    cat ${fastqs} > ${sample}.fastq.gz
    """
}
