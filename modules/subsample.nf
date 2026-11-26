process SUBSAMPLE {
    input:
    tuple val(sample), path(fq), val(depth)

    output:
    tuple val(sample), path("${sample}_${depth}.fastq.gz"), val(depth)

    script:
    """
    seqtk sample -2 -s 42 ${fq} ${depth} | gzip > ${sample}_${depth}.fastq.gz
    """
}
