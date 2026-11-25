process FASTQC {
    input:
    tuple val(sample), path(input_fastq)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_reports  // glob: both files FastQC writes

    script:
    """
    fastqc -t ${task.cpus} -q --outdir . ${input_fastq}
    """
}
