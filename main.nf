// MCF7 ± estradiol: more replicates or deeper sequencing?
// Subsample each sample to several depths, quantify with kallisto.

params.samplesheet = 'samplesheet.csv'
params.index       = 'ref/ensembl_r116.idx'
params.depths      = [2500000, 5000000, 10000000]   // reads per sample
params.outdir      = 'results'

// Single-end only (see the hint in chat). Delete these if the data is paired-end.
params.frag_len    = 200
params.frag_sd     = 20

process SUBSAMPLE {
    tag "${sample}_${depth}"

    input:
    tuple val(sample), path(fq), val(depth)

    output:
    tuple val(sample), val(depth), path("${sample}_${depth}.fastq.gz")

    script:
    """
    seqtk sample -2 -s 42 ${fq} ${depth} | gzip > ${sample}_${depth}.fastq.gz
    """
}

process KALLISTO_QUANT {
    tag "${sample}_${depth}"
    publishDir "${params.outdir}/kallisto", mode: 'copy'

    input:
    tuple val(sample), val(depth), path(fq)
    path index

    output:
    tuple val(sample), val(depth), path("${sample}_${depth}")

    script:
    """
    kallisto quant -i ${index} -o ${sample}_${depth} -t ${task.cpus} \\
        --single -l ${params.frag_len} -s ${params.frag_sd} ${fq}
    """
}

// TODO 1: process FASTQC  (raw reads only)
// TODO 2: process MULTIQC (runs once, over everything)

workflow {
    // One item per sample: (sample, fastq)
    reads = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fastq)) }

    depths = channel.fromList(params.depths)

    // Every sample x every depth
    subsampled = SUBSAMPLE(reads.combine(depths))

    // Full-depth reads skip subsampling
    full = reads.map { sample, fq -> tuple(sample, 'full', fq) }

    KALLISTO_QUANT(subsampled.mix(full), file(params.index))
}
