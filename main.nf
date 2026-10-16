// MCF7 ± estradiol: more replicates or deeper sequencing?
// QC each lane, merge lanes per sample, subsample to several depths, quantify with kallisto.
// Tuple convention everywhere: (sample, fastq, depth)

params.samplesheet = 'sample_sheets/samplesheet_all.csv'
params.index       = 'ref/ensembl_r116.idx'
params.depths = [2500000, 5000000, 10000000, 15000000, 20000000, 25000000]
params.outdir      = 'results'

// Single-end reads: kallisto needs the fragment length distribution
params.frag_len    = 200
params.frag_sd     = 20

process FASTQC {
    tag "${fq.simpleName}"                       // e.g. SRR1012917
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(fq)

    output:
    path "*_fastqc.{zip,html}"                   // glob: both files FastQC writes

    script:
    """
    fastqc -t ${task.cpus} -q --outdir . ${fq}
    """
}

process MERGE_LANES {
    tag "${sample}"

    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("${sample}.fastq.gz")

    script:
    """
    cat ${fastqs} > ${sample}.fastq.gz
    """
}

process SUBSAMPLE {
    tag "${sample}_${depth}"

    input:
    tuple val(sample), path(fq), val(depth)

    output:
    tuple val(sample), path("${sample}_${depth}.fastq.gz"), val(depth)

    script:
    """
    seqtk sample -2 -s 42 ${fq} ${depth} | gzip > ${sample}_${depth}.fastq.gz
    """
}

process KALLISTO_QUANT {
    tag "${sample}_${depth}"
    publishDir "${params.outdir}/kallisto", mode: 'copy'

    input:
    tuple val(sample), path(fq), val(depth)
    path index

    output:
    tuple val(sample), path("${sample}_${depth}"), val(depth), emit: quant
    path "${sample}_${depth}.log",                           emit: log

    script:
    """
    kallisto quant -i ${index} -o ${sample}_${depth} -t ${task.cpus} \\
        --single -l ${params.frag_len} -s ${params.frag_sd} ${fq} \\
        2> ${sample}_${depth}.log
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports                                 // a list of every report file

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

workflow {
    // One item per run (lane): (sample, fastq)
    runs = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fastq)) }

    // QC every lane individually
    FASTQC(runs)

    // Group lanes by sample, then merge: (sample, merged fastq)
    merged = MERGE_LANES(runs.groupTuple())

    depths     = channel.fromList(params.depths)
    subsampled = SUBSAMPLE(merged.combine(depths))
    full       = merged.map { sample, fq -> tuple(sample, fq, 'full') }

    KALLISTO_QUANT(subsampled.mix(full), file(params.index))

    // One MultiQC report over all FastQC reports + kallisto logs
    MULTIQC( FASTQC.out.mix(KALLISTO_QUANT.out.log).collect() )
}
