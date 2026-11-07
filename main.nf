// MCF7 ± estradiol: more replicates or deeper sequencing?
// QC each lane, merge lanes per sample, subsample to several depths, quantify with kallisto.

params {
    samplesheet: Path
    index: Path
    depths: List<Integer>
    outdir: Path
    frag_len: Integer
    frag_sd: Integer
}

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

process KALLISTO_QUANT {
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

workflow {

    main:
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
    all_reports = FASTQC.out.fastqc_reports.mix(KALLISTO_QUANT.out.log)
    MULTIQC( all_reports.collect() )

    publish:
    fastqc_reports = FASTQC.out.fastqc_reports
    multiqc_data = MULTIQC.out.data
    multiqc_html = MULTIQC.out.html
    kallisto_quant = KALLISTO_QUANT.out.quant
    kallisto_log = KALLISTO_QUANT.out.log
}

output {
    fastqc_reports {
        path 'fastqc'
    }
    multiqc_data {
        path 'multiqc'
    }
    multiqc_html{
        path 'multiqc'
    }
    kallisto_quant {
        path 'kallisto'
    }
    kallisto_log {
        path 'kallisto'
    }
}
