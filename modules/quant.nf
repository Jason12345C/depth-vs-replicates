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
