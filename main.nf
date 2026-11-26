// MCF7 ± estradiol: more replicates or deeper sequencing?
// QC each lane, merge lanes per sample, subsample to several depths, quantify with kallisto.
include { DOWNLOAD_FASTQ } from './modules/download_fastq.nf'
include { FASTQC } from './modules/fastqc.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { KALLISTO_INDEX } from './modules/index.nf'
include { KALLISTO_QUANT } from './modules/quant.nf'
include { SUBSAMPLE } from './modules/subsample.nf'
include { MERGE_LANES } from './modules/merge_lanes.nf'
include { MAKE_SAMPLESHEET } from './modules/make_samplesheet.nf'
include { DE_FULL } from './modules/de_full.nf'
include { TX2GENE } from './modules/tx2gene.nf'
include {RUN_TRUTH} from './modules/run_truth.nf'
include {RUN_JOB} from './modules/run_job.nf'
include {RUN_ROC} from './modules/run_roc.nf'

params {
    samplesheet: Path
    sra_run_table: Path
    urls: Path
    md5: Path
    index: Path
    depths: List<Integer>
    outdir: Path
    frag_len: Integer
    frag_sd: Integer
    cdna: Path
    ncrna: Path 
    n_draws: Integer
    replicates: List<Integer>
}

workflow {
    main:
    // One item per run (lane): (sample, fastq)
    samplesheet = MAKE_SAMPLESHEET(
    file(params.sra_run_table, checkIfExists: true),
    file(params.urls, checkIfExists: true),
    file(params.md5, checkIfExists: true)
)
    rows = samplesheet.splitCsv(header: true)

    split = rows.branch { r ->
        local:  file("fastq/${r.run_accession}.fastq.gz").exists()
        remote: true                                  // everything else
    }

    local_fastqs = split.local
        .map { r -> tuple(r.sample, file("fastq/${r.run_accession}.fastq.gz")) }

    downloaded = DOWNLOAD_FASTQ(
        split.remote.map { r -> tuple(r.sample, r.run_accession, r.url, r.fastq_md5) })

    fastqs = local_fastqs.mix(downloaded)            // (sample, fastq) either way

    FASTQC(fastqs)                            // QC the files that actually exist

    // Group lanes by sample, then merge: (sample, merged fastq)
    merged = MERGE_LANES(fastqs.groupTuple())

    depths     = channel.fromList(params.depths)
    subsampled = SUBSAMPLE(merged.combine(depths))
    full       = merged.map { sample, fq -> tuple(sample, fq, 'full') }

    KALLISTO_INDEX(file(params.cdna), file(params.ncrna))
    KALLISTO_QUANT(subsampled.mix(full), KALLISTO_INDEX.out)

    // One MultiQC report over all FastQC reports + kallisto logs
    all_reports = FASTQC.out.fastqc_reports.mix(KALLISTO_QUANT.out.log)
    MULTIQC( all_reports.collect() )

    // Analysis
    TX2GENE(file(params.cdna), file(params.ncrna))
    full_quants = KALLISTO_QUANT.out.quant
    .filter  { sample, dir, depth -> depth == 'full' }   // keep only full depth
    .map     { sample, dir, depth -> dir }               // keep just the folder
    .collect()                                           // all 14 as one list
    DE_FULL(samplesheet, TX2GENE.out, full_quants)

    quants_25M = KALLISTO_QUANT.out.quant
    .filter { sample, dir, depth -> depth == 25000000}
    .map    { sample, dir, depth -> dir}
    .collect()

    sub_quants = KALLISTO_QUANT.out.quant
    .filter { sample, dir, depth -> depth.toString() != 'full' }   // all six subsampled depths
    .map    { sample, dir, depth -> dir }                          // just the folder
    .collect()      
                                                   // all 84 as one list
    RUN_TRUTH(samplesheet, TX2GENE.out, quants_25M)

    draws = channel.of(1..params.n_draws)                // 1, 2, ..., 20
    reps  = channel.fromList(params.replicates)          // 2, 3, 4, 5, 6, 7
    jobs = draws
    .combine(reps)                                   // every draw with every n: 120 pairs
    .filter { draw, n -> !(n == 7 && draw > 1) }     // keep n = 7 only for draw 1: 101 pairs

    RUN_JOB(samplesheet, TX2GENE.out, jobs,
        RUN_TRUTH.out.truth, RUN_TRUTH.out.negatives, sub_quants)    

    power_runs = RUN_JOB.out.table
    .collectFile(name: 'power_runs.csv', keepHeader: true, sort: { it.name })
    
    RUN_ROC(RUN_TRUTH.out.truth, RUN_TRUTH.out.negatives, RUN_JOB.out.padj.collect())

    publish:
    sample_sheet = MAKE_SAMPLESHEET.out
    fastqc_reports = FASTQC.out.fastqc_reports
    multiqc_data = MULTIQC.out.data
    multiqc_html = MULTIQC.out.html
    kallisto_quant = KALLISTO_QUANT.out.quant
    kallisto_log = KALLISTO_QUANT.out.log
    index = KALLISTO_INDEX.out
    de_results = DE_FULL.out.results
    de_ma = DE_FULL.out.ma_plot
    de_genes = DE_FULL.out.de_genes
    truth = RUN_TRUTH.out.truth
    negatives =  RUN_TRUTH.out.negatives
    roc = RUN_ROC.out.roc
    power_runs = RUN_JOB.out.table
    }

output {
    sample_sheet {
        path 'samplesheet'
    }
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
    index {
        path 'ref'
    }
    de_results {
        path 'analysis/deseq'
    }
    de_ma {
        path 'analysis/deseq'
    }
    de_genes {
        path 'analysis/deseq'
    }
    truth {
        path 'analysis'
    }
    negatives {
        path 'analysis'
    }
    power_runs {
        path 'analysis'
    }
    roc {
        path 'analysis'
    }



}
