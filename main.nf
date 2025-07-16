#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = null
params.input_glob = null
params.input_filelist = null
params.output_folder = './results'
params.cpus = 8
params.mem_gbs = 64

workflow {

    samples_ch = Channel.empty()
    if (params.input) {
        log.info("Multi-sample mode via CSV: ${params.input}")
        samples_ch = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
    }
    else if (params.input_glob) {
        log.info("Glob-input mode: ${params.input_glob}")
        samples_ch = Channel.fromFilePairs(params.input_glob, flat: true)
            .map { sample, reads -> tuple(sample, reads[0], reads[1]) }
    }
    else if (params.input_filelist) {
        log.info("Filelist mode: ${params.input_filelist}")
        samples_ch = Channel.fromPath(params.input_filelist)
            .splitText()
            .map { f ->
                def sample = file(f).baseName.replaceAll(/_R[12].*$/, '')
                tuple(sample, file(f), null)
            }
    }
    else {
        error("Specify one of --input, --input_glob or --input_filelist")
    }

    samples_ch
        | scrub
        | parseStats
        | collect
        | report
        | view { "Summary report written to: ${it}" }
}

process scrub {
    container "${params.container__scrubber}:${params.scrubber__version}"
    tag "${sample}"
    publishDir "${params.output_folder}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}.stats.txt")

    script:
    """
    zcat ${r1} ${r2 ? " ${r2}" : ""} | scrub.sh \
      -i - \
      -o ${sample}_clean_R1.fastq.gz \
      -p ${params.cpus} \
      2>&1 | tee ${sample}.stats.txt
    """
}

process parseStats {
    tag "${sample}"
    publishDir "${params.output_folder}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(stats_txt)

    output:
    path("*.metrics.csv")

    script:
    """
    total=\$(grep -m1 'total read count' ${stats_txt} | cut -d: -f4 | tr -d ' ')
    removed=\$(grep -m1 'spot(s) masked or removed' ${stats_txt} | awk '{print \$1}')
    total=\${total}
    removed=\${removed}
    remaining=\$(( total - removed ))
    echo "${sample},\${total},\${removed},\${remaining}" > ${sample}.metrics.csv
    """
}

process report {
    publishDir("${params.output_folder}/summary_report.csv", mode: 'copy', overwrite: true)

    input:
    file metrics

    output:
    path "summary_report.csv"

    script:
    """
    echo "sample,total_reads,reads_removed,reads_remaining" > summary_report.csv
    cat ${metrics.join(' ')} >> summary_report.csv
    """
}
