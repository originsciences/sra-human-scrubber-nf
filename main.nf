#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input          = null            // new: CSV sample sheet
params.input_glob     = null            // existing
params.input_filelist = null            // existing
params.output_folder  = './results'
params.cpus           = 8
params.mem_gbs        = 64

workflow {

    Channel
        .empty()
        .ifEmpty {
            if ( params.input ) {
                log.info "Launching multi-sample mode, using samplesheet: ${params.input}"
                Channel
                  .fromPath(params.input)
                  .splitCsv(header: true)
                  .map { row ->
                      tuple(row.sample, file(row.fastq_1), file(row.fastq_2))
                  }
            }
            else if ( params.input_glob ) {
                log.info "Launching glob-input mode: ${params.input_glob}"
                Channel.fromFilePairs(params.input_glob, flat: true)
                       .map { sample, reads -> tuple(sample, reads[0], reads[1]) }
            }
            else if ( params.input_filelist ) {
                log.info "Launching filelist mode: ${params.input_filelist}"
                Channel.fromPath(params.input_filelist)
                       .splitText()
                       .map { f -> def sample = file(f).baseName.replaceAll(/_R[12].*$/, '')
                                   tuple(sample, file(f), null) }
            }
            else {
                error "No input specified! Provide --input, --input_glob, or --input_filelist"
            }
        }
        .set{ samples_ch }

    scrub_workflow(samples_ch)
}

workflow scrub_workflow {

    input:
    tuple val(sample), path(r1), path(r2)

    main:
    scrub_ch = Channel.of(tuple(sample, r1, r2))

    process scrub {
        tag "$sample"

        input:
        tuple val(sample), path(r1), path(r2)

        output:
        path "${sample}*.fastq*" into cleaned_ch

        script:
        """
        scrubber \\
          -i ${r1} ${ r2 ? "-i ${r2}" : "" } \\
          -o ${sample}_clean_R1.fastq.gz \\
          -p ${params.cpus}
        """
    }

    cleaned_ch.view { "Finished scrub on: ${it}" }
}
