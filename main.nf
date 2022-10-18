#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Scrub human sequences from a single FASTQ file
process scrub {
    container "${params.container__scrubber}:${params.scrubber__version}"
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
    memory "${params.mem_gbs}.GB"
    cpus "${params.cpus}"
    
    input:
    path "input/"

    output:
    path "*"

    script:
    template 'scrub.sh'
}

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/sra-human-scrubber-nf <ARGUMENTS>

Required Arguments:

  Input Data:
  --input_glob          Wildcard expression indicating the input files to be processed
    AND/OR
  --input_filelist      Path to file containing the list of files to process, one per line

  Output Location:
  --output_folder       Folder for output files

  Version:
  --container__scrubber The Docker image being used for the SRA Human Scrubber
                        (default: ${params.container__scrubber})
  --scrubber__version   The specific tag (version) of the Docker image being used
                        (default: ${params.scrubber__version})

  Resources:
  --mem_gbs             Amount of memory used per process (in gigabytes)
  --cpus                Number of CPUs used per process
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.output_folder == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }
    

    // Show help message if the user does not specify any inputs
    if ( params.input_glob == false && params.input_filelist == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }
    
    // Construct the inputs from both of the optional input sources
    if ( params.input_glob ){
        input_glob_ch = Channel.fromPath(
            "${params.input_glob}",
            checkIfExists: true,
            glob: true
        )
    } else {
        input_glob_ch = Channel.empty()
    }

    if ( params.input_filelist ){
        input_filelist_ch = Channel.fromPath(
            "${params.input_filelist}",
            checkIfExists: true,
            glob: false
        )
        .splitText()
        .map {
            it -> file(
                // Remove the newline character from each line
                it.substring(0, it.length() - 1),
                checkIfExists: true
            )
        }
    } else {
        input_filelist_ch = Channel.empty()
    }

    // Analyze data from both sources
    scrub(
        input_filelist_ch
        .mix(input_glob_ch)
    )

}