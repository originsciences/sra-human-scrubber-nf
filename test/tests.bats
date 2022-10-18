#!/usr/bin/env bats

@test "Minimal test" {
    nextflow run \
        ../main.nf \
        --input_glob scrubber_test.fastq \
        --output_folder output \
        -profile docker

    cmp output/scrubber_test.fastq scrubber_expected_output.fastq
}

