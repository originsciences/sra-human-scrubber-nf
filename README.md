# sra-human-scrubber-nf
Apply the NCBI human read removal tool (HRRT) to mask human reads from metagenomic data

## Background

When submitting genomic sequences to public repositories, it is important to
remove any human sequences which may have been inadvertently included.
This is particularly needed for specimens which are obtained from a human
source, but for which the primary organisms of interest are non-human
(for example, when studying the human microbiome).

This workflow will use the NCBI-approved tool for masking all human sequences
with N's in the raw FASTQ data.
While this can be used to scrub previously-analyzed datasets in preparation
for submission to public repositories (as is required for the Sequence Read
Archive), it could also be used to scrub datasets at the start of a project
prior to running any analyses.

## Usage:

```
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
  ```
