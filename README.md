# sra-human-scrubber-nf
Apply the NCBI human read removal tool (HRRT) to mask human reads from metagenomic FASTQ data, enabling safe submission to public repositories like NCBI SRA.

[![Test Workflow](https://github.com/FredHutch/sra-human-scrubber-nf/actions/workflows/test.yaml/badge.svg)](https://github.com/FredHutch/sra-human-scrubber-nf/actions/workflows/test.yaml)

## Multi-sample CSV input

Running on multiple samples using a `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
S1,S1_R1.fastq.gz,S1_R2.fastq.gz
S2,S2_R1.fastq.gz,S2_R2.fastq.gz
```

You can run:

```bash
nextflow run originsciences/sra-human-scrubber-nf \
  --input samplesheet.csv \
  --output_folder scrub_results \
  --cpus 8 \
  --mem_gbs 64
```

The pipeline will process each sample in parallel.

---

## Single-sample modes

Single samples:

1. **Filename patterns**

   ```bash
   nextflow run . \
     --input_glob "*_R{1,2}.fastq.gz" \
     --output_folder scrub_results
   ```

2. **File list**

   ```bash
   nextflow run . \
     --input_filelist files.txt \
     --output_folder scrub_results
   ```

---

## Arguments

* `--input` : `samplesheet.csv` (new multi-sample mode)
* `--input_glob` : Wildcard pattern (e.g., `*_R{1,2}.fastq.gz`)
* `--input_filelist` : Plain file list (one FASTQ per line)
* `--output_folder` : Directory for scrubbed FASTQs
* `--container__scrubber` : Docker image (default from config)
* `--scrubber__version` : Tag for the image
* `--cpus` : CPU count per job (default: 8)
* `--mem_gbs` : Memory (GB) per job (default: 64)

---

## Examples

**Multi-sample (CSV)**:

```bash
nextflow run . \
  --input my_samplesheet.csv \
  --output_folder out_dir \
  --cpus 12 \
  --mem_gbs 96
```

**Glob pattern (paired-end)**:

```bash
nextflow run . \
  --input_glob "*_R{1,2}.fastq.gz"
```

**File list (single-end)**:

```bash
nextflow run . \
  --input_filelist my_samples.txt
```

---

## Details

* Supports **single-end** and **paired-end** reads.
* Utilizes Docker/Singularity container for `ncbi/sra-human-scrubber`.
* Supports per-job resource settings (`cpus`, `mem_gbs`).
* Results are written under `<output_folder>/<sample>/…` with scrubbed FASTQ files and logs.

---

## Background

The NCBI Human Read Removal Tool (HRRT) masks reads with k‑mers from a human-derived reference database. It preserves pathogen sequences while sanitizing human contamination ([nf-co.re][1], [github.com][2], [github.com][3], [ncbiinsights.ncbi.nlm.nih.gov][4], [bactopia.github.io][5]).

---

## Integration Tips

* Sample sheets follow the format used in nf-core pipelines (e.g., RNA-seq) .
* For advanced uses, you can incorporate additional columns (e.g. sample metadata) by extending the `tuple(...)` channel mapping logic in `main.nf`.

---

## License

MIT

---

[1]: https://nf-co.re/rnaseq/2.0/docs/usage "rnaseq: Usage - nf-core"
[2]: https://github.com/ncbi/sra-human-scrubber "ncbi/sra-human-scrubber - GitHub"
[3]: https://github.com/FredHutch/sra-human-scrubber-nf "FredHutch/sra-human-scrubber-nf - GitHub"
[4]: https://ncbiinsights.ncbi.nlm.nih.gov/2023/02/02/scrubbing-human-sequences-sra-submissions "Scrubbing human sequence contamination from ... - NCBI Insights"
[5]: https://bactopia.github.io/v2.2.0/enhancements "Enhancements to OSS - Bactopia"
