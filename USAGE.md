# mhap_to_methyl Usage Guide

This document walks through building and running the `mhap_to_methyl` converter on large collections of mHap files, with storage-minded output guidance that keeps subsequent dataset-level aggregation straightforward.

## Build

```bash
make
```

The Makefile links against `zlib` and `pthread`. Rebuild after code updates with the same command. Clean the binary with `make clean`.

## Input data

- `--cpg`: bgzip-compressed and tabix-indexed CpG coordinate index (e.g. `hg19_CpG.gz`). Only the first two columns (chromosome, CpG position) are consumed.
- `--mhap-list`: optional text file listing mHap files to process, one per line. Lines beginning with `#` or blank lines are ignored.
- Positional arguments: additional `.mhap.gz` files. When mixing `--mhap-list` and positional filenames, all inputs are processed.

Each mHap file is decoded independently. Records with invalid counts, unknown chromosomes, absent CpGs, or mismatched haplotype lengths are skipped and tallied in the warning summary that prints to `stderr` once the run finishes.

## Output organisation

- `--output <path>`: write a single sample summary to this file. Use when processing exactly one input file.
- `--output-dir <dir>`: write one file per input. Each output is named after the input basename (stripping `.gz`/`.mhap`) plus `--output-suffix` (default `.tsv.gz`).
- `--output-suffix .ext`: optional override for the per-file suffix (e.g. `.tsv` for an uncompressed TSV). Files whose names end in `.gz` are bgzip-compatible.
- `--min-cov C`: filter out CpGs with total coverage below `C` before writing.

Every output contains `chrom`, `pos`, `cov`, `meth`, and `beta` columns. Gzipped TSVs keep per-sample files compact, easy to ship between systems, and trivial to stream into aggregation pipelines or convert to parquet later.

## Memory and threading

- `--threads N` strictly bounds the number of worker threads consuming mHap files. The program never spawns more threads than requested (default 1).
- Each worker allocates two dense 64-bit arrays per chromosome (coverage and methylated counts) tied to the shared CpG coordinate list. Memory usage therefore scales with `#CpGs × 16 bytes × active_threads`, independent of the number of mHap records.

Example: one mHap file to one output

```bash
./mhap_to_methyl \
  --cpg /data/hg19_CpG.gz \
  --output /results/sample123.tsv.gz \
  --threads 10 \
  sample123.mhap.gz
```

Example: thousands of inputs listed in a manifest

```bash
./mhap_to_methyl \
  --cpg /data/hg19_CpG.gz \
  --output-dir /results/studyA \
  --threads 10 \
  --mhap-list studyA_mhap_files.txt
```

Each line of `studyA_mhap_files.txt` points to `sampleXYZ.mhap.gz`; the converter produces `/results/studyA/sampleXYZ.tsv.gz` for downstream aggregation.

## Aggregation strategy guidance

1. Generate per-sample gzipped TSVs using `--output-dir` grouped by dataset.
2. Persist the gzipped TSVs; they compress extremely well and can be streamed (`zcat file.tsv.gz`) when building dataset-level summaries in Python/R/Polars.
3. Convert the TSVs to parquet in a separate lightweight step if needed for analytics platforms—no converter changes required.

Warnings printed at the end flag unusual input issues—capture them in batch logs to spot problematic sources early.
