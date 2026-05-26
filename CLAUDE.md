# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Snakemake workflow for RNA-seq analysis: FASTQ → trimmed reads → HISAT2 alignment → featureCounts quantification → count matrix → DESeq2 differential expression.

## How to Run

```bash
# Edit config files first, then:
sh run_snakemake_workflow.sh
# Equivalent to:
nohup snakemake --cores 20 --keep-going --use-conda 1>snak.log 2>&1 &
```

## Configuration

All config in `config/`:

- **`config.yaml`** — points to reference genome, annotation GTF, units table, contrast table, and sets `workdir: result`
- **`units.tsv`** — tab-separated: `group`, `sample`, `fastq1`, `fastq2`. Maps sample IDs to FASTQ filenames and defines experimental groups for DE analysis.
- **`de_contrast.tsv`** — pairwise contrasts for DESeq2, one `groupA\tgroupB` per line.
- **`samples.tsv`** — auto-generated from `units.tsv` by `rule make_samples_file`. Two columns: group, sample.

## Pipeline Stages

1. **`trim_fastq`** — fastp: adapter/quality trimming of raw FASTQ → `01cleandata/`
2. **`fastqc_multiqc`** — FastQC + MultiQC report → `02fastqc/`
3. **`unzip_genome_annotation`** — decompress .gz reference/GTF if needed → `data/ref/`
4. **`build_hisat2_index`** — HISAT2 index from reference genome
5. **`align_hisat2`** — paired-end alignment with `--rna-strandness RF` → SAM
6. **`samtools_sort_idx`** — SAM → sorted BAM + index → `03aligned/`
7. **`quantification`** — featureCounts (gene-level, exon, stranded mode 2) via `workflow/script/run-featurecounts.R` → `04quantification/`
8. **`abundance_estimates_to_matrix`** — merge per-sample counts into gene × sample matrices (counts, FPKM, TPM, TMM, TPM+TMM) via `workflow/script/abundance_to_matrices.R` → `05expression/`
9. **`de_analysis`** — runs `run_DE_analysis.pl` (from Trinity toolkit) with DESeq2, then intelligently fixes missing gene column headers via `workflow/script/add_gene_header.py` → `06de_analysis/DESeq2_Result/`
10. **`generate_report`** — generates `report.html`: a lightweight summary with sample QC, alignment stats, expression summary, and DEG counts (pure Python, no extra deps).
11. **`generate_qmd_report`** — renders `workflow/script/report.qmd` via Quarto into `report_qmd.html`: a rich report with TOC, embedded figures (correlation heatmap, PCA, volcano plot, alignment bar charts), and callout blocks for user customization. Self-contained for offline viewing.

## Key Files

- `workflow/snakefile` — primary workflow definition with all rules inline
- `workflow/rules/*.smk` — older/alternative rule definitions (the main snakefile is authoritative; these contain some duplicate rules with different input patterns)
- `workflow/script/run-featurecounts.R` — wrapper for Rsubread::featureCounts, outputs gene counts + FPKM + TPM per sample
- `workflow/script/abundance_to_matrices.R` — merges per-sample quant files into expression matrices with TMM normalization (via edgeR)
- `workflow/script/abundance_to_matrices.py` — Python port of the same merge logic (not used by the pipeline currently)
- `workflow/script/add_gene_header.py` — intelligently detects and fixes missing gene column headers in DE_results and count_matrix files (compares field counts, not string matching)
- `workflow/script/make_report.py` — generates HTML summary report from pipeline outputs (QC, alignment, expression, DE); pure stdlib, no extra dependencies
- `workflow/script/report.qmd` — Quarto template for rich HTML report with embedded plots (PCA, heatmap, volcano, bar charts) and callout blocks for user customization
- `workflow/script/corr.R` — standalone script for sample correlation heatmaps and PCA plots (expects `data/samples.txt` and pre-computed matrix, not wired into the pipeline)
- `workflow/env/rnaseq_env.yaml` — conda environment spec with bioinformatics tools (HISAT2, samtools, fastp, fastqc, multiqc, R/Bioconductor packages, Trinity)

## Data Directory

- `data/rawdata/` — input FASTQ files
- `data/ref/` — reference genome FASTA and GTF annotation
- `result/` — all pipeline output (working directory set in config.yaml)

## Note: Rule Duplication

Some rules in `workflow/rules/*.smk` overlap with the main `workflow/snakefile`. The snakefile is the authoritative pipeline definition. The rule files use different input patterns (e.g., `filter.smk` expects `{sample}_1.fastq.gz` directly from rawdata, while the snakefile uses the units.tsv lookup).
