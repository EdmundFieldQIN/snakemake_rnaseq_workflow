# RNA-seq Analysis Pipeline

A Snakemake workflow for RNA-seq data analysis, from raw FASTQ files to differential expression results and a comprehensive HTML report.

## Pipeline Overview

```
FASTQ → fastp (trimming) → FastQC/MultiQC → HISAT2 (alignment)
→ samtools (sort/index) → featureCounts (quantification)
→ expression matrix → DESeq2 (DE analysis) → HTML report
```

## Requirements

- [Conda](https://docs.conda.io/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Snakemake](https://snakemake.readthedocs.io/) (installed automatically via conda environment)

## Quick Start

### 1. Prepare your data

Place your input files in the appropriate directories:

```
data/
├── rawdata/          # FASTQ files (paired-end: *_R1.fq.gz, *_R2.fq.gz)
└── ref/              # Reference genome (FASTA) and annotation (GTF)
```

### 2. Configure the pipeline

Edit the configuration files in `config/`:

**`config.yaml`** — Main configuration:
```yaml
workdir: result                           # Output directory
units: ../config/units.tsv                # Sample table
contrast: ../config/de_contrast.tsv       # Contrast table
PE: true                                  # Paired-end mode

ref:
  genome: ../data/ref/genome.fa           # Reference genome (supports .gz)
  annotation: ../data/ref/genes.gtf       # Gene annotation (supports .gz)
```

**`units.tsv`** — Sample-to-FASTQ mapping (tab-separated):
```
group   sample   fastq1           fastq2
WT      WT-1     WT-1_R1.fq.gz    WT-1_R2.fq.gz
WT      WT-2     WT-2_R1.fq.gz    WT-2_R2.fq.gz
Treat   T-1      T-1_R1.fq.gz     T-1_R2.fq.gz
Treat   T-2      T-2_R1.fq.gz     T-2_R2.fq.gz
```

| Column | Description |
|--------|-------------|
| `group` | Experimental group for DE analysis |
| `sample` | Unique sample identifier |
| `fastq1` | Read 1 FASTQ filename |
| `fastq2` | Read 2 FASTQ filename |

**`de_contrast.tsv`** — Pairwise comparisons (tab-separated, one per line): Use Group name
```
Treatment1      Control
Treatment2      Control
```

### 3. Run the pipeline

```bash
sh run_snakemake_workflow.sh
```

Or manually:
```bash
nohup snakemake --cores 20 --keep-going --use-conda 1>snak.log 2>&1 &
```

Monitor progress:
```bash
tail -f snak.log
```

## Output Structure

```
result/
├── 01cleandata/                     # Trimmed FASTQ + fastp reports
├── 02fastqc/                        # FastQC + MultiQC report
├── 03aligned/                       # Sorted BAM + index files
├── 04quantification/                # Per-sample gene counts
├── 05expression/                    # Merged expression matrices
│   ├── genes.counts.matrix          # Raw counts
│   ├── genes.FPKM.EXPR.matrix       # FPKM
│   ├── genes.TPM.EXPR.matrix        # TPM
│   ├── genes.TMM.EXPR.matrix        # TMM-normalized counts
│   └── genes.TPM.TMM.EXPR.matrix    # TMM-normalized TPM
├── 06de_analysis/
│   └── DESeq2_Result/               # DE results + volcano/MA plots
├── report.html                      # Lightweight HTML summary
├── report_qmd.html                  # Rich Quarto report with figures
└── logs/                            # Rule execution logs
```

## Pipeline Stages

| # | Rule | Tool | Description |
|---|------|------|-------------|
| 1 | `trim_fastq` | fastp | Adapter/quality trimming, auto-detection |
| 2 | `fastqc_multiqc` | FastQC + MultiQC | Per-sample QC + aggregated report |
| 3 | `unzip_genome_annotation` | gzip | Decompress .gz reference/GTF if needed |
| 4 | `build_hisat2_index` | HISAT2 | Build genome index |
| 5 | `align_hisat2` | HISAT2 | Paired-end alignment (`--rna-strandness RF`) |
| 6 | `samtools_sort_idx` | samtools | SAM → sorted BAM + index |
| 7 | `quantification` | featureCounts | Gene-level counts, exon features, stranded mode 2 |
| 8 | `abundance_estimates_to_matrix` | edgeR | Merge per-sample counts + TMM normalization |
| 9 | `de_analysis` | DESeq2 (Trinity) | Differential expression for all contrasts |
| 10 | `generate_report` | Python | Lightweight HTML summary (no extra deps) |
| 11 | `generate_qmd_report` | Quarto | Rich HTML report with plots, TOC, and customization |

## Reports

Two HTML reports are generated automatically:

### Lightweight Report (`report.html`)
- Sample QC summary from fastp (read counts, Q20/Q30, GC, duplication)
- HISAT2 alignment statistics per sample
- Expression matrix dimensions and per-sample gene counts
- DEG counts per contrast
- **Zero extra dependencies** — pure Python stdlib

### Rich Report (`report_qmd.html`)
- **Table of Contents** with left sidebar navigation
- **Embedded figures:** QC bar charts, alignment composition, correlation heatmap, PCA, volcano plots, DEG bar charts
- **Tabbed volcano plots** — one tab per contrast
- **Customizable callout blocks** — each section has a template for adding personal notes
- **Software version table** — all tool and R package versions
- **Self-contained** — single HTML file, viewable offline
- **Template editable** — modify `workflow/script/report.qmd` to add custom content

To customize the Quarto report, edit `workflow/script/report.qmd` and look for `<!-- USER: ... -->` comments. Re-render with:
```bash
quarto render workflow/script/report.qmd \
  -P workdir:result -P samples_file:config/samples.tsv
```

## Configuration Reference

### `config.yaml`

| Key | Description | Example |
|-----|-------------|---------|
| `workdir` | Output directory | `result` |
| `units` | Path to units table | `../config/units.tsv` |
| `contrast` | Path to contrast table | `../config/de_contrast.tsv` |
| `PE` | Paired-end mode | `true` |
| `ref.genome` | Reference genome FASTA (.gz supported) | `../data/ref/genome.fa` |
| `ref.annotation` | Gene annotation GTF (.gz supported) | `../data/ref/genes.gtf` |

### `units.tsv`

Tab-separated file mapping samples to FASTQ files. The `group` column defines experimental groups for differential expression analysis. Each group should have at least 2 biological replicates for DESeq2.

### `de_contrast.tsv`

Tab-separated file with one contrast per line: `groupA    groupB`. DESeq2 compares groupA vs groupB (groupB is the reference/baseline).

## Advanced Usage

### Running specific rules

```bash
# Run only trimming
snakemake --cores 20 trim_fastq

# Run only DE analysis (skipping upstream if outputs exist)
snakemake --cores 20 de_analysis

# Force re-run a specific rule
snakemake --cores 20 --force de_analysis
```

### Dry run (check what will execute)

```bash
snakemake --cores 20 --dryrun
```

### Visualize the DAG

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### Using a different conda environment

```bash
snakemake --cores 20 --use-conda --conda-frontend mamba
```

## Key Scripts

| Script | Description |
|--------|-------------|
| `workflow/script/run-featurecounts.R` | featureCounts wrapper, outputs gene counts + FPKM + TPM |
| `workflow/script/abundance_to_matrices.R` | Merges per-sample counts into matrices with TMM normalization |
| `workflow/script/add_gene_header.py` | Fixes missing gene column headers in DE output files |
| `workflow/script/make_report.py` | Generates lightweight HTML summary report |
| `workflow/script/report.qmd` | Quarto template for rich HTML report |

## Troubleshooting

**Protected BAM file error:** If Snakemake refuses to overwrite protected BAM files, remove them first:
```bash
rm result/03aligned/*.bam result/03aligned/*.bam.bai
```

**Missing conda environment:** The first run with `--use-conda` will create the environment. This may take 10-20 minutes. Subsequent runs use the cached environment.

**Out of memory during alignment:** Reduce the thread count for `align_hisat2` in the snakefile, or reduce `--cores`.

**Quarto not found:** Ensure the conda environment is properly created. Quarto is included in the environment spec.
