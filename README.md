# SARS-CoV-2 SNP Frequency Analyzer

A desktop application for identifying and interpreting SNPs and indels across SARS-CoV-2 sequencing samples from NCBI SRA, using the Wuhan-Hu-1 reference genome (NC_045512.2) as the baseline.

## Purpose

This tool parses FASTQ sequencing data to calculate mutation frequencies across the SARS-CoV-2 genome and interprets their biological significance — identifying whether each mutation is synonymous (silent) or non-synonymous (amino acid change), and which viral protein it affects. Comparing multiple samples reveals which mutations are conserved across strains and which are unique to individual samples.

## What it does

1. Fetches SARS-CoV-2 sequencing runs from NCBI SRA linked to NC_045512.2
2. Streams reads directly through alignment — nothing written to disk
3. Calls SNP and indel frequencies using pysam pileup against the Wuhan-Hu-1 reference
4. Annotates each variant with its gene (S, N, ORF1ab, etc.) and amino acid effect
5. Displays results as a frequency plot and three tables

## Why SARS-CoV-2

- **Largest public dataset** — thousands of SRA runs from patients, time points, and variants of concern (Alpha, Delta, Omicron)
- **Small genome (29,903 bp)** — fast runtimes; 3–5 samples run in minutes
- **Well-annotated genes** — every mutation maps to a known protein with understood function
- **Known mutations with documented significance** — spike protein mutations like D614G, N501Y, and E484K are in the literature and directly relatable to immune evasion and vaccine escape
- **Ideal for comparison** — selecting samples from different variants of concern shows how the virus evolved across the pandemic

## Features

- **Streaming pipeline** — `fastq-dump --stdout | minimap2 | samtools sort` runs concurrently with no intermediate FASTQ files on disk
- **Platform auto-detection** — reads the sequencing platform (Illumina, Nanopore, PacBio) per run and applies the correct minimap2 preset automatically
- **Max reads cap** (default 500,000) — limits download per sample; reduces typical runtime to ~2 minutes
- **Configurable thresholds** — minimum read depth and minimum allele frequency
- **Synonymous / non-synonymous classification** — for every SNP inside a coding sequence, computes the codon change and amino acid effect on both plus- and minus-strand genes
- **Multi-sample** — select multiple SRA runs; results are overlaid on the same plot and combined in the tables

## Prerequisites

| Tool | Purpose |
|------|---------|
| Python 3.8 via conda | Runtime (see note below) |
| `minimap2` | Read alignment |
| `samtools` | BAM sorting |
| `fastq-dump` (SRA Toolkit) | Read streaming from NCBI |

> **macOS note:** The app uses PyQt5, which requires a Python build linked against the system Qt libraries. The conda base environment (Python 3.8) works. Python 3.10+ from a standard venv breaks the Qt plugin on macOS and will show a "cocoa platform plugin" error at launch.

### Install bioinformatics tools (conda)

```bash
conda install -c bioconda minimap2 samtools sra-tools
```

### Install Python dependencies

```bash
pip install PyQt5 pysam pandas matplotlib seaborn biopython
```

## Running the app

```bash
# Make sure the conda environment is active (not a venv)
deactivate 2>/dev/null; true
python desktop_app.py
```

## Workflow

### 1. Enter your NCBI email
NCBI requires an email for Entrez API access. Enter it in the Email field and click **Update Email**. Any valid email works.

### 2. Fetch sequencing runs
- Set the **Max results** limit (how many SRA runs to list, default 10)
- Click **Fetch Runs**
- Each checkbox shows the run ID, title, date, total bases, and sequencing platform
- Runs without size data are dimmed — they are often inaccessible
- For best results, prefer **SRR-prefixed** runs (NCBI) over ERR (ENA)

### 3. Configure analysis settings

| Setting | Default | Notes |
|---------|---------|-------|
| Max reads/sample | 500,000 | Set to 0 for all reads (warns if run >500 Mbp) |
| Min depth | 20 | Minimum reads at a position to call a variant |
| Min frequency | 5% | Minimum allele frequency to report a variant |

### 4. Select samples — more is better

The app is designed for **multiple samples**. Selecting runs from different variants of concern (e.g. Alpha, Delta, Omicron) makes every part of the output more informative:

- The **plot** overlays each sample with a different marker so you can see which mutations are shared vs unique
- The **Fixed Mutations table** `Samples` column shows how many runs carry each mutation — a mutation in 5/5 samples is a reliable strain marker; one in 1/5 may be noise
- **Minority variants** appearing at the same position across independent samples are strong evidence they are real rather than sequencing artefacts

**Recommended:** select 3–5 samples from different patients or variant backgrounds.

### 5. Run analysis
Click **Run Analysis**. The status bar shows live progress and estimated time remaining. Click **Abort** to cancel at any time.

### 6. Interpret results

#### SNP / INDEL Profile plot
- X-axis: genome position (bp) across the 29,903 bp SARS-CoV-2 genome
- Y-axis: allele frequency (proportion of reads carrying that base)
- Red = SNPs, blue = INDELs; marker shape varies per sample
- Dashed line at 0.8 = 80% threshold — variants above this are effectively fixed
- Gene track below shows only SARS-CoV-2 genes containing called variants (S, N, M, E, ORF1ab, etc.)

#### SNP / INDEL Count by Gene
Overview of which viral genes carry the most mutations across all selected samples. High counts in the spike (S) gene are typical of variants of concern.

#### Fixed Mutations (>80%)

Mutations present in more than 80% of reads — effectively fixed in the sample and likely defining characteristics of that variant relative to Wuhan-Hu-1.

| Column | Meaning |
|--------|---------|
| Position | Genome position (1-based) |
| Gene | Viral gene (S, N, ORF1ab, etc.) or Intergenic |
| Base | Alternate nucleotide |
| Variant_Type | SNP or INDEL |
| Mutation_Effect | Biological effect (see table below) |
| AA_Change | Amino acid change e.g. D614G, N501Y; "—" for non-coding |
| Avg_Frequency | Mean frequency across all samples carrying this mutation |
| Samples | Number of selected samples in which this mutation was found |

`Avg_Frequency` is averaged across samples. `Samples` tells you how consistent the mutation is — high frequency across many samples indicates a conserved variant-defining change.

#### Minority Variants (20%–80%)

Mutations present in 20–80% of reads. These represent within-sample diversity — mixed infection, intra-host evolution, or quasi-species dynamics. Variants below 20% are excluded from display (indistinguishable from sequencing error at typical coverage) but are included in the CSV export.

| Column | Meaning |
|--------|---------|
| Sample_ID | The SRR run this variant was found in |
| Position | Genome position (1-based) |
| Gene | Viral gene or Intergenic |
| Variant_Type | SNP or INDEL |
| Base | Alternate nucleotide |
| Mutation_Effect | Biological effect |
| AA_Change | Amino acid change e.g. E484K; "—" for non-coding |
| Count | Reads supporting this base |
| Depth | Total reads at this position |
| Frequency | Count / Depth — raw per-sample value, not averaged |

## Biological interpretation

| Mutation_Effect | Meaning |
|-----------------|---------|
| **Synonymous** | Codon changed but amino acid unchanged — silent, no effect on protein sequence |
| **Non-synonymous** | Amino acid changed — may affect protein structure, antibody binding, or drug susceptibility |
| **Stop gained** | Premature stop codon — likely truncates the viral protein |
| **Stop lost** | Stop codon mutated — protein may extend beyond its normal terminus |
| **Start lost** | Start codon disrupted — protein may not be translated |
| **Frameshift (INDEL)** | Reading frame shifted — disrupts all downstream codons |
| **Intergenic** | Outside any annotated gene — may affect regulatory elements |

### Notable SARS-CoV-2 mutations to look for

| Mutation | Gene | Effect | Associated with |
|----------|------|--------|----------------|
| D614G | S | Non-synonymous | Increased transmissibility, dominant globally from mid-2020 |
| N501Y | S | Non-synonymous | Increased ACE2 binding; Alpha, Beta, Gamma, Omicron |
| E484K | S | Non-synonymous | Immune evasion, reduced antibody neutralisation; Beta, Gamma |
| K417N | S | Non-synonymous | Immune evasion; Beta, Omicron |
| P681H/R | S | Non-synonymous | Furin cleavage site; increased cell entry; Alpha, Delta |
| del69-70 | S | Frameshift | Spike deletion; Alpha, Omicron |

## Known limitations

- **ERR (ENA) runs** are less reliable than SRR (NCBI) — `fastq-dump` can silently return 0 reads for some ENA accessions; prefer runs that show a file size
- **Amplicon sequencing** concentrates reads in specific regions and may produce 0 variants elsewhere — whole-genome sequencing runs give the most complete profile
- **Single sample** — the app works with one sample but comparison features only become meaningful with two or more
- **500 K read cap** is sufficient for reliable frequency estimates on the 30 kb genome; raising it increases runtime without significantly improving results for most samples

## File summary

| File | Purpose |
|------|---------|
| `desktop_app.py` | Entire application |
| `requirements.txt` | Python package list |
| `packages.txt` | System-level dependencies (Linux) |
| `install_and_run.sh` | Installs dependencies and launches the app |
| `README.md` | This file |

> A multi-organism version supporting 16 viral and bacterial presets is preserved on the `multi-organism` branch.

## Dependencies credited

- [BioPython](https://biopython.org/) — sequence I/O, Entrez API, codon translation (PMID: 19304878)
- [minimap2](https://github.com/lh3/minimap2) — read alignment (PMID: 29750223)
- [SAMtools](http://www.htslib.org/) — BAM processing (PMID: 21320865)
- [pysam](https://pysam.readthedocs.io/) — pileup-based variant calling
- [NCBI SRA / Entrez](https://www.ncbi.nlm.nih.gov/sra) — sequencing data source
