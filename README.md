# SNP Frequency Analyzer

A desktop application for identifying and interpreting SNPs and indels across viral and bacterial genomes using publicly available sequencing data from NCBI SRA.

## What it does

1. Selects a reference genome (from 16 presets or a manual accession)
2. Fetches linked sequencing runs from NCBI SRA
3. Streams reads directly through alignment — nothing written to disk
4. Calls SNP and indel frequencies using pysam pileup
5. Annotates each variant with its gene and whether it is synonymous (silent) or non-synonymous (amino acid change)
6. Displays results as an interactive plot and three tables

## Features

- **16 preset organisms** (SARS-CoV-2, HIV-1, Ebola, E. coli, M. tuberculosis, S. pneumoniae D39, and more) plus manual accession entry
- **Platform auto-detection** — reads the sequencing platform (Illumina, Nanopore, PacBio) per run and applies the correct minimap2 preset automatically
- **Streaming pipeline** — `fastq-dump --stdout | minimap2 | samtools sort` runs concurrently with no intermediate FASTQ files on disk
- **Max reads cap** (default 500,000) — limits download per sample; reduces typical runtime to ~2 minutes for viral genomes
- **Configurable thresholds** — minimum read depth and minimum allele frequency before a variant is reported
- **Gene map** — fetches CDS annotations from GenBank; handles RefSeq CON records (e.g. all bacterial NC_ accessions) by falling back to the primary submission automatically
- **Synonymous / non-synonymous classification** — for every SNP inside a coding sequence, computes the codon change and amino acid effect; handles both plus- and minus-strand genes
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

### 2. Select a reference genome
- Choose from the **Preset Organisms** list (radio button), or
- Switch to **Manual Accession** and type any NCBI nucleotide accession (e.g. `NC_045512.2`)

### 3. Fetch sequencing runs
- Set the **Max runs** limit (how many SRA runs to list)
- Click **Fetch Runs**
- Each checkbox shows the run ID, title, date, total bases, and detected sequencing platform
- Runs without size data are dimmed — they are often inaccessible

### 4. Configure analysis settings

| Setting | Default | Notes |
|---------|---------|-------|
| Max reads/sample | 500,000 | Set to 0 for all reads (warns if run >500 Mbp) |
| Min depth | 20 | Lower to 5–10 for large bacterial genomes |
| Min frequency | 5% | Minimum allele frequency to report a variant |

### 5. Select samples — more is better

The app is designed to be run with **multiple samples selected**. Each table and the plot become more informative with more than one run:

- The plot overlays each sample with a different marker shape so you can see which mutations are shared vs unique to one sample
- The `Samples` column in the Fixed Mutations table shows how many runs carry each mutation — a mutation present in 5/5 samples is a reliable strain characteristic; one present in 1/5 may be noise or a unique subvariant
- Minority variants appearing at the same position across multiple independent samples are strong evidence they are real and not sequencing artefacts

For viral genomes, **3–5 samples** from different patients or time points is the recommended starting point. For bacterial genomes, **2–3 samples** given the longer runtimes.

### 6. Run analysis
Click **Run Analysis**. The status bar shows live progress and estimated time remaining. Click **Abort** to cancel at any time.

### 7. Interpret results

#### SNP / INDEL Profile plot
- X-axis: genome position (bp)
- Y-axis: allele frequency (proportion of reads carrying that base)
- Red = SNPs, blue = INDELs; marker shape varies by sample when multiple are selected
- Dashed line at 0.8 = 80% threshold
- Gene track below shows only genes that contain at least one called variant

#### SNP / INDEL Count by Gene
Overview of where mutations concentrate. Counts of SNPs and INDELs per gene across all selected samples, sorted by count. Useful for identifying genes under selection pressure.

#### Fixed Mutations (>80%)

Mutations present in more than 80% of reads — effectively fixed in the sample and likely defining characteristics of that strain relative to the reference.

| Column | Meaning |
|--------|---------|
| Position | Genome position (1-based) |
| Gene | Annotated gene name, or "Intergenic" |
| Base | Alternate base |
| Variant_Type | SNP or INDEL |
| Mutation_Effect | Biological effect (see below) |
| AA_Change | Amino acid change, e.g. K47R; "—" for non-coding |
| Avg_Frequency | Mean frequency across all samples that carry this mutation |
| Samples | Number of selected samples in which this mutation was found |

`Avg_Frequency` is the average across all samples that carried the mutation. `Samples` tells you how consistent it is — a mutation seen in all samples at high frequency is a reliable strain marker.

#### Minority Variants (20%–80%)

Mutations present in 20–80% of reads. These represent genuine sub-population variation — mixed infection, within-host evolution, or heterozygosity. Variants below 20% are excluded from display as they are ambiguous (indistinguishable from sequencing error without very high coverage) but are included in the CSV export.

| Column | Meaning |
|--------|---------|
| Sample_ID | The SRR run this variant was found in |
| Position | Genome position (1-based) |
| Gene | Annotated gene name, or "Intergenic" |
| Variant_Type | SNP or INDEL |
| Base | Alternate base |
| Mutation_Effect | Biological effect (see below) |
| AA_Change | Amino acid change, e.g. K47R; "—" for non-coding |
| Count | Number of reads supporting this base |
| Depth | Total reads at this position |
| Frequency | Count / Depth (raw per-sample value, not averaged) |

`Frequency` here is a raw per-sample value — it tells you exactly what fraction of reads in that one run carried this base. Unlike the Fixed Mutations table it is not averaged across samples.

## Biological interpretation

| Mutation_Effect | Meaning |
|-----------------|---------|
| **Synonymous** | Codon changed but amino acid unchanged — usually neutral, no effect on protein |
| **Non-synonymous** | Amino acid changed — may affect protein structure, function, drug binding, or immune evasion |
| **Stop gained** | Premature stop codon introduced — likely truncates the protein |
| **Stop lost** | Stop codon mutated — protein may be extended beyond its normal end |
| **Start lost** | Start codon disrupted — protein may not be translated |
| **Frameshift (INDEL)** | Reading frame shifted — typically disrupts all downstream codons |
| **Intergenic** | Outside any annotated gene — may affect regulatory elements |

## Tips for specific genome types

| Genome type | Recommended settings |
|-------------|---------------------|
| Viral (~10–30 kb) | Min depth 20, Max reads 500 K, Min freq 5% |
| Small bacterial (~2 Mb, e.g. S. pneumoniae) | Min depth 10, Max reads 1–2 M, Min freq 5% |
| Large bacterial (~4–5 Mb, e.g. E. coli) | Min depth 5, Max reads 2 M, Min freq 5% |

## Known limitations

- **ERR (ENA) runs** are less reliable than SRR (NCBI) — `fastq-dump` can silently return 0 reads for some ENA accessions; prefer runs that show a file size
- **Monkeypox (NC_063383.1)** has no genome-linked SRA runs in NCBI; organism search works but prefer SRR-prefixed results
- **Amplicon sequencing** (e.g. Influenza SRR3165632) concentrates reads in a few regions, which may produce 0 variants elsewhere in the genome
- **Single sample** — the app works with one sample but the comparison features (Samples column, per-sample plot markers) only become meaningful with two or more
- Minus-strand gene annotation requires the reference FASTA to match the sequence used for alignment; custom accessions work as long as they are in NCBI Nucleotide

## File summary

| File | Purpose |
|------|---------|
| `desktop_app.py` | Entire application |
| `requirements.txt` | Python package list |
| `packages.txt` | System-level dependencies (Linux) |
| `install_and_run.sh` | Installs dependencies and launches the app |
| `README.md` | This file |

## Dependencies credited

- [BioPython](https://biopython.org/) — sequence I/O, Entrez API, codon translation (PMID: 19304878)
- [minimap2](https://github.com/lh3/minimap2) — read alignment (PMID: 29750223)
- [SAMtools](http://www.htslib.org/) — BAM processing (PMID: 21320865)
- [pysam](https://pysam.readthedocs.io/) — pileup-based variant calling
- [NCBI SRA / Entrez](https://www.ncbi.nlm.nih.gov/sra) — sequencing data source
