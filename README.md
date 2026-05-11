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

### 5. Run analysis
Select one or more runs and click **Run Analysis**. The status bar shows live progress and estimated time remaining. Click **Abort** to cancel at any time.

### 6. Interpret results

**SNP / INDEL Profile plot**
- X-axis: genome position (bp)
- Y-axis: allele frequency (proportion of reads carrying that base)
- Red circles = SNPs, blue squares = INDELs (shape varies by sample)
- Dashed line at 0.8 = 80% threshold; variants above this are likely fixed in the sample
- Gene track below the plot shows only genes that contain at least one called variant

**SNP / INDEL Count by Gene table**
Counts of SNPs and INDELs per gene, sorted by frequency. Gives a quick view of which genes are most variable.

**High-Frequency Mutations (>80%) table**
Variants present in >80% of reads — likely fixed differences from the reference. Includes `Mutation_Effect` and `AA_Change` for biological interpretation.

**All Variants table**
Every variant passing the depth and frequency thresholds, with:

| Column | Meaning |
|--------|---------|
| Position | Genome position (1-based) |
| Gene | Annotated gene name, or "Intergenic" |
| Variant_Type | SNP or INDEL |
| Base | Alternate base (or "Indel") |
| Mutation_Effect | Synonymous, Non-synonymous, Stop gained, Frameshift (INDEL), Intergenic, Reference allele |
| AA_Change | Amino acid change in single-letter code (e.g. K47R); "—" for non-coding |
| Count | Number of reads supporting this base |
| Depth | Total reads at this position |
| Frequency | Count / Depth |

## Biological interpretation

- **Synonymous** mutations change the DNA codon but not the amino acid — they are usually neutral and do not affect protein function
- **Non-synonymous** mutations change the amino acid — they may affect protein structure, function, drug binding, or immune evasion
- **Stop gained** mutations introduce a premature stop codon, likely truncating the protein
- **Frameshift (INDEL)** mutations shift the reading frame, typically disrupting all downstream codons
- Variants in **Intergenic** regions may affect regulatory elements but do not directly alter protein sequence

## Tips for specific genome types

| Genome type | Recommended settings |
|-------------|---------------------|
| Viral (~10–30 kb) | Min depth 20, Max reads 500 K, Min freq 5% |
| Small bacterial (~2 Mb, e.g. S. pneumoniae) | Min depth 10, Max reads 1–2 M, Min freq 5% |
| Large bacterial (~4–5 Mb, e.g. E. coli) | Min depth 5, Max reads 2 M, Min freq 5% |

## Known limitations

- **ERR (ENA) runs** are less reliable than SRR (NCBI) — `fastq-dump` can silently return 0 reads for some ENA accessions; prefer runs that show a file size
- **Monkeypox (NC_063383.1)** has no genome-linked SRA runs in NCBI; organism search works but prefer SRR-prefixed results
- **Amplicon sequencing** (e.g. Influenza SRR3165632) concentrates reads in a few regions, which may produce 0 variants in the rest of the genome
- Minus-strand gene annotation requires the reference FASTA to be the same sequence used for alignment; custom accessions work as long as they are in NCBI Nucleotide

## File summary

| File | Purpose |
|------|---------|
| `desktop_app.py` | Entire application |
| `requirements.txt` | Python package list |
| `packages.txt` | System-level dependencies (Linux) |
| `README.md` | This file |

## Dependencies credited

- [BioPython](https://biopython.org/) — sequence I/O, Entrez API, codon translation (PMID: 19304878)
- [minimap2](https://github.com/lh3/minimap2) — read alignment (PMID: 29750223)
- [SAMtools](http://www.htslib.org/) — BAM processing (PMID: 21320865)
- [pysam](https://pysam.readthedocs.io/) — pileup-based variant calling
- [NCBI SRA / Entrez](https://www.ncbi.nlm.nih.gov/sra) — sequencing data source
