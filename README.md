# SARS-CoV-2 SNP Frequency Analyzer

A bioinformatics tool that downloads real SARS-CoV-2 sequencing data from NCBI, aligns it to the Wuhan-Hu-1 reference genome, and calculates how frequently each mutation appears — identifying which viral proteins are affected and whether the mutation changes the protein sequence.

---

## Background

**SNPs (Single Nucleotide Polymorphisms)** are positions in the genome where a single DNA base differs from the reference. For example, at position 23,403 in SARS-CoV-2, an A→G change produces the D614G mutation in the spike protein — a change that made the virus significantly more transmissible.

**Indels** are insertions or deletions of one or more bases. Even a single base deletion shifts the entire reading frame downstream, typically disrupting the protein completely.

**Allele frequency** is the proportion of sequencing reads that carry a particular base at a position. A frequency of 0.95 means 95% of reads show that base — the mutation is effectively fixed in the sample. A frequency of 0.30 means only 30% of reads carry it, suggesting a mixed population or within-host variation.

---

## What it does

1. Fetches SARS-CoV-2 sequencing runs from NCBI SRA (linked to NC_045512.2)
2. Streams reads through alignment with no intermediate files written to disk
3. Calculates SNP and indel frequencies at every position using pysam pileup
4. Annotates each variant with its gene and whether it is **synonymous** (amino acid unchanged — silent) or **non-synonymous** (amino acid changed — potentially affects protein function)
5. Displays results as a frequency plot and three tables

**The tool is designed for multiple samples.** Selecting runs from different variants of concern (Alpha, Delta, Omicron) lets you see which mutations are conserved across strains and which are unique to one sample — this is where the comparison becomes biologically meaningful.

---

## Prerequisites

| Tool | Install via |
|------|------------|
| Python 3.8 (conda) | `conda` — do not use a standard venv on macOS |
| `minimap2` | `conda install -c bioconda minimap2` |
| `samtools` | `conda install -c bioconda samtools` |
| `fastq-dump` | `conda install -c bioconda sra-tools` |

```bash
# Install all bioinformatics tools at once
conda install -c bioconda minimap2 samtools sra-tools

# Install Python dependencies
pip install PyQt5 pysam pandas matplotlib seaborn biopython
```

## Running

```bash
python desktop_app.py
```

---

## Workflow

1. **Enter your NCBI email** — required for Entrez API access
2. **Fetch Runs** — retrieves sequencing runs linked to the SARS-CoV-2 reference; prefer SRR-prefixed runs (NCBI) over ERR (ENA)
3. **Select 3–5 samples** from different patients or variant backgrounds for meaningful comparison
4. **Run Analysis** — streams, aligns, and calls variants; takes ~2 minutes per sample

---

## Results

### SNP / INDEL Profile Plot
Shows where mutations fall across the 29,903 bp genome and how frequently they appear. Each sample gets a different marker shape so shared vs unique mutations are immediately visible. The gene track below highlights only genes with called variants. The dashed line at 80% separates fixed mutations from minority variants.

### Fixed Mutations (>80%)
Mutations present in more than 80% of reads — effectively fixed in that sample relative to Wuhan-Hu-1. These are the defining characteristics of a strain. The `Samples` column shows how many of your selected runs carry each mutation; a mutation seen in all samples is a reliable strain marker.

### Minority Variants (20–80%)
Mutations present in 20–80% of reads, indicating a mixed population — possibly co-infection, within-host evolution, or a variant emerging under selection. Variants below 20% are excluded from display as they are indistinguishable from sequencing error at typical coverage.

### SNP / INDEL Count by Gene
A summary of how many mutations fall in each viral gene across all samples. High counts in the spike (S) gene are characteristic of variants of concern.

---

## Mutation effects explained

| Effect | Meaning |
|--------|---------|
| **Synonymous** | DNA changed but amino acid unchanged — silent, no effect on protein |
| **Non-synonymous** | Amino acid changed — may affect protein folding, antibody binding, or drug susceptibility |
| **Stop gained** | Premature stop codon introduced — truncates the protein |
| **Frameshift (INDEL)** | Reading frame shifted — disrupts all codons downstream |
| **Intergenic** | Falls outside any annotated gene |

---

## Notable SARS-CoV-2 mutations

| Mutation | Gene | Associated with |
|----------|------|----------------|
| D614G | Spike | Increased transmissibility — dominant globally from mid-2020 |
| N501Y | Spike | Increased ACE2 binding — Alpha, Beta, Gamma, Omicron |
| E484K | Spike | Reduced antibody neutralisation — Beta, Gamma |
| K417N | Spike | Immune evasion — Beta, Omicron |
| P681H/R | Spike | Enhanced cell entry via furin cleavage site — Alpha, Delta |

If your analysis calls any of these, the `AA_Change` column will show the exact notation (e.g. `D614G`) and `Mutation_Effect` will confirm it as non-synonymous.

---

## Known limitations

- ERR (ENA) runs can silently return 0 reads — prefer SRR runs with a listed file size
- Amplicon sequencing data concentrates reads in specific regions and may miss variants elsewhere; whole-genome sequencing runs give the most complete profile
- The 500K read cap is sufficient for the 30kb SARS-CoV-2 genome; raising it increases runtime without meaningfully improving results

---

## Dependencies

- [BioPython](https://biopython.org/) — sequence I/O, Entrez API, codon translation (PMID: 19304878)
- [minimap2](https://github.com/lh3/minimap2) — read alignment (PMID: 29750223)
- [SAMtools](http://www.htslib.org/) — BAM processing (PMID: 21320865)
- [pysam](https://pysam.readthedocs.io/) — pileup-based variant calling
- [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) — sequencing data source

> A multi-organism version supporting 16 viral and bacterial presets is available on the `multi-organism` branch.
