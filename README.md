# Nucleotide SNP & Variant Analyzer

A local desktop application for analyzing SNP and indel frequencies across any NCBI nucleotide accession and its linked sequencing samples.

## Overview

This tool:

- Downloads a reference sequence from NCBI using an accession ID
- Finds linked SRA (Sequence Read Archive) sequencing samples
- Aligns reads to the reference using minimap2
- Computes SNP and indel frequencies across the genome
- Displays variant counts by gene and high-frequency mutation summaries
- **Automatically finds alternative accessions with SRA data** if the provided accession has no linked samples

## Who should use it

Any researcher or analyst who wants to compare sample-level SNP and variant information for NCBI nucleotide sequences. It supports any valid NCBI accession with a FASTA and GenBank record.

### What accessions work?

- **Direct accessions**: Accessions with direct SRA sample links (e.g., `NC_045512.2` for SARS-CoV-2)
- **RefSeq sequences**: Reference sequences like HIV-1 (`NC_001802.1`) are automatically mapped to available strain data
- **Generic sequences**: Any organism sequence — the app finds related samples

## Prerequisites

- Python 3.10 or newer
- `fasterq-dump`, `minimap2`, and `samtools` installed and available in your PATH
- Internet access for NCBI and SRA
- Valid email address (for NCBI API access)

## Install required tools

### macOS (Homebrew)
```bash
brew install sratoolkit minimap2 samtools
```

### Linux
Use your package manager or install from project websites.

## Set up the Python environment

From the application folder:

```bash
cd /path/to/SNPApp
python3 -m venv snp_venv
source snp_venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

If `requirements.txt` is missing, install packages directly:

```bash
pip install PyQt5 pysam pandas matplotlib seaborn biopython
```

## Run the app

```bash
cd /path/to/SNPApp
source snp_venv/bin/activate
python desktop_app.py
```

## Using the app

### Setup
1. Enter your email address (required for NCBI API access) and click **Update Email**.
2. This email is used to track API requests; it can be a generic or shared email.

### Finding & Selecting an Accession
1. In the **"Find Accessions with Available Samples"** section:
   - Enter an organism or keyword (e.g., `SARS-CoV-2`, `HIV`, `Influenza`, `Monkeypox`)
   - Click **Search Accessions**
2. The app will:
   - Verify that sequencing samples are available
   - Find complete genome accessions for that organism
   - Display options with the total number of available samples
3. **Select an accession** from the dropdown
   - The selected accession appears in the "Selected Reference Sequence" field

### Analysis workflow
1. Optionally adjust the **"Max Samples"** limit (default: 10)
2. Click **Fetch Sequencing Samples**
   - The app downloads the reference and finds matching SRA samples
3. Select one or more samples from the list
4. Click **Run Analysis**
5. Review the results:
   - Mutation profile plot
   - SNP count by gene
   - High-frequency variant summary (>80% frequency)

## Notes

- **Accession Selection**: Use the "Find Accessions with Available Samples" feature to search for valid sequences. This ensures you only select accessions that have linked SRA samples.
- The app automatically extracts the reference sequence name from the downloaded FASTA, so it works with different accessions.
- **Complete Genomes**: The search prioritizes complete genome sequences; partial sequences are included as fallback options.
- If a GenBank record contains CDS annotations, the SNP counts are assigned to genes; otherwise, they're labeled "Intergenic."
- Temporary files (FASTQ, SAM, BAM) are cleaned up after processing.

## Handling Sample Download Failures

If one or more sequencing samples fail to download:
- The app will **automatically skip** the failed sample and continue analysis with successful ones
- Failed samples are retried up to 3 times with exponential backoff
- A summary message shows which samples succeeded and which failed
- Possible reasons for failure:
  - Sample is from DDBJ (Japan, DRR prefix) or ENA (Europe, ERR prefix) and may have temporary unavailability
  - Network issues or NCBI downtime
  - Sample has been removed or restricted

## Filtering by Sample Source

You can filter which sample sources to include:
- **SRR (NCBI)**: Most reliable ✓ (enabled by default)
- **ERR (ENA Europe)**: Generally reliable ✓ (enabled by default)
- **DRR (DDBJ Japan)**: Less reliable ⚠️ (disabled by default due to frequent availability issues)

To adjust filters:
1. Use the **"Filter Sample Sources"** checkboxes after fetching samples
2. Uncheck DRR to prevent unreliable Japanese samples
3. The sample list updates immediately to show only selected sources

### Why Disable DRR?
DDBJ samples (DRR prefix) frequently fail with "exit status 3" because:
- Metadata is less reliable than NCBI/SRA
- Network connectivity to DDBJ from outside Japan is slower
- Samples may have restricted availability

For more reliable analysis, **disable DRR** unless you specifically need data from Japanese repositories.

## Troubleshooting

- **Search returns no results**: Try a broader search term (e.g., use "SARS" instead of "SARS-CoV-2 strain XYZ").
- **Organism not found**: Verify the organism name spelling. Use common names (e.g., "HIV", "COVID", "Flu").
- **No samples found for selected accession**: This shouldn't happen if you used the search feature. If it does, try a different accession from the search results.
- **Sample download failed (Exit status 3)**: The sample is unavailable. The app will skip it and continue with others. Try re-running or selecting different samples.
- **Command not found**: Install missing command-line tools and ensure they are on PATH.
- **Missing Python packages**: Activate the virtual environment and run `pip install -r requirements.txt`.
- **Network issues**: Confirm internet access and NCBI availability.
- **Email error**: Ensure you've entered a valid email format and clicked **Update Email**.


## File summary

- `desktop_app.py` — Main application code
- `requirements.txt` — Python dependencies
- `packages.txt` — System dependencies for Linux environments
- `README.md` — This file

## Recommended workflow

Use the app locally on your machine. This avoids cloud restrictions and supports the required external bioinformatics tools (fasterq-dump, minimap2, samtools).

## Citation

If you use this tool in your research, please cite the underlying tools:
- BioPython (PMID: 19304878)
- minimap2 (PMID: 29750223)
- SAMtools (PMID: 21320865)
- NCBI SRA/Entrez (NCBI resource)

