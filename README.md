# SNP Frequency Analyzer - Desktop Application

A professional desktop application for analyzing SNP (Single Nucleotide Polymorphism) frequencies across genomic samples.

## Features

✅ **Professional GUI** - Built with PyQt5 for a native desktop experience  
✅ **Real-time Progress Tracking** - Visual progress bar with elapsed time counter  
✅ **Threading** - Analysis runs in background without freezing the UI  
✅ **Interactive Results** - View mutation profiles, SNP counts by gene, and high-frequency mutations  
✅ **Responsive UI** - Scroll through SRR codes with checkboxes for easy selection  
✅ **Error Handling** - Graceful error messages for network or processing issues  

## Installation

### 1. Prerequisites
Ensure you have Homebrew-installed bioinformatics tools:
```bash
brew install sratoolkit minimap2 samtools
```

### 2. Install Python Dependencies
Navigate to the SNPApp folder and install dependencies:
```bash
cd /Users/benedictmugar/Documents/SNPApp
pip install -r requirements.txt
```

Or install manually:
```bash
pip install PyQt5 pandas matplotlib seaborn biopython pysam
```

## Running the Application

### From Terminal
```bash
cd /Users/benedictmugar/Documents/SNPApp
python desktop_app.py
```

### From VS Code
- Open `desktop_app.py`
- Press `Ctrl+F5` (or `Cmd+F5` on Mac) to run

## How to Use

1. **Enter Nucleotide ID** - Type the NCBI nucleotide accession (e.g., `NC_045512.2` for SARS-CoV-2)
2. **Set Result Limit** - Choose how many SRR codes to fetch (1-50)
3. **Fetch SRR Codes** - Click to download available sequencing datasets
4. **Select Samples** - Check the boxes for samples you want to compare
5. **Run Analysis** - Click to start processing
   - Watch the progress bar and elapsed time
   - Analysis runs in background (UI stays responsive)
6. **View Results**
   - **Mutation Profile Plot** - Scatter plot showing SNP frequencies across the genome
   - **SNP Count by Gene** - Table showing total SNPs per gene
   - **High Frequency Mutations** - Top 20 mutations present in >80% of reads

## Advantages Over Streamlit Cloud

| Feature | Streamlit Cloud | Desktop App |
|---------|-----------------|-------------|
| External Tools | ❌ Not supported | ✅ Full support |
| Download/Upload Files | ⚠️ Limited | ✅ Full support |
| UI Responsiveness | ⚠️ Reloads on interact | ✅ True threading |
| Installation | Cloud-based | Local control |
| Deployment | Easy but limited | Requires local setup |

## Troubleshooting

**"Command not found: fasterq-dump"**
- Install SRA toolkit: `brew install sratoolkit`

**"ModuleNotFoundError: No module named 'PyQt5'"**
- Install PyQt5: `pip install PyQt5`

**"Entrez error: connection timeout"**
- Network issue accessing NCBI. Try again or use a different nucleotide ID.

**Progress bar stuck**
- Analysis is still running. Large datasets can take 10+ minutes. Check terminal for subprocess output.

## System Requirements

- Python 3.10+
- macOS (or Linux/Windows with appropriate bioinformatics tools)
- ~2GB RAM minimum
- Internet connection (for NCBI access and SRA downloads)

## Files

- `desktop_app.py` - Main PyQt5 application
- `requirements.txt` - Python package dependencies
- `packages.txt` - System-level dependencies

## Notes

- Temporary files (FASTQ, SAM, BAM) are automatically cleaned up after analysis
- Reference FASTA files are saved as `ref.fasta` in the working directory
- Large SRR datasets (>5GB) may take significant time to download and process
