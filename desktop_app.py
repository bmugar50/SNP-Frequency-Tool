import sys
import os
import time

# Fix "Could not find Qt platform plugin cocoa" on macOS when running from a venv.
# Must be set before any PyQt5 import.
def _fix_qt_plugin_path():
    try:
        import PyQt5 as _qt
        _base = os.path.dirname(_qt.__file__)
        for _candidate in (
            os.path.join(_base, "Qt5", "plugins", "platforms"),
            os.path.join(_base, "Qt",  "plugins", "platforms"),
        ):
            if os.path.isdir(_candidate):
                os.environ.setdefault("QT_QPA_PLATFORM_PLUGIN_PATH", _candidate)
                break
    except Exception:
        pass

_fix_qt_plugin_path()
import pysam
import pandas as pd
import matplotlib.patches as mpatches
import seaborn as sns
import re
import xml.etree.ElementTree as ET
import multiprocessing
import bisect
from collections import OrderedDict
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
import subprocess

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLineEdit, QLabel, QTableWidget, QTableWidgetItem,
    QProgressBar, QCheckBox, QScrollArea, QFrame, QGroupBox,
    QMessageBox, QSpinBox, QRadioButton, QButtonGroup,
    QComboBox, QSplitter, QSizePolicy, QHeaderView, QFileDialog, QDoubleSpinBox
)
from PyQt5.QtCore import Qt, pyqtSignal, QThread, QTimer
from PyQt5.QtGui import QFont, QColor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

STYLESHEET = """
QMainWindow, QWidget {
    background-color: #f5f6fa;
    color: #2c3e50;
    font-family: -apple-system, 'Segoe UI', Arial, sans-serif;
    font-size: 13px;
}
QGroupBox {
    background-color: #ffffff;
    border: 1px solid #dde1e7;
    border-radius: 6px;
    margin-top: 10px;
    padding: 8px;
    font-weight: 600;
    color: #2c3e50;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
}
QPushButton {
    background-color: #3498db;
    color: white;
    border: none;
    border-radius: 4px;
    padding: 6px 14px;
    font-weight: 500;
}
QPushButton:hover  { background-color: #2980b9; }
QPushButton:pressed { background-color: #2471a3; }
QPushButton:disabled { background-color: #bdc3c7; color: #ecf0f1; }
QPushButton#run_btn {
    background-color: #27ae60;
    font-size: 13px;
    font-weight: bold;
    padding: 8px 20px;
}
QPushButton#run_btn:hover    { background-color: #229954; }
QPushButton#run_btn:disabled { background-color: #a9cbb8; color: #ecf0f1; }
QPushButton#abort_btn {
    background-color: #e74c3c;
    font-size: 13px;
    font-weight: bold;
    padding: 8px 20px;
    color: white;
}
QPushButton#abort_btn:hover    { background-color: #c0392b; }
QPushButton#abort_btn:disabled { background-color: #e8a9a3; color: #f9ebea; }
QLineEdit, QComboBox, QSpinBox {
    background-color: #ffffff;
    border: 1px solid #ced4da;
    border-radius: 4px;
    padding: 5px 8px;
    color: #2c3e50;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus { border-color: #3498db; }
QProgressBar {
    background-color: #ecf0f1;
    border: 1px solid #dde1e7;
    border-radius: 4px;
    text-align: center;
    min-height: 18px;
    color: #2c3e50;
}
QProgressBar::chunk {
    background-color: #3498db;
    border-radius: 3px;
}
QTableWidget {
    background-color: #ffffff;
    alternate-background-color: #f4f6f8;
    border: 1px solid #dde1e7;
    gridline-color: #ecf0f1;
}
QTableWidget QHeaderView::section {
    background-color: #2c3e50;
    color: #ffffff;
    padding: 6px;
    border: none;
    font-weight: 600;
    font-size: 12px;
}
QScrollArea {
    border: 1px solid #dde1e7;
    border-radius: 4px;
    background-color: #ffffff;
}
QCheckBox { spacing: 6px; padding: 3px 6px; }
QCheckBox:hover { background-color: #ebf5fb; border-radius: 3px; }
QRadioButton { spacing: 6px; }
QSplitter::handle { background-color: #dde1e7; height: 2px; }
"""


# ---------------------------------------------------------------------------
# Preset organisms: display label -> NCBI accession
# ---------------------------------------------------------------------------
PRESET_ORGANISMS = OrderedDict([
    ("-- Select Organism --", ""),
    # ── Viruses ──────────────────────────────────────────────────────────
    ("SARS-CoV-2 (COVID-19)          NC_045512.2", "NC_045512.2"),
    ("Influenza A H1N1 2009 (HA seg) NC_026431.1", "NC_026431.1"),
    ("HIV-1                          NC_001802.1", "NC_001802.1"),
    ("Ebola virus (Zaire)            NC_002549.1", "NC_002549.1"),
    ("Dengue virus 1                 NC_001477.1", "NC_001477.1"),
    ("Hepatitis B virus              NC_003977.2", "NC_003977.2"),
    ("Hepatitis C virus              NC_004102.1", "NC_004102.1"),
    ("Measles virus                  NC_001498.1", "NC_001498.1"),
    ("Monkeypox virus                NC_063383.1", "NC_063383.1"),
    ("Norovirus GI                   NC_001959.1", "NC_001959.1"),
    # ── Bacteria ─────────────────────────────────────────────────────────
    ("E. coli K-12 MG1655            NC_000913.3", "NC_000913.3"),
    ("M. tuberculosis H37Rv          NC_000962.3", "NC_000962.3"),
    ("S. aureus MRSA252              NC_002952.2", "NC_002952.2"),
    ("Salmonella Typhimurium LT2     NC_003197.2", "NC_003197.2"),
    ("K. pneumoniae NTUH-K2044       NC_012731.1", "NC_012731.1"),
    ("S. pneumoniae D39               NC_008533.1", "NC_008533.1"),
])

# SRR prefixes known to be reliably accessible from NCBI
ACCESSIBLE_PREFIXES = ("SRR", "ERR")

# Fallback organism-name search terms — tried only when accession[Genome] returns nothing.
# Verified 2025-05 against NCBI SRA. NC_026431.1 (Influenza HA segment) intentionally
# omitted — organism search returns nothing for individual segments; use curated list.
ORGANISM_SEARCH_TERMS = {
    "NC_001802.1": "Human immunodeficiency virus 1[Organism]",
    "NC_002549.1": "Zaire ebolavirus[Organism]",
    "NC_001477.1": "Dengue virus 1[Organism]",
    "NC_004102.1": "Hepatitis C virus[Organism]",
    "NC_001498.1": "Measles morbillivirus[Organism]",
    "NC_063383.1": "Monkeypox virus[Organism]",
    "NC_001959.1": "Norovirus[Organism]",
    "NC_002952.2": "Staphylococcus aureus[Organism]",
    "NC_003197.2": "Salmonella enterica[Organism]",
    "NC_012731.1": "Klebsiella pneumoniae[Organism]",
    "NC_008533.1": "Streptococcus pneumoniae[Organism]",
}

# Pre-verified SRR/ERR runs — used instantly (no network call) if both searches fail.
# Influenza H1N1 always uses this list since the segment accession can't be found
# via organism search. All other entries are last-resort fallbacks.
CURATED_RUNS = {
    "NC_026431.1": [
        ("SRR1048819", "Influenza H1N1 2009 pandemic"),
        ("SRR3165632", "Influenza H1N1 amplicon"),
        ("SRR1562345", "Influenza H1N1 HA sequencing"),
    ],
}


# Maps SRA platform strings to minimap2 -x presets and short display labels.
PLATFORM_PRESETS = {
    "ILLUMINA":       "sr",
    "ION_TORRENT":    "sr",
    "LS454":          "sr",
    "ABI_SOLID":      "sr",
    "BGISEQ":         "sr",
    "DNBSEQ":         "sr",
    "OXFORD_NANOPORE":"map-ont",
    "PACBIO_SMRT":    "map-pb",
}
PLATFORM_LABELS = {
    "ILLUMINA":       "Illumina",
    "ION_TORRENT":    "IonTorrent",
    "LS454":          "454",
    "OXFORD_NANOPORE":"Nanopore",
    "PACBIO_SMRT":    "PacBio",
    "BGISEQ":         "BGI",
    "DNBSEQ":         "DNB",
}
DEFAULT_PLATFORM = "ILLUMINA"


def _parse_platform(exp_xml: str) -> str:
    """Extract the SRA platform string from an ExpXml summary blob."""
    m = re.search(r'<Platform[^>]*>\s*([A-Z_0-9]+)\s*</Platform>', exp_xml, re.I)
    if m:
        return m.group(1).upper()
    # Fallback: instrument_model attribute sometimes encodes the brand
    m2 = re.search(r'instrument_model="([^"]+)"', exp_xml, re.I)
    if m2:
        model = m2.group(1).upper()
        if "NANOPORE" in model or "MINION" in model or "PROMETHION" in model:
            return "OXFORD_NANOPORE"
        if "PACBIO" in model or "SEQUEL" in model or "SMRT" in model:
            return "PACBIO_SMRT"
    return DEFAULT_PLATFORM


# ---------------------------------------------------------------------------
# Worker: fetch SRR codes
# ---------------------------------------------------------------------------
class SRRFetchWorker(QThread):
    finished = pyqtSignal(pd.DataFrame)
    error_occurred = pyqtSignal(str)
    status_updated = pyqtSignal(str)

    def __init__(self, accession_id, limit, email):
        super().__init__()
        self.accession_id = accession_id
        self.limit = limit
        self.email = email

    def run(self):
        Entrez.email = self.email
        try:
            result = self._fetch_sra_runs(self.accession_id, limit=self.limit)
            self.finished.emit(result)
        except Exception as e:
            self.error_occurred.emit(str(e))

    def _fetch_sra_runs(self, accession_id, limit=10, retries=3):
        # 1. Try primary search: accession linked as reference genome
        df = self._search_ncbi(f"{accession_id}[Genome]", limit, retries)

        # 2. If empty, find an organism-name search term
        if df.empty:
            if accession_id in ORGANISM_SEARCH_TERMS:
                # Preset organism — use pre-verified term (no extra API call)
                org_term = ORGANISM_SEARCH_TERMS[accession_id]
            else:
                # Manual accession — derive organism name from GenBank record
                org_term = self._organism_term_from_genbank(accession_id)

            if org_term:
                self.status_updated.emit("No genome-linked runs — trying organism search...")
                df = self._search_ncbi(org_term, limit, retries)

        # 3. If still empty, use curated list instantly (no network call)
        if df.empty and accession_id in CURATED_RUNS:
            self.status_updated.emit("Using curated known-good runs...")
            rows = [
                {"SRR_ID": run_id, "Title": desc, "Date": "curated",
                 "Bases": "N/A", "Platform": DEFAULT_PLATFORM}
                for run_id, desc in CURATED_RUNS[accession_id]
            ]
            df = pd.DataFrame(rows)

        return df

    def _organism_term_from_genbank(self, accession_id):
        """Fetch organism name from GenBank — only called for unknown/manual accessions."""
        try:
            self.status_updated.emit("Fetching organism name from GenBank...")
            with Entrez.efetch(db="nucleotide", id=accession_id,
                               rettype="gb", retmode="text") as h:
                record = SeqIO.read(h, "genbank")
            organism = record.annotations.get("organism", "")
            if organism:
                return f"{organism}[Organism]"
        except Exception:
            pass
        return None

    def _search_ncbi(self, term, limit, retries):
        for attempt in range(retries):
            try:
                self.status_updated.emit(f"Searching SRA: {term[:60]}...")
                handle = Entrez.esearch(db="sra", term=term, retmax=limit * 3)
                uids = Entrez.read(handle)["IdList"]
                handle.close()

                if not uids:
                    return pd.DataFrame()

                self.status_updated.emit(f"Fetching summaries for {len(uids)} SRA records...")
                sum_handle = Entrez.esummary(db="sra", id=",".join(uids))
                summaries = Entrez.read(sum_handle)
                sum_handle.close()

                rows = []
                for s in summaries:
                    try:
                        platform = _parse_platform(s.get("ExpXml", ""))
                        root = ET.fromstring(f"<root>{s['Runs']}</root>")
                        for run_el in root.findall(".//Run"):
                            acc = run_el.get("acc", "")
                            if not acc.startswith(ACCESSIBLE_PREFIXES):
                                continue
                            rows.append({
                                "SRR_ID":   acc,
                                "Title":    s.get("Title", "No Title")[:60],
                                "Date":     s.get("CreateDate", "Unknown"),
                                "Bases":    run_el.get("total_bases", "N/A"),
                                "Platform": platform,
                            })
                    except Exception:
                        continue
                    if len(rows) >= limit:
                        break

                return pd.DataFrame(rows[:limit])

            except Exception as e:
                if attempt < retries - 1:
                    time.sleep(2)
                else:
                    raise e

        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Worker: full SNP analysis pipeline
# ---------------------------------------------------------------------------
class AnalysisWorker(QThread):
    progress_updated = pyqtSignal(int)
    status_updated = pyqtSignal(str)
    time_updated = pyqtSignal(str, str)   # elapsed, eta
    finished = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)
    sample_error = pyqtSignal(str)        # non-fatal per-sample errors
    aborted = pyqtSignal()

    def __init__(self, accession_id, selected_srr, email, max_reads=0,
                 platform_map=None, min_depth=20, min_freq=0.05):
        super().__init__()
        self.accession_id = accession_id
        self.selected_srr = selected_srr        # list of SRR ID strings
        self.email = email
        self.max_reads = max_reads
        self.platform_map = platform_map or {}  # srr_id -> platform string
        self.min_depth = min_depth
        self.min_freq = min_freq
        self._start_time = None
        self._aborted = False
        self._current_proc = None          # currently running subprocess
        self._sample_times = []            # seconds each completed sample took

    def abort(self):
        self._aborted = True
        if self._current_proc and self._current_proc.poll() is None:
            self._current_proc.kill()

    # ── helpers ─────────────────────────────────────────────────────────

    def _emit_progress(self, pct):
        self.progress_updated.emit(pct)
        elapsed = time.time() - self._start_time

        if self._sample_times:
            avg = sum(self._sample_times) / len(self._sample_times)
            remaining_samples = len(self.selected_srr) - len(self._sample_times)
            eta_str = self._fmt_seconds(avg * remaining_samples) if remaining_samples > 0 else "almost done"
        elif self._setup_elapsed > 0 and pct > 13:
            # Linear interpolation over the sample phase (pct 13 → 95)
            samples_elapsed = max(elapsed - self._setup_elapsed, 0.1)
            frac = max(min((pct - 13) / 82, 0.99), 0.01)
            remaining = max((samples_elapsed / frac) * (1 - frac), 0)
            eta_str = self._fmt_seconds(remaining)
        else:
            eta_str = "calculating..."

        self.time_updated.emit(self._fmt_seconds(elapsed), eta_str)

    @staticmethod
    def _fmt_seconds(secs):
        if secs < 60:
            return f"{secs:.0f}s"
        m, s = divmod(int(secs), 60)
        return f"{m}m {s:02d}s"

    @staticmethod
    def _check_tools():
        """Return list of required CLI tools that are not on PATH."""
        missing = []
        for tool in ("minimap2", "samtools", "fastq-dump"):
            r = subprocess.run(f"which {tool}", shell=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if r.returncode != 0:
                missing.append(tool)
        return missing

    def _run_proc(self, cmd):
        """Run a shell command; capture stderr so failures produce useful messages."""
        self._current_proc = subprocess.Popen(
            cmd, shell=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
        )
        self._current_proc.wait()
        rc = self._current_proc.returncode
        stderr = (self._current_proc.stderr.read()
                  .decode(errors="replace").strip()) if rc != 0 else ""
        self._current_proc = None
        if self._aborted:
            raise InterruptedError("Aborted")
        if rc != 0:
            short_cmd = cmd.split("|")[0].strip()[:80]
            tail = stderr[-400:] if len(stderr) > 400 else stderr
            raise RuntimeError(
                f"Command failed (exit {rc}):\n  {short_cmd}"
                + (f"\n\nTool output:\n{tail}" if tail else "")
            )

    # ── main pipeline ────────────────────────────────────────────────────

    def run(self):
        Entrez.email = self.email
        self._start_time = time.time()
        self._setup_elapsed = 0
        self._sample_errors = []

        try:
            # Pre-flight: check all required tools are on PATH
            missing = self._check_tools()
            if missing:
                self.error_occurred.emit(
                    f"Required tools not found on PATH:\n\n  {', '.join(missing)}\n\n"
                    "Install via conda:\n"
                    "  conda install -c bioconda minimap2 samtools sra-tools"
                )
                return

            n = len(self.selected_srr)

            # 1. Download reference FASTA
            self.status_updated.emit(f"[1/4] Downloading reference sequence ({self.accession_id})...")
            with Entrez.efetch(db="nucleotide", id=self.accession_id,
                               rettype="fasta", retmode="text") as h:
                with open("ref.fasta", "w") as f:
                    f.write(h.read())
            if self._aborted:
                raise InterruptedError("Aborted")
            self._emit_progress(4)

            # 2. Index reference for minimap2
            self.status_updated.emit("[2/4] Building minimap2 index for reference...")
            self._run_proc("minimap2 -d ref.mmi ref.fasta")
            self._emit_progress(8)

            # 2b. Read reference sequence for codon annotation
            ref_record = SeqIO.read("ref.fasta", "fasta")
            ref_seq = str(ref_record.seq)

            # 3. Build gene map from GenBank record
            self.status_updated.emit("[3/4] Fetching gene annotations from GenBank...")
            gene_map, cds_map = self._build_gene_map(self.accession_id)
            n_genes = len(gene_map)
            self.status_updated.emit(f"[3/4] Found {n_genes} annotated gene{'s' if n_genes != 1 else ''}.")
            if self._aborted:
                raise InterruptedError("Aborted")
            self._emit_progress(13)
            self._setup_elapsed = time.time() - self._start_time

            # Build a sorted interval index for O(log n) gene lookup per variant.
            # Intervals sorted by start position; bisect finds the candidate in one step.
            _gene_intervals = sorted(
                [(s, e, g) for g, (s, e) in gene_map.items()], key=lambda x: x[0]
            )
            _gene_starts = [iv[0] for iv in _gene_intervals]

            def _lookup_gene(pos):
                i = bisect.bisect_right(_gene_starts, pos) - 1
                if i >= 0 and _gene_intervals[i][0] <= pos <= _gene_intervals[i][1]:
                    return _gene_intervals[i][2]
                return "Intergenic"

            # 4. Process each SRR sample
            self.status_updated.emit(f"[4/4] Processing {n} sample{'s' if n > 1 else ''}...")
            master_list = []
            per_sample = 82 / n   # progress points allocated per sample

            for i, srr in enumerate(self.selected_srr):
                if self._aborted:
                    raise InterruptedError("Aborted")

                base = 13 + int(i * per_sample)
                sample_start = time.time()
                self._emit_progress(base + 1)  # kick ETA update before download starts

                bam = f"{srr}_sorted.bam"  # pre-compute so finally can always clean up
                try:
                    self._download_and_align(srr, i + 1, n, base, per_sample)
                    self.status_updated.emit(f"[4/4] Sample {i+1}/{n} — calling SNP variants...")
                    df, diag = self._call_snps(bam, self.accession_id,
                                               min_depth=self.min_depth,
                                               min_freq=self.min_freq)
                    df["Sample_ID"] = srr
                    if df.empty:
                        mapped = diag["mapped_reads"]
                        if mapped == 0:
                            reason = ("0 reads mapped to the reference — "
                                      "the sample may be from a different organism or strain")
                        elif diag["positions_above_depth"] == 0:
                            reason = (f"{mapped:,} reads mapped but max depth was only "
                                      f"{diag['max_depth']}× across "
                                      f"{diag['positions_seen']:,} positions "
                                      f"(need ≥ {self.min_depth}× to call variants) — "
                                      "try increasing Max reads or lowering Min depth")
                        else:
                            reason = (f"{mapped:,} reads mapped, "
                                      f"{diag['positions_above_depth']:,} positions "
                                      f"above depth threshold, but no allele frequency "
                                      f"≥ {self.min_freq:.0%} — "
                                      "sample may be nearly identical to the reference")
                        msg = f"{srr}: {reason}"
                        self._sample_errors.append(msg)
                        self.sample_error.emit(f"Warning: {msg}")
                    else:
                        df["Gene"] = df["Position"].apply(_lookup_gene)
                        df = self._annotate_variants(df, cds_map, ref_seq)
                        n_vars = len(df)
                        self.status_updated.emit(
                            f"[4/4] Sample {i+1}/{n} — {n_vars} variant call{'s' if n_vars != 1 else ''} found."
                        )
                        master_list.append(df)
                except InterruptedError:
                    raise
                except Exception as e:
                    msg = f"{srr}: {e}"
                    self._sample_errors.append(msg)
                    self.sample_error.emit(f"Skipped {msg}")
                finally:
                    for tmp in [bam, bam + ".bai"]:
                        if os.path.exists(tmp):
                            os.remove(tmp)

                self._sample_times.append(time.time() - sample_start)
                self._emit_progress(min(95, base + int(per_sample)))

            if not master_list:
                detail = "\n\n".join(self._sample_errors) if self._sample_errors else "No error details available."
                self.error_occurred.emit(
                    f"No usable variant calls produced from {len(self.selected_srr)} sample(s).\n\n"
                    f"{detail}\n\n"
                    f"Try increasing Max reads or lowering min_depth in _call_snps."
                )
                return

            # 5. Summarise
            self.status_updated.emit("Compiling and summarising all results...")
            master_df = pd.concat(master_list, ignore_index=True)
            self._emit_progress(98)

            self.finished.emit({
                "master_df": master_df,
                "gene_map": gene_map,
                "accession_id": self.accession_id,
            })
            self._emit_progress(100)
            total = self._fmt_seconds(time.time() - self._start_time)
            self.status_updated.emit(f"Analysis complete! ({len(master_list)} sample(s), {total} total)")

        except InterruptedError:
            self.status_updated.emit("Analysis aborted.")
            self.progress_updated.emit(0)
            self.aborted.emit()
        except Exception as e:
            self.error_occurred.emit(str(e))
        finally:
            for f in ["ref.mmi", "ref.fasta"]:
                if os.path.exists(f):
                    os.remove(f)

    # ── gene map ─────────────────────────────────────────────────────────

    def _build_gene_map(self, accession_id):
        with Entrez.efetch(db="nucleotide", id=accession_id,
                           rettype="gb", retmode="text") as h:
            record = SeqIO.read(h, "genbank")

        # RefSeq bacterial chromosomes are CON records with no features of their
        # own — annotations live in the primary GenBank submission.  Parse the
        # primary accession from the CONTIG field and re-fetch.
        has_features = any(f.type in ("CDS", "gene") for f in record.features)
        if not has_features:
            contig = record.annotations.get("contig", "")
            m = re.search(r'([A-Z]{1,2}_?\d{5,9}\.\d+)', contig)
            if m:
                primary = m.group(1)
                with Entrez.efetch(db="nucleotide", id=primary,
                                   rettype="gb", retmode="text") as h:
                    record = SeqIO.read(h, "genbank")

        gene_map = {}   # name -> (start, end)  1-based, used for region lookup
        cds_map  = {}   # name -> {start, end, strand, codon_start}  CDS only
        for feat in record.features:
            if feat.type not in ("CDS", "gene"):
                continue
            quals = feat.qualifiers
            name = (
                quals.get("gene", quals.get("product", [None]))[0]
                or quals.get("product", [None])[0]
            )
            if not name:
                continue
            start = int(feat.location.start) + 1   # 1-based inclusive
            end   = int(feat.location.end)          # 1-based inclusive (= 0-based exclusive)
            strand = feat.location.strand if feat.location.strand is not None else 1
            # Prefer CDS entries; avoid duplicating if gene + CDS both present
            if name not in gene_map or feat.type == "CDS":
                gene_map[name] = (start, end)
            if feat.type == "CDS":
                codon_start = int(quals.get("codon_start", [1])[0])
                cds_map[name] = {
                    "start": start, "end": end,
                    "strand": strand, "codon_start": codon_start,
                }
        return gene_map, cds_map

    # ── download + align ─────────────────────────────────────────────────

    def _annotate_variants(self, df, cds_map, ref_seq):
        """Add Mutation_Effect and AA_Change columns to the variant DataFrame.

        For each SNP inside a CDS, extracts the affected codon from the reference
        sequence, substitutes the alternate base (taking strand into account), and
        translates both codons to classify the change.
        """
        _comp = str.maketrans("ACGTacgt", "TGCAtgca")

        effects, aa_changes = [], []

        for _, row in df.iterrows():
            pos  = int(row["Position"])      # 1-based genome position
            base = str(row["Base"])
            vtype = row["Variant_Type"]
            gene  = row.get("Gene", "Intergenic")

            # Reference base at this genome position
            ref_base = ref_seq[pos - 1].upper() if pos - 1 < len(ref_seq) else "N"

            if ref_base == base.upper() and vtype == "SNP":
                effects.append("Reference allele")
                aa_changes.append("—")
                continue

            if gene == "Intergenic" or gene not in cds_map:
                effects.append("Intergenic" if gene == "Intergenic" else "Non-coding")
                aa_changes.append("—")
                continue

            if vtype == "INDEL":
                effects.append("Frameshift (INDEL)")
                aa_changes.append("—")
                continue

            cds   = cds_map[gene]
            start = cds["start"]        # 1-based inclusive
            end   = cds["end"]          # 1-based inclusive (= 0-based exclusive)
            strand      = cds["strand"]
            codon_start = cds.get("codon_start", 1)

            try:
                if strand == 1:
                    cds_offset = (pos - start) - (codon_start - 1)
                else:
                    cds_offset = (end - pos) - (codon_start - 1)

                if cds_offset < 0:
                    effects.append("Non-coding")
                    aa_changes.append("—")
                    continue

                codon_idx = cds_offset // 3
                pos_in_codon = cds_offset % 3

                if strand == 1:
                    cs = start - 1 + (codon_start - 1) + codon_idx * 3   # 0-based
                    ref_codon = ref_seq[cs:cs + 3].upper()
                    alt_codon = ref_codon[:pos_in_codon] + base.upper() + ref_codon[pos_in_codon + 1:]
                else:
                    # Codon on the coding strand sits at the high end of the genome
                    ge = end - (codon_start - 1) - codon_idx * 3         # 0-based exclusive
                    ref_codon_fwd = ref_seq[ge - 3:ge].upper()
                    ref_codon = ref_codon_fwd.translate(_comp)[::-1]
                    alt_base_coding = base.upper().translate(_comp)
                    alt_codon = ref_codon[:pos_in_codon] + alt_base_coding + ref_codon[pos_in_codon + 1:]

                if len(ref_codon) != 3 or len(alt_codon) != 3:
                    effects.append("Unknown")
                    aa_changes.append("—")
                    continue

                ref_aa = str(Seq(ref_codon).translate())
                alt_aa = str(Seq(alt_codon).translate())
                aa_pos = codon_idx + 1

                if ref_aa == alt_aa:
                    effects.append("Synonymous")
                    aa_changes.append(f"{ref_aa}{aa_pos}{alt_aa} (silent)")
                elif alt_aa == "*":
                    effects.append("Stop gained")
                    aa_changes.append(f"{ref_aa}{aa_pos}*")
                elif ref_aa == "*":
                    effects.append("Stop lost")
                    aa_changes.append(f"*{aa_pos}{alt_aa}")
                elif aa_pos == 1 and ref_aa == "M":
                    effects.append("Start lost")
                    aa_changes.append(f"M1{alt_aa}")
                else:
                    effects.append("Non-synonymous")
                    aa_changes.append(f"{ref_aa}{aa_pos}{alt_aa}")

            except Exception:
                effects.append("Unknown")
                aa_changes.append("—")

        out = df.copy()
        out["Mutation_Effect"] = effects
        out["AA_Change"]       = aa_changes
        return out

    def _download_and_align(self, srr_id, sample_num, sample_total, base_pct, span):
        bam_sorted = f"{srr_id}_sorted.bam"
        prefix = f"[4/4] Sample {sample_num}/{sample_total} ({srr_id})"
        threads = min(multiprocessing.cpu_count(), 8)
        ref = "ref.mmi" if os.path.exists("ref.mmi") else "ref.fasta"

        platform = self.platform_map.get(srr_id, DEFAULT_PLATFORM)
        preset = PLATFORM_PRESETS.get(platform, "sr")
        plat_label = PLATFORM_LABELS.get(platform, platform.capitalize())

        reads_flag = f"-X {self.max_reads}" if self.max_reads > 0 else ""
        cap_note = f" (capped at {self.max_reads:,} reads)" if self.max_reads > 0 else ""
        cmd = (
            f"fastq-dump {reads_flag} --stdout --skip-technical {srr_id} | "
            f"minimap2 -ax {preset} -t {threads} {ref} - | "
            f"samtools view -bS -F 4 | "
            f"samtools sort -o {bam_sorted}"
        )

        self.status_updated.emit(
            f"{prefix} — streaming {plat_label} reads → align [{preset}]{cap_note}..."
        )
        self._run_proc(cmd)
        self._emit_progress(base_pct + int(span * 0.90))

        self.status_updated.emit(f"{prefix} — indexing BAM...")
        self._run_proc(f"samtools index {bam_sorted}")
        self._emit_progress(base_pct + int(span * 0.95))

    # ── SNP calling ──────────────────────────────────────────────────────

    def _call_snps(self, bam_path, accession_id, min_depth=20, min_freq=0.05):
        bam = pysam.AlignmentFile(bam_path, "rb")

        mapped_reads = bam.mapped  # from BAM index

        # Resolve reference name from BAM header
        ref_names = list(bam.references)
        if not ref_names:
            bam.close()
            raise ValueError(f"BAM {bam_path} has no mapped references — check alignment.")
        matching = [r for r in ref_names if accession_id in r or r in accession_id]
        ref_name = matching[0] if matching else ref_names[0]

        results = []
        positions_seen = 0
        positions_above_depth = 0
        max_depth = 0

        for col in bam.pileup(ref_name, min_mapping_quality=20):
            positions_seen += 1
            depth = col.nsegments
            max_depth = max(max_depth, depth)
            if depth < min_depth:
                continue
            positions_above_depth += 1
            pos = col.pos + 1

            counts = {"A": 0, "T": 0, "C": 0, "G": 0, "Indel": 0}
            for read in col.pileups:
                if read.is_del or read.is_refskip:
                    counts["Indel"] += 1
                elif read.indel != 0:
                    counts["Indel"] += 1
                else:
                    base = read.alignment.query_sequence[read.query_position]
                    if base in counts:
                        counts[base] += 1

            for var, count in counts.items():
                freq = count / depth
                if freq >= min_freq:
                    results.append({
                        "Position": pos,
                        "Variant_Type": "SNP" if var != "Indel" else "INDEL",
                        "Base": var,
                        "Count": count,
                        "Depth": depth,
                        "Frequency": round(freq, 4),
                    })
        bam.close()

        diag = {
            "mapped_reads":         mapped_reads,
            "positions_seen":       positions_seen,
            "positions_above_depth": positions_above_depth,
            "max_depth":            max_depth,
        }
        _cols = ["Position", "Variant_Type", "Base", "Count", "Depth", "Frequency"]
        df = pd.DataFrame(results, columns=_cols) if results else pd.DataFrame(columns=_cols)
        return df, diag


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------
class SNPAnalyzerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SNP Frequency Analyzer")
        self.setGeometry(100, 100, 1300, 950)
        self.df_candidates = pd.DataFrame()
        self.srr_checkboxes = []
        self.current_results = None
        self._worker = None
        self._fetch_worker = None
        self._analysis_start = None

        self._elapsed_timer = QTimer(self)
        self._elapsed_timer.setInterval(1000)
        self._elapsed_timer.timeout.connect(self._tick_elapsed)

        self._init_ui()

    # ── UI construction ──────────────────────────────────────────────────

    def _init_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setSpacing(8)
        outer.setContentsMargins(10, 10, 10, 10)

        splitter = QSplitter(Qt.Vertical)
        outer.addWidget(splitter)

        # ── TOP PANEL (inputs) ──────────────────────────────────────────
        top = QWidget()
        top_layout = QVBoxLayout(top)
        top_layout.setSpacing(6)
        top_layout.setContentsMargins(0, 0, 0, 0)

        # Email row
        email_group = QGroupBox("NCBI Email (required for Entrez searches)")
        email_layout = QHBoxLayout()
        self.email_input = QLineEdit()
        self.email_input.setPlaceholderText("Enter your email address (e.g. you@example.com)")
        email_layout.addWidget(self.email_input)
        email_group.setLayout(email_layout)
        top_layout.addWidget(email_group)

        # Accession input group
        acc_group = QGroupBox("Reference Accession")
        acc_outer = QVBoxLayout()

        # Mode selector
        mode_row = QHBoxLayout()
        self._mode_group = QButtonGroup()
        self._rb_preset = QRadioButton("Select from preset organisms")
        self._rb_manual = QRadioButton("Enter accession manually")
        self._rb_preset.setChecked(True)
        self._mode_group.addButton(self._rb_preset, 0)
        self._mode_group.addButton(self._rb_manual, 1)
        mode_row.addWidget(self._rb_preset)
        mode_row.addWidget(self._rb_manual)
        mode_row.addStretch()
        acc_outer.addLayout(mode_row)

        # Preset dropdown
        self._preset_widget = QWidget()
        preset_row = QHBoxLayout(self._preset_widget)
        preset_row.setContentsMargins(0, 0, 0, 0)
        self._organism_combo = QComboBox()
        self._organism_combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        for label in PRESET_ORGANISMS:
            self._organism_combo.addItem(label)
        self._organism_combo.currentIndexChanged.connect(self._on_preset_changed)
        preset_row.addWidget(QLabel("Organism:"))
        preset_row.addWidget(self._organism_combo)
        acc_outer.addWidget(self._preset_widget)

        # Manual input
        self._manual_widget = QWidget()
        manual_row = QHBoxLayout(self._manual_widget)
        manual_row.setContentsMargins(0, 0, 0, 0)
        self._manual_input = QLineEdit()
        self._manual_input.setPlaceholderText("e.g. NC_045512.2")
        manual_row.addWidget(QLabel("Accession ID:"))
        manual_row.addWidget(self._manual_input)
        self._manual_widget.setVisible(False)
        acc_outer.addWidget(self._manual_widget)

        self._mode_group.buttonClicked.connect(self._on_mode_changed)

        # Resolved accession display
        resolved_row = QHBoxLayout()
        resolved_row.addWidget(QLabel("Resolved accession:"))
        self._resolved_label = QLabel("—")
        self._resolved_label.setFont(QFont("Monospace", 10))
        resolved_row.addWidget(self._resolved_label)
        resolved_row.addStretch()
        acc_outer.addLayout(resolved_row)

        acc_group.setLayout(acc_outer)
        top_layout.addWidget(acc_group)

        # SRR fetch controls
        fetch_group = QGroupBox("Sequencing Run Discovery")
        fetch_outer = QVBoxLayout()

        fetch_row = QHBoxLayout()
        fetch_row.addWidget(QLabel("Max results:"))
        self._limit_spin = QSpinBox()
        self._limit_spin.setRange(1, 50)
        self._limit_spin.setValue(10)
        fetch_row.addWidget(self._limit_spin)
        self._fetch_btn = QPushButton("Fetch Runs")
        self._fetch_btn.clicked.connect(self._fetch_srr)
        fetch_row.addWidget(self._fetch_btn)

        fetch_row.addSpacing(20)
        fetch_row.addWidget(QLabel("Max reads/sample:"))
        self._max_reads_spin = QSpinBox()
        self._max_reads_spin.setRange(0, 10_000_000)
        self._max_reads_spin.setSingleStep(100_000)
        self._max_reads_spin.setValue(500_000)
        self._max_reads_spin.setSpecialValueText("all reads")
        self._max_reads_spin.setToolTip(
            "Limit reads downloaded per sample.\n"
            "500 K is usually enough for reliable frequency estimates.\n"
            "Set to 0 to download all reads."
        )
        fetch_row.addWidget(self._max_reads_spin)
        fetch_row.addStretch()
        fetch_outer.addLayout(fetch_row)

        # Note about DRR filtering
        note = QLabel(
            "Note: Only SRR (NCBI) and ERR (ENA) runs are shown — DRR (DDBJ) runs are excluded. "
            "Runs without a size are dimmed and may be inaccessible."
        )
        note.setStyleSheet("color: #888; font-size: 11px;")
        note.setWordWrap(True)
        fetch_outer.addWidget(note)

        # SRR checkbox list
        self._srr_scroll = QScrollArea()
        self._srr_scroll.setWidgetResizable(True)
        self._srr_scroll.setMaximumHeight(180)
        self._srr_container = QWidget()
        self._srr_vbox = QVBoxLayout(self._srr_container)
        self._srr_scroll.setWidget(self._srr_container)
        fetch_outer.addWidget(self._srr_scroll)

        # Select all / none row
        sel_row = QHBoxLayout()
        self._select_all_btn = QPushButton("Select All")
        self._select_all_btn.clicked.connect(lambda: self._set_all_checked(True))
        self._select_none_btn = QPushButton("Deselect All")
        self._select_none_btn.clicked.connect(lambda: self._set_all_checked(False))
        sel_row.addWidget(self._select_all_btn)
        sel_row.addWidget(self._select_none_btn)
        sel_row.addStretch()
        fetch_outer.addLayout(sel_row)

        fetch_group.setLayout(fetch_outer)
        top_layout.addWidget(fetch_group)

        # Progress
        prog_group = QGroupBox("Analysis Progress")
        prog_layout = QVBoxLayout()

        self._progress_bar = QProgressBar()
        self._progress_bar.setTextVisible(True)
        prog_layout.addWidget(self._progress_bar)

        time_row = QHBoxLayout()
        self._status_label = QLabel("Ready")
        self._elapsed_label = QLabel("Elapsed: —")
        self._eta_label = QLabel("ETA: —")
        time_row.addWidget(self._status_label, 3)
        time_row.addWidget(self._elapsed_label, 1)
        time_row.addWidget(self._eta_label, 1)
        prog_layout.addLayout(time_row)

        prog_group.setLayout(prog_layout)
        top_layout.addWidget(prog_group)

        # Variant calling thresholds
        thresh_group = QGroupBox("Variant Calling Thresholds")
        thresh_row = QHBoxLayout()

        thresh_row.addWidget(QLabel("Min depth:"))
        self._min_depth_spin = QSpinBox()
        self._min_depth_spin.setRange(1, 500)
        self._min_depth_spin.setValue(20)
        self._min_depth_spin.setToolTip(
            "Minimum read depth at a position to call a variant.\n"
            "Lower this (e.g. 5–10) for large bacterial genomes\n"
            "or when using a small read cap."
        )
        thresh_row.addWidget(self._min_depth_spin)

        thresh_row.addSpacing(16)
        thresh_row.addWidget(QLabel("Min frequency:"))
        self._min_freq_spin = QDoubleSpinBox()
        self._min_freq_spin.setRange(0.01, 0.5)
        self._min_freq_spin.setSingleStep(0.01)
        self._min_freq_spin.setDecimals(2)
        self._min_freq_spin.setValue(0.05)
        self._min_freq_spin.setToolTip(
            "Minimum allele frequency (0–1) to report a variant.\n"
            "0.05 = 5%.  Lower for rare-variant detection."
        )
        thresh_row.addWidget(self._min_freq_spin)
        thresh_row.addStretch()
        thresh_group.setLayout(thresh_row)
        top_layout.addWidget(thresh_group)

        # Run / Abort buttons
        run_row = QHBoxLayout()
        self._run_btn = QPushButton("Run Analysis")
        self._run_btn.setObjectName("run_btn")
        self._run_btn.setEnabled(False)
        self._run_btn.setMinimumHeight(36)
        self._run_btn.clicked.connect(self._run_analysis)
        run_row.addWidget(self._run_btn)

        self._abort_btn = QPushButton("Abort")
        self._abort_btn.setObjectName("abort_btn")
        self._abort_btn.setEnabled(False)
        self._abort_btn.setMinimumHeight(36)
        self._abort_btn.clicked.connect(self._abort_analysis)
        run_row.addWidget(self._abort_btn)
        top_layout.addLayout(run_row)

        splitter.addWidget(top)

        # ── BOTTOM PANEL (results) ──────────────────────────────────────
        results_scroll = QScrollArea()
        results_scroll.setWidgetResizable(True)
        results_inner = QWidget()
        results_layout = QVBoxLayout(results_inner)
        results_layout.setSpacing(8)

        # Plot
        plot_group = QGroupBox("SNP Frequency Plot")
        plot_vbox = QVBoxLayout()
        self._fig = Figure(figsize=(12, 5), dpi=100)
        self._canvas = FigureCanvas(self._fig)
        self._canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._canvas.setMinimumHeight(320)
        plot_vbox.addWidget(self._canvas)

        # Export button
        export_row = QHBoxLayout()
        self._export_plot_btn = QPushButton("Save Plot as PNG")
        self._export_plot_btn.clicked.connect(self._export_plot)
        self._export_csv_btn = QPushButton("Export Results as CSV")
        self._export_csv_btn.clicked.connect(self._export_csv)
        export_row.addWidget(self._export_plot_btn)
        export_row.addWidget(self._export_csv_btn)
        export_row.addStretch()
        plot_vbox.addLayout(export_row)
        plot_group.setLayout(plot_vbox)
        results_layout.addWidget(plot_group)

        # SNP by Gene table
        gene_group = QGroupBox("SNP / INDEL Count by Gene")
        gene_vbox = QVBoxLayout()
        self._gene_table = QTableWidget()
        self._gene_table.setMinimumHeight(120)
        self._gene_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        gene_vbox.addWidget(self._gene_table)
        gene_group.setLayout(gene_vbox)
        results_layout.addWidget(gene_group)

        # High frequency mutations table
        hf_group = QGroupBox("High-Frequency Mutations (>80%)")
        hf_vbox = QVBoxLayout()
        self._hf_table = QTableWidget()
        self._hf_table.setMinimumHeight(120)
        self._hf_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        hf_vbox.addWidget(self._hf_table)
        hf_group.setLayout(hf_vbox)
        results_layout.addWidget(hf_group)

        # All variants table
        all_group = QGroupBox("All Variant Calls")
        all_vbox = QVBoxLayout()
        self._all_table = QTableWidget()
        self._all_table.setMinimumHeight(150)
        self._all_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        all_vbox.addWidget(self._all_table)
        all_group.setLayout(all_vbox)
        results_layout.addWidget(all_group)

        results_layout.addStretch()
        results_scroll.setWidget(results_inner)
        splitter.addWidget(results_scroll)

        splitter.setSizes([420, 680])

    # ── UI event handlers ────────────────────────────────────────────────

    def _on_mode_changed(self, btn):
        is_preset = self._mode_group.checkedId() == 0
        self._preset_widget.setVisible(is_preset)
        self._manual_widget.setVisible(not is_preset)
        if is_preset:
            self._on_preset_changed()
        else:
            self._resolved_label.setText("—")

    def _on_preset_changed(self):
        label = self._organism_combo.currentText()
        acc = PRESET_ORGANISMS.get(label, "")
        self._resolved_label.setText(acc if acc else "—")

    def _get_accession(self):
        if self._mode_group.checkedId() == 0:
            label = self._organism_combo.currentText()
            return PRESET_ORGANISMS.get(label, "").strip()
        return self._manual_input.text().strip()

    def _get_email(self):
        return self.email_input.text().strip()

    def _validate_email(self):
        email = self._get_email()
        if not email or "@" not in email:
            QMessageBox.warning(self, "Email Required",
                                "Please enter a valid email address.\n"
                                "NCBI requires an email for Entrez API usage.")
            return False
        return True

    def _set_all_checked(self, state):
        for cb in self.srr_checkboxes:
            cb.setChecked(state)

    # ── SRR fetch ────────────────────────────────────────────────────────

    def _fetch_srr(self):
        if not self._validate_email():
            return
        acc = self._get_accession()
        if not acc:
            QMessageBox.warning(self, "No Accession",
                                "Please select an organism or enter an accession ID.")
            return

        self._fetch_btn.setEnabled(False)
        self._status_label.setText(f"Searching SRA for {acc}...")

        self._fetch_worker = SRRFetchWorker(acc, self._limit_spin.value(), self._get_email())
        self._fetch_worker.finished.connect(self._on_srr_fetched)
        self._fetch_worker.error_occurred.connect(self._on_fetch_error)
        self._fetch_worker.status_updated.connect(self._status_label.setText)
        self._fetch_worker.start()

    def _on_srr_fetched(self, df):
        self.df_candidates = df
        self._fetch_btn.setEnabled(True)
        if df.empty:
            self._status_label.setText("No accessible SRR/ERR runs found for this accession.")
            self._run_btn.setEnabled(False)
            QMessageBox.information(self, "No Results",
                                    "No SRR or ERR runs were found linked to this accession.\n\n"
                                    "Possible reasons:\n"
                                    "• The accession may not have associated SRA data\n"
                                    "• Try a more commonly sequenced reference genome")
        else:
            self._status_label.setText(f"Found {len(df)} accessible run(s).")
            self._populate_srr_checkboxes(df)
            self._run_btn.setEnabled(True)

    def _on_fetch_error(self, msg):
        self._fetch_btn.setEnabled(True)
        self._status_label.setText("Fetch error — see dialog.")
        QMessageBox.critical(self, "Fetch Error", msg)

    def _populate_srr_checkboxes(self, df):
        # Clear old checkboxes
        while self._srr_vbox.count():
            item = self._srr_vbox.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.srr_checkboxes = []

        for _, row in df.iterrows():
            srr_id = str(row["SRR_ID"])
            title = str(row.get("Title", ""))[:55]
            date = str(row.get("Date", ""))
            platform_raw = str(row.get("Platform", DEFAULT_PLATFORM))
            platform_label = PLATFORM_LABELS.get(platform_raw, platform_raw.capitalize())
            raw_bases = row.get("Bases", "N/A")
            try:
                n = int(raw_bases)
                if n >= 1_000_000_000:
                    bases_str = f"{n / 1e9:.1f} Gbp"
                elif n >= 1_000_000:
                    bases_str = f"{n / 1e6:.1f} Mbp"
                elif n >= 1_000:
                    bases_str = f"{n / 1e3:.1f} Kbp"
                else:
                    bases_str = f"{n} bp"
                has_size = True
            except (ValueError, TypeError):
                bases_str = "size unknown"
                has_size = False
            label = f"{srr_id}   {title}   ({date}, {bases_str}, {platform_label})"
            cb = QCheckBox(label)
            cb.setProperty("srr_platform", platform_raw)
            try:
                cb.setProperty("srr_bases", int(raw_bases))
            except (ValueError, TypeError):
                cb.setProperty("srr_bases", -1)
            if not has_size:
                cb.setStyleSheet("color: #999; font-style: italic;")
                cb.setToolTip("No size data from NCBI — this run may be inaccessible or incompletely deposited")
            self.srr_checkboxes.append(cb)
            self._srr_vbox.addWidget(cb)

        self._srr_vbox.addStretch()

    # ── Analysis ─────────────────────────────────────────────────────────

    def _run_analysis(self):
        if not self._validate_email():
            return
        acc = self._get_accession()
        if not acc:
            QMessageBox.warning(self, "No Accession", "Please select or enter an accession ID.")
            return

        checked = [cb for cb in self.srr_checkboxes if cb.isChecked()]
        if not checked:
            QMessageBox.warning(self, "No Selection", "Please select at least one SRR run.")
            return
        selected = [cb.text().split()[0] for cb in checked]

        if self._max_reads_spin.value() == 0:
            large = [
                cb.text().split()[0]
                for cb in checked
                if (cb.property("srr_bases") or -1) > 500_000_000
            ]
            if large:
                names = ", ".join(large)
                reply = QMessageBox.warning(
                    self, "Large Run — No Read Cap",
                    f"Max reads is set to <b>all reads</b> and the following run(s) "
                    f"exceed 500 Mbp:<br><br><b>{names}</b><br><br>"
                    f"This may take <b>10–30 minutes per sample</b> depending on "
                    f"network speed and CPU.<br><br>"
                    f"Consider setting Max reads to 500,000–1,000,000 for a faster result. "
                    f"Continue anyway?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.No,
                )
                if reply != QMessageBox.Yes:
                    return

        platform_map = {
            cb.text().split()[0]: (cb.property("srr_platform") or DEFAULT_PLATFORM)
            for cb in checked
        }

        self._run_btn.setEnabled(False)
        self._abort_btn.setEnabled(True)
        self._fetch_btn.setEnabled(False)
        self._progress_bar.setValue(0)
        self._status_label.setText("Starting analysis...")

        self._analysis_start = time.time()
        self._elapsed_timer.start()

        self._worker = AnalysisWorker(
            acc, selected, self._get_email(),
            max_reads=self._max_reads_spin.value(),
            platform_map=platform_map,
            min_depth=self._min_depth_spin.value(),
            min_freq=self._min_freq_spin.value(),
        )
        self._worker.progress_updated.connect(self._progress_bar.setValue)
        self._worker.status_updated.connect(self._status_label.setText)
        self._worker.time_updated.connect(self._update_time_labels)
        self._worker.finished.connect(self._display_results)
        self._worker.error_occurred.connect(self._on_analysis_error)
        self._worker.aborted.connect(self._on_analysis_aborted)
        self._worker.sample_error.connect(
            lambda msg: self._status_label.setText(f"Warning: {msg}")
        )
        self._worker.start()

    def _abort_analysis(self):
        if self._worker and self._worker.isRunning():
            self._abort_btn.setEnabled(False)
            self._abort_btn.setText("Aborting...")
            self._worker.abort()

    def _reset_controls(self):
        self._elapsed_timer.stop()
        self._analysis_start = None
        self._run_btn.setEnabled(True)
        self._abort_btn.setEnabled(False)
        self._abort_btn.setText("Abort")
        self._fetch_btn.setEnabled(True)
        self._elapsed_label.setText("Elapsed: —")
        self._eta_label.setText("ETA: —")

    def _tick_elapsed(self):
        if self._analysis_start is not None:
            secs = time.time() - self._analysis_start
            m, s = divmod(int(secs), 60)
            self._elapsed_label.setText(f"Elapsed: {m}m {s:02d}s" if m else f"Elapsed: {s}s")

    def _on_analysis_aborted(self):
        self._reset_controls()

    def _update_time_labels(self, _elapsed, eta):
        # elapsed is driven by the 1-second QTimer; only update ETA here
        self._eta_label.setText(f"ETA: {eta}")

    def _on_analysis_error(self, msg):
        self._reset_controls()
        self._status_label.setText("Analysis failed.")
        QMessageBox.critical(self, "Analysis Error", msg)

    # ── Results display ──────────────────────────────────────────────────

    def _display_results(self, results):
        self._reset_controls()
        self.current_results = results

        master_df = results["master_df"]
        gene_map = results["gene_map"]
        acc = results["accession_id"]

        # Exclude reference alleles from display — they are not mutations and
        # add noise to every table. Full data is still available via CSV export.
        variants_df = (
            master_df[master_df["Mutation_Effect"] != "Reference allele"]
            if "Mutation_Effect" in master_df.columns
            else master_df
        )

        self._draw_plot(variants_df, gene_map, acc)
        self._fill_gene_table(variants_df)
        self._fill_hf_table(variants_df)
        self._fill_all_table(variants_df)

    def _draw_plot(self, df, gene_map, acc):
        self._fig.clear()

        # Two axes: top = scatter, bottom = gene map track
        ax_snp = self._fig.add_axes([0.07, 0.30, 0.88, 0.60])
        ax_gene = self._fig.add_axes([0.07, 0.05, 0.88, 0.18], sharex=ax_snp)

        # ── scatter plot ────────────────────────────────────────────────
        palette = {"SNP": "#E74C3C", "INDEL": "#3498DB"}
        markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
        samples = df["Sample_ID"].unique()

        for i, sample in enumerate(samples):
            sub = df[df["Sample_ID"] == sample]
            snps = sub[sub["Variant_Type"] == "SNP"]
            indels = sub[sub["Variant_Type"] == "INDEL"]
            m = markers[i % len(markers)]
            if not snps.empty:
                ax_snp.scatter(snps["Position"], snps["Frequency"],
                               color=palette["SNP"], marker=m, alpha=0.6,
                               s=28, label=f"{sample} SNP" if i == 0 else "_")
            if not indels.empty:
                ax_snp.scatter(indels["Position"], indels["Frequency"],
                               color=palette["INDEL"], marker=m, alpha=0.6,
                               s=28, label=f"{sample} INDEL" if i == 0 else "_")

        # Restrict gene map to genes that actually contain variants so large
        # bacterial genomes (1000+ genes) don't flood the plot with shading.
        variant_genes = set(df["Gene"].unique()) - {"Intergenic"}
        plot_gene_map = (
            {g: v for g, v in gene_map.items() if g in variant_genes}
            if variant_genes else gene_map
        )

        # Gene region shading on SNP axis
        colors = sns.color_palette("Set2", max(len(plot_gene_map), 1))
        for i, (gene, (start, end)) in enumerate(plot_gene_map.items()):
            ax_snp.axvspan(start, end, alpha=0.08, color=colors[i % len(colors)])

        ax_snp.set_ylabel("Allele Frequency", fontsize=10)
        ax_snp.set_title(f"SNP / INDEL Profile — {acc}  ({len(samples)} sample(s))",
                         fontsize=11, fontweight="bold", pad=8)
        ax_snp.set_ylim(-0.05, 1.1)
        ax_snp.axhline(0.8, color="#7f8c8d", linestyle="--", linewidth=0.8, alpha=0.8,
                       label="80% threshold")
        ax_snp.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6, color="#bdc3c7")
        ax_snp.set_axisbelow(True)
        ax_snp.tick_params(labelbottom=False)
        ax_snp.spines["top"].set_visible(False)
        ax_snp.spines["right"].set_visible(False)

        # Legend for SNP/INDEL types
        legend_handles = [
            mpatches.Patch(color=palette["SNP"], label="SNP"),
            mpatches.Patch(color=palette["INDEL"], label="INDEL"),
        ]
        ax_snp.legend(handles=legend_handles, loc="upper right", fontsize=8)

        # ── gene map track ──────────────────────────────────────────────
        ax_gene.set_ylim(0, 1)
        ax_gene.set_yticks([])
        ax_gene.set_xlabel("Genome Position (bp)", fontsize=10)
        ax_gene.spines["top"].set_visible(False)
        ax_gene.spines["right"].set_visible(False)
        ax_gene.spines["left"].set_visible(False)

        if plot_gene_map:
            for i, (gene, (start, end)) in enumerate(plot_gene_map.items()):
                color = colors[i % len(colors)]
                rect = mpatches.FancyBboxPatch(
                    (start, 0.15), end - start, 0.70,
                    boxstyle="round,pad=0.01",
                    facecolor=color, edgecolor="white", linewidth=0.5, alpha=0.85
                )
                ax_gene.add_patch(rect)
                mid = (start + end) / 2
                gene_label = gene if len(gene) <= 6 else gene[:5] + "…"
                ax_gene.text(mid, 0.5, gene_label,
                             ha="center", va="center",
                             fontsize=6, fontweight="bold", color="white",
                             clip_on=True)
            label = "Genes (with variants)" if variant_genes else "Genes"
            ax_gene.set_ylabel(label, fontsize=8)
        else:
            ax_gene.text(0.5, 0.5, "No gene annotations available",
                         transform=ax_gene.transAxes, ha="center", va="center",
                         color="grey", fontsize=9)

        self._canvas.draw()

    def _fill_gene_table(self, df):
        summary = (
            df.groupby(["Gene", "Variant_Type"])
            .agg(Count=("Position", "size"), Avg_Freq=("Frequency", "mean"))
            .reset_index()
            .sort_values("Count", ascending=False)
        )
        self._populate_table(self._gene_table, summary)
        self._fit_table(self._gene_table)

    def _fill_hf_table(self, df):
        hf = df[df["Frequency"] > 0.8].copy()
        if hf.empty:
            self._hf_table.setRowCount(1)
            self._hf_table.setColumnCount(1)
            self._hf_table.setHorizontalHeaderLabels(["Result"])
            self._hf_table.setItem(0, 0, QTableWidgetItem("No high-frequency mutations (>80%) found."))
            self._fit_table(self._hf_table)
            return
        ann_cols = [c for c in ["Mutation_Effect", "AA_Change"] if c in hf.columns]
        group_cols = ["Position", "Gene", "Base", "Variant_Type"] + ann_cols
        summary = (
            hf.groupby(group_cols)
            .agg(Avg_Frequency=("Frequency", "mean"), Samples=("Sample_ID", "nunique"))
            .reset_index()
            .sort_values("Avg_Frequency", ascending=False)
            .head(30)
        )
        self._populate_table(self._hf_table, summary)
        self._fit_table(self._hf_table)

    def _fill_all_table(self, df):
        display_cols = ["Sample_ID", "Position", "Gene", "Variant_Type",
                        "Base", "Mutation_Effect", "AA_Change",
                        "Count", "Depth", "Frequency"]
        cols = [c for c in display_cols if c in df.columns]
        self._populate_table(self._all_table, df[cols].sort_values("Frequency", ascending=False))
        self._fit_table(self._all_table, max_visible=12)

    @staticmethod
    def _populate_table(widget, df):
        widget.setRowCount(len(df))
        widget.setColumnCount(len(df.columns))
        widget.setHorizontalHeaderLabels(list(df.columns))
        widget.setAlternatingRowColors(True)
        widget.setShowGrid(True)
        widget.verticalHeader().setDefaultSectionSize(26)
        widget.verticalHeader().setVisible(False)
        for r, (_, row) in enumerate(df.iterrows()):
            for c, val in enumerate(row):
                txt = f"{val:.4f}" if isinstance(val, float) else str(val)
                item = QTableWidgetItem(txt)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                widget.setItem(r, c, item)

    @staticmethod
    def _fit_table(table, max_visible=8):
        """Resize table height to show up to max_visible rows, then enable internal scroll."""
        table.resizeRowsToContents()
        header_h = table.horizontalHeader().height()
        row_h = sum(table.rowHeight(r) for r in range(min(table.rowCount(), max_visible)))
        target = header_h + row_h + 4
        table.setMinimumHeight(target)
        table.setMaximumHeight(target)

    # ── Export ───────────────────────────────────────────────────────────

    def _export_plot(self):
        if self.current_results is None:
            QMessageBox.warning(self, "No Results", "Run an analysis first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Plot", "snp_plot.png",
                                              "PNG Images (*.png)")
        if path:
            self._fig.savefig(path, dpi=150, bbox_inches="tight")
            QMessageBox.information(self, "Saved", f"Plot saved to:\n{path}")

    def _export_csv(self):
        if self.current_results is None:
            QMessageBox.warning(self, "No Results", "Run an analysis first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", "snp_results.csv",
                                              "CSV Files (*.csv)")
        if path:
            self.current_results["master_df"].to_csv(path, index=False)
            QMessageBox.information(self, "Exported", f"Results exported to:\n{path}")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(STYLESHEET)
    window = SNPAnalyzerApp()
    window.show()
    sys.exit(app.exec_())
