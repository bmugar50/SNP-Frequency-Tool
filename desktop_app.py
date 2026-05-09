import sys
import os
import time
import pysam
import pandas as pd
import matplotlib.patches as mpatches
import seaborn as sns
import xml.etree.ElementTree as ET
from collections import OrderedDict
from Bio import SeqIO, Entrez
import subprocess

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLineEdit, QLabel, QTableWidget, QTableWidgetItem,
    QProgressBar, QCheckBox, QScrollArea, QFrame, QGroupBox,
    QMessageBox, QSpinBox, QRadioButton, QButtonGroup,
    QComboBox, QSplitter, QSizePolicy, QHeaderView
)
from PyQt5.QtCore import Qt, pyqtSignal, QThread
from PyQt5.QtGui import QFont, QColor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


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
    ("S. pneumoniae R6               NC_003098.1", "NC_003098.1"),
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
    "NC_003098.1": "Streptococcus pneumoniae[Organism]",
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
                {"SRR_ID": run_id, "Title": desc, "Date": "curated", "Bases": "N/A"}
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
                        root = ET.fromstring(f"<root>{s['Runs']}</root>")
                        for run_el in root.findall(".//Run"):
                            acc = run_el.get("acc", "")
                            if not acc.startswith(ACCESSIBLE_PREFIXES):
                                continue
                            rows.append({
                                "SRR_ID": acc,
                                "Title": s.get("Title", "No Title")[:60],
                                "Date": s.get("CreateDate", "Unknown"),
                                "Bases": run_el.get("total_bases", "N/A"),
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

    def __init__(self, accession_id, selected_srr, email):
        super().__init__()
        self.accession_id = accession_id
        self.selected_srr = selected_srr
        self.email = email
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

        # Use average completed-sample time for ETA once we have at least one sample
        if self._sample_times:
            avg = sum(self._sample_times) / len(self._sample_times)
            remaining_samples = len(self.selected_srr) - len(self._sample_times)
            eta_str = self._fmt_seconds(avg * remaining_samples) if remaining_samples > 0 else "almost done"
        elif pct >= 15:
            # Past the fast setup phase — linear estimate is more meaningful now
            samples_elapsed = elapsed - self._setup_elapsed
            samples_done_frac = max((pct - 15) / 75, 0.01)
            remaining = (samples_elapsed / samples_done_frac) * (1 - samples_done_frac)
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

    def _run_proc(self, cmd):
        """Run a shell command via Popen so it can be killed on abort."""
        self._current_proc = subprocess.Popen(
            cmd, shell=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        self._current_proc.wait()
        rc = self._current_proc.returncode
        self._current_proc = None
        if self._aborted:
            raise InterruptedError("Aborted")
        if rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)

    # ── main pipeline ────────────────────────────────────────────────────

    def run(self):
        Entrez.email = self.email
        self._start_time = time.time()
        self._setup_elapsed = 0

        try:
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

            # 3. Build gene map from GenBank record
            self.status_updated.emit("[3/4] Fetching gene annotations from GenBank...")
            gene_map = self._build_gene_map(self.accession_id)
            n_genes = len(gene_map)
            self.status_updated.emit(f"[3/4] Found {n_genes} annotated gene{'s' if n_genes != 1 else ''}.")
            if self._aborted:
                raise InterruptedError("Aborted")
            self._emit_progress(13)
            self._setup_elapsed = time.time() - self._start_time

            # 4. Process each SRR sample
            self.status_updated.emit(f"[4/4] Processing {n} sample{'s' if n > 1 else ''}...")
            master_list = []
            per_sample = 82 / n   # progress points allocated per sample

            for i, srr in enumerate(self.selected_srr):
                if self._aborted:
                    raise InterruptedError("Aborted")

                base = 13 + int(i * per_sample)
                sample_start = time.time()

                try:
                    bam = self._download_and_align(srr, i + 1, n, base, per_sample)
                    self.status_updated.emit(f"[4/4] Sample {i+1}/{n} — calling SNP variants...")
                    df = self._call_snps(bam, self.accession_id)
                    df["Sample_ID"] = srr
                    df["Gene"] = df["Position"].apply(
                        lambda x: next(
                            (g for g, (s, e) in gene_map.items() if s <= x <= e),
                            "Intergenic"
                        )
                    )
                    n_vars = len(df)
                    self.status_updated.emit(
                        f"[4/4] Sample {i+1}/{n} — {n_vars} variant call{'s' if n_vars != 1 else ''} found."
                    )
                    master_list.append(df)
                    for f in [bam, bam + ".bai"]:
                        if os.path.exists(f):
                            os.remove(f)
                except InterruptedError:
                    raise
                except Exception as e:
                    self.sample_error.emit(f"Skipped {srr}: {e}")
                    continue

                self._sample_times.append(time.time() - sample_start)
                self._emit_progress(min(95, base + int(per_sample)))

            if not master_list:
                self.error_occurred.emit(
                    "No samples could be processed.\n"
                    "Check that the selected SRR runs contain reads that align to this reference."
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
            for f in ["ref.mmi"]:
                if os.path.exists(f):
                    os.remove(f)

    # ── gene map ─────────────────────────────────────────────────────────

    def _build_gene_map(self, accession_id):
        with Entrez.efetch(db="nucleotide", id=accession_id,
                           rettype="gb", retmode="text") as h:
            record = SeqIO.read(h, "genbank")
        gene_map = {}
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
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            # Prefer CDS entries; avoid duplicating if gene + CDS both present
            if name not in gene_map or feat.type == "CDS":
                gene_map[name] = (start, end)
        return gene_map

    # ── download + align ─────────────────────────────────────────────────

    def _download_and_align(self, srr_id, sample_num, sample_total, base_pct, span):
        bam_sorted = f"{srr_id}_sorted.bam"
        prefix = f"[4/4] Sample {sample_num}/{sample_total} ({srr_id})"

        # Step A — download reads (~55% of per-sample time)
        self.status_updated.emit(f"{prefix} — downloading reads via fasterq-dump...")
        self._run_proc(f"fasterq-dump --split-files {srr_id}")
        self._emit_progress(base_pct + int(span * 0.55))

        # Find FASTQ file(s)
        single = f"{srr_id}.fastq"
        r1 = f"{srr_id}_1.fastq"
        r2 = f"{srr_id}_2.fastq"
        if os.path.exists(single):
            fastq_input = single
            layout = "single-end"
        elif os.path.exists(r1) and os.path.exists(r2):
            fastq_input = f"{r1} {r2}"
            layout = "paired-end"
        elif os.path.exists(r1):
            fastq_input = r1
            layout = "single-end (R1 only)"
        else:
            raise FileNotFoundError(f"No FASTQ files found after downloading {srr_id}")

        # Step B — align with minimap2 (~35% of per-sample time)
        ref = "ref.mmi" if os.path.exists("ref.mmi") else "ref.fasta"
        self.status_updated.emit(f"{prefix} — aligning {layout} reads with minimap2...")
        self._run_proc(
            f"minimap2 -a {ref} {fastq_input} | "
            f"samtools view -bS -F 4 | "
            f"samtools sort -o {bam_sorted}"
        )
        self._emit_progress(base_pct + int(span * 0.90))

        # Step C — index BAM (~10% of per-sample time)
        self.status_updated.emit(f"{prefix} — indexing BAM file...")
        self._run_proc(f"samtools index {bam_sorted}")
        self._emit_progress(base_pct + int(span * 0.95))

        # Clean FASTQ
        for fq in [single, r1, r2, f"{srr_id}.sam"]:
            if os.path.exists(fq):
                os.remove(fq)

        return bam_sorted

    # ── SNP calling ──────────────────────────────────────────────────────

    def _call_snps(self, bam_path, accession_id, min_depth=20, min_freq=0.05):
        bam = pysam.AlignmentFile(bam_path, "rb")

        # Resolve reference name from BAM header
        ref_names = list(bam.references)
        if not ref_names:
            bam.close()
            raise ValueError(f"BAM {bam_path} has no mapped references — check alignment.")
        matching = [r for r in ref_names if accession_id in r or r in accession_id]
        ref_name = matching[0] if matching else ref_names[0]

        results = []
        for col in bam.pileup(ref_name, min_mapping_quality=20):
            pos = col.pos + 1
            depth = col.nsegments
            if depth < min_depth:
                continue

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
        return pd.DataFrame(results)


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
        fetch_group = QGroupBox("SRR Sample Discovery")
        fetch_outer = QVBoxLayout()

        fetch_row = QHBoxLayout()
        fetch_row.addWidget(QLabel("Max results:"))
        self._limit_spin = QSpinBox()
        self._limit_spin.setRange(1, 50)
        self._limit_spin.setValue(10)
        fetch_row.addWidget(self._limit_spin)
        self._fetch_btn = QPushButton("Fetch SRR Codes")
        self._fetch_btn.clicked.connect(self._fetch_srr)
        fetch_row.addWidget(self._fetch_btn)
        fetch_row.addStretch()
        fetch_outer.addLayout(fetch_row)

        # Note about DRR filtering
        note = QLabel(
            "Note: DRR-prefixed runs (DDBJ) are excluded as they are typically inaccessible. "
            "Only SRR and ERR runs are shown."
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

        # Run / Abort buttons
        run_row = QHBoxLayout()
        self._run_btn = QPushButton("Run Analysis")
        self._run_btn.setEnabled(False)
        self._run_btn.setMinimumHeight(36)
        self._run_btn.setFont(QFont("Arial", 11, QFont.Bold))
        self._run_btn.clicked.connect(self._run_analysis)
        run_row.addWidget(self._run_btn)

        self._abort_btn = QPushButton("Abort")
        self._abort_btn.setEnabled(False)
        self._abort_btn.setMinimumHeight(36)
        self._abort_btn.setFont(QFont("Arial", 11, QFont.Bold))
        self._abort_btn.setStyleSheet("QPushButton { color: white; background-color: #c0392b; }"
                                      "QPushButton:disabled { background-color: #888; }")
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
        self._canvas.setMinimumHeight(420)
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
        self._gene_table.setMaximumHeight(200)
        self._gene_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        gene_vbox.addWidget(self._gene_table)
        gene_group.setLayout(gene_vbox)
        results_layout.addWidget(gene_group)

        # High frequency mutations table
        hf_group = QGroupBox("High-Frequency Mutations (>80%)")
        hf_vbox = QVBoxLayout()
        self._hf_table = QTableWidget()
        self._hf_table.setMaximumHeight(200)
        self._hf_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        hf_vbox.addWidget(self._hf_table)
        hf_group.setLayout(hf_vbox)
        results_layout.addWidget(hf_group)

        # All variants table
        all_group = QGroupBox("All Variant Calls")
        all_vbox = QVBoxLayout()
        self._all_table = QTableWidget()
        self._all_table.setMaximumHeight(250)
        self._all_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        all_vbox.addWidget(self._all_table)
        all_group.setLayout(all_vbox)
        results_layout.addWidget(all_group)

        results_layout.addStretch()
        results_scroll.setWidget(results_inner)
        splitter.addWidget(results_scroll)

        splitter.setSizes([520, 430])

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
            title = str(row.get("Title", ""))[:60]
            date = str(row.get("Date", ""))
            bases = str(row.get("Bases", ""))
            label = f"{srr_id}  |  {title}  |  {date}  |  {bases} bases"
            cb = QCheckBox(label)
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

        selected = [cb.text().split("|")[0].strip() for cb in self.srr_checkboxes if cb.isChecked()]
        if not selected:
            QMessageBox.warning(self, "No Selection", "Please select at least one SRR run.")
            return

        self._run_btn.setEnabled(False)
        self._abort_btn.setEnabled(True)
        self._fetch_btn.setEnabled(False)
        self._progress_bar.setValue(0)
        self._status_label.setText("Starting analysis...")

        self._worker = AnalysisWorker(acc, selected, self._get_email())
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

    def _on_analysis_aborted(self):
        self._run_btn.setEnabled(True)
        self._abort_btn.setEnabled(False)
        self._abort_btn.setText("Abort")
        self._fetch_btn.setEnabled(True)
        self._elapsed_label.setText("Elapsed: —")
        self._eta_label.setText("ETA: —")

    def _update_time_labels(self, elapsed, eta):
        self._elapsed_label.setText(f"Elapsed: {elapsed}")
        self._eta_label.setText(f"ETA: {eta}")

    def _on_analysis_error(self, msg):
        self._run_btn.setEnabled(True)
        self._abort_btn.setEnabled(False)
        self._abort_btn.setText("Abort")
        self._fetch_btn.setEnabled(True)
        self._status_label.setText("Analysis failed.")
        QMessageBox.critical(self, "Analysis Error", msg)

    # ── Results display ──────────────────────────────────────────────────

    def _display_results(self, results):
        self._run_btn.setEnabled(True)
        self._abort_btn.setEnabled(False)
        self._abort_btn.setText("Abort")
        self._fetch_btn.setEnabled(True)
        self.current_results = results

        master_df = results["master_df"]
        gene_map = results["gene_map"]
        acc = results["accession_id"]

        self._draw_plot(master_df, gene_map, acc)
        self._fill_gene_table(master_df)
        self._fill_hf_table(master_df)
        self._fill_all_table(master_df)

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
                               color=palette["SNP"], marker=m, alpha=0.55,
                               s=25, label=f"{sample} SNP" if i == 0 else "_")
            if not indels.empty:
                ax_snp.scatter(indels["Position"], indels["Frequency"],
                               color=palette["INDEL"], marker=m, alpha=0.55,
                               s=25, label=f"{sample} INDEL" if i == 0 else "_")

        # Gene region shading on SNP axis
        colors = sns.color_palette("Set2", max(len(gene_map), 1))
        for i, (gene, (start, end)) in enumerate(gene_map.items()):
            ax_snp.axvspan(start, end, alpha=0.07, color=colors[i % len(colors)])

        ax_snp.set_ylabel("Allele Frequency")
        ax_snp.set_title(f"SNP / INDEL Profile — {acc}  ({len(samples)} sample(s))")
        ax_snp.set_ylim(-0.05, 1.1)
        ax_snp.axhline(0.8, color="grey", linestyle="--", linewidth=0.7, alpha=0.7)
        ax_snp.tick_params(labelbottom=False)

        # Legend for SNP/INDEL types
        legend_handles = [
            mpatches.Patch(color=palette["SNP"], label="SNP"),
            mpatches.Patch(color=palette["INDEL"], label="INDEL"),
        ]
        ax_snp.legend(handles=legend_handles, loc="upper right", fontsize=8)

        # ── gene map track ──────────────────────────────────────────────
        ax_gene.set_ylim(0, 1)
        ax_gene.set_yticks([])
        ax_gene.set_xlabel("Genome Position (bp)")

        if gene_map:
            for i, (gene, (start, end)) in enumerate(gene_map.items()):
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
            ax_gene.set_ylabel("Genes", fontsize=8)
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

    def _fill_hf_table(self, df):
        hf = df[df["Frequency"] > 0.8].copy()
        if hf.empty:
            self._hf_table.setRowCount(1)
            self._hf_table.setColumnCount(1)
            self._hf_table.setHorizontalHeaderLabels(["Result"])
            self._hf_table.setItem(0, 0, QTableWidgetItem("No high-frequency mutations (>80%) found."))
            return
        summary = (
            hf.groupby(["Position", "Gene", "Base", "Variant_Type"])
            .agg(Avg_Frequency=("Frequency", "mean"), Samples=("Sample_ID", "nunique"))
            .reset_index()
            .sort_values("Avg_Frequency", ascending=False)
            .head(30)
        )
        self._populate_table(self._hf_table, summary)

    def _fill_all_table(self, df):
        display_cols = ["Sample_ID", "Position", "Gene", "Variant_Type",
                        "Base", "Count", "Depth", "Frequency"]
        cols = [c for c in display_cols if c in df.columns]
        self._populate_table(self._all_table, df[cols].sort_values("Frequency", ascending=False))

    @staticmethod
    def _populate_table(widget, df):
        widget.setRowCount(len(df))
        widget.setColumnCount(len(df.columns))
        widget.setHorizontalHeaderLabels(list(df.columns))
        for r, (_, row) in enumerate(df.iterrows()):
            for c, val in enumerate(row):
                if isinstance(val, float):
                    txt = f"{val:.4f}"
                else:
                    txt = str(val)
                item = QTableWidgetItem(txt)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                widget.setItem(r, c, item)

    # ── Export ───────────────────────────────────────────────────────────

    def _export_plot(self):
        if self.current_results is None:
            QMessageBox.warning(self, "No Results", "Run an analysis first.")
            return
        from PyQt5.QtWidgets import QFileDialog
        path, _ = QFileDialog.getSaveFileName(self, "Save Plot", "snp_plot.png",
                                              "PNG Images (*.png)")
        if path:
            self._fig.savefig(path, dpi=150, bbox_inches="tight")
            QMessageBox.information(self, "Saved", f"Plot saved to:\n{path}")

    def _export_csv(self):
        if self.current_results is None:
            QMessageBox.warning(self, "No Results", "Run an analysis first.")
            return
        from PyQt5.QtWidgets import QFileDialog
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", "snp_results.csv",
                                              "CSV Files (*.csv)")
        if path:
            self.current_results["master_df"].to_csv(path, index=False)
            QMessageBox.information(self, "Exported", f"Results exported to:\n{path}")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = SNPAnalyzerApp()
    window.show()
    sys.exit(app.exec_())
