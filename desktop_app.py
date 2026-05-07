import sys
import os
import time
import threading
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import xml.etree.ElementTree as ET
from Bio import SeqIO, Entrez
import subprocess

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLineEdit, QLabel, QTableWidget, QTableWidgetItem,
    QProgressBar, QCheckBox, QScrollArea, QFrame, QGroupBox,
    QMessageBox, QSpinBox, QFileDialog
)
from PyQt5.QtCore import Qt, pyqtSignal, QThread, QTimer
from PyQt5.QtGui import QFont
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

Entrez.email = "bmugar50@gmail.com"


class SRRFetchWorker(QThread):
    """Worker thread for fetching SRR codes"""
    finished = pyqtSignal(pd.DataFrame)
    error_occurred = pyqtSignal(str)

    def __init__(self, nucleotide_id, limit):
        super().__init__()
        self.nucleotide_id = nucleotide_id
        self.limit = limit

    def run(self):
        try:
            result = self._automate_sra_discovery(self.nucleotide_id, limit=self.limit)
            self.finished.emit(result)
        except Exception as e:
            self.error_occurred.emit(f"Error fetching SRR codes: {str(e)}")

    def _automate_sra_discovery(self, fasta_id, limit=5, retries=3):
        for i in range(retries):
            try:
                search_handle = Entrez.esearch(db="sra", term=f"{fasta_id}[Genome]", retmax=limit)
                uids = Entrez.read(search_handle)['IdList']
                search_handle.close()
                if not uids:
                    return pd.DataFrame()
                sum_handle = Entrez.esummary(db="sra", id=",".join(uids[:limit]))
                summaries = Entrez.read(sum_handle)
                sum_handle.close()
                rows = []
                for s in summaries:
                    run_acc = ET.fromstring(f"<root>{s['Runs']}</root>").find(".//Run").get("acc")
                    rows.append({
                        "SRR_ID": run_acc,
                        "Title": s.get('Title', 'No Title'),
                        "Date": s.get('CreateDate')
                    })
                return pd.DataFrame(rows)
            except Exception as e:
                if i < retries - 1:
                    time.sleep(2)
        return pd.DataFrame()


class AnalysisWorker(QThread):
    """Worker thread for running SNP analysis without blocking UI"""
    progress_updated = pyqtSignal(int)
    status_updated = pyqtSignal(str)
    elapsed_updated = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)

    def __init__(self, nucleotide_id, selected_srr):
        super().__init__()
        self.nucleotide_id = nucleotide_id
        self.selected_srr = selected_srr
        self.start_time = time.time()

    def run(self):
        try:
            # 1. Fetch reference
            self.status_updated.emit("Fetching reference...")
            with Entrez.efetch(db="nucleotide", id=self.nucleotide_id, rettype="fasta", retmode="text") as h:
                with open("ref.fasta", "w") as f:
                    f.write(h.read())
            self.progress_updated.emit(10)
            self._update_elapsed()

            # 2. Build gene map
            self.status_updated.emit("Building gene map...")
            GENE_MAP = self._automate_gene_map(self.nucleotide_id)
            self.progress_updated.emit(20)
            self._update_elapsed()

            # 3. Process samples
            master_list = []
            for i, srr in enumerate(self.selected_srr):
                self.status_updated.emit(f"Processing {srr}...")
                try:
                    bam = self._align_sra_to_fasta("ref.fasta", srr)
                    df = self._calculate_snp_frequencies(bam, self.nucleotide_id)
                    df['Sample_ID'] = srr
                    df['Gene'] = df['Position'].apply(
                        lambda x: next((g for g, (s, e) in GENE_MAP.items() if s <= x <= e), "Intergenic")
                    )
                    master_list.append(df)
                    if os.path.exists(bam):
                        os.remove(bam)
                        if os.path.exists(bam + ".bai"):
                            os.remove(bam + ".bai")
                except Exception as e:
                    self.error_occurred.emit(f"Error processing {srr}: {str(e)}")
                    continue
                self.progress_updated.emit(20 + int((i + 1) / len(self.selected_srr) * 60))
                self._update_elapsed()

            if not master_list:
                self.error_occurred.emit("No samples processed successfully")
                return

            # 4. Analyze
            self.status_updated.emit("Analyzing results...")
            master_df = pd.concat(master_list, ignore_index=True)
            self.progress_updated.emit(80)
            self._update_elapsed()

            # 5. Generate results
            results = {
                'master_df': master_df,
                'gene_map': GENE_MAP,
                'nucleotide_id': self.nucleotide_id
            }

            self.progress_updated.emit(100)
            self.status_updated.emit("Analysis Complete!")
            self.finished.emit(results)

        except Exception as e:
            self.error_occurred.emit(f"Error during analysis: {str(e)}")

    def _automate_gene_map(self, accession_id):
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
        auto_map = {}
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("product"))[0]
                auto_map[gene_name] = (int(feature.location.start) + 1, int(feature.location.end))
        return auto_map

    def _automate_sra_discovery(self, fasta_id, limit=5, retries=3):
        for i in range(retries):
            try:
                search_handle = Entrez.esearch(db="sra", term=f"{fasta_id}[Genome]", retmax=limit)
                uids = Entrez.read(search_handle)['IdList']
                search_handle.close()
                if not uids:
                    return pd.DataFrame()
                sum_handle = Entrez.esummary(db="sra", id=",".join(uids[:limit]))
                summaries = Entrez.read(sum_handle)
                sum_handle.close()
                rows = []
                for s in summaries:
                    run_acc = ET.fromstring(f"<root>{s['Runs']}</root>").find(".//Run").get("acc")
                    rows.append({
                        "SRR_ID": run_acc,
                        "Title": s.get('Title', 'No Title'),
                        "Date": s.get('CreateDate')
                    })
                return pd.DataFrame(rows)
            except Exception as e:
                if i < retries - 1:
                    time.sleep(2)
        return pd.DataFrame()

    def _align_sra_to_fasta(self, fasta_path, srr_id):
        bam_sorted = f"{srr_id}_sorted.bam"
        subprocess.run(f"fasterq-dump {srr_id}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        fastq = f"{srr_id}.fastq" if os.path.exists(f"{srr_id}.fastq") else f"{srr_id}_1.fastq"
        subprocess.run(f"minimap2 -a {fasta_path} {fastq} > {srr_id}.sam", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(f"samtools view -bS {srr_id}.sam | samtools sort -o {bam_sorted}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(f"samtools index {bam_sorted}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for f in [f"{srr_id}.sam", f"{srr_id}.fastq", f"{srr_id}_1.fastq", f"{srr_id}_2.fastq"]:
            if os.path.exists(f):
                os.remove(f)
        return bam_sorted

    def _calculate_snp_frequencies(self, bam_path, ref_id):
        bam = pysam.AlignmentFile(bam_path, "rb")
        results = []
        for col in bam.pileup(ref_id):
            pos, depth = col.pos + 1, col.nsegments
            if depth < 20:
                continue
            counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'Indel': 0}
            for read in col.pileups:
                if read.is_del or read.is_refskip or read.indel != 0:
                    counts['Indel'] += 1
                elif not read.is_del:
                    base = read.alignment.query_sequence[read.query_position]
                    if base in counts:
                        counts[base] += 1
            for var, count in counts.items():
                if count / depth > 0.05:
                    results.append({
                        "Position": pos,
                        "Variant_Type": "SNP" if var != 'Indel' else "INDEL",
                        "Base": var,
                        "Count": count,
                        "Depth": depth,
                        "Frequency": round(count / depth, 4)
                    })
        bam.close()
        return pd.DataFrame(results)

    def _update_elapsed(self):
        elapsed = time.time() - self.start_time
        self.elapsed_updated.emit(f"Elapsed: {elapsed:.1f}s")


class SNPAnalyzerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SNP Frequency Analyzer")
        self.setGeometry(100, 100, 1200, 900)
        self.df_candidates = pd.DataFrame()
        self.current_results = None
        self.worker_thread = None
        self.srr_fetch_worker = None

        self.init_ui()

    def init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Input section
        input_group = QGroupBox("Input Parameters")
        input_layout = QHBoxLayout()
        input_layout.addWidget(QLabel("Nucleotide ID:"))
        self.nucleotide_input = QLineEdit("NC_045512.2")
        input_layout.addWidget(self.nucleotide_input)
        input_layout.addWidget(QLabel("Limit SRR results:"))
        self.limit_spinbox = QSpinBox()
        self.limit_spinbox.setValue(10)
        self.limit_spinbox.setMinimum(1)
        self.limit_spinbox.setMaximum(50)
        input_layout.addWidget(self.limit_spinbox)
        input_group.setLayout(input_layout)
        main_layout.addWidget(input_group)

        # Fetch button
        self.fetch_button = QPushButton("Fetch SRR Codes")
        self.fetch_button.clicked.connect(self.fetch_srr_codes)
        main_layout.addWidget(self.fetch_button)

        # SRR selection section
        srr_group = QGroupBox("Select SRR IDs to Compare")
        srr_layout = QVBoxLayout()
        self.srr_scroll = QScrollArea()
        self.srr_scroll.setWidgetResizable(True)
        self.srr_container = QWidget()
        self.srr_container_layout = QVBoxLayout(self.srr_container)
        self.srr_scroll.setWidget(self.srr_container)
        srr_layout.addWidget(self.srr_scroll)
        srr_group.setLayout(srr_layout)
        main_layout.addWidget(srr_group)

        # Progress section
        progress_group = QGroupBox("Analysis Progress")
        progress_layout = QVBoxLayout()
        self.progress_bar = QProgressBar()
        progress_layout.addWidget(self.progress_bar)
        self.status_label = QLabel("Ready")
        progress_layout.addWidget(self.status_label)
        self.elapsed_label = QLabel("Elapsed: 0.0s")
        progress_layout.addWidget(self.elapsed_label)
        progress_group.setLayout(progress_layout)
        main_layout.addWidget(progress_group)

        # Run button
        self.run_button = QPushButton("Run Analysis")
        self.run_button.clicked.connect(self.run_analysis)
        self.run_button.setEnabled(False)
        main_layout.addWidget(self.run_button)

        # Results section
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout()

        # Plot area
        self.plot_figure = Figure(figsize=(10, 4), dpi=100)
        self.plot_canvas = FigureCanvas(self.plot_figure)
        results_layout.addWidget(self.plot_canvas)

        # SNP by Gene table
        self.snp_gene_table = QTableWidget()
        results_layout.addWidget(QLabel("SNP Count by Gene:"))
        results_layout.addWidget(self.snp_gene_table)

        # High Frequency Mutations table
        self.high_freq_table = QTableWidget()
        results_layout.addWidget(QLabel("High Frequency Mutations (>80%):"))
        results_layout.addWidget(self.high_freq_table)

        results_group.setLayout(results_layout)
        main_layout.addWidget(results_group)

    def fetch_srr_codes(self):
        nucleotide_id = self.nucleotide_input.text().strip()
        if not nucleotide_id:
            QMessageBox.warning(self, "Input Error", "Please enter a nucleotide ID")
            return

        self.fetch_button.setEnabled(False)
        self.status_label.setText("Fetching SRR codes...")

        self.srr_fetch_worker = SRRFetchWorker(nucleotide_id, self.limit_spinbox.value())
        self.srr_fetch_worker.finished.connect(self.on_srr_codes_fetched)
        self.srr_fetch_worker.error_occurred.connect(self.handle_fetch_error)
        self.srr_fetch_worker.start()

    def on_srr_codes_fetched(self, df):
        self.df_candidates = df
        if self.df_candidates.empty:
            self.status_label.setText("No SRR codes found")
            self.run_button.setEnabled(False)
        else:
            self.status_label.setText(f"Found {len(self.df_candidates)} SRR codes")
            self.populate_srr_checkboxes()
        self.fetch_button.setEnabled(True)

    def handle_fetch_error(self, error_msg):
        self.status_label.setText(f"Error: {error_msg}")
        self.fetch_button.setEnabled(True)
        QMessageBox.critical(self, "Fetch Error", error_msg)

    def populate_srr_checkboxes(self):
        # Clear existing checkboxes
        while self.srr_container_layout.count():
            item = self.srr_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        self.srr_checkboxes = []
        if self.df_candidates.empty:
            label = QLabel("No SRR codes found")
            self.srr_container_layout.addWidget(label)
            self.run_button.setEnabled(False)
            return

        for idx, row in self.df_candidates.iterrows():
            try:
                srr_id = str(row['SRR_ID']).strip()
                title = str(row.get('Title', 'No Title')).strip()[:50]
                checkbox = QCheckBox(f"{srr_id} - {title}")
                self.srr_checkboxes.append(checkbox)
                self.srr_container_layout.addWidget(checkbox)
            except Exception as e:
                print(f"Error creating checkbox: {e}")

        self.srr_container_layout.addStretch()
        self.run_button.setEnabled(True)
        self.srr_scroll.setVisible(True)

    def run_analysis(self):
        selected_srr = [self.srr_checkboxes[i].text().split(' - ')[0] 
                       for i, cb in enumerate(self.srr_checkboxes) if cb.isChecked()]

        if not selected_srr:
            QMessageBox.warning(self, "Selection Error", "Please select at least one SRR ID")
            return

        nucleotide_id = self.nucleotide_input.text().strip()
        self.run_button.setEnabled(False)
        self.progress_bar.setValue(0)

        self.worker_thread = AnalysisWorker(nucleotide_id, selected_srr)
        self.worker_thread.progress_updated.connect(self.progress_bar.setValue)
        self.worker_thread.status_updated.connect(self.status_label.setText)
        self.worker_thread.elapsed_updated.connect(self.elapsed_label.setText)
        self.worker_thread.finished.connect(self.display_results)
        self.worker_thread.error_occurred.connect(self.handle_error)
        self.worker_thread.start()

    def display_results(self, results):
        master_df = results['master_df']
        gene_map = results['gene_map']
        nucleotide_id = results['nucleotide_id']

        # Plot
        self.plot_figure.clear()
        ax = self.plot_figure.add_subplot(111)
        colors = sns.color_palette("husl", len(gene_map))
        for i, (gene, (start, end)) in enumerate(gene_map.items()):
            ax.axvspan(start, end, alpha=0.15, color=colors[i])
            ax.text((start + end) / 2, 1.05, gene, rotation=45, ha='center', fontsize=8)
        sns.scatterplot(data=master_df, x='Position', y='Frequency', hue='Variant_Type', 
                       style='Sample_ID', alpha=0.6, ax=ax)
        ax.set_title(f"Mutation Profile for {nucleotide_id}")
        ax.set_ylim(-0.05, 1.2)
        self.plot_canvas.draw()

        # SNP by Gene table
        snp_by_gene = master_df.groupby('Gene').size().reset_index(name='Total SNPs').sort_values('Total SNPs', ascending=False)
        self.populate_table(self.snp_gene_table, snp_by_gene)

        # High Frequency Mutations table
        high_freq = master_df[master_df['Frequency'] > 0.8].copy()
        if not high_freq.empty:
            summary_table = high_freq.groupby(['Position', 'Gene', 'Base']).agg({
                'Frequency': 'mean',
                'Sample_ID': 'count'
            }).rename(columns={'Frequency': 'Avg Frequency', 'Sample_ID': 'Samples'}).reset_index()
            summary_table = summary_table.sort_values('Avg Frequency', ascending=False).head(20)
            self.populate_table(self.high_freq_table, summary_table)
        else:
            self.high_freq_table.setRowCount(0)
            self.high_freq_table.setColumnCount(0)

        self.run_button.setEnabled(True)
        self.current_results = results

    def populate_table(self, table_widget, df):
        table_widget.setRowCount(len(df))
        table_widget.setColumnCount(len(df.columns))
        table_widget.setHorizontalHeaderLabels(list(df.columns))

        for row_idx, (idx, row) in enumerate(df.iterrows()):
            for col_idx, value in enumerate(row):
                try:
                    # Format numeric values nicely
                    if isinstance(value, float):
                        display_value = f"{value:.4f}"
                    else:
                        display_value = str(value)
                    item = QTableWidgetItem(display_value)
                    table_widget.setItem(row_idx, col_idx, item)
                except Exception as e:
                    print(f"Error setting table item: {e}")

        table_widget.resizeColumnsToContents()

    def handle_error(self, error_msg):
        self.status_label.setText(f"Error: {error_msg}")
        self.run_button.setEnabled(True)
        QMessageBox.critical(self, "Analysis Error", error_msg)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SNPAnalyzerApp()
    window.show()
    sys.exit(app.exec_())
