import streamlit as st
import os, shutil, subprocess, time, pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import xml.etree.ElementTree as ET
from Bio import SeqIO, Entrez

Entrez.email = "bmugar50@gmail.com"

def automate_gene_map(accession_id):
    print(f"Fetching annotations for {accession_id}...")
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
    auto_map = {}
    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("product"))[0]
            auto_map[gene_name] = (int(feature.location.start) + 1, int(feature.location.end))
    return auto_map

def automate_sra_discovery(fasta_id, limit=5, retries=3):
    for i in range(retries):
        try:
            search_handle = Entrez.esearch(db="sra", term=f"{fasta_id}[Genome]", retmax=limit)
            uids = Entrez.read(search_handle)['IdList']
            search_handle.close()
            if not uids: return pd.DataFrame()
            sum_handle = Entrez.esummary(db="sra", id=",".join(uids[:limit]))
            summaries = Entrez.read(sum_handle)
            sum_handle.close()
            rows = []
            for s in summaries:
                run_acc = ET.fromstring(f"<root>{s['Runs']}</root>").find(".//Run").get("acc")
                rows.append({"SRR_ID": run_acc, "Title": s.get('Title', 'No Title'), "Date": s.get('CreateDate')})
            return pd.DataFrame(rows)
        except Exception as e:
            print(f"Retry {i+1} due to: {e}"); time.sleep(3)
    return pd.DataFrame()

def align_sra_to_fasta(fasta_path, srr_id):
    bam_sorted = f"{srr_id}_sorted.bam"
    print(f"📥 Downloading and Aligning {srr_id}...")
    subprocess.run(f"fasterq-dump {srr_id}", shell=True, check=True)
    fastq = f"{srr_id}.fastq" if os.path.exists(f"{srr_id}.fastq") else f"{srr_id}_1.fastq"
    subprocess.run(f"minimap2 -a {fasta_path} {fastq} > {srr_id}.sam", shell=True, check=True)
    subprocess.run(f"samtools view -bS {srr_id}.sam | samtools sort -o {bam_sorted}", shell=True, check=True)
    subprocess.run(f"samtools index {bam_sorted}", shell=True, check=True)
    for f in [f"{srr_id}.sam", f"{srr_id}.fastq", f"{srr_id}_1.fastq", f"{srr_id}_2.fastq"]:
        if os.path.exists(f): os.remove(f)
    return bam_sorted

def calculate_snp_frequencies(bam_path, ref_id):
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = []
    for col in bam.pileup(ref_id):
        pos, depth = col.pos + 1, col.nsegments
        if depth < 20: continue
        counts = {'A':0,'T':0,'C':0,'G':0,'Indel':0}
        for read in col.pileups:
            if read.is_del or read.is_refskip or read.indel != 0: counts['Indel'] += 1
            elif not read.is_del:
                base = read.alignment.query_sequence[read.query_position]
                if base in counts: counts[base] += 1
        for var, count in counts.items():
            if count/depth > 0.05:
                results.append({"Position": pos, "Variant_Type": "SNP" if var != 'Indel' else "INDEL",
                                "Base": var, "Count": count, "Depth": depth, "Frequency": round(count/depth, 4)})
    bam.close()
    return pd.DataFrame(results)

def plot_genome_variants(df, gene_map, title):
    plt.figure(figsize=(15, 5))
    colors = sns.color_palette("husl", len(gene_map))
    for i, (gene, (start, end)) in enumerate(gene_map.items()):
        plt.axvspan(start, end, alpha=0.15, color=colors[i])
        plt.text((start+end)/2, 1.05, gene, rotation=45, ha='center', fontsize=8)
    sns.scatterplot(data=df, x='Position', y='Frequency', hue='Variant_Type', style='Sample_ID', alpha=0.6)
    plt.title(title); plt.ylim(-0.05, 1.2); plt.show()

st.title("SNP Frequency Analyzer")

nucleotide_id = st.text_input("Enter Nucleotide Accession ID", value="NC_045512.2")

if 'df_candidates' not in st.session_state:
    st.session_state.df_candidates = pd.DataFrame()

if st.button("Fetch SRR Codes"):
    with st.spinner("Fetching SRR codes..."):
        st.session_state.df_candidates = automate_sra_discovery(nucleotide_id, limit=10)

if not st.session_state.df_candidates.empty:
    st.dataframe(st.session_state.df_candidates)
    selected_srr = st.multiselect("Select SRR IDs to compare", st.session_state.df_candidates['SRR_ID'].tolist())

    if st.button("Run Analysis"):
        if not selected_srr:
            st.error("Please select at least one SRR ID")
        else:
            progress_bar = st.progress(0)
            status_text = st.empty()
            start_time = time.time()

            # 1. Fetch ref
            status_text.text(f"Fetching reference... Elapsed: {time.time() - start_time:.1f}s")
            with Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text") as h:
                with open("ref.fasta", "w") as f: f.write(h.read())
            progress_bar.progress(10)

            # Gene map
            status_text.text(f"Building gene map... Elapsed: {time.time() - start_time:.1f}s")
            GENE_MAP = automate_gene_map(nucleotide_id)
            progress_bar.progress(20)

            # Process samples
            master_list = []
            for i, srr in enumerate(selected_srr):
                status_text.text(f"Processing {srr}... Elapsed: {time.time() - start_time:.1f}s")
                bam = align_sra_to_fasta("ref.fasta", srr)
                df = calculate_snp_frequencies(bam, nucleotide_id)
                df['Sample_ID'] = srr
                df['Gene'] = df['Position'].apply(lambda x: next((g for g, (s, e) in GENE_MAP.items() if s <= x <= e), "Intergenic"))
                master_list.append(df)
                if os.path.exists(bam): os.remove(bam); os.remove(bam+".bai")
                progress_bar.progress(20 + int((i+1)/len(selected_srr)*60))

            status_text.text(f"Analyzing... Elapsed: {time.time() - start_time:.1f}s")
            master_df = pd.concat(master_list)
            progress_bar.progress(80)

            # Plot
            fig, ax = plt.subplots(figsize=(15, 5))
            colors = sns.color_palette("husl", len(GENE_MAP))
            for i, (gene, (start, end)) in enumerate(GENE_MAP.items()):
                ax.axvspan(start, end, alpha=0.15, color=colors[i])
                ax.text((start+end)/2, 1.05, gene, rotation=45, ha='center', fontsize=8)
            sns.scatterplot(data=master_df, x='Position', y='Frequency', hue='Variant_Type', style='Sample_ID', alpha=0.6, ax=ax)
            ax.set_title(f"Mutation Profile for {nucleotide_id}"); ax.set_ylim(-0.05, 1.2)
            st.pyplot(fig)
            progress_bar.progress(90)

            # SNP Count by Gene
            snp_by_gene = master_df.groupby('Gene').size().reset_index(name='Total SNPs').sort_values('Total SNPs', ascending=False)
            st.subheader("SNP Count by Gene")
            st.dataframe(snp_by_gene, use_container_width=True)

            # Table - High Frequency Mutations Summary
            high_freq = master_df[master_df['Frequency'] > 0.8].copy()
            if not high_freq.empty:
                summary_table = high_freq.groupby(['Position', 'Gene', 'Base']).agg({
                    'Frequency': 'mean',
                    'Sample_ID': 'count'
                }).rename(columns={'Frequency': 'Avg Frequency', 'Sample_ID': 'Samples'}).reset_index()
                summary_table = summary_table.sort_values('Avg Frequency', ascending=False).head(20)
                st.subheader("High Frequency Mutations (>80%)")
                st.dataframe(summary_table, use_container_width=True)
            else:
                st.info("No high-frequency mutations (>80%) found")
            progress_bar.progress(100)
            status_text.text(f"Analysis Complete! Total time: {time.time() - start_time:.1f}s")