"""
prepare_data.py
---------------
Downloads TCGA RNA-seq expression and clinical survival data from UCSC Xena,
filters primary tumor (01) vs solid-tissue normal (11) samples, runs differential
expression analysis using PyDESeq2, and exports lightweight parquet files.
"""

import os
import argparse
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Dictionary of available TCGA Cancers (Expandable list)
CANCERS = {
    "BRCA": "Breast Invasive Carcinoma",
    "KIRC": "Kidney Renal Clear Cell Carcinoma",
    "LUAD": "Lung Adenocarcinoma",
    "COAD": "Colon Adenocarcinoma",
    "PRAD": "Prostate Adenocarcinoma",
    "LIHC": "Liver Hepatocellular Carcinoma",
    "SKCM": "Skin Cutaneous Melanoma",
    "UCEC": "Uterine Corpus Endometrial Carcinoma",
    "THCA": "Thyroid Carcinoma"
}

UCSC_XENA_HUB = "https://gdc.xenahubs.net/download/"

import os
import subprocess

def download_tcga_data(cancer_code, raw_dir="raw_downloads"):
    """Downloads expression counts and survival metadata using the correct TCGA star_counts naming."""
    os.makedirs(raw_dir, exist_ok=True)
    
    counts_file = f"TCGA-{cancer_code}.star_counts.tsv.gz"
    surv_file = f"TCGA-{cancer_code}.survival.tsv.gz"
    
    counts_path = os.path.join(raw_dir, counts_file)
    surv_path = os.path.join(raw_dir, surv_file)
    
    counts_url = f"https://gdc.xenahubs.net/download/TCGA-{cancer_code}.star_counts.tsv.gz"
    surv_url = f"https://gdc.xenahubs.net/download/TCGA-{cancer_code}.survival.tsv.gz"
    
    if not os.path.exists(counts_path) or os.path.getsize(counts_path) < 1000:
        print(f"Downloading RNA-seq counts for {cancer_code}...")
        subprocess.run(["curl", "-L", "-o", counts_path, counts_url], check=True)
        
    if not os.path.exists(surv_path) or os.path.getsize(surv_path) < 1000:
        print(f"Downloading survival data for {cancer_code}...")
        subprocess.run(["curl", "-L", "-o", surv_path, surv_url], check=True)

    return counts_path, surv_path

def process_cancer(cancer_code, data_dir="data", max_tumor_samples=200):
    """Processes raw counts into differential expression and survival parquet files."""
    cancer_out_dir = os.path.join(data_dir, cancer_code.lower())
    os.makedirs(cancer_out_dir, exist_ok=True)
    
    counts_path, surv_path = download_tcga_data(cancer_code)
    
    # Load RNA-seq counts (Log2(count+1) converted back to raw counts)
    df_counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts = np.rint(2**df_counts - 1).astype(int)
    
    # Classify TCGA samples: 01 = Primary Tumor, 11 = Solid Tissue Normal
    sample_ids = counts.columns
    tumors = [s for s in sample_ids if s.split("-")[3].startswith("01")]
    normals = [s for s in sample_ids if s.split("-")[3].startswith("11")]
    
    if len(normals) < 5:
        print(f"Warning: {cancer_code} only has {len(normals)} normal samples. DE analysis might be underpowered.")
    
    # Subsample tumor samples to limit memory usage while keeping all available normals
    if len(tumors) > max_tumor_samples:
        np.random.seed(42)
        tumors = list(np.random.choice(tumors, size=max_tumor_samples, replace=False))
        
    selected_samples = tumors + normals
    counts = counts[selected_samples].T # PyDESeq2 expects samples in rows and genes in columns
    
    # Metadata dataframe for PyDESeq2
    sample_metadata = pd.DataFrame({
        "condition": ["Tumor" if s in tumors else "Normal" for s in selected_samples]
    }, index=selected_samples)
    
    # Differential Expression Analysis with PyDESeq2
    print(f"Running PyDESeq2 for {cancer_code} ({len(tumors)} Tumors vs. {len(normals)} Normals)...")
    dds = DeseqDataSet(counts=counts, metadata=sample_metadata, design_factors="condition", refit_cooks=True)
    dds.deseq2()
    
    stat_res = DeseqStats(dds, contrast=["condition", "Tumor", "Normal"])
    stat_res.summary()
    
    df_de = stat_res.results_df.reset_index().rename(columns={"index": "gene"})
    df_de.to_parquet(os.path.join(cancer_out_dir, "de_results.parquet"))
    
    # Prepare Expression Data Matrix for App Visualization
    counts_norm = np.log2(counts + 1) # Normalization for plotting
    counts_norm.to_parquet(os.path.join(cancer_out_dir, "expression.parquet"))
    
    # Process Survival Data
    df_surv = pd.read_csv(surv_path, sep="\t")
    df_surv = df_surv[df_surv["sample"].isin(selected_samples)]
    df_surv["type"] = df_surv["sample"].apply(lambda s: "Tumor" if s in tumors else "Normal")
    df_surv.rename(columns={"OS": "OS", "OS.time": "OS.time"}, inplace=True)
    df_surv.to_parquet(os.path.join(cancer_out_dir, "samples.parquet"))

    print(f"Successfully processed {cancer_code}!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TCGA datasets for Cancer Gene Explorer.")
    parser.add_argument("--cancers", nargs="+", default=list(CANCERS.keys()), help="TCGA cancer codes to process.")
    args = parser.parse_args()
    
    for c in args.cancers:
        if c.upper() in CANCERS:
            process_cancer(c.upper())
        else:
            print(f"Skipping unknown cancer code: {c}")
