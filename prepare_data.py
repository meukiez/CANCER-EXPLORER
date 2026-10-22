"""
prepare_data.py
---------------
Downloads TCGA RNA-seq expression and clinical survival data from UCSC Xena,
filters primary tumor (01) vs solid-tissue normal (11) samples, runs differential
expression analysis using PyDESeq2, and exports lightweight parquet files.
"""

import os
import io
import json
import argparse
import subprocess
import requests
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

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

def setup_app_metadata(data_dir="data"):
    """Saves metadata and the gene map required by the Streamlit app."""
    os.makedirs(data_dir, exist_ok=True)
    
    # Save cancers.json
    cancers_json_path = os.path.join(data_dir, "cancers.json")
    with open(cancers_json_path, "w") as f:
        json.dump(CANCERS, f, indent=2)
    print(f"Updated {cancers_json_path}")

    # Save the Ensembl-to-symbol map used to label the prepared expression data.
    gene_map_path = os.path.join(data_dir, "gene_map.parquet")
    if not os.path.exists(gene_map_path):
        probemap_path = os.path.join("raw_downloads", "gencode_v36_probemap.tsv")
        if os.path.exists(probemap_path):
            gene_map = pd.read_csv(probemap_path, sep="\t", usecols=["id", "gene"])
        else:
            response = requests.get(
                f"{UCSC_XENA_HUB}gencode_v36_probemap.tsv", timeout=60
            )
            response.raise_for_status()
            gene_map = pd.read_csv(
                io.BytesIO(response.content), sep="\t", usecols=["id", "gene"]
            )
        gene_map = (gene_map.rename(columns={"id": "Ensembl_ID"})[
            ["Ensembl_ID", "gene"]
        ].dropna())
        if gene_map.empty:
            raise ValueError("The GENCODE probemap did not contain any gene identifiers.")
        gene_map.to_parquet(gene_map_path, index=False)
        print(f"Created {gene_map_path}")

    # Download MSigDB Hallmark Gene Sets for pathway analysis
    gene_sets_path = os.path.join(data_dir, "gene_sets.json")
    if not os.path.exists(gene_sets_path):
        print("Downloading MSigDB Hallmark gene sets for pathway analysis...")
        url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
        try:
            res = requests.get(url, timeout=30)
            if res.status_code == 200:
                gene_sets = {}
                for line in res.text.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) > 2:
                        pathway_name = parts[0]
                        genes = [g.split(",")[0].strip() for g in parts[2:] if g.strip()]
                        gene_sets[pathway_name] = genes
                
                with open(gene_sets_path, "w") as f:
                    json.dump({"MSigDB Hallmark": gene_sets}, f)
                print(f"Successfully created {gene_sets_path}")
        except Exception as e:
            print(f"Warning: Could not download gene sets: {e}")

def download_tcga_data(cancer_code, raw_dir="raw_downloads"):
    """Downloads expression counts and survival metadata from UCSC Xena GDC hub."""
    os.makedirs(raw_dir, exist_ok=True)
    
    counts_file = f"TCGA-{cancer_code}.star_counts.tsv.gz"
    surv_file = f"TCGA-{cancer_code}.survival.tsv.gz"
    
    counts_path = os.path.join(raw_dir, counts_file)
    surv_path = os.path.join(raw_dir, surv_file)
    
    counts_url = f"{UCSC_XENA_HUB}TCGA-{cancer_code}.star_counts.tsv.gz"
    surv_url = f"{UCSC_XENA_HUB}TCGA-{cancer_code}.survival.tsv.gz"
    
    if not os.path.exists(counts_path) or os.path.getsize(counts_path) < 1000:
        print(f"Downloading RNA-seq counts for {cancer_code}...")
        subprocess.run(["curl", "-L", "-o", counts_path, counts_url], check=True)
        
    if not os.path.exists(surv_path) or os.path.getsize(surv_path) < 1000:
        print(f"Downloading survival data for {cancer_code}...")
        subprocess.run(["curl", "-L", "-o", surv_path, surv_url], check=True)

    return counts_path, surv_path

def process_cancer(cancer_code, data_dir="data", max_tumor_samples=50):
    """Processes raw counts into differential expression, expression matrix, and survival files."""
    cancer_out_dir = os.path.join(data_dir, cancer_code.lower())
    os.makedirs(cancer_out_dir, exist_ok=True)
    
    counts_path, surv_path = download_tcga_data(cancer_code)
    
    # Load RNA-seq counts (log2(count+1) converted back to un-logged counts)
    df_counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts = np.rint(2**df_counts - 1).astype(int)
    
    sample_ids = counts.columns
    tumors = [s for s in sample_ids if s.split("-")[3].startswith("01")]
    normals = [s for s in sample_ids if s.split("-")[3].startswith("11")]
    
    if len(normals) < 3:
        print(f"Skipping {cancer_code}: Only {len(normals)} normal sample(s) available.")
        return

    if len(tumors) > max_tumor_samples:
        np.random.seed(42)
        tumors = list(np.random.choice(tumors, size=max_tumor_samples, replace=False))
        
    selected_samples = tumors + normals
    counts_sub = counts[selected_samples].T  # Rows = Samples, Cols = Genes
    
    sample_metadata = pd.DataFrame({
        "condition": ["Tumor" if s in tumors else "Normal" for s in selected_samples]
    }, index=selected_samples)
    
    # Differential Expression Analysis
    print(f"Running PyDESeq2 for {cancer_code} ({len(tumors)} Tumors vs. {len(normals)} Normals)...")
    dds = DeseqDataSet(counts=counts_sub, metadata=sample_metadata, design_factors="condition", refit_cooks=True)
    dds.deseq2()
    
    stat_res = DeseqStats(dds, contrast=["condition", "Tumor", "Normal"])
    stat_res.summary()
    
    df_de = stat_res.results_df.reset_index().rename(columns={"index": "gene"})
    df_de.to_parquet(os.path.join(cancer_out_dir, "de_results.parquet"))
    
    # Save Normalized Expression Matrix (Rows = Samples, Cols = Genes)
    counts_norm = np.log2(counts_sub + 1)
    counts_norm.reset_index().rename(columns={"index": "sample"}).to_parquet(os.path.join(cancer_out_dir, "expression.parquet"))
    
    # Process Survival Data
    df_surv = pd.read_csv(surv_path, sep="\t")
    df_surv = df_surv[df_surv["sample"].isin(selected_samples)]
    df_surv["type"] = df_surv["sample"].apply(lambda s: "Tumor" if s in tumors else "Normal")
    df_surv.to_parquet(os.path.join(cancer_out_dir, "samples.parquet"))

    print(f"Successfully processed {cancer_code}!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TCGA datasets for Cancer Gene Explorer.")
    parser.add_argument("--cancers", nargs="+", default=list(CANCERS.keys()), help="TCGA cancer codes to process.")
    args = parser.parse_args()
    
    setup_app_metadata()
    for c in args.cancers:
        if c.upper() in CANCERS:
            try:
                process_cancer(c.upper())
            except Exception as err:
                print(f"Error processing {c}: {err}")