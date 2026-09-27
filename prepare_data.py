"""
prepare_data.py - run this ONCE on your own computer before starting the app.

What it does, for each cancer type:
  1. Downloads real TCGA RNA-seq counts + survival data from the UCSC Xena GDC hub
  2. Keeps primary tumors (barcode code 01) and solid-tissue normals (code 11)
  3. Runs differential expression (tumor vs normal) with PyDESeq2
  4. Saves small, app-ready files into the data/ folder

It also downloads pathway gene-set libraries from Enrichr so the app can do
pathway analysis offline.

Usage:
    python prepare_data.py                  # all cancers in CANCERS
    python prepare_data.py --cancers KIRC   # just one (fastest, good for a first test)

Expect roughly 5-20 minutes per cancer type depending on your computer.
"""

import argparse
import io
import json
import os

import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------- settings
CANCERS = {
    "BRCA": "Breast invasive carcinoma",
    "KIRC": "Kidney renal clear cell carcinoma",
    "LUAD": "Lung adenocarcinoma",
}
MAX_TUMORS = 200        # random subset of tumors keeps the app fast and small
RANDOM_SEED = 42
MIN_COUNT = 10          # a gene must have >= MIN_COUNT reads...
MIN_FRACTION = 0.2      # ...in at least this fraction of samples to be kept

XENA = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP_URL = f"{XENA}/gencode.v36.annotation.gtf.gene.probemap"
ENRICHR_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"
GENE_SET_LIBRARIES = {
    "Hallmark (MSigDB)": "MSigDB_Hallmark_2020",
    "KEGG": "KEGG_2021_Human",
    "Reactome": "Reactome_2022",
    "GO Biological Process": "GO_Biological_Process_2023",
}

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RAW = os.path.join(HERE, "raw_downloads")   # big files; NOT pushed to GitHub


# ---------------------------------------------------------------- helpers
def download(url, dest):
    """Download a file once and reuse it on later runs."""
    if os.path.exists(dest):
        print(f"   (already downloaded) {os.path.basename(dest)}")
        return dest
    print(f"   downloading {url}")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = dest + ".part"
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
        os.replace(tmp, dest)
    return dest


def sample_type(barcode):
    """TCGA barcode characters 14-15 give the sample type: 01 = primary tumor, 11 = normal."""
    code = barcode[13:15]
    return {"01": "Tumor", "11": "Normal"}.get(code)


def load_gene_symbols():
    path = download(PROBEMAP_URL, os.path.join(RAW, "gencode_v36_probemap.tsv"))
    pm = pd.read_csv(path, sep="\t")
    return dict(zip(pm.iloc[:, 0], pm.iloc[:, 1]))   # Ensembl ID -> gene symbol


def choose_samples(barcodes):
    """Keep primary tumors + normals; randomly subsample tumors to MAX_TUMORS."""
    types = pd.Series({s: sample_type(s) for s in barcodes}).dropna()
    tumors = types.index[types == "Tumor"]
    normals = types.index[types == "Normal"]
    rng = np.random.default_rng(RANDOM_SEED)
    if len(tumors) > MAX_TUMORS:
        tumors = pd.Index(sorted(rng.choice(tumors, MAX_TUMORS, replace=False)))
    keep = tumors.append(normals)
    return pd.DataFrame({"type": types.loc[keep]})


def load_counts(code, id_to_symbol):
    path = download(f"{XENA}/TCGA-{code}.star_counts.tsv.gz",
                    os.path.join(RAW, f"TCGA-{code}.star_counts.tsv.gz"))
    print("   reading counts (this can take a minute)...")
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    samples = choose_samples(header[1:])
    # only read the columns we need - saves a lot of memory
    df = pd.read_csv(path, sep="\t", index_col=0,
                     usecols=[header[0], *samples.index], dtype={s: np.float32 for s in samples.index})
    # Xena stores log2(count + 1); undo that to get raw integer counts for DESeq2
    counts = np.round(np.power(2.0, df.values) - 1).clip(min=0)
    counts = pd.DataFrame(counts, index=df.index, columns=df.columns)

    # Ensembl IDs -> gene symbols; if a symbol appears twice, keep the more expressed one
    counts["symbol"] = counts.index.map(id_to_symbol)
    counts = counts.dropna(subset=["symbol"])
    counts = counts.loc[counts.drop(columns="symbol").mean(axis=1)
                        .groupby(counts["symbol"]).idxmax()]
    counts = counts.set_index("symbol")
    counts.index.name = "gene"
    return counts.T.loc[samples.index], samples   # rows = samples, columns = genes


def load_survival(code):
    try:
        path = download(f"{XENA}/TCGA-{code}.survival.tsv.gz",
                        os.path.join(RAW, f"TCGA-{code}.survival.tsv.gz"))
        s = pd.read_csv(path, sep="\t")
    except Exception as err:   # the app still works without survival data
        print(f"   WARNING: could not load survival data ({err})")
        return pd.DataFrame(columns=["OS", "OS.time"])
    s = s.rename(columns={s.columns[0]: "sample"})
    return s[["sample", "OS", "OS.time"]].drop_duplicates("sample").set_index("sample")


def run_deseq2(counts, samples):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    meta = samples[["type"]].rename(columns={"type": "condition"})
    dds = DeseqDataSet(counts=counts.astype(int), metadata=meta,
                       design="~condition", refit_cooks=True, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "Tumor", "Normal"], quiet=True)
    stats.summary()
    res = stats.results_df.rename(columns={"log2FoldChange": "log2FC"})
    res.index.name = "gene"

    # normalized expression for plotting: log2(DESeq2-normalized counts + 1)
    size_factors = dds.obs["size_factors"] if "size_factors" in dds.obs else dds.obsm["size_factors"]
    norm = counts.div(np.asarray(size_factors), axis=0)
    return res, np.log2(norm + 1)


# ---------------------------------------------------------------- main steps
def prepare_cancer(code, id_to_symbol):
    print(f"\n=== {code}: {CANCERS[code]} ===")
    counts, samples = load_counts(code, id_to_symbol)
    n = samples["type"].value_counts()
    print(f"   {n.get('Tumor', 0)} tumors, {n.get('Normal', 0)} normals")

    # drop genes that are barely expressed (too noisy to test)
    expressed = (counts >= MIN_COUNT).mean(axis=0) >= MIN_FRACTION
    counts = counts.loc[:, expressed]
    print(f"   {counts.shape[1]} expressed genes kept")

    surv = load_survival(code)
    samples = samples.join(surv, how="left")

    print("   running DESeq2 (the slow part)...")
    res, logexpr = run_deseq2(counts, samples)

    out = os.path.join(DATA, code)
    os.makedirs(out, exist_ok=True)
    res.reset_index().to_parquet(os.path.join(out, "de_results.parquet"), index=False)
    logexpr.astype(np.float16).reset_index(names="sample") \
        .to_parquet(os.path.join(out, "expression.parquet"), index=False)
    samples.reset_index(names="sample").to_parquet(os.path.join(out, "samples.parquet"), index=False)
    print(f"   saved to data/{code}/")


def prepare_gene_sets():
    print("\n=== Pathway gene sets (Enrichr) ===")
    libraries = {}
    for label, name in GENE_SET_LIBRARIES.items():
        print(f"   {name}")
        txt = requests.get(ENRICHR_URL.format(name), timeout=120).text
        sets = {}
        for line in txt.splitlines():
            parts = line.split("\t")
            if len(parts) > 2:
                # Enrichr lines: name <tab> description <tab> gene1 <tab> gene2 ...
                genes = sorted({g.split(",")[0].strip() for g in parts[2:] if g.strip()})
                sets[parts[0]] = genes
        libraries[label] = sets
    with open(os.path.join(DATA, "gene_sets.json"), "w") as f:
        json.dump(libraries, f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancers", nargs="+", default=list(CANCERS), choices=list(CANCERS))
    ap.add_argument("--skip-gene-sets", action="store_true")
    args = ap.parse_args()

    os.makedirs(DATA, exist_ok=True)
    os.makedirs(RAW, exist_ok=True)
    id_to_symbol = load_gene_symbols()
    for code in args.cancers:
        prepare_cancer(code, id_to_symbol)
    if not args.skip_gene_sets:
        prepare_gene_sets()

    with open(os.path.join(DATA, "cancers.json"), "w") as f:
        done = {c: n for c, n in CANCERS.items() if os.path.isdir(os.path.join(DATA, c))}
        json.dump(done, f, indent=2)
    print("\nAll done! Now run:  streamlit run app.py")


if __name__ == "__main__":
    main()
