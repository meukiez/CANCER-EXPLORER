# Cancer Gene Explorer 🧬

An interactive web app for exploring how gene activity differs between **tumors and normal
tissue**, using real patient data from The Cancer Genome Atlas (TCGA).

Pick a cancer type and you can:

- **Look up any gene**: tumor vs. normal box plot, survival curve, and a plain-language explanation
- **See every gene at once** in an interactive volcano plot and heatmap
- **Find pathways**: which biological processes are switched on or off in the tumor
- **Compare cancers**: which changes are shared, and which are unique to one cancer
- **Learn**: a built-in guide to the biology and statistics, written for beginners

Included cancers: Breast (BRCA), Kidney clear cell (KIRC), Lung adenocarcinoma (LUAD).

> Screenshot: add one here once the app is running, e.g. `![app](screenshot.png)`

---

## Setup (about 30-60 minutes the first time)

### 1. Install Python and the packages
Install Python 3.10 or newer, then in a terminal inside this folder:

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# Mac/Linux: source .venv/bin/activate
pip install -r requirements.txt
pip install pydeseq2          # only needed for step 2
```

### 2. Download and prepare the data (once)
```bash
python prepare_data.py --cancers KIRC     # start with one cancer to test (~5-10 min)
python prepare_data.py                    # then all three
```
This downloads real TCGA data (about 1 GB in total) into `raw_downloads/`, runs the
differential-expression analysis, and saves small app-ready files into `data/`.

### 3. Run the app
```bash
streamlit run app.py
```
It opens in your browser at http://localhost:8501.

---

## Put it online for free (Streamlit Community Cloud)

1. Create a GitHub repository and upload this folder **including `data/`**
   (but not `raw_downloads/`, which `.gitignore` already skips).
2. Go to https://share.streamlit.io, sign in with GitHub, click **New app**, pick the repo
   and `app.py`.
3. You get a public link like `https://your-name-cancer-explorer.streamlit.app` to use in
   your application and portfolio.

---

## How it works

| Step | Method | File |
|---|---|---|
| Data | TCGA STAR gene counts + overall survival, UCSC Xena GDC hub | `prepare_data.py` |
| Samples | Primary tumors (barcode `01`, up to 200 random) vs. solid-tissue normals (`11`) | `prepare_data.py` |
| Differential expression | PyDESeq2 (negative binomial model, Benjamini-Hochberg correction) | `prepare_data.py` |
| Pathways | Over-representation (hypergeometric) test on MSigDB Hallmark, KEGG, Reactome and GO sets from Enrichr | `analysis.py` |
| Survival | Kaplan-Meier curves + log-rank test, split at median expression | `analysis.py` |
| Overview | PCA on the 1,000 most variable genes | `analysis.py` |

## Ideas to make it your own

- Add another cancer: add its TCGA code to `CANCERS` in `prepare_data.py` (e.g. `"COAD"`,
  `"PRAD"`, `"LIHC"`) and add example genes in `genes_info.py`
- Write notes for more genes in `genes_info.py`
- Compare breast cancer **subtypes** (HER2+, triple-negative) instead of tumor vs. normal
- Add a "gene of the week" story page with your own findings

## Limitations

TCGA normal samples are tissue next to the tumor, not from healthy people. Gene changes show
association, not cause. Survival plots don't adjust for age, stage or treatment. **This app is
for education only, not medical advice.**

## Credits

Built by Charlene Siawira. Data from the TCGA Research Network (https://www.cancer.gov/tcga),
accessed via UCSC Xena (Goldman et al., *Nature Biotechnology*, 2020). Gene sets from Enrichr and
MSigDB. Analysis with PyDESeq2 (Muzellec et al., *Bioinformatics*, 2023). Inspired by GEPIA2
and UALCAN.
