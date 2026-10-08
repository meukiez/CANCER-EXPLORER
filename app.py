"""
Cancer Gene Explorer - an interactive tool for exploring how gene activity
differs between tumors and normal tissue, using real TCGA data.

Run with:  streamlit run app.py
"""

import json
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import analysis as an
from genes_info import EXAMPLES, GENE_NOTES

st.set_page_config(page_title="Cancer Gene Explorer", page_icon="🧬", layout="wide")

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
COLORS = {"Up in tumor": "#d1495b", "Down in tumor": "#2e86ab", "Not significant": "#c8c8c8",
          "Tumor": "#d1495b", "Normal": "#2e86ab"}


# ---------------------------------------------------------------- data loading
@st.cache_data
def load_cancers():
    path = os.path.join(DATA, "cancers.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


@st.cache_data
def load_de(code):
    de = pd.read_parquet(os.path.join(DATA, code, "de_results.parquet")).dropna(subset=["padj"])
    de["neg_log10_padj"] = -np.log10(de["padj"].clip(lower=1e-300))
    return de.set_index("gene")


@st.cache_resource
def load_expression(code):
    e = pd.read_parquet(os.path.join(DATA, code, "expression.parquet")).set_index("sample")
    return e.astype(np.float32)


@st.cache_data
def load_samples(code):
    return pd.read_parquet(os.path.join(DATA, code, "samples.parquet")).set_index("sample")


@st.cache_data
def load_gene_sets():
    path = os.path.join(DATA, "gene_sets.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


@st.cache_data
def cached_pca(code):
    expr, samples = load_expression(code), load_samples(code)
    return an.pca(expr, samples.loc[expr.index, "type"])


cancers = load_cancers()
if not cancers:
    st.title("Cancer Gene Explorer")
    st.error("No data found yet. Run `python prepare_data.py` first (see README.md).")
    st.stop()


# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.title("🧬 Cancer Gene Explorer")
    code = st.selectbox("Cancer type", list(cancers), format_func=lambda c: f"{cancers[c]} ({c})")
    st.markdown("**What counts as a 'changed' gene?**")
    padj_cut = st.select_slider("Significance cutoff (adjusted p-value)",
                                options=[0.001, 0.01, 0.05, 0.1], value=0.05,
                                help="Smaller = stricter. 0.05 means we accept about a 5% "
                                     "false-discovery rate among the genes we call significant.")
    lfc_cut = st.slider("Minimum change (|log2 fold change|)", 0.0, 4.0, 1.0, 0.5,
                        help="1 means at least 2x higher or lower in tumors; 2 means 4x.")
    st.divider()
    st.caption("Data: The Cancer Genome Atlas (TCGA), via the UCSC Xena GDC hub. "
               "For learning and exploration only, not medical advice.")

de = load_de(code)
de["status"] = an.classify(de, padj_cut, lfc_cut)
samples = load_samples(code)
n_up = int((de["status"] == "Up in tumor").sum())
n_down = int((de["status"] == "Down in tumor").sum())

if "gene" not in st.session_state:
    st.session_state.gene = next((g for g in EXAMPLES.get(code, []) if g in de.index), de.index[0])
if st.session_state.get("pending_gene"):          # set by clicking a dot in the volcano plot
    st.session_state.gene = st.session_state.pop("pending_gene")


def set_gene(g):
    st.session_state.gene = g

tabs = st.tabs(["Overview", "Gene lookup", "Volcano plot", "Pathways", "Compare cancers", "Learn"])


# ---------------------------------------------------------------- Overview
with tabs[0]:
    st.header(cancers[code])
    c = st.columns(4)
    c[0].metric("Tumor samples", int((samples["type"] == "Tumor").sum()))
    c[1].metric("Normal samples", int((samples["type"] == "Normal").sum()))
    c[2].metric("Genes up in tumor", n_up)
    c[3].metric("Genes down in tumor", n_down)

    left, right = st.columns([3, 2])
    with left:
        st.subheader("Do tumors and normal tissue look different overall?")
        pcs, var = cached_pca(code)
        fig = px.scatter(pcs, x="PC1", y="PC2", color="type", color_discrete_map=COLORS,
                         hover_name=pcs.index, opacity=0.75,
                         labels={"PC1": f"PC1 ({var[0]:.0f}% of variation)",
                                 "PC2": f"PC2 ({var[1]:.0f}% of variation)"})
        fig.update_layout(height=430, legend_title=None, margin=dict(t=10))
        st.plotly_chart(fig, width="stretch")
        with st.expander("What am I looking at?"):
            st.markdown(
                "Each dot is one person's sample. **PCA** squeezes thousands of genes into two "
                "axes that capture the biggest differences between samples. If tumor and normal "
                "dots form separate clouds, their overall gene activity is very different.")
    with right:
        st.subheader("Biggest changes")
        sig = de[de["status"] != "Not significant"]
        show = ["log2FC", "padj"]
        st.markdown("**Most increased in tumors**")
        st.dataframe(sig.nlargest(8, "log2FC")[show].style.format({"log2FC": "{:.2f}", "padj": "{:.1e}"}),
                     width="stretch")
        st.markdown("**Most decreased in tumors**")
        st.dataframe(sig.nsmallest(8, "log2FC")[show].style.format({"log2FC": "{:.2f}", "padj": "{:.1e}"}),
                     width="stretch")


# ---------------------------------------------------------------- Gene lookup
with tabs[1]:
    st.header("Look up a gene")
    examples = [g for g in EXAMPLES.get(code, []) if g in de.index]
    if examples:
        st.caption("Try one of these:")
        cols = st.columns(len(examples))
        for col, g in zip(cols, examples):
            col.button(g, key=f"ex_{g}", width="stretch", on_click=set_gene, args=(g,))

    genes = sorted(de.index)
    if st.session_state.gene not in genes:
        st.session_state.gene = genes[0]
    gene = st.selectbox("Search for any gene (type to search)", genes, key="gene")

    row = de.loc[gene]
    fold = 2 ** row["log2FC"]
    direction = "higher" if fold >= 1 else "lower"
    times = fold if fold >= 1 else 1 / fold

    info, links = st.columns([3, 1])
    with info:
        st.subheader(gene)
        if gene in GENE_NOTES:
            st.info(GENE_NOTES[gene])
        verdict = {"Up in tumor": "significantly **higher** in tumors",
                   "Down in tumor": "significantly **lower** in tumors",
                   "Not significant": "**not significantly changed** at your current cutoffs"}[row["status"]]
        st.markdown(f"On average this gene is **{times:.1f}x {direction}** in tumors "
                    f"(log2FC = {row['log2FC']:.2f}, adjusted p = {row['padj']:.1e}). "
                    f"That makes it {verdict}.")
    with links:
        st.markdown("**Learn more**")
        st.markdown(f"- [GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene})\n"
                    f"- [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/?term={gene}%5Bsym%5D+AND+human%5Borgn%5D)\n"
                    f"- [Human Protein Atlas](https://www.proteinatlas.org/search/{gene})")

    expr = load_expression(code)
    left, right = st.columns(2)
    with left:
        st.markdown("**Tumor vs. normal**")
        df = pd.DataFrame({"expression": expr[gene],
                           "type": samples.loc[expr.index, "type"]})
        fig = px.box(df, x="type", y="expression", color="type", points="all",
                     color_discrete_map=COLORS, category_orders={"type": ["Normal", "Tumor"]},
                     labels={"expression": "Expression (log2 normalized counts)", "type": ""})
        fig.update_traces(marker=dict(size=4, opacity=0.5), jitter=0.4)
        fig.update_layout(showlegend=False, height=400, margin=dict(t=10))
        st.plotly_chart(fig, width="stretch")

    with right:
        st.markdown("**Does this gene relate to survival?** (tumor patients only)")
        tum = samples[(samples["type"] == "Tumor")].dropna(subset=["OS", "OS.time"])
        tum = tum[tum.index.isin(expr.index)]
        if len(tum) < 20:
            st.warning("Not enough survival data for this cancer type.")
        else:
            level = expr.loc[tum.index, gene]
            high = level > level.median()
            fig = go.Figure()
            groups = {"High expression": tum[high], "Low expression": tum[~high]}
            for (label, g), color in zip(groups.items(), ["#d1495b", "#2e86ab"]):
                t, s = an.kaplan_meier(g["OS.time"].values / 365.25, g["OS"].values)
                fig.add_trace(go.Scatter(x=t, y=s, mode="lines", line_shape="hv",
                                         name=f"{label} (n={len(g)})", line=dict(color=color, width=2.5)))
            p = an.logrank_p(tum[high]["OS.time"].values, tum[high]["OS"].values,
                             tum[~high]["OS.time"].values, tum[~high]["OS"].values)
            fig.update_layout(height=400, margin=dict(t=10), yaxis_range=[0, 1.02],
                              xaxis_title="Years after diagnosis", yaxis_title="Fraction still alive",
                              legend=dict(orientation="h", y=-0.25))
            st.plotly_chart(fig, width="stretch")
            st.caption(f"Log-rank p = {p:.3g}. Patients are split at the median expression. "
                       + ("The curves differ more than expected by chance."
                          if p < 0.05 else "No clear difference between the groups."))

    # the same gene across every cancer in the app
    across = []
    for c2 in cancers:
        d2 = load_de(c2)
        if gene in d2.index:
            across.append({"cancer": c2, "log2FC": d2.loc[gene, "log2FC"], "padj": d2.loc[gene, "padj"]})
    if len(across) > 1:
        st.markdown(f"**{gene} across all cancer types in this app**")
        a = pd.DataFrame(across)
        a["status"] = an.classify(a, padj_cut, lfc_cut)
        fig = px.bar(a, x="cancer", y="log2FC", color="status", color_discrete_map=COLORS,
                     labels={"log2FC": "log2 fold change (tumor vs normal)", "cancer": ""})
        fig.add_hline(y=0, line_color="gray")
        fig.update_layout(height=300, margin=dict(t=10), legend_title=None)
        st.plotly_chart(fig, width="stretch")


# ---------------------------------------------------------------- Volcano
with tabs[2]:
    st.header("Volcano plot: every gene at once")
    st.markdown(f"Right side = **higher in tumors**, left side = **lower in tumors**, "
                f"higher up = **more confident**. At your cutoffs: "
                f"**{n_up}** genes up, **{n_down}** down. Click a dot to select that gene.")
    n_labels = st.slider("How many top genes to label?", 0, 30, 10)
    plot_df = de.reset_index()
    fig = px.scatter(plot_df, x="log2FC", y="neg_log10_padj", color="status",
                     color_discrete_map=COLORS, hover_name="gene",
                     hover_data={"log2FC": ":.2f", "padj": ":.1e", "neg_log10_padj": False, "status": False},
                     labels={"log2FC": "log2 fold change (tumor vs normal)",
                             "neg_log10_padj": "-log10(adjusted p-value)"},
                     category_orders={"status": ["Up in tumor", "Down in tumor", "Not significant"]})
    fig.update_traces(marker=dict(size=5, opacity=0.7))
    fig.add_vline(x=lfc_cut, line_dash="dot", line_color="gray")
    fig.add_vline(x=-lfc_cut, line_dash="dot", line_color="gray")
    fig.add_hline(y=-np.log10(padj_cut), line_dash="dot", line_color="gray")
    sig = plot_df[plot_df["status"] != "Not significant"].copy()
    sig["score"] = sig["neg_log10_padj"] * sig["log2FC"].abs()
    for _, r in sig.nlargest(n_labels, "score").iterrows():
        fig.add_annotation(x=r["log2FC"], y=r["neg_log10_padj"], text=r["gene"],
                           showarrow=True, arrowhead=0, ax=0, ay=-18, font=dict(size=11))
    fig.update_layout(height=560, legend_title=None, margin=dict(t=10))
    event = st.plotly_chart(fig, width="stretch", on_select="rerun",
                            selection_mode="points", key="volcano")
    points = event.selection.points if event and event.selection else []
    if points:
        pt = points[0]
        picked = fig.data[pt["curve_number"]].hovertext[pt["point_index"]]
        st.success(f"Selected **{picked}**. Open the *Gene lookup* tab to see its details.")
        if st.session_state.get("last_pick") != picked:   # only react to NEW clicks
            st.session_state.last_pick = picked
            st.session_state.pending_gene = picked
            st.rerun()

    st.subheader("Heatmap of the top changed genes")
    n_heat = st.slider("Genes per direction", 5, 30, 15)
    top = list(sig.nlargest(n_heat, "log2FC")["gene"]) + list(sig.nsmallest(n_heat, "log2FC")["gene"])
    if top:
        expr = load_expression(code)
        order = samples.loc[expr.index].sort_values("type").index
        z = expr.loc[order, top]
        z = ((z - z.mean()) / z.std()).clip(-3, 3).T
        fig = px.imshow(z, aspect="auto", color_continuous_scale="RdBu_r", zmin=-3, zmax=3,
                        labels=dict(color="z-score"))
        n_norm = int((samples.loc[order, "type"] == "Normal").sum())
        fig.add_vline(x=n_norm - 0.5, line_color="black", line_width=2)
        fig.update_xaxes(showticklabels=False, title=f"← Normal ({n_norm})  |  Tumor ({len(order) - n_norm}) →")
        fig.update_layout(height=max(350, 16 * len(top)), margin=dict(t=10))
        st.plotly_chart(fig, width="stretch")
        st.caption("Each row is a gene, each column a sample. Red = higher than that gene's "
                   "average, blue = lower.")

    st.subheader("Full results table")
    table = de[["baseMean", "log2FC", "pvalue", "padj", "status"]].sort_values("padj")
    only_sig = st.checkbox("Show only significant genes", value=True)
    if only_sig:
        table = table[table["status"] != "Not significant"]
    st.dataframe(table, width="stretch", height=300)
    st.download_button("Download table (CSV)", table.to_csv().encode(),
                       file_name=f"{code}_tumor_vs_normal.csv", mime="text/csv")


# ---------------------------------------------------------------- Pathways
with tabs[3]:
    st.header("Pathways: what are these genes doing together?")
    st.markdown("A single gene rarely acts alone. Here we ask: **are genes from a known biological "
                "pathway over-represented among the changed genes?** The results update when you "
                "move the cutoffs in the sidebar.")
    gene_sets = load_gene_sets()
    if not gene_sets:
        st.warning("No gene sets found. Run `python prepare_data.py` to download them.")
    else:
        lib = st.selectbox("Pathway collection", list(gene_sets),
                           help="Hallmark is the easiest to start with: 50 well-defined processes.")
        universe = de.index
        cols = st.columns(2)
        for col, direction in zip(cols, ["Up in tumor", "Down in tumor"]):
            with col:
                hits = de.index[de["status"] == direction]
                res = an.enrichment(hits, universe, gene_sets[lib])
                st.subheader(f"{direction} ({len(hits)} genes)")
                if res.empty or (res["padj"] >= 0.05).all():
                    st.info("No pathways reach significance. Try loosening the cutoffs.")
                    continue
                show = res[res["padj"] < 0.05].head(15).iloc[::-1]
                fig = px.bar(show, x=-np.log10(show["padj"]), y="pathway", orientation="h",
                             color="fold_enrichment", color_continuous_scale="Reds"
                             if direction.startswith("Up") else "Blues",
                             hover_data={"overlap": True, "pathway_size": True},
                             labels={"x": "-log10(adjusted p-value)", "pathway": "",
                                     "fold_enrichment": "Fold<br>enrichment"})
                fig.update_layout(height=max(300, 32 * len(show)), margin=dict(t=10, l=10))
                st.plotly_chart(fig, width="stretch")
                with st.expander("Table and genes"):
                    st.dataframe(res[["pathway", "overlap", "pathway_size", "fold_enrichment", "padj", "genes"]],
                                 width="stretch", height=250)
        with st.expander("How does this work?"):
            st.markdown(
                "For each pathway we count how many of its genes are in our 'changed' list, then "
                "use a **hypergeometric test** to ask how likely that overlap would be if we had "
                "picked genes at random. Because we test hundreds of pathways, p-values are "
                "adjusted (Benjamini-Hochberg) to avoid false alarms. This is the same idea used "
                "by tools like Enrichr and DAVID.")


# ---------------------------------------------------------------- Compare
with tabs[4]:
    st.header("Compare two cancer types")
    others = [c for c in cancers if c != code]
    if not others:
        st.info("Add at least two cancer types in prepare_data.py to use this tab.")
    else:
        other = st.selectbox("Compare with", others, format_func=lambda c: f"{cancers[c]} ({c})")
        d2 = load_de(other)
        d2["status"] = an.classify(d2, padj_cut, lfc_cut)
        both = de[["log2FC", "status"]].join(d2[["log2FC", "status"]], lsuffix=f"_{code}",
                                             rsuffix=f"_{other}", how="inner")
        s1, s2 = both[f"status_{code}"], both[f"status_{other}"]
        both["group"] = np.select(
            [(s1 == "Up in tumor") & (s2 == "Up in tumor"),
             (s1 == "Down in tumor") & (s2 == "Down in tumor"),
             (s1 != "Not significant") & (s2 != "Not significant"),
             s1 != "Not significant", s2 != "Not significant"],
            ["Up in both", "Down in both", "Opposite directions", f"Only {code}", f"Only {other}"],
            "Neither")
        counts = both["group"].value_counts()
        c = st.columns(4)
        c[0].metric("Up in both", int(counts.get("Up in both", 0)))
        c[1].metric("Down in both", int(counts.get("Down in both", 0)))
        c[2].metric(f"Only {code}", int(counts.get(f"Only {code}", 0)))
        c[3].metric(f"Only {other}", int(counts.get(f"Only {other}", 0)))

        gcolors = {"Up in both": "#d1495b", "Down in both": "#2e86ab", "Opposite directions": "#edae49",
                   f"Only {code}": "#8e6c8a", f"Only {other}": "#66a182", "Neither": "#dddddd"}
        fig = px.scatter(both.reset_index(), x=f"log2FC_{code}", y=f"log2FC_{other}", color="group",
                         hover_name="gene", color_discrete_map=gcolors, opacity=0.7,
                         labels={f"log2FC_{code}": f"log2FC in {code}", f"log2FC_{other}": f"log2FC in {other}"})
        fig.update_traces(marker=dict(size=5))
        fig.add_hline(y=0, line_color="gray")
        fig.add_vline(x=0, line_color="gray")
        fig.update_layout(height=520, legend_title=None, margin=dict(t=10))
        st.plotly_chart(fig, width="stretch")
        r = np.corrcoef(both[f"log2FC_{code}"], both[f"log2FC_{other}"])[0, 1]
        st.caption(f"Correlation of fold changes: r = {r:.2f}. Genes in the top-right corner go up in "
                   "both cancers and may be part of a shared 'cancer program', like cell division.")
        l, rt = st.columns(2)
        for col, grp in zip([l, rt], ["Up in both", "Down in both"]):
            genes_grp = both[both["group"] == grp]
            genes_grp = genes_grp.reindex(
                (genes_grp[f"log2FC_{code}"] + genes_grp[f"log2FC_{other}"]).abs()
                .sort_values(ascending=False).index)
            col.markdown(f"**{grp}** (strongest first)")
            col.dataframe(genes_grp[[f"log2FC_{code}", f"log2FC_{other}"]].head(50),
                          width="stretch", height=250)


# ---------------------------------------------------------------- Learn
with tabs[5]:
    st.header("Learn: how to read this app")
    st.markdown("""
### The big idea
Every cell in your body has the same DNA, but each cell **switches genes on and off**
differently. A lung cell switches on lung genes; a dividing cell switches on cell-division genes.
**RNA sequencing** measures how "switched on" each of ~20,000 genes is by counting its RNA
messages. Cancer cells change which genes are on, and comparing tumors to normal tissue shows
us *what changed*.

### Key terms
- **Differential expression**: a gene is used more (or less) in tumors than in normal tissue.
- **log2 fold change (log2FC)**: how big the change is. `1` = 2x higher, `2` = 4x higher,
  `-1` = 2x lower. `0` = no change.
- **Adjusted p-value**: how confident we are the change is real and not random noise. We test
  ~20,000 genes, so p-values are corrected to limit false positives.
- **Pathway**: a group of genes that work together on one job, like copying DNA or responding
  to low oxygen.
- **Kaplan-Meier curve**: shows what fraction of patients are still alive over time.

### Important limitations
- TCGA "normal" samples come from tissue **next to** the tumor, not from healthy people.
- A gene that changes is not necessarily a **cause** of the cancer; it could be an effect.
- Tumors are a mix of cancer cells, immune cells and support cells, and all of them contribute RNA.
- Some key cancer genes (like **TP53**) are damaged by **mutation**, which RNA levels can miss.
- Survival plots show associations and don't account for age, stage or treatment.
- This app is for learning, **not** for medical decisions.

### Methods
Data: TCGA RNA-seq (STAR gene counts) and overall survival from the UCSC Xena GDC hub. Primary
tumors (up to 200, randomly chosen) were compared with solid-tissue normals using
**PyDESeq2** (a Python version of DESeq2). Pathway analysis is an over-representation test on
gene sets from Enrichr (MSigDB Hallmark, KEGG, Reactome, GO).

### About
Built by Charlene Siawira as an independent project. Code on
[GitHub](https://github.com/your-username/cancer-gene-explorer). Inspired by professional tools
such as GEPIA2 and UALCAN, but designed for learners.
""")
