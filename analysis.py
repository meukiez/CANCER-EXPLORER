"""The math behind the app, kept separate so app.py stays readable."""

import numpy as np
import pandas as pd
from scipy import stats


def classify(de, padj_cut, lfc_cut):
    """Label each gene Up / Down / Not significant."""
    sig = (de["padj"] < padj_cut) & (de["log2FC"].abs() >= lfc_cut)
    return np.select([sig & (de["log2FC"] > 0), sig & (de["log2FC"] < 0)],
                     ["Up in tumor", "Down in tumor"], "Not significant")


def pca(expr, labels, n_genes=1000):
    """PCA on the most variable genes. Returns a DataFrame with PC1, PC2 and % variance."""
    top = expr.var().sort_values(ascending=False).index[:n_genes]
    x = expr[top].values.astype(float)
    x = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    var = s ** 2 / np.sum(s ** 2) * 100
    out = pd.DataFrame({"PC1": u[:, 0] * s[0], "PC2": u[:, 1] * s[1],
                        "type": labels.values}, index=expr.index)
    return out, var[:2]


def kaplan_meier(time, event):
    """Kaplan-Meier survival curve. Returns (times, survival probabilities) as step points."""
    order = np.argsort(time)
    time, event = np.asarray(time)[order], np.asarray(event)[order]
    t_out, s_out, surv = [0.0], [1.0], 1.0
    for t in np.unique(time[event == 1]):
        at_risk = np.sum(time >= t)
        died = np.sum((time == t) & (event == 1))
        surv *= 1 - died / at_risk
        t_out.append(t)
        s_out.append(surv)
    t_out.append(time.max())
    s_out.append(surv)
    return np.array(t_out), np.array(s_out)


def logrank_p(t1, e1, t2, e2):
    """Log-rank test: are two survival curves different? Returns a p-value."""
    t = np.concatenate([t1, t2])
    e = np.concatenate([e1, e2])
    g = np.concatenate([np.zeros(len(t1)), np.ones(len(t2))])
    obs_minus_exp, var = 0.0, 0.0
    for time in np.unique(t[e == 1]):
        at_risk = t >= time
        n, n1 = at_risk.sum(), (at_risk & (g == 0)).sum()
        d = ((t == time) & (e == 1)).sum()
        d1 = ((t == time) & (e == 1) & (g == 0)).sum()
        obs_minus_exp += d1 - d * n1 / n
        if n > 1:
            var += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if var == 0:
        return np.nan
    return float(stats.chi2.sf(obs_minus_exp ** 2 / var, df=1))


def bh_adjust(p):
    """Benjamini-Hochberg correction for testing many things at once."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def enrichment(hits, universe, gene_sets, min_size=10, max_size=500):
    """Over-representation analysis (a hypergeometric test, like Enrichr/DAVID).

    Question for each pathway: do my significant genes include MORE pathway
    members than we'd expect by chance?
    """
    universe = set(universe)
    hits = set(hits) & universe
    N, n = len(universe), len(hits)
    rows = []
    if n == 0:
        return pd.DataFrame()
    for name, genes in gene_sets.items():
        members = set(genes) & universe
        K = len(members)
        if not (min_size <= K <= max_size):
            continue
        overlap = hits & members
        k = len(overlap)
        if k == 0:
            continue
        p = stats.hypergeom.sf(k - 1, N, K, n)
        rows.append({"pathway": name, "overlap": k, "pathway_size": K,
                     "gene_ratio": k / K, "fold_enrichment": (k / n) / (K / N),
                     "pvalue": p, "genes": ", ".join(sorted(overlap))})
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["padj"] = bh_adjust(res["pvalue"])
    return res.sort_values("pvalue").reset_index(drop=True)
