"""
genes_info.py
-------------
Stores biological summaries and curated key driver genes for each cancer type.
"""

CAN_INFO = {
    "BRCA": {
        "name": "Breast Invasive Carcinoma",
        "description": "The most common cancer among women worldwide. Driven by estrogen receptor (ER), progesterone receptor (PR), and HER2 expression subtypes.",
        "key_genes": ["TP53", "ERBB2", "ESR1", "BRCA1", "BRCA2", "PIK3CA"]
    },
    "KIRC": {
        "name": "Kidney Renal Clear Cell Carcinoma",
        "description": "The most common type of kidney cancer. Characterized by loss-of-function mutations in VHL.",
        "key_genes": ["VHL", "PBRM1", "SETD2", "BAP1", "VEGFA", "EGLN1"]
    },
    "LUAD": {
        "name": "Lung Adenocarcinoma",
        "description": "A prevalent form of non-small cell lung cancer (NSCLC). Frequently exhibits EGFR and KRAS mutations.",
        "key_genes": ["EGFR", "KRAS", "ALK", "TP53", "NKX2-1", "ROS1"]
    },
    "COAD": {
        "name": "Colon Adenocarcinoma",
        "description": "A prominent gastrointestinal cancer driven by alterations in Wnt/β-catenin, EGFR, and MMR pathways.",
        "key_genes": ["APC", "KRAS", "TP53", "SMAD4", "BRAF", "PIK3CA"]
    },
    "PRAD": {
        "name": "Prostate Adenocarcinoma",
        "description": "Strongly influenced by androgen receptor (AR) activation and PTEN deletions.",
        "key_genes": ["AR", "PTEN", "ERG", "TP53", "SPOP", "FOXA1"]
    },
    "LIHC": {
        "name": "Liver Hepatocellular Carcinoma",
        "description": "Primary liver cancer associated with chronic liver injury and metabolic dysfunction.",
        "key_genes": ["TP53", "CTNNB1", "ALB", "AFP", "AXIN1", "ARID1A"]
    },
    "SKCM": {
        "name": "Skin Cutaneous Melanoma",
        "description": "Aggressive skin cancer driven primarily by UV-induced somatic mutations in RAF/RAS pathways.",
        "key_genes": ["BRAF", "NRAS", "CDKN2A", "TP53", "PTEN", "KIT"]
    },
    "UCEC": {
        "name": "Uterine Corpus Endometrial Carcinoma",
        "description": "Features high rates of PTEN mutation, MSI, and PI3K pathway hyperactivation.",
        "key_genes": ["PTEN", "PIK3CA", "ARID1A", "TP53", "CTNNB1", "KRAS"]
    },
    "THCA": {
        "name": "Thyroid Carcinoma",
        "description": "Characterized by mutations in BRAF V600E, RET/PTC rearrangements, or RAS alterations.",
        "key_genes": ["BRAF", "NRAS", "HRAS", "KRAS", "TG", "RET"]
    }
}

GENE_DESCRIPTIONS = {
    "TP53": "Known as the 'Guardian of the Genome', TP53 encodes a tumor suppressor protein regulating apoptosis and DNA repair.",
    "ERBB2": "Also called HER2; encodes a receptor tyrosine kinase overexpressed in aggressive breast cancers.",
    "ESR1": "Encodes Estrogen Receptor Alpha. Crucial for hormone-dependent signaling in breast cancer.",
    "BRCA1": "A critical tumor suppressor involved in repairing double-strand DNA breaks.",
    "BRCA2": "Works alongside BRCA1 to maintain genomic stability through DNA repair.",
    "VHL": "Von Hippel-Lindau tumor suppressor. Inactivation drives clear cell renal carcinoma.",
    "PBRM1": "A chromatin remodeling gene frequently mutated in renal cell carcinoma.",
    "EGFR": "Epidermal Growth Factor Receptor. Driver mutations spur unchecked cell proliferation.",
    "KRAS": "A major GTPase oncogene. Mutations lock KRAS in an active signaling state.",
    "APC": "Adenomatous Polyposis Coli. A key negative regulator of the Wnt/β-catenin pathway.",
    "AR": "Androgen Receptor. The master transcriptomic driver in prostate adenocarcinoma.",
    "PTEN": "A lipid phosphatase tumor suppressor regulating the PI3K/AKT/mTOR axis.",
    "BRAF": "A serine/threonine protein kinase. The BRAF V600E mutation hyperactivates the MAPK pathway.",
    "CTNNB1": "Encodes β-Catenin, involved in cell adhesion and transcriptomic activation in Wnt signaling.",
    "CDKN2A": "Encodes p16(INK4a) and p14(ARF) tumor suppressors regulating CDK4/6 progression."
}

EXAMPLES = {
    "BRCA": ["TP53", "ERBB2", "ESR1", "BRCA1"],
    "KIRC": ["VHL", "PBRM1", "SETD2"],
    "LUAD": ["EGFR", "KRAS", "ALK"],
    "COAD": ["APC", "KRAS", "TP53"],
    "PRAD": ["AR", "PTEN", "ERG"],
    "LIHC": ["TP53", "CTNNB1", "ALB"],
    "SKCM": ["BRAF", "NRAS", "CDKN2A"],
    "UCEC": ["PTEN", "PIK3CA", "ARID1A"],
    "THCA": ["BRAF", "NRAS", "RET"]
}

GENE_NOTES = GENE_DESCRIPTIONS

def get_gene_description(gene_symbol):
    """Returns plain-language summary for known genes or a fallback default message."""
    return GENE_DESCRIPTIONS.get(
        gene_symbol.upper(),
        f"**{gene_symbol}** is an active human gene monitored in TCGA differential expression datasets."
    )