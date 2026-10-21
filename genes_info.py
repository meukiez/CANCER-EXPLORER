"""
genes_info.py
-------------
Stores biological summaries and curated key driver genes for each cancer type
to display beginner-friendly gene explanations in the Streamlit app.
"""

CAN_INFO = {
    "BRCA": {
        "name": "Breast Invasive Carcinoma",
        "description": "The most common cancer among women worldwide. Driven by estrogen receptor (ER), progesterone receptor (PR), and HER2 expression subtypes.",
        "key_genes": ["TP53", "ERBB2", "ESR1", "BRCA1", "BRCA2", "PIK3CA"]
    },
    "KIRC": {
        "name": "Kidney Renal Clear Cell Carcinoma",
        "description": "The most common type of kidney cancer. Characterized by loss-of-function mutations in the VHL tumor suppressor and alterations in hypoxia-inducible factors (HIF).",
        "key_genes": ["VHL", "PBRM1", "SETD2", "BAP1", "VEGFA", "EGLN1"]
    },
    "LUAD": {
        "name": "Lung Adenocarcinoma",
        "description": "A prevalent form of non-small cell lung cancer (NSCLC). Frequently exhibits oncogenic driver mutations in EGFR and KRAS pathways.",
        "key_genes": ["EGFR", "KRAS", "ALK", "TP53", "NKX2-1", "ROS1"]
    },
    "COAD": {
        "name": "Colon Adenocarcinoma",
        "description": "A prominent gastrointestinal cancer. Driven by key alterations in Wnt/β-catenin, EGFR, and Mismatch Repair (MMR) pathway signaling.",
        "key_genes": ["APC", "KRAS", "TP53", "SMAD4", "BRAF", "PIK3CA"]
    },
    "PRAD": {
        "name": "Prostate Adenocarcinoma",
        "description": "The second most common cancer in men. Strongly influenced by androgen receptor (AR) activation, TMPRSS2-ERG fusions, and PTEN deletions.",
        "key_genes": ["AR", "PTEN", "ERG", "TP53", "SPOP", "FOXA1"]
    },
    "LIHC": {
        "name": "Liver Hepatocellular Carcinoma",
        "description": "The major form of primary liver cancer, strongly associated with chronic liver injury, viral hepatitis (HBV/HCV), and metabolic dysfunction.",
        "key_genes": ["TP53", "CTNNB1", "ALB", "AFP", "AXIN1", "ARID1A"]
    },
    "SKCM": {
        "name": "Skin Cutaneous Melanoma",
        "description": "An aggressive skin cancer originating from melanocytes. Driven primarily by UV-induced somatic mutations in RAF/RAS signaling pathways.",
        "key_genes": ["BRAF", "NRAS", "CDKN2A", "TP53", "PTEN", "KIT"]
    },
    "UCEC": {
        "name": "Uterine Corpus Endometrial Carcinoma",
        "description": "A primary cancer of the uterine lining, often featuring high rates of PTEN mutation, microsatellite instability (MSI), and PI3K pathway hyperactivation.",
        "key_genes": ["PTEN", "PIK3CA", "ARID1A", "TP53", "CTNNB1", "KRAS"]
    },
    "THCA": {
        "name": "Thyroid Carcinoma",
        "description": "The most frequent endocrine malignancy. Commonly characterized by targeted mutations in BRAF V600E, RET/PTC rearrangements, or RAS alterations.",
        "key_genes": ["BRAF", "NRAS", "HRAS", "KRAS", "TG", "RET"]
    }
}

GENE_DESCRIPTIONS = {
    "TP53": "Known as the 'Guardian of the Genome', TP53 encodes a tumor suppressor protein that regulates cell cycle arrest, DNA repair, and apoptosis.",
    "ERBB2": "Also called HER2; encodes a receptor tyrosine kinase. Overexpressed or amplified in aggressive subtypes of breast and gastric cancers.",
    "ESR1": "Encodes Estrogen Receptor Alpha. Crucial for hormone-dependent breast cancer signaling and targeted hormone therapy response.",
    "BRCA1": "A critical tumor suppressor involved in homologous recombination repair of double-strand DNA breaks.",
    "BRCA2": "Works alongside BRCA1 to maintain genomic stability through error-free double-strand DNA break repair.",
    "VHL": "Von Hippel-Lindau tumor suppressor. Inactivation drives clear cell renal carcinoma via HIF accumulation and VEGF-mediated blood vessel proliferation.",
    "PBRM1": "A chromatin remodeling gene frequently mutated in clear cell renal cell carcinoma, playing a key role in epigenetic regulation.",
    "EGFR": "Epidermal Growth Factor Receptor. Driver mutations spur unchecked cell proliferation, commonly targeted by tyrosine kinase inhibitors in lung adenocarcinoma.",
    "KRAS": "A major GTPase oncogene. Mutations lock KRAS in an active state, signaling continuous cellular proliferation across colorectal, lung, and pancreatic cancers.",
    "APC": "Adenomatous Polyposis Coli. A key negative regulator of the Wnt/β-catenin pathway; its loss initiates most colorectal adenomas and carcinomas.",
    "AR": "Androgen Receptor. The master transcriptomic driver in prostate epithelial development and prostate adenocarcinoma progression.",
    "PTEN": "A lipid phosphatase tumor suppressor that negatively regulates the pro-survival PI3K/AKT/mTOR signaling axis.",
    "BRAF": "A serine/threonine protein kinase. The BRAF V600E mutation hyperactivates the MAPK pathway in melanomas, thyroid, and colorectal cancers.",
    "CTNNB1": "Encodes β-Catenin, a dual-function protein involved in cell-cell adhesion and transcriptomic activation within the Wnt pathway.",
    "CDKN2A": "Encodes the p16(INK4a) and p14(ARF) tumor suppressors, regulating CDK4/6-dependent cell cycle progression and p53 stability."
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

def get_gene_description(gene_symbol):
    """Returns plain-language summary for known genes or a fallback default message."""
    return GENE_DESCRIPTIONS.get(
        gene_symbol.upper(),
        f"**{gene_symbol}** is an active human gene monitored in TCGA differential expression datasets."
    )

EXAMPLES = {} 
GENE_NOTES = GENE_DESCRIPTIONS