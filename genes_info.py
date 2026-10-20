<<<<<<< HEAD
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

def get_gene_description(gene_symbol):
    """Returns plain-language summary for known genes or a fallback default message."""
    return GENE_DESCRIPTIONS.get(
        gene_symbol.upper(),
        f"**{gene_symbol}** is an active human gene monitored in TCGA differential expression datasets."
    )
=======
"""Plain-language notes for genes with a well-known cancer story.
Add your own! Keep each note to 1-2 sentences a classmate could understand."""

GENE_NOTES = {
    "MKI67": "Makes the Ki-67 protein, which is only present in cells that are actively dividing. "
             "Doctors use it as a 'speedometer' for how fast a tumor is growing.",
    "TOP2A": "An enzyme that untangles DNA so a cell can divide. Several chemotherapy drugs "
             "(like doxorubicin) work by blocking it.",
    "CCNB1": "Cyclin B1, a 'go' signal that pushes a cell into division.",
    "CDK1": "A key engine of the cell cycle; it teams up with cyclin B1 to trigger cell division.",
    "ERBB2": "Codes for HER2, a growth-signal receptor. In some breast cancers extra copies of this "
             "gene make cells grow too fast; the drug trastuzumab (Herceptin) targets it.",
    "ESR1": "The estrogen receptor. Most breast cancers are 'ER-positive' and grow in response to "
            "estrogen, which is why hormone therapies like tamoxifen work.",
    "BRCA1": "Helps repair broken DNA. Inherited mutations in BRCA1 raise the risk of breast and "
             "ovarian cancer.",
    "CA9": "Carbonic anhydrase IX, normally switched on when cells lack oxygen. In clear cell kidney "
           "cancer the VHL gene is often broken, so cells act as if they are always low on oxygen.",
    "NDUFA4L2": "Another 'low-oxygen' gene that is strongly switched on in clear cell kidney cancer.",
    "VEGFA": "Sends signals to grow new blood vessels. Tumors use it to get a blood supply; "
             "several drugs block this pathway.",
    "SFTPC": "Surfactant protein C, made by healthy lung air-sac cells to keep the lungs open. "
             "Tumor cells often lose this 'lung identity'.",
    "EGFR": "A growth-signal receptor. Some lung adenocarcinomas have EGFR mutations and can be "
            "treated with targeted pills such as osimertinib.",
    "AMACR": "An enzyme that is much higher in prostate cancer cells; pathologists use it to help "
             "spot cancer in biopsies.",
    "TP53": "The 'guardian of the genome' and the most commonly mutated gene in cancer. It is usually "
            "damaged by mutation rather than changed in amount, so expression data can miss it.",
    "MYC": "A master switch for cell growth. It is overactive in many cancers.",
    "COL11A1": "A collagen gene made by 'cancer-associated fibroblasts', the support cells that "
               "tumors recruit to remodel the tissue around them.",
    "MMP11": "An enzyme that cuts up the scaffolding between cells, helping tumors remodel and invade tissue.",
    "GAPDH": "A 'housekeeping' gene used for basic energy metabolism in almost every cell. A good "
             "control: it should look similar in tumor and normal.",
    "ACTB": "Beta-actin, part of every cell's skeleton. Another housekeeping control gene.",
}

# Suggested genes to try for each cancer (only the ones present in the data are shown)
EXAMPLES = {
    "BRCA": ["ERBB2", "ESR1", "MKI67", "COL11A1", "GAPDH"],
    "KIRC": ["CA9", "NDUFA4L2", "VEGFA", "MKI67", "GAPDH"],
    "LUAD": ["SFTPC", "EGFR", "TOP2A", "MKI67", "GAPDH"],
}
>>>>>>> e6747d96e5d679f63505dd0d562c7b7dc818bedb
