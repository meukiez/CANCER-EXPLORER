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
