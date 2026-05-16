"""
Centralized configuration for the unified PGx Reporter pipeline.

Holds paths, JAR/Java auto-detection, default timeouts, the gene
functional categorization, drug therapeutic categories, action
symbols (colorblind-safe), phenotype-to-risk mapping, and the
disclaimer text used by both reports.
"""

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# src/ directory (this file lives in src/config/)
SRC_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SRC_DIR.parent

# PharmCAT JAR — checked in this order:
#   1) explicit CLI flag
#   2) lib/pharmcat*.jar (any version, sorted descending)
#   3) lib/pharmcat.jar (Marco's hardcoded path, kept as fallback)
LIB_DIR = PROJECT_ROOT / "lib"
PHARMCAT_JAR_DEFAULT = LIB_DIR / "pharmcat.jar"

# Templates directory (Marco's location: reports/templates/)
TEMPLATES_DIR = SRC_DIR / "reports" / "templates"

# Default output directory
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output"


# ---------------------------------------------------------------------------
# Java auto-detection
# ---------------------------------------------------------------------------

def find_java_executable() -> str:
    """Locate a Java binary. Prefers a bundled JDK in lib/jdk/, then Eclipse
    Adoptium on Windows, then `java` on PATH."""
    # Bundled JDK under lib/jdk/ (Adib's pattern)
    for candidate in (
        LIB_DIR / "jdk" / "Contents" / "Home" / "bin" / "java",
        LIB_DIR / "jdk" / "bin" / "java",
        LIB_DIR / "jdk" / "bin" / "java.exe",
    ):
        if candidate.exists():
            return str(candidate)

    # Eclipse Adoptium on Windows (Marco's pattern)
    adoptium = Path(os.environ.get("ProgramFiles", r"C:\Program Files")) / "Eclipse Adoptium"
    if adoptium.is_dir():
        for entry in sorted(os.listdir(adoptium), reverse=True):
            cand = adoptium / entry / "bin" / "java.exe"
            if cand.is_file():
                return str(cand)

    return "java"


JAVA_EXECUTABLE = find_java_executable()


# ---------------------------------------------------------------------------
# PharmCAT execution defaults
# ---------------------------------------------------------------------------

# Default timeout (seconds). Whole-genome VCFs can take 30–60 min without
# preprocessing; with --preprocess the run is normally a few seconds.
PHARMCAT_TIMEOUT_DEFAULT = 3600

# Default research-mode flags. Empty by default: PharmCAT 3.2.0+ disables the
# full reporter (report.json / report.html) when ANY research mode is on, so
# turning research on costs you all CPIC/DPWG drug recommendations. Pass
# --research cyp2d6 explicitly to opt into best-effort CYP2D6 calling at that cost.
PHARMCAT_RESEARCH_DEFAULT: list[str] = []

# Max VCF size accepted by the validator (MB). Sized for HG005 WGS.
MAX_VCF_SIZE_MB = 3000


# ---------------------------------------------------------------------------
# Gene functional categories (Adib)
# ---------------------------------------------------------------------------

GENE_CATEGORIES: dict[str, dict] = {
    "Phase I Metabolism (CYP Enzymes)": {
        "genes": {"CYP2B6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2"},
        "description": (
            "Cytochrome P450 enzymes are the body's primary drug-metabolizing system. "
            "They break down many common medications in the liver."
        ),
    },
    "Phase II Metabolism": {
        "genes": {"DPYD", "NAT2", "TPMT", "NUDT15", "UGT1A1"},
        "description": (
            "These enzymes modify drugs and their breakdown products to make them easier "
            "to eliminate. Variations can affect how the body handles certain chemotherapy "
            "agents and other medications."
        ),
    },
    "Drug Transporters": {
        "genes": {"ABCG2", "SLCO1B1"},
        "description": (
            "Transporter proteins move drugs into and out of cells. Variations can alter "
            "drug absorption and distribution, affecting how much medication reaches its target."
        ),
    },
    "Immune Markers (HLA)": {
        "genes": {"HLA-A", "HLA-B"},
        "description": (
            "HLA genes help the immune system distinguish the body's own cells from "
            "foreign substances. Certain HLA variants are associated with severe drug "
            "hypersensitivity reactions."
        ),
    },
    "Other Pharmacogenes": {
        "genes": {"CACNA1S", "CFTR", "F2", "F5", "G6PD", "IFNL3", "MT-RNR1", "RYR1", "VKORC1"},
        "description": (
            "Additional genes that influence drug response through various mechanisms "
            "including drug targets, enzyme deficiencies, and receptor sensitivity."
        ),
    },
}

CATEGORY_ORDER: list[str] = list(GENE_CATEGORIES.keys())


def gene_category(gene_symbol: str) -> tuple[str, str]:
    """Return (category_name, category_description) for a gene symbol."""
    for cat_name, info in GENE_CATEGORIES.items():
        if gene_symbol in info["genes"]:
            return cat_name, info["description"]
    return "Other Pharmacogenes", GENE_CATEGORIES["Other Pharmacogenes"]["description"]


# ---------------------------------------------------------------------------
# Per-gene plain-language metadata (for the patient report)
# ---------------------------------------------------------------------------

GENE_DESCRIPTIONS: dict[str, str] = {
    "CYP2C19": "helps the body break down many medications including blood thinners, antidepressants, and stomach acid drugs",
    "CYP2C9": "helps the body process pain medications, blood thinners like warfarin, and anti-inflammatory drugs",
    "CYP2D6": "is involved in processing about 25% of all medications, including many antidepressants and pain medications",
    "CYP3A4": "is the most important drug-metabolizing enzyme, involved in processing nearly half of all medications",
    "CYP3A5": "helps process immunosuppressants like tacrolimus and some blood pressure medications",
    "CYP2B6": "helps process certain HIV medications and the smoking cessation drug bupropion",
    "CYP4F2": "affects how the body uses vitamin K, which influences warfarin dosing",
    "DPYD": "determines how the body handles fluoropyrimidine chemotherapy drugs",
    "TPMT": "affects how the body processes thiopurine medications used for cancer and autoimmune conditions",
    "NUDT15": "also affects thiopurine medication processing, similar to TPMT",
    "VKORC1": "controls sensitivity to warfarin, a common blood thinner",
    "SLCO1B1": "affects how the body transports statin medications used to lower cholesterol",
    "UGT1A1": "helps the body process certain medications and bilirubin (a waste product)",
    "ABCG2": "is a transporter protein that affects how the body absorbs certain medications",
    "CACNA1S": "is related to response to certain anesthetics used during surgery",
    "CFTR": "is related to cystic fibrosis and affects certain medication responses",
    "G6PD": "affects red blood cells' sensitivity to certain medications",
    "HLA-A": "is part of the immune system and affects risk of severe drug reactions",
    "HLA-B": "is part of the immune system and affects risk of severe drug reactions to specific drugs",
    "IFNL3": "affects response to hepatitis C treatments",
    "MT-RNR1": "affects sensitivity to aminoglycoside antibiotics that can cause hearing loss",
    "NAT2": "helps the body process certain medications including some antibiotics",
    "RYR1": "is related to risk of malignant hyperthermia, a rare reaction to anesthetics",
    "F2": "affects blood clotting and may influence response to anticoagulant medications",
    "F5": "affects blood clotting and may influence response to anticoagulant medications",
}

GENE_PROTEIN_TYPE: dict[str, str] = {
    "CYP2C19": "Liver enzyme", "CYP2C9": "Liver enzyme", "CYP2D6": "Liver enzyme",
    "CYP3A4": "Liver enzyme", "CYP3A5": "Liver enzyme", "CYP2B6": "Liver enzyme",
    "CYP4F2": "Liver enzyme",
    "DPYD": "Metabolic enzyme", "TPMT": "Metabolic enzyme", "NUDT15": "Metabolic enzyme",
    "NAT2": "Metabolic enzyme", "UGT1A1": "Metabolic enzyme",
    "G6PD": "Red-blood-cell enzyme",
    "VKORC1": "Drug target",
    "SLCO1B1": "Drug transporter", "ABCG2": "Drug transporter",
    "CACNA1S": "Ion channel", "RYR1": "Ion channel", "CFTR": "Ion channel",
    "HLA-A": "Immune marker", "HLA-B": "Immune marker",
    "IFNL3": "Immune signaling",
    "MT-RNR1": "Mitochondrial RNA",
    "F2": "Clotting factor", "F5": "Clotting factor",
}


# ---------------------------------------------------------------------------
# Per-gene primary expression tissue (for patient-report icons)
# ---------------------------------------------------------------------------

# Primary site of expression for each pharmacogene. Some genes are widely
# expressed; the value here is the tissue most relevant to PGx (i.e. where
# the drug-metabolism / drug-target activity happens). ABCG2 is dual-site
# (liver + intestine) — we tag it as liver to match the other transporters.
GENE_TISSUE: dict[str, str] = {
    "CYP2B6": "liver", "CYP2C9": "liver", "CYP2C19": "liver",
    "CYP2D6": "liver", "CYP3A4": "liver", "CYP3A5": "liver",
    "CYP4F2": "liver",
    "DPYD": "liver", "NAT2": "liver", "TPMT": "liver", "UGT1A1": "liver",
    "F2": "liver", "F5": "liver", "IFNL3": "liver", "VKORC1": "liver",
    "SLCO1B1": "liver", "ABCG2": "liver",
    "NUDT15": "immune",
    "HLA-A": "immune", "HLA-B": "immune",
    "G6PD": "rbc",
    "CACNA1S": "muscle", "RYR1": "muscle",
    "CFTR": "lung",
    "MT-RNR1": "mitochondria",
}

TISSUE_LABELS: dict[str, str] = {
    "liver": "Liver",
    "immune": "Immune cells",
    "rbc": "Red blood cells",
    "muscle": "Skeletal muscle",
    "lung": "Lung / airway",
    "mitochondria": "Mitochondria",
}

# Tissue icons. We use Unicode emoji where Unicode has one (RBC, muscle,
# lung, immune-as-shield). For liver and mitochondria, no Unicode emoji
# exists, so we fall back to small inline SVG silhouettes — smooth vector
# paths (not pixel art) so they read at small sizes.

_LIVER_SVG = (
    '<svg viewBox="0 0 24 24" width="22" height="22" aria-hidden="true" '
    'style="vertical-align:middle">'
    # Triangular liver silhouette: wide rounded right base (right lobe)
    # tapering to a point on the left (left lobe tip), with falciform line
    # dividing the lobes.
    '<path d="M3 12 L20 6 Q22 6 22 9 L22 15 Q22 18 20 18 Z" '
    'fill="#A0522D" stroke="#5C2C0A" stroke-width="1.2" '
    'stroke-linejoin="round"/>'
    # Falciform ligament.
    '<line x1="13" y1="9" x2="13" y2="15" '
    'stroke="#5C2C0A" stroke-width="1" stroke-linecap="round"/>'
    "</svg>"
)

_MITO_SVG = (
    '<svg viewBox="0 0 24 24" width="22" height="22" aria-hidden="true" '
    'style="vertical-align:middle">'
    # Outer membrane: smooth brown capsule.
    '<rect x="2" y="8" width="20" height="8" rx="4" '
    'fill="#A0522D" stroke="#5C2C0A" stroke-width="1.2"/>'
    # Inner membrane: continuous wavy serpentine line representing the
    # folded cristae running through the matrix.
    '<path d="M4 12 Q5 9.5 6 12 Q7 14.5 8 12 Q9 9.5 10 12 Q11 14.5 12 12 '
    'Q13 9.5 14 12 Q15 14.5 16 12 Q17 9.5 18 12 Q19 14.5 20 12" '
    'stroke="#3E2723" stroke-width="1.2" fill="none" stroke-linecap="round"/>'
    "</svg>"
)

TISSUE_ICONS: dict[str, str] = {
    "liver":        _LIVER_SVG,
    "rbc":          "🩸",
    "muscle":       "💪",
    "lung":         "🫁",
    "immune":       "🛡️",
    "mitochondria": _MITO_SVG,
}


def gene_tissue_icon(gene_symbol: str) -> tuple[str, str, str]:
    """Return (tissue_key, tissue_label, emoji) for a gene, or ('', '', '')
    if no tissue is mapped."""
    tissue = GENE_TISSUE.get(gene_symbol, "")
    if not tissue:
        return "", "", ""
    return tissue, TISSUE_LABELS.get(tissue, tissue), TISSUE_ICONS.get(tissue, "")


# ---------------------------------------------------------------------------
# Drug therapeutic categories (Marco)
# ---------------------------------------------------------------------------

DRUG_CATEGORIES: dict[str, str] = {
    # Antidepressants
    "amitriptyline": "Antidepressants", "citalopram": "Antidepressants",
    "clomipramine": "Antidepressants", "doxepin": "Antidepressants",
    "escitalopram": "Antidepressants", "imipramine": "Antidepressants",
    "sertraline": "Antidepressants", "trimipramine": "Antidepressants",
    "nortriptyline": "Antidepressants", "paroxetine": "Antidepressants",
    "fluvoxamine": "Antidepressants", "venlafaxine": "Antidepressants",
    "fluoxetine": "Antidepressants", "desipramine": "Antidepressants",
    "protriptyline": "Antidepressants", "amoxapine": "Antidepressants",
    "vortioxetine": "Antidepressants", "flibanserin": "Antidepressants",
    # Antiplatelet / Anticoagulants
    "clopidogrel": "Antiplatelet / Anticoagulant agents",
    "warfarin": "Antiplatelet / Anticoagulant agents",
    "acenocoumarol": "Antiplatelet / Anticoagulant agents",
    "phenprocoumon": "Antiplatelet / Anticoagulant agents",
    # Statins
    "simvastatin": "Statins (cholesterol)", "atorvastatin": "Statins (cholesterol)",
    "lovastatin": "Statins (cholesterol)", "rosuvastatin": "Statins (cholesterol)",
    "fluvastatin": "Statins (cholesterol)", "pitavastatin": "Statins (cholesterol)",
    "pravastatin": "Statins (cholesterol)",
    # Proton pump inhibitors
    "omeprazole": "Proton pump inhibitors (stomach acid)",
    "lansoprazole": "Proton pump inhibitors (stomach acid)",
    "pantoprazole": "Proton pump inhibitors (stomach acid)",
    "dexlansoprazole": "Proton pump inhibitors (stomach acid)",
    "esomeprazole": "Proton pump inhibitors (stomach acid)",
    "rabeprazole": "Proton pump inhibitors (stomach acid)",
    # Antifungals
    "voriconazole": "Antifungals",
    # Cardiovascular
    "mavacamten": "Cardiovascular", "metoprolol": "Cardiovascular",
    "propranolol": "Cardiovascular", "carvedilol": "Cardiovascular",
    "nebivolol": "Cardiovascular", "flecainide": "Cardiovascular",
    "propafenone": "Cardiovascular",
    "hydralazine": "Cardiovascular",
    # Immunosuppressants
    "tacrolimus": "Immunosuppressants",
    # Chemotherapy / Immunosuppressants
    "fluorouracil": "Chemotherapy", "capecitabine": "Chemotherapy",
    "mercaptopurine": "Chemotherapy / Immunosuppressants",
    "azathioprine": "Chemotherapy / Immunosuppressants",
    "thioguanine": "Chemotherapy / Immunosuppressants",
    "tamoxifen": "Chemotherapy",
    "irinotecan": "Chemotherapy",
    # Pain medications
    "codeine": "Pain medications", "tramadol": "Pain medications",
    "hydrocodone": "Pain medications",
    # Antivirals (HIV)
    "abacavir": "Antivirals (HIV)", "atazanavir": "Antivirals (HIV)",
    # Antiepileptics
    "carbamazepine": "Antiepileptics", "phenytoin": "Antiepileptics",
    "oxcarbazepine": "Antiepileptics", "brivaracetam": "Antiepileptics",
    "clobazam": "Antiepileptics", "diazepam": "Antiepileptics",
    # Antipsychotics
    "aripiprazole": "Antipsychotics", "haloperidol": "Antipsychotics",
    "pimozide": "Antipsychotics", "clozapine": "Antipsychotics",
    "risperidone": "Antipsychotics", "brexpiprazole": "Antipsychotics",
    "iloperidone": "Antipsychotics", "perphenazine": "Antipsychotics",
    "thioridazine": "Antipsychotics", "zuclopenthixol": "Antipsychotics",
    # Gout
    "allopurinol": "Gout medications",
    # ADHD
    "atomoxetine": "ADHD medications", "amphetamine": "ADHD medications",
    "viloxazine": "ADHD medications",
    # Anesthetics
    "desflurane": "Anesthetics", "enflurane": "Anesthetics",
    "halothane": "Anesthetics", "isoflurane": "Anesthetics",
    "sevoflurane": "Anesthetics", "succinylcholine": "Anesthetics",
    # Anti-nausea
    "ondansetron": "Anti-nausea medications", "tropisetron": "Anti-nausea medications",
    # Other specific
    "eliglustat": "Gaucher disease medications",
    "gefitinib": "Targeted cancer therapy",
    "donepezil": "Alzheimer medications",
    "galantamine": "Alzheimer medications",
    "metoclopramide": "GI motility agents",
    "lofexidine": "Opioid withdrawal agents",
    "pitolisant": "Narcolepsy medications",
    "belzutifan": "Targeted cancer therapy",
    "abrocitinib": "Dermatology",
    # CFTR modulators (effectiveness depends on the patient's CFTR genotype)
    "ivacaftor": "Cystic Fibrosis Modulators",
    "lumacaftor/ivacaftor": "Cystic Fibrosis Modulators",
    "tezacaftor/ivacaftor": "Cystic Fibrosis Modulators",
    "elexacaftor/tezacaftor/ivacaftor": "Cystic Fibrosis Modulators",
}


def drug_category(drug_name: str) -> str:
    """Return the therapeutic category for a drug, or 'Other medications'."""
    return DRUG_CATEGORIES.get(drug_name.strip().lower(), "Other medications")


# ---------------------------------------------------------------------------
# Risk levels — colorblind-safe, 4-level (Adib)
# ---------------------------------------------------------------------------

# Order matters: action > review > normal > nodata for sorting.
RISK_LEVELS: tuple[str, ...] = ("action", "review", "normal", "nodata")
RISK_PRIORITY: dict[str, int] = {level: i for i, level in enumerate(RISK_LEVELS)}

ACTION_SYMBOLS: dict[str, str] = {
    "action": "▲",   # ▲
    "review": "◆",   # ◆
    "normal": "✓",   # ✓
    "nodata": "—",   # —
}

ACTION_LABELS: dict[str, str] = {
    "action": "Action",
    "review": "Review",
    "normal": "Normal",
    "nodata": "No Data",
}

# Phenotype string -> 4-level risk level
PHENOTYPE_RISK: dict[str, str] = {
    "poor metabolizer": "action",
    "ultrarapid metabolizer": "action",
    "likely poor metabolizer": "action",
    "increased function": "action",
    "poor function": "action",
    "intermediate metabolizer": "review",
    "rapid metabolizer": "review",
    "likely intermediate metabolizer": "review",
    "possible intermediate metabolizer": "review",
    "decreased function": "review",
    "indeterminate": "review",
    "normal metabolizer": "normal",
    "extensive metabolizer": "normal",
    "normal function": "normal",
    "normal": "normal",
    "uncertain susceptibility": "normal",
    "no result": "nodata",
    "n/a": "nodata",
}


def phenotype_to_risk(phenotype: str | None) -> str:
    """Map a free-form phenotype string to a 4-level risk level."""
    if not phenotype:
        return "nodata"
    lower = phenotype.lower().strip()
    for key, risk in PHENOTYPE_RISK.items():
        if key in lower:
            return risk
    if "normal" in lower:
        return "normal"
    return "review"


# ---------------------------------------------------------------------------
# CYP2D6 caveat (best-effort calling from short-read VCFs)
# ---------------------------------------------------------------------------

CYP2D6_CAVEAT = (
    "CYP2D6 results may be unreliable from short-read whole-genome sequencing "
    "due to the gene's complex structural variation (deletions, duplications, "
    "hybrid alleles). Clinical CYP2D6 testing is recommended for actionable decisions."
)


# ---------------------------------------------------------------------------
# Disclaimers
# ---------------------------------------------------------------------------

PATIENT_DISCLAIMER = (
    "This report is for educational and informational purposes only. "
    "It is NOT a substitute for professional medical advice, diagnosis, or treatment. "
    "Always consult a healthcare provider before making any changes to medications. "
    "Results shown are based on computational analysis of genetic data and should be "
    "confirmed by a certified clinical laboratory before any clinical decisions are made."
)

CLINICIAN_DISCLAIMER = (
    "This report was generated using PharmCAT (Pharmacogenomics Clinical Annotation "
    "Tool) developed by PharmGKB / Stanford University. Recommendations are based on "
    "CPIC and DPWG guidelines and FDA-approved drug labels. Results should be "
    "interpreted in the context of the patient's complete clinical picture. Confirm "
    "genotyping results with a CLIA-certified laboratory before making prescribing decisions."
)


# ---------------------------------------------------------------------------
# App identity
# ---------------------------------------------------------------------------

APP_TITLE = "PGx Reporter"
APP_SUBTITLE = "Personalized Pharmacogenomics Report"
