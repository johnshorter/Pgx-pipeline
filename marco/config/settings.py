"""
Centralized configuration for the PGx Reporter application.

All paths, gene lists, and application settings are defined here
so they can be easily modified without changing code throughout the project.
"""

import os

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# Root directory of the pgx-reporter project
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# PharmCAT JAR location
PHARMCAT_JAR = os.path.join(PROJECT_ROOT, "lib", "pharmcat.jar")

# Java executable - auto-detect from PATH, fall back to known install location
JAVA_EXECUTABLE = "java"
_ADOPTIUM_PATH = os.path.join(
    os.environ.get("ProgramFiles", r"C:\Program Files"),
    "Eclipse Adoptium",
)
if os.path.isdir(_ADOPTIUM_PATH):
    # Pick the first JDK directory found
    for _entry in sorted(os.listdir(_ADOPTIUM_PATH), reverse=True):
        _candidate = os.path.join(_ADOPTIUM_PATH, _entry, "bin", "java.exe")
        if os.path.isfile(_candidate):
            JAVA_EXECUTABLE = _candidate
            break

# Data files directory
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

# Temporary directory for PharmCAT output
TEMP_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "temp_output")

# Report templates directory
TEMPLATES_DIR = os.path.join(PROJECT_ROOT, "reports", "templates")

# ---------------------------------------------------------------------------
# PharmCAT Settings
# ---------------------------------------------------------------------------

# Maximum time (seconds) to wait for PharmCAT to finish processing.
# Whole-genome VCFs (1-3 GB) can take 30-60 minutes on typical laptops
# because PharmCAT reads the entire file looking for PGx positions.
# Override with the --timeout CLI flag for specific runs.
PHARMCAT_TIMEOUT = 3600  # 60 minutes default

# Maximum VCF file size allowed for upload (in MB)
MAX_VCF_SIZE_MB = 3000  # 3 GB to accommodate the HG005 benchmark file

# ---------------------------------------------------------------------------
# Target Pharmacogenes
# ---------------------------------------------------------------------------

# Genes we report on, with GRCh38 coordinates for reference
# These are the most clinically important pharmacogenes
TARGET_GENES = {
    "CYP2D6": {
        "chromosome": "chr22",
        "start": 42126499,
        "end": 42130881,
        "description": "Drug metabolism enzyme",
        "flag": "no_data",  # No variants in HG005 VCF for this gene
    },
    "CYP2C19": {
        "chromosome": "chr10",
        "start": 94762681,
        "end": 94855547,
        "description": "Drug metabolism enzyme",
    },
    "CYP2C9": {
        "chromosome": "chr10",
        "start": 94938683,
        "end": 94990091,
        "description": "Drug metabolism enzyme",
    },
    "CYP3A4": {
        "chromosome": "chr7",
        "start": 99756960,
        "end": 99784247,
        "description": "Major drug metabolism enzyme",
    },
    "CYP3A5": {
        "chromosome": "chr7",
        "start": 99648194,
        "end": 99680437,
        "description": "Drug metabolism enzyme",
    },
    "DPYD": {
        "chromosome": "chr1",
        "start": 97543299,
        "end": 98386615,
        "description": "Fluoropyrimidine metabolism",
    },
    "TPMT": {
        "chromosome": "chr6",
        "start": 18128542,
        "end": 18155374,
        "description": "Thiopurine metabolism",
    },
    "VKORC1": {
        "chromosome": "chr16",
        "start": 31096068,
        "end": 31107301,
        "description": "Warfarin sensitivity",
    },
    "SLCO1B1": {
        "chromosome": "chr12",
        "start": 21130388,
        "end": 21241524,
        "description": "Statin transporter",
    },
    "UGT1A1": {
        "chromosome": "chr2",
        "start": 233757013,
        "end": 233773299,
        "description": "Bilirubin/drug metabolism",
    },
}

# ---------------------------------------------------------------------------
# Metabolizer Phenotype Colors (Traffic-Light System)
# ---------------------------------------------------------------------------

PHENOTYPE_COLORS = {
    # Normal function - no action needed
    "Normal Metabolizer": "#4CAF50",          # Green
    "Extensive Metabolizer": "#4CAF50",       # Green (legacy term)
    "Normal Function": "#4CAF50",             # Green

    # Altered function - caution, discuss with clinician
    "Intermediate Metabolizer": "#FF9800",    # Orange/Yellow
    "Rapid Metabolizer": "#FF9800",           # Orange/Yellow
    "Possible Intermediate Metabolizer": "#FF9800",
    "Decreased Function": "#FF9800",          # Orange

    # Significantly altered - action likely needed
    "Poor Metabolizer": "#F44336",            # Red
    "Ultrarapid Metabolizer": "#F44336",      # Red
    "Increased Function": "#F44336",          # Red

    # Unknown or no data
    "Indeterminate": "#9E9E9E",              # Grey
    "No Data": "#9E9E9E",                    # Grey
}

# ---------------------------------------------------------------------------
# Streamlit App Settings
# ---------------------------------------------------------------------------

APP_TITLE = "PGx Reporter"
APP_SUBTITLE = "Personalized Pharmacogenomics Report"
APP_ICON = "\U0001f9ec"  # DNA emoji

# Page layout
PAGE_LAYOUT = "wide"

# ---------------------------------------------------------------------------
# Disclaimer Text
# ---------------------------------------------------------------------------

PATIENT_DISCLAIMER = (
    "This report is for educational and informational purposes only. "
    "It is NOT a substitute for professional medical advice, diagnosis, or treatment. "
    "Always consult your healthcare provider before making any changes to your medications. "
    "The results shown are based on computational analysis of genetic data and should be "
    "confirmed by a certified clinical laboratory before any clinical decisions are made."
)

CLINICIAN_DISCLAIMER = (
    "This report was generated using PharmCAT (Pharmacogenomics Clinical Annotation Tool) "
    "developed by PharmGKB/Stanford University. Recommendations are based on CPIC guidelines, "
    "DPWG guidelines, and FDA-approved drug labels. Results should be interpreted in the context "
    "of the patient's complete clinical picture. Confirm genotyping results with a CLIA-certified "
    "laboratory before making prescribing decisions."
)
