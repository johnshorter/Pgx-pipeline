"""
Generate a patient-friendly pharmacogenomics report as HTML and PDF.

Uses the parsed PharmCAT output to render a Jinja2 template with
traffic-light gene results, inconclusive warnings, and plain-language
explanations suitable for patients.
"""

import os
from datetime import datetime

from jinja2 import Environment, FileSystemLoader

from config.settings import PATIENT_DISCLAIMER, TEMPLATES_DIR
from pharmcat_wrapper.output_parser import parse_full_report

# ---------------------------------------------------------------------------
# Plain-language gene descriptions (what the gene does)
# ---------------------------------------------------------------------------

GENE_DESCRIPTIONS = {
    "CYP2C19": "helps your body break down many medications including blood thinners, antidepressants, and stomach acid drugs",
    "CYP2C9": "helps your body process pain medications, blood thinners like warfarin, and anti-inflammatory drugs",
    "CYP2D6": "is involved in processing about 25% of all medications, including many antidepressants and pain medications",
    "CYP3A4": "is the most important drug-metabolizing enzyme, involved in processing nearly half of all medications",
    "CYP3A5": "helps process immunosuppressants like tacrolimus and some blood pressure medications",
    "CYP2B6": "helps process certain HIV medications and the smoking cessation drug bupropion",
    "CYP4F2": "affects how your body uses vitamin K, which influences warfarin dosing",
    "DPYD": "determines how your body handles fluoropyrimidine chemotherapy drugs",
    "TPMT": "affects how your body processes thiopurine medications used for cancer and autoimmune conditions",
    "NUDT15": "also affects thiopurine medication processing, similar to TPMT",
    "VKORC1": "controls your sensitivity to warfarin, a common blood thinner",
    "SLCO1B1": "affects how your body transports statin medications used to lower cholesterol",
    "UGT1A1": "helps your body process certain medications and bilirubin (a waste product)",
    "ABCG2": "is a transporter protein that affects how your body absorbs certain medications",
    "CACNA1S": "is related to your response to certain anesthetics used during surgery",
    "CFTR": "is related to cystic fibrosis and affects certain medication responses",
    "G6PD": "affects your red blood cells' sensitivity to certain medications",
    "HLA-A": "is part of your immune system and affects risk of severe drug reactions",
    "HLA-B": "is part of your immune system and affects risk of severe drug reactions to specific drugs",
    "IFNL3": "affects your response to hepatitis C treatments",
    "MT-RNR1": "affects your sensitivity to aminoglycoside antibiotics that can cause hearing loss",
    "NAT2": "helps your body process certain medications including some antibiotics",
    "RYR1": "is related to your risk of a rare but serious reaction to anesthetics called malignant hyperthermia",
    "F2": "affects blood clotting and may influence response to anticoagulant medications",
    "F5": "affects blood clotting and may influence response to anticoagulant medications",
}

# ---------------------------------------------------------------------------
# Plain-language phenotype explanations
# ---------------------------------------------------------------------------

PHENOTYPE_PLAIN = {
    "Normal Metabolizer": "Your body processes this medication normally. Standard doses should work as expected.",
    "Intermediate Metabolizer": "Your body processes this medication somewhat differently. Your doctor may need to adjust your dose.",
    "Poor Metabolizer": "Your body has difficulty processing this medication. Alternative drugs or significant dose changes may be needed.",
    "Ultrarapid Metabolizer": "Your body processes this medication much faster than normal. Standard doses may be too low or the drug may not work well.",
    "Rapid Metabolizer": "Your body processes this medication faster than normal. Dose adjustments may be considered.",
    "Increased Function": "This gene's activity is higher than normal, which may change how certain medications affect you.",
    "Decreased Function": "This gene's activity is lower than normal, which may change how certain medications affect you.",
    "Normal Function": "This gene works normally. Standard medication doses should be appropriate.",
    "Poor Function": "This gene has significantly reduced activity. Medications affected by this gene may need alternatives or dose changes.",
    "-1639 AA": "You have increased sensitivity to warfarin. If prescribed warfarin, you may need a lower dose.",
}

# ---------------------------------------------------------------------------
# Short protein-type labels for the overview cards
# ---------------------------------------------------------------------------

GENE_PROTEIN_TYPE = {
    "CYP2C19": "Liver enzyme",
    "CYP2C9": "Liver enzyme",
    "CYP2D6": "Liver enzyme",
    "CYP3A4": "Liver enzyme",
    "CYP3A5": "Liver enzyme",
    "CYP2B6": "Liver enzyme",
    "CYP4F2": "Liver enzyme",
    "DPYD": "Metabolic enzyme",
    "TPMT": "Metabolic enzyme",
    "NUDT15": "Metabolic enzyme",
    "NAT2": "Metabolic enzyme",
    "UGT1A1": "Metabolic enzyme",
    "G6PD": "Red-blood-cell enzyme",
    "VKORC1": "Drug target",
    "SLCO1B1": "Drug transporter",
    "ABCG2": "Drug transporter",
    "CACNA1S": "Ion channel",
    "RYR1": "Ion channel",
    "CFTR": "Ion channel",
    "HLA-A": "Immune marker",
    "HLA-B": "Immune marker",
    "IFNL3": "Immune signaling",
    "MT-RNR1": "Mitochondrial RNA",
    "F2": "Clotting factor",
    "F5": "Clotting factor",
}

# ---------------------------------------------------------------------------
# Accessibility symbols for color-blind users
# ---------------------------------------------------------------------------

# Each status gets a Unicode symbol that conveys meaning without color.
_STATUS_SYMBOLS = {
    "green": "\u2714",   # ✔ Heavy check mark
    "orange": "\u25B2",  # ▲ Triangle (caution)
    "red": "\u2716",     # ✖ Heavy X mark
    "grey": "\u2014",    # — Em dash (unknown)
}

# Traffic-light color codes
_GREEN = "#4CAF50"
_ORANGE = "#FF9800"
_RED = "#F44336"
_GREY = "#9E9E9E"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_patient_report(
    report_json_path: str,
    match_json_path: str,
    output_dir: str,
) -> tuple[str, str]:
    """Generate patient-friendly report as HTML and PDF.

    Args:
        report_json_path: Path to PharmCAT .report.json file.
        match_json_path:  Path to PharmCAT .match.json file.
        output_dir:       Directory where output files will be written.

    Returns:
        (html_path, pdf_path) -- paths to the generated files.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Parse PharmCAT output
    parsed = parse_full_report(report_json_path, match_json_path)

    # Build template context
    context = _build_context(parsed)

    # Render HTML from Jinja2 template
    env = Environment(
        loader=FileSystemLoader(TEMPLATES_DIR),
        autoescape=True,
    )
    template = env.get_template("patient_report.html")
    html_content = template.render(**context)

    # Write HTML
    html_path = os.path.join(output_dir, "patient_report.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    # Convert to PDF using WeasyPrint
    pdf_path = os.path.join(output_dir, "patient_report.pdf")
    try:
        from weasyprint import HTML as WeasyprintHTML
        WeasyprintHTML(string=html_content).write_pdf(pdf_path)
    except (ImportError, OSError):
        # WeasyPrint not installed or missing system libraries (GTK/Pango)
        pdf_path = _write_pdf_fallback(pdf_path)

    return html_path, pdf_path


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def _build_context(parsed: dict) -> dict:
    """Transform parsed PharmCAT data into template context variables."""
    definitive = parsed.get("definitive_genes", [])
    ambiguous = parsed.get("ambiguous_genes", [])
    no_call = parsed.get("no_call_genes", [])

    # Filter definitive genes: only those with real results
    display_genes = [
        g for g in definitive
        if g.get("phenotype") not in ("No Result", "Uncertain Susceptibility")
    ]

    # Traffic-light counts
    green_count = sum(1 for g in display_genes if g.get("color") == _GREEN)
    orange_count = sum(1 for g in display_genes if g.get("color") == _ORANGE)
    red_count = sum(1 for g in display_genes if g.get("color") == _RED)
    grey_count = len(ambiguous) + len(no_call)

    # Summary sentence
    total = green_count + orange_count + red_count + grey_count
    summary_parts = []
    if green_count:
        summary_parts.append(
            f"{green_count} gene{'s' if green_count != 1 else ''} "
            f"show{'s' if green_count == 1 else ''} normal results"
        )
    if orange_count or red_count:
        action = orange_count + red_count
        summary_parts.append(
            f"{action} gene{'s' if action != 1 else ''} "
            f"may affect your medication plan"
        )
    if grey_count:
        summary_parts.append(
            f"{grey_count} gene{'s' if grey_count != 1 else ''} "
            f"could not be fully determined"
        )
    summary_sentence = (
        ". ".join(summary_parts) + "."
        if summary_parts
        else f"Your test analyzed {total} genes."
    )

    # Filter no-call genes to those with affected drugs
    no_call_with_drugs = [
        g for g in no_call if g.get("affected_drugs")
    ]

    # Build overview cards for ALL genes
    overview_cards = _build_overview_cards(definitive, ambiguous, no_call)

    return {
        "date_generated": datetime.now().strftime("%B %d, %Y"),
        "disclaimer": PATIENT_DISCLAIMER,
        "metadata": parsed.get("metadata", {}),
        "counts": {
            "green": green_count,
            "orange": orange_count,
            "red": red_count,
            "grey": grey_count,
        },
        "summary_sentence": summary_sentence,
        "overview_cards": overview_cards,
        "definitive_genes": display_genes,
        "ambiguous_genes": ambiguous,
        "no_call_genes_with_drugs": no_call_with_drugs,
        "gene_descriptions": GENE_DESCRIPTIONS,
        "phenotype_plain": PHENOTYPE_PLAIN,
    }


# ---------------------------------------------------------------------------
# Overview card builder
# ---------------------------------------------------------------------------

def _color_to_status(color: str) -> str:
    """Map a hex color to a status key."""
    return {
        _GREEN: "green",
        _ORANGE: "orange",
        _RED: "red",
    }.get(color, "grey")


def _build_overview_cards(
    definitive: list, ambiguous: list, no_call: list,
) -> list:
    """
    Build a flat list of overview-card dicts for every gene, suitable for
    the card-grid in the template.

    Each card has: gene, status, color, symbol, status_label, protein_type,
    phenotype_short, med_count.
    """
    cards = []

    for g in definitive:
        status = _color_to_status(g.get("color", _GREY))
        phenotype = g.get("phenotype", "No Result")
        med_count = len(g.get("related_drugs", []))

        cards.append({
            "gene": g["gene"],
            "status": status,
            "color": g.get("color", _GREY),
            "symbol": _STATUS_SYMBOLS.get(status, "\u2014"),
            "status_label": _status_label(status),
            "protein_type": GENE_PROTEIN_TYPE.get(g["gene"], "Gene"),
            "phenotype_short": _shorten_phenotype(phenotype),
            "med_count": med_count,
        })

    for g in ambiguous:
        med_count = len(g.get("related_drugs", []))
        cards.append({
            "gene": g["gene"],
            "status": "orange",
            "color": _ORANGE,
            "symbol": "?",
            "status_label": "Inconclusive",
            "protein_type": GENE_PROTEIN_TYPE.get(g["gene"], "Gene"),
            "phenotype_short": "Inconclusive",
            "med_count": med_count,
        })

    for g in no_call:
        med_count = len(g.get("related_drugs", []))
        cards.append({
            "gene": g["gene"],
            "status": "grey",
            "color": _GREY,
            "symbol": _STATUS_SYMBOLS["grey"],
            "status_label": "Not tested",
            "protein_type": GENE_PROTEIN_TYPE.get(g["gene"], "Gene"),
            "phenotype_short": "Not tested",
            "med_count": med_count,
        })

    # Sort: red first, then orange, then green, then grey; alphabetically within
    priority = {"red": 0, "orange": 1, "green": 2, "grey": 3}
    cards.sort(key=lambda c: (priority.get(c["status"], 3), c["gene"]))

    return cards


def _status_label(status: str) -> str:
    """Human-readable label for the card status badge."""
    return {
        "green": "Normal",
        "orange": "Caution",
        "red": "Action needed",
        "grey": "Not tested",
    }.get(status, "Unknown")


def _shorten_phenotype(phenotype: str) -> str:
    """Shorten a phenotype for the overview card."""
    short = {
        "Normal Metabolizer": "Normal",
        "Intermediate Metabolizer": "Intermediate",
        "Poor Metabolizer": "Poor",
        "Ultrarapid Metabolizer": "Ultrarapid",
        "Rapid Metabolizer": "Rapid",
        "Normal Function": "Normal",
        "Increased Function": "Increased",
        "Decreased Function": "Decreased",
        "Poor Function": "Poor",
        "No Result": "No result",
        "Uncertain Susceptibility": "Uncertain",
        "-1639 AA": "High sensitivity",
    }
    return short.get(phenotype, phenotype)


# ---------------------------------------------------------------------------
# PDF fallback when WeasyPrint is not installed
# ---------------------------------------------------------------------------

def _write_pdf_fallback(pdf_path: str) -> str:
    """Write a minimal text file explaining WeasyPrint is needed."""
    fallback_path = pdf_path.replace(".pdf", ".pdf_not_generated.txt")
    with open(fallback_path, "w", encoding="utf-8") as f:
        f.write(
            "PDF generation requires WeasyPrint.\n"
            "Install it with: pip install weasyprint\n"
            "Then re-run the report generation.\n"
        )
    return fallback_path
