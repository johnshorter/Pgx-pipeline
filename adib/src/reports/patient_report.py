"""
Patient report generator.

Translates PGx data into plain language with a colorblind-safe design system.
Uses shapes/symbols alongside colors and groups genes by functional category.
"""

from collections import OrderedDict
from datetime import date
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from src.pharmcat.output_parser import (
    ParsedResults, GeneResult, DrugRecommendation,
    GENE_CATEGORIES, ACTION_SYMBOLS,
)

TEMPLATES_DIR = Path(__file__).parent.parent / "templates"

# Plain-language explanations for phenotypes
PHENOTYPE_EXPLANATIONS = {
    "normal metabolizer": {
        "brief": "Your body processes medications involving this gene at a typical rate.",
        "detail": (
            "You carry two normally functioning copies of this gene. Your body is expected "
            "to process medications involving this gene at a standard rate. Standard dosing "
            "recommendations typically apply."
        ),
    },
    "intermediate metabolizer": {
        "brief": "Your body may process some medications slightly slower than average.",
        "detail": (
            "You carry gene variants that may result in somewhat reduced enzyme activity. "
            "This means your body might process certain medications a bit more slowly than "
            "average. Your doctor may want to consider this when prescribing affected medications."
        ),
    },
    "poor metabolizer": {
        "brief": "Your body processes some medications much slower than average — discuss with your doctor.",
        "detail": (
            "You carry gene variants associated with significantly reduced or absent enzyme "
            "activity. This means certain medications may build up in your body more than "
            "expected, potentially increasing the risk of side effects. Your doctor may need "
            "to adjust doses or choose alternative medications."
        ),
    },
    "rapid metabolizer": {
        "brief": "Your body processes some medications faster than average.",
        "detail": (
            "You carry gene variants associated with increased enzyme activity. This means "
            "your body may break down certain medications faster than average, which could "
            "reduce their effectiveness at standard doses. Your doctor may need to consider "
            "dose adjustments."
        ),
    },
    "ultrarapid metabolizer": {
        "brief": "Your body processes some medications much faster than average — discuss with your doctor.",
        "detail": (
            "You carry gene variants associated with significantly increased enzyme activity. "
            "Standard doses of certain medications may be broken down too quickly to be "
            "effective, or in some cases, faster metabolism can lead to dangerously "
            "high levels of active drug. This is an important finding to discuss with your "
            "healthcare provider."
        ),
    },
    "indeterminate": {
        "brief": "Your result could not be clearly determined — discuss with your doctor.",
        "detail": (
            "The analysis was unable to clearly determine your metabolizer status for this gene. "
            "This may be due to uncommon genetic variants or limitations of the testing method. "
            "Your healthcare provider can help determine if additional testing is warranted."
        ),
    },
    "normal function": {
        "brief": "This gene is functioning normally. No medication adjustments expected.",
        "detail": (
            "Your genetic result indicates normal function for this gene. Drug transport "
            "and processing related to this gene are expected to occur at standard rates."
        ),
    },
    "normal": {
        "brief": "Your result for this gene is within the normal range.",
        "detail": (
            "Your genetic result indicates normal activity for this gene. No special "
            "medication considerations are expected based on this finding."
        ),
    },
    "decreased function": {
        "brief": "This gene may have reduced activity — discuss with your doctor.",
        "detail": (
            "Your genetic result suggests decreased function of this gene. This may affect "
            "how certain medications are transported or processed in your body. Your doctor "
            "may want to consider this when prescribing affected medications."
        ),
    },
    "increased function": {
        "brief": "This gene may have increased activity — discuss with your doctor.",
        "detail": (
            "Your genetic result suggests increased function of this gene. This may affect "
            "how certain medications are transported or processed in your body."
        ),
    },
    "uncertain susceptibility": {
        "brief": "Your susceptibility status for this gene is within the typical range.",
        "detail": (
            "Your genetic result does not indicate a known increased susceptibility "
            "associated with this gene. Standard precautions apply."
        ),
    },
    "no result": {
        "brief": "No result could be determined for this gene from the available data.",
        "detail": (
            "This gene could not be analyzed from the available genetic data. This may be "
            "due to limitations of the sequencing method or the gene requiring specialized "
            "analysis tools. Consider clinical testing if results for this gene are needed."
        ),
    },
}

# Patient-friendly drug guidance by action level
DRUG_GUIDANCE = {
    "action": "Important: Talk to your doctor before taking this medication. Your genetics suggest it may need special consideration.",
    "review": "Note: Discuss this medication with your doctor. A dose adjustment or monitoring may be helpful.",
    "normal": "Standard use is expected to be appropriate based on your genetics.",
    "nodata": "Insufficient data to provide guidance for this medication.",
}

# Category display order
CATEGORY_ORDER = list(GENE_CATEGORIES.keys())


def generate_patient_report(results: ParsedResults, output_dir: str | Path) -> dict[str, Path]:
    """Generate the patient-friendly report in HTML and PDF formats."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_data = [_enrich_gene_for_patient(g) for g in results.gene_results]

    # Group by category, maintaining category order
    categories = _group_by_category(gene_data)

    drug_data = [_enrich_drug_for_patient(d) for d in results.drug_recommendations]
    action_order = {"action": 0, "review": 1, "normal": 2, "nodata": 3}
    drug_data.sort(key=lambda d: action_order.get(d["risk_level"], 2))

    has_action = any(g["risk_level"] == "action" for g in gene_data)

    # Flat list for summary table (sorted: action first within each category)
    all_genes_flat = []
    for cat_name, genes in categories.items():
        sorted_genes = sorted(genes, key=lambda g: action_order.get(g["risk_level"], 2))
        all_genes_flat.extend(sorted_genes)

    context = {
        "report_date": date.today().strftime("%B %d, %Y"),
        "sample_id": results.sample_id,
        "gene_results": all_genes_flat,
        "categories": categories,
        "drug_recommendations": drug_data,
        "has_action_results": has_action,
        "action_symbols": ACTION_SYMBOLS,
    }

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=True)
    template = env.get_template("patient_report.html")
    html_content = template.render(**context)

    html_path = output_dir / "patient_report.html"
    html_path.write_text(html_content, encoding="utf-8")

    output_files = {"html": html_path}

    pdf_path = output_dir / "patient_report.pdf"
    try:
        from weasyprint import HTML
        HTML(string=html_content).write_pdf(str(pdf_path))
        output_files["pdf"] = pdf_path
    except ImportError:
        pass
    except Exception as e:
        import logging
        logging.getLogger(__name__).warning("PDF generation failed: %s", e)

    return output_files


def _group_by_category(gene_data: list[dict]) -> OrderedDict:
    """Group gene results by their functional category."""
    categories: OrderedDict[str, list[dict]] = OrderedDict()
    for cat_name in CATEGORY_ORDER:
        categories[cat_name] = []

    for gene in gene_data:
        cat = gene.get("category", "Other Pharmacogenes")
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(gene)

    # Remove empty categories
    return OrderedDict((k, v) for k, v in categories.items() if v)


def _enrich_gene_for_patient(gene: GeneResult) -> dict:
    """Add plain-language fields to a gene result for the patient template."""
    phenotype_lower = gene.phenotype.lower().strip()

    explanations = None
    for key, val in PHENOTYPE_EXPLANATIONS.items():
        if key in phenotype_lower:
            explanations = val
            break

    if explanations is None:
        explanations = PHENOTYPE_EXPLANATIONS.get("indeterminate", {
            "brief": "Discuss this result with your healthcare provider.",
            "detail": "Your healthcare provider can help interpret this result.",
        })

    return {
        "gene": gene.gene,
        "diplotype": gene.diplotype,
        "phenotype": gene.phenotype,
        "activity_score": gene.activity_score,
        "star_alleles": gene.star_alleles,
        "risk_level": gene.risk_level,
        "caveat": gene.caveat,
        "plain_language": explanations["brief"],
        "detailed_explanation": explanations["detail"],
        "category": gene.category,
        "category_desc": gene.category_desc,
        "symbol": ACTION_SYMBOLS.get(gene.risk_level, "\u2014"),
    }


def _enrich_drug_for_patient(drug: DrugRecommendation) -> dict:
    """Add patient-friendly guidance to a drug recommendation."""
    return {
        "drug": drug.drug,
        "gene": drug.gene,
        "risk_level": drug.risk_level,
        "patient_guidance": DRUG_GUIDANCE.get(drug.risk_level, DRUG_GUIDANCE["review"]),
        "symbol": ACTION_SYMBOLS.get(drug.risk_level, "\u2014"),
    }
