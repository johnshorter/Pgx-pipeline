"""
Clinician report generator.

Presents full technical detail: star alleles, diplotypes, CPIC guideline
citations, evidence strength, dosing recommendations, and publication links.
Grouped by gene functional category with colorblind-safe design.
"""

from collections import OrderedDict
from datetime import date
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from src.pharmcat.output_parser import (
    ParsedResults, GENE_CATEGORIES, ACTION_SYMBOLS,
)

TEMPLATES_DIR = Path(__file__).parent.parent / "templates"

CATEGORY_ORDER = list(GENE_CATEGORIES.keys())


def generate_clinician_report(results: ParsedResults, output_dir: str | Path) -> dict[str, Path]:
    """Generate the clinician-facing technical report in HTML and PDF."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    action_order = {"action": 0, "review": 1, "normal": 2, "nodata": 3}

    gene_data = []
    for g in results.gene_results:
        gene_data.append({
            "gene": g.gene,
            "diplotype": g.diplotype,
            "phenotype": g.phenotype,
            "activity_score": g.activity_score,
            "star_alleles": g.star_alleles,
            "risk_level": g.risk_level,
            "caveat": g.caveat,
            "has_data": g.has_data,
            "category": g.category,
            "category_desc": g.category_desc,
            "symbol": ACTION_SYMBOLS.get(g.risk_level, "\u2014"),
        })

    # Group by category
    categories: OrderedDict[str, list[dict]] = OrderedDict()
    for cat_name in CATEGORY_ORDER:
        categories[cat_name] = []
    for gene in gene_data:
        cat = gene.get("category", "Other Pharmacogenes")
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(gene)
    categories = OrderedDict((k, v) for k, v in categories.items() if v)

    # Sort within each category: action first
    for cat_name in categories:
        categories[cat_name].sort(key=lambda g: action_order.get(g["risk_level"], 2))

    drug_data = []
    for d in results.drug_recommendations:
        drug_data.append({
            "drug": d.drug,
            "gene": d.gene,
            "guideline_source": d.guideline_source,
            "recommendation": d.recommendation,
            "classification": d.classification,
            "implications": d.implications,
            "phenotype": d.phenotype,
            "risk_level": d.risk_level,
            "url": d.url,
            "symbol": ACTION_SYMBOLS.get(d.risk_level, "\u2014"),
        })
    drug_data.sort(key=lambda d: action_order.get(d["risk_level"], 2))

    # Flat sorted list for summary table
    all_genes_flat = []
    for cat_name, genes in categories.items():
        all_genes_flat.extend(genes)

    context = {
        "report_date": date.today().strftime("%B %d, %Y"),
        "sample_id": results.sample_id,
        "gene_results": all_genes_flat,
        "categories": categories,
        "category_descs": {k: v["description"] for k, v in GENE_CATEGORIES.items()},
        "drug_recommendations": drug_data,
        "missing_genes": results.missing_genes,
        "action_symbols": ACTION_SYMBOLS,
    }

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=True)
    template = env.get_template("clinician_report.html")
    html_content = template.render(**context)

    html_path = output_dir / "clinician_report.html"
    html_path.write_text(html_content, encoding="utf-8")

    output_files = {"html": html_path}

    pdf_path = output_dir / "clinician_report.pdf"
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
