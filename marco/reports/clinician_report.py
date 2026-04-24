"""
Generate a clinician/researcher pharmacogenomics report as HTML and PDF.

Uses the parsed PharmCAT output to render a Jinja2 template with
full technical detail: executive summary, definitive results with drug
annotations, ambiguous gene analysis, no-call genes, all drug annotations,
methodology notes, references, and disclaimers.
"""

import os
from datetime import datetime

from jinja2 import Environment, FileSystemLoader

from config.settings import CLINICIAN_DISCLAIMER, TEMPLATES_DIR
from pharmcat_wrapper.output_parser import parse_full_report


# ---------------------------------------------------------------------------
# Urgency color mapping for template
# ---------------------------------------------------------------------------

URGENCY_COLORS = {
    "red": "#c62828",
    "yellow": "#f9a825",
    "green": "#2e7d32",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_clinician_report(
    report_json_path: str,
    match_json_path: str,
    output_dir: str,
) -> tuple[str, str]:
    """Generate clinician/researcher report as HTML and PDF.

    Args:
        report_json_path: Path to PharmCAT .report.json file.
        match_json_path:  Path to PharmCAT .match.json file.
        output_dir:       Directory where output files will be written.

    Returns:
        (html_path, pdf_path) -- paths to the generated files.
        pdf_path may be None if WeasyPrint is unavailable.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Parse PharmCAT output
    parsed = parse_full_report(report_json_path, match_json_path)

    # Build template context
    context = _build_context(parsed, report_json_path)

    # Render HTML from Jinja2 template
    env = Environment(
        loader=FileSystemLoader(TEMPLATES_DIR),
        autoescape=True,
    )
    template = env.get_template("clinician_report.html")
    html_content = template.render(**context)

    # Write HTML
    html_path = os.path.join(output_dir, "clinician_report.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    # Convert to PDF using WeasyPrint
    pdf_path = os.path.join(output_dir, "clinician_report.pdf")
    try:
        from weasyprint import HTML as WeasyprintHTML
        WeasyprintHTML(string=html_content).write_pdf(pdf_path)
    except (ImportError, OSError):
        pdf_path = _write_pdf_fallback(pdf_path)

    return html_path, pdf_path


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def _build_context(parsed: dict, report_json_path: str) -> dict:
    """Transform parsed PharmCAT data into template context variables."""
    definitive = parsed.get("definitive_genes", [])
    ambiguous = parsed.get("ambiguous_genes", [])
    no_call = parsed.get("no_call_genes", [])
    all_drugs = parsed.get("drugs", [])
    metadata = parsed.get("metadata", {})
    messages = parsed.get("messages", [])

    # Derive sample ID from filename
    sample_id = _derive_sample_id(report_json_path)

    # Filter drugs to CPIC and DPWG sources only (exclude FDA)
    cpic_dpwg_drugs = [
        d for d in all_drugs
        if _is_cpic_or_dpwg(d.get("source", ""))
    ]

    # Build drug-to-gene mapping for the drug annotations table
    for drug in cpic_dpwg_drugs:
        drug["gene_list"] = ", ".join(drug.get("affected_genes", []))
        phenotypes = drug.get("phenotypes", {})
        drug["phenotype_list"] = "; ".join(
            f"{gene}: {pheno}" for gene, pheno in sorted(phenotypes.items())
        )
        drug["urgency_color"] = URGENCY_COLORS.get(drug.get("urgency", "green"), "#2e7d32")

    # Sort drugs by source then drug name
    cpic_dpwg_drugs.sort(key=lambda d: (d.get("source", ""), d.get("drug", "")))

    # Build per-gene drug recommendations for definitive genes
    gene_drug_map = _build_gene_drug_map(cpic_dpwg_drugs)

    # Deduplicate citations
    all_citations = _collect_citations(all_drugs)

    # Coverage summary
    total_positions_found = sum(
        g.get("positions_found", 0) for g in definitive
    ) + sum(
        g.get("positions_found", 0) for g in ambiguous
    )
    total_positions_missing = sum(
        g.get("positions_missing", 0) for g in definitive
    ) + sum(
        g.get("positions_missing", 0) for g in ambiguous
    ) + sum(
        g.get("positions_missing", 0) for g in no_call
    )

    coverage_summary = {
        "total_genes": len(definitive) + len(ambiguous) + len(no_call),
        "definitive_count": len(definitive),
        "ambiguous_count": len(ambiguous),
        "no_call_count": len(no_call),
        "total_positions_found": total_positions_found,
        "total_positions_missing": total_positions_missing,
    }

    return {
        "title": "Pharmacogenomic Analysis Report",
        "date": datetime.now().strftime("%B %d, %Y"),
        "sample_id": sample_id,
        "disclaimer": CLINICIAN_DISCLAIMER,
        "metadata": metadata,
        "definitive_genes": definitive,
        "ambiguous_genes": ambiguous,
        "no_call_genes": no_call,
        "all_drugs": cpic_dpwg_drugs,
        "all_citations": all_citations,
        "messages": messages,
        "coverage_summary": coverage_summary,
        "gene_drug_map": gene_drug_map,
        "urgency_colors": URGENCY_COLORS,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _derive_sample_id(path: str) -> str:
    """Extract sample ID from the report JSON filename."""
    basename = os.path.basename(path)
    # Remove common suffixes
    for suffix in (".report.json", ".match.json", ".json"):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    return basename


def _is_cpic_or_dpwg(source: str) -> bool:
    """Check if a drug source is CPIC or DPWG (not FDA)."""
    source_lower = source.lower()
    return "cpic" in source_lower or "dpwg" in source_lower


def _build_gene_drug_map(drugs: list) -> dict:
    """Build a mapping of gene -> list of drug recommendations."""
    gene_map = {}
    for drug in drugs:
        for gene in drug.get("affected_genes", []):
            if gene not in gene_map:
                gene_map[gene] = []
            gene_map[gene].append(drug)
    return gene_map


def _collect_citations(drugs: list) -> list:
    """Collect and deduplicate all citations from drug recommendations."""
    seen_pmids = set()
    citations = []
    for drug in drugs:
        for cite in drug.get("citations", []):
            pmid = cite.get("pmid", "")
            if pmid and pmid not in seen_pmids:
                seen_pmids.add(pmid)
                citations.append(cite)
            elif not pmid:
                # Include citations without PMID but deduplicate by title
                title = cite.get("title", "")
                if title and title not in seen_pmids:
                    seen_pmids.add(title)
                    citations.append(cite)
    # Sort by year descending, then title
    citations.sort(key=lambda c: (-(c.get("year") or 0), c.get("title", "")))
    return citations


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
