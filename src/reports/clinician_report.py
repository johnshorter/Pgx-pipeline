"""
Clinician report — full technical detail.

Combines:
- Marco's three-bucket layout (definitive / ambiguous / no-call)
- Marco's CPIC+DPWG drug filter (FDA excluded)
- Marco's coverage_summary, gene→drug map, citation collection
- Adib's gene functional categories (each bucket grouped by category)
- Adib's 4-level action symbols + colorblind-safe palette
"""

from collections import OrderedDict, defaultdict
from datetime import datetime
from pathlib import Path

from config.settings import (
    ACTION_LABELS, ACTION_SYMBOLS, APP_TITLE, CATEGORY_ORDER,
    CLINICIAN_DISCLAIMER, GENE_CATEGORIES, RISK_PRIORITY,
)
from pharmcat.output_parser import ParsedResults
from reports._render import html_and_pdf


def generate_clinician_report(
    parsed: ParsedResults,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Render the clinician report (HTML + best-effort PDF)."""
    context = _build_context(parsed)
    return html_and_pdf("clinician_report.html", context, Path(output_dir), "clinician_report")


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def _build_context(parsed: ParsedResults) -> dict:
    coverage = _coverage_summary(parsed)
    cpic_dpwg = _filter_cpic_dpwg(parsed)
    gene_drug_map = _build_gene_drug_map(cpic_dpwg)

    # Bucket-by-category groupings so the template can show all three
    # buckets organised by functional category.
    definitive_by_cat = _group_by_category(
        [_enrich_drugs(_enrich_definitive(g), gene_drug_map) for g in parsed.definitive_genes]
    )
    ambiguous_by_cat = _group_by_category(
        [_enrich_ambiguous(g) for g in parsed.ambiguous_genes]
    )
    no_call_by_cat = _group_by_category(
        [_enrich_no_call(g) for g in parsed.no_call_genes]
    )

    drugs_by_category = _group_drugs_by_category(cpic_dpwg)

    return {
        "title": "Pharmacogenomic Analysis Report",
        "app_title": APP_TITLE,
        "report_date": datetime.now().strftime("%B %d, %Y"),
        "sample_id": parsed.sample_id or "Unknown",
        "metadata": parsed.metadata,
        "disclaimer": CLINICIAN_DISCLAIMER,
        "messages": parsed.messages,
        "action_symbols": ACTION_SYMBOLS,
        "action_labels": ACTION_LABELS,
        "coverage_summary": coverage,
        "definitive_by_category": definitive_by_cat,
        "ambiguous_by_category": ambiguous_by_cat,
        "no_call_by_category": no_call_by_cat,
        "drugs_by_category": drugs_by_category,
        "all_drugs": cpic_dpwg,
        "citations": parsed.citations,
        "gene_drug_map": gene_drug_map,
        "category_order": CATEGORY_ORDER,
        "category_descriptions": {
            k: v["description"] for k, v in GENE_CATEGORIES.items()
        },
    }


def _coverage_summary(parsed: ParsedResults) -> dict:
    total_found = (
        sum(g.positions_found for g in parsed.definitive_genes)
        + sum(g.positions_found for g in parsed.ambiguous_genes)
    )
    total_missing = (
        sum(g.positions_missing for g in parsed.definitive_genes)
        + sum(g.positions_missing for g in parsed.ambiguous_genes)
        + sum(g.positions_missing for g in parsed.no_call_genes)
    )
    return {
        "total_genes": (
            len(parsed.definitive_genes)
            + len(parsed.ambiguous_genes)
            + len(parsed.no_call_genes)
        ),
        "definitive_count": len(parsed.definitive_genes),
        "ambiguous_count": len(parsed.ambiguous_genes),
        "no_call_count": len(parsed.no_call_genes),
        "total_positions_found": total_found,
        "total_positions_missing": total_missing,
    }


# ---------------------------------------------------------------------------
# Gene enrichment
# ---------------------------------------------------------------------------

def _enrich_definitive(g) -> dict:
    return {
        "gene": g.gene,
        "diplotype": g.diplotype,
        "phenotype": g.phenotype,
        "activity_score": g.activity_score,
        "star_alleles": g.star_alleles,
        "allele1_name": g.allele1_name,
        "allele1_function": g.allele1_function,
        "allele2_name": g.allele2_name,
        "allele2_function": g.allele2_function,
        "risk_level": g.risk_level,
        "symbol": ACTION_SYMBOLS[g.risk_level],
        "label": ACTION_LABELS[g.risk_level],
        "category": g.category,
        "category_desc": g.category_desc,
        "positions_found": g.positions_found,
        "positions_missing": g.positions_missing,
        "related_drugs": g.related_drugs,
        "caveat": g.caveat,
        "drug_recs": [],
    }


def _enrich_ambiguous(g) -> dict:
    return {
        "gene": g.gene,
        "diplotype_count": g.diplotype_count,
        "phenotype_range": g.phenotype_range,
        "harmful_phenotypes": g.harmful_phenotypes,
        "risk_level": g.risk_level,
        "symbol": ACTION_SYMBOLS[g.risk_level],
        "label": ACTION_LABELS[g.risk_level],
        "category": g.category,
        "category_desc": g.category_desc,
        "positions_found": g.positions_found,
        "positions_missing": g.positions_missing,
        "actionable_drugs": g.actionable_drugs,
        "related_drugs": g.related_drugs,
        "caveat": g.caveat,
    }


def _enrich_no_call(g) -> dict:
    return {
        "gene": g.gene,
        "risk_level": "nodata",
        "symbol": ACTION_SYMBOLS["nodata"],
        "label": ACTION_LABELS["nodata"],
        "category": g.category,
        "category_desc": g.category_desc,
        "positions_missing": g.positions_missing,
        "affected_drugs": g.affected_drugs,
        "related_drugs": g.related_drugs,
        "caveat": g.caveat,
    }


def _enrich_drugs(gene_dict: dict, gene_drug_map: dict) -> dict:
    """Attach the matching CPIC/DPWG drug recommendations to a gene dict."""
    gene_dict["drug_recs"] = gene_drug_map.get(gene_dict["gene"], [])
    return gene_dict


# ---------------------------------------------------------------------------
# Drug helpers
# ---------------------------------------------------------------------------

def _filter_cpic_dpwg(parsed: ParsedResults) -> list[dict]:
    out: list[dict] = []
    for d in parsed.drugs:
        source_lower = d.source.lower()
        if "cpic" not in source_lower and "dpwg" not in source_lower:
            continue
        phenotype_list = "; ".join(
            f"{gene}: {pheno}" for gene, pheno in sorted(d.phenotypes.items())
        )
        out.append({
            "drug": d.drug,
            "source": d.source,
            "recommendation": d.recommendation,
            "classification": d.classification,
            "implications": d.implications,
            "affected_genes": d.affected_genes,
            "gene_list": ", ".join(d.affected_genes),
            "phenotypes": d.phenotypes,
            "phenotype_list": phenotype_list,
            "population": d.population,
            "risk_level": d.risk_level,
            "symbol": ACTION_SYMBOLS[d.risk_level],
            "label": ACTION_LABELS[d.risk_level],
            "therapeutic_category": d.therapeutic_category,
            "urls": d.urls,
            "citations": d.citations,
        })
    out.sort(key=lambda d: (RISK_PRIORITY.get(d["risk_level"], 3), d["source"], d["drug"]))
    return out


def _build_gene_drug_map(drugs: list[dict]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    for d in drugs:
        for gene in d["affected_genes"]:
            out[gene].append(d)
    return dict(out)


def _group_drugs_by_category(drugs: list[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for d in drugs:
        grouped[d["therapeutic_category"]].append(d)

    out: list[dict] = []
    for cat in sorted(grouped.keys()):
        items = sorted(
            grouped[cat],
            key=lambda d: (RISK_PRIORITY.get(d["risk_level"], 3), d["drug"]),
        )
        out.append({"category": cat, "drugs": items})

    out.sort(
        key=lambda group: min(
            RISK_PRIORITY.get(d["risk_level"], 3) for d in group["drugs"]
        )
    )
    return out


# ---------------------------------------------------------------------------
# Functional-category grouping
# ---------------------------------------------------------------------------

def _group_by_category(genes: list[dict]) -> "OrderedDict[str, dict]":
    cats: OrderedDict[str, dict] = OrderedDict()
    for cat_name in CATEGORY_ORDER:
        cats[cat_name] = {
            "name": cat_name,
            "description": GENE_CATEGORIES[cat_name]["description"],
            "genes": [],
        }
    for g in genes:
        cats.setdefault(g["category"], {
            "name": g["category"],
            "description": "",
            "genes": [],
        })["genes"].append(g)

    for cat in cats.values():
        cat["genes"].sort(
            key=lambda g: (RISK_PRIORITY.get(g["risk_level"], 3), g["gene"])
        )
    return OrderedDict((k, v) for k, v in cats.items() if v["genes"])
