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


# Phenotype strings considered "normal-equivalent" — drugs whose recommendation
# lists the patient's gene at one of these phenotypes don't drive the gene's
# risk classification. Mirrors the same set used by patient_report_v2 so both
# reports bucket genes the same way.
_NORMAL_PHENOTYPES_LOWER = {
    "normal metabolizer",
    "extensive metabolizer",
    "normal function",
    "normal",
    "uncertain susceptibility",
    "n/a",
    "",
    "no result",
}


def _is_normal_pheno(p: str | None) -> bool:
    return (p or "").lower().strip() in _NORMAL_PHENOTYPES_LOWER


def _effective_gene_risk(parsed: ParsedResults) -> dict[str, str]:
    """Each gene's effective risk = max risk among drugs whose recommendation
    this gene's phenotype actually drives. Falls back to the gene's own
    phenotype-derived risk_level when no drugs contribute. Same logic as in
    patient_report_v2 so the two reports agree on which bucket a gene
    belongs to."""
    drug_derived: dict[str, str] = {}
    for d in parsed.drugs:
        for gene in d.affected_genes:
            if _is_normal_pheno(d.phenotypes.get(gene, "")):
                continue
            current = drug_derived.get(gene)
            if (
                current is None
                or RISK_PRIORITY.get(d.risk_level, 3)
                < RISK_PRIORITY.get(current, 3)
            ):
                drug_derived[gene] = d.risk_level

    effective: dict[str, str] = {}
    for g in parsed.definitive_genes:
        effective[g.gene] = drug_derived.get(g.gene, g.risk_level)
    for g in parsed.ambiguous_genes:
        effective[g.gene] = drug_derived.get(g.gene, g.risk_level)
    return effective


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
    effective_risk = _effective_gene_risk(parsed)

    # Bucket-by-category groupings, with per-gene risk now derived from drug
    # recommendations rather than the phenotype name alone.
    definitive_by_cat = _group_by_category([
        _enrich_drugs(_enrich_definitive(g, effective_risk), gene_drug_map)
        for g in parsed.definitive_genes
    ])
    ambiguous_by_cat = _group_by_category(
        [_enrich_ambiguous(g, effective_risk) for g in parsed.ambiguous_genes]
    )
    no_call_by_cat = _group_by_category(
        [_enrich_no_call(g) for g in parsed.no_call_genes]
    )

    drugs_by_category = _group_drugs_by_category(cpic_dpwg)

    clinical_priority = _build_clinical_priority(
        parsed, effective_risk, gene_drug_map,
    )
    consolidated_caveats = _consolidated_caveats(parsed)

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
        "clinical_priority": clinical_priority,
        "definitive_by_category": definitive_by_cat,
        "ambiguous_by_category": ambiguous_by_cat,
        "no_call_by_category": no_call_by_cat,
        "drugs_by_category": drugs_by_category,
        "all_drugs": cpic_dpwg,
        "citations": parsed.citations,
        "gene_drug_map": gene_drug_map,
        "consolidated_caveats": consolidated_caveats,
        "category_order": CATEGORY_ORDER,
        "category_descriptions": {
            k: v["description"] for k, v in GENE_CATEGORIES.items()
        },
    }


def _build_clinical_priority(
    parsed: ParsedResults,
    effective_risk: dict[str, str],
    gene_drug_map: dict[str, list[dict]],
) -> dict:
    """Build the per-risk priority summary at the top of the report: Action
    and Review genes with a one-line headline (phenotype + driven-drug count)
    so a clinician can see at-a-glance which gene results need attention."""

    def _entry(gene: str, phenotype: str) -> dict:
        drugs = gene_drug_map.get(gene, [])
        action_count = sum(1 for d in drugs if d["risk_level"] == "action")
        review_count = sum(1 for d in drugs if d["risk_level"] == "review")
        return {
            "gene": gene,
            "phenotype": phenotype or "—",
            "action_drug_count": action_count,
            "review_drug_count": review_count,
        }

    action: list[dict] = []
    review: list[dict] = []

    for g in parsed.definitive_genes:
        risk = effective_risk.get(g.gene, g.risk_level)
        if risk == "action":
            action.append(_entry(g.gene, g.phenotype))
        elif risk == "review":
            review.append(_entry(g.gene, g.phenotype))

    for g in parsed.ambiguous_genes:
        risk = effective_risk.get(g.gene, g.risk_level)
        pheno_summary = "Inconclusive (" + str(g.diplotype_count) + " candidates)"
        if risk == "action":
            action.append(_entry(g.gene, pheno_summary))
        elif risk == "review":
            review.append(_entry(g.gene, pheno_summary))

    action.sort(key=lambda e: e["gene"])
    review.sort(key=lambda e: e["gene"])

    return {
        "action": action,
        "review": review,
    }


def _consolidated_caveats(parsed: ParsedResults) -> list[str]:
    """Collect unique caveat strings across all gene buckets. The only
    caveat currently used is the CYP2D6 short-read warning, but the
    mechanism generalises to any per-gene caveat the parser sets."""
    seen: set[str] = set()
    out: list[str] = []
    for bucket in (
        parsed.definitive_genes,
        parsed.ambiguous_genes,
        parsed.no_call_genes,
    ):
        for g in bucket:
            cav = getattr(g, "caveat", None)
            if cav and cav not in seen:
                seen.add(cav)
                out.append(cav)
    return out


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
    total_positions = total_found + total_missing
    coverage_percent = (
        int(100 * total_found / total_positions)
        if total_positions else 0
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
        "total_positions": total_positions,
        "coverage_percent": coverage_percent,
    }


# ---------------------------------------------------------------------------
# Gene enrichment
# ---------------------------------------------------------------------------

def _enrich_definitive(g, effective_risk: dict[str, str]) -> dict:
    risk = effective_risk.get(g.gene, g.risk_level)
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
        "risk_level": risk,
        "symbol": ACTION_SYMBOLS[risk],
        "label": ACTION_LABELS[risk],
        "category": g.category,
        "category_desc": g.category_desc,
        "positions_found": g.positions_found,
        "positions_missing": g.positions_missing,
        "related_drugs": g.related_drugs,
        # Per-gene caveat is suppressed in the template (hoisted to one
        # consolidated block at the end of the report). Keep on the dict
        # for downstream code that may still want it.
        "caveat": g.caveat,
        "drug_recs": [],
    }


def _enrich_ambiguous(g, effective_risk: dict[str, str]) -> dict:
    risk = effective_risk.get(g.gene, g.risk_level)
    return {
        "gene": g.gene,
        "diplotype_count": g.diplotype_count,
        "phenotype_range": g.phenotype_range,
        "harmful_phenotypes": g.harmful_phenotypes,
        "risk_level": risk,
        "symbol": ACTION_SYMBOLS[risk],
        "label": ACTION_LABELS[risk],
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
