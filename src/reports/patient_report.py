"""
Patient report — plain-language, colorblind-safe.

Combines:
- Adib's gene functional categories with category descriptions
- Adib's 4-level action symbols (▲ ◆ ✓ —)
- Marco's overview cards with traffic-light counts
- Marco's drug therapeutic categories
- Plain-language phenotype explanations + per-gene descriptions
"""

from collections import OrderedDict, defaultdict
from datetime import date
from pathlib import Path

from config.settings import (
    ACTION_LABELS, ACTION_SYMBOLS, APP_TITLE, CATEGORY_ORDER,
    GENE_CATEGORIES, PATIENT_DISCLAIMER, RISK_PRIORITY,
    drug_category, gene_tissue_icon,
)
from pharmcat.output_parser import (
    AmbiguousGene, DefinitiveGene, NoCallGene, ParsedResults,
)
from reports._render import html_and_pdf

# ---------------------------------------------------------------------------
# Plain-language explanations (4-level keyed)
# ---------------------------------------------------------------------------

PHENOTYPE_PLAIN: dict[str, dict] = {
    "normal metabolizer": {
        "brief": "Your body processes related medications at a typical rate — standard dosing applies.",
        "detail": "",
    },
    "intermediate metabolizer": {
        "brief": "You carry variants with somewhat reduced enzyme activity. Affected medications or their byproducts may build up to higher-than-usual levels, raising the risk of side effects — a reduced starting dose with extra monitoring is typically recommended.",
        "detail": "",
    },
    "poor metabolizer": {
        "brief": "You carry variants with significantly reduced enzyme activity — affected drugs or their byproducts can build up to dangerously high levels, raising the risk of severe side effects. Substantial dose reduction (or an alternative drug) is typically needed.",
        "detail": "",
    },
    "rapid metabolizer": {
        "brief": "You carry variants with increased enzyme activity — affected medications may break down faster than usual, potentially reducing their effect.",
        "detail": "",
    },
    "ultrarapid metabolizer": {
        "brief": "You carry variants with significantly increased enzyme activity — standard doses may not be effective, or for some drugs may produce dangerously high active-drug levels.",
        "detail": "",
    },
    "decreased function": {
        "brief": "Reduced gene activity — may affect how certain medications are processed or transported in your body.",
        "detail": "",
    },
    "increased function": {
        "brief": "Elevated gene activity — may affect how certain medications are processed or transported in your body.",
        "detail": "",
    },
    "normal function": {
        "brief": "Standard function — typical dosing applies for related medications.",
        "detail": "",
    },
    "indeterminate": {
        "brief": "Your result for this gene could not be clearly determined — discuss with your doctor.",
        "detail": "",
    },
    "no result": {
        "brief": "No result could be determined for this gene from the available data.",
        "detail": "",
    },
    # ── Phenotypes that need direct lookup (not substring) ──────────────
    # PharmCAT reports these as standalone strings; without an exact-match
    # entry they would fall through to "indeterminate".
    "normal": {  # G6PD reports just "Normal" — identical text to "normal function" so they merge
        "brief": "Standard function — typical dosing applies for related medications.",
        "detail": "",
    },
    "n/a": {  # IFNL3 and some others when PharmCAT doesn't emit a labelled phenotype
        "brief": "No specific result was reported for this gene from your data.",
        "detail": "",
    },
    "uncertain susceptibility": {  # CACNA1S, RYR1 — reference / negative state for MH
        "brief": "No known risk variants were found — your test result falls within the expected range.",
        "detail": "",
    },
    # CFTR — drug-specific phenotype strings emitted by PharmCAT
    "ivacaftor non-responsive in cf patients": {
        "brief": "Your CFTR variant is not the type ivacaftor (Kalydeco) is designed to treat — this drug is unlikely to be effective.",
        "detail": "",
    },
    "ivacaftor responsive in cf patients": {
        "brief": "Your CFTR variant is one of the types ivacaftor (Kalydeco) is designed to treat — this drug is likely to be effective.",
        "detail": "",
    },
    # VKORC1 — reported as a literal SNP genotype string
    "-1639 aa": {
        "brief": "You carry the variant that increases warfarin sensitivity — lower-than-usual doses are typically needed to control blood clotting safely.",
        "detail": "",
    },
    "-1639 ag": {
        "brief": "You carry one copy of the warfarin-sensitivity variant — your prescriber may adjust the starting dose accordingly.",
        "detail": "",
    },
    "-1639 ga": {  # same as AG depending on allele order
        "brief": "You carry one copy of the warfarin-sensitivity variant — your prescriber may adjust the starting dose accordingly.",
        "detail": "",
    },
    "-1639 gg": {
        "brief": "Standard warfarin sensitivity — typical dosing applies.",
        "detail": "",
    },
    # Less common metabolizer / function labels
    "extensive metabolizer": {
        "brief": "Your body processes related medications at a typical rate — standard dosing applies.",
        "detail": "",
    },
    "poor function": {
        "brief": "You carry variants associated with reduced gene activity — may affect how related medications are processed or how the drug target responds.",
        "detail": "",
    },
}

DRUG_GUIDANCE: dict[str, str] = {
    "action": "Important: Talk to your doctor before taking this medication. Your genetics suggest it may need special consideration.",
    "review": "Note: Discuss this medication with your doctor. A dose adjustment or extra monitoring may be helpful.",
    "normal": "Standard use is expected to be appropriate based on your genetics.",
    "nodata": "Insufficient data to provide guidance for this medication.",
}

PHENOTYPE_SHORT: dict[str, str] = {
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
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_patient_report(
    parsed: ParsedResults,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Render the patient report (HTML + best-effort PDF)."""
    context = _build_context(parsed)
    return html_and_pdf("patient_report.html", context, Path(output_dir), "patient_report")


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def _build_context(parsed: ParsedResults) -> dict:
    overview_cards = _build_overview_cards(parsed)
    counts = _count_by_risk(overview_cards)

    # Functional category groupings (with all three buckets)
    categories = _group_by_category(parsed)

    # Drug guidance grouped by therapeutic category, action-first
    drugs_by_category = _group_drugs_for_patient(parsed)

    summary_sentence = _summary_sentence(counts, parsed)

    return {
        "title": "Your Personal Pharmacogenetics Report",
        "app_title": APP_TITLE,
        "report_date": date.today().strftime("%B %d, %Y"),
        "sample_id": parsed.sample_id,
        "metadata": parsed.metadata,
        "disclaimer": PATIENT_DISCLAIMER,
        "action_symbols": ACTION_SYMBOLS,
        "action_labels": ACTION_LABELS,
        "summary_sentence": summary_sentence,
        "counts": counts,
        "overview_cards": overview_cards,
        "categories": categories,
        "drugs_by_category": drugs_by_category,
        "phenotype_plain": PHENOTYPE_PLAIN,
        "drug_guidance": DRUG_GUIDANCE,
        "category_order": CATEGORY_ORDER,
    }


# ---------------------------------------------------------------------------
# Overview cards (Marco)
# ---------------------------------------------------------------------------

def _build_overview_cards(parsed: ParsedResults) -> list[dict]:
    cards: list[dict] = []

    def _tissue(gene: str) -> dict:
        key, label, icon = gene_tissue_icon(gene)
        return {"tissue": key, "tissue_label": label, "tissue_icon": icon}

    for g in parsed.definitive_genes:
        cards.append({
            "gene": g.gene,
            "risk_level": g.risk_level,
            "symbol": ACTION_SYMBOLS[g.risk_level],
            "label": ACTION_LABELS[g.risk_level],
            "protein_type": g.protein_type,
            "phenotype_short": PHENOTYPE_SHORT.get(g.phenotype, g.phenotype),
            "med_count": len(g.related_drugs),
            "bucket": "definitive",
            **_tissue(g.gene),
        })

    for g in parsed.ambiguous_genes:
        cards.append({
            "gene": g.gene,
            "risk_level": g.risk_level,
            "symbol": ACTION_SYMBOLS[g.risk_level],
            "label": ACTION_LABELS[g.risk_level],
            "protein_type": g.protein_type,
            "phenotype_short": "Inconclusive",
            "med_count": len(g.related_drugs),
            "bucket": "ambiguous",
            **_tissue(g.gene),
        })

    for g in parsed.no_call_genes:
        cards.append({
            "gene": g.gene,
            "risk_level": "nodata",
            "symbol": ACTION_SYMBOLS["nodata"],
            "label": ACTION_LABELS["nodata"],
            "protein_type": g.protein_type,
            "phenotype_short": "Not tested",
            "med_count": len(g.related_drugs),
            "bucket": "no_call",
            **_tissue(g.gene),
        })

    cards.sort(key=lambda c: (RISK_PRIORITY.get(c["risk_level"], 3), c["gene"]))
    return cards


def _count_by_risk(cards: list[dict]) -> dict[str, int]:
    counts = {"action": 0, "review": 0, "normal": 0, "nodata": 0}
    for c in cards:
        counts[c["risk_level"]] = counts.get(c["risk_level"], 0) + 1
    counts["total"] = sum(counts.values())
    return counts


def _summary_sentence(counts: dict[str, int], parsed: ParsedResults) -> str:
    total = counts["total"]
    if total == 0:
        return "No genes were analyzed."

    parts: list[str] = []
    if counts["normal"]:
        n = counts["normal"]
        parts.append(f"{n} gene{'s' if n != 1 else ''} show{'s' if n == 1 else ''} normal results")
    if counts["action"] or counts["review"]:
        n = counts["action"] + counts["review"]
        parts.append(f"{n} gene{'s' if n != 1 else ''} may affect your medication plan")
    if counts["nodata"]:
        n = counts["nodata"]
        parts.append(f"{n} gene{'s' if n != 1 else ''} could not be fully determined")

    return ". ".join(parts) + "." if parts else f"Your test analyzed {total} genes."


# ---------------------------------------------------------------------------
# Functional category grouping
# ---------------------------------------------------------------------------

def _group_by_category(parsed: ParsedResults) -> "OrderedDict[str, dict]":
    """Group every gene (in any bucket) by functional category, then within
    each category split into definitive / ambiguous / no-call sub-lists."""
    cats: OrderedDict[str, dict] = OrderedDict()
    for cat_name in CATEGORY_ORDER:
        cats[cat_name] = {
            "name": cat_name,
            "description": GENE_CATEGORIES[cat_name]["description"],
            "definitive": [],
            "ambiguous": [],
            "no_call": [],
        }

    for g in parsed.definitive_genes:
        cats[g.category]["definitive"].append(_enrich_definitive(g))
    for g in parsed.ambiguous_genes:
        cats[g.category]["ambiguous"].append(_enrich_ambiguous(g))
    for g in parsed.no_call_genes:
        cats[g.category]["no_call"].append(_enrich_no_call(g))

    # Drop empty categories
    return OrderedDict(
        (k, v) for k, v in cats.items()
        if v["definitive"] or v["ambiguous"] or v["no_call"]
    )


def _enrich_definitive(g: DefinitiveGene) -> dict:
    explanations = _explain_phenotype(g.phenotype)
    return {
        "gene": g.gene,
        "diplotype": g.diplotype,
        "phenotype": g.phenotype,
        "phenotype_short": PHENOTYPE_SHORT.get(g.phenotype, g.phenotype),
        "activity_score": g.activity_score,
        "star_alleles": g.star_alleles,
        "risk_level": g.risk_level,
        "symbol": ACTION_SYMBOLS[g.risk_level],
        "label": ACTION_LABELS[g.risk_level],
        "description": g.description,
        "protein_type": g.protein_type,
        "plain_language": explanations["brief"],
        "detail": explanations["detail"],
        "caveat": g.caveat,
    }


def _enrich_ambiguous(g: AmbiguousGene) -> dict:
    return {
        "gene": g.gene,
        "diplotype_count": g.diplotype_count,
        "phenotype_range": g.phenotype_range,
        "harmful_phenotypes": g.harmful_phenotypes,
        "risk_level": g.risk_level,
        "symbol": ACTION_SYMBOLS[g.risk_level],
        "label": ACTION_LABELS[g.risk_level],
        "description": g.description,
        "protein_type": g.protein_type,
        "actionable_drugs": g.actionable_drugs,
        "caveat": g.caveat,
    }


def _enrich_no_call(g: NoCallGene) -> dict:
    return {
        "gene": g.gene,
        "risk_level": "nodata",
        "symbol": ACTION_SYMBOLS["nodata"],
        "label": ACTION_LABELS["nodata"],
        "description": g.description,
        "protein_type": g.protein_type,
        "affected_drugs": g.affected_drugs,
        "caveat": g.caveat,
    }


def _explain_phenotype(phenotype: str) -> dict[str, str]:
    if not phenotype:
        return PHENOTYPE_PLAIN["no result"]
    lower = phenotype.lower().strip()
    # Exact match first so standalone strings like "Normal", "n/a", or
    # "Uncertain Susceptibility" land on their specific entries instead of
    # falling through the substring loop to "indeterminate".
    if lower in PHENOTYPE_PLAIN:
        return PHENOTYPE_PLAIN[lower]
    # Substring fallback for compound phenotypes ("Normal Metabolizer",
    # "Likely Poor Metabolizer", etc.).
    for key, val in PHENOTYPE_PLAIN.items():
        if key in lower:
            return val
    return PHENOTYPE_PLAIN["indeterminate"]


# ---------------------------------------------------------------------------
# Drug grouping for the patient view
# ---------------------------------------------------------------------------

def _group_drugs_for_patient(parsed: ParsedResults) -> list[dict]:
    """Group drug recommendations by therapeutic category, action-first
    within each category. One row per (drug, gene-or-blank) — we deduplicate
    multi-source recommendations to the highest-risk one. Drugs whose final
    risk level is 'normal' (no need to worry) are dropped from the patient view."""
    by_key: dict[tuple[str, str], dict] = {}

    for d in parsed.drugs:
        gene_label = ", ".join(d.affected_genes) if d.affected_genes else ""
        key = (d.drug, gene_label)
        guidance = DRUG_GUIDANCE.get(d.risk_level, DRUG_GUIDANCE["review"])

        candidate = {
            "drug": d.drug,
            "gene": gene_label,
            "risk_level": d.risk_level,
            "symbol": ACTION_SYMBOLS[d.risk_level],
            "label": ACTION_LABELS[d.risk_level],
            "patient_guidance": guidance,
            "therapeutic_category": d.therapeutic_category,
        }
        existing = by_key.get(key)
        if (
            existing is None
            or RISK_PRIORITY.get(d.risk_level, 3)
            < RISK_PRIORITY.get(existing["risk_level"], 3)
        ):
            by_key[key] = candidate

    # Drop drugs whose final (highest-risk) recommendation is 'normal'.
    # We keep 'nodata' so the patient knows we couldn't analyze the gene.
    grouped: dict[str, list[dict]] = defaultdict(list)
    for d in by_key.values():
        if d["risk_level"] == "normal":
            continue
        grouped[d["therapeutic_category"]].append(d)

    out: list[dict] = []
    for cat in sorted(grouped.keys()):
        drugs = sorted(
            grouped[cat],
            key=lambda d: (RISK_PRIORITY.get(d["risk_level"], 3), d["drug"]),
        )
        out.append({"category": cat, "drugs": drugs})

    out.sort(
        key=lambda group: min(
            RISK_PRIORITY.get(d["risk_level"], 3) for d in group["drugs"]
        )
    )
    return out
