"""
Parse PharmCAT report JSON and match JSON into clean Python data structures
for the patient and clinician report views.

Handles three categories of gene results:
1. Definitive  - single diplotype identified (match.json has 1 diplotype candidate)
2. Ambiguous   - multiple diplotype candidates (match.json has >1 diplotype candidates)
3. No Call     - no variants found in VCF (match.json has 0 variants, 0 diplotype candidates)
"""

import json
import os
from collections import defaultdict

from config.settings import PHENOTYPE_COLORS, TARGET_GENES

# ---------------------------------------------------------------------------
# Drug therapeutic category mapping
# ---------------------------------------------------------------------------

DRUG_CATEGORIES = {
    # Antidepressants
    "amitriptyline": "Antidepressants",
    "citalopram": "Antidepressants",
    "clomipramine": "Antidepressants",
    "doxepin": "Antidepressants",
    "escitalopram": "Antidepressants",
    "imipramine": "Antidepressants",
    "sertraline": "Antidepressants",
    "trimipramine": "Antidepressants",
    "nortriptyline": "Antidepressants",
    "paroxetine": "Antidepressants",
    "fluvoxamine": "Antidepressants",
    "venlafaxine": "Antidepressants",
    "fluoxetine": "Antidepressants",
    "desipramine": "Antidepressants",
    "protriptyline": "Antidepressants",
    "amoxapine": "Antidepressants",
    "vortioxetine": "Antidepressants",
    "flibanserin": "Antidepressants",
    # Antiplatelet / Anticoagulants
    "clopidogrel": "Antiplatelet / Anticoagulant agents",
    "warfarin": "Antiplatelet / Anticoagulant agents",
    "acenocoumarol": "Antiplatelet / Anticoagulant agents",
    # Statins
    "simvastatin": "Statins (cholesterol)",
    "atorvastatin": "Statins (cholesterol)",
    "lovastatin": "Statins (cholesterol)",
    "rosuvastatin": "Statins (cholesterol)",
    "fluvastatin": "Statins (cholesterol)",
    "pitavastatin": "Statins (cholesterol)",
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
    "mavacamten": "Cardiovascular",
    "metoprolol": "Cardiovascular",
    "propranolol": "Cardiovascular",
    "carvedilol": "Cardiovascular",
    "nebivolol": "Cardiovascular",
    "flecainide": "Cardiovascular",
    "propafenone": "Cardiovascular",
    # Immunosuppressants
    "tacrolimus": "Immunosuppressants",
    # Chemotherapy / Immunosuppressants
    "fluorouracil": "Chemotherapy",
    "capecitabine": "Chemotherapy",
    "mercaptopurine": "Chemotherapy / Immunosuppressants",
    "azathioprine": "Chemotherapy / Immunosuppressants",
    "thioguanine": "Chemotherapy / Immunosuppressants",
    "tamoxifen": "Chemotherapy",
    # Pain medications
    "codeine": "Pain medications",
    "tramadol": "Pain medications",
    "hydrocodone": "Pain medications",
    # Antivirals (HIV)
    "abacavir": "Antivirals (HIV)",
    "atazanavir": "Antivirals (HIV)",
    # Antiepileptics
    "carbamazepine": "Antiepileptics",
    "phenytoin": "Antiepileptics",
    "oxcarbazepine": "Antiepileptics",
    "brivaracetam": "Antiepileptics",
    "clobazam": "Antiepileptics",
    "diazepam": "Antiepileptics",
    # Antipsychotics
    "aripiprazole": "Antipsychotics",
    "haloperidol": "Antipsychotics",
    "pimozide": "Antipsychotics",
    "clozapine": "Antipsychotics",
    "risperidone": "Antipsychotics",
    "brexpiprazole": "Antipsychotics",
    "iloperidone": "Antipsychotics",
    "perphenazine": "Antipsychotics",
    "thioridazine": "Antipsychotics",
    "zuclopenthixol": "Antipsychotics",
    # Gout
    "allopurinol": "Gout medications",
    # ADHD
    "atomoxetine": "ADHD medications",
    "amphetamine": "ADHD medications",
    "viloxazine": "ADHD medications",
    # Anesthetics
    "desflurane": "Anesthetics",
    "enflurane": "Anesthetics",
    "halothane": "Anesthetics",
    "isoflurane": "Anesthetics",
    "sevoflurane": "Anesthetics",
    "succinylcholine": "Anesthetics",
    # Anti-nausea
    "ondansetron": "Anti-nausea medications",
    "tropisetron": "Anti-nausea medications",
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
}

# Phenotypes that indicate potential clinical concern
_HARMFUL_PHENOTYPES = {
    "Poor Metabolizer",
    "Ultrarapid Metabolizer",
    "Increased Function",
    "Poor Function",
    "Likely Poor Metabolizer",
}

# Normal/benign phenotypes (no action needed)
_NORMAL_PHENOTYPES = {
    "Normal Metabolizer",
    "Normal Function",
    "Extensive Metabolizer",
}

# Keywords in recommendation text that indicate actionable guidance
_ACTIONABLE_KEYWORDS = [
    "avoid", "alternative", "reduce", "increase", "adjust",
    "lower", "caution", "consider", "contraindicated",
    "not recommended", "do not use",
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_full_report(report_json_path: str, match_json_path: str = None) -> dict:
    """
    Parse PharmCAT report.json and match.json into structured results
    with three gene categories: definitive, ambiguous, and no-call.

    Args:
        report_json_path: Path to the PharmCAT .report.json file.
        match_json_path:  Path to the PharmCAT .match.json file.
                          If None, derived from report_json_path by replacing
                          '.report.json' with '.match.json'.

    Returns:
        dict with keys: metadata, definitive_genes, ambiguous_genes,
        no_call_genes, drugs, messages.
    """
    # Load report JSON
    with open(report_json_path, encoding="utf-8") as f:
        report = json.load(f)

    # Resolve and load match JSON
    if match_json_path is None:
        match_json_path = report_json_path.replace(".report.json", ".match.json")

    match_data = None
    if os.path.isfile(match_json_path):
        with open(match_json_path, encoding="utf-8") as f:
            match_data = json.load(f)

    # Build a lookup of match results keyed by gene symbol
    match_by_gene = {}
    if match_data:
        for result in match_data.get("results", []):
            gene = result.get("gene")
            if gene:
                match_by_gene[gene] = result

    # Parse drugs first so we can reference them from gene results
    drugs = _parse_drugs(report)

    # Classify genes
    definitive_genes = []
    ambiguous_genes = []
    no_call_genes = []

    genes_raw = report.get("genes", {})
    for gene_symbol, gene_data in genes_raw.items():
        match_result = match_by_gene.get(gene_symbol)
        _classify_gene(
            gene_symbol, gene_data, match_result, drugs,
            definitive_genes, ambiguous_genes, no_call_genes,
        )

    # Sort definitive: actionable (red, orange) first, then green
    _color_priority = {"#F44336": 0, "#FF9800": 1, "#4CAF50": 2, "#9E9E9E": 3}
    definitive_genes.sort(
        key=lambda g: (_color_priority.get(g["color"], 3), g["gene"])
    )

    # Sort ambiguous by number of actionable drugs (most first)
    ambiguous_genes.sort(
        key=lambda g: (-len(g.get("actionable_drugs", {})), g["gene"])
    )

    # Sort no-call alphabetically
    no_call_genes.sort(key=lambda g: g["gene"])

    return {
        "metadata": _parse_metadata(report),
        "definitive_genes": definitive_genes,
        "ambiguous_genes": ambiguous_genes,
        "no_call_genes": no_call_genes,
        "drugs": drugs,
        "messages": [m.get("message", "") for m in report.get("messages", [])],
    }


def parse_report(report_json_path: str) -> dict:
    """
    Backward-compatible parser that returns the old flat gene list format.

    Returns:
        dict with keys: metadata, genes, drugs, messages.
    """
    full = parse_full_report(report_json_path)

    # Reshape the three gene categories back into a flat list
    genes = []
    for g in full["definitive_genes"]:
        genes.append(_to_legacy_gene(g, has_result=True))
    for g in full["ambiguous_genes"]:
        genes.append(_to_legacy_gene_ambiguous(g))
    for g in full["no_call_genes"]:
        genes.append(_to_legacy_gene_no_call(g))

    # Sort: target genes first, then alphabetically
    genes.sort(key=lambda g: (not g["is_target"], g["gene"]))

    return {
        "metadata": full["metadata"],
        "genes": genes,
        "drugs": full["drugs"],
        "messages": full["messages"],
    }


def get_summary_stats(parsed: dict) -> dict:
    """
    Generate summary statistics from parsed results.
    Works with both old format (genes key) and new format
    (definitive_genes / ambiguous_genes / no_call_genes keys).

    Returns:
        dict with: total_genes, genes_with_results, genes_needing_testing,
                   total_drug_recommendations, actionable_recommendations.
    """
    # Support both old and new formats
    if "genes" in parsed:
        genes = parsed["genes"]
        genes_with_results = sum(1 for g in genes if g.get("has_result", False))
        genes_needing_testing = sum(1 for g in genes if g.get("needs_testing", False))
        total_genes = len(genes)
    else:
        n_definitive = len(parsed.get("definitive_genes", []))
        n_ambiguous = len(parsed.get("ambiguous_genes", []))
        n_no_call = len(parsed.get("no_call_genes", []))
        total_genes = n_definitive + n_ambiguous + n_no_call
        genes_with_results = n_definitive
        genes_needing_testing = n_no_call

    drugs = parsed.get("drugs", [])

    return {
        "total_genes": total_genes,
        "genes_with_results": genes_with_results,
        "genes_needing_testing": genes_needing_testing,
        "total_drug_recommendations": len(drugs),
        "actionable_recommendations": sum(
            1 for d in drugs if d["urgency"] in ("red", "yellow")
        ),
    }


# ---------------------------------------------------------------------------
# Gene classification
# ---------------------------------------------------------------------------

def _classify_gene(
    gene_symbol, gene_data, match_result, drugs,
    definitive_out, ambiguous_out, no_call_out,
):
    """
    Classify a single gene into definitive, ambiguous, or no-call
    and append it to the appropriate output list.
    """
    diplotypes_report = gene_data.get("sourceDiplotypes", [])
    related_drugs_raw = gene_data.get("relatedDrugs", [])
    related_drugs = [d.get("name", d) if isinstance(d, dict) else d
                     for d in related_drugs_raw]

    # Coverage stats from match.json
    if match_result:
        match_data_section = match_result.get("matchData", {})
        n_variants = len(match_result.get("variants", []))
        n_diplotypes_match = len(match_result.get("diplotypes", []))
        n_missing = len(match_data_section.get("missingPositions", []))
        n_found = n_variants
    else:
        # Gene not present in match.json (e.g. CYP2D6 when callCyp2d=false)
        n_variants = 0
        n_diplotypes_match = 0
        n_missing = 0
        n_found = 0

    # Determine category based on match.json diplotype candidate count
    # No Call: 0 variants AND 0 diplotype candidates in match.json,
    #   OR gene missing from match.json with Unknown/Unknown in report
    is_no_call = (
        n_diplotypes_match == 0
        and n_variants == 0
        and _is_unknown_diplotype(diplotypes_report)
    )

    is_ambiguous = n_diplotypes_match > 1

    if is_no_call:
        _build_no_call(gene_symbol, n_missing, related_drugs, no_call_out)
    elif is_ambiguous:
        _build_ambiguous(
            gene_symbol, gene_data, match_result,
            n_diplotypes_match, n_found, n_missing,
            related_drugs, drugs, ambiguous_out,
        )
    else:
        _build_definitive(
            gene_symbol, gene_data,
            n_found, n_missing, related_drugs, definitive_out,
        )


def _is_unknown_diplotype(diplotypes_report):
    """Check if all report diplotypes are Unknown/Unknown with No Result."""
    if not diplotypes_report:
        return True
    return all(
        sd.get("label", "") in ("Unknown/Unknown", "Unknown", "")
        and sd.get("phenotypes", ["No Result"]) == ["No Result"]
        for sd in diplotypes_report
    )


# ---------------------------------------------------------------------------
# Definitive gene builder
# ---------------------------------------------------------------------------

def _build_definitive(gene_symbol, gene_data, n_found, n_missing,
                      related_drugs, out):
    """Build a definitive gene result dict."""
    diplotypes = gene_data.get("sourceDiplotypes", [])
    top = diplotypes[0] if diplotypes else {}

    label = top.get("label", "Unknown")
    phenotypes = top.get("phenotypes", [])
    phenotype = phenotypes[0] if phenotypes else "No Result"
    activity_score = top.get("activityScore")

    allele1 = top.get("allele1") or {}
    allele2 = top.get("allele2") or {}

    color = PHENOTYPE_COLORS.get(
        phenotype, PHENOTYPE_COLORS.get("Indeterminate", "#9E9E9E")
    )

    out.append({
        "gene": gene_symbol,
        "diplotype": label,
        "phenotype": phenotype,
        "activity_score": activity_score,
        "color": color,
        "positions_found": n_found,
        "positions_missing": n_missing,
        "allele1_name": allele1.get("name", "Unknown"),
        "allele1_function": allele1.get("function", "Unknown"),
        "allele2_name": allele2.get("name", "Unknown"),
        "allele2_function": allele2.get("function", "Unknown"),
        "related_drugs": related_drugs,
    })


# ---------------------------------------------------------------------------
# Ambiguous gene builder
# ---------------------------------------------------------------------------

def _build_ambiguous(gene_symbol, gene_data, match_result,
                     diplotype_count, n_found, n_missing,
                     related_drugs, all_drugs, out):
    """Build an ambiguous gene result dict."""
    diplotypes = gene_data.get("sourceDiplotypes", [])

    # Collect all unique phenotypes across candidate diplotypes
    phenotype_set = set()
    for sd in diplotypes:
        for p in sd.get("phenotypes", []):
            if p and p != "n/a":
                phenotype_set.add(p)
    phenotype_range = sorted(phenotype_set)

    # Identify harmful phenotypes in the range
    harmful = sorted(phenotype_set & _HARMFUL_PHENOTYPES)
    has_harmful = len(harmful) > 0

    # Find actionable drugs for this ambiguous gene
    actionable_drugs = _find_actionable_drugs(gene_symbol, all_drugs)

    out.append({
        "gene": gene_symbol,
        "diplotype_count": diplotype_count,
        "positions_found": n_found,
        "positions_missing": n_missing,
        "phenotype_range": phenotype_range,
        "has_harmful_phenotypes": has_harmful,
        "harmful_phenotypes": harmful,
        "actionable_drugs": actionable_drugs,
        "related_drugs": related_drugs,
    })


def _find_actionable_drugs(gene_symbol, all_drugs):
    """
    Scan all drug recommendations to find actionable ones for a gene.

    A drug is actionable if:
    - The recommendation involves the gene
    - The phenotype is NOT a normal/benign phenotype
    - The recommendation text contains action keywords
    """
    categorized = defaultdict(set)

    for drug_rec in all_drugs:
        if gene_symbol not in drug_rec.get("affected_genes", []):
            continue

        phenotypes = drug_rec.get("phenotypes", {})
        gene_phenotype = phenotypes.get(gene_symbol, "")

        # Skip normal phenotypes
        if gene_phenotype in _NORMAL_PHENOTYPES:
            continue

        rec_text = drug_rec.get("recommendation", "").lower()
        if any(kw in rec_text for kw in _ACTIONABLE_KEYWORDS):
            drug_name = drug_rec["drug"]
            category = DRUG_CATEGORIES.get(drug_name, "Other medications")
            categorized[category].add(drug_name)

    # Convert sets to sorted lists
    return {cat: sorted(drugs) for cat, drugs in sorted(categorized.items())}


# ---------------------------------------------------------------------------
# No-call gene builder
# ---------------------------------------------------------------------------

def _build_no_call(gene_symbol, n_missing, related_drugs, out):
    """Build a no-call gene result dict."""
    # Group related drugs by category
    affected_drugs = defaultdict(list)
    for drug_name in related_drugs:
        # Normalize: related_drugs may contain compound names; use base name
        base = drug_name.strip().lower()
        category = DRUG_CATEGORIES.get(base, "Other medications")
        affected_drugs[category].append(drug_name)

    # Sort drugs within each category
    affected_sorted = {
        cat: sorted(drugs)
        for cat, drugs in sorted(affected_drugs.items())
    }

    out.append({
        "gene": gene_symbol,
        "positions_missing": n_missing,
        "affected_drugs": affected_sorted,
        "related_drugs": related_drugs,
    })


# ---------------------------------------------------------------------------
# Metadata parser
# ---------------------------------------------------------------------------

def _parse_metadata(raw: dict) -> dict:
    return {
        "pharmcat_version": raw.get("pharmcatVersion", "Unknown"),
        "data_version": raw.get("dataVersion", "Unknown"),
        "timestamp": raw.get("timestamp", "Unknown"),
    }


# ---------------------------------------------------------------------------
# Drug recommendation parser
# ---------------------------------------------------------------------------

def _parse_drugs(raw: dict) -> list:
    """Extract per-drug recommendations from all guideline sources."""
    drugs_section = raw.get("drugs", {})
    results = []

    for source_name, drugs_by_source in drugs_section.items():
        for drug_name, drug_data in drugs_by_source.items():
            guidelines = drug_data.get("guidelines", [])
            drug_messages = [
                m.get("message", "") for m in drug_data.get("messages", [])
            ]
            urls = drug_data.get("urls", [])
            citations = drug_data.get("citations", [])

            for guideline in guidelines:
                annotations = guideline.get("annotations", [])

                for ann in annotations:
                    recommendation = ann.get("drugRecommendation", "")
                    if not recommendation:
                        continue

                    classification = ann.get("classification", "")
                    implications = ann.get("implications", [])
                    phenotype_map = ann.get("phenotypes", {})
                    population = ann.get("population", "")

                    affected_genes = list(phenotype_map.keys())
                    urgency = _classify_urgency(recommendation, classification)

                    results.append({
                        "drug": drug_name,
                        "source": source_name,
                        "recommendation": recommendation,
                        "classification": classification,
                        "implications": implications,
                        "affected_genes": affected_genes,
                        "phenotypes": phenotype_map,
                        "population": population,
                        "urgency": urgency,
                        "urls": urls,
                        "citations": [
                            {
                                "pmid": c.get("pmid", ""),
                                "title": c.get("title", ""),
                                "journal": c.get("journal", ""),
                                "year": c.get("year"),
                            }
                            for c in citations
                        ],
                        "messages": drug_messages,
                    })

    # Sort: actionable recommendations first, then alphabetically by drug
    results.sort(
        key=lambda d: (
            d["urgency"] != "red",
            d["urgency"] != "yellow",
            d["drug"],
        )
    )
    return results


# ---------------------------------------------------------------------------
# Urgency classifier
# ---------------------------------------------------------------------------

def _classify_urgency(recommendation: str, classification: str) -> str:
    """
    Classify a drug recommendation into traffic-light urgency.

    Returns: "green", "yellow", or "red"
    """
    rec_lower = recommendation.lower()

    # Red: avoid, contraindicated, do not use
    red_keywords = [
        "avoid", "contraindicated", "do not use",
        "not recommended", "consider alternative",
    ]
    if any(kw in rec_lower for kw in red_keywords):
        return "red"

    # Yellow: dose adjustment, caution, reduced dose, increased dose
    yellow_keywords = [
        "reduce", "increase", "adjust", "lower dose", "higher dose",
        "caution", "monitor", "consider", "decreased dose",
    ]
    if any(kw in rec_lower for kw in yellow_keywords):
        return "yellow"

    # Green: standard dose, use as directed
    return "green"


# ---------------------------------------------------------------------------
# Legacy format helpers (for backward compatibility)
# ---------------------------------------------------------------------------

def _to_legacy_gene(definitive_gene: dict, has_result: bool = True) -> dict:
    """Convert a definitive gene dict to the legacy flat format."""
    gene_symbol = definitive_gene["gene"]
    is_target = gene_symbol in TARGET_GENES
    phenotype = definitive_gene.get("phenotype", "No Result")
    return {
        "gene": gene_symbol,
        "diplotype": definitive_gene.get("diplotype", "Unknown"),
        "phenotype": phenotype,
        "has_result": has_result and phenotype != "No Result",
        "is_target": is_target,
        "needs_testing": False,
        "color": definitive_gene.get("color", "#9E9E9E"),
        "activity_score": definitive_gene.get("activity_score"),
        "allele1_name": definitive_gene.get("allele1_name", "Unknown"),
        "allele1_function": definitive_gene.get("allele1_function", "Unknown"),
        "allele2_name": definitive_gene.get("allele2_name", "Unknown"),
        "allele2_function": definitive_gene.get("allele2_function", "Unknown"),
        "related_drugs": definitive_gene.get("related_drugs", []),
    }


def _to_legacy_gene_ambiguous(ambiguous_gene: dict) -> dict:
    """Convert an ambiguous gene dict to the legacy flat format."""
    gene_symbol = ambiguous_gene["gene"]
    is_target = gene_symbol in TARGET_GENES
    # Use first phenotype in range or "Indeterminate"
    phenotype_range = ambiguous_gene.get("phenotype_range", [])
    phenotype = phenotype_range[0] if phenotype_range else "Indeterminate"
    color = PHENOTYPE_COLORS.get(
        phenotype, PHENOTYPE_COLORS.get("Indeterminate", "#9E9E9E")
    )
    return {
        "gene": gene_symbol,
        "diplotype": f"{ambiguous_gene.get('diplotype_count', '?')} possible diplotypes",
        "phenotype": phenotype,
        "has_result": False,
        "is_target": is_target,
        "needs_testing": is_target,
        "color": color,
        "activity_score": None,
        "allele1_name": "Unknown",
        "allele1_function": "Unknown",
        "allele2_name": "Unknown",
        "allele2_function": "Unknown",
        "related_drugs": ambiguous_gene.get("related_drugs", []),
    }


def _to_legacy_gene_no_call(no_call_gene: dict) -> dict:
    """Convert a no-call gene dict to the legacy flat format."""
    gene_symbol = no_call_gene["gene"]
    is_target = gene_symbol in TARGET_GENES
    return {
        "gene": gene_symbol,
        "diplotype": "Unknown/Unknown",
        "phenotype": "No Result",
        "has_result": False,
        "is_target": is_target,
        "needs_testing": is_target,
        "color": PHENOTYPE_COLORS.get("Indeterminate", "#9E9E9E"),
        "activity_score": None,
        "allele1_name": "Unknown",
        "allele1_function": "Unknown",
        "allele2_name": "Unknown",
        "allele2_function": "Unknown",
        "related_drugs": no_call_gene.get("related_drugs", []),
    }
