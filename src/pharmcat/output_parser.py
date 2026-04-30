"""
PharmCAT output parser — unified model.

Reads phenotype.json + report.json + match.json (Marco's JSON-based approach,
no HTML scraping) and produces the three-bucket gene model (definitive /
ambiguous / no-call) with the four-level colorblind-safe risk vocabulary
(action / review / normal / nodata).

Each gene is annotated with:
    - Functional category (Phase I / II / Transporters / HLA / Other)
    - Plain-language description and protein type (for the patient report)
    - CYP2D6 caveat where applicable

Each drug recommendation is annotated with:
    - 4-level risk_level
    - Therapeutic category (Antidepressants, Statins, etc.)
    - Citation list (collected and dedup'd at the parsed-results level)
"""

import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from config.settings import (
    CYP2D6_CAVEAT,
    GENE_DESCRIPTIONS,
    GENE_PROTEIN_TYPE,
    drug_category,
    gene_category,
    phenotype_to_risk,
)

logger = logging.getLogger(__name__)


# Phenotypes that clearly indicate clinical concern in the ambiguous bucket
_HARMFUL_PHENOTYPES = {
    "Poor Metabolizer",
    "Ultrarapid Metabolizer",
    "Increased Function",
    "Poor Function",
    "Likely Poor Metabolizer",
}

_NORMAL_PHENOTYPES = {
    "Normal Metabolizer",
    "Normal Function",
    "Extensive Metabolizer",
}

# Keywords used to refine drug-recommendation risk
_ACTION_KEYWORDS = (
    "avoid", "contraindicated", "do not use", "not recommended",
    "consider alternative",
)
_REVIEW_KEYWORDS = (
    "reduce", "increase", "adjust", "lower dose", "higher dose",
    "caution", "monitor", "consider", "decreased dose", "titrate",
    "may need", "dose adjustment",
)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class DefinitiveGene:
    gene: str
    diplotype: str
    phenotype: str
    activity_score: str | None
    star_alleles: list[str]
    risk_level: str
    allele1_name: str
    allele1_function: str
    allele2_name: str
    allele2_function: str
    positions_found: int
    positions_missing: int
    related_drugs: list[str]
    category: str
    category_desc: str
    description: str          # plain-language gene description
    protein_type: str
    caveat: str | None = None


@dataclass
class AmbiguousGene:
    gene: str
    diplotype_count: int
    phenotype_range: list[str]
    has_harmful_phenotypes: bool
    harmful_phenotypes: list[str]
    risk_level: str
    actionable_drugs: dict[str, list[str]]   # therapeutic category -> [drug names]
    related_drugs: list[str]
    positions_found: int
    positions_missing: int
    category: str
    category_desc: str
    description: str
    protein_type: str
    caveat: str | None = None


@dataclass
class NoCallGene:
    gene: str
    positions_missing: int
    affected_drugs: dict[str, list[str]]
    related_drugs: list[str]
    risk_level: str = "nodata"
    category: str = ""
    category_desc: str = ""
    description: str = ""
    protein_type: str = ""
    caveat: str | None = None


@dataclass
class Citation:
    pmid: str
    title: str
    journal: str
    year: int | None


@dataclass
class DrugRecommendation:
    drug: str
    source: str               # CPIC / DPWG / FDA / etc.
    recommendation: str
    classification: str
    implications: list[str]
    affected_genes: list[str]
    phenotypes: dict[str, str]
    population: str
    risk_level: str           # 4-level: action / review / normal / nodata
    therapeutic_category: str
    urls: list[str]
    citations: list[Citation]
    messages: list[str]


@dataclass
class ParsedResults:
    sample_id: str | None = None
    metadata: dict = field(default_factory=dict)
    definitive_genes: list[DefinitiveGene] = field(default_factory=list)
    ambiguous_genes: list[AmbiguousGene] = field(default_factory=list)
    no_call_genes: list[NoCallGene] = field(default_factory=list)
    drugs: list[DrugRecommendation] = field(default_factory=list)
    citations: list[Citation] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_pharmcat_output(
    phenotype_json_path: str | Path | None = None,
    report_json_path: str | Path | None = None,
    match_json_path: str | Path | None = None,
    sample_id: str | None = None,
) -> ParsedResults:
    """Parse PharmCAT output into the unified model.

    All paths are optional but `report_json_path` is needed for drug
    recommendations and `match_json_path` is needed to classify
    definitive vs. ambiguous vs. no-call.
    """
    results = ParsedResults(sample_id=sample_id)

    # Load JSONs
    report = _safe_load(report_json_path)
    match_data = _safe_load(match_json_path)

    if report:
        results.metadata = _parse_metadata(report)
        results.messages = [m.get("message", "") for m in report.get("messages", []) if m.get("message")]

    # Drug recs first so we can reference them when classifying genes
    if report:
        results.drugs = _parse_drugs(report)
        results.citations = _collect_citations(results.drugs)

    # Build a lookup of match results by gene
    match_by_gene: dict[str, dict] = {}
    if match_data:
        for entry in match_data.get("results", []):
            gene = entry.get("gene")
            if gene:
                match_by_gene[gene] = entry

    # Classify each gene
    if report:
        for gene_symbol, gene_data in report.get("genes", {}).items():
            _classify_gene(
                gene_symbol, gene_data, match_by_gene.get(gene_symbol),
                results,
            )

    # Sort each bucket
    _sort_buckets(results)

    # Attach CYP2D6 caveat where applicable
    for bucket in (results.definitive_genes, results.ambiguous_genes, results.no_call_genes):
        for g in bucket:
            if g.gene == "CYP2D6":
                g.caveat = CYP2D6_CAVEAT

    # Sample ID from phenotype.json metadata if not provided
    if not results.sample_id:
        pheno = _safe_load(phenotype_json_path)
        if pheno:
            meta = pheno.get("matcherMetadata", {})
            input_file = meta.get("inputFilename", "")
            if input_file:
                stem = Path(input_file).name
                for suffix in (".vcf.gz", ".vcf"):
                    if stem.lower().endswith(suffix):
                        stem = stem[: -len(suffix)]
                        break
                results.sample_id = stem

    return results


# ---------------------------------------------------------------------------
# Gene classification
# ---------------------------------------------------------------------------

def _classify_gene(
    gene_symbol: str,
    gene_data: dict,
    match_result: dict | None,
    results: ParsedResults,
) -> None:
    diplotypes_report = gene_data.get("sourceDiplotypes", [])
    related_drugs = [
        d.get("name", d) if isinstance(d, dict) else d
        for d in gene_data.get("relatedDrugs", []) or []
    ]

    if match_result:
        match_section = match_result.get("matchData", {})
        n_variants = len(match_result.get("variants", []) or [])
        n_diplotypes_match = len(match_result.get("diplotypes", []) or [])
        n_missing = len(match_section.get("missingPositions", []) or [])
        n_found = n_variants
    else:
        n_variants = n_diplotypes_match = n_missing = n_found = 0

    is_no_call = (
        n_diplotypes_match == 0
        and n_variants == 0
        and _is_unknown_diplotype(diplotypes_report)
    )
    is_ambiguous = n_diplotypes_match > 1

    cat_name, cat_desc = gene_category(gene_symbol)
    description = GENE_DESCRIPTIONS.get(gene_symbol, "")
    protein_type = GENE_PROTEIN_TYPE.get(gene_symbol, "Gene")

    if is_no_call:
        results.no_call_genes.append(_build_no_call(
            gene_symbol, n_missing, related_drugs,
            cat_name, cat_desc, description, protein_type,
        ))
    elif is_ambiguous:
        results.ambiguous_genes.append(_build_ambiguous(
            gene_symbol, gene_data,
            n_diplotypes_match, n_found, n_missing,
            related_drugs, results.drugs,
            cat_name, cat_desc, description, protein_type,
        ))
    else:
        results.definitive_genes.append(_build_definitive(
            gene_symbol, gene_data,
            n_found, n_missing, related_drugs,
            cat_name, cat_desc, description, protein_type,
        ))


def _is_unknown_diplotype(diplotypes_report: list[dict]) -> bool:
    if not diplotypes_report:
        return True
    return all(
        sd.get("label", "") in ("Unknown/Unknown", "Unknown", "")
        and sd.get("phenotypes", ["No Result"]) == ["No Result"]
        for sd in diplotypes_report
    )


def _build_definitive(
    gene_symbol, gene_data, n_found, n_missing, related_drugs,
    cat_name, cat_desc, description, protein_type,
) -> DefinitiveGene:
    diplotypes = gene_data.get("sourceDiplotypes", []) or []
    top = diplotypes[0] if diplotypes else {}
    label = top.get("label", "Unknown")
    phenotypes = top.get("phenotypes", []) or []
    phenotype = phenotypes[0] if phenotypes else "No Result"
    activity_score = top.get("activityScore")
    activity_score = str(activity_score) if activity_score is not None else None

    a1 = top.get("allele1") or {}
    a2 = top.get("allele2") or {}
    a1_name = a1.get("name", "Unknown")
    a2_name = a2.get("name", "Unknown")
    star_alleles = [a for a in (a1_name, a2_name) if a and a != "Unknown"]

    return DefinitiveGene(
        gene=gene_symbol,
        diplotype=label,
        phenotype=phenotype,
        activity_score=activity_score,
        star_alleles=star_alleles,
        risk_level=phenotype_to_risk(phenotype),
        allele1_name=a1_name,
        allele1_function=a1.get("function", "Unknown"),
        allele2_name=a2_name,
        allele2_function=a2.get("function", "Unknown"),
        positions_found=n_found,
        positions_missing=n_missing,
        related_drugs=related_drugs,
        category=cat_name,
        category_desc=cat_desc,
        description=description,
        protein_type=protein_type,
    )


def _build_ambiguous(
    gene_symbol, gene_data, diplotype_count, n_found, n_missing,
    related_drugs, all_drugs,
    cat_name, cat_desc, description, protein_type,
) -> AmbiguousGene:
    diplotypes = gene_data.get("sourceDiplotypes", []) or []
    phenotype_set: set[str] = set()
    for sd in diplotypes:
        for p in sd.get("phenotypes", []) or []:
            if p and p != "n/a":
                phenotype_set.add(p)

    harmful = sorted(phenotype_set & _HARMFUL_PHENOTYPES)
    has_harmful = bool(harmful)

    # Risk: action if any harmful phenotype is in the range, otherwise review
    # (ambiguity itself warrants clinician review).
    risk = "action" if has_harmful else "review"

    actionable = _find_actionable_drugs(gene_symbol, all_drugs)

    return AmbiguousGene(
        gene=gene_symbol,
        diplotype_count=diplotype_count,
        phenotype_range=sorted(phenotype_set),
        has_harmful_phenotypes=has_harmful,
        harmful_phenotypes=harmful,
        risk_level=risk,
        actionable_drugs=actionable,
        related_drugs=related_drugs,
        positions_found=n_found,
        positions_missing=n_missing,
        category=cat_name,
        category_desc=cat_desc,
        description=description,
        protein_type=protein_type,
    )


def _build_no_call(
    gene_symbol, n_missing, related_drugs,
    cat_name, cat_desc, description, protein_type,
) -> NoCallGene:
    affected: dict[str, list[str]] = defaultdict(list)
    for drug_name in related_drugs:
        affected[drug_category(drug_name)].append(drug_name)
    affected_sorted = {
        cat: sorted(set(drugs))
        for cat, drugs in sorted(affected.items())
    }

    return NoCallGene(
        gene=gene_symbol,
        positions_missing=n_missing,
        affected_drugs=affected_sorted,
        related_drugs=related_drugs,
        category=cat_name,
        category_desc=cat_desc,
        description=description,
        protein_type=protein_type,
    )


def _find_actionable_drugs(
    gene_symbol: str, drugs: list[DrugRecommendation],
) -> dict[str, list[str]]:
    """For an ambiguous gene, list drugs whose recommendation hinges on a
    non-normal phenotype of this gene. Grouped by therapeutic category."""
    categorized: dict[str, set[str]] = defaultdict(set)

    for d in drugs:
        if gene_symbol not in d.affected_genes:
            continue
        gene_phenotype = d.phenotypes.get(gene_symbol, "")
        if gene_phenotype in _NORMAL_PHENOTYPES:
            continue
        rec_lower = d.recommendation.lower()
        if any(kw in rec_lower for kw in _ACTION_KEYWORDS + _REVIEW_KEYWORDS):
            categorized[d.therapeutic_category].add(d.drug)

    return {cat: sorted(s) for cat, s in sorted(categorized.items())}


# ---------------------------------------------------------------------------
# Drug parser
# ---------------------------------------------------------------------------

def _parse_drugs(report: dict) -> list[DrugRecommendation]:
    out: list[DrugRecommendation] = []

    for source_name, drugs_by_source in (report.get("drugs", {}) or {}).items():
        for drug_name, drug_data in (drugs_by_source or {}).items():
            urls = drug_data.get("urls", []) or []
            citations_raw = drug_data.get("citations", []) or []
            messages = [
                m.get("message", "")
                for m in drug_data.get("messages", []) or []
                if m.get("message")
            ]

            citations = [
                Citation(
                    pmid=c.get("pmid", "") or "",
                    title=c.get("title", "") or "",
                    journal=c.get("journal", "") or "",
                    year=c.get("year"),
                )
                for c in citations_raw
            ]

            for guideline in drug_data.get("guidelines", []) or []:
                for ann in guideline.get("annotations", []) or []:
                    recommendation = ann.get("drugRecommendation", "") or ""
                    if not recommendation:
                        continue
                    classification = ann.get("classification", "") or ""
                    implications_raw = ann.get("implications", []) or []
                    if isinstance(implications_raw, str):
                        implications = [implications_raw]
                    else:
                        implications = list(implications_raw)
                    phenotype_map = ann.get("phenotypes", {}) or {}
                    population = ann.get("population", "") or ""
                    affected_genes = list(phenotype_map.keys())
                    risk_level = _drug_risk_level(recommendation, classification)

                    out.append(DrugRecommendation(
                        drug=drug_name,
                        source=source_name,
                        recommendation=recommendation,
                        classification=classification,
                        implications=implications,
                        affected_genes=affected_genes,
                        phenotypes=phenotype_map,
                        population=population,
                        risk_level=risk_level,
                        therapeutic_category=drug_category(drug_name),
                        urls=urls,
                        citations=citations,
                        messages=messages,
                    ))

    return out


def _drug_risk_level(recommendation: str, classification: str) -> str:
    """Map a drug recommendation to the 4-level risk vocabulary."""
    cls = classification.lower().strip()
    rec = (recommendation or "").lower()

    if cls in ("no action", "no action needed"):
        return "normal"
    if cls in ("strong", "actionable"):
        return "action"
    if cls in ("moderate", "informative"):
        return "review"

    if any(kw in rec for kw in _ACTION_KEYWORDS):
        return "action"
    if any(kw in rec for kw in _REVIEW_KEYWORDS):
        return "review"
    return "normal"


# ---------------------------------------------------------------------------
# Citation collection
# ---------------------------------------------------------------------------

def _collect_citations(drugs: list[DrugRecommendation]) -> list[Citation]:
    seen: set[str] = set()
    out: list[Citation] = []
    for d in drugs:
        for c in d.citations:
            key = c.pmid or c.title
            if not key or key in seen:
                continue
            seen.add(key)
            out.append(c)
    out.sort(key=lambda c: (-(c.year or 0), c.title))
    return out


# ---------------------------------------------------------------------------
# Sorting
# ---------------------------------------------------------------------------

def _sort_buckets(results: ParsedResults) -> None:
    risk_order = {"action": 0, "review": 1, "normal": 2, "nodata": 3}

    results.definitive_genes.sort(
        key=lambda g: (risk_order.get(g.risk_level, 3), g.gene)
    )
    results.ambiguous_genes.sort(
        key=lambda g: (-len(g.actionable_drugs), g.gene)
    )
    results.no_call_genes.sort(key=lambda g: g.gene)
    results.drugs.sort(
        key=lambda d: (risk_order.get(d.risk_level, 3), d.drug)
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_metadata(raw: dict) -> dict:
    return {
        "pharmcat_version": raw.get("pharmcatVersion", "Unknown"),
        "data_version": raw.get("dataVersion", "Unknown"),
        "timestamp": raw.get("timestamp", "Unknown"),
    }


def _safe_load(path: str | Path | None) -> dict | None:
    if path is None:
        return None
    p = Path(path)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        logger.warning("Failed to parse %s: %s", p, e)
        return None
