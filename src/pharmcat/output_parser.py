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
# Phrases that explicitly indicate "no change needed". CPIC frequently pairs
# a 'Strong' evidence classification with a "use standard dose" recommendation
# for normal metabolizers — those are NOT action items for the patient.
_NORMAL_PHRASES = (
    "standard dose",
    "standard dosing",
    "recommended starting dose",
    "label-recommended",
    "label recommended",
    "no indication to change",
    "no change in dose",
    "usual dose",
    "use as labeled",
    "no genotype-related",
)
# Negated-action phrases: text that contains an action keyword ("avoid",
# "contraindicated") but is actually telling the prescriber NOT to avoid /
# NOT to contraindicate. E.g. "No reason to avoid based on G6PD status".
# Checked BEFORE the action-keyword scan so these don't get flagged Action.
_FALSE_ACTION_PHRASES = (
    "no need to avoid",
    "no reason to avoid",
    "not contraindicated",
    "is not contraindicated",
    "may be used",
    "can be used",
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


_BENIGN_FUNCTION_KEYWORDS = ("normal function", "favorable")
_REFERENCE_ALLELE_NAMES = {"*1", "reference"}


def _is_benign_allele(name: str | None, function: str | None) -> bool:
    """True if a single allele looks unambiguously normal — either the
    canonical reference star allele (`*1`) or a function string PharmCAT
    explicitly tags as benign ('Normal function', 'Favorable response
    allele', etc.)."""
    if name and any(ref in name.lower() for ref in _REFERENCE_ALLELE_NAMES):
        return True
    if function:
        lowered = function.lower()
        return any(kw in lowered for kw in _BENIGN_FUNCTION_KEYWORDS)
    return False


def _benign_call(
    a1_name: str | None, a1_function: str | None,
    a2_name: str | None, a2_function: str | None,
) -> bool:
    """True when both alleles of a called diplotype look benign. Used to
    reclassify phenotype-less but definitively-called genes (CYP4F2 *1/*1,
    IFNL3 rs12979860 C/C) from 'nodata' to 'normal'."""
    return (
        _is_benign_allele(a1_name, a1_function)
        and _is_benign_allele(a2_name, a2_function)
    )


def _is_unknown_diplotype(diplotypes_report: list[dict]) -> bool:
    """True when PharmCAT couldn't determine a diplotype. The label is
    `Unknown/Unknown` (or `Unknown` for haploid genes like MT-RNR1), and
    the phenotypes field is either `["No Result"]` (most genes) or `[]`
    (HLA-A/HLA-B emit an empty list)."""
    if not diplotypes_report:
        return True
    return all(
        sd.get("label", "") in ("Unknown/Unknown", "Unknown", "")
        and (sd.get("phenotypes") or ["No Result"]) == ["No Result"]
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
    a1_function = a1.get("function", "Unknown")
    a2_function = a2.get("function", "Unknown")
    star_alleles = [a for a in (a1_name, a2_name) if a and a != "Unknown"]

    # Some genes (CYP4F2, IFNL3) have a definitive diplotype call but CPIC
    # doesn't define a diplotype-level phenotype — CYP4F2 only contributes
    # to the warfarin dosing algorithm, IFNL3 is interpreted at the allele
    # level. If both alleles look benign (reference *1, or function strings
    # tagged 'Normal'/'Favorable'), classify as normal rather than nodata
    # so the patient isn't mis-told their gene is unknown.
    risk_level = phenotype_to_risk(phenotype)
    if risk_level == "nodata" and _benign_call(a1_name, a1_function, a2_name, a2_function):
        risk_level = "normal"

    return DefinitiveGene(
        gene=gene_symbol,
        diplotype=label,
        phenotype=phenotype,
        activity_score=activity_score,
        star_alleles=star_alleles,
        risk_level=risk_level,
        allele1_name=a1_name,
        allele1_function=a1_function,
        allele2_name=a2_name,
        allele2_function=a2_function,
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
                    risk_level = _drug_risk_level(
                        recommendation, classification,
                        list(phenotype_map.values()),
                    )

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


# Phenotype phrases that may appear in a recommendation as a "gate" — they
# restrict an action keyword to a specific population (e.g. "consider
# alternative in poor metabolizers"). If the patient doesn't have that
# phenotype, the action keyword shouldn't apply to them.
_PHENOTYPE_GATE_KEYWORDS = (
    "normal metabolizer", "normal metabolizers",
    "intermediate metabolizer", "intermediate metabolizers",
    "poor metabolizer", "poor metabolizers",
    "rapid metabolizer", "rapid metabolizers",
    "ultrarapid metabolizer", "ultrarapid metabolizers",
    "extensive metabolizer", "extensive metabolizers",
    "normal function", "increased function", "decreased function",
    "poor function",
)


def _action_kw_gated_to_other_phenotype(
    rec_lower: str, kw_pos: int, patient_phenotypes: list[str],
) -> bool:
    """True if the action keyword starting at `kw_pos` is followed (within
    the same sentence) by an "in/for/with [phenotype]" clause whose phenotype
    is NOT in the patient's phenotype list. This catches cases like FDA-label
    text saying "Consider alternative therapy in poor metabolizers" when the
    patient is an intermediate metabolizer — the action recommendation is
    for a different population.

    If `patient_phenotypes` is empty, the heuristic can't decide and returns
    False (no suppression) — safer to keep the action trigger than to silently
    hide it."""
    if not patient_phenotypes:
        return False

    # Bound the search to the current sentence — avoid bleed-over to the
    # next sentence's gating phrase.
    sentence_end = rec_lower.find(". ", kw_pos)
    if sentence_end == -1:
        sentence_end = len(rec_lower)
    window = rec_lower[kw_pos:sentence_end]
    patient_text = " ".join(p.lower() for p in patient_phenotypes)

    for ph in _PHENOTYPE_GATE_KEYWORDS:
        for prep in ("in ", "for ", "with "):
            if prep + ph in window:
                # Found gating. If the patient matches it, keep the trigger.
                if ph.rstrip("s") in patient_text or ph in patient_text:
                    return False
                return True
    return False


def _drug_risk_level(
    recommendation: str,
    classification: str,
    patient_phenotypes: list[str] | None = None,
) -> str:
    """Map a drug recommendation to the 4-level risk vocabulary.

    The CPIC `classification` field describes the strength of the underlying
    *guideline evidence* — NOT whether this patient needs to take action.
    A 'Strong' classification paired with "use standard dose" means "strong
    evidence that the normal dose is appropriate", which is patient-normal.
    So we read the recommendation TEXT first and fall back to classification
    only when the text is silent.

    `patient_phenotypes` carries the per-gene phenotypes that this annotation
    applies to. They're used to suppress action keywords that are gated to a
    different phenotype than the patient has (e.g. "consider alternative in
    poor metabolizers" should not flag Action for an intermediate metabolizer).
    """
    cls = classification.lower().strip()
    rec = (recommendation or "").lower()
    patient_phenotypes = patient_phenotypes or []

    # 0. Catch CPIC reassurance phrasings that contain an action keyword but
    #    are negating it ("No reason to avoid", "not contraindicated", …).
    if any(phrase in rec for phrase in _FALSE_ACTION_PHRASES):
        return "normal"

    # 1. Explicit "avoid / contraindicated / consider alternative" → action,
    #    but suppress action keywords that are gated to a phenotype the
    #    patient doesn't have.
    for kw in _ACTION_KEYWORDS:
        pos = 0
        while True:
            idx = rec.find(kw, pos)
            if idx == -1:
                break
            if not _action_kw_gated_to_other_phenotype(rec, idx, patient_phenotypes):
                return "action"
            pos = idx + 1

    # 2. Explicit "use standard / label-recommended / usual dose" → normal,
    #    even when CPIC classification is "Strong".
    if any(phrase in rec for phrase in _NORMAL_PHRASES):
        return "normal"

    # 3. Dose-adjustment / monitoring language → review
    if any(kw in rec for kw in _REVIEW_KEYWORDS):
        return "review"

    # 4. Fall back to classification when the text gives no clear signal.
    #    Default Strong/Moderate/Actionable/Informative to "review", not
    #    "action" — safer than over-warning when the text is ambiguous.
    if cls in ("no action", "no action needed"):
        return "normal"
    if cls in ("strong", "actionable", "moderate", "informative"):
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
