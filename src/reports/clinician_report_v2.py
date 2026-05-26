"""
Clinician report v2 — risk-grouped redesign.

Differences from v1:
- Genes are grouped by risk (Action / Review / Normal / No-Data) like the
  v2 patient report, instead of by functional category. Functional category
  (Phase I / Phase II / Transporters / HLA / Other) becomes a tag next to
  each gene name within the risk bucket.
- Per-gene CPIC/DPWG drug recommendations are rendered visibly (no toggle).
- All technical fields (diplotype, allele names + functions, coverage stats)
  preserved per gene.
- Clinical Priority Summary up top.
- Consolidated caveats at the end (same as v1).
- The full Drug Recommendations table at the bottom is kept as a
  cross-reference but rendered inside a single collapsed `<details>` so
  it doesn't dominate the page (the per-gene rendering already covers it).

The shared infrastructure (drug-derived gene risk, filter to CPIC+DPWG,
gene→drug map, coverage summary, consolidated caveats) is imported from
the v1 module rather than copied.
"""

import html
import re
from datetime import datetime
from pathlib import Path

from markupsafe import Markup

from config.settings import (
    ACTION_LABELS, ACTION_SYMBOLS as _DEFAULT_SYMBOLS, APP_TITLE,
    CATEGORY_ORDER, CLINICIAN_DISCLAIMER, GENE_CATEGORIES,
    GENE_DESCRIPTIONS, RISK_PRIORITY,
)

# Match the patient v2 visual vocabulary: exclamation for Action, monochrome
# eye SVG for Review. The eye is sized in `em` and strokes in `currentColor`
# so it inherits the surrounding text's size and color (which is itself
# risk-level-tinted via the parent CSS class). v1 keeps the original symbols
# for side-by-side comparison.
_EYE_SVG = Markup(
    '<svg viewBox="0 0 24 24" width="1em" height="1em" fill="none" '
    'stroke="currentColor" stroke-width="2.2" stroke-linecap="round" '
    'stroke-linejoin="round" style="vertical-align:-0.12em">'
    '<path d="M2 12s3.5-7 10-7 10 7 10 7-3.5 7-10 7-10-7-10-7z"/>'
    '<circle cx="12" cy="12" r="3"/>'
    '</svg>'
)
ACTION_SYMBOLS = {**_DEFAULT_SYMBOLS, "action": "!", "review": _EYE_SVG}
from pharmcat.output_parser import ParsedResults
from reports._render import html_and_pdf
from reports.clinician_report import (
    _build_gene_drug_map,
    _consolidated_caveats,
    _coverage_summary,
    _effective_gene_risk,
    _filter_cpic_dpwg,
    _group_drugs_by_category,
    _is_normal_pheno,
)


# Role tag for each gene comes from `GENE_PROTEIN_TYPE` (already on each
# parsed gene as `g.protein_type`). Same labels as the v2 patient report:
# "Drug breakdown", "Drug transporter", "Drug target",
# "Cystic-fibrosis drug response", "Anesthesia response", "Immune response",
# "Hepatitis-C drug response", "Antibiotic sensitivity",
# "Red-blood-cell drug response", "Blood clotting".


# ---------------------------------------------------------------------------
# Mechanism briefs — short technical sentence linking gene+phenotype to drug
# effect. Surfaced inside each Medication Guidance card so the prescriber sees
# the *why* without having to bounce to the citing guideline.
#
# Keys are matched on (gene, normalized_phenotype). The normalizer lowercases
# the phenotype and strips surrounding whitespace; callers should keep keys
# in lower case.
# ---------------------------------------------------------------------------

MECHANISM_BRIEFS: dict[tuple[str, str], str] = {
    # --- CFTR -----------------------------------------------------------
    ("CFTR", "ivacaftor non-responsive in cf patients"):
        "Variant not in CFTR ivacaftor-responsive set (gating / residual-function "
        "mutations); potentiator therapy unlikely to improve channel activity.",
    ("CFTR", "ivacaftor responsive in cf patients"):
        "Gating- or residual-function CFTR mutation; ivacaftor potentiates "
        "channel open-probability and improves conductance.",

    # --- CYP3A5 ---------------------------------------------------------
    ("CYP3A5", "poor metabolizer"):
        "Non-expresser (*3/*3 or *6/*6 or *7/*7); CYP3A5-mediated tacrolimus "
        "clearance is ~50% of expressers — standard mg/kg risks supratherapeutic exposure.",
    ("CYP3A5", "intermediate metabolizer"):
        "Heterozygous expresser; intermediate CYP3A5 activity — typically slightly "
        "higher tacrolimus exposure than full expressers.",
    ("CYP3A5", "normal metabolizer"):
        "Expresser; CYP3A5-mediated clearance is preserved.",

    # --- CYP2C9 ---------------------------------------------------------
    ("CYP2C9", "poor metabolizer"):
        "Severely reduced CYP2C9 activity; substrate clearance impaired — "
        "supratherapeutic levels at standard dose for narrow-index drugs (phenytoin, warfarin, NSAIDs).",
    ("CYP2C9", "intermediate metabolizer"):
        "Reduced CYP2C9 activity; modestly elevated exposure to CYP2C9 substrates.",
    ("CYP2C9", "normal metabolizer"):
        "Preserved CYP2C9 activity; standard dosing applies for CYP2C9 substrates.",

    # --- CYP2C19 --------------------------------------------------------
    ("CYP2C19", "poor metabolizer"):
        "Loss-of-function CYP2C19; prodrugs (clopidogrel) under-activated, "
        "parent-drug exposure (PPIs, voriconazole, SSRIs) elevated.",
    ("CYP2C19", "intermediate metabolizer"):
        "Reduced CYP2C19 activity; partial under-activation of prodrugs and "
        "modest accumulation of parent-drug substrates.",
    ("CYP2C19", "rapid metabolizer"):
        "Increased CYP2C19 activity (*17 carrier); accelerated metabolism — "
        "lower parent-drug exposure for PPIs/voriconazole/SSRIs.",
    ("CYP2C19", "ultrarapid metabolizer"):
        "Very high CYP2C19 activity (*17/*17); markedly reduced parent-drug "
        "exposure for PPIs and voriconazole — therapeutic failure risk.",
    ("CYP2C19", "normal metabolizer"):
        "Preserved CYP2C19 activity; standard dosing applies.",

    # --- CYP2D6 ---------------------------------------------------------
    ("CYP2D6", "poor metabolizer"):
        "Loss-of-function CYP2D6; codeine/tramadol under-activated, "
        "TCA / atomoxetine / certain antipsychotics over-exposed at standard dose.",
    ("CYP2D6", "intermediate metabolizer"):
        "Reduced CYP2D6 activity; modest under-activation of prodrugs and "
        "increased exposure to active-drug substrates.",
    ("CYP2D6", "ultrarapid metabolizer"):
        "Gene duplication / high-activity alleles; codeine/tramadol "
        "over-activation (opioid-toxicity risk in some patients).",

    # --- DPYD -----------------------------------------------------------
    ("DPYD", "poor metabolizer"):
        "Severe DPD deficiency; fluoropyrimidines (5-FU, capecitabine) "
        "contraindicated — life-threatening toxicity risk.",
    ("DPYD", "intermediate metabolizer"):
        "Partial DPD deficiency; fluoropyrimidine dose reduction (typically 50%) advised.",

    # --- NAT2 (acetylator vocabulary; PharmCAT sometimes uses "metabolizer") ---
    ("NAT2", "slow acetylator"):
        "Reduced N-acetylation; prolonged exposure to hydralazine / isoniazid "
        "→ lupus-like syndrome and peripheral neuropathy risk.",
    ("NAT2", "slow metabolizer"):
        "Reduced N-acetylation; prolonged exposure to hydralazine / isoniazid "
        "→ lupus-like syndrome and peripheral neuropathy risk.",
    ("NAT2", "poor metabolizer"):
        "Severely reduced N-acetylation; markedly prolonged exposure to "
        "hydralazine / isoniazid — high lupus-like syndrome and neuropathy risk.",
    ("NAT2", "intermediate acetylator"):
        "Intermediate N-acetylation; modestly elevated parent-drug exposure.",
    ("NAT2", "intermediate metabolizer"):
        "Intermediate N-acetylation; modestly elevated parent-drug exposure.",
    ("NAT2", "rapid acetylator"):
        "Accelerated N-acetylation; risk of subtherapeutic hydralazine "
        "exposure at standard doses.",
    ("NAT2", "rapid metabolizer"):
        "Accelerated N-acetylation; risk of subtherapeutic hydralazine "
        "exposure at standard doses.",

    # --- NUDT15 / TPMT (thiopurine pathway) -----------------------------
    ("NUDT15", "poor metabolizer"):
        "Severely reduced NUDT15; cannot inactivate cytotoxic thiopurine "
        "triphosphates — high risk of fatal myelosuppression at standard doses.",
    ("NUDT15", "intermediate metabolizer"):
        "Reduced NUDT15 activity; thiopurine triphosphates accumulate — "
        "myelosuppression risk; CPIC suggests 30–50% starting-dose reduction.",
    ("NUDT15", "normal metabolizer"):
        "Preserved NUDT15 activity; standard thiopurine dosing applies on the NUDT15 axis.",
    ("TPMT", "poor metabolizer"):
        "Severely reduced TPMT; thiopurine S-methylation impaired — "
        "TGN accumulation and high myelosuppression risk at standard doses.",
    ("TPMT", "intermediate metabolizer"):
        "Reduced TPMT activity; modest TGN accumulation — CPIC suggests "
        "30–80% starting-dose reduction.",
    ("TPMT", "normal metabolizer"):
        "Preserved TPMT activity; standard thiopurine dosing applies on the TPMT axis.",

    # --- UGT1A1 ---------------------------------------------------------
    ("UGT1A1", "poor metabolizer"):
        "Severely reduced UGT1A1 glucuronidation; SN-38 (active irinotecan "
        "metabolite) accumulates — high risk of neutropenia / severe diarrhea.",
    ("UGT1A1", "intermediate metabolizer"):
        "Reduced UGT1A1 glucuronidation (*28/*28 or similar); dose-dependent "
        "SN-38 accumulation — CPIC advises reducing irinotecan starting dose ≥150 mg/m².",
    ("UGT1A1", "normal metabolizer"):
        "Preserved UGT1A1 glucuronidation; standard irinotecan dosing applies.",

    # --- VKORC1 ---------------------------------------------------------
    ("VKORC1", "-1639 aa"):
        "High-sensitivity VKORC1 haplotype (low VKORC1 expression); "
        "lower warfarin/acenocoumarol dose requirement — standard dose risks supratherapeutic INR.",
    ("VKORC1", "-1639 ag"):
        "Intermediate VKORC1 expression; modest dose reduction per "
        "genotype-guided dosing algorithms (e.g., IWPC).",
    ("VKORC1", "-1639 ga"):
        "Intermediate VKORC1 expression; modest dose reduction per "
        "genotype-guided dosing algorithms (e.g., IWPC).",
    ("VKORC1", "-1639 gg"):
        "Wild-type VKORC1 expression; standard or higher coumarin dose requirement.",

    # --- SLCO1B1 / ABCG2 (statin transporters) --------------------------
    ("SLCO1B1", "poor function"):
        "Severely reduced hepatic OATP1B1 uptake (*5/*5 etc.); statin "
        "(simvastatin > others) plasma exposure markedly increased — myopathy risk.",
    ("SLCO1B1", "decreased function"):
        "Reduced hepatic OATP1B1 uptake; statin exposure increased — "
        "modest myopathy risk; consider lower-dose or non-OATP1B1 statin.",
    ("SLCO1B1", "increased function"):
        "Increased OATP1B1 uptake; statin exposure may be reduced.",
    ("SLCO1B1", "possible decreased function"):
        "Possibly reduced OATP1B1 uptake; modest statin exposure increase possible.",
    ("ABCG2", "poor function"):
        "Severely reduced ABCG2 efflux; rosuvastatin plasma exposure markedly "
        "increased — myopathy risk.",
    ("ABCG2", "decreased function"):
        "Reduced ABCG2 efflux; rosuvastatin exposure modestly increased.",

    # --- CACNA1S / RYR1 (anesthesia) ------------------------------------
    ("CACNA1S", "malignant hyperthermia susceptibility"):
        "Pathogenic CACNA1S variant; malignant-hyperthermia-susceptible — "
        "volatile anesthetics and succinylcholine are contraindicated.",
    ("RYR1", "malignant hyperthermia susceptibility"):
        "Pathogenic RYR1 variant; malignant-hyperthermia-susceptible — "
        "volatile anesthetics and succinylcholine are contraindicated.",

    # --- HLA-B / HLA-A --------------------------------------------------
    ("HLA-B", "positive"):
        "Carries the HLA-B risk allele; severe cutaneous adverse reaction "
        "(SJS/TEN, DRESS) risk with the associated drug.",
    ("HLA-A", "positive"):
        "Carries the HLA-A risk allele; severe cutaneous adverse reaction "
        "risk with the associated drug.",

    # --- IFNL3 ----------------------------------------------------------
    ("IFNL3", "favorable response genotype"):
        "Favorable IFNL3 (rs12979860 CC) genotype; higher sustained virologic "
        "response to interferon-based hepatitis C therapy.",
    ("IFNL3", "unfavorable response genotype"):
        "Unfavorable IFNL3 (CT/TT) genotype; reduced sustained virologic "
        "response to interferon-based hepatitis C therapy.",

    # --- MT-RNR1 --------------------------------------------------------
    ("MT-RNR1", "increased risk of aminoglycoside-induced hearing loss"):
        "Mitochondrial m.1555A>G or m.1494C>T variant; aminoglycosides "
        "cause irreversible ototoxicity even at therapeutic doses — avoid.",

    # --- G6PD -----------------------------------------------------------
    ("G6PD", "deficient"):
        "G6PD-deficient erythrocytes; oxidant drugs (rasburicase, primaquine, "
        "dapsone) trigger acute hemolytic anemia — contraindicated.",
    ("G6PD", "deficient with cnshb"):
        "Deficient with chronic non-spherocytic hemolytic anemia; same "
        "oxidant-drug contraindication, baseline hemolysis present.",
    ("G6PD", "variable"):
        "Variable G6PD activity (heterozygous female); X-inactivation may "
        "produce a deficient cell population — caution with oxidant drugs.",
}


def _norm_pheno_key(pheno: str) -> str:
    return (pheno or "").strip().lower()


def _mechanism_for(gene: str, phenotype: str) -> str | None:
    """Return the brief mechanism sentence for (gene, phenotype), or None
    if no specific entry — caller can fall back to a generic phrase."""
    return MECHANISM_BRIEFS.get((gene, _norm_pheno_key(phenotype)))


def _generic_mechanism(phenotype: str) -> str:
    """Fallback when no gene-specific brief is registered. Tries to convey
    the directionality of the phenotype without claiming gene-specific
    pharmacology that we haven't curated."""
    p = _norm_pheno_key(phenotype)
    if "poor" in p:
        return "Severely reduced activity; substrate handling impaired at standard dose."
    if "intermediate" in p:
        return "Reduced activity; exposure to substrates of this gene is modestly altered."
    if "rapid" in p or "ultrarapid" in p:
        return "Increased activity; exposure to substrates of this gene may be reduced."
    if "normal" in p:
        return "Preserved activity on this axis; no genotype-driven adjustment."
    return "Phenotype contributes to the recommendation per the citing guideline."


def _number_citations_in_body(
    parsed: ParsedResults,
    grouped: dict[str, list[dict]],
) -> list[dict]:
    """Walk every source entry inside grouped's Action/Review gene
    monographs, assign sequential citation numbers in order of first
    appearance, and attach `cite_numbers` to each source entry. Returns
    the ordered References list — `[{"number": N, "citation": Citation},
    ...]` — containing only citations actually used by surviving body
    content. Citations attached to drug-source rows that were filtered
    out (e.g. Normal-phenotype guideline arms) are excluded.

    Side-effect: mutates source entries inside `grouped` to add the
    `cite_numbers` field that the template renders as `[3,7]` markers."""
    canonical = {c.pmid: c for c in parsed.citations if c.pmid}
    key_to_num: dict[str, int] = {}
    order: list[str] = []

    for bucket in ("action", "review"):
        for g in _iter_genes(grouped.get(bucket, [])):
            for d in (g.get("drugs") or []):
                for s in d.get("sources", []):
                    nums = []
                    for k in s.get("citation_keys") or []:
                        if k not in key_to_num:
                            key_to_num[k] = len(order) + 1
                            order.append(k)
                        nums.append(key_to_num[k])
                    s["cite_numbers"] = nums

    return [
        {"number": key_to_num[k], "citation": canonical[k]}
        for k in order if k in canonical
    ]


def _build_drug_index(parsed: ParsedResults) -> list[dict]:
    """Compact drug-first index. One row per (drug, contributing-gene)
    pair, sorted by risk then category. Used at the top of the report so
    the prescriber can scan affected medications and jump to the gene
    monograph below. Iterates all parsed.drugs (every source) so that a
    drug whose only Action/Review annotation comes from an FDA label
    still surfaces here."""
    rows: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for d in parsed.drugs:
        if d.risk_level not in ("action", "review"):
            continue
        for gene in d.affected_genes:
            pheno = (d.phenotypes or {}).get(gene, "")
            if _is_normal_pheno(pheno):
                continue
            key = (d.drug, gene)
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "drug": d.drug,
                "risk_level": d.risk_level,
                "symbol": ACTION_SYMBOLS[d.risk_level],
                "label": ACTION_LABELS[d.risk_level],
                "therapeutic_category": d.therapeutic_category,
                "gene": gene,
                "phenotype": _display_phenotype(gene, pheno),
                "anchor": f"gene-{gene}",
            })
    rows.sort(key=lambda r: (
        RISK_PRIORITY.get(r["risk_level"], 3),
        r["therapeutic_category"] or "",
        r["drug"],
        r["gene"],
    ))
    return rows


_WHITESPACE_RE = re.compile(r"\s+")


def _clean_recommendation(text: str) -> Markup:
    """Normalize a PharmCAT recommendation string for safe inline HTML
    rendering. PharmCAT's report.json mixes two encodings:

    - FDA Label / FDA PGx Association entries pre-encode quotes as
      `&quot;` HTML entities (and similar). Without unescaping, Jinja
      autoescape double-encodes the entity and the reader sees the raw
      `&quot;` text.
    - DPWG Guideline entries embed real HTML markup (`<ul>`, `<li>`,
      `<br />`) plus literal newlines for paragraph structure. With
      autoescape on, the markup renders as literal text instead of a
      bullet list.

    We unescape entities, collapse runs of whitespace (newlines included)
    to a single space, and wrap the result as `Markup` so Jinja preserves
    the structural tags. PharmCAT is a trusted source — the only HTML it
    emits is the block-level list/break markup above — so this is safe."""
    if not text:
        return Markup("")
    cleaned = html.unescape(text)
    cleaned = _WHITESPACE_RE.sub(" ", cleaned).strip()
    return Markup(cleaned)


def _drugs_for_gene_block(parsed: ParsedResults, gene: str) -> list[dict]:
    """All drugs whose Action/Review classification is driven by this
    gene's non-normal phenotype. Grouped by drug name; each drug's
    `sources` list stacks every guideline annotation (CPIC, DPWG, FDA
    Label, etc.) so the clinician sees source-specific recommendation
    text inline within the gene monograph."""
    by_drug: dict[str, dict] = {}
    for d in parsed.drugs:
        if d.risk_level not in ("action", "review"):
            continue
        if gene not in d.affected_genes:
            continue
        if _is_normal_pheno((d.phenotypes or {}).get(gene, "")):
            continue
        entry = by_drug.setdefault(d.drug, {
            "drug": d.drug,
            "therapeutic_category": d.therapeutic_category,
            "risk_level": d.risk_level,
            "symbol": ACTION_SYMBOLS[d.risk_level],
            "sources": [],
        })
        if (
            RISK_PRIORITY.get(d.risk_level, 3)
            < RISK_PRIORITY.get(entry["risk_level"], 3)
        ):
            entry["risk_level"] = d.risk_level
            entry["symbol"] = ACTION_SYMBOLS[d.risk_level]
        entry["sources"].append({
            "source": d.source,
            "classification": d.classification,
            "recommendation": _clean_recommendation(d.recommendation),
            # PMID-bearing citations only. FDA-label / Drugs@FDA "citations"
            # have no PMID and no clickable URL, so they offer no
            # navigability in the References section — their provenance is
            # already covered by the visible source tag on the drug row.
            "citation_keys": [c.pmid for c in d.citations if c.pmid],
        })
    out = list(by_drug.values())
    for e in out:
        e["sources"].sort(key=lambda s: s["source"])
        # Partition into CPIC/DPWG (actionable, shown by default) and FDA
        # (legally binding / regulatory framing, often vague — collapsed
        # behind a <details> toggle in the template). If a drug has only
        # FDA sources (e.g. irinotecan, warfarin), promote them so the
        # drug block always shows something at the top level.
        primary = [s for s in e["sources"] if "FDA" not in s["source"]]
        fda = [s for s in e["sources"] if "FDA" in s["source"]]
        if not primary:
            primary, fda = fda, []
        e["primary_sources"] = primary
        e["fda_sources"] = fda
    out.sort(key=lambda d: (RISK_PRIORITY.get(d["risk_level"], 3), d["drug"]))
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_clinician_report_v2(
    parsed: ParsedResults,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Render the v2 clinician report (HTML + best-effort PDF)."""
    context = _build_context(parsed)
    return html_and_pdf(
        "clinician_report_v2.html", context,
        Path(output_dir), "clinician_report_v2",
    )


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def _build_context(parsed: ParsedResults) -> dict:
    coverage = _coverage_summary(parsed)
    cpic_dpwg = _filter_cpic_dpwg(parsed)
    gene_drug_map = _build_gene_drug_map(cpic_dpwg)
    effective_risk = _effective_gene_risk(parsed)

    # Drug Index — compact navigation table. One row per (drug, contributing
    # gene) pair, with an anchor linking to the gene monograph in the body.
    # Built first so the Action/Review gene buckets can be ordered to match
    # the drug-index order (the clinician reads the index, then jumps to
    # gene monographs in the same sequence).
    drug_index = _build_drug_index(parsed)

    grouped = _group_by_risk(parsed, effective_risk, drug_index)
    # Filter + number citations to only those backing surviving body content.
    # Mutates source entries inside `grouped` so the template can stamp
    # `[N]` markers next to each cited recommendation.
    numbered_citations = _number_citations_in_body(parsed, grouped)
    consolidated_caveats = _consolidated_caveats(parsed)
    drugs_by_category = _group_drugs_by_category(cpic_dpwg)

    # Gene-overview summary: just the gene names per risk bucket, in a flat
    # alphabetical list. Used by the four single-line "Gene Overview" blocks
    # at the top of the report so the prescriber sees the full panel at a
    # glance before diving into the monographs.
    gene_summary_by_risk = {
        bucket: sorted(g["gene"] for g in _iter_genes(grouped[bucket]))
        for bucket in ("action", "review", "normal", "no_data")
    }

    return {
        "title": "Clinician report",
        "app_title": APP_TITLE,
        "report_date": datetime.now().strftime("%B %d, %Y"),
        "sample_id": parsed.sample_id or "Unknown",
        "metadata": parsed.metadata,
        "disclaimer": CLINICIAN_DISCLAIMER,
        "messages": parsed.messages,
        "action_symbols": ACTION_SYMBOLS,
        "action_labels": ACTION_LABELS,
        "coverage_summary": coverage,
        "gene_summary_by_risk": gene_summary_by_risk,
        "grouped_by_risk": grouped,
        "drug_index": drug_index,
        "drugs_by_category": drugs_by_category,
        "all_drugs": cpic_dpwg,
        "numbered_citations": numbered_citations,
        "consolidated_caveats": consolidated_caveats,
        "category_order": CATEGORY_ORDER,
    }


# ---------------------------------------------------------------------------
# Risk-bucket grouping (clinician-flavored: keeps all technical fields)
# ---------------------------------------------------------------------------

def _group_by_risk(
    parsed: ParsedResults,
    effective_risk: dict[str, str],
    drug_index: list[dict] | None = None,
) -> dict[str, list[dict]]:
    """Group genes by risk bucket.

    Action/Review come back as **flat** gene lists ordered by first appearance
    in `drug_index`, so the gene monograph order matches the order the
    clinician scans in the Drug Index above. Normal/No-Data come back as
    protein-type subgroups (`[{protein_type, genes}, ...]`) since drug order
    doesn't apply there.
    """
    action: list[dict] = []
    review: list[dict] = []
    normal: list[dict] = []
    no_data: list[dict] = []

    for g in parsed.definitive_genes:
        risk = effective_risk.get(g.gene, g.risk_level)
        entry = _enrich_definitive(g, risk, parsed)
        {"action": action, "review": review, "normal": normal,
         "nodata": no_data}.get(risk, no_data).append(entry)

    for g in parsed.ambiguous_genes:
        risk = effective_risk.get(g.gene, g.risk_level)
        entry = _enrich_ambiguous(g, risk, parsed)
        {"action": action, "review": review, "normal": normal,
         "nodata": no_data}.get(risk, no_data).append(entry)

    for g in parsed.no_call_genes:
        no_data.append(_enrich_no_call(g))

    # Sort Action/Review by first appearance in drug_index (a gene with no
    # drug-index row sinks to the end, then sorted by gene name as tiebreak).
    action_order = _gene_order_from_drug_index(drug_index or [], "action")
    review_order = _gene_order_from_drug_index(drug_index or [], "review")
    action.sort(key=lambda g: (action_order.get(g["gene"], 10_000), g["gene"]))
    review.sort(key=lambda g: (review_order.get(g["gene"], 10_000), g["gene"]))

    return {
        "action": action,
        "review": review,
        "normal": _subgroup_by_protein_type(normal),
        "no_data": _subgroup_by_protein_type(no_data),
    }


def _gene_order_from_drug_index(drug_index: list[dict], risk_level: str) -> dict[str, int]:
    """Map each gene to its first-row position among `drug_index` rows of
    the given risk level. The drug_index is already sorted by therapeutic
    category and drug name, so this preserves that ordering."""
    out: dict[str, int] = {}
    for i, r in enumerate(drug_index):
        if r.get("risk_level") != risk_level:
            continue
        out.setdefault(r["gene"], i)
    return out


def _iter_genes(bucket_value: list[dict]):
    """Iterate gene dicts from a risk bucket, transparently handling both
    flat lists (Action/Review) and protein-type subgroup lists
    (Normal/No-Data)."""
    if not bucket_value:
        return
    if isinstance(bucket_value[0], dict) and "genes" in bucket_value[0]:
        for sub in bucket_value:
            yield from sub["genes"]
    else:
        yield from bucket_value


def _subgroup_by_protein_type(genes: list[dict]) -> list[dict]:
    """Group genes inside a risk bucket by their protein_type (role), so the
    clinician sees genes clustered by function within each risk category.
    Returns a list of {protein_type, genes} dicts, sub-groups sorted
    alphabetically by protein_type label, genes within each sub-group sorted
    by gene name."""
    from collections import defaultdict
    buckets: dict[str, list[dict]] = defaultdict(list)
    for g in genes:
        buckets[g.get("protein_type") or "Other"].append(g)
    out: list[dict] = []
    for pt in sorted(buckets.keys()):
        out.append({
            "protein_type": pt,
            "genes": sorted(buckets[pt], key=lambda d: d["gene"]),
        })
    return out


# ---------------------------------------------------------------------------
# Gene enrichment (keep all technical clinician-facing fields)
# ---------------------------------------------------------------------------

# Genes whose raw PharmCAT "phenotype" is a genotype string rather than a
# functional label. Translated to a functional phrase so the clinician's
# phenotype pill reads as the clinical meaning, not the SNP shorthand.
_PHENOTYPE_DISPLAY_OVERRIDES: dict[tuple[str, str], str] = {
    ("VKORC1", "-1639 AA"): "Highly increased coumarin sensitivity",
    ("VKORC1", "-1639 GA"): "Increased coumarin sensitivity",
    ("VKORC1", "-1639 AG"): "Increased coumarin sensitivity",
    ("VKORC1", "-1639 GG"): "Normal coumarin sensitivity",
    # CYP4F2 *1/*1 and IFNL3 rs12979860 C/C are reclassified to Normal in
    # the parser (see _benign_call). The raw phenotype is "No Result" /
    # "n/a" because CPIC has no diplotype-level phenotype label — these
    # overrides give the pill a meaningful label that matches the bucket.
    ("CYP4F2", "No Result"): "Normal function",
    ("IFNL3", "n/a"): "Favorable variant",
}


def _display_phenotype(gene: str, raw: str | None) -> str:
    """Override raw phenotype strings for genes where PharmCAT returns a
    genotype-style label that isn't clinician-readable on its own. Anything
    not in the override map passes through unchanged."""
    if not raw:
        return raw or ""
    return _PHENOTYPE_DISPLAY_OVERRIDES.get((gene, raw.strip()), raw)


def _enrich_definitive(g, risk: str, parsed: ParsedResults) -> dict:
    mechanism = (
        _mechanism_for(g.gene, g.phenotype)
        or _generic_mechanism(g.phenotype)
    )
    return {
        "gene": g.gene,
        "anchor": f"gene-{g.gene}",
        "kind": "definitive",
        "protein_type": g.protein_type,
        "description": GENE_DESCRIPTIONS.get(g.gene),
        "diplotype": g.diplotype,
        "phenotype": _display_phenotype(g.gene, g.phenotype),
        "star_alleles": g.star_alleles,
        "allele1_name": g.allele1_name,
        "allele1_function": g.allele1_function,
        "allele2_name": g.allele2_name,
        "allele2_function": g.allele2_function,
        "risk_level": risk,
        "symbol": ACTION_SYMBOLS[risk],
        "label": ACTION_LABELS[risk],
        "positions_found": g.positions_found,
        "positions_missing": g.positions_missing,
        "related_drugs": g.related_drugs,
        "mechanism": mechanism,
        "drugs": (
            _drugs_for_gene_block(parsed, g.gene)
            if risk in ("action", "review") else []
        ),
    }


def _enrich_ambiguous(g, risk: str, parsed: ParsedResults) -> dict:
    return {
        "gene": g.gene,
        "anchor": f"gene-{g.gene}",
        "kind": "ambiguous",
        "protein_type": g.protein_type,
        "description": GENE_DESCRIPTIONS.get(g.gene),
        "diplotype_count": g.diplotype_count,
        "phenotype_range": g.phenotype_range,
        "harmful_phenotypes": g.harmful_phenotypes,
        "risk_level": risk,
        "symbol": ACTION_SYMBOLS[risk],
        "label": ACTION_LABELS[risk],
        "positions_found": g.positions_found,
        "positions_missing": g.positions_missing,
        "actionable_drugs": g.actionable_drugs,
        "related_drugs": g.related_drugs,
        "mechanism": None,
        "drugs": (
            _drugs_for_gene_block(parsed, g.gene)
            if risk in ("action", "review") else []
        ),
    }


def _enrich_no_call(g) -> dict:
    return {
        "gene": g.gene,
        "anchor": f"gene-{g.gene}",
        "kind": "no_call",
        "protein_type": g.protein_type,
        "description": GENE_DESCRIPTIONS.get(g.gene),
        "risk_level": "nodata",
        "symbol": ACTION_SYMBOLS["nodata"],
        "label": ACTION_LABELS["nodata"],
        "positions_missing": g.positions_missing,
        "affected_drugs": g.affected_drugs,
        "related_drugs": g.related_drugs,
        "mechanism": None,
        "drugs": [],
    }
