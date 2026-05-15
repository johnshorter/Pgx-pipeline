"""
Patient report v2 — proposed redesign based on patient-readability feedback.

Differences from v1:
- Genetics jargon stripped from patient-facing strings (no star alleles, no
  activity scores, no diplotype counts, no raw harmful-phenotype names, no
  PharmCAT version in the header).
- Phenotype subtitles in overview cards are sanitized to the PHENOTYPE_SHORT
  enum; anything else falls back to "See details below".
- "Your next steps" block right after the dashboard with personalized
  bullets driven by the actual results.
- Medication Guidance is rendered BEFORE Genes-by-Function (action precedes
  education).
- Per-action drug gets a one-line "what might happen" consequence.
- Common drugs gain brand names in parentheses.
- New "What this report does NOT tell you" section.
- New "Learn more" external resources footer block.
- No-call affected-drug lists are collapsed into <details>.

The original `patient_report.py` is unchanged.
"""

from collections import OrderedDict, defaultdict
from datetime import date
from pathlib import Path

from markupsafe import Markup

from config.settings import (
    ACTION_LABELS, ACTION_SYMBOLS as _DEFAULT_SYMBOLS, APP_TITLE,
    CATEGORY_ORDER, GENE_CATEGORIES, PATIENT_DISCLAIMER, RISK_PRIORITY,
    drug_category, gene_tissue_icon,
)

# Phenotype strings that are "normal-equivalent" — when a drug recommendation
# lists a gene with one of these phenotypes for the patient, that gene is NOT
# driving the recommendation's risk classification, so it shouldn't get its
# own risk bumped from the drug.
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

# v2 uses an exclamation mark for Action and a monochrome eye SVG for Review
# instead of the default ▲ / ◆ shapes. The eye is sized in `em` units and
# strokes in `currentColor`, so it scales and recolors with the surrounding
# text. v1 keeps the original symbols for side-by-side comparison.
_EYE_SVG = Markup(
    '<svg viewBox="0 0 24 24" width="1em" height="1em" fill="none" '
    'stroke="currentColor" stroke-width="2.2" stroke-linecap="round" '
    'stroke-linejoin="round" style="vertical-align:-0.12em">'
    '<path d="M2 12s3.5-7 10-7 10 7 10 7-3.5 7-10 7-10-7-10-7z"/>'
    '<circle cx="12" cy="12" r="3"/>'
    '</svg>'
)
ACTION_SYMBOLS = {**_DEFAULT_SYMBOLS, "action": "!", "review": _EYE_SVG}
from pharmcat.output_parser import (
    AmbiguousGene, DefinitiveGene, NoCallGene, ParsedResults,
)
from reports._render import html_and_pdf
from reports.patient_report import (
    DRUG_GUIDANCE, PHENOTYPE_PLAIN, PHENOTYPE_SHORT, _explain_phenotype,
)

# ---------------------------------------------------------------------------
# v2-only content
# ---------------------------------------------------------------------------

# Trimmed, single-sentence category descriptions for the deep-dive section.
CATEGORY_DESC_SHORT: dict[str, str] = {
    "Phase I Metabolism (CYP Enzymes)":
        "Liver enzymes that break down many common medications.",
    "Phase II Metabolism":
        "Enzymes that make drugs easier for the body to eliminate.",
    "Drug Transporters":
        "Proteins that move drugs into and out of cells.",
    "Immune Markers (HLA)":
        "Genes that can trigger severe allergic-type reactions to specific drugs.",
    "Other Pharmacogenes":
        "Other genes that influence how the body responds to certain drugs.",
}

# Common brand names for the drugs PharmCAT reports on. Not exhaustive —
# only the drugs a patient is likely to recognize from a pill bottle.
DRUG_BRAND_NAMES: dict[str, str] = {
    "warfarin": "Coumadin, Jantoven",
    "clopidogrel": "Plavix",
    "simvastatin": "Zocor",
    "atorvastatin": "Lipitor",
    "rosuvastatin": "Crestor",
    "pravastatin": "Pravachol",
    "fluvastatin": "Lescol",
    "lovastatin": "Mevacor, Altoprev",
    "omeprazole": "Prilosec",
    "pantoprazole": "Protonix",
    "esomeprazole": "Nexium",
    "lansoprazole": "Prevacid",
    "codeine": "(in Tylenol with Codeine)",
    "tramadol": "Ultram, ConZip",
    "hydrocodone": "Vicodin, Norco",
    "tacrolimus": "Prograf, Astagraf",
    "carbamazepine": "Tegretol, Carbatrol",
    "phenytoin": "Dilantin",
    "oxcarbazepine": "Trileptal",
    "amitriptyline": "Elavil",
    "citalopram": "Celexa",
    "escitalopram": "Lexapro",
    "fluoxetine": "Prozac",
    "paroxetine": "Paxil",
    "sertraline": "Zoloft",
    "venlafaxine": "Effexor",
    "fluvoxamine": "Luvox",
    "voriconazole": "Vfend",
    "allopurinol": "Zyloprim",
    "abacavir": "Ziagen",
    "azathioprine": "Imuran, Azasan",
    "mercaptopurine": "Purinethol, Purixan",
    "tamoxifen": "Soltamox",
    "ondansetron": "Zofran",
    "atomoxetine": "Strattera",
    "amphetamine": "Adderall",
    "metoprolol": "Lopressor, Toprol-XL",
    "propranolol": "Inderal",
    "carvedilol": "Coreg",
    "aripiprazole": "Abilify",
    "haloperidol": "Haldol",
    "clozapine": "Clozaril",
    "risperidone": "Risperdal",
    "donepezil": "Aricept",
    "galantamine": "Razadyne",
    "fluorouracil": "Adrucil, 5-FU",
    "capecitabine": "Xeloda",
}

# Per-therapeutic-category plain-English consequence for the Action bucket.
# Falls back to a generic statement if the category isn't listed.
# Each Action consequence is a (directive, reason_template) pair.
# - The directive leads the rationale and starts with an actionable verb
#   (Don't, Do not, Before, Talk). It's what the patient should *do*.
# - The reason explains *why* and includes a `{genes}` placeholder where the
#   patient's relevant gene names get interpolated (already wrapped in <strong>).
ACTION_CONSEQUENCE_BY_CATEGORY: dict[str, tuple[str, str]] = {
    "Antidepressants": (
        "Don't start at the usual dose without talking to your doctor first.",
        "{genes}: the usual dose may cause stronger side effects or may not "
        "relieve symptoms for you. Ask about a lower starting dose or a "
        "different antidepressant.",
    ),
    "Antiplatelet / Anticoagulant agents": (
        "Don't start without a dose-adjustment plan — discuss with your doctor.",
        "{genes}: the usual dose may not protect against blood clots for you. "
        "Ask about a different drug or an adjusted dosing strategy.",
    ),
    "Statins (cholesterol)": (
        "Don't start at the usual dose — ask about a lower dose or different statin.",
        "{genes}: you have a higher-than-usual risk of muscle pain or muscle "
        "damage at the standard dose.",
    ),
    "Proton pump inhibitors (stomach acid)": (
        "Don't expect the usual dose to fully work — ask about adjustments.",
        "{genes}: standard doses may not fully control stomach acid for you; "
        "a higher dose or a different acid-reducer may be needed.",
    ),
    "Antifungals": (
        "Don't start without a plan for blood-level monitoring.",
        "{genes}: drug levels can swing too high or too low at the standard "
        "dose for you.",
    ),
    "Cardiovascular": (
        "Don't start at the usual dose without discussing — ask about a lower dose.",
        "{genes}: the usual dose may be too strong for you, raising the risk of "
        "side effects like dizziness, fatigue, or slow heart rate.",
    ),
    "Immunosuppressants": (
        "Don't start without a plan for close blood-level monitoring.",
        "{genes}: drug levels are hard to predict at the usual dose for you; "
        "careful monitoring and dose adjustment are essential.",
    ),
    "Chemotherapy": (
        "Don't start at the usual dose — ask your oncologist about a reduced dose.",
        "{genes}: the standard dose can cause severe, potentially life-threatening "
        "side effects for you.",
    ),
    "Chemotherapy / Immunosuppressants": (
        "Don't start without dose reduction — ask your prescriber.",
        "{genes}: standard doses can cause severe drops in blood counts for you. "
        "Ask about a reduced starting dose or an alternative drug.",
    ),
    "Pain medications": (
        "Don't rely on this for pain control — ask about alternatives.",
        "{genes}: this drug may not relieve pain for you, or in rare cases may "
        "produce dangerously high levels of active drug; a different pain "
        "medication may be safer.",
    ),
    "Antivirals (HIV)": (
        "Do not take this drug — ask for an alternative.",
        "{genes}: you're at serious risk of a hypersensitivity (allergic-type) "
        "reaction.",
    ),
    "Antiepileptics": (
        "Do not take this drug — ask for an alternative.",
        "{genes}: you're at increased risk of a severe, potentially "
        "life-threatening skin reaction.",
    ),
    "Antipsychotics": (
        "Don't start at the usual dose — ask about a lower starting dose.",
        "{genes}: the usual dose may be too strong for you, with more sedation, "
        "weight gain, or movement-related side effects.",
    ),
    "Gout medications": (
        "Do not take this drug — ask for an alternative.",
        "{genes}: you're at serious risk of a severe, potentially "
        "life-threatening skin reaction.",
    ),
    "ADHD medications": (
        "Don't start at the usual dose — ask about a lower starting dose.",
        "{genes}: the usual dose may be too strong for you, with more side "
        "effects like jitteriness, trouble sleeping, or reduced appetite.",
    ),
    "Anesthetics": (
        "Before any procedure, tell your anesthesiologist about this result.",
        "{genes}: there is a potentially life-threatening risk of malignant "
        "hyperthermia with these drugs.",
    ),
    "Anti-nausea medications": (
        "Don't expect the usual dose to fully work — ask about alternatives.",
        "{genes}: the usual dose may not control nausea, or may affect heart "
        "rhythm.",
    ),
    "Cystic Fibrosis Modulators": (
        "Don't start this drug — ask your doctor about alternative CF therapies.",
        "{genes}: this drug is unlikely to work for your specific CFTR variant.",
    ),
}

ACTION_CONSEQUENCE_FALLBACK: tuple[str, str] = (
    "Talk to your prescriber before starting this drug.",
    "{genes}: the usual prescribing approach may not be right for you. "
    "Ask whether an alternative dose or drug would be better.",
)


# Patient-friendly display names for therapeutic categories — used as the
# 'Type' column in the Medication Guidance tables. Singular form, since
# each row is a single drug.
CATEGORY_DISPLAY: dict[str, str] = {
    "Antiplatelet / Anticoagulant agents": "Blood thinner",
    "Statins (cholesterol)": "Cholesterol medication",
    "Proton pump inhibitors (stomach acid)": "Stomach-acid reducer",
    "Antidepressants": "Antidepressant",
    "Antifungals": "Antifungal",
    "Cardiovascular": "Heart / blood-pressure medication",
    "Immunosuppressants": "Immunosuppressant",
    "Chemotherapy": "Chemotherapy",
    "Chemotherapy / Immunosuppressants": "Chemotherapy / Immunosuppressant",
    "Pain medications": "Pain reliever",
    "Antivirals (HIV)": "HIV medication",
    "Antiepileptics": "Seizure medication",
    "Antipsychotics": "Antipsychotic",
    "Gout medications": "Gout medication",
    "ADHD medications": "ADHD medication",
    "Anesthetics": "Anesthetic",
    "Anti-nausea medications": "Anti-nausea medication",
    "Cystic Fibrosis Modulators": "Cystic-fibrosis medication",
    "Gaucher disease medications": "Gaucher-disease medication",
    "Targeted cancer therapy": "Targeted cancer therapy",
    "Alzheimer medications": "Alzheimer medication",
    "GI motility agents": "Gut-motility medication",
    "Opioid withdrawal agents": "Opioid-withdrawal medication",
    "Narcolepsy medications": "Narcolepsy medication",
    "Dermatology": "Skin medication",
    "Other medications": "Other",
}


def _category_display(category: str) -> str:
    return CATEGORY_DISPLAY.get(category, category)



LIMITATIONS: list[str] = [
    "It is a computational prediction — not a clinical diagnosis. Results "
    "should be confirmed by a clinical lab before being used to prescribe.",
    "It does not predict drug allergies, drug-drug interactions, or how your "
    "kidney/liver function, age, diet, and other medications affect a dose.",
    "It only covers the genes PharmCAT analyzes. Other genes that affect drug "
    "response are not included.",
    "It is not a substitute for genetic counseling. If a result surprises you, "
    "a genetic counselor or pharmacist can help interpret it.",
]

LEARN_MORE: list[dict] = [
    {"label": "CPIC — Clinical Pharmacogenetics Implementation Consortium",
     "url": "https://cpicpgx.org/", "blurb": "Plain-language gene/drug summaries written by the same group that publishes the prescribing guidelines."},
    {"label": "FDA Table of Pharmacogenomic Biomarkers in Drug Labeling",
     "url": "https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling",
     "blurb": "Official list of drugs whose FDA label mentions genetic factors."},
    {"label": "MedlinePlus Genetics",
     "url": "https://medlineplus.gov/genetics/",
     "blurb": "U.S. National Library of Medicine — plain-language gene encyclopedia."},
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_patient_report_v2(
    parsed: ParsedResults,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Render the v2 patient report (HTML + best-effort PDF)."""
    context = _build_context(parsed)
    return html_and_pdf(
        "patient_report_v2.html", context, Path(output_dir), "patient_report_v2"
    )


# ---------------------------------------------------------------------------
# Context builder
# ---------------------------------------------------------------------------

def _build_context(parsed: ParsedResults) -> dict:
    drugs_by_gene = _drugs_by_gene(parsed)
    effective_risk = _effective_gene_risk(parsed)
    overview_cards = _build_overview_cards(parsed, effective_risk, drugs_by_gene)
    counts = _count_by_risk(overview_cards)
    # Contributing-gene counts (a gene that drives drugs in both Action and
    # Review is counted in BOTH dashboard sentences). The bucket counts in
    # `counts` stay mutually exclusive — those drive the at-a-glance grid
    # and deep-dive sections where each gene appears exactly once.
    counts["action_contributing"] = sum(
        1 for drugs in drugs_by_gene.values()
        if any(d["risk_level"] == "action" for d in drugs)
    )
    counts["review_contributing"] = sum(
        1 for drugs in drugs_by_gene.values()
        if any(d["risk_level"] == "review" for d in drugs)
    )
    grouped_by_risk = _group_by_risk(parsed, drugs_by_gene, effective_risk)
    drug_view = _build_drug_view(parsed)
    summary_sentence = _summary_sentence(counts)
    next_steps = _build_next_steps(counts, drug_view["action_drugs"])

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
        "next_steps": next_steps,
        "counts": counts,
        "overview_cards": overview_cards,
        "grouped_by_risk": grouped_by_risk,
        "action_drugs": drug_view["action_drugs"],
        "review_drugs": drug_view["review_drugs"],
        "action_total": len(drug_view["action_drugs"]),
        "review_total": drug_view["review_total"],
        "phenotype_plain": PHENOTYPE_PLAIN,
        "drug_guidance": DRUG_GUIDANCE,
        "category_order": CATEGORY_ORDER,
        "limitations": LIMITATIONS,
        "learn_more": LEARN_MORE,
    }


# ---------------------------------------------------------------------------
# Overview cards (jargon-pruned)
# ---------------------------------------------------------------------------

# Patient-facing phenotype labels used in the at-a-glance grid. Standardised
# so the same CPIC concept always renders the same way (no more "Rapid" vs
# "Rapid Metabolizer" vs "Not tested" vs "No result" inconsistency). Word
# choices: "Slow" instead of "Poor" (avoids the judgmental tone), "Very rapid"
# instead of "Ultrarapid", and "No known risk variants" for the negative MH
# state which was previously rendered as the confusing "Uncertain".
PHENOTYPE_DISPLAY: dict[str, str] = {
    # Metabolizer status
    "Normal Metabolizer": "Normal metabolizer",
    "Extensive Metabolizer": "Normal metabolizer",
    "Intermediate Metabolizer": "Intermediate metabolizer",
    "Likely Intermediate Metabolizer": "Likely intermediate metabolizer",
    "Possible Intermediate Metabolizer": "Possibly intermediate metabolizer",
    "Poor Metabolizer": "Slow metabolizer",
    "Likely Poor Metabolizer": "Likely slow metabolizer",
    "Rapid Metabolizer": "Rapid metabolizer",
    "Ultrarapid Metabolizer": "Very rapid metabolizer",
    # Function status
    "Normal Function": "Normal function",
    "Normal": "Normal function",
    "Increased Function": "Increased function",
    "Decreased Function": "Decreased function",
    "Poor Function": "Low function",
    # Susceptibility (MH genes)
    "Uncertain Susceptibility": "No known risk variants",
    "Malignant Hyperthermia Susceptible": "Risk variant present",
    # CFTR-specific
    "ivacaftor non-responsive in CF patients": "Non-responsive to ivacaftor",
    "ivacaftor responsive in CF patients": "Responsive to ivacaftor",
    # VKORC1 (literal genotype labels — translated to functional meaning)
    "-1639 AA": "Warfarin-sensitive variant",
    "-1639 GA": "Intermediate warfarin sensitivity",
    "-1639 AG": "Intermediate warfarin sensitivity",
    "-1639 GG": "Standard warfarin response",
    # No result / not tested — unified label
    "No Result": "Not tested",
    "n/a": "Not tested",
    "N/A": "Not tested",
}


def _display_phenotype(raw: str | None) -> str:
    """Map a raw CPIC phenotype string to the patient-facing display label.
    Falls back to 'Not tested' for empty/unknown inputs (rather than leaking
    raw clinical strings into the overview card)."""
    if not raw:
        return "Not tested"
    canonical = PHENOTYPE_DISPLAY.get(raw.strip())
    if canonical is not None:
        return canonical
    # Try case-insensitive fallback
    lower = raw.strip().lower()
    for k, v in PHENOTYPE_DISPLAY.items():
        if k.lower() == lower:
            return v
    return "Not tested"


def _build_overview_cards(
    parsed: ParsedResults,
    effective_risk: dict[str, str],
    drugs_by_gene: dict[str, list[dict]],
) -> list[dict]:
    """Build the at-a-glance grid cards.

    `med_count` is the number of Action/Review drugs *driven by this gene's
    non-normal phenotype* (via `_drugs_by_gene`) — i.e. the medications the
    patient should actually pay attention to because of this gene. Not the
    raw size of PharmCAT's `related_drugs` catalogue, which counts every
    drug whose guideline references the gene regardless of the patient's
    own phenotype."""
    cards: list[dict] = []

    def _tissue(gene: str) -> dict:
        key, label, icon = gene_tissue_icon(gene)
        return {"tissue": key, "tissue_label": label, "tissue_icon": icon}

    def _risk(g) -> str:
        return effective_risk.get(g.gene, g.risk_level)

    def _affected_count(gene: str) -> int:
        return len(drugs_by_gene.get(gene, []))

    for g in parsed.definitive_genes:
        risk = _risk(g)
        cards.append({
            "gene": g.gene,
            "risk_level": risk,
            "symbol": ACTION_SYMBOLS[risk],
            "label": ACTION_LABELS[risk],
            "phenotype_short": _display_phenotype(g.phenotype),
            "med_count": _affected_count(g.gene),
            **_tissue(g.gene),
        })

    for g in parsed.ambiguous_genes:
        risk = _risk(g)
        cards.append({
            "gene": g.gene,
            "risk_level": risk,
            "symbol": ACTION_SYMBOLS[risk],
            "label": ACTION_LABELS[risk],
            "phenotype_short": "Inconclusive",
            "med_count": _affected_count(g.gene),
            **_tissue(g.gene),
        })

    for g in parsed.no_call_genes:
        cards.append({
            "gene": g.gene,
            "risk_level": "nodata",
            "symbol": ACTION_SYMBOLS["nodata"],
            "label": ACTION_LABELS["nodata"],
            "phenotype_short": "Not tested",
            "med_count": _affected_count(g.gene),
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


def _summary_sentence(counts: dict[str, int]) -> str:
    if counts["total"] == 0:
        return "No genes were analyzed."
    if counts["action"] == 0 and counts["review"] == 0:
        if counts["nodata"]:
            return (
                f"Good news — none of the tested genes require special attention. "
                f"{counts['nodata']} gene{'s were' if counts['nodata'] != 1 else ' was'} "
                f"not testable from this sample."
            )
        return "Good news — none of the tested genes require special attention."
    pieces: list[str] = []
    if counts["action"]:
        pieces.append(
            f"{counts['action']} gene{'s' if counts['action'] != 1 else ''} need"
            f"{'' if counts['action'] != 1 else 's'} attention before certain prescriptions"
        )
    if counts["review"]:
        pieces.append(
            f"{counts['review']} should be reviewed at your next visit"
        )
    if counts["nodata"]:
        pieces.append(
            f"{counts['nodata']} could not be tested from this sample"
        )
    return _join_clauses(pieces).capitalize() + "."


def _join_clauses(parts: list[str]) -> str:
    if not parts:
        return ""
    if len(parts) == 1:
        return parts[0]
    return ", ".join(parts[:-1]) + " and " + parts[-1]


# ---------------------------------------------------------------------------
# Next-steps block
# ---------------------------------------------------------------------------

def _build_next_steps(
    counts: dict[str, int],
    action_drugs: list[dict],
) -> list[dict]:
    steps: list[dict] = []

    # Always: bring this to the doctor
    steps.append({
        "icon": "🩺",
        "text": "Bring this report to your healthcare provider — especially before starting a new prescription.",
    })

    # Action-bucket drugs to name explicitly (cap at 5 + "and others")
    if action_drugs:
        names: list[str] = []
        seen: set[str] = set()
        for d in action_drugs:
            n = d["drug"]
            if n.lower() in seen:
                continue
            seen.add(n.lower())
            names.append(n)
            if len(names) >= 5:
                break
        suffix = " and others" if len(names) == 5 and len(action_drugs) > 5 else ""
        steps.append({
            "icon": "⚠️",
            "text": f"Pay special attention to: {', '.join(names)}{suffix}. "
                    f"If any of these are prescribed, mention this report first.",
        })

    # No-call genes
    if counts.get("nodata"):
        steps.append({
            "icon": "🔬",
            "text": f"{counts['nodata']} gene{'s were' if counts['nodata'] != 1 else ' was'} "
                    "not testable from this sample. If a doctor needs a result for "
                    "one of these, ask about clinical pharmacogenomic testing.",
        })

    # Always: confirm before prescribing
    steps.append({
        "icon": "✅",
        "text": "Results should be confirmed by a CLIA-certified lab before any "
                "prescribing decisions are made.",
    })

    return steps


# ---------------------------------------------------------------------------
# Functional category grouping (deep-dive, jargon stripped)
# ---------------------------------------------------------------------------

def _effective_gene_risk(parsed: ParsedResults) -> dict[str, str]:
    """Each gene's effective risk = max risk among drugs whose recommendation
    this gene's phenotype actually drives (gene's phenotype in the rec is
    non-normal). If a gene has NO such contributing drugs, fall back to its
    phenotype-derived risk so genes with notable phenotypes but no covered
    drugs still surface.

    Resolves two symmetric attribution bugs:
      - A Normal-phenotype gene (TPMT, ABCG2) getting bumped just because a
        Review/Action drug happens to list it among `affected_genes`.
      - A non-normal-phenotype gene (CYP3A5 PM) staying in Action just
        because the phenotype name sounds severe — when its only covered
        drug parses as Review or Normal, the gene should follow the drug.
    """
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


def _drugs_by_gene(parsed: ParsedResults) -> dict[str, list[dict]]:
    """Index of gene → [non-normal drugs THIS gene's phenotype contributes to].
    A drug is attributed to a gene only if the gene's phenotype in the
    recommendation is non-normal, mirroring the rule used in
    `_effective_gene_risk` so the gene-card affected-meds list and the
    gene's risk bucket are derived from the same evidence."""
    by_gene: dict[str, dict[str, dict]] = defaultdict(dict)
    for d in parsed.drugs:
        if d.risk_level not in ("action", "review"):
            continue
        entry = {
            "drug": d.drug,
            "brand": _brand(d.drug),
            "risk_level": d.risk_level,
            "symbol": ACTION_SYMBOLS[d.risk_level],
        }
        for gene in d.affected_genes:
            if _is_normal_pheno(d.phenotypes.get(gene, "")):
                continue
            existing = by_gene[gene].get(d.drug)
            if (
                existing is None
                or RISK_PRIORITY.get(d.risk_level, 3)
                < RISK_PRIORITY.get(existing["risk_level"], 3)
            ):
                by_gene[gene][d.drug] = entry
    return {
        gene: sorted(
            entries.values(),
            key=lambda e: (RISK_PRIORITY.get(e["risk_level"], 3), e["drug"]),
        )
        for gene, entries in by_gene.items()
    }


def _affected_drugs_by_category(related_drugs: list[str]) -> dict[str, list[str]]:
    """Group a flat list of drug names by therapeutic category. Used for
    definitive-but-no-result genes (HLA-A/B/CYP4F2) so their no-data card
    shows the same 'Medications that depend on X' list as parser-side
    no-calls."""
    buckets: dict[str, list[str]] = defaultdict(list)
    for drug_name in related_drugs or []:
        buckets[drug_category(drug_name)].append(drug_name)
    return {
        cat: sorted(set(drugs))
        for cat, drugs in sorted(buckets.items())
    }


def _group_normal_definitives(normals: list[dict]) -> list[dict]:
    """Merge Normal-risk gene cards that share the same plain-language /
    detail text into a single card. Each merged card lists the genes it
    covers, with the per-gene description and protein type preserved so
    they can be shown inline."""
    groups: "OrderedDict[tuple[str, str], dict]" = OrderedDict()
    for g in normals:
        key = (g["plain_language"], g["detail"])
        if key not in groups:
            groups[key] = {
                "plain_language": g["plain_language"],
                "detail": g["detail"],
                "genes": [],
            }
        groups[key]["genes"].append({
            "gene": g["gene"],
            "description": g["description"],
            "protein_type": g["protein_type"],
            "note": g.get("note"),
        })
    for grp in groups.values():
        grp["genes"].sort(key=lambda x: (x["protein_type"], x["gene"]))
    return list(groups.values())


def _group_by_risk(
    parsed: ParsedResults,
    drugs_by_gene: dict[str, list[dict]],
    effective_risk: dict[str, str],
) -> dict:
    """Group enriched gene cards by risk bucket (action / review / normal /
    no_call) instead of by functional category. Sorting within each bucket
    is (protein_type, gene) so similar protein types cluster naturally,
    with no sub-headings between them."""
    action: list[dict] = []
    review: list[dict] = []
    normal_raw: list[dict] = []
    extra_no_data: list[dict] = []

    for g in parsed.definitive_genes:
        risk = effective_risk.get(g.gene, g.risk_level)
        enriched = _enrich_definitive(g, drugs_by_gene.get(g.gene, []), risk)
        if risk == "action":
            action.append(enriched)
        elif risk == "review":
            review.append(enriched)
        elif risk == "normal":
            normal_raw.append(enriched)
        else:  # 'nodata' (e.g. HLA-A/B/CYP4F2 with No Result phenotype)
            extra_no_data.append({
                "gene": g.gene,
                "risk_level": "nodata",
                "symbol": ACTION_SYMBOLS["nodata"],
                "label": ACTION_LABELS["nodata"],
                "description": g.description,
                "protein_type": g.protein_type,
                "affected_drugs": _affected_drugs_by_category(g.related_drugs),
                "caveat": g.caveat,
            })

    for g in parsed.ambiguous_genes:
        risk = effective_risk.get(g.gene, g.risk_level)
        enriched = _enrich_ambiguous(g, drugs_by_gene.get(g.gene, []), risk)
        if risk == "action":
            action.append(enriched)
        elif risk == "review":
            review.append(enriched)
        elif risk == "normal":
            normal_raw.append(enriched)

    no_call: list[dict] = [_enrich_no_call(g) for g in parsed.no_call_genes]
    no_call.extend(extra_no_data)

    def _sort_key(d: dict) -> tuple[str, str]:
        return (d.get("protein_type") or "", d["gene"])

    action.sort(key=_sort_key)
    review.sort(key=_sort_key)
    normal_raw.sort(key=_sort_key)
    no_call.sort(key=_sort_key)

    return {
        "action": action,
        "review": review,
        "normal_groups": _group_normal_definitives(normal_raw),
        "no_call": no_call,
    }


def _enrich_definitive(
    g: DefinitiveGene, affected_drugs: list[dict], risk_level: str,
) -> dict:
    explanations = _explain_phenotype(g.phenotype)
    return {
        "gene": g.gene,
        "risk_level": risk_level,
        "symbol": ACTION_SYMBOLS[risk_level],
        "label": ACTION_LABELS[risk_level],
        "description": g.description,
        "protein_type": g.protein_type,
        "plain_language": explanations["brief"],
        "detail": explanations["detail"],
        "caveat": g.caveat,
        "affected_drugs": affected_drugs,
        "is_ambiguous": False,
        "note": _gene_note(g.gene, g.phenotype),
    }


def _gene_note(gene: str, phenotype: str | None) -> str | None:
    """Per-gene contextual notes for cases where the patient's genetic test
    result has clinical relevance beyond what the drug list captures. Used
    today for CACNA1S/RYR1 when 'Uncertain Susceptibility' — the standard
    'negative' result that an anesthesiologist may still want to know about
    before any procedure."""
    pheno = (phenotype or "").lower().strip()
    if gene in ("CACNA1S", "RYR1") and "uncertain susceptibility" in pheno:
        return (
            "No known malignant hyperthermia variants were found in this gene. "
            "Worth mentioning to an anesthesiologist before any procedure — "
            "the test cannot rule out every rare variant, but no catalogued "
            "susceptibility variant is present."
        )
    return None


def _enrich_ambiguous(
    g: AmbiguousGene, affected_drugs: list[dict], risk_level: str,
) -> dict:
    return {
        "gene": g.gene,
        "risk_level": risk_level,
        "symbol": ACTION_SYMBOLS[risk_level],
        "label": ACTION_LABELS[risk_level],
        "description": g.description,
        "protein_type": g.protein_type,
        "plain_language": None,
        "detail": None,
        "caveat": g.caveat,
        "affected_drugs": affected_drugs,
        "is_ambiguous": True,
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


# ---------------------------------------------------------------------------
# Drug grouping for the patient view (with brand names + consequences)
# ---------------------------------------------------------------------------

def _brand(drug_name: str) -> str:
    return DRUG_BRAND_NAMES.get(drug_name.strip().lower(), "")


def _action_consequence(therapeutic_category: str) -> tuple[str, str]:
    return ACTION_CONSEQUENCE_BY_CATEGORY.get(
        therapeutic_category, ACTION_CONSEQUENCE_FALLBACK
    )


def _build_drug_view(parsed: ParsedResults) -> dict:
    """Patient-facing drug view as two flat tables (Action, Review).

    Each entry is deduplicated by drug name (multiple source/gene rows
    for the same drug merge into one — the gene list becomes the union
    and the risk level is the most concerning). Normal-risk drugs are
    hidden.

    Each entry carries a pre-rendered `rationale` HTML sentence that
    names the patient's gene(s) inline — so the column reads on its own
    without bouncing between sections."""
    by_drug: dict[str, dict] = {}

    for d in parsed.drugs:
        existing = by_drug.get(d.drug)
        if existing is None:
            by_drug[d.drug] = {
                "drug": d.drug,
                "brand": _brand(d.drug),
                "phenotypes": dict(d.phenotypes or {}),
                "risk_level": d.risk_level,
                "therapeutic_category": d.therapeutic_category,
            }
            continue
        # Merge per-gene phenotypes — prefer a non-normal one if any exists.
        for gene, pheno in (d.phenotypes or {}).items():
            current = existing["phenotypes"].get(gene)
            if current is None or _is_normal_pheno(current):
                existing["phenotypes"][gene] = pheno
        if (
            RISK_PRIORITY.get(d.risk_level, 3)
            < RISK_PRIORITY.get(existing["risk_level"], 3)
        ):
            existing["risk_level"] = d.risk_level
            existing["therapeutic_category"] = d.therapeutic_category

    action_drugs: list[dict] = []
    review_drugs: list[dict] = []

    for entry in by_drug.values():
        if entry["risk_level"] not in ("action", "review"):
            continue
        contributing = sorted([
            gene for gene, p in entry["phenotypes"].items()
            if not _is_normal_pheno(p)
        ])
        if not contributing:
            continue

        gene_text = ", ".join(contributing)

        if entry["risk_level"] == "action":
            plural = "s" if len(contributing) > 1 else ""
            genes_phrase = f"Your <strong>{gene_text}</strong> result{plural}"
            directive, reason_template = _action_consequence(
                entry["therapeutic_category"]
            )
            reason = reason_template.format(genes=genes_phrase)
            rationale = (
                f'<span class="rationale-directive">{directive}</span>'
                f"{reason}"
            )
            action_drugs.append({
                "drug": entry["drug"],
                "brand": entry["brand"],
                "gene": gene_text,
                "rationale": rationale,
                "therapeutic_category": entry["therapeutic_category"],
                "type_display": _category_display(entry["therapeutic_category"]),
            })
        else:
            review_drugs.append({
                "drug": entry["drug"],
                "brand": entry["brand"],
                "gene": gene_text,
                "therapeutic_category": entry["therapeutic_category"],
                "type_display": _category_display(entry["therapeutic_category"]),
            })

    action_drugs.sort(key=lambda d: (d["therapeutic_category"], d["drug"]))
    review_drugs.sort(key=lambda d: d["drug"])

    return {
        "action_drugs": action_drugs,
        "review_drugs": review_drugs,
        "review_total": len(review_drugs),
    }
