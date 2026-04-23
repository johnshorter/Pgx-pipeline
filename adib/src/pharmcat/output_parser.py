"""
PharmCAT output parser (v3 compatible).

Reads PharmCAT's phenotype JSON and report HTML to extract:
- Per-gene results (star alleles, diplotype, phenotype, activity score)
- Per-drug recommendations (CPIC/DPWG guidance, evidence strength, dosing)
- Flags CYP2D6 as "no data" caveat.
"""

import json
import logging
import re
from pathlib import Path
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# ── Gene categories and descriptions ──────────────────────────────────────
GENE_CATEGORIES: dict[str, dict] = {
    "Phase I Metabolism (CYP Enzymes)": {
        "genes": {"CYP2B6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2"},
        "description": (
            "Cytochrome P450 enzymes are the body's primary drug-metabolizing system. "
            "They break down many common medications in the liver."
        ),
    },
    "Phase II Metabolism": {
        "genes": {"DPYD", "NAT2", "TPMT", "NUDT15", "UGT1A1"},
        "description": (
            "These enzymes modify drugs and their breakdown products to make them "
            "easier to eliminate. Variations can affect how your body handles "
            "certain chemotherapy agents and other medications."
        ),
    },
    "Drug Transporters": {
        "genes": {"ABCG2", "SLCO1B1"},
        "description": (
            "Transporter proteins move drugs into and out of cells. Variations can "
            "alter drug absorption and distribution, affecting how much medication "
            "reaches its target."
        ),
    },
    "Immune Markers (HLA)": {
        "genes": {"HLA-A", "HLA-B"},
        "description": (
            "HLA genes help the immune system distinguish the body's own cells from "
            "foreign substances. Certain HLA variants are associated with severe drug "
            "hypersensitivity reactions."
        ),
    },
    "Other Pharmacogenes": {
        "genes": {
            "CACNA1S", "CFTR", "F2", "F5", "G6PD", "IFNL3",
            "MT-RNR1", "RYR1", "VKORC1",
        },
        "description": (
            "Additional genes that influence drug response through various mechanisms "
            "including drug targets, enzyme deficiencies, and receptor sensitivity."
        ),
    },
}


def gene_category(gene_symbol: str) -> tuple[str, str]:
    """Return (category_name, category_description) for a gene."""
    for cat_name, cat_info in GENE_CATEGORIES.items():
        if gene_symbol in cat_info["genes"]:
            return cat_name, cat_info["description"]
    return "Other Pharmacogenes", GENE_CATEGORIES["Other Pharmacogenes"]["description"]


# ── Action levels (colorblind-safe) ──────────────────────────────────────
# "action"  = needs clinical attention (was red)   — symbol: ▲
# "review"  = monitor / discuss          (was yellow) — symbol: ◆
# "normal"  = standard / no action       (was green)  — symbol: ✓
# "nodata"  = could not be determined                  — symbol: —
ACTION_SYMBOLS = {
    "action": "\u25B2",   # ▲
    "review": "\u25C6",   # ◆
    "normal": "\u2713",   # ✓
    "nodata": "\u2014",   # —
}

# CYP2D6 is not reliably called from short-read sequencing
CYP2D6_CAVEAT = (
    "CYP2D6 results may be unreliable from short-read whole genome sequencing "
    "due to the gene's complex structural variation (deletions, duplications, "
    "hybrid alleles). Clinical CYP2D6 testing is recommended for actionable decisions."
)

# Risk level mapping for reports
# Uses colorblind-safe action levels: action / review / normal / nodata
PHENOTYPE_RISK = {
    "poor metabolizer": "action",
    "intermediate metabolizer": "review",
    "normal metabolizer": "normal",
    "rapid metabolizer": "review",
    "ultrarapid metabolizer": "action",
    "likely poor metabolizer": "action",
    "likely intermediate metabolizer": "review",
    "possible intermediate metabolizer": "review",
    "increased function": "review",
    "decreased function": "review",
    "normal function": "normal",
    "normal": "normal",
    "indeterminate": "review",
    "uncertain susceptibility": "normal",
    "no result": "nodata",
    "n/a": "normal",
}


@dataclass
class GeneResult:
    gene: str
    diplotype: str
    phenotype: str
    activity_score: str | None
    star_alleles: list[str]
    has_data: bool
    risk_level: str  # "action", "review", "normal", "nodata"
    caveat: str | None = None
    category: str = ""          # e.g. "Phase I Metabolism (CYP Enzymes)"
    category_desc: str = ""     # plain-language description of the category


@dataclass
class DrugRecommendation:
    drug: str
    gene: str
    guideline_source: str  # e.g. "CPIC", "DPWG"
    recommendation: str
    classification: str  # e.g. "Strong", "Moderate", "No Action"
    implications: str
    phenotype: str
    risk_level: str  # "green", "yellow", "red"
    url: str | None = None


@dataclass
class ParsedResults:
    gene_results: list[GeneResult] = field(default_factory=list)
    drug_recommendations: list[DrugRecommendation] = field(default_factory=list)
    sample_id: str | None = None
    missing_genes: list[str] = field(default_factory=list)


def parse_pharmcat_output(
    phenotyper_json_path: str | Path | None = None,
    reporter_json_path: str | Path | None = None,
    reporter_html_path: str | Path | None = None,
) -> ParsedResults:
    """Parse PharmCAT output files into structured results.

    Args:
        phenotyper_json_path: Path to the .phenotype.json file (gene calls).
        reporter_json_path: Path to the .report.json file (unused in v3).
        reporter_html_path: Path to the .report.html file (drug recommendations).

    Returns:
        ParsedResults with gene results and drug recommendations.
    """
    results = ParsedResults()

    # Parse phenotype JSON (gene-level results)
    if phenotyper_json_path:
        phenotyper_json_path = Path(phenotyper_json_path)
        if phenotyper_json_path.exists():
            _parse_phenotype_json(phenotyper_json_path, results)
        else:
            logger.warning("Phenotype JSON not found: %s", phenotyper_json_path)

    # Parse reporter HTML (drug recommendations)
    if reporter_html_path:
        reporter_html_path = Path(reporter_html_path)
        if reporter_html_path.exists():
            _parse_reporter_html(reporter_html_path, results)
        else:
            logger.warning("Reporter HTML not found: %s", reporter_html_path)

    # Deduplicate drug recommendations: keep one per drug+gene+source,
    # preferring higher risk level and non-empty recommendations
    results.drug_recommendations = _deduplicate_drugs(results.drug_recommendations)

    # Ensure CYP2D6 caveat
    for gene_result in results.gene_results:
        if gene_result.gene == "CYP2D6":
            gene_result.caveat = CYP2D6_CAVEAT

    return results


def _parse_phenotype_json(path: Path, results: ParsedResults) -> None:
    """Parse the PharmCAT v3 phenotype JSON for gene-level results."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        logger.error("Failed to parse phenotype JSON: %s", e)
        return

    # Extract sample ID from metadata
    metadata = data.get("matcherMetadata", {})
    input_file = metadata.get("inputFilename", "")
    if input_file:
        # Use filename stem as sample ID
        results.sample_id = Path(input_file).stem

    # PharmCAT v3: geneReports is a dict keyed by gene symbol
    gene_reports = data.get("geneReports", {})

    if isinstance(gene_reports, list):
        # Fallback for older formats
        for entry in gene_reports:
            gene = entry.get("gene", entry.get("geneSymbol", ""))
            if gene:
                gene_reports_dict = {gene: entry}
        gene_reports = gene_reports_dict if gene_reports else {}

    for gene_symbol, gene_data in gene_reports.items():
        # Get best diplotype from recommendationDiplotypes (highest matchScore)
        rec_diplotypes = gene_data.get("recommendationDiplotypes", [])

        cat_name, cat_desc = gene_category(gene_symbol)

        if not rec_diplotypes:
            # Gene was in panel but no calls made
            results.gene_results.append(
                GeneResult(
                    gene=gene_symbol,
                    diplotype="Not determined",
                    phenotype="No result",
                    activity_score=None,
                    star_alleles=[],
                    has_data=False,
                    risk_level="nodata",
                    category=cat_name,
                    category_desc=cat_desc,
                )
            )
            continue

        # Get the top-scoring diplotype
        top = max(rec_diplotypes, key=lambda x: x.get("matchScore", 0))

        # Extract alleles (may be None for some genes)
        allele1_obj = top.get("allele1")
        allele2_obj = top.get("allele2")
        a1_name = allele1_obj.get("name", "Unknown") if allele1_obj else "Unknown"
        a2_name = allele2_obj.get("name", "Unknown") if allele2_obj else "Unknown"

        diplotype = f"{a1_name}/{a2_name}"
        star_alleles = [a for a in [a1_name, a2_name] if a != "Unknown"]

        # Extract phenotype
        phenotypes = top.get("phenotypes", [])
        phenotype = phenotypes[0] if phenotypes else "No result"

        # Activity score
        activity_score = top.get("activityScore")
        if activity_score is not None:
            activity_score = str(activity_score)

        # Determine if there's meaningful data
        has_data = a1_name != "Unknown" or a2_name != "Unknown"

        # Determine risk level from phenotype
        risk_level = _phenotype_to_risk(phenotype)

        results.gene_results.append(
            GeneResult(
                gene=gene_symbol,
                diplotype=diplotype if has_data else "Not determined",
                phenotype=phenotype,
                activity_score=activity_score,
                star_alleles=star_alleles,
                has_data=has_data,
                risk_level=risk_level,
                category=cat_name,
                category_desc=cat_desc,
            )
        )

    # Identify genes with no data
    results.missing_genes = [
        g.gene for g in results.gene_results if not g.has_data
    ]


def _parse_reporter_html(path: Path, results: ParsedResults) -> None:
    """Parse the PharmCAT v3 reporter HTML for drug recommendations.

    PharmCAT v3 HTML structure: each drug is in a <section class="guideline drugname">
    containing a <table> with <tbody>. Rows use <tr class="top-aligned ..."> but may
    not have closing </tr> tags, so we split by <tr markers instead.
    """
    try:
        content = path.read_text(encoding="utf-8")
    except Exception as e:
        logger.error("Failed to read reporter HTML: %s", e)
        return

    # Find all drug h3 IDs
    drug_ids = re.findall(r'<h3[^>]*id="([^"]+)"[^>]*>[^<]*</h3>', content)

    for drug_name in drug_ids:
        h3_pos = content.find(f'id="{drug_name}"')
        if h3_pos == -1:
            continue

        # Find <tbody> after this drug's h3
        tbody_pos = content.find("<tbody>", h3_pos)
        if tbody_pos == -1:
            continue

        tbody_end = content.find("</tbody>", tbody_pos)
        if tbody_end == -1:
            continue

        tbody_html = content[tbody_pos:tbody_end]

        # Split tbody into row chunks (since </tr> may be missing)
        tr_chunks = re.split(r"(?=<tr\s)", tbody_html)
        tr_chunks = [c for c in tr_chunks if c.startswith("<tr")]

        for chunk in tr_chunks:
            # Extract 5 cells: Source | Genes | Implications | Recommendation | Classification
            cells = re.findall(r"<td[^>]*>(.*?)(?:</td>|$)", chunk, re.DOTALL)
            if len(cells) < 4:
                continue

            # Cell 0: Source (guideline name, URL, classification tag)
            source_cell = cells[0]
            source_match = re.search(
                r'<a\s+href="([^"]+)"[^>]*>([^<]+)</a>', source_cell
            )
            url = source_match.group(1) if source_match else None
            source_text = source_match.group(2) if source_match else ""

            if "cpic" in source_text.lower():
                guideline_source = "CPIC"
            elif "dpwg" in source_text.lower():
                guideline_source = "DPWG"
            elif "fda" in source_text.lower():
                guideline_source = "FDA"
            else:
                guideline_source = _strip_html(source_text) or "Unknown"

            # Classification tag
            class_match = re.search(
                r'<div\s+class="tag\s+[^"]*"[^>]*>([^<]+)</div>', source_cell
            )
            classification = class_match.group(1).strip() if class_match else ""

            # Cell 1: Genes (genotype + phenotype)
            gene_cell = cells[1]
            gene_match = re.search(
                r'<a\s+href="#([^"]+)">([^<]+)</a>', gene_cell
            )
            gene = gene_match.group(2) if gene_match else ""

            # Phenotype from cell 1
            pheno_match = re.search(
                r'Phenotype.*?</div>\s*<p[^>]*>(.*?)</p>', gene_cell, re.DOTALL
            )
            phenotype_text = _strip_html(pheno_match.group(1)) if pheno_match else ""

            # Cell 2: Implications
            implications = _strip_html(cells[2])

            # Cell 3: Recommendation
            recommendation = _strip_html(cells[3])

            # Determine risk level
            risk_level = _classification_to_risk(classification, recommendation)

            if drug_name and (recommendation or classification):
                results.drug_recommendations.append(
                    DrugRecommendation(
                        drug=drug_name,
                        gene=gene,
                        guideline_source=guideline_source,
                        recommendation=recommendation,
                        classification=classification,
                        implications=implications,
                        phenotype=phenotype_text,
                        risk_level=risk_level,
                        url=url,
                    )
                )


def _strip_html(text: str) -> str:
    """Remove HTML tags and normalize whitespace."""
    text = re.sub(r"<[^>]+>", " ", text)
    text = text.replace("&nbsp;", " ").replace("&amp;", "&")
    text = text.replace("&quot;", '"').replace("&lt;", "<").replace("&gt;", ">")
    return re.sub(r"\s+", " ", text).strip()


def _deduplicate_drugs(drugs: list[DrugRecommendation]) -> list[DrugRecommendation]:
    """Keep one recommendation per drug+gene+source, preferring higher risk."""
    risk_priority = {"red": 0, "yellow": 1, "green": 2}
    best: dict[tuple[str, str, str], DrugRecommendation] = {}

    for d in drugs:
        key = (d.drug, d.gene, d.guideline_source)
        if key not in best:
            best[key] = d
        else:
            existing = best[key]
            # Prefer higher risk, or longer recommendation text
            if risk_priority.get(d.risk_level, 2) < risk_priority.get(existing.risk_level, 2):
                best[key] = d
            elif d.risk_level == existing.risk_level and len(d.recommendation) > len(existing.recommendation):
                best[key] = d

    return list(best.values())


def _phenotype_to_risk(phenotype: str) -> str:
    """Map a phenotype string to an action level."""
    if not phenotype:
        return "nodata"
    lower = phenotype.lower().strip()
    for key, risk in PHENOTYPE_RISK.items():
        if key in lower:
            return risk
    if "normal" in lower:
        return "normal"
    return "review"


def _classification_to_risk(classification: str, recommendation: str) -> str:
    """Map PharmCAT classification tag + recommendation text to action level."""
    cls_lower = classification.lower().strip()
    rec_lower = recommendation.lower() if recommendation else ""

    if cls_lower in ("no action", "no action needed"):
        return "normal"

    if cls_lower in ("strong", "actionable"):
        return "action"

    if cls_lower in ("moderate",):
        return "review"

    if cls_lower in ("informative",):
        return "review"

    action_keywords = [
        "avoid", "contraindicated", "do not use", "not recommended",
        "consider alternative", "use with caution",
    ]
    if any(kw in rec_lower for kw in action_keywords):
        return "action"

    review_keywords = [
        "caution", "monitor", "consider", "may need", "dose adjustment",
        "reduced", "increased", "titrate",
    ]
    if any(kw in rec_lower for kw in review_keywords):
        return "review"

    return "normal"
