#!/usr/bin/env python3
"""
PGx Reporter — Standalone pharmacogenomic report generator.

Given a VCF file as input, this tool produces two reports:

  1. A patient-friendly report (plain language, colorblind-safe design).
  2. A clinician report (full technical detail, PharmCAT-backed).

Both are produced as HTML and (when WeasyPrint is installed) PDF.

Pipeline
--------
  1. Validate the input VCF.
  2. Discover PGx positions (via an initial PharmCAT pass, unless a reference
     phenotype.json is supplied via --reference-phenotype).
  3. Screen VCF coverage + preprocess VCF (fill reference calls at missing
     PGx positions — required for GIAB-style benchmark VCFs).
  4. Final PharmCAT run on the preprocessed VCF (with -research cyp2d6
     by default to also call CYP2D6).
  5. Parse PharmCAT JSON output into a structured model.
  6. Render patient + clinician reports (HTML + optional PDF).

Requirements
------------
  * Python 3.10+
  * Java 17+ on PATH (required by PharmCAT)
  * PharmCAT JAR (download with: python download_pharmcat.py)
  * pip install -r requirements.txt   # jinja2, weasyprint (optional for PDFs)

Typical usage
-------------
    # One-off — runs the two-pass pipeline end to end.
    python pgx_report.py sample.vcf --output-dir ./output

    # Skip PDF, output HTML only.
    python pgx_report.py sample.vcf --no-pdf

    # Batch-process multiple samples sharing a phenotype.json cache
    # (much faster — avoids redoing position discovery per sample).
    python pgx_report.py HG001.vcf -o out/HG001
    python pgx_report.py HG002.vcf -o out/HG002 \\
        --reference-phenotype out/HG001/pharmcat_raw/HG001.phenotype.json
"""

import argparse
import logging
import sys
from pathlib import Path

# Make the bundled `src` package importable regardless of how the script is invoked.
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.pharmcat.vcf_validator import validate_vcf
from src.pharmcat.runner import run_pharmcat, find_pharmcat_jar
from src.pharmcat.output_parser import parse_pharmcat_output
from src.reports.patient_report import generate_patient_report
from src.reports.clinician_report import generate_clinician_report
from src.screening.screen_pharmacogenes import screen_vcf, print_screening_report
from src.preprocessing.preprocess_vcf import preprocess_vcf, clean_vcf


# ──────────────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────────────

def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pgx_report",
        description=(
            "Run the standalone PGx reporting pipeline on a VCF file. "
            "Outputs patient + clinician reports (HTML, optionally PDF)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python pgx_report.py sample.vcf\n"
            "  python pgx_report.py sample.vcf -o ./results\n"
            "  python pgx_report.py sample.vcf --no-pdf\n"
            "  python pgx_report.py sample.vcf --skip-preprocess\n"
            "  python pgx_report.py sample.vcf --reference-phenotype cached.json\n"
        ),
    )
    p.add_argument("vcf", type=Path, help="Input VCF file (GRCh38).")
    p.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=PROJECT_ROOT / "output",
        help="Directory to write results into (default: ./output).",
    )
    p.add_argument(
        "--jar",
        type=Path,
        default=None,
        help="Path to the PharmCAT JAR (default: auto-detect in lib/).",
    )
    p.add_argument(
        "--no-pdf",
        action="store_true",
        help="Do not generate PDF reports (HTML only).",
    )
    p.add_argument(
        "--skip-preprocess",
        action="store_true",
        help=(
            "Skip the reference-fill preprocessing step. "
            "Use this only if your VCF already contains genotype calls at all PGx positions."
        ),
    )
    p.add_argument(
        "--reference-phenotype",
        type=Path,
        default=None,
        help=(
            "Path to an existing PharmCAT phenotype.json to use for PGx position discovery. "
            "When provided, the initial PharmCAT pass is skipped. "
            "Useful for batch-processing multiple samples."
        ),
    )
    p.add_argument(
        "--research",
        type=str,
        default="cyp2d6",
        help=(
            "Comma-separated PharmCAT research-mode flags "
            "(default: cyp2d6 — required to call CYP2D6). "
            "Pass an empty string to disable."
        ),
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="PharmCAT execution timeout in seconds (default: 300).",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )
    return p


# ──────────────────────────────────────────────────────────────────────────────
# Pipeline
# ──────────────────────────────────────────────────────────────────────────────

def main() -> int:
    args = build_parser().parse_args()
    setup_logging(args.verbose)
    log = logging.getLogger("pgx-reporter")

    if not args.vcf.exists():
        log.error("Input VCF not found: %s", args.vcf)
        return 1

    research_mode = (
        [r.strip() for r in args.research.split(",") if r.strip()]
        if args.research else None
    )

    # ── Step 1: Validate ─────────────────────────────────────────────────────
    log.info("Step 1/6: Validating VCF: %s", args.vcf)
    validation = validate_vcf(args.vcf)
    for w in validation.warnings:
        log.warning("  VCF warning: %s", w)
    if not validation.valid:
        for e in validation.errors:
            log.error("  VCF error: %s", e)
        log.error("VCF validation failed. Aborting.")
        return 1
    log.info("  VCF validation passed.")

    # Locate PharmCAT JAR
    jar_path = args.jar or find_pharmcat_jar(PROJECT_ROOT / "lib")
    if jar_path is None:
        log.error(
            "PharmCAT JAR not found. Download with:\n"
            "    python download_pharmcat.py\n"
            "…or pass --jar /path/to/pharmcat.jar."
        )
        return 1

    pharmcat_raw_dir = args.output_dir / "pharmcat_raw"
    final_vcf = args.vcf

    if not args.skip_preprocess:
        # ── Step 2: Position discovery ──────────────────────────────────────
        phenotype_json_path = args.reference_phenotype
        if phenotype_json_path:
            log.info("Step 2/6: Using reference phenotype.json: %s", phenotype_json_path)
        else:
            log.info("Step 2/6: Initial PharmCAT pass (discovering PGx positions)...")
            cleaned_vcf = args.output_dir / "cleaned" / (args.vcf.stem + "_cleaned.vcf")
            log.info("  Cleaning VCF (fixing FORMAT/sample field mismatches)...")
            clean_vcf(args.vcf, cleaned_vcf)

            initial = run_pharmcat(
                vcf_path=cleaned_vcf,
                output_dir=pharmcat_raw_dir,
                jar_path=jar_path,
                timeout=args.timeout,
            )
            if not initial.success:
                log.error("Initial PharmCAT run failed: %s", initial.error_message)
                return 1
            if not initial.phenotyper_json:
                log.error("Initial run produced no phenotype.json — cannot determine missing positions.")
                return 1
            phenotype_json_path = initial.phenotyper_json
            log.info("  Got %s", phenotype_json_path.name)

        # ── Step 3: Screen + preprocess ─────────────────────────────────────
        log.info("Step 3/6: Screening PGx coverage and preprocessing VCF...")
        screening_report = screen_vcf(args.vcf, phenotype_json_path)
        print_screening_report(screening_report)

        preprocessed_vcf = args.output_dir / "preprocessed" / "preprocessed.vcf"
        pre = preprocess_vcf(
            input_vcf=args.vcf,
            phenotype_json=phenotype_json_path,
            output_vcf=preprocessed_vcf,
        )
        log.info(
            "  Preprocessor: added %d reference-fill records across %d genes",
            pre.added_records, len(pre.genes_filled),
        )
        for gene, count in sorted(pre.genes_filled.items()):
            log.info("    %s: %d positions filled", gene, count)

        final_vcf = preprocessed_vcf

    # ── Step 4: Final PharmCAT run ──────────────────────────────────────────
    step = "Step 4/6" if not args.skip_preprocess else "Step 2/4"
    log.info("%s: Running PharmCAT%s (JAR: %s)...",
             step, " with research mode" if research_mode else "", jar_path.name)
    final_out = args.output_dir / "pharmcat_final" if not args.skip_preprocess else pharmcat_raw_dir
    result = run_pharmcat(
        vcf_path=final_vcf,
        output_dir=final_out,
        jar_path=jar_path,
        timeout=args.timeout,
        research_mode=research_mode,
    )
    if not result.success:
        log.error("PharmCAT failed: %s", result.error_message)
        if result.stderr:
            log.error("stderr: %s", result.stderr[:1000])
        return 1
    log.info("  PharmCAT completed.")

    # ── Step 5: Parse ───────────────────────────────────────────────────────
    step = "Step 5/6" if not args.skip_preprocess else "Step 3/4"
    log.info("%s: Parsing PharmCAT output...", step)
    parsed = parse_pharmcat_output(
        phenotyper_json_path=result.phenotyper_json,
        reporter_json_path=result.reporter_json,
        reporter_html_path=result.reporter_html,
    )
    log.info("  %d gene results, %d drug recommendations.",
             len(parsed.gene_results), len(parsed.drug_recommendations))
    if parsed.missing_genes:
        log.info("  Genes not reported: %s", ", ".join(parsed.missing_genes))

    # ── Step 6: Reports ─────────────────────────────────────────────────────
    step = "Step 6/6" if not args.skip_preprocess else "Step 4/4"
    log.info("%s: Generating reports...", step)
    report_dir = args.output_dir / "reports"

    want_pdf = not args.no_pdf
    patient_files = generate_patient_report(parsed, report_dir)
    clinician_files = generate_clinician_report(parsed, report_dir)

    log.info("─" * 60)
    log.info("Reports written to: %s", report_dir)
    log.info("  patient_report.html       → %s", patient_files["html"])
    if want_pdf and "pdf" in patient_files:
        log.info("  patient_report.pdf        → %s", patient_files["pdf"])
    log.info("  clinician_report.html     → %s", clinician_files["html"])
    if want_pdf and "pdf" in clinician_files:
        log.info("  clinician_report.pdf      → %s", clinician_files["pdf"])
    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
