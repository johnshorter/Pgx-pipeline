#!/usr/bin/env python3
"""
PGx Reporter — unified pharmacogenomics CLI.

Takes a VCF, runs PharmCAT, and produces both a patient-friendly and a
clinician/researcher report (HTML + best-effort PDF).

Pipeline stages (configurable via flags):

  1. Validate VCF.
  2. Filter (optional, --filter)            — Marco's WGS fast-path.
  3. Position discovery (optional)          — initial PharmCAT pass for ref-fill.
  4. Reference-fill (optional, --ref-fill)  — Adib's GIAB rescue.
  5. Final PharmCAT run.
  6. Parse output → unified 3-bucket / 4-level model.
  7. Render patient + clinician reports.

Examples:
    # GIAB benchmark VCF (HG001/HG002/HG005): use ref-fill
    python src/pgx_report.py HG005.vcf --ref-fill -o out/

    # Whole-genome VCF: filter first to PharmCAT positions
    python src/pgx_report.py wgs.vcf.gz --filter --java-memory 4g -o out/

    # Apply both (filter first, then fill the small filtered VCF)
    python src/pgx_report.py wgs.vcf.gz --filter --ref-fill -o out/

    # Batch: reuse phenotype.json across samples (skips position discovery)
    python src/pgx_report.py HG002.vcf --ref-fill \\
        --reference-phenotype out/HG001/pharmcat_raw/HG001.phenotype.json
"""

import argparse
import logging
import os
import shutil
import sys
import time
from pathlib import Path

# Make the unified `src/` package importable as flat modules.
SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from config.settings import (
    DEFAULT_OUTPUT_DIR, JAVA_EXECUTABLE, LIB_DIR,
    PHARMCAT_RESEARCH_DEFAULT, PHARMCAT_TIMEOUT_DEFAULT,
)
from pharmcat.runner import find_pharmcat_jar, run_pharmcat
from pharmcat.vcf_validator import validate_vcf
from pharmcat.output_parser import parse_pharmcat_output
from preprocessing.clean_vcf import clean_vcf
from preprocessing.filter_vcf import FilterError, filter_vcf
from preprocessing.reference_fill import reference_fill_vcf
from reports.clinician_report import generate_clinician_report
from reports.patient_report import generate_patient_report
from screening.screen_pharmacogenes import print_screening_report, screen_vcf


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def derive_sample_id(vcf_path: Path) -> str:
    name = vcf_path.name
    for suffix in (".vcf.gz", ".vcf"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return vcf_path.stem


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pgx_report",
        description=(
            "Run the unified PGx reporting pipeline on a VCF file. Outputs "
            "patient + clinician reports (HTML + best-effort PDF)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("vcf", type=Path, help="Input VCF file (GRCh38).")
    p.add_argument(
        "-o", "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help="Output root directory (default: ./output).",
    )
    p.add_argument(
        "--sample-id", default=None,
        help="Sample ID (used for the per-sample output subfolder). "
             "Defaults to the VCF filename without extension.",
    )
    p.add_argument(
        "--jar", type=Path, default=None,
        help="Path to the PharmCAT JAR (default: auto-detect in lib/).",
    )
    p.add_argument(
        "--filter", action="store_true",
        help="Pre-filter the VCF to PharmCAT-relevant positions (Marco's WGS "
             "fast-path). Strongly recommended for whole-genome VCFs.",
    )
    p.add_argument(
        "--ref-fill", action="store_true",
        help="Reference-fill missing PGx positions (Adib's GIAB rescue). Enables "
             "the two-pass discovery flow unless --reference-phenotype is given.",
    )
    p.add_argument(
        "--skip-validation", action="store_true",
        help="Skip the VCF validation step.",
    )
    p.add_argument(
        "--reference-phenotype", type=Path, default=None,
        help="Existing phenotype.json to use for PGx position discovery. When "
             "provided with --ref-fill, the initial PharmCAT pass is skipped.",
    )
    p.add_argument(
        "--research", default=",".join(PHARMCAT_RESEARCH_DEFAULT),
        help="Comma-separated PharmCAT research-mode flags (default: none). "
             "Pass 'cyp2d6' for best-effort CYP2D6 calls — but PharmCAT 3.2.0+ "
             "disables the full reporter in research mode, which means you lose "
             "all CPIC/DPWG drug recommendations.",
    )
    p.add_argument(
        "--java-memory", default=None,
        help="Java -Xmx heap size (e.g. '4g'). Recommended for large VCFs "
             "without --filter.",
    )
    p.add_argument(
        "--java-executable", default=JAVA_EXECUTABLE,
        help="Java binary to invoke (default: auto-detected).",
    )
    p.add_argument(
        "--timeout", type=int, default=PHARMCAT_TIMEOUT_DEFAULT,
        help=f"PharmCAT subprocess timeout in seconds (default: {PHARMCAT_TIMEOUT_DEFAULT}).",
    )
    p.add_argument(
        "--no-pdf", action="store_true",
        help="Skip PDF generation (HTML only). Note: PDFs are already best-effort.",
    )
    p.add_argument(
        "--keep-intermediate", action="store_true",
        help="Keep PharmCAT raw output directories. Default: deleted after reports.",
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable DEBUG logging.",
    )
    return p


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def main() -> int:
    args = build_parser().parse_args()
    setup_logging(args.verbose)
    log = logging.getLogger("pgx-reporter")

    if not args.vcf.is_file():
        log.error("Input VCF not found: %s", args.vcf)
        return 1

    sample_id = args.sample_id or derive_sample_id(args.vcf)
    sample_out = args.output_dir / sample_id
    sample_out.mkdir(parents=True, exist_ok=True)

    research_mode = (
        [r.strip() for r in args.research.split(",") if r.strip()]
        if args.research else None
    )

    # ── Step 1: Validate ────────────────────────────────────────────────
    if not args.skip_validation:
        log.info("Step: Validating VCF: %s", args.vcf)
        v = validate_vcf(args.vcf)
        for w in v.warnings:
            log.warning("  warning: %s", w)
        if not v.valid:
            for e in v.errors:
                log.error("  error: %s", e)
            log.error("VCF validation failed. Aborting.")
            return 2
        log.info("  VCF passed validation. (size: %s MB, build: %s, samples: %d)",
                 v.info.get("size_mb"), v.info.get("build"), v.info.get("sample_count", 0))

    # Locate PharmCAT JAR
    jar_path = args.jar or find_pharmcat_jar(LIB_DIR)
    if jar_path is None or not jar_path.is_file():
        log.error(
            "PharmCAT JAR not found in %s. Download with:\n"
            "    python src/download_pharmcat.py\n"
            "…or pass --jar /path/to/pharmcat.jar.", LIB_DIR,
        )
        return 3

    pharmcat_dir = sample_out / "pharmcat"
    pharmcat_raw_dir = sample_out / "pharmcat_raw"  # for the position-discovery pass

    current_vcf = args.vcf

    # ── Step 2: Filter (Marco's WGS fast-path) ─────────────────────────
    if args.filter:
        log.info("Step: Filtering VCF to PharmCAT positions...")
        filtered_path = sample_out / "intermediates" / "filtered.vcf.gz"
        try:
            stats = filter_vcf(current_vcf, filtered_path, jar_path)
        except FilterError as e:
            log.error("Filter failed: %s", e)
            return 4
        log.info(
            "  Read %d lines, kept %d PGx variants. %.1f MB → %.3f MB (%.1fx). %.1fs.",
            stats.lines_read, stats.variants_kept,
            stats.input_size_mb, stats.output_size_mb,
            stats.reduction_ratio, stats.elapsed_seconds,
        )
        if stats.variants_kept == 0:
            log.error(
                "Filter kept 0 variants. The VCF may use an incompatible "
                "chromosome naming or build. Re-run without --filter."
            )
            return 4
        current_vcf = filtered_path

    # ── Step 3 + 4: Reference-fill (Adib's GIAB rescue) ────────────────
    if args.ref_fill:
        # Position discovery (initial PharmCAT pass) — unless cached
        if args.reference_phenotype:
            phenotype_json = args.reference_phenotype
            log.info("Step: Using reference phenotype.json: %s", phenotype_json)
        else:
            log.info("Step: Initial PharmCAT pass (discovering PGx positions)...")
            cleaned_path = sample_out / "intermediates" / "cleaned.vcf"
            clean_vcf(current_vcf, cleaned_path)
            initial = run_pharmcat(
                vcf_path=cleaned_path,
                output_dir=pharmcat_raw_dir,
                jar_path=jar_path,
                timeout=args.timeout,
                java_memory=args.java_memory,
                java_executable=args.java_executable,
            )
            if not initial.success or not initial.phenotype_json:
                log.error("Initial PharmCAT pass failed: %s",
                          initial.error_message or "no phenotype.json produced")
                if initial.stderr:
                    log.error("stderr: %s", initial.stderr[:500])
                return 5
            phenotype_json = initial.phenotype_json

        # Coverage screening (informational)
        log.info("Step: Screening PGx coverage...")
        report = screen_vcf(current_vcf, phenotype_json)
        print_screening_report(report)

        # Reference-fill
        log.info("Step: Reference-filling missing PGx positions...")
        filled_path = sample_out / "intermediates" / "preprocessed.vcf"
        fill = reference_fill_vcf(current_vcf, phenotype_json, filled_path,
                                  sample_id=sample_id)
        log.info("  Added %d records across %d genes → %s",
                 fill.added_records, len(fill.genes_filled), fill.output_path)
        current_vcf = filled_path

    # ── Step 5: Final PharmCAT run ─────────────────────────────────────
    if research_mode:
        log.warning(
            "Research mode is on (%s). PharmCAT 3.2.0+ disables the full "
            "reporter in research mode, so the report will contain gene calls "
            "only — drug recommendations and citations will be empty.",
            ",".join(research_mode),
        )
    log.info("Step: Running PharmCAT (final pass)%s...",
             " with research mode" if research_mode else "")
    t0 = time.time()
    result = run_pharmcat(
        vcf_path=current_vcf,
        output_dir=pharmcat_dir,
        jar_path=jar_path,
        timeout=args.timeout,
        research_mode=research_mode,
        java_memory=args.java_memory,
        java_executable=args.java_executable,
    )
    if not result.success:
        log.error("PharmCAT failed: %s", result.error_message)
        if result.stderr:
            log.error("stderr: %s", result.stderr[:1000])
        return 6
    log.info("  PharmCAT finished in %.1fs.", time.time() - t0)

    # ── Step 6: Parse ──────────────────────────────────────────────────
    log.info("Step: Parsing PharmCAT output...")
    parsed = parse_pharmcat_output(
        phenotype_json_path=result.phenotype_json,
        report_json_path=result.report_json,
        match_json_path=result.match_json,
        sample_id=sample_id,
    )
    log.info(
        "  %d definitive, %d ambiguous, %d no-call genes; %d drug recs.",
        len(parsed.definitive_genes), len(parsed.ambiguous_genes),
        len(parsed.no_call_genes), len(parsed.drugs),
    )

    # ── Step 7: Reports ────────────────────────────────────────────────
    log.info("Step: Generating reports...")
    patient_files = generate_patient_report(parsed, sample_out)
    clinician_files = generate_clinician_report(parsed, sample_out)

    # Cleanup intermediates if requested
    if not args.keep_intermediate:
        for path in (pharmcat_dir, pharmcat_raw_dir, sample_out / "intermediates"):
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)

    log.info("─" * 60)
    log.info("Sample:        %s", sample_id)
    log.info("Output dir:    %s", sample_out)
    log.info("  patient_report.html       → %s", patient_files["html"])
    if "pdf" in patient_files:
        log.info("  patient_report.pdf        → %s", patient_files["pdf"])
    log.info("  clinician_report.html     → %s", clinician_files["html"])
    if "pdf" in clinician_files:
        log.info("  clinician_report.pdf      → %s", clinician_files["pdf"])
    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
