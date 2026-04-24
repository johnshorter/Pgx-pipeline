"""
PGx Reporter -- End-to-end command-line pipeline.

Takes a VCF file and produces both the patient-friendly and the
clinician/researcher pharmacogenomic reports (HTML + PDF).

Pipeline:
    VCF validation -> PharmCAT -> Output parser -> Patient + Clinician reports

Usage:
    python generate_reports.py path/to/sample.vcf
    python generate_reports.py path/to/sample.vcf.gz --output my_results
    python generate_reports.py sample.vcf --sample-id HG005 --keep-intermediate
    python generate_reports.py huge_wgs.vcf.gz --preprocess --java-memory 4g

Outputs:
    <output_dir>/
        <sample-id>/
            patient_report.html     <- patient-friendly report
            patient_report.pdf      (if WeasyPrint available)
            clinician_report.html   <- clinician/researcher report
            clinician_report.pdf    (if WeasyPrint available)
            pharmcat/               (PharmCAT raw JSON output, optional)
"""

import argparse
import os
import shutil
import sys
import time


def _print_header(title: str) -> None:
    line = "=" * 70
    print(f"\n{line}\n  {title}\n{line}")


def _print_step(step_num: int, total_steps: int, title: str) -> None:
    print(f"\n[{step_num}/{total_steps}] {title}")


def _derive_sample_id(vcf_path: str) -> str:
    """Get a clean sample ID from the VCF filename (drop .vcf / .vcf.gz)."""
    basename = os.path.basename(vcf_path)
    for suffix in (".vcf.gz", ".vcf"):
        if basename.lower().endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    return basename


def run_pipeline(
    vcf_path: str,
    output_dir: str,
    sample_id: str | None = None,
    keep_intermediate: bool = False,
    skip_validation: bool = False,
    timeout: int | None = None,
    java_memory: str | None = None,
    preprocess: bool = False,
) -> dict:
    """
    Execute the full PGx Reporter pipeline on a VCF file.

    Args:
        vcf_path: Path to the input VCF (.vcf or .vcf.gz) file.
        output_dir: Directory where per-sample output folders will be written.
        sample_id: Optional sample identifier. Defaults to the VCF basename.
        keep_intermediate: If True, keep PharmCAT's raw JSON files in a
            pharmcat/ subdirectory; otherwise remove them after the reports
            have been generated.
        skip_validation: If True, skip the VCF validation step (useful when
            you already know the file is compatible).
        timeout: Max seconds to wait for PharmCAT. Defaults to PHARMCAT_TIMEOUT.
        java_memory: Java max heap size passed as -Xmx (e.g. "4g", "8g").
        preprocess: If True, pre-filter the VCF down to PharmCAT-relevant
            positions before invoking PharmCAT. Recommended for
            whole-genome VCFs.

    Returns:
        Dict with keys: sample_id, sample_output_dir, patient_html,
        patient_pdf, clinician_html, clinician_pdf.

    Raises:
        FileNotFoundError: VCF file not found.
        ValueError: VCF failed validation.
        pharmcat_wrapper.runner.PharmCATError: PharmCAT execution failed.
    """
    # Lazy imports so the script can display --help without the config loading
    from pharmcat_wrapper.vcf_validator import validate_vcf
    from pharmcat_wrapper.runner import run_pharmcat
    from pharmcat_wrapper.preprocessor import preprocess_vcf
    from reports.patient_report import generate_patient_report
    from reports.clinician_report import generate_clinician_report

    if sample_id is None:
        sample_id = _derive_sample_id(vcf_path)

    # Per-sample output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    pharmcat_dir = os.path.join(sample_output_dir, "pharmcat")

    # Total steps: validation (optional) + optional preprocess + pharmcat
    # + patient + clinician + optional cleanup
    total_steps = 3  # pharmcat + patient + clinician
    if not skip_validation:
        total_steps += 1
    if preprocess:
        total_steps += 1
    if not keep_intermediate:
        total_steps += 1
    step = 0

    # ---------------------------------------------------------------
    # Step 1: Validate VCF
    # ---------------------------------------------------------------
    if not skip_validation:
        step += 1
        _print_step(step, total_steps, "Validating VCF file")

        if not os.path.isfile(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        result = validate_vcf(vcf_path)
        info = result.get("info", {})
        if info.get("size_mb") is not None:
            print(f"      File size:    {info['size_mb']} MB")
        if info.get("build"):
            print(f"      Genome build: {info['build']}")
        if info.get("sample_count") is not None:
            print(f"      Samples:      {info['sample_count']}")

        for warn in result.get("warnings", []):
            print(f"      [warning] {warn}")

        if not result.get("valid"):
            print("      [error] VCF validation failed:")
            for err in result.get("errors", []):
                print(f"        - {err}")
            raise ValueError("VCF validation failed. See errors above.")

        print("      VCF is valid.")

    # ---------------------------------------------------------------
    # Optional step: Preprocess VCF (keep only PharmCAT-relevant positions)
    # ---------------------------------------------------------------
    vcf_for_pharmcat = vcf_path
    if preprocess:
        step += 1
        _print_step(step, total_steps, "Preprocessing VCF (fast-path filter)")

        os.makedirs(pharmcat_dir, exist_ok=True)
        filtered_vcf = os.path.join(pharmcat_dir, "prefiltered.vcf.gz")
        print(f"      Output: {filtered_vcf}")

        try:
            stats = preprocess_vcf(vcf_path, filtered_vcf)
        except Exception as e:
            raise RuntimeError(f"VCF preprocessing failed: {e}") from e

        print(
            f"      Loaded {stats['positions_loaded']} PGx positions "
            f"from PharmCAT JAR."
        )
        print(
            f"      Read {stats['lines_read']:,} lines, "
            f"kept {stats['variants_kept']:,} PGx variants "
            f"({stats['header_lines']:,} header lines)."
        )
        print(
            f"      Size: {stats['input_size_mb']} MB -> "
            f"{stats['output_size_mb']} MB "
            f"(reduction {stats['reduction_ratio']}x)."
        )
        print(f"      Elapsed: {stats['elapsed_seconds']:.1f}s.")

        if stats["variants_kept"] == 0:
            raise RuntimeError(
                "Preprocessor kept 0 variants -- the VCF may use an "
                "incompatible chromosome naming convention (e.g. '1' vs "
                "'chr1') or a different genome build. Re-run without "
                "--preprocess to let PharmCAT handle the raw VCF."
            )

        vcf_for_pharmcat = filtered_vcf

    # ---------------------------------------------------------------
    # Step 2: Run PharmCAT
    # ---------------------------------------------------------------
    step += 1
    _print_step(step, total_steps, "Running PharmCAT")
    if preprocess:
        print(f"      Using prefiltered VCF (fast).")
    else:
        print(f"      (Whole-genome VCFs can take 30-60 minutes.)")
    if java_memory:
        print(f"      Java heap: -Xmx{java_memory}")
    if timeout:
        print(f"      Timeout:   {timeout}s")

    t0 = time.time()
    pharmcat_result = run_pharmcat(
        vcf_for_pharmcat,
        output_dir=pharmcat_dir,
        timeout=timeout,
        java_memory=java_memory,
    )
    elapsed = time.time() - t0
    print(f"      PharmCAT finished in {elapsed:.1f}s.")
    print(f"      Output directory: {pharmcat_dir}")

    report_json = pharmcat_result["report_json"]
    match_json = pharmcat_result["match_json"]

    # ---------------------------------------------------------------
    # Step 3: Patient report
    # ---------------------------------------------------------------
    step += 1
    _print_step(step, total_steps, "Generating patient report")

    patient_html, patient_pdf = generate_patient_report(
        report_json_path=report_json,
        match_json_path=match_json,
        output_dir=sample_output_dir,
    )
    print(f"      HTML: {patient_html}")
    print(f"      PDF:  {patient_pdf}")

    # ---------------------------------------------------------------
    # Step 4: Clinician report
    # ---------------------------------------------------------------
    step += 1
    _print_step(step, total_steps, "Generating clinician report")

    clinician_html, clinician_pdf = generate_clinician_report(
        report_json_path=report_json,
        match_json_path=match_json,
        output_dir=sample_output_dir,
    )
    print(f"      HTML: {clinician_html}")
    print(f"      PDF:  {clinician_pdf}")

    # ---------------------------------------------------------------
    # Step 5: Cleanup intermediate files
    # ---------------------------------------------------------------
    if not keep_intermediate and os.path.isdir(pharmcat_dir):
        step += 1
        _print_step(step, total_steps, "Cleaning up intermediate files")
        shutil.rmtree(pharmcat_dir, ignore_errors=True)
        print(f"      Removed: {pharmcat_dir}")

    return {
        "sample_id": sample_id,
        "sample_output_dir": sample_output_dir,
        "patient_html": patient_html,
        "patient_pdf": patient_pdf,
        "clinician_html": clinician_html,
        "clinician_pdf": clinician_pdf,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="generate_reports.py",
        description=(
            "PGx Reporter -- generate patient and clinician "
            "pharmacogenomic reports from a VCF file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python generate_reports.py sample.vcf\n"
            "  python generate_reports.py sample.vcf.gz --output my_results\n"
            "  python generate_reports.py sample.vcf --sample-id HG005 "
            "--keep-intermediate\n"
        ),
    )
    parser.add_argument(
        "vcf",
        help="Path to the input VCF file (.vcf or .vcf.gz, GRCh38).",
    )
    parser.add_argument(
        "-o", "--output",
        default="output",
        help="Directory where the per-sample reports folder will be written "
             "(default: ./output).",
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        help="Sample identifier (used for the output subfolder name). "
             "Defaults to the VCF filename without extension.",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep PharmCAT's raw JSON output files (in the pharmcat/ "
             "subfolder) after report generation. Useful for debugging.",
    )
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip the VCF validation step.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Max seconds to wait for PharmCAT. Default: 3600 (60 min). "
             "Increase for very large VCFs on slow machines.",
    )
    parser.add_argument(
        "--java-memory",
        default=None,
        help="Java max heap size passed as -Xmx (e.g. '4g', '8g'). "
             "Recommended for large VCFs. Default: JVM default.",
    )
    parser.add_argument(
        "--preprocess",
        action="store_true",
        help="Pre-filter the VCF down to PharmCAT-relevant positions "
             "before running PharmCAT. Recommended for whole-genome VCFs "
             "(reduces PharmCAT runtime from ~30-60 min to ~1 min).",
    )

    args = parser.parse_args()

    _print_header("PGx Reporter -- End-to-end report generation")
    print(f"Input VCF:    {args.vcf}")
    print(f"Output root:  {os.path.abspath(args.output)}")

    try:
        result = run_pipeline(
            vcf_path=args.vcf,
            output_dir=args.output,
            sample_id=args.sample_id,
            keep_intermediate=args.keep_intermediate,
            skip_validation=args.skip_validation,
            timeout=args.timeout,
            java_memory=args.java_memory,
            preprocess=args.preprocess,
        )
    except FileNotFoundError as e:
        print(f"\n[error] {e}", file=sys.stderr)
        return 2
    except ValueError as e:
        print(f"\n[error] {e}", file=sys.stderr)
        return 3
    except Exception as e:
        print(f"\n[error] Pipeline failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    _print_header("All done!")
    print(f"Sample:                {result['sample_id']}")
    print(f"Reports directory:     {result['sample_output_dir']}")
    print(f"  Patient report:      {os.path.basename(result['patient_html'])}")
    print(f"  Clinician report:    {os.path.basename(result['clinician_html'])}")
    print()
    print("Open the HTML files in your browser to view the reports.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
