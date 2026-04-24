"""
PharmCAT runner — executes PharmCAT as a subprocess and returns output paths.
"""

import os
import subprocess
import tempfile
import shutil

from config.settings import PHARMCAT_JAR, JAVA_EXECUTABLE, PHARMCAT_TIMEOUT, TEMP_OUTPUT_DIR


class PharmCATError(Exception):
    """Raised when PharmCAT execution fails."""


def run_pharmcat(
    vcf_path: str,
    output_dir: str | None = None,
    timeout: int | None = None,
    java_memory: str | None = None,
) -> dict:
    """
    Run PharmCAT on a VCF file.

    Args:
        vcf_path: Path to the input VCF file.
        output_dir: Directory for output files. Defaults to TEMP_OUTPUT_DIR.
        timeout: Max seconds to wait. Defaults to PHARMCAT_TIMEOUT.
        java_memory: Java max heap size passed as -Xmx (e.g. "4g", "8g").
            If None, the JVM default is used.

    Returns:
        dict with keys:
            success (bool)
            output_dir (str): Directory containing output files.
            report_json (str): Path to the .report.json file.
            phenotype_json (str): Path to the .phenotype.json file.
            match_json (str): Path to the .match.json file.
            stderr (str): Any stderr output from PharmCAT.
    """
    if not os.path.isfile(PHARMCAT_JAR):
        raise PharmCATError(
            f"PharmCAT JAR not found at {PHARMCAT_JAR}. "
            "Download it from https://github.com/PharmGKB/PharmCAT/releases"
        )

    if timeout is None:
        timeout = PHARMCAT_TIMEOUT

    if output_dir is None:
        output_dir = TEMP_OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    # Build command
    cmd = [JAVA_EXECUTABLE]
    if java_memory:
        cmd.append(f"-Xmx{java_memory}")
    cmd += [
        "-jar", PHARMCAT_JAR,
        "-vcf", vcf_path,
        "-o", output_dir,
        "-reporterJson",  # We need the JSON for parsing
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        raise PharmCATError(
            f"PharmCAT did not finish within {timeout} seconds. "
            "The VCF file may be too large or Java may be running out of memory."
        )
    except FileNotFoundError:
        raise PharmCATError(
            f"Java not found at '{JAVA_EXECUTABLE}'. "
            "Please install Java 17+ and ensure it is on your PATH."
        )

    # PharmCAT writes output files based on the input VCF basename
    base_name = os.path.splitext(os.path.basename(vcf_path))[0]
    if base_name.endswith(".vcf"):
        base_name = base_name[:-4]  # Handle .vcf.gz case

    report_json = os.path.join(output_dir, f"{base_name}.report.json")
    phenotype_json = os.path.join(output_dir, f"{base_name}.phenotype.json")
    match_json = os.path.join(output_dir, f"{base_name}.match.json")

    # Check for success
    if result.returncode != 0:
        raise PharmCATError(
            f"PharmCAT exited with code {result.returncode}.\n"
            f"stderr: {result.stderr[:2000]}"
        )

    if not os.path.isfile(report_json):
        raise PharmCATError(
            f"PharmCAT completed but report JSON not found at {report_json}.\n"
            f"stdout: {result.stdout[:1000]}\n"
            f"stderr: {result.stderr[:1000]}"
        )

    return {
        "success": True,
        "output_dir": output_dir,
        "report_json": report_json,
        "phenotype_json": phenotype_json,
        "match_json": match_json,
        "stderr": result.stderr,
    }


def cleanup_output(output_dir: str) -> None:
    """Remove a PharmCAT output directory and all its contents."""
    if os.path.isdir(output_dir) and output_dir != TEMP_OUTPUT_DIR:
        shutil.rmtree(output_dir, ignore_errors=True)
