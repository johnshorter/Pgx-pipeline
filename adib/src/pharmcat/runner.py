"""
PharmCAT runner module.

Uses Python's subprocess to execute PharmCAT's Java JAR file.
Constructs the command, executes with a timeout, and captures output/errors.
"""

import subprocess
import logging
from pathlib import Path
from dataclasses import dataclass

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT_SECONDS = 300  # 5 minutes


@dataclass
class PharmCATResult:
    success: bool
    output_dir: Path
    phenotyper_json: Path | None
    reporter_json: Path | None
    reporter_html: Path | None
    stdout: str
    stderr: str
    error_message: str | None


def find_pharmcat_jar(lib_dir: str | Path) -> Path | None:
    """Find the PharmCAT JAR file in the lib directory."""
    lib_dir = Path(lib_dir)
    if not lib_dir.exists():
        return None

    jars = sorted(lib_dir.glob("pharmcat*.jar"), reverse=True)
    return jars[0] if jars else None


def run_pharmcat(
    vcf_path: str | Path,
    output_dir: str | Path,
    jar_path: str | Path,
    timeout: int = DEFAULT_TIMEOUT_SECONDS,
    research_mode: list[str] | None = None,
) -> PharmCATResult:
    """Run PharmCAT on a VCF file.

    Args:
        vcf_path: Path to the input VCF file.
        output_dir: Directory for PharmCAT output files.
        jar_path: Path to the PharmCAT JAR file.
        timeout: Maximum execution time in seconds.
        research_mode: List of research mode flags (e.g., ["cyp2d6"]).

    Returns:
        PharmCATResult with paths to output files and any error info.
    """
    vcf_path = Path(vcf_path).resolve()
    output_dir = Path(output_dir).resolve()
    jar_path = Path(jar_path).resolve()

    # Validate inputs
    if not vcf_path.exists():
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotyper_json=None,
            reporter_json=None,
            reporter_html=None,
            stdout="",
            stderr="",
            error_message=f"VCF file not found: {vcf_path}",
        )

    if not jar_path.exists():
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotyper_json=None,
            reporter_json=None,
            reporter_html=None,
            stdout="",
            stderr="",
            error_message=f"PharmCAT JAR not found: {jar_path}",
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    # Find Java binary: check for local JDK in lib/jdk first, then system PATH
    java_bin = "java"
    local_jdk = jar_path.parent / "jdk" / "Contents" / "Home" / "bin" / "java"
    if not local_jdk.exists():
        local_jdk = jar_path.parent / "jdk" / "bin" / "java"
    if local_jdk.exists():
        java_bin = str(local_jdk)

    # Build PharmCAT command
    # PharmCAT expects: java -jar pharmcat.jar -vcf <input> -o <output_dir>
    cmd = [
        java_bin,
        "-jar",
        str(jar_path),
        "-vcf",
        str(vcf_path),
        "-o",
        str(output_dir),
    ]

    if research_mode:
        cmd.extend(["-research", ",".join(research_mode)])

    logger.info("Running PharmCAT: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        stdout = result.stdout
        stderr = result.stderr

        if result.returncode != 0:
            logger.error("PharmCAT failed (exit code %d): %s", result.returncode, stderr)
            return PharmCATResult(
                success=False,
                output_dir=output_dir,
                phenotyper_json=None,
                reporter_json=None,
                reporter_html=None,
                stdout=stdout,
                stderr=stderr,
                error_message=f"PharmCAT exited with code {result.returncode}. {stderr[:500]}",
            )

    except subprocess.TimeoutExpired:
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotyper_json=None,
            reporter_json=None,
            reporter_html=None,
            stdout="",
            stderr="",
            error_message=f"PharmCAT timed out after {timeout} seconds.",
        )
    except FileNotFoundError:
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotyper_json=None,
            reporter_json=None,
            reporter_html=None,
            stdout="",
            stderr="",
            error_message="Java not found. Please install Java 17+ and ensure it is on PATH.",
        )

    # Find output files - PharmCAT names them based on input VCF stem
    vcf_stem = vcf_path.stem
    if vcf_stem.endswith(".vcf"):
        vcf_stem = vcf_stem[:-4]

    # PharmCAT v3 uses .phenotype.json (not .phenotyper.json) and .report.html (no .report.json)
    phenotyper_json = _find_output(output_dir, vcf_stem, ".phenotype.json")
    if not phenotyper_json:
        phenotyper_json = _find_output(output_dir, vcf_stem, ".phenotyper.json")
    reporter_json = _find_output(output_dir, vcf_stem, ".report.json")
    reporter_html = _find_output(output_dir, vcf_stem, ".report.html")

    if not phenotyper_json and not reporter_json:
        # Try finding any json files PharmCAT may have produced
        all_json = list(output_dir.glob("*.json"))
        logger.warning(
            "Expected output files not found. Found: %s",
            [f.name for f in all_json],
        )
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotyper_json=None,
            reporter_json=None,
            reporter_html=None,
            stdout=stdout,
            stderr=stderr,
            error_message="PharmCAT ran but expected output files were not found.",
        )

    logger.info("PharmCAT completed successfully.")
    return PharmCATResult(
        success=True,
        output_dir=output_dir,
        phenotyper_json=phenotyper_json,
        reporter_json=reporter_json,
        reporter_html=reporter_html,
        stdout=stdout,
        stderr=stderr,
        error_message=None,
    )


def _find_output(output_dir: Path, stem: str, suffix: str) -> Path | None:
    """Find a PharmCAT output file by stem and suffix."""
    # Try exact match first
    exact = output_dir / f"{stem}{suffix}"
    if exact.exists():
        return exact

    # Try glob for any file with that suffix
    matches = list(output_dir.glob(f"*{suffix}"))
    return matches[0] if matches else None
