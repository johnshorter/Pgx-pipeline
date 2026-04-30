"""
PharmCAT runner — invokes the PharmCAT JAR via subprocess.

Combines Adib's auto-detect / dataclass result with Marco's heap-size
control. Always passes -reporterJson so report.json is produced
(needed by the unified parser). Locates all four output artefacts:
phenotype.json, match.json, report.json, report.html.
"""

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

from config.settings import (
    JAVA_EXECUTABLE, LIB_DIR, PHARMCAT_JAR_DEFAULT, PHARMCAT_TIMEOUT_DEFAULT,
)

logger = logging.getLogger(__name__)


@dataclass
class PharmCATResult:
    success: bool
    output_dir: Path
    phenotype_json: Path | None
    match_json: Path | None
    report_json: Path | None
    report_html: Path | None
    stdout: str
    stderr: str
    error_message: str | None


def find_pharmcat_jar(lib_dir: str | Path = LIB_DIR) -> Path | None:
    """Find the PharmCAT JAR. Checks `pharmcat*.jar` (newest first), then
    falls back to a hardcoded `pharmcat.jar` in the lib directory."""
    lib_dir = Path(lib_dir)
    if not lib_dir.is_dir():
        return None
    matches = sorted(lib_dir.glob("pharmcat*.jar"), reverse=True)
    if matches:
        return matches[0]
    if PHARMCAT_JAR_DEFAULT.is_file():
        return PHARMCAT_JAR_DEFAULT
    return None


def run_pharmcat(
    vcf_path: str | Path,
    output_dir: str | Path,
    jar_path: str | Path,
    timeout: int = PHARMCAT_TIMEOUT_DEFAULT,
    research_mode: list[str] | None = None,
    java_memory: str | None = None,
    java_executable: str = JAVA_EXECUTABLE,
) -> PharmCATResult:
    """Run PharmCAT on a VCF.

    Args:
        vcf_path:        Input VCF (absolute path preferred).
        output_dir:      Destination directory for PharmCAT output.
        jar_path:        Path to pharmcat.jar.
        timeout:         Subprocess timeout in seconds.
        research_mode:   List of `-research` flags (e.g. ["cyp2d6"]).
        java_memory:     `-Xmx` heap size (e.g. "4g") or None for JVM default.
        java_executable: Java binary to invoke.

    Returns:
        PharmCATResult — always returned (no exceptions raised on failure;
        check `success` and `error_message`).
    """
    vcf_path = Path(vcf_path).resolve()
    output_dir = Path(output_dir).resolve()
    jar_path = Path(jar_path).resolve()

    if not vcf_path.is_file():
        return _fail(output_dir, f"VCF file not found: {vcf_path}")
    if not jar_path.is_file():
        return _fail(output_dir, f"PharmCAT JAR not found: {jar_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [java_executable]
    if java_memory:
        cmd.append(f"-Xmx{java_memory}")
    cmd += [
        "-jar", str(jar_path),
        "-vcf", str(vcf_path),
        "-o", str(output_dir),
        "-reporterJson",  # ensure report.json is written
    ]
    if research_mode:
        cmd += ["-research", ",".join(research_mode)]

    logger.info("Running PharmCAT: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return _fail(output_dir, f"PharmCAT timed out after {timeout} seconds.")
    except FileNotFoundError:
        return _fail(
            output_dir,
            f"Java not found at '{java_executable}'. Install Java 17+ and ensure "
            "it is on PATH or pass --java-executable.",
        )

    stdout, stderr = result.stdout, result.stderr

    if result.returncode != 0:
        logger.error("PharmCAT failed (exit code %d): %s", result.returncode, stderr)
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotype_json=None, match_json=None,
            report_json=None, report_html=None,
            stdout=stdout, stderr=stderr,
            error_message=(
                f"PharmCAT exited with code {result.returncode}. {stderr[:500]}"
            ),
        )

    # Locate output artefacts. PharmCAT v3 names them by VCF stem.
    stem = vcf_path.name
    if stem.endswith(".vcf.gz"):
        stem = stem[: -len(".vcf.gz")]
    elif stem.endswith(".vcf"):
        stem = stem[: -len(".vcf")]

    phenotype_json = _find_output(output_dir, stem, ".phenotype.json")
    if not phenotype_json:
        # Older PharmCAT versions use .phenotyper.json
        phenotype_json = _find_output(output_dir, stem, ".phenotyper.json")
    match_json = _find_output(output_dir, stem, ".match.json")
    report_json = _find_output(output_dir, stem, ".report.json")
    report_html = _find_output(output_dir, stem, ".report.html")

    if not (phenotype_json or report_json):
        return PharmCATResult(
            success=False,
            output_dir=output_dir,
            phenotype_json=None, match_json=None,
            report_json=None, report_html=None,
            stdout=stdout, stderr=stderr,
            error_message="PharmCAT ran but expected output files were not found.",
        )

    return PharmCATResult(
        success=True,
        output_dir=output_dir,
        phenotype_json=phenotype_json,
        match_json=match_json,
        report_json=report_json,
        report_html=report_html,
        stdout=stdout, stderr=stderr,
        error_message=None,
    )


def _find_output(output_dir: Path, stem: str, suffix: str) -> Path | None:
    exact = output_dir / f"{stem}{suffix}"
    if exact.is_file():
        return exact
    matches = list(output_dir.glob(f"*{suffix}"))
    return matches[0] if matches else None


def _fail(output_dir: Path, message: str) -> PharmCATResult:
    return PharmCATResult(
        success=False,
        output_dir=output_dir,
        phenotype_json=None, match_json=None,
        report_json=None, report_html=None,
        stdout="", stderr="",
        error_message=message,
    )
