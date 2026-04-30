"""
VCF validator — combines Adib's strict header checks with Marco's richer
metadata extraction (build, sample names, contig sniffing).

Returns a ValidationResult dataclass with errors, warnings, and an info dict.
"""

import gzip
import os
from dataclasses import dataclass, field
from pathlib import Path

from config.settings import MAX_VCF_SIZE_MB


@dataclass
class ValidationResult:
    valid: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    info: dict = field(default_factory=dict)


def validate_vcf(vcf_path: str | Path, max_size_mb: int = MAX_VCF_SIZE_MB) -> ValidationResult:
    """Validate a VCF for PharmCAT compatibility.

    Checks:
        - File exists and has a .vcf or .vcf.gz extension
        - File size is within the configured limit and non-empty
        - VCF has a ##fileformat= header and a #CHROM column line
        - All required columns are present
        - Genome build is GRCh38 (rejects GRCh37/hg19)
    """
    path = Path(vcf_path)
    errors: list[str] = []
    warnings: list[str] = []
    info: dict = {}

    if not path.is_file():
        return ValidationResult(False, [f"File not found: {path}"], [], {})

    name_lower = path.name.lower()
    if name_lower.endswith(".vcf.gz"):
        is_gz = True
    elif name_lower.endswith(".vcf"):
        is_gz = False
    else:
        return ValidationResult(
            False,
            [f"Invalid file extension on {path.name}. Expected .vcf or .vcf.gz."],
            [], {},
        )

    size_bytes = path.stat().st_size
    size_mb = size_bytes / (1024 * 1024)
    info["size_mb"] = round(size_mb, 1)
    if size_bytes == 0:
        return ValidationResult(False, ["File is empty."], [], info)
    if size_mb > max_size_mb:
        return ValidationResult(
            False,
            [f"File is {size_mb:.0f} MB, exceeds maximum of {max_size_mb} MB."],
            [], info,
        )

    opener = gzip.open if is_gz else open

    try:
        has_fileformat = False
        has_chrom_line = False
        sample_names: list[str] = []
        contig_lines: list[str] = []
        column_header: str | None = None
        line_count = 0

        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line_count += 1
                stripped = line.rstrip("\n")

                if stripped.startswith("##fileformat="):
                    has_fileformat = True
                    info["fileformat"] = stripped.split("=", 1)[1]

                elif stripped.startswith("##contig="):
                    contig_lines.append(stripped)

                elif stripped.startswith("#CHROM"):
                    has_chrom_line = True
                    column_header = stripped
                    parts = stripped.split("\t")
                    if len(parts) < 8:
                        errors.append("#CHROM line has fewer than 8 required columns.")
                    if len(parts) >= 10:
                        sample_names = parts[9:]
                    else:
                        warnings.append("No sample columns found in header.")
                    break

                # Safety bound: bail if we never find the #CHROM line
                if line_count > 10_000 and not has_chrom_line:
                    errors.append(
                        "Could not find #CHROM header line in the first 10,000 lines."
                    )
                    break

        info["sample_names"] = sample_names
        info["sample_count"] = len(sample_names)

        if not has_fileformat:
            errors.append("Missing ##fileformat= header — file may not be a valid VCF.")
        if not has_chrom_line:
            errors.append("Missing #CHROM column header line.")

        # Required column check
        if column_header is not None:
            cols = column_header.split("\t")
            for required in ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"):
                if required not in cols:
                    errors.append(f"Missing required column: {required}")

        # Genome build detection (header contigs first, fall back to chr-prefix sniff)
        info["build"] = _detect_build(contig_lines, path, is_gz)
        build = info["build"]
        if build and ("GRCh37" in build or "hg19" in build.lower() or "37" in build):
            if "GRCh37" in build or "hg19" in build.lower():
                errors.append(
                    f"VCF appears to use genome build {build}. PharmCAT requires "
                    "GRCh38 (hg38). Please liftover the VCF first."
                )
        elif not build or build.startswith("Unknown"):
            warnings.append(
                "Could not confirm GRCh38 genome build from headers. Ensure the VCF "
                "uses GRCh38 coordinates — PharmCAT will fail otherwise."
            )

    except UnicodeDecodeError:
        errors.append("File encoding error. Ensure the VCF is UTF-8 / ASCII.")
    except Exception as e:
        errors.append(f"Error reading file: {e}")

    return ValidationResult(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        info=info,
    )


def _detect_build(contig_lines: list[str], path: Path, is_gz: bool) -> str:
    """Detect the genome build from contig headers, then chromosome prefix."""
    for line in contig_lines:
        low = line.lower()
        if "grch38" in low or "hg38" in low:
            return "GRCh38"
        if "grch37" in low or "hg19" in low:
            return "GRCh37"

    opener = gzip.open if is_gz else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if line.startswith("chr"):
                    return "GRCh38 (inferred from chr prefix)"
                return "Unknown (no chr prefix, could be GRCh37)"
    except Exception:
        pass
    return "Unknown"
