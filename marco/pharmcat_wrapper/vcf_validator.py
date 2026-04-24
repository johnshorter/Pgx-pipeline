"""
VCF file validation before PharmCAT processing.

Checks file format, required headers, and basic integrity
so we can give clear error messages before running PharmCAT.
"""

import os


class VCFValidationError(Exception):
    """Raised when a VCF file fails validation."""


def validate_vcf(file_path: str, max_size_mb: int = 3000) -> dict:
    """
    Validate a VCF file for PharmCAT compatibility.

    Returns a dict with:
        valid (bool): Whether the file passed all checks.
        errors (list[str]): List of error messages (empty if valid).
        warnings (list[str]): Non-fatal issues.
        info (dict): Metadata extracted (sample count, contig build, etc.)
    """
    errors = []
    warnings = []
    info = {}

    # --- Check file exists ---
    if not os.path.isfile(file_path):
        return {"valid": False, "errors": [f"File not found: {file_path}"], "warnings": [], "info": {}}

    # --- Check extension ---
    lower = file_path.lower()
    if not (lower.endswith(".vcf") or lower.endswith(".vcf.gz")):
        errors.append("File must have a .vcf or .vcf.gz extension.")

    # --- Check file size ---
    size_bytes = os.path.getsize(file_path)
    size_mb = size_bytes / (1024 * 1024)
    info["size_mb"] = round(size_mb, 1)
    if size_mb > max_size_mb:
        errors.append(f"File is {size_mb:.0f} MB, exceeds maximum of {max_size_mb} MB.")
    if size_bytes == 0:
        errors.append("File is empty.")
        return {"valid": False, "errors": errors, "warnings": warnings, "info": info}

    # --- Read and check header lines ---
    try:
        opener = open
        mode = "r"
        if lower.endswith(".vcf.gz"):
            import gzip
            opener = gzip.open
            mode = "rt"

        with opener(file_path, mode, encoding="utf-8", errors="replace") as f:
            has_fileformat = False
            has_chrom_line = False
            sample_names = []
            contig_lines = []
            line_count = 0

            for line in f:
                line_count += 1

                if line.startswith("##fileformat="):
                    has_fileformat = True
                    info["fileformat"] = line.strip().split("=", 1)[1]

                elif line.startswith("##contig="):
                    contig_lines.append(line.strip())

                elif line.startswith("#CHROM"):
                    has_chrom_line = True
                    parts = line.strip().split("\t")
                    # Standard VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE...
                    if len(parts) < 8:
                        errors.append("Header line has fewer than 8 required columns.")
                    if len(parts) >= 10:
                        sample_names = parts[9:]
                    else:
                        warnings.append("No sample columns found in header.")
                    info["sample_names"] = sample_names
                    info["sample_count"] = len(sample_names)
                    break  # Stop after header

                # Safety: don't read the entire file just for headers
                if line_count > 10000 and not has_chrom_line:
                    errors.append("Could not find #CHROM header line in the first 10,000 lines.")
                    break

            if not has_fileformat:
                errors.append("Missing ##fileformat= header. This may not be a valid VCF file.")

            if not has_chrom_line:
                errors.append("Missing #CHROM header line.")

            # --- Check for GRCh38 build ---
            info["build"] = _detect_build(contig_lines, file_path, lower.endswith(".vcf.gz"))
            if info["build"] and "37" in info["build"]:
                errors.append(
                    f"VCF appears to use genome build {info['build']}. "
                    "PharmCAT requires GRCh38 (hg38). Please lift over your VCF."
                )

    except Exception as e:
        errors.append(f"Error reading file: {e}")

    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "info": info,
    }


def _detect_build(contig_lines: list, file_path: str, is_gzipped: bool) -> str:
    """Try to detect the genome build from contig headers or chromosome naming."""
    for line in contig_lines:
        low = line.lower()
        if "grch38" in low or "hg38" in low:
            return "GRCh38"
        if "grch37" in low or "hg19" in low:
            return "GRCh37"

    # Fall back: check if chromosomes use "chr" prefix (typical of GRCh38)
    try:
        opener = open
        mode = "r"
        if is_gzipped:
            import gzip
            opener = gzip.open
            mode = "rt"

        with opener(file_path, mode, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                # First data line
                if line.startswith("chr"):
                    return "GRCh38 (inferred from chr prefix)"
                else:
                    return "Unknown (no chr prefix, could be GRCh37)"
    except Exception:
        pass

    return "Unknown"
