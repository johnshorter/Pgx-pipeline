"""
VCF file validator for PharmCAT input.

Checks file format, header lines, genome build (GRCh38), and file size
before sending to PharmCAT. Prevents confusing downstream errors.
"""

from pathlib import Path
from dataclasses import dataclass

MAX_FILE_SIZE_MB = 500


@dataclass
class ValidationResult:
    valid: bool
    errors: list[str]
    warnings: list[str]


def validate_vcf(vcf_path: str | Path) -> ValidationResult:
    """Validate a VCF file for PharmCAT compatibility.

    Checks:
    - File exists and has .vcf or .vcf.gz extension
    - File size is within limits
    - Has valid VCF header (##fileformat=VCFv4.x)
    - Contains required header columns
    - Genome build is GRCh38 (required by PharmCAT)

    Returns a ValidationResult with errors/warnings.
    """
    vcf_path = Path(vcf_path)
    errors: list[str] = []
    warnings: list[str] = []

    # Check file exists
    if not vcf_path.exists():
        return ValidationResult(False, [f"File not found: {vcf_path}"], [])

    # Check extension — use the last 1-2 suffixes to handle complex filenames
    name_lower = vcf_path.name.lower()
    if name_lower.endswith(".vcf.gz"):
        suffixes = ".vcf.gz"
    elif name_lower.endswith(".vcf"):
        suffixes = ".vcf"
    else:
        suffixes = "".join(vcf_path.suffixes).lower()
    if suffixes not in (".vcf", ".vcf.gz"):
        errors.append(
            f"Invalid file extension '{suffixes}'. Expected .vcf or .vcf.gz"
        )
        return ValidationResult(False, errors, warnings)

    # Check file size
    size_mb = vcf_path.stat().st_size / (1024 * 1024)
    if size_mb > MAX_FILE_SIZE_MB:
        errors.append(
            f"File too large ({size_mb:.1f} MB). Maximum is {MAX_FILE_SIZE_MB} MB."
        )
        return ValidationResult(False, errors, warnings)

    if size_mb == 0:
        errors.append("File is empty.")
        return ValidationResult(False, errors, warnings)

    # For .vcf.gz files, we can't easily read headers without gzip
    if suffixes == ".vcf.gz":
        import gzip

        opener = gzip.open
    else:
        opener = open

    try:
        with opener(vcf_path, "rt") as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.rstrip("\n"))
                if i > 500:
                    break

            if not lines:
                errors.append("File is empty.")
                return ValidationResult(False, errors, warnings)

            # Check VCF format header
            if not lines[0].startswith("##fileformat=VCF"):
                errors.append(
                    f"Missing VCF format header. First line: '{lines[0][:60]}...'"
                )
                return ValidationResult(False, errors, warnings)

            # Parse meta-information lines and find column header
            meta_lines = []
            column_header = None
            has_grch38 = False

            for line in lines:
                if line.startswith("##"):
                    meta_lines.append(line)
                    lower = line.lower()
                    # Check for genome build references
                    if "grch38" in lower or "hg38" in lower:
                        has_grch38 = True
                    # Check for GRCh37/hg19 (wrong build)
                    if "grch37" in lower or "hg19" in lower:
                        errors.append(
                            "VCF appears to use GRCh37/hg19 genome build. "
                            "PharmCAT requires GRCh38/hg38. "
                            "Please use a liftover tool to convert."
                        )
                elif line.startswith("#CHROM"):
                    column_header = line
                    break

            if column_header is None:
                errors.append(
                    "Missing column header line (#CHROM POS ID REF ALT ...)."
                )
                return ValidationResult(False, errors, warnings)

            # Check required columns
            columns = column_header.split("\t")
            required = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
            for col in required:
                if col not in columns:
                    errors.append(f"Missing required column: {col}")

            # Check for sample column(s)
            if len(columns) < 10:
                warnings.append(
                    "No sample columns found. PharmCAT may need at least one sample."
                )

            # Warn if genome build is not explicitly stated
            if not has_grch38 and not errors:
                warnings.append(
                    "Could not confirm GRCh38 genome build from headers. "
                    "PharmCAT requires GRCh38. Ensure your VCF uses GRCh38 coordinates."
                )

    except UnicodeDecodeError:
        errors.append("File encoding error. Ensure the VCF is UTF-8 or ASCII encoded.")
    except Exception as e:
        errors.append(f"Error reading file: {e}")

    return ValidationResult(valid=len(errors) == 0, errors=errors, warnings=warnings)
