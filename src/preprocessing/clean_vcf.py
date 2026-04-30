"""
VCF cleaner — fixes FORMAT/sample field-count mismatches.

Some GIAB benchmark VCFs (e.g. HG001/HG002/HG005) have lines where the
number of FORMAT fields doesn't match the number of sample data fields,
which PharmCAT's strict parser rejects. We pad missing trailing sample
fields with '.' to match the FORMAT count.
"""

import gzip
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def fix_format_mismatch(line: str) -> str:
    """Pad sample fields to match the FORMAT field count (if shorter)."""
    cols = line.split("\t")
    if len(cols) < 10:
        return line
    fmt_fields = cols[8].split(":")
    n_fmt = len(fmt_fields)
    for i in range(9, len(cols)):
        sample = cols[i].split(":")
        if len(sample) < n_fmt:
            sample.extend(["."] * (n_fmt - len(sample)))
            cols[i] = ":".join(sample)
    return "\t".join(cols)


def clean_vcf(input_vcf: str | Path, output_vcf: str | Path) -> Path:
    """Clean a VCF by fixing FORMAT/sample mismatches. Output is uncompressed."""
    in_path = Path(input_vcf)
    out_path = Path(output_vcf)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if in_path.name.endswith(".gz") else open
    fixed = 0

    with opener(in_path, "rt", encoding="utf-8", errors="replace") as fin, \
            open(out_path, "w", encoding="utf-8") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line.startswith("#"):
                original = line
                line = fix_format_mismatch(line)
                if line != original:
                    fixed += 1
            fout.write(line + "\n")

    logger.info("Cleaned VCF: fixed %d FORMAT/sample mismatches → %s", fixed, out_path)
    return out_path
