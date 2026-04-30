"""
WGS filter preprocessor (Marco).

Stream-filters a VCF down to PharmCAT-relevant positions only. Cuts a
30–60 min PharmCAT run on a whole-genome VCF down to seconds.

PharmCAT only inspects ~500 positions across 23 genes. We extract those
positions directly from the JAR's allele-definition JSON files (no
bcftools dependency), then keep only matching VCF lines. Probes
pos-1/pos/pos+1 to forgive indel-anchor differences between conventions.
"""

import gzip
import json
import logging
import os
import time
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class FilterResult:
    output_path: Path
    lines_read: int = 0
    variants_kept: int = 0
    header_lines: int = 0
    positions_loaded: int = 0
    elapsed_seconds: float = 0.0
    input_size_mb: float = 0.0
    output_size_mb: float = 0.0
    reduction_ratio: float = 0.0
    genes_filled: dict[str, int] = field(default_factory=dict)


class FilterError(Exception):
    """Raised when filtering fails."""


def filter_vcf(
    input_vcf: str | Path,
    output_vcf: str | Path,
    jar_path: str | Path,
    progress_every: int = 5_000_000,
) -> FilterResult:
    """Filter a VCF to PharmCAT-relevant positions only.

    Args:
        input_vcf:      Input VCF (.vcf or .vcf.gz).
        output_vcf:     Output path. Extension determines compression
                        (.gz writes gzipped output).
        jar_path:       Path to pharmcat.jar (positions extracted from inside).
        progress_every: Print progress line every N input records.
    """
    in_path = Path(input_vcf)
    out_path = Path(output_vcf)
    jar_path = Path(jar_path)

    if not in_path.is_file():
        raise FilterError(f"Input VCF not found: {in_path}")
    if not jar_path.is_file():
        raise FilterError(f"PharmCAT JAR not found: {jar_path}")

    positions = extract_pgx_positions(jar_path)
    if not positions:
        raise FilterError(
            "No PGx positions extracted from the PharmCAT JAR — "
            "the JAR may be corrupt or from an unsupported version."
        )
    chroms_of_interest = {chrom for chrom, _ in positions}

    in_is_gz = in_path.name.lower().endswith(".gz")
    out_is_gz = out_path.name.lower().endswith(".gz")
    in_opener = gzip.open if in_is_gz else open
    out_opener = gzip.open if out_is_gz else open
    in_mode = "rt" if in_is_gz else "r"
    out_mode = "wt" if out_is_gz else "w"

    out_path.parent.mkdir(parents=True, exist_ok=True)

    input_size = in_path.stat().st_size
    t0 = time.time()
    lines_read = variants_kept = header_lines = 0

    with in_opener(in_path, in_mode, encoding="utf-8", errors="replace") as inf, \
            out_opener(out_path, out_mode, encoding="utf-8") as outf:
        for line in inf:
            lines_read += 1

            if line.startswith("#"):
                outf.write(line)
                header_lines += 1
                continue

            tab1 = line.find("\t")
            if tab1 == -1:
                continue
            chrom = line[:tab1]
            if chrom not in chroms_of_interest:
                continue
            tab2 = line.find("\t", tab1 + 1)
            if tab2 == -1:
                continue
            try:
                pos = int(line[tab1 + 1:tab2])
            except ValueError:
                continue

            if (
                (chrom, pos) in positions
                or (chrom, pos - 1) in positions
                or (chrom, pos + 1) in positions
            ):
                outf.write(line)
                variants_kept += 1

            if progress_every and lines_read % progress_every == 0:
                elapsed = time.time() - t0
                rate = lines_read / elapsed / 1e6 if elapsed else 0
                logger.info(
                    "%d M lines read, %d PGx variants kept (%.1f M/s, %.0fs elapsed)",
                    lines_read // 1_000_000, variants_kept, rate, elapsed,
                )

    elapsed = time.time() - t0
    output_size = out_path.stat().st_size

    return FilterResult(
        output_path=out_path,
        lines_read=lines_read,
        variants_kept=variants_kept,
        header_lines=header_lines,
        positions_loaded=len(positions),
        elapsed_seconds=elapsed,
        input_size_mb=round(input_size / (1024 * 1024), 1),
        output_size_mb=round(output_size / (1024 * 1024), 3),
        reduction_ratio=(
            round(input_size / output_size, 1) if output_size > 0 else 0
        ),
    )


def extract_pgx_positions(jar_path: str | Path) -> set[tuple[str, int]]:
    """Read all gene allele-definition JSON files from the PharmCAT JAR and
    return the set of (chromosome, position) PharmCAT scans for."""
    positions: set[tuple[str, int]] = set()
    prefix = "org/pharmgkb/pharmcat/definition/alleles/"

    with zipfile.ZipFile(os.fspath(jar_path)) as z:
        gene_files = [
            n for n in z.namelist()
            if n.startswith(prefix) and n.endswith("_translation.json")
        ]
        for name in gene_files:
            with z.open(name) as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    continue
            for variant in data.get("variants", []) or []:
                chrom = variant.get("chromosome")
                if not chrom:
                    continue
                for key in ("position", "cpicPosition"):
                    pos = variant.get(key)
                    if isinstance(pos, int):
                        positions.add((chrom, pos))
    return positions
