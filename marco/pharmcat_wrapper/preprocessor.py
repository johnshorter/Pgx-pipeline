"""
VCF preprocessor -- fast-path filter that keeps only PharmCAT-relevant positions.

PharmCAT itself only looks at a small set of pharmacogene positions (~500
across 23 genes), but when given a whole-genome VCF it still scans every line.
This preprocessor extracts those positions from the PharmCAT JAR and uses them
to stream-filter the input VCF, producing a tiny (<1 MB) VCF that PharmCAT can
process in seconds instead of minutes.

The official PharmCAT Python preprocessor depends on bcftools (painful to
install on Windows). This pure-Python implementation is simpler and has no
external dependencies.
"""

import gzip
import json
import os
import time
import zipfile

from config.settings import PHARMCAT_JAR


class PreprocessorError(Exception):
    """Raised when preprocessing fails."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def preprocess_vcf(
    input_path: str,
    output_path: str,
    jar_path: str | None = None,
    progress_every: int = 5_000_000,
) -> dict:
    """
    Stream-filter a VCF file down to PharmCAT-relevant positions only.

    Args:
        input_path:  Path to the input VCF (.vcf or .vcf.gz).
        output_path: Path for the filtered VCF. Extension determines
                     compression (.gz for gzipped output).
        jar_path:    Path to pharmcat.jar. Defaults to PHARMCAT_JAR.
        progress_every: Print a progress line every N input records.

    Returns:
        dict with stats: lines_read, variants_kept, header_lines,
            positions_loaded, elapsed_seconds, input_size_mb,
            output_size_mb, reduction_ratio.

    Raises:
        PreprocessorError: On I/O errors or malformed input.
    """
    if jar_path is None:
        jar_path = PHARMCAT_JAR

    if not os.path.isfile(input_path):
        raise PreprocessorError(f"Input VCF not found: {input_path}")
    if not os.path.isfile(jar_path):
        raise PreprocessorError(f"PharmCAT JAR not found: {jar_path}")

    # Step 1: Load PGx positions from the PharmCAT JAR
    positions = extract_pgx_positions(jar_path)
    if not positions:
        raise PreprocessorError(
            "No PGx positions extracted from the PharmCAT JAR. "
            "The JAR may be corrupt or from an unsupported version."
        )

    chroms_of_interest = {chrom for chrom, _ in positions}

    # Step 2: Stream filter
    input_size = os.path.getsize(input_path)
    t0 = time.time()

    in_is_gz = input_path.lower().endswith(".gz")
    out_is_gz = output_path.lower().endswith(".gz")

    in_opener = gzip.open if in_is_gz else open
    out_opener = gzip.open if out_is_gz else open
    in_mode = "rt" if in_is_gz else "r"
    out_mode = "wt" if out_is_gz else "w"

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)

    lines_read = 0
    variants_kept = 0
    header_lines = 0

    with in_opener(
        input_path, in_mode, encoding="utf-8", errors="replace"
    ) as inf, out_opener(
        output_path, out_mode, encoding="utf-8"
    ) as outf:

        for line in inf:
            lines_read += 1

            # Headers: keep all of them verbatim
            if line.startswith("#"):
                outf.write(line)
                header_lines += 1
                continue

            # Parse just the first two fields (CHROM \t POS) for speed
            tab1 = line.find("\t")
            if tab1 == -1:
                continue
            chrom = line[:tab1]

            # Early reject: chromosome isn't one PharmCAT cares about
            if chrom not in chroms_of_interest:
                continue

            tab2 = line.find("\t", tab1 + 1)
            if tab2 == -1:
                continue

            pos_str = line[tab1 + 1 : tab2]
            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # Keep the line if (chrom, pos) matches a PGx position.
            # We also probe pos-1 and pos+1 to forgive off-by-one
            # situations between VCF conventions (e.g. indels that start
            # at a different base than PharmCAT's reference position).
            if (
                (chrom, pos) in positions
                or (chrom, pos - 1) in positions
                or (chrom, pos + 1) in positions
            ):
                outf.write(line)
                variants_kept += 1

            # Progress reporting
            if progress_every and lines_read % progress_every == 0:
                elapsed = time.time() - t0
                rate = lines_read / elapsed / 1e6 if elapsed else 0
                print(
                    f"      ...{lines_read // 1_000_000} M lines read, "
                    f"{variants_kept} PGx variants kept "
                    f"({rate:.1f} M/s, {elapsed:.0f}s elapsed)"
                )

    elapsed = time.time() - t0
    output_size = os.path.getsize(output_path)

    return {
        "lines_read": lines_read,
        "variants_kept": variants_kept,
        "header_lines": header_lines,
        "positions_loaded": len(positions),
        "elapsed_seconds": elapsed,
        "input_size_mb": round(input_size / (1024 * 1024), 1),
        "output_size_mb": round(output_size / (1024 * 1024), 3),
        "reduction_ratio": (
            round(input_size / output_size, 1) if output_size > 0 else 0
        ),
    }


# ---------------------------------------------------------------------------
# PharmCAT JAR resource extraction
# ---------------------------------------------------------------------------

def extract_pgx_positions(jar_path: str) -> set[tuple[str, int]]:
    """
    Read all gene allele-definition JSON files from the PharmCAT JAR and
    return the set of (chromosome, position) tuples PharmCAT scans for.
    """
    positions: set[tuple[str, int]] = set()
    prefix = "org/pharmgkb/pharmcat/definition/alleles/"

    with zipfile.ZipFile(jar_path) as z:
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
