"""
VCF preprocessor for pharmacogenomics.

Fills missing PGx positions with reference genotype calls (0/0) so that
PharmCAT can make diplotype assignments.  This is valid for GIAB benchmark
VCFs where absence of a position means the sample matches the reference.

The position list comes from a PharmCAT phenotype.json file, which records
every position PharmCAT checked including those with no genotype call.
"""

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path

from src.screening.pgx_positions import PgxPosition, load_pgx_positions

logger = logging.getLogger(__name__)

# Chromosome sort order
CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})

# Genes to skip filling (need specialized tools)
SKIP_GENES = {"HLA-A", "HLA-B", "MT-RNR1"}


@dataclass
class PreprocessResult:
    output_path: Path
    original_records: int
    added_records: int
    total_records: int
    genes_filled: dict[str, int] = field(default_factory=dict)


def _fix_format_mismatch(line: str) -> str:
    """Fix VCF lines where FORMAT field count doesn't match sample field count.

    Some GIAB benchmark VCFs have lines like:
        GT:AD:PS    1/0:1,1   (FORMAT has 3 fields, sample has 2)
    PharmCAT's VCF parser rejects these. We pad with '.' to match.
    """
    cols = line.split("\t")
    if len(cols) < 10:
        return line

    fmt_fields = cols[8].split(":")
    for i in range(9, len(cols)):
        sample_fields = cols[i].split(":")
        if len(sample_fields) < len(fmt_fields):
            sample_fields.extend(["."] * (len(fmt_fields) - len(sample_fields)))
            cols[i] = ":".join(sample_fields)

    return "\t".join(cols)


def clean_vcf(
    input_vcf: str | Path,
    output_vcf: str | Path,
) -> Path:
    """Clean a VCF file by fixing FORMAT/sample field mismatches.

    Some GIAB benchmark VCFs have lines where the number of FORMAT fields
    doesn't match the number of sample data fields.  PharmCAT's strict
    parser rejects these.

    Returns:
        Path to the cleaned VCF.
    """
    input_path = Path(input_vcf)
    output_path = Path(output_vcf)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if input_path.name.endswith(".gz") else open
    fixed_count = 0

    with opener(input_path, "rt") as fin, open(output_path, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line.startswith("#"):
                original = line
                line = _fix_format_mismatch(line)
                if line != original:
                    fixed_count += 1
            fout.write(line + "\n")

    logger.info("Cleaned VCF: fixed %d FORMAT/sample mismatches -> %s", fixed_count, output_path)
    return output_path


def preprocess_vcf(
    input_vcf: str | Path,
    phenotype_json: str | Path,
    output_vcf: str | Path,
    sample_id: str | None = None,
) -> PreprocessResult:
    """Preprocess a VCF by filling missing PGx positions with reference calls.

    Args:
        input_vcf: Path to input VCF or VCF.GZ.
        phenotype_json: Path to a PharmCAT phenotype.json.
        output_vcf: Path for the output VCF (uncompressed).
        sample_id: Override sample ID. If None, uses the one from the input VCF.

    Returns:
        PreprocessResult with statistics.
    """
    input_path = Path(input_vcf)
    output_path = Path(output_vcf)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load PGx positions
    gene_positions = load_pgx_positions(phenotype_json)

    # Collect all missing positions to fill (excluding HLA/mito genes)
    positions_to_fill: dict[tuple[str, int], PgxPosition] = {}
    for gene_symbol, gp in gene_positions.items():
        if gene_symbol in SKIP_GENES:
            continue
        for pos in gp.positions:
            key = (pos.chromosome, pos.position)
            positions_to_fill[key] = pos

    logger.info("Loaded %d PGx positions to potentially fill", len(positions_to_fill))

    # Read existing VCF
    opener = gzip.open if input_path.name.endswith(".gz") else open
    header_lines: list[str] = []
    data_lines: list[str] = []
    existing_positions: set[tuple[str, int]] = set()
    detected_sample_id = None

    with opener(input_path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
            elif line.startswith("#CHROM"):
                header_lines.append(line)
                cols = line.split("\t")
                if len(cols) > 9:
                    detected_sample_id = cols[9]
            else:
                line = _fix_format_mismatch(line)
                data_lines.append(line)
                parts = line.split("\t", 3)
                if len(parts) >= 2:
                    try:
                        existing_positions.add((parts[0], int(parts[1])))
                    except ValueError:
                        pass

    if sample_id is None:
        sample_id = detected_sample_id or "SAMPLE"

    original_count = len(data_lines)
    logger.info("Read %d existing VCF records (sample: %s)", original_count, sample_id)

    # Generate synthetic records for missing positions
    added_count = 0
    genes_filled: dict[str, int] = {}

    for (chrom, pos), pgx_pos in positions_to_fill.items():
        if (chrom, pos) in existing_positions:
            continue

        rsid = pgx_pos.rsid or "."
        ref = pgx_pos.reference_allele
        record = f"{chrom}\t{pos}\t{rsid}\t{ref}\t.\t50\tPASS\tPGX_REF_FILL\tGT\t0/0"
        data_lines.append(record)
        added_count += 1

        gene = pgx_pos.gene
        genes_filled[gene] = genes_filled.get(gene, 0) + 1

    logger.info("Added %d reference-fill records across %d genes", added_count, len(genes_filled))

    # Sort all data lines by chromosome + position
    def sort_key(line: str) -> tuple[int, int]:
        parts = line.split("\t", 3)
        chrom_num = CHROM_ORDER.get(parts[0], 99)
        try:
            pos_num = int(parts[1])
        except (ValueError, IndexError):
            pos_num = 0
        return (chrom_num, pos_num)

    data_lines.sort(key=sort_key)

    # Add PGX_REF_FILL INFO header before #CHROM line
    info_header = '##INFO=<ID=PGX_REF_FILL,Number=0,Type=Flag,Description="Position filled with reference genotype for PharmCAT PGx analysis">'

    with open(output_path, "w") as f:
        for hline in header_lines:
            if hline.startswith("#CHROM"):
                f.write(info_header + "\n")
            f.write(hline + "\n")
        for dline in data_lines:
            f.write(dline + "\n")

    total = original_count + added_count
    logger.info("Wrote %d total records to %s", total, output_path)

    return PreprocessResult(
        output_path=output_path,
        original_records=original_count,
        added_records=added_count,
        total_records=total,
        genes_filled=genes_filled,
    )
