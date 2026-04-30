"""
Reference-fill preprocessor (Adib).

Fills missing PGx positions with synthetic 0/0 reference calls so that
PharmCAT can make diplotype assignments. Valid for GIAB benchmark VCFs
where absence of a position means the sample matches the reference.

Position list comes from a PharmCAT phenotype.json (which records every
position PharmCAT checked, including positions with no genotype call).
"""

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path

from preprocessing.clean_vcf import fix_format_mismatch
from screening.pgx_positions import PgxPosition, load_pgx_positions

logger = logging.getLogger(__name__)

CHROM_ORDER: dict[str, int] = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})

# Genes that need specialized tools — never fill these.
SKIP_GENES = {"HLA-A", "HLA-B", "MT-RNR1"}


@dataclass
class RefFillResult:
    output_path: Path
    original_records: int
    added_records: int
    total_records: int
    genes_filled: dict[str, int] = field(default_factory=dict)


def reference_fill_vcf(
    input_vcf: str | Path,
    phenotype_json: str | Path,
    output_vcf: str | Path,
    sample_id: str | None = None,
) -> RefFillResult:
    """Fill missing PGx positions with reference-genotype records."""
    in_path = Path(input_vcf)
    out_path = Path(output_vcf)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    gene_positions = load_pgx_positions(phenotype_json)

    positions_to_fill: dict[tuple[str, int], PgxPosition] = {}
    for gene_symbol, gp in gene_positions.items():
        if gene_symbol in SKIP_GENES:
            continue
        for pos in gp.positions:
            positions_to_fill[(pos.chromosome, pos.position)] = pos

    logger.info("Loaded %d PGx positions to potentially fill", len(positions_to_fill))

    opener = gzip.open if in_path.name.endswith(".gz") else open
    header_lines: list[str] = []
    data_lines: list[str] = []
    existing_positions: set[tuple[str, int]] = set()
    detected_sample_id: str | None = None

    with opener(in_path, "rt", encoding="utf-8", errors="replace") as f:
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
                line = fix_format_mismatch(line)
                data_lines.append(line)
                parts = line.split("\t", 3)
                if len(parts) >= 2:
                    try:
                        existing_positions.add((parts[0], int(parts[1])))
                    except ValueError:
                        pass

    sid = sample_id or detected_sample_id or "SAMPLE"
    original = len(data_lines)
    logger.info("Read %d existing VCF records (sample: %s)", original, sid)

    added = 0
    genes_filled: dict[str, int] = {}
    for (chrom, pos), pgx in positions_to_fill.items():
        if (chrom, pos) in existing_positions:
            continue
        rsid = pgx.rsid or "."
        ref = pgx.reference_allele
        record = f"{chrom}\t{pos}\t{rsid}\t{ref}\t.\t50\tPASS\tPGX_REF_FILL\tGT\t0/0"
        data_lines.append(record)
        added += 1
        genes_filled[pgx.gene] = genes_filled.get(pgx.gene, 0) + 1

    logger.info("Added %d reference-fill records across %d genes", added, len(genes_filled))

    def _sort_key(line: str) -> tuple[int, int]:
        parts = line.split("\t", 3)
        chrom_num = CHROM_ORDER.get(parts[0], 99)
        try:
            pos_num = int(parts[1])
        except (ValueError, IndexError):
            pos_num = 0
        return (chrom_num, pos_num)

    data_lines.sort(key=_sort_key)

    info_header = (
        '##INFO=<ID=PGX_REF_FILL,Number=0,Type=Flag,'
        'Description="Position filled with reference genotype for PharmCAT PGx analysis">'
    )

    with open(out_path, "w", encoding="utf-8") as f:
        for hline in header_lines:
            if hline.startswith("#CHROM"):
                f.write(info_header + "\n")
            f.write(hline + "\n")
        for dline in data_lines:
            f.write(dline + "\n")

    total = original + added
    logger.info("Wrote %d total records to %s", total, out_path)

    return RefFillResult(
        output_path=out_path,
        original_records=original,
        added_records=added,
        total_records=total,
        genes_filled=genes_filled,
    )
