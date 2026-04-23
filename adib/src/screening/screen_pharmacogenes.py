"""
Pharmacogene screening tool.

Checks a VCF file for coverage at PharmCAT-required PGx positions
and classifies each gene by its rescue strategy.
"""

import gzip
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path

from src.screening.pgx_positions import GenePositions, load_pgx_positions

logger = logging.getLogger(__name__)

# Genes that need specialized tools beyond standard variant calling
HLA_GENES = {"HLA-A", "HLA-B"}
MITO_GENES = {"MT-RNR1"}
RESEARCH_GENES = {"CYP2D6"}


@dataclass
class GeneScreening:
    gene: str
    total_positions: int
    found_positions: int
    missing_positions: int
    category: str  # CALLED, RESCUABLE, NEEDS_RESEARCH_MODE, NEEDS_HLA_TYPING, NEEDS_CHRM_DATA, NO_POSITIONS


@dataclass
class ScreeningReport:
    sample_id: str | None
    vcf_path: str
    total_genes: int
    gene_results: list[GeneScreening] = field(default_factory=list)

    @property
    def rescuable_genes(self) -> list[GeneScreening]:
        return [g for g in self.gene_results if g.category == "RESCUABLE"]

    @property
    def called_genes(self) -> list[GeneScreening]:
        return [g for g in self.gene_results if g.category == "CALLED"]

    @property
    def uncallable_genes(self) -> list[GeneScreening]:
        return [g for g in self.gene_results
                if g.category in ("NEEDS_HLA_TYPING", "NEEDS_CHRM_DATA", "NO_POSITIONS")]


def _index_vcf_positions(vcf_path: str | Path) -> tuple[set[tuple[str, int]], str | None]:
    """Stream through VCF and collect all (chrom, pos) pairs.

    Returns:
        Tuple of (position_set, sample_id).
    """
    path = Path(vcf_path)
    opener = gzip.open if path.name.endswith(".gz") else open
    positions: set[tuple[str, int]] = set()
    sample_id = None

    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                if len(cols) > 9:
                    sample_id = cols[9]
                continue
            parts = line.split("\t", 3)
            if len(parts) >= 2:
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                    positions.add((chrom, pos))
                except ValueError:
                    pass

    return positions, sample_id


def screen_vcf(
    vcf_path: str | Path,
    phenotype_json_path: str | Path,
) -> ScreeningReport:
    """Screen a VCF file for PGx variant coverage.

    Args:
        vcf_path: Path to VCF or VCF.GZ file.
        phenotype_json_path: Path to a PharmCAT phenotype.json (from any prior run).

    Returns:
        ScreeningReport with per-gene coverage details.
    """
    logger.info("Loading PGx positions from %s", phenotype_json_path)
    gene_positions = load_pgx_positions(phenotype_json_path)

    logger.info("Indexing VCF positions from %s", vcf_path)
    vcf_positions, sample_id = _index_vcf_positions(vcf_path)
    logger.info("Found %d positions in VCF (sample: %s)", len(vcf_positions), sample_id)

    report = ScreeningReport(
        sample_id=sample_id,
        vcf_path=str(vcf_path),
        total_genes=len(gene_positions),
    )

    for gene_symbol, gp in sorted(gene_positions.items()):
        total = len(gp.positions)
        found = sum(
            1 for p in gp.positions
            if (p.chromosome, p.position) in vcf_positions
        )
        missing = total - found

        category = _classify_gene(gene_symbol, gp, found, missing)

        report.gene_results.append(GeneScreening(
            gene=gene_symbol,
            total_positions=total,
            found_positions=found,
            missing_positions=missing,
            category=category,
        ))

    return report


def _classify_gene(
    gene: str, gp: GenePositions, found: int, missing: int
) -> str:
    if gene in HLA_GENES:
        return "NEEDS_HLA_TYPING"
    if gene in MITO_GENES:
        return "NEEDS_CHRM_DATA"
    if gene in RESEARCH_GENES:
        return "NEEDS_RESEARCH_MODE"
    if gp.positions and missing == 0:
        return "CALLED"
    if gp.positions and missing > 0:
        return "RESCUABLE"
    return "NO_POSITIONS"


def print_screening_report(report: ScreeningReport) -> None:
    """Print a formatted screening report to stdout."""
    print(f"\n{'=' * 80}")
    print(f"  PGx Screening Report — {report.vcf_path}")
    print(f"  Sample: {report.sample_id or 'Unknown'}")
    print(f"{'=' * 80}\n")

    cat_order = {
        "CALLED": 0,
        "RESCUABLE": 1,
        "NEEDS_RESEARCH_MODE": 2,
        "NEEDS_HLA_TYPING": 3,
        "NEEDS_CHRM_DATA": 4,
        "NO_POSITIONS": 5,
    }
    results = sorted(report.gene_results, key=lambda g: (cat_order.get(g.category, 9), g.gene))

    print(f"  {'Gene':<15} {'Positions':>10} {'Found':>7} {'Missing':>9} {'Category'}")
    print(f"  {'-' * 15} {'-' * 10} {'-' * 7} {'-' * 9} {'-' * 25}")

    for g in results:
        pct = f"({g.found_positions}/{g.total_positions})" if g.total_positions > 0 else "(n/a)"
        print(f"  {g.gene:<15} {g.total_positions:>10} {g.found_positions:>7} {g.missing_positions:>9} {g.category} {pct}")

    called = len(report.called_genes)
    rescuable = len(report.rescuable_genes)
    uncallable = len(report.uncallable_genes)
    research = sum(1 for g in report.gene_results if g.category == "NEEDS_RESEARCH_MODE")

    print(f"\n  Summary:")
    print(f"    Already called:       {called}")
    print(f"    Rescuable (ref-fill): {rescuable}")
    print(f"    Research mode:        {research}")
    print(f"    Uncallable:           {uncallable}")
    print(f"    Total positions to fill: {sum(g.missing_positions for g in report.gene_results if g.category == 'RESCUABLE')}")
    print()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: python -m src.screening.screen_pharmacogenes <vcf_path> <phenotype_json>")
        sys.exit(1)

    logging.basicConfig(level=logging.INFO)
    report = screen_vcf(sys.argv[1], sys.argv[2])
    print_screening_report(report)
