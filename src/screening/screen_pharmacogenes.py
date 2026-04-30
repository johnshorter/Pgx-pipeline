"""
Pharmacogene coverage screening.

Reports per-gene coverage of a VCF against the PharmCAT-required PGx
positions (loaded from a phenotype.json) and classifies each gene by
its rescue strategy:
    CALLED              — all positions present
    RESCUABLE           — partially missing, fixable by reference-fill
    NEEDS_RESEARCH_MODE — CYP2D6 (best-effort via -research cyp2d6)
    NEEDS_HLA_TYPING    — HLA-A / HLA-B (specialized typing required)
    NEEDS_CHRM_DATA     — MT-RNR1 (mitochondrial; absent from chr1-22 VCFs)
    NO_POSITIONS        — gene has no positions in the panel (rare)
"""

import gzip
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path

from screening.pgx_positions import GenePositions, load_pgx_positions

logger = logging.getLogger(__name__)

HLA_GENES = {"HLA-A", "HLA-B"}
MITO_GENES = {"MT-RNR1"}
RESEARCH_GENES = {"CYP2D6"}


@dataclass
class GeneScreening:
    gene: str
    total_positions: int
    found_positions: int
    missing_positions: int
    category: str


@dataclass
class ScreeningReport:
    sample_id: str | None
    vcf_path: str
    total_genes: int
    gene_results: list[GeneScreening] = field(default_factory=list)

    @property
    def called_genes(self) -> list[GeneScreening]:
        return [g for g in self.gene_results if g.category == "CALLED"]

    @property
    def rescuable_genes(self) -> list[GeneScreening]:
        return [g for g in self.gene_results if g.category == "RESCUABLE"]

    @property
    def uncallable_genes(self) -> list[GeneScreening]:
        return [
            g for g in self.gene_results
            if g.category in ("NEEDS_HLA_TYPING", "NEEDS_CHRM_DATA", "NO_POSITIONS")
        ]


def _index_vcf_positions(vcf_path: str | Path) -> tuple[set[tuple[str, int]], str | None]:
    path = Path(vcf_path)
    opener = gzip.open if path.name.endswith(".gz") else open
    positions: set[tuple[str, int]] = set()
    sample_id: str | None = None

    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 9:
                    sample_id = cols[9]
                continue
            parts = line.split("\t", 3)
            if len(parts) >= 2:
                try:
                    positions.add((parts[0], int(parts[1])))
                except ValueError:
                    pass
    return positions, sample_id


def screen_vcf(
    vcf_path: str | Path,
    phenotype_json_path: str | Path,
) -> ScreeningReport:
    """Screen a VCF for PGx variant coverage."""
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
        report.gene_results.append(GeneScreening(
            gene=gene_symbol,
            total_positions=total,
            found_positions=found,
            missing_positions=missing,
            category=_classify(gene_symbol, gp, missing),
        ))
    return report


def _classify(gene: str, gp: GenePositions, missing: int) -> str:
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
    print(f"\n{'=' * 80}")
    print(f"  PGx Screening Report — {report.vcf_path}")
    print(f"  Sample: {report.sample_id or 'Unknown'}")
    print(f"{'=' * 80}\n")

    cat_order = {
        "CALLED": 0, "RESCUABLE": 1, "NEEDS_RESEARCH_MODE": 2,
        "NEEDS_HLA_TYPING": 3, "NEEDS_CHRM_DATA": 4, "NO_POSITIONS": 5,
    }
    rows = sorted(report.gene_results, key=lambda g: (cat_order.get(g.category, 9), g.gene))

    print(f"  {'Gene':<15} {'Positions':>10} {'Found':>7} {'Missing':>9} {'Category'}")
    print(f"  {'-' * 15} {'-' * 10} {'-' * 7} {'-' * 9} {'-' * 25}")
    for g in rows:
        print(f"  {g.gene:<15} {g.total_positions:>10} {g.found_positions:>7} "
              f"{g.missing_positions:>9} {g.category}")

    print(f"\n  Summary:")
    print(f"    Already called:       {len(report.called_genes)}")
    print(f"    Rescuable (ref-fill): {len(report.rescuable_genes)}")
    print(f"    Research mode:        {sum(1 for g in report.gene_results if g.category == 'NEEDS_RESEARCH_MODE')}")
    print(f"    Uncallable:           {len(report.uncallable_genes)}")
    print(f"    Total positions to fill: "
          f"{sum(g.missing_positions for g in report.rescuable_genes)}")
    print()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python -m screening.screen_pharmacogenes <vcf_path> <phenotype_json>")
        sys.exit(1)
    logging.basicConfig(level=logging.INFO)
    print_screening_report(screen_vcf(sys.argv[1], sys.argv[2]))
