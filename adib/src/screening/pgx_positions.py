"""
PGx position database.

Extracts the complete set of pharmacogenomic positions that PharmCAT
checks from a PharmCAT phenotype.json output file.  Each position records
the gene, chromosome, genomic coordinate, rsID, and reference allele.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PgxPosition:
    gene: str
    chromosome: str
    position: int
    rsid: str | None
    reference_allele: str


@dataclass
class GenePositions:
    gene: str
    chromosome: str | None
    positions: list[PgxPosition] = field(default_factory=list)
    call_source: str | None = None  # MATCHER, NONE, etc.


def load_pgx_positions(phenotype_json_path: str | Path) -> dict[str, GenePositions]:
    """Load all PGx positions from a PharmCAT phenotype.json.

    The phenotype.json contains every position PharmCAT checked for each gene,
    including positions where no genotype call was found (call=None).

    Returns:
        Dict mapping gene symbol to GenePositions.
    """
    path = Path(phenotype_json_path)
    with open(path) as f:
        data = json.load(f)

    gene_positions: dict[str, GenePositions] = {}

    for gene_symbol, report in data["geneReports"].items():
        chrom = report.get("chr")
        call_source = report.get("callSource")

        gp = GenePositions(
            gene=gene_symbol,
            chromosome=chrom,
            call_source=call_source,
        )

        for variant in report.get("variants", []):
            gp.positions.append(PgxPosition(
                gene=gene_symbol,
                chromosome=variant["chromosome"],
                position=variant["position"],
                rsid=variant.get("dbSnpId"),
                reference_allele=variant["referenceAllele"],
            ))

        gene_positions[gene_symbol] = gp

    return gene_positions
