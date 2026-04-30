"""
PGx position database — extracts the full set of pharmacogenomic positions
PharmCAT checks, from a phenotype.json file.

Used by the reference-fill preprocessor to know exactly which positions
to synthesize 0/0 records for.
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
    call_source: str | None = None


def load_pgx_positions(phenotype_json_path: str | Path) -> dict[str, GenePositions]:
    """Load all PGx positions from a PharmCAT phenotype.json.

    The phenotype.json records every position PharmCAT checked for each
    gene, including positions where no genotype call was found.
    """
    path = Path(phenotype_json_path)
    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    out: dict[str, GenePositions] = {}
    for gene_symbol, report in (data.get("geneReports", {}) or {}).items():
        gp = GenePositions(
            gene=gene_symbol,
            chromosome=report.get("chr"),
            call_source=report.get("callSource"),
        )
        for variant in report.get("variants", []) or []:
            gp.positions.append(PgxPosition(
                gene=gene_symbol,
                chromosome=variant["chromosome"],
                position=variant["position"],
                rsid=variant.get("dbSnpId"),
                reference_allele=variant["referenceAllele"],
            ))
        out[gene_symbol] = gp
    return out
