"""Look up where a fly gene sits in the dm6 genome.

`resources/dm6_genes.tsv` was built from UCSC's ncbiRefSeqCurated table for
dm6 (downloaded 2026-09-22), one row per gene symbol, spanning all its
transcripts. 17,184 genes.

    from scripts.genes import gene_window
    chrom, start, end = gene_window("ple", flank=10_000)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

TABLE = Path(__file__).resolve().parent.parent / "resources" / "dm6_genes.tsv"

_genes: pd.DataFrame | None = None


def gene_table() -> pd.DataFrame:
    """The whole table, loaded once."""
    global _genes
    if _genes is None:
        _genes = pd.read_csv(TABLE, sep="\t").set_index("gene")
    return _genes


def gene_window(gene: str, flank: int = 10_000) -> tuple[str, int, int]:
    """Where is this gene, plus some space either side?

    The extra space (flank) matters: enhancers often sit outside the gene
    itself, sometimes tens of kilobases away.
    """
    t = gene_table()
    if gene not in t.index:
        hint = [g for g in t.index if g.lower() == gene.lower()]
        raise KeyError(
            f"No gene called '{gene}' in dm6."
            + (f" Did you mean '{hint[0]}'? (names are case-sensitive)" if hint else "")
        )
    row = t.loc[gene]
    return str(row.chrom), max(0, int(row.start) - flank), int(row.end) + flank


def region_string(gene: str, flank: int = 10_000) -> str:
    """Same thing as 'chr3L:6703355-6729525', the way genome browsers write it."""
    chrom, start, end = gene_window(gene, flank)
    return f"{chrom}:{start}-{end}"


def genes_in_window(chrom: str, start: int, end: int) -> pd.DataFrame:
    """Which genes overlap this window? Useful after you spot a peak."""
    t = gene_table()
    hit = t[(t.chrom == chrom) & (t.end > start) & (t.start < end)]
    return hit.sort_values("start")
