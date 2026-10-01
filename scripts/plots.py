"""Draw ATAC coverage the way a genome browser would.

One row per cell type, all rows on the same y-axis so you can compare heights
honestly, with the genes in that window drawn underneath.

    from scripts.plots import plot_tracks
    plot_tracks(bigwigs, "ple", flank=20_000)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

from .genes import gene_window, genes_in_window


def coverage(bigwig: str | Path, chrom: str, start: int, end: int, bins: int = 800):
    """Average signal in each of `bins` equal slices of the window."""
    with pyBigWig.open(str(bigwig)) as bw:
        if chrom not in bw.chroms():
            raise KeyError(f"{Path(bigwig).name} has no chromosome '{chrom}'. "
                           f"It has e.g. {list(bw.chroms())[:5]}")
        values = bw.stats(chrom, start, end, nBins=bins, type="mean")
    return np.array([0.0 if v is None else v for v in values])


def read_peaks(bed: str | Path, chrom: str, start: int, end: int):
    """The peaks from a BED file that fall in this window, as (start, end) pairs."""
    out = []
    with open(bed) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if f[0] == chrom and int(f[2]) > start and int(f[1]) < end:
                out.append((int(f[1]), int(f[2])))
    return out


def plot_tracks(
    tracks: dict[str, str | Path],
    gene: str | None = None,
    flank: int = 20_000,
    region: tuple[str, int, int] | None = None,
    peaks: dict[str, str | Path] | None = None,
    highlight: tuple[int, int] | None = None,
    colors: dict[str, str] | None = None,
    title: str | None = None,
):
    """Plot one coverage track per cell type over a gene or a region.

    tracks     {"KC_g": "path/to/KC_g.bw", ...}
    gene       gene symbol, e.g. "ple" (case-sensitive) - or use `region`
    flank      how much space to show either side of the gene
    region     ("chr3L", 6700000, 6730000) instead of a gene name
    peaks      optional {"KC_g": "path/to/KC_g.bed"} to mark called peaks
    highlight  (start, end) to shade one candidate enhancer
    """
    if region is not None:
        chrom, start, end = region
    elif gene is not None:
        chrom, start, end = gene_window(gene, flank)
    else:
        raise ValueError("Give either a gene name or a region.")

    names = list(tracks)
    n = len(names)
    fig, axes = plt.subplots(
        n + 1, 1, figsize=(11, 1.15 * n + 1.4), sharex=True,
        gridspec_kw={"height_ratios": [1] * n + [0.55], "hspace": 0.15},
    )
    axes = np.atleast_1d(axes)
    x = np.linspace(start, end, 800)

    signals = {name: coverage(path, chrom, start, end) for name, path in tracks.items()}
    ymax = max((s.max() for s in signals.values()), default=1.0) * 1.1 or 1.0

    for ax, name in zip(axes, names):
        color = (colors or {}).get(name, "#3b6ea5")
        ax.fill_between(x, signals[name], color=color, linewidth=0)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=9)
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        if highlight:
            ax.axvspan(highlight[0], highlight[1], color="#f0a202", alpha=0.25, zorder=0)
        if peaks and name in peaks:
            for ps, pe in read_peaks(peaks[name], chrom, start, end):
                ax.axvspan(ps, pe, ymin=0.0, ymax=0.06, color="#222222", alpha=0.85)

    # gene models along the bottom
    gax = axes[-1]
    hits = genes_in_window(chrom, start, end)
    for i, (name, row) in enumerate(hits.iterrows()):
        y = -(i % 3)
        gax.plot([max(row.start, start), min(row.end, end)], [y, y], lw=3, color="#444444",
                 solid_capstyle="butt")
        gax.text((max(row.start, start) + min(row.end, end)) / 2, y + 0.25,
                 f"{name} {'▶' if row.strand == '+' else '◀'}",
                 ha="center", va="bottom", fontsize=8)
    gax.set_ylim(-3.2, 1.0)
    gax.set_yticks([])
    gax.set_ylabel("genes", rotation=0, ha="right", va="center", fontsize=9)
    for side in ("top", "right", "left"):
        gax.spines[side].set_visible(False)
    gax.set_xlim(start, end)
    # Plain positions in kb - matplotlib's default "+1.82e7" offset confuses people.
    gax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v / 1000:,.0f} kb")
    )
    gax.set_xlabel(f"{chrom}:{start:,}-{end:,}  (dm6)")
    if highlight:
        gax.axvspan(highlight[0], highlight[1], color="#f0a202", alpha=0.25, zorder=0)

    axes[0].set_title(title or (f"Chromatin accessibility around {gene}" if gene
                                else f"Chromatin accessibility at {chrom}:{start:,}-{end:,}"),
                      fontsize=11, pad=8)
    return fig
