"""Find known transcription factor binding sites in a sequence.

CREsted ships a motif database of 7,446 position weight matrices with a table
mapping them to fly transcription factors. This module pulls out the motifs for
a named TF and scans a sequence for matches, so a contribution score plot can be
annotated with the sites the model might be reading.

    from scripts.motifs import scan
    hits = scan(sequence, ["Mef2", "ey", "sr", "onecut"])

This is deliberately simple: a log-odds scan against a uniform background,
scoring both strands, keeping hits above a fraction of each motif's best
possible score. It is not a replacement for tf-modisco or TOMTOM, which discover
patterns rather than look for known ones.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# The motifs used in Figure 4 of Janssens et al. 2022, which dissects two
# enhancers. Fly names on the left; the database's motif id on the right.
# Two of these are human matrices, for reasons given at each entry: the fly
# database has no Mef2 motif at all, and its only Ey motif is the wrong one of
# that factor's two binding modes.
FIGURE4_MOTIFS = {
    "Mef2": "metacluster_179.4.taipale__MEF2A_DBD_KCTAWAAATAGM_repr",
    # Ey is Pax6, which binds DNA through a paired domain and a homeodomain.
    # The only fly `ey` entry in the database is a homeodomain-type site,
    # TAATCGAATTA, and the site Figure 4 annotates in enh47 is the paired-domain
    # one, TCACGCTTAATTAAGC. Scanning with the fly entry puts its best match at
    # position 381, which carries 0.68x the average contribution: a site the
    # model ignores. The human PAX6 matrix finds the real one at 332 at 3.14x.
    # Same reasoning as Mef2 above, and the paired domain is conserved.
    "ey": "metacluster_154.2.taipale_cyt_meth__PAX6_NYACGCNTSANYGNNYN_eDBD",
    "sr": "metacluster_131.7.flyfactorsurvey__sr_SOLEXA_5_FBgn0003499",
    # The database holds two Onecut motifs. Use the CGATCGATA one, which is the
    # cut-homeobox site shown in Figure 4. The flyfactorsurvey entry
    # (`flyfactorsurvey__onecut_SOLEXA_FBgn0028996`) is a TAATTAA-like variant
    # that matches homeodomain sites all over the genome and is misleading here.
    "onecut": "metacluster_165.4.nitta__onecut_TTAACT20NCAA_KG_NNCGATCGATANNN_m1_c2",
    "TfAP-2": "metacluster_125.2.nitta__AP-2_TGTCCG20NGA_KH_NGCCNNNNNGGCN_m0_c3",
    "fkh": "idmmpmm__fkh",
    "acj6": "metacluster_145.2.flyfactorsurvey__acj6_SOLEXA_5_FBgn0000028",
}

_cache: dict = {}


def _load_meme() -> dict[str, np.ndarray]:
    """Every motif in the CREsted database, as a position weight matrix."""
    if "meme" not in _cache:
        import crested

        path, _ = crested.get_motif_db()
        motifs, name, rows, reading = {}, None, [], False
        for line in open(path):
            line = line.strip()
            if line.startswith("MOTIF"):
                if name and rows:
                    motifs[name] = np.array(rows, dtype=float)
                name, rows, reading = line.split()[1], [], False
            elif line.startswith("letter-probability"):
                reading = True
            elif reading:
                parts = line.split()
                if len(parts) == 4:
                    rows.append([float(x) for x in parts])
                else:
                    reading = False
        if name and rows:
            motifs[name] = np.array(rows, dtype=float)
        _cache["meme"] = motifs
    return _cache["meme"]


def _one_hot(seq: str) -> np.ndarray:
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    out = np.zeros((len(seq), 4))
    for i, b in enumerate(seq.upper()):
        if b in idx:
            out[i, idx[b]] = 1
    return out


def _reverse_complement(pwm: np.ndarray) -> np.ndarray:
    return pwm[::-1, ::-1]


def _trim(pwm: np.ndarray, min_bits: float = 0.4) -> np.ndarray:
    """Drop uninformative positions from both ends of a matrix.

    Several database entries are padded with near-uniform columns, written as N
    in the motif name: the Onecut entry is `NNCGATCGATANNN`, so two of its
    fourteen columns on the left and three on the right match anything. Keeping
    them makes a reported site wider than the real binding site and puts the
    left edge several bases before the sequence actually starts matching, which
    is visible as a highlight box offset from the nucleotides it should mark.

    Information content per column is 2 + sum(p * log2 p): 2 bits for a fixed
    position, 0 for a uniform one.
    """
    p = np.clip(pwm, 1e-9, None)
    bits = 2 + (p * np.log2(p)).sum(axis=1)
    keep = np.where(bits >= min_bits)[0]
    if len(keep) == 0:
        return pwm
    return pwm[keep[0]: keep[-1] + 1]


def motifs_for(tf: str) -> list[str]:
    """Every motif in the database annotated to this fly transcription factor.

    The database ships a table mapping motifs to the factors that bind them,
    which lets any fly factor be scanned for by name rather than only the
    handful recorded in `FIGURE4_MOTIFS`. 672 motifs covering 295 factors can be
    resolved this way, including `gcm` and `repo`.

    The two files disagree on naming: the MEME file prefixes clustered motifs
    with `metacluster_<n>.<m>.`, the annotation table does not, so the prefix is
    stripped before matching. Without that only 247 of 1,368 resolve.
    """
    if "by_tf" not in _cache:
        import re

        import crested

        _, annotation = crested.get_motif_db()
        table = pd.read_csv(annotation, sep="\t").dropna(subset=["Fly_Direct_annot"])

        stripped = {}
        for key in _load_meme():
            stripped.setdefault(re.sub(r"^metacluster_[0-9.]+\.", "", key), key)

        index: dict[str, list[str]] = {}
        for row in table.itertuples():
            key = stripped.get(row.Motif_name)
            if key is not None:
                index.setdefault(row.Fly_Direct_annot, []).append(key)
        _cache["by_tf"] = index
    return _cache["by_tf"].get(tf, [])


def fimo_scan(sequence: str, tfs=None, p_value: float = 1e-3,
              motif_ids: dict | None = None) -> pd.DataFrame:
    """The same question asked with FIMO, the standard tool.

    `scan` above scores a match as a fraction of the range between the worst and
    best possible score for that matrix. That is easy to explain but it is not
    calibrated: the number means something different for every motif, and there
    is no way to say how often a score that high would arise by chance.

    FIMO, from the MEME suite, builds the null distribution of scores for each
    matrix and reports a p-value per hit, so one threshold means the same thing
    for every motif. It is used here through `memelite`, a fast reimplementation
    that `crested[motif]` already installs, so there is nothing extra to fetch.

    The two disagree on real data, and the disagreement is instructive rather
    than a bug. See section 6 of Practical 6.

    Returns one row per hit: tf, start, end, strand, score, p_value.
    """
    from memelite import fimo

    motif_ids = motif_ids or FIGURE4_MOTIFS
    tfs = list(tfs) if tfs is not None else list(motif_ids)
    db = _load_meme()

    # A factor may have several matrices in the database. All of them are
    # scanned and the name is recovered afterwards from the part before the bar.
    motifs = {}
    for tf in tfs:
        keys = [motif_ids[tf]] if tf in motif_ids else motifs_for(tf)
        if not keys:
            raise KeyError(
                f"No motif in the database is annotated to the fly factor '{tf}'. "
                f"Check the spelling against FlyBase; names are case-sensitive."
            )
        for n, key in enumerate(keys):
            motifs[f"{tf}|{n}"] = _trim(db[key]).T

    # memelite wants each matrix as (alphabet, width) and the sequence as
    # (n_sequences, alphabet, length), which is the transpose of what we use.
    encoded = _one_hot(sequence).T[None].astype("float32")

    found = fimo(motifs, encoded, threshold=p_value, dim=1)
    if not len(found) or not len(found[0]):
        return pd.DataFrame(columns=["tf", "start", "end", "strand", "score", "p_value"])

    hits = found[0].rename(columns={"motif_name": "tf", "p-value": "p_value"})
    hits["tf"] = hits["tf"].astype(str).str.split("|").str[0]
    return (hits[["tf", "start", "end", "strand", "score", "p_value"]]
            .sort_values("start").reset_index(drop=True))


def scan(sequence: str, tfs=None, threshold: float = 0.82,
         motif_ids: dict | None = None) -> pd.DataFrame:
    """Where in this sequence does each transcription factor have a site?

    threshold runs from 0 (the worst possible match) to 1 (a perfect one), so
    0.82 keeps reasonably good sites. Lower it if nothing is found.

    Returns one row per hit: tf, start, end, strand, sequence, score.
    """
    motif_ids = motif_ids or FIGURE4_MOTIFS
    tfs = list(tfs) if tfs is not None else list(motif_ids)
    db = _load_meme()
    hot = _one_hot(sequence)

    rows = []
    for tf in tfs:
        if tf not in motif_ids:
            raise KeyError(f"No motif id recorded for '{tf}'. Known: {list(motif_ids)}")
        pwm = db.get(motif_ids[tf])
        if pwm is None:
            raise KeyError(f"Motif {motif_ids[tf]} is not in the database.")
        pwm = _trim(pwm)

        # log-odds against a uniform background
        logodds = np.log2(np.clip(pwm, 1e-6, None) / 0.25)
        w = len(logodds)
        for strand, mat in (("+", logodds), ("-", _reverse_complement(logodds))):
            # Score on the full achievable range, not as a fraction of the best
            # possible score. A mismatch contributes a large negative log-odds,
            # so for a long matrix `score / best` stays low even for a good site
            # and long motifs are silently never reported. Rescaling so that 1 is
            # the best possible match and 0 the worst makes motifs of different
            # lengths comparable.
            best, worst = mat.max(axis=1).sum(), mat.min(axis=1).sum()
            span = best - worst
            for i in range(len(sequence) - w + 1):
                score = float((hot[i:i + w] * mat).sum())
                fraction = (score - worst) / span
                if fraction >= threshold:
                    rows.append({
                        "tf": tf, "start": i, "end": i + w, "strand": strand,
                        "sequence": sequence[i:i + w],
                        "score": round(fraction, 3),
                    })

    hits = pd.DataFrame(rows, columns=["tf", "start", "end", "strand",
                                       "sequence", "score"])
    if len(hits):
        hits = hits.sort_values(["start", "score"]).reset_index(drop=True)
    return hits


def best_per_tf(sequence: str, tfs=None, threshold: float = 0.82) -> pd.DataFrame:
    """The single best site for each transcription factor."""
    hits = scan(sequence, tfs, threshold=threshold)
    if not len(hits):
        return hits
    return (hits.sort_values("score", ascending=False)
                .groupby("tf", as_index=False)
                .first()
                .sort_values("start")
                .reset_index(drop=True))
