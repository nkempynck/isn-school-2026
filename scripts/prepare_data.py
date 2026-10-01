"""Turn the downloaded Fly Cell Atlas head file into a small brain-only h5ad.

Why this exists
---------------
The published FCA head file is 2.5 GB. Most of that is material the practicals
do not use: twenty stored Leiden clusterings with their marker gene tables, a
PCA, and a highly-variable-gene matrix that omits most genes. Reading it with
`sc.read_h5ad` takes about 50 seconds and 2.5 GB, and subsetting it afterwards
with `.copy()` needs more memory than an 8 GB laptop has. It was killed by the
operating system on a 17 GB machine during testing.

So we do the reduction once, here, reading the file with h5py rather than
anndata so that nothing is held in memory that we are about to discard. The
result is a normal h5ad of about 350 MB that `scanpy` opens in a few seconds,
and every analysis in the practicals is then plain scanpy on that file.

Nothing is recomputed. The values, the UMAP coordinates and the annotations are
the authors' own; we are only selecting rows and dropping unused fields.

Note on normalisation, which differs between the two files:

- `fca_brain.h5ad` holds the atlas `raw` matrix, which is **already**
  `normalize_total(target_sum=1e4)` followed by `log1p`. Verified: taking
  `expm1` gives exactly 10,000 per cell, and dividing by the smallest nonzero
  value per cell returns whole numbers. Do not normalise it again.
- `larval_brain.h5ad` holds **raw integer counts**, per-cell totals between
  2,015 and 43,644. It does need normalising before use.

Run once:

    python scripts/prepare_data.py
"""

from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "data" / "raw" / "fca_head_stringent.h5ad"
DEST = ROOT / "data" / "processed" / "fca_brain.h5ad"

# Building the small files from the 2.5 GB original takes a few minutes and
# 2.5 GB of memory. That is fine once on a laptop, but wasteful when thirty
# people each do it at the start of a one hour practical, and on Colab it would
# be repeated every time a session times out. So the prepared files are also
# kept on the lab's resource server and fetched from there when available.
#
# Upload fca_brain.h5ad and larval_brain.h5ad to this directory before the
# school. Point ISN_PREPARED_URL somewhere else, or set it empty, to build
# everything from the authors' original files instead.
PREPARED_URL = os.environ.get(
    "ISN_PREPARED_URL", "https://resources.aertslab.org/isn_school_2026"
).rstrip("/")


def _try_mirror(dest: Path) -> bool:
    """Fetch an already prepared file, if a mirror is configured."""
    if not PREPARED_URL:
        return False
    from .download_data import fetch

    try:
        fetch(f"{PREPARED_URL}/{dest.name}", dest)
        return True
    except Exception as err:                       # fall back to building it
        print(f"  mirror unavailable ({err}), building from the original")
        return False

# The head sample includes eye, cuticle, muscle and fat. These are the two
# broad classes that make up the brain.
BRAIN_CLASSES = ["neuron", "glial cell"]


def _decode(a) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def make_fca_brain_h5ad(overwrite: bool = False) -> Path:
    if DEST.exists() and not overwrite:
        print(f"  already prepared: {DEST.name} ({DEST.stat().st_size / 1e6:.0f} MB)")
        return DEST
    if _try_mirror(DEST):
        return DEST
    if not SOURCE.exists():
        raise FileNotFoundError(
            f"{SOURCE} not found.\n"
            f"  Run: from scripts.download_data import fca_head; fca_head()"
        )

    print(f"Preparing {DEST.name} from {SOURCE.name}")
    with h5py.File(SOURCE, "r") as f:
        cat = lambda c: _decode(f[f"uns/{c}_categories"][:])[f["obs"][c][:]]
        annotation = cat("annotation")
        broad = cat("annotation_broad")
        broad_filled = cat("annotation_broad_extrapolated")

        keep = np.isin(broad, BRAIN_CLASSES)
        print(f"  keeping {keep.sum():,} brain cells of {len(keep):,}")

        symbols = _decode(f["raw.cat/Symbol_categories"][:])[f["raw.var"]["Symbol"][:]]
        is_tf = f["raw.var"]["flybase_mg__transcription_factors"][:]

        # raw.X holds all 13,056 genes and is stored gene by gene, so we build
        # it in that orientation and transpose once.
        counts = sp.csr_matrix(
            (f["raw.X/data"][:], f["raw.X/indices"][:], f["raw.X/indptr"][:]),
            shape=(len(symbols), len(annotation)),
        )
        X = counts.T.tocsr()[keep].copy()
        del counts

        obs = pd.DataFrame({
            "annotation": pd.Categorical(annotation[keep]),
            "annotation_broad": pd.Categorical(broad[keep]),
            "annotation_broad_extrapolated": pd.Categorical(broad_filled[keep]),
            "sex": pd.Categorical(cat("sex")[keep]),
            "n_genes": f["obs"]["n_genes"][:][keep],
        })
        obs.index = _decode(f["obs"]["index"][:])[keep]

        obsm = {
            "X_umap": f["obsm/X_umap"][:][keep].astype("float32"),
            "X_tsne": f["obsm/X_tsne"][:][keep].astype("float32"),
        }

    var = pd.DataFrame({"is_TF": is_tf.astype(bool)}, index=pd.Index(symbols, name=None))
    adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm)
    adata.var_names_make_unique()
    adata.uns["source"] = (
        "Fly Cell Atlas, 10x head, stringent annotation (Li, Janssens et al. 2022). "
        "Brain cell types only. https://flycellatlas.org"
    )

    DEST.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(DEST, compression="gzip")
    print(f"  wrote {DEST} ({DEST.stat().st_size / 1e6:.0f} MB), shape {adata.shape}")
    return DEST


LARVAL_SOURCE = ROOT / "data" / "raw" / "Ravenscroft_2019_LarvalBrain.loom"
LARVAL_DEST = ROOT / "data" / "processed" / "larval_brain.h5ad"


def make_larval_h5ad(overwrite: bool = False) -> Path:
    """Convert the larval brain loom to h5ad.

    `sc.read_loom` cannot open this file. It is a SCope-style loom, and it
    stores `Embedding` and `Clusterings` as compound dtypes, which the Loom 2.0
    specification does not allow, so loompy rejects it during validation.
    Reading it with h5py and writing an h5ad avoids the problem and means the
    practicals only ever call `sc.read_h5ad`.
    """
    if LARVAL_DEST.exists() and not overwrite:
        print(f"  already prepared: {LARVAL_DEST.name}")
        return LARVAL_DEST
    if _try_mirror(LARVAL_DEST):
        return LARVAL_DEST
    if not LARVAL_SOURCE.exists():
        raise FileNotFoundError(
            f"{LARVAL_SOURCE} not found.\n"
            f"  Run: from scripts.download_data import larval_brain_loom; larval_brain_loom()"
        )

    print(f"Preparing {LARVAL_DEST.name} from {LARVAL_SOURCE.name}")
    with h5py.File(LARVAL_SOURCE, "r") as f:
        genes = _decode(f["row_attrs/Gene"][:])
        counts = sp.csr_matrix(f["matrix"][:])        # genes x cells, 9853 x 5056
        emb = f["col_attrs/Embedding"][:]
        annotation = _decode(f["col_attrs/annotation"][:])
        obs = pd.DataFrame({
            "annotation": pd.Categorical(annotation),
            # Clusters the authors left as bare numbers were never named.
            "is_annotated": [not str(v).isdigit() for v in annotation],
            "genotype": pd.Categorical(_decode(f["col_attrs/genotype"][:])),
            "n_genes": f["col_attrs/nFeature_RNA"][:],
        })
        obs.index = _decode(f["col_attrs/CellID"][:])
        obsm = {"X_umap": np.column_stack([emb["_X"], emb["_Y"]]).astype("float32")}

    adata = ad.AnnData(X=counts.T.tocsr(), obs=obs,
                       var=pd.DataFrame(index=pd.Index(genes)), obsm=obsm)
    adata.var_names_make_unique()
    adata.uns["source"] = (
        "Ravenscroft et al. 2019, larval brain. https://flybrain.aertslab.org"
    )
    LARVAL_DEST.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(LARVAL_DEST, compression="gzip")
    print(f"  wrote {LARVAL_DEST} ({LARVAL_DEST.stat().st_size / 1e6:.0f} MB), "
          f"shape {adata.shape}")
    return LARVAL_DEST


if __name__ == "__main__":
    make_fca_brain_h5ad()
    make_larval_h5ad()
