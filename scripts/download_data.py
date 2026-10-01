"""Download the public datasets used in the ISN School 2026 practicals.

Every file here is the *author-deposited processed file* from the original
source (the flybrain portal or GEO). Nothing is re-processed by us.

All URLs below were checked and returned HTTP 200 on 2026-09-22; the recorded
sizes are the Content-Length the server reported.

Usage from a notebook:

    from scripts.download_data import janssens_adult_bigwigs
    paths = janssens_adult_bigwigs(["KC_g", "T1", "Astrocyte_like"])

Offline fallback: if the environment variable ISN_USB points at the USB stick
(e.g. /Volumes/ISN2026), files are copied from there instead of downloaded.
"""

from __future__ import annotations

import os
import shutil
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------- locations

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "data" / "raw"
GENOME = REPO_ROOT / "genome"

FLYBRAIN = "https://flybrain.aertslab.org/downloads"
GEO_MERRILL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197760/suppl"
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips"


def _usb_root() -> Path | None:
    p = os.environ.get("ISN_USB")
    return Path(p) if p else None


# ------------------------------------------------------------------ fetcher


def fetch(url: str, dest: Path, expected_bytes: int | None = None) -> Path:
    """Get one file. Skips it if it is already there, tries the USB first."""
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size > 0:
        print(f"  already here: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest

    usb = _usb_root()
    if usb is not None:
        candidate = usb / "data" / "raw" / dest.name
        if candidate.exists():
            print(f"  copying from USB: {dest.name}")
            shutil.copy2(candidate, dest)
            return dest

    size = f" ({expected_bytes / 1e6:.0f} MB)" if expected_bytes else ""
    print(f"  downloading {dest.name}{size}\n    from {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh, length=1024 * 512)
    except urllib.error.HTTPError as err:
        tmp.unlink(missing_ok=True)
        if err.code == 404:
            raise FileNotFoundError(
                f"The server does not have {dest.name}.\n"
                f"  Tried: {url}\n"
                f"  If you asked for a named cell type, check the spelling: the "
                f"portal's names are exact."
            ) from None
        raise
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    tmp.replace(dest)
    print(f"    done: {dest.stat().st_size / 1e6:.1f} MB")
    return dest


# ------------------------------------- Janssens et al. 2022 (GSE163697)
# Nature 601:630. Portal: https://flybrain.aertslab.org
# Adult bigwigs are ~58 MB each; there are 82 adult cell types.

ADULT_BIGWIG = f"{FLYBRAIN}/bigwig/AdultCellTypes"
ADULT_PEAKS = f"{FLYBRAIN}/peaks/AdultCellTypes"
DEVEL_BIGWIG = f"{FLYBRAIN}/bigwig/EarlyDevelCellTypes"

# A small default set: mushroom body neurons, an optic lobe neuron, one glia.
PRACTICAL5_CELL_TYPES = ["KC_g", "KC_ab", "T1", "Astrocyte_like"]

# Careful: the portal does NOT use the same file name for a cell type's bigwig
# and its peaks. Checked 2026-09-22; these are the only disagreements, and
# T4/T5 subtypes have a bigwig but no peak file at all.
PEAK_NAME = {
    "KC_ab": "KC_a_b",            # alpha/beta Kenyon cells
    "KC_abPrime": "KC_a__b_",     # alpha'/beta' Kenyon cells
}
NO_PEAK_FILE = {"T4_a_b", "T4_c_d", "T5_a_b", "T5_c_d"}


def janssens_adult_bigwigs(cell_types=None) -> dict[str, Path]:
    """Per-cell-type accessibility tracks, adult brain. ~58 MB each."""
    cell_types = cell_types or PRACTICAL5_CELL_TYPES
    print("Janssens et al. 2022 - adult brain scATAC bigwigs")
    return {
        ct: fetch(f"{ADULT_BIGWIG}/{ct}.bw", RAW / "janssens_adult_bw" / f"{ct}.bw")
        for ct in cell_types
    }


def janssens_adult_peaks(cell_types=None) -> dict[str, Path]:
    """Per-cell-type peak calls, adult brain. ~0.4 MB each.

    Keys are the bigwig-style cell type names, so peaks and bigwigs line up
    even where the portal spells them differently.
    """
    cell_types = cell_types or PRACTICAL5_CELL_TYPES
    print("Janssens et al. 2022 - adult brain peaks per cell type")
    out = {}
    for ct in cell_types:
        if ct in NO_PEAK_FILE:
            print(f"  skipping {ct}: the portal has a bigwig but no peak file for it")
            continue
        remote = PEAK_NAME.get(ct, ct)
        out[ct] = fetch(
            f"{ADULT_PEAKS}/{remote}.bed", RAW / "janssens_adult_peaks" / f"{ct}.bed"
        )
    return out


def janssens_devel_bigwigs(cell_types) -> dict[str, Path]:
    """Early development (larva to 12h APF) tracks. Project 2."""
    print("Janssens et al. 2022 - early development bigwigs")
    return {
        ct: fetch(f"{DEVEL_BIGWIG}/{ct}.bw", RAW / "janssens_devel_bw" / f"{ct}.bw")
        for ct in cell_types
    }


def janssens_brain_peaks() -> Path:
    """Whole-brain peak set, resized to max 500 bp. 4.8 MB."""
    print("Janssens et al. 2022 - whole-brain peak set")
    return fetch(
        f"{FLYBRAIN}/other/BrainPeaks_ResizedToMax500.bed",
        RAW / "BrainPeaks_ResizedToMax500.bed",
        4_774_826,
    )


def janssens_dars() -> Path:
    """Differentially accessible regions per cell type. 18 MB zip."""
    print("Janssens et al. 2022 - differentially accessible regions")
    return fetch(f"{FLYBRAIN}/other/tbl_DARs.txt.zip", RAW / "tbl_DARs.txt.zip", 18_105_749)


def janssens_tested_enhancers() -> Path:
    """The 54 enhancers Janssens et al. tested in vivo with reporter lines.

    BED with a name per region (enh14, enh47, ...). These are worth more than an
    arbitrary accessible region, because their activity was measured in flies
    rather than inferred. Figure 4 of the paper dissects two of them:

      enh47  chr2L:20040501-20041072  near sNPF,   gamma Kenyon cells, Ey + Mef2
      enh44  chr3R:21945795-21946260  near Eip93F, T4 neurons,         Fkh + TfAP-2

    Widths run from 299 to 1,731 bp, so centre them before giving them to a
    model that expects 500 bp.
    """
    print("Janssens et al. 2022 - enhancers tested in vivo")
    return fetch(f"{FLYBRAIN}/other/enhancersTested.bed", RAW / "enhancersTested.bed")


def janssens_celltype_info() -> Path:
    """Cell type names, abbreviations and cell counts. 7 KB."""
    return fetch(f"{FLYBRAIN}/cellTypeInfo.tsv", RAW / "cellTypeInfo.tsv", 7_241)


# --------------------------------------------- expression (scanpy) looms


def larval_brain_loom() -> Path:
    """Ravenscroft et al. 2019 larval brain scRNA-seq. 21 MB."""
    print("Ravenscroft et al. 2019 - larval brain scRNA-seq")
    return fetch(
        f"{FLYBRAIN}/loom/Ravenscroft_et_al_2019_LarvalBrain.loom",
        RAW / "Ravenscroft_2019_LarvalBrain.loom",
        20_945_628,
    )


def fca_head() -> Path:
    """Fly Cell Atlas, 10x head, stringent annotation. 2.5 GB.

    Li, Janssens et al. 2022, Science 375:eabk2432. The community reference
    atlas for Drosophila. The head sample covers brain, eye, cuticle and some
    muscle and fat; `scripts/prepare_data.py` reduces it to the brain cell types.

    The atlas is served from a Nextcloud share, which does not report a file
    size up front, so this download shows no total. It is 2.5 GB.
    """
    print("Fly Cell Atlas - 10x head, stringent (2.5 GB)")
    print("  Tissue list and other formats: https://flycellatlas.org")
    return fetch(
        "https://cloud.flycellatlas.org/index.php/s/LAEybPc2HZnpzKs/download",
        RAW / "fca_head_stringent.h5ad",
    )


def davie_adult_loom() -> Path:
    """Davie et al. 2018 adult brain atlas, 118,687 cells. 519 MB - big!

    Only call this if you really need the full adult atlas. On an 8 GB laptop
    open it with backed mode or subset it first.
    """
    print("Davie et al. 2018 - adult brain atlas (519 MB)")
    return fetch(
        f"{FLYBRAIN}/loom/Davie_Janssens_Koldere_et_al_2018_AdultBrain.loom",
        RAW / "Davie_2018_AdultBrain.loom",
        519_015_198,
    )


# ------------------------------------- Merrill et al. 2022 (GSE197760)
# Sorted dopaminergic (TH) vs GABAergic neurons, bulk ATAC. Project 3.

MERRILL_FILES = {
    "TH_bw": ("GSE197760_TH-merged.bw", 15_414_605),
    "GABA_bw": ("GSE197760_GABA-merged.bw", 16_482_910),
    "differential": ("GSE197760_differential.csv.gz", 4_951_817),
    "peaks_tar": ("GSE197760_RAW.tar", 7_536_640),
}


def merrill_atac(which=None) -> dict[str, Path]:
    """TH vs GABA bigwigs, differential table and per-sample narrowPeaks."""
    which = which or list(MERRILL_FILES)
    print("Merrill et al. 2022 - TH vs GABA neuron bulk ATAC (GEO GSE197760)")
    out = {}
    for key in which:
        name, size = MERRILL_FILES[key]
        out[key] = fetch(f"{GEO_MERRILL}/{name}", RAW / "merrill" / name, size)
    return out


# ------------------------------------------------------------------ genome


def dm6_genome() -> Path:
    """dm6 FASTA from UCSC, unzipped. ~44 MB gzipped, ~170 MB unzipped."""
    import gzip

    fa = GENOME / "dm6.fa"
    if fa.exists():
        print(f"  already here: dm6.fa ({fa.stat().st_size / 1e6:.0f} MB)")
        return fa
    gz = fetch(f"{UCSC}/dm6.fa.gz", GENOME / "dm6.fa.gz")
    print("  unzipping dm6.fa.gz")
    with gzip.open(gz, "rb") as src, open(fa, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 512)
    gz.unlink()
    print(f"    done: {fa.stat().st_size / 1e6:.0f} MB")
    return fa


# --------------------------- Dickmanken et al. 2026 embryo scATAC (GSE293575)
# Nature Communications, "Evaluating single-cell ATAC-seq atlasing technologies
# using sequence-to-function modeling". DOI 10.1038/s41467-026-68742-4.
# Drosophila embryo, 16-20 h after egg laying, one stage.
# Matches the CREsted model "Embryo10x".

EMBRYO_BW = ("https://ucsctracks.aertslab.org/papers/hydrop_v2_paper/dm6/bw/"
             "drosophila_embryo_10x_v2_only")

# The model's class names and the bigwig file names disagree for three cell
# types, two of them because the file names are misspelled. Checked 2026-09-23.
# Left: the model output name. Right: the file on the server.
EMBRYO_FILE_NAME = {
    "Malpighian_Tubule": "Malpighian_tubule",
    "Pharynx": "Pharnyx",
    "muscle_attachment_Stripe": "muscle_attachement_Stripe",
}

# Sizes on the server, in MB, read 2026-09-23. These tracks are much larger
# than the Janssens ones, so choose a few rather than taking all twenty.
EMBRYO_SIZES_MB = {
    "Midgut_acidification": 49, "Visceral_muscles": 56, "Yolk": 62,
    "Head_Ectoderm": 65, "Malpighian_Tubule": 71, "Salivary_gland": 75,
    "Hemocytes": 78, "PNS_sens_neurons": 83, "Tracheal_system": 95,
    "Glia": 97, "muscle_attachment_Stripe": 114, "Primordium_all": 120,
    "Hindgut": 126, "Neuroblasts": 146, "Somatic_muscles": 148,
    "Pharynx": 152, "Fat_body": 203, "Epidermis": 263, "Midgut": 265,
    "Neuronal": 343,
}


def dickmanken_embryo_bigwigs(cell_types=("Neuronal", "Glia")) -> dict[str, Path]:
    """Per-cell-type embryo accessibility tracks. Large: 49 to 343 MB each.

    Keys are the Embryo10x model's class names, so predictions and tracks line
    up even where the server spells a name differently.

    Note this is a single timepoint (16-20 h AEL), not a developmental series.
    For a series within one dataset, use `janssens_devel_bigwigs`.
    """
    print("Dickmanken et al. 2026 - embryo scATAC bigwigs (16-20 h AEL)")
    total = sum(EMBRYO_SIZES_MB.get(ct, 0) for ct in cell_types)
    print(f"  {len(cell_types)} tracks, about {total} MB in total")
    out = {}
    for ct in cell_types:
        remote = EMBRYO_FILE_NAME.get(ct, ct)
        out[ct] = fetch(
            f"{EMBRYO_BW}/{remote}.bw",
            RAW / "dickmanken_embryo_bw" / f"{ct}.bw",
            EMBRYO_SIZES_MB.get(ct, 0) * 1_000_000 or None,
        )
    return out


if __name__ == "__main__":
    janssens_celltype_info()
    janssens_adult_peaks()
    janssens_brain_peaks()
    janssens_adult_bigwigs()
    merrill_atac()
    dm6_genome()
