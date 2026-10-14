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

# Some files are also kept on the lab's own server, because the original host is
# slow. UCSC serves dm6 at roughly 1.7 MB/s against 8 MB/s here, and the genome
# is the single slowest download in the practicals. Anything missing from the
# mirror falls back to the original source automatically.
MIRROR = os.environ.get(
    "ISN_MIRROR", "https://resources.aertslab.org/workshops/isn_course_oct2026"
).rstrip("/")


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
            declared = r.headers.get("Content-Length")
        # A connection that drops midway still leaves a readable file, so the
        # only way to notice is to compare what arrived against what the server
        # said it was sending.
        if declared is not None and tmp.stat().st_size != int(declared):
            got, want = tmp.stat().st_size, int(declared)
            tmp.unlink(missing_ok=True)
            raise IOError(
                f"{dest.name} arrived incomplete: {got:,} of {want:,} bytes. "
                f"The connection probably dropped. Run the cell again."
            )
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


def _streaming_works(url: str) -> bool:
    """Can pyBigWig read this file over the network, without downloading it?

    The bigwig format is indexed, so a reader can ask the server for just the
    bytes covering one locus instead of fetching the whole file. That needs the
    server to answer HTTP range requests, which the flybrain portal does, and it
    needs nothing in between to interfere. Some institutional proxies strip
    range requests and return the whole file, which breaks the reader.

    Checked once per session rather than per file.
    """
    try:
        import pyBigWig
    except ImportError:
        return False

    # When this fails, the C library underneath pyBigWig writes its own errors
    # straight to file descriptor 2, where Python cannot intercept them. They
    # look alarming and are followed by a perfectly good fallback, so they are
    # silenced for the duration of the check.
    saved = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 2)
        bw = pyBigWig.open(url)
        ok = bool(bw.chroms())
        bw.close()
        return ok
    except Exception:
        return False
    finally:
        os.dup2(saved, 2)
        os.close(saved)
        os.close(devnull)


def _resolve_tracks(entries: dict, stream: bool, sizes: dict | None = None) -> dict:
    """Decide, per track, whether to read it over the network or fetch it.

    `entries` maps a name to (url, local path). The order of preference is:

    1. a copy already on this machine, so nothing is repeated
    2. the USB stick, if ISN_USB is set
    3. streaming, which for a bigwig means a few range requests per locus
       instead of tens or hundreds of megabytes
    4. downloading, if streaming does not work here

    Whether streaming works is checked once per call rather than per file.
    """
    have = {n: p for n, (_, p) in entries.items() if p.exists() and p.stat().st_size > 0}
    missing = [n for n in entries if n not in have]

    use_stream = False
    if missing and stream and _usb_root() is None:
        use_stream = _streaming_works(entries[missing[0]][0])
        if use_stream:
            print(f"  reading {len(missing)} track(s) over the network, not downloading")
        else:
            mb = sum((sizes or {}).get(n, 0) for n in missing)
            extra = f", about {mb} MB" if mb else ""
            print(f"  cannot read over the network here, downloading instead{extra}")

    out = {}
    for name, (url, path) in entries.items():
        if name in have:
            print(f"  already here: {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
            out[name] = path
        elif use_stream:
            out[name] = url
        else:
            out[name] = fetch(url, path,
                              (sizes or {}).get(name, 0) * 1_000_000 or None)
    return out


def janssens_adult_bigwigs(cell_types=None, stream: bool = True):
    """Per-cell-type accessibility tracks, adult brain.

    Returns a dict of cell type -> URL or local path, either of which can be
    handed straight to `pyBigWig.open` or to `scripts.plots.plot_tracks`.

    Each file is about 58 MB, and the practical only ever looks at a few loci in
    them, so by default they are read over the network rather than downloaded.
    That turns 182 MB of downloading into a few hundred kilobytes of range
    requests. See `_resolve_tracks` for the order of preference.

    Pass `stream=False` to force the download, which is what you want if you are
    about to go offline.
    """
    cell_types = cell_types or PRACTICAL5_CELL_TYPES
    print("Janssens et al. 2022 - adult brain scATAC bigwigs")
    return _resolve_tracks(
        {ct: (f"{ADULT_BIGWIG}/{ct}.bw", RAW / "janssens_adult_bw" / f"{ct}.bw")
         for ct in cell_types},
        stream,
    )


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


def merrill_bigwigs(stream: bool = True):
    """The two accessibility tracks, dopaminergic and GABAergic.

    Returns a dict ready for `scripts.plots.plot_tracks`. Like the fly brain
    tracks these are read over the network by default rather than downloaded.

    The two files are already normalised against each other: their `sumData`
    totals are 1.9687e9 and 1.9734e9, 0.24 per cent apart, so peak heights can
    be compared directly with no scaling. That is not true of every pair of
    bigwigs and is worth checking before trusting a side-by-side plot.
    """
    print("Merrill et al. 2022 - TH vs GABA neuron bulk ATAC (GEO GSE197760)")
    return _resolve_tracks(
        {label: (f"{GEO_MERRILL}/{MERRILL_FILES[key][0]}",
                 RAW / "merrill" / MERRILL_FILES[key][0])
         for key, label in (("TH_bw", "TH (dopaminergic)"), ("GABA_bw", "GABA"))},
        stream,
        {"TH (dopaminergic)": 15, "GABA": 16},
    )


def merrill_differential() -> Path:
    """The authors' own differential accessibility table, 5 MB gzipped.

    80,590 regions, each already tested for a difference between the two cell
    types with DESeq2 and annotated with its nearest gene by ChIPseeker. Read it
    with `pd.read_csv(path, skiprows=1, low_memory=False)`: there is a one-line
    comment above the header. Positive `log2FoldChange` means more open in
    dopaminergic neurons.
    """
    name, size = MERRILL_FILES["differential"]
    print("Merrill et al. 2022 - differential accessibility table")
    return fetch(f"{GEO_MERRILL}/{name}", RAW / "merrill" / name, size)


# --------------------------------- Pech et al. 2025 (eLife 98348, GSE228843)
# Five fly models of familial Parkinsonism, whole-brain scRNA-seq of young
# flies, plus human snRNA-seq. The group project uses the authors' own result
# tables rather than the raw data: the 24 per-sample h5ads on GEO carry no cell
# type annotations, so using them would mean clustering and annotating 75,000
# cells, which is a week of work rather than an evening.

ELIFE_PECH = "https://cdn.elifesciences.org/articles/98348"

PECH_FILES = {
    # Supplementary file 3: how many cells of each cluster were captured, per
    # genotype. This is what tells you whether a result rests on 90 cells or
    # 900, and there are three separate control genotypes.
    "counts": ("elife-98348-supp3-v1.xlsx", 65_919),
    # Supplementary file 5: how deregulated each cell type is in each model,
    # plus the authors' own confidence-interval call on a second sheet.
    "celltypes": ("elife-98348-supp5-v1.xlsx", 43_284),
    # Supplementary file 7: differential expression, including a fly gene list
    # ranked by how deregulated it is across all five models.
    "genes": ("elife-98348-supp7-v1.xlsx", 11_632_233),
}


def pech_tables(which=None) -> dict[str, Path]:
    """The Pech et al. 2025 supplementary tables, as local .xlsx paths.

    Read them with `pd.read_excel`, which needs `openpyxl`.

    Both sheets need tidying before use; the project notebook does it in the
    open rather than hiding it here, because dealing with a spreadsheet laid
    out for human eyes is most of what using published tables involves. The
    two gotchas: in the cell type sheet the model columns are followed by
    unnamed spacer columns and summary columns, and in the gene sheet the
    column holding the gene names is called `Unnamed: 10`, because
    `sorted by median` is a heading spanning the block rather than a column.
    """
    which = which or list(PECH_FILES)
    print("Pech et al. 2025 - five fly Parkinsonism models (eLife 98348)")
    out = {}
    for key in which:
        name, size = PECH_FILES[key]
        out[key] = fetch(f"{ELIFE_PECH}/{name}", RAW / "pech" / name, size)
    return out


# ------------------------------------------------------------------ genome


def dm6_genome() -> Path:
    """dm6 FASTA from UCSC, unzipped. ~44 MB gzipped, ~170 MB unzipped."""
    import gzip

    fa = GENOME / "dm6.fa"
    if fa.exists():
        print(f"  already here: dm6.fa ({fa.stat().st_size / 1e6:.0f} MB)")
        return fa

    # The archive is unzipped to a temporary name and only moved into place once
    # it is complete. Without that, a truncated download leaves a partial
    # dm6.fa, which the check above would accept on the next run: the practical
    # would then be reading a genome that silently stops part way through.
    gz = GENOME / "dm6.fa.gz"
    part = GENOME / "dm6.fa.part"
    last_error = None
    for source in (f"{MIRROR}/dm6.fa.gz", f"{UCSC}/dm6.fa.gz"):
        try:
            fetch(source, gz)
            print("  unzipping dm6.fa.gz")
            with gzip.open(gz, "rb") as src, open(part, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 512)
            part.replace(fa)
            gz.unlink(missing_ok=True)
            print(f"    done: {fa.stat().st_size / 1e6:.0f} MB")
            return fa
        except Exception as err:
            last_error = err
            part.unlink(missing_ok=True)
            gz.unlink(missing_ok=True)
            print(f"  {source} did not work ({type(err).__name__}: {err})")

    raise RuntimeError(f"Could not get the dm6 genome. Last error: {last_error}")


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


EMBRYO_CELL_TYPES = tuple(EMBRYO_SIZES_MB)


def dickmanken_embryo_bigwigs(cell_types=("Neuronal", "Glia"), stream: bool = True):
    """Per-cell-type embryo accessibility tracks, the data behind Embryo10x.

    Keys are the Embryo10x model's class names, so predictions and tracks line
    up even where the server spells a name differently.

    These files are large, 49 to 343 MB each and about 2.4 GB for all twenty,
    which used to make them impractical for a course. They are read over the
    network instead, so asking for all twenty costs no more than asking for two.
    Pass `cell_types=EMBRYO_CELL_TYPES` for the lot, or `stream=False` to force
    a download.

    Note this is a single timepoint (16-20 h AEL), not a developmental series.
    For a series within one dataset, use `janssens_devel_bigwigs`.
    """
    cell_types = list(cell_types)
    print("Dickmanken et al. 2026 - embryo scATAC bigwigs (16-20 h AEL)")
    return _resolve_tracks(
        {ct: (f"{EMBRYO_BW}/{EMBRYO_FILE_NAME.get(ct, ct)}.bw",
              RAW / "dickmanken_embryo_bw" / f"{ct}.bw")
         for ct in cell_types},
        stream,
        EMBRYO_SIZES_MB,
    )


if __name__ == "__main__":
    janssens_celltype_info()
    janssens_adult_peaks()
    janssens_brain_peaks()
    janssens_adult_bigwigs()
    merrill_atac()
    dm6_genome()
