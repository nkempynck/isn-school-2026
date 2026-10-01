# Setup: ISN School 2026 practicals (local environment)

This sets up a local folder and Python environment on your MacBook Pro (M1) that mirrors what participants will run on their own laptops. Everything runs on CPU, with no GPU and no cloud.

## Design decisions

- **Environment manager: `uv`.** CREsted recommends it, it is fast, and it works the same on macOS, Windows and Linux. A `uv.lock` file pins exact versions, so every participant gets the same environment.
- **Python 3.12.** CREsted needs Python ≥ 3.11, and TOMTOM motif matching (via memelite) only works on Python ≤ 3.12. 3.12 satisfies both.
- **Deep learning backend: PyTorch (CPU).** CREsted runs on Keras 3 with either TensorFlow or PyTorch. PyTorch installs cleanly on Apple Silicon and on Windows without GPU drivers, which matters for participants' laptops. We only run inference with pretrained models (e.g. DeepFlyBrain, 12 MB), so CPU speed is fine.
- **Jupyter notebooks in JupyterLab**, run locally.

## 1. Install uv (once)

macOS / Linux:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Windows (PowerShell):
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Close and reopen the terminal, then check with `uv --version`.

## 2. Create the project folder

```bash
mkdir isn-school-2026 && cd isn-school-2026
uv init --python 3.12 --no-readme
mkdir -p data/raw data/processed models genome notebooks/practicals notebooks/projects scripts docs
```

Put `CLAUDE.md` (handover notes) and this `SETUP.md` in the folder root.

## 3. Add dependencies

```bash
uv add "crested[motif]" torch scanpy jupyterlab ipykernel pandas matplotlib seaborn pyBigWig
```

Notes:
- `crested[motif]` adds tf-modisco-lite and motif tools. Drop the extra if install problems arise; basic CREsted (import, preprocessing, predictions) works without it.
- `pyBigWig` is for simple bigwig reading and plotting in notebooks. Check whether CREsted already pulls in a bigwig reader and remove it if redundant.
- Optional genome browser view in notebooks: try `igv-notebook`. Test that it works in local JupyterLab before relying on it.

After this, `pyproject.toml` and `uv.lock` define the environment. Commit both.

## 4. Set the Keras backend

CREsted uses Keras 3, which reads the `KERAS_BACKEND` environment variable. Put this at the top of every notebook, before importing crested:

```python
import os
os.environ["KERAS_BACKEND"] = "torch"
```

Doing it in the notebook rather than in shell config means it works identically on every participant's machine.

## 5. Launch

```bash
uv run jupyter lab
```

## 6. Smoke test (becomes `notebooks/00_setup_check.ipynb`)

```python
import os
os.environ["KERAS_BACKEND"] = "torch"

import crested, scanpy as sc, keras
print("crested", crested.__version__, "| scanpy", sc.__version__, "| keras backend", keras.backend.backend())

model_path, output_names = crested.get_model("DeepFlyBrain")
model = crested.utils.load_model(model_path)
pred = crested.tl.predict("A" * 500, model)
print("DeepFlyBrain OK, output shape:", pred.shape, "(expect 81 topics)")
```

If this runs, the core stack works.

## 7. Genome (needed to fetch sequences for predictions)

DeepFlyBrain uses dm6. Download once into `genome/`:

```bash
curl -L -o genome/dm6.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip genome/dm6.fa.gz
```

Check how `crested.Genome` / `crested.register_genome` expect the FASTA (and whether a `.fai` index or chrom sizes file is needed).

## 8. Participant laptops: what to test before October

The school's value depends on this working on their machines, not just yours.

1. **Test on Windows**, ideally a modest laptop with 8 GB RAM. Most participants will likely be on Windows. Check that `uv sync` works, that `torch` installs, that `pyBigWig` installs (it has had build issues on Windows; if so, fall back to CREsted's own bigwig reading or another reader), and memory use when loading the scRNA-seq files.
2. **Participant install should be three commands:** install uv, then `uv sync` in the course folder, then `uv run jupyter lab`. Package the course folder as a zip (or GitHub repo) with `pyproject.toml`, `uv.lock`, notebooks and scripts, but no data.
3. **Pre-school check:** send the zip about 2 weeks before, and ask everyone to run `00_setup_check.ipynb` and send a screenshot. This catches broken installs early.
4. **Offline fallback on USB sticks:** the course folder, the downloaded data files, the dm6 genome, the DeepFlyBrain model file, and ideally a pre-downloaded wheel cache for Windows and macOS (`uv` can install from a local cache / `--find-links` folder). Test installing from USB with Wi-Fi switched off.
5. **Memory:** the adult brain scRNA-seq loom is about 495 MB on disk and can be several GB in memory. Test on an 8 GB machine; consider `backed="r"` mode in scanpy or subsetting.