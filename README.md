# ISN Neurochemistry School 2026: genomics practicals

Working with public single-cell data from the *Drosophila* brain, and reading
enhancers with a deep learning model.

ISN Neurochemistry School 2026, "Neurogenetics for Brain Health", University of
Port Harcourt, Nigeria, 25-31 October 2026. Taught by Niklas Kempynck, Aerts
lab, VIB-KU Leuven.

## For participants

Click a badge. The notebook opens in Google Colab. Run the first cell, wait a
minute while it installs, then work down the page.

You need a Google account. Nothing is installed on your own computer.

| | |
|---|---|
| **Practical 5**: public single-cell RNA and ATAC data | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nkempynck/isn-school-2026/blob/main/notebooks/practicals/05_public_atac.ipynb) |
| **Practical 6**: enhancer modelling with CREsted | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nkempynck/isn-school-2026/blob/main/notebooks/practicals/06_crested_deepflybrain.ipynb) |

Two things that are normal and not your fault:

- after the install, Colab may show a **Restart session** button. Click it, then
  run the first cell again.
- if you leave the page for a long time Colab disconnects and forgets
  everything. Reconnect and run the cells from the top again.

Each notebook already contains its outputs, so if something fails you can still
read what should have happened.

## What is in them

**Practical 5** uses the Fly Cell Atlas to find which cell types express a gene,
runs a differential expression test to find marker genes, looks at an earlier
developmental stage in the larval brain, and then asks whether the chromatin
near those genes is accessible in the same cell types, using ATAC data from
Janssens et al. 2022.

**Practical 6** gives 500 bp of DNA to two pretrained networks, DeepFlyBrain and
Embryo10x, and asks which cell type the sequence belongs to and which
nucleotides the model used to decide. It ends on an enhancer that the authors
tested in flies, where the answer is known in advance.

Everything runs on a CPU. No model is trained.

## Running it on your own machine instead

See `Setup.md`. In short:

```bash
uv sync
uv run jupyter lab
```

The notebooks detect whether they are on Colab and skip the install step when
they are not.

## Data

Nothing in this repository is data. Every file is downloaded from the original
source by `scripts/download_data.py`, which prints each URL as it goes:

- Fly Cell Atlas, https://flycellatlas.org
- Janssens et al. 2022 fly brain scATAC, https://flybrain.aertslab.org and GEO GSE163697
- Merrill et al. 2022 TH vs GABA ATAC, GEO GSE197760
- Dickmanken et al. 2026 embryo scATAC, GEO GSE293575
- dm6 genome, UCSC

`scripts/prepare_data.py` reduces the 2.5 GB Fly Cell Atlas head file to a 42 MB
brain-only file once. The reduced files are mirrored so that Colab sessions do
not have to repeat it.

## Credits

Course materials by Niklas Kempynck. The data belong to the groups that
published it; please cite them, not this repository.
