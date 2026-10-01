"""Cell type names, and the fact that the two datasets do not agree on them.

The Janssens ATAC files and the Fly Cell Atlas describe the same cells but use
different names. This table is the bridge. Checked against the Fly Cell Atlas
annotation list on 2026-09-22; every value here exists in the atlas.
"""

ATAC_TO_RNA = {
    "KC_g": "gamma Kenyon cell",
    "KC_ab": "alpha/beta Kenyon cell",
    "KC_abPrime": "alpha'/beta' Kenyon cell",
    "T1": "columnar neuron T1",
    "T2": "T neuron T2",
    "T2a": "T neuron T2a",
    "T3": "T neuron T3",
    "Astrocyte_like": "adult reticular neuropil associated glial cell",
    "Ensheathing_glia": "ensheathing glial cell",
    "Cortex_glia": "optic-lobe-associated cortex glial cell",
    "Perineurial_glia": "adult brain perineurial glial cell",
    "Subperineurial_glia": "subperineurial glial cell",
    "Chiasm_glia": "adult optic chiasma glial cell",
}

RNA_TO_ATAC = {v: k for k, v in ATAC_TO_RNA.items()}
