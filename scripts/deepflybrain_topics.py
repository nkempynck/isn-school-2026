"""Which cell type does each DeepFlyBrain topic belong to?

DeepFlyBrain has 81 outputs. They are called Topic_1 ... Topic_81, which tells
you nothing on its own. A "topic" is a set of regions that tend to be open in
the same cells - cisTopic found these groups, and the model learned to predict
them from DNA sequence. Some topics line up neatly with one cell type; most
were never annotated.

Source: the "Topic - Cell Type correspondence" legend on
https://flybrain.aertslab.org (Janssens et al. 2022, Extended Data Fig. 10).
Read on 2026-09-22. 26 of the 81 topics are annotated; the rest are left as
"not annotated" in the paper, and we keep them that way rather than guessing.
"""

TOPIC_CELL_TYPE = {
    1: "Pan-neuron",
    3: "BEAF-32, Pan-glia",
    8: "Pan-neuron",
    9: "BEAF-32",
    10: "T3",
    18: "T4",
    20: "T2",
    21: "Alpha/Beta Kenyon cells",
    22: "Chiasm glia",
    23: "T1",
    25: "Cortex glia",
    32: "T4/T5",
    34: "Perineurial glia",
    35: "Gamma Kenyon cells",
    36: "Subperineurial glia",
    40: "Pan-glia",
    43: "BEAF-32",
    44: "T2a",
    56: "Ensheathing glia",
    59: "Pan-glia",
    60: "Pan-neuron",
    65: "Pan-neuron",
    68: "Astrocyte-like",
    72: "Ensheathing glia (partly)",
    77: "Alpha'/Beta' Kenyon cells",
}

# BEAF-32 is an insulator protein, not a cell type. Topics 3, 9 and 43 are
# dominated by its binding sites, so a high score there usually means
# "this looks like an insulator", not "this is a specific neuron".
INSULATOR_TOPICS = {3, 9, 43}


def label(topic) -> str:
    """Turn 'Topic_35' (or 35) into 'Topic_35 (Gamma Kenyon cells)'."""
    n = int(str(topic).replace("Topic_", ""))
    name = TOPIC_CELL_TYPE.get(n)
    return f"Topic_{n} ({name})" if name else f"Topic_{n} (not annotated)"


def annotated_topics() -> list[str]:
    """The 26 topic names that have a cell type, in numeric order."""
    return [f"Topic_{n}" for n in sorted(TOPIC_CELL_TYPE)]


# Handy groupings for the practicals.
KENYON_CELL_TOPICS = ["Topic_21", "Topic_35", "Topic_77"]
T_NEURON_TOPICS = ["Topic_10", "Topic_18", "Topic_20", "Topic_23", "Topic_32", "Topic_44"]
GLIA_TOPICS = ["Topic_22", "Topic_25", "Topic_34", "Topic_36", "Topic_56", "Topic_68"]
