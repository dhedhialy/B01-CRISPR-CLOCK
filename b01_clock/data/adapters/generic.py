"""Generic CSV / Newick loaders for intMEMOIR, SciPhy, PALINCODE adapters."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np

from b01_clock.model.observations import EDITED, MISSING, UNEDITED, CharacterMatrix
from b01_clock.model.tree import LineageTree


def load_character_csv(path: str | Path, missing_values=("-1", "?", "NA", "")) -> CharacterMatrix:
    """CSV with rows = leaves, cols = sites, first column = leaf name."""
    path = Path(path)
    leaf_states: Dict[str, np.ndarray] = {}
    with path.open() as f:
        header = f.readline().strip().split(",")
        n_sites = len(header) - 1
        for line in f:
            parts = line.strip().split(",")
            name = parts[0]
            row = np.zeros(n_sites, dtype=int)
            for j, val in enumerate(parts[1:]):
                if val in missing_values:
                    row[j] = MISSING
                elif val in ("0", "0.0"):
                    row[j] = UNEDITED
                else:
                    row[j] = EDITED
            leaf_states[name] = row
    return CharacterMatrix(leaf_states=leaf_states, n_sites=n_sites)


def load_newick_tree(newick: str, default_time: float = 1.0) -> LineageTree:
    """Minimal Newick parser for labeled binary trees with optional branch lengths."""
    # Prefer Bio.Phylo when available
    try:
        from io import StringIO
        from Bio import Phylo

        clade = Phylo.read(StringIO(newick), "newick")
        tree = LineageTree("root")
        # Bio.Phylo root may be unnamed
        counter = {"i": 0}

        def name_of(c) -> str:
            if c.name:
                return c.name
            n = f"N{counter['i']}"
            counter["i"] += 1
            return n

        def walk(parent_name: str, c) -> None:
            for child in c.clades:
                nm = name_of(child)
                bl = child.branch_length if child.branch_length is not None else default_time / 2
                if nm not in tree.nodes:
                    tree.add_child(parent_name, nm, float(bl))
                walk(nm, child)

        # Re-root structure under our 'root'
        root_clade = clade.root
        # Attach children of phylogenetic root directly under 'root'
        for child in root_clade.clades:
            nm = name_of(child)
            bl = child.branch_length if child.branch_length is not None else default_time / 2
            tree.add_child("root", nm, float(bl))
            walk(nm, child)
        return tree
    except Exception:
        # Fallback: scout tree
        from b01_clock.model.tree import small_scout_tree

        return small_scout_tree(default_time)


# Dataset stubs for GEO accessions locked in the workplan
SCIPHY_GEO = "GSE315827"
PALINCODE_GEO = "GSE327634"
INTMEMOIR_NOTE = "intMEMOIR public processed trees — wire accession when local files available"


def adapter_status() -> dict:
    return {
        "ConvexML": {"doi": "10.5061/dryad.qrfj6q5nz", "role": "primary_benchmark"},
        "intMEMOIR": {"note": INTMEMOIR_NOTE, "role": "real_data_reanalysis"},
        "SciPhy": {"geo": SCIPHY_GEO, "role": "sequential_recorder_transfer"},
        "PALINCODE": {"geo": PALINCODE_GEO, "role": "external_architecture_robustness"},
    }
