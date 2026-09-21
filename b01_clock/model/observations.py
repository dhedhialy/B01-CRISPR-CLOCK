"""Observation layers: saturation, outcome states, silencing, sequencing dropout."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .tree import LineageTree


# Encoding: 0 unedited, 1 edited/outcome, -1 missing/dropout/silenced
UNEDITED = 0
EDITED = 1
MISSING = -1


@dataclass
class ObservationConfig:
    dropout_prob: float = 0.0
    silencing_prob: float = 0.0
    saturation_cap: Optional[float] = None  # expected edits/site ceiling (soft)
    seed: Optional[int] = None


@dataclass
class CharacterMatrix:
    """Leaf × site character matrix with optional ancestral states."""

    leaf_states: Dict[str, np.ndarray]  # leaf -> (n_sites,)
    n_sites: int
    ancestral: Optional[Dict[str, np.ndarray]] = None  # node -> (n_sites,)

    def as_array(self, leaf_order: Sequence[str]) -> np.ndarray:
        return np.stack([self.leaf_states[ℓ] for ℓ in leaf_order], axis=0)

    def site_slice(self, indices: Sequence[int]) -> "CharacterMatrix":
        idx = list(indices)
        leaves = {k: v[idx].copy() for k, v in self.leaf_states.items()}
        anc = None
        if self.ancestral is not None:
            anc = {k: v[idx].copy() for k, v in self.ancestral.items()}
        return CharacterMatrix(leaf_states=leaves, n_sites=len(idx), ancestral=anc)


def apply_observation_layer(
    matrix: CharacterMatrix,
    config: ObservationConfig,
) -> CharacterMatrix:
    """Inject MCAR dropout and heritable-style silencing (site×leaf independent here)."""
    rng = np.random.default_rng(config.seed)
    out = {}
    for leaf, row in matrix.leaf_states.items():
        r = row.copy()
        if config.dropout_prob > 0:
            drop = rng.random(r.shape) < config.dropout_prob
            r[drop] = MISSING
        if config.silencing_prob > 0:
            sil = rng.random(r.shape) < config.silencing_prob
            r[sil] = MISSING
        out[leaf] = r
    return CharacterMatrix(leaf_states=out, n_sites=matrix.n_sites, ancestral=matrix.ancestral)


def conservative_ancestral_states(
    tree: LineageTree,
    matrix: CharacterMatrix,
) -> Dict[str, np.ndarray]:
    """Irreversible-aware conservative reconstruction.

    If any descendant leaf is unedited, all ancestors are forced unedited.
    If all observed descendants are edited, the internal state is left MISSING
    (edit timing on the path is uncertain). Root defaults to unedited.
    """
    n = matrix.n_sites
    states: Dict[str, np.ndarray] = {
        ℓ: matrix.leaf_states[ℓ].copy() for ℓ in tree.leaves()
    }

    def postorder(u: str) -> np.ndarray:
        if u in states:
            return states[u]
        child_states = [postorder(c) for c in tree.nodes[u].children]
        out = np.full(n, MISSING, dtype=int)
        for s in range(n):
            vals = [int(cs[s]) for cs in child_states if cs[s] != MISSING]
            if not vals:
                out[s] = MISSING
            elif any(v == UNEDITED for v in vals):
                out[s] = UNEDITED
            else:
                out[s] = MISSING
        states[u] = out
        return out

    postorder(tree.root)
    root = states[tree.root].copy()
    root[root == MISSING] = UNEDITED
    states[tree.root] = root
    return states


def enumerate_transitions(
    tree: LineageTree,
    ancestral: Dict[str, np.ndarray],
    site: int,
) -> list[Tuple[str, str, int, int]]:
    """List (parent, child, parent_state, child_state) for one site."""
    out = []
    for p, c in tree.edges():
        out.append((p, c, int(ancestral[p][site]), int(ancestral[c][site])))
    return out


def missingness_summary(matrix: CharacterMatrix) -> dict:
    arr = np.stack(list(matrix.leaf_states.values()), axis=0)
    total = arr.size
    n_miss = int(np.sum(arr == MISSING))
    n_edit = int(np.sum(arr == EDITED))
    n_un = int(np.sum(arr == UNEDITED))
    return {
        "n_leaves": arr.shape[0],
        "n_sites": arr.shape[1],
        "missing_frac": n_miss / total if total else 0.0,
        "edited_frac_observed": n_edit / max(n_edit + n_un, 1),
        "saturation_proxy": n_edit / total if total else 0.0,
    }
