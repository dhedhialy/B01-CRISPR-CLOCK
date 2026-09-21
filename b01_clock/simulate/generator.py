"""Ground-truth simulation under the irreversible exposure model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from b01_clock.model.observations import (
    EDITED,
    MISSING,
    UNEDITED,
    CharacterMatrix,
    ObservationConfig,
    apply_observation_layer,
)
from b01_clock.model.tree import LineageTree


@dataclass
class SimulationTruth:
    tree: LineageTree
    rates: np.ndarray
    branch_times: np.ndarray
    total_time: float
    matrix: CharacterMatrix
    ancestral: Dict[str, np.ndarray]
    meta: dict = field(default_factory=dict)


def simulate_recorder(
    tree: LineageTree,
    rates: np.ndarray,
    total_time: Optional[float] = None,
    branch_rate_scales: Optional[np.ndarray] = None,
    observation: Optional[ObservationConfig] = None,
    seed: int = 0,
) -> SimulationTruth:
    """Simulate irreversible edits site-by-site along branches; root starts unedited."""
    rng = np.random.default_rng(seed)
    tree = tree.copy()
    times = tree.branch_times_vector()
    if branch_rate_scales is None:
        branch_rate_scales = np.ones(len(times))
    n_sites = len(rates)
    ancestral: Dict[str, np.ndarray] = {
        tree.root: np.zeros(n_sites, dtype=int),
    }

    def dfs(u: str) -> None:
        for v in tree.nodes[u].children:
            t = tree.nodes[v].time_from_parent
            eidx = tree.edge_index()[(u, v)]
            scale = float(branch_rate_scales[eidx])
            parent = ancestral[u]
            child = parent.copy()
            for s in range(n_sites):
                if parent[s] == EDITED:
                    child[s] = EDITED
                    continue
                if parent[s] == MISSING:
                    child[s] = MISSING
                    continue
                p_edit = 1.0 - np.exp(-rates[s] * scale * t)
                child[s] = EDITED if rng.random() < p_edit else UNEDITED
            ancestral[v] = child
            dfs(v)

    dfs(tree.root)
    leaf_states = {ℓ: ancestral[ℓ].copy() for ℓ in tree.leaves()}
    matrix = CharacterMatrix(leaf_states=leaf_states, n_sites=n_sites, ancestral=None)
    if observation is not None:
        matrix = apply_observation_layer(matrix, observation)
    T = total_time if total_time is not None else tree.total_depth()
    return SimulationTruth(
        tree=tree,
        rates=rates.copy(),
        branch_times=times.copy(),
        total_time=float(T),
        matrix=matrix,
        ancestral=ancestral,
        meta={
            "seed": seed,
            "branch_rate_scales": branch_rate_scales.tolist(),
            "observation": None if observation is None else observation.__dict__,
        },
    )


def draw_site_rates(
    n_sites: int,
    mean: float = 1.0,
    heterogeneity: float = 0.0,
    seed: int = 0,
    clip: Tuple[float, float] = (1e-3, 50.0),
) -> np.ndarray:
    """Log-normal site rates; heterogeneity = sd of log rates."""
    rng = np.random.default_rng(seed)
    if heterogeneity <= 0:
        return np.full(n_sites, mean, dtype=float)
    logs = rng.normal(np.log(mean), heterogeneity, size=n_sites)
    rates = np.exp(logs)
    return np.clip(rates, clip[0], clip[1])
