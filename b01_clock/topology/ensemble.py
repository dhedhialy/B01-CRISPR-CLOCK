"""Topology-ensemble timing: union of feasible sets across audited topologies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np

from b01_clock.model.observations import CharacterMatrix
from b01_clock.model.tree import LineageTree
from b01_clock.solver.exact import ExactSolverResult, exact_profile_bounds


@dataclass
class TopologyEnsembleResult:
    target: str
    union_lower: float
    union_upper: float
    robust_order: bool
    per_topology: list
    note: str


def leaf_swap_ensemble(tree: LineageTree, max_topologies: int = 4) -> List[Tuple[str, LineageTree]]:
    """Finite audited ensemble: identity + leaf swaps across clades."""
    leaves = tree.leaves()
    out = [("identity", tree.copy())]
    if len(leaves) < 4:
        return out
    for i, j in [(0, 2), (1, 3), (0, 3)]:
        if len(out) >= max_topologies:
            break
        t = tree.copy()
        a, b = leaves[i], leaves[j]
        pa, pb = t.nodes[a].parent, t.nodes[b].parent
        if pa is None or pb is None or pa == pb:
            continue
        t.nodes[pa].children = [b if c == a else c for c in t.nodes[pa].children]
        t.nodes[pb].children = [a if c == b else c for c in t.nodes[pb].children]
        t.nodes[a].parent = pb
        t.nodes[b].parent = pa
        t._recompute_ages()
        out.append((f"swap:{a}<->{b}", t))
    return out


def ensemble_timing_union(
    tree: LineageTree,
    matrix: CharacterMatrix,
    rate_bounds: Sequence[Tuple[float, float]],
    total_time: float,
    targets: Sequence[str],
    confidence: float = 0.95,
    grid_size: int = 31,
    max_topologies: int = 4,
) -> List[TopologyEnsembleResult]:
    """Compute union of compatible timing intervals across topology ensemble.

    Temporal orders that remain true across all admissible topologies are flagged robust.
    """
    ensemble = leaf_swap_ensemble(tree, max_topologies=max_topologies)
    per_target = {t: [] for t in targets}

    for topo_id, topo in ensemble:
        # Remap matrix leaf keys if swaps renamed — leaf names stay, only parents change
        res = exact_profile_bounds(
            topo,
            matrix,
            rate_bounds,
            total_time,
            targets=list(targets),
            confidence=confidence,
            grid_size=grid_size,
        )
        for r in res:
            per_target[r.target].append({"topology": topo_id, "result": r})

    out = []
    for target, rows in per_target.items():
        lowers = [row["result"].lower for row in rows]
        uppers = [row["result"].upper for row in rows]
        u_lo, u_hi = float(min(lowers)), float(max(uppers))
        # robust order for order_gap targets: all intervals strictly positive or all negative
        robust = False
        if target.startswith("order_gap:"):
            if all(row["result"].lower > 0 for row in rows):
                robust = True
            if all(row["result"].upper < 0 for row in rows):
                robust = True
        out.append(
            TopologyEnsembleResult(
                target=target,
                union_lower=u_lo,
                union_upper=u_hi,
                robust_order=robust,
                per_topology=[
                    {
                        "topology": row["topology"],
                        "lower": row["result"].lower,
                        "upper": row["result"].upper,
                        "status": row["result"].status,
                    }
                    for row in rows
                ],
                note=(
                    "Union of profile LR sets across audited topologies; "
                    "robust_order means every topology agrees on sign of the gap."
                ),
            )
        )
    return out
