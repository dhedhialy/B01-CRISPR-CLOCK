#!/usr/bin/env python3
"""Score `age:` certificate coverage against known internal-node times.

The ConvexML Dryad deposit is a Cassiopeia-generated benchmark (see its
README): the newick branch lengths ARE the true node times. This is therefore
a genuine, leakage-free coverage check (analysis sites only) — it is NOT a
claim about wet-lab data.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from b01_clock.data.adapters.convexml import load_convexml_example
from b01_clock.model.observations import CharacterMatrix
from b01_clock.solver.optimize import scalable_timing_certificate
from b01_clock.validation.calibration import held_out_rate_bounds, split_sites

T = 1.0
OUT = ROOT / "B01_CRISPR_CLOCK" / "07_REAL_DATA"


def main(rep: int = 0):
    bundle = load_convexml_example(ROOT / "data_store" / "convexml", total_time=T, rep=rep)
    tree, matrix, meta = bundle.tree, bundle.matrix, bundle.meta
    calib, analysis = split_sites(matrix.n_sites, calib_frac=0.4, seed=7)
    bounds = held_out_rate_bounds(tree, matrix, calib, T)
    ana = matrix.site_slice(analysis)

    rng = np.random.default_rng(13)  # same draw as run_convexml_reanalysis.py
    keep = sorted(rng.choice(sorted(tree.leaves()), size=64, replace=False).tolist())
    ctree = tree.induced_subtree(keep)
    ana_mat = CharacterMatrix({l: ana.leaf_states[l] for l in keep}, ana.n_sites)

    cand = [u for u in ctree.internal_nodes() if u != ctree.root]
    ranked = sorted((ctree.node_age(u), u) for u in cand)
    n = len(ranked)
    picks = []
    for q in (0.05, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 0.98):
        u = ranked[min(n - 1, int(q * n))][1]
        if u not in picks:
            picks.append(u)

    # canonical headline target: lexicographically smallest child of the basal split
    u0, a0, b0 = ctree.basal_split()
    n1 = sorted([a0, b0])[0]
    if n1 not in picks:
        picks.insert(0, n1)

    rows = []
    for u in picks:
        true = ctree.node_age(u)
        res = scalable_timing_certificate(
            ctree, ana_mat, bounds, T, targets=[f"age:{u}"],
            confidence=0.95, n_starts=2, seed=7,
        )[0]
        rows.append({
            "target": f"age:{u}",
            "true_age": true,
            "lower": res.lower,
            "upper": res.upper,
            "mle": res.mle,
            "status": res.status,
            "covered": res.lower <= true <= res.upper,
            "width_frac": (res.upper - res.lower) / T,
        })

    covered = sum(r["covered"] for r in rows)
    out = {
        "dataset": "convexml_official_default",
        "dryad_doi": "10.5061/dryad.qrfj6q5nz",
        "note": "Deposit is a Cassiopeia-simulated benchmark; newick branch lengths are known truth.",
        "n_leaves": len(tree.leaves()),
        "n_leaves_certified": len(keep),
        "n_sites": matrix.n_sites,
        "transferred_bounds": list(bounds),
        "targets": len(rows),
        "coverage_frac": covered / len(rows),
        "rows": rows,
    }
    (OUT / "convexml_ground_truth_coverage.json").write_text(
        json.dumps(out, indent=2, default=str)
    )
    print(json.dumps({
        "coverage_frac": out["coverage_frac"],
        "n_targets": len(rows),
        "bounds": list(bounds),
    }, indent=2))


if __name__ == "__main__":
    main()