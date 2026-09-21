#!/usr/bin/env python3
"""ConvexML reanalysis: holdout rate bounds + timing certificates on the official Dryad bundle."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from b01_clock.data.adapters.convexml import load_convexml_example
from b01_clock.data.audit import audit_dataset
from b01_clock.solver.optimize import scalable_timing_certificate
from b01_clock.validation.calibration import (
    held_out_rate_bounds,
    lock_split_record,
    split_sites,
)

T = 1.0
OUT = ROOT / "B01_CRISPR_CLOCK" / "07_REAL_DATA"
OUT.mkdir(parents=True, exist_ok=True)


def main(rep: int = 0):
    bundle = load_convexml_example(ROOT / "data_store" / "convexml", total_time=T, rep=rep)
    tree, matrix, meta = bundle.tree, bundle.matrix, bundle.meta
    from b01_clock.model.observations import CharacterMatrix

    calib, analysis = split_sites(matrix.n_sites, calib_frac=0.4, seed=7)
    split = lock_split_record(calib, analysis, 7, note=f"official Dryad tree_{rep} default regime")
    bounds = held_out_rate_bounds(tree, matrix, calib, T)
    ana = matrix.site_slice(analysis)
    # Full 400-leaf tree leaves ~399 free internal ages; the profile solver scales
    # poorly there (O(d) evals per sweep at ~150ms each). Certify on the induced
    # subtree over a fixed random 64-leaf subset of the official tree — all official
    # Dryad edit states, free-parameter count ~O(64).
    rng = np.random.default_rng(13)
    keep = sorted(rng.choice(sorted(tree.leaves()), size=64, replace=False).tolist())
    ctree = tree.induced_subtree(keep)

    ana_mat = CharacterMatrix(
        leaf_states={ℓ: ana.leaf_states[ℓ] for ℓ in keep},
        n_sites=ana.n_sites,
    )
    _, i, j = ctree.basal_split()
    targets = [f"age:{i}", f"order_gap:{i},{j}"]
    res = scalable_timing_certificate(
        ctree, ana_mat, bounds, T, targets=targets, confidence=0.95, n_starts=2, seed=7
    )
    prov = bundle.provenance
    audit = audit_dataset(
        "ConvexML-official", tree, matrix, prov,
        calibration_sites=list(calib), analysis_sites=list(analysis),
    )
    headline = [
        {
            "target": r.target,
            "lower": r.lower,
            "upper": r.upper,
            "mle": r.mle,
            "status": r.status,
            "width_frac": (r.upper - r.lower) / T,
        }
        for r in res
    ]
    out = {
        "dryad_doi": "10.5061/dryad.qrfj6q5nz",
        "mode": "official_dryad",
        "regime": meta.get("regime"),
        "rep": rep,
        "n_leaves": len(tree.leaves()),
        "n_internals": len(tree.internal_nodes()),
        "n_sites": matrix.n_sites,
        "n_leaves_certified": len(keep),
        "certified_tree_n_internals": len(ctree.internal_nodes()),
        "split": split,
        "transferred_bounds": list(bounds),
        "headline": headline,
        "pass_narrow": any(h["width_frac"] < 0.95 and h["status"] != "unidentified" for h in headline),
        "audit": audit,
    }
    (OUT / "convexml_protocol_reanalysis.json").write_text(json.dumps(out, indent=2, default=str))
    (OUT / "convexml_protocol_certificate.json").write_text(json.dumps({
        "dataset": "convexml_official_default",
        "dryad_doi": "10.5061/dryad.qrfj6q5nz",
        "regime": meta.get("regime"),
        "estimand": targets[0],
        "rate_bounds": list(bounds),
        "anchors": {"total_time": T},
        "validation_split": split,
        "result": headline[0],
        "solver": "scalable_timing_certificate",
    }, indent=2))
    (ROOT / "B01_CRISPR_CLOCK" / "03_DATA_AUDIT" / "convexml_audit.json").write_text(
        json.dumps(audit, indent=2, default=str)
    )
    print(json.dumps({"pass_narrow": out["pass_narrow"], "bounds": list(bounds), "headline": headline}, indent=2))


if __name__ == "__main__":
    main()