#!/usr/bin/env python3
"""Conditional (mle-bin) conformal recalibration on the saved recalibration study rows.

Profile-LR certs over all internals are expensive and already saved in
experiments/recalibration_study.json. This analyzes whether conditioning the
conformal width on the node's MLE (depth proxy) keeps nominal coverage with
much narrower intervals, LOOCV and cross-draw.

Writes experiments/recalibration_conditional.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from b01_clock.validation.conformal import N_BINS, bin_of, conditional_loo, conditional_split

STUDY = json.loads((ROOT / "experiments" / "recalibration_study.json").read_text())


def main():
    out = {}
    for name, rows in [("canonical_draw13", STUDY["canonical_draw13"]["rows"]),
                       ("alt_draw42", STUDY["alt_draw42"]["rows"])]:
        cov, mean_w = conditional_loo(rows)
        out[name] = {"conditional_loo_coverage": cov,
                     "conditional_mean_width": mean_w,
                     "profile_coverage": float(np.mean([x["profile_covered"] for x in rows]))}
        print(f"{name}: cond true-LOOCV cov={cov:.3f}  mean_width={mean_w:.3f}  (profile cov={out[name]['profile_coverage']:.3f})")

    a = STUDY["canonical_draw13"]["rows"]
    b = STUDY["alt_draw42"]["rows"]
    for fname, train, test in [("13->42", a, b), ("42->13", b, a)]:
        cov, mean_w, _ = conditional_split(train, test)
        out[f"transfer_{fname}"] = {"coverage": cov, "mean_width": mean_w}
        print(f"transfer {fname}: cov={cov:.3f}  mean_width={mean_w:.3f}")

    # panel-of-9 must be DISJOINT from the fit: load full-fidelity rows from the
    # frozen ground-truth JSON, and fit bins on draw13 rows whose nodes are NOT
    # in the panel (same draw, excluded-set split).
    gt_path = ROOT / "B01_CRISPR_CLOCK" / "07_REAL_DATA" / "convexml_ground_truth_coverage.json"
    gt = json.loads(gt_path.read_text())
    panel = [{"node": r["target"].split(":", 1)[1], "mle": r["mle"], "true": r["true_age"]}
             for r in gt["rows"]]
    panel_ids = {p["node"] for p in panel}
    train = [x for x in a if x["node"] not in panel_ids]
    cov, mean_w, bins = conditional_split(train, panel)
    res = []
    for p in panel:
        w = bins[bin_of(p["mle"])]
        res.append({"node": p["node"], "mle": p["mle"], "true": p["true"],
                    "width": w, "covered": bool(p["mle"] - w <= p["true"] <= p["mle"] + w)})
    out["canonical_panel_9_nodes_disjoint"] = {
        "coverage": cov, "mean_width": mean_w,
        "train_n_nodes": len(train), "panel_n_nodes": len(panel),
        "rows": res}
    print(f"panel-of-9 (bins fit on {len(train)} NON-panel draw13 nodes): cov={cov:.3f} mean_width={mean_w:.3f}")
    for r in res:
        print(f"  {r['node']:5s} mle={r['mle']:.3f} true={r['true']:.3f} w={r['width']:.3f} covered={r['covered']}")

    (ROOT / "experiments" / "recalibration_conditional.json").write_text(json.dumps(out, indent=2))
    print("wrote experiments/recalibration_conditional.json")


if __name__ == "__main__":
    main()