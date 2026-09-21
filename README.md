# B01 — When Is a CRISPR Lineage Recorder a Clock?

Sharp timing sets under heterogeneous editing rates — theory, computation, and validation.

## Start here

**Research workspace:** [`B01_CRISPR_CLOCK/`](B01_CRISPR_CLOCK/README.md)  
**Software package:** `b01_clock/`  
**Reproduce figures + tables:**

```bash
source .venv/bin/activate
export MPLCONFIGDIR=$PWD/.mplconfig
python scripts/run_full_research.py
python -m b01_clock scout --config configs/scout.yaml
python scripts/validate_convexml_ground_truth.py
```

## Status

Decisive synthetic scout **PASS**. Exact small-tree intervals narrower than `[0,T]`, synthetic coverage near nominal, holdout transfer works.

External-benchmark signal is **negative and reported as such** (open science, not a headline limitation-burial): the official ConvexML Dryad bundle is a **simulated Cassiopeia deposit**, so its newick branch lengths give known ground truth. Baseline certificates are sharp — `age:N1 = [0.073, 0.146]` — but ground-truth coverage on the deposit is **1/9 (11%)**: deep splits are systematically over-aged or pinned at \(T\) (`convexml_ground_truth_coverage.json`). The primary binary/global-rate lock does **not** transfer to this generated process; no wet-lab recorder source is certified yet. A scalable-solver bound bug that had produced spurious `[1,1]` certificates was found during this check and fixed (`profile_grid_target` in `b01_clock.solver.optimize`).

## Framing

Not a better clock estimator. The estimand is the **sharp compatible interval** for node ages / branch times given locked rate bounds and anchors, plus exact conditions for partial vs point identification.
