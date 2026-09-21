# B01 research status

**Portfolio:** LAUNCH — Independent high-ceiling theory  
**Framing:** Sharp partial-identification of compatible timing sets (not a better clock estimator)

## Go / no-go (decisive scout)

| Criterion | Status |
|-----------|--------|
| Exact small-tree feasible interval under rate bounds + dated T | PASS — age(A) materially narrower than [0,T] |
| Synthetic coverage near nominal | PASS — 0.95 at σ=0; false-ID=0 |
| Holdout rate-bound transfer with narrowing | PASS — width/T ≈ 0.46 |
| Rate–time invariance + guardrail | PASS |
| ConvexML Dryad primary tree | Official `trees.tgz` downloaded (403 MB) & sha256-certified; 400 leaves × 39 sites |

**Verdict: PASS for the synthetic gate; FAIL for the external benchmark — reported honestly.** External holdout on official Dryad tree_0 (a **simulated Cassiopeia deposit**, not wet-lab data; fixed 64-leaf induced certification): with the corrected solver, basal-split `age:N1` = [0.073, 0.146] partially identified (scout draw [0.048, 0.067]). Ground-truth coverage check: **1/9 (11%)** — sharp intervals are **miscalibrated** on this generated process (late-shifted deep splits; several pinned at T). Primary lock does not transfer to the benchmark; no wet-lab source is certified yet. Anubis 'fast' PoW solved in `b01_clock.data.adapters.convexml.anubis_get`; idempotent `download_convexml` re-verifies sha256; provenance tagged `topology_provenance: simulated`. Target naming is deterministic (`basal_split` canonicalizes children).

**Real-data consistency (GESTALT, McKenna 2016, GSE81713):** see `b01_clock/data/adapters/gestalt.py` + `scripts/run_gestalt_real_consistency.py`. Real fish ADR1 MIX matrix (601 taxa x 1156 chars, edited_frac 0.004) reproduces the same late-shift the benchmark showed: all age MLEs hug the 5 dpf harvest (inputs are nearly uninformative for timing). Benchmark-trained conformal bins transfer in-domain; the boundary-pinned degenerate cert is flagged, not reported as a point. Real lineages have **no internal-age ground truth**, so this is a consistency check with **no coverage claim** -> `07_REAL_DATA/gestalt_real_consistency.json`. Empirical recalibration evidence lives in `../experiments/recalibration_study.json` and `recalibration_conditional.json`.

## Workspace map

See folders `00_LOCK` … `12_CERTIFICATES`. Reproduce with:

```bash
source .venv/bin/activate
export MPLCONFIGDIR=$PWD/.mplconfig
python scripts/run_full_research.py
python -m b01_clock scout --config configs/scout.yaml
python scripts/validate_convexml_ground_truth.py   # external-benchmark coverage check
```

## Role split (locked)

- Theory/solver/adversarial sims: identification side  
- Data audit/adapters/empirical validation: coworker side  
- Joint: estimand, source lock, manuscript
