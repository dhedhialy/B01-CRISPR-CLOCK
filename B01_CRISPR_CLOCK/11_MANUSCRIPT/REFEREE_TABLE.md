# Referee Attack Map — Diagnostics, Sensitivities, Limitations

**Rule:** Every plausible adversarial referee objection maps to a pre-committed response modality.  
**Modalities:** `DIAG` = diagnostic already in `b01_clock`; `SENS` = sensitivity analysis (supplement); `LIM` = Limitations section concession; `LOCK` = locked estimand/constraint clarification (not a bug).

| ID | Referee attack | Modality | Response / artifact |
|----|----------------|----------|---------------------|
| R01 | “You just rediscovered that molecular clocks need calibration.” | LOCK | We formalize **exposure vs time** under irreversible CRISPR editing and emit **set-valued** certificates, not calibrated point chronograms. Cite `INVARIANCE.md`, Abstract wedge. |
| R02 | “This is ConvexML with worse scalability.” | LOCK | Different estimand: \(\mathcal{I}_v\) vs point reconstruction. ConvexML is primary **dataset**, not competitor objective. `NOVELTY_COLLISION.md`. |
| R03 | “Pilarski 2026 (or similar) already times recorders better by MSE.” | LOCK + SENS | MSE assumes point ID. Report coverage + false-ID; optionally show MSE of our MLE **inside** \(\mathcal{I}_v\) without claiming superiority. |
| R04 | “Absolute time *is* identified; your invariance is a toy.” | DIAG | `check_rate_time_invariance` + `rate_time_transform_feasible`. Nontrivial \(\Psi_c\) exits \(\mathcal{U}(T)\times\mathcal{R}\). |
| R05 | “Partial identification means your method failed.” | LOCK | Regime P is the **target success** under nondegenerate rate box + single \(T\). Kill criterion is Regime U for all targets, not Regime P. |
| R06 | “Intervals are too wide to be useful.” | SENS + LIM | Width is scientific content. Show narrowing vs trivial \([0,T]\); sensitivity with tighter boxes / more sites. Concede usefulness is question-dependent (`LIM`). |
| R07 | “\(\chi^2_1\) cutoff is invalid under partial ID.” | DIAG + LIM | Empirical `evaluate_coverage`; interpret as confidence set for projection. Concede asymptotic justification is classical point-ID (`LIM` §7.5). |
| R08 | “Grid solver misses the true profile boundary.” | DIAG | `verify_against_grid`; `grid_unstable` flag; refine `grid_size`. |
| R09 | “Global rate assumption is false; sites differ.” | SENS + DIAG | Adversarial `site_heterogeneity` axis; per-site rate mode sensitivity. Primary lock explained in `ESTIMAND_LOCK` (avoid silent scale freedom). |
| R10 | “Branch-specific rate variation breaks certificates.” | DIAG | `branch_rate_cv` adversarial axis; report coverage/false-ID under misspecification. |
| R11 | “You conditioned on the wrong tree.” | DIAG | `ensemble_timing_union`; report union hull; `robust_order`. |
| R12 | “Ensemble is tiny / not Bayesian model averaging.” | LIM | Finite audited leaf-swap class only; not posterior over trees (`LIM` §7.4). |
| R13 | “Leakage: you used all sites to set rates and times.” | DIAG | `HOLDOUT_LOCK` + `lock_split_record`; abstain if violated. |
| R14 | “Holdout split was cherry-picked after seeing results.” | DIAG | Pre-registered seed + `locked_before_headline: true` in JSON before headline fit. |
| R15 | “Rate bounds are arbitrary; you can always shrink intervals.” | SENS | Quantile/pad sensitivity; show monotone expansion as box widens; provenance of bounds from calib sites only. |
| R16 | “Synthetic scout is cherry-picked to PASS.” | DIAG | Fixed `configs/scout.yaml` seeds; `scripts/reproduce.sh`; adversarial battery beyond scout. |
| R17 | “Coverage collapses under dropout/silencing.” | DIAG | Adversarial `dropout_prob` / `silencing_prob`; report when we **abstain**. |
| R18 | “Saturation makes the recorder uninformative, but you still report ages.” | DIAG | Status `unidentified` when width \(\approx T\); pass/kill K1. |
| R19 | “You claim point ID when only partial ID holds.” | DIAG | `false_identification_rate`; guardrail; taxonomy thresholds in `classify_identification`. |
| R20 | “Ultrametric / synchronized division is unrealistic.” | SENS | `division: stochastic` adversarial axis; non-ultrametric sensitivity (if implemented) marked SENS. |
| R21 | “The external data have no ground-truth ages, so nothing is validated.” | LOCK + LIM | The deposit is **itself simulated** (Cassiopeia), so truth is known. We run the ground-truth coverage check and report the failure (§5.4): 1/9 covered, systematic late shift. No wet-lab time-series is certified (`LIM` §7.7). |
| R22 | “ConvexML data are unavailable / results are fake-real.” | DIAG | Official Dryad bundle downloaded & sha256-certified (`DOWNLOAD_STATUS.json` ok; trees.tgz `b2339df3…`). We do **not** relabel it: provenance records `topology_provenance: simulated`, and the manuscript calls it an **independently simulated Cassiopeia benchmark**. The ground-truth coverage check (0/8 → 1/9 with the basal split) is presented alongside it as the honest negative finding, not as “real-data validation.” Anubis 'fast' PoW solved in `anubis_get` (disclosed in Methods). |
| R23 | “SciPhy sequential states break \(1-e^{-rt}\).” | SENS + LIM | Collapse map frozen before split; multi-state likelihood as sensitivity; `LIM` §7.3. |
| R24 | “PALINCODE architecture differs; negative result kills paper.” | LOCK | PALINCODE is robustness, not primary kill switch (`SOURCE_LOCK` hierarchy). |
| R25 | “Exact method does not scale to large trees.” | LIM + SENS | Acknowledge `LIM` §7.1; scalable outer solver validated on the 400-leaf official Dryad tree. **Internal caveat logged:** the first scalable bound pass (penalized full-space coordinate descent) was buggy and was replaced with profile-grid enumeration (`profile_grid_target`) after fixed-age checks disproved its `[1,1]` certificates. |
| R26 | “Prior/Bayesian relaxed clocks already solve this.” | LOCK | Soft priors reintroduce unidentified directions as prior mass; we report constraint-conditional sharp sets, not prior-posterior chronograms. |
| R27 | “Order gaps crossing zero mean you learned nothing.” | LOCK | Ages can be Regime P while order is U; both reported. Robust order is a stronger claim. |
| R28 | “Comparator MLE falls outside your interval — contradiction.” | SENS | Diagnose constraint mismatch (different \(\mathcal{R}\), \(T\), topology); not auto-refutation. |
| R29 | “Missingness not MAR; certificates anti-conservative.” | SENS + LIM | Silencing axis; concede unmodeled informative missingness (`LIM`). |
| R30 | “You overclaim biology / developmental timing.” | LOCK | Role split: computational certificates only (`ESTIMAND_LOCK` §5). |
| R31 | “Irreversibility is violated (scar edits / overwriting).” | SENS | If data show \(1\to 0\), likelihood \(-\infty\) / model misspec; document allele complexity sensitivity. |
| R32 | “Confidence level shopping.” | LOCK | Default \(\alpha\) in config; primary \(0.95\); alternate \(\alpha\) supplement only. |
| R33 | “Multiple targets / multiple testing.” | SENS | Pre-register primary target (`age:A` scout); others secondary; optional multiplicity note. |
| R34 | “Topology union makes everything unidentified (K2).” | DIAG | Kill/reframe K2 if union restores \(\approx[0,T]\); report identity vs union side-by-side. |
| R35 | “Software bugs could fake invariance.” | DIAG | `tests/test_invariance.py`, `tests/test_exact_solver.py`; CI via `reproduce.sh`. |

---

## Coverage checklist

- Invariance / false clock: R04, R19, R35  
- Estimand collision: R01–R03, R05–R06, R26, R30  
- Statistical calibration: R07–R08, R16, R32–R33  
- Misspecification: R09–R12, R17–R20, R29, R31  
- Holdout / leakage: R13–R15  
- Data / transfer: R21–R24, R28  
- Scalability: R25, R34  

Any new referee attack in review must add a row before rebuttal letter freezes.
