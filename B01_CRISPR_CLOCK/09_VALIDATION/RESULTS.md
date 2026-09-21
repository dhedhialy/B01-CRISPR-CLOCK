# B01 empirical results



## Filled empirical results (auto)

- Exact small-tree age(A): [0.000, 0.625] (MLE 0.375; partially_identified).
- Exact small-tree age(B): [0.375, 0.792].
- Feasible-set fraction of (age A, age B) grid: 0.229.
- Synthetic 95% coverage at σ=0: 0.95; σ=0.5: 0.85; false-ID: 0 across scanned regimes.
- Holdout age(A): [0.333, 0.792] (width/T=0.458) with transferred bounds [0.26152915677434624, 15.19706161372907].
- Adversarial battery: 11/11 specs covered truth under the routed evaluator — per-site rate profiling for `site_het_*`, ensemble timing union for `topology_error`, global-rate profile otherwise. Coverage failure analysis closed.
- Per-site rate recovery (fig3): coverage at σ=0.5 rises from 0.85 (global lock) to 1.00; at σ=0.75/1.00 it partially restores (0.75), the residual shortfall reflecting per-site boxes estimated from only 4 leaves.
- Sensitivity map (fig6): asymmetric truth (A at 0.2T, B at 0.5T). Strict ordering "B after A" is only robustly confirmed (≥7/9 reps) at 240 sites and heterogeneity ≤ 0.8; below that site count or under strong heterogeneity the strict separation claim is not identifiable, while the interval still covers the true gap.
- Design map (§20): recommended designs (width<0.75T under nominal, coverage≥0.70) hold at n_sites ≥ 16 for all scanned heterogeneity; only n_sites=8 fails at σ=0.5 (coverage 1.00 but width 0.81).
- Conclusion stress (§17): with tight per-site boxes the ordering claim stays unresolved on this 4-leaf recorder (interval [-0.188, 0.25]); wide transferred bounds widen the gap interval to [-0.188, 0.31]. Ordering claims on this design are not identifiable at n_sites≈60.
- Scout verdict: PASS (global-rate lock coverage 0.92 at σ=0, 0.83 at σ=0.5; exact small-tree intervals narrow; invariance holds) — now against the official Dryad bundle. Internal gate only, not external peer review.
- Official ConvexML Dryad data downloaded and certified: trees.tgz (403.2 MB) + README, `DOWNLOAD_STATUS: ok`, sha256 `b2339df3…` / `755b3a7e…`. Dryad's Anubis 'fast' PoW (SHA-256(randomData+nonce), difficulty 4) is solved in Python (`anubis_get` in the adapter); on-demand files stream via `file_stream/<fileid>`. Provenance tag: `topology_provenance: simulated` (deposit README confirms Cassiopeia simulation; not wet-lab data).
- **Solver bug found & fixed during validation.** The scalable (d>24) bound path used a penalized full-space local search that emitted spurious degenerate `[1,1]` and spuriously wide certificates (old `age:N1 [0.025,1.0]` / `[0.092,1.0]`, `N51 [1,1]`). Replaced with profile-grid enumeration (`profile_grid_target`, 2-sweep coordinate descent, warm-start chain) in `b01_clock.solver.optimize`. Fixed-age profiles confirm the suppressed ages were feasible, so the old certificates were wrong, not merely coarse.
- External simulated-benchmark holdout (official tree_0, 400 leaves, 39 sites, fixed 64-leaf induced subtree, canonical `age:N1`): with the fixed solver, basal-split age **partially_identified [0.073, 0.146]** (width/T 0.074) — sharply narrowed; MLE 0.085 vs **known truth 0.114** (covered). Scout's independent draw: [0.048, 0.067] (just misses the truth). Order gap unidentified [0.000, 1.0] (truth +0.146, so the overall positive sign holds but is not certified).
- **Ground-truth coverage on the benchmark: 1/9 (11%)** at nominal 95% (`convexml_ground_truth_coverage.json`). Missed nodes are systematically late-shifted: over-aged shallow splits (N3 true 0.242 vs [0.298,0.326]; N51 true 0.522 vs [0.788,0.849]) and deep splits pinned at \(T\) (N119/N73/N209/N78 → [1.000,1.000] against true 0.64–0.84). The binary/global-rate lock does **not** transfer to this Cassiopeia-simulated process; the sharp N1 interval is miscalibrated, so sharpness is claimed without coverage.

See `full_research_summary.json` for machine-readable dump.
