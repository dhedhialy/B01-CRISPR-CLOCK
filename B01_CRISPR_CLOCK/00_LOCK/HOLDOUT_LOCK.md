# Holdout Lock — Leakage-Free Site Splits

**Status:** LOCKED (pre-registration)  
**Package mirror:** `b01_clock.validation.calibration`  
**Purpose:** Prevent rate-bound estimation from contaminating headline timing certificates.

---

## 1. Unit of split

Splits are over **recording sites** (columns of the character matrix), not over cells/leaves and not over trees.

- Calibration sites \(\mathcal{S}_{\mathrm{cal}}\) estimate transferable global rate bounds \([r_{\min},r_{\max}]\).
- Analysis sites \(\mathcal{S}_{\mathrm{an}}\) alone enter the profile likelihood for \(\mathcal{I}_v(D)\).
- \(\mathcal{S}_{\mathrm{cal}}\cap\mathcal{S}_{\mathrm{an}}=\emptyset\), \(\mathcal{S}_{\mathrm{cal}}\cup\mathcal{S}_{\mathrm{an}}=\{0,\ldots,S-1\}\) (up to documented drop of empty sites).

**Rationale.** Site rates are the nuisance parameters that trade with time. Using the same sites to learn \(r\)-bounds and to time the tree creates circular narrowing (leakage).

---

## 2. Pre-registered split rule

### 2.1 Default fractions

| Parameter | Locked default | Allowed sensitivity |
|-----------|----------------|---------------------|
| `calib_frac` | **0.40** | \(\{0.30,0.50\}\) only in supplement |
| RNG | `numpy.random.default_rng(seed)` | seed recorded in split JSON |
| Headline seed (scout) | config `seed` + 99 (as in `run_decisive_scout`) | must appear in certificate |

Algorithm (`split_sites`):

1. Draw a uniform random permutation \(\pi\) of \(\{0,\ldots,S-1\}\).
2. \(n_{\mathrm{cal}}=\max\bigl(1,\mathrm{round}(\texttt{calib\_frac}\cdot S)\bigr)\).
3. \(\mathcal{S}_{\mathrm{cal}}=\mathrm{sort}(\pi[:n_{\mathrm{cal}}])\), \(\mathcal{S}_{\mathrm{an}}=\mathrm{sort}(\pi[n_{\mathrm{cal}}:])\).
4. If \(\mathcal{S}_{\mathrm{an}}=\emptyset\) (degenerate tiny \(S\)), set \(\mathcal{S}_{\mathrm{an}}\leftarrow\mathcal{S}_{\mathrm{cal}}\) and **flag** `degenerate_split: true` (abstain from headline claims).

### 2.2 Lock-before-headline

Before any call to `exact_profile_bounds` / scalable certificates on analysis sites, write:

```json
{
  "calib_sites": [...],
  "analysis_sites": [...],
  "seed": <int>,
  "locked_before_headline": true,
  "note": "..."
}
```

via `lock_split_record(...)`.  
**Forbidden:** re-drawing the split after inspecting \(\mathcal{I}_v\); cherry-picking sites with high edit diversity post hoc.

---

## 3. Calibration → rate bounds (only)

From calibration sites only, estimate a shared box for global-rate profiling:

\[
\hat r_s \;=\; \frac{-\log(1-\hat p_s)}{T},
\qquad
\hat p_s \;=\; \frac{\#\{\text{edited leaves at }s\}}{\#\{\text{observed leaves at }s\}}.
\]

Locked transfer (`held_out_rate_bounds`):

\[
r_{\min}
=\max\Bigl(r_{\mathrm{floor}},\;
Q_{0.1}(\{\hat r_s\}_{s\in\mathcal{S}_{\mathrm{cal}}})/\texttt{pad}\Bigr),
\quad
r_{\max}
=\max\bigl(r_{\min}\cdot 1.01,\;
Q_{0.9}(\{\hat r_s\})\cdot\texttt{pad}\bigr),
\]

with defaults \(\texttt{pad}=1.1\), \(r_{\mathrm{floor}}=10^{-3}\), quantiles \((0.1,0.9)\).

**These bounds may not be recomputed on \(\mathcal{S}_{\mathrm{an}}\).** Sensitivity may recompute under alternate quantiles/pads and report interval expansion.

---

## 4. Headline analysis (only)

On \(\mathcal{S}_{\mathrm{an}}\) with frozen \([r_{\min},r_{\max}]\) and dated \(T\):

1. Profile LR feasible set \(\widehat{\Theta}_\alpha\).
2. Project to \(\mathcal{I}_v\) for pre-registered targets (scout: `age:A`, `order_gap:A,B`; manuscript: list in certificate).
3. Emit identification status + topology ensemble union if topology uncertainty is on.

---

## 5. Multi-dataset rules

| Dataset | Split inventory |
|---------|-----------------|
| ConvexML | Sites = columns after adapter encoding; split per analyzed matrix (certified on fixed 64-leaf induced subtree of official Dryad tree_0) |
| SciPhy / PALINCODE | Same; if multi-state, define irreversible collapse **before** split and freeze collapse map |
| intMEMOIR | Same; do not split leaves |
| Synthetic scout | Same protocol on simulated matrix |

Cross-dataset: **never** transfer \([r_{\min},r_{\max}]\) from ConvexML calibration to SciPhy analysis without an explicit “external bound transfer” sensitivity label.

---

## 6. Leakage checklist (must all pass)

- [ ] Split JSON exists with `locked_before_headline: true` and timestamp/seed.
- [ ] Rate bounds provenance cites **only** `calib_sites`.
- [ ] Likelihood / certificate cites **only** `analysis_sites`.
- [ ] No site appears in both lists.
- [ ] Grid refinement / optimizer restarts do not reshuffle sites.
- [ ] Topology ensemble uses the same locked split and bounds.

Failure of any item → **abstain** from that headline interval (`status: abstain`).

---

## 7. What holdout does *not* claim

Site holdout is not a substitute for biological replicates, nor a cross-validation of tree topology. It only blocks **nuisance-rate leakage** into timing certificates.
