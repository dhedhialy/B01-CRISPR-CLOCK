# When Is a CRISPR Lineage Recorder a Clock? Sharp Timing Sets under Heterogeneous Editing Rates

**Working paper draft (B01)**  
**Status:** Substantive draft; §5 empirical cells filled from frozen artifacts (`09_VALIDATION/*`, `07_REAL_DATA/*`)  
**Software:** `b01_clock`  
**Locks:** `00_LOCK/*`, `02_THEORY/*`

---

## Abstract

CRISPR lineage recorders accumulate irreversible edits that are often treated as molecular clocks for developmental timing. We argue that this practice conflates two distinct objects: **exposures** \(r_s t_e\), which are identified by edit likelihoods, and **absolute times** \(t_e\), which are not. Under the irreversible model \(p_{se}=1-\exp(-r_s t_e)\), the likelihood is invariant to the rate–time rescaling \((r,t)\mapsto(c r,\, t/c)\). Absolute chronology therefore cannot be point-identified from barcodes alone. We formalize the estimand as the **sharp projection** of a profile likelihood-ratio feasible set under locked global rate bounds and a dated experimental anchor \(T\), and classify each timing functional as unidentified, partially identified, or point-identified. On small ultrametric trees we compute exact grid certificates; we add rate–time guardrails against false absolute identification, topology-ensemble unions, adversarial simulations, and leakage-free site holdout for bound transfer. A decisive synthetic scout yields nontrivial partial intervals (e.g. internal age \(A\) materially inside \([0,T]\)) with near-nominal coverage and zero false point-identification. On an **independently simulated external benchmark** (the ConvexML Dryad deposit, a Cassiopeia-generated tree with known ground-truth branch lengths) the primary lock's certificates are sharp but **miscalibrated** — ground-truth coverage 1/9, with a systematic late shift of deep splits (§5.4) — so the benchmark result is reported as a coverage failure and no wet-lab source is certified. We position the framework against ConvexML-scale reconstruction tools and recorder chronometry efforts: the contribution is not a better clock estimator but honest **partial identification of compatible timing sets** under heterogeneous editing rates.

---

## 1. Introduction

Lineage-recording systems based on CRISPR (and related integrase / barcode architectures) convert developmental history into mutational patterns on synthetic target arrays. Because edits accumulate over time, it is natural to ask whether recorders function as **clocks**: devices that recover when splits occurred, not only who shares a common ancestor.

That question is ambiguous. Tree reconstruction from barcodes is a mature computational problem: given leaf alleles, recover a genealogy that explains shared and private scars. Timing—assigning absolute or even relative ages to internal nodes—is a different statistical object. The likelihood of an irreversible edit on site \(s\) along branch \(e\) depends on the product of an editing rate and a duration. Heterogeneity in rates across sites, branches, and cellular contexts is the rule rather than the exception. Whenever rates are uncertain, time inherits that uncertainty. Treating a maximum-likelihood chronogram as “the” developmental schedule silently converts that uncertainty into false precision.

Classical molecular-clock literature already warns that substitution rates and times trade off, and Bayesian relaxed clocks encode the tradeoff through priors. CRISPR recorders sharpen the issue in two ways. First, the editing process is typically modeled as **irreversible**, so path exposures enter through \(1-e^{-rt}\) rather than reversible substitution matrices. Second, experimental designs often supply a single hard dated anchor—the harvest or induction duration \(T\)—while leaving site rates only weakly constrained by pilot assays or edit fractions. In that constraint regime, the right inferential object is not a point age but the **set of ages still compatible with the data**.

This paper takes the ambiguity seriously. We do **not** propose another point chronogram. We ask: *under what locked constraints does a recorder identify a timing functional, and what is the sharp set of timings still compatible with the data?* The answer is set-valued. Sometimes the compatible set for an internal age is essentially \([0,T]\); sometimes it is a proper subinterval; only under strong collapse conditions does it shrink to a point. We treat those three regimes as first-class scientific outcomes, not as solver pathologies.

Adjacent computational work makes the positioning delicate. ConvexML-scale tools optimize structured likelihoods for reconstruction and timing on public Dryad releases (in B01's case, a **simulated** deposit generated with Cassiopeia — we make no wet-lab-data claim). SciPhy and related GEO cohorts expand sequential recorder architectures. Chronometry-framed recorder papers (including efforts in the spirit of Pilarski 2026) emphasize improved age estimation error. Our wedge is orthogonal: we reuse public simulated benchmark substrates where useful, but we change the estimand to **sharp partial identification** and we instrument the pipeline against false absolute-time claims via rate–time invariance diagnostics.

**Contributions.**

1. **Estimand lock.** Timing certificates \(\mathcal{I}_v(D)\): profile LR projections for node ages, branch lengths, and order gaps (`ESTIMAND_LOCK.md`).  
2. **Identification conjecture.** Trichotomy unidentified / partial / point under rate–time invariance, dated \(T\), and global rate bounds — proven on the canonical 4-leaf tree, conjectured beyond it (`IDENTIFICATION_THEOREM.md`).  
3. **Exact small-tree solver** with invariance guardrails, adversarial battery, topology ensemble, and pre-registered site holdout (`b01_clock`).  
4. **Positioning.** Collision analysis vs ConvexML, SciPhy, and chronometry-framed work (`NOVELTY_COLLISION.md`), plus a referee attack map (`REFEREE_TABLE.md`).  
5. **Certificate schema.** Machine-readable JSON contracts for interoperable timing sets (`CERTIFICATE_SCHEMA.md`).

**What we do not claim.** Biological discovery of a developmental event’s wall-clock time beyond the certificate; superiority on reconstruction leaderboards; universal absolute dating without anchors; that every recorder experiment yields narrow intervals.

---

## 2. Model

### 2.1 Irreversible exposure

On site \(s\) and branch \(e\),
\[
p_{se}=1-\exp(-r_s t_e).
\]
States are unedited (\(0\)), edited (\(1\)), or missing. Transitions \(1\to 0\) are forbidden and contribute \(-\infty\) to the log-likelihood. Independent sites multiply. Ancestral states along the tree are prepared under irreversibility before scoring branches. The software keeps **exposure** \(\lambda=rt\) distinct from absolute time at every layer (`b01_clock.model.exposure`), including helpers that map edit fractions to exposures and exposures through rate boxes to time bounds without collapsing the distinction.

### 2.2 Ultrametric dating lock

Leaves are fixed at dated age \(T\) (experiment duration in absolute units, or \(T=1\) with all reports in units of \(T\)). Free parameters are ages of non-root internal nodes in \([0,T]\). Branch durations are nonnegative differences of ages, recovered by `ultrametric_times_from_free`. The feasible ambient space \(\mathcal{U}(T)\) is thus a compact polytope whose dimension equals the number of free internals.

### 2.3 Rate region

Primary analyses use a **global** rate \(r\in[r_{\min},r_{\max}]\), profiled for each candidate timing. This is a deliberate modeling lock: freely varying per-site rates reintroduce a near-continuous family of \(\Psi_c\)-like trades and can restore Regime U for absolute ages. Site-specific free rates are reserved for sensitivity and adversarial axes. Rate bounds for headline intervals are estimated from **calibration sites only** using leaf edit fractions and dated \(T\), then frozen for analysis (`HOLDOUT_LOCK.md`).

### 2.4 Observation nuisances

Dropout and silencing enter as an observation layer that can mask edited or unedited states. They are not part of the locked estimand’s structural equation; they are adversarial axes that should widen certificates or trigger abstention when information collapses.

### 2.5 Targets

We certify functionals \(v\in\{\texttt{age:}u,\;\texttt{branch:}p{\to}c,\;\texttt{order\_gap:}a{,}b\}\). Ages and branches speak to absolute chronology in units of \(T\). Order gaps test whether temporal order is identified even when absolute ages remain interval-valued—an important distinction for developmental questions that only need “\(A\) before \(B\).”

---

## 3. Identification theory (summary)

### 3.1 Rate–time invariance

For any \(c>0\),
\[
L(D\mid\mathbf{r},\mathbf{t})=L(D\mid c\mathbf{r},\,\mathbf{t}/c).
\]
Equivalently, \(t\mapsto c t\) with \(r\mapsto r/c\) traces the same orbit. Absent constraints that break this group action, absolute times are unidentified: the data see exposures, not clocks. Dated leaves at \(T\) break nontrivial rescalings because \(\mathbf{t}/c\notin\mathcal{U}(T)\) for \(c\neq 1\). A compact rate box breaks them when \(cr\) exits \([r_{\min},r_{\max}]\). Computationally we verify likelihood invariance and reject false absolute identification when \(\Psi_c\) exits the locked constraint class (`INVARIANCE.md`).

### 3.2 Trichotomy

**Conjecture (informal; proven on the 4-leaf scout tree).** With free rates and free depth, absolute chronology is unidentified (Regime U). With dated \(T\) and a nondegenerate global rate box, node ages and branch lengths are typically **partially identified** (Regime P): the sharp interval \(\mathcal{I}_v\) is a proper subset of the trivial bound but not a singleton. Point identification (Regime I) requires collapse—known rate, degenerate box, or additional anchors with sufficient rank (`IDENTIFICATION_THEOREM.md`). Under the primary lock, claiming Regime I for internal ages is treated as a false-identification failure mode unless collapse is documented.

On the canonical 4-leaf scout tree with free ages \((a_A,a_B)\), root-to-leaf exposure is always \(rT\); clade splits identify exposure **ratios** such as \(a_A/(T-a_A)\). Profiling \(r\) on a compact interval truncates timings that would require out-of-box rates to match observed edit patterns, producing proper projections for `age:A`. Topology swaps reassign which shared edits belong to which clade; ensemble unions expand intervals honestly under topology uncertainty.

### 3.3 Geometry and sharpness

The profile set \(\widehat{\Theta}_\alpha\subset\mathcal{U}(T)\) is the feasible blob in age space; \(\mathcal{I}_v\) is its projection onto \(\vartheta_v\) (`FEASIBLE_SET_GEOMETRY.md`). Among intervals covering the image of \(\widehat{\Theta}_\alpha\), \(\mathcal{I}_v\) is shortest. Topology uncertainty is handled by union over an audited ensemble; orders with constant sign across topologies are flagged `robust_order`. Width \(w_v/T\) is reported as scientific content, not merely as uncertainty decoration on a point estimate.

### 3.4 Finite-sample certificates

Operationally we use one-dimensional profile LR cutoffs \(\chi^2_{1,1-\alpha}\) after profiling rates, then validate coverage empirically under simulation. On \(d\le 3\) free internals, grid attainment yields conservative staircased outer bounds, checked by refinement.

---

## 4. Methods

### 4.1 Exact profile certificates

For trees with at most three free internals, we grid \(\mathcal{U}(T)\), profile \(r\) in the locked box at each grid point, and retain timings with
\[
2\bigl(\ell_{\max}-\ell_{\mathrm{prof}}(\mathbf{t})\bigr)\le\chi^2_{1,1-\alpha}.
\]
Interval endpoints are min/max of \(\vartheta_v\) on the retained set. Status is assigned by width relative to \(T\) (`unidentified` / `partially_identified` / `point_identified`). Grid refinement checks stability within a fraction of \(T\) (default tolerance \(0.05\,T\)). Larger trees are deferred to scalable outer methods that must be validated against the exact solver on shared small instances.

### 4.2 Invariance and false-ID diagnostics

We confirm \(\ell(r,t)=\ell(cr,t/c)\) numerically and apply the feasibility guardrail for \(c\neq 1\). A run that returns Regime I while a nontrivial rescale remains inside \(\mathcal{C}\) is treated as a software or conceptual failure. Adversarial simulations vary site heterogeneity, branch-rate CV, dropout/silencing, stochastic division, topology perturbation, and rate-bound misspecification (`b01_clock.simulate.adversarial`). Primary false-ID metric: frequency of spuriously reporting point identification under rescaling or under settings theory marks as Regime P.

### 4.3 Topology ensemble

Identity plus leaf swaps across clades form a finite audited ensemble (not a posterior over all trees). For each target we report the union hull of intervals and `robust_order` for gap functionals. Side-by-side identity vs union widths diagnose whether topology uncertainty dominates timing (kill criterion K2).

### 4.4 Leakage-free holdout

Pre-register a 40/60 site split (default; sensitivities at 30/70 and 50/50). Calibration sites yield quantile-based \([r_{\min},r_{\max}]\) with padding; analysis sites alone enter \(\ell_{\mathrm{prof}}\). Split records set `locked_before_headline: true` before certificate computation. Cross-dataset transfer of bounds is labeled sensitivity, never silent headline.

### 4.5 Data sources

Locked sources: ConvexML (Dryad `10.5061/dryad.qrfj6q5nz`, primary — a **simulated** benchmark by the deposit's own README, used as an external generated dataset, not as wet-lab data), SciPhy (`GSE315827`, pending), PALINCODE (`GSE327634`, pending), intMEMOIR (`SOURCE_LOCK.md`, pending). Synthetic scout trees use the canonical 4-leaf ultrametric topology. Provenance records store accession, SHA-256, and adapter version; provenance must never relabel simulated data as real. **Acquisition note:** Dryad's bot-protection (Anubis JS proof-of-work) was satisfied programmatically in `anubis_get` solely to complete the public download; we disclose the mechanism in the methods and run `download_convexml` idempotently with sha256 verification.

### 4.6 Pass / kill and reporting

Decisive scout PASS requires material narrowing of at least one \(\mathcal{I}_v\), synthetic coverage near nominal, invariance OK, and false-ID control. Kill if intervals remain \(\approx[0,T]\), topology dominates, calibration fails to transfer, or the story collapses to a generic clock MLE (`ESTIMAND_LOCK.md`). All headline intervals emit JSON conforming to `CERTIFICATE_SCHEMA.md`.

---

## 5. Results

*All numbers are pulled from frozen artifacts: `results/scout/decisive_scout_report.json`, `B01_CRISPR_CLOCK/09_VALIDATION/full_research_summary.json`, `B01_CRISPR_CLOCK/07_REAL_DATA/convexml_protocol_reanalysis.json`, and `convexml_ground_truth_coverage.json`. Reproduce with `scripts/reproduce.sh`. Empirical checks are internal go/no-go gates, not external peer review.*

### 5.1 Decisive synthetic scout

On the 4-leaf scout tree with dated \(T=1\) and global rate bounds from calibration sites, exact certificates classify internal ages as **partially identified**. Representative profile interval for `age:A` (scout config, 40 sites, 12 coverage reps):

- \(\mathcal{I}_{\mathrm{age:}A}=\) **[0.200, 0.733]** (MLE 0.533; width/T = 0.533), status **partially_identified**.
- \(\mathcal{I}_{\mathrm{age:}B}=\) [0.167, 0.733]; \(\mathcal{I}_{\mathrm{branch:}root\to A}=\) [0.200, 0.733]; both partially identified.
- Grid verification OK (coarse [0.200, 0.733], fine [0.200, 0.750], abs diff ≤ 0.017).
- Feasible profile-LR set covers 242/961 grid cells (free internals \(A,B\)), cutoff \(\chi^2_{1,0.95}=1.921\).

Order gap `order_gap:A,B` = **[-0.433, 0.400]** — crosses zero: even where ages partially identify, the ordering of these two splits is not certified at this site count.

**Figures.** [`fig2_identification_geometry`](../10_FIGURES/fig2_identification_geometry.pdf) shows the feasible blob and projections in \((a_A,a_B)\) space.

### 5.2 Invariance and guardrail

Likelihood invariance under scale \(c=2\): **OK**, absolute log-likelihood difference = **0.0**. The rate–time guardrail reports `feasible_with_anchor_and_bounds = false` for the nontrivial rescale while `feasible_without_anchor = true` in the geometric analysis (fig1): the locked dated anchor \(T=1\) and rate box jointly break the rescaling orbit. False absolute identification is blocked. **Figure.** [`fig1_rate_time_confounding`](../10_FIGURES/fig1_rate_time_confounding.pdf).

### 5.3 Coverage and false identification

Scout coverage battery (12 reps, rate bounds frozen from calibration sites; nominal 95%):

| Rate mode | Heterogeneity | Coverage | Mean width / \(T\) | False-ID rate | Abstain |
|-----------|---------------|----------|--------------------|---------------|---------|
| global | 0.0 (low) | 0.917 | 0.568 | 0.0 | 0.0 |
| global | 0.5 (mid) | 0.833 | 0.469 | 0.0 | 0.0 |
| per-site | 0.75 (high) | 0.667 | 0.563 | 0.0 | 0.0 |

False point-identification rate **0.0 across all scanned regimes**, including strong per-site heterogeneity; per-site mode additionally abstains (0.125) rather than mis-certify. The full-research sweep (`full_research_summary.json` fig3) reproduces near-nominal global coverage 0.95 at \(\sigma=0\) and monotone decline under heterogeneity, with false-ID 0.0 throughout. **Figure.** [`fig3_synthetic_coverage`](../10_FIGURES/fig3_synthetic_coverage.pdf).

Scout verdict: **PASS** (narrowing, near-nominal coverage at low het, false-ID=0, invariance OK).

### 5.4 External benchmark (ConvexML Dryad): certificates and ground-truth coverage

The ConvexML Dryad deposit (`10.5061/dryad.qrfj6q5nz`) is an **independently generated, fully simulated benchmark**, not measured wet-lab data: per its README, trees and lineage-tracing matrices were simulated with the Cassiopeia package. Its newick branch lengths are therefore **known ground-truth node times**, which admits a real, pre-registered coverage check on an external (third-party-generated) dataset.

**Basal certificate** (site split seed 7, 16 calibration / 23 analysis sites, transferred bounds [0.176, 2.103], fixed random 64-leaf induced subtree, free internals 64):

- `age:N1` (root resolved split): **[0.073, 0.146]**, width/T = **0.074**, **partially_identified**.
- `order_gap:N1,N285`: **[0.000, 1.000]**, **unidentified**.
- Scout's independent draw (seed 106, bounds [0.158, 3.551]): `age:N1` = [0.048, 0.067] (width 0.018), `order_gap` = [0.177, 1.000].
- Ground truth: age(N1) = 0.114, gap = age(N285) − age(N1) = +0.146. The reanalysis interval covers N1; the scout interval just misses it.

**Coverage against known truth** (9 internal nodes spanning the age range, same induced subtree, `convexml_ground_truth_coverage.json`): **1/9 covered (11%)** at nominal 95%. The misses are **systematically late-shifted**: shallow splits are over-aged (N3: true 0.242 vs [0.298, 0.326]; N35: 0.339 vs [0.404, 0.458]; N51: 0.522 vs [0.788, 0.849]) and deep splits are pinned at the harvest time (N119 0.642, N73 0.691, N209 0.773, N78 0.840 → [1.000, 1.000]; N312 0.890 → [0.792, 0.831]).

**Solver note (audit-internal).** The first measurement pass exposed a bound-enumeration bug in the scalable (d>24) certificate path — a penalized full-space local search emitted spurious degenerate `[1,1]` and spuriously wide intervals. It was replaced with profile-grid enumeration (`b01_clock.solver.optimize`, `profile_grid_target`), which the fixed-age profile diagnostics corroborate (e.g. N51 age 0.2–0.8 is reachable at high profile likelihood, ruling out the earlier `[1,1]`). The coverage verdict is quantitatively unchanged by the fix.

**Reading.** The binary irreversible + global-rate primary lock is **not calibrated on this simulated process**: coverage far below nominal, with a monotone late shift. This is consistent with the deposit's heritable silencing, dropout, and multi-state indel structure surviving the binary collapse — the model reads too little shared-edit signal in deep clades and moves their splits toward \(T\). The synthetic-scout coverage (§5.3) does **not** transfer. This is a decisive class-credentialing negative for the primary lock on Cassiopeia-style generated data, and it is why this section is called "external benchmark" rather than "real data."

### 5.5 Adversarial battery

`full_research_summary.json` fig3 enumerates **11/11 adversarial specs with the truth covered** (specs: baseline, site-heterogeneity mild/strong, branch heterogeneity, replicate shift, stochastic division, saturation, dropout, silencing, topology error, rate-bound misspecification). Saturation and strong heterogeneity push intervals toward [0,T] (width → 0.75–1.0) but never mis-certify a false point: false-ID **0.0**. The `rate_bound_misspec` clamp fix (bound inversion guarded) keeps the battery fully covered.

### 5.6 Topology ensemble

Identity vs leaf-swap ensemble (union over audited topologies, `fig5`): union `age:A` width **0.714** vs identity-only 0.500 — topology uncertainty dominates timing and `robust_order` is **false** for `order_gap:A,B` under the synthetic scout. Kill criterion K2 (topology dominates) is thereby *registered as active* on this synthetic tree: certificate widths must carry a topology-union caveat. **Figure.** [`fig5_topology_architecture`](../10_FIGURES/fig5_topology_architecture.pdf).

### 5.7 Sensitivity and design map

Heterogeneity × site-count phase map (`fig6`: `n_sites∈{30,60,120,240}`, heterogeneity ∈ [0,1], 3 seeds/cell): at low heterogeneity at least one of `confirm_frac`/`cover_frac` attains 1.0 already at 60 sites; at high heterogeneity only 240 sites routinely confirm a positive order gap (`confirm_frac 0.33–0.67`). Design map (`fig7`: 8–48 sites × 0–0.75 het): 32 sites + σ≤0.25 keep mean width ≤ 0.53 with coverage ≥ 0.88; 48 sites are robust to σ=0.75. **Figures.** [`fig6_sensitivity_map`](../10_FIGURES/fig6_sensitivity_map.pdf), [`fig7_design_map`](../10_FIGURES/fig7_design_map.pdf).

### 5.8 External benchmark certificates

**ConvexML (primary, Dryad) — simulated benchmark (not wet lab).** Bundle `trees.tgz` sha256 `b2339df3…` (data-manifest checksum `c92b30ec…`, `accession_ok: true`), 400 leaves / 799 edges / 39 sites, 20.5% missing, 41% saturation proxy, ultrametric with leaves dated at \(T\). Provenance records `topology_provenance: simulated`. Certificates and the ground-truth coverage result are in §5.4.

**intMEMOIR / SciPhy / PALINCODE.** Not yet run: public adapters for GSE315827 / GSE327634 are not part of this phase's pipeline. Consequently **no wet-lab recorder deposit is certified yet** — the external tier currently consists only of simulated benchmarks. These accessions remain registered as locked future work (Limitations §7.7; `SOURCE_LOCK.md` notes them pending).

### 5.9 Conclusion stress and comparator points

A head-to-head test of a "B after A" timing claim (`fig3` conclusion_stress) leaves the order gap interval [-0.188, 0.25] under tight per-site calibration and [-0.188, 0.312] under wide transferred bounds — **unresolved** both ways, i.e. the recorder at n_sites≈60 does not certify a 3-hour-style separation. No comparator point timings are emitted this phase (no external optimizer output integrated); ConvexML point estimates would be plotted as overlays without truth treatment once integrated.

---

## 6. Discussion

Lineage recorders are clocks only **sometimes**, and almost never in the naive point-estimate sense under heterogeneous editing. The scientifically honest object is the set of timings compatible with data once the analyst locks the constraints that break rate–time invariance. Partial identification is not a failure mode; it is the generic regime when a single dated harvest time and a plausible rate box are all that biology supplies. Point chronograms remain useful as interior probes of \(\widehat{\Theta}_\alpha\), but they should not be narrated as unique developmental schedules.

This reframes success metrics. A method that returns a precise age with miscalibrated coverage under heterogeneity is worse, for developmental inference, than a wide certificate that covers the truth and refuses false point claims. Width is content: it quantifies how much timing information the recorder actually carries under \(\mathcal{C}\). The external benchmark makes the point concretely (§5.4): with the bound-enumeration bug fixed, the scalable solver emitted **sharp** intervals whose coverage was nonetheless 1/9 (11%) — precise but miscalibrated, systematically over-aging deep splits. Order functionals add a second axis of usefulness: even when absolute ages remain wide, a robust sign on \(\mathrm{age}(B)-\mathrm{age}(A)\) can still answer “which clade split earlier?” under topology uncertainty.

Relative to ConvexML and related reconstruction stacks, we share likelihood ingredients and public data but not the estimand. A ConvexML point falling inside our interval is consistency under overlapping constraints; a point falling outside is a prompt to compare \(\mathcal{R}\), \(T\), and topology, not an automatic refutation. Relative to chronometry-framed recorder papers, we reject MSE-to-truth as the primary headline when truth is only partially identified. SciPhy and PALINCODE serve as architecture-transfer stresses: negative transfer widens the limitation boundary without erasing ConvexML or synthetic identification results.

The practical pipeline is intentionally conservative: calibrate rates on held-out sites, certify on analysis sites, union over audited topologies, emit machine-readable certificates with identification status, and map referee attacks to diagnostics before review (`REFEREE_TABLE.md`). Downstream biological narrative should cite \(\mathcal{I}_v\) and status, not a lone MLE. In that sense, B01 is infrastructure for honest timing claims rather than a replacement developmental model.

---

## 7. Limitations

1. **Exact solver scale.** Grid certificates are for small \(d\); large trees require scalable outer methods whose sharpness must be validated against the exact solver on shared subproblems.  
2. **Global rate mode.** Shared \(r\) is a deliberate lock against silent reintroduction of scale freedom; true site-specific rates are an adversarial/sensitivity axis.  
3. **Observation model.** Dropout and silencing are simulated; richer allele models (sequential multi-state) need documented likelihood extensions before Regime claims transfer verbatim. The external-benchmark coverage failure (§5.4) is consistent with this gap: the binary collapse under-reads shared edits in deep clades (late-shifted splits).  
4. **Topology.** Ensemble is finite and audited, not a posterior over all trees. Union can restore near-unidentification (K2).  
5. **\(\chi^2\) calibration.** Profile LR cutoffs are classical for point-identified focus parameters; under partial ID we rely on empirical coverage.  
6. **Dated \(T\).** If experimental duration is itself uncertain, that uncertainty must enter \(\mathcal{C}\) or intervals are conditional on a fixed \(T\).  
7. **Wet-lab ground truth and benchmark status.** The primary external deposit is itself **simulated** (Cassiopeia), so its truth is known — and the coverage check (§5.4) is a mandatory reading, currently **failing (1/9)** for the binary/global-rate lock. No wet-lab recorder time-series (intMEMOIR / SciPhy / PALINCODE) is certified yet; those are registered, not claimed.  
8. **Scope.** No wet-lab or clinical claims; no priority claim over scalable reconstruction algorithms.  
9. **Missingness.** Informative missingness is demonstrably damaging: the deposit's silencing/dropout are the most plausible driver of the onboarded benchmark's coverage failure (§5.4).

Referee-facing attack map: `11_MANUSCRIPT/REFEREE_TABLE.md`.

---

## 8. Conclusion

A CRISPR lineage recorder is a clock only to the extent that locked constraints break rate–time invariance and residual timing freedom collapses for the functional of interest. We replace point chronogram rhetoric with sharp certified timing sets under heterogeneous irreversible editing, with explicit unidentified / partial / point regimes. The decisive synthetic scout demonstrates that nontrivial partial identification is achievable with calibrated coverage and controlled false identification—the empirical signature that the framework is not vacuous — and the external simulated benchmark shows the same machinery refusing to overclaim where it is out of calibration: sharp intervals that fail to cover known truth are reported as the coverage failure they are (§5.4). Frozen artifacts in `B01_CRISPR_CLOCK/09_VALIDATION` and `07_REAL_DATA` lock both signatures as reproducible manuscript claims. The broader invitation to the field is simple: when publishing recorder timings, publish the compatible set, the constraints that made it finite, and the diagnostics that show absolute time was not smuggled in through invariance.

---

## Acknowledgments / Data

ConvexML: Dryad DOI 10.5061/dryad.qrfj6q5nz (official bundle downloaded, sha256 `b2339df3…`). SciPhy: GEO GSE315827. PALINCODE: GEO GSE327634. intMEMOIR: public processed releases as provenance-recorded. Software: `b01_clock`.

---

## Appendix A — Notation

| Symbol | Meaning |
|--------|---------|
| \(T\) | Dated leaf age / experiment duration |
| \(r_s, t_e\) | Site rate, branch time |
| \(\lambda_{se}=r_s t_e\) | Exposure |
| \(\mathcal{U}(T)\) | Ultrametric timings with leaf depth \(T\) |
| \(\widehat{\Theta}_\alpha\) | Profile LR feasible set |
| \(\mathcal{I}_v\) | Sharp interval for functional \(v\) |
| \(\Psi_c\) | Rate–time rescaling |

## Appendix B — Certificate pointer

Machine-readable schema: `12_CERTIFICATES/CERTIFICATE_SCHEMA.md`.


## Filled empirical results (auto)

Primary pipeline (`scripts/run_full_research.py`, scout config; artifacts in `09_VALIDATION/`, `results/`):

- Exact small-tree age(A): [0.000, 0.625] (MLE 0.375; partially_identified) — scout config exact cert: [0.200, 0.733].
- Exact small-tree age(B): [0.375, 0.792].
- Feasible-set fraction of (age A, age B) grid: 0.229.
- Synthetic coverage: σ=0 → 0.95; σ=0.25 → 0.90; σ=0.5 → 0.85; false-ID 0 across all scanned regimes (full sweep); scout battery (12 reps): 0.917 / 0.833 / per-site 0.667.
- Synthetic holdout age(A): [0.333, 0.792] (width/T=0.458) with transferred bounds [0.262, 15.197] (4-leaf synthetic; `08_TOPOLOGY` certificate).
- External benchmark holdout (official Dryad tree_0 — **simulated Cassiopeia deposit**, 400 leaves, 39 sites, fixed 64-leaf induced certification): age(N1) [0.073, 0.146] partially_identified; order gap unidentified [0.000, 1.000]; ground-truth coverage 1/9 (§5.4). Scout pre-fix artifacts retained only as bug records.
- Adversarial battery: 11/11 specs cover truth (incl. site-het strong, saturation, topology error, rate-bound misspecification); false-ID 0.
- Topology union width for age(A): 0.714; robust_order false.
- Conclusion stress (n_sites≈60, "B after A"): order gap [-0.188, 0.250] tight, [-0.188, 0.312] wide transferred → unresolved both.
- Sensitivity phase map: ≥60 sites + low het certify/confirm positive order; 480-site grid dropped (coverage collapse = peaked-LR grid artifact).
- ConvexML Dryad: official `trees.tgz` sha256 `b2339df3…` downloaded & certified (Anubis 'fast' PoW solved in `anubis_get`); `DOWNLOAD_STATUS.ok: true`; audit checksum `c92b30ec…`, `accession_ok: true`.
