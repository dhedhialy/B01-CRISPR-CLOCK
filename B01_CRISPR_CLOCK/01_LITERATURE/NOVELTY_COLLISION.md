# Novelty Collision — What Adjacent Work Owns vs Our Wedge

**Status:** Literature lock for positioning  
**Wedge (one line):** Sharp **partial identification of compatible timing sets** under heterogeneous irreversible editing rates — not a better clock point estimator.

---

## 1. Collision map

| Work | Core contribution | Assumes / optimizes | Collision risk | Our differentiation |
|------|-------------------|---------------------|----------------|---------------------|
| **ConvexML** (Dryad `10.5061/dryad.qrfj6q5nz`) | Scalable convex / structured ML for lineage reconstruction & timing under recorder likelihoods | Often seeks **point** timings / scalable MAP–MLE; rate heterogeneity handled as model feature for fit quality | High on data + likelihood form \(1-e^{-rt}\) | We reuse ConvexML **data**, not their estimand. We output **certified sets** \(\mathcal{I}_v\) with explicit unidentified / partial / point regimes and rate–time **false-ID guardrails** |
| **SciPhy** (GEO `GSE315827`) | Phylogenetic / sequential lineage-recorder analysis pipelines | Tree + substitution / sequential edit models; emphasis on reconstruction accuracy | Medium: sequential states vs binary irreversible | Transfer test only; our theory is exposure-scale identification under locked \(\mathcal{R},T\). Multi-state extensions are sensitivity, not the locked estimand |
| **Pilarski et al. 2026** (lineage-recorder chronometry / related clock framing — cite final venue when frozen) | Improved molecular / recorder **clock estimation** or calibration under biological constraints | Point chronograms, calibration tricks, or variance reduction for age estimates | **Highest narrative collision:** “we also time CRISPR recorders” | We **refuse** the “better clock” frame. Heterogeneous rates make absolute time a **set-valued** object; when the set is wide we say so. Success = honest \(\mathcal{I}_v\), not smaller MSE vs truth when truth is unidentified |
| Classical molecular clocks (Zuckerkandl–Pauling → BEAST-style) | Substitution clocks with relaxed rates | Prior-driven dating; often continuous-time Markov substitutions | Low on CRISPR barcode mechanics; high on “relaxed clock” language | CRISPR irreversibility + barcode saturation change the likelihood; our invariance is \(r\leftrightarrow t\) on **exposure**, not CIR/UCLN priors |
| Cassiopeia / PRESCIENT / related lineage tools | Tree reconstruction from barcodes | Topology first; timing secondary or heuristic | Medium if readers equate trees with clocks | Topology is an **input constraint class**; timing certificates may **union** over ensembles rather than pick one chronogram |

---

## 2. ConvexML — detailed wedge

**They own:** Computational scalability; convex relaxations; empirical reconstruction benchmarks on Dryad release.

**We own:**  
1. Formal three-regime identification taxonomy under global rate bounds + dated \(T\).  
2. Profile LR **sharp projections** as the estimand.  
3. Explicit test that rate–time rescaling does not yield false absolute identification.  
4. Leakage-free site holdout for bound transfer.

**Non-claim:** B01 is not “ConvexML but sharper ages.” If ConvexML reports a point age \(a^\star\) outside our \(\mathcal{I}_v\), that is a **model/constraint mismatch** to diagnose, not automatic refutation.

---

## 3. SciPhy — detailed wedge

**They own:** Sequential recorder biology and GEO-scale processing.

**We own:** Whether sequential architectures, after irreversible collapse or documented multi-state likelihood, still yield **nontrivial** \(\mathcal{I}_v\) under locked \(\mathcal{C}\).

**Collision avoidance text (for intro):** “SciPhy provides an external sequential-recorder cohort; we ask identification questions, not reconstruction leaderboards.”

---

## 4. Pilarski 2026 — detailed wedge

Treat Pilarski 2026 as the nearest **rhetorical** neighbor: papers that sell CRISPR lineage recorders as clocks with improved estimators.

| Their typical success metric | Our success metric |
|------------------------------|--------------------|
| Low error of \(\hat t\) vs simulated truth under homogeneous or mildly heterogeneous rates | Coverage of true \(\vartheta_v\) by \(\mathcal{I}_v\); **false-ID rate** under adversarial heterogeneity; width as scientific content |
| Tighter posterior / CI under a point-identified model | Correct **expansion** of intervals when rates are free or bounds weak |
| “Recorder \(X\) enables developmental timing” | “Recorder \(X\) identifies timing **only when** rank/anchor/topology conditions hold; otherwise abstain” |

**If Pilarski reports point IDs where our theory predicts partial ID:** either their constraints are stronger than advertised (must map into \(\mathcal{C}\)), or they are estimating a different functional (e.g., relative order with extra anchors).

---

## 5. Novelty claims we are allowed to make

1. **Estimand:** \(\mathcal{I}_v(D)=\mathrm{proj}_v\widehat{\Theta}_\alpha(D)\) under irreversible exposure + locked rate box + dated \(T\).  
2. **Theorem:** trichotomy unidentified / partial / point with explicit collapse conditions (global rate bound + \(T\) + topology).  
3. **Invariance guardrail** as a first-class computational diagnostic against false clocks.  
4. **Topology-ensemble union** for robust order claims.  
5. **Pre-registered site holdout** for rate-bound transfer.

## 6. Novelty claims we are forbidden to make

- “State-of-the-art accuracy” chronogram MSE.  
- Universal absolute dating without anchors.  
- Biological discovery of a developmental event time beyond \(\mathcal{I}_v\).  
- Priority over ConvexML for scalable reconstruction algorithms.

---

## 7. One-paragraph positioning (copy-ready)

CRISPR lineage recorders accumulate irreversible edits whose likelihood depends on exposures \(r_s t_e\). Prior tools—including ConvexML-scale optimizers and recorder chronometry efforts such as Pilarski 2026—primarily deliver point trees or point timings. Because exposures are invariant to rate–time rescaling, absolute chronology is not a classical point parameter under heterogeneous rates. We therefore target the **sharp set of timings compatible with data under explicit rate bounds and a dated anchor**, classify each functional as unidentified, partially identified, or point-identified, and certify intervals with profile likelihood ratios. Adjacent datasets (ConvexML, SciPhy, PALINCODE, intMEMOIR) stress-transfer this identification framework rather than re-litigate reconstruction accuracy.
