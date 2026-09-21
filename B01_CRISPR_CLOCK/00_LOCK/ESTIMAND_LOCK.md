# Estimand Lock — B01 CRISPR Clock

**Status:** LOCKED  
**Package mirror:** `b01_clock.solver.exact`, `b01_clock.solver.certificates`  
**Date lock:** 2026-09-08

---

## 1. Scientific claim (one sentence)

Under an irreversible recorder model \(p_{se}=1-\exp(-r_s t_e)\) with heterogeneous site rates, we do **not** claim a better chronogram MLE; we claim **sharp partial identification** of the set of timings compatible with data, global rate bounds, a dated anchor \(T\), and (when used) an audited topology ensemble.

---

## 2. Locked estimand

### 2.1 Data and parameters

Let \(\mathcal{T}=(V,E)\) be a rooted ultrametric tree with leaf set \(\mathcal{L}\) and root \(\rho\). Let \(T>0\) be a **dated anchor** (experiment duration, harvest time, or other absolute scale lock) with all leaves at age \(T\).

Site rates \(\mathbf{r}=(r_s)_{s=1}^{S}\), \(r_s>0\). Branch durations \(\mathbf{t}=(t_e)_{e\in E}\), \(t_e\ge 0\), consistent with ultrametricity to \(T\).

Observation layer \(D\): irreversible \(0/1\) (or missing) leaf states for each site, with ancestral states reconstructed under irreversibility (parent edited \(\Rightarrow\) child edited).

Exposure on site \(s\), edge \(e\):
\[
\lambda_{se}=r_s t_e,\qquad
p_{se}=1-\exp(-\lambda_{se}).
\]

### 2.2 Functional of interest

For a named target functional \(v\) (node age, branch length, or order gap), define
\[
\vartheta_v(\mathbf{t})\;=\;
\begin{cases}
\mathrm{age}(u), & v=\texttt{age:}u,\\
t_{\mathrm{parent}\to\mathrm{child}}, & v=\texttt{branch:}p{\to}c,\\
\mathrm{age}(b)-\mathrm{age}(a), & v=\texttt{order\_gap:}a{,}b.
\end{cases}
\]

### 2.3 Feasible timing set (population)

Given locked constraint class \(\mathcal{C}\) (rate box, dated \(T\), topology class), the **identified set** is
\[
\Theta(D;\mathcal{C})
\;=\;
\Bigl\{\,\mathbf{t}\;:\;
\exists\,\mathbf{r}\in\mathcal{R}\ \text{s.t.}\ 
\mathbf{t}\in\mathcal{U}(T),\ 
L(D\mid\mathbf{r},\mathbf{t})\ \text{is maximal under }\mathcal{C}
\,\Bigr\},
\]
where \(\mathcal{U}(T)\) is the ultrametric simplex of branch times with leaf depth \(T\), and \(\mathcal{R}\) is the locked rate region (default: global rate \(r\in[r_{\min},r_{\max}]\)).

### 2.4 Sharp timing interval (estimand)

The **locked estimand** for target \(v\) is the **projection** of the profile likelihood-ratio set:
\[
\boxed{
\mathcal{I}_v(D)
\;=\;
\Bigl[\,
\inf\bigl\{\vartheta_v(\mathbf{t}):\mathbf{t}\in\widehat{\Theta}_\alpha(D)\bigr\},\;
\sup\bigl\{\vartheta_v(\mathbf{t}):\mathbf{t}\in\widehat{\Theta}_\alpha(D)\bigr\}
\,\Bigr]
}
\]
where the finite-sample certificate set is
\[
\widehat{\Theta}_\alpha(D)
\;=\;
\Bigl\{\mathbf{t}\in\mathcal{U}(T):
2\bigl(\ell_{\max}-\ell_{\mathrm{prof}}(\mathbf{t})\bigr)
\le \chi^2_{1,1-\alpha}
\Bigr\},
\]
and \(\ell_{\mathrm{prof}}(\mathbf{t})=\sup_{\mathbf{r}\in\mathcal{R}} \ell(D\mid\mathbf{r},\mathbf{t})\) profiles rates inside locked bounds.

**Interpretation:** \(\mathcal{I}_v(D)\) is a **sharp certified interval** for \(\vartheta_v\) under \(\mathcal{C}\), not a classical Wald CI for a point-identified clock parameter.

### 2.5 Identification regimes (locked taxonomy)

| Status | Definition on \(\mathcal{I}_v\) | Code enum |
|--------|----------------------------------|-----------|
| **Unidentified** | \(\mathrm{width}(\mathcal{I}_v)\ge 0.99\,T\) (essentially \([0,T]\) up to sign for gaps) | `unidentified` |
| **Partially identified** | \(0 < \mathrm{width}(\mathcal{I}_v) < 0.99\,T\) | `partially_identified` |
| **Point identified** | \(\mathrm{width}(\mathcal{I}_v)\le \varepsilon_{\mathrm{pt}}\,T\) with \(\varepsilon_{\mathrm{pt}}=10^{-6}\) | `point_identified` |
| **Infeasible / abstain** | empty feasible set or pre-registered abstention rule | `infeasible` / `abstain` |

---

## 3. Locked constraints \(\mathcal{C}\) (default)

1. **Irreversible Poisson editing** on each site; independent sites given \((\mathbf{r},\mathbf{t})\).
2. **Global rate mode** (primary): one shared \(r\in[r_{\min},r_{\max}]\) profiled; site heterogeneity absorbed into adversarial/sensitivity axes, not into silent free rates that re-open scale.
3. **Dated anchor** \(T\): all leaves fixed at age \(T\); free parameters = ages of non-root internals in \([0,T]\).
4. **Rate–time guardrail:** nontrivial maps \((r,\mathbf{t})\mapsto(c r,\mathbf{t}/c)\) that exit \(\mathcal{R}\) or break \(T\) are **rejected as false absolute-time identification**.
5. **Topology:** either a fixed audited tree, or **union** of \(\mathcal{I}_v\) over a finite audited ensemble (leaf-swap class); robust order claims require sign agreement on all ensemble members.

Anything outside \(\mathcal{C}\) is sensitivity, not headline.

---

## 4. Pass / kill criteria (decisive scout + manuscript gate)

### 4.1 PASS (all required)

| Gate | Criterion |
|------|-----------|
| **P1 Narrowing** | \(\ge 1\) target with \(\mathrm{width}(\mathcal{I}_v)/T < 0.95\) and status \(\neq\) unidentified |
| **P2 Coverage** | Synthetic coverage of \(\mathcal{I}_v\) for true \(\vartheta_v\) within \(\approx\) nominal \(1-\alpha\) (scout tolerance: \(\ge (1-\alpha)-0.25\) at small \(n_{\mathrm{rep}}\)) |
| **P3 False ID** | False point-identification rate \(=0\) under rate–time rescaling adversarial checks |
| **P4 Invariance** | Likelihood invariance under \(r\mapsto c r\), \(t\mapsto t/c\) holds numerically; guardrail blocks nontrivial feasible rescale when \(T\) + bounds locked |
| **P5 Holdout hygiene** | Headline \(\mathcal{I}_v\) uses analysis sites only; rate bounds from calibration sites only |

**Scout reference (synthetic, \(T=1\)):** age:\(A\) partial interval materially inside \([0,1]\); coverage near nominal; false-ID \(=0\). Fill manuscript numbers from `{{RESULT:scout_ageA}}`, `{{RESULT:scout_coverage}}`.

### 4.2 KILL / REFRAME (any one)

| Kill | Meaning |
|------|---------|
| **K1** | All headline intervals \(\approx[0,T]\) after locking \(\mathcal{C}\) |
| **K2** | Topology ensemble union restores near-full width (topology dominates timing) |
| **K3** | Coverage fails to transfer across heterogeneity / observation axes |
| **K4** | Contribution collapses to “another clock MLE” (point estimates without certified sets) |
| **K5** | Rate–time invariance broken in code, or guardrail fails to block false absolute ID |

---

## 5. Role split (what B01 owns vs does not)

| Owns | Does **not** own |
|------|------------------|
| Mathematical identification of \(\mathcal{I}_v(D)\) | Biological narrative of developmental “when” beyond certified intervals |
| Exact small-tree profile certificates | Wet-lab recorder design / CRISPR chemistry optimization |
| Invariance + adversarial diagnostics | Claims of universal molecular clocks |
| Leakage-free site holdout protocol | Choosing which biological question is scientifically interesting |
| Topology-ensemble unions | Exhaustive search over all tree space |
| Software (`b01_clock`) and reproducible artifacts | Clinical or therapeutic timing claims |

**Role one-liner:** B01 is a **certificate engine for compatible timing sets**, not a chronogram product.

---

## 6. Forbidden silent defaults

- No hard-coded \([r_{\min},r_{\max}]\) in headline runs without provenance in config / certificate JSON.
- No reporting a point MLE as “the age” without \(\mathcal{I}_v\) and identification status.
- No mixing calibration and analysis sites for bounds and headline.
- No claiming absolute time from exposures alone (violates rate–time invariance).

---

## 7. Change control

Amendments to this lock require: (i) updated theorem statement in `02_THEORY/IDENTIFICATION_THEOREM.md`, (ii) matching tests in `tests/`, (iii) scout re-run, (iv) manuscript Methods sync. Soft preferences (grid size, \(\alpha\)) may move in `configs/` without reopening the estimand.
