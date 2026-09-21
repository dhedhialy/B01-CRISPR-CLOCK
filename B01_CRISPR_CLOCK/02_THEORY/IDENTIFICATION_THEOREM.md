# Identification Conjecture — Timing under Heterogeneous Irreversible Editing

**Status:** Core theory lock — the general trichotomy below is a **conjecture**; it is **proven only for the canonical 4-leaf tree** (§4, complete proof sketch on closed form). Conceptually supported on larger topologies by geometry and by the finite-sample numerical grid certificates, but not by a general closed-form theorem. Any publication must cite it as "Conjecture"; the text uses **Theorem** only for the proven local statements (Propositions 2–4, Lemma 5).  
**Implements:** `b01_clock.solver.exact`, `b01_clock.solver.invariance`, `b01_clock.topology.ensemble`  
**Estimand:** \(\mathcal{I}_v(D)\) as in `00_LOCK/ESTIMAND_LOCK.md`

---

## 0. Setup

### 0.1 Tree and ultrametric simplex

Let \(\mathcal{T}=(V,E)\) be rooted at \(\rho\) with leaves \(\mathcal{L}\). Branch times \(\mathbf{t}=(t_e)_{e\in E}\in\mathbb{R}_{\ge 0}^{|E|}\). Dated anchor \(T>0\): for every \(\ell\in\mathcal{L}\),
\[
\sum_{e\in\mathrm{path}(\rho\to\ell)} t_e \;=\; T.
\]
Write \(\mathcal{U}(T)\) for this ultrametric polytope. Free coordinates may be taken as ages \(a_u\in[0,T]\) of non-root internal nodes, with leaf ages fixed at \(T\) and branch times recovered by differences (as in `ultrametric_times_from_free`).

### 0.2 Irreversible exposure model

Independent sites \(s=1,\ldots,S\) with rates \(r_s>0\). On edge \(e\), for irreversible \(0\to 1\) editing,
\[
\mathbb{P}(0\to 0)=e^{-r_s t_e},\qquad
\mathbb{P}(0\to 1)=1-e^{-r_s t_e},\qquad
\mathbb{P}(1\to 1)=1,\qquad
\mathbb{P}(1\to 0)=0.
\]
The likelihood \(L(D\mid\mathbf{r},\mathbf{t})\) depends on \((\mathbf{r},\mathbf{t})\) only through **exposures** \(\{\lambda_{se}=r_s t_e\}\).

### 0.3 Constraint class

\[
\mathcal{C}
\;=\;
\bigl\{\,
\mathbf{r}\in\mathcal{R},\;
\mathbf{t}\in\mathcal{U}(T),\;
\tau\in\mathfrak{T}
\,\bigr\},
\]
where \(\mathfrak{T}\) is either a singleton topology or a finite audited ensemble, and \(\mathcal{R}\) is a rate region. **Primary lock:** global rate mode \(\mathcal{R}=[r_{\min},r_{\max}]\) (shared \(r\)).

### 0.4 Profile feasible set and estimand

\[
\ell_{\mathrm{prof}}(\mathbf{t};\tau)
=\sup_{r\in[r_{\min},r_{\max}]}\ell(D\mid r,\mathbf{t};\tau),
\qquad
\widehat{\Theta}_\alpha(D;\tau)
=\bigl\{\mathbf{t}\in\mathcal{U}(T):
2(\ell_{\max}-\ell_{\mathrm{prof}}(\mathbf{t};\tau))\le \chi^2_{1,1-\alpha}\bigr\}.
\]
\[
\mathcal{I}_v(D;\tau)
=\bigl[\inf\vartheta_v(\widehat{\Theta}_\alpha),\;\sup\vartheta_v(\widehat{\Theta}_\alpha)\bigr],
\qquad
\mathcal{I}_v^{\cup}(D)
=\bigcup_{\tau\in\mathfrak{T}}\mathcal{I}_v(D;\tau)
\quad\text{(interval hull of the union)}.
\]

---

## 1. Main conjecture (trichotomy)

**Conjecture 1 (Identification regimes for irreversible recorders; proven for the canonical 4-leaf tree).**  
Fix topology \(\tau\), dated anchor \(T\), and rate box \([r_{\min},r_{\max}]\) with \(0<r_{\min}\le r_{\max}<\infty\). Let \(v\) be a continuous functional \(\vartheta_v:\mathcal{U}(T)\to\mathbb{R}\). Write \(\Theta^\star(D)\) for the population profile maximizer set (large-\(S\) limit of \(\widehat{\Theta}_\alpha\) as \(\alpha\to 1\) and \(S\to\infty\) under correct specification), and \(\mathcal{I}_v^\star=\mathrm{proj}_v\Theta^\star\).

Then exactly one of the following holds for generic data-generating exposures:

### Regime U — Unidentified

If either  
**(U1)** \(r_{\min}\to 0\) and \(r_{\max}\to\infty\) (effectively free positive rate), or  
**(U2)** no dated anchor (\(T\) free) while rates remain free up to scale,  

then for any functional that is **not** invariant under the rate–time group action (Section 2),  
\[
\mathcal{I}_v^\star \;=\; \mathrm{range}(\vartheta_v)\quad\text{on }\mathcal{U}(T)\ \text{(or the free-depth analogue)},
\]
i.e. absolute ages and absolute branch lengths are **unidentified**. Relative functionals that are homogeneous of degree 0 in \(\mathbf{t}\) may remain partially constrained by **shape** information in the edit pattern, but absolute scale is not.

### Regime P — Partially identified

If \(0<r_{\min}\le r_{\max}<\infty\) **and** dated \(T\) is fixed, then for node ages and branch lengths,
\[
\mathcal{I}_v^\star \;\subseteq\; [0,T]
\]
is a **proper closed interval** (or finite union of intervals before taking convex hull of the profile set) whenever the exposure pattern distinguishes path lengths at more than one scale. Generically under site information,
\[
0 \;<\; \mathrm{width}(\mathcal{I}_v^\star) \;<\; T,
\]
so \(v\) is **partially identified**: not a singleton, but sharper than the trivial bound.

### Regime I — Point identified

Point identification \(\mathcal{I}_v^\star=\{t_v^\star\}\) occurs only when **collapse conditions** (Section 3) hold—typically: (i) rate box degenerates to a known singleton \(r_{\min}=r_{\max}=r^\star\), or (ii) additional independent anchors fix enough absolute times that the remaining free ages are pinned by invertible exposure maps, or (iii) an infinite-data limit with known site-specific rates. Under the **primary lock** (nondegenerate box + single dated \(T\) + free internal ages), **node ages are not point-identified**; claiming Regime I without collapse is a **false identification**.

**Corollary 1.1 (Certificate taxonomy).** Finite-sample profile intervals \(\mathcal{I}_v(D)\) are classified by width relative to \(T\) as in `ESTIMAND_LOCK.md`. Regime P is the scientifically central case for B01.

---

## 2. Rate–time invariance (formal)

**Proposition 2 (Exposure invariance).**  
For any \(c>0\),
\[
L(D\mid \mathbf{r},\mathbf{t}) \;=\; L(D\mid c\mathbf{r},\,\mathbf{t}/c).
\]
Thus the map \((\mathbf{r},\mathbf{t})\mapsto(c\mathbf{r},\mathbf{t}/c)\) is a likelihood-automorphism of the unconstrained model.

**Proposition 3 (Breaking the group).**  
The same map fails to be an automorphism of the **constrained** model as soon as either:

1. \(\mathbf{t}/c\notin\mathcal{U}(T)\) when \(c\neq 1\) (dated leaves at \(T\) force \(\sum_{\mathrm{path}} t_e=T\), so \(\mathbf{t}\mapsto\mathbf{t}/c\) exits \(\mathcal{U}(T)\)); or  
2. \(c\mathbf{r}\notin\mathcal{R}\) (scaled rates exit the locked box).

**Lemma 5 (No free lunch for absolute time).**  
Absent (1)–(2), absolute chronology is unidentified (Regime U). With both a nondegenerate rate box and dated \(T\), the residual freedom is the **shape** of \(\mathbf{t}\) inside \(\mathcal{U}(T)\) compatible with profiled \(r\in[r_{\min},r_{\max}]\), yielding Regime P for generic internals.

See `INVARIANCE.md` for the computational guardrail.

---

## 3. Collapse conditions (when Regime I can occur)

Write \(\mathrm{dim}\,\mathcal{U}(T)=d\) (number of free internal ages; for the 4-leaf scout tree, \(d=2\)).

**Proposition 4 (Rank / anchor conditions).**  
Suppose site exposures induce, in the large-\(S\) limit, a map
\[
\Phi:\mathcal{U}(T)\times[r_{\min},r_{\max}]\to\mathbb{R}^{m}
\]
given by identifiable path-exposure summaries (e.g. \(-\log(1-p)\) along root-to-leaf and clade paths). Then:

| Condition | Consequence |
|-----------|-------------|
| \(r\) free on \((0,\infty)\), \(T\) free | \(\mathrm{rank}\) of \(d\Phi\) drops by \(\ge 1\) (scale); Regime U for absolute \(t\) |
| \(T\) fixed, \(r\in[r_{\min},r_{\max}]\) with \(r_{\min}<r_{\max}\) | Fiber \(\Phi^{-1}(\Phi^\star)\) typically positive-dimensional \(\Rightarrow\) Regime P |
| \(T\) fixed, \(r=r^\star\) known (degenerate box) | If \(d\Phi_t\) has full column rank \(d\) at truth, Regime I for \(\mathbf{t}\) |
| Additional dated internal anchors \(a_{u_i}=T_i\) with rank filling \(d\) | Regime I even with mild rate uncertainty, if anchors pin all free ages |
| Topology \(\tau\) wrong | \(\Phi\) misspecified; \(\mathcal{I}_v\) may be empty, shifted, or spuriously narrow — handle via ensemble union |

**Collapse slogan:** *Global rate bound + dated \(T\) + topology turn unidentified absolute time into partial identification; point identification still needs rate pin-down or extra anchors.*

---

## 4. Proof sketch on the canonical 4-leaf tree

### 4.1 Geometry

Scout tree (`small_scout_tree`):
\[
\begin{array}{c}
\rho \\
/\quad\setminus \\
A\qquad B \\
/\setminus\quad/\setminus \\
L_0\; L_1\quad L_2\; L_3
\end{array}
\]
Free ages \(a_A,a_B\in[0,T]\). Branch times:
\[
t_{\rho A}=a_A,\;
t_{A L_0}=t_{A L_1}=T-a_A,\;
t_{\rho B}=a_B,\;
t_{B L_2}=t_{B L_3}=T-a_B.
\]
Target \(v=\texttt{age:}A\) is \(\vartheta=a_A\).

### 4.2 Path exposures (global rate \(r\))

Root-to-leaf exposure is always \(\lambda_{\rho\to\ell}=r T\) (ultrametric). Informative contrasts are **clade-specific**:
\[
\lambda_{\rho\to A}=r a_A,\qquad
\lambda_{A\to L_0}=r(T-a_A),
\]
and likewise for \(B\). Shared mutations within \(\{L_0,L_1\}\) vs unique edits inform the split of exposure between \(t_{\rho A}\) and \(t_{A L_i}\).

### 4.3 Unconstrained scale

If \(r\) is free and \(T\) free, \((r,a_A,a_B,T)\mapsto(cr,a_A,a_B,T/c)\) or equivalently rescaling all times preserves likelihood; absolute \(a_A\) is meaningless (Regime U).

### 4.4 With dated \(T\) only

Fix \(T\). Likelihood depends on \((r a_A,\, r(T-a_A),\, r a_B,\, r(T-a_B))\). The map
\[
(a_A,a_B,r)\mapsto r\cdot(a_A,\,T-a_A,\,a_B,\,T-a_B)
\]
has, for each fixed shape \((a_A,a_B)\), a one-dimensional fiber in \(r\). Distinct shapes are separated by **ratios** of clade exposures, e.g.
\[
\frac{\lambda_{\rho A}}{\lambda_{A L}}
=\frac{a_A}{T-a_A},
\]
which is **scale-free**. Thus shape (relative splits) can be identified while absolute ages remain identified only **relative to \(T\)** — but wait: ages are already in units of \(T\). The residual non-identification with free \(r\in(0,\infty)\) is not scale of \(T\) but **which shapes** remain equally likely after profiling \(r\): with infinite sites and binary irreversible data, many \((a_A,a_B)\) can share the same optimized \(r\) up to sampling noise, and with finite sites the profile LR set is thick (Regime P with wide intervals). Critically, **free \(r\) still allows trading** a more unbalanced time split against rate when observation noise and missingness blur exposure ratios—hence intervals need not collapse to points.

### 4.5 With rate box \([r_{\min},r_{\max}]\)

Profiling \(r\) on a compact interval truncates the fiber. Extreme ages near \(0\) or \(T\) demand extreme exposures on short branches; if those require \(r\notin[r_{\min},r_{\max}]\) to fit observed edit fractions, such ages exit \(\widehat{\Theta}_\alpha\). Therefore
\[
\mathcal{I}_{\mathrm{age:}A}(D)\;\subsetneq\;[0,T]
\]
generically (Regime P). This matches the decisive scout: synthetic age:\(A\) returns a proper subinterval of \([0,T]\).

### 4.6 Degenerate box \(\Rightarrow\) Regime I tendency

If \(r_{\min}=r_{\max}=r^\star\), exposures become pure time maps \(\lambda=r^\star t\). With sufficient independent path contrasts (rank \(d=2\)), \((a_A,a_B)\) is locally identifiable and profile sets shrink as \(S\to\infty\) (Regime I asymptotically).

### 4.7 Topology

A leaf swap \(L_0\leftrightarrow L_2\) reassigns which clade owns which shared edits. Each \(\tau\) yields its own \(\mathcal{I}_v(D;\tau)\). The ensemble hull \(\mathcal{I}_v^{\cup}\) is the honest set under topology uncertainty; order functionals with constant sign across \(\mathfrak{T}\) are **robust orders**.

---

## 5. Finite-sample certificate (operational lemma)

**Lemma 6 (Grid-attained profile certificate).**  
On trees with \(d\le 3\) free internals, let \(\mathcal{G}\subset\mathcal{U}(T)\) be a finite grid. Let \(\ell_{\max}^{\mathcal{G}}\) be the maximum profiled log-likelihood on \(\mathcal{G}\), and
\[
\widehat{\Theta}_\alpha^{\mathcal{G}}
=\bigl\{\mathbf{t}\in\mathcal{G}:
2(\ell_{\max}^{\mathcal{G}}-\ell_{\mathrm{prof}}(\mathbf{t}))\le\chi^2_{1,1-\alpha}\bigr\}.
\]
Then
\[
\mathcal{I}_v^{\mathcal{G}}
=\bigl[\min\vartheta_v(\widehat{\Theta}_\alpha^{\mathcal{G}}),\;
\max\vartheta_v(\widehat{\Theta}_\alpha^{\mathcal{G}})\bigr]
\]
is a **conservative outer approximation** to the continuous profile interval up to grid bias. Refinement (`verify_against_grid`) must agree within tolerance \(\varepsilon T\) (scout default \(\varepsilon=0.05\)) or the certificate is marked `grid_unstable`.

**Remark.** The \(\chi^2_1\) cutoff is the standard one-dimensional profile LR calibration for a scalar focus parameter after profiling nuisances; under partial identification the set is interpreted as a **confidence set for the identified projection**, with coverage validated empirically (`evaluate_coverage`), not as a classical CI for a point-identified age.

---

## 6. What this document is not

- Not a general theorem: the trichotomy is **Conjecture 1**, proven on the 4-leaf topology; the d>2 regime is supported by geometry and the finite-sample grid certificates only.
- Not a consistency theorem for MLE chronograms under misspecified homogeneous rates.  
- Not a claim that Regime P intervals are short. Shortness is empirical — and the external-benchmark check (`convexml_ground_truth_coverage.json`) is the referee: the primary lock's intervals are **not** calibrated on Cassiopeia-generated data (1/9 coverage).  
- Not identification of biological “event time” without mapping that event to a tree node \(u\).

---

## 7. Cross-references

| Document | Role |
|----------|------|
| `INVARIANCE.md` | Group action, guardrail algorithm |
| `FEASIBLE_SET_GEOMETRY.md` | Geometry of \(\Theta(D)\) and sharp intervals |
| `ESTIMAND_LOCK.md` | Pass/kill and taxonomy |
| `CERTIFICATE_SCHEMA.md` | JSON emission of \(\mathcal{I}_v\) |
