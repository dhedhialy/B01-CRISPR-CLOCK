# Rate–Time Invariance

**Status:** Theory + computational guardrail lock  
**Package:** `b01_clock.solver.invariance`

---

## 1. Formal statement

### 1.1 Group action

Let \(c>0\). Define
\[
\Psi_c:\quad
(\mathbf{r},\mathbf{t})\;\longmapsto\;
\bigl(c\,\mathbf{r},\;\mathbf{t}/c\bigr).
\]
Componentwise, \(r_s\mapsto c r_s\), \(t_e\mapsto t_e/c\).

### 1.2 Invariance of the unconstrained likelihood

**Proposition (Exposure invariance).**  
Under the irreversible model \(p_{se}=1-\exp(-r_s t_e)\), every transition probability depends on the product \(\lambda_{se}=r_s t_e\) only. Hence
\[
\lambda_{se}\;\longmapsto\; (c r_s)\,(t_e/c)\;=\;\lambda_{se},
\]
and for all \(c>0\),
\[
\boxed{L(D\mid\mathbf{r},\mathbf{t})\;=\;L\bigl(D\mid\Psi_c(\mathbf{r},\mathbf{t})\bigr)}.
\]
Equivalently in log space: \(\ell(\mathbf{r},\mathbf{t})=\ell(c\mathbf{r},\mathbf{t}/c)\).

**Dual form (as in user-facing docs).** Writing times scaled up and rates scaled down,
\[
t\mapsto c\,t,\qquad r\mapsto r/c
\]
is the same orbit (replace \(c\) by \(1/c\) in \(\Psi_c\)).

### 1.3 What is invariant vs what is not

| Object | Invariant under \(\Psi_c\)? |
|--------|----------------------------|
| Exposures \(\lambda_{se}\) | Yes |
| Likelihood \(L(D\mid\cdot)\) | Yes (unconstrained) |
| Edit probabilities \(p_{se}\) | Yes |
| Absolute times \(t_e\), ages \(a_u\) | **No** |
| Absolute rates \(r_s\) | **No** |
| Ratios \(t_e/t_{e'}\), \(a_u/T\) when \(T\) also scales | Yes if \(T\mapsto T/c\) |
| Ratios \(a_u/T\) with **fixed** dated \(T\) | Ages change under \(\Psi_c\) \(\Rightarrow\) exit \(\mathcal{U}(T)\) |

---

## 2. What breaks the invariance (constrained model)

The scientific clock question is whether \(\Psi_c\) remains feasible under locked constraints \(\mathcal{C}\).

### 2.1 Dated anchor \(T\)

Ultrametric lock: \(\sum_{e\in\mathrm{path}(\rho\to\ell)} t_e = T\) for all leaves. Under \(\Psi_c\),
\[
\sum t_e/c \;=\; T/c.
\]
For \(c\neq 1\), \(\mathbf{t}/c\notin\mathcal{U}(T)\).  
**Break mechanism:** absolute leaf dating forbids free rescaling of all branch times.

### 2.2 Rate box \(\mathcal{R}=[r_{\min},r_{\max}]\)

Require \(r_s\in[r_{\min},r_{\max}]\). Then \(c r_s\) must remain inside the box for all sites (global mode: single \(r\)).  
**Break mechanism:** if \(c>1\) and \(r\) near \(r_{\max}\), or \(c<1\) and \(r\) near \(r_{\min}\), \(\Psi_c\) exits \(\mathcal{R}\).

### 2.3 Joint lock (primary B01 setting)

**Proposition (False absolute-ID blocker).**  
If \(\mathbf{t}\in\mathcal{U}(T)\) and \(r\in[r_{\min},r_{\max}]\), then for all \(c\neq 1\),
\[
\Psi_c(r,\mathbf{t})\;\notin\;\mathcal{U}(T)\times[r_{\min},r_{\max}]
\]
whenever either the dated constraint is enforced strictly or the scaled rate exits the box. In the package guardrail, nontrivial scales are marked
`feasible_with_anchor_and_bounds: false`.

**Interpretation:** A procedure that returns a single absolute chronology from exposures alone, without using constraints that break \(\Psi_c\), is committing a **false identification**.

### 2.4 Additional breakers (sensitivity axes)

| Constraint | Effect on \(\Psi_c\) |
|------------|----------------------|
| Site-specific known rates \(r_s=r_s^\star\) | Only \(c=1\) survives |
| Multiple dated internals | Overconstrains scaled ages |
| Branch-specific rate multipliers with fixed biological scale | Can break or reintroduce pseudo-invariance if multipliers are free |
| Soft priors on \(r\) (Bayesian clocks) | Softly break; posterior may still be prior-dominated — not B01’s estimand |

---

## 3. Computational checks (locked)

### 3.1 Numerical invariance (`check_rate_time_invariance`)

Given fitted/true \((\mathbf{r},\mathbf{t})\) and scale \(c\) (default \(2\)):

1. Compute \(\ell_0=\ell(\mathbf{r},\mathbf{t})\), \(\ell_1=\ell(c\mathbf{r},\mathbf{t}/c)\).  
2. Require \(|\ell_0-\ell_1|\le \texttt{atol}\,(1+|\ell_0|)\).  
3. Failure \(\Rightarrow\) implementation bug (kill gate P4).

### 3.2 Feasibility guardrail (`rate_time_transform_feasible`)

1. Form \(r'=c r\), \(\mathbf{t}'=\mathbf{t}/c\).  
2. Test \(r'\in\mathcal{R}\).  
3. Test whether dated depth / \(\mathcal{U}(T)\) is preserved (it is not for \(c\neq 1\)).  
4. Set `feasible_with_anchor_and_bounds` true only if the transform remains inside locked \(\mathcal{C}\) (only \(c=1\) in the primary lock).

**Expected scout outcome:** invariance `ok: true`, guardrail `feasible_with_anchor_and_bounds: false` for \(c=2\).

---

## 4. Implications for estimators

1. **MLE without \(\mathcal{C}\):** at best identifies exposures / shapes up to \(\Psi_c\).  
2. **MLE with \(T\) + rate box:** selects a profile optimum inside a partially identified set; the optimum is **not** “the true age,” only a point in \(\Theta(D)\).  
3. **Certificates:** report \(\mathcal{I}_v\), not only \(\hat\vartheta_v\).  
4. **Adversarial sims:** rescale true rates/times along \(\Psi_c\) and verify certificates do not claim Regime I (`false_identification_rate`).

---

## 5. Minimal counterexample (verbal)

Take one branch, one site, observe a single \(0\to 1\) with probability \(1-e^{-rt}\). Any pair \((r,t)\) with the same product fits equally. Dating \(t=T\) identifies \(r=-\log(1-p)/T\); bounding \(r\) without dating leaves a \(t\)-interval \([\lambda/r_{\max},\lambda/r_{\min}]\). Both constraints together sandwich time; neither alone yields a clock without remainder.
