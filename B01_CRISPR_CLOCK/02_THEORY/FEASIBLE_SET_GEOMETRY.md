# Feasible-Set Geometry — \(\Theta(D)\) and Sharp Intervals

**Status:** Theory lock  
**Implements:** projection step in `exact_profile_bounds`; union in `ensemble_timing_union`

---

## 1. Ambient space

### 1.1 Ultrametric polytope

With dated leaf age \(T\) and free internal ages \(\mathbf{a}=(a_u)_{u\in V_{\mathrm{int}}\setminus\{\rho\}}\),
\[
\mathcal{U}(T)
=\bigl\{\mathbf{a}: 0\le a_{\mathrm{parent}(u)}\le a_u\le T\ \text{for all }u;\ 
a_\ell=T\ \forall\ell\in\mathcal{L}\bigr\}
\]
(expressed in age coordinates). Equivalently, branch times are nonnegative differences along edges. \(\mathcal{U}(T)\) is a compact convex polytope of dimension \(d=\#\{u:u\text{ free}\}\).

**Scout:** \(d=2\), coordinates \((a_A,a_B)\in[0,T]^2\) (with parent age \(0\)).

### 1.2 Extended parameter space

\[
\Xi \;=\; \mathcal{U}(T)\times[r_{\min},r_{\max}].
\]
Likelihood level sets live in \(\Xi\); after profiling \(r\), geometry collapses to subsets of \(\mathcal{U}(T)\).

---

## 2. Population feasible set

Under correct specification with true \(( \mathbf{t}^\star, r^\star)\), the **population identified set** for times is
\[
\Theta^\star(D_\infty)
=\bigl\{\mathbf{t}\in\mathcal{U}(T):
\exists r\in[r_{\min},r_{\max}]\ \text{s.t.}\ 
\mathbb{P}_{r,\mathbf{t}}=\mathbb{P}_{r^\star,\mathbf{t}^\star}\bigr\}.
\]
Because probabilities depend on exposures, this is the set of ultrametric timings whose exposure vectors (as \(r\) varies in the box) can match the true exposure vector.

**Geometry facts:**

1. \(\Theta^\star\) is compact (closed subset of compact \(\mathcal{U}(T)\)).  
2. \(\Theta^\star\) need not be convex: profile LR sublevel sets in finite samples can be nonconvex in age coordinates when clade exposures interact.  
3. \(\mathbf{t}^\star\in\Theta^\star\) whenever \(r^\star\in[r_{\min},r_{\max}]\).  
4. If \(r_{\min}<r_{\max}\), \(\Theta^\star\) typically has positive \(d\)-dimensional measure or positive \((d-1)\)-dimensional measure (partial ID).

---

## 3. Finite-sample certificate set

\[
\widehat{\Theta}_\alpha(D)
=\bigl\{\mathbf{t}\in\mathcal{U}(T):
\ell_{\mathrm{prof}}(\mathbf{t})\;\ge\;\ell_{\max}-\tfrac12\chi^2_{1,1-\alpha}\bigr\}.
\]

**Properties:**

| Property | Statement |
|----------|-----------|
| Nestedness | \(\alpha<\alpha'\Rightarrow\widehat{\Theta}_{\alpha}\subseteq\widehat{\Theta}_{\alpha'}\) (larger sets for higher confidence) |
| Nestedness in rate box | Widening \([r_{\min},r_{\max}]\) weakly expands \(\widehat{\Theta}_\alpha\) |
| Nestedness in \(S\) | More analysis sites typically shrink \(\widehat{\Theta}_\alpha\) (information), all else equal |
| Grid outer approx | \(\widehat{\Theta}_\alpha^{\mathcal{G}}\) on a grid yields interval hulls that are staircased outer bounds |

---

## 4. Sharp timing intervals (projections)

For functional \(\vartheta_v:\mathcal{U}(T)\to\mathbb{R}\),
\[
\mathcal{I}_v(D)
=\bigl[\,\underline{\vartheta}_v,\,\overline{\vartheta}_v\,\bigr],
\quad
\underline{\vartheta}_v=\inf_{\mathbf{t}\in\widehat{\Theta}_\alpha}\vartheta_v(\mathbf{t}),
\quad
\overline{\vartheta}_v=\sup_{\mathbf{t}\in\widehat{\Theta}_\alpha}\vartheta_v(\mathbf{t}).
\]

**Sharpness.** Among all intervals guaranteed to cover \(\{\vartheta_v(\mathbf{t}):\mathbf{t}\in\widehat{\Theta}_\alpha\}\), \(\mathcal{I}_v\) is the shortest. It is **not** sharp for the true \(\vartheta_v(\mathbf{t}^\star)\) beyond what \(\widehat{\Theta}_\alpha\) covers; coverage is a frequentist property of the set estimator, validated by simulation.

### 4.1 Linear functionals (ages, branches, gaps)

If \(\vartheta_v\) is linear in age coordinates (all locked targets are), then \(\mathcal{I}_v\) is the range of a linear form over \(\widehat{\Theta}_\alpha\). If \(\widehat{\Theta}_\alpha\) were convex, extrema would lie on the boundary; with grid approximation we enumerate feasible grid points.

### 4.2 Width as scientific content

\[
w_v \;=\; \overline{\vartheta}_v-\underline{\vartheta}_v,
\qquad
w_v/T\in[0,1]\ \text{(ages/branches)}.
\]

| Width | Geometry | Status |
|-------|----------|--------|
| \(w_v/T\approx 1\) | Projection fills the ambient range | Unidentified |
| \(0<w_v/T<1\) | Proper slab / band in \(\mathcal{U}(T)\) | Partially identified |
| \(w_v/T\approx 0\) | Fiber collapses in direction of \(\vartheta_v\) | Point identified |

---

## 5. Pictures in the 4-leaf age plane

Coordinates \((a_A,a_B)\in[0,T]^2\).

1. **Unconstrained profile ridges:** level sets of \(\ell_{\mathrm{prof}}\) form ridges along directions trading clade balance when \(r\) is free on a wide box.  
2. **Rate box:** cuts the ends of ridges where required \(r\) exits \([r_{\min},r_{\max}]\), leaving a compact blob \(\widehat{\Theta}_\alpha\).  
3. **Projection to age:\(A\):** horizontal span of the blob \(\Rightarrow\mathcal{I}_{\mathrm{age:}A}\).  
4. **Projection to order_gap:\(A,B\):** span of \(a_B-a_A\) over the blob; may still cover both signs (order not identified) even when each age is partially identified.  
5. **Topology ensemble:** each \(\tau\) yields a blob \(\widehat{\Theta}_\alpha(\tau)\); union of projections expands \(\mathcal{I}_v\). Robust order = all blobs lie in \(\{a_B>a_A\}\) or all in \(\{a_B<a_A\}\).

*(Manuscript figures: `{{FIG:feasible_blob}}`, `{{FIG:projection_ageA}}`, `{{FIG:ensemble_union}}`.)*

---

## 6. Operations that enlarge or shrink \(\Theta\)

| Operation | Effect on \(\widehat{\Theta}_\alpha\) / \(\mathcal{I}_v\) |
|-----------|------------------------------------------------------|
| Increase \(S_{\mathrm{an}}\) | Shrink (typically) |
| Widen rate box | Enlarge |
| Tighten rate box toward \(r^\star\) | Shrink; may approach Regime I |
| Raise confidence \(1-\alpha\) | Enlarge |
| Add topology ensemble members | Enlarge union |
| Dropout / silencing | Enlarge (information loss) |
| Branch-rate heterogeneity unmodeled | Possible miscoverage (adversarial axis) — width may shrink **spuriously** |

Spurious narrowing under misspecification is why adversarial coverage and false-ID rates are first-class, not optional.

---

## 7. Relationship to convex relaxations (ConvexML collision)

ConvexML-style methods optimize a point \(\hat{\mathbf{t}}\) (possibly unique under a strongly convex surrogate). Geometrically that selects a **single point** inside or near \(\Theta(D)\). B01’s object is the **set**. A point estimator can sit inside a wide \(\widehat{\Theta}_\alpha\) while the sharp interval remains large — the clock is still only partially identified.

---

## 8. Certificate geometry fields (for JSON)

Emit with each target:

- `lower`, `upper`, `width_frac = (upper-lower)/T`  
- `status` \(\in\) `{unidentified, partially_identified, point_identified, infeasible, abstain}`  
- `n_feasible_grid` (volume proxy on grid)  
- optional `ensemble_union_lower/upper`, `robust_order`

These fields are the operational geometry of \(\mathcal{I}_v\).
