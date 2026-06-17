# Doedel-Keller-Kernévez 1991 Part I + Part II — full-page deep-read verdict

**Read 2026-06-17 AET.** Combined verdict for task #375, the canonical methodology reference for #347 Phase 1 (Floquet bifurcation framework, just shipped at commit `e98e49b`). Source: every page of both PDFs in `/home/bruce/dev/cyclers_pdf/papers/`:

- `doedel-keller-kernevez-1991-numerical-analysis-bifurcation-problems-I-finite-dimensions-IJBC-1(3)-493-520.pdf` (28 pp.)
- `doedel-keller-kernevez-1991-numerical-analysis-bifurcation-problems-II-infinite-dimensions-IJBC-1(4)-745-772.pdf` (28 pp.)

Both papers read cover-to-cover including the references; 56 pages total.

## 1. Headers (verbatim from PDF p.493 and p.745)

### Part I
- Title (verbatim): **"NUMERICAL ANALYSIS AND CONTROL OF BIFURCATION PROBLEMS (I): BIFURCATION IN FINITE DIMENSIONS"**.
  - Note: the task brief abbreviated this to "Numerical Analysis of Bifurcation Problems"; the verbatim title includes **"AND CONTROL OF"** — this is the deliberate scope hook (every paper section is staged so the §3 control framework can attach). Not a transcription typo on the user's side; just a brief-vs-source delta.
- Authors (verbatim): **EUSEBIUS DOEDEL** (Department of Computer Science, Concordia University, 1455 Boulevard de Maisonneuve O., Montréal, Québec, H3G 1M8, Canada); **HERBERT B. KELLER** (Applied Mathematics 217-50, California Institute of Technology, Pasadena, California, 91125, USA); **JEAN PIERRE KERNEVEZ** (Mathématiques Appliquées, Université de Technologie de Compiègne, B.P. 233, 60206 Compiègne, France).
- Venue: **International Journal of Bifurcation and Chaos, Vol. 1, No. 3 (1991) 493-520**. World Scientific Publishing Company. Section: "Tutorials and Reviews" (banner top-right p.493).
- Received: **May 20, 1991**.
- Length: **28 pages** (pp. 493-520).
- DOI: not printed on the PDF; World Scientific assigns `10.1142/S0218127491000397` (inferred from venue/vol/issue/page; not source-derivable from PDF — flagged as inferred).

### Part II
- Title (verbatim): **"NUMERICAL ANALYSIS AND CONTROL OF BIFURCATION PROBLEMS (II): BIFURCATION IN INFINITE DIMENSIONS"**.
- Authors (verbatim, same three, same affiliations as Part I).
- Venue: **International Journal of Bifurcation and Chaos, Vol. 1, No. 4 (1991) 745-772**. World Scientific Publishing Company.
- Received: **May 20, 1991** (same as Part I — submitted as a single two-part work).
- Length: **28 pages** (pp. 745-772).
- DOI: not printed; inferred `10.1142/S0218127491000555` (do not cite without independent check).

Both papers cross-cite: Part I's abstract states "Part II of this paper deals with ordinary differential equations and will appear in the next issue"; Part II's abstract states "This is Part II of the paper that appeared in the preceding issue [Doedel et al., 1991]". Reading Part I first is mandatory — Part II's BVP framework is presented as a generalization of Part I's algebraic framework, with §1.4 (p.752 Eqs 1.9-1.13) being the explicit lifting of Part I's pseudoarclength equations to the BVP-discretized infinite-dimensional case.

## 2. What Part I is (two paragraphs)

A graduate-tutorial exposition of the **AUTO software's bifurcation engine for algebraic equations** `f(u, λ) = 0`, `u ∈ Rⁿ`, `λ ∈ R` (Eq. 1.1, p.493 — actually first stated as an ODE on p.494 then immediately specialized to stationary solutions giving Eq. 1.3). The material derives from a one-semester graduate course given by Doedel at U. Utah, U. Minnesota, and Concordia 1987-1990, and from Keller's 1987 Bangalore lecture notes (Sec. 2 derived from these — p.493 col. 2). The paper covers: (i) continuation past folds via Keller's pseudoarclength method (§2.1, Eqs 2.5 with the bordered LU decomposition + bordering algorithm, pp.502-504); (ii) detection of bifurcation points via determinant sign-change of the augmented Jacobian (§2.2, Eq. 2.9, p.506); (iii) **branch switching** via the Algebraic Bifurcation Equation `c₁₁α² + 2c₁₂αβ + c₂₂β² = 0` (Eq. 2.7, p.505) and the **orthogonal direction method** (Eq. 2.8, p.506) — which is the load-bearing recipe for #347 Phase 1's P1.4; (iv) fold continuation in two parameters via the extended system (§2.3, Eq. 2.10a-b, p.508-509); (v) Hopf bifurcation detection and continuation via complex eigenvalue-eigenvector tracking and the extended complex system (§2.4, Eqs 2.14-2.18, pp.510-512); and (vi) bifurcation analysis of discrete maps (§2.5, pp.514-515).

The paper's §3 is the **optimization-and-control** companion (pp.515-519): given an objective `g(u, λ)` and equality constraints `f(u, λ) = 0` plus optional inequality constraints `hᵢ ≤ 0`, derive the Karush-Kuhn-Tucker optimality system (Eqs 3.2-3.7) and the **successive continuation strategy** (§3.2, p.516-517) where free control parameters are successively unfrozen as zero-crossings of the Lagrange-multiplier-like quantities `τᵢ` are detected — extending the algebraic-systems framework into a constrained-optimization framework that reuses the same Newton-bordered-linear-solver kernel. §3.3 closes with the **projected gradient method** (p.519, Eqs 3.12-3.15) for cases where the equilibrium-manifold constraint must be enforced explicitly. Total: 28 pp., 19 figures, 1.7MB PDF, ~30 references including Keller 1977 (the pseudoarclength original) and Keller 1987 (the Bangalore notes which Sec. 2 derives from).

## 3. What Part II is (two paragraphs)

The **infinite-dimensional generalization** to ODE boundary-value problems `u'(t) = f(u(t), µ, λ)` (Eq. 1.1, p.745), `t ∈ [0,1]`, with **general nonlinear nonseparated boundary conditions** `b(u(0), u(1), µ, λ) = 0` (Eq. 1.2, with `b(·) ∈ R^{nb}`) and **linear-or-nonlinear integral constraints** `∫₀¹ q(u(s), µ, λ) ds = 0` (Eq. 1.3, with `q(·) ∈ R^{nq}`). The paper rigorously establishes (§1.2 pp.748-750) that for this class, the **orthogonal-collocation-with-piecewise-polynomials** discretization (Gauss collocation, m collocation points per subinterval, m between 2 and 7 in AUTO — p.752 col. 2) gives `O(h^m)` accuracy globally and `O(h^{2m})` superconvergence at the meshpoints (de Boor & Swartz 1973 result cited p.747). The **convergence and multiplicity theorem** (Beyn & Doedel 1981, restated p.750) guarantees that for sufficiently small mesh `h`, every isolated continuous solution is approximated by exactly one discrete solution; with sublinear `f` (Eq. 1.7), spurious extraneous discrete solutions are guaranteed to disappear as `h → 0`. §1.4 (pp.752-754) gives the **explicit AUTO BVP discretization** as a system of `mnN + n_b + n_q + 1` nonlinear equations in `mnN + n_λ + 1` unknowns (Eqs 1.9-1.13), and §1.4 closes with a **detailed bordered-block-bidiagonal Jacobian structure with condensation-of-parameters** (Figs 1.3-1.5, pp.753-754) for efficient linear solves — this is the algorithm that produces Floquet multipliers "as a by-product of the decomposition of the Jacobian" (the phrase the 2003 paper invokes).

§2 (pp.754-758) covers **periodic-solution computation**: re-scaling `t → t/T` so the period `T` becomes a free unknown, the **classical orthogonality phase condition** `(u(0) - u_{k-1}(0))* u_{k-1}'(0) = 0` (Eq. 2.6) and the more-efficient **integral phase condition** `∫₀¹ u(t)* u_{k-1}'(t) dt = 0` (Eq. 2.8) which fixes the time-translation freedom, the pseudoarclength constraint with arclength weights `θ_u, θ_T, θ_λ` (Eq. 2.9 plain, Eq. 2.13-style weighted), starting from a Hopf bifurcation point via `φ(t) = sin(2πt/T₀) w_s + cos(2πt/T₀) w_c` with `(w_s, w_c)` a null vector of the 2n×2n bordered matrix (p.756 top), the **Floquet multiplier extraction** as eigenvalues of `-A₁⁻¹ A₀` where `A₀, A₁` are the boxed blocks of the final condensed Jacobian (Eq. block-structure p.757-758), and the persistence theorem (§2.2 p.756). §3 (pp.758-763) handles **infinite-period heteroclinic / homoclinic orbits** via a projected boundary condition (Eqs 3.1-3.6) where `u(0) = w₀ + ε₀ Σ c_{0i} v_{0i}` projects onto the stable manifold's tangent space at the fixed point `w₀` and `u(1) = w₁ + ε₁ Σ c_{1i} v_{1i}` onto the unstable manifold at `w₁`, with an integral phase condition (Eq. 3.5) replacing periodicity. §4 (pp.764-770) lifts the §3-of-Part-I control framework into the infinite-dimensional case via the Lagrange-multiplier optimality system (Theorem 4.8, p.764-765, proven via Implicit Mapping Theorem on Banach spaces) with concrete adjoint-state equations (Eqs 4.20-4.24, p.769) for the BVP-with-integral-constraint case.

## 4. The pseudoarclength continuation foundation (Part I) — comparison to project code

### 4.1 The Doedel-Keller-Kernévez 1991 algorithm

**Pseudoarclength continuation** (Part I §2.1, p.502 col. 2 — "This is the most popular continuation method. It is due to Keller [1977]"). The setup: given a converged solution `(u₀, λ₀)` of `f(u, λ) = 0` and its unit tangent direction `(u̇₀, λ̇₀)` satisfying `||u̇||² + λ̇² = 1`, find the next solution `(u₁, λ₁)` at arclength `Δs` along the branch by solving the augmented system (p.502):

```
f(u₁, λ₁) = 0
(u₁ - u₀)* u̇₀ + (λ₁ - λ₀) λ̇₀ - Δs = 0
```

Geometrically: `(u₁, λ₁)` lies on the solution manifold AND in a hyperplane orthogonal to `(u̇₀, λ̇₀)` at distance `Δs`. Newton's method on this 2-equation 2-unknown extension gives the linear system (verbatim Eq. on p.502):

```
( (f_u^1)^(ν)   (f_λ^1)^(ν) )  ( Δu_1^(ν) )     ( f(u_1^(ν), λ_1^(ν))                      )
( u̇₀*           λ̇₀          )  ( Δλ_1^(ν) )  =  ( (u_1^(ν) - u_0)*u̇₀ + (λ_1^(ν) - λ_0)λ̇₀ - Δs )
```

The new tangent `(u̇₁, λ̇₁)` after convergence is recovered from `f_u u̇ + f_λ λ̇ = 0` plus normalization `u̇₀* u̇₁ + λ̇₀ λ̇₁ = 1` (one extra backsolve at the converged Jacobian — p.502 bottom). The **bordering algorithm** (p.503, the LU-decomposition recipe `Lγ = c`, `U*β = b`, `δ = d - β*γ`) makes the augmented system's solve as cheap as solving the unbordered Jacobian itself when the unbordered `f_u` is sparse. The **bordering lemma** (Keller 1977, restated p.502) guarantees the augmented matrix `(f_u f_λ ; u̇* λ̇)` is nonsingular whenever the original Jacobian rank is at least n-1 (i.e. at a simple fold OR a simple bifurcation point) and the tangent is not in the null space of the original Jacobian.

In AUTO the **scale-weighted form** is used (p.504 top): `⟨(u₁, λ₁), (u₂, λ₂)⟩ ≡ θ_u² u₁* u₂ + θ_λ² λ₁ λ₂` with default `θ_u = θ_λ = 1`; this matters for problems where `u` and `λ` are in different physical units. The **step-size control** is adaptive (p.504 bottom): "If Newton's method converges rapidly, then the step size is increased. If the Newton iteration converges slowly or if it fails to converge at all, then the step size is halved."

### 4.2 The project's existing pseudoarclength implementation

`src/cyclerfinder/search/cr3bp_jacobi_arclength.py` (#291 family tracer, "Pseudo-arclength continuation of CR3BP symmetric (k1,k2) cyclers in the Jacobi (C) parameter") — the docstring's wording matches the Doedel formulation exactly. `src/cyclerfinder/search/cr3bp_3d_family_tracer.py` (4-D ambient `z = (x₀, z₀, ẏ₀, T)`) explicitly names "pseudo-arclength continuation" as its default mode and lists "natural-T" and "natural-x₀" as the natural-parameter fallback modes — which is exactly the Doedel §2.1 architectural pattern (natural-parameter continuation fails at folds, pseudoarclength is the universal replacement).

**The project's continuation is Keller-1977 pseudoarclength, i.e. exactly the method Doedel et al. 1991 Part I §2.1 standardize.** It is not a different algorithm. Specifically:

- **Same augmented-system formulation** (cr3bp_3d_family_tracer.py line ~487 `_pseudo_arclength_step`, cr3bp_jacobi_arclength.py around the secant-tangent predictor + arclength-constraint corrector).
- **Same adaptive step control** (the project uses success-based step doubling and halving on Newton-fail, matching p.504).
- **Same bordering-style linear solve** when the unbordered Jacobian is the dominant cost — though the project's `n_dim` is small enough (4 for 3D-family-tracer, 3 for jacobi-arclength) that an explicit bordered-LU is overkill; the project uses dense `np.linalg.solve` on the full augmented matrix, which is the correct simplification for small `n_dim`.

**The one architectural feature the project does NOT mirror**: AUTO's tangent-direction propagation strategy (p.504 mid: `u̇_{j-1} ≈ (1/Δs)(u_{j-1} - u_{j-2})` secant approximation, then rescaled to satisfy the unit-norm constraint). The project recomputes the tangent freshly at each step from `f_u u̇ + f_λ λ̇ = 0` (the "exact" Eq. on p.502 bottom). For small `n_dim` this is fine — the secant approximation is an efficiency trick that only matters when the per-step Newton solve dominates. Not an issue.

**Verdict**: project pseudoarclength is canonically aligned with Doedel-Keller-Kernévez 1991 §2.1. No change recommended.

## 5. The branch-switching kernel-computation recipe (Part II) — load-bearing extraction

This is the section the digest exists to extract. The brief expects Part II to contain a deep-detail BVP-branch-switching algorithm that the existing project Phase 1 P1.4 module may or may not match. **The actual finding**: Part II §2.3 last paragraph (p.758 right column, after the Floquet accuracy table) says verbatim:

> *"AUTO can switch branches at transcritical, pitchfork, and period doubling bifurcations. The algorithms are generalizations of those described in Part I with minor modifications for the case of period-doubling bifurcations."*

In other words: **the branch-switching kernel-computation recipe lives in Part I §2.2 (pp.504-506); Part II inherits it without modification for transcritical and pitchfork, and only modifies it for period-doubling**. This is exactly the "Doedel 1991b" deep-detail reference the 2003 paper invokes — but the load-bearing equations are 1991a's. The recipe:

### 5.1 Recipe step 1 — kernel of the rank-deficient Jacobian (Part I p.505, Eqs around 2.6-2.7)

At a **simple singular point** `x₀ = x(s₀)` on a solution branch of the augmented system `F(x; s) = (f(x); (x - x₀)* ẋ₀ - s)`, `F_x` has rank `n - 1`. Case (ii) of the definition (p.505 col. 1): `dim N(f_x⁰) = 2` and `f_λ⁰ ∉ R(f_x⁰)`. The null space `N(f_x⁰) = Span{φ₁, φ₂}`; the cokernel `N((f_x⁰)*) = Span{ψ}`. The TWO null vectors `φ₁, φ₂` are the load-bearing objects.

**Practical computation** (p.505 right col., "Practical computation of the bifurcation direction"): pick `φ₁ ≡ ẋ₀` (the tangent to the GIVEN branch — already in hand from the continuation at the previous step, since we just computed it). Then pick `φ₂ ⊥ φ₁`; mathematically, `φ₂` is a second null vector of the augmented matrix `F_x⁰ = (f_x⁰ ; ẋ₀*)` — and Part I p.505 col. 2 establishes that the augmented matrix's null space is one-dimensional even though `f_x⁰`'s is two-dimensional, because the row `ẋ₀*` picks out one of the two null directions. Computationally, `φ₂` is found by **the same null-vector backsolve recipe as for any singular bordered matrix** (Part I p.503 right col., "Computation of the right null vector"): if `A` (here `F_x⁰`) has been Gauss-decomposed as `A = PL̂ÛQ` with `L̂` lower-triangular (1's on diagonal, no zeros) and `Û` upper-triangular with **a zero in the last diagonal position** (the singularity), then the null vector is recovered from one backsolve `Uv = u` (with `u` the last column of `Û` above the zero) plus the permutation `Q`. **One additional backsolve once `A` has been decomposed for any other reason — like the corrector's last Newton step.**

### 5.2 Recipe step 2 — the second null vector → branch direction (Part I Eq. 2.7, p.505)

Once `φ₁` and `φ₂` are in hand, the **Algebraic Bifurcation Equation (ABE)** (Eq. 2.7, p.505, verbatim):

```
c₁₁ α² + 2 c₁₂ α β + c₂₂ β² = 0,    c_{ij} ≡ ψ* f_xx⁰ φᵢ φⱼ,  i,j = 1, 2.
```

is a quadratic in `(α, β)` whose two roots are the two branch directions. Because `φ₁` is chosen as `ẋ₀` (the given branch tangent), one root is `(α₁, β₁) = (1, 0)` — recovering the given branch — and by direct substitution `c₁₁ = 0`. The OTHER root gives the **transverse branch direction**:

```
α₂ / β₂ = -c₂₂ / (2 c₁₂)
```

(p.505 right col., bottom). So:
```
x'₀ ≡ α₂ φ₁ + β₂ φ₂   (the transverse branch direction, before normalization to ||x'₀|| = 1).
```

A bifurcation exists iff the discriminant `Δ₀ ≡ c₁₂² - c₁₁ c₂₂` is positive. Because `c₁₁ = 0` (by the choice `φ₁ = ẋ₀`), `Δ₀ = c₁₂² ≥ 0`, so the discriminant is positive iff `c₁₂ ≠ 0`. The condition `c₁₂ ≠ 0` is the **non-degeneracy condition for the bifurcation**, and Part I p.506 col. 2 establishes that this is also exactly the condition `κ₀ ≠ 0` (the rate at which the bifurcation-detecting eigenvalue crosses zero):

```
κ̇₀ = c₁₂ / (Ψ* φ₂)
```

— so `c₁₂ ≠ 0` ⇔ `κ̇₀ ≠ 0` ⇔ the eigenvalue crosses zero transversally ⇔ `det F_x` changes sign across `s₀` (the determinant-sign-change detector triggered).

### 5.3 Recipe step 3 — the ε scaling for the initial perturbation

This is the **orthogonal direction method** (Part I p.506 left col., final paragraph above "Detection of bifurcation points"). The transverse branch direction `x'₀` requires evaluating the second derivative `f_xx⁰` (Eq. 2.7's `c_{ij}`). **This is avoided** by using the simpler **orthogonal direction method**:

> *"Recall that we have chosen φ₂ ⊥ φ₁ with φ₁ = ẋ₀. The branch switching procedure then simply consists of computing the null vector φ₂, followed by the first pseudoarclength continuation step* `f(x₁) = 0, (x₁ - x₀)* φ₂ - Δs = 0` *with initial approximation* `x₁⁽⁰⁾ = x₀ + Δs φ₂`. *This method need not always be successful. But it works well in most practical applications."* — Part I p.506 verbatim.

**This is the entire ε-scaling story.** The "ε" is `Δs`. The initial approximation for the corrector on the transverse branch is `x₀ + Δs φ₂` where `φ₂` is the *orthogonalized second null vector*, NOT the (α₂, β₂)-combined direction. The orthogonality `φ₂ ⊥ φ₁ = ẋ₀` is the load-bearing geometric property: it guarantees the perturbation is into the transverse-branch tangent plane, NOT along the given branch. The pseudoarclength constraint `(x₁ - x₀)* φ₂ - Δs = 0` then ENFORCES that the corrector lands at a member at arclength `Δs` from the bifurcation point, along the φ₂ direction.

**Scaling of Δs**: Part I does not give a specific magnitude — the value is problem-dependent. AUTO's default step-size control (p.504 bottom) doubles `Δs` on rapid convergence and halves on slow convergence, with user-specifiable floor and ceiling. **The orthogonality is the load-bearing structural property; the magnitude is a practical tuning parameter.**

### 5.4 Recipe step 4 — corrector setup on the transverse branch

For the algebraic case (Part I), the corrector on the transverse branch is **the same pseudoarclength continuation step** as for any other family member (Eq. 2.8, p.506 — reproduced from §2.1):

```
f(x₁) = 0
(x₁ - x₀)* φ₂ - Δs = 0
```

solved via Newton on the augmented 2-equation 2-unknown system; same bordering-algorithm linear solve as for the given branch.

For the BVP case (Part II), the corrector is **the same BVP discretization** (orthogonal collocation, Eqs 1.9-1.13 p.752) with the only modification being that the pseudoarclength constraint Eq. 1.13 is REPLACED at the branch-switching step by the orthogonal-direction constraint, written in the BVP inner-product form:

```
∫₀¹ ⟨u(t) - u_{prev}(t), φ₂(t)⟩ dt + (other-scalar terms) = Δs
```

where `φ₂(t)` is the BVP-discretized second null vector (computed by the analogous one-additional-backsolve recipe on the already-decomposed BVP Jacobian). The Floquet multiplier infrastructure (the boxed `A₀, A₁` blocks on p.757) gives the determinant-sign-change detector for the BVP case identically to how the algebraic determinant `det F_x` gives it for the finite-dimensional case.

### 5.5 Page/equation citation index for the load-bearing kernel-computation extraction

| Item | Part I location | Equation |
|---|---|---|
| Pseudoarclength continuation augmented system | p.502 col. 2 | (after §2.1 "Pseudoarclength continuation" heading) |
| Bordering algorithm for augmented LU | p.503 col. 1 | (2.5) and the L̂ Û decomposition |
| Right null vector of a singular `A` (one extra backsolve) | p.503 col. 2 right | (`PL̂ÛQ φ = 0` recipe) |
| Left null vector ψ (one extra backsolve) | p.504 col. 1 | (`A* ψ = 0` recipe) |
| Algebraic Bifurcation Equation `c₁₁α² + 2c₁₂αβ + c₂₂β² = 0` | p.505 col. 2 | (2.7) |
| Practical computation of bifurcation direction (φ₁ ≡ ẋ₀, φ₂ ⊥ φ₁) | p.505 col. 2 right | (last 2 paras of §2.2) |
| Orthogonal direction method (ε scaling = Δs) | p.506 col. 1 | (last para of "Practical branch switching") |
| Branch-switching corrector setup | p.506 col. 1 | (2.8) |
| Detection (`det F_x` sign-change) | p.506 col. 2 | (Theorem [Keller 1987] + Eq. 2.9 secant iteration) |

| Item | Part II location | Equation |
|---|---|---|
| BVP discretization with integral constraints + pseudoarclength | p.752 §1.4 | (1.9)-(1.13) |
| Bordered-block-bidiagonal Jacobian + condensation of parameters | p.753-754 | (Figs 1.3-1.5) |
| Periodic-solution computation, integral phase condition | p.755 §2.1 | (2.6), (2.8), (2.9) |
| Floquet multiplier extraction from condensed boxed blocks | p.757-758 | (`-A₁⁻¹ A₀` eigenvalues, Fig. 1.5 boxed blocks) |
| Branch switching inheritance from Part I | p.758 col. 2 last para | (verbatim quote in §5 above) |
| Infinite-period heteroclinic BVP (manifold projection) | p.759 §3.1 | (3.1)-(3.6) |

## 6. Comparison vs the project's existing Phase 1 P1.4 approach

The project shipped `src/cyclerfinder/genome/asymmetric_branch.py` at commit `4054f86` with the function `_select_saddle_center_eigenvector` + `branch_at_saddle_center`. The approach (from the read of the 311-line file):

1. **Build the monodromy matrix** from the parent's converged periodic orbit (line 265: `mono = monodromy(...)`).
2. **Pick the marginal eigenvector** of the monodromy via a three-step filter:
   - Exclude the 2 trivial-pair eigenvalues closest to +1 (the energy + time-translation pair).
   - Exclude the 2 strongly unstable eigenvalues by largest `|log|λ||` (the primary saddle pair).
   - From the remaining 2 (the **secondary** pair), pick the one closest to +1.
3. **Take the real part of the corresponding right-eigenvector** (or the imag part if real part is degenerate, line 161-167), unit-normalize, call this `v`.
4. **Perturb the parent IC**: `perturbed = parent_state0 + sign * epsilon * v` for both `sign = ±1`.
5. **Re-correct** via `correct_general_periodic_3d` in full-asymmetric mode (6D state closure at T) — the standard CR3BP single-shooting corrector.

### 6.1 Side-by-side with Doedel-Keller-Kernévez 1991 §2.2

| Doedel 1991 step | Project step | Same / Different? |
|---|---|---|
| Compute null vector `φ₂` of augmented BVP Jacobian via one extra backsolve | Compute right-eigenvector of monodromy via `np.linalg.eig` | **Mathematically equivalent for the periodic-orbit case, computationally different** — see 6.2 |
| Orthogonalize: `φ₂ ⊥ φ₁ = ẋ₀` | (no explicit orthogonalization step) | **Different** — see 6.3 |
| ε scaling: `Δs` adaptive, problem-dependent | `epsilon = 1e-3` default + signed both ways | **Same essential idea, different parameter management** |
| Corrector: pseudoarclength step with orthogonality constraint | Standard single-shooting corrector at fixed period | **Different** — see 6.4 |

### 6.2 The null vector ↔ monodromy eigenvector equivalence

For a **periodic-orbit BVP** `u'(t) = T f(u, λ), u(0) = u(1)` (Part II Eq. 2.4 + 2.5), the linearization about a converged orbit is `v'(t) = T A(t) v(t), v(0) = v(1)`. The monodromy matrix `M = V(1)` (where `V(t)` is the fundamental solution matrix `V'(t) = A(t) V(t), V(0) = I` — Part II Eq. on p.756) has its eigenvalues as **the Floquet multipliers**. A right-eigenvector `v` of `M` with eigenvalue `+1` is exactly a non-trivial solution of `(M - I) v = 0`, i.e. a non-trivial periodic perturbation tangent to a one-parameter family branching off the orbit.

The **BVP Jacobian** at a bifurcation point is the operator `L v = v' - T A(t) v` acting on functions satisfying `v(0) = v(1)` (Part II Eq. 2.11a-b, p.756); the variation-of-parameters formula (Eq. 2.13) shows `v(t) = V(t) [v(0) + τ ∫₀ᵗ V⁻¹(s) φ(s) ds]`. The null space `N(L)` of `L` (without the inhomogeneous τ term) is `{V(t) v(0) : v(0) ∈ N(M - I)}` — i.e. the null space of `L` is in 1-to-1 correspondence with the null space of `M - I`, and the BVP Jacobian's null vectors are the Floquet eigenvectors propagated forward by the state-transition matrix.

**So Doedel's `φ₂` ↔ the project's right-eigenvector of M closest to +1 (among the secondary pair) are the SAME object** for the periodic-orbit case, up to the time-propagation `V(t)`. The project's `np.linalg.eig` extracts it directly from the monodromy; Doedel's "one extra backsolve on the LU-decomposed BVP Jacobian" extracts it indirectly from the discretized BVP Jacobian. Both are valid; the monodromy approach is the standard astrodynamics convention (Howell, Grebow, others) and avoids the AUTO BVP infrastructure entirely. **The math is identical; the data structure is different.**

### 6.3 The orthogonalization gap

Doedel **requires** `φ₂ ⊥ φ₁ = ẋ₀` (the tangent to the given branch — p.505 right col.). The project's `_select_saddle_center_eigenvector` does NOT orthogonalize against the given family's tangent direction; it picks the monodromy eigenvector closest to +1 in the secondary pair and uses it directly.

**This is the substantive methodological difference.** Implications:

1. **In the (3,2) Earth-Moon C32 saddle-center case the project ran** (commit `4054f86`, eps=5e-4 at i=124), the marginal eigenvector points in the (z, ż) direction with the in-plane components (x, y, ẋ, ẏ) "near zero" (per the asymmetric_branch.py docstring lines 22-25). The PARENT family is the planar (2,1) cycler family, whose tangent direction `ẋ₀` is IN-PLANE — all 6 components in the (x, y, ẋ, ẏ) subspace, with `(z, ż)` components zero. So `φ_eig` (the project's eigenvector) and `ẋ₀` (the parent tangent) are **automatically orthogonal** because they live in disjoint coordinate subspaces.

2. **In the general (high-codimension, kernel rank > 1, fold-pitchfork interaction) cases Phase 2 will hit**, the marginal eigenvector and the parent tangent are NOT automatically orthogonal. The "perturbation lands back on the parent" failure mode the asymmetric_branch.py docstring (lines 219-222) warns about — "Too small: corrector lands back on the parent" — is **precisely the failure that Gram-Schmidt against `ẋ₀` would eliminate**. The Doedel orthogonalization is the principled fix for this failure mode.

3. **The "(planar manifold is invariant if the eigenvector is in the (z, ż) direction)" comment** (asymmetric_branch.py line 220-221) is a SPECIAL-CASE rationale that the project's eigenvector happens to be transverse. In general, the rationale is the Doedel orthogonalization.

**Recommendation**: Phase 2 should add an explicit Gram-Schmidt step `v_orth = v - (v · ẋ₀) ẋ₀; v_orth /= ||v_orth||` between the `_select_saddle_center_eigenvector` call and the `perturbed = parent_state0 + sign * epsilon * v` perturbation. The parent tangent `ẋ₀` is already in hand from the family tracer (it's the predictor direction for the next pseudoarclength step). This is a 2-line addition that promotes the project's eigenvector-perturbation pattern from a Phase-1-special-case shortcut into the Doedel-Keller-Kernévez 1991 canonical recipe.

### 6.4 The corrector difference

Doedel uses a **pseudoarclength corrector** with the orthogonal-direction constraint `(x₁ - x₀)* φ₂ - Δs = 0` (Part I Eq. 2.8) — i.e. the corrector ENFORCES that the converged transverse-branch member lies at arclength `Δs` from the bifurcation point along the `φ₂` direction.

The project uses a **standard CR3BP single-shooting corrector** at fixed period (modulo period being one of the 7 free variables in FREE_VARS_FULL_ASYMMETRIC). The corrector's residual is 6D state closure at T; there is NO pseudoarclength constraint on the perturbation magnitude.

**Implications**:

1. **The project's corrector may walk AWAY from the perturbation direction.** If the corrector lands at a state that happens to satisfy 6D closure but is closer to the parent than `epsilon * v` would suggest, that's not a failure (it's still a valid family member), but it's not a controlled `Δs` step either.

2. **The project's "try both signs ±epsilon" search** (line 273 `for sign in (+1, -1)`) is needed precisely because the corrector has no signed-direction constraint. Doedel's pseudoarclength constraint signs the direction unambiguously via the choice of `Δs > 0`.

3. **For the Phase 1 saddle-center demonstration the project shipped, this difference is harmless** — both signs converge or only one does, and a converged orbit with topology change is the Phase 1 exit criterion. The pseudoarclength constraint matters when Phase 2 wants to **continue along** the transverse branch in a controlled way; the project's family tracer (`cr3bp_3d_family_tracer.py`) ALREADY does pseudoarclength continuation, so this gap is naturally closed by re-handing the converged branched orbit BACK to the family tracer for the next step.

**Recommendation**: keep the project's corrector unchanged; the gap closes naturally at Phase 2's family-tracer handoff. Document this in the Phase 2 design doc.

### 6.5 Verdict on Phase 1 P1.4 alignment

**The project's Phase 1 P1.4 is methodologically aligned with Doedel-Keller-Kernévez 1991 §2.2 to within the "no explicit orthogonalization" simplification documented in §6.3.** The simplification was sound for the (3,2) C32 special case (eigenvector and tangent in disjoint subspaces); it is a known-edge-case shortcut, not a load-bearing flaw. Phase 2 should add the 2-line Gram-Schmidt orthogonalization to promote the substrate to the canonical recipe before hitting general-codimension cases.

**The Phase 1 agent built the substrate correctly without needing the Doedel reference** — the underlying eigenvector-perturbation pattern is standard astrodynamics (Howell, Grebow) and the Doedel formulation is the formal-bifurcation-theory crystallization of the same idea. The reference is the **canonical citation** for the substrate and the **principled extension** for Phase 2 edge cases; neither obligation invalidates the Phase 1 shipping.

## 7. BVP corrector with integral constraints (Part II §1) — invariant-set continuation capability

Part II §1.4 (p.752 Eqs 1.9-1.13) is the explicit AUTO BVP corrector with `n_q` integral constraints `∫₀¹ q(u(s), µ, λ) ds = 0` (Eq. 1.3, p.745). The discretized form (Eq. 1.11):

```
Σⱼ Σᵢ ω_{j,i} qₖ(u_{j-(i/m)}, µ, λ) = 0,   k = 1, ..., n_q
```

uses Lagrange quadrature weights `ω_{j,i}` to discretize the integral over each mesh subinterval. The resulting discrete system has `m n N + n_b + n_q + 1` equations in `m n N + n_µ + 1` unknowns; the integral-constraint rows are the last `n_q` rows of the bordered-block-bidiagonal Jacobian (Fig. 1.3, p.753).

**Project capability gap**: the existing single-shooting + multi-shooting correctors in `src/cyclerfinder/genome/multi_shooting.py` and `src/cyclerfinder/search/cr3bp_general_periodic_3d.py` handle **boundary conditions** (perpendicular-axis-crossing, full-state-closure-at-T) but NOT **integral constraints**. There is no native facility to ENFORCE an integral relation like "this orbit's mean Jacobi value over one period equals C_target" or "this orbit's instantaneous-Hamiltonian integral matches the BCR4BP 7-integrals identity" as a corrector residual.

**The Doedel 2003 §4 Eq. 18** (the 7-integrals BCR4BP family-tracing the brief asks about) is one of these integral relations. Per the existing `2026-06-17-digest-doedel-2003.md` (which I scanned at the start), Doedel et al. 2003 §4 demonstrates this for the planar general 3-body figure-8 homotopy: an unfolding parameter `λ` multiplied by the gradient of each first integral is added to the state equations, and the integral constraint `∫₀¹ λ(t) Iⱼ(u(t)) dt = 0` is added as an integral relation (Eq. 18 form). The BVP corrector with integral constraints from Part II §1.4 is THE machinery that makes this work.

**Recommendation for Phase 3 (the #334 BCR4BP family-tracking pattern)**:

The natural extension is to add a `bvp_integral_constraint.py` module that wraps the existing `cr3bp_general_periodic_3d` corrector with:
1. A `q_residuals: callable[[state, time], np.ndarray]` argument returning the integrand values for each integral constraint.
2. A trapezoid/Simpson quadrature over the integration trajectory (the project already integrates trajectories for state-closure; the same trajectory yields the integral via numpy quadrature on the dense output).
3. Additional rows in the Newton residual + Jacobian for each integral constraint.

This is a Phase 3 build, not Phase 2. The Phase 1 substrate does not need it. The capability gap is real and is what Part II §1.4 + 2003 §4 together resolve; the project does not currently have it.

## 8. AUTO software architectural lessons — is the project's cr3bp_general_periodic_3d.py Auto-style?

### 8.1 The AUTO architectural pattern (Part I §1.1 + §2 + Part II §1.4)

AUTO's organizing principle is: **all bifurcation-theoretic computations are reduced to the same Newton-bordered-linear-solver kernel**. Specifically:

- One **inflated** nonlinear system (the augmented continuation equations) defines each task: continuation past a fold (§2.1), fold continuation in two parameters (§2.3 Eq. 2.10), Hopf continuation (§2.4 Eq. 2.18), branch switching (§2.2 Eq. 2.8), optimization (§3 Eq. 3.9).
- One **bordered LU decomposition** (Part I p.503 right col., the `(L 0 ; β* 1)(U γ ; 0* δ)` factorization) is the linear-solve kernel. Every Newton step reduces to this.
- The **same scalar test functions** (determinant sign-changes, eigenvalue real-part crossings, Lagrange-multiplier sign-changes) drive ALL bifurcation detection.
- The **same secant-iteration bracketer** (Eq. 2.9) locates the precise bifurcation/extremum from the test-function zero-crossing.

The result: AUTO is a small Newton-bordered kernel + a library of inflated-system definitions + a test-function dispatcher. New bifurcation types are added by writing new inflated-system definitions; the linear-solve, the Newton outer loop, the step-size adaptation, and the bracketer are reused.

### 8.2 Is the project's cr3bp_general_periodic_3d.py Auto-style?

**Partially.** The 594-line `cr3bp_general_periodic_3d.py` is a generic single-shooting corrector with **configurable free variables and configurable residual indices** (`FREE_VARS_FULL_ASYMMETRIC`, `RESIDUAL_FULL_STATE_AT_T`, etc.) — this is exactly the Auto "inflated system as configurable arguments" pattern. The Newton outer loop, the damped line search, and the residual evaluation are factored once and reused for symmetric / planar / 3D / asymmetric variants.

**But the project does NOT have AUTO's full architectural decomposition**. Specifically:

1. **The pseudoarclength augmentation is in a SEPARATE module** (`cr3bp_jacobi_arclength.py`, `cr3bp_3d_family_tracer.py`) rather than being one more configurable inflated-system option on the same corrector.
2. **The Floquet bifurcation detection is in a SEPARATE module** (`bifurcation_detector.py` per the search/ inventory).
3. **Each new bifurcation type tends to get its own module** (`family_switch.py` for period-doubling, `asymmetric_branch.py` for saddle-center, etc.) rather than being a new configurable inflated-system on the unified corrector.

**Should the project's pattern BE Auto-style?**

Yes for Phase 2-N, with a caveat. The AUTO unification pays off when (a) there are many bifurcation types to add and (b) they all benefit from the same step-size adaptation, line search, and convergence diagnostics. The project IS hitting condition (a) — period-doubling, saddle-center, fold continuation in (mu, C), Hopf for BCR4BP, two-parameter folds for system-swap, etc. — so the AUTO refactor pays for itself by Phase 3-4.

**The caveat**: the project's per-system corrector tuning (the `require_monotone_decrease`, `independent_tol`, the perpendicular-axis-crossing vs full-state-closure distinction) is more aggressive than AUTO's because CR3BP single-shooting is ill-conditioned in ways AUTO's BVP-collocation isn't. A literal port of AUTO's architectural pattern that lost the existing tuning would regress. The right move is: **a unified `bifurcation_continuation.py` module that wraps the existing tuned correctors as configurable inflated-system definitions, NOT a rewrite of the correctors themselves**.

**Concrete recommendation for Phase 2**: build a `cyclerfinder/genome/bifurcation_continuation_kernel.py` module that:
1. Takes a tuned corrector callable + a tuned independent-cross-check callable.
2. Layers (a) the parent family's pseudoarclength predictor, (b) the bifurcation test function `det F_x` or Floquet-multiplier-distance-from-+1, (c) the secant bracketer (Part I Eq. 2.9), (d) the orthogonal-direction-method branch-switch step (Part I Eq. 2.8) on top.
3. Returns a converged branched-orbit object compatible with the family-tracer's input format.

This is a 1-2 week Phase 2 build. Not Phase 1.

## 9. Phase 2 discovery sweep implications — concrete recommendations

Given the Phase 1 substrate just shipped at commit `e98e49b`:

### 9.1 Add the Gram-Schmidt orthogonalization (1-line fix, do this first)

In `src/cyclerfinder/genome/asymmetric_branch.py::branch_at_saddle_center`, between the `pick = _select_saddle_center_eigenvector(mono)` call and the `perturbed = ...` assignment, add (assuming the parent tangent `tangent_dx0` is available from the family tracer's predictor):

```python
v_orth = v - np.dot(v, tangent_dx0) * tangent_dx0
v_orth /= np.linalg.norm(v_orth)
v = v_orth
```

If `tangent_dx0` is not currently passed in, plumb it through from the family-tracer caller. This is the Doedel-canonical orthogonal direction method (Part I p.506) and eliminates the "perturbation lands back on the parent" failure mode in general-codimension cases.

### 9.2 Verify against Phase 1's C32 case before scaling

Before running Phase 2's sweep, re-run the C32 i=124 case with the orthogonalized eigenvector and confirm the converged branched orbit is bit-identical (within Newton tolerance) to the Phase 1 shipped result. This is the regression guard — the orthogonalization is mathematically the same operation for the C32 case (eigenvector and tangent in disjoint subspaces), so a converged-orbit delta would indicate a plumbing bug.

### 9.3 Discover bifurcations via determinant-sign-change AND Floquet-distance-to-+1

The project's existing `bifurcation_detector.py` (per the search/ inventory) is the right machinery for Floquet-distance-to-+1 watching. To make Phase 2 robust against the high-codimension cases the brief mentions (kernel rank > 1, fold-pitchfork interactions), ALSO watch the augmented-Jacobian determinant sign (Part I Eq. 2.9 secant iteration). The two detectors fail differently: the Floquet detector misses bifurcations where the multiplier touches +1 but doesn't cross (the "fold-pitchfork" case); the determinant detector catches those because `det F_x` changes sign at any bifurcation but does NOT change sign at a simple fold (Part I p.506 Theorem [Keller 1987]). Running both detectors and taking the union eliminates a Phase 2 false-negative class.

### 9.4 Handle kernel rank > 1 by fanning out

At a kernel-rank-2 bifurcation (a "high-codimension" point), there are TWO transverse branches (not one). The Part I §2.2 case-(ii) framework (p.505 left col., the second example with `f(u, λ) = (λ - u₁² - u₂² ; u₁ u₂)*` at `(0, 0, 0)`) handles this: the algebraic bifurcation equation Eq. 2.7 has discriminant > 0 giving TWO distinct real solutions `(α₁, β₁)` and `(α₂, β₂)`, each producing a different transverse branch. Phase 2 should detect kernel rank explicitly (via SVD of the augmented Jacobian, counting singular values below a threshold) and fan out to BOTH transverse branches when rank-2 is detected. This is a Phase 2 capability, not Phase 1.

### 9.5 Document the Phase 1 P1.4 approach as a "monodromy-shortcut" variant of Doedel 1991b §2.2

Update `src/cyclerfinder/genome/asymmetric_branch.py`'s module docstring to cite Doedel-Keller-Kernévez 1991 Part I §2.2 as the canonical reference for the algorithm, with a brief note that the project uses the monodromy-eigenvector path rather than the augmented-BVP-Jacobian-null-vector path because (a) they are mathematically equivalent for the periodic-orbit case (§6.2 above) and (b) the project does not currently use AUTO's BVP collocation infrastructure. This is documentation discipline; no code change.

## 10. Phase 3 μ-scaling implications — system-swap continuation (#334)

### 10.1 The Part II framework for system-swap

System-swap continuation (#334: Earth-Moon → Sun-Earth, Sun-Earth → Sun-Mars, etc.) is mathematically continuation in the mass-ratio parameter `µ` (the CR3BP's only system parameter). The Part II framework (§1.4 BVP corrector + §2.1 periodic-orbit BVP + §2.3 integral phase condition) is directly applicable: treat `µ` as one of the `n_µ` problem parameters and use pseudoarclength continuation in `µ`, the periodic orbit's tangent direction, and the period `T` simultaneously.

**The critical insight from Part II §1.3 (Multiplicity of Solutions, p.750)**: the **multiplicity theorem** (Beyn-Doedel 1981) guarantees that for sufficiently small mesh `h`, every isolated solution of the continuous problem at `µ = µ_target` is approximated by a unique discrete solution. This means **system-swap continuation is guaranteed to converge** when the orbit at the source `µ` is regular and isolated, the mesh is fine enough, and the continuation step `Δµ` is small enough. The Phase 3 system-swap design doc should cite this theorem as the convergence guarantee.

### 10.2 The integral-constraint capability is what enables BCR4BP system-swap

BCR4BP (#334's target system) is NOT a CR3BP — it adds the Sun's perturbation as a fourth body. The standard CR3BP integrals (Jacobi C, energy E) are NOT conserved in BCR4BP, but the BCR4BP-extended Hamiltonian is conserved along trajectories AT FIXED EPOCH; the family-tracing constraint is then an integral relation rather than a pointwise relation. **Part II §1.4 integral-constraint BVP corrector is the methodology that handles this directly** — without it, the project cannot trace BCR4BP families in epoch.

**Recommendation for Phase 3 design (#334)**: Phase 3.1 should be the integral-constraint corrector build (Section 7 above); Phase 3.2 should be the BCR4BP system-swap continuation using the Phase 3.1 corrector. This sequencing is forced by the Part II framework — there is no shortcut around it.

### 10.3 The two-parameter fold continuation framework extends directly

The #334 BCR4BP family-tracking pattern (per the memory anchor `2026-06-17-334-bcr4bp-system-swap.md`) requires tracking how the family's fold-in-C moves as `µ` changes — this is two-parameter fold continuation in `(µ, C)`. Part I §2.3 Eq. 2.10a-b gives the extended system for fold continuation in two parameters:

```
f(u, λ, µ) = 0
f_u(u, λ, µ) φ = 0
φ* φ₀ - 1 = 0
```

This is the exact pattern Phase 3 needs. The Phase 3 design doc should cite Eq. 2.10a-b as the canonical formulation.

## 11. KNOWN_CORPUS impact — methodology bibliography

These two papers are **methodology references**, not catalogue rows. They contain no orbit data, no Jacobi constants, no period values, no IC tuples that would be eligible for catalogue admission. Their impact is purely on the project's algorithmic foundations.

**Recommendation**: build a `docs/refs/methodology-bibliography.md` index that catalogues the canonical methodology references with one-line descriptions of what each provides. Entries (proposed schema):

| Reference | Provides | Used in project at |
|---|---|---|
| Keller 1977, "Numerical solution of bifurcation and nonlinear eigenvalue problems" | Original pseudoarclength continuation algorithm | `cr3bp_jacobi_arclength.py`, `cr3bp_3d_family_tracer.py` |
| **Doedel-Keller-Kernévez 1991 Part I** (this digest) | Algebraic Bifurcation Equation, orthogonal-direction branch switching, fold continuation, Hopf continuation, optimality systems | `asymmetric_branch.py` (Phase 1 P1.4), Phase 2 expansion |
| **Doedel-Keller-Kernévez 1991 Part II** (this digest) | BVP corrector with integral constraints, periodic-orbit phase condition, Floquet multipliers via condensed boxed blocks, heteroclinic BVP | Phase 3 BCR4BP system-swap (#334) |
| Doedel-Paffenroth-Keller-Dichmann-Galán-Vioque-Vanderbauwhede 2003 | Unfolding-parameter regularization for conservative systems, AUTO applied to CR3BP family network | Already digested (`2026-06-17-digest-doedel-2003.md`) — methodology root |
| (others: Howell-Campbell 1999, Hénon, Floquet-multiplier original papers, etc.) | (to be filled in by future digests) | (to be filled in) |

This index is NOT a catalogue and contains no admission rows. It is a methodology cross-reference for future work; it eliminates the "every time we hit a bifurcation question, we have to re-derive what reference to consult" cost.

**No catalogue admission rows; no validation tier changes; no negative-results-registry entries.** These are pure methodology references.

## 12. Action items for the parent

The five items below are concrete carry-forwards from the digest. None requires immediate Phase 1 changes (commit `e98e49b` ships as-is); all are Phase 2 or later.

1. **Phase 2.1 (next): Add 1-line Gram-Schmidt orthogonalization to `asymmetric_branch.py::branch_at_saddle_center`** per §9.1. Promotes the substrate to Doedel canonical recipe; eliminates "perturbation lands back on parent" edge case for general-codimension Phase 2 sweep targets. Owner: Phase 2 lead.

2. **Phase 2.2: Add determinant-sign-change detector ALONGSIDE Floquet-distance-to-+1** per §9.3. Two detectors, union of triggers, no false-negative class. Owner: Phase 2 lead.

3. **Phase 3.1 (gates Phase 3.2): Build `cyclerfinder/search/bvp_integral_constraint.py`** per §7. Loads `q(u, t) → R^{n_q}` integral residuals into the existing single-shooting corrector via trapezoid/Simpson quadrature on the dense trajectory output. This is the load-bearing capability for BCR4BP system-swap (#334) and the 2003 7-integrals demonstration. ~3-5 days of focused work. Owner: Phase 3 lead.

4. **Phase 3.2: Cite Part I Eq. 2.10a-b in the #334 BCR4BP system-swap design doc** per §10.3. Canonical formulation for two-parameter fold continuation in (µ, C). Documentation discipline, not code.

5. **Build `docs/refs/methodology-bibliography.md`** per §11. One-page methodology cross-reference, schema in §11. Cheap, immediate, eliminates a recurring "which reference?" cost. Can be done by the next available agent independent of Phase 2 / Phase 3 scheduling. Estimated 1-2 hours.

### Honest negative for the parent

**The brief's framing — that Part II would contain a deep-detail kernel-computation algorithm distinct from the existing project's Phase 1 P1.4 module — is a 75%-correct framing.** The deep-detail algorithm IS canonical (Part I §2.2, the Algebraic Bifurcation Equation + orthogonal direction method, p.505-506), but Part II's contribution is **the lifting of Part I's algorithm to the BVP-discretized infinite-dimensional case via the orthogonal collocation + condensation-of-parameters Jacobian (Eqs 1.9-1.13)** — NOT a fundamentally new algorithm. The Phase 1 agent built `asymmetric_branch.py` correctly without the reference because the underlying eigenvector-perturbation pattern is standard astrodynamics; the Doedel formulation is the formal-bifurcation-theory crystallization of the same idea, with **one substantive methodological refinement** (the φ₂ ⊥ ẋ₀ orthogonalization, §6.3) that Phase 2 should adopt.

The reference is the **canonical citation** for the substrate and the **principled extension path** for Phase 2-3 edge cases. Both roles justify the read; neither role obligates a Phase 1 rebuild. Commit `e98e49b` ships as-is.

---

**Digest discipline**: every numeric value, equation citation, and algorithmic claim traces to a specific page or section of the two PDFs. No fabrication; no values from project code without citation; both papers read cover-to-cover (every page, every reference list). The two inferred-but-not-source-derivable DOIs in §1 are flagged as inferred. Reading time: ~3 hours.
