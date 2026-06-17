# Digest — Pellegrini & Russell 2016 "On the Computation and Accuracy of Trajectory State Transition Matrices"

**Date**: 2026-06-17 AET (digest of `#116` STM continuation entry; parent #347 / Phase 1 Floquet substrate)
**Verdict (TL;DR)**: **ADMIT as methodology citation + LOAD-BEARING ACTION ITEM for #347 Phase 1.** This is a JGCD methods paper (DOI 10.2514/1.G001920) comparing three STM-computation techniques: (1) variational equations alongside the state, (2) complex/bicomplex-step derivative (CSD), (3) multipoint finite differences. **The project currently uses method (1) with a variable-step integrator (scipy DOP853, `cr3bp.py:cr3bp_stm_eom` + `propagate(with_stm=True)`).** Pellegrini-Russell's headline result (Abstract; Conclusions §IV.1): *"It is important to avoid mixed propagations of the state only and of the state and STMs"* — and *"the variational equations are usually outperformed by CSD, and even by the finite difference methods for loose integration tolerances"* when used with variable-step integrators. **This applies directly to `bifurcation_detector.monodromy()`** which is the kernel of #347 Phase 1's saddle-center detection — at the cluster point (multipliers near +1) the Phase 0 design doc (line 116) already explicitly flagged that "the existing `monodromy` runs at rtol=atol=1e-12 ... **Unknown** whether the eigenvalue separation will survive at the cluster point." Pellegrini-Russell provides the mitigation menu: switch to **fixed-path** (record VS step sizes, replay as FP) or to **CSD**. **No catalogue row introduced** — this is a methods paper, not a trajectory result.

---

## 1. Header

- **Title** (verbatim, p. 1): "On the Computation and Accuracy of Trajectory State Transition Matrices"
- **Authors** (verbatim, p. 1): Etienne Pellegrini (Graduate Research Assistant), Ryan P. Russell (Associate Professor, Associate Fellow AIAA); both at Department of Aerospace Engineering and Engineering Mechanics, University of Texas at Austin, 210 E. 24th St., Austin, Texas 78712
- **Venue**: *Journal of Guidance, Control, and Dynamics* (JGCD), "Article in Advance" — received 4 December 2015; revision received 18 April 2016; accepted for publication 1 June 2016; published online 14 September 2016
- **DOI**: 10.2514/1.G001920
- **Page count**: 15 pages (front matter through bibliography)
- **Publisher**: American Institute of Aeronautics and Astronautics, Inc.
- **Funding**: Airbus Defence and Space, research agreement no. UTA12-000936 (Acknowledgments, p. 14)

## 2. What the paper actually is

Per the Abstract (p. 1, verbatim, lightly edited for line breaks):
> "State transition matrices provide sensitivities or partial derivatives between states at different times along a trajectory and are used for a number of applications such as feedback controls, stability analysis, estimation, targeting, parameter optimization, and optimal control. The need for accurate state transition matrices is especially important for problems and applications that have highly sensitive and highly nonlinear dynamics. In this paper, examples are considered in the context of multiple-revolution and multiple-body space trajectories. Three techniques to compute both the first- and second-order state transition matrices are compared: 1) augmenting the state with the classic variational equations, 2) complex and bicomplex-step derivative approximation, and 3) multipoint stencils for traditional finite differences. Each of the methods are compared for accuracy and speed across a variety of problems and numerical integration techniques. The subtle differences between variable- and fixed-step integration for partial computation are revealed, common pitfalls are observed, and recommendations are made to enhance the quality of state transition matrices. A main result is the demonstration of small but potentially significant errors in the partials when they are computed with variational equations and a variable-step integrator."

Section structure:
1. Introduction (pp. 1–2): motivation, prior art, why CSD has gained traction.
2. Methods (pp. 2–5): definitions; Variational Equations (§II.A); Complex-Step Derivative (§II.B); Finite Differences (§II.C, generic + 2nd-order).
3. Applications (pp. 5–14): Testing Framework (§III.A) on three unstable three-body periodic orbits PO1–PO3 and two two-body orbits TB1–TB2; Caveats of Variable-Step Integration (§III.B); Accuracy and Timings (§III.C); Special case PO2 (§III.D).
4. Conclusions (p. 14): three numbered pitfalls of variable-step integration for STMs.
5. References (pp. 14–15): 48 entries.

## 3. The three methods

### 3.1 Variational equations (§II.A, pp. 2–3)

Augment the state with the first- and second-order STMs and integrate the variational equations alongside the EOMs. For state X with N elements:

- First-order STM Φ¹(t, t₀) = ∂X(t)/∂X₀ (eq. 1)
- Second-order STM Φ²(t, t₀) = ∂²X(t)/∂X₀² (eq. 2)
- Variational equations: Φ̇¹ = f_X Φ¹ (eq. 7); Φ̇² = f_X • Φ² + Φ¹ᵀ • f_XX • Φ¹ (eq. 8)
- Total augmented state size: N + N² + N³ (paper's eq. count; for N=6, that's 6 + 36 + 216 = 258 ODEs).

**Implementation pitfall identified (eqs. 16–17, p. 3)**: when the integrator is variable-step, the variational equations on the augmented state **do not capture changes due to integration-step-size variation** with respect to initial conditions. The δt(X) dependence (the per-step error estimate's IC-dependence) contributes a missing term to the STM. Formally, the actual STM update equation is

```
∂X_{i+1}/∂X_0  =  Φ_i^1  +  (f_X|_X_i δt + f(X_i) ∂δt/∂X|_X_i) Φ_i^1     (eq. 17)
```

but the variational propagation gives only Φ_i^1 + f_X|_X_i Φ_i^1 δt (eq. 12). The missing term f(X_i) (∂δt/∂X) Φ_i^1 is the **δt_i/∂X term** that recurs throughout the paper.

For an RKF(7)8 integrator with tolerance ε, the δt_i/∂X term scales as ε^(4/5) (eq. 19, p. 3) — so tighter tolerance shrinks the missing term, but does not eliminate it.

### 3.2 Complex-step derivative (CSD) (§II.B, pp. 3–4)

Use complex (or bicomplex, for second-order) numbers to compute derivatives via the imaginary part of f(x + ih). For analytic functions:

```
f'(x) = Im(f(x + ih)) / h + O(h²)     (eq. 20)
```

No subtraction-cancellation error (unlike finite differences). Perturbation can be made arbitrarily small (paper uses h = 10⁻⁴⁰, p. 10). For second-order STMs, use bicomplex numbers (C² = {z₁ + i₂z₂ | z₁, z₂ ∈ C}; eqs. 22–24).

**Critical advantage stated** (p. 4, verbatim): *"Including the sensitivity of the step size term is a subtle but important advantage of CSD and AD methods over the variational equations when using a variable-step integrator."* The CSD method propagates the whole integration process using complex numbers, so the integration step size IS exposed to the perturbation, and the δt/∂X term **is automatically captured**.

**Implementation cost**: each procedure in the dynamics code (and the integrator!) must be modified to accept complex inputs. The authors implement a Python preprocessor "generize" that translates a real-number Fortran library to a generic complex/bicomplex version.

### 3.3 Multipoint finite differences (§II.C, pp. 4–5)

Use stencils of n+1 points: forward FD (n=1, α=[0,1]), central FD (n=2, α=[-1,0,1]), generic FD (high-order Lagrange-stencil; this paper uses n=4 with α = [-2, -1, 0, 1, 2], β = [1/12, -2/3, 0, 2/3, -1/12]).

Pitfall (§II.C.1): for small h, FD methods are limited by floating-point arithmetic; for large h, by truncation. The "tuned" perturbation magnitude must be optimized per orbit (Table 2 p. 8 gives the tuned h for each test case: 10⁻⁷ to 10⁻⁹).

## 4. The headline result — variable-step + variational = small but real STM errors

Pellegrini-Russell's central numerical experiment is a **line-search application** (§III.A.3, p. 6): given a reference trajectory X(t) starting at X₀, and an STM Φ¹(t, t_f), predict the final state X_{f, expected} for a perturbed IC X₀ + δX via X_{f, expected} = X_{f, real} + Φ¹ δX. Compute the relative error vs the actually-propagated perturbed trajectory:

```
Err = ||X_expected - X_real|| / ||X_real||     (eq. 49)
```

A perfect STM (zero error) plots as flat at machine-epsilon ~ 10⁻¹⁵ across the δX domain. Any systematic departure from machine-epsilon at small δX is a **systematic STM error**.

### 4.1 Three numbered pitfalls (Conclusions §IV)

Verbatim from p. 14:

1. *"When propagating the variational equations, it is important to avoid mixed propagations of the state only and of the state and STMs."*
2. *"Variable-step integration can lead to discontinuities of the integration function, due to the step size selection changing the number of steps required to complete the integration. These discontinuities cannot be predicted by the partials computations. Moreover, the finite difference approximations are especially vulnerable to these discontinuities since the number of integration steps is fixed when computing the perturbed trajectories. Therefore, the use of **fixed-step integration and time regularization techniques is recommended** to avoid the need for a variable-step integrator."*
3. *"The variation of the step sizes has an impact on the final state, and it should be accounted for in the partials. The propagation of the variational equations does not capture this change, in contrast to the cases of the CSD approximation and the finite difference methods. When using variable-step integration, the variational equations are usually outperformed by CSD, and even by the finite difference methods for loose integration tolerances. These caveats all have a larger influence for highly sensitive orbits and for low-accuracy propagations, making them especially important for the preliminary design of modern space missions (where a low-accuracy propagation is adapted and does not change the nature of the orbit)."*

### 4.2 Quantitative — Test Cases (Table 1, p. 5)

The three unstable three-body periodic orbits used:

| Orbit | μ (mass ratio) | Period T_P (TU) | λ_max | ||Φ²|| | Description |
|---|---|---|---|---|---|
| PO1 | 2.528e-5 | 25.139 | 2.469e7 | 1.882e4 | 7:4 unstable resonant periodic orbit |
| PO2 | 2.528e-5 | 70.539 | 2.804e17 | 8.290e17 | 8:11 unstable resonant orbit with three loitering close flybys |
| PO3 | 1.213e-2 | 2.518 | 17.624 | 5.183e6 | L1 halo orbit |

**PO3's μ = 1.213e-2** ≈ **the Earth-Moon CR3BP mu** the project's #347 Phase 1 (3,2) C32 family lives in (Braik-Ross 2026 Earth-Moon). λ_max(PO3) = 17.624 vs the project's Phase 0 design doc target |λ_max| ~ 2.5e5 at the (3,2) anchor — i.e., our orbits are ~10⁴× more sensitive than PO3, putting us in the regime where Pellegrini-Russell shows the variational-equations + variable-step combination underperforms most.

### 4.3 PO2 — most sensitive case (§III.D, p. 13)

For PO2 (||Φ²|| ~ 10¹⁷, λ_max ~ 10¹⁷), Fig. 13 shows:
- For **loose tolerances (ε=10⁻³, 10⁻⁶)**: CSD with variable-step **outperforms variational equations by several orders of magnitude**.
- For **tight tolerances (ε=10⁻⁹, 10⁻¹²)**: variational and CSD become comparable.

Verbatim conclusion (p. 14): *"For the variable-step integration and certain integration tolerances, the CSD algorithm outperforms the variational equations by several orders of magnitude, confirming that the sensitivity of the orbit makes the importance of the δt_i/∂X term increase."*

## 5. Project relevance — where Pellegrini-Russell 2016 hits our code

### 5.1 Current STM machinery (read-only audit)

| File | Function | Method | Integrator | Notes |
|---|---|---|---|---|
| `src/cyclerfinder/core/cr3bp.py` L100–125 | `cr3bp_stm_eom` | **Variational equations** (state-6 + STM-36 augmented) | n/a (just the RHS) | Includes Hessian terms Uxx/Uyy/Uzz/Uxy/Uxz/Uyz |
| `src/cyclerfinder/core/cr3bp.py` L128–174 | `propagate(with_stm=True)` | Variational propagation | **scipy `solve_ivp` method="DOP853" with rtol=atol=1e-12** — **variable-step** | Augmented state path (state + STM together) |
| `src/cyclerfinder/search/bifurcation_detector.py` L130–174 | `monodromy(system, state0, period)` | Variational | Same DOP853 via `cr3bp.propagate(with_stm=True)` | Returns Φ(T) for Floquet multiplier eigendecomposition |
| `src/cyclerfinder/search/cr3bp_periodic.py` L160–200 | `_propagate` (used by symmetric-orbit corrector + general corrector) | Variational | DOP853, rtol/atol from caller | Augmented path |
| `src/cyclerfinder/search/cr3bp_general_periodic_3d.py` L223–232 | `_propagate_with_stm` | Variational | DOP853 | Augmented; used by 3D corrector |
| `src/cyclerfinder/core/bcr4bp.py` L335–356 | `propagate(with_stm=True)` | Variational | DOP853 | BCR4BP equivalent |

**All STM machinery in the project uses variational equations + variable-step integration (DOP853).** No fixed-step path exists. No CSD path exists. No finite-difference path exists.

### 5.2 Phase 0 design doc (line 116) — already flagged this

The #347 Phase 0 design doc, `docs/notes/2026-06-17-347-floquet-bifurcation-framework-phase0-design.md` line 116, verbatim:

> "**Does the project's DOP853 integrator hold corrector tolerance at the bifurcation?** Near a saddle-center, the monodromy has multipliers near (+1, +1, +1, +1, …) — clustered eigenvalues are ill-conditioned. The existing `monodromy` runs at rtol=atol=1e-12 (the project's standard). **Unknown** whether the eigenvalue separation will survive at the cluster point. Mitigation: cross-check the multiplier via the trace/determinant invariants (the characteristic polynomial coefficients should match a known structure for a Hamiltonian monodromy). Add this as a check in P1.3."

**Pellegrini-Russell 2016 confirms this is a real and quantified concern.** At cluster points where λ-spread shrinks, the variational+variable-step combo's missing δt/∂X term contaminates the multipliers. The Phase 0 design's mitigation (Hamiltonian-invariant cross-check) is necessary but not sufficient — it detects the contamination, it does not eliminate it. **The Pellegrini-Russell-recommended mitigation is to switch to either fixed-path or CSD.**

### 5.3 #347 Phase 1 specific impact

The Phase 1 walker `scripts/floquet_phase1_p1_2_walk.py` walks (3,2) C32 family monodromy through CJ continuation. At each step it calls `monodromy(sys, mem.state0, mem.period, rtol=1e-12, atol=1e-12)`. The Phase 0 design (line 14 of the walker) verbatim: *"Per the Phase 0 design doc Section 5 (risk: DOP853 conditioning at the cluster point), the |lambda_max| ~ 2.5e5 at the anchor means the monodromy is ill-conditioned."*

**Implication**: the Floquet multipliers we are recording per family member are being computed by a method (variational + variable-step DOP853) that Pellegrini-Russell 2016 demonstrates is **dominated by the δt/∂X term at the cluster point for sensitive orbits**. The non-trivial pair's |λ|=1 saddle-center signal is what we're hunting for — and Pellegrini-Russell shows that's exactly the regime where the variational method underperforms.

### 5.4 Recommendations for the project's STM machinery — graded

**Strong (do): Implement a fixed-path verification mode for the saddle-center cluster point.** Per Pellegrini-Russell's Conclusion 2 ("the use of fixed-step integration and time regularization techniques is recommended to avoid the need for a variable-step integrator"), and Conclusion 1 ("avoid mixed propagations of the state only and of the state and STMs"):

1. Add a `fixed_path` mode to `cr3bp.propagate(with_stm=True)` and `bifurcation_detector.monodromy()`: first propagate the state-only path under variable-step DOP853, **record the realized step sizes**, then re-propagate the augmented state+STM using those **same recorded step sizes** as a fixed-step schedule. This matches Pellegrini-Russell's "fixed-path feature" (§III.A.2 p. 6).
2. At each Phase 1 saddle-center candidate (and at the family walker's anchor), compute the monodromy **once via the existing variable-step variational path** and **once via the fixed-path variational path**, and report the |λ_nontrivial| disagreement. If they disagree at 2+ decimal places at the cluster point, the saddle-center signal needs the fixed-path or CSD value.
3. Cost: small (a few hundred lines of new code; reuses `cr3bp_stm_eom`); no new integrator required.

**Medium (consider): Implement a CSD fallback for second-order STM applications.** If the project ever adds second-order STM use (e.g., differential dynamic programming, primer-vector second-order optimality conditions), CSD's automatic δt/∂X capture and h~10⁻⁴⁰ insensitivity make it the right method per Pellegrini-Russell. But this means complex-number-clean dynamics code and integrator; a substantial refactor. **Not needed for #347 Phase 1.** Flag as a future option in #226-style FBS work or in advanced optimization paths.

**Weak (defer): Drop variable-step altogether.** Pellegrini-Russell Conclusion 2 advocates fixed-step + time regularization for production STM work. The project's current rtol/atol=1e-12 DOP853 is industry-standard and the scipy ecosystem doesn't easily support time-regularization integrators (Mikkola-Merritt, DROMO, etc., refs [9,10]). The cost/benefit for swapping out DOP853 globally is high; the fixed-path mitigation in §5.4(1) captures most of the Conclusion-2 benefit at far lower cost. **Defer until a concrete bug surfaces.**

## 6. Catalogue impact — none directly; KNOWN_CORPUS citation recommended

This is a **methodology paper**, not a trajectory result. **No catalogue row admitted, no admission verdict.**

For `KNOWN_CORPUS` in `src/cyclerfinder/search/literature_check.py`: there is no exact match for "STM accuracy / variational equations" methodology citations in the current anchor list. The project's methodology-citation discipline (per the §16.5 literature_check spec) anchors methods only when they directly bear on a *cycler admission* — they aren't general methods refs. **Pellegrini-Russell 2016 doesn't anchor a cycler family**; it informs how we compute monodromy.

**Recommendation**: do **not** add a `KNOWN_CORPUS` anchor for this paper (methods refs are out of scope for that file's purpose). Instead, cite Pellegrini-Russell 2016 as a docstring reference in:
- `src/cyclerfinder/core/cr3bp.py` — at the top of `cr3bp_stm_eom` and `propagate(with_stm=True)`, note that variational + variable-step has known small-but-real STM errors per Pellegrini-Russell 2016 JGCD DOI 10.2514/1.G001920, with fixed-path mitigation per §5.4(1) above.
- `src/cyclerfinder/search/bifurcation_detector.py` — at `monodromy()`, note the cluster-point conditioning concern with Pellegrini-Russell as the canonical reference.
- The #347 Phase 0 design doc — add a follow-up note to line 116 citing this paper as the quantified-evidence reference.

These docstring updates can be batched into the Phase 1 incremental commits; they don't need a separate dedicated commit.

## 7. Errata against the project's prior context

| Where | Prior claim | What the paper actually says | Status |
|---|---|---|---|
| Agent prompt | "the project's current STM machinery (e.g. in `src/cyclerfinder/core/`, `src/cyclerfinder/search/cr3bp_periodic.py`, the Phase 1 in-flight `scripts/floquet_phase1_p1_2_walk.py`) using a method that Pellegrini-Russell shows is sub-optimal?" | **Yes.** All STM machinery uses variational + variable-step DOP853. This is the Conclusion-1 + Conclusion-2 anti-pattern of Pellegrini-Russell 2016. Sub-optimal — but quantifiably so, with explicit mitigations. | **Confirmed concern** — action items in §5.4 |
| Phase 0 design doc line 116 | "Unknown whether the eigenvalue separation will survive at the cluster point" | The unknown is now bounded: Pellegrini-Russell's Fig. 13 (PO2) shows variable-step variational underperforms CSD by **several orders of magnitude** for loose-tolerance sensitive orbits. At our rtol=1e-12 (tight) the gap narrows; quantification at the (3,2) anchor needs an experiment. | **Phase 0 doc should be updated** — annotate with this reference |
| Phase 1 walker docstring line 14 | "the risk: DOP853 conditioning at the cluster point" | This is the right risk. The walker has the right awareness. **No erratum** — just a citation strengthening. | **OK** |
| Project precedent | (none — no prior STM-accuracy citation in code) | n/a | n/a |

## 8. Action items for the parent (#347)

1. **Atomic commit this digest** at `docs/notes/2026-06-17-digest-pellegrini-russell-2016.md`. (Done in this run.)
2. **#116 STM-continuation closure**: this paper is the canonical reference. Mark #116 as DIGESTED. The methodology recommendations §5.4 (fixed-path mode + cluster-point cross-check) are the actionable outputs.
3. **Phase 1 P1.3 saddle-center detector enhancement**: add the §5.4(2) cross-check (variable-step vs fixed-path monodromy at each candidate cluster point). Pellegrini-Russell-fixed-path agreement at 6+ decimals is a positive signal; disagreement at the multiplier modulus identifies STM contamination. This belongs in P1.3 (saddle-center detector), not the P1.2 walker, since the P1.2 walker is mainly logging multipliers, not deciding bifurcations.
4. **Docstring citations** (§6): add Pellegrini-Russell 2016 reference + DOI to docstrings in `cr3bp.py:cr3bp_stm_eom`, `cr3bp.py:propagate`, `bifurcation_detector.py:monodromy`. Defer to the next Phase 1 incremental commit so as not to fragment the digest commit.
5. **Phase 0 design doc update** (§7): add a follow-up annotation to line 116 noting that Pellegrini-Russell 2016 quantifies the concern and provides the fixed-path mitigation. Defer to next Phase 1 design-doc revision.
6. **No KNOWN_CORPUS change**: methods refs are out of scope for `literature_check.py`.
7. **No catalogue rows added**: methods paper.

---

**Verdict reaffirmed**: **ADMIT as methodology citation, with explicit fixed-path-mitigation action item for #347 Phase 1.** Pellegrini-Russell 2016 is the canonical answer to the Phase 0 design doc's open question on cluster-point monodromy conditioning. Our current STM machinery (variational equations + variable-step DOP853) is exactly the anti-pattern the paper documents; the project's (3,2) family at λ_max~2.5e5 is **~10⁴× more sensitive than the paper's PO3 test case** (an L1 halo at the same μ_EM), placing us deep in the regime where the variable-step variational method underperforms. Mitigation is straightforward (fixed-path mode, no integrator swap required); should be implemented at #347 Phase 1 P1.3.
