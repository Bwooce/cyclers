# Digest — D'Amario, Byrnes & Stanford 1981 "A New Method for Optimizing Multiple-Flyby Trajectories" (AIAA 80-1676R)

**Date**: 2026-06-19 AET (Agent dispatched by parent #384; corpus-document-policy digest — METHODOLOGY PROVENANCE, lighter digest)
**Verdict (TL;DR)**: **METHODOLOGY PROVENANCE for the project's corrector/optimizer lane.** This is the canonical JPL multiple-flyby trajectory optimizer (the "MOSES"/MULCON lineage): a **bounds-constrained Newton parameter-optimization** using **analytic first and second derivatives** (gradient + Hessian) of a **weighted sum-of-squares-of-ΔV cost function**, with **multiconic propagation + analytic state-transition matrices** for each trajectory segment. NOT a cycler paper — no catalogue admission. Directly relevant as a citable, sourced reference for the project's BVP corrector (#380) and STM work (#347): it documents (a) why sum-of-squares (not sum-of-magnitudes) gives a smooth gradient and a Newton-amenable Hessian, (b) the restart-reweighting trick to recover the sum-of-magnitudes optimum, and (c) the analytic STM-based partial-derivative chain for multi-segment flyby trajectories.

---

## 1. Header

- **Title** (verbatim, p.591): "A New Method for Optimizing Multiple-Flyby Trajectories"
- **Authors** (verbatim, p.591 byline + footnotes): **L.A. D'Amario** (Member Technical Staff, Member AIAA), **D.V. Byrnes** (Consultant, Member AIAA), **R.H. Stanford** (Member Technical Staff, Member AIAA)
- **Affiliation**: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, Calif.
- **Journal**: J. Guidance and Control, Vol. 4, No. 6, Nov.-Dec. 1981, pp. 591-596. DOI 10.2514/3.56115
- **Originally**: AIAA Paper 80-1676 / AIAA-80-1676R, AIAA/AAS Astrodynamics Conference, Danvers Mass., Aug 11-13, 1980 (submitted Sept 16, 1980) (p.591 footnote)
- **cyclers_pdf path**: `papers/damario-byrnes-stanford-1981-new-method-optimizing-multiple-flyby-trajectories-jgc-doi-10.2514-3.56115.pdf`
- **Page count**: 6 pages (591-596), text layer present

## 2. What the paper is

A new procedure to **minimize total impulsive ΔV for multiple-flyby trajectories** subject to constraints on flyby parameters and maneuver times (Abstract, p.591). It supersedes an earlier penalty-function quasi-Newton method (ref 7, D'Amario-Byrnes-Sackett-Stanford 1979) that was "exceedingly slow" with terminal-convergence and segment-fitting problems (p.591). Demonstrated by optimizing **Galileo satellite tours containing up to eleven flybys** (Abstract).

## 3. The method (sourced, with page citations)

### 3.1 Cost function (Section II.A, p.592)
- **F = Σ f_i · ΔV_i^T ΔV_i** (Eq. 1) — **weighted sum of the SQUARES of the ΔV magnitudes**, in contrast to the sum-of-magnitudes of ref 7.
- Advantages stated (p.591-592): (1) eliminates the gradient discontinuity that the sum-of-magnitudes form has whenever any ΔV → 0; (2) a good Hessian approximation can be built from **first derivatives of the velocity vectors only** (the second-derivative term in Eq. 3 is much smaller than the first, p.592); ⇒ enables a **Newton algorithm** with quadratic convergence ("spectacular improvement in terminal convergence and computing time," p.591).
- Gradient (Eq. 2): ḡ^T = ∂F/∂X̄ = 2 Σ f_i ΔV_i^T ∂ΔV_i/∂X̄. Hessian (Eq. 3): two terms, first dominant.
- Newton step (Eq. 8): X̄ = X̄_0 − H_0^{-1} ḡ_0, used as a trial step then refined by a 1-D line search (p.593).

### 3.2 Independent variables (Section II, p.592; Eqs 4-5)
X̄ = stacked vector of, for each flyby i: **altitude h_i, B-plane angle θ_i, periapse time T_{p,i}** (vector P̄_i), plus the **maneuver times T_j**. The B-plane angle θ defines the orientation of the approach hyperbola (Fig 2, B-plane coordinate system, p.592).

### 3.3 Constraint handling (p.591-593)
- Any independent variable may be **bounded by simple range constraints** (upper/lower; equality when upper=lower). This **eliminates penalty functions** — the key practical win over ref 7. During the line search a variable that hits a limit is fixed; freed later if a lower minimum is indicated (p.593).
- The bounded Newton uses an NPL (National Physical Laboratory, England) algorithm based on Gill & Murray Ch. 2 (ref 8) (p.593).

### 3.4 Weight-factor restart-reweighting (Section II.B, p.593)
- Setting f_i = 1/|ΔV_i*| makes the weighted sum-of-squares solution equal the sum-of-magnitudes solution (Eq. 9). In practice: solve, then **restart with new weights computed from the just-found ΔV magnitudes**; "with several restarts of a few iterations each, the desired sum-of-magnitudes solution is easily found" (p.593). This is the mechanism that gets the true minimum-total-ΔV (sum-of-magnitudes) answer while keeping the smooth Newton-friendly sum-of-squares machinery.

### 3.5 Trajectory generation + STM (Section III, p.593-594)
- Each segment generated in two steps (Fig 3, p.593): (1) target a trajectory from the preceding maneuver point to the upcoming flyby's periapse conditions (h, θ, T_p) by varying the 3 components of initial velocity V̄_0 via an **ordinary Newton iteration** ΔV̄_0 = −K^{-1} ΔP̄ (Eqs 10-12), where K = (∂P̄/∂X̄_p)(∂X̄_p/∂V̄_0) uses the **state transition matrix** for the segment; (2) propagate to the next maneuver time.
- **Two propagation models** (p.593-594):
  - **1STEP** — three-body dynamic model (central body + flyby satellite + spacecraft), one-step multiconic via Byrnes' pseudostate theory (ref 9) built on Wilson (ref 10). Eliminates 80-90% of simple-conic error; fast/cheap; good for preliminary design + tradeoffs.
  - **MULCON** — multistep multiconic (central body, arbitrary number of satellites, solar perturbations, central-body oblateness), based on Kwok-Nacozy (ref 11). At least an order of magnitude slower than 1STEP, used when more precise dynamics needed.
- **Multiconic STM is computed ANALYTICALLY without significant additional cost** — the property that makes multiconic "especially well-suited" here, since the STM drives every segment's targeting and the gradient/Hessian (p.593-594).

### 3.6 Partial derivatives (Appendix, p.596)
The ∂ΔV̄_i/∂P̄_j and ∂ΔV̄_i/∂T_j chains (Eqs A1-A6) are cascading: ΔV̄_i depends on every independent variable preceding it plus the periapse parameters immediately following. Built from segment STMs Φ(T_i, T_{i-1}) (Eq A4) and local Jacobians W_i, U_i (Eqs A5-A6).

## 4. Results (Section IV, p.594-595) — the worked Galileo 10-flyby tour
- A 10-flyby Galileo satellite tour (1 Callisto, 5 Ganymede, 4 Europa; 15-month duration), optimized in four stages (sum-of-squares → reweighted sum-of-magnitudes → maneuver-time activation → flyby-variable reactivation) (Table 2-4, p.594-595).
- 1STEP total ΔV reduced from conic 74 m/s → 32.3 m/s across the 4 stages (Table 2). MULCON total 44.0 m/s (Table 5) — 11.7 m/s higher than 1STEP, primarily from solar-perturbation modeling absent in 1STEP and from the suboptimal 4-flyby-segmenting forced by core-memory limits (p.595).
- 1STEP vs MULCON altitudes agree to within ~1 km, B-plane angles within ~0.1° (p.596); against a precision integrated trajectory the MULCON ΔV's agree to within a few tenths of a m/s (p.596).
- Physical insight (p.595-596): most maneuvers become ~zero at the optimum (ballistic where possible); residual nonzero ΔV's are tied to a flyby parameter on a limit; maneuvers are most efficiently placed near apoapse where velocity is smallest.

## 5. References of note (p.596)
- Ref 7: D'Amario, Byrnes, Sackett, Stanford, "Optimization of Multiple Flyby Trajectories," AAS 79-162, 1979 — the **superseded** penalty-function predecessor.
- Ref 8: Gill & Murray, *Numerical Methods for Constrained Optimization*, Academic Press 1974 — the NPL Newton basis (same ref appears in the Galileo 1983 paper).
- Ref 9: Byrnes, "Applications of the Pseudostate Theory to the Three-Body Lambert Problem," AAS 79-163, 1979.
- Ref 11: Kwok & Nacozy, MULCON Final Report / User's Guide, Univ. of Texas at Austin.
- Ref 14: D'Amario, "Minimum Impulse Three Body Trajectories," PhD diss., MIT, 1973.

## 6. Catalogue / KNOWN_CORPUS relevance

### 6.1 Catalogue admission — NONE
Not a cycler, not a trajectory family. It is a numerical method. No row.

### 6.2 Corrector/optimizer lane provenance (the value — #380 BVP corrector, #347 STM)
This is a **sourced, citable foundation** for several choices the project either makes or could make:
- **Sum-of-squares-of-ΔV cost** with restart-reweighting to recover sum-of-magnitudes (§3.1, §3.4): a smooth, Newton-friendly objective that avoids the |ΔV|→0 gradient kink. If the project's optimizer lane ever needs to minimize total impulsive ΔV across a multi-flyby/cycler establishment chain, this is the canonical justification for the objective form.
- **Bounds-constrained variables instead of penalty functions** (§3.3): altitude, B-plane angle, periapse time, maneuver time — all box-constrained, handled inside the Newton solver. Relevant to how the corrector's free variables / constraints should be posed.
- **Analytic STM-driven targeting + gradient/Hessian** (§3.5, §3.6): the multiconic STM is computed analytically per segment and reused for both targeting and the optimization derivatives. This is the same STM-centric structure the project's #347 STM work and #380 BVP corrector rest on; the paper is a precedent that the STM chain across multiple flyby segments is analytically tractable.
- **1STEP vs MULCON two-fidelity strategy** (§3.5): fast pseudostate-multiconic for preliminary design, precise multistep-multiconic for final — a fidelity-ladder pattern directly analogous to the project's circular-coplanar → analytic-ephemeris → CR3BP/integrated ladder.

### 6.3 KNOWN_CORPUS
Methodology paper; not a cycler-novelty datapoint. Register as the corrector/optimizer-lane provenance source. It is the *house optimizer* behind the D'Amario-Byrnes Galileo design work (the PLATO program in AIAA-83-0099 is the same Newton + analytic-derivative + bounds-constrained family), tying the two 1983 Galileo papers and this 1981 method into one JPL design-and-optimize lineage.
