# Digest: Singh, Anderson, Taheri & Junkins 2021, "Low-Thrust Transfers to Southern L2 Near-Rectilinear Halo Orbits Facilitated by Invariant Manifolds"

**Task**: #353 (low-thrust digest Agent C, paper 3 of 4).
**Read**: 2026-06-17, all 28 pages.

## Bibliographic

- **Authors**: Sandeep K. Singh (TAMU), Brian D. Anderson (JPL), Ehsan Taheri
  (Auburn), John L. Junkins (TAMU). Same author quartet as paper 2.
- **Venue**: Journal of Optimization Theory and Applications (JOTA), Springer.
- **DOI**: 10.1007/s10957-021-01898-9.
- **Received**: 12 September 2020; accepted 16 June 2021; online 1 July 2021.
- **Communicated by**: Mauro Pontani.
- **Preliminary version**: AAS 20-565 at the 2020 AAS/AIAA Astrodynamics
  Specialist Virtual Lake Tahoe Conference.

## What the paper actually is

A piecewise minimum-time and minimum-fuel low-thrust transfer study connecting
a **super-GTO (sGTO)** parking orbit with **three candidate southern L2
Near-Rectilinear Halo Orbits (NRHOs)** in the Earth-Moon system. The three
NRHOs are picked at synodic resonances 9:2 (the Lunar Gateway operational
candidate), 24:5, and 4:1. Stable-manifold legs serve as the long terminal
coast arcs for **NRHO-bound** transfers; unstable-manifold legs serve the same
role for **Earth-bound** return. The paper compares:

1. **Direct** transfers sGTO → NRHO (no manifold leverage).
2. **Manifold-coupled** transfers: sGTO → stable-manifold IC → coast on the
   manifold → NRHO insertion.

For both minimum-time and minimum-fuel cost functionals, in **both** CR3BP
(IC computation) and **High-Fidelity Model (HFM)** point-mass N-body Earth +
Moon + Sun + Jupiter with JPL DE436 ephemerides (manifold propagation).

The novelty over Singh et al. 2021 Acta Astronautica (paper 2) is:
- NRHO target instead of L1 halo + Lunar Polar Orbit.
- HFM ephemeris point-mass perturbation analysis instead of BCP solar-only.
- Explicit minimum-fuel HTC bang-off-bang solution with switching-function
  diagnostics (Figs 13, 16).
- Three NRHOs compared (table-style trade study), not one.

## Method digest

- **Optimisation framework**: indirect, PMP-based, single-shooting TPBVP.
  MATLAB `fsolve` (Levenberg-Marquardt), `OptimalityTolerance` and
  `StepTolerance` 1e-12. Homotopy with ~20 (min-time) and ~30 (min-fuel)
  steps. State-costate dynamics propagated with mex-compiled `ode45`.
- **Coordinates**: MEE (Modified Equinoctial Elements) `X = [P, e_x, e_y,
  h_x, h_y, L]`. Eq 3 gives `Ẋ = b + (T_max/(m_0 m g_0)) M δ_T σ`,
  `ṁ = −(T_max/(m_0 c)) σ`. Eq 4 gives the M matrix. Pure 6-MEE +
  mass state.
- **PMP**: Hamiltonian H_MT = Γ^T Ẋ (Eq 6). Optimal thrust direction
  δ_T* = −M^T Λ / ‖M^T Λ‖ (Eq 7) — same primer-vector form as paper 1.
- **Min-fuel cost**: J_MF = m_p = (T_max/(m_0 c)) ∫ σ dt (Eq 8). Switching
  function SF = (c A^T M δ_T*)/m + λ_m − 1 (Eq 13). σ*(SF) = 1 if SF > 0;
  0 if SF < 0 (Eq 14).
- **Smoothing**: hyperbolic-tangent (HTS), Eq 15-16:
  `σ*(SF, ρ) = (1/2)(1 + tanh(SF/ρ))`. Same recipe as paper 1
  (Junkins-Taheri 2019); paper 1 is the canonical HTS reference here
  (cited as ref [33] in this paper).
- **HFM dynamics**: Eq 2, point-mass N-body relative-to-Earth with Moon, Sun,
  Jupiter as perturbers, using JPL DE436 (ref [11]).
- **Spacecraft parameters**: m_0 = 1000 kg, I_sp = 1500 s, T_max = 1 N
  (much larger thruster than paper 2's 0.096 N).

## **Sourced NRHO ICs (Table 1, p.5) — KNOWN_CORPUS candidates**

Earth-Moon CR3BP, EM synodic, μ = 0.012150586632602, DU = 386 274.56245094 km,
TU = 375 699.8074372613 s:

| Case | Synodic resonance | Period (days) | Perilune (km) | Jacobi C | Stability idx 1 | Stability idx 2 |
|------|------------------|---------------|---------------|----------|-----------------|-----------------|
| A    | 9:2 (Gateway)    | 6.573         | 3269          | 3.0446   | −1.3753         | 0.6626          |
| B    | 24:5             | 6.074         | 2047          | 3.0540   | −1.1119         | 0.7451          |
| C    | 4:1              | 7.322         | 5542          | 3.0332   | −1.6434         | 0.4837          |

These are "mildly unstable" NRHOs — stability indices nearly within
the |λ| ≤ 1 box per Howell's classification (ref [9]). The paper notes
manifolds exist for all three (necessary condition for transport).

The paper's Table 2 (p.10) and Table 3 (p.11) further give **15-digit
synodic-frame Cartesian IC/FC states** for the selected stable manifold IC and
unstable manifold FC of each case — directly usable as anchors. Examples
(Case A, 9:2 NRHO stable manifold IC):

`(x, y, z) = (−0.01215, −0.56415, 0.05927) DU`
`(v_x, v_y, v_z) = (0.81277, 0.05946, 0.33727) VU`
`TOF = 127.5501 days` (on manifold, CR3BP).

Case A unstable manifold FC (Earth side):
`(x, y, z) = (−0.01215, 0.56415, 0.05926) DU`
`(v_x, v_y, v_z) = (−0.81277, 0.05945, −0.33727) VU`
`TOF = 127.5505 days`. Mirror-image of stable (as expected by CR3BP
y-symmetry).

## Key end-to-end results

**Table 7 (p.18) Min-time transfers** (sGTO → NRHO):

| Param         | 9:2 D | 9:2 M | 24:5 D | 24:5 M | 4:1 D | 4:1 M |
|---------------|-------|-------|--------|--------|-------|-------|
| m_f (kg)      | 829.49| 876.63| 829.85 | 879.97 | 828.89| 862.58|
| Total TOF (d) | 29.03 | 78.73 | 28.96  | 77.85  | 29.13 | 82.25 |
| Coast time (d)| —     | 57.73 | —      | 57.42  | —     | 58.86 |
| Thrust time (d)| 29.03| 21.00 | 28.96  | 20.43  | 29.13 | 23.39 |
| Total ΔV (m/s)|2749.5 |1936.5 |2743.1  |1880.6  |2760.1 |2174.1 |

D = direct, M = manifold. Manifold saves 586-860 m/s of ΔV in trade for
49-58 extra days.

**Table 8 (p.21) Min-fuel transfers** (sGTO → NRHO via manifold):

| Param        | 9:2 MF | 24:5 MF | 4:1 MF |
|--------------|--------|---------|--------|
| m_f (kg)     | 900.51 | 902.65  | 889.59 |
| Total TOF (d)| 93.20  | 92.89   | 94.33  |
| Coast (d)    | 75.18  | 74.58   | 75.84  |
| Thrust (d)   | 18.02  | 18.31   | 18.50  |
| Total ΔV (m/s)| 1540.8| 1506.4  | 1720.7 |

Min-fuel further saves 300-350 m/s over min-time-via-manifold, at the cost
of ~15 days.

**Conclusion**: manifold-aided min-fuel saves ~1040-1236 m/s ΔV vs direct
min-time, at the cost of ~64-65 extra days. For uncrewed cargo/resupply,
this is the recommended mode.

## Relevance to this project

### To #309 (Sims-Flanagan low-thrust substrate) — **modest, not direct**

- Same caveats as paper 2: this is **indirect-method** TPBVP, not Sims-Flanagan
  direct. The HTS smoothing trick does not port to SF's already-piecewise-
  linear thrust bounds.
- However, the paper's clear demonstration that **manifold-coupled coast arcs
  save ~30-45 % ΔV vs direct low-thrust transfers** is a quantitative anchor
  the cyclers project should respect. Any future powered-cycler-precursor
  optimisation that has access to invariant manifolds of a CR3BP orbit
  (e.g. a powered Mars-cycler precursor with EM L1 NRHO drop-off) should
  budget for this saving.
- **Recipe to potentially port**: the workflow is
  1. Pre-compute manifold IC/FC on a periodic orbit (CR3BP).
  2. Sort piercing points on an `e-i-r` phase plot, pick favourable cluster.
  3. Re-validate in HFM with date-sliding (epoch as free variable). Select
     best epoch by minimum-radius criterion.
  4. Solve TPBVP for the body-bound piece only.

  This recipe is **the standard procedure** for low-thrust manifold transport
  and is essentially what #309 should mimic if it ever adds a "deliver via L1
  halo manifold" precursor task.

### To #347 (Floquet bifurcation framework Phase 1) — **MODERATE / VALUABLE**

- Table 1 gives **three sourced NRHO state-vector ICs** with their **stability
  indices** computed. This is an immediate Floquet cross-check target.
- Stability index 1 = (λ₁ + 1/λ₁)/2, stability index 2 = (λ₃ + 1/λ₃)/2 in
  standard Howell convention. The paper's reported values are dimensional
  ground truth for #347's monodromy code.
- Specifically: for Case A (9:2 NRHO), stability indices are −1.3753 and
  0.6626. This means λ₁ + 1/λ₁ = −2.7506 (negative ⇒ real-pair with sign
  flip, i.e. period-doubling neighbourhood) and λ₃ + 1/λ₃ = 1.3252.
- A period-doubling Floquet bifurcation is *exactly* the RTR2026 saddle-
  center mechanism the #347 substrate is meant to detect. **These three
  NRHO ICs and their stability data are perfect spec-validation anchors for
  Phase 1's bifurcation classifier.**
- Singh et al. do not actually run a bifurcation analysis — they just report
  the indices and use the manifolds — but the indices indicate a near-
  period-doubling regime that #347's substrate should classify correctly
  when fed these ICs.

### Catalogue / KNOWN_CORPUS impact — **VALUABLE**

Three sourced NRHO ICs ready for KNOWN_CORPUS addition as
**earth_moon_libration** anchors:

1. Southern L2 NRHO at 9:2 synodic resonance, T = 6.573 days, C = 3.0446,
   stability (−1.3753, 0.6626) — Gateway target.
2. Southern L2 NRHO at 24:5, T = 6.074 days, C = 3.0540,
   stability (−1.1119, 0.7451).
3. Southern L2 NRHO at 4:1, T = 7.322 days, C = 3.0332,
   stability (−1.6434, 0.4837).

Plus the **manifold IC/FC** Tables 2 and 3 (six rows total) as precursor_mga
transit-state anchors if/when that schema admits "manifold piercing points"
as anchored objects.

## Catalogue impact

This paper does **not directly admit catalogue rows** in any class
(cycler / quasi_cycler / precursor_mga / mga_tour). The sGTO → NRHO transfer
is a one-shot delivery, not a cycler. However:

- The three NRHOs themselves could become **precursor_mga** rows once the
  cross-system framework (#316) admits "EM L2 NRHO → ..." chains.
- The 9:2 NRHO is the operational Lunar Gateway orbit and is a high-
  visibility anchor for the catalogue.

## Action items

1. **High priority: stage three NRHO KNOWN_CORPUS anchors** (9:2, 24:5, 4:1
   southern L2) for a future literature_check.py refresh. Cite this paper
   (DOI 10.1007/s10957-021-01898-9). Recommendation only; do not edit corpus
   in this task.
2. **Spec-validation cross-check for #347 Phase 1**: feed the three NRHO ICs
   (after converting period to dimensionless and re-deriving the IC vector
   from the perilune state via standard NRHO-correction) into the Floquet
   monodromy code; expect to recover stability indices to ~3-4 digits.
   Becomes a sourced regression test.
3. **Note for #309 V2 powered-precursor mode**: budget for ~40 % ΔV
   improvement when manifold coupling is allowed. This is the ground truth
   number for "is the manifold worth the wait" trade studies.
4. **Bibliographic logging**: pair with paper 2 (same author quartet,
   different target orbit and dynamics fidelity). Cite both as the canonical
   **manifold-coupled low-thrust** references for EM cislunar transport.

## Verdict

**Methodology + sourced ICs paper.** Off-scope for #309 SF substrate
(indirect, not direct). **VALUABLE for #347 Floquet Phase 1**: three sourced
NRHO ICs with stability indices — immediate spec-validation anchors. **HIGH
VALUE for KNOWN_CORPUS**: three southern L2 NRHO anchors (incl. the Gateway
9:2) plus six manifold IC/FC piercing states. Recommend three primary anchor
additions in a future curator pass.
