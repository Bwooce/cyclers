# Digest: Woollands, Taheri & Junkins 2019, "Efficient Computation of Optimal Low Thrust Gravity Perturbed Orbit Transfers"

**Task**: #353 (low-thrust digest Agent C, paper 4 of 4).
**Read**: 2026-06-17, all 27 pages.

## Bibliographic

- **Authors**: Robyn Woollands (JPL + TAMU), Ehsan Taheri (TAMU), John L.
  Junkins (TAMU).
- **Venue**: The Journal of the Astronautical Sciences.
- **DOI**: 10.1007/s40295-019-00152-9.
- **Published online**: 14 May 2019.
- **Funding**: AFOSR (Staci Williams), AFRL (Alok Das); JPL portion under
  NASA Spontaneous R&T Grant (Fred Hadaegh).

## What the paper actually is

A computational-efficiency paper for indirect-method low-thrust fuel-optimal
orbit transfers around a single large body (planet or asteroid) with a
**high-fidelity (70×70) spherical-harmonic gravity model**. The paper does
not invent a new optimal-control formulation; it speeds up an existing one
by combining three known tools in a novel cocktail:

1. **Hyperbolic-tangent smoothing** on the throttle switching function (same
   as paper 1 Junkins-Taheri 2019 and paper 3 Singh 2021 JOTA; the canonical
   reference is paper 1 [ref 42-43 in this paper]).
2. **Method of Particular Solutions (MPS)** as the shooting-method TPBVP
   solver — iteratively updates initial costates from n = 7 neighbour
   particular-solution propagations rather than from STM-Jacobian inversion.
3. **Picard-Chebyshev (PC) numerical integration** as the path-approximation
   propagator inside MPS — converges all nodes of a Chebyshev-spaced grid
   simultaneously (vs step-by-step Gauss-Jackson / RK), and reuses the
   reference-trajectory gravity evaluation across the n + 1 particular
   solutions.

The bench problems are **four 9-11-hour Earth-bound LEO-class transfers**
in a 70×70 spherical-harmonic field, with the spacecraft of m₀ = 100 kg,
T_max = 2 N, I_sp = 3000 s. The fidelity comparison studies whether the
**costate** equations need J₂, J₃, J₄, or just two-body, and concludes
that two-body costates produce engineering-precision results (mass
differences at the milligram level, Table 2 p.20).

The headline efficiency result: PC-MPS uses ~1.2-1.4 million full
70×70 force-model evaluations to converge each of the 4 cases; MATLAB
`fsolve` + `ode45` uses ~6.4-10.2 million (Fig 15, p.20). So ~5-8x
speedup with equivalent accuracy.

## Method digest

- **Optimisation framework**: indirect, PMP-based, single-shooting TPBVP.
  Same indirect family as papers 1, 2, 3.
- **Coordinates**: MEE (Modified Equinoctial Elements) `[p, f, g, h, k, l]`
  + mass `m`. Costate `λ = [λ_p, λ_f, λ_g, λ_h, λ_k, λ_l, λ_m]`. 7-d
  initial-costate unknown vector.
- **Dynamics**: Eq 4-5, MEE dynamics with high-fidelity 70×70 spherical
  harmonic gravity in the **state** equations. Costate equations use a
  **reduced** force model (two-body + J₂ + J₃ + J₄ at most) — the paper's
  key efficiency insight.
- **PMP / Hamiltonian**: Eq 7-10. Switching function S = c‖B^T λ‖/m + λ_m
  − 1 (Eq 15). σ*(S) = (1+sign(S))/2 (Eq 16).
- **HTS smoothing**: Eq 18: σ*(S, ρ) = (1/2)[1 + tanh(S/ρ)]. Continuation
  from ρ = 0.5 down to ρ < 1e-5 (deeper than paper 1's 1e-2). Fig 1
  shows the switching-function shape sharpen with ρ.
- **MPS shooter** (Sec "Method of Particular Solutions", p.10-12):
  - Reference trajectory x_ref(t), λ_ref(t) propagated from initial guess.
  - n = 7 particular solutions x_j(t), λ_j(t) with perturbed initial
    costates λ_ref(t₀) + Δλ_j (Eq 29).
  - Linear combination αⱼ solved (Eq 33) from the 7×7 system that drives
    Δx(t_f) = x_f − x_ref(t_f).
  - New initial costate estimate λ_new(t₀) = λ_ref(t₀) + Σ αⱼ Δλ_j(t₀)
    (Eq 36). Iterate.
  - Crucial detail: each particular solution uses the SAME Chebyshev
    grid as the reference, and the gravity gradient is "constant to >9
    digits within 50 m of a converged node" (cite [81-83]), so a **local
    force approximation** (Eq 40) computes the perturbed acceleration as
    `a_full,ref + (a_low,particular − a_low,ref)` — i.e. compute the full
    70×70 force only ONCE on the reference and add the J₂-J₆ delta
    between the perturbed and reference. This is the heart of the speed
    win.
- **Picard-Chebyshev iteration** (Sec PC, p.8-9): rearrange ẋ = f(t, x)
  to its integral form x(t) = x(t₀) + ∫ f(τ, x(τ)) dτ (Eq 22). Picard
  sequence (Eq 23) converges geometrically inside `|t_f − t₀| < d` (large
  for Cartesian; "up to 3 LEO orbits" or ">10 in MEEs" per ref [67]).
  Integrand approximated by Chebyshev polynomials. Segment-break times
  placed at switching-function zeros to handle thrust discontinuity.
- **Costate fidelity table** (Table 2, p.20): for the 4 test cases,
  final-mass differences between two-body, +J₂, +J₂+J₃, +J₂+J₃+J₄
  costate models are sub-mg. Verdict: **two-body costates suffice**.

## Relevance to this project

### To #309 (Sims-Flanagan low-thrust substrate) — **modest, NOT direct port**

- This paper is **indirect** method (PMP/TPBVP/MPS/PC); #309 is **direct**
  (Sims-Flanagan + differential-evolution + SLSQP). The two solve the same
  *problem* but they don't share solver machinery.
- **What COULD partially port**: the **local force approximation** trick
  (paper's Eq 39-40 — compute full force once on a reference, then
  cheaply propagate perturbations as low-fidelity-plus-delta) is generic
  and applies to ANY ensemble propagation. Sims-Flanagan with thousands
  of differential-evolution candidate vectors per generation is exactly
  the workload where this trick saves CPU. If #309 ever spends serious
  wall-time on multi-revolution heliocentric SF legs with non-trivial
  perturbations (J₂ of Earth, asteroid mascons), this approximation could
  speed up the population evaluation phase.
- **Picard-Chebyshev itself** is also a candidate accelerator for SF coast
  arcs that span many revolutions, BUT SF coast arcs are short single-rev
  Kepler segments and current Kepler propagator is already cheap. PC is
  unlikely to help current substrate; it would help only an expanded
  variant with explicit perturbations.
- **Practical verdict**: catalogue the local-force-approximation idea for a
  future #309 efficiency pass; do not act now.

### To #347 (Floquet bifurcation framework Phase 1) — **off-scope**

- Paper is body-bound LEO-class with no CR3BP / multi-body content. No
  Floquet content. Off-scope for #347.

### Catalogue / KNOWN_CORPUS impact — **NONE**

- The 4 bench cases (Table 1, p.15) are unnamed Earth-orbit transfers
  parameterised by initial/final MEE. They are not anchored to any named
  mission, asteroid, or named cycler/periodic-orbit. No corpus value.
- No periodic orbit, no cycler IC, no resonance condition is reported.
- The 70×70 spherical harmonic model is well-known (Earth gravity); not
  novel.

### To wider project — **methodology archive value**

- This is the canonical reference if the cyclers project ever needs to
  add **near-asteroid** or **near-comet** low-thrust precursor optimization
  with mascon / spherical-harmonic gravity. The paper explicitly mentions
  small-asteroid rendezvous and orbital-debris peer-to-peer
  ("traveling-salesman") tours as motivating applications. Neither of
  these is in our current scope.

## KNOWN_CORPUS impact

**No recommended anchors.** Pure body-bound LEO methodology; no periodic
orbit, no cycler structure, no resonance.

## Catalogue impact

**None.** This paper admits zero catalogue rows across all four classes.

## Action items

1. **No code changes required.** Defer.
2. **Capture two efficiency ideas in the project's "future-references"
   ledger**:
   - Local force approximation (paper Eq 39-40) for ensemble propagation.
     Candidate for #309 V2 if mascon / J₂ perturbations are added.
   - Picard-Chebyshev path integration as an alternative to step-by-step
     RK / Gauss-Jackson for long propagations. Candidate for any future
     multi-rev heliocentric SF lane.
3. **Bibliographic logging**: file under "indirect method efficiency"
   alongside paper 1 (which is the HTS reference). Cite paper 4 as the
   canonical MPS-PC fusion reference.

## Verdict

**Pure methodology paper**, off-scope for both #309 (Sims-Flanagan direct
substrate is a different solver family) and #347 (Floquet bifurcation is
multi-body, this is body-bound LEO). **NO catalogue rows, NO KNOWN_CORPUS
anchors, NO near-term action.** Keep on file as the canonical reference for
"MPS + Picard-Chebyshev + HTS" indirect-method efficiency stack, in case the
project ever needs a near-asteroid / small-body low-thrust precursor lane.
