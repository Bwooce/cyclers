# Digest: Junkins & Taheri 2019, "Exploration of Alternative State Vector Choices for Low-Thrust Trajectory Optimization"

**Task**: #353 (low-thrust digest Agent C, paper 1 of 4).
**Read**: 2026-06-17, all 18 pages.

## Bibliographic

- **Authors**: John L. Junkins, Ehsan Taheri (Texas A&M).
- **Venue**: Journal of Guidance, Control, and Dynamics 42(1), 2019; article-in-advance pagination 1-18.
- **DOI**: 10.2514/1.G003686.
- **Received**: 7 March 2018; accepted 19 June 2018; published online 4 October 2018.
- **Preliminary version**: This paper is the journal extension of the conference Taheri-Junkins ideas in refs [17], [21], [22] — i.e. it consolidates prior bang-off-bang smoothing and co-state initialisation work into a single coordinate-choice study.

## What the paper actually is

A systematic 8-coordinate-set comparison study for indirect minimum-fuel
low-thrust optimal-control TPBVPs. The paper does NOT propose a single new
method; it asks **which state representation gives the most robust / fast
indirect TPBVP** for two canonical bench problems:

1. Heliocentric interplanetary transfer Earth → asteroid Dionysus (e = 0.542,
   i = 13.6 deg, multi-revolution, mostly-coast, near-bang-bang structure).
2. Earth-centred GTO → GEO orbit raising (8 revolutions in 6 days, high force-to-
   thrust ratio, very nonlinear).

The eight coordinate sets are:

- Four **hybrid sets** combining orbital angles (Ω, i, θ) or (α, δ, ψ) with
  radial/transversal/normal velocity components in a rotating osculating triad
  — explicitly novel to this paper (Sec II.A-D, Eqs 21-37).
- Modified equinoctial elements (MEE, Sec II.E, Eq 39).
- The h-e set (specific-angular-momentum + eccentricity-vector components, Sec
  II.F, Eq 41).
- Classical Cartesian (Sec II.G, Eq 43).
- Spherical (Sec II.H, Eq 45).

Each set's dynamics is cast as `ẋ = A + B u` and fed into a Pontryagin minimum-
principle (PMP) indirect formulation with hyperbolic-tangent smoothing of the
bang-off-bang throttle (Eq 61-62). The smoothing parameter ρ is reduced by
homotopy from ρ = 1 (smooth) down to ρ_min ≈ 0.0078 (near bang-off-bang).

## Method digest

- **Optimisation framework**: indirect, PMP-based, single-shooting TPBVP. The
  Hamiltonian is the standard `H = (T_max/c) σ + λ_x^T A + (T_max/m) λ_x^T B û
  σ − λ_m (T_max/c) σ` (Eq 54). PMP gives `û* = −B^T λ_x / ‖B^T λ_x‖` (Eq 56)
  and `σ*` from the sign of switching function `S = c‖B^T λ_x‖/m + λ_m − 1`
  (Eqs 58-60).
- **Smoothing**: hyperbolic-tangent (HTSL) on the throttle switching function,
  Eq 62: `σ* ≈ (1/2)[1 + tanh(S/ρ)]`. Continuation in ρ from 1 → 0.0078
  (interplanetary) or → 0.0343 (GTO-to-GEO).
- **Integrator**: MATLAB `ode45` with abs+rel tol 1e-10. No Picard-Chebyshev
  here — that's paper 4 (Woollands-Taheri-Junkins 2019 JAS).
- **Solvers**: MATLAB `fsolve` (gradient/Newton with finite-difference Jacobian)
  and a Method of Particular Solutions (MPS, ref [31] Miele-Iyer 1970) solver.
  Random costate initialisation η₀ ∈ [0, 0.1]^7, 50 simulations per set per
  problem for statistical robustness measurement.
- **Performance metric**: % converged out of 50 random costate seeds, final
  mass, position/velocity defects, number of iterations, per-simulation CPU
  time.

## Key findings (Tables 4, 9, 11)

**Interplanetary (Earth → Dionysus, m_f ≈ 2718 kg)**:
- All hybrid sets, MEE, h-e, spherical: 100 % convergence with fsolve.
- Cartesian: only 10 % convergence (43/50 nominally converged but to suboptimal
  local minima; only 5 hit the true optimum).
- MPS fails on hybrid set 4 and Cartesian.
- Hybrid set 1 fastest (2.8 s/sim); MEE next (2.9 s/sim).
- Propagation regularity: MEE = 110 ms, h-e = 120 ms beat hybrid sets
  (Table 6).

**GTO → GEO (m_f ≈ 94.15 kg, 8 revs, 6 days)**:
- This problem distinguishes the sets sharply.
- MEE, h-e: 70 % and 76 % convergence (fsolve).
- Hybrid sets 1-4: 10-24 % convergence.
- Cartesian: 0 % convergence (NC across all 50 seeds).
- Order of merit for high force-to-thrust ratio: h-e > MEE > hybrid sets >
  spherical > Cartesian.

**Conclusion order** (paper's Sec VI):
- Interplanetary: hybrid 1 > MEE > h-e > hybrid 4 > hybrid 2 > hybrid 3 >
  spherical > Cartesian.
- Geocentric: h-e > MEE > hybrid 1 > hybrid 4 > hybrid 2 > hybrid 3 >
  spherical > Cartesian.

## Relevance to this project

This is a **methodology** paper about coordinate choice for indirect-method
low-thrust optimal control. The cyclers project's #309 low-thrust substrate
(`src/cyclerfinder/search/lowthrust.py`,
`src/cyclerfinder/search/lowthrust_maintenance.py`,
`src/cyclerfinder/search/low_thrust_cycler_search.py` and the underlying
`cyclerfinder.core.sims_flanagan` leg model) uses **Sims-Flanagan
transcription** — a *direct* method with per-segment impulses — solved with
scipy `differential_evolution` (global) + SLSQP (local polish). That is a
fundamentally different family from Junkins-Taheri's indirect / PMP / TPBVP
approach.

**To #309 (Sims-Flanagan substrate)**: **Limited direct value**.
- This paper does not improve the SF transcription itself.
- The HTSL smoothing trick (Eq 62) is a *throttle* smoothing for indirect
  bang-off-bang structure; SF's per-segment ΔV bounds are already a piecewise-
  linear convex constraint in the direct method and don't need smoothing.
- Honest verdict: a future *indirect* lane (if we ever build one) would
  benefit from this paper's coordinate ranking. The current direct lane does
  not.

**Possible refinement to #309 propagation/state choice**:
- The MEE coordinate set IS relevant if the SF segments' coasting propagation
  ever switches from Cartesian Kepler to MEE for multi-revolution heliocentric
  transfers. MEE is 5-slow + 1-fast and well-behaved over many revolutions —
  the same regularity that lets MEE/h-e dominate for indirect would also reduce
  numerical drift across SF coast arcs that span >1 revolution. Worth noting
  for a future #309 efficiency pass.
- The h-e set is novel and lacks broad community adoption; ignore for now.

**To #347 (Floquet bifurcation Phase 1)**: **Off-scope**. Floquet bifurcation
of CR3BP periodic orbits has no immediate need for indirect low-thrust
optimal-control coordinate choice. No overlap.

**KNOWN_CORPUS anchors**: **None**. The two test problems are
- Earth → Dionysus (asteroid rendezvous, not a cycler).
- GTO → GEO (planet-bound transfer).
Neither produces a periodic / quasi-periodic orbit suitable as a corpus
anchor.

## KNOWN_CORPUS impact

**No recommended anchors.** This is pure methodology; no cycler IC, no Floquet
multipliers, no periodic-orbit constants. The Dionysus rendezvous final mass
2718.16 kg (Table 4) is interesting as a community benchmark but is irrelevant
to the cycler corpus.

## Catalogue impact

**None.** This paper admits zero catalogue rows in any class (cycler /
quasi_cycler / precursor_mga / mga_tour). The methodology is body-bound
rendezvous and orbit-raising, not transport between bodies in resonance.

## Action items

1. **No code changes required.** Defer to a hypothetical future indirect-method
   lane.
2. **Coordinate-set note for #309 V2 (if a propagation-regularity pass is ever
   scheduled)**: the MEE result is the cleanest take-away — if SF segments grow
   long enough that Cartesian propagation drift becomes a defect-budget
   problem, MEE is the canonical fix. Capture in a future plan, not now.
3. **Bibliographic logging**: mark the paper as **read & off-scope for the
   active substrates**. Cite only if a future indirect lane is opened.

## Verdict

**Pure methodology paper, off-scope for both #309 (Sims-Flanagan direct
substrate) and #347 (Floquet bifurcation). No catalogue rows, no KNOWN_CORPUS
anchors, no near-term action.** Keep on file as the canonical reference for
indirect-method coordinate choice in case the project ever opens an indirect
low-thrust lane.
