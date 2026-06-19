# Gurfil (ed.) 2007 — Modern Astrodynamics (Elsevier Astrodynamics Series, vol. 1)

**Deep-read verdict note, 2026-06-19 AET.** Multi-author edited volume,
textbook-summary scope. TOC + foreword + editor's introduction read in
full; deep-read of the four project-relevant chapters
(Ch3 Guibout-Scheeres BVP, Ch4 Belbruno low-energy, Ch5 Junge-Dellnitz
set-oriented, Ch6 Ross optimization); section-header sampling of the
other four (Ch1 Vallado, Ch2 Efroimsky, Ch7 McInnes-Cartmell, Ch8
formation flying). Text layer present (machine-readable PDF), so the
size discipline was easy here — extraction via pdftotext -layout, no
page-image OCR needed.

## Header

- **Title:** *Modern Astrodynamics* (Elsevier Astrodynamics Series,
  Volume 1)
- **Editor:** Pini Gurfil (Technion — Israel Institute of Technology)
- **Series foreword:** David A. Vallado (Analytical Graphics Inc.)
- **Publisher:** Elsevier / Butterworth-Heinemann (Academic Press
  imprint)
- **Year:** 2007 (1st edition)
- **Pages:** 293 (PDF); 8 chapters
- **Format:** Invited-contribution research monograph; each chapter is
  a self-contained survey by a leading researcher, with its own
  reference list.
- **PDF:** machine-readable text layer present (unlike Szebehely 1967).

## Chapter map (PDF page = chapter start; running headers confirm
numbering)

| Ch | Author(s) | Title | PDF p | Project relevance |
|----|-----------|-------|-------|-------------------|
| 1 | David Vallado | Perturbed Motion | 7 | low (Earth-satellite force models) |
| 2 | Michael Efroimsky | Gauge Freedom in Astrodynamics | 29 | low-moderate (osculating-element subtlety) |
| 3 | V. Guibout & D. Scheeres | Two-Point Boundary Value Problems | 59 | **HIGH — #380** |
| 4 | Edward Belbruno | Low-Energy Transfers and Applications | 108 | moderate — #347 P3 (subsumed by #376) |
| 5 | O. Junge & M. Dellnitz | Set Oriented Numerical Methods | 132 | **HIGH — #314** |
| 6 | I. Michael Ross | Space Trajectory Optimization | 163 | moderate — optimization framing |
| 7 | C. McInnes & M. Cartmell | Orbital Mechanics of Propellantless Propulsion (solar sails + tethers) | 199 | low |
| 8 | (formation-flying authors) | Cooperative Spacecraft Formation: Open/Closed-Loop Robustness | 248 | low |

**Scope-task correction (honest):** the task brief named a "Scheeres —
strongly-perturbed orbital mechanics" chapter. There is **no
standalone Scheeres chapter** in this volume. Scheeres appears only as
**co-author of Chapter 3** (with Guibout, on two-point BVPs via
generating functions). The brief's "Scheeres" target therefore folds
into the Chapter 3 deep-read below. Chapter 7 ("Orbital Mechanics of
Propellantless Propulsion") is McInnes-Cartmell on solar sails and
tethers, not Scheeres perturbation theory.

---

## Chapter 3 (Guibout & Scheeres) — Two-Point Boundary Value Problems
### → directly feeds #380 (BVP integral-constraint corrector)

This is the highest-value chapter for the project. It develops a
**generating-function / Hamilton-Jacobi method for solving two-point
boundary value problems (2PBVPs)** that is an alternative to
shooting/differential-correction.

### Core theory (Sections 3.2.1-3.2.2, pp 62-71)

- A canonical transformation between two phase-space coordinate sets is
  characterized by a **generating function** F. The Legendre
  transformation (eq 3.14) converts one generating-function type into
  another (F1, F2, ... families indexed by which variables are held).
- The key object is the transformation **φ induced by the phase flow**
  (Prop 3.2.4, p 66): φ maps the current state to the initial state.
  Its inverse maps the system to equilibrium. The generating functions
  of φ satisfy the **Hamilton-Jacobi equation** (eq 3.27, p 70).
- **Lemma 3.2.6 (p 69): "Generating functions solve two-point boundary
  value problems."** This is the load-bearing result. Once you have the
  generating function F1(q, q0, t) for the phase flow, the relations
  p = ∂F1/∂q and p0 = -∂F1/∂q0 (eqs 3.28-3.29, p 71) let you read off
  the momenta from a *prescribed pair of positions* (q0 at t0, q at
  t1) — i.e., they solve the position-to-position BVP directly, no
  iteration.

### Caustics = solution non-uniqueness (Section 3.2.4, pp 80-90)

- **Prop 3.2.10 / Thm 3.2.12-3.2.13:** a generating function is
  *singular* exactly when the local projection of the Lagrangian
  submanifold fails to be a diffeomorphism — geometrically a
  **caustic**. A caustic corresponds to **multiple solutions of the
  BVP** (the number k of solutions is the multiplicity). Prop 3.2.3
  guarantees at least one well-defined generating function exists at
  every instant, so when one type becomes singular you switch to
  another via the Legendre transform.
- This is the rigorous analogue of "Lambert's problem can have multiple
  revs / multiple branches" — singularities of F enumerate the
  multiple targeting solutions.

### The computational algorithm (Section 3.4, pp 91-99 + Appendix A)

- For **relative motion about a reference trajectory** (worked example:
  relative to the L2 point of a three-body problem), expand the
  Hamiltonian Hh as a **Taylor series of order N** in the relative
  state. Seek the generating function as a polynomial of order N.
  Substituting into the Hamilton-Jacobi equation reduces it to a set of
  **ODEs for the Taylor coefficients** of F, integrated alongside the
  reference trajectory. Appendix A gives the order-3 expansion
  explicitly.
- Caveat stated by the authors (p 93): a finite-order Hamiltonian does
  NOT have a finite-order generating function — truncating F at order N
  is always an approximation, and one must verify the Taylor series
  converges (radius-of-convergence warning, p 95).

### #380 recommendation (BVP integral-constraint corrector)

The Guibout-Scheeres method is a **constraint-respecting, derivative-
exact 2PBVP solver** that differs from the project's current
shooting/Newton correctors in three useful ways:

1. **It solves the position-to-position BVP without an initial guess
   for the velocity** — the generating function delivers p and p0
   analytically once F is computed. For #380's integral-constraint
   corrector this means the boundary constraints (endpoint positions /
   periodic-closure positions) enter as data, not as residuals to be
   nulled by iteration.
2. **Caustic detection = multiplicity enumeration.** Where #380 might
   silently converge to one of several closure solutions, the
   singularity structure of F tells you *how many* solutions exist and
   *where* they merge. This is a principled diagnostic for the
   "it closed — but is this the family member I wanted?" failure mode
   (cf. memory: orbit-closure discipline).
3. **The relative-motion Taylor-coefficient ODE machinery** maps onto
   the project's existing variational/STM propagation: the order-2
   coefficients are essentially the STM, higher orders are the
   state-transition tensors. #380 could lift the order-N generating
   function from the same integration that already produces the STM.

Cost note: implementing the order-N generating-function corrector is a
**multi-week build** (Taylor-tensor bookkeeping + Legendre-switch logic
+ convergence guards). The order-2/order-3 version is tractable and
maps directly onto existing STM infrastructure. Cite Guibout PhD thesis
(Michigan 2004, ref [12]) and Guibout-Scheeres JGCD 27(4):693-704
(2003, ref [16]) for the full algorithm.

---

## Chapter 5 (Junge & Dellnitz) — Set Oriented Numerical Methods
### → directly feeds #314 (heteroclinic mass-transport framework)

This chapter is the **set-oriented (box-covering) computational
machinery** for invariant manifolds and connecting orbits — the GAIO
approach. It is the natural numerical backbone for #314's transport
framework.

### The subdivision / box-covering core (Section 5.3, pp 135-139)

- Phase space is covered by **boxes (rectangles)** organized in a tree;
  a **multilevel subdivision algorithm** (5.3.1) repeatedly bisects
  boxes and keeps only those whose image under the map intersects the
  current collection. The sequence of coverings converges to the
  **relative global attractor** / all invariant sets in the region Q,
  *together with their unstable manifolds* (worked on the Hénon map,
  Figs 5.3-5.4). The reference software is **GAIO** (MATLAB; a script
  is printed in the chapter).
- Key property (vs direct simulation, p 139): subdivision covers ALL
  invariant sets and their unstable manifolds in Q, including ones a
  forward simulation would never land on.

### Invariant-manifold continuation (Section 5.4, pp 140-143)

- A **continuation algorithm** (5.4.1) grows an initial box covering of
  the *local* unstable manifold of a hyperbolic fixed point into the
  *global* unstable manifold (Prop 5.4.1, Hausdorff convergence).
- **Prop 5.4.2 + the modified continuation step:** naive continuation
  blows up error under weak transversal contraction; the fix is to
  continue with a single application of a **time-T map** per step
  (perform one continuation step while computing the time-T image),
  bounding the covering size.
- **Worked application: the NASA Genesis halo orbit** — they compute the
  global unstable manifold of an unstable halo orbit near L1 in the
  Sun-Earth CR3BP (μ = 3.040423398444176e-6), Fig 5.5. (A flythrough
  movie was hosted at the Paderborn group's site.)

### Connecting orbits + controlled transport (Sections 5.5-5.6, pp 144-156)

- **Connecting (heteroclinic/homoclinic) orbit detection** (5.5): the
  "hat algorithm" intersects box coverings of W^u(x*) and W^s(y*) to
  locate connecting orbits between steady states / periodic orbits —
  exactly the **energy-efficient transport channels** #314 is built on.
- **Patched three-body + controlled extension** (5.6): the patched-3BP
  manifold-tube method (replacing each tube intersection by **reachable
  sets** of a low-thrust controlled CR3BP) is demonstrated on a
  **low-thrust mission to Venus** (Earth-Venus transfer via L1/L2
  Lyapunov "gateways," reachable-set intersection in a common section,
  Color plate 4). Claimed payoff: ~3x flight time buys ~1/3 the fuel of
  the Hohmann VenusExpress baseline.

### #314 recommendation (heteroclinic mass-transport)

Junge-Dellnitz is the **canonical set-oriented methodology anchor** for
#314. Concretely it supplies:

1. **The box-covering subdivision algorithm** as a global,
   simulation-independent way to enumerate the invariant sets and
   unstable manifolds bounding a transport region — complements the
   project's trajectory-following manifold tools, which can miss
   measure-positive structures.
2. **The hat-algorithm connecting-orbit detector** (W^u ∩ W^s box
   intersection) as a constructive heteroclinic-connection finder. This
   is the operational core of a mass-transport calculation.
3. **The reachable-set / controlled-CR3BP extension** for the
   low-thrust variant of transport — relevant if #314 ever moves beyond
   ballistic heteroclinic connections.
4. The **almost-invariant-set / transition-rate** lineage (the
   Dellnitz-Junge-Padberg "transport" body of work that this chapter
   sits in) is the rigorous framework for *quantifying* mass transport
   (transition probabilities between almost-invariant regions) — the
   literal "mass-transport" in #314's name. The chapter itself focuses
   on the geometric (manifold/connecting-orbit) side; the
   transition-matrix side is in the cited Dellnitz-Junge references.

Cite: Dellnitz-Junge set-oriented methods; the GAIO toolbox; and for
the transport-matrix extension follow the chapter's Dellnitz/Junge
reference crumbs (refs to almost-invariant-set and Pareto-covering
work). This is a methodology anchor for KNOWN_CORPUS, not a cycler
family.

---

## Chapter 4 (Belbruno) — Low-Energy Transfers and Applications
### → #347 Phase 3; SUBSUMED by the Belbruno-2004 digest (#376)

This is a **condensed 19-page restatement** of the WSB/ballistic-capture
material already digested at full depth from Belbruno's 2004 book
(`docs/notes/2026-06-17-digest-belbruno-2004.md`, #376). Content
parity:

- Same **planar elliptic restricted four-body** model (E=P1, M=P2,
  sc=P3, Sun=P4), reducing to PCR3BP when m4=0 (pp 108-109).
- Same **capture problem** A1-A4 (ΔV0, ΔV1, ΔVC; capture = two-body
  Kepler energy E2 ≤ 0 at QF), eqs 4.7-4.8, pp 109-111.
- Same **Hill regions / zero-velocity curves / Lagrange-point Jacobi
  values** C4=C5=3 < C3 < C1 < C2, with the 3-4-digit approximations for
  C1, C2 (eq 4.4, p 116) — these are the standard CR3BP results
  (Szebehely territory; see the Szebehely digest).
- Same **weak stability boundary** W = {J=C, C<C1, E2≤0, ṙ2=0}
  (Def 4.3.3, p 119) and its 2-D annular structure r2 = f(θ2, e2).
- Same **chaos result** (Thm A / Thm B, pp 121-122): a hyperbolic
  network Λ exists on the extended WSB W̃H, so weak capture is chaotic
  (this is Theorem 3.58 of the 2004 book in compressed form).

### What Chapter 4 adds beyond the 2004-book digest

- **Section 4.5 "Origin of the Moon" (pp 123-125)** — the
  **Belbruno-Gott (2005)** application: the hypothetical Mars-sized
  Moon-forming impactor may have grown at a **Sun-Earth L4/L5 Trojan
  point**, then bifurcated out of a horseshoe orbit via a low-energy
  (WSB-type) escape to collide with Earth. This is a *scientific*
  (planetary-formation) application of low-energy dynamics, not a
  mission-design one — adjacent to, but not a, cycler. Worth a one-line
  KNOWN_CORPUS crumb only.
- Mentions **SMART-1** (ESA) as a second flown ballistic-capture
  example alongside Hiten (refs [16,17]).

### #347 Phase 3 recommendation

No new methodology beyond #376. If #347 Phase 3 fires, cite the 2004
book (deeper) and this chapter as the compact companion. The
Belbruno-Gott origin-of-Moon result (ref [5]) is a citable curiosity,
not a corpus entry.

---

## Chapter 6 (Ross) — Space Trajectory Optimization
### moderate relevance — optimization framing, not a cycler method

A **functional-analysis perspective on low-thrust trajectory
optimization**. Useful framing for any future ΔV/fuel-optimal cycler or
transfer work, but no direct cycler-family content.

Key takeaways:

- **Fuel cost is an L1-norm, not L2/quadratic** (Sections 6.3.1-6.3.2,
  pp 169-173). The proper minimum-fuel cost family is indexed by the
  Lp-norm of the thrust function; the *physically correct* fuel measure
  is the **L1-norm** of ‖T(t)‖. The ubiquitous quadratic cost JQ =
  ‖T‖²_L2 is **not** a fuel measure — using it gives non-optimal or
  "non-bang" controls (p 169). This is a clean, citable caution for any
  project optimization objective: minimum-ΔV ≠ minimum-quadratic.
- **Control-space geometry matters** (Section 6.2, pp 166-169): for real
  electric thrusters the admissible control set Ω can be **non-convex,
  disjoint** (power-limited, discrete thruster combinations), and the
  mismatch between control-space geometry (an l∞-ball for ungimbaled
  thrusters) and mass-flow geometry produces counter-intuitive
  time-optimal maneuvers.
- **Pseudospectral methods + the Covector Mapping Principle**
  (Section 6.x, pp 185-190): direct pseudospectral discretization (the
  DIDO lineage) with a covector-mapping theorem linking the discretized
  duals to the continuous costates; argues most "difficulties in
  solving optimal control" are discretization-method artifacts (some
  Runge-Kutta direct methods fail the covector map).

### Recommendation

Cite Ross Ch6 if the project adds a fuel-optimal-cycler or
low-thrust-transfer optimization track. The L1-vs-L2 cost result is the
single most reusable takeaway. No catalogue impact.

---

## Chapters 1, 2, 7, 8 — sampled, low project relevance (honest negative)

- **Ch1 Vallado — Perturbed Motion (p7):** Earth-satellite force-model
  survey (geopotential / drag / third-body / SRP / tides; secular vs
  short/long-periodic). Foundational orbit-determination material,
  well outside the multi-body cycler regime. No project use.
- **Ch2 Efroimsky — Gauge Freedom (p29):** the osculating-vs-
  non-osculating-element subtlety; Lagrange-constraint relaxation gives
  a gauge family of element sets; "osculation and canonicity are
  incompatible" for Delaunay/Andoyer elements (p 31). Mathematically
  elegant and *potentially* relevant if the project ever does
  perturbation-element-based averaging, but no current hook. Flag as a
  background reference for any future variation-of-parameters work.
- **Ch7 McInnes & Cartmell — Propellantless Propulsion (p199):** solar
  sails (orbital mechanics, artificial three-body equilibria, mission
  applications) and tethers. Solar-sail artificial-equilibria
  (Section 7.4) is a mild CR3BP-adjacent curiosity (sail-displaced
  Lagrange points) but not a cycler topic. Low relevance.
- **Ch8 — Cooperative Spacecraft Formation (p248):** model-predictive-
  control robustness for formation flying (open/closed-loop, virtual
  center, replan frequency). Control-systems material, no orbital-
  mechanics-cycler content. No project use.

---

## KNOWN_CORPUS methodology-anchor recommendations (no catalogue
writeback)

Recommend three methodology-bibliography anchors (substrate methods,
not cycler families):

1. **Guibout & Scheeres 2007 (Gurfil Ch3) — generating-function 2PBVP
   solver** → anchor for #380. Companion primary refs: Guibout PhD
   thesis (Michigan 2004); Guibout-Scheeres JGCD 27(4):693-704 (2003).

2. **Junge & Dellnitz 2007 (Gurfil Ch5) — set-oriented box-covering
   manifold/connecting-orbit method (GAIO)** → anchor for #314.
   Companion lineage: Dellnitz-Junge almost-invariant-set / transport
   work for the transition-matrix (true "mass-transport") side.

3. **Ross 2007 (Gurfil Ch6) — L1-cost / pseudospectral optimal control**
   → anchor for any future fuel-optimal track (the "fuel cost is L1,
   not quadratic" result).

The full citation for the volume:

  Gurfil, P. (ed.) (2007). *Modern Astrodynamics*. Elsevier
  Astrodynamics Series, Vol. 1. Butterworth-Heinemann / Elsevier.

## Action items for the parent

1. **#380 (BVP corrector):** Guibout-Scheeres generating-function 2PBVP
   method is the strongest new substrate this volume offers. The
   caustic = solution-multiplicity result is a principled
   "which-family-did-I-close-on" diagnostic. Order-2/3 version is
   tractable on existing STM infra; full order-N is multi-week.
2. **#314 (mass transport):** Junge-Dellnitz set-oriented / GAIO method
   is the canonical box-covering backbone (subdivision + hat-algorithm
   connecting-orbit detection + reachable-set controlled extension).
   Follow the Dellnitz-Junge transport-matrix reference crumbs for the
   transition-probability side.
3. **#347 Phase 3:** Belbruno Ch4 adds nothing over #376 except the
   Belbruno-Gott origin-of-Moon curiosity (Sun-Earth L4/L5 Trojan
   growth + horseshoe-bifurcation escape). Use #376 as the deep source.
4. **Scope correction:** no standalone Scheeres "strongly-perturbed"
   chapter exists in this volume — Scheeres is the Ch3 co-author.
   Ch7 is solar sails + tethers (McInnes-Cartmell). Recorded honestly.
5. **Honest negative on cycler corpus:** zero new cycler families. This
   is a methods volume; its value is three substrate-method anchors
   (#380, #314, fuel-optimization), not catalogue members.
