# Belbruno 2004 — Capture Dynamics and Chaotic Motions in Celestial Mechanics: With Applications to the Construction of Low Energy Transfers

**Deep-read verdict note, 2026-06-17 AET.** Textbook scope discipline:
TOC + foreword + preface + chapter abstracts pass, with deep-read of
Chapter 3 (sections 3.1 through 3.5) which is the WSB/BCT meat and
the only section with concrete cycler-relevant content.

## Header

- **Title:** *Capture Dynamics and Chaotic Motions in Celestial
  Mechanics: With Applications to the Construction of Low Energy
  Transfers*
- **Author:** Edward Belbruno (Princeton University; PhD student of
  Jurgen Moser; former JPL trajectory analyst)
- **Publisher:** Princeton University Press, Princeton and Oxford
- **Year:** 2004 (1st edition; camera-ready copy supplied by author per
  copyright page)
- **ISBN:** 0-691-09480-2 (cl, alk. paper)
- **LoC catalog:** QB362.M3B45 2004 (Many-body problem / chaotic
  behavior / celestial mechanics)
- **Pages:** 231 (xvi prefatory + 191 body + ~17 bibliography + ~3
  index); 3 chapters total
- **Foreword:** Jerry Marsden (Pasadena, Spring 2003)
- **Page size:** 444.6 x 653.3 pts (small octavo)
- **Format:** Graduate / research mathematics textbook with new
  theorems; Chapter 1 suitable for upper undergrad, Chapter 2-3 for
  graduate, sections 3.1.5 / 3.3.2 / 3.4 are the only purely-applied
  parts.

## Table of contents (chapter abstracts)

### Chapter 1 — Introduction to the N-Body Problem (pp 1-48)

Foundational mathematical setup. Sections:

- **1.1 The N-Body Problem** (p 1) — Newtonian n-body equations, n>=2;
  conservation laws (linear momentum, angular momentum, energy);
  centre-of-mass reduction. Brief discussion of Sundman's theorems and
  the Painleve conjecture (proven by Xia).
- **1.2 Planar Three-Body Problem** (p 9) — Jacobi coordinate
  formulation.
- **1.3 Two-Body Problem** (p 11) — classical solution; polar
  coordinates, momenta, sectorial area, Kepler ellipse geometry.
- **1.4 Regularization of Collision** (p 16) — Levi-Civita
  transformation; Kustaanheimo-Stiefel transformation.
- **1.5 The Restricted Three-Body Problem: Formulations** (p 24) —
  CR3BP / ER3BP in inertial and rotating coordinates; Jacobi integral;
  Lagrange points L1-L5; barycentric and primary-centred frames. THE
  KEY SECTION FOR PROJECT CR3BP-LITERACY CROSS-REFERENCE.
- **1.6 The Kepler Problem and Equivalent Geodesic Flows** (p 35) —
  geodesic equivalence on spaces of constant curvature; stereographic
  projection from S^2; Lobachevsky disc; NEW SHORT PROOF of the Kepler
  equivalence (Theorem 1.39). Mathematical-physics adjacent; no direct
  cycler application.

### Chapter 2 — Bounded Motion, Cantor Sets, and Twist Maps (pp 49-102)

KAM theory and area-preserving maps. Sections:

- **2.1 Quasi-Periodicity and the KAM Theorem** (p 50) — KAM applied
  to the restricted three-body problem; precessing elliptical motion;
  flow on torus.
- **2.2 The Moser Twist Theorem, Cantor Sets** (p 59) — twist maps;
  twist maps induced on a section of KAM tori; Cantor set
  construction.
- **2.3 Area-Preserving Maps, Fixed Points, Hyperbolicity** (p 66) —
  hyperbolic points; homoclinic and heteroclinic points; transverse
  intersections; hyperbolic networks; invariant manifolds near
  hyperbolic orbits.
- **2.4 Periodic Orbits and Elliptic Fixed Points** (p 78) —
  complicated dynamics near elliptic points; Poincare continuation
  method; periodic orbits of the first kind; surface of section;
  periodic point of type (p,q).
- **2.5 Aubrey-Mather Sets and the Restricted Three-Body Problem**
  (p 91) — NEW THEOREM 2.33 on Aubrey-Mather set existence in the
  restricted problem; lift maps on the covering space.

### Chapter 3 — Capture (pp 103-191)

Half the book. The cycler / cislunar / WSB / BCT material. Sections:

- **3.1 Introduction to Capture** (p 105) — permanent capture (Chazy,
  Hopf, Sitnikov, Alekseev); definitions of unbounded oscillatory
  motion, hyperbolic / parabolic orbits, ejection; temporary capture;
  BALLISTIC CAPTURE analytically defined (Def 3.11) via Kepler two-
  body energy sign; THE WEAK STABILITY BOUNDARY W introduced
  (Def 3.12, eq 3.9); the capture problem in the elliptic restricted
  four-body problem ER4BP-2D.
- **3.2 The Weak Stability Boundary** (p 120) — numerical algorithm
  definition of W via iterating periapsis-state stability under
  primary interchange; analytic W approximation (Lemma 3.21,
  eq 3.29); the C-problem (Def 3.23, restricted to C < C_1);
  visualization for Earth-Moon (mu=0.0123) at theta_2 = pi (Fig 3.7)
  showing the (r_2, e_2, -C) surface; C range
  [2.22, C_1=3.184] across the surface.
- **3.3 Existence of Primary Interchange Capture and an Application**
  (p 131) — Conley's transit orbits in the neck region near L_2;
  Lyapunov orbits and their invariant manifolds (Fig 3.8); proof
  that a transit orbit IS a ballistic capture transfer (Lemma 3.30);
  primary interchange capture (Lemma 3.32, infinitely many); applied
  example subsection 3.3.2 — first known realistic three-dimensional
  ballistic capture transfer to the Moon (r_13 = 154,089 km, 16 days
  TOF, near L_2 neck, ~75,000 km out of E-M plane). This subsection
  3.3.2 example is the INTERIOR WSB-TRANSFER.
- **3.4 A Low Energy Lunar Transfer Using Ballistic Capture**
  (p 144) — THE HITEN RESCUE. PR4BP-3D four-body model; backwards
  integration method from W-tilde at r_M + 100 km; matching arcs
  I (Earth phasing, 14 m/s burn at Q_0) and II (apoapsis-anchored
  arc to W-tilde) at Q_a (apoapsis ~1.5e6 km from Earth, 30 m/s
  delta-V_a mismatch); total delta-V = 44 m/s for capture into lunar
  orbit (capture delta-V_C = 0 by ballistic capture definition); Hiten
  arrival 2 Oct 1991. Forward method (differential-correction) added
  later in Remark 4. Subsection 3.4.3 historical perspective: SMART1
  (ESA, originally based on the interior transfer of subsection 3.3.2),
  Lunar A (Japan, planned), Planet B (Japan Mars mission using
  ballistic-ejection variant), HGS-1 1998 idea, Hughes-1992 ITER
  studies, Europa ballistic capture proposal (Sweetser); references
  Genesis, ISEE-3 / ICE halo orbit.
- **3.5 Parabolic Motion, Hyperbolic Extension of W** (p 156) — Easton-
  Robinson parabolic-orbit existence theorem; extension of W to W-tilde
  including slight hyperbolic motion (e_2 > 1 slightly); Jacobi-integral
  value on W-tilde for mu=0 parabolic orbits with respect to P_1 is
  C = +/- sqrt(2) (Lemma 3.34, eq 3.39).
- **3.6 Existence of a Hyperbolic Network on W-tilde-H** (p 165) — the
  book's main theorem (Theorem 3.58): a hyperbolic invariant set Lambda
  on the hyperbolic extension W-tilde-H of the weak stability boundary
  in the planar circular restricted problem with mu << 1, proving
  ballistic capture on W-tilde is a CHAOTIC process. Summaries of
  Moser's proof for the Sitnikov problem and Xia's proof for the
  planar circular restricted problem (refs [175], [227]).

### Bibliography (pp 193-208), Index (pp 209-211)

231 numbered references including 23 Belbruno entries (refs [21]-[42]).
Marsden / Lo / Koon / Ross "Low Energy Transfer to the Moon" (Celestial
Mechanics 2001, ref [126]) and their "Heteroclinic Connections" (Chaos
2000, ref [127]) appear. NO entries for: Galileo VEEGA, D'Amario,
Diehl, Roberts, Byrnes, Aldrin, cyclers. See section 8 below.

## What the textbook actually is

A graduate-to-research mathematics textbook focused on **chaotic
dynamics in the three-body problem**, with the construction of
**ballistic-capture lunar transfers** as the central applied motivator
(Chapter 3 is half the book). Belbruno positions the book around two
goals (preface, p xiii-xiv):

1. **Unify two capture definitions** — *permanent capture* (geometric,
   global, set of measure zero, chaotic per Alekseev / Moser / Xia) and
   *ballistic capture* (analytic, local, Kepler-energy-sign-based,
   used for spacecraft applications). The book's main result
   (Theorem 3.58) proves the hyperbolic-network set Lambda intersects
   the hyperbolic-extended weak stability boundary W-tilde, making the
   two capture types equivalent on that intersection and making
   ballistic capture on W-tilde a chaotic dynamical process.

2. **Extend Aubrey-Mather theory to the restricted three-body
   problem** — Theorem 2.33 (new) on existence of Aubrey-Mather sets
   near periodic orbits of the first kind.

A secondary goal is the new short proof (Theorem 1.39) of the
two-body / geodesic-flow equivalence.

The applied content sits in subsections 3.1.5, 3.3.2, and section 3.4
only. The rest is rigorous-mathematics: definitions, lemmas, theorems,
proofs. Marsden's foreword places the book in the lineage Poincare ->
Conley, McGehee -> Belbruno's WSB / Hiten work -> Lo / Koon / Marsden /
Ross's tube-manifold and Multi-Moon Orbiter program (the "1963 invariant
manifold tube ideas of Conley and McGehee," foreword p xii). The
*operational* deliverable of the book is the WSB / BCT methodology — used
in Hiten (1991), based on which SMART-1 was originally designed (per
subsection 3.4.3, p 155), and intended for Lunar A / Planet B / Europa /
HGS-1 / future lunar base missions.

## Identified project-relevant chapters

Chapter 3 sections 3.1, 3.2, 3.3, 3.4, and 3.5 (deep-read). Chapter 1
section 1.5 reviewed but not deep-read (foundational CR3BP /
formulations; covered in higher-fidelity form in our existing CR3BP
module and in Hintz 2023 sections 4.4-4.5). Chapter 2 KAM and Aubrey-
Mather theory NOT deep-read — relevant to #347 Floquet bifurcation
framework only at the level of "the same chaotic-dynamics machinery
applies to area-preserving maps near elliptic fixed points," which our
existing literature corpus (Doedel 1991, Restrepo-Russell, Wittal IAC-22)
already covers more directly. Chapter 3 sections 3.6 (hyperbolic-network
existence proof) skimmed — pure-mathematics adversarial-symbolic-
dynamics proof; the project does not yet have a use case requiring this
result at proof-detail level.

**Deep-read justification:**

- **3.1 (Introduction to Capture)** — gives the formal definition of
  ballistic capture (Def 3.11) and weak stability boundary W (Def 3.12,
  eq 3.9). Required for any future cislunar BCT work.
- **3.2 (The Weak Stability Boundary)** — analytic approximation of W
  (Lemma 3.21, eq 3.29); the C-problem (Def 3.23); the Earth-Moon W
  surface visualization. Foundational for #347 Phase 3 if it fires.
- **3.3 (Primary Interchange Capture and an Application)** — Conley
  transit orbits and the Lyapunov-orbit invariant-manifold construction
  near L_2; Lemma 3.30 that transit orbits ARE ballistic capture
  transfers. This is the conceptual bridge from CR3BP periodic orbits
  to BCT — directly relevant to the #316 Sun-Earth <-> Earth-Moon
  manifold transit framework.
- **3.4 (A Low Energy Lunar Transfer Using Ballistic Capture)** —
  the operational Hiten BCT recipe: PR4BP-3D model, backwards
  integration from W-tilde at r_M + 100 km, two-arc matching at the
  ~1.5e6 km apoapsis, ~44 m/s total delta-V vs ~200 m/s Hohmann.
  Canonical methodology for #347 Phase 3 cislunar BCT.
- **3.5 (Parabolic Motion, Hyperbolic Extension of W)** — the W-tilde
  extension with E_2 slightly positive; Jacobi constant C = +/- sqrt(2)
  for parabolic orbits with mu=0 (Lemma 3.34, eq 3.39). Background for
  Theorem 3.58 (the book's main result) and for any
  hyperbolic-tube-network analysis.

## Deep-read: Chapter 3 sections 3.1-3.5

### Section 3.1 (pp 105-120) — Introduction to Capture

Defines five capture / motion types in the three-dimensional elliptic
restricted three-body problem (E, M, S spacecraft P_3):

| Definition | Type | Page |
|---|---|---|
| Def 3.1 | Permanent capture: bounded forward, unbounded backward | 105 |
| Def 3.4 | Unbounded oscillatory motion: limsup r = inf, liminf r < inf | 109 |
| Def 3.6 | Parabolic orbit: |Q|->inf, |Q-dot|->0 | 109 |
| Def 3.7 | Hyperbolic orbit: |Q|->inf, |Q-dot|->const>0 | 110 |
| Def 3.9 | Temporary capture: bounded at t_0, unbounded at +/- inf | 110 |
| Def 3.11 | **BALLISTIC CAPTURE**: two-body Kepler energy E_2 <= 0 at t_1 | 111 |

**The two-body Kepler energy** (Def 3.10, eq 3.6) in P_2-centred
inertial coordinates X = (X_1, X_2, X_3):

  E_2(X, X-dot) = (1/2) |X-dot|^2 - mu / r_23

where r_23 = |X|, 0 <= mu < 1/2. Ballistic capture means E_2 <= 0 (P_3
is two-body-bound to P_2) — not permanent, just instantaneous-energy-
negative.

Permanent capture is geometric (whole-trajectory). Ballistic capture is
analytic (point-in-time energy sign). Belbruno's central insight: they
are equivalent on the hyperbolic-network-intersection set W-tilde ∩
Lambda (Theorem 3.58, the book's main result).

The **capture problem** (Section 3.1.5, p 114): given that the
spacecraft starts at periapsis r_13 = r_E + 200 km about Earth and ends
at periapsis r_23 = r_M + 100 km about Moon with E_2 <= 0, find the
trajectory. Modelled as a restricted four-body problem.

### Section 3.2 (pp 120-131) — The Weak Stability Boundary

The book's central mathematical object. The numerical-algorithm
definition iterates the stability of P_3 as it leaves P_2 along a
periapsis-state trajectory, classifying it stable (returns to
periapsis after one revolution) or unstable (escapes to P_1 or crashes
to P_2). W is the boundary in (r_2, theta_2, e_2) space at which
stability flips.

**The analytic approximation W** (Lemma 3.21, eq 3.29) restricts the
Jacobi-integral surface J^-1(C) to the periapsis subset (r-dot_23 = 0)
and the capture subset (E_2 <= 0):

  W = J^-1(C) ∩ Σ ∩ σ                  (eq 3.9, p 111)
  Σ = {x, x-dot | E-tilde_2 <= 0}      (Kepler-bound subset)
  σ = {x, x-dot | r-dot_23 = 0}        (periapsis subset)

With theta_2 in [0, 2pi] and e_2 in [0, 1], the explicit Jacobi-energy
expression on W is (eq 3.29):

  C = -r_2 * (+/- 2*sqrt(mu*(1+e_2)/r_2) + r_2)
      + mu * (1-e_2) / r_2 + A(r_2, theta_2)

where A is the residual J-tilde term, +/- gives direct/retrograde.

For Earth-Moon (mu = 0.0123), the Lagrange-point distances from P_2 are
d_1, d_2 ≈ 0.169 (eq 3.31 approximation valid to 4 digits for
mu <= 0.001) and C_1 ≈ 3.184. Numerical and applied evidence shows W
exists for C in [2.22, C_1] across the (r_2, e_2) surface (Figure 3.7
at theta_2 = pi).

### Section 3.3 (pp 131-144) — Primary Interchange Capture; the 1986 interior WSB-transfer

The Hill's regions of the CR3BP (Figure 3.6) evolve as C decreases past
C_2, C_1, 3 in succession (Theorem 3.27, Conley). For C just below C_2,
the Hill's region near L_2 has a "neck" through which the spacecraft
can transit from H_E to H_M. **Transit orbits** (Def 3.28) cross this
neck; Lemma 3.29 (Conley) proves they exist; Lemma 3.30 proves they
ARE ballistic-capture transfers (E_2(phi(t_F)) < 0 by continuity from
E_2(L_2) = -3^(-1/3) + (1/2)*3^(-2/3)*mu^(2/3) + O(mu^b), which
evaluates to -1.20187 < 0 for mu << 1).

**Asymptotic capture to a Lyapunov orbit** (Def 3.31) and **primary
interchange capture** (Lemma 3.32) extend this: infinitely many
transfers exist within a thin annular region A of width O(mu^(1/3))
near the stable manifold W_E^s of the Lyapunov orbit ζ(t), and they
spiral asymptotically to ζ(t) as t -> inf.

The applied example (subsection 3.3.2, p 142, the **interior
WSB-transfer**) is the first known realistic 3D ballistic-capture
transfer (citation [25] = Belbruno 1987, AIAA/DGLR/JSASS Electric
Propulsion conference proc. 87-1054): r_13 = 154,089 km initial,
delta-t = 16 days, trajectory goes ~75,000 km out of the E-M plane,
sharply down to ~100,000 km below it near L_2, then sharply up over
the north lunar pole to Q_F. C just below C_2. This is the interior
WSB-transfer prototype. SMART-1 was originally designed on this
template (p 155).

### Section 3.4 (pp 144-156) — Hiten BCT; the operational recipe

The exterior WSB-transfer, used by Hiten / MUSES-A on 2 Oct 1991.
This is the only operational ballistic-capture transfer flown when the
book was written. The setup:

- Spacecraft: MUSES-A renamed Hiten after MUSES-B comm failure;
  delta-V budget ~100 m/s (designed for Earth orbit, not lunar).
- Designed by Belbruno and Miller [38] = JPL IOM 312/90.4-1731-EAB,
  June 1990. Published in detail in [39] = J. Guidance, Control, and
  Dynamics 16(4):770-775, Jul-Aug 1993.
- Model: pseudo-circular restricted four-body problem PR4BP-3D
  (E, M, S, P_3; Sun gravity via DE403 ephemeris).
- **Backwards integration recipe** (p 145-146, Figure 3.14):
  1. Pick a target QF on W-tilde at r_23 = r_M + 100 km with
     e_2 ≈ 0.95 (i.e., on the W-tilde surface at the Moon).
  2. Integrate the four-body equations backward from QF for ~45 days.
     The trajectory reaches an apoapsis Qa at ~1.5e6 km (~4 x E-M
     distance) from Earth. This arc is **arc II** (ballistic-capture
     arc).
  3. From Hiten's actual Earth-orbit periapsis Q_0 (r_13 = 8,900 km),
     burn ~14 m/s and integrate forward ~100 days to reach near Qa
     **(arc I)**.
  4. Match at Qa: residual mismatch delta-V_a = 30 m/s (achieved by
     adjustment).
  5. Total delta-V = 14 + 30 = 44 m/s for ballistic capture at the
     Moon. **delta-V_C = 0** (capture is ballistic by definition of
     W-tilde).
- TOF: ~150 days vs ~5 days for Hohmann.
- delta-V comparison (p 147-148): Hohmann to Moon at r_23 = r_M +
  100 km circular requires ~200 m/s capture + ~648 m/s circularization
  = 848 m/s total. BCT requires 44 + 0 + 648 = 692 m/s total. **20%
  savings**, growing to 25% with the forward-method improvements
  (refs [33], [34]).
- Actual flight: delta-V_a = 34 m/s (vs 30 m/s designed). Hiten
  arrived 2 Oct 1991, modified to a 72,422 km altitude periapsis for
  L_4 / L_5 science. Returned to Moon 15 Feb 1992 for proper lunar
  orbit (delta-V_C = 82 m/s). Crashed 10 Apr 1993.

**Forward method** (Remark 4, p 148): later replaced the backward
method. A 2x2 differential-correction algorithm |V_0|, gamma_0 ->
r_23, i_M targets r_23 and lunar inclination i_M directly from Q_0.
No arc-matching needed. References [33], [34] for the algorithm
detail.

**The four-body decomposition** (subsection 3.4.2, p 150-154): the
CR4BP can be approximated as a switch between two CR3BPs based on the
distance r_13:

  CR4BP = CR3BEM (E,M,P_3, mu=0.012)    if r_13(Phi) <= rho
        = CR3BES (E,S,P_3, mu=3e-6)     if r_13(Phi) > rho

where rho = d(E,M) + rho* and rho* = 0.368 * d(E,M) is the maximum
extent of W in r_2 (at theta_2 = pi/2 or 3pi/2, e_2 = 1, C = C_1).
Table 3.2 (p 152) gives Jacobi-constant values at four checkpoints
along Phi (start, apoapsis, 7-days-pre-capture, capture). At
checkpoint E_3 (7 days pre-capture), C(E,M) = 3.17466 just below
C_1(E,M) = 3.184077, meaning the Moon's L_1 neck JUST opens, and Phi
threads it via the stable / unstable manifolds of the Lyapunov orbit.
This is the **dynamical mechanism**: the Sun shapes the apoapsis
geometry to drop the spacecraft into a near-stationary state relative
to the Moon, then the Moon-system Lyapunov-manifold network funnels
it through L_1 to W-tilde.

The Sun's role is energy-removal during the fall-back from Qa:
"prior to capture the point P_3 moves in approximate parallel
formation with M for about a week" (p 149). This is the
characteristic-signature of an exterior WSB-transfer.

### Section 3.5 (pp 156-165) — Parabolic motion, the W-tilde extension

For the mu=0 limit, parabolic orbits with respect to P_1 have Jacobi
constant C = +/- sqrt(2) (Lemma 3.34, eq 3.39 — direct + sign,
retrograde - sign). For mu > 0, the value is C = +/- sqrt(2) + O(mu).
The W-tilde extension includes E_2 slightly positive (pseudo-ballistic
capture, e_2 slightly > 1) and r-dot_23 slightly nonzero (not exactly
at periapsis).

This extension is the substrate for Theorem 3.58 (the book's main
result, section 3.6): the hyperbolic invariant set Lambda intersects
W-tilde-H non-trivially, making ballistic capture chaotic.

## WSB formal definition (canonical Belbruno)

**Weak Stability Boundary W** (Def 3.12, p 112, with eq 3.9 on p 111):

  W = J^-1(C) ∩ Σ ∩ σ

where for the planar circular restricted three-body problem in
barycentric rotating coordinates x = (x_1, x_2):

  J^-1(C) = the Jacobi-integral surface at constant C (Hill region
            boundary)
  Σ = {(x, x-dot) | E-tilde_2 <= 0}    (Kepler two-body bound to P_2)
  σ = {(x, x-dot) | r-dot_23 = 0}      (periapsis with respect to P_2)

The Kepler energy E_2 is computed in P_2-centred inertial coordinates:

  E_2(X, X-dot) = (1/2) |X-dot|^2 - mu / r_23

with r_23 = |X|, mu the smaller-primary mass parameter.

**Explicit W expression on the (r_2, theta_2, e_2) surface**
(eq 3.29, Lemma 3.21):

  C = -r_2 * (+/- 2 * sqrt(mu*(1+e_2)/r_2) + r_2)
      + mu * (1-e_2) / r_2
      + A(r_2, theta_2)

where A is the residual Jacobi term, theta_2 in [0, 2pi], e_2 in
[0, 1], r_2 >= 0, mu in (0, 1/2), +/- is direct/retrograde.

**Validity domain (Def 3.22):** W is a valid approximation of the
numerical-algorithm W when C < C_1 (the L_1 Jacobi constant).
For Earth-Moon (mu = 0.0123), C_1 = 3.184 and W exists across
C in [2.22, C_1] (verified by Figure 3.7).

**Hyperbolic extension W-tilde** (section 3.5): same construction
but with E_2 slightly > 0 (pseudo-ballistic capture) and r-dot_23
slightly nonzero.

**Generalization-of-Lagrange-points view** (p 128): on or near W,
F + G ≈ 0 with x-dot ≠ 0 — i.e., gravity + centrifugal-force balance
with non-zero velocity. The Lagrange points are the F + G = 0,
x-dot = 0 special case. WSB is the generalization to non-static
balance.

## BCT construction methodology (canonical Belbruno recipe)

Two recipes are given. The **backwards integration method** (used
operationally for Hiten 1991) and the **forward differential-
correction method** (the modern improvement).

### Backwards integration method (sections 3.4.1, Figure 3.14)

For an exterior WSB-transfer in PR4BP-3D (CR4BP including Sun gravity
via planetary ephemeris, e.g. DE403):

1. **Target state** at QF on W-tilde at the Moon:
   r_23 = r_M + r_capture (e.g., r_M + 100 km),
   e_2 prescribed (e.g., 0.95), theta_2, alpha_2, beta_2 chosen
   (orientation degrees of freedom).
2. **Backwards-integrate** the four-body equations from QF for
   ~45 days. The trajectory reaches an apoapsis Qa at ~3.5-4x the
   Earth-Moon distance (~1.5e6 km). Stop at Qa. Call this **arc II**.
3. **Forward-integrate** from spacecraft's Earth-orbit periapsis Q_0
   (any r_13) with a small initial burn delta-V_0 (e.g., 14 m/s for
   Hiten from r_13 = 8,900 km) for ~100 days to reach near Qa.
   Call this **arc I**.
4. **Match the arcs at Qa**: residual mismatch delta-V_a (e.g.,
   30 m/s for Hiten) is the apoapsis-correction burn.
5. **Total transfer delta-V** = delta-V_0 + delta-V_a + delta-V_C,
   where delta-V_C = 0 by definition (ballistic capture at QF).

This method is "unwieldy" (Remark 4, p 148) — convergence and
parameter-sweep are hard because matching r_13, i_E, t_0 from arc I
to W-tilde target via arc II requires manual iteration.

### Forward differential-correction method (Remark 4, p 148-149)

A 2x2 targeting algorithm:

  |V_0|, gamma_0  ->  r_23, i_M

Control variables at Q_0: initial-velocity magnitude |V_0|, flight-path
angle gamma_0 (angle between V_0 and the normal direction N to Q_0 in
the (Q_0, V_0) plane). Target variables at QF: r_23 (capture distance)
and i_M (lunar-equator inclination). Other parameters (t_0, r_13, i_E)
are prescribed and held fixed.

**Heuristic**: choose |V_0| so the trajectory has an apoapsis near
1.5e6 km from Earth. Adjust t_0 to align E-M-S phasing. P_3 then
"naturally" reaches W-tilde with e_2 < 1 (no need to prescribe e_2).
Method is described in refs [33] Belbruno-Carrico 2000 and [34]
Belbruno-Humble-Coil 1997 (USAF Academy Blue Moon Mission).

### Alternative: invariant-manifold method (Koon-Lo-Marsden-Ross)

Reference [126] = Koon-Lo-Marsden-Ross, "Low energy transfer to the
Moon", Celest. Mech. Dyn. Astron. 81:63-73, 2001. Uses the global
stable/unstable invariant manifolds of the L_1 / L_2 Lyapunov orbits
in the CR3BP. Belbruno notes this method requires C just below C_1
and therefore omits "a significant portion of the orbit space" where
the exterior WSB-transfer actually lives (C substantially less than 3).

### Genetic-algorithm method

Reference [42] = Bello-Mora et al. IAF-00-A.6.03, 2000. Mentioned as
an alternative search method.

## Galileo attribution check

**Result: zero hits across the entire textbook**, including the 231-
entry bibliography and the index.

The full PDF was extracted to plain text and grep'd for: `Galileo`,
`VEEGA`, `D'Amario`, `Diehl`, `Byrnes`, `Roberts`, `cycler`, `Aldrin`.
Only one hit: the literal string "gravity assists" appears once, in
the title of bibliography entry [124] = Kawaguchi et al.
"On making use of lunar and solar gravity assists in lunar a, planet
b missions", Acta Astronautica 35:633-642, 1995 — concerning Japanese
missions, not Galileo.

The author's own publication record (refs [21] through [42], 22
entries) covers:

- Kuiper Belt 2:3 resonance (1990s)
- Two-body geodesic flows (Celestial Mechanics 1977)
- Restricted-three-body regularization (Celestial Mechanics 1981, two
  papers)
- **The 1987 interior WSB-transfer** (ref [25] = AIAA 87-1054, May 1987
  — Lunar GAS mission proposal, the original SMART-1-template paper)
- **The 1990 Hiten BCT paper** (ref [38] = JPL IOM 312/90.4-1731-EAB,
  15 Jun 1990, Belbruno + Miller; published as ref [39] = JGCD 16(4),
  Jul-Aug 1993)
- 1990 Earth-Moon ballistic-capture/escape (AIAA 90-2896)
- 1994 four-body invariant-manifold analysis (CRM Research Report 270,
  Barcelona)
- 1997 Comet resonance-hopping (with Marsden, Astronomical Journal)
- 1997 Fast-resonance-shifting (Annals NYAS 822)
- 1997 USAF Blue Moon (with Humble + Coil)
- 2000 Belbruno-Carrico WSB ballistic lunar capture (AIAA 2000-4142)
- 2000 single-parameter capture prediction (JPL Contract 1213585)
- 2002 Analytic estimation of WSB (Contemp. Math. 292:17-47)
- 2002 Hill's problem periodic orbits (Kluwer 2003)
- 2003 HGS-1 retrospective (with Ridenoure and Ocampo, submitted)

Belbruno's bibliography includes **zero papers on the Galileo mission**,
**zero papers from before 1977**, and **zero papers in any
cycler-relevant or VEEGA-trajectory line**. His career profile in
this book is entirely: WSB / BCT / chaotic-capture mathematics +
operational implementation on Hiten + planning involvement on
SMART-1, USAF Blue Moon, Lunar A, Planet B, Europa, HGS-1.

**Verdict for the project's KNOWN_CORPUS "Diehl-Belbruno-Roberts /
Galileo 1986" attribution:**

This 231-page textbook with a comprehensive author-publication
bibliography (including JPL IOMs and obscure conference proceedings)
does NOT claim Galileo VEEGA contribution. Belbruno was at JPL during
the Galileo-development era (refs [25] / [38] both list JPL
affiliation), but his book — written by him, on his own work, with a
self-attribution-friendly chapter 3 acknowledgments page — is silent
on Galileo. Section 3.4.3 historical perspective (p 154-156) is the
natural place where a Galileo contribution would be claimed if it
existed, and it is not mentioned.

**This is corroborative evidence for Agent D's D'Amario 1992 finding**
(Belbruno does NOT appear in the D'Amario 1992 reference list). The
"Diehl-Belbruno-Roberts" KNOWN_CORPUS attribution for Galileo 1986
remains SUSPECT-INCORRECT. The triangulation now has TWO independent
sources: (a) D'Amario 1992 (the canonical Galileo VEEGA-design paper)
omits Belbruno from its references, and (b) Belbruno 2004 (the
canonical Belbruno-self-attribution publication) is silent on Galileo.
The likely correction is that the historical Galileo VEEGA design team
is D'Amario, Byrnes, et al. (per D'Amario 1992) and the Diehl /
Belbruno / Roberts triple is a transcription error somewhere in the
project's historical-reading chain.

A possible adjacent attribution worth checking separately: Belbruno's
1986 work was on the **Lunar GAS** (gravity-assist surveys) mission
proposal (ref [25] is dated May 1987 AIAA Electric Propulsion paper),
not on Galileo. The "1986 GAS work" and "1986 VEEGA work" are
plausibly confused upstream of the project.

## #347 Phase 3 implications (if Phase 3 fires)

If #347 graduates from current Phase 1 (Floquet bifurcation mining in
Earth-Moon CR3BP) through Phase 2 (mu-scaling to other systems) to
**Phase 3 (cislunar BCT integration)**, Belbruno 2004 supplies:

1. **The analytic-W formula (eq 3.29, Lemma 3.21)** as the surface on
   which Phase 3 capture-targeting must occur. This is mu-parameterized
   (single mu input), so the mu-scaling lineage from Phase 2 carries
   through naturally.

2. **The numerical-algorithm W definition** (section 3.2.1, p 122-125) —
   stability-class iteration from periapsis state, classified into one
   of {stable, unstable, capture, escape, primary-interchange} based on
   future trajectory. Implementation: O(N_periapsis x N_orbits)
   trajectory propagations in a CR3BP shooter with categorical labels.
   Phase 3 should test against this numerical definition as the ground
   truth.

3. **The forward differential-correction method** (Remark 4 of
   section 3.4.1, p 148-149) as the operational targeting algorithm:
   2x2 (|V_0|, gamma_0) -> (r_23, i_M). Maps to a finite-difference
   Jacobian + Newton iteration around any nominal trajectory. This is
   a small, tractable implementation.

4. **The four-body CR4BP decomposition** (eq on p 151) as a
   coarse-grained CR3BP-pair stand-in for full ephemeris when ephemeris
   integration is too expensive: CR3BEM when r_13 < rho = d(E,M) +
   0.368 d(E,M); CR3BES otherwise. This is a Phase-3-ready model
   reduction.

5. **The hyperbolic-network result (Theorem 3.58)** as a **theoretical
   anchor** that the cislunar BCT region admits chaotic ballistic
   capture. Phase 3 doesn't need the proof — it needs the *existence
   guarantee*, which justifies investing in cislunar-BCT search at all.

6. **The Hiten delta-V signature** as a benchmark: 44 m/s total
   delta-V (vs 200 m/s Hohmann capture), 150 days TOF (vs 5 days
   Hohmann), Q_a apoapsis at ~1.5e6 km. Any Phase 3 cislunar BCT
   result should be in this delta-V / TOF / apoapsis ballpark or
   present a defensible reason to differ.

**Phase 3 cost estimate (rough)**: 2-4 weeks for a Belbruno-style W
implementation in the existing CR3BP module, plus 1-2 weeks for a
PR4BP-3D / ephemeris bridge. The forward method is the right entry
point (much easier than backward) and maps cleanly onto the project's
existing differential-correction infrastructure.

**What Belbruno 2004 does NOT supply for Phase 3:**

- Modern porkchop-plot tooling (better: Hintz 2023 Appendix D for
  Mars; need equivalent for cislunar BCT).
- Three-body Lambert-like patched algorithms (better: Sanchez et al.
  / FBS line per #243).
- Real-ephemeris validation against a flown mission (Hiten is the only
  reference data point; SMART-1, ARTEMIS could supplement).
- Floquet-bifurcation analysis at all (this is #347 Phase 1's
  contribution).

## KNOWN_CORPUS impact

**No new cycler-family anchor.** WSB orbits are not cyclers — they are
single-leg low-energy lunar-transfer trajectories whose hallmark is
ballistic capture (zero capture delta-V) at the target Moon-region
W-tilde. The mathematical machinery is the same chaotic-dynamics /
invariant-manifold lineage that informs cislunar low-energy work, but
no specific cycler member is named in the book. Even the closest-
adjacent concept (the Multi-Moon Orbiter, mentioned in Marsden's
foreword) is not a cycler either — it is a phase-locked-orbit
hierarchy.

**Methodology bibliography entry recommended (do not edit catalogue).**
The book should be cited in any future cislunar BCT / WSB methodology
context, e.g., from #347 Phase 3, from #316 Sun-Earth <-> Earth-Moon
manifold transit reporting, or from a future cislunar-BCT track. The
recommended citation:

  Belbruno, E.A. (2004). Capture Dynamics and Chaotic Motions in
  Celestial Mechanics: With Applications to the Construction of Low
  Energy Transfers. Princeton University Press. ISBN 0-691-09480-2.

Specific equations and theorems to cite:

- Def 3.10 / eq 3.6 (Kepler two-body energy E_2)
- Def 3.11 (ballistic capture definition)
- Def 3.12 + eq 3.9 (WSB definition W = J^-1(C) ∩ Σ ∩ σ)
- Lemma 3.21 / eq 3.29 (analytic-W approximation)
- Def 3.22 (validity domain C < C_1)
- Lemma 3.30 (transit-orbits-are-BCTs)
- Section 3.4.1 + Figure 3.14 (Hiten backwards-integration recipe)
- Remark 4, section 3.4.1 (forward differential-correction method)
- Section 3.4.2 + eq on p 151 (four-body CR3BP-pair decomposition)
- Theorem 3.58 (hyperbolic network on W-tilde, ballistic capture is
  chaotic)

**Galileo attribution correction recommended.** Update the
"Diehl-Belbruno-Roberts / Galileo 1986" KNOWN_CORPUS entry to
SUSPECT-INCORRECT-CONFIRMED, with the two corroborating sources:
D'Amario 1992 omits Belbruno, and Belbruno 2004 makes no Galileo
claim. Suggest investigating whether the actual attribution should be
to D'Amario / Byrnes / et al. and whether the "Belbruno 1986" entry
in upstream sources actually refers to Belbruno's Lunar GAS work
(ref [25], May 1987 AIAA proceedings).

## Action items for the parent

1. **Galileo / KNOWN_CORPUS correction (independent corroboration).**
   Belbruno 2004 makes zero Galileo claim across 231 pages including a
   detailed self-bibliography. Combined with Agent D's D'Amario 1992
   finding, the "Diehl-Belbruno-Roberts" attribution is now
   double-corroborated as SUSPECT-INCORRECT. Recommend marking the
   KNOWN_CORPUS entry accordingly and opening a correction-research
   task. Plausible upstream confusion: Belbruno's 1986/1987 **Lunar
   GAS** (gravity-assist surveys) proposal -> historian conflation
   with VEEGA-era Galileo work.

2. **Cislunar BCT methodology anchor lodged.** Belbruno 2004 is now the
   canonical citation for any future #347 Phase 3, #316 followup, or
   net-new cislunar BCT track. The specific equations / theorems /
   recipe steps to cite are listed in the KNOWN_CORPUS impact section
   above. NO catalogue edit recommended at this time (no cycler
   member).

3. **Higher-priority-than-Parker-2007 verdict CONFIRMED.** The recovery
   agent's flag was correct. Belbruno 2004 is the canonical WSB / BCT
   reference; Parker 2007 (482-page PhD thesis on low-energy transfers)
   is the next layer of operational detail and only needed if #347
   Phase 3 actually fires. Belbruno 2004 supplies the foundational
   definitions; Parker 2007 would supply the operational depth on
   specific lunar missions.

4. **Phase 3 capture-readiness assessment.** If #347 Phase 1 graduates
   to Phase 3, the implementation cost using Belbruno 2004 alone is
   ~3-6 weeks (W-implementation + forward-method-targeting + four-body
   CR3BP-pair decomposition). The book's analytic W formula (eq 3.29)
   is implementable directly in the existing CR3BP module.

5. **Honest negative on cycler corpus.** Belbruno 2004 supplies no new
   cycler-family anchors. WSB / BCT is a single-leg low-energy capture
   technique, not a periodic-encounter trajectory. The book is a
   methodology reference for FUTURE work, not a corpus contribution
   for CURRENT mining.

6. **Recursive bibliography crumbs.** Several entries in Belbruno's
   bibliography may be worth follow-up:
   - **Ref [218] = Sweetser et al. AAS 97-174 (1997)** — Europa BCT
     proposal. Adjacent to project's Jupiter-system work.
   - **Ref [126] = Koon-Lo-Marsden-Ross "Low energy transfer to the
     Moon" (CMDA 81:63-73, 2001)** — likely already in project corpus.
   - **Ref [127] = Koon-Lo-Marsden-Ross "Heteroclinic Connections
     and resonance transitions" (Chaos 10:427-469, 2000)** — confirms
     hyperbolic network for Sun-Jupiter system numerically; may
     supply a separate methodology anchor.
   - **Ref [42] = Bello-Mora et al. IAF-00-A.6.03 (2000)** — early
     genetic-algorithm WSB search; possibly relevant if the project's
     Forge ML-surrogate line (#226 / #242 / #243) explores
     evolutionary methods.
   - **Ref [157] = Mendell, Space Policy 17:13-17, 2001** — "Gateway
     for human exploration: the weak stability boundary". Adjacent to
     human-cycler scope expansion.
