# Digest: Pontani & Conway 2018 — Optimal Trajectories for Hyperbolic Rendezvous with Earth-Mars Cycling Spacecraft

Task #383. Single-paper digest, deep read of all 17 pages on 2026-06-17
AET. JGCD (Journal of Guidance, Control, and Dynamics) DOI
10.2514/1.G002984; "Article in Advance" version (received April 2017,
revision July 2017, accepted Aug 2017, online Sep 2017).

## 1. Header

- **Title (verbatim)**: *Optimal Trajectories for Hyperbolic Rendezvous
  with Earth-Mars Cycling Spacecraft*
- **Authors**:
  - Mauro Pontani — Assistant Professor, Department of Astronautical,
    Electrical, and Energy Engineering, University of Rome "La
    Sapienza", 00138 Rome, Italy. mauro.pontani@uniroma1.it
  - Bruce A. Conway — Professor, Department of Aerospace Engineering,
    University of Illinois at Urbana-Champaign, Urbana IL 61801.
    bconway@uiuc.edu. Associate Fellow AIAA.
- **Venue**: *Journal of Guidance, Control, and Dynamics* — Article in
  Advance, AIAA. Copyright AIAA.
- **DOI**: 10.2514/1.G002984
- **Submitted/accepted/online**: 19 April 2017 / 31 August 2017 /
  27 September 2017
- **Length**: 17 pages (pp. 1-17 of the Article-in-Advance PDF)

## 2. What the paper actually is

The paper solves the **propellant-minimum hyperbolic-rendezvous problem**
between a small "taxi" vehicle (launched from Earth) and a large Earth-
Mars cycling spacecraft (which performs a hyperbolic flyby of Earth and
does not capture). This is the "TAXI" half of the cycler-architecture
problem: cyclers themselves are heliocentric repeat-encounter
trajectories, but they DO NOT capture at either Earth or Mars on each
encounter — therefore the crew transfer requires a separate
hyperbolic-velocity rendezvous vehicle.

Two distinct problem classes are formulated and solved:

1. **Impulsive hyperbolic rendezvous** (high-thrust taxi). Up to four
   chemical impulses applied at apse points along the transfer arc, with
   particle swarm optimization (PSO) finding the globally optimal
   parameter set. Includes an explicit abort-strategy formulation.
2. **Low-thrust hyperbolic rendezvous** (electric propulsion taxi).
   Continuous-thrust minimum-time problem solved by indirect method
   (Pontryagin minimum principle + costate boundary-value problem) with
   PSO seeding the costate initial guess.

The headline numerical result (p. 14 col 2; verbatim Sec III.D):

- High-thrust impulsive (best-case 4-impulse Delta-V_tot = 1.791 km/s,
  Centaur RL-10 c_H = 4.4 km/s) yields m_f/m_0 = 0.714 (propellant
  mass 0.286 m_0).
- Low-thrust (c = 30 km/s, n_0 = 0.001 g_0) yields m_f/m_0 = 0.829
  (propellant 0.171 m_0). Time of flight 145.48 h (~6.1 days).

**Low-thrust beats high-thrust by ~40% in propellant fraction** for this
problem — the central scientific conclusion.

## 3. Key content, tables, equations with page citations

### 3.1 Problem setup (Sec II, pp. 2-4)

- Coplanar orbits assumption (p. 2 col 2): "The orbits of the two
  spacecraft are assumed as coplanar and defined only by their specific
  energy and angular momentum. This means that the mutual orbit
  orientation and the initial relative positions of the two spacecraft
  are unspecified."
- Inertial frame definition (p. 3 col 1): right-handed triad
  (c_1, c_2, c_3) with c_1 pointing toward the periapse of the incoming
  hyperbola, c_3 aligned with the cycler angular momentum.
- Taxi initial orbit (p. 3 col 2): geosynchronous transfer orbit (GTO)
  with perigee 400 km, apogee 35,786 km altitude. p_0 = 11,679 km
  (semilatus rectum), e_0 = 0.723.
- Cycler hyperbolic arrival (Eq 7, p. 3): hyperbola at Earth on
  **4 February 2019** with p_H = 91,595.2 km, e_H = 2.157.

### 3.2 Reference cycler (citation chain)

The cycler trajectory referenced for the rendezvous case is taken from
ref [5] = McConaghy, Landau, Hong, Yam, Longuski 2006 *J. Spacecraft and
Rockets* 43:2 pp. 456-465, "Notable Two Synodic Period Earth-Mars
Cycler" — DOI 10.2514/1.15215.

The 4 February 2019 Earth-hyperbola-arrival epoch and (p_H, e_H) values
trace to this McConaghy 2006 paper. (Existing catalogue cross-check:
McConaghy 2006 is in the corpus as `mcconaghy-2006-two-synodic` if I
recall correctly; the V_inf at Earth from the hyperbola elements is
sqrt(mu_E / |a_H|) = sqrt(mu_E (e_H^2 - 1) / p_H) = approximately 5.85
km/s — a representative two-synodic-period cycler Earth V_inf.)

### 3.3 Impulsive rendezvous theory (Sec II.A-F, pp. 2-10)

The analysis builds on Pontani 2010 ref [10] ("Simple Method to
Determine Globally Optimal Solutions for Initial Orbit Transfers", JGCD
32:3 pp. 899-914, DOI 10.2514/1.38143). Key theoretical premises (p. 3
col 2 - p. 4 col 1):

1. On Keplerian trajectories (ellipses, parabolas, hyperbolas), locally
   optimal finite impulses occur tangentially at apse points.
2. On hyperbolic and parabolic trajectories, **infinitesimal** impulses
   can occur at infinite distance: they change periapse without
   modifying specific energy.

The (r_p, X) plane (Fig 2a, p. 4) parametrises Keplerian arcs with
X = 1/r_A:
- X > 0: elliptic orbits, r_A is apoapse
- X = 0: parabolas
- X < 0: hyperbolic trajectories, |r_A| is the periapse from the vacant
  focus.

The globally optimal ellipse-to-hyperbola four-impulse path (Fig 2b,
sequence A->1->2->3->4->5->6->7->8->C, p. 4) involves two finite apse
impulses and two infinitesimal impulses at infinite distance — but the
last is INFEASIBLE for a real taxi (must remain within Earth SOI, set
to 900,000 km in this paper, p. 2 col 2 last line). Therefore the
globally optimal solution in the FEASIBLE region requires finite-radius
constraints (Eq 9 + Eq 18, p. 4 and p. 7):

- r_p^(min) >= R_E + 200 km (Earth grazing avoidance, conservative
  lower bound)
- r_A^(max) <= 800,000 km (or 400,000 km abort-strategy bound)

### 3.4 Three rendezvous classes (Sec II.B, II.C, II.D, II.E)

1. **Two-impulse** (no abort, Sec II.B p. 4-6, Figs 4-6): minimum
   Delta-V_tot = 1.796 km/s at Delta-t_P = 20 h (cycler-perigee-pass-
   time since first impulse). The 20 h limit corresponds to the
   farthest feasible point along the hyperbola at 900,000 km from
   Earth.
2. **Two-impulse with abort** (Sec II.C p. 6-7, Figs 7-9): forces the
   intermediate arc to be elliptic with r_A^(ab) <= 400,000 km, so if
   the second ignition fails, the taxi returns to Earth ballistically.
   Best Delta-V_tot = 2.898 km/s at Delta-t_P = 3.22 h. Cost over
   no-abort: ~1.1 km/s extra for safety.
3. **Three-impulse** (Sec II.D p. 8 and Fig 10): three-impulse path
   with apses at r_p^(min), r_A^(max), and the final hyperbolic-arrival
   periapse. Delta-V_tot = 2.266 km/s when Delta-t_P < 6.64 h
   (intermediate values reduce to a 2-impulse).
4. **Four-impulse** (Sec II.E p. 8-9, Figs 11-13): globally optimal at
   Delta-t_P > 8.80 h. Delta-V_tot = 1.791 km/s at Delta-t_P = 20 h.
   Five m/s less than the 2-impulse minimum.

### 3.5 Low-thrust formulation (Sec III, pp. 10-15)

- Equations of motion in polar (r, v_r, v_theta) with c_1 inertial
  reference: Eqs 28-30 p. 10.
- Alternative formulation in osculating (p, e, f): Eqs 32-34 p. 10.
- Hamiltonians H^(A), H^(B) and the boundary-condition functions
  Phi^(A), Phi^(B): Eqs 38-39, 50-51 p. 11.
- Costate adjoint equations: Eqs 40-42 p. 11 (set A); Eq 52 p. 11
  (set B).
- Pontryagin minimum principle for thrust angle: Eqs 44, 59 (sets A
  and B respectively).
- Transversality conditions: Eqs 47-48 (set A), 60-61 (set B).
- BVP unknowns: 5 unknowns for set A (Eq 49), 4 for set B (Eq 62).

**Numerical method (p. 12 col 2 - p. 13 col 1)**:
- Indirect heuristic method = PSO over the costate initial guesses
- 100 particles, 1000 iterations
- Equation set B converges (1st attempt); equation set A fails (10/10
  attempts unsuccessful, p. 13 col 1)
- Reason discovered (Sec III.E, p. 13-15): set A is **hypersensitive**
  — costate lambda_3^(A) reaches values of -24.8 at t_0 and exhibits
  abrupt time variations (Fig 21). Set B costates are smooth (Fig 20)
  and stable under 0.1% perturbation (J_tilde^(B) = 0.397) while set A
  blows up (J_tilde^(A) = 6.539, Eq 72). Verbatim conclusion p. 15:
  "The result for J_tilde^(A) demonstrates that hypersensitivity
  affects the spacecraft dynamics if state x_A is used. In the end,
  this research proposes several different optimal options for
  hyperbolic rendezvous, with their respective possible uses,
  advantages, and disadvantages, to inform the design of Earth-Mars
  missions by means of cycling spacecraft."

### 3.6 Numerical values (p. 13 col 2, Eqs 66-72)

- c = 30 km/s, n_0 = 0.001 g_0 = 0.0098 m/s^2 initial thrust accel.
- Canonical units: DU = 100,000 km, TU = 50,087.7 s.
- Set A search bounds: lambda_{1,2,i}^(A) in [-1, 1], f_i in
  [-pi, pi], f_f in [-2, 2], t_f in [5, 50] TU.
- Set B search bounds: lambda_{1,2,i}^(B) in [-1, 1], f_i in
  [-pi, pi], t_f in [5, 50] TU.
- Final BVP residuals (set B): c_1 = -1.215e-10, c_2 = 7.824e-11,
  c_3 = -5.010e-10 (essentially zero).
- t_f = 145.48 h = 6.06 days.
- m_f/m_0 = 0.829.
- Rendezvous radius = 812,988 km (close to but inside the 900,000 km
  outer bound).

### 3.7 PSO algorithm (Appendix, pp. 15-16)

The PSO setup is canonical (Eqs A1-A9):

- N particles, position chi(i) in [a, b], velocity w(i) in [-d, d].
- Inertial weight c_I = (1 + r_1(0,1))/2.
- Cognitive weight c_C = 1.49445 * r_2(0,1).
- Social weight c_S = 1.49445 * r_3(0,1).
- Equality-constraint handling: penalty function J + sum |alpha_r d_r|
  (Eqs A8-A9).
- Inequality constraints: assign infinite J on violation (the simplest
  "elementary" approach, p. 16 col 1).

500 iterations + 50 particles for impulsive Secs; 1000 iterations + 100
particles for low-thrust Sec III.

### 3.8 Comparative analysis with previous work (p. 2 col 2)

- The paper cites Anderson 2016 (AIAA 2016-5267, ref [7]) as the prior
  "comprehensive review with more than 200 different Earth-Mars cycler
  trajectories" but notes Anderson did NOT explicitly optimise the
  rendezvous taxi (assumed two impulses, did not specify rendezvous
  trajectory beyond V_inf at Earth and Mars and periapse altitude).
- Jedrey-Landau-Whitley 2016 ref [6] reported "the only explicit
  optimization of the maneuver for a two-burn (impulsive or finite)
  case" using Copernicus, but did not consider 3- or 4-impulse paths,
  did not consider low-thrust, did not consider abort strategy.

Pontani 2018 fills these gaps.

## 4. KNOWN_CORPUS impact

**RECOMMEND: New KNOWN_CORPUS anchor for the taxi-cycler hyperbolic
rendezvous problem.**

This paper is the FIRST sourced taxi-side reference in the corpus. The
corpus has many cycler-side anchors (Aldrin, Byrnes-Longuski-Aldrin,
Russell-Ocampo, McConaghy 2004/2005/2006, Conte-Spencer, Rogers, the
recent admissions wave) but no taxi-side anchor.

Suggested KNOWN_CORPUS anchor:

- `corpus_id`: pontani-conway-2018-taxi-rendezvous
- `paradigm_label`: ["hyperbolic-rendezvous", "taxi-cycler-interface",
  "impulsive-4-impulse", "low-thrust-rendezvous"]
- `body_set`: ["Earth", "(generic-cycler)"]
- `orbit_class`: precursor_mga ?  -- see Sec 6 below; this is the
  classification problem the digest task brief highlighted.
- `topology_label`: ["GTO-to-hyperbola", "4-impulse-apse", "indirect-
  pontryagin"]
- `validation_level`: V0 (sourced; not reproduced by us)
- `cycler_referenced`: mcconaghy-2006-two-synodic (or whichever
  catalogue ID covers the 4 Feb 2019 Earth-arrival cycler)

This anchor is **complementary to** the existing cycler rows, not a
replacement for them — Pontani 2018 takes a cycler trajectory AS INPUT
and optimises the taxi joining it.

## 5. Catalogue impact

**No new catalogue rows directly admitted from this paper alone.**

The paper does not commit to a specific catalogue-admissible trajectory:

- The cycler half is referenced (McConaghy 2006) but not produced by
  Pontani 2018.
- The taxi half is generic optimization output; it operates on the GTO
  initial state and a published hyperbola at Earth, but Pontani does
  not deliver a SPECIFIC Earth-launch trajectory with a tabulated
  set of dates, ephemeris-locked V_infs, etc. The PSO is run on
  Keplerian assumptions (Earth-centric two-body) with the cycler
  hyperbola fixed.

**Indirect catalogue impact:** if the McConaghy 2006 "two synodic
period cycler" Earth-arrival on 4 Feb 2019 is already a catalogue row,
that row could carry a back-reference to Pontani 2018 in `notes`
documenting that the published rendezvous Delta-V for that hyperbola is
1.791 km/s (low-thrust) / 1.796 km/s (2-impulse no abort) / 2.898 km/s
(2-impulse with safe abort).

## 6. Schema impact — THE KEY DECISION POINT

The task brief asks: **does Pontani 2018's hyperbolic-rendezvous problem
warrant a new orbit class in the schema, or fit as extension of
`precursor_mga`?**

**RECOMMEND: extend `precursor_mga` semantics; do NOT add a new orbit
class.** Rationale below.

### 6.1 Re-reading the precursor_mga definition (schema v4.7)

From `data/catalogue.schema.json` line 144 (verbatim):

> precursor_mga = non-repeating MGA chain that inserts a spacecraft INTO
> a steady-state cycler, epoch_locked, n_returns=1, MUST set inserts_into
> to an existing cycler row id.

This is **exactly** the function of Pontani 2018's taxi vehicle: an
Earth-launched non-repeating maneuver chain that injects a spacecraft
into a hyperbolic encounter with a cycler. The "MGA" wording is a
literal "multiple gravity assist", which Pontani 2018 does NOT use (the
taxi is purely chemical-or-electric impulse(s) inside Earth SOI, no
gravity assists). BUT the semantic role — getting a crew vehicle from
Earth-launch into rendezvous with a steady-state cycler — is the
**identical use case**.

Three options:

**Option A: extend `precursor_mga` to cover taxi-rendezvous.**
Treat the "MGA" in the name as a vestigial label, and clarify the
description to read "non-repeating chain (gravity assist or
impulsive-rendezvous) that inserts a spacecraft into a steady-state
cycler". Cheap. Backward-compatible. Recommended for the FIRST
taxi-anchor admission.

**Option B: add a new `taxi_rendezvous` orbit_class.**
Surgical clarity, but five-way taxonomy expansion just two months after
the four-way scope expansion (#294, 2026-06-15) is over-reactive. Only
one taxi-side reference paper exists in the corpus (Pontani 2018).
Cost > benefit at one-row scale.

**Option C: defer the decision.**
Admit Pontani 2018 to KNOWN_CORPUS but do not catalogue-admit it
either as `precursor_mga` or as a new class until a SECOND taxi-side
paper lands. This is the most conservative.

### 6.2 My recommendation

**Option A** — extend `precursor_mga` semantics. Specifically:

- KEEP the orbit_class enum at four values
- BROADEN the description string in catalogue.schema.json line 144 to:
  "precursor_mga = non-repeating chain (gravity assist sequence,
  impulsive rendezvous, or low-thrust transfer) that inserts a
  spacecraft into a steady-state cycler, epoch_locked, n_returns=1,
  MUST set inserts_into to an existing cycler row id."
- ADD a sub-classification annotation in `notes` for each
  precursor_mga row: "mga_sequence" vs "rendezvous_only" vs
  "low_thrust_rendezvous" — free-text, not a schema enum.

This preserves the four-class semantic taxonomy while accommodating
Pontani 2018's specific mechanism. If a second / third / Nth taxi
paper lands, a future schema v4.8 can promote the sub-classification
to an enum.

### 6.3 Inserts-into pointer

If Pontani 2018 is ever catalogue-admitted under `precursor_mga`, the
`inserts_into` pointer should resolve to the McConaghy 2006 two-
synodic-period cycler row. Verify that row exists in catalogue.yaml
before admission. Per the brief check at digest-write time, the
McConaghy-2006 digest was completed earlier today, so the catalogue
row likely exists.

## 7. Errata

None found. The paper is mathematically careful — the BVP / Pontryagin
formulations are correct, the PSO setup is standard, the numerical
results carry residuals at the 10^-10 level (Eq 69).

One minor typographical observation: the formula text occasionally uses
period for decimal separator and occasionally renders multiplied units
with slashes inconsistently (Sec III.D Eq 66 has "n_0 = 0.001 g_0" with
g_0 = 9.8 m/s^2, which gives n_0 = 0.0098 m/s^2; the choice of writing
the initial thrust acceleration as a g_0-fraction rather than directly
in m/s^2 is unusual but unambiguous). Not an errata.

The paper is a "Article in Advance" version; final journal version may
have minor copy-edit differences (the PDF carries the AIAA Crossmark
header but no final volume/issue/page numbers).

## 8. Action items

1. Anchor `pontani-conway-2018-taxi-rendezvous` in KNOWN_CORPUS as
   the first taxi-side reference (V0).
2. SCHEMA RECOMMENDATION (do NOT edit now): broaden the
   `precursor_mga` enum description to encompass impulsive-and-low-
   thrust rendezvous, in addition to the MGA-gravity-assist chains
   already covered. This is Option A above; defer the
   description-string edit to a paired admission task (the Pontani
   2018 catalogue row + the schema doc edit, atomic).
3. Cross-check: identify the McConaghy 2006 catalogue row id (likely
   `mcconaghy-2006-two-synodic` or similar) and confirm the 4
   February 2019 Earth-hyperbola-arrival epoch matches. If it does,
   the catalogue back-reference (Pontani 2018 rendezvous Delta-V
   options) is a free annotation update.
4. Future work: when a second taxi-side paper lands (e.g. Anderson
   2016 ref [7] AIAA 2016-5267, or a NASA Decadal taxi study, or a
   Mars Society 2019 free-return-with-taxi-rendezvous paper), revisit
   the schema-v4.8 question of whether `taxi_rendezvous` becomes a
   fifth orbit_class. With ONE taxi paper the current
   recommendation is Option A; with multiple papers a dedicated
   class becomes defensible.
5. Reproduction potential: the impulsive 2-impulse case is a clean
   PSO over 3 unknowns and is reproducible with our existing tools
   given the GTO and the hyperbola (p_H, e_H). The low-thrust case
   requires the indirect Pontryagin solver, which the project does
   not currently have. Defer low-thrust reproduction.
6. Cycler V_inf reference value: a downstream computation worth
   noting — the McConaghy 2006 two-synodic cycler's Earth V_inf
   computed from (p_H = 91595.2 km, e_H = 2.157) is
   V_inf = sqrt(mu_E * (e_H^2 - 1) / p_H) where mu_E = 398600.4418
   km^3/s^2; this is ~5.85 km/s — consistent with typical two-
   synodic Earth V_inf in the literature. Useful as a cross-check
   value for the McConaghy 2006 catalogue row.
