# Mining Note: Shakouri, Kiani & Pourtakdoust (2019), SBGM multiple-impulse coplanar maneuvers

**Source:** Amir Shakouri, Maryam Kiani, Seid H. Pourtakdoust, "A New Shape-Based
Multiple-Impulse Strategy for Coplanar Orbital Maneuvers," arXiv:1905.04543v1 [math.OC],
11 May 2019 (preprint submitted to Acta Astronautica).

**Mined:** 2026-06-10 from full preprint (26 pp., 40 refs).

**Mining lens:** (1) shape-based parameterisation as a SEED generator for our
multi-impulse/DSM lane (`search/dsm_leg.py` Takao η-transcription + MBH); (2) published
numeric transfers usable as independent goldens; (3) validity limits vs our eccentric
coplanar Earth–Mars legs.

**Verdict: USEFUL — as a golden source and as a closed-form 2-impulse seed primitive;
NOT a drop-in seeder for our time-anchored DSM legs (no time-of-flight control).**

---

## 1. The method (SBGM) in one paragraph

The transfer trajectory is a chain of **confocal elliptic arcs** (all with the central
body at a focus, Assumption 2.3) joined so adjacent arcs intersect AND share a tangent
direction — G0 + G1 continuity ("smooth", Definition 2.1, p.4). Because the velocity
direction is continuous across each junction, **every impulse is tangent by
construction**; only the speed jumps. The geometry is two-body, coplanar, elliptic-only.
Hohmann (Proposition 3.2, p.10) and bi-elliptic (Proposition 3.3, p.11) transfers are
proved to be special cases.

## 2. Parameterisation, free parameters, equation map

- Planar arc description: `p = [a, e, θ, ω]ᵀ`, polar conic
  `r = a(1−e²)/(1+e cos(θ+ω))` — **Eq. (1), p.4**; ω is the apse-line rotation in the
  inertial frame, θ+ω the true anomaly (Fig. 1, p.5).
- Element extraction from (r, v): **Eqs. (2)–(5), pp.4–5**; tangent slope
  `dr/dθ = e sin(θ+ω) r / (1+e cos(θ+ω))` — **Eq. (6), p.5**.
- Junction conditions at intersection i (same radius, same polar angle, same slope):
  **Eqs. (7)–(9), p.6**, reduced to two scalar residuals per junction:
  - `f_{i1}`: **Eq. (10), p.6** (tangency, in e_i, e_{i+1}, ω_i, ω_{i+1}, θ_{i(i+1)})
  - `f_{i2}`: **Eq. (11), p.6** (radius match, brings in a_i, a_{i+1})
- **Counting (the headline result):** an N-impulse maneuver has N−1 designed arcs, N
  junctions ⇒ 2N equations vs 4N−5 unknowns {a_i, e_i, ω_i (i=2..N), θ_{i(i+1)}
  (i=2..N−1)}, so **M = 2N−5 adjustable parameters ε** (p.6). Table 1 (p.8) compares
  against 3N−6 … 3N−1 design variables for conventional fixed/free-endpoint impulsive
  formulations. For N=3 the entire transfer has **one** free parameter (they use ω₂,
  p.12); for N=2 it has **zero** (fully determined by the smoothness equations).
- Algorithm 1 (p.8): assemble f = [f_{11} … f_{N2}]ᵀ, solve by relaxed Newton (or GA).
  **Remark 2.6 (p.7):** the Jacobian is singular if any intermediate arc is circular —
  never seed Newton with e_i = 0 (Proposition 2.5, p.6: intermediate arcs can only be
  circular if impulses sit exactly at apses).
- **Impulse recovery — Eq. (17), p.9:** because impulses are tangent,
  `Δv_i = v_i⁺ − v_i⁻` with vis-viva speeds
  `v⁻ = sqrt(μ(2/r_i − 1/a_i))`, `v⁺ = sqrt(μ(2/r_i − 1/a_{i+1}))`. Scalar subtraction —
  no flight-path-angle bookkeeping needed.
- **Closed-form 2-impulse special case (Remark 2.7, p.7):** first impulse at perigee of
  the initial orbit ⇒ system collapses to **Eqs. (14)–(15), p.7** in (e₂, θ₂₃); if the
  target is circular it is fully closed-form: `θ₂₃ = π`,
  `e₂ = (a₃ − a₁(1−e₁))/(a₃ + a₁(1−e₁))` — **Eq. (16), p.7**. This is a generalized
  Hohmann from an eccentric departure orbit, recovering Eq. (25) (p.11) when e₁ = 0.
- Section 5 (pp.20–21) is a small-variation/continuous-thrust extension
  (**Eqs. (28)–(29), p.20**) — tangential low-thrust shaping, not relevant to our
  impulsive lane.

## 3. Published numeric examples (golden candidates)

Both case studies are Earth-centred two-body, coplanar. μ is **not stated** in the
paper; the element-extraction reference is Curtis [their ref. 40], so the standard
μ_Earth ≈ 398 600 km³/s² is the presumed value — a golden test must confirm which μ
reproduces the tabulated numbers before locking the tolerance (sourced-expected rule
still holds: the EXPECTED side is the published table; μ is a documented model input).

### Case study 1 — eccentric LEO → circular LEO (Sec. 4.1, p.13; Table 2, p.14)

Inputs (p.13):

| | a (km) | e | ω | departure/arrival anomaly |
|---|---|---|---|---|
| initial (i=1) | 13 756 | 0.5 | 10° | θ₁₂ = 270° |
| final (i=4) | 13 756 | 0 | 60° | θ₃₄ = 30° |

Outputs (Table 2, p.14; all costs km/s, times s):

| Method | J_c (sum ‖Δv‖) | J_m (max ‖Δv‖) | t_f (s) |
|---|---|---|---|
| 1-impulse | 2.6305 | — | 2631 or 3463 (two solutions) |
| 2-imp Lambert (a) (t_f free) | 4.4539 | 2.2989 | 3750 / 3184 |
| 2-imp Lambert (b) (t_f, θ_f free) | 1.4677 | 0.7831 | 11 640 / 12 090 |
| 2-imp perigee-start (Remark 2.7, closed-form) | 1.5210 | 0.9878 | 2315 |
| 2-impulse SBGM | 1.5746 | 0.9487 | 25 415 |
| 3-impulse SBGM (opt. ω₂) | J_c* = 1.5746 | J_m* = 0.9471 | 24 581 / 23 156 |

### Case study 2 — near-circular LEO → Molniya (Sec. 4.2, p.16; Table 3, p.18)

Inputs (p.16):

| | a (km) | e | ω | anomaly |
|---|---|---|---|---|
| initial (i=1) | 6 644.4 | 0.01 | 60° | θ₁₂ = 45° |
| final (i=4) | 26 562 | 0.74105 | 30° | θ₃₄ = 15° |

Outputs (Table 3, p.18):

| Method | J_c | J_m | t_f (s) |
|---|---|---|---|
| 2-imp Lambert (a) | 5.1176 | 7.9455 | 2894 / 3570 |
| 2-imp Lambert (b) | 1.3344 | 2.5604 | 826 / 764 |
| 2-impulse SBGM | 2.3263 | 2.5659 | 5009 |
| 3-impulse SBGM (opt. ω₂) | J_c* = 1.3815 | J_m* = 2.5659 | 4560 / 5009 |

Golden value: case study 2 exercises an **eccentric (e = 0.74) oblique target** —
exactly the terminal class of our heliocentric Earth–Mars legs — and the 2-impulse SBGM
row needs **no optimizer at all** (zero adjustable parameters), so it is the cleanest
independently-reproducible row. The Remark 2.7 closed-form row of case study 1
(J_c = 1.5210, J_m = 0.9878, t_f = 2315 s) is reproducible from Eq. (16) by hand.
Hohmann/bi-elliptic special-case propositions give additional analytic sanity anchors.

## 4. Limits of validity (precise)

1. **Two-body, single focus** (Assumption 2.3, p.4). Fine for our heliocentric
   patched-conic legs; the focus becomes the Sun.
2. **Coplanar only** — the whole formulation is polar-planar (Eq. 1). Matches our
   coplanar Earth–Mars leg lane; useless for inclined legs.
3. **Elliptic arcs only** (Assumption 2.3). No hyperbolic or parabolic transfer arcs.
   Heliocentric Earth–Mars arcs are elliptic, so OK; it could NOT shape an escape leg.
4. **Eccentric oblique terminals supported** — this is the paper's selling point
   (abstract; case study 2 e = 0.74105). Our eccentric coplanar terminals are in scope.
5. **Tangent impulses only.** G1 continuity forces every interior AND terminal impulse
   to be velocity-aligned. The reachable transfer family is a strict subset of general
   multi-impulse transfers; their own Table 2/3 discussion (item 6, p.20) concedes the
   unconstrained-direction 2-impulse solutions (Lambert b, Remark 2.7) beat SBGM on
   cost. Our Takao η-DSM lane allows arbitrary impulse vectors — SBGM seeds would land
   in the tangent-impulse sub-basin only.
6. **No time-of-flight control — the critical mismatch.** Maneuver time t_f is an
   *output* of the geometry (Tables 2–3 report it; it is never constrained). SBGM is an
   orbit-to-orbit method with free phase. Our DSM legs are **point-to-point,
   epoch-anchored** (Lambert from DSM position to a target position over a fixed
   remaining ToF): the planet must be at the arrival point at the arrival epoch. SBGM
   has no mechanism to hit a (position, time) pair, so it cannot directly seed a
   time-anchored leg or any cycler phasing problem.
7. **Numerical hygiene constraint:** Newton solution must avoid e_i = 0 intermediate
   arcs (Remark 2.6, p.7) — relevant if we ever implement it: clamp e seeds away from
   zero, the same pattern as our η-endpoint clamping in `dsm_leg.py`.

## 5. Seeding-utility assessment for our solver

- **Direct seeding of `evaluate_dsm_chain`: NO.** Our genome coordinates per leg are
  (V∞, α, β, τ, η) with fixed epochs at planetary encounters. SBGM produces
  free-time orbit-to-orbit arc chains with tangent impulses; there is no map from its
  output to a leg whose ToF is pinned by planetary ephemerides. Bolting a time
  constraint onto SBGM destroys its variable-count advantage (the 2N−5 economy comes
  precisely from leaving time free).
- **Indirect use 1 — cheap geometric warm start for the DSM interior.** For a leg whose
  endpoints are *nearly* apse-aligned, the closed-form generalized-Hohmann Eq. (16) /
  Remark 2.7 system (Eqs. 14–15) gives, at negligible cost, an intermediate-ellipse
  (a₂, e₂, ω₂) whose junction anomaly θ₂₃ locates a near-tangent DSM point. Converting
  that junction's (position, epoch-fraction) into an (η, τ) pair would seed MBH inside
  the tangent-impulse basin — plausibly better than descriptor-based seeding for
  low-energy legs where the optimum is near-tangential, but unproven, and only worth
  trying if basin statistics show DSM solutions clustering near tangency.
- **Indirect use 2 — basin pruning.** The N=3 one-parameter family (sweep ω₂, Fig. 2
  p.13 / Fig. 6 p.17) is a 1-D scan that exposes where two-impulse degenerations sit
  (|Δv₁| → 0 or |Δv₃| → 0 points). As a pre-scan it is far cheaper than MBH restarts,
  but again only inside the free-time tangent family.
- **Primary value: independent goldens.** Two complete input/output case studies
  (Sec. 3 above) plus closed-form Eq. (16) give a sourced expected side for a
  multi-impulse coplanar solver test — something our DSM lane currently lacks at the
  2-/3-impulse two-body level. Per repo rule, the expected values trace to Tables 2–3
  of the publication, not to our own code.

## 6. What this does NOT give us

- No heliocentric or Earth–Mars numeric example (both case studies are Earth-orbital).
- No V∞ / flyby concept anywhere — nothing to check free-return anchors against.
- No time-fixed or rendezvous variant; phasing is explicitly out of scope.
- No 3D extension; no perturbations.
- μ used in the case studies is unstated (must be pinned down empirically before a
  golden test tolerance is set; record the resolved value in the test docstring).
