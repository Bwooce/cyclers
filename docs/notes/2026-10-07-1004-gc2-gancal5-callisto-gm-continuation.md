# #1004: gc-2 and GanCal#5 in continuous gravity, a continuation in Callisto's GM

Status: PRE-REGISTRATION (written and committed before any #1004 run). No catalogue writes. Owner
ruling to respect: gc-2 is a "GanCal-family relative", NOT novel; this work refines that relation,
it does not reopen novelty.

## 1. Question and why the continuation starts from gc-2

In the patched conic no Callisto-mass path joins gc-2 and GanCal#5 (generator note 6.20-6.21: the
date residual contains no GM; gc-2 needs at least 0.214 of Callisto's GM for its 200 km floor). In a
continuous-gravity model the number of encounters is not fixed, so the open test is whether gc-2's
branch, continued in Callisto's GM, joins GanCal#5.

GanCal#5 cannot be the starting point: it passes through Callisto's centre (Callisto massless in
R-S), so for any s_C > 0 it is a collision orbit. The continuation starts from gc-2 at full mass and
goes down in s_C (Callisto GM scale).

## 2. Model and the comparison endpoint

- The #968 lane model, validated on GanCal#5 in the ideal continuous model (#968 note sec. 7):
  R-S 2009 Table 2 circular coplanar ephemeris (`run_942_enumerate.rs_moon_system`), Jupiter central
  GM from the lane, Ganymede AND Callisto point masses with the Jupiter-frame indirect term
  (`jovian_stm.propagate_with_stm`, DOP853 + analytic STM; GM/radius overrides, f5662f1a).
  Ganymede stays at full mass throughout.
- The endpoint "GanCal#5" in this model is the #968 orbit, NOT the published numbers:
  `data/968_control/a_state.json`, Ganymede V_inf 3.0403 km/s, Ganymede altitude 576.5 km, Callisto
  speed 3.2513 km/s at its centre pass (s_C = 0). Comparing with 3.24 / 328 km would repeat the
  criterion-4 model mismatch of #968.
- Reference conic elements: GanCal#5's G->C leg a 1,683,097 km, e 0.4205; gc-2's C-C leg
  a 1,686,915 km, e 0.3856 (patched conic, note 6.21).

## 3. Method (`scripts/run_1004_gc2.py`)

- Corrector: forward-backward multiple shooting as #968 corrector A, generalised to four massive-body
  nodes in time order B (Ganymede, 11.68 d), C1 (Callisto, 14.29 d), C2 (Callisto, 35.80 d),
  A (Ganymede, 38.41 d), and B' = (Q x_B, t_B + T), T = 3 S_GC = 37.5697 d, Q the rotation by the
  common advance (90.42 deg). Each node: state and epoch (7 unknowns), periapsis gauge relative to
  its own moon. 28 unknowns, 24 match rows + 4 gauges = 28 residuals. With two moons at different
  rates the problem is non-autonomous (no Jacobi integral, no time-shift symmetry), so the system
  is square and should be full rank.
- Seed: the patched-conic gc-2 (`data/943_cell_gc_gauntlet.json`, cell "gc", both moons massive),
  nodes from `periapsis_node` at each flyby.
- Stage 1 (identity of the s_C = 1 orbit): joint GM scale sigma on both moons (radii scaled by
  sigma) from sigma = 0.05 up to 1. At sigma = 0.05 the orbit must be near the patched-conic gc-2
  (V_inf G within 0.02 of 3.617, C within 0.02 of 3.039); the branch is followed to sigma = 1. No
  agreement with the patched-conic numbers is required at sigma = 1.
- Stage 2 (the question): pseudo-arclength continuation in (z, s_C) from the sigma = 1 orbit, with
  s_C decreasing, Ganymede at full mass. Variables scaled by the column norms at the start; tangent =
  the null vector of the column-scaled [J_z | J_s]; step from a target Delta s_C. At each point a
  clear gap between the smallest and second-smallest singular values of [J_z | J_s] is required
  (ratio > 10), else the step is rejected and halved.
- Two variants of Stage 2:
  - (S) scaled radius: Callisto softening radius s_C x 2410.3 km (point-mass limit; the mathematical
    branch).
  - (P) physical radius 2410.3 km (registry). Every accepted point must have both Callisto periapses
    above 2408 km (R-S radius); the run stops at the first point where one is not, and records where
    the 200 km floor (2608 km) is crossed. No converged point inside the softened core is accepted.
- Recorded per point: s_C; Ganymede V_inf at B and A; Callisto V_inf, turn, periapsis and
  periapsis / SOI (SOI 37,681 km) at C1 and C2; a and e of the C1->C2 arc at its mid-time; the
  smallest two singular values of J_z and of [J_z | J_s]; the sign of d s_C / d l.

## 4. Stopping rules

- s_C <= 1e-3; or
- 6 step halvings without convergence: a NUMERICAL stop (stiffness grows about as 1/s_C), not a
  physical end; or
- after a fold: 60 more points, or s_C back at 1; or
- variant (P): a Callisto periapsis at or below 2408 km;
- each call under 8 minutes, checkpointed in `data/1004_gc2/`; live logs only in
  `data/1004_gc2/live/` (gitignored).

## 5. Fold versus limit degeneracy (fixed now)

- Fold: sigma_min(J_z) -> 0 at finite s_C with the sign of d s_C / d l flipping.
- Limit degeneracy: sigma_min(J_z) shrinking roughly in proportion to s_C as s_C -> 0, no sign flip.

## 6. What "a continuous path to GanCal#5" would look like (fixed now)

Along the gc-2 branch, or a branch reached through a fold, at s_C -> 0:
- one Callisto encounter's periapsis moves out past the SOI as its turn goes to zero;
- the other becomes a centre pass with zero turn;
- Ganymede V_inf moves from about 3.6 to the #968 value 3.04 (within 0.05 km/s), and the long leg's
  (a, e) moves to GanCal#5's.
Anything else is "no path".

## 7. Expected outcome (pre-registered)

NO PATH. Reason: the patched-conic date residual contains no GM, so with the radius scaled by s_C
the gc-2 branch should reach s_C -> 0 with both Callisto encounters still near 6.9 deg each and
periapses proportional to s_C; its limit is gc-2's own generating orbit. Only a fold could change
that, and a fold would be the surprising outcome. Variant (P), patched-conic estimates (r_p at
s_C = 1 about 12,208 km, R-S radius 2408): surface reached near s_C = 0.197, 200 km floor near
0.214 (matching note 6.21); the continuous values will differ and are recorded.

## 8. Results

(pending)

### 8.1 Stage 1 so far: a fold in the joint mass scale (2026-10-08)

- Natural continuation in sigma from 0.05 (`data/1004_gc2/sigma.json`) converged at sigma = 0.05 with
  V_inf G 3.6042, C 3.0248 (patched conic 3.617 / 3.039: within the 0.02 of sec. 3, identity OK),
  then at 0.0625 ... 0.1615. Callisto r_p / sigma fell 12,073 -> 11,320 km at an increasing rate;
  the step to 0.1696 failed.
- Pseudo-arclength from 0.1615 (`cont_J.json`): sigma rose to 0.16395, then the tangent's sigma
  component changed sign and the branch came back DOWN in sigma: 0.15953, 0.147, ..., 0.0571, with
  Callisto r_p / sigma 10,996 -> 10,224 km, Callisto turns 8.0 -> 8.9 deg, V_inf C 2.944 -> 2.892.
  At sigma about 0.098 the two branches differ: V_inf C 3.009 vs 2.906 km/s, r_p / sigma 11,875 vs
  10,446 km. Reading (not yet verified): the gc-2 branch from the patched-conic limit FOLDS at sigma
  about 0.164 and does not reach the physical masses along this path.
- Rule violation, disclosed: sec. 3 required a singular-value gap > 10 in [J_z | J_s] before each
  step; the code did not enforce it, and the logged gap was about 1.0-1.1 at every step. J_z
  (column-scaled) has four singular values of 6e-7 to 2.6e-6, against 3e-4 for the next. The
  converged points are solutions (each meets the floors), but the tangent was not well defined by
  the registered test.

### 8.2 AMENDMENT 1 (before the runs it covers)

The fold claim must not rest on the tangent. Checks, each with its pass rule:
1. Two solutions at one sigma: at sigma = 0.150 and 0.160, damped Newton at FIXED sigma from the
   nearest lower-branch point and from the nearest upper-branch point. Fold supported if both
   converge to the floors and differ by more than 0.02 km/s in V_inf C and 100 km in Callisto
   r_p / sigma.
2. Fold bracket: natural continuation on the lower branch from 0.1615 upward in steps of 0.0005 (then
   0.0001 near failure); record the last converged sigma and that steps down to 1e-5 fail beyond it.
3. Independent cross-check: IAS15 re-fly of the eight half-arcs of each of the two sigma = 0.150
   solutions (agreement < 1e-2 km, 1e-7 km/s).
4. Direct attempt at the physical masses: damped Newton at sigma = 1 from the patched-conic seed and
   from each branch's last point (offsets scaled); record convergence and, if it converges, the
   geometry (a different orbit is reported as such, not as gc-2).
5. The weak J_z directions: recompute J_z at rtol 1e-13 at one point; if the four small singular
   values move by more than a factor 3, they are integration noise; otherwise structural.
