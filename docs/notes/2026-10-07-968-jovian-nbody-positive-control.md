# #968: Jovian n-body lane positive control (GanCal#5), and the periodicity-wrap fix

Status: PRE-REGISTRATION (sec. 3) committed before any control run. No catalogue writes. Nothing here
is called novel.

## 1. Why the lane had no control (history read before designing)

- Member D (`#223`, `docs/notes/2026-06-13-liang-member-d-nbody.md`): `shoot_cycle` (no periodicity
  wrap; nodes at the patched-conic flyby periapsis; FD Jacobian, 60 nfev; a boundary V_inf pin).
  Position continuity reached 1-11 km, velocity stalled at 135-216 m/s. The note's own sec. 5 names
  the periapsis-node seed as the source of the 1e4-1e6 km seed gap.
- EGGIE (`#480`, notes `2026-06-29-480-eggie-*` and `2026-06-30-480-eggie-level3-nbody.md`):
  `ideal_eggie_shoot` / `jovian_shoot` (FD, analytic STM, epoch-free and sub-arc variants) all
  plateaued at 0.1-0.4 km/s velocity continuity, wrap block about 0.19 km/s.

## 2. Shared-code defect found and fixed (lead ruling 2026-10-07)

The periodicity wrap in `jovian_defect_residual` and `subarc_defect_residual` compared the end node
and the start node relative to the home moon by translation only. A cycler repeats rotated by the
home moon's advance over the period, so that wrap can be met only when the advance is 0 mod 360 deg.

- Evidence (computed from `resonant_conic.ideal_moon_smas` / `ideal_t_syn`, lane mu): over the EGGIE
  ideal cycle (4 T_syn = 28.02 d) Ganymede advances 0.00 deg, Europa 339.50 deg (-20.50) and Io
  298.49 deg (-61.51). The home moon is Europa, so the old wrap asked for an orbit that repeats
  un-rotated while Europa moved 20.5 deg. Also the three moons do not advance by a common angle, so the
  ideal three-moon configuration does not repeat at all: no exact periodic orbit exists in that model.
- Test first: `8d79b3fa` (`tests/nbody/test_968_wrap_rotation.py`; fixture exactly periodic in the
  moon's rotating frame with a 120-deg advance, coplanar and 25-deg inclined; the old wrap leaves
  1.33e6 km and 16.6 km/s, the legs close to 2.6e-7).
- Fix: `b9c27f77`. `wrap_rotation()` maps the moon's instantaneous orbital frame
  (r_hat, h_hat x r_hat, h_hat) at t0 onto the one at tn; the wrap is `reln - Q rel0`. Circular
  coplanar model: rotation about the normal by the advance. Real ephemeris: the actual advance in the
  moon's orbital plane. The STM Jacobians carry `-R_W blockdiag(Q, Q)` on node 0. Q is the identity
  when the advance is 0 mod 360 deg, so such results do not change. tests/nbody (99) pass.

### 2.1 Results computed with the old wrap: now VOID

A bug fix voids past negatives produced by the buggy code. VOID (the lead registers the re-run):
- `#480` EGGIE Stage 2 (`2026-06-29-480-eggie-stage2-nbody-verdict.md`), Stage 3 STM
  (`...-stage3-stm-verdict.md`), Stage 4 sub-arc (`...-stage4-subarc-verdict.md`), the level-3
  real-ephemeris run and its CORRECTION (`2026-06-30-480-eggie-level3-nbody.md`) and
  `2026-06-30-480-eggie-realeph-stm.md`: every plateau there was measured against the old wrap. The
  skip reason in `tests/verify/test_ieg_reproduction_golden.py` ("a basin/model feature, not a
  numerics one") rests on these runs. Note the ideal three-moon model has no exact periodic orbit
  at all (sec. 2), so a re-run needs a criterion other than exact periodicity there.
- `#318` (`scripts/scan_318_joint_sobol_smoke.py`) and `#501`
  (`scripts/scan_501_broadened_joint_search.py`): their n-body stage is `jovian_shoot(jacobian="stm")`,
  so "26 shot, 0 closed" is VOID, and so is the n-body part of the empty-region stamps
  `jovian-cgcec-sobol-smoke-318-2026-06-30` and the six `jovian-*-sobol-broadened-501-2026-06-30`
  rows in `data/empty_regions.jsonl`. Their "positive control (Liang Member D) PASSED" was the
  patched-conic prefilter only (`scan_501_broadened_joint_search.py` lines 772-785); the n-body
  stage never had a control.
- NOT affected: Member D (`#223`). `shoot_cycle` has no wrap (`_cycle_residual` blocks are leg
  continuity, encounter hinges and a boundary V_inf pin). Its failure has other causes, named in
  sec. 1.

## 3. PRE-REGISTRATION (written and committed before any control run)

### 3.1 Control and model

- Control: Russell & Strange 2009 GanCal#5 (published V_inf G 3.24 / C 3.34 km/s, minimum Ganymede
  altitude 328 km, r from 821,915 to 2,390,844 km, period 3 G-C synodic periods = 37.5697 d). The
  patched-conic member is the #943 production enumerator's LITERAL recall
  (`data/943_cell_gc_gauntlet.json`, key `k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s`),
  rebuilt by `scripts/run_942_enumerate.py:rs_moon_system` (R-S 2009 Table 2 constants).
- Model (rung a): the published model made continuous. Ganymede and Callisto on circular coplanar
  rails (R-S periods 618,153 s and 1,441,931 s; radii from R-S mu = 126,686,535); Jupiter central GM
  from the lane (126,686,534, 8e-9 relative difference, recorded); Ganymede a point mass (lane GM,
  9887.834; R-S 9887.83) with the Jupiter-frame indirect term (the lane's force); Callisto MASSLESS,
  as in R-S (craft passes its centre). Altitudes use the R-S Ganymede radius 2634 km.
- What a pass on rung (a) validates: the propagators and correctors of the lane in a continuous
  model. It does NOT validate the jup365 rails path. Rung (b) (sec. 3.6) is the real-ephemeris step.

### 3.2 Corrector A (independent; `scripts/run_968_control.py`)

Forward-backward multiple shooting. Nodes: Ganymede flyby A (state 6 + epoch), Ganymede flyby B
(6 + epoch), Callisto node C (velocity 3 + epoch; position pinned to Callisto at the node epoch, the
massless-target hit). The end of the last leg is A's periodic image: state `Q x_A` at `t_A + T`,
`Q` the rotation by Ganymede's advance over T (90.42 deg mod 360 for GanCal#5, the same for Callisto since T is 3 synodic periods). Each leg is
matched at its mid-time by a forward arc from its start node and a backward arc from its end node.
Gauge: (r - r_G).(v - v_G) = 0 at A and B (nodes at the Ganymede periapsis). Count: 18 unknowns,
20 residuals, one cyclic Jacobi-constant redundancy expected; zero residual is achievable. Arcs:
`jovian_stm.propagate_with_stm` (DOP853 + analytic STM, rtol 1e-12, atol 1e-10); Jacobian analytic
(STM plus epoch columns), checked once against finite differences before use. Seed: the patched-conic
member, A and B from `periapsis_node` geometry (moon-centred hyperbola), C at Callisto with
V_inf_in = V_inf_out. Velocity rows weighted 1e3 (lane convention).

### 3.3 Corrector B (the lane, fixed wrap)

`jovian_defect_residual` (REBOUND IAS15 on rails) with `jovian_stm_jacobian`, via `least_squares`
(trf), nodes at the patched-conic encounters (A, B, C, A'), epochs fixed (the lane's documented
choice), moons = ("Ganymede",), the ideal R-S ephemeris. The lane has no Callisto-hit residual;
its C node may slide. Lane PASS is judged on its own floors and on agreement with corrector A's
orbit (sec. 3.4 (5)).

### 3.4 PASS criteria (all required, on each corrector where stated)

1. Continuity: every match or leg defect below the lane floors, 1e-3 km and 1e-6 km/s (A: every
   mid-leg match and the wrapped end; B: `jovian_shoot`-style leg block below the floor and the
   rotated wrap below the same floors). Gauge residuals of A below 1e-9 (normalised).
2. Massless-target hit (A): C node within 1e-3 km of Callisto's centre (pinned by construction;
   checked). (B): distance of the C node to Callisto reported, not judged.
3. Sequence: the converged orbit makes exactly the encounters G, G, C per period: Ganymede closest
   approaches below Ganymede's Hill radius (31,715 km) only at A and B; altitude at A and B >= 100 km (project
   floor), on a 0.01-d sampling plus local minimisation.
4. Agreement with the published member: Ganymede V_inf (two-body energy relative to Ganymede at the
   periapsis) within 0.10 km/s of 3.24 at both flybys; Callisto relative speed within 0.10 km/s of
   3.34; minimum Ganymede altitude within 150 km of 328 km; r_min and r_max within 2 % of
   821,915 and 2,390,844 km.
5. Lane B agrees with A on the same physical orbit: Ganymede V_inf within 1e-3 km/s and periapsis
   altitude within 1 km at both flybys.
6. Independent cross-check of A's orbit with a second integrator: REBOUND IAS15
   (`JovianRestrictedNBody`, exact circular rails, `max_wall_sec` large; a wall-clock stop is
   reported as "timeout", never as a defect) re-flies every half-arc from its node: agreement with
   DOP853 at the match point < 1e-2 km and < 1e-7 km/s. The full-period IAS15 single shot from A is
   reported (descriptive; flyby amplification is large).

Expected outcome (stated now): corrector A converges; V_inf within 0.05 km/s and altitude within
100 km of the published values; corrector B with the fixed wrap converges to the same orbit.

### 3.5 If it fails

Diagnose which side failed with evidence: the seed (seed defects, basin), the lane (propagator
agreement, Jacobian against FD), or the criterion (e.g. a published value that is a patched-conic
display value). Two distinct, diagnosed approaches that fail go to the lead (stop and ask).

### 3.6 Rung (b), real ephemeris (after (a) passes; criteria amended here before it runs)

GanEur#316 at R-S's 2019 epoch on jup365 rails, a 1-2 cycle open chain seeded
from the #943 real-ephemeris patched-conic chain (`data/943_ganeur316_realeph/n10_rs2019_rel/`), no
periodicity wrap. The force model keeps R-S's massless target massless (Europa off, its hit
pinned as in 3.2); whether the other moons act as perturbers is fixed in the amendment. Criteria to be written into this note before the run.

## 4. Results

(pending)

### 4.1 First attempt (pre-registered corrector A, hand LM loop): plateau, not converged

- Seed match defects (km, km/s): leg B->C 287,856 / 1.94, leg C->A 8,066 / 0.069, leg A->B' 54,867 /
  0.391. Analytic Jacobian FD check at the seed: relative column errors 2e-7 to 1e-4.
- The hand LM loop (normal equations) stopped at |r| 164.75 (match dv 0.108 km/s, dr 0.85 km); a
  rerun with column-scaled augmented least squares reached |r| 98.6 (dv 0.097 km/s) and crawled.
  States kept: `data/968_control/a_state_normaleq_plateau.json`, `a_state_scaledlm_plateau.json`.
- Diagnosis so far: the stop is a stationary point that is not a root. The residual lies along the
  second-weakest scaled singular direction (sv 2e-6); the Gauss-Newton step there is 3e4 km and a
  line search raises the residual from alpha = 0.01. The FD check of all 18 columns at the plateau
  passes (worst 3e-3, the out-of-plane columns). The node Jacobi constants differ
  (B -172.692, C -172.629, A -172.606 km^2/s^2); Jacobi is conserved along each arc to 1e-11.
- Found while diagnosing: the published patched-conic GanCal#5 B->C leg passes GANYMEDE at 94,989 km
  (3.0 Hill radii) at t = 25.15 d with v_rel 2.70 km/s. The R-S model does not see this encounter;
  its impulse estimate 2 mu/(b v) is 0.077 km/s, the scale of the plateau. INFERRED link, not shown.
- A bug in my own script, fixed: the FD step for the gauge rows was chosen by column index mod 6,
  which is right only for node B. The gauge derivative is now analytic (it is a small row at weight
  1e4; not believed to be the plateau cause).

### 4.2 AMENDMENT 1 (before the runs it covers; criteria of sec. 3.4 unchanged)

- Solver for corrector A: scipy `least_squares` trf with the analytic Jacobian (a real trust region)
  replaces the hand LM loop. Restart from the patched-conic seed.
- Diagnostic run (does not replace the control): the same three-leg shooting with T FREE, the
  Callisto pin DROPPED (C a free node at a fixed epoch), and the Jacobi constant at B FIXED (dJ = 0
  first, from the plateau state). Outcomes and their reading, fixed now:
  - It converges to the floors: the corrector works; map T(J) over a small J range. If 3 S_GC is
    inside the range, a seed is interpolated and the pinned, fixed-T control is rerun (amendment 2).
    If T(J) turns back before 3 S_GC, there is no fixed-period continuation of GanCal#5 near the
    seed: a diagnosed negative on the published-orbit side, not on the lane.
  - It stalls near 0.1 km/s with T free: the corrector or formulation is at fault; go to approach 2,
    a continuation in Ganymede's GM from s_G about 1e-4 (the patched-conic limit) to 1, with the
    softening radius scaled.

### 4.3 T-free diagnostic result (amendment 1), and AMENDMENT 2 (before the runs it covers)

- From the plateau state, with T free, Callisto unpinned and J at B fixed (dJ = 0), trf converged:
  |r| 3.2e-6, match dr 5.6e-6 km, dv 2.8e-11 km/s, gauge 6e-14. So corrector A is NOT the cause of the
  plateau: it closes a Ganymede-only periodic orbit near GanCal#5 to far below the lane floors.
  That orbit has T = 37.518130 d, against 3 S_GC = 37.569705 d (short by 0.0516 d, 74 min), Jacobi
  -172.691866 km^2/s^2. State: `data/968_control/free_state_J0.json`.
- A J step of +0.02 stalled at J = -172.6746 (T 37.5184): the family may fold in J there.
- AMENDMENT 2: map the family by pseudo-arclength continuation (T and J both free; tangent from the
  smallest right singular vector of the column-scaled Jacobian; trf corrector with the arclength
  row), both directions from `free_state_J0.json`, recording T, J, V_inf and altitudes per point.
  Reading, fixed now: if T reaches 3 S_GC on the family, bracket it and rerun the pinned, fixed-T
  control from the bracketed point (amendment 3, criteria of 3.4 unchanged). If T turns back below
  3 S_GC in both directions within the explored range, GanCal#5 has no fixed-period continuation
  on this family: a diagnosed negative on the published-orbit side, with the T(J) curve as evidence.

### 4.4 Family result (amendment 2) and AMENDMENT 3 (before the runs it covers)

- The pseudo-arclength attempt took steps along an exactly null direction (the node-A slide; sv
  1e-16) and did not move; replaced by natural-parameter continuation in T (J free), recorded in
  `data/968_control/tcont.json`. Every step converged to |r| < 1e-5 (lane floors met by 2-3 orders):

| T (d) | J (km^2/s^2) | V_inf B / A (km/s) | alt B / A (km, R-S radius) |
|---|---|---|---|
| 37.518130 | -172.6919 | 3.1528 / 3.1532 | 536.7 / 573.8 |
| 37.5300 | -172.7658 | 3.1293 / 3.1293 | 556.5 / 562.9 |
| 37.5500 | -172.9046 | 3.0843 / 3.0843 | 568.3 / 568.3 |
| 37.569705 (= 3 S_GC) | -173.0387 | 3.0403 / 3.0403 | 576.5 / 576.5 |

- So a Ganymede-only periodic orbit with GanCal#5's topology and EXACTLY its period exists in
  continuous gravity, but with V_inf 3.040 km/s (published 3.24) and Ganymede altitude 576.5 km
  (published 328). Along this family dV_inf/dT is about -2.2 km/s per day; the member with
  V_inf = 3.24 would have T about 37.475 d. From T = 37.5325 d on, the two flybys are equal (a
  symmetric orbit).
- AMENDMENT 3 (criteria of 3.4 unchanged, and still judged as registered):
  1. Pinned control: seed corrector A (fixed T, Callisto pin) from the T = 3 S_GC member, time-shifted
     so that Callisto sits at the orbit's inbound crossing of Callisto's radius nearest the old C node
     (`--stage pinseed`), then `--stage a`, `--stage check`.
  2. Identity check (is this orbit the continuation of the published one?): continuation in
     Ganymede's GM, s_G from 1 down toward the patched-conic limit, at fixed T = 3 S_GC with the
     Callisto pin, softening radius scaled as s_G^(1/3) * R. Expected: V_inf tends to 3.24 and the
     altitude, scaled by s_G, to the patched-conic r_p as s_G -> 0. Needs a GM/radius override on
     `jovian_stm` (additive kwarg, pinned by a fast test; #1004 needs the same knob).

### 4.5 Pinned control result (amendment 3, step 1) and lane corrector B, first run; AMENDMENT 4

Corrector A, pinned layout, seeded per amendment 3 (`--stage pinseed`: inbound Callisto-radius
crossing at 35.9997 d, time shift -0.04134 d): converged in 11 evaluations. `data/968_control/a_check.json`:

| Criterion (3.4) | Value | Verdict |
|---|---|---|
| 1. continuity (floors 1e-3 km, 1e-6 km/s) | max match 3.2e-6 km, 1.6e-11 km/s; gauge 3e-13 | PASS |
| 2. Callisto hit | 0 km (pinned) | PASS |
| 3. sequence G, G, C; alt >= 100 km | Ganymede inside its Hill radius only at A and B; 576.5 km both | PASS |
| 4. V_inf G within 0.10 of 3.24 | 3.0403 (diff 0.200) | FAIL |
| 4. Callisto speed within 0.10 of 3.34 | 3.2513 (diff 0.089) | pass |
| 4. min Ganymede alt within 150 km of 328 | 576.5 (diff 248) | FAIL |
| 4. r_min within 2 % of 821,915 | 841,395 (+2.37 %) | FAIL |
| 4. r_max within 2 % of 2,390,844 | 2,376,185 (-0.61 %) | pass |
| 6. IAS15 re-fly of the six half-arcs | max 2.2e-5 km, 1.4e-10 km/s | PASS |
| 6. IAS15 full period from B (descriptive) | return miss 2.4e-3 km, 8.7e-7 km/s | (DOP853 single shot: 1.5 km) |

So, as registered, the control FAILS criterion 4: the continuous orbit with GanCal#5's topology,
period and Callisto hit exists and closes to the floors on two integrators, but it sits 0.20 km/s
below the printed V_inf, with the flyby 248 km higher.

Lane corrector B from the patched-conic periapsis seed (the Member D / EGGIE seed type): seed defect
5.4e5 km; after 40 evaluations still 2.4e5 km (leg dv up to 2.95 km/s); not converged
(`data/968_control/b_state_pcseed.json`). Forward node-to-node legs from periapsis nodes put the
24-day B->C leg 2.3e5 km off; this is the seed weakness the Member D note (sec. 5) already named.

AMENDMENT 4 (before running): (a) lane corrector B is also run from corrector A's converged orbit
(nodes and epochs from `a_state.json`): criterion 5 then asks whether the lane's own residual, with
the fixed wrap, holds that orbit at its floors and agrees with A; (b) the patched-conic-seed lane run
gets up to 3 more chunks of 40 evaluations; if still unconverged it is recorded as a seed failure.

### 4.6 Lane B from orbit A (amendment 4) and the GM continuation so far; AMENDMENT 5

- Lane B seeded from orbit A (`b_state_aseed.json`): the lane residual starts at 0.068 km /
  2.4e-5 km/s (its rails are a 0.02-d spline of Ganymede; A uses the exact circle) and converges to
  3.7e-6 km / 5.2e-9 km/s, wrap 3.8e-7 km. V_inf at both Ganymede nodes agrees with A to 6e-8 km/s and
  the node distances to 4e-5 km. Criterion 5 PASS. The old translation-only wrap could not hold this
  orbit (Ganymede advances 90.4 deg per period). From the patched-conic seed, 160 evaluations stall
  at leg dv 1.5-3.5 km/s (`b_state_pcseed.json`): a seed failure, as with Member D.
- GM continuation (amendment 3.2), `data/968_control/gm.json`: converged at s_G = 1, 0.95, 0.855,
  0.7695 (V_inf 3.0403, 3.0506, 3.0698, 3.0867 km/s), then failed at 0.73 and 0.6925. Cause: at
  s_G = 0.7695 the periapsis is 2430 km, against a softening radius of 0.7695^(1/3) x 2631 = 2412 km.
  The pre-registered s^(1/3) radius scaling shrinks slower than r_p, which scales with s_G, so the
  flyby runs into the softened core. My error in the pre-registration.
- AMENDMENT 5 (before the further runs): the softening radius scales linearly with s_G (r_surf =
  s_G x 2631 km), so r_p / r_surf is constant along the continuation, as in the patched-conic limit.
  Points already converged are unaffected (their periapses lie above both radii). Expectation
  unchanged: V_inf -> 3.24 and r_p / s_G -> about 2962 km as s_G -> 0.

### 4.7 Identity check result (amendments 3.2 and 5): the s_G = 1 orbit IS GanCal#5's continuation

`data/968_control/gm.json`: 35 converged points of the pinned, fixed-T problem from s_G = 1 down to
s_G = 0.036 (Ganymede GM scaled by s_G, softening radius by s_G after amendment 5). The run stopped
at s_G = 0.032 (|r| 4.6e-2 within 80 evaluations; the flyby time scale shrinks with s_G and the
arcs get stiff; not pursued further). Selected points (V_inf at both Ganymede flybys is equal):

| s_G | V_inf G (km/s) | r_p / s_G (km) | Callisto speed (km/s) |
|---|---|---|---|
| 1 | 3.0403 | 3210.5 | 3.2513 |
| 0.5072 | 3.1371 | 3098.0 | |
| 0.2191 | 3.1917 | 3027.7 | |
| 0.1048 | 3.2144 | 2996.7 | 3.3291 |
| 0.0526 | 3.2255 | 2981.1 | 3.3340 |
| 0.036 | 3.2292 | 2975.7 | 3.3356 |
| s_G -> 0 (quadratic fit, s_G <= 0.12, 12 points) | 3.2375 | 2963.6 | 3.3392 |
| patched conic (#943 enumerator) / published | 3.2383 / 3.24 | 2962.0 (alt 328) | 3.3395 / 3.34 |

Reading: the continuous orbit at full Ganymede mass is joined, by a smooth branch with no fold, to
an orbit whose s_G -> 0 limit reproduces the published GanCal#5 V_inf (to 0.001 km/s), periapsis
(r_p / s_G to 1.6 km of the patched-conic 2962 km) and Callisto speed (to 0.0003 km/s). So the
published orbit closes in the lane's continuous model in the patched-conic limit, and its
full-mass continuation sits 0.20 km/s lower in V_inf. The failed criterion 4 compared a
continuous-gravity orbit with a patched-conic display value; the tolerance (0.10 km/s, 150 km) was
set too tight for a 0.94 turn-ratio Ganymede flyby plus an unmodelled 95,000 km Ganymede pass.
This is my reading; the verdict as registered stays FAIL on criterion 4, and the lead rules on it.

## 5. Lead ruling 2026-10-07 and AMENDMENT 6 (written before the runs it covers)

Ruling: #968 is a PARTIAL. "The lane closes a GanCal#5-topology orbit at 3 S_GC in the continuous
R-S model; the published patched-conic values are not reproduced; criterion 4 FAILED as
pre-registered." Criterion 4 is not loosened. To make the shift explanation SHOWN, and to validate
the lane, three things, each pre-registered here.

### 5.1 GM continuation: tolerances (stated after the data in 4.7 exist; said plainly)

The 35 points of 4.7 were computed before this tolerance was written, so judging them now is
post hoc. Tolerance, as asked by the lead: the s_G -> 0 limit (quadratic fit, s_G <= 0.12) must give
V_inf within 0.01 km/s of 3.24 and altitude (r_p / s_G - 2634 km) within 10 km of 328 km; Callisto
speed within 0.01 km/s of 3.34. On the 4.7 fit these read 3.2375, 329.6 km, 3.3392 (post hoc).
A genuine prediction, written now from that fit, before the runs: continuing to s_G = 0.030 and
0.026, the pinned solution must give V_inf 3.2305 and 3.2314 km/s (within 0.002) and r_p / s_G
2973.7 and 2972.4 km (within 3 km).

### 5.2 Is the 0.2 km/s shift the unscheduled Ganymede pass? (independent patched-conic route)

- Done before this amendment, diagnostic only: the lead's suggested zero-SOI insertion of the pass as a
  scheduled Ganymede encounter (`two_working_body.correct_dates`, legs LG>G/1l | LG>G/0s | LG>C/0s |
  LC>G/0s, seeded at the GanCal#5 dates plus the pass at 25.15 d) converges, but to a DIFFERENT
  orbit: V_inf 2.143 km/s, turns 53.6 and 67.5 deg. The other three rev/branch splits of the B->C leg
  do not converge. A zero-SOI model forces the craft through Ganymede's centre, a 95,000 km move,
  so it cannot represent a distant pass. Not used as evidence either way.
- Pre-registered route: an impulsive-kick patched conic (`scripts/run_968_kick.py`). Kepler legs
  about Jupiter; the scheduled Ganymede flybys stay zero-SOI, as in R-S; the B->C leg is propagated
  by Kepler steps, and at its closest approach to Ganymede the Ganymede-relative velocity is rotated
  toward Ganymede by delta = 2 atan(mu_G / (d v_rel^2)), with d the Kepler closest-approach distance,
  in the plane of the relative position and velocity. Unknowns: the three dates and the B->C
  departure velocity (planar). Residuals: arrival at Callisto (position) at its date, V_inf magnitude
  match at B and A, vector match at Callisto. Solved by least squares from the GanCal#5 dates; the
  kick is recomputed inside every residual call.
- Reading, fixed now: Ganymede V_inf within 0.05 km/s of 3.040 and Ganymede altitude within 100 km
  of 576 km: the shift is the pass, SHOWN. V_inf within 0.05 km/s of 3.24: the pass is NOT the
  cause. Anything between: the pass explains a stated fraction (V_inf shift / 0.198).

### 5.3 Second published control: GanEur#43 (Russell & Strange 2009)

- Choice: GanEur#43 against GanEur#5. #43 is in the #943 ge gauntlet as a LITERAL recall (key
  `k2|LGanymede>Europa/1h|LEuropa>Ganymede/1l`, R-S Table 2 constants, Europa massless); #5 would
  need a new solve. Screen for unscheduled Ganymede passes (Kepler legs, 0.002-d grid, local minima of
  the Ganymede distance): #43 has none inside 51 Hill radii apart from its one scheduled flyby
  (GanCal#5, same screen: the 94,988 km = 3.00 R_H pass). Ganymede turn ratio 0.37 (22.7 of 61.0 deg).
- Published values: V_inf G 1.87, E 3.89 km/s; Ganymede altitude 8861 km; r 564,558 to 1,072,330 km;
  period 14.10 d (2 G-E synodic periods).
- Model and criteria: as 3.1 and 3.4, with Europa as the massless target (its centre hit pinned) and
  one Ganymede flyby node per period (one periapsis gauge). Criterion 3: Ganymede inside its Hill
  radius only at the scheduled flyby, altitude >= 100 km. Criterion 4: V_inf G within 0.10 of 1.87,
  Europa speed within 0.10 of 3.89, altitude within 150 km of 8861, r_min and r_max within 2 %.
  Criterion 5: the lane (fixed wrap) seeded from corrector A's solution (not from periapsis states)
  holds it and agrees on V_inf (1e-3 km/s) and altitude (1 km). Criterion 6: IAS15 half-arc re-fly.
- Method: corrector A, pinned layout, trf from the patched-conic seed. If it does not converge in
  200 evaluations, the fallback is the 4.3-4.4 route (T free, then continuation in T to 2 S_GE, then
  the Callisto-style pin by time shift), recorded as such.
- Expected (stated now): converges; V_inf G within 0.03 km/s and altitude within 100 km of the
  published values (moderate turn, no unscheduled pass); criteria 1-6 pass.

### 5.4 Context for #1023 (the VOID re-run): seed type matters

The lane's residual, with the fixed wrap, holds a correct orbit (4.6), but from patched-conic
PERIAPSIS seeds with node-to-node forward legs it stalls at km/s velocity defects (GanCal#5: 1.5-3.5
km/s after 160 evaluations; Member D and EGGIE used the same seed type). A re-run of the voided
EGGIE / #318 / #501 results should seed from a forward-backward (mid-leg match) solution, not from
periapsis states.

## 6. Results of amendment 6

### 6.1 GM continuation (5.1)

- Post-hoc tolerances on the 4.7 fit: V_inf 3.2375 (|diff| 0.0025 <= 0.01), altitude 329.6 km
  (|diff| 1.6 <= 10), Callisto speed 3.3392 (|diff| 0.0008 <= 0.01): within, but judged after the
  data.
- Fresh prediction (written in 5.1 before the runs): s_G = 0.030 gave V_inf 3.2306 (predicted 3.2305)
  and r_p / s_G 2973.7 km (2973.7); s_G = 0.026 gave 3.2315 (3.2314) and 2972.3 km (2972.4). Both
  within 0.002 km/s and 3 km: PASS. The branch from the full-mass orbit to the published patched-conic
  limit is shown.

### 6.2 The unscheduled pass, impulsive-kick patched conic (5.2), `data/968_control/kick.json`

- Tool check: with the kick switched off, the solver reproduces GanCal#5 exactly (V_inf 3.2383,
  altitude 328.0 km, Callisto 3.3395; residual 7e-12).
- With the kick: the pass is at 91,414 km, v_rel 2.670 km/s, turn 1.739 deg. Ganymede V_inf 3.1937
  km/s, altitudes 331.9 / 332.6 km, Callisto speed 3.3217 km/s (residual 6e-5; the closest-approach
  search tolerance; not pushed further).
- Reading fixed in 5.2: neither branch. The pass explains 0.045 of the 0.198 km/s V_inf shift (23 %)
  and 4 km of the 248 km altitude shift (about 2 %). My earlier inference (4.1, 4.5) that the pass
  is the cause is therefore WRONG for most of the shift. What 6.1 does show is that the whole shift
  is a smooth finite-Ganymede-mass effect on one branch. Which part of continuous gravity carries the
  rest (the near-limit flyby itself, at turn ratio 0.94, or Ganymede's distant pull along the legs)
  is not decomposed: an open item, not needed for the lane verdict.

### 6.3 Second control GanEur#43 (5.3), `data/968_control2/`

Corrector A converged straight from the patched-conic seed (61 evaluations; FD check of all 11
columns at the seed, worst 1.5e-7). `a_check.json`, `lane_state.json`:

| Criterion | Value | Verdict |
|---|---|---|
| 1. continuity | 5.4e-8 km, 1.8e-11 km/s; gauge 2e-15 | PASS |
| 2. Europa hit | 0 km (pinned) | PASS |
| 3. sequence; alt >= 100 km | no Ganymede approach inside its Hill radius other than the flyby; 10,036 km | PASS |
| 4. V_inf G within 0.10 of 1.87 | 1.9993 (diff 0.129) | FAIL |
| 4. Europa speed within 0.10 of 3.89 | 4.1755 (diff 0.286) | FAIL |
| 4. alt within 150 km of 8861 | 10,036 (diff 1175) | FAIL |
| 4. r_min within 2 % of 564,558 | 550,879 (-2.4 %) | FAIL |
| 4. r_max within 2 % of 1,072,330 | 1,085,176 (+1.2 %) | pass |
| 5. lane (fixed wrap) seeded from A | seed 7.6e-3, converged 4.8e-8; V_inf agrees to 6e-9 km/s, node distance to 1e-5 km | PASS |
| 6. IAS15 half-arcs | max 3.6e-6 km, 4.1e-11 km/s; full period 4.1e-6 km | PASS |

So the second control, chosen for having no unscheduled pass and a moderate turn ratio (0.37),
ALSO fails criterion 4, with a V_inf shift of +0.13 km/s at Ganymede and +0.29 km/s at Europa. My
expectation in 5.3 (within 0.03 km/s) was wrong. Two controls now show that the patched-conic
values of these low-V_inf (1.9-3.2 km/s) Ganymede cyclers move by 0.1-0.3 km/s when Ganymede's
gravity is continuous. That suggests criterion 4 at full mass tests the patched-conic
approximation, not the lane. This is my reading, for the lead.

### 6.4 AMENDMENT 7 (before the run): GanEur#43 GM continuation

Same method as 4.7/5.1 (pinned, fixed T, Ganymede GM scaled by s_G, softening radius by s_G,
secant predictor in log s). Pass rule, fixed now: the quadratic fit over the converged points with
s_G <= 0.12 must extrapolate to V_inf G within 0.01 km/s of 1.87, Europa speed within 0.01 km/s of
3.89, and altitude (r_p / s_G - 2634 km) within 20 km of 8861 km (the R-S table value 8861; our
patched-conic solve gives 8862.2).

### 6.5 GanEur#43 GM continuation result (amendment 7), `data/968_control2/gm.json`

58 converged points from s_G = 1 to 0.0055 (stopped at 0.0047, |r| 0.21 within 60 evaluations).
V_inf G falls monotonically 1.9993 -> 1.8722 km/s, Europa speed 4.1755 -> 3.8981 (at s_G = 0.02),
r_p / s_G 12,669.6 -> 11,541.2 km.

Judged as registered (quadratic fit over the 17 points with s_G <= 0.12):
- V_inf G 1.8716 km/s (|diff| from 1.87: 0.0016 <= 0.01): PASS.
- Europa speed 3.8930 km/s (|diff| from 3.89: 0.003 <= 0.01): PASS.
- altitude 8895.8 km (|diff| from 8861: 34.8 > 20): FAIL as registered.
Descriptive, not judged: r_p / s_G is still falling at the smallest s_G (8907 km altitude at s_G = 0.0055,
about 6 km per 15 % step), so the limit is not polynomial in s_G over this range; a quadratic fit in
sqrt(s_G) gives 8796 km. The two fits bracket 8861 km; the altitude limit is fit-dependent here, and
the pre-registered rule is the verdict.

## 7. Lead verdict on GanCal#5 (2026-10-08, ideal continuous model)

PASSED, recorded as: criterion 4 FAIL as registered (it compared a full-mass continuous orbit with
patched-conic numbers: a model mismatch in the criterion, now documented); identity shown by the
pre-registered GM continuation: the branch is GanCal#5's and reproduces the published values in the
patched-conic limit (V_inf 3.2375 vs 3.2383 patched conic vs 3.24 published; r_p / s 2963.6 vs
2962.0; Callisto 3.3392 vs 3.3395). Amendment 5 (my softening-scaling error, fixed before the
rerun) stays as written.

"Lane validated" needs the second control. GanEur#43 was run (6.3, 6.5) before this verdict
arrived, with criteria 1-6 pre-registered (5.3) and the continuation limit as amendment 7 (in
effect criterion 4b): criteria 1, 2, 3, 5, 6 PASS; criterion 4 FAIL as registered; 4b V_inf and
Europa speed PASS, altitude FAIL as registered (34.8 km against a 20 km rule; the curve is not
converged at s_G = 0.0055). The lead rules on whether that completes "lane validated".

## 8. VERDICT (lead ruling 2026-10-08): lane validated, ideal continuous R-S model

#968 = "lane validated (ideal continuous R-S model)". Two published controls (GanCal#5, GanEur#43):
- close in the lane with the fixed wrap;
- are held by the lane's own residual;
- agree with an independent integrator (REBOUND IAS15);
- continue on one smooth branch to the published V_inf and target speed in the patched-conic limit
  (both PASS as pre-registered).

Criterion 4 at full Ganymede mass FAILS for both. It tests the patched-conic approximation, not the
lane: 0.1-0.3 km/s shifts for 1.9-3.2 km/s Ganymede cyclers. The 95,000 km pass explains 23 % of
GanCal#5's shift; my earlier inference that it caused the shift is withdrawn (6.2). The GanEur#43
altitude limit FAILS as registered (quadratic fit 8895.8 km against 8861, tolerance 20 km); the sqrt
fit brackets it (8796 km) and the limit is not converged at s = 0.0055. This is recorded as a
fit-rule limitation and was not re-registered after the fact. No deeper altitude run and no GanEur#5
(option (i)). Real ephemeris (rung (b)) is not covered by this verdict.

### 8.1 Criterion 4c for any FUTURE control (pre-registered now)

The continuation in the working body's GM must reach s <= 0.01. The last computed point must lie
within a tolerance stated before the run of the published value (V_inf, target speed, and altitude
via r_p / s minus the radius). Any extrapolation rule (form of the fit, range of s) is fixed before
the run, and the verdict is the last computed point unless the registered fit rule says otherwise.

### 8.2 For the catalogue and the n-body rung of candidate rows

R-S patched-conic values shift by 0.1-0.3 km/s in continuous gravity at these V_inf (1.9-3.2 km/s at
Ganymede). A candidate's n-body standing must be judged by closure and by continuation to its own
patched-conic limit, not by agreement with its patched-conic numbers at full mass. (This is the
wording the n-body rung of #1025 and of the #942/#943 rows needs.) See also #1004 (gc-2's continuous
branch folds at 0.164 of the real moon masses): a patched-conic real-ephemeris PASS does not imply
that the orbit exists in continuous gravity.

## 9. Rung (b), GanEur#316 on jup365 at 2019: AMENDMENT 8 to sec. 3.6 (before any rung-(b) run)

What is published (AAS 07-118 p.15, Fig. 10(a) title; generator note 6.27): "10 cycles in ephemeris
model, 40 G. & 10 E. flybys, start=4-24-2019, TOF=493.5 days, Delta-v_TOTAL=0 m/s". Topology,
epoch and ballistic status only; no per-flyby numbers. The per-flyby reference values below are
therefore OUR real-ephemeris patched-conic reconstruction (#943 note 6.38,
`data/943_ganeur316_realeph/n10_rs2019_rel/`), not published numbers.

Sec. 3.6 is amended (it said Europa stays massless): the chain was solved in the both-massive cell
"ge" (the chain tool refuses massless targets in 3-D) and Europa bends in it (gate ratio 0.569), so a
massless Europa would contradict the seed.

1. Seed reconstruction, its own positive control. The chain tool is copied from commit 1e4b7aeb (the
   version that produced the 6.38 run) into my scratch directory, run once with `chain_eval` at
   lambda = 1 in the chain's own system, and the flyby list (body, epoch, V_inf in/out) plus the tool
   hash written to `data/968_rungb/seed_chain.json`; it is not imported again. PASS only if the
   reconstruction reproduces the stored result: max residual < 1e-6, gate worst 0.804 (to 1e-3),
   minimum required altitude 982 km (to 1 km).
2. Time and frame: the chain's times (seconds past JD 2440000, `core.Ephemeris` spice backend) are
   converted to the lane's TDB seconds past J2000; Ganymede's and Europa's positions from both must
   agree to < 1 km at one seed flyby epoch, else stop.
3. Model: jup365 rails (a spline wrapper with `.state()`, validated once against spkezr: position
   < 0.1 km and velocity < 1e-6 km/s at 20 epochs), Ganymede and Europa massive with GMs and radii
   scaled by sigma; Io and Callisto OFF during the continuation (the chain does not model them).
4. Open chain, N = 1 cycle (5 flybys) for the verdict, then N = 2 (descriptive if N = 1 passes);
   10 cycles only as a descriptive extension within 8-minute calls. Nodes at each flyby periapsis
   (state, epoch, periapsis gauge relative to its moon), forward-backward mid-leg matches, no wrap.
   The FIRST node is anchored (state and epoch fixed at the sigma-scaled patched-conic periapsis of
   the first flyby), so the system is square: 7N - 7 unknowns, 6(N-1) + (N-1) residuals.
5. Continuation in sigma: start sigma = 0.005 (seed offsets scaled), natural continuation upward to
   1 (factor 1.2, halving to 1e-5); a fold, if any, is checked as in #1004 amendment 1.
6. Acceptance at every point (binding constraints): each node is a periapsis of its moon (gauge),
   within 0.5 SOI of it, and above the sigma-scaled floor (Ganymede 100 km, Europa 100 km, times
   sigma, above sigma x radius); an unscheduled pass inside any of the four moons' Hill radii
   (sigma-scaled for G and E, unscaled for Io and Callisto, which are off) is recorded and stops the
   run.
7. Criteria (4c form): (i) closure at sigma = 1 at the lane floors (1e-3 km, 1e-6 km/s); (ii) the
   continuation reaches sigma <= 0.01 and the last computed point (sigma = 0.005) is within 0.01 km/s
   of the reconstructed chain's per-flyby V_inf; (iii) IAS15 re-fly of every half-arc at sigma = 1 on
   the lane's own `JovianRailsCache(JovianEphemeris)` and `JovianRestrictedNBody` with
   moons = (Ganymede, Europa) at registry GMs, agreement < 1e-2 km, 1e-7 km/s (a wall-clock stop
   is reported as "timeout"). Agreement with the chain's values at full mass is NOT a criterion.
8. Extra step at sigma = 1 (descriptive): switch Io and Callisto on at full mass and re-converge;
   record whether it converges and the V_inf shifts.
9. Expected outcome: closes at sigma = 1 for N = 1 and N = 2 (GanEur#316 has turn ratios up to 0.80
   at 3.2 km/s and the ideal-model controls closed); full-mass V_inf shifted by 0.1-0.3 km/s.
10. Meaning of a pass: the lane works on jup365 rails, in 3-D, with two massive moons, as an open
   chain. It does not validate periodicity (no wrap) or any candidate.

### 9.1 Rung (b) so far, and AMENDMENT 9 (before the runs it covers)

- Step 1 (seed reconstruction, `scripts/run_968_rungb_seed.py`, tool pinned at 1e4b7aeb): PASS.
  Max residual 1.4e-8, gate worst 0.80384, minimum required altitude 982.05 km. The flyby list (50
  encounters over 10 cycles; 5 per cycle: G, G, G, E, G) is in `data/968_rungb/seed_chain.json`.
- Step 2 (time and frame): chain time + (2440000 - 2451545) x 86400 s = lane TDB seconds past J2000;
  Ganymede and Europa positions from the chain's `Ephemeris("spice")` and the lane's
  `JovianEphemeris` agree to 0.0 km at three seed epochs (same kernel, same J2000 frame). PASS.
- Step 3 (spline ephemeris): at a 0.01-d grid the spline velocity error was 6.2e-6 km/s (fails the
  registered 1e-6); at a 0.004-d grid 4.6e-5 km and 4.1e-7 km/s: PASS (the grid is an implementation
  choice; the threshold is unchanged).
- Node coordinates changed to moon-relative (r_rel, v_rel, epoch): with absolute coordinates the
  sigma = 0.005 problem was so stiff that no Newton step descended (an epoch change of 1 s moved a
  20-km periapsis node by 10 km relative to its moon). A reparametrisation, not a criterion change.
- The anchored chain (node 0 a fixed state one day after the first departure) is WRONG as
  designed: a fixed full state at a fixed time fixes the whole downstream trajectory, so the
  encounters are no longer constraints, and the Jacobian has two exact null directions (singular
  values 3e-16; right null vectors on the last node's state, left null vectors on the anchor leg's
  rows). Newton stalled at |r| 4.18 (anchor-leg dv 3.9e-3 km/s). My design error, after the advisor
  suggestion "anchor the first node"; no result is drawn from it.
- AMENDMENT 9: no anchor. Nodes 1..M (the chain's flybys) are all free, with periapsis gauges and
  M - 1 forward-backward legs: 7M unknowns, 7M - 6 residuals, so six free directions (the open
  chain's end freedom). Steps are MINIMUM-NORM Gauss-Newton (column-scaled lstsq). Fold detection is
  NOT claimed under this choice; if the continuation stops, the result is reported as a numerical
  stop. Criteria 7 (i)-(iii), the acceptance checks (item 6) and the meaning (item 10) are unchanged.
  With N = 1 the chain is the five flybys of cycle 1 (G 9.39 d, G 20.12, G 29.51, E 37.10, G 49.34).
