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
