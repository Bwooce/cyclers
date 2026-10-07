# #1041: the continuous-gravity EGGIE in the paper's model, quasi-periodic criterion

Status: PRE-REGISTRATION (written and committed before any #1041 run). No catalogue writes; nothing
is called novel. Lead ruling 2026-10-08 (after #1023 sec. 7).

## 1. Starting point

In the paper's ideal model (Hernandez et al. 2017; `resonant_conic.ideal_moon_smas` with
T_syn = `ideal_t_syn()`, the ideal Ganymede period, 7.0042 d; Laplace angle 180 deg) the #1023
date corrector reproduces Table 4: root (iii) (`data/1023_eggie/pc_roots_coded.json`, legs E>G 0s,
G>G 1h, G>I 1l, I>E 1h; V_inf E 9.068, G 7.082, I 8.207 km/s; leg times 1.57 / 8.42 / 7.27 / 10.76 d;
gate pass). Over one cycle T = 4 T_syn = 28.017 d the moons advance by different angles (Ganymede 0,
Europa -20.50, Io -61.51 deg), so no orbit repeats exactly: the published EGGIE is quasi-periodic.
The date corrector accepts it because it matches only V_inf MAGNITUDES at the junctions.

## 2. Method (`scripts/run_1041_eggie.py`, data `data/1041_eggie/`, live log `data/1041_eggie/live/`)

1. Patched-conic n-cycle chain, n = 3. The date corrector on the 3-cycle cycle (the four legs x 3,
   period 3T), seeded from root (iii) tiled over three cycles. Accepted: max residual < 1e-9 km/s,
   Kepler re-fly miss < 1 km, gate pass at every flyby (25 km floor). Its 12 flybys are the seed. The
   wrap junction is a magnitude match, as for the paper.
2. Continuous open chain. Nodes at the 13 encounters (E G G I x 3, then E), in moon-relative
   coordinates with each node's epoch, a periapsis gauge relative to its own moon, and forward-backward
   mid-time matches on the 12 legs. Like #968 amendment 12, the first node's inbound and the last
   node's outbound V_inf vectors are pinned to the patched-conic chain's (6 rows); the system is
   square (91 x 91). Io, Europa and Ganymede are point masses, with GM and softening radius scaled by
   sigma. DOP853 + analytic STM, rtol 1e-13.
3. Continuation in sigma from 0.02 up to 1 (factor 1.2, halving to 1e-5), damped Newton with the
   half-distance step cap. Intermediate points may be noise-floor-limited (dr < 1e-2 km,
   dv < 1e-6 km/s; #968 amendment 11); verdict points may not. Then down to 0.01 and 0.005 for the
   identity.
4. Fallback, if the n = 3 chain does not converge at sigma = 0.02 within 25 Newton iterations: the
   same with n = 1 and n = 2 (diagnostic; the verdict is then stated for the largest n that works).

## 3. Criteria

EGGIE CLOSES as a quasi-periodic continuous-gravity object (n cycles) if:
1. at sigma = 1 every leg match meets the lane floors (1e-3 km, 1e-6 km/s), the gauges are < 1e-9
   and the end pins are < 1e-6 km/s;
2. every node altitude is >= 25 km (the paper's floor), and there is no unscheduled pass inside any
   moon's Hill radius;
3. identity: at the smallest sigma <= 0.01 converged at the floors, every flyby V_inf is within
   0.02 km/s of the patched-conic chain (criterion 4c form); full-mass agreement with Table 4 is not a
   criterion;
4. IAS15 re-fly of every half-arc at sigma = 1 (`JovianRestrictedNBody`, exact circular rails,
   registry GMs): < 1e-2 km, 1e-7 km/s.

Quasi-periodicity, MEASURED, not imposed, and reported per moon at sigma = 1: for each encounter of
cycle k+1 against cycle k, the change in V_inf and in periapsis altitude, and the change in the
moon-relative periapsis state after rotating cycle k's by THAT moon's own advance over T (Ganymede 0,
Europa -20.50, Io -61.51 deg).

## 4. Outcomes and meaning (fixed now)

- CLOSES: EGGIE exists as a quasi-periodic ballistic object in continuous gravity in the paper's
  model, over n = 3 cycles, with the measured cycle-to-cycle drift.
- DOES NOT CLOSE: it stops before sigma = 1 on a fold (#1004 checks), on an impact or acceptance
  failure (the node and the growing defect named), or because the n = 3 chain does not start (the
  fallback is then reported). The defect that grows is identified: which leg's match, which node's
  altitude, or which end pin.
- NUMERICAL STOP: undecided, with the last point and the failure mode.

Expected outcome: CLOSES with probability about 0.5. EGGIE's V_inf of 7-9 km/s favour small
patched-conic to continuous shifts and the turns are small (1-6 deg), but #968 rung (b) showed that a
pinned 3-cycle chain can fail to start.
