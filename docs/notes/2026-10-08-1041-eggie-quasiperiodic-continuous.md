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

## 5. Result: the patched-conic n-cycle seed fails the gate (n = 3 and n = 2); #1041 stops at step 1

`data/1041_eggie/pc_chain_n3.json`, `pc_chain_n2.json`. The date corrector, seeded from root (iii)
tiled over n cycles, converges exactly in both cases (residual 2e-13 and 1.6e-13 km/s; re-fly miss
5e-6 and 1e-6 km). Neither root passes the gate, so neither passes the pre-registered acceptance of
step 1:

| n | flybys in time order: body, V_inf (km/s), demanded turn (deg) |
|---|---|
| 3 | G 6.818 5.9, G 6.818 3.5, I 7.827 0.2, E 8.917 3.4, G 6.224 4.4, G 6.224 3.0, **I 5.418 86.7**, E 8.065 0.1, G 6.320 3.8, G 6.320 0.04, **I 7.354 131.2**, E 8.622 2.4 |
| 2 | G 5.938 0.2, G 5.938 0.4, I 5.840 2.4, E 8.166 3.2, G 5.836 3.5, G 5.836 0.6, **I 5.410 90.9**, E 7.940 1.0 |

- The defect that grows is the Io turn in the second cycle (87-91 deg) and in the third (131 deg).
  The first cycle stays small-turn. This is consistent with the model: Io's position shifts by
  -61.5 deg per cycle relative to Ganymede (and Europa by -20.5), so the cycle-1 geometry does not
  carry over at Io.
- The V_inf also drift away from Table 4 when n cycles are solved together (G 5.9-6.8 km/s against
  7.07): the closed n-cycle chain distributes the mismatch over all cycles.
- By sec. 4 this is DOES NOT CLOSE, with the growing defect identified (the cycle-2 Io turn), at the
  patched-conic seed stage; the continuous steps 2-4 were not run. The n = 1 object is root (iii)
  itself (one cycle, a gate pass).
- This matches the paper (p.10): ballistic repeatability "in general will only last for a few
  cycles", then maintenance Delta-V.
- Not tested (a more faithful quasi-periodic test, not pre-registered here): MARCHING the open chain
  cycle by cycle from root (iii), with each new cycle's four dates solved from the junction magnitude
  matches given the previous cycle's end. That would show after how many cycles the gate first fails
  without spreading the mismatch over all cycles. Proposed as a follow-up.

## 6. Cycle-by-cycle MARCH: PRE-REGISTRATION (lead go 2026-10-08; before running)

Script `scripts/run_1041_eggie.py --stage march`; data `data/1041_eggie/march_<start>.json`.

- Model: the paper's (sec. 1). Legs per cycle are E>G 0s, G>G 1h, G>I 1l, I>E 1h (root (iii)'s
  topology).
- March (open chain, no spreading): cycle k starts at the Europa encounter E_k (its date fixed by the
  previous cycle, its inbound V_inf the previous I>E arrival). The four unknown dates G1, G2, I and
  E_{k+1} are solved from the four junction magnitude matches at E_k, G1, G2 and I (least squares
  from the previous cycle's dates + T). Exact: max residual < 1e-9 km/s.
- Start A, root (iii): cycle 1 is root (iii) itself (its closing magnitude match gives E_1's inbound).
- Start B, the Table 4 values: E_0 at root (iii)'s E date. Cycle 1 has its G1, G2 and I dates solved
  from the three interior junction matches, seeded at the Table 4 leg times (1.59, 8.60, 7.34 d), with
  E_1 fixed at E_0 + 28.22 d (the Table 4 total). There is no inbound at E_0, so that junction is
  not imposed. It is then marched as in A.
- Recorded per cycle: V_inf at every node, demanded turns, gate ratios and status at the paper's
  25 km floor and at the project floors (`gate_cycle`). The march runs for up to 10 cycles, or until
  the corrector fails (recorded as such).
- Reported: the first cycle at which the gate fails, at each floor, for each start.
- Expected outcome: the gate fails by cycle 2 (probability 0.7) or cycle 3 (0.2), with Io as the
  failing flyby, from Io's -61.5 deg per-cycle shift.
- Conclusion rule (lead): if cycle 1 stands (a gate pass) and the march fails by cycle 2-3, record
  "the published EGGIE is a one-cycle ballistic object in its own model; ballistic repeatability is
  limited to that, consistent with the paper's own statement", then stop.
