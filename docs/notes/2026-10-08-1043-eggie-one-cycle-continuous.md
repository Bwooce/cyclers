# #1043: does the published one-cycle EGGIE exist in continuous gravity?

Status: PRE-REGISTRATION (written and committed before any #1043 run). No catalogue writes; nothing
is called novel. Lead ruling 2026-10-08, after #1041.

## 1. Object and model

Root (iii) of #1023 (`data/1023_eggie/pc_roots_coded.json`[0]) in the paper's ideal model
(`resonant_conic.ideal_moon_smas`, T_syn = `ideal_t_syn()`, Laplace angle 180 deg). Legs E>G 0s,
G>G 1h, G>I 1l, I>E 1h; V_inf E 9.068, G 7.082, I 8.207 km/s; leg times 1.57 / 8.42 / 7.27 /
10.76 d; it is the Table 4 EGGIE reproduced at 0.173 km/s, and it is one-cycle ballistic in its own
model (#1041).

## 2. Method (`scripts/run_1043_eggie.py`, a thin wrapper of `scripts/run_1041_eggie.py` with n = 1;
data `data/1043_eggie/`, live log `data/1043_eggie/live/`)

- The patched-conic one-cycle chain is root (iii) re-solved by the date corrector with n = 1. It must
  reproduce root (iii): V_inf within 1e-6 km/s, gate pass.
- Open chain of 5 nodes: E0, G1, G2, I, E1, in moon-relative coordinates with epochs, periapsis gauges
  and forward-backward legs. The ends are pinned as in #968 amendment 12 and #1041 sec. 2: E0's inbound
  is the closing I>E arrival rotated back by Europa's own advance; E1's outbound is the first departure
  rotated forward by it. Square 35 x 35.
- Io, Europa and Ganymede are massive, with GM and softening radius scaled by sigma. DOP853 + STM,
  rtol 1e-13.
- Continuation from sigma = 0.02 up to 1 (factor 1.2, halving to 1e-5, half-distance step cap;
  intermediate points may be noise-floor-limited), then down to 0.01 and 0.005.

## 3. Criteria (CLOSES requires all)

1. At sigma = 1, every leg match meets the lane floors (1e-3 km, 1e-6 km/s), gauges < 1e-9, end pins
   < 1e-6 km/s.
2. Every node altitude >= 25 km (the paper's floor), no unscheduled pass inside any moon's Hill
   radius.
3. Identity: at the smallest sigma <= 0.01 converged at the floors, every node V_inf is within
   0.02 km/s of the patched-conic chain.
4. IAS15 re-fly of every half-arc at sigma = 1 (exact circular rails, registry GMs): < 1e-2 km,
   1e-7 km/s.

Reported at sigma = 1: per-node V_inf, turns (from the moon-relative hyperbola: 2 asin(1/e)) and
altitudes. Full-mass agreement with Table 4 is NOT a criterion (#968 sec. 8.2).

## 4. Outcomes and meaning (fixed now)

- CLOSES: the published one-cycle EGGIE exists in continuous gravity (ideal model). This would be the
  first continuous-gravity confirmation of any published Jovian triple cycler in this project.
- DOES NOT CLOSE: the node whose defect grows is named (an acceptance failure such as an impact, a
  fold with the #1004 two-solution check, or a stalled leg).
- NUMERICAL STOP: undecided, with the last point and the failure mode.

Expected outcome: CLOSES with probability about 0.7. The V_inf are high (7-9 km/s), the turns small
(1-6 deg), and the analogous one-cycle pinned control (#968 rung (b), N = 1) converged straight from
its seed; that one failed only at sigma 0.62, through a Europa flyby with a 7.2-deg turn.

### 2.1 AMENDMENT 1 (seed construction, before the continuation that counts)

- The first sigma = 0.02 attempt stalled at 4.1e4 km and 1.42 km/s. Diagnosis: the I>E leg's seed
  defect did not shrink with sigma (2.7e5 km at both sigma = 0.001 and 0.02). `periapsis_node` clamps
  the periapsis distance to 0.6 x the moon's SOI by default. Io's 0.9-deg turn at 8.2 km/s needs
  r_p = 11,170 km, beyond the 6,330 km cap, so the clamped seed node turned by about 1.6 deg instead
  of 0.9: about 0.1 km/s of error carried over a 10.8-day leg.
- Fix: the seed nodes use `max_offset_km = 1e12` (no clamp). Seed defects now scale with sigma
  (sigma 0.001: <= 419 km, 0.014 km/s; sigma 0.02: <= 9,383 km, 0.32 km/s). This is a seed
  construction fix only; criteria unchanged.
- The same clamp did not affect #968 rung (b): its turns of 10-24 deg at Ganymede and 7.2 deg at
  Europa give periapses inside the cap.

### 2.2 AMENDMENT 2 (code aligned with the registered criteria, before continuing)

- The continuation converged at the floors at every step from sigma = 0.02 up to 0.178. It then
  stopped at sigma = 0.214 on an acceptance check that #1043 never registered: "periapsis within
  0.5 SOI". That check was inherited from the #968 rung-(b) code (amendment 8 item 6).
- The Io node tripped it (r_p / SOI = 0.51). Its 0.9-deg turn at 8.2 km/s needs r_p of about
  11,000 km at full mass, about 1.06 Io SOI even in the patched conic. So for this object the bound
  is wrong, not the orbit.
- The #1043 criteria (sec. 3) are the 25 km floor and no unscheduled pass inside a Hill radius, and
  they are unchanged. The inherited SOI bound is switched off for #1043 (`SOI_LIMIT = inf` in the
  wrapper), and the continuation resumes from its checkpoint at sigma = 0.178.

## 5. Result: DOES NOT CLOSE, a fold at sigma = 0.5050 driven by the pinned end Europa nodes

Data: `data/1043_eggie/sigma_n1.json`, `monitor_n1.json`, `two_n1.json` and `direct_sigma1.json`.

- Identity: at sigma = 0.02 every node V_inf matches the patched-conic chain (E 9.0679, G 7.082,
  I 8.198, against 9.068 / 7.082 / 8.207). The continuation is on root (iii)'s branch.
- Natural continuation met the lane floors at every accepted point from sigma = 0.02 to 0.505009 (no
  noise-floor-limited points). Steps beyond 0.50502 fail down to a 1e-5 step.
- Fold check (sec. 4, #1004 checks):
  - The monitor-variable continuation passed the turning point. Sigma peaked at 0.505011 and the
    second branch runs back down (0.4552 after 120 points).
  - Two solutions at fixed sigma: at 0.48 they differ by 0.082 km/s in V_inf and 1,634 km in r_p; at
    0.46 by 0.109 km/s and 2,192 km. Both above the 0.02 km/s and 100 sigma km thresholds: FOLD
    SUPPORTED.
- Direct Newton at sigma = 1, from the seed and from the lower branch's end: no convergence (dr 1.2e4
  and 7.2e4 km).
- The growing defect is at the END Europa nodes, which carry the pinned asymptotes:

  | sigma | r_p / sigma (km): E0 | G1 | G2 | I | E1 |
  |---|---|---|---|---|---|
  | 0.02 | 2,850 | 3,648 | 4,340 | 10,549 | 2,946 |
  | 0.37 | 2,871 | 3,532 | 4,267 | 9,689 | 4,368 |
  | 0.505 (fold) | 3,792 | 3,357 | 4,172 | 9,198 | 6,904 |
  | 0.455, upper branch | 5,259 | 3,238 | 4,101 | 9,002 | 10,220 |

  The interior nodes (G1, G2, I) drift by 8-13 %, while E1's scaled periapsis grows 2.3x to the fold.
  All altitudes stay above the floors; there is no impact.
- Reading (INFERRED): as in #968 rung (b), holding the end asymptotes at their patched-conic values
  forces the end flybys to absorb the full-mass shift, and the branch folds. The interior of the
  one-cycle EGGIE stays close to its patched-conic geometry up to half mass. So the fold says the
  pinned-end formulation has no solution beyond sigma = 0.505. It does not say the one-cycle EGGIE
  fails to exist; that question needs an end condition that lets the ends move (the #1039
  formulation question).
- Outcome by sec. 4: DOES NOT CLOSE (fold at sigma = 0.505, with the end Europa nodes as the growing
  defect). The continuous existence of the published one-cycle EGGIE is undecided under this
  formulation.

Amendments made before the runs they cover: 1 (unclamped seed periapses), 2 (the inherited SOI bound
removed). Expected outcome (CLOSES, probability 0.7) was not met.
