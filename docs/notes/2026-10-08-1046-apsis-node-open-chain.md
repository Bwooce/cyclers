# #1046: a full-revolution-aware open chain (apojove nodes)

Status: PRE-REGISTRATION (written and committed before any #1046 corrector run). Lead ruling
2026-10-08. Time box: 4 hours of calls from 16:05 AEDT (ends 20:05 AEDT); I stop with a report
either way. No catalogue writes; nothing is called novel.

## 1. Problem (#1044 sec. 5)

The two real-ephemeris objects with a FULL-revolution resonant leg do not start at sigma = 0.02 under
the periapsis-node, mid-leg-match open chain:
- gc-1 e2 stalls at a stationary point that is not a root. The residual sits on the velocity rows of
  the 16.7-day Callisto-Callisto 1:1 return, along the weakest left singular vector (1.5e-6).
- GanCal#1 at 2013: the Newton steps ask for node moves larger than the node-moon distances, on either
  side of its 2:1 return. The line search crawls.
The EGGIE (multi-rev Lambert legs, no resonant return) and GanEur#316 (half-revolution leg) start.

## 2. Formulation: apojove nodes

The lead offered two forms: a node at the resonant leg's apojove, or the return parametrised by its
V_inf direction on the resonant circle. I take the first, in a general form:
- One extra multiple-shooting node at EVERY apojove passage of the seed trajectory between two flyby
  nodes, kept only if it is at least 1 day from both flyby nodes. The seed apojove is the
  Jupiter-centred two-body conic leaving flyby node k with the seed's outbound asymptote (moon
  velocity + V_inf out from the seed hyperbola at the full moon GM), propagated by `kepler_step`.
- Unknowns per apsis node: its absolute Jupiter-centred state and epoch (7). Row: the
  Jupiter-centred apsis gauge r.v / (|r||v|) = 0 (weight 1e4, as the periapsis gauges). Its leg on
  each side has the usual forward-backward mid-time match (6 rows). The square structure is kept: 7
  new unknowns and 7 new rows per apsis node.
- Implementation: the apsis node is a pseudo-moon "APO" whose ephemeris state is zero, so the rungb
  node, gauge, leg, Jacobian and Newton code is reused unchanged. Apsis offsets are not scaled with
  sigma. `describe` and the acceptance skip them, so every criterion acts on the flyby nodes as before.
- Everything else as the object's previous run: end rows, solver, floors, noise-floor rule, sigma
  schedule (0.02 up by 1.2, halving to 1e-5; down to 0.01 and 0.005), acceptance, IAS15.

Inserted nodes (build stage, which only constructs the chain and evaluates the seed residual):
- EGGIE one cycle: E, G, APO, APO, G, APO, I, APO, APO, E. Apojoves in the G-G leg (1.57 and
  6.85 d after G), the G-I leg (4.08 d) and the I-E leg (2.68 and 8.40 d). The control exercises
  the new code five times.
- GanCal#1 at 2013: G, APO, G, APO, C, G. One apojove in the G-G leg (6.49 d) and one in the 2:1
  return (6.29 d after G).
- gc-1 e2: G, APO, C, APO, C, G. One apojove in the G-C leg (5.17 d) and one in the 1:1 Callisto
  return (11.91 d after C).
- Seed residual at sigma 0.02 (scaled offsets): about 2.1e3 km on every object, the same scale as
  #1044. The apsis nodes do not add seed error.

## 3. Order and criteria

1. Control, one-cycle EGGIE, paper's ideal model, #1039 formulation (b), damped Newton (the #1039
   PASS configuration). PASS = the #1039 sec. 3 criteria: sigma = 1 at the lane floors; the #1043
   acceptance (25 km floor, no unscheduled Hill-radius pass); identity going down (sigma 0.01 and
   0.005 converged at the floors, max |dV_inf| < 0.01 km/s against the chain); IAS15 < 1e-2 km,
   < 1e-7 km/s. Also it must reproduce the #1039 (b) full-mass flyby nodes: same orbit, so V_inf per
   node within 1e-3 km/s and r_p within 1 km of `data/1039_b_eggie/sigma_n1.json`. The apsis nodes
   change the parametrisation, not the problem.
   If the control fails, I stop, and I do not run the real-ephemeris objects.
2. GanCal#1 at 2013 on jup365, #1044 formulation (b') with LM then Newton (#1044 amendments 1-2),
   the #1044 sec. 1 PASS criteria (i)-(iv).
3. gc-1 e2 on jup365, the same.

Each call runs under `timeout 470` (8-minute rule). Every stage resumes from its checkpoint (the
sigma JSON and the LM checkpoint).

## 4. What each outcome means (fixed now)

- Control PASS, and a real-ephemeris object STARTS at sigma = 0.02: the #1044 non-start was the
  formulation, as #1044 inferred. That object's run then goes up as far as it goes, and its result
  is a real verdict under criteria (i)-(iv).
- Control PASS, and an object still does not start: the apojove split is not enough. Report the
  residual structure (which rows, smallest singular values, compared with #1044). The second form
  (the return parametrised by its V_inf direction on the resonant circle) would be the next step,
  and it is not run in this time box unless time remains and the lead has ruled.
- GanCal#1 PASS: the jup365 lane is validated by a published control. FAIL at an interior node:
  recorded as in #1044 sec. 2.

## 5. Expectations (fixed now)

- EGGIE control: PASS with probability about 0.8. Extra shooting nodes on a problem that already
  converged should not move the root.
- gc-1 e2: starts with probability about 0.5. The apsis node halves the propagation span either side
  of the 1:1 return's match. I am not sure that this removes the near-degenerate direction, which is
  physical in the patched-conic limit (the crank on the resonant circle) and only weakly fixed by the
  flybys at small sigma. Reaches sigma = 1, given a start: about 0.6.
- GanCal#1: starts with probability about 0.35. Its weak flybys (4.0 and 1.0 deg turns) were part of
  the #1044 problem, and apojove nodes do not touch them.

## 6. Results

### 6.1 Control: one-cycle EGGIE with five apsis nodes (`data/1046_eggie/`): PASS

- Sigma 0.02 to 1 in 23 points, every point at the lane floors (no noise-floor point). No
  unscheduled Hill-radius pass at sigma = 1.
- Full mass: V_inf 9.0323 / 7.0457 / 7.0411 / 7.8799 / 9.0323 km/s. Every flyby node equals the #1039
  (b) result (`data/1039_b_eggie/sigma_n1.json`) to < 1e-6 km/s in V_inf and < 1e-3 km in r_p. Same
  orbit, as expected. Altitudes 1,082 / 1,028 / 1,764 / 7,726 / 1,664 km (floor 25 km).
- Identity going down: sigma 0.01 max |dV_inf| 0.0044, sigma 0.005 0.0023 km/s, both at the floors.
  PASS.
- IAS15 re-fly of every half-arc at sigma = 1: max 1.3e-7 km. PASS.
- Wall time 6.4 min for the whole control (machine load 8-70).

### 6.2 AMENDMENT 1 (solver bookkeeping only; before the real-ephemeris calls that follow it)

- The #1044 LM wrapper checkpoints the LAST evaluated point. That can be a rejected LM trial, so a
  call stopped by `timeout 470` resumes from a worse point. #1046 uses its own copy of the wrapper,
  which keeps the BEST evaluated point (`_install_lm_best`). Criteria and everything else are
  unchanged. The first GanCal#1 call (16:10-16:18) and the first gc-1 call (16:19-16:27) ran before
  this change. Their final checkpoints were each that call's best point, so nothing is lost.
- Cost under this load (40-70): a residual-plus-Jacobian evaluation takes about 17 s for GanCal#1
  and about 110 s for gc-1. Every arc that starts at a flyby node takes 15-35 s at sigma 0.02, while
  arcs from an apsis node take 0.1 s. So gc-1 gets about 4 LM evaluations per call.

### 6.3 Real-ephemeris objects at sigma = 0.02 (`data/1046_gc1/`, `data/1046_gancal1/`)

- gc-1 e2 (21 LM evaluations, 4 calls): |r| 2.8e3 -> 0.107 -> 0.0275, then flat (0.0278-0.0275 over
  the last 6 evaluations). DOES NOT START. Structure (`diag.json`): 99.9 % of the residual lies on the
  weakest left singular vector (singular value 2.3e-6; #1044 without apsis nodes: 1.5e-6). It is the
  velocity rows of leg 3, the APO-to-Callisto second half of the 1:1 return (dv 2.7e-5 km/s). Every
  other leg is at or near its floor. Same structure as #1044. The apsis node did not change the
  conditioning.
- GanCal#1 at 2013 (about 200 LM evaluations, 8 calls): |r| 3.3e5 -> 4.24 in 9 evaluations, then a
  slow crawl: 4.24 -> 2.53, slowing to about 5 % per 6 minutes (the last 75 evaluations). DID NOT
  START within about 200 evaluations; still descending, so no stationary point is shown. 99.7 % of
  the residual lies on the weakest left vector (singular value 1.1e-5), on the velocity rows of leg
  0, the G-to-APO first half of the 1-rev G-G leg (dv 2.5e-3 km/s).
- Frame of the weak direction (`frame_shares` in `diag.json`; each block projected on the orbit
  frame (r_hat, h_hat x r_hat, h_hat) of its node's moon; the lane frame is the ephemeris frame, not
  Jupiter-equatorial, so raw z components were not used):
  - gc-1: the stalled leg-3 velocity residual is 100.0 % orbit-normal. The weakest right vector is
    dominated by the APO node in the 1:1 return (block norm 0.44; position and velocity 100 %
    orbit-normal). The Callisto nodes on either side are 96-99 % orbit-normal.
  - GanCal#1: the leg-0 velocity residual is 99.3 % orbit-normal. The weakest right vector is
    dominated by the APO node in the G-G leg (block 0.20; 94 % / 96 % orbit-normal) and the APO node
    in the 2:1 return (block 0.02; 98-99 %).
  So in both objects the near-degenerate direction is the OUT-OF-PLANE motion across a
  full-revolution leg, and the residual that will not go is an out-of-plane velocity mismatch there.
  This fits the Kepler limit, where the out-of-plane variational map over a full period is close to
  the identity: an out-of-plane change at one end of a full-revolution leg comes back almost
  unchanged at the other end, so the matching rows hardly see it (INFERRED, not tested).

### 6.4 AMENDMENT 2: free-end DIAGNOSTIC (registered before it runs; not a verdict)

- Question: is the obstruction the (b') end pins, or the interior full-revolution leg? The answer
  decides which form comes next: the resonant-circle parametrisation (interior) or a different end
  form (pins).
- Method (`--stage freeend`): from each object's sigma-0.02 LM checkpoint, the six end rows are set
  to zero rows (the system is underdetermined) and the plain damped Gauss-Newton of rungb
  (min-norm lstsq steps, step cap, backtracking) runs for 20 iterations. LM is not used: it needs at
  least as many rows as unknowns. Recorded: convergence of matches and gauges; the drift of the end
  V_inf directions and of the magnitude gap from the pins; node and epoch moves.
- Readings, fixed now:
  - matches and gauges reach the lane floors with free ends -> the pins are the obstruction (the
    pinned problem has no root near the chain), and the drift is the size of the inconsistency;
  - they stall at the same out-of-plane velocity residual -> the obstruction is interior (the
    full-revolution leg), and the end form is not the issue.
- Order: GanCal#1 first (about 5 s per evaluation now), gc-1 if the time box allows (35-110 s per
  evaluation).

### 6.5 Amendment-2 diagnostic results: the obstruction is INTERIOR (the full-revolution leg)

- GanCal#1 (`data/1046_gancal1/freeend.json`): with the end rows removed, the first min-norm
  Gauss-Newton step does not lower the residual (2.532 -> 2.532), and the line search fails, so the
  run stops after one step. The leg-0 velocity mismatch stays at 2.47e-3 km/s. The end directions
  drift only 3.8e-5 and 2.2e-6 rad, the magnitude gap 2.1e-5 km/s and the span 0.11 s. Freeing the
  ends gives the solver nothing to use.
- gc-1 e2: the call reached its 470 s limit after two Gauss-Newton evaluations (about 3.5 min each
  under load 19), so there is no JSON. From the runlog: |r| 0.02750 -> 0.02746, with the leg-3
  velocity mismatch unchanged at 2.74e-5 km/s. The same stall with free ends.
- Reading (registered in 6.4): the (b') end pins are NOT the obstruction. The residual that will not
  go is interior: an out-of-plane velocity mismatch across the full-revolution leg (gc-1: the 1:1
  Callisto return; GanCal#1: the 1-rev G-G leg), along a near-null out-of-plane direction. The
  apojove nodes do not remove it. The diagnostic is short (one and two steps), so this is a
  diagnosis, not a proof.

## 7. Verdict of #1046 (time used: 16:04-17:21 AEDT, about 1 h 20 min of the 4 h box)

- Formulation: apojove nodes are sound. The ideal EGGIE control reproduces #1039 (b) exactly, with
  five apsis nodes.
- Real ephemeris: NEITHER object starts at sigma = 0.02.
  - gc-1 e2: stationary at |r| 0.0275, the same structure as #1044.
  - GanCal#1: still crawling after about 200 evaluations; no stationary point shown.
  NO VERDICT on either object. The jup365 lane is still not validated by a published control. The gc-1
  row's real-ephemeris continuous-gravity standing remains UNDECIDED (#1044 wording stands).
- New: the obstruction is located. It is an out-of-plane velocity mismatch across the
  full-revolution leg, along a near-null out-of-plane direction (moon orbit frame, 94-100 %), and it
  is not caused by the end pins. A next formulation must give that leg an out-of-plane control that
  the matching rows can see. Candidates, not run (the lead rules):
  - the resonant-circle parametrisation with the crank (out-of-plane rotation of V_inf about the moon
    velocity) as an explicit unknown;
  - or a node at the full-revolution leg's line of nodes (where an out-of-plane velocity change is
    most effective) instead of its apojove.
