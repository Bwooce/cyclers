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

(pending)
