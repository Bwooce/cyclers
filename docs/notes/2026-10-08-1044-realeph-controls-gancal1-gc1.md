# #1044: second real-ephemeris control (GanCal#1 at 2013) and gc-1 at one epoch, formulation (b')

Status: PRE-REGISTRATION (written and committed before any #1044 run). Lead ruling 2026-10-08. No
catalogue writes; nothing is called novel.

## 1. Method (common)

The #968 rung-(b) / #1039 machinery (`scripts/run_968_rungb.py`, END_MODE = "direction_seedmag"
(b')), on jup365 rails (spline wrapper, rtol 1e-13), one cycle as an open chain:
- moon-relative periapsis nodes at each encounter;
- end V_inf directions pinned, the end magnitude gap held at the seed's, and the span held;
- the encounter moons (Ganymede and Callisto) massive and scaled by sigma, Io and Europa off;
- continuation from sigma = 0.02 to 1 (intermediate points may be noise-floor-limited; verdict
  points may not), then down to 0.01 and 0.005.

Seeds are reconstructed from the stored #943 chain results by the chain tool, PINNED by `git show` at
the commit that produced each run. The reconstruction is its own positive control: the stored max
residual (< 1e-6), gate worst ratio (to 1e-3) and minimum required altitude (to 1 km) must be
reproduced.

PASS (each object): (i) sigma = 1 at the lane floors (1e-3 km, 1e-6 km/s), gauges < 1e-9, end rows
< 1e-6; (ii) every node altitude >= the sigma-scaled project floor (Ganymede 100 km, Callisto
200 km) and within 0.5 SOI (the rung-(b) acceptance), no unscheduled pass inside a Hill radius;
(iii) identity: at the smallest sigma <= 0.01 converged at the floors, every node V_inf within
0.01 km/s of the reconstructed chain; (iv) IAS15 re-fly of every half-arc at sigma = 1 on the
lane's `JovianRailsCache(JovianEphemeris)` (registry GMs): < 1e-2 km, 1e-7 km/s.

## 2. Control: GanCal#1 (R-S 2009, C4 path) at R-S's 2013 epoch

- Source: `data/943_c4_rs2013/n10_grow/` (note 6.44): cell gc, key
  `k3|LGanymede>Ganymede/1l|RGanymede/2:1|LGanymede>Callisto/0s|LCallisto>Ganymede/0s`, 10 cycles,
  epoch JD 2456562.889, `--direct --grow-chain --shoot-rel-time`. Tool commit 7d753140 (the run started
  at 16:53:46, after that commit at 16:53:36). Stored: gate indeterminate, worst ratio 0.99170
  (Ganymede), minimum required altitude 130.47 km.
- One cycle = its first four encounters (G after the G-G leg, G after the 2:1 return, C, G).
- Meaning: PASS validates the jup365 lane with a published control, and #316's interior Europa
  failure is specific to #316. Record what differs: #316's Europa periapsis at 982 km required
  altitude against GanCal#1's Ganymede at 130 km. FAIL at an interior node: the chain tool's
  patched-conic real-ephemeris passes do not generally survive continuous gravity at these
  altitudes, and that statement goes into every Jovian row's data_gaps.
- Expected: FAIL with probability about 0.55. GanCal#1's Ganymede flyby is near its turn limit
  (ratio 0.96-0.99, 130 km required altitude), and continuous gravity moves near-limit flybys
  (#968 GanCal#5: +248 km at full mass, though upward there).

## 3. gc-1 at one epoch

- Source: `data/943_gc1_realeph/e2/` (the best of the five #943 rung-(d) epochs: worst ratio 0.75512,
  minimum required altitude 1,675.3 km; epoch JD 2467182.832). Cell gc, key
  `k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s`, 10 cycles,
  5 epochs (index 2). Tool commit 4735f609 (the lead launched at 14:22; the next tool commit, c4ff9a41
  at 14:50, came after). If the reconstruction does not reproduce the stored result, the run stops
  there (an honest NOT RE-RUNNABLE, with the mismatch).
- One cycle = its first four encounters (G, C, C after the 1:1 return, G).
- Meaning: the gc-1 row's real-ephemeris continuous-gravity statement, either way. Run only after the
  control, whatever the control's outcome; the control's outcome conditions how it is read.
- Expected: PASS with probability about 0.5 (moderate turn ratios 0.65-0.76; required altitudes
  1,650-2,440 km in the patched conic; the ideal-model gc-1 closed in continuous gravity, #1034).

## 4. Reconstructions (step 1) and AMENDMENT 1 (before the counted runs)

- GanCal#1 at 2013 (tool 7d753140): max residual 5.6e-8; worst ratio 0.99170 (Ganymede; Callisto
  0.408); minimum required altitude 130.47 km. PASS. `data/1044_gancal1/seed_chain.json`.
- gc-1 e2 (tool 4735f609): max residual 1.6e-9; worst ratio 0.75512 (Callisto; Ganymede 0.661);
  minimum required altitude 1,675.28 km. PASS. `data/1044_gc1/seed_chain.json`.
- The first GanCal#1 start at sigma = 0.02 stalled at 2e4 km. The seed defects did not scale with
  sigma (3.1e4 km at both 0.001 and 0.02). This is the #1043 cause again: the second Ganymede flyby
  (turn 4.0 deg) and the Callisto flyby (0.97 deg) need periapses of 28,327 km and 79,393 km, beyond
  the lane's default 0.6-SOI clamp (19,030 and 30,084 km).
- AMENDMENT 1: the seed periapsis nodes are unclamped (`SEED_CAP_KM = 1e12`, as in #1043). This is a
  seed construction fix only; criteria unchanged. The rung-(b) GanEur#316 runs were not affected:
  their periapses lie inside the clamp.

### 4.1 AMENDMENT 2 (solver only; criteria unchanged; before the counted runs)

- With amendment 1 both chains start at sigma = 0.02 with seed defects that scale with sigma
  (GanCal#1 up to 3,464 km and 0.0085 km/s; gc-1 up to 1,867 km and 0.013 km/s). But the damped
  line-search Newton crawls on both. A full step RAISES the residual: GanCal#1 4.0e3 -> 6.6e5 at a
  full step and 5.8e3 at 0.1; gc-1 2.4e3 -> 1.1e4 at a full step and 2.4e3 at 0.25.
- What drives it:
  - GanCal#1: its weak flybys (G2 4.0 deg, Callisto 1.0 deg; patched-conic periapses 28,000 and
    79,000 km) make the Newton step ask for node moves larger than the nodes' moon distances.
  - gc-1: the steps ask for epoch moves of 500-840 s at its Callisto nodes, on either side of the 1:1
    full-revolution return.
  - Both objects carry a FULL-revolution resonant leg (2:1 and 1:1); GanEur#316, which converged
    directly, has a half-revolution one.
- AMENDMENT 2: Levenberg-Marquardt (scipy `least_squares` method lm, analytic Jacobian, checkpointed
  so that a call can stop and the next resume), then the damped-Newton polish.
- A diagnostic LM run on gc-1 at sigma = 0.02 took the residual from 2.4e3 to 0.17 in 30 evaluations
  (410 s at machine load 42). The continuation reuses that state.

## 5. Result: neither chain STARTS at sigma = 0.02. No verdict on either; time box reached.

- gc-1 e2 (LM, then Newton): the residual fell from 2.4e3 to 0.0155, and then stalled over two
  further calls. It sits at a stationary point that is not a root:
  - the residual is concentrated on the VELOCITY rows of leg 2, the 16.7-day Callisto-Callisto 1:1
    full-revolution return (up to 1.4e-5 km/s);
  - it lies along the weakest left singular vector of the column-scaled Jacobian (singular value
    1.5e-6, against 3.8e-5 for the next);
  - every other leg, gauge and end row is at or near its floor.
  The full-revolution return is the near-degenerate part of the problem: a 1:1 return leaves a
  one-parameter family of resonant directions, and the #943 chain tool had to shoot that leg
  separately for the same reason. Reading (INFERRED): the periapsis-node, mid-leg-match formulation is
  ill-conditioned across a full-revolution return. This is not a property of gc-1.
- GanCal#1 at 2013: the Newton step at the start asks for node moves larger than the nodes' moon
  distances at its weak flybys (G2 4.0 deg, Callisto 1.0 deg), on either side of its 2:1
  full-revolution return. The line search crawls; a 200-evaluation LM diagnostic exceeded the 8-minute
  call. No converged start.
- So #1044 gives NO verdict on either object. The jup365 lane is still not validated by a published
  control. The gc-1 row's real-ephemeris continuous-gravity standing is UNDECIDED, limited by the
  formulation rather than tested.
- What the three real-ephemeris attempts now share: the closed-chain EGGIE (no full-revolution leg)
  passes under (b); GanEur#316 (a half-revolution leg) converges but loses its interior Europa flyby;
  the two full-revolution objects (GanCal#1, gc-1) do not start. Proposed follow-up, not run: a
  full-revolution-aware formulation. Put a node at the resonant leg's apojove with its own gauge, or
  parametrise the 1:1 / 2:1 return by its V_inf direction on the resonant circle, as the chain tool
  does.
- Suggested data_gaps wording for the gc rows (lead's decision): "real-ephemeris continuous-gravity
  standing untested: the open-chain n-body formulation does not converge across a full-revolution
  return (#1044); the ideal-model continuous-gravity closure stands (#1034)".
