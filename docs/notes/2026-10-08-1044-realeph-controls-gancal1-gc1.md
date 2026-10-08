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
