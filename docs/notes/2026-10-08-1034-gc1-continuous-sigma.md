# #1034: does gc-1 exist in continuous gravity? Joint-sigma continuation

Status: PRE-REGISTRATION (written and committed before any #1034 run). Gating item for the
#942/#943 writeback (lead ruling 2026-10-08). No catalogue writes; nothing is called novel (gc-1's
owner ruling, candidate-novel, stands or falls with the owner, not here).

## 1. Question

gc-1 (k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s; V_inf G
2.397, C 1.807 km/s; Ganymede 2 x 29.35 deg, Callisto 2 x 40.15 deg) passed the patched-conic
real-ephemeris chain (generator note 6.26). #1004 showed that a patched-conic pass does not imply
existence in continuous gravity: gc-2's branch folds at sigma = 0.164 of the real moon masses. Does
gc-1's branch reach the real masses?

## 2. Model and method (`scripts/run_1034_gc1.py`, reusing `scripts/run_1004_gc2.py`)

- The #968-validated lane model (ideal continuous R-S model; #968 note sec. 8): R-S 2009 Table 2
  circular coplanar ephemeris, Ganymede and Callisto point masses with the indirect term, DOP853 +
  analytic STM, GM and softening radius both scaled by sigma (sigma = 1: real GMs and radii).
- Corrector: the #1004 forward-backward shooter, four massive-body nodes in time order
  G (11.76 d), C (18.27 d), C (34.96 d), G (38.75 d) and the rotated wrap; periapsis gauge at each
  node; square 28 x 28.
- Seed: the patched-conic gc-1 (`data/943_cell_gc_gauntlet.json`), nodes from `periapsis_node`,
  moon-relative offsets scaled by sigma.
- Identity check at sigma = 0.05: V_inf G within 0.02 km/s of 2.397 and C within 0.02 of 1.807.
  If it fails, the run stops there (wrong branch or no generating orbit).
- Continuation: natural parameter sigma upward (secant predictor in log sigma, damped Newton at
  fixed sigma) to sigma = 1. Step factor 1.2, halved on failure down to a 1e-5 absolute step.
- Tangent rule, ENFORCED in code: a pseudo-arclength step on an SVD tangent is taken only if the
  singular-value gap (second smallest / smallest of the column-scaled [J_z | J_s]) exceeds 10; in
  #1004 the gap was about 1, because of weak J_z directions of tiny physical size (#1004 note sec. 9).
  When the gap test fails, a fold is passed by MONITOR-VARIABLE continuation instead: sigma becomes
  an unknown and the state component that changed most over the last secant (in column-scaled
  units) is stepped and fixed by an extra row. No step ever uses an ambiguous tangent.

## 3. Checks if natural continuation stops before sigma = 1 (amendment-1 checks of #1004)

1. Fold bracket: refine the step to 1e-5 on the branch; record the last converged sigma.
2. Sibling: pass the turn by monitor-variable continuation; then two solutions at one sigma (two
   values below the fold): converge at FIXED sigma from each branch; the fold is supported if both
   converge and differ by more than 0.02 km/s in V_inf C and 100 km in Callisto r_p / sigma.
3. Direct Newton at sigma = 1 from the patched-conic seed and from each branch's last point.

## 4. End-point checks

- IAS15 re-fly (second integrator, scaled GMs as in #1004 check 3) of every half-arc at the end point
  (sigma = 1 if reached, else the two sigma-below-fold solutions): agreement < 1e-2 km, 1e-7 km/s.
- Impact: every accepted point must have each periapsis above that body's scaled radius
  (sigma x R); a branch that reaches a surface stops there (an impact end).
- At sigma = 1 (if reached) record V_inf at each node, altitudes against R-S radii (Ganymede 2634,
  Callisto 2408 km) and the project floors (100 and 200 km), turns, and r range.

## 5. Stopping rules

sigma = 1 reached; or a fold confirmed by checks 1-2; or an impact; or six halvings without
convergence (a NUMERICAL stop, not a physical end). Each call under 8 minutes, checkpointed in
`data/1034_gc1/`; live log in `data/1034_gc1/live/` (gitignored).

## 6. What each outcome means (fixed now)

- Reaches sigma = 1 on the branch that passed the identity check: gc-1 EXISTS in the continuous ideal
  model; its full-mass V_inf, altitudes and turns are recorded. It is judged by closure and by
  continuation to its own patched-conic limit (the #968 sec. 8.2 rule), not by agreement with its
  patched-conic numbers at full mass. (Real ephemeris is a separate rung.)
- Folds before sigma = 1 (checks 1-2 pass): gc-1 is a patched-conic object only, on this branch;
  existence elsewhere (isola, other seed) untested.
- Impact end: as a fold for the writeback: no physical full-mass orbit on this branch.
- Numerical stop: undecided; reported with the last point and the failure mode.

## 7. Expected outcome (stated now)

Uncertain. gc-1 has a lower Callisto V_inf (1.81 km/s) than gc-2 (3.04) and a full-revolution 1:1
Callisto return, both of which make continuous Callisto gravity matter more. My expectation: a fold
before sigma = 1 (probability about 0.6).

## 8. Results

(pending)

### 8.1 Identity check at sigma = 0.05: FAIL as registered (by 0.0007 km/s)

Converged at sigma = 0.05 from the patched-conic seed: V_inf G 2.3763 / 2.3763, C 1.8070 / 1.8097
(the two Callisto nodes differ by 0.0027). Ganymede is 0.0207 km/s from 2.397, against a 0.02
tolerance; Callisto passes. By sec. 2 the run stopped. The likely reason is the finite-mass shift
already at sigma = 0.05 (GanCal#5 shifted 0.013 km/s at s = 0.05 with V_inf 3.24; gc-1's V_inf is
lower), but that is a reading, not a check.

### 8.2 AMENDMENT 1 (before the runs it covers)

Identity by the limit, as in #968 sec. 8.1 (criterion 4c style): converge at sigma = 0.02, 0.01 and
0.005 (from the patched-conic seed, offsets scaled, then each from the previous). Pass: at the last
computed point (sigma = 0.005) V_inf G within 0.005 km/s of 2.3972 and C within 0.005 of 1.8067
(the #943 patched-conic values), and the G shift falls monotonically with sigma. If it passes, the
upward continuation starts from the sigma = 0.05 point as registered; if not, the run stops.

## 9. Results (`data/1034_gc1/`)

- Identity by the limit (amendment 1, `identity.json`): sigma = 0.02, 0.01, 0.005 gave V_inf G 2.3877,
  2.3921, 2.3944 and C 1.8069/1.8080, 1.8068/1.8074, 1.8068/1.8070. At sigma = 0.005, G is 0.0028
  and C 0.0003 from the patched-conic 2.3972 / 1.8067 (tolerance 0.005), and the G shift falls
  monotonically: PASS. The branch is gc-1's.
- Natural continuation from sigma = 0.05 (`sigma.json`): every step converged at factor 1.2, no
  halving, no impact, through sigma = 0.06, 0.072, ..., 0.770, 0.924, 1. No fold: the monitor-variable
  stage and the fold checks of sec. 3 were not needed.
- **gc-1 REACHES sigma = 1: it exists in the continuous ideal model** (R-S constants, both moons point
  masses with real GMs and radii). Full-mass orbit:

| Node | t (d) | V_inf (km/s) | periapsis (km) | altitude (R-S radius) | turn (deg) |
|---|---|---|---|---|---|
| Ganymede | 11.7 | 2.1298 | 5598.3 | 2964 km | 32.55 |
| Callisto | 18.5 | 1.7845 | 4373.0 | 1965 km | 39.77 |
| Callisto | 35.1 | 1.8428 | 5249.7 | 2842 km | 33.37 |
| Ganymede | 38.8 | 2.1294 | 5595.9 | 2962 km | 32.57 |

  r from 916,320 to 2,294,978 km; the C-C arc at mid-time has a 1,881,915 km, e 0.2199. All
  altitudes are far above the floors (100 / 200 km). The two Callisto flybys are no longer
  symmetric. Against the patched conic (V_inf G 2.397, C 1.807; 2,437 / 1,801 km required altitudes)
  the full-mass shift is -0.27 km/s at Ganymede, the same order as #968's controls; per #968 sec. 8.2
  that is not a failure criterion.
- Sequence (descriptive, DOP853 over one period at 0.002 d): Ganymede distance minima at the node
  (5,599 km) and at 280,861 and 957,807 km; Callisto minima only at the two nodes. No unscheduled
  pass inside either Hill radius.
- IAS15 re-fly of the eight half-arcs at sigma = 1: max 2.4e-6 km, 1.1e-11 km/s: PASS. (Full-period
  DOP853 single shot from the first node: 17 km return miss, descriptive; flyby amplification.)
- Tangent rule: no pseudo-arclength step was taken, so the gap test never applied.

Meaning (fixed in sec. 6): gc-1 EXISTS in the continuous ideal model on its own branch, with
full-mass V_inf 2.130 (G) and 1.785 / 1.843 (C). Real ephemeris in continuous gravity is a separate
rung, not run.
