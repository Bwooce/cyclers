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
