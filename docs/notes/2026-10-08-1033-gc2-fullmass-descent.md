# #1033: descending continuation from the full-mass G-C-C-G orbit found by the gc-2 seed

Status: PRE-REGISTRATION (committed before the run). Lead ruling 2026-10-08. No catalogue writes;
nothing is called novel.

## 1. Starting orbit

#1025 control gc-2, check 6 (`docs/notes/2026-10-08-1025-sigma-batch.md` sec. 7): damped Newton at
sigma = 1 from gc-2's unclamped patched-conic seed converges to a symmetric G-C-C-G orbit. V_inf G
3.507, C 2.932; periapses 4,075 km (G) and 11,258 km (C); nodes at 11.66, 14.32, 35.77, 38.43 d;
IAS15 re-fly 6.2e-6 km. Ideal model: R-S circular, Ganymede and Callisto massive, joint sigma.

## 2. Method

`scripts/run_1033_gc2_descent.py`, which reuses the #1025 driver's corrector and checks unchanged:
- reconverge the orbit at sigma = 1 from the seed (the floors are required);
- natural continuation DOWN: sigma divided by 1.2 per step, secant predictor in log sigma, halving
  on failure to a relative step of 1e-5; noise-floor rule for intermediate points; stop at sigma =
  0.005;
- impact: a periapsis at or inside sigma x R ends the run as an impact;
- if the step falls below 1e-5: monitor-variable continuation (the #1025 `monitor_step`, no SVD
  tangent), up to 30 points. If sigma turns back up, it is a fold on the descent. If sigma keeps
  falling, the run resumes natural descent;
- at each point: V_inf per node, r_p / sigma, and the maximum |V_inf - gc-2 patched conic (3.617 /
  3.039)|; where the gc-2 #1025 branches have a point within 2 % in sigma, the V_inf difference to
  the nearest lower-branch and upper-branch point;
- IAS15 re-fly of every half-arc at sigma = 1 and at the lowest sigma reached (< 1e-2 km,
  < 1e-7 km/s).

## 3. Outcomes (fixed now)

- (i) It reaches sigma = 0.005 with every node within 0.005 km/s of gc-2's patched conic, and the
  shift falls monotonically over the last three points: gc-2's generating orbit. Then gc-2 EXISTS at
  full mass on a second branch, and the "folds at 0.164" statement is branch-specific; the row needs
  an addendum (lead, next window).
- (ii) It reaches sigma = 0.005 at a DIFFERENT limit (some node more than 0.005 km/s from gc-2's
  patched conic): a full-mass G-C-C-G object that is not gc-2. I record its sigma -> 0 V_inf and node
  epochs (its patched-conic identity). The #973 screens and the literature checks are a follow-up,
  run only after this result is reported.
- (iii) It folds or impacts before sigma = 0.005: an isola or a disconnected branch.
- Numerical stop (six monitor halvings without convergence): reported as such, no outcome.

## 4. Expectation

The C r_p / sigma at full mass (11,258 km) lies between gc-2's lower (11,356 at sigma 0.16) and upper
(11,007) branches, and V_inf C 2.932 is close to the upper branch (2.945 at 0.16). I expect (iii)
or (i) through the upper branch: the descent may run into the gc-2 fold region from above and turn
there (probability about 0.4 for a fold near sigma 0.16-0.5), or reach gc-2 (about 0.35), or reach
another limit (about 0.25).

## 5. Results

(pending)
