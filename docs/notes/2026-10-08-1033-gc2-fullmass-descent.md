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

### 5.1 Run (`data/1033_gc2_descent/`; two calls, 17:22 and 17:26 AEDT)

- Sigma = 1 reconverged at the floors: V_inf 3.507 / 2.9323 / 2.9323 / 3.507. IAS15 6.2e-6 km,
  9.6e-9 km/s.
- The descent converged at the floors to sigma = 0.8208, with almost no change (V_inf 3.5096 /
  2.9221; C r_p / sigma 11,090 km). It then could not step below 0.82072: steps down to a relative
  1e-5 failed, and that last point met only the noise floor (dr 9.9e-3 km).
- The monitor stage crashed on its first attempt: a NaN sigma-derivative made lstsq fail. That is a
  bug in the shared #1025 `monitor_step`, fixed in 95ce7550 (#1025 amendment 1). After the fix it
  ended with six halvings without convergence: a numerical stop.

### 5.2 Why: the orbit flies through Ganymede (diagnostic, `scan` in the scratch dir; DOP853 at
0.0005 d, rtol 1e-12)

- On leg 1 (Callisto 14.33 d to Callisto 35.76 d) the trajectory passes Ganymede at t = 25.05 d at a
  distance of 83 km from Ganymede's CENTRE at sigma = 1 (94 km at sigma 0.8208). Ganymede's radius is
  2,631 km, which is also the softening radius at sigma = 1. The path crosses the moon's interior,
  where the softened force is harmonic. No other pass inside 3 Hill radii.
- The residual is also extremely sensitive to sigma there: a relative change of 1e-6 in sigma moves
  the leg-1 match by 3.5 km. That is why the descent stalls.
- The #1025 driver checks impact only at the flyby NODES (r_p <= sigma R). This pass is not a node,
  so it was missed. By the impact rule of #1034 sec. 4 ("every periapsis above sigma R"), this
  closest approach to Ganymede is a periapsis inside the body.

## 6. Outcome (iii): IMPACT. The starting orbit is not a physical orbit.

- The full-mass G-C-C-G "orbit" from the #1025 control's check 6 passes through Ganymede (83 km from
  the centre). It exists only because the force model softens the moon's interior. It is not a
  second gc-2 branch and not a new object. The #1025 sec. 7 observation is withdrawn: gc-2's "folds
  at 0.164" statement stands, and no addendum to the row is needed.
- Consequence for #1025 (recorded there as amendment 2): every EXISTS member must also pass an
  unscheduled-pass scan at sigma = 1, with no closest approach to either moon inside sigma x R,
  before it counts as EXISTS. I will run that scan when I analyse the batch. #1034's gc-1 had this
  scan and passed it (its Ganymede minima were at 5,599, 280,861 and 957,807 km).
