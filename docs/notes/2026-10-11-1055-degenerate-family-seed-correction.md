# #1055: seed correction along the degenerate full-revolution family (Lyapunov-Schmidt)

Status: PRE-REGISTRATION (written and committed before any #1055 corrector run). Lead ruling
2026-10-11. No catalogue writes; nothing is called novel.

## 1. What #1046 left, and what the premise check found (before any design)

#1046 (secs. 6.3-9) located the stall of GanCal#1@2013 and gc-1 e2 on jup365: an out-of-plane
velocity mismatch across the full-revolution leg, along the weakest singular direction. It survives
apojove nodes, an explicit crank and free ends. My #1046 report proposed "a seed corrected with the
out-of-plane rows included". Before designing that, I checked its premise and the scaling. These
are three diagnostics, run 2026-10-11 11:05-11:15 AEDT, with scratch scripts; none is a #1055 run.

1. The seeds are NOT out-of-plane inconsistent. I propagated each seed leg (Kepler, the lane's
   Jupiter GM, the lane's jup365 moon states) from moon + V_inf out to the next flyby epoch:
   - GanCal#1: every leg closes to < 1e-3 km and < 1e-10 km/s;
   - gc-1 e2: 0.03-0.19 km and < 2e-6 km/s, of which the orbit-normal part is < 0.002 km and
     < 3e-9 km/s;
   - GanEur#316: < 1e-3 km.
   The chain tool's Jupiter GM equals the lane's to 8e-9. So the patched-conic seeds are consistent
   in every direction, and there is no out-of-plane seed error to correct. The premise of the
   #1046 proposal is false.
2. The near-null direction is a degeneracy of the patched-conic limit. Seed Jacobian, column-scaled,
   #1044 formulation (b'), GanCal#1 on jup365:

   | sigma | sv_min | next sv | u_min . r | weak-direction Newton step (scaled) |
   |---|---|---|---|---|
   | 0.02 | 1.27e-5 | 6.5e-5 | 4.00 | 3.2e5 |
   | 0.005 | 3.18e-6 | 6.4e-5 | 1.18 | 3.7e5 |
   | 0.001 | 4.85e-7 | 6.3e-5 | 0.248 | 5.1e5 |

   sv_min is proportional to sigma, so at sigma = 0 the direction is exactly singular: a one-parameter
   FAMILY of patched-conic chains. The residual's projection on it is also proportional to sigma. So
   the Newton move needed along the family does NOT shrink as sigma -> 0: it stays O(1). Starting at
   a smaller sigma cannot help.
3. The same diagnostic in the IDEAL model (the #1025 closed chain, R-S circular, coplanar):
   - gc-1 has the same sigma-proportional sv_min (3.2e-6, 9.9e-7, 3.0e-7 at 0.02, 0.005, 0.00125),
     but u_min . r is 1e-13. The forcing along the family vanishes by the coplanar symmetry, which is
     why gc-1 starts in the ideal model (#1034) and not on jup365.
   - gc4-19 (#1025 NUMERICAL STOP at the start) has TWO sigma-proportional singular values (3.7e-6
     and 7.4e-6 at 0.02), with u_min . r = 9.5e-3 (proportional to sigma). The same obstruction, in
     the ideal model, explains why it never started.
   - gc-2 has no sigma-proportional singular value (1.6e-4 at every sigma).

Reading: a full-revolution (resonant) leg leaves a family of patched-conic chains. Continuous
gravity selects members of that family by a solvability condition at first order in sigma: the
residual's projection onto the family's left null vector(s) must vanish. The seed is one member,
chosen by the chain tool through the ephemeris's weak symmetry breaking, and not in general the
member that continuous gravity selects. The correction needed is a move ALONG the family to the
member where the projected residual vanishes. This is the Lyapunov-Schmidt reduction. It is the
corrected seed #1055 builds, in place of "out-of-plane rows".

## 2. Seed correction (the method)

Inputs: a formulation's residual and Jacobian (the #1044 open chain for jup365, the #1025 closed
chain for the ideal model), the seed z0 scaled to sigma0 = 0.02, and sigma0.

1. Detection. Compute the column-scaled SVD at sigma0 and at sigma0 / 4. A singular direction is
   DEGENERATE if its singular value falls by more than 2.5 (proportional to sigma would be 4), while
   the next non-degenerate one changes by less than 25 %. q = the number of degenerate directions
   (0, 1 or 2). If q = 0 the correction is the identity (the seed is returned unchanged).
2. Family continuation at fixed sigma0 (q = 1; q = 2 runs one parameter at a time in the order of
   |u . r|, then both refined together). The family parameter is the arclength alpha along the null
   right vector v. At each step:
   - the predictor is z + h v;
   - the corrector is Gauss-Newton on the bordered system: all residual rows except their projection
     onto the left null vector u (the rows are multiplied by (I - u u^T)), plus one arclength row
     v . (z - z_pred) = 0;
   - u and v are recomputed at the converged point.
   Monitored: g(alpha) = u . r(z(alpha)) and the complement residual |(I - u u^T) r|.
   Step h: start at 1 % of |Newton step along v| (from the table), double after a convergence in
   4 iterations or fewer, halve on failure, at least 1e-4 of it. Both directions (+v and -v) are
   scanned, up to 3 times the Newton-step length each way.
3. Root: at the first sign change of g, a secant on alpha (bordered corrector at each) until
   |g| < 1e-3 |g(0)|. The corrected seed is z* = z(alpha*).
4. Hand-off: plain damped Newton (#1044 polish) on the full system at sigma0 from z*. Then the
   object's own continuation (#1044 for jup365, #1025 for the ideal model) from that start.
- The patched-conic member that the corrected seed represents is recorded: z* rescaled to
  sigma -> 0 (offsets / sigma0 x sigma) and its V_inf per node. Identity is judged against this
  member (sec. 4), and its distance from the chain tool's member is reported (V_inf per node,
  crank angle).

Code: `src/cyclerfinder/nbody/degenerate_family.py` (additive, pinned by fast tests at production
scales: km, km/s, the Jovian GM), driver `scripts/run_1055_seedfix.py` (task 1055, preflight).
Tests (fast):
- an analytic two-parameter toy with a known family and a known solvability root (a
  sigma-proportional singular value), in km and km/s magnitudes: the root found to 1e-8;
- the detection returns q = 0 on a well-conditioned toy;
- the bordered corrector keeps the arclength row to 1e-12.

## 3. Controls and targets, in order

1. Positive control A (ideal, q = 1, zero forcing): gc-1 in the #1025 closed chain. Expected: q = 1
   detected, g identically below 1e-9, so alpha* = 0. The corrected start must reproduce the #1025 /
   #1034 sigma = 0.02 point (node V_inf within 1e-4 km/s), and the continuation must reach sigma = 1
   with full-mass V_inf 2.130 / 1.785 / 1.843 (#1034) within 2e-3 km/s.
2. Positive control B (the lane already closes it, as the lead asked): GanEur#316 one cycle on
   jup365 under #1039 (b') (it starts at sigma 0.02 from the uncorrected seed, #1039). Expected:
   q = 0, or q >= 1 with alpha* close to 0. The corrected seed must still start, at the same
   sigma-0.02 point (node V_inf within 1e-4 km/s). (GanCal#5 and GanEur#43 were closed by corrector
   A in the ideal model, a different code path; GanEur#316 is the control the open-chain lane
   closes. My first GanCal#5 build attempt through the #1025 builder gave a seed residual of 1e10:
   it is the wrong candidate match, and I do not use it.)
3. Negative control (the mechanism): for each target whose corrected seed starts, a start from
   z(alpha* + 0.5 L) and from z(alpha* - 0.5 L), where L is the Newton-step length along v: the same
   family moved off the root. Both must FAIL to start within the same budget (40 Newton iterations).
   If either starts, the start does not depend on the correction, and the correction is not shown
   to be the cause. This guards against "it closed!".
4. Ideal-model target: gc4-19 (#1025 NUMERICAL STOP at the start, q = 2 expected). Pass: the
   corrected seed starts, identity holds, and the #1025 continuation then gives an outcome
   (EXISTS / FOLDS / IMPACT) under #1025's rules, including the amendment-2 scan.
5. jup365 targets: GanCal#1@2013, then gc-1 e2, #1044 formulation (b') with the #1044 criteria
   (i)-(iv), the amendment-2 unscheduled-pass scan, and identity against the corrected member.

## 4. Pass criteria per target

- START: Newton from z* converges at sigma0 to the lane floors (1e-3 km, 1e-6 km/s, gauges 1e-9,
  ends 1e-6).
- IDENTITY: converge going down at sigma 0.01 and 0.005. Every node V_inf within 0.01 km/s of the
  corrected patched-conic member (sec. 2), and within 0.05 km/s of the chain tool's chain (reported
  either way). The family move may change V_inf slightly, and the 0.05 band is the bound for it to
  still count as the same chain.
- The object's continuation and end checks as in #1044 / #1025.
- An object whose corrected seed starts and passes identity, then reaches sigma = 1 with every check
  passing: the jup365 continuous-gravity lane is validated on a full-revolution chain (GanCal#1 is
  the published control), and gc-1's real-ephemeris continuous-gravity standing is decided (EXISTS),
  with its corrected-member geometry recorded.

## 5. Outcomes and meanings (fixed now)

- (a) Control A and B pass, the negative controls fail as required, and a target starts: the #1046
  obstruction was the family selection, and the method fixes it. The target's continuation result
  is its verdict.
- (b) The controls pass, but g has no sign change within 3 Newton-step lengths for a target: no
  solvable member near the seed. No verdict on that object. Report g(alpha) over the scanned range.
- (c) g has a root, but the start from z* fails: the reduction is not enough (second-order terms).
  No verdict. Report.
- (d) A control fails: stop and report. No target runs.

## 6. Expected runtimes (measured cost per residual-plus-Jacobian, 2026-10-11 at load 11-15)

- Ideal model (gc-1, gc4-19): about 1 s per evaluation. The family scan plus root plus start take
  about 50-150 evaluations, so 2-4 minutes. The continuation is #1025's (2-12 minutes).
- GanEur#316 and GanCal#1 on jup365: 3-17 s per evaluation, so 5-40 minutes for the family stage.
  The stage is checkpointed per family point. Calls of up to 8 minutes; if more than 3 calls are
  needed, I hand the exact command to the lead.
- gc-1 e2: 35-110 s per evaluation, so 1-4 hours for the family stage. That is the lead's launch: I
  give the command.

## 7. Expectations (fixed now)

- Control A passes: 0.9. Control B passes: 0.85.
- gc4-19 (ideal) starts after correction: 0.6. A root exists: 0.75.
- GanCal#1 starts: 0.45.
- gc-1 e2 starts: 0.5.
- A target that starts reaches sigma = 1: about 0.6 each.

## 8. Results

(pending)

### Amendment 1 (implementation, before any run)

- The null vectors u and v are taken at the BASE point (the seed with its complement converged) and
  held fixed along the scan, instead of being recomputed at each point (sec. 2 item 2). With fixed
  vectors the bordered system and g are the same function along the whole scan, so a sign change of
  g is meaningful. The family is traced as a graph over the fixed arclength coordinate, which is
  valid while it does not fold with respect to that coordinate within the scanned range. The
  predictor is the secant along the traced family.
- Root at base: if |g(0)| is already within the complement tolerance, the seed is returned
  unchanged (positive control A is expected to take this path: its g is 1e-13).
- q >= 2 (gc4-19): damped Newton on the q reduced equations, with the Jacobian of g by forward
  differences (1 % of L) and steps capped at L, instead of coordinate scans followed by a refine.
- Tests (`tests/nbody/test_1055_degenerate_family.py`, 7 fast tests, km and km/s magnitudes, the
  lane's 1e3 velocity weight). Toys whose roots are set by construction:
  - detection gives q = 1 and q = 2, and 0 for a well-conditioned toy;
  - the bordered corrector holds the arclength row to 1e-12;
  - the 1-d root is within 1e-3 of the constructed root (4e5 km, O(1) from the seed), and plain
    Newton from it reaches the exact root to 1e-9;
  - the 2-d root is within 2e-3;
  - a seed already at the root is returned unchanged.

### Amendment 2 (identity wording, before any run)

The corrected member's patched-conic V_inf is not known independently: it is the sigma -> 0 limit
of the corrected branch. So "within 0.01 km/s of the corrected member" (sec. 4) would be circular.
Identity is instead:
- (1) the branch converges going down at sigma 0.01, 0.005 and 0.0025 (the #1025 schedule for the
  ideal model; the #1044 `stage_down` 0.01 / 0.005 for jup365), with the successive V_inf
  differences shrinking (a Cauchy limit);
- (2) the Richardson limit V(0) = 2 V(sigma_last) - V(2 sigma_last) is within 0.05 km/s of the chain
  tool's V_inf at every node.
I report the #1025 / #1044 identity numbers against the chain as they come out. Control A, gc-1
with alpha* = 0, must also pass #1025's own identity gate unchanged.

Driver: `scripts/run_1055_seedfix.py` (stages fix, start, neg, cont, down, ias15).

### Amendment 3 (start criterion on jup365, before any counted start; 11:42 AEDT)

On jup365 the START at sigma0 accepts the noise-floor rule (dr < 1e-2 km, dv < 1e-6 km/s, ends
< 1e-6). That is the rule rungb's own start uses (#1039 amendment 11, #1044), and it is the rule
under which the uncorrected GanEur#316 started in #1039, so control B is judged by the same
criterion as the original. The first control-B start call (11:33) hit the 470 s limit at dr
2.0e-3 km, dv 3.3e-8 km/s under a machine load of 35-146 (a CI runner with 6 workers). Verdict
points at sigma = 1 still need the floors. The ideal-model start (#1025 closed chain) keeps the
floors: control A met them.

### Results so far

- Control A (gc-1, ideal): PASS.
  - Detection q = 1: sv_min 3.18e-6 at sigma 0.02 and 9.88e-7 at 0.005.
  - g(0) = 7.0e-14, so the root is at the base (alpha* = 0) and the seed is unchanged.
  - The start converges at the floors, and its V_inf 2.38775 / 1.80691 / 1.80797 / 2.38774 equals
    the #1025 gc-1 sigma-0.02 point to 1e-5.
  - The #1025 continuation repeats exactly: identity shifts 0.0095 / 0.0052 / 0.0028 / 0.0015, sigma 1
    reached, full-mass V_inf 2.1298 / 1.7845 / 1.8428 / 2.1294 (#1034), IAS15 pass.
  - `data/1055_gc-1/`.
- Control B (GanEur#316, jup365): detection q = 0 (sv_min 1.53e-5 at 0.02, 1.25e-5 at 0.005), so the
  corrected seed is the seed. The start runs under amendment 3.
- Control B (GanEur#316, jup365): PASS. q = 0, so the seed is unchanged. The start (amendment 3)
  converges at dr 3.0e-3 km and dv 5.1e-8 km/s, with V_inf 3.19029 / 3.19107 / 3.19554 / 3.77569 /
  3.17503, equal to the #1039 (b') sigma-0.02 start to 1e-5. (Three calls: the first two hit the
  470 s limit under machine load 35-146, then amendment 3's per-iteration checkpoint resumed.)

### Amendment 4: implementation BUG found by the first gc4-19 run, and fixed (12:00-12:25 AEDT)

- Run (12:00): gc4-19, detection q = 2 (sv 7.37e-6 and 3.70e-6 at 0.02; 1.84e-6 and 1.00e-6 at
  0.005). The q-d Newton moved 4 steps of L = 2.6e3, and |g| fell only from 8.15e-3 to 7.72e-3:
  "no root within the scanned range".
- Cause: the code took the null vectors, the scales and L from the RAW seed. Amendment 1 registered
  them at the BASE (the seed with its complement converged). The raw seed is far from the family
  (|r| about 3.7e3), so its vectors and L do not describe the family the corrector then traces.
  Fixed in `family_frame`: converge the complement with the seed's vectors, then re-linearise at
  that point and converge again with the base's own vectors; the arclength is measured from the
  base. The driver's negative-control stage uses the same frame.
- Test toys: with base vectors, the original toy's curvature (x shifted by 1e-7 w^2 km, which is
  16,000 km at the root) put the root 10 L away in the fixed-v coordinate. The toy curvature is now
  1e-9 w^2 (160 km at the root), and the family row's vy coupling vanishes on the family, so the
  root is exactly 4e5 km by construction. Added test: an off-family seed gives a converged base
  with its residual along u. 8 tests pass, and `mypy src tests` is clean.
- Re-runs: gc4-19 fix again (its first run is VOID). gc-1's fix is re-run for the record; its
  result cannot change (g = 7e-14, root at base). Control B is q = 0 and does not use the frame.
