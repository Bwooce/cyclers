# #1055: seed correction along the degenerate full-revolution family (Lyapunov-Schmidt)

Status: PRE-REGISTRATION (written and committed before any #1055 corrector run). Lead ruling
2026-10-11. No catalogue writes; nothing is called novel.

How to read this driver's results: judge each stage by its output file in `data/1055_<target>/`,
never by the exit code. fix -> `fix.json` (a paused call writes only `fix_state.json` and exits 0),
start -> `start.json` (`converged`), neg -> `neg.json` (`started` must be false), ideal cont ->
`<target>.json` (`outcome`), jup365 cont -> `sigma_n1.json`, limit -> `limit.json`. A stage that
crashes leaves its file absent or unchanged.

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

Summary at 12:55 AEDT. The details are in the dated subsections below, which are interleaved with
the amendments in the order they happened.

| Item | Pre-registered criterion | Result |
|---|---|---|
| Control A, gc-1 ideal | q = 1, alpha* = 0; start V_inf equals #1025 at sigma 0.02 to 1e-4; sigma = 1 with V_inf equal to #1034 to 2e-3 | PASS. q = 1, g(0) = 7.0e-14, root at the base, seed unchanged. Start at the floors, V_inf 2.38775 / 1.80691 / 1.80797 / 2.38774 (equal to 1e-5). #1025 identity 0.0095 / 0.0052 / 0.0028 / 0.0015 passes. Sigma = 1: 2.1298 / 1.7845 / 1.8428 / 2.1294 (equal). IAS15 passes. |
| Control B, GanEur#316 jup365 (b') | q = 0, or a root near 0; the start at the same sigma-0.02 point to 1e-4 | PASS. q = 0 (sv_min 1.53e-5 at 0.02, 1.25e-5 at 0.005), seed unchanged. Start (amendment 3 rule): dr 3.0e-3 km, dv 5.1e-8 km/s, V_inf 3.19029 / 3.19107 / 3.19554 / 3.77569 / 3.17503 (equal to #1039 (b') to 1e-5). |
| Negative control (mechanism), on gc4-19 | every start from the family 0.5 L off the root fails (40 iterations) | PASS. 4 of 4 fail: dv 5.3e-4, 5.2e-4, 9.4e-6, 5.6e-6 km/s (above 1e-6). |
| gc4-19 ideal, start | START at the floors; identity (amendment 2) | STARTS. q = 2, root at 1.0 L (|g| 7.7e-3 -> 5.3e-6). dr 8.3e-4 km, dv 2.3e-9 km/s. Identity shifts 0.0009 / 0.0004 / 0.0002 / 0.0001 against the chain tool's V_inf, so the limit is the chain's own V_inf. PASS. |
| gc4-19 continuation | #1025 outcome rules plus the amendment-2 scan | RUNNING (lead launch 12:37). The first step, sigma 0.024, converged at 12:44 after about 30 min at machine load about 45. |
| GanCal#1@2013 | fix, start, negative control | RUNNING (lead launch 12:37). Detection q = 1 (sv 1.27e-5 -> 3.18e-6). |
| gc-1 e2 | fix, start, negative control | Queued after GanCal#1 (lead). |

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

### Amendment 5: noise floor of the family stage (12:20 AEDT, before the gc4-19 re-run)

The diagnostic at gc4-19's base (scratch, 12:05-12:15) showed the following:
- With the seed's vectors, the complement converges to 6.3e-4 (weighted), and the whole residual is
  then 8.2e-3, along u: g = (-1.2e-3, 8.1e-3).
- Re-linearised there, the complement stalls at 1.45e-3. That is integration noise at sigma 0.02,
  the level the #968 noise-floor rule was made for.
- At the base the smallest singular values are 1.3e-5, 6.7e-6 and 3.9e-8.

Changes:
1. `bordered_solve` projects every iterate exactly onto the linear arclength rows, and its line
   search acts on the complement alone. Before, a reduced line-search step left the arclength rows
   just above their 1e-10 test, so a converged complement was reported as a failure.
2. Complement tolerance TOL = 3e-3 (weighted), above that noise floor. Root criterion: |g| <
   max(1e-3 |g(0)|, 1e-3), i.e. below the start's noise-floor rule in weighted units (1e-3 km,
   1e-6 km/s). The pre-registered 1e-3 |g(0)| alone would ask for g below the noise when g(0) is
   already small. The START criteria (sec. 4, amendment 3) are unchanged.
The 8 tests pass.

### gc4-19 (ideal model), after amendments 4-5 (12:06-12:24 AEDT)

- Fix: q = 2. Base |g(0)| 7.7e-3 (weighted), L = 1.96e5 (scaled). Two Newton steps on the reduced
  equations give |g| 5.3e-6: root at |alpha*| = 1.92e5, about 1 L from the base, so an O(1) move
  along the family, as sec. 1 predicted. Corrected-seed V_inf 4.1099 / 4.1082 / 2.3562 / 2.3564
  (raw seed at sigma 0.02: 4.109 / 4.109 / 2.3565 / 2.3565).
- START: converges at the FLOORS at sigma 0.02 (dr 8.3e-4 km, dv 2.3e-9 km/s). #1025 never got a
  start for gc4-19.
- NEGATIVE CONTROL: from the family moved 0.5 L off the root along each null direction (4 points),
  Newton fails every time (dv 5.3e-4, 5.2e-4, 9.4e-6 and 5.6e-6 km/s, all above 1e-6). The start
  depends on the correction: PASS.
- IDENTITY (#1025 gate): shifts 0.0009 / 0.0004 / 0.0002 / 0.0001 km/s at sigma 0.02 / 0.01 / 0.005 /
  0.0025, monotone, PASS. The corrected member's limit equals the chain tool's patched-conic V_inf
  to 1e-4 (amendment 2's Richardson limit is within 0.05 by a wide margin). The family move changes
  the geometry along the degenerate direction, not the V_inf.
- Continuation (#1025 rules) is slow: a factor-1.2 step does not converge within the 25 + 15
  budget, and a 2.5 % step needs 8-12 Newton iterations (diagnostic, scratch). The run needs more
  than three 8-minute calls, so it goes to the lead (`--max-wall`).

### Amendment 6: resumable family stage (12:45 AEDT; lead ruling: under-8-minute calls, checkpointed)

`family_root_1d` keeps its whole state (frame vectors, scan position, bracket, secant) in a
JSON-serialisable dict, saved after every step. `bordered_solve` takes a deadline and pauses after
at least one step. Every call therefore makes progress, even one started past its deadline. The
driver keeps the detection and the family state in `fix_state.json`, so `--stage fix` resumes. A
test runs the toy with an already-past deadline (one step per call, state through JSON) and
reaches the same root as an uninterrupted run: 9 tests pass. The method and the criteria are
unchanged. The first GanCal#1 fix call (12:28-12:36, pre-amendment code) reached detection, q = 1
(sv 1.27e-5 at 0.02, 3.18e-6 at 0.005), then hit the call limit with no state saved.

Launch overlap (recorded): the lead launched gc4-19 cont and the GanCal#1 fix/start/neg loop at
12:36. The GanCal#1 fix process loaded the driver before the amendment-6 edit: it is the
uninterrupted, non-resumable family stage, and it writes `fix.json` when it ends. My own resumable
fix call ran on the same target from 12:39 to 12:44, concurrently with it. It wrote only
`fix_state.json` (paused in phase frame2), which the lead's process does not read. I make no further
GanCal#1 calls while the lead's run is alive. `fix_state.json` stays: a valid partial state of the
same method, used only if the stage is re-run.

## 9. Prepared jup365 continuation commands (run only after a target's start and negative control pass)

From /Users/bruce/dev/cyclers, for `<t>` = gancal1 or gc1e2:

    uv run python scripts/run_1055_seedfix.py --target <t> --stage cont && \
    uv run python scripts/run_1055_seedfix.py --target <t> --stage down && \
    uv run python scripts/run_1055_seedfix.py --target <t> --stage ias15

- cont: writes `sigma_n1.json` with the corrected start as its first point, then rungb's
  `stage_sigma` (natural continuation 0.02 -> 1, factor 1.2, halving to 1e-5, the noise-floor rule
  for intermediate points, the #1044 acceptance at every point, the unscheduled Hill-radius scan at
  sigma = 1). Resumable per point; no time limit of its own.
- down: rungb's `stage_down` (sigma 0.01 and 0.005 from the start). The amendment-2 identity
  (Cauchy limit plus a Richardson limit within 0.05 km/s) is computed from it in the analysis.
- ias15: rungb's `stage_ias15` at sigma = 1, every half-arc, on `JovianRailsCache` with registry GMs
  (#1044 criterion iv).
- After the run I apply the amendment-2 inside-body scan to the sigma = 1 point (at 0.0005-day
  sampling, as #1025 sec. 10).
- Runtime: about 23 natural steps if none fail. Each step costs 3-10 evaluations at 5-40 s
  (GanCal#1) or 35-110 s (gc-1 e2), depending on the machine load. That is about 0.5-2.5 h for
  GanCal#1 and 2-8 h for gc-1 e2, plus the down and ias15 stages (minutes for GanCal#1; up to 1 h
  for gc-1 e2).

## 10. gc4-19 continuation result (lead launch 12:37; finished 12:57 AEDT, exit 0)

Checked against the #1025 outcome rules and #1025 amendment 2 (14:00 AEDT check brought forward):
- Outcome EXISTS. 23 points, factor 1.2, sigma 0.02 -> 1, no passed stops, no impact at any point,
  every node above the scaled floors at every point. Sigma = 1 is at the floors (dr 1.4e-5 km).
- IAS15 re-fly of every half-arc at sigma = 1: max 1.5e-5 km, 5.4e-11 km/s, PASS (< 1e-2 km, < 1e-7
  km/s).
- Amendment-2 scan at sigma = 1 (DOP853, 0.0005 d, both moons, all legs): no minimum within 3 Hill
  radii, so none inside a body. The positive control of the scan again flags the #1033 orbit
  (Ganymede, 83 km, inside the body).
- Identity (#1025 gate, sec. 8): 0.0009 / 0.0004 / 0.0002 / 0.0001 km/s, PASS. Amendment 2's
  Richardson limit equals the chain tool's patched-conic V_inf to 1e-4.

Full mass against the patched conic:

| Node | t (d) | V_inf km/s (patched conic) | altitude km (patched conic) | turn deg (patched conic) |
|---|---|---|---|---|
| Callisto | 19.78 | 4.1445 (4.1090) | 7,263 (6,903) | 4.75 (5.01) |
| Callisto | 36.48 | 4.0616 (4.1090) | 6,730 (6,903) | 5.21 (5.01) |
| Ganymede | 39.26 | 2.3356 (2.3565) | 697 (583) | 41.29 (41.77) |
| Ganymede | 53.58 | 2.3352 (2.3565) | 568 (583) | 42.42 (41.77) |

What the correction changed, at the sigma-0.02 start against the raw seed scaled to 0.02:
- the Ganymede node epochs moved by -600 s and -692 s, and the Callisto node epochs by -141 s and
  -139 s;
- the Ganymede nodes' moon-relative positions moved by 9 km (of 64 km) and their velocities by
  0.17 km/s; the Callisto nodes by 0.5 km and 2.7e-3 km/s.
That is the move along the two degenerate directions. The node V_inf is unchanged in the limit
(identity to 1e-4).

Evidence standing against the 15 EXISTS members of #1025 (c90568bf):
- The same model, corrector, continuation rules, identity gate, IAS15 criterion, floors and
  amendment-2 scan. All pass.
- One difference: the start is the #1055-corrected seed, not the raw seed. The correction method
  passed its pre-registered controls (A and B), and the negative control shows the start depends on
  the correction.
- The sigma -> 0 limit of the branch is a member of gc4-19's degenerate patched-conic family, with
  the chain tool's V_inf at every node (to 1e-4). Its node epochs differ from the chain tool's
  member by up to about 11.5 min, which is the family move. That member was not separately run
  through the #973 patched-conic gate: its V_inf is the same, and its full-mass altitudes are far
  above the floors.

Statement: gc4-19 EXISTS in the continuous ideal model (R-S constants, Ganymede and Callisto
massive, real GMs), with the same checks as the 15 #1025 EXISTS members. It is the 16th of the 20
clean #973 gc members to do so. Caveat for the row: it was reached from the #1055 corrected seed (a
different member of the same degenerate patched-conic family; same V_inf). Nothing here is a
novelty claim, and no catalogue write is made.

## 11. GanCal#1 fix stage CRASHED (lead launch 12:37-12:58): diagnosis and AMENDMENT 7

What ran: the lead's process loaded the pre-amendment-6 driver, i.e. a non-resumable fix call with
no state saved. It used the module's fixed-vector 1-d scan (amendments 1 and 4). From its runlog:
- Base: L = 4.59e5 (scaled), g(0) = -4.03 (weighted).
- The + direction: g = -3.99, -3.91, -3.76, -3.48, -2.99, -2.28, -1.85, -1.71 at alpha = 4.6e3 ... 5.09e5
  (1.1 L). The complements converged at every point (4e-4 to 2.5e-3).
- At the next predictor, `propagate_with_stm` raised "Required step size is less than spacing
  between numbers" inside `bordered_solve`'s FIRST residual evaluation. That evaluation was not
  guarded, so the whole stage died. The - direction was never scanned, and no fix.json was written.
  start and neg then stopped on the missing file. (Python exits 1 on a traceback; I checked that
  `uv run python -c "raise RuntimeError"` returns 1. The "exit 0" in the lead's loop is the loop
  wrapper's own status, so judge by the files.)

Diagnosis:
1. The degenerate direction is the WEAK Ganymede flyby, G2 (turn 4.0 deg). At the base, the
   family's right null vector moves G2's moon-relative periapsis by 512 km per 5.09e5 of alpha,
   against an offset of 432 km at sigma 0.02. It moves G1 by 13 km, the Callisto flyby by 184 km,
   and G3 by 0.1 km. The family is the position of G2's periapsis around Ganymede at nearly fixed
   V_inf: a CIRCLE-like family. A fixed-vector (chord) scan cuts across it towards the moon. The
   linear prediction puts G2's offset at a minimum of 96 km near alpha = 4.2e5. That is 1.8 times
   the sigma-scaled radius (52.6 km) and 1.75 times radius plus floor (54.6 km). So the scan was
   skimming Ganymede when the integration failed. INFERRED: the crash point was not recorded, and
   the nonlinear path may have passed closer.
2. g rose monotonically, with a slope falling by about 0.6 per 1.47e5 of alpha. On the chord, an
   extrapolation gives no zero before about 3 L. That is an extrapolation along the wrong (straight)
   path, not a result.
3. So the pre-registered procedure did NOT complete: + reached 1.1 L of 3 L and crashed, and - was
   never scanned. GanCal#1 is NOT a clean negative under sec. 5 (b). It is a numerical failure of
   the fix stage, caused by the fixed-vector scan on a curved family plus an unguarded predictor
   evaluation.

AMENDMENT 7 (method changes; registered before any further GanCal#1 or gc-1 e2 run):
- (a) Curved-family tracing. The scan is pseudo-arclength: at each converged point the local null
  pair (u, v) is recomputed (scales held at the base) and oriented continuously. The next step's
  arclength row is v . (d (z - z_last)) = sgn h about the last point, and the secant predictor is
  along the traced family. g on a step is the step's own u . r, so a sign change is detected on one
  continuous function per segment, and the secant refines within that segment in that segment's
  frame. This restores the pre-registration's own sec. 2 ("u and v recomputed at the converged
  point"), which amendment 1 had replaced with fixed vectors. Fixed vectors were adequate for the
  gentle families of the controls and of gc4-19, but not for GanCal#1's circle-like one.
- (b) An unintegrable predictor (a `RuntimeError` in the corrector's first evaluation) counts as a
  failed step, and the step is halved. It no longer crashes the stage.
- (c) Domain bound: a converged family point with any flyby node's moon-relative distance below
  sigma0 (R + floor) (Ganymede 100 km, Callisto 200 km, Europa 100 km floors) ends the scan in that
  direction ("domain"). A root beyond it would be a below-floor (impact-side) member, outside the
  physical domain.
- (d) Both directions are always scanned (+ then -) to 3 L of path length or a bound, unchanged.
- Tests: 11 pass. New: a circle family (R = 1e5 km; root 115 deg round from the seed, which no
  chord in a fixed coordinate reaches) is followed to its root to 1e-3 rad; a domain bound ends a
  scan cleanly with no root. gc4-19's result (sec. 10) used the fixed-vector scan, which found
  its root at 1 L on a gentle family. It does not depend on (a), and it is not re-run.

Corrected command (lead launch; resumable, each call continues from fix_state.json):

    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage fix --max-wall 3600

Repeat until `data/1055_gancal1/fix.json` exists (a paused call logs "fix: paused" and exits 0).
Then:

    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage start && \
    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage neg

`fix_state.json` holds my 12:39 partial base frame (phase frame2, the same method), which the
first call resumes. The outcome rules are sec. 5: a root, then start and negative control; or "no
sign change within 3 L" in both directions, or domain bounds on both sides, which is the clean
negative (b) for this target.

### 10.1 CORRECTION to sec. 10: the limit member IS the chain tool's member (13:04-13:12 AEDT)

Sec. 10 said that the branch's sigma -> 0 limit is "a different member of the same family", with
epochs shifted by up to about 11.5 min. That was WRONG. I read the shifts at the sigma-0.02 start
as if they were the limit. The new driver stage `limit` (`data/1055_gc4-19/limit.json`) converges
the corrected branch going down from the start at the floors (sigma 0.01, 0.005 and 0.0025; 0.00125
does not converge within 40 + 60 iterations) and stores each point's node epochs.

Node-epoch offsets from the #973 gauntlet's flyby epochs (s), nodes C, C, G, G:

| sigma | C | C | G | G | V_inf (C, C, G, G) |
|---|---|---|---|---|---|
| 0.02 | -138 | -143 | -600 | -692 | 4.10987, 4.10822, 2.35619, 2.35643 |
| 0.01 | -84 | -92 | -375 | -431 | 4.10945, 4.10863, 2.35634, 2.35647 |
| 0.005 | -43 | -51 | -204 | -238 | 4.10923, 4.10882, 2.35641, 2.35649 |
| 0.0025 | -20 | -28 | -106 | -126 | 4.10912, 4.10892, 2.35646, 2.35649 |
| Richardson limit | -0.5 | -1.2 | -8.8 | -14.0 | 4.10901, 4.10901, 2.35650, 2.35650 |

The shifts fall about in proportion to sigma, and the limit equals the gauntlet member to within
14 s, which is the Richardson error. In #973 convention the limit's start dates are x_days =
(3.498195, 36.485538), against the gauntlet's (3.498357, 36.485552).

So the correction does not pick a different patched-conic member here. It moves the sigma-0.02
start along the near-null direction by an amount proportional to sigma: about 10 min at the
Ganymede nodes at sigma 0.02. In the limit that move vanishes. The full-mass orbit therefore
continues gc4-19's own #973 member. The caveat in sec. 10 reduces to "reached from the #1055
corrected seed", and no separate #973 gate run of a new member is needed.

Independent check (twobody-ext-opus, `docs/notes/2026-10-11-1055-gc4-19-973-gate.md`): from the
limit x_days (3.498195, 36.485538), the #973 `correct_dates` lands on the gc4-19 gauntlet zero within
1e-8 s (date residual 7.4e-5 km/s before correction), and it passes the full #973 gate (worst ratio
0.899, both screens, DOP853 5.3e-5 km).

## 12. GanCal#1@2013 on jup365 under amendment 7 (lead launch 13:04 AEDT)

- Fix (`fix.json`, 13:33; one call, 1,786 s): q = 1, L = 4.59e5. Curved-family scan:
  - + direction: no sign change within 3 L;
  - - direction: g = -4.03 (base), -6.66, -8.07, -7.02, -1.65 and +3.67 at path -2.9e5, -5.8e5, -8.8e5,
    -1.17e6 and -1.46e6. The secant root is at path -1.26e6 (about 2.7 L, inside the 3 L range),
    |g| = 1.9e-3.
  - No domain bound was reached.
- START (`start.json`, 13:34): CONVERGED. Two Newton iterations from the corrected seed, dr 1.56e-3
  km, dv 9.3e-7 km/s, ends within 1e-6 (amendment 3 rule). V_inf 3.0976 / 3.0761 / 3.2533 / 3.1626 km/s
  (raw seed at sigma 0.02: 3.0994 / 3.0958 / 3.2513 / 3.1644). This is the first start of GanCal#1
  under any formulation: #1044, #1046 apsis and #1046 crank all failed.
- Negative control: running at 14:15. The +0.5 L point is at |r| 5.6 with dv 5.4e-3 km/s after
  about 33 of its 40 iterations, which is failing to start as required. Result pending.

### 12.1 Negative control as run (`neg.json`, 14:37 AEDT) and AMENDMENT 8 (registered before it runs)

As run (the sec. 3 item 3 procedure, base null frame, offsets +-0.5 L = +-2.29e5 about the root):
- v0_plus: a family point was constructed, and Newton (40) FAILED to start: dr 1.29 km, dv 5.3e-3
  km/s (floors 1e-2 km noise rule, 1e-6 km/s). As required.
- v0_minus: the bordered corrector could not construct a family point, so no start was attempted.
  This is not a test.

How the pre-registration reads it: sec. 3 item 3 asks for BOTH points to fail. One point is a valid
failure; the other is missing. So the negative control is INCOMPLETE, neither satisfied nor
violated. Two weaknesses of the as-run procedure also matter:
- The offsets were taken along the BASE frame's null vector, but the root sits 2.7 L along a curved
  family (amendment 7). At the root, the base vector is not the family tangent. The "minus" point
  then lies off the family, which plausibly explains why it could not be constructed.
- The SVD sign of the base vector is not tied to the scan's direction of travel. So which of
  v0_plus / v0_minus lies beyond the root and which lies back towards the base is not recorded.

AMENDMENT 8, registered now at 14:45 AEDT, before the run:
- `--stage neg2`: the null pair is recomputed at the root (the local frame, as amendment 7 traces
  the family). It is oriented along the scan's direction of travel, from `fix_state.json`'s secant
  segment.
- On EACH side, beyond the root (+0.5 L, then +0.25 L if the first cannot be constructed) and back
  towards the base (-0.5 L, then -0.25 L), the first constructed point is started with Newton (40).
- SATISFIED if the first constructed point on each side fails to start. The original v0_plus failure
  stands as supporting evidence. If a side has no constructible point at 0.5 L or 0.25 L, the
  control is one-sided, and that is said plainly with the GanCal#1 verdict.
- VIOLATED if either side's point starts. Then the start is not shown to depend on the correction,
  and the GanCal#1 continuation result (launched by the lead at 14:38) is not used as #1055
  evidence.
- Expected runtime: about 2 x (corrector + up to 40 Newton iterations) at 0.5-1.5 min each, i.e.
  0.5-1.5 h, running alongside the lead's continuation.

### 12.2 AMENDMENT 9 (registered 14:55 AEDT, BEFORE the trigger; reactive to a foreseen trigger, and said so)

- What I saw (`sigma_n1.json` of the lead's 14:38 continuation, points at sigma 0.02 and 0.024): the
  corrected member's weak G2 flyby has its periapsis at 2,123 km at sigma 0.024. That is rp/SOI =
  0.388, and 88,000 km at full scale, against the chain tool's 28,327 km. rp/SOI grows about as
  sigma^0.6, so rungb's acceptance (rp/SOI > 0.5, from #968 rung (b), inherited through #1044
  criterion (ii)) will stop the run within a step or two, near sigma 0.03. It has not fired yet.
- The precedent: #1043 met the same issue with the EGGIE's 0.9-deg Io flyby, whose periapsis lies
  beyond the SOI even in the patched conic. It dropped the SOI bound and kept the floors and the
  unscheduled-pass check.
- AMENDMENT 9 (lead to accept or reject): for #1055's jup365 targets, the 0.5-SOI acceptance is
  dropped (`--soi-limit inf`). Unchanged: the floors at every point, IAS15 (#1044 iv), identity
  (amendment 2), and the inside-body scan at sigma = 1 (`--stage scan`: DOP853, rtol 1e-12,
  0.0005 d, all four moons; the same algorithm as #1025 amendment 2, whose positive control flagged
  the #1033 interior pass). rp/SOI per node is recorded as a flag. If the lead rejects it, GanCal#1
  fails #1044 (ii) at about sigma 0.03, and the full-mass question stays open.
- The relaunch resumes from the last accepted point of `sigma_n1.json` (rungb's resume), so nothing
  computed is lost.

### 12.3 Results (pending: neg2 launched 14:40; continuation launched 14:38)

To be filled against:
- neg2: SATISFIED if the first constructed point on each side fails to start (amendment 8).
- continuation: sigma = 1 at the floors (#1044 i); the floors at every point (ii, with amendment 9
  if accepted); identity going down (amendment 2: Cauchy limit, and a Richardson limit within
  0.05 km/s of the chain); IAS15 < 1e-2 km, < 1e-7 km/s (iv); the inside-body scan clean.

### 12.2.1 Amendment 9 accepted by the lead (option (a)), on these terms, recorded before the relaunch

1. Two readings, both kept here and in any row text:
   - UNDER THE ORIGINAL PRE-REGISTRATION: the stop on the 0.5-SOI acceptance near sigma 0.03 means
     GanCal#1 FAILS #1044 criterion (ii). The current run is left to reach that stop, so the
     result is on disk (`sigma_n1.json` `acceptance_fail`).
   - UNDER AMENDMENT 9: the continued run (0.5-SOI acceptance dropped, as #1043 did for weak
     flybys) is labelled so in every statement of its result.
2. Binding under amendment 9: the floors at every point, IAS15 (< 1e-2 km, < 1e-7 km/s), identity
   (a Cauchy limit going down and a Richardson limit within 0.05 km/s of the chain tool's V_inf), and
   the inside-body scan at sigma = 1. rp/SOI per node is recorded as a flag.
3. The identity gate decides one more question, and I will state the answer: is the corrected member,
   with G2's periapsis near 88,000 km at full scale, the chain tool's object (G2 at 28,327 km), or
   a DIFFERENT member of the degenerate family? V_inf alone does not answer it, so I will also
   compare the limit's node epochs and G2 periapsis direction with the chain tool's, as was done for
   gc4-19 (sec. 10.1).

### 12.2.2 Result UNDER THE ORIGINAL PRE-REGISTRATION: GanCal#1 FAILS #1044 criterion (ii)

From `data/1055_gancal1/sigma_n1.json` (the lead's 14:38 continuation, stopped as foreseen):
- Accepted points: sigma 0.02, 0.024, 0.0288, 0.03456, each at the floors or under the noise-floor
  rule.
- At sigma 0.041472, rungb's acceptance FAILED on the second Ganymede node (G2, t = 25.03 d from x0):
  rp 3,646 km, rp/SOI 0.535 > 0.5 (altitude 3,537 km, floor ok, V_inf 3.0558 km/s).
- Per point, as (V_inf km/s, rp km, rp/SOI):

| sigma | G1 | G2 | C | G3 |
|---|---|---|---|---|
| 0.02 | 3.0976, 57, 0.011 | 3.0761, 1,772, 0.348 | 3.2533, 1,062, 0.135 | 3.1626, 62, 0.012 |
| 0.024 | 3.0973, 69, 0.013 | 3.0723, 2,123, 0.388 | 3.2536, 1,215, 0.143 | 3.1623, 74, 0.014 |
| 0.0288 | 3.0969, 83, 0.014 | 3.0677, 2,543, 0.432 | 3.2538, 1,386, 0.152 | 3.1619, 89, 0.015 |
| 0.03456 | 3.0965, 99, 0.016 | 3.0623, 3,045, 0.480 | 3.2542, 1,576, 0.161 | 3.1615, 107, 0.017 |

- G2's rp / sigma stays near 88,000 km (88,600 to 88,100 km), about 3.1 times the chain tool's G2
  periapsis (28,327 km). G2's V_inf drifts down with sigma (3.0761 to 3.0623), while the other nodes
  move by less than 1.2e-3.
- Under the original pre-registration this is the GanCal#1 result: it fails (ii) at sigma 0.0415. It
  says nothing about full-mass existence. The continued run under amendment 9 (lead relaunch 14:42:
  cont --soi-limit inf, down, ias15, scan) is reported separately and labelled as such.

### 12.4 Continuation UNDER AMENDMENT 9 (lead relaunch 14:42-15:32 AEDT): classification and verdict

From `data/1055_gancal1/sigma_n1.json` and the runlog:
- 30 accepted points, sigma 0.02 to 0.737400, every node above its floor, no noise-floor-limited
  point after the start region.
- The natural step then collapsed: 0.737418 failed with the factor at 1.00001, "step below 1e-5 at
  0.737400; stop" (`stopped_at` 0.737400). No sigma = 1 point exists, so the ias15 stage had nothing
  to fly and the scan did not run (the chain used &&).
- Last converged nodes (sigma 0.7374):
  - G1: V_inf 3.0706, altitude 152 km (scaled floor 74 km);
  - G2: 2.5660, rp 74,841 km, rp/SOI 3.47;
  - C: 3.2908, rp 11,179 km;
  - G3: 3.1356, altitude 339 km.
  G2's rp is 2.36 times Ganymede's Hill radius (31,717 km): the "flyby" is now a distant encounter.

Is the stop at 0.7374 a fold?
- G2's rp against sigma over the last 10 points fits a square-root law, rp = rp_f - b sqrt(sigma_f -
  sigma), with sigma_f = 0.73777, rp_f = 75,592 km and residuals 13-148 km. A straight line leaves
  220-880 km. The derivative d rp / d sigma grows 2.8e5, 4.2e5, 6.5e5, 9.2e5, 1.3e6, 2.1e6, 4.4e6 km
  over the last 8 points, which is the divergence a turning point produces.
- No monitor-variable continuation and no two-solution check were run. rungb's `stage_sigma` has
  neither; #1025 has both.
- So under the #1025 outcome rules the classification is NUMERICAL STOP, with a FOLD SIGNATURE
  (inferred, not shown).

Identity (amendment 2), from the down stage (0.01 and 0.005 at the floors):
- max |dV_inf| against the chain tool is 0.0100 at 0.01 and 0.0051 at 0.005: it halves with sigma,
  a Cauchy limit. The Richardson limits are G1 3.0992, G2 3.0956, C 3.2515 and G3 3.1642, against
  the chain's 3.0994 / 3.0958 / 3.2513 / 3.1644 (within 2e-4). PASS (0.05 band).
- The node epochs converge to the chain's. Offsets at sigma 0.02, 0.01 and 0.005 are G2 -591, -293,
  -146 s and C +118, +72, +41 s, falling in proportion to sigma.
- BUT G2's periapsis does NOT converge to the chain's: rp / sigma = 88,617, 89,091 and 89,399 km
  against the chain's 28,327 km. Its limit turn is about 1.3 deg against the chain's 4.03 deg. (The
  Callisto node, 53,084 -> 61,742 -> 68,528 km, is still moving towards the chain's 79,393 km.)
  For comparison, gc4-19's limit converges to its chain in every respect: rp / sigma 9,315 / 9,310 /
  3,222 / 3,203 against 9,313 / 9,313 / 3,214 / 3,214, and turns 5.00 / 5.01 / 41.70 / 41.87 against
  5.01 / 5.01 / 41.77 / 41.77.

The answer to the lead's question (3): the corrected GanCal#1 branch is a DIFFERENT member of the
degenerate patched-conic family. It has the chain tool's flyby epochs and V_inf magnitudes, but
G2's periapsis is near 89,000 km (turn 1.3 deg) instead of 28,327 km (turn 4.03 deg), with the
difference absorbed by the crank of the full-revolution return. It is not the chain tool's
GanCal#1@2013 object. The identity gate as worded (V_inf only) passes. It does not by itself decide
same-or-different, and the periapsis comparison does.

VERDICT for #1055's GanCal#1 target, with both readings:
- Under the original pre-registration: FAILS #1044 criterion (ii) at sigma 0.0415 (G2 rp/SOI 0.535,
  sec. 12.2.2).
- Under amendment 9: the corrected seed STARTS. The branch continues at the floors to 74 % of the
  moon masses and stops there in a NUMERICAL STOP with a fold signature. The branch is a different
  family member from the chain tool's object (G2's periapsis three times farther). No full-mass
  point, so IAS15 and the scan do not apply.
- So: the jup365 lane is NOT validated on a full-revolution chain at full mass by GanCal#1. The
  method change itself is shown to work on jup365: a start where no formulation had one, identity
  in V_inf and epochs, and a negative control pending (neg2).
- Negative control (amendment 8): neg2 still running; the verdict above is conditional on it.

### 12.5 Single most informative follow-up (PRE-REGISTERED here, NOT run)

The fold test at the stop, by the #1004/#1025 rules:
- Monitor-variable continuation from the last two points (sigma becomes an unknown; the
  fastest-changing scaled state component is stepped), up to 30 points, using `monitor_step`
  generalised to the Target face. Then, if sigma turns back below 0.9 sigma_f, two solutions at
  fixed sigma 0.675 and 0.720 (0.915 and 0.976 x 0.7374), from the lower and upper branches.
- FOLD supported if both converge at the floors and some node differs by > 0.02 km/s in V_inf and
  > 100 km in r_p / sigma (the #1025 rule). FOLD then means "no full-mass orbit on this branch",
  as for gc-2. If the monitor carries sigma past 0.7374 upward instead: NOT a fold; natural
  continuation resumes towards sigma = 1 (which is then the jup365 verdict).
- Runtime: about 15-40 s per evaluation at the current load, 30 monitor points of about 5-10
  evaluations each, plus two solutions: 1.5-3.5 h (lead launch). Code: about 1 h to port
  `monitor_step` and `two_solutions` to the jup365 Target, with a fast test on the ideal gc-2 fold
  (it must reproduce sigma_f 0.1639 and the two-solution numbers of #1025 sec. 7).

### 12.5.1 Fold-test code and its control (registered before the GanCal#1 fold run)

- Code: `src/cyclerfinder/nbody/fold_tools.py` (`monitor_step`, `sigma_derivative`,
  `fold_supported`), generalised from #1025's `monitor_step` and the #1004 two-solution rule. Two
  fast tests (`tests/nbody/test_1055_fold_tools.py`):
  - an analytic fold at sigma_f = 0.1639 with branches x = +-2e5 sqrt(sigma_f - sigma) km. The
    monitor passes the turn, never exceeds sigma_f, and lands exactly on the other branch;
  - the two-solution rule on the #1004 sec. 8.3 numbers.
  Driver stage `--stage fold` (resumable per point, `fold.json`).
- CONTROL, ideal gc-2 (`data/1055_gc-2/fold.json`, 15:36-15:41 AEDT), from #1025's natural points:
  - monitor sigmas 0.163948, 0.163909, 0.163617, 0.160595 and 0.121001, identical to #1025 sec. 7;
  - two solutions at 0.15001: max dV_inf 0.0475 km/s, d r_p / sigma 657 km;
  - at 0.16001: 0.0252 km/s and 349 km (#1004: 0.047 / 657 and 0.025 / 350);
  - status FOLDS. PASS.
- AMENDMENT 11 (jup365 only): the monitor points and the two fixed-sigma solutions accept the
  noise-floor rule (dr < 1e-2 km, dv < 1e-6 km/s, ends < 1e-6), as every intermediate jup365 point
  of #1044 / #1055 does. The ideal control kept the floors.

### 12.5.2 What the different-member finding means for the method (AMENDMENT 10, registered now)

GanCal#1 (sec. 12.4) shows that a corrected seed can lie on a DIFFERENT member of the degenerate
patched-conic family. The node epochs and V_inf magnitudes are the chain's, but a flyby's periapsis
and turn are not (G2: 89,000 km and 1.3 deg against 28,327 km and 4.03 deg). An identity gate on
V_inf alone (amendment 2) cannot see this, because the family is degenerate in V_inf by
construction.

AMENDMENT 10, part of the method's acceptance from now on, for every #1055-corrected seed used on any
row (the gc rows, gc-1 e2, and any later use):
- MEMBER IDENTITY. At the smallest converged sigma going down, and in the Richardson limit, every
  flyby node must match the chain tool's patched-conic member within:
  - V_inf: 0.01 km/s;
  - r_p / sigma: 5 %;
  - turn: 0.2 deg or 5 %, whichever is larger;
  - node epoch: 60 s (the limit; the sigma-0.02 offsets are O(sigma) and are reported).
- PASS: the closure is the chain tool's object. FAIL: the closure belongs to a different member of
  the same family. It is then reported as such and may not stand for the chain tool's object.
- Applied to the two targets so far:
  - gc4-19: PASS. At sigma 0.0025, r_p / sigma is 9,315 / 9,310 / 3,222 / 3,203 km against 9,313 /
    9,313 / 3,214 / 3,214 (within 0.4 %); the turns are within 0.1 deg; the Richardson limit of the
    epochs is within 14 s; V_inf within 1e-4.
  - GanCal#1: FAIL at G2 (r_p / sigma 89,399 km against 28,327 km; turn about 1.3 deg against
    4.03 deg). The Callisto node (68,528 km at sigma 0.005 against 79,393 km) is still moving
    towards the chain's value. As sec. 12.4 already said, the GanCal#1 closure attempt concerns a
    different member.
- Why this belongs to acceptance, not to analysis: the correction selects the member by the
  first-order solvability condition. When that condition's root differs from the chain tool's
  member, a "successful" closure validates a different orbit. Only a member-identity check makes
  "the row's object exists in continuous gravity" a true statement.

### 12.5.3 Prepared (not run): resume toward sigma = 1 if the fold test reads "NOT A FOLD"

`--stage resume` appends the monitor points beyond the natural stop to `sigma_n1.json` and resumes
rungb's natural continuation (factor 1.01, regrowing), under amendment 9 (`--soi-limit inf`).
The chain then runs down, ias15 and scan, with the amendment-10 member-identity check applied in
the analysis:

    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage resume --soi-limit inf && \
    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage down && \
    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage ias15 && \
    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage scan

## 12.6 Fold test on GanCal#1 (lead launch 15:43; ended 15:46 AEDT): result and reading

From `data/1055_gancal1/fold.json` and `live/lead_fold_1.out`:
- The monitor (state component 11, the stepped one: G2's moon-relative velocity) went from sigma
  0.737364 and 0.737400 (the last two natural points) through 0.737395, 0.737307, 0.736748 and
  0.730763 to 0.607634. sigma DECREASED after the turn, by 0.13 in total: the branch turns back.
  Status "two" was reached at 0.6076 < 0.9 sigma_stop.
- Two fixed-sigma solves, all four converged (amendment 11 noise rule), lower start = nearest
  natural point, upper start = nearest monitor point after the turn:
  - sigma 0.67472: two DIFFERENT orbits. G2 rp 62,924 km (lower) against 81,027 km (upper),
    i.e. 26,830 km in rp / sigma, with V_inf 2.5840 against 2.5748. The Callisto node differs by
    0.0203 km/s and 477 km rp / sigma, which meets the #1025 rule (> 0.02 and > 100 at one node):
    SUPPORTED.
  - sigma 0.71970: two DIFFERENT orbits again. G2 rp 69,423 against 78,710 km (12,904 km in
    rp / sigma), but the largest V_inf difference at any node is 0.0108 km/s (Callisto), below
    the rule's 0.02: NOT SUPPORTED.
- So the status "NUMERICAL STOP: turn passed, two solutions not supported" does not mean the two
  solves collapsed onto one point or failed. They are two distinct converged orbits at each sigma.
  At 0.976 sigma_f their V_inf separation is smaller than the fixed 0.02 km/s threshold. That
  threshold was set on gc-2, where the branches separate mainly in V_inf. GanCal#1's branches
  separate mainly in G2's periapsis. Near a fold the separation shrinks as sqrt(sigma_f - sigma).
  From 0.67472 to 0.71970 the predicted ratio is about 0.53. Observed: 0.48 in G2 rp / sigma and
  0.53 in the Callisto V_inf difference.

VERDICT for GanCal#1 under the pre-registration (unchanged by this reading):
- Original pre-registration: FAILS #1044 (ii) at sigma 0.0415.
- Amendment 9: NUMERICAL STOP at sigma 0.7374. The fold is UNDECIDED by the pre-registered rule
  (supported at 0.915 sigma_f, not at 0.976 sigma_f). No full-mass orbit is shown.
- Amendment 10: the branch is a DIFFERENT member of the degenerate family (G2 periapsis), so even a
  decided fold or a full-mass closure would concern that member, not the chain tool's
  GanCal#1@2013 object.
- Evidence beyond the rule (reported, not a verdict): the monitor turns back by 0.13 in sigma, two
  distinct converged orbits exist at both test sigmas, and their separation follows the square-root
  scaling of a fold. All of this points to a fold at about sigma 0.7377 on this member's branch.

Would one further step decide it cheaply? PRE-REGISTERED here as AMENDMENT 12, post-hoc: written
after the 15:46 result, and labelled so in any statement of its outcome. It is not run without the
lead's launch.
- `--stage foldsqrt`: the same two fixed-sigma solves at sigma 0.65, 0.69 and 0.72. FOLD CONFIRMED
  if all six converge, the two solutions are distinct at every sigma (some node > 1,000 km in
  rp / sigma), and the square of the most separated node's rp / sigma separation is linear in sigma
  with R^2 > 0.95, vanishing at sigma_f within 0.01 of 0.73777 (the natural-stop fit, sec. 12.4).
  The tolerance is wide because the two existing points alone extrapolate to 0.7332: higher-order
  terms matter about 0.06 below the fold. Runtime: 6 solves at about 20-40 s each plus the
  integrations, i.e. 5-10 minutes.
- What it would add: a decided fold for this member's branch, so "no full-mass orbit on this
  branch", as for gc-2. It would NOT test the chain tool's object (amendment 10), and so it does not
  change whether the jup365 lane is validated on a full-revolution chain at full mass. My
  recommendation: run it only if a decided GanCal#1 statement is wanted for the record. The
  method's open question is now amendment 10's member selection, not this fold.

Command (lead launch):

    uv run python scripts/run_1055_seedfix.py --target gancal1 --stage foldsqrt

## 12.7 Member selection on GanCal#1: why the root is not the chain tool's member (analysis only, data on disk)

Data: `data/1055_gancal1/fix_state.json` (the amendment-7 scan: path, g, complement), `fix.json`
(the root) and the seed. Units: path in scaled arclength, g in weighted residual units, L = 4.59e5.

1. WHERE THE CHAIN TOOL'S MEMBER SITS. The scan's base (path 0) is the chain tool's seed with its
   complement rows converged at sigma 0.02. The base has G2 rp / sigma = 21,646 km (the seed
   28,327 km). Converging the complement moves the non-degenerate coordinates too, at O(sigma).
   At the base, g = -4.03. That is the same size as g anywhere on the scan (|g| 1.6 to 8.1 away
   from the root): the chain tool's member is a clear NON-ROOT of the solvability function.
2. THE SCAN, BOTH SIDES:
   - + side: g = -4.03, -3.99, -3.91, -3.76, -3.47, -2.95, -2.19 and -1.65 at path 0 ... 4.36e5
     (0.95 L), rising and roughly linear in path. The scan then STOPPED AT THE DOMAIN BOUND
     (amendment 7c): a flyby node dropped below sigma (R + floor). On this side the family drives
     G2's periapsis INTO Ganymede (sec. 11), so the bound is G2 reaching its floor. Extrapolating
     the slope, g would cross zero near path 9e5 (about 2 L), past the bound. That root, if it
     exists, is an impact-side member with G2 below its floor: outside the physical domain.
   - - side: g = -4.07, -4.15, -4.32, -4.66, -5.36, -6.66, -8.07, -7.02, -1.65 and +3.67 at path
     -4.6e3 ... -1.46e6. Exactly ONE sign change, between -1.17e6 and -1.46e6, refined to -1.26e6
     (2.7 L). That root's member has G2 rp / sigma = 88,617 km and Callisto 53,086 km (base:
     21,646 and 42,964).
   - Resolution caveat: the samples double in spacing up to 2.9e5. A pair of roots between two
     samples of one sign cannot be excluded. g is smooth and of one sign over each sampled
     stretch, so nothing points to one.
3. SO, within the scanned domain there is ONE physical solvability root, at -2.7 L. The chain
   tool's member (path 0) is not a root. The other candidate root, on the + side, lies past the
   G2 floor.

What it means:
- The chain tool chooses its member of the degenerate family through the ephemeris's weak symmetry
  breaking (the fixed-leg crank solved for arrival at the moon, #942). Continuous gravity selects
  members through the first-order solvability condition, which depends on the moons' gravity along
  the legs. Nothing ties the two choices together. On GanCal#1@2013 they differ by 2.7 L, the
  difference between G2 at about 28,000 km and at about 89,000 km.
- This EXPLAINS the #1044 / #1046 negatives: every one of those formulations started from the chain
  tool's member, which is not the limit of any continuous-gravity orbit. A solver can only crawl or
  stall there. The out-of-plane mismatch they all met (#1046) is g, the residual's projection on
  the null direction, which no formulation can remove at a non-root.
- The chain tool's GanCal#1@2013 object, as such, is not the sigma -> 0 limit of a continuous-gravity
  orbit within the scanned domain. Its family has one physical continuous-gravity member nearby,
  which is not it.

For the other gc rows:
- The 16 rows reached from their RAW seeds (the 15 #1025 members, plus gc-1) did not use the
  correction: their raw seed started and continued, with identity to the chain. Nothing here
  changes them. In the coplanar ideal model, g vanishes by symmetry for gc-1 (sec. 1), and the
  raw-seed starts show the chain members are roots (or so near one that Newton absorbs the
  difference).
- gc4-19 (corrected) passed amendment 10: there the root and the chain's member coincide in the
  limit (sec. 10.1).
- For any FUTURE row that needs the correction (the real-ephemeris rung of the gc rows, gc-1 e2):
  amendment 10 decides. A FAIL is a finding about the row: its patched-conic object does not
  continue into continuous gravity, and a nearby family member does instead. It then cannot be
  reported as "the row's object exists". GanCal#1 shows the real-ephemeris chain tool can pick a
  non-root, so expect amendment-10 failures on jup365 to be possible, not rare.

Cheapest informative next step (PRE-REGISTERED here, not run): gc-1 e2 fix only. Run
`--stage fix` (amendment 7) on gc-1 e2 and record where its chain-tool member sits relative to the
roots: g at the base, the sign changes within 3 L on each side, the domain bounds, and the root
member's periapses against the chain's. No start, neg or continuation unless the root passes
amendment 10 against the chain tool's member. Runtime: the fix stage at about 5-7 times
GanCal#1's cost per evaluation, so 3-6 h at the current load (lead launch). This answers the
question that matters for the gc rows' real-ephemeris standing: does the gc-1 chain tool's member
coincide with a continuous-gravity root on jup365?

## 13. gc-1 e2 on jup365: the fix stage (lead launch 15:50 AEDT)

For the record: gc-1's catalogue row (ganymede-callisto-two-working-body-cycler-gc1-2026) carries,
in its data_gaps, "Real-ephemeris continuous-gravity standing untested: the open-chain n-body
formulation does not converge across a full-revolution return (#1044) ... obstruction: out-of-plane
mismatch across the full-revolution leg (#1046)". This fix stage is the FIRST DIRECT TEST of that
standing. Sec. 12.7 explains the #1044/#1046 obstruction as a start from a non-root. The test
records:
- g at the base (gc-1's chain-tool member);
- the sign changes and domain bounds within 3 L on each side;
- any root member's periapses, turns and epochs against the chain's (amendment 10).
Start, neg and continuation follow ONLY if the root passes amendment 10, on my call. A FAIL means
gc-1's real-ephemeris chain-tool member does not continue into continuous gravity on jup365, and the
row's data_gaps wording would then change (a lead / owner decision).

### 12.8 foldsqrt (amendment 12, POST-HOC; lead launch 15:52 AEDT on fa6c247f): FOLD CONFIRMED

From `data/1055_gancal1/foldsqrt.json` (all six fixed-sigma solves converged, amendment 11 rule):

| sigma | d (r_p / sigma), nodes G1, G2, C, G3 (km) | d V_inf, same nodes (km/s) |
|---|---|---|
| 0.65 | 47, 33,494, 561, 37 | 0.0096, 0.0171, 0.0240, 0.0096 |
| 0.69 | 35, 22,550, 416, 27 | 0.0068, 0.0052, 0.0177, 0.0068 |
| 0.72 | 21, 12,785, 254, 16 | 0.0040, 0.0001, 0.0107, 0.0040 |

- The two solutions are distinct at every sigma (G2 differs by 12,785-33,494 km in r_p / sigma).
- The squared G2 separation is linear in sigma: R^2 = 0.9939. It vanishes at sigma_f = 0.73006,
  0.0077 below the natural-stop square-root fit (0.73777), inside the registered 0.01 tolerance.
- Under amendment 12: FOLD CONFIRMED (labelled post-hoc: the criterion was written after the
  15:46 fold-test result). The fold location is bracketed by 0.730 (separation fit) and 0.7377
  (natural-stop fit). The natural continuation's last point is 0.7374.

FINAL GanCal#1@2013 verdict (#1055 jup365 target), all readings:
- Under the ORIGINAL PRE-REGISTRATION, GanCal#1 FAILS #1044 criterion (ii) at sigma 0.0415 (G2
  rp/SOI 0.535).
- Under AMENDMENT 9 (0.5-SOI acceptance dropped, as #1043 did for weak flybys):
  - the #1055-corrected seed STARTS, the first start under any formulation;
  - it tracks at the floors to sigma 0.7374;
  - the branch then FOLDS at sigma 0.730-0.738: the pre-registered two-solution rule leaves it
    undecided, and the post-hoc amendment 12 confirms it;
  - so there is NO full-mass orbit on the corrected branch.
- Under AMENDMENT 10, that branch is a DIFFERENT member of the degenerate family (G2 periapsis
  near 89,000 km against the chain tool's 28,327 km). The chain tool's own GanCal#1@2013 object is
  therefore UNTESTED at full mass by this run. Per sec. 12.7 it is a NON-ROOT of the first-order
  solvability condition, so it is not the small-mass limit of any continuous-gravity orbit within
  the scanned domain.
- Negative control (amendment 8, neg2): the "beyond" point fails as required; the "back" point is
  still running. The start's dependence on the correction is conditional on it (sec. 12.3).
- The jup365 lane is NOT validated on a full-revolution chain at full mass by GanCal#1.

### 13.1 gc-1 e2 fix result (15:50-15:57 AEDT): q = 0, a weak direction that does NOT scale with sigma

From `data/1055_gc1e2/fix.json`:
- Column-scaled singular values (smallest four):
  - at sigma 0.02: 4.55e-4, 2.47e-4, 3.81e-5, 2.75e-7;
  - at sigma 0.005: 4.54e-4, 2.47e-4, 3.80e-5, 2.66e-7.
- The detection rule (sec. 2 item 1: a trailing singular value is degenerate if it drops by more
  than 2.5 between sigma0 and sigma0 / 4, proportional-to-sigma being 4) sees a ratio of 1.04, so
  q = 0, the seed is unchanged, and amendment 10 is trivially passed (the seed is the chain's member).

Against my expectation. I expected gc-1 e2 to show the sigma-proportional direction that the
ideal-model gc-1 has (3.2e-6, 9.9e-7 and 3.0e-7 at 0.02, 0.005 and 0.00125, sec. 1). It does not.
Instead it has a VERY weak direction, 2.75e-7 (140 times below the next singular value), whose size
is INDEPENDENT of sigma. In the ideal coplanar model the 1:1 Callisto return's family is exactly
degenerate at sigma = 0. On jup365, Callisto's real eccentricity and inclination break that
degeneracy by a fixed amount (INFERRED). The family is then nearly degenerate, with a strength set
by the ephemeris, not by sigma. This is below the detection threshold by construction: the rule
looks for a singular value that FALLS with sigma, and this one does not fall at all. A run at sigma
0.0025 or 0.00125 would show the same 2.7e-7, by the evidence of 0.02 against 0.005.

What it does to the #1046 / sec. 12.7 explanation for gc-1:
- The sec. 12.7 mechanism (the chain member is a non-root of an exactly degenerate family) does NOT
  apply to gc-1 e2 as stated: there is no exactly degenerate family, so there is no solvability
  root to miss.
- The mechanism that fits is its near-degenerate version. The ephemeris selects the chain's member
  with a stiffness of about 2.7e-7. Continuous gravity forces the residual along that direction by
  an amount proportional to sigma (u . r ~ c sigma). The Newton move needed along it is then about
  c sigma / 2.7e-7. That is large at sigma 0.02, but it shrinks IN PROPORTION TO SIGMA, unlike
  GanCal#1's O(1) move. So for gc-1 e2 a start at smaller sigma should help, which is the opposite
  of GanCal#1.
- So, for gc-1, the #1044/#1046 start failures are explained (INFERRED) by a near-degenerate,
  ephemeris-stiffened direction with an O(sigma) forcing, not by a non-root. The out-of-plane
  character of the stuck residual (#1046) fits: Callisto's real inclination is what breaks the
  out-of-plane degeneracy.

PRE-REGISTERED reading of the start stage (lead launch 15:57; written before it ends):
- The start is plain rungb Newton from the UNCHANGED seed at sigma 0.02 (q = 0). Up to 40
  iterations, with the noise-floor rule.
- CONVERGED: the #1044 non-start was specific to its solver path (LM, then a stalled polish) or to
  #1044's iteration budget, not to the problem. The continuation and the amendment-10 identity (to
  the chain, trivial here) follow; I will give the commands.
- STALL (no convergence in 40 iterations, residual on the weak direction): the same obstruction as
  #1044, as the near-degenerate mechanism above predicts. The pre-registered next step is then
  `--stage weakscale` (not yet built; about 30 min of code): the residual's projection on the
  smallest singular direction and the implied Newton move at sigma 0.02, 0.005 and 0.00125. If
  the move falls in proportion to sigma, the start is moved to the sigma at which the move is
  below 10 % of the node offsets, and the continuation goes up from there. Runtime: 3 Jacobians
  at 35-110 s each, so under 10 minutes.
