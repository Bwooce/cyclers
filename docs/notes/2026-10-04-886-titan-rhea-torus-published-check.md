# #886 part 1 -- the #882 invariant-circle corrector against Kumar's Saturn-Titan-Rhea tori

**Date:** 2026-10-04. **Task:** `#886`, first part (tori only; no manifolds or connections).
**Script:** `scripts/screen_886_titan_rhea_torus_check.py` (staged). **Tests:**
`tests/search/test_titan_rhea_torus_886.py`. **Outputs:** `data/found/886_titan_rhea_torus/`.

**Published result tested.** B. Kumar, "Analysis of Unstable Resonant Orbits for Saturn Tour
Design: Between Titan and Rhea", IAC-25-C1.9.6, 2025, filed in the private paper corpus as
`kumar-2025-analysis-unstable-resonant-orbits-saturn-tour-design-titan-rhea-IAC-25-C1.9.6-doi-10.52202-083087-0076.pdf`
(digest `docs/notes/2026-10-03-digest-kumar-2025-iac-titan-rhea-unstable-resonant-orbits.md`).
Section 6: the Titan 3:2 low/mid-e unstable family, continued in Rhea's mass to
mu3 = 4.05746e-6 as invariant circles of the stroboscopic map, persists "from a wide range of
Jacobi constants", and "many mid-e orbits having periods near the aforementioned secondary
resonances [4/21, 5/26, 6/31, 7/36, 8/41, 9/47] failed to continue as tori"; the resonance
regions "do not overlap". The paper prints no initial conditions, no per-orbit outcome and no
tolerances.

**Machinery under test (read-only).** `cyclerfinder.search.ccr4bp_strob_connection`
(`seed_circle_from_periodic_orbit`, `correct_invariant_circle`, `circle_residual`,
`hyperbolic_bundles`, `symmetric_periodic_orbit`) on `cyclerfinder.core.ccr4bp.CCR4BPSystem`.
The continuation driver (adaptive steps, secant predictor, retry from the last converged
circle) is a wrapper in the script; the module's own `continue_circle_in_mass` uses fixed
steps and stops at the first failure, which cannot separate a corrector failure from a torus
failure.

## 1. Model comparison (paper against `core/ccr4bp.py`)

Checked against the paper's Section 2.2, Eq. [4]-[5] (p. 2-3), read in full.

| Item | Paper | `core/ccr4bp.py` | Effect |
|---|---|---|---|
| Frame, units | m1-m2 synodic frame, barycentre origin, `G(m1+m2) = r12 = Omega2 = 1` | same | none |
| Base problem | PCRTBP, Saturn at `(-mu, 0)`, Titan at `(1-mu, 0)` | same (`cr3bp_eom`) | none |
| Variables | position-momentum `(x, y, px, py)` | position-velocity `(x, y, vx, vy)`, `px = vx - y`, `py = vy + x` | none: a change of variables; the stroboscopic circles are compared in velocity form throughout |
| mu3 | `m3/(m1+m2)` (typeset fraction garbled in extraction; denominator `m1+m2` for both mu and mu3) | `mu_gan = m3/(m1+m2)` | same definition |
| Perturber rate | `theta3 = (Omega3 - 1) t + theta3,0`, `Omega3 = sqrt(G(m1+m3)/r13^3)` | `theta = theta_gan0 + omega_gan t`; `two_body_synodic_rate` is the same Kepler law | same sign convention: Rhea is inside Titan's orbit, so `omega_gan = +2.5297 > 0` (Rhea advances in Titan's frame) |
| **Perturber circle centre** | **Saturn**: `(-mu + r13 cos theta3, r13 sin theta3)`; indirect term `-mu3 (cos, sin)(theta3) / r13^2`, i.e. Saturn's acceleration toward Rhea | **barycentre**: `(a cos theta, a sin theta)`; indirect term `-mu_gan r_gan / a^3`, the barycentre's | **the one real difference**: the perturber is displaced by `mu = 2.4e-4` in length. The force difference is of order `mu * mu3`; quantified on converged circles in section 5 by the invariance residual under a Saturn-centred right-hand side written in the script (Kumar Eq. [4] in velocity form) |
| Mass continuation | `r13` varied with mu3 so that `Omega3`, hence the rotation number, stays fixed | `omega_gan` and `a_gan` fixed while `mu_gan` varies | rotation number fixed in both; the radius difference is `O(mu3)` relative (about 1e-6 in `a_gan`), negligible |
| Stroboscopic map | time-`2 pi/|Omega3 - 1|` map at `theta3,f = 0` | `strob_iterates`, period `2 pi/|omega_gan|`, `theta_gan0 = 0`, `t0 = 0` | same section (phase 0 measured about Saturn in the paper and about the barycentre here; the same up to the centring above) |

The project model is the paper's model up to the perturber's centre. That changes no
qualitative behaviour (the secondary-resonance structure depends on Rhea's mass and frequency,
both identical), so a persistence/failure comparison is meaningful; it does prevent
digit-level comparison of four-body quantities, of which the paper prints none anyway.

## 2. Constants

All from `data/found/886_titan_rhea_torus/constants.json`.

- Used exactly as printed: `mu = 2.36639e-4`, `mu3 = 4.05746e-6`, `Tp = 2.48376`. The model's
  `omega_gan = 2 pi / Tp = 2.529707100` (so the model's Tp is the printed value exactly) and
  `a_gan = ((1 - mu + mu3)/(1 + omega)^2)^(1/3) = 0.431327565` (Kepler, as the paper defines
  `Omega3`).
- The printed "about 0.4315" for Rhea's radius does NOT reproduce the printed Tp: it gives
  `Tp = 2.485839` (0.08 per cent off). The radius consistent with Tp is 0.431328. So 0.4315 is
  loose rounding and Tp is the defining constant; this script lets Tp define the rate.
- Project ephemeris values (JPL SSD registry in `core/satellites.py`; Saturn SYSTEM GM
  3.7931207e7 km^3/s^2, Titan 8978.14, Rhea 153.94; semi-major axes 1221870 and 527070 km):
  `mu = 2.366953e-4` (+2.4e-4 relative to the paper), `mu3 = 4.058400e-6` (+2.3e-4),
  `a_Rhea/a_Titan = 0.431363` (+8.3e-5), which gives `Tp = 2.484192` (+1.7e-4). The paper's
  Saturn-Rhea two-body value 4.05841e-6 (its p. 2) matches the registry's Rhea/system ratio, and
  its 4.05746e-6 is that value times `(1 - mu)`, consistent with the `m3/(m1+m2)` definition.
  The paper's mu for Titan differs from the registry by 2.4e-4 relative; the source of its GM
  values is not stated.
- Farey check (independent of any orbit): the reduced fractions with `q < 50` strictly inside
  the printed range (0.1898475, 0.1952735) are exactly the paper's six, 4/21, 9/47, 5/26, 6/31,
  7/36, 8/41. With `q <= 100` there are more (listed in `constants.json`); they are used only to
  name any failure found away from the six.

## 3. The three-body family

Computed with the module's `symmetric_periodic_orbit` (x-axis-symmetric orbits, perpendicular
crossing at fixed `x`), continued in the crossing point `x` on the negative x-axis (the
periapse; Fig. 2 plots exactly this "x-intercept at symmetric point" against C). 119 members,
`data/found/886_titan_rhea_torus/family.json`. Each row has period `T`, Jacobi constant `C`,
`Tp/T`, the monodromy multiplier, and the distances to Titan, to Saturn and to Rhea's orbit
circle.

Against what the paper prints:

| Quantity | Paper | This family | Verdict |
|---|---|---|---|
| Upper end of Tp/T | 0.1952735 | **0.195273282** (interior minimum of T = 12.719405191 at x = -0.6822600, C = 3.0467074) | agrees to 2.2e-7, inside the rounding of the printed Tp: the Tp that would give 0.1952735 exactly is 2.4837643, which rounds to the printed 2.48376. Independent of where the paper cut its family, because it is an interior extremum |
| Lower end of Tp/T | 0.1898475 | reached at x = -0.49063, C = 2.950663, 1671 km above Titan's surface (minimum distance 0.003475 Titan units) | Fig. 2's family ends at C about 2.951 by eye: agrees. The family continues to Titan's surface at x = -0.4695, C = 2.9345, Tp/T = 0.18805; the paper's end is a stopping choice (perhaps an altitude floor), not a feature |
| Near-circular end | Fig. 2 ends at x about -0.745, C about 3.0555, with a vertical tangent | stability boundary (multiplier 1) at x = -0.7427208, C = 3.0547515, which is also the maximum of C (the fold) | agrees by eye; beyond it the orbits are stable and not part of the unstable family |
| Stability | "unstable" family | real multiplier from 1.0 (boundary) to 483 (paper's end), 930 at Titan's surface | agrees |
| Fig. 11 shape (T against C) | falls from about 13.1 near C = 2.95, local minimum near 2.99, local maximum near 3.01, minimum near 3.05, steep rise at the C fold | minima T = 12.797621 (C = 2.99229) and 12.719405 (C = 3.04671), maximum 12.841712 (C = 3.01280), steep rise to T = 12.853 at the fold | agrees by eye |
| Fig. 1/12 shape | three-lobed rotating-frame orbit, apoapse loops near Titan's orbit, periapse inside Titan's orbit and outside Rhea's | periapse radius 0.742 (boundary) to 0.490 (paper's end), never inside Rhea's orbit (smallest gap to Rhea's circle 0.058 at the paper's end) | agrees with "none ... intersect that of Rhea" |

**Monotone segments and crossings.** The three extrema of T split the family into four segments
on which Tp/T is monotone. The six ratios are crossed 11 times:

| Segment | x range | Tp/T range | crossings (x) |
|---|---|---|---|
| 0 (near-circular) | -0.74272 to -0.68226 | 0.193242 to 0.195273 | 6/31 (-0.740092), 7/36 (-0.728129), 8/41 (-0.705120) |
| 1 | -0.68226 to -0.58936 | 0.195273 to 0.193413 | 8/41 (-0.657045), 7/36 (-0.624974), 6/31 (-0.599104) |
| 2 | -0.58936 to -0.55190 | 0.193413 to 0.194080 | 6/31 (-0.578996) |
| 3 (high energy) | -0.55190 to -0.46951 | 0.194080 to 0.188053 | 6/31 (-0.532048), 5/26 (-0.515756), 9/47 (-0.507294), 4/21 (-0.497136) |

4/21, 9/47 and 5/26 are crossed only in segment 3, the high-energy part where the paper puts its
secondary-resonant orbits ("plotted in the higher-energy part of the family"). The three T
extrema are twistless points (`dTp/T` along the family is zero): Tp/T = 0.195273 (1.5e-4 above
8/41), 0.193413 (1.35e-4 below 6/31) and 0.194080.

## 4. Pre-registration (written and committed before the four-body stage was run)

Written after the three-body family and the mu3 = 0 seed resolution survey, BEFORE any
mu3 > 0 computation, and committed on its own.

**4.1 A resolution limit found before the four-body stage.** The invariant circle of the
stroboscopic map is parameterised by the dynamics (rigid rotation), so at mu3 = 0 it is the
periodic orbit sampled UNIFORMLY IN TIME, and its Fourier content is the orbit's content in time.
In segment 3 the orbit passes Titan at 0.003 to 0.03 Titan units, a fraction of about 1e-3 of the
period. Seed residuals (the exact periodic orbit sampled and checked against the map, i.e. pure
interpolation error) measured with the module: x = -0.55: 701 nodes give 9e-13; x = -0.52:
701 nodes give 3e-5; x = -0.50 and -0.49: 1e-2 at 701 nodes. One Gauss-Newton iteration of the
module's corrector (dense least squares, `O((4N)^3)`) took 18 s at N = 1001 and 26 s at N = 1201
(4 BLAS threads), with seed residuals 5e-9 and 3.5e-10 at the segment-3 6/31 crossing
(x = -0.532). The paper's method (Kumar, Anderson & de la Llave, CMDA 2022, its ref. [20]) is a
quasi-Newton solve in `O(N log N)` time per step, so the paper can use node counts this
corrector cannot.

**4.2 Two regions.**
- **Region A:** members whose mu3 = 0 seed passes the gate below at some
  N in {151, 201, 251, 321, 401, 501, 601, 701}. By the survey this is segments 0 to 2 and the
  start of segment 3 (x up to about -0.55, C >= about 2.99). It contains 7 of the 11 crossings:
  6/31, 7/36 and 8/41 only, all 0.13 to 0.31 Titan units from Rhea's orbit.
- **Region B:** members that fail the gate at 701 nodes. They are reported as NOT TESTABLE BY
  THIS CORRECTOR, with an estimate of the node count they would need (log-linear fit of seed
  residual against N). They are never counted as torus failures. Region B contains every
  crossing of 4/21, 9/47 and 5/26, and the segment-3 crossing of 6/31.

**Consequence, on record before any result:** Region A can test the paper's PERSISTENCE claim
and the absence of failures away from the ratios, but it cannot test the FAILURE claim at the
ratios the paper locates in the high-energy part. The best verdict this run can reach is
**"partly reproduced"**. A full test needs either a structured (Fourier-diagonal) linear solve
or a long run at N of order 1300 and more for the segment-3 6/31 crossing only (4/21, 9/47 and
5/26 need more than that).

**4.3 Per-member procedure.**
1. Seed gate (positive control of the representation): the mu3 = 0 circle sampled from the
   periodic orbit must have node invariance residual <= 1e-10 at the chosen N (smallest level
   that passes).
2. Fixed rotation number `rho = 2 pi Tp/T mod 2 pi`, perturber phase 0 at t = 0 (the paper's
   `theta3,f = 0`). Natural-parameter continuation in mu3 from 0 to 4.05746e-6 with a secant
   predictor: first step 0.05 of the target, growth x1.5 per success up to 0.25, halving on a
   failed step and retry from the last converged circle, floor 1/1024 of the target. A step
   succeeds if the module's `correct_invariant_circle` reaches node residual <= 1e-10 (the
   module's tolerance; the plan's 1e-9 would be looser) within 12 iterations.
3. **PERSIST**: reached the target mass; final node residual <= 1e-10; residual at the
   half-node angles (not collocated) <= 1e-8; top-quarter Fourier tail <= 1e-6.
4. **FAIL**: not PERSIST (either the floor was hit, or the final circle fails the off-node or
   tail test), AND the refinement rerun (N raised to the next odd integer >= 1.5 N, first step and
   floor divided by 4) is also not PERSIST. If the rerun persists, the member is a
   **corrector failure (rescued)**, counted as persisting and reported against the corrector.
5. Recorded per member: Tp/T, C, N, mass fraction reached, final and off-node residuals, tail,
   deformation from the mu3 = 0 circle and its amplitude at harmonic q of the nearest listed
   ratio, the per-map unstable multiplier (module's `hyperbolic_bundles`) and the periodic
   orbit's monodromy multiplier.

**4.4 Windows.** A failure agrees with the paper if `|Tp/T - p/q| <= 1.5e-4` for a listed
p/q (a quarter of the smallest gap between listed ratios, 6.8e-4). Twistless bands:
`|Tp/T - (Tp/T)_extremum| <= 2e-5` around each of the three extrema; a failure there is
AMBIGUOUS (the fixed-rho problem is degenerate at a twistless point, independent of resonance)
and is not counted either way. A failure outside both is a DISAGREEMENT unless it lies within
1.5e-4 of a rational with q <= 100 (reported as "unlisted resonance", ambiguous, since the paper
says "including").

**4.5 Members.** A uniform grid in x (step 0.004) over the whole unstable family, the three
extrema, and at every crossing the members with Tp/T = p/q + delta for
delta in {0, +-1e-6, +-3e-6, +-1e-5, +-3e-5, +-1e-4, +-3e-4} that exist on that segment
(207 members, `members.json`).

**4.6 Agreement test (Region A).** (a) no FAIL outside the windows, apart from twistless bands;
(b) for each listed ratio with a crossing in Region A, whether a FAIL occurs inside its window;
(c) FAIL windows of different ratios separated by PERSIST members (the "do not overlap" claim).

**4.7 Predictions (Region A, written before running).** Every Region A member PERSISTS,
including delta = 0. Reason: the crossings in Region A are 0.13 to 0.31 from Rhea's orbit,
where the q-th (q = 31 to 41) harmonic of Rhea's forcing along the orbit should be far below the
1e-10 tolerance, so the island chains are narrower than anything this test can see; a delta = 0
"persist" is then tolerance-limited, not a contradiction, and the deformation at harmonic q is
recorded to show it. This agrees with the paper's "many MID-e orbits ... failed", i.e. not the
low-e ones. Twistless members: expected to persist too if the mass shift of the extremum is
smaller than the distance in rho, otherwise to fail (no prediction made).

**4.8 Method control (not a paper check).** The failure detector has not been shown to detect a
resonance. Control: the 8/41 crossing on segment 1 (x = -0.65704) with all offsets, at
100 x and 1000 x the physical mass. Resonance widths grow like `sqrt(mu3)`. Expected: FAIL at
delta = 0 and small offsets, PERSIST at the largest, and the harmonic-q deformation growing
toward the window edge.

**4.9 Independent checks on a sample of PERSIST members.** (i) closure: each of 16 nodes mapped
5 forcing periods forward in ONE integration of the core 6-state `ccr4bp_eom` (not the module's
batched planar right-hand side), compared with the interpolated circle rotated by `k rho`;
expected error growth bounded by the per-map multiplier to the power k times the residual.
(ii) node refinement: the same member at 1.5 N, curve-to-curve distance <= 1e-8.
(iii) Saturn-centred residual (section 1).

## 5. Results

Run 2026-10-04, outputs `scan.jsonl` (207 members), `refine.jsonl` (73), `control.jsonl`
(24), `posthoc.jsonl` (8), `closure.jsonl` (13) and `summary.json` in
`data/found/886_titan_rhea_torus/`. The pre-registered classes come first; post-hoc
additions are labelled as such.

### 5.1 Counts (pre-registered classification)

| Class | Members |
|---|---|
| persist (scan) | 61 |
| rescued (scan failed, 1.5 N rerun persisted) | 29 |
| fail-floor (step floor hit in both runs) | 32 |
| fail-unresolved (target mass reached but off-node or tail test failed in both runs) | 12 |
| not refined (scan at N = 601 or 701; the 1.5 N rerun exceeds the per-command budget) | 3 |
| Region B (seed gate failed at 701 nodes; not testable) | 70 |

Failures by location (44 in total): 40 inside a window, 2 in a twistless band, 2 outside.

### 5.2 Region A crossings: failure windows

Offsets delta = Tp/T - p/q at which the member did not persist (pre-registered classes; a
`*` marks fail-unresolved members that persisted at 3 N in the post-hoc rerun, section 5.4).

| Ratio, segment (C) | Failed | Persisted | Width of failure |
|---|---|---|---|
| 6/31, seg 0 (3.0547) | 0, +-1e-6, +-3e-6 | +-1e-5 and wider | +-3e-6 |
| 7/36, seg 0 (3.0541) | 0 only | +-1e-6 and wider | below 1e-6 |
| 8/41, seg 0 (3.0513) | 0 only (+-1e-6 rescued at 1.5 N) | the rest | below 1e-6 |
| 8/41, seg 1 (3.0398) | 0, +-1e-6, +-3e-6`*` | +-1e-5 (rescued) and wider | +-1e-6 to 3e-6 |
| 7/36, seg 1 (3.0283) | 0, +-1e-6, +-3e-6, +-1e-5`*`, -3e-4`*` | +-3e-5 and wider (rescued) | +-3e-6 to 1e-5 |
| 6/31, seg 1 (3.0173) | 0 to +-1e-5, +-3e-5`*` | +-1e-4, +3e-4 | +-1e-5 to 3e-5 |
| 6/31, seg 2 (3.0077) | 0 to +-1e-5, +-3e-5 (not re-run at 3 N) | +-1e-4, +-3e-4 | +-1e-5 to 3e-5 |

The three seg-0 failures are TOLERANCE-EDGE failures: the corrector's residual stagnates at
1.0e-10 to 1.2e-10, just above the 1e-10 tolerance, at mass fractions 0.03 to 0.17. They are
real in the sense that the residual cannot be driven below the resonant forcing at exactly
rational Tp/T, but the resonance there is about as weak as the tolerance; a 1e-9 tolerance
would have called them persistent. In segments 1 and 2 (C 3.01 to 3.04) the failures are
not marginal: at |delta| <= 3e-6 the corrector cannot take even the first 1/1024 step, and
the deformation of the circle at harmonic q grows like 1/delta toward the ratio (8/41 seg 1:
1.0e-8 at delta = 1e-4, 9.0e-8 at 3e-5, 3.4e-7 at 1e-5, 1.2e-6 at 3e-6), the small-divisor
signature the method control below also shows.

Twistless points: the T minimum at Tp/T = 0.195273 (grid member at delta_ext = 5.5e-7 and
the extremum member) fails, as the fixed-rho formulation is degenerate there; AMBIGUOUS as
pre-registered. The T maximum at 0.193413 persists (at 1.5 N). The third extremum (0.194080,
N = 601) failed in the scan and was not refined.

### 5.3 Agreement test against the paper

- (a) No failure outside the windows except twistless bands: **not met as pre-registered**,
  2 failures outside: member 127 (7/36 seg 1, delta = -3e-4, fail-unresolved; post-hoc at
  3 N it PERSISTS, so it was a resolution failure of the corrector setting, not of the torus)
  and grid member 44 (x = -0.5637, Tp/T = 0.193939, 3.9e-4 above 6/31; fail-unresolved at 401
  and 603 nodes, off-node residual 1e-6 to 4e-6; nearest q <= 100 rational 19/98 at 6.2e-5;
  not re-run at 3 N because 1203 nodes exceeds the budget). Post-hoc, member 44 is the only
  open item.
- (b) A failure inside the window at every Region A crossing of 6/31, 7/36 and 8/41: **met**
  (7 of 7 crossings). 4/21, 9/47 and 5/26 are crossed only in Region B: **not tested**.
- (c) Failure windows separated by persisting members: **met**. Every failure window is at
  most +-3e-5 wide, while neighbouring listed ratios are 6.8e-4 or more apart, and persisting
  members sit at +-1e-4 around every crossing (except where the segment ends).
- Prediction 4.7 ("every Region A member persists") was **wrong**: the low-e and mid-e
  crossings at C 3.01 to 3.05, 0.13 to 0.31 from Rhea's orbit, do fail inside very narrow
  windows. The paper's text ("many mid-e orbits ... failed") does not say whether such narrow
  failures occurred in its own computation at these energies; it prints no per-orbit list.

### 5.4 Corrector failure against torus failure

- Refinement at 1.5 N and a quarter of the step rescued 29 members: in this problem the
  perturbed circles near a resonance need several times the node count of the mu3 = 0 seed
  (the response is amplified at harmonics near q, 2q, ...), so the seed gate alone does not
  pick a sufficient N.
- Post-hoc (not pre-registered): the 12 fail-unresolved members were rerun at 3 N where that
  was at most 1001 nodes. All 8 that fitted the budget (N = 453) PERSIST, including the
  outside-window member 127 and the in-window members at 8/41 +-3e-6, 7/36 +-1e-5,
  6/31 +-3e-5. So "fail-unresolved" was a resolution failure in every case checked. The
  remaining 4 (6/31 seg 2 at -3e-5, +3e-5 and grid 40, needing 753 nodes; grid 44 needing
  1203) were not run: a 753-node continuation exceeded 470 s per member on the loaded machine.
- The fail-floor members (32) did not converge at 1.5 N with a quarter of the step and a
  floor of 1/4096 of the mass either, and in segments 1 and 2 they fail at the very first
  step: these are the torus failures, all inside |delta| <= 1e-5 except the twistless pair.

### 5.5 Method control (8/41 crossing, segment 1, mass raised)

| Mass | Failed | Persisted |
|---|---|---|
| 100 x | -1e-5 to +3e-6 (fail), +1e-5 (unresolved), +1e-4 (fail, 5e-5 from the twistless T minimum), -3e-4 (unresolved) | -1e-4, -3e-5, +3e-5 |
| 1000 x | every offset | none |

At 100 x the failure window widens from about +-3e-6 (physical mass) to about +-1e-5, as
expected for a width growing like sqrt(mu3) (factor 10). At 1000 x nothing persists (no
refinement was run on the control). The detector does see a resonance and its growth with
mass; the 100 x members at -3e-4 and +1e-5 were not refined, so their status is open.

### 5.6 Independent checks on persisting circles (13 sampled across Region A)

- Closure with the core 6-state `ccr4bp_eom` (one integration per node over 5 forcing
  periods, 16 nodes): maximum errors against the rotated circle 2e-12 to 2e-8 at k = 5;
  growth is consistent with the per-map multipliers (1.0 to 2.5) acting on the off-node
  interpolation error. Pass.
- Node refinement at 1.5 N: 12 of 12 persist; curve-to-curve distance 5e-12 to 5e-9 (the
  pre-registered 1e-8 met). The 13th (x = -0.5477, N = 701) was not refined (1053 nodes).
  The module's `distance_to_circle` could not be used for this: its bounded Brent search in
  theta has a relative tolerance of about 1.5e-8 theta and floors distances at about 2e-8;
  the script uses a Newton projection instead.
- Saturn-centred residual (the paper's perturber centring, Kumar Eq. [4], written in the
  script): 9e-8 to 3e-7 for low/mid-e members, 1.5e-6 at x = -0.548. The circles' deformation
  by Rhea is 1.6e-5 to 1.7e-3, so the centring changes the circle by roughly 1 per cent of
  Rhea's effect; it moves no resonance (frequencies are unchanged) and is not a factor in the
  pattern above.
- Per-map unstable multipliers of persisting circles (module `hyperbolic_bundles`): 1.0008 at
  the stability boundary to 2.50 at the edge of Region A, consistent with the periodic-orbit
  multiplier raised to Tp/T.

### 5.7 Region B (not testable)

70 members (x >= -0.5397, C <= 2.985), including every crossing of 4/21, 9/47, 5/26 and the
segment-3 crossing of 6/31. Estimated node counts for a 1e-10 seed (log-linear fit): about
900 at C = 2.985, 1250 at 2.979 (6/31 crossing), 2000 to 5400 at 2.97 to 2.963 (5/26, 9/47),
8000 to 12000 near 4/21 (C = 2.955). One Gauss-Newton iteration of the dense solve takes
18 s at 1001 nodes and 26 s at 1201 (4 BLAS threads), so only the segment-3 6/31 crossing is
within reach of this corrector at all (estimated 15 to 30 min per member, about 4 to 7 hours
for its 13 offsets). The paper's FFT-based quasi-Newton method is `O(N log N)` per step.

## 6. Verdict

**Partly reproduced.** What was run shows: (1) the three-body family matches every printed
number and figure (max Tp/T to the rounding of the printed Tp); (2) in the part of the family
the corrector can represent (C from 2.99 to 3.055, ratios 6/31, 7/36, 8/41), invariant circles
persist to Rhea's printed mass away from the ratios and fail only within +-3e-5 of them (or at
a twistless point), with non-overlapping failure windows, which is the published pattern;
(3) the failure detector responds to resonance as theory requires (1/delta growth of the
resonant harmonic, wider windows at higher mass). Not shown: the paper's failures in the
high-energy part (4/21, 9/47, 5/26 and the C = 2.98 crossing of 6/31), because the circles
there need 10^3 to 10^4 nodes that a dense least-squares corrector cannot handle.

**Trust in the corrector:** when it converges and passes the off-node and tail tests, the
circle is right (13 of 13 independent closures, refinement agreement to 5e-9). Its failures
need care: 37 of the 76 Region A members that did not persist in the scan were resolution or
step failures (29 rescued at 1.5 N, 8 more at 3 N), and 3 more sit at twistless points or
were not refined, so a single run's
"fail" must not be read as "no torus" without the refinement ladder used here. Its reach is
limited to orbits without close flybys of the base moon (here, distance to Titan above about
0.025, i.e. 30 000 km).

## 7. Open items and the long run

- Post-hoc 3 N on the three 6/31 seg-2 members (753 nodes, about 30 min each; 3 in parallel):
  `uv run python scripts/screen_886_titan_rhea_torus_check.py --stage posthoc --budget-s 7200`
  (member 44 is skipped by the 1001-node cap).
- The segment-3 6/31 crossing at about 1301 to 1401 nodes: needs either several hours of the
  dense solve or a structured (Fourier-diagonal) Newton step, which is also what 4/21, 9/47
  and 5/26 need.
