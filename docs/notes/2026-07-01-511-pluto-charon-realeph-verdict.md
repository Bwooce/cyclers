# Task #511: Pluto-Charon (3,2) real-ephemeris DIFFERENTIAL-CORRECTION verdict

**Date:** 2026-07-01
**Row:** `ross-rt-pc-cycler-32-2026` (V2-ballistic since #505)
**Task:** build the PROPER real-eph lever #506 scoped but did not attempt --
differentially correct the CR3BP (3,2) periodic orbit into a real-ephemeris
analog (NOT the naive propagation test #506 already rejected).
**Status:** BUILT AND RUN. Clean, well-characterized negative. **Row stays
V2-ballistic.** No catalogue edit (out of scope; task is characterization).

---

## 0. TL;DR

**No strictly-periodic real-ephemeris analog of the PC (3,2) orbit exists,
for a STRUCTURAL reason, not merely because Charon's real orbit is slightly
eccentric.** The orbit's period is `T = 11.8334625 TU = 1.8834` Charon
orbital periods -- not an integer. Under real (non-circular) Pluto-Charon
dynamics, exact periodicity in a frame that tracks Charon's true anomaly
requires the spacecraft period to be an integer multiple of Charon's own
period (the "period_f trap", `#441`). Because `1.8834` is nowhere near an
integer, the differential corrector hits this trap: the symmetric
half-period residual "converges" (`2.5e-12`, looks like success) but an
**independent** full-orbit re-propagation shows the orbit does **not**
actually close -- residual `1.726e-2` nd = **338 km**, 74x the radial
excursion (4.6 km) Charon's real eccentricity alone would produce. This is a
genuine physics/topology mismatch between this orbit's `(3,2)` encounter
cadence and Charon's real periodicity, confirmed numerically on real
SPICE-sourced input, not an infra gap.

---

## 1. What #506 asked for, and what was actually missing

#506 (same day) characterized two gates:
- **Gate (a) kernel availability** -- CLEARED by #510 (`plu060.bsp` fetched
  and verified, same NAIF PLU060 source `satellites.py` already cites for
  Charon's GM).
- **Gate (b) methodology** -- #506 explicitly said a naive real-eph
  propagation is the WRONG test (it only measures Charon's eccentricity
  model-mismatch oscillation, ~756 m, 5800x the V2 drift signal -- noise, not
  discrimination) and named the RIGHT test: *"differentially correcting the
  CR3BP IC in the real-eph model to find the real-eph analog of the periodic
  orbit ... a new, harder computation."* This task builds exactly that.

## 2. Real Charon eccentricity from SPICE (plu060.bsp)

`ensure_pluto_kernel()` (new, `verify/spice_kernels.py`, mirrors
`ensure_jup365_kernel` exactly -- local-only 129 MB, `RuntimeError` with the
NAIF URL if absent, no auto-download) resolves
`~/dev/references/kernels/plu060.bsp`.

`charon_osculating_elements()` (new, `verify/pluto_charon_realeph.py`) reads
Charon (NAIF 901) relative to Pluto (999), J2000 frame, and computes classical
two-body osculating elements (vis-viva energy + Laplace-Runge-Lenz
eccentricity vector) using `GM = 975.5 km^3/s^2` (the Pluto+Charon SYSTEM GM
`satellites.py` already uses for `mu = 0.10876473603280369` -- same-source,
zero model-mismatch risk per #510's mass-consistency check).

| Epoch | a (km) | e (osculating) | T (days) |
|---|---|---|---|
| 2000-01-01 | 19594.301 | 2.2807e-4 | 6.386268 |
| 2015-07-14 | 19594.304 | 2.2651e-4 | 6.386270 |
| **2026-07-01** | **19594.305** | **2.3306e-4** | **6.386270** |
| 2050-01-01 | 19594.307 | 2.1652e-4 | 6.386271 |

This is the FIRST time this project has measured Charon's real eccentricity
directly from ephemeris rather than citing a published bound. It is
epoch-stable (spread `1.65e-5` across 50 years) and **~4.6x #506's
back-of-envelope `e < 5e-5`** (Brozovic et al. 2015's MEAN/forced
eccentricity) -- the osculating value folds in short-period terms from
Nix/Hydra/Kerberos/Styx and solar perturbation that the mean-element bound
excludes. Both are tiny; `a = 19594.3` km matches the catalogue's `19600` km
CR3BP length unit to 0.03%.

## 3. The differential correction: ER3BP e-continuation

`differential_correct_pc32_to_eccentricity()` (new,
`verify/pluto_charon_realeph.py`) bridges the CR3BP (3,2) seed (`#494`
C-sweep nu=0 midpoint) into the ER3BP pulsating frame at `e=0` (exact per
`#441` Sec. 1 -- ER3BP IS CR3BP there), then walks the existing,
independently-validated (`#293`, Fitzgerald positive control + Broucke
e-continuation) secant e-continuator
(`genome/er3bp_continuation.continue_er3bp_family_in_e`) from `e=0` to the
real `e=2.3306e-4`.

**Actual run** (`scripts/pc_v3_realeph_correction.py`, 2026-07-01):

```
Charon real osculating orbit (2026-07-01T00:00:00):
  a=19594.30451 km, e=2.330574e-04, T=6.386270 d
CR3BP seed: T_nd=11.8334625171, T/(2*pi)=1.883354
e=0 bridge: corrector_res=9.699e-12, independent_res=3.121e-09
e=2.331e-04 target: corrector_res=2.494e-12, independent_res=1.726e-02 (gate 1e-08)
CONVERGED (real-eph analog found): False
```

## 4. Why: the period_f trap, confirmed on real data

The ER3BP pulsating frame is `2*pi`-periodic in true anomaly `f` (the scale
factor `1/(1+e*cos f)` depends on where Charon sits in its own orbit). A
strictly periodic orbit in this frame therefore requires its period, in `f`,
to be an integer multiple of `2*pi` -- i.e. the spacecraft orbit must be
commensurate with Charon's own orbital period. This is a documented, durable
finding from `#441` ("Phase 2 bridge spike"), previously demonstrated only on
synthetic Phase-1 CR3BP seeds (Broucke family, up to `e=0.30`, residual grew
to `1.44`).

PC (3,2)'s catalogue-documented topology is **"3 Pluto + 2 Charon passes per
period"** -- a statement about encounter counts within ONE spacecraft period,
not a claim that the spacecraft period is commensurate with Charon's orbital
period. Measuring it: `T_nd / (2*pi) = 1.8834` -- **0.117 short of the
nearest integer (2)**, i.e. genuinely incommensurate.

The numerical run reproduces the exact #441 signature on this real,
catalogue-sourced orbit for the first time:
- **Corrector residual converges** (`2.5e-12`): the symmetric half-period
  crossing condition `(y=0, xdot=0)` is a LOCAL condition and is satisfiable
  at any nearby `e` regardless of commensurability -- this is why a naive
  "did the corrector converge?" check would give a FALSE POSITIVE.
- **Independent full-orbit residual does NOT converge** (`1.726e-2` nd =
  **338 km**): re-propagating over the full period with an independent
  (Radau) integrator shows `X(T) != X(0)` by 338 km.
- **338 km is 74x the 4.6 km radial excursion** Charon's real eccentricity
  alone would produce (`a*e = 19594*2.33e-4`) -- confirming this is the
  STRUCTURAL incommensurability effect, not merely an amplified version of
  #506's eccentricity-oscillation concern. Even a perfectly circular-but-
  displaced-phase model would show a comparable gap; the effect is
  topological, not perturbative.

## 5. V2 -> V3 recommendation

**STAYS V2-ballistic. Real-eph analog orbit at this exact topology does not
exist** (to the precision tested; `independent_tol=1e-8`, actual `1.7e-2`,
9+ orders above tolerance -- not a marginal near-miss).

This is a genuine, positive scientific finding, not an infra gap:
1. **Both #506 gates are now fully resolved.** Gate (a) cleared (#510 kernel
   fetch). Gate (b) is resolved by actually building and running the harder
   computation #506 named, rather than leaving it as a scope note.
2. **The blocker is physics/topology, not a missing capability.** The
   differential corrector, the e-continuation infrastructure, and the
   independent-closure gate are all landed, tested, and reusable for any
   FUTURE Pluto-Charon row whose period happens to be commensurate with
   Charon's own period (e.g. a `k:1` resonant PC cycler, which #442-style
   connected-family continuation could target directly).
3. **The result generalizes the #441 period_f-trap finding** from a
   synthetic Phase-1 probe to a real catalogue row with real SPICE input --
   strengthening confidence that the trap is a genuine ER3BP structural
   feature, not a corrector-convention artifact specific to one seed.
4. **V2-ballistic remains the strongest applicable evidence** for this row:
   100-period REBOUND/IAS15 bounded-drift in the row's DEFINING model
   (CR3BP), independent integrator + frame (#505). The real-eph attempt does
   not add a stronger positive claim; it adds a precisely-quantified boundary
   on where the CR3BP idealization stops being physically realizable exactly.

## 6. Deliverables

- `src/cyclerfinder/verify/spice_kernels.py::ensure_pluto_kernel()` -- kernel
  resolver, mirrors `ensure_jup365_kernel` exactly.
- `src/cyclerfinder/verify/pluto_charon_realeph.py` -- new module:
  `charon_osculating_elements()` (SPICE osculating-element extraction) +
  `differential_correct_pc32_to_eccentricity()` (ER3BP e-continuation
  wrapper with the independent-closure gate surfaced, not hidden).
- `scripts/pc_v3_realeph_correction.py` -- runnable driver (mirrors
  `scripts/pc_v2_longspan.py` conventions), reuses the sourced CR3BP
  constants from `pc_v2_longspan.py` (no duplication).
- `tests/verify/test_511_pc_realeph_correction.py` -- 4 tests: 2 kernel-free
  (structural period-ratio regression pin + exact e=0 bridge), 2 kernel-gated
  (skip cleanly without `plu060.bsp`, per the just-fixed CI convention):
  real-eccentricity sanity band + the headline non-closure result. **4/4
  pass** (`171 s`, dominated by the `n_steps=20` continuation's many
  full-period STM integrations at `T_nd~11.83`).
- No catalogue edit.

## 7. Discipline anchors

- `feedback_check_dont_guess` -- every number in this note (osculating
  elements, residuals, ratios) came from an actual run, captured in
  `scripts/pc_v3_realeph_correction.py`'s output, not estimated.
- `feedback_verify_gauntlet_with_positive_control` -- the e=0 bridge
  (`independent_res=3.1e-9`) is the positive control proving the corrector
  +gate correctly recognizes a genuinely-closing orbit; the e=2.33e-4 result
  is judged by the SAME independent-closure criterion, not a weaker one.
- `feedback_bugfix_invalidates_past_searches` -- N/A here (no bug fixed);
  this is a NEW capability applied for the first time, not a re-run.
- `project_negative_results_registry` -- "no real-eph analog at this
  topology" is method-versioned: it is conditional on the `(3,2)` orbit's
  specific (non-commensurate) period, not a general Pluto-Charon real-eph
  verdict. A future `k:1`-resonant PC row is NOT pre-judged by this result.
- `feedback_orbit_closure_discipline` -- judged on the INDEPENDENT residual,
  never the corrector's own (locally-blind) residual; "it converged!" was
  the exact danger signal this note's Section 4 dissects rather than trusts.
