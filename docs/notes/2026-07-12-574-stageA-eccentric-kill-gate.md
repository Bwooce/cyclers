# #574 Stage A -- Titan-Iapetus eccentric-Keplerian 3D-closure kill gate

Date: 2026-07-12. Script: `scripts/run_574_titan_iapetus_eccentric_kill_gate.py`.
Data: `data/probe_574_titan_iapetus_eccentric_kill_gate.jsonl`. Git sha at run
time: `b987b81`.

## What this tests

#573 found 22 distinct 3D-closure branches (Titan-Iapetus, Saturn) under an
IDEALIZED CIRCULAR model. The one unresolved risk: both moons have real,
non-negligible eccentricity (Titan e~=0.0288, Iapetus e~=0.028 -- JPL SSD
mean orbital elements, per the #574 Stage-A spec text) that the circular
model cannot address. This Stage A gate answers: does the family survive
when both moons are placed on real eccentric Keplerian orbits?

Method: continuation in eccentricity, per the spec's own recommendation over
a fresh grid sweep. Each of the 22 circular branches is its own known
solution; e is stepped 0 -> real value over 6 stages (fractions 0, 0.2, 0.4,
0.6, 0.8, 1.0), refining the 4D free-parameter space (Omega, tof_scale,
M0_Titan, M0_Iapetus -- mean anomaly at epoch for each moon, per C1: no free
per-encounter phase, both moons' every later state is Kepler-propagated by
mean motion from these epoch M0 values) with a small bounded Nelder-Mead
window at each stage, seeded from the previous stage's converged point.

## C2 positive control (mandatory, checked before crediting any verdict)

Two independent checks, both passed:

1. **Analytic reduction smoke test** (`_smoke_test_kepler_reduction`): the
   new eccentric Kepler propagator (`kepler_state_3d`) reduces EXACTLY
   (`dr < 1e-6` km, `dv < 1e-9` km/s) to `_moon_state` (Titan) and the
   circular `iapetus_state_3d` (Iapetus) at e=0, at a grid of M0/Omega/u
   test points. PASS.
2. **Per-branch e=0 continuation stage 0**: all 22 branches, run through the
   SAME 4D-refinement machinery at e=0, reproduce their #573 circular
   residual to well inside the 0.05 km/s gate (all 22 stage-0 residuals are
   at or near machine precision, `< 1e-11` km/s for all but branch 21 which
   nonetheless converges to `2.9e-3` -> refines cleanly). 0/22 smoke
   failures.

**C2 verdict: PASS.** The eccentric machinery is trusted.

## Main-population result (17 eccentricity-robust branches)

| branch | n_rev | circ V_inf(Iap) | final V_inf(Iap) | drift | circ bend | final bend | survives |
|---|---|---|---|---|---|---|---|
| 1  | (1,1) | 1.378 | 1.339 | -0.039 | 8.11  | 8.56  | YES |
| 2  | (0,0) | 1.380 | 1.448 | +0.069 | 8.09  | 7.39  | YES |
| 3  | (0,0) | 1.590 | 1.539 | -0.051 | 6.20  | 6.59  | YES |
| 4  | (0,0) | 1.574 | 1.335 | -0.239 | 6.32  | 8.60  | YES |
| 5  | (0,0) | 1.546 | 1.499 | -0.047 | 6.53  | 6.93  | YES |
| 6  | (0,0) | 1.204 | 1.177 | -0.027 | 10.38 | 10.83 | YES |
| 7  | (0,0) | 1.193 | 1.156 | -0.037 | 10.57 | 11.20 | YES |
| 8  | (0,0) | 1.617 | 1.655 | +0.038 | 6.00  | 5.74  | YES |
| 9  | (1,1) | 1.420 | 1.677 | +0.257 | 7.66  | 5.60  | YES |
| 10 | (2,2) | 1.335 | 1.245 | -0.091 | 8.60  | 9.78  | YES |
| 12 | (1,1) | 1.341 | 1.295 | -0.046 | 8.53  | 9.09  | YES |
| 13 | (0,0) | 1.475 | 1.436 | -0.039 | 7.18  | 7.51  | YES |
| 16 | (2,2) | 1.482 | 1.449 | -0.033 | 7.08  | 7.39  | YES |
| 19 | (0,0) | 1.380 | 1.360 | -0.020 | 8.09  | 8.31  | YES |
| 20 | (0,0) | 1.450 | 1.329 | -0.121 | 7.37  | 8.67  | YES |
| 0  | (2,2) | 1.547 | 2.740 | +1.194 | 6.53  | 2.16  | **NO** |
| 14 | (0,0) | 1.574 | 1.839 | +0.265 | 6.32  | 4.69  | **NO** |

**Raw survivors: 15/17. Deduped survivors: 15 (no two branches converged to
the same eccentric point -- every survivor stayed in its own dedup
singleton).** N_rev classes spanned: 3 (all of (0,0)/(1,1)/(2,2) have >=1
surviving deduped branch). 15 >= the pre-registered PASS bar of 5 with wide
margin.

Branches 0 and 14 die on the #324 bend gate at real eccentricity (residual
stays near machine precision -- these ARE genuine closures, just with the
Iapetus asymptote bend falling below the 5 deg floor). Both show the
LARGEST V_inf drift under the continuation (0: +1.19 km/s; 14: +0.27 km/s),
and both drift SMOOTHLY and monotonically stage-by-stage (not a single-step
jump), consistent with a real physical effect rather than the optimizer
hopping Lambert branches: branch 0 has the longest TOF in the whole
population (tof_scale=2.80, n_rev=(2,2)); branch 14's tof_scale (0.560) is
around the middle of its n_rev=(0,0) class but sits at a basin (mirror-paired
with branch 4, see below) that also shows above-median drift. The pattern is
consistent with more elapsed TOF exposing more mean-anomaly-dependent
orbital-radius variation, i.e. a real "eccentricity hurts longer-TOF/higher-
n_rev branches more" signal, not an artifact -- but this is a 2-branch signal,
not independently re-verified past the base gate check.

## Floor-hugger control (ids 11, 15, 17, 18, 21; bends 5.28/5.37/5.02/5.32/5.02 deg)

| branch | n_rev | circ bend | final bend | survives |
|---|---|---|---|---|
| 11 | (1,1) | 5.28 | 5.00 | YES |
| 15 | (1,1) | 5.37 | 5.64 | YES |
| 17 | (1,1) | 5.02 | 5.60 | YES |
| 18 | (0,0) | 5.32 | 2.90 | **NO** |
| 21 | (1,1) | 5.02 | 6.17 | YES |

**4/5 floor-huggers survived** (deduped: 4, all distinct). This is ABOVE the
pre-registered "if >=3 survive, the ecc-robust proxy is non-discriminating"
flag threshold. **FLAGGED EXPLICITLY AND PROMINENTLY per the spec's own
mandatory control-check instruction: the 6.0 deg eccentricity-robust proxy
used to build the #573 17/5 split is NOT a reliable discriminator of which
branches survive real eccentricity** -- 4 of the 5 branches predicted most
likely to die (all sitting right at the circular 5 deg gate floor) instead
comfortably cleared the gate under real eccentricity, in one case (21) even
GAINING over a degree of bend. This does not overturn the PASS verdict per
the pre-registered primary threshold (>=5 deduped survivors among the 17
ecc-robust branches, independent of the control-check outcome), but it
materially discounts confidence that "ecc-robust bend >= 6 deg" is doing
useful predictive work -- the main-population 15/17 survival rate should be
read as "most closures near this circular family survive a real, but still
idealized, eccentricity perturbation," not as evidence the 6 deg proxy
specifically identified the survivors.

## Mirror-pair check (C3 dedup extension to (M0_Titan, M0_Iapetus))

The two known near-mirror pairs from #573 (identical circular residual/
V_inf, related by an exact Omega+180/u+180 transform):

- **{branch 2, branch 19}: both survive, but DID NOT MERGE** under the
  extended (M0_Titan, M0_Iapetus) dedup criterion (final Omega 189.6 deg vs
  6.6 deg -- still ~180 deg apart, so NOT literally the same point; final
  V_inf 1.448 vs 1.360 km/s, a real ~6% difference that opened up under
  eccentricity).
- **{branch 4, branch 14}: NOT_BOTH_SURVIVING** -- 4 survives (bend rose to
  8.60 deg), 14 died (bend fell to 4.69 deg).

Both outcomes match the a priori physical prediction written into the
script docstring: the (Omega+180, u+180) transform is an EXACT degeneracy
only for a constant-radius (circular) orbit, since it leaves a circular
orbit's (x, y) position unchanged while flipping z. Under eccentricity,
`r(nu) != r(nu+180)` in general (nu=0 is fixed at periapsis under this
script's omega=0 convention, nu=180 at apoapsis), so the transform is no
longer an exact symmetry of the eccentric problem -- the two "mirror"
branches are free to (and did) diverge once continuation is run. This is a
positive, mechanistically-understood result, not a surprise or a bug.

## Verdict

**Deduped robust survivors: 15 (>= 5 PASS bar), spanning 3/3 n_rev classes
(>= 2 PASS bar). VERDICT: PASS**, carried with the explicit floor-hugger
control caveat above (the discriminating power of the >=6.0 deg circular
bend-margin proxy specifically is not established by this run, even though
the population-level survival rate clears the bar by a wide margin).

Per #574 Stage A scope: this is quasi-cycler-class evidence only, same
standing as #312's own family (V2 fails on drift by design) -- NOT a
ballistic-cycler finding and NOT a novelty claim, an internal fact about our
own enumerated (now eccentric-Keplerian, still non-ephemeris) search space.

**Explicitly NOT done here (out of scope):** Stage B (productized corrector,
SPICE SAT441 fetch, V-gauntlet, catalogue.yaml writeback), the Opus/Fable
A->B transition adjudication (C6) that would be a prerequisite before any
Stage B dispatch, and resolving/adjudicating the floor-hugger non-
discriminating flag beyond reporting it.
