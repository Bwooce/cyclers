# #574 Stage B -- Titan-Iapetus corrector productization + Saturn V2->V3->V4->V4-strict gauntlet

Date: 2026-07-12. Scripts: `scripts/run_574_stageB_saturn_gauntlet.py` (gauntlet
runner). Modules built: `src/cyclerfinder/genome/titan_iapetus_corrector.py`,
`src/cyclerfinder/data/validation/{v2_saturn_3d,v3_saturn_3d,v4_saturn,
v4_saturn_strict}.py`. Data: `data/gauntlet_574_saturn_stageB.jsonl`. Git sha at run
time: see the jsonl's `_meta` record.

## What this covers

Stage B's remaining scope from `data/OUTSTANDING.md` #574: (1) productize the #572/
#573/#574-Stage-A idealized 3D-closure engine into a reusable, tested corrector module,
and (3) generalize the Uranian V2->V3->V4->V4-strict validation chain
(`v2_moontour.py`/`v3_3d.py`/`v4_uranus.py`/`v4_uranus_strict.py`) to Saturn, running it
on the 15 Stage-A eccentricity-robust survivors. Kernel fetch (item 2) was already done
in a prior run (see the OUTSTANDING #574 "STAGE B KERNEL FETCH" note).

## Part 1: the corrector

`titan_iapetus_corrector.py` holds: `TitanIapetusClosureParams` (the stable 4-free-
parameter input contract -- `omega_deg`, `tof_scale`, `m0_titan_deg`, `m0_iapetus_deg`,
plus fixed `n_rev`/`e_titan`/`e_iapetus`/`inclination_deg`), `kepler_state_3d` (the
eccentric-Keplerian 3D Kepler propagator, C1-compliant -- every later encounter state is
derived by mean-motion propagation from the SAME epoch M0, never a free
re-specification), `evaluate_closure`/`ClosureResult` (the Lambert-closure evaluation),
and `closure_passes_gate` (wraps `candidate_passes_physical_gate` -- the #324 gate --
verbatim). `tests/genome/test_titan_iapetus_corrector.py` pins: the C2 e=0 reduction
positive control (both the Titan in-plane case against `_moon_state` and the Iapetus
inclined case against an independently-hand-written circular-inclined reference, NOT
imported from any throwaway script), a known-branch reproduction against the
already-committed #574 Stage-A jsonl, and structural/input-contract checks.

## Part 2: the Saturn validation chain

Built as genuinely Saturn-specific modules (NOT a reuse of the generic
`v2_moontour.py`/`v3_3d.py`), because the whole point of #574 is that the Titan-Iapetus
candidates' defining model is the ECCENTRIC 3D Kepler model, not the circular-coplanar
`_moon_state` the generic Uranian driver uses -- feeding these candidates through the
generic driver would silently drop the exact fidelity axis #574 was built to test:

* `v2_saturn_3d.py` -- re-solves the same 2 Lambert legs at `k * cycle_period_days`
  advances using the corrector's eccentric 3D states (mirrors `v2_moontour`'s "cycle
  ToFs held fixed" discipline verbatim).
* `v3_saturn_3d.py` -- REBOUND IAS15 (or scipy LSODA fallback) cross-check; reuses
  `v3_3d._ias15_propagate_planet_frame` VERBATIM (it is bare two-body Kepler
  propagation, already fully system-agnostic).
* `v4_saturn.py` -- J2 + 8-moon (Mimas/Enceladus/Tethys/Dione/Rhea/Titan/Iapetus/
  Hyperion) third-body scipy DOP853 fallback. Saturn J2 = 16290.573e-6, R_eq = 60268 km
  (Iess et al. 2019, Science 364(6445), Cassini Grand Finale gravity science). Reuses
  `v4_uranus._j2_acceleration_kms2`/`_third_body_acceleration_kms2`/`_hill_radius_km`
  VERBATIM (fully body-agnostic). Titan/Iapetus perturber positions use the corrector's
  eccentric 3D model; the other 6 moons use circular-coplanar `_moon_state` at a fixed
  zero phase (mirrors the Uranian V4's own documented non-tour-perturber treatment).
* `v4_saturn_strict.py` -- real SAT441 SPICE ephemeris for all 8 moons.
  **Written FRESH with both #567 fixes inherited from the start** (verified, not just
  claimed): (a) continuous Lambert branch selection by ACTUAL propagated terminal
  offset, never a pre-propagation departure-velocity proxy; (b) planet-crossing legs
  pre-screened and TAGGED `FAILURE_MODE_PLANET_CROSSING` with periapsis recorded, never
  silently excluded from the pass-rate denominator -- confirmed firing on a real
  trigger (branch 2 at epoch 2000-09-15, periapsis 32,409-45,000ish km inside Saturn's
  60,268 km R_eq; pinned as a regression test, not a synthetic case).
  `ensure_sat441_kernel()` added to `verify/spice_kernels.py`, mirroring
  `ensure_jup365_kernel`/`ensure_pluto_kernel` exactly.

All four new validation modules + the corrector have tests in
`tests/genome/test_titan_iapetus_corrector.py` (6 tests) and
`tests/data/test_saturn_v2v3v4_gauntlet.py` (8 tests, 2 of which skip gracefully if
SAT441 is absent) -- all pass. `ruff check`/`ruff format --check` clean; `mypy` clean.

## Part 3: the 15-candidate gauntlet -- load-bearing finding

**Every one of the 15 Stage-A survivors FAILS the governing V2 gate.** V2 was run at
`n_cycles` in {3, 5, 10} (the #566 grid). Three candidates (branch ids 1, 10, 16 -- all
`n_rev != (0,0)`) fail to even complete 3 cycles: the SAME `n_rev` Lambert transfer
physically ceases to exist past cycle 0 (a genuine multi-revolution feasibility-window
violation, confirmed at e=0 too -- not an eccentricity effect). The other 12 complete
all `n_cycles` but with a per-cycle V_inf-continuity residual that grows to several km/s
(0.3-6.4 km/s) and a rendezvous drift of ~0.7-2.4 million km -- both far past the
50,000 km drift floor and, when checked across the full {3,5,10} grid, past even the
#566-style 0.5 km/s "quasi-bounded" floor (only branch 6 is quasi-bounded at nc=3 AND
nc=5, but not nc=10). **Root cause (confirmed by direct instrumentation, not assumed):**
unlike the Uranian #558-#569 family -- whose `tof_scale`/`rel_offset_deg` were found by
#563's DEDICATED symmetric/commensurate-closure enumeration specifically so a fixed-TOF
multi-cycle repeat reproduces the SAME encounter geometry -- the Titan-Iapetus #571-#574
closures were found by a free (Omega, tof_scale, rel_offset) search for a SINGLE
V_inf-continuity closure, with NO periodicity/commensurability constraint. Branch 1's
leg-0 transfer angle measured 117.6 deg at cycle 0, 78.0 deg at cycle 1, 38.3 deg at
cycle 2 (checked at e=0, so this is not an eccentricity artifact either). **This means
most of Stage A's 15 "closures" are single V_inf-continuity transfers, not repeating
cyclers, under a literal multi-cycle test** -- a materially different (and more
fundamental) finding than "the family dies on real eccentricity"; it is a
periodicity-formulation gap in how the #571-574 discovery search was originally posed,
orthogonal to the eccentricity question Stage A itself answered.

V3/V4/V4-strict were still computed for the 12 candidates that complete 3 V2 cycles
(the spec's "run the chain on all 15" instruction, and because it's cheap -- <1s/stage/
candidate), for completeness and because they are informative on their own terms (V3
agrees with V2 to near machine precision for all 12, as expected for two integrators
sharing the same analytic model; V4 J2+n-body passes for 9/12; V4-strict at the single
reference epoch 2000-06-21 passes for only 3/12 -- branches 2, 8, 19). **None of these
downstream results change the governing verdict**: since V2 itself never clears PASS or
FAIL_QUASI_BOUNDED for any candidate, the full-chain verdict for all 15 is
`FAIL_AT_V2_*` -- **0/15 reach `PASS`/`PASS_AS_QUASI_CYCLER`.**

### Full per-candidate table

| branch | n_rev | V2 status (3/5/10-cycle grid) | V3 | V4 (J2+8-moon) | V4-strict (2000-06-21) | chain verdict |
|---|---|---|---|---|---|---|
| 1  | (1,1) | FAIL_UNBOUNDED (Lambert infeasible past cycle 0) | -- | -- | -- | FAIL_AT_V2_LAMBERT_INFEASIBLE |
| 2  | (0,0) | FAIL_UNBOUNDED | PASS | PASS | PASS | FAIL_AT_V2_UNBOUNDED |
| 3  | (0,0) | FAIL_UNBOUNDED | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |
| 4  | (0,0) | FAIL_UNBOUNDED | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |
| 5  | (0,0) | FAIL_UNBOUNDED | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |
| 6  | (0,0) | FAIL_UNBOUNDED (quasi-bounded at nc=3,5 only) | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |
| 7  | (0,0) | FAIL_UNBOUNDED | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |
| 8  | (0,0) | FAIL_UNBOUNDED | PASS | FAIL | PASS | FAIL_AT_V2_UNBOUNDED |
| 9  | (1,1) | FAIL_UNBOUNDED | PASS | FAIL | FAIL (converged but agreement/bounded-drift floor) | FAIL_AT_V2_UNBOUNDED |
| 10 | (2,2) | FAIL_UNBOUNDED (Lambert infeasible past cycle 0) | -- | -- | -- | FAIL_AT_V2_LAMBERT_INFEASIBLE |
| 12 | (1,1) | FAIL_UNBOUNDED | PASS | FAIL | FAIL (lambert_no_solution) | FAIL_AT_V2_UNBOUNDED |
| 13 | (0,0) | FAIL_UNBOUNDED | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |
| 16 | (2,2) | FAIL_UNBOUNDED (Lambert infeasible past cycle 0) | -- | -- | -- | FAIL_AT_V2_LAMBERT_INFEASIBLE |
| 19 | (0,0) | FAIL_UNBOUNDED | PASS | PASS | PASS | FAIL_AT_V2_UNBOUNDED |
| 20 | (0,0) | FAIL_UNBOUNDED | PASS | PASS | FAIL (planet-crossing) | FAIL_AT_V2_UNBOUNDED |

**0/15 PASS. 0/15 PASS_AS_QUASI_CYCLER.** Full per-cycle numeric detail:
`data/gauntlet_574_saturn_stageB.jsonl`.

## Epoch-sensitivity spot check (not a full sweep, per scope)

For the 3 candidates whose V4-strict PASSED at the single reference epoch (2, 8, 19), a
6-point spot check (5 months across 2000 + one point in 2015) was run to sanity-check
this wasn't a knife-edge single point (per the #568 duty-cycle framing -- a raw
single-epoch PASS/FAIL is not the final word):

| branch | 2000-01-15 | 2000-04-15 | 2000-06-21 | 2000-09-15 | 2000-12-15 | 2015-06-21 |
|---|---|---|---|---|---|---|
| 2  | PASS | PASS | PASS | FAIL (planet-crossing) | FAIL (planet-crossing) | PASS |
| 8  | FAIL (planet-crossing) | FAIL (planet-crossing) | PASS | FAIL (planet-crossing) | FAIL (planet-crossing) | PASS |
| 19 | PASS | PASS | PASS | PASS | FAIL (planet-crossing) | FAIL (planet-crossing) |

This is genuine epoch/synodic-boundary structure (not a single knife-edge point -- each
candidate PASSes at multiple, non-adjacent sample epochs), consistent with the
#567/#568-established `planet_crossing_infeasible` synodic-aliasing pattern found for
the Uranian family. **This does NOT change the headline result**: all three candidates
already fail the governing V2 multi-cycle gate regardless of V4-strict epoch behavior,
so this spot check is reported for completeness/transparency only, not as grounds for a
different verdict. A full annual/daily epoch-robustness sweep (per #567/#568-style
follow-on work) was explicitly out of scope for this dispatch and was not run.

## Framing (mandatory)

This is quasi-cycler-CLASS evidence about our own idealized + real-ephemeris-tested
search space, same standing as #312's own Uranian family (`FAIL_QUASI_BOUNDED` is a
legitimate class member, not a disqualifier) -- but here the result is a genuinely
different (harsher) verdict: **not one of the 15 Stage-A survivors clears even the V2
multi-cycle gate**, let alone the full V4-strict real-ephemeris chain. This is NOT
evidence the underlying Stage-A eccentric closures are computationally wrong (branch
6's cycle-0 residual is machine precision, reproduced independently by this dispatch's
own test suite) -- it is evidence that the #571-#574 discovery search's parameterization
never targeted multi-cycle-repeat commensurability the way the Uranian #563 search did,
so most "closures" found are single V_inf-continuity transfers rather than periodic
tours. NOT a ballistic-cycler finding and NOT a novelty claim.

## Explicitly NOT done here (out of scope per the #574 Stage B dispatch)

* Opus + Fable adjudication of this result (separate next step, per dispatch scope).
* `data/catalogue.yaml` / `data/empty_regions.jsonl` writeback.
* A full annual/daily V4-strict epoch-robustness sweep (the 6-point spot check above is
  informational only).
* Re-scoping the discovery search to target multi-cycle commensurability (a defensible
  follow-up given this run's root-cause finding, but a new task, not part of #574).
