# 2026-06-17 - #344 Phase 2 Stage B verdict on the Saturn Titan-Rhea-Titan SILVER

## Stage B scope

#344 Phase 1 (commit `781e2ee`) identified a 0.0102 km/s closure on the Saturn
Titan-Rhea-Titan (1,1) repeated-moon cycle at `phase0=273.75°`,
`rel_offset=288.75°`, `tof_scale=2.0`, `V_inf=[1.7375, 1.6463, 1.7273] km/s`.
Stage A (commit `63809ec`) verified the IC, the post-#346 anchor count = 0
(`literature_check`), the #324 physical-sanity gate, and the #256 ML flagger,
and returned `PASS_PROCEED_TO_STAGE_B`.

Stage B runs the same V2-moontour gauntlet
(`cyclerfinder.data.validation.v2_moontour.run_v2_moontour`) that #330 ran on
the #327 Umbriel-Oberon-Umbriel SILVER: re-solve the SILVER's Lambert legs
over `n_cycles` consecutive cycles, advancing the moon ephemerides through
their natural Keplerian motion across each cycle. Strict gates: every cycle's
Lambert converged AND per-cycle V_inf-continuity residual <= 0.05 km/s AND
inter-cycle rendezvous drift <= 50,000 km.

## Phase convention translation

The Stage A IC was mined with `_sweep_one_cycle`'s phasing convention
(`scripts/scan_344_saturn_titan_rhea_finer.py`):
`theta = {anchor: phase0, intermediate: phase0 + rel_offset}`. The `v2_moontour`
driver uses `theta_base = {sorted-first: phase0, sorted-second: phase0 + rel_off}`.

For `("Titan", "Rhea", "Titan")`:

* Stage A places Titan (anchor) at `273.75°` and Rhea (intermediate) at
  `(273.75 + 288.75) mod 360 = 202.5°`.
* `sorted({"Titan", "Rhea"}) = ("Rhea", "Titan")`, so the v2 driver places
  Rhea at `phase0_v2` and Titan at `phase0_v2 + rel_off_v2`.
* Conversion to reproduce Stage A geometry under v2 convention:
  `phase0_v2 = 202.5°` (Rhea), `rel_off_v2 = 71.25°` (Titan minus Rhea).

The #330 Umbriel-Oberon SILVER had `rel_offset = 180°` which is palindromic
under the swap, so the convention difference vanished there; for #344 the
swap is load-bearing. Cycle-0 closure under the v2 driver reproduces the
Stage A residual `0.0101881 km/s` exactly when the conversion is applied
(test `test_silver_344_v2_cycle_zero_reproduces_stage_a_closure`).

## SILVER V2 verdict at `n_cycles = {3, 5, 10}`

Stored fields driving the run (sourced from `data/silver_344_verified.jsonl`,
commit `63809ec`):

```
sequence        = ('Titan', 'Rhea', 'Titan')
V_inf (km/s)    = (1.7375055995850324, 1.6462740278228238, 1.7273175030110421)
leg ToFs (days) = (16.977266455394638, 16.977266455394638)   # cycle = 33.954 d
rel_offset (deg, Stage A) = 288.75    # gate-passing basin floor (#344 Part A.2)
phase0 (deg, Stage A)     = 273.75
rel_offset (deg, v2 convention) = 71.25
phase0 (deg, v2 convention)     = 202.5
n_rev           = (1, 1)
```

V2-moontour verdict at the three `n_cycles` scales (sourced from
`data/silver_344_moontour_v2_verdicts.jsonl`):

| n_cycles | passes_v2 | n_completed | max_drift (km) | max_closure (km/s) |
|----------|-----------|-------------|----------------|--------------------|
| 3        | **False** | 3 / 3       | 1.772e+06      | 5.444              |
| 5        | **False** | 5 / 5       | 2.441e+06      | 5.444              |
| 10       | **False** | 10 / 10     | 2.441e+06      | 5.444              |

Per-cycle trace at `n_cycles = 10`:

```
cycle 0: drift=0.000e+00 km, residual=1.019e-02 km/s, converged_legs=2/2
cycle 1: drift=9.641e+05 km, residual=5.444e+00 km/s, converged_legs=2/2
cycle 2: drift=1.772e+06 km, residual=2.663e+00 km/s, converged_legs=2/2
cycle 3: drift=2.292e+06 km, residual=9.833e-01 km/s, converged_legs=2/2
cycle 4: drift=2.441e+06 km, residual=4.262e+00 km/s, converged_legs=2/2
cycle 5: drift=2.193e+06 km, residual=1.657e-01 km/s, converged_legs=2/2
cycle 6: drift=1.590e+06 km, residual=4.000e+00 km/s, converged_legs=2/2
cycle 7: drift=7.290e+05 km, residual=3.880e+00 km/s, converged_legs=2/2
cycle 8: drift=2.503e+05 km, residual=1.866e-01 km/s, converged_legs=2/2
cycle 9: drift=1.189e+06 km, residual=4.286e+00 km/s, converged_legs=2/2
```

Verdict label: **`FAIL`** - the SILVER fails the strict V2 gates and ALSO
fails the v4.7 quasi-cycler envelope. The cycle-0 closure (0.0102 km/s) is
the Stage A SILVER residual reproduced exactly, but by cycle 1 the closure
residual has jumped to 5.44 km/s (109x the strict floor of 0.05 km/s and
10.9x the v4.7 quasi envelope of 0.5 km/s). The closure oscillates - cycles
5 and 8 dip back below 0.5 km/s (0.166 and 0.187 km/s respectively) - but
the intervening cycles peak at 4-5 km/s, so the trajectory is structurally a
**closure-divergent oscillator**, not a bounded near-resonant quasi-cycler.

### Reading the trace

The drift hits ~9.6e5 km at cycle 1 (19x the 50,000 km strict floor) and
peaks at ~2.4e6 km by cycle 4, then oscillates between 2.5e5 km and 2.4e6
km without ever returning below the strict floor. The pattern resembles a
2:1-ish near-commensurable phasing - drifts and closures briefly return
near the cycle-0 phasing at cycles 5 and 8 - but the recovery is partial
and the closure residual at peaks is ~30x what the #330 Umbriel-Oberon
SILVER showed (#330 max 0.349 km/s vs #344 max 5.44 km/s, a 15.6x ratio).

The cycle period is 33.95 d. With Titan SMA 1,221,870 km and Rhea SMA 527,070
km (from `src/cyclerfinder/core/satellites.py`) and Saturn GM 3.7931207e7
km^3/s^2, the moons' sidereal periods are 15.948 d (Titan) and 4.518 d
(Rhea), giving a Titan-Rhea synodic period of
`1 / (1/4.518 - 1/15.948) = 6.304` d. So 33.95 d / 6.304 d = 5.386 x synodic
- close to a 5:1 near-resonance but at 7.72% offset from exact 5:1, far
from the #330 Umbriel-Oberon 0.18% offset from exact 5:1
(cycle 29.88 d / synodic 5.987 d = 4.991, `(5 - 4.991)/5 = 0.18%`).

The Lambert legs continue to converge at every cycle - the geometry is
solvable - but the V_inf-continuity envelope diverges quickly and never
heals to the quasi-cycler floor.

## Comparison vs the #330 Umbriel-Oberon precedent

| Metric (n_cycles = 10)           | #330 Umbriel-Oberon | #344 Titan-Rhea-Titan |
|----------------------------------|---------------------|-----------------------|
| Stage A IC residual (km/s)       | 0.025               | 0.0102                |
| Cycle-0 reproduction (km/s)      | 0.025               | 0.0102                |
| Max drift (km)                   | 5.30e+05            | 2.44e+06              |
| Max closure (km/s)               | 0.349               | 5.444                 |
| Cycles with closure > 0.5 km/s   | 0                   | 7 of 10               |
| Drift returns near seed by ?     | cycle 5 (~86k km)   | cycle 8 (~2.5e5 km)   |
| Resonance offset from N:1        | 0.18% off 5:1       | 7.7% off 5:1          |
| Verdict label                    | FAIL_QUASI_BOUNDED  | FAIL                  |
| Admission slot                   | quasi_cycler (v4.7) | none                  |

Key structural difference: the #330 SILVER kept the V_inf-continuity envelope
inside 0.5 km/s for ALL 10 cycles (Lambert geometry "wobbled" but the cycler's
identity persisted); the #344 SILVER's V_inf-continuity envelope ruptures
from cycle 1 - the trajectory may briefly recover at cycles 5 and 8, but
the recovery is partial and the peaks are 10-15x the quasi envelope ceiling.

## Stage C recommendation

**HALT** - the candidate does NOT proceed to Stage C (V3 REBOUND IAS15). The
honest verdict is that the 0.0102 km/s closure at Stage A is real but
**fragile**: it does NOT generate a near-resonant tour that survives even a
single cycle of Keplerian moon drift in the planet-frame patched-conic model.
A V3 real-ephemeris IAS15 run would inherit a far worse starting point than
either the #330 SILVER OR a healthy near-resonant tour, and would compound
the geometry's already-broken V_inf-continuity envelope.

The Stage A SILVER closure looks like a single-cycle resonance accident:
the specific Titan/Rhea longitudes at `phase0=273.75°` give a deep planar
Lambert closure, but the geometry's underlying near-resonance (`5.385x`
Titan-Rhea synodic) is far enough off exact 5:1 that the phasing drifts
out of the closing basin within one synodic cycle.

Per the `feedback_orbit_closure_discipline` "clean negative is success" rule
and the `project_negative_results_registry` method-versioned anti-catalogue
discipline, the candidate retires to the negative-results registry:

* Method: V2-moontour (`v2_moontour.run_v2_moontour` at the spec-§14 strict
  floors).
* Conditional empty: at this phasing under the planar circular-coplanar
  patched-conic model, the SILVER does NOT support a 3-cycle bounded
  trajectory.

A future re-sweep is warranted IF:

1. A method with **richer ephemeris** (BCR4BP including Saturn oblateness;
   real-eph DOP853 / IAS15 with SPICE Saturn kernels) finds the trajectory
   stabilises - this would subsume the planar-CR3BP / patched-conic empty.
   Cost is high enough to not chase speculatively from a closure-divergent
   V2 verdict.
2. A finer Stage A phasing sweep finds a phase0/rel_offset basin with
   **bounded** V2 trajectory at the strict floors - the current Stage A
   record is the deepest of the post-#346 ps=96 scan; a ps=192 sweep or
   continuation off this seed would test whether the closure-divergent
   pattern is robust to phasing.

## Discipline anchors

* NO catalogue writeback. V2 FAIL does NOT admit, and Stage E
  `_LEVEL_EVIDENCE` registration is deferred until/unless a future Stage
  passes.
* The V2 verdict is whatever the math says. The SILVER's closure residual
  blows past the strict floor at cycle 1 (109x), past the v4.7 quasi envelope
  at cycle 1 (10.9x), and peaks at 5.44 km/s. The honest verdict is FAIL.
* `feedback_orbit_closure_discipline`: clean negative is success.
* Lambert geometry uses the same kernel Stage A closed under (cycle-0
  reproduces the stored 0.0102 km/s tightly to ~1e-9 km/s).

## Artifacts

* `scripts/run_344_stage_b_v2_moontour.py` - the Stage B driver
* `data/silver_344_moontour_v2_verdicts.jsonl` - the verdict JSONL
* `tests/verify/test_silver_344_v2_quasi_cycler.py` - the frozen-gate pytest
  (5 tests, all passing)
* this note

Run the driver via:

```
uv run python scripts/run_344_stage_b_v2_moontour.py
```

Verify with:

```
uv run pytest tests/verify/test_silver_344_v2_quasi_cycler.py -v
```
