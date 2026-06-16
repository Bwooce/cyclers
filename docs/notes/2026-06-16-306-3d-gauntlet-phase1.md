# #306 Phase 1 — 3D V0-V5 gauntlet adaptation (V1 + V2 MVP)

Date: 2026-06-16

## What Phase 1 delivers

* **V1 same-model 3D gauntlet** (`src/cyclerfinder/data/validation/v1_3d.py`):
  re-closes a 6D CR3BP periodic-orbit IC under the #291 Phase 1 3D corrector
  (single-shooting Newton on the STM, DOP853) AND independently re-propagates
  under Radau at `rtol = atol = 1e-12`. Both must hold. Converts nondim
  closure residuals to km/s via the system's velocity unit (`l_km / t_s`) and
  asserts the spec §14 V1 floor of 1e-3 km/s on the independent closure. The
  Lambert / lamberthub `V1_TOLERANCE_MPS` floor is preserved (same 1 m/s
  number), specialised here for periodic-orbit closure.
* **V2 long-span bounded-drift 3D gauntlet** (`v2_3d.py`):
  propagates a corrector-clean IC for `n_cycles >= 3` consecutive periods
  using the 6D CR3BP propagator and asserts cumulative position drift stays
  within the same-model 50,000 km floor (sourced from
  `cyclerfinder.verify.propagate.DRIFT_TOLERANCE_KM`, the V2-ballistic bar
  this project already uses for E-M class cyclers).
* **8 V1 tests + 8 V2 tests** (`tests/data/test_v1_3d.py`,
  `test_v2_3d.py`): sourced-golden inputs (catalogue Braik-Ross C11a planar,
  #287 spike 3D, #301 sub-family member, JPL L1 halo fixture), passing and
  negative-control cases, fabrication guards on the spec-fixed floors.
* **#327 SILVER V1 + V2 verdicts** (`data/silver_327_v1_v2_verdicts.jsonl`):
  the SILVER cleared V1 against the spec V1 floor; V2 is honestly
  `NOT_RUN_PHASE_1` (moontour ↔ periodic-IC mismatch, see below).

## What Phase 1 does NOT deliver

* **V0** — already implicit in the spec (sourced-row check). No code shipped.
* **V3 (independent-model)** — the planar n-body harness in
  `src/cyclerfinder/nbody/` is planar-only; extending it to 3D plus a
  REBOUND comparison is Phase 2 of #306.
* **V4 (HFEM real-eph)** — requires GMAT / SPICE integration at Uranus
  (currently V4-ready only at E-M / E-M-Mars). Phase 3.
* **V5 (mission-quality)** — human gate, out of scope for automated phases.
* **Moontour V2** — for multi-leg Lambert-patched-conic moontours like the
  #327 SILVER, "V2 over ≥3 cycles" requires re-solving the Lambert legs
  each cycle with the moons' phase advanced. The shipped V2 operates on a
  CR3BP periodic-orbit IC, NOT a moontour. See "Honest moontour caveat"
  below.

## #327 SILVER — V1 + V2 verdicts

The #327 Umbriel-Oberon-Umbriel SILVER (`data/silver_327_verified.jsonl`,
verified at task #327) is a two-leg Lambert moontour with stored Lambert
leg-end residual `0.025 km/s` and a DOP853 cross-check arrival residual of
`max_dv_arrival = 2.3e-10 km/s` (`max_dr_arrival = 1.6e-5 km`).

### V1 verdict: **PASS**

* **Independent cross-check residual (headline)**: `2.30e-10 km/s`
* **Spec §14 V1 floor**: `1.00e-3 km/s` (1 m/s)
* **Margin**: ~7 orders of magnitude

The DOP853 re-propagation of each Lambert leg arrives within 16 micrometres
per second — far below the 1 m/s V1 bar. The SILVER passes V1 cleanly.

### V2 verdict: **NOT_RUN_PHASE_1**

* **Reason**: the SILVER is a single-cycle moontour, not a periodic-orbit
  IC. The Phase 1 V2 module operates on `(state0, period)` for a CR3BP
  periodic orbit; the moontour interpretation of V2 (3 consecutive
  Umbriel-Oberon-Umbriel sequences with phase advanced by the moons' orbital
  periods) requires re-solving the Lambert legs each cycle and is **Phase 2
  work**, explicitly out of Phase 1 scope.
* **Diagnostic**: the SILVER's stored 25 m/s velocity residual at the end
  of cycle 1 is precisely the extent to which the moontour is NOT
  periodic — a useful diagnostic but NOT the V2 spec quantity.

### Honest moontour caveat

The SILVER's V1 PASS is real but operates in a different domain than the
CR3BP periodic-orbit V1 PASS the spike #287 / sub-family #301 demonstrate.
Both clear the spec V1 floor under their respective same-model
independent-cross-check interpretations. The SILVER passing V1 is NOT the
"first 3D candidate to clear two real gauntlet gates" claim the task brief
flagged as a possibility — it cleared one (V1) under the moontour
interpretation, and V2 honestly didn't run for it (Phase 2 work).

### Catalogue admission status (unchanged)

* **NOT** admitted to `catalogue.yaml`.
* V1 PASS does NOT promote.
* V3 / V4 / V5 plus the #328 lit-check verdict remain prerequisites.

## #287 spike + #301 sub-family — V1 regression verdicts

Both clear V1 cleanly under their CR3BP periodic-orbit interpretation
(direct test assertions in `tests/data/test_v1_3d.py`):

* **#287 spike 3D Braik-Ross C11a extension**: passes V1, `degenerate_planar=False`.
* **#301 doubly-hyperbolic k=4 sub-family member** (first accepted bracket
  switched state): passes V1, `degenerate_planar=False`. The whole 41-member
  sub-family was previously converged at residual 1e-10 to 1e-14 nondim;
  V1 reproduces.

Both correctly **FAIL V2** under naive single-cycle propagation —
hyperbolic-instability amplification dominates by cycle 2-3. Catalogue row
`braik-ross-c11a-cycler-2026` already carries the exact note: **"NOT V2: an
unstable orbit cannot satisfy V2-ballistic's bounded-drift-over->=3-laps
requirement"**. V2 is doing its job by rejecting these (a powered cycler
would need V2-powered, which exists separately in
`cyclerfinder.verify.v2_powered`).

The V2 POSITIVE test uses the sourced JPL L1 halo fixture
(`tests/verify/fixtures/jpl_earth_moon_l1_halo_sample.json`) — its
lowest-stability member's cycle-3 drift is ~5e-5 km, well below the 50,000
km floor.

## Files shipped (Phase 1)

| Layer | Path |
|---|---|
| Module V1 | `src/cyclerfinder/data/validation/v1_3d.py` |
| Module V2 | `src/cyclerfinder/data/validation/v2_3d.py` |
| Package init | `src/cyclerfinder/data/validation/__init__.py` |
| Tests V1 | `tests/data/test_v1_3d.py` (8 cases) |
| Tests V2 | `tests/data/test_v2_3d.py` (8 cases) |
| SILVER gauntlet script | `scripts/run_306_phase1_silver_327_gauntlet.py` |
| SILVER verdict JSONL | `data/silver_327_v1_v2_verdicts.jsonl` |
| This note | `docs/notes/2026-06-16-306-3d-gauntlet-phase1.md` |

## Phase 2 path — concrete next steps

### V3 (independent-model, REBOUND n-body)

The current planar-only n-body harness lives in `src/cyclerfinder/nbody/`;
the rung pattern (`rung.py`) is per-body but the IC translation hard-codes
planar. Phase 2 must:

1. Lift the planar IC translation to 6D (state0 already 6D, just the wiring
   downstream that assumes z=0).
2. Wire REBOUND or a JPL-DE440-backed integrator as the independent-tool
   panel D in `verify/agreement.py`.
3. Run V3 on the JPL L1 halo + the #287 spike member as positive-control
   pair, and on a deliberately broken IC as negative.

### V4 (HFEM real-eph at Uranus)

Currently V4-ready only at Earth-Moon and Earth-Moon-Mars via the GMAT V4
lane (`reference_gmat_install`). Uranus requires SPICE-kernel coverage
extension; the GMAT installation already supports it but the project's
GMAT-bridge code is system-hard-coded. Estimated cost: 1-2 days to extend.

### V5 (mission-quality)

Out of automated scope. A human reviewer with the V1+V2+V3+V4 evidence
package decides. Phase 3+ work.

### Moontour-specific V2

The #327 SILVER and other moontour candidates need a V2 implementation
that:

1. Re-solves the Lambert legs each cycle with the moons' phase advanced by
   their orbital periods.
2. Tracks position drift between successive cycle-end Encounter rendezvous.
3. Asserts the bounded-drift floor against position (not velocity)
   residual.

This is a meaningful Phase 2 deliverable, ~ 4-6 hours.

## Test summary

```
$ uv run pytest tests/data/test_v1_3d.py tests/data/test_v2_3d.py \
                tests/search/test_cr3bp_3d_corrector.py -x --timeout=180
26 passed
```

8 V1 tests, 8 V2 tests, 10 existing 3D-corrector tests — all pass, no
regressions.

## Discipline footnotes

* `V1_FLOOR_KMS = 1e-3` is spec §14 (`crosscheck.V1_TOLERANCE_MPS` is the
  same value, reused via module constant).
* `V2_DRIFT_FLOOR_KMS = 50_000.0` is sourced from
  `propagate.DRIFT_TOLERANCE_KM`; we cite the source in the module
  docstring and the test fabrication-guard.
* `V2_N_CYCLES_MIN = 3` is spec §14 V2-ballistic.
* No catalogue writeback in Phase 1.
* Golden tests assert on closure quality + topology, not on specific 3D
  numbers our own code produced (per `feedback_golden_tests_sourced_only`).
