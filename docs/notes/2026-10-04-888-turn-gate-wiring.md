# `#888`: wiring the demanded-turn gate into the validation lane

Date: 2026-10-04. Follows `docs/notes/2026-10-04-888-demanded-turn-gate.md` (section 6 listed the
places to wire). Commits `214315ec` (helper) and `f91c957a` (wiring and tests). The gate modules
(`verify/turn_gate.py`, `verify/turn_gate_closures.py`) are unchanged.

## 1. What was wired where

The glue is a new module, `src/cyclerfinder/data/validation/moontour_turn.py`. Each tier's
`_cycle_*` function takes an optional `CycleTurnRecord` and, when one is passed, stores every leg's
departure and arrival V-infinity vectors (Lambert velocity minus the moon's velocity, in the tier's
own frame) and the departure of the next cycle's first leg. `chain_encounters` builds one encounter
per intermediate flyby (`"e1"`, ...) plus the anchor wrap (`"wrap"`: the last leg's arrival against
the next cycle's first departure, at the same epoch and the same moon state, so no rotation is
needed). The gate is `demanded_turn_gate` with registry floors, and the verdict is `turn_feasible`,
not `ballistic` (the magnitude continuity stays with each tier's own closure residual). A
missing wrap leg or a zero V-infinity is recorded as an error and counts as not feasible; it never
raises. Callers that pass no record see no change.

| file, function | change |
|---|---|
| `data/validation/v2_moontour.py`, `_cycle_residual` / `run_v2_moontour` | the cycle keeps the vectors and solves the wrap leg (same branch rule); `passes_v2` requires `turn_feasible`; new `turn_alt_floor_km` keyword |
| `data/validation/v3_3d.py`, `run_v3_3d` | `passes_v3` requires `v2_verdict.turn_feasible` (V3 re-propagates the same legs leg by leg, so integrator agreement says nothing about the flybys). Not on the original list, added because the six rows passed V3 too |
| `data/validation/v4_uranus.py`, `_cycle_v4` / `run_v4_uranus` | as V2; the wrap leg is a Lambert solve only |
| `data/validation/v4_uranus_strict.py`, `_select_leg_transfer` / `_cycle_v4_strict` / `run_v4_uranus_strict` | `_LegOutcome` carries the chosen branch's Lambert `v1`/`v2`; the cycle stores them against the SPICE moon velocities; the run loop takes each cycle's wrap from the next cycle's first leg and selects one extra leg after the final cycle (no extra propagation for the other cycles); the per-cycle verdicts are rebuilt with the turn fields |
| `data/validation/v2_saturn_3d.py`, `v3_saturn_3d.py`, `v4_saturn.py` (`_cycle_v4`), `v4_saturn_strict.py` (`_select_leg_transfer`, `_cycle_v4_strict`) | the same, for the Titan-Iapetus copies (Titan wrap at the 1500 km registry floor). `v2_saturn_3d` and `v3_saturn_3d` were not on the list; they are the same lane |
| `search/correct.py`, `_bend_feasible` | new keyword `wrap_rotation_rad` adds the periodicity wrap (`_wrap_bend_feasible`, the `turn_ratio_check.wrap_node_turn` construction); `ballistic_correct` always reports `wrap_bend_feasible` and folds it into `bend_feasible` only with `gate_wrap=True` |
| `scripts/scan_558_uranus_all_pairs_offset_sweep.py`, `gate_candidate` | `all_gates_passed` now requires the turn gate on the rebuilt closure (`symmetric_closure`) instead of the `#324` capacity gate. `physical_gate_passed` (capacity) is still reported because `refine_562` and `compare_576` read it |

Every verdict object now carries `turn_feasible`, `worst_turn_ratio` and `turn_failure_reason`
(the binding flyby: body, label, cycle, demanded and available angle, floor, ratio, required
altitude); per-cycle verdicts carry `turn_feasible` and `turn_encounters` (the gate's
`EncounterTurn` rows). The strict JSON writers add the same fields. All new fields have defaults
(`turn_feasible=False` when not evaluated), so existing constructors keep working.

Deviation, deliberate: in `ballistic_correct` the wrap is opt-in. Making it the default would change
`bend_feasible` on every heliocentric row whose published turn ratio is reproduced on the
intermediate nodes only, and gating the heliocentric rows is plan item 2(b), which was not run.

Powered V2 (`releg`): the gate is applied to the ballistic Lambert legs, because the powered close
pins V-infinity magnitudes only and does not expose its leg vectors. This is conservative (a DSM
could fix a direction) and is stated on the verdict field.

## 2. What now fails that used to pass (all measured in the test runs)

* The six withdrawn Uranian rows, fed exactly as the `#566` gauntlet and the `#330` run fed them:
  V2 (with drift and closure floors opened, so the turn is the only criterion left), V3, V4-scipy
  and V4-strict all reject, with the turn named. The rebuilt angles match section 4 of the gate note
  to the printed digits (for example Titania-Oberon 152.6 / 7.0 and 177.2 / 6.3 degrees).
* V4-strict at the canonical epoch 2000-06-21: Umbriel-Titania no longer reaches the turn gate,
  because since `#567` every branch of its first leg is planet-crossing there; the stored `#566`
  PASS predates `#567`. One day later all cycles converge and the turn rejects it (Umbriel wrap,
  34.8 against 7.7 degrees).
* All 17 of the 22 `#574` Stage-A Titan-Iapetus branches that complete three V2 cycles are
  turn-infeasible (worst ratios 7.9 to 44.6, from a run of all branches); `#574` Stage B had
  already passed 0 of 15, so no claim changes.
* `search.correct`: an Earth-to-Earth `nPr` chain has no intermediate flyby, so `_bend_feasible`
  without the wrap accepts every McConaghy cycler; with the wrap it accepts exactly 6S7, 6S8, 6S9
  and rejects 1L1, as Table 4 footnote e says.

Acceptance (so the wiring does not reject everything): the chain builder accepts all eight
Russell-Strange generic-leg cyclers that clear the project floors and reproduces their published
minimum flyby altitude within 4 km (GanIo#403: 20 km). The two-moon lane cannot ingest those cyclers
(their target moon is massless and unphased), so the lane-level control is the `#890`
Titania-Oberon chain, a project computation: V2 passes it (ratio 0.57) under the same open drift
floor that the six rows fail under. At the default 50,000 km floor it fails on drift (the V2 drift
is an inertial offset of the closing encounter, about 7e5 km here), not on the turn. A wrap-only
case, Titania-Umbriel-Titania (legs of 3 synodic periods, 1 revolution each, also a project
computation), passes every pre-`#888` criterion once the drift floor is opened (magnitude
closure at its default floor) and has a feasible Umbriel flyby (ratio 0.3); it is now rejected at
the Titania wrap (ratio 22.8). A check of the intermediate encounters alone
would still pass it.

## 3. Tests

New: `tests/verify/test_888_turn_gate_wiring.py` (41 tests, all pass, URA111 kernel present). The
first run, before any wiring, failed 32 of them (log kept outside the repo).

Changed, none deleted:

* `tests/data/test_v2_moontour.py::test_drift_floor_correctly_passes_small_drift`: asserted
  `passes_v2 is True` with opened floors; now asserts the floors accept and `passes_v2` is False on
  the turn.
* `tests/data/test_v3_3d.py::test_silver_v3_ias15_passes_agreement_floor`: agreement assertions
  kept; `passes_v3` is now asserted False on the inherited turn verdict.
* `tests/data/test_v4_uranus_strict.py::test_560_silver_312_canonical_epoch_unchanged`:
  convergence, branch-selection and drift-agreement assertions kept (they still guard the `#560`
  fixes); the headline is now `passes_v4_strict is False`, on the turn.
* Frozen-gate tests `tests/verify/test_566_five_representatives_v4.py` and
  `tests/verify/test_silver_327_{v1_passes,v2_quasi_cycler,v3_nbody,v4_strict_passes}.py`: the
  stored-data checks are untouched, each docstring now says the stored verdicts come from an
  insufficient gauntlet and back no claim, and each gains a test that the gated lane rejects the
  same inputs on the turn.
* `tests/data/test_saturn_v2v3v4_gauntlet.py`: a new class checking the Saturn tiers' flyby
  structure, registry floors and the pass identity.

## 4. What is still leg by leg

Gating the turn makes the verdicts honest about the flybys, but no tier flies them. V2 and V3 solve
or re-propagate each Lambert leg separately. V4-scipy and V4-strict start each leg from the moon's
position with that leg's own Lambert departure velocity, and the next leg starts again from the
next moon, whatever velocity and position the previous leg arrived with (the arrival miss is
reported as `v4_terminal_offset_vs_moon_kms`). The V4 drift is therefore a sequence of independent
legs, not a trajectory. The turn gate uses the Lambert-conic vectors of the legs the tier targeted,
not the propagated arrival velocity.

What a continuous propagation would need (not built):

1. The tour moons as unsoftened point masses. V4 zeroes or softens each moon's third-body term
   inside its Hill sphere to avoid the patched-conic singularity, so even an unbroken propagation
   in the current force model cannot bend at a moon.
2. Step control through periapsis (a flyby at about 50 km above a 500 to 800 km moon at under
   1 km/s).
3. Per-encounter targeting: B-plane or periapsis-state multiple shooting so that each leg's
   arrival hyperbola leaves on the next leg, with full-state continuity across encounters and
   across cycles (the wrap), and a powered correction only where the gate says the turn is short.
4. For the real-ephemeris tier, the same with SPICE moon states, and a check that the converged
   flyby altitudes stay above the floors.

`verify/flyby_integrate.py` (one flyby, from several Hill radii in to several out) is the natural
building block for item 3.

## 5. Not wired (still capacity-only or magnitude-only)

* `genome/titan_iapetus_corrector.py` and the scripts that call `candidate_passes_physical_gate`
  directly (`enumerate_600` 3-moon gate `gate_candidate_3moon`, `scan_816`, `run_573`, `run_574`,
  `probe_575*`, `probe_655*`, `verify_576`, the older `scan_3xx`/`scan_433` scripts and the
  empty-region stamping scripts) still use the capacity gate. Only `scan_558.gate_candidate` and
  its generic callers (`enumerate_563`, `scan_571`, `refine_562`, `compare_576`, `_apply_609`,
  `certify_610`) now see the turn. No script that writes `data/` was re-run.
* The bounded-drift labels live in scripts, not in `src`. `FAIL_QUASI_BOUNDED` is stamped on a
  strict-V2 FAIL with a bounded drift shape, so gating `passes_v2` alone would NOT have removed it:
  a turn-infeasible chain still fails V2 and still qualified. Its conditions in
  `scripts/run_330_silver_moontour_v2.py`, `run_566_gauntlet_five_representatives.py` and
  `run_574_stageB_saturn_gauntlet.py` now also require `turn_feasible` (edited, not re-run).
  `PASS_AS_QUASI_CYCLER` additionally needs V3, V4 and V4-strict passes, which now require the
  turn. `EFFECTIVELY_CYCLIC` (`analyze_338_boundary.py`) reads stored epoch-sweep output; a new
  sweep would inherit the gate through `passes_v4_strict`. None of these scripts was re-run.
* `search.correct.ballistic_correct` gates the wrap only with `gate_wrap=True` (section 1). Its
  ephemeris-driven wrap is tested against `turn_ratio_check.wrap_node_turn` on
  `russell-ch4-9.353Gg2` (`tests/search/test_turn_ratio_check.py`).
