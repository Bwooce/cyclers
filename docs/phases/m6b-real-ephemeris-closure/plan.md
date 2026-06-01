# M6b — Real-ephemeris closure verification

**Spec reference:** spec.md §4 (architecture — `verify/propagate.py` extends from M6a; `search/phase_match.py` already shipped via the M6 slice `9b2611d`; `search/optimize.py::optimise_cell_ephemeris` is the locked stub M6b fills), §5 step 6 (refine on ephemeris) + step 7 (verify — extended from M6a's idealized verification to real DE440 closure over N cycles), §6 (top-level `find_cyclers` interfaces unchanged at M6b; M6b adds `verify_real_closure` as a thin wrapper around the extended pipeline), §8 (M6 milestone, real-ephemeris half), §9.1 (Aldrin and 2-synodic E-M anchors that still bind in real ephemeris), §10 (closure-frame correctness risk now applied across multi-cycle horizons; multi-rev Lambert risk reactivated as a hard blocker), §12(a) (idealised vs ephemeris-mode optimisation — **M6b is the ephemeris-mode half**; M6a wrote the verification machinery, M6b uses it), §12(c) (the dynamic ephemeris frame M6a built — M6b reuses it without changes), **§12.1 (idealised→ephemeris bridge — `phase_match.find_real_windows` is the entry point M6b composes; the M6 slice already shipped it)**, §12.2 (the three-representation framework — **M6b promotes a catalogue entry from "idealised form" V1 to "real-ephemeris instance" V2** but does not flatten one representation into another), §13.5 (signature-based dedup — M6b's outputs share signatures with M6a's), §14 (V0–V5 validation gauntlet — **M6b stands up the V2 → V3 promotion gate**: V2 was multi-lap *idealised* periodicity; V3 is real-ephemeris horizon TCM bounded), §16.1 (catalogue schema v2; `validation.gates.V3.horizon_tcm_mps` field exists already — M6b populates it), §16.2 (canonical signatures — M6b's real-ephemeris instances inherit their *idealised parent's* signature, do not generate their own).

**Purpose:** stand up the **real-ephemeris closure verification gate** — the gate that promotes a catalogue entry from V1 (reproduced in circular-coplanar; M6a's binding state) to V2-real (closes ≥ 1 full cycle within drift tolerance on JPL DE440; **M6b's contribution**). M6a's `verify_long_term_stability` already handles the multi-lap drift measurement against a passed-in `Ephemeris`; what M6b adds is (1) the construction path that produces a *real-ephemeris instance* `Cycler` from an idealised V1 entry, (2) the multi-cycle drift tolerance and rationale, (3) the catalogue loader that turns the V1-passed catalogue into M6b inputs, (4) the regression set of literature-anchored real-ephemeris closures and the EXPECTED_SKIPS registry of entries known not to close, and (5) the V3 hand-off shape that fills `validation.gates.V3.horizon_tcm_mps`. **TCM ΔV budgeting itself is M7's concern**; M6b measures positional drift and reports it. The boundary between "real-ephemeris closure (geometric)" and "TCM budget (operational)" is one of the load-bearing decisions of this plan — see §3.5.

Pascarella et al. 2024 (Solar System Pony Express) is the architectural template per `docs/v2-future-references.md`: their patched-conic → medium-fidelity → high-fidelity pipeline directly maps onto M6a (medium-fidelity stage) and M6b (the *output of medium-fidelity, before high-fidelity GMAT*). M6b stays inside v1's patched-conic + impulsive-flyby scope; their N-body stage is M7's optional V4 gate.

**Gate (definition of done):**

1. `tests/verify/test_real_closure.py::test_aldrin_cycler_periodic_over_2_cycles_astropy` (**M6b BINDING GATE — spec §8 M6 anchor, real-ephemeris extension**) asserts the Aldrin classic E-M k=1 outbound cycler, loaded from catalogue entry `aldrin-classic-em-k1-outbound`, **constructed at a phase-matched real launch date via `find_real_windows`**, fed to `verify_real_closure(cycler, n_cycles=2, ephem=Ephemeris("astropy"))` returns `result.closes == True` with `result.max_drift_km < REAL_DRIFT_TOLERANCE_KM` (200,000 km per §4.3) and `result.n_cycles_propagated == 2`. **The cycler is the fixture, not the tolerance** — no cherry-picking, no per-test widening. The Aldrin gate is the M6b binding choice (vs the catalogue's 2-syn `s1l1-2syn-em-cpom` entry) because (a) Aldrin is the canonical free-return reference, (b) the 2-syn S1L1 construction requires multi-rev Lambert (an outstanding hard blocker), and (c) Aldrin closes at k=1 synodic ≈ 2.135 yr so a 2-cycle verification covers ≈ 4.3 yr — long enough for real eccentricity to bite but short enough to keep CI runtime under the 60-second budget per §6.
2. `tests/verify/test_real_closure.py::test_2syn_em_cpom_periodic_over_2_cycles_astropy` (**M6b ASPIRATIONAL GATE; xfail at M6b time pending multi-rev Lambert**) asserts the 2-synodic E-M S1L1 cycler from `s1l1-2syn-em-cpom` closes over 2 cycles on astropy. Marked `pytest.mark.xfail(strict=False, reason="multi-rev Lambert blocker — see §5 risk #2")`; **flips to passing once multi-rev Lambert lands** (M-future or stretch). Documented in the EXPECTED_SKIPS registry per §4.4.
3. `tests/verify/test_real_closure.py::test_real_drift_rejects_open_trajectory` asserts that a deliberately-perturbed Aldrin cycler (whose lap-0 V∞ has been rotated by 5° at one encounter, breaking ballistic continuity) returns `result.closes == False` with `result.max_drift_km > 5 × REAL_DRIFT_TOLERANCE_KM` over 2 cycles — i.e. the gate has rejection power, not just acceptance.
4. `tests/verify/test_real_closure.py::test_real_closure_regression_set` runs the M6b regression set (§4.2: 5 catalogue entries known to close in literature on real ephemeris) and asserts every entry not in `EXPECTED_SKIPS` returns `closes == True`. Skipped entries are documented per §4.4.
5. `tests/verify/test_real_closure.py::test_real_closure_result_frozen_and_v3_fields_locked` asserts `RealClosureResult` raises `FrozenInstanceError` on attribute assignment; that `horizon_tcm_mps == 0.0` and `per_cycle_tcm_mps == (0.0,) * n_cycles` (M6b populates **only positional drift**, not TCMs — those are M7); and that `result.v3_status in ("v3-real-closure-pass", "v3-real-closure-fail", "v3-skipped-multirev")`.
6. `tests/verify/test_real_closure.py::test_real_closure_uses_m6a_machinery` asserts that `verify_real_closure` internally invokes `multi_lap_propagation` + `lap_to_lap_drift` from M6a; it does not duplicate the propagator. This is a binding **composition** assertion — the failure mode it catches is "M6b reimplemented propagation and the two implementations silently diverge."
7. `tests/search/test_phase_match.py::test_find_real_windows_for_aldrin_signature_within_priority_window` asserts that `find_real_windows` returns at least one window within ±5 years of the catalogue's `priority_date` for the Aldrin entry. This is the *plumbing* test for M6b's construction path: if `find_real_windows` returns no window for a literature-anchored cycler, M6b cannot construct anything to verify.
8. `tests/data/test_catalogue_loader_m6b.py::test_loader_filters_v1_pass_circular_coplanar_ballistic_sun_only` asserts the catalogue loader (§3.2) returns exactly the entries whose `model_assumption == "circular-coplanar"` AND `trajectory_regime == "ballistic"` AND `primary == "Sun"` (or absent, which defaults to Sun per the schema-v2 backfill rules) AND who passed V1 in the prior catalogue run. The count should be ~213 (per the catalogue audit) minus the 4 `cr3bp` entries minus the 2 `analytic-ephemeris` entries minus any non-Sun entries, less the V1-fail subset.
9. `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests` all clean. CI green on the M6b commit.

The `REAL_DRIFT_TOLERANCE_KM = 200,000.0` is M6b's **third binding architectural choice** (after M6a's dynamic-frame and 50,000 km idealised tolerance) and is justified in §4.3 below. It is 4× M6a's tolerance because (a) M6b verifies over N=2 cycles rather than M6a's 3 laps, doubling the integrated breathing amplitude, (b) the Aldrin k=1 cycler closes only every 2.135 yr (vs the 2-syn cycler's 4.27 yr) so 2 cycles ≈ 4.3 yr matches M6a's 3-lap horizon, and (c) the real-ephemeris case absorbs ~3% Earth eccentricity + ~9% Mars eccentricity in the geometric breathing, vs the idealised case's near-zero breathing in M6a's circular fallback. See §4.3 for the derivation.

---

## 0. Naming and scope conventions used in this plan

To keep the M6a and M6b machinery distinguishable in code and prose:

- **"Lap"** is M6a's concept: one revolution of the cycler's encoded geometry, propagated from leg start to leg end without re-anchoring to the real planet. M6a measures lap-to-lap drift in the dynamic rotating frame.
- **"Cycle"** is M6b's concept: one full cycler period propagated on **real** DE440 ephemeris, with planet positions read at the actual cycle epoch. A cycle and a lap are the same length of time (≈ one cycler period); the difference is which ephemeris drives the planet states.
- **"Real-ephemeris closure"** is the M6b gate: the spacecraft state at the end of cycle N matches the spacecraft state at the start of cycle 0 to within `REAL_DRIFT_TOLERANCE_KM` when both are projected into the dynamic rotating frame.
- **"V1"** (per spec §14) is the multi-lap idealised periodicity gate — passed by `verify_long_term_stability` returning `stable=True` on `Ephemeris("circular")`. **M6a's binding gate is V1.**
- **"V2"** is multi-lap *bounded drift on real ephemeris* — passed by `verify_real_closure` returning `closes=True` on `Ephemeris("astropy")`. **M6b's binding gate is V2.**
- **"V3"** is real-ephemeris **horizon-TCM** bounded — TCMs computed and summed against a published budget. **M7 stands up V3**; M6b ships the V3 placeholder fields on `RealClosureResult` so M7 only populates them rather than reshaping the dataclass.

This nomenclature aligns with spec §14 (gauntlet V0–V5) and §12.2 (the idealised vs real-ephemeris representation distinction). The "V2-real" shorthand is used where ambiguity matters.

---

## 1. What this milestone delivers

Four new pieces of source code and their test files; one new module (`verify/real_closure.py`); two extensions to existing modules. M6b is **additive** for everything except `search/optimize.py::optimise_cell_ephemeris` (whose locked stub from M5 may be partly filled — see §3.4 for the policy).

### 1.1 `src/cyclerfinder/verify/real_closure.py` (NEW)

The headline M6b module. Public surface, in dependency order:

- `REAL_DRIFT_TOLERANCE_KM: Final[float] = 200_000.0` — the binding tolerance per §4.3.
- `N_CYCLES_DEFAULT: Final[int] = 2` — the default cycle count for the binding gate. Documented as the trade-off between sensitivity (more cycles = more accumulated drift visible) and CI budget (more cycles = longer wall-clock).
- `RealClosureResult` — frozen dataclass per §3.3 below; the locked spec §14 V2-real verification interface.
- `verify_real_closure(cycler, n_cycles, ephem, *, t_start=None, frame_bodies=None, cycler_id=None, signature_priority_date=None) -> RealClosureResult` — the public M6b entry point. Wraps `verify_long_term_stability` with:
  - real-ephemeris epoch defaulting via `find_real_windows` if `t_start` is `None`,
  - drift tolerance comparison against `REAL_DRIFT_TOLERANCE_KM`,
  - V3 placeholder fields (`horizon_tcm_mps=0.0`, `per_cycle_tcm_mps=()`, `v3_status="v3-real-closure-pass" | "v3-real-closure-fail" | "v3-skipped-multirev"`) populated for M7 hand-off.
- `construct_real_ephemeris_cycler(catalogue_entry, ephem, launch_window) -> Cycler` — given a catalogue YAML entry + a `LaunchWindow` from `find_real_windows`, build a `Cycler` whose `encounters[].r` are the real planet positions at the launched epoch and whose `legs[]` are Lambert-solved across those real positions.
- `EXPECTED_SKIPS: dict[str, str]` — module-level registry of `catalogue_id → reason` for entries known not to close on real ephemeris (e.g. multi-rev-Lambert blocker, near-grazing flybys outside Kepler tolerance). Documented per §4.4.

Module-internal helpers (private):

- `_resolve_real_t_start(signature, ephem, signature_priority_date) -> float | None` — pick the best `find_real_windows` epoch closest to `signature_priority_date`; return `None` if no window is below the M6b mismatch cap (the open-trajectory case).
- `_construct_cycler_from_lambert_chain(catalogue_entry, ephem, t_start) -> Cycler` — Lambert-solve each leg in order using real planet positions. Single-rev only at M6b; raises `MultiRevLambertRequiredError` (a new exception subclass) when a leg's published `n_revs > 0`. The catch-and-skip pathway in `verify_real_closure` uses this exception to populate `v3_status="v3-skipped-multirev"`.
- `_check_vinf_continuity(cycler, ephem, tolerance_kms=0.5) -> tuple[float, ...]` — at each interior encounter, compute the |V∞_in| vs |V∞_out| mismatch. Pure ballistic flybys preserve |V∞|; M6b reports the per-encounter mismatch as a *diagnostic*, not a gate (drift is the gate). Used in `RealClosureResult.per_encounter_vinf_mismatch_kms`.

### 1.2 `src/cyclerfinder/verify/__init__.py` (EXTEND — additive only)

Re-export the new public surface. After M6b:

```python
from cyclerfinder.verify.propagate import (  # M6a (unchanged)
    DRIFT_TOLERANCE_KM,
    StabilityReport,
    lap_to_lap_drift,
    multi_lap_propagation,
    propagate_lap,
    verify_long_term_stability,
)
from cyclerfinder.verify.real_closure import (  # M6b (new)
    EXPECTED_SKIPS,
    N_CYCLES_DEFAULT,
    REAL_DRIFT_TOLERANCE_KM,
    RealClosureResult,
    construct_real_ephemeris_cycler,
    verify_real_closure,
)
```

### 1.3 `tests/data/_catalogue_loader_m6b.py` (NEW test helper)

A test-only catalogue loader filtered to M6b scope. Public surface:

- `load_m6b_entries() -> list[dict]` — read `data/catalogue.yaml`, filter to `model_assumption == "circular-coplanar"` AND `trajectory_regime == "ballistic"` AND `primary in (None, "Sun")` AND the entry passed V1 in the M6a regression record. Returns the raw dict per entry (caller turns it into a `Cycler` via `construct_real_ephemeris_cycler`).
- `M6B_REGRESSION_IDS: tuple[str, ...]` — the 5-entry regression set documented in §4.2.

Lives under `tests/data/` (NOT `src/`) because the v1 spec defers a full `data/catalog.py` module to M7; M6b's loader is test infrastructure that M7's loader supersedes.

### 1.4 `tests/verify/test_real_closure.py` (NEW)

The M6b gate tests + helper tests live here. Structure per §4:

- Gate 1: `test_aldrin_cycler_periodic_over_2_cycles_astropy` (M6b binding).
- Gate 2: `test_2syn_em_cpom_periodic_over_2_cycles_astropy` (xfail under multi-rev Lambert blocker).
- Gate 3: `test_real_drift_rejects_open_trajectory`.
- Gate 4: `test_real_closure_regression_set`.
- Gate 5: `test_real_closure_result_frozen_and_v3_fields_locked`.
- Gate 6: `test_real_closure_uses_m6a_machinery`.
- Helper tests: `test_construct_real_ephemeris_cycler_aldrin`, `test_construct_raises_on_multi_rev_leg`, `test_resolve_real_t_start_prefers_priority_window`, `test_check_vinf_continuity_diagnostic`.

### 1.5 `tests/data/test_catalogue_loader_m6b.py` (NEW)

Loader filter tests per gate #8.

### 1.6 `tests/search/test_phase_match.py` (EXTEND)

Add `test_find_real_windows_for_aldrin_signature_within_priority_window` per gate #7. The existing M6-slice tests remain unchanged; M6b extends.

### 1.7 `src/cyclerfinder/search/optimize.py::optimise_cell_ephemeris` (STUB POLICY — leave RAISING)

Spec §12(a) calls for an ephemeris-mode optimiser that minimises summed TCM ΔV. **M6b explicitly does NOT fill this in.** The body of `optimise_cell_ephemeris` stays raising `NotImplementedError("requires M6 ephemeris backend (multi-lap propagator + TCM budget). Stub locked in M5; body lands in M6b.")` because:

- The TCM budget is an **M7 concern** (spec §16.1 catalogue field `validation.gates.V3.horizon_tcm_mps`). M6b's job is the geometric closure gate (V2-real), not the operational maintenance budget (V3).
- The ephemeris-mode optimiser is a *consumer* of `verify_real_closure` (it would call M6b's drift measurement inside its objective), so wiring it requires M6b's surface to exist first.
- Filling the body without TCM-budget tests creates a half-implemented optimiser that future M7 work has to rewrite. Better to leave the stub explicit until M7's plan defines the TCM budget shape.

The error message is **updated** to reflect M6b's new state: `"requires M6b real-ephemeris closure (shipped) + M7 TCM budget machinery (not yet shipped). M6b's verify_real_closure is the right drift-feasibility check; wiring lands in M7."` This is the only edit to `optimize.py` in M6b.

### 1.8 Explicit non-goals (M6b boundaries)

These belong to M7 or later; **do not stub or partially implement them in M6b** beyond what the locked dataclass shape requires:

| Out of M6b | Where it lands | Why deferred |
|---|---|---|
| TCM ΔV budget computation (per-cycle ΔV, summed horizon TCM) | **M7** | M6b populates `horizon_tcm_mps=0.0` and `per_cycle_tcm_mps=()` as locked placeholders. M7 fills the bodies. Computing TCMs requires the M7 catalogue ingest writer that consumes `RealClosureResult` and writes back into `validation.gates.V3.horizon_tcm_mps` — that flow doesn't exist yet. |
| `optimise_cell_ephemeris` body fill | **M7** | Per §1.7 above. M6b updates the error message; body stays raising. |
| Multi-revolution Lambert solver | **M-future (stretch) or M7 prerequisite** | The `lambert(..., max_revs=N>0)` body is a known blocker per spec §10 and the M1 ports. M6b documents the blocker; entries needing multi-rev Lambert go into `EXPECTED_SKIPS`. |
| `verify/crosscheck.py` (V1 Lambert cross-check on real ephemeris) | **M7** | Spec §14 V1 is *Lambert-vs-Lamberthub cross-check*; that's an M7 gate, structurally independent of M6b's V2-real propagation. |
| GMAT bridge (V4) | **Stretch** | Spec §14 V4, explicit stretch goal per spec §7. |
| Catalogue YAML write-back (`validation.gates.V2.*` field population) | **M7** | M7's batch-validate runner consumes `RealClosureResult` and writes the YAML; M6b only produces the result. |
| `find_real_windows` for non-Sun primaries | **M-future** | The M6-slice loader raises `NotImplementedError` for `primary != "Sun"`. M6b inherits that scope; lunar/Jovian/Saturnian cyclers (4 entries) sit in `EXPECTED_SKIPS` with reason `"non-Sun primary"`. |
| Real-ephemeris instance writer (`windows.json` derived view) | **M7** | spec §12.2 / overview note that `cyclers.space/src/data/windows.json` is the public-facing real-ephemeris instance store. M7 fills it from `RealClosureResult` outputs. M6b does not write any data files. |
| 3+ cycle horizons (per spec §12(a) "3–5 lap horizon") | **M7** | M6b defaults to `n_cycles=2` to stay inside CI budget. 5-cycle runs are M7 batch concern; the dataclass supports `n_cycles >= 2` so M7 can extend without reshape. |
| VEM (3+ body) real-closure | **M8** | The VEM entries (`jones-2017-vem-triple-family`, `vem-emeeve-3syn`) go into `EXPECTED_SKIPS` with reason `"VEM real-closure is M8 scope"`. M6b's frame helper from M6a (`_resolve_frame_bodies`) handles VEM in principle, but the catalogue's VEM entries' `priority_date` may not yield a `find_real_windows` hit and the multi-rev Lambert constraint hits harder for 3-body chains. |

---

## 2. File tree after M6b

```
cyclers/
├── … (M0/M1/M2/M3/M4/M5/M6-slice/M6a layout preserved unchanged)
├── src/cyclerfinder/
│   ├── core/                              # M0–M6a — unchanged
│   ├── search/
│   │   ├── … (M3/M4/M5/M6-slice — unchanged)
│   │   └── optimize.py                    # M5 + M6b (one-line error message edit)
│   ├── model/                             # M3 — unchanged
│   └── verify/
│       ├── __init__.py                    # M6a + M6b additive re-exports
│       ├── propagate.py                   # M6a — UNCHANGED
│       └── real_closure.py                # NEW (M6b)
└── tests/
    ├── core/                              # M0–M6a — unchanged
    ├── search/
    │   └── test_phase_match.py            # M6-slice + M6b (one new test added)
    ├── verify/
    │   ├── test_propagate.py              # M6a — unchanged; xfail astropy test stays xfail
    │   └── test_real_closure.py           # NEW (M6b — gate + regression + helper tests)
    └── data/                              # NEW test sub-package
        ├── __init__.py                    # NEW
        ├── _catalogue_loader_m6b.py       # NEW
        └── test_catalogue_loader_m6b.py   # NEW
```

Subpackages `data/` (under `src/`) and `viz/` still remain uncreated — M7/M8 territory. The only edits to pre-existing source files are (1) `verify/__init__.py` adding three re-exports, (2) `search/optimize.py::optimise_cell_ephemeris` updating its error-message string. No edits to `verify/propagate.py`, `core/frames.py`, `core/ephemeris.py`, `core/kepler.py`, `model/cycler.py`, or any other M0–M6a module body.

---

## 3. Module designs

This section walks the four pieces of M6b in dependency order: (3.1) the construction path that turns a catalogue entry into a real-ephemeris `Cycler`; (3.2) the catalogue loader and its V1-pass filter; (3.3) the `RealClosureResult` dataclass + `verify_real_closure` API; (3.4) the relationship to `optimise_cell_ephemeris`; (3.5) the explicit TCM-budget deferral; (3.6) the imports / dependency graph.

### 3.1 Construction path — primary: Lambert-chain across real planet positions

**Binding architectural decision.** There are two candidate construction paths per the M6a hand-off summary:

- **Path A — Real-ephemeris optimiser pass via `optimise_cell_ephemeris`.** M5 left the stub raising; the natural shape is "given a cell + an ephemeris, search free-parameter timing to minimise drift over N cycles." This is the spec §12(a) ephemeris-mode optimiser. **Currently blocked**: the M5 idealised optimiser's binding gate (`test_find_cyclers_em_top_level`, task #54 per the M6a hand-off) is broken, so calling `optimise_cell_ephemeris` even with a body would consume a broken upstream optimiser state.
- **Path B — Phase-matched single-window Lambert chain across encounters.** Use `find_real_windows` to pick a launch date that matches the catalogue's signature; then for each leg, Lambert-solve over real planet positions at the matched epoch + cumulative leg ToF. This is what `tests/verify/test_propagate.py::test_2syn_em_cycler_periodic_over_3_laps_astropy` (currently xfailed) attempts.

**M6b's primary path is B (Lambert-chain).** Reasons in priority order:

1. **Path A's broken precondition.** Wiring `optimise_cell_ephemeris` requires fixing the M5 binding gate first — work that's neither M6a's nor M6b's scope. Building M6b on a known-broken M5 path would couple our progress to a separate fix.
2. **Path B's simpler test surface.** A Lambert chain is deterministic given (catalogue_entry, launch_epoch); it has no optimiser-internal seed, no convergence-failure mode beyond Lambert's own (already handled in `find_real_windows`), and produces a `Cycler` that's directly consumable by `verify_long_term_stability`. The same code path can be unit-tested for the Aldrin case without invoking any optimisation.
3. **Pascarella 2024 alignment.** The Pascarella patched-conic → medium-fidelity pipeline uses Lambert chaining for the patched-conic stage and only invokes optimisation at the medium-fidelity-impulsive stage. M6b is the *output of medium-fidelity*; Lambert-chain construction matches the architectural template.
4. **Path A as future enhancement.** Once M5's binding gate is fixed (task #54) and M7 defines the TCM-budget shape, `optimise_cell_ephemeris` can be wired by passing `verify_real_closure` as its drift-feasibility check. The Lambert-chain construction then becomes the initial-guess generator for the optimiser. This is documented as M7 follow-up.

**Path B's hard dependency: multi-rev Lambert.** The catalogue's 2-syn S1L1 entry has legs published with `n_revs > 0` (specifically `n_revs=1` per the catalogue YAML for the E-E and M-M intermediate legs). The M1 Lambert solver raises `LambertGeometryError` for `n_revs > 0` because the body is not yet implemented (per `core/lambert.py:252`). M6b's response:

- The Lambert-chain constructor catches `LambertGeometryError` on any leg with `n_revs > 0` and raises a new `MultiRevLambertRequiredError(catalogue_id, leg_index)`.
- `verify_real_closure` catches `MultiRevLambertRequiredError` and returns a `RealClosureResult` with `closes=False`, `v3_status="v3-skipped-multirev"`, `max_drift_km=float("inf")`.
- The catalogue entry is added to `EXPECTED_SKIPS` so the regression test marks it as known-skip (not failure).

This keeps M6b honest: a multi-rev-needing entry is *not silently treated as failing*; it's surfaced with a precise reason that M7 (or M-future-multi-rev) can act on.

#### 3.1.1 Algorithm — `construct_real_ephemeris_cycler`

```text
INPUT: catalogue_entry (dict from YAML), ephem (Ephemeris("astropy")),
       launch_window (LaunchWindow from find_real_windows)

ALGORITHM:
  1. bodies = tuple(catalogue_entry["bodies"])     # e.g. ("E", "M", "E", "M") for 2-syn
  2. legs_meta = catalogue_entry["legs"]            # list of {tof_days, n_revs, branch, ...}
  3. t_start_sec = _dt_to_t_sec(launch_window.departure_date)

  4. # Construct encounter epochs
     encounter_times = [t_start_sec]
     for leg in legs_meta:
       encounter_times.append(encounter_times[-1] + leg["tof_days"] * SECONDS_PER_DAY)

  5. # For each leg, Lambert-solve with REAL planet positions at the encounter epochs
     legs = []
     encounters = []
     for j, leg in enumerate(legs_meta):
       body_dep = bodies[j]
       body_arr = bodies[j + 1]
       t_dep    = encounter_times[j]
       t_arr    = encounter_times[j + 1]

       r1, v1_planet = ephem.state(body_dep, t_dep)
       r2, v2_planet = ephem.state(body_arr, t_arr)

       if leg.get("n_revs", 0) > 0:
         raise MultiRevLambertRequiredError(
           catalogue_entry["id"], leg_index=j
         )

       try:
         sols = lambert(r1, r2, t_arr - t_dep, max_revs=0)
       except (LambertConvergenceError, LambertGeometryError) as e:
         raise RealClosureConstructionError(
           catalogue_entry["id"], leg_index=j, cause=e
         )

       sol = sols[0]                                 # single-rev only
       vinf_out_dep = sol.v1 - v1_planet            # 3-vector
       vinf_in_arr  = sol.v2 - v2_planet            # 3-vector

       # Build Encounter for j-th departure body
       encounters.append(Encounter(
         body=body_dep,
         t=t_dep,
         r=r1,
         vinf_out=vinf_out_dep,
         vinf_in=ENCOUNTERS[j-1].vinf_in if j > 0 else None,
       ))
       legs.append(Leg(...))                         # per Cycler dataclass shape

     # Append final encounter
     encounters.append(Encounter(
       body=bodies[-1], t=encounter_times[-1], r=r_last_planet,
       vinf_in=last_vinf_in, vinf_out=None,
     ))

  6. period_sec = catalogue_entry["period_years"] * SECONDS_PER_YEAR
     return Cycler(bodies=bodies, encounters=encounters, legs=legs,
                   period=period_sec, ...)
```

**Notes on the algorithm:**

- **No V∞-magnitude enforcement.** Real-ephemeris Lambert solutions will not perfectly preserve |V∞| at intermediate encounters (the published catalogue values assume circular-coplanar). The mismatch is computed in `_check_vinf_continuity` and reported as a diagnostic; it's not enforced as a hard constraint. A *truly* ballistic cycler should have small mismatches (≲ 0.5 km/s); a degenerate-closure or open trajectory will have large mismatches.
- **No flyby bend re-application.** M6b's construction assumes the catalogue's published encounter sequence is geometrically achievable on real ephemeris. The actual flyby mechanics (turning angle, periapsis radius) are computed implicitly by the Lambert solutions at each intermediate encounter. M7's V3 gate will check flyby feasibility explicitly; M6b leaves that to the V0 (internal consistency) gate already shipped in M5.
- **Period from catalogue.** M6b reads `period_years` from the catalogue rather than re-deriving from `k_synodic * T_synodic`. This is the schema v2 convention; the two should agree to round-trip precision but the catalogue value is the canonical one.

#### 3.1.2 Algorithm — `verify_real_closure`

```text
INPUT: cycler (Cycler or catalogue_entry), n_cycles (int >= 2),
       ephem (Ephemeris), t_start (float | None),
       frame_bodies (tuple[str, ...] | None),
       cycler_id (str | None),
       signature_priority_date (datetime | None)

ALGORITHM:
  1. # Resolve t_start if not provided
     if t_start is None:
       if signature_priority_date is None:
         raise ValueError("t_start or signature_priority_date required")
       signature = phase_signature(cycler)
       t_start = _resolve_real_t_start(signature, ephem, signature_priority_date)
       if t_start is None:
         return RealClosureResult(
           cycler_id=cycler_id, n_cycles_propagated=0,
           max_drift_km=inf, closes=False,
           v3_status="v3-no-real-window", ...
         )

  2. # Optionally reconstruct the cycler from real ephemeris at t_start
     if isinstance(cycler, dict):                        # catalogue-entry path
       try:
         cycler = construct_real_ephemeris_cycler(cycler, ephem, t_start)
       except MultiRevLambertRequiredError as e:
         return RealClosureResult(
           cycler_id=cycler_id, n_cycles_propagated=0,
           max_drift_km=inf, closes=False,
           v3_status="v3-skipped-multirev",
           per_encounter_vinf_mismatch_kms=(),
           ...
         )

  3. # Delegate the propagation + drift to the M6a entry point.
     # This is the binding composition assertion (gate #6).
     stability = verify_long_term_stability(
       cycler, n_laps=n_cycles, ephem=ephem,
       t_start=t_start, frame_bodies=frame_bodies,
       cycler_id=cycler_id, use_uniform_frame=False,
     )

  4. # Translate StabilityReport into RealClosureResult.
     # The drift tolerance is the M6b binding REAL_DRIFT_TOLERANCE_KM,
     # NOT M6a's DRIFT_TOLERANCE_KM.
     closes = stability.max_drift_km < REAL_DRIFT_TOLERANCE_KM
     v3_status = "v3-real-closure-pass" if closes else "v3-real-closure-fail"
     vinf_mismatch = _check_vinf_continuity(cycler, ephem)

  5. return RealClosureResult(
       cycler_id=cycler_id,
       n_cycles_propagated=stability.n_laps_propagated,
       max_drift_km=stability.max_drift_km,
       per_cycle_drift_km=stability.per_lap_drift_km,
       per_encounter_vinf_mismatch_kms=vinf_mismatch,
       closes=closes,
       v3_status=v3_status,
       # M7 placeholders — zeros at M6b, populated at M7
       horizon_tcm_mps=0.0,
       per_cycle_tcm_mps=(0.0,) * n_cycles,
       frame_used=stability.frame_used,
       t_start_sec=t_start,
     )
```

**Notes:**

- **Step 3 is the binding composition.** `verify_real_closure` reuses `verify_long_term_stability` instead of duplicating the multi-lap propagation logic. This is asserted in gate #6 (`test_real_closure_uses_m6a_machinery`) so an accidental reimplementation breaks the test.
- **Step 4's tolerance choice.** The drift comes from M6a's `StabilityReport.max_drift_km` (which was computed against M6a's 50,000 km tolerance), but M6b compares it against `REAL_DRIFT_TOLERANCE_KM = 200,000 km`. The two reports therefore disagree on `stable` vs `closes` for cyclers whose drift sits between 50k and 200k km. **This is intentional**: M6a is the idealised gate (where real-eccentricity breathing is absent and 50k is a tight bound); M6b is the real-ephemeris gate (where breathing is present and 200k absorbs it). The shared `max_drift_km` value is meaningful; the two `stable`/`closes` flags reflect different tolerances.

### 3.2 Catalogue loader filter (`tests/data/_catalogue_loader_m6b.py`)

The loader is **test infrastructure** at M6b — see §1.3. Spec §16.1 schema v2 carries `model_assumption` and `trajectory_regime` fields per entry; M6b filters to:

- `model_assumption in (None, "circular-coplanar")` — accepts entries that omit the field (per schema v2 default rule) and explicit circular-coplanar entries. Rejects `cr3bp` (4 entries) and `analytic-ephemeris` (2 entries). Per spec §12.2: "M6b real-ephemeris optimisation consumes `circular-coplanar` as the *seed* and produces a `real-ephemeris-instance` view stored separately."
- `trajectory_regime in (None, "ballistic")` — v1 catalogue is ballistic-only per the 2026-06-01 schema-v2 backfill; the field is forward-looking for low-thrust. M6b respects the convention even though the field is currently 1-valued.
- `primary in (None, "Sun")` — heliocentric only. Rejects the 3 Earth-primary, 2 Jupiter-primary, 1 Saturn-primary entries (Earth-Moon Arenstorf-class, Jovian / Saturnian CR3BP). The `find_real_windows` slice already raises `NotImplementedError` for non-Sun; M6b filters at load time so the gate test doesn't trip the loader.
- **V1 pass** — entries whose `validation.gates.V1.pass == True` in the catalogue. M6a's V1 gate populates this field (per the M6a hand-off note). Entries not yet V1-checked OR V1-fail are excluded from M6b's regression set.

Expected entry counts:

- Total catalogue entries: 219
- After `model_assumption` filter: 213 (drop 4 `cr3bp` + 2 `analytic-ephemeris`)
- After `primary` filter: ≈ 207 (drop 6 non-Sun entries)
- After V1-pass filter: depends on M6a's reproduction success. The 5-entry regression set is hand-selected from this pool.

**M6B_REGRESSION_IDS** (the 5-entry regression set):

```python
M6B_REGRESSION_IDS: Final[tuple[str, ...]] = (
    "aldrin-classic-em-k1-outbound",     # Aldrin classic, k=1 outbound (M6b binding)
    "aldrin-classic-em-k1-inbound",      # Aldrin classic, k=1 inbound (counterpart)
    "mcconaghy-2006-em-k2",              # McConaghy 2-syn S2L1 — single-rev legs
    "russell-ocampo-2.1.1+2-case2",      # Russell-Ocampo strict-ARMIN/TRMIN, 2-syn
    "russell-ocampo-2.5.1+0",            # Russell-Ocampo wider net, 2-syn
)
```

These are the published-in-literature E-M cyclers most likely to close on real ephemeris with single-rev Lambert legs. The 2-syn S1L1 entry `s1l1-2syn-em-cpom` is **not** in this set because it requires multi-rev Lambert; it lives in `EXPECTED_SKIPS`. See §4.2 for the per-entry expected drift bands and citations.

### 3.3 `RealClosureResult` dataclass + `verify_real_closure` — the M6b locked interface

```python
from dataclasses import dataclass
from datetime import datetime
from typing import Final, Literal

REAL_DRIFT_TOLERANCE_KM: Final[float] = 200_000.0
"""Maximum permissible cycle-to-cycle drift on real DE440 ephemeris over
N >= 2 cycles. Derivation in plan §4.3:

  - 4x M6a's idealised 50,000 km tolerance to absorb real eccentricity
    breathing on E-M trajectories (Earth e ~0.017, Mars e ~0.093).
  - Calibrated against the Aldrin classic k=1 cycler: Pascarella 2024 and
    Russell 2004 report cycle-to-cycle drift on real ephemeris in the
    100,000-200,000 km range before TCMs.
  - Tight enough to reject the spec §10 degenerate-closure basin
    (V_inf > ~11 km/s), which produces AU/cycle drift.

M7 may tighten this if its TCM-budget computation reports post-TCM drift
consistently below 20,000 km."""

N_CYCLES_DEFAULT: Final[int] = 2
"""Default cycle count for the binding gate. Trade-off:
  - n_cycles=2 covers ~4.3 yr for k=1 cyclers, ~8.5 yr for k=2.
  - n_cycles=3 (M6a's lap count) doubles CI runtime per gate test.
  - n_cycles=5 (spec §12(a) horizon) is M7's batch concern.
Plan §6 documents the runtime budget that drives this choice."""


@dataclass(frozen=True)
class RealClosureResult:
    """Result of :func:`verify_real_closure`.

    All fields are immutable. The ``horizon_tcm_mps`` and
    ``per_cycle_tcm_mps`` fields are zeros at M6b; M7 populates them
    when the catalogue ingest writer consumes this result.

    Spec references: §12.1 (the idealised->ephemeris bridge), §12.2
    (real-ephemeris instance representation), §14 V2-real (multi-cycle
    real-ephemeris closure), §16.1 (catalogue ``validation.gates.V2``).

    Attributes
    ----------
    cycler_id:
        Catalogue entry id (e.g. ``"aldrin-classic-em-k1-outbound"``)
        if the result was produced from a catalogue cycler, else
        ``None``. Pass-through; M7's batch-validate runner sets it.
    n_cycles_propagated:
        Number of cycles the propagation actually completed. Equals the
        ``n_cycles`` argument unless an early-termination tripped
        (multi-rev-Lambert blocker, propagator divergence, no real
        launch window).
    max_drift_km:
        Maximum consecutive-cycle-pair drift across the
        ``range(n_cycles - 1)`` pairs in the dynamic rotating frame.
        Sourced from :class:`StabilityReport.max_drift_km`. The basis
        for ``closes``.
    per_cycle_drift_km:
        Cumulative drift at each cycle boundary, sourced from
        :class:`StabilityReport.per_lap_drift_km`. Length
        ``n_cycles_propagated``.
    per_encounter_vinf_mismatch_kms:
        At each interior encounter, ``||V_inf_in| - |V_inf_out||`` in
        km/s — diagnostic measure of how far the real-ephemeris
        construction departs from the published ballistic ideal. Empty
        tuple if construction skipped (multi-rev case).
    closes:
        ``max_drift_km < REAL_DRIFT_TOLERANCE_KM``. The headline
        boolean for spec §14 V2-real (the M6b binding gate).
    v3_status:
        One of ``"v3-real-closure-pass"`` (passes),
        ``"v3-real-closure-fail"`` (drift > tolerance),
        ``"v3-skipped-multirev"`` (multi-rev Lambert blocker hit),
        ``"v3-no-real-window"`` (``find_real_windows`` returned no
        match), or ``"v3-construction-error"`` (Lambert convergence /
        geometry error on a single-rev leg). M7's catalogue writer
        records this in the schema-v2 record's
        ``validation.gates.V2.status`` field.
    horizon_tcm_mps:
        **M6b: 0.0**. M7 populates with the summed TCM ΔV over the
        horizon (m/s). Locked here so M7 doesn't reshape the dataclass.
    per_cycle_tcm_mps:
        **M6b: zero-tuple of length** ``n_cycles_propagated``. M7
        populates with the per-cycle TCM ΔV (m/s). Locked here so M7
        doesn't reshape the dataclass.
    frame_used:
        ``"dynamic"`` always at M6b — the M6a-locked dynamic frame.
        Field exists for forward compatibility with M7 batch runs that
        may exercise the uniform frame for diagnostic comparison.
    t_start_sec:
        The inertial-frame launch epoch (seconds since J2000) the
        result was computed against. Either passed in by the caller
        or derived from :func:`_resolve_real_t_start`.
    """

    cycler_id: str | None
    n_cycles_propagated: int
    max_drift_km: float
    per_cycle_drift_km: tuple[float, ...]
    per_encounter_vinf_mismatch_kms: tuple[float, ...]
    closes: bool
    v3_status: str
    horizon_tcm_mps: float
    per_cycle_tcm_mps: tuple[float, ...]
    frame_used: str
    t_start_sec: float | None


def verify_real_closure(
    cycler: Cycler | dict,
    n_cycles: int,
    ephem: Ephemeris,
    *,
    t_start: float | None = None,
    frame_bodies: tuple[str, ...] | None = None,
    cycler_id: str | None = None,
    signature_priority_date: datetime | None = None,
) -> RealClosureResult:
    """Spec §14 V2-real gate machinery; M6b's binding entry point.

    Pipeline:
      1. Resolve ``t_start`` if not provided
         (``_resolve_real_t_start(signature, ephem, signature_priority_date)``).
      2. If ``cycler`` is a catalogue dict, construct it from real
         ephemeris at ``t_start`` via
         :func:`construct_real_ephemeris_cycler`.
      3. Delegate propagation + drift to
         :func:`cyclerfinder.verify.propagate.verify_long_term_stability`
         with the dynamic frame.
      4. Compare drift against :data:`REAL_DRIFT_TOLERANCE_KM` (not
         M6a's :data:`DRIFT_TOLERANCE_KM`).
      5. Compute V_inf continuity diagnostics.
      6. Return :class:`RealClosureResult` with V3-placeholder fields
         (``horizon_tcm_mps=0.0``, ``per_cycle_tcm_mps=()``) for M7.

    Parameters and exceptions documented per the algorithm in plan
    §3.1.2.
    """
```

**Design notes:**

- **`cycler` can be a `Cycler` or a `dict`.** The `dict` path is the catalogue-driven path (caller passes the YAML entry directly). The `Cycler` path is for callers who built one externally (e.g. M5's `find_cyclers`, when fixed). Type-union accepted in the signature; runtime dispatch on `isinstance`.
- **`closes` is the *binding* boolean.** M7's catalogue writer consumes this directly into `validation.gates.V2.pass`. The four `v3_status` values communicate *why* a failure occurred — multi-rev blocker, no real window, drift exceeded, or construction error — so M7 can mark entries appropriately rather than just "failed."
- **TCM fields stay zero at M6b.** Documented contract: a `RealClosureResult` produced by `verify_real_closure` always has `horizon_tcm_mps == 0.0`. M7's TCM-budget machinery (not built yet) will produce its own `RealClosureResult` variants with populated TCM fields. The frozen dataclass guarantees M7 cannot mutate the M6b-produced result; M7 must construct a new one.

### 3.4 The `optimise_cell_ephemeris` policy

Per §1.7, M6b does **not** fill the body. The only change is the error message:

```python
def optimise_cell_ephemeris(
    cell: Cell,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    n_laps: int = 5,
    n_starts: int = 5,
    seed: int = 0,
    rp_factors: dict[str, float] | None = None,
) -> OptimisationResult:
    """..."""
    del cell, ephem, vinf_cap, n_laps, n_starts, seed, rp_factors
    raise NotImplementedError(
        "optimise_cell_ephemeris requires M6b real-ephemeris closure "
        "(shipped: verify.real_closure.verify_real_closure) AND M7 "
        "TCM budget machinery (not yet shipped). M6b's verify_real_closure "
        "is the right drift-feasibility check; full wiring lands in M7."
    )
```

**Why this is its own §:** the error-message text is what callers (including future automated agents) read to know what to do next. Keeping it precise is load-bearing.

### 3.5 The TCM ΔV deferral — drawing the M6b/M7 boundary

Pascarella et al. 2024 (the architectural template per `docs/v2-future-references.md`) explicitly separates two concerns:

- **Patched-conic / medium-fidelity ballistic closure** — does the trajectory geometrically close on real ephemeris, ignoring corrective ΔV? This is M6b's V2-real gate.
- **TCM budgeting** — given the residual mismatch from medium-fidelity, what is the propellant budget needed to keep the spacecraft on the cycler trajectory over a 20–30 year horizon? This is V3 of the §14 gauntlet.

**M6b ships only the first.** Concretely:

- `RealClosureResult.max_drift_km` is the positional drift after `n_cycles` cycles in the dynamic rotating frame. **No corrective ΔV is applied** — the spacecraft is propagated purely on its constructed initial state, without TCMs at cycle boundaries.
- The TCM-budget computation would require: (a) deciding when (start-of-cycle? mid-cycle? at flybys?) to apply correctives, (b) what cost function to minimise (`sum(|ΔV|)` vs `max(|ΔV|)`), (c) which guidance/reachability formulation to use (impulsive Lambert-target-state vs continuous-feedback control). These are M7 + M-future architectural decisions.
- The V3 placeholder fields on `RealClosureResult` (`horizon_tcm_mps`, `per_cycle_tcm_mps`) are *zeros* at M6b. M7's TCM-budget machinery produces new `RealClosureResult` instances with populated values; the frozen dataclass shape is the contract.

**What this means for the catalogue.** Spec §16.1 carries both `validation.gates.V2.max_drift_km` and `validation.gates.V3.horizon_tcm_mps` fields per entry. M6b populates only the V2 field (via M7's writer consuming `RealClosureResult.max_drift_km`). The V3 field stays null in the catalogue until M7 + future TCM-budget work fills it.

**Justification for the split.** The geometric-closure question ("does this idealised cycler have a real-ephemeris instance at all?") is *prior* to the operational-cost question ("how much propellant does maintaining it cost?"). A cycler that fails V2-real cannot have a meaningful V3 — there's nothing to budget TCMs for. Splitting the milestones means M6b can ship a gate that classifies the catalogue cleanly into closes / open / multi-rev-blocked, and M7 can extend only the closes set without revisiting the others.

### 3.6 Composition with existing code — explicit reuse map

M6b is intentionally a thin layer over existing infrastructure. The composition map:

| Existing function (location) | Used by M6b for |
|---|---|
| `cyclerfinder.verify.propagate.verify_long_term_stability` | The body of `verify_real_closure` (gate #6 binding composition). |
| `cyclerfinder.verify.propagate.multi_lap_propagation` | Called transitively via `verify_long_term_stability`. |
| `cyclerfinder.verify.propagate.lap_to_lap_drift` | Called transitively. |
| `cyclerfinder.verify.propagate._resolve_frame_bodies` | Used to default `frame_bodies` if caller omits — same policy as M6a. |
| `cyclerfinder.core.frames.to_rotating_dynamic` | Called transitively from `lap_to_lap_drift`. |
| `cyclerfinder.core.frames.from_rotating_dynamic` | Not used directly by M6b; used inside M6a's propagator helpers. |
| `cyclerfinder.core.frames.synodic_omega_dynamic` | Called transitively from M6a's machinery. |
| `cyclerfinder.core.ephemeris.Ephemeris("astropy")` | The default ephemeris for `verify_real_closure`; passed through to M6a. |
| `cyclerfinder.core.ephemeris.Ephemeris.state` | Called by `construct_real_ephemeris_cycler` to read real planet positions at encounter epochs. |
| `cyclerfinder.core.kepler.propagate` | Called transitively from `propagate_lap`. |
| `cyclerfinder.core.lambert.lambert` (single-rev only) | Called by `construct_real_ephemeris_cycler` to build each leg's velocity vectors. Multi-rev (`max_revs > 0`) is the documented blocker. |
| `cyclerfinder.search.phase_match.find_real_windows` | Called by `_resolve_real_t_start` to pick a launch epoch matching the cycler's signature near the catalogue's `priority_date`. |
| `cyclerfinder.search.phase_match.phase_signature` | Called to extract the cycler's geometric fingerprint for the `find_real_windows` query. |
| `cyclerfinder.search.phase_match.phase_signature_from_catalogue_entry` | Used in the catalogue-dict path to build the signature without first constructing a `Cycler`. |
| `cyclerfinder.model.cycler.Cycler`, `Encounter`, `Leg` | Constructed by `construct_real_ephemeris_cycler`; consumed by `verify_long_term_stability`. |
| `cyclerfinder.search.optimize.optimise_cell_ephemeris` (stub) | Not called by M6b. Documented as the M7+ wiring target. |
| `cyclerfinder.search.optimize.find_cyclers` | Not called by M6b. M5's binding gate is broken; M6b's Lambert-chain path doesn't need it. |

**New M6b symbols (full list):**

- `verify_real_closure`, `construct_real_ephemeris_cycler` — public functions
- `RealClosureResult` — public frozen dataclass
- `REAL_DRIFT_TOLERANCE_KM`, `N_CYCLES_DEFAULT` — public constants
- `EXPECTED_SKIPS` — public registry dict
- `MultiRevLambertRequiredError`, `RealClosureConstructionError` — public exception types
- `_resolve_real_t_start`, `_construct_cycler_from_lambert_chain`, `_check_vinf_continuity` — module-internal helpers

### 3.7 Imports / dependency graph after M6b

```
constants.py             (M0)
ephemeris.py             (M1 + M6 slice — astropy backend)
lambert.py               (M1 — single-rev only)
kepler.py                (M1)
flyby.py                 (M2)
frames.py                (M3 + M6a)
tisserand.py             (M2)
resonance.py             (M2)
model/cycler.py          (M3)
search/construct.py      (M3)
search/sequence.py       (M4)
model/score.py           (M4)
search/optimize.py       (M5 + M6b one-line edit)
search/phase_match.py    (M6 slice)
verify/propagate.py      (M6a — UNCHANGED)
verify/real_closure.py   (M6b)  <-- frames, ephemeris, kepler, lambert,
                                     verify/propagate, search/phase_match,
                                     model/cycler
```

No new cyclic dependencies. `verify/real_closure.py` sits at the top of the import hierarchy (depends on M0–M6a), so M6b cannot leak into any earlier module's import path.

---

## 4. Tests + gate

Tests live under `tests/verify/test_real_closure.py`, `tests/data/test_catalogue_loader_m6b.py`, and a single addition to `tests/search/test_phase_match.py`. Tolerances are named at the module level (`REAL_DRIFT_TOLERANCE_KM`, `N_CYCLES_DEFAULT`) or hardcoded with cited rationale.

### 4.1 Gate tests (spec §8 binding)

| Test | Assertion | Tolerance |
|---|---|---|
| `test_aldrin_cycler_periodic_over_2_cycles_astropy` (**M6b BINDING GATE — spec §8 M6, real-ephemeris half**) | Load `aldrin-classic-em-k1-outbound` from catalogue; call `verify_real_closure(entry, n_cycles=2, ephem=Ephemeris("astropy"), signature_priority_date=date(1985, 10, 28))`. Assert `result.closes == True`, `result.max_drift_km < REAL_DRIFT_TOLERANCE_KM`, `result.n_cycles_propagated == 2`, `result.v3_status == "v3-real-closure-pass"`, `result.horizon_tcm_mps == 0.0`, `result.per_cycle_tcm_mps == (0.0, 0.0)`. | `max_drift_km < 200_000` |
| `test_2syn_em_cpom_periodic_over_2_cycles_astropy` (**xfail at M6b**) | Same shape, on `s1l1-2syn-em-cpom`. `pytest.mark.xfail(strict=False, reason="multi-rev Lambert blocker — see §5 risk #2")`. The test runs; if it unexpectedly passes (i.e. someone implemented multi-rev Lambert), the xfail is non-strict so the test is treated as xpass rather than failure — the M6b regression record is not destabilised. | n/a (xfail) |
| `test_real_drift_rejects_open_trajectory` | Construct a perturbed Aldrin cycler by rotating one `vinf_out` by 5 degrees; assert the resulting trajectory has `result.closes == False` and `result.max_drift_km > 5 * REAL_DRIFT_TOLERANCE_KM` over 2 cycles. | `max_drift_km > 1_000_000` |
| `test_real_closure_uses_m6a_machinery` (**spec §10 BINDING composition test**) | Mock-patch `cyclerfinder.verify.real_closure.verify_long_term_stability` to count calls; assert that `verify_real_closure(cycler, n_cycles=2, ...)` calls it exactly once with `n_laps=2`. | Call-count exact |
| `test_real_closure_result_frozen_and_v3_fields_locked` | `result.closes = ...` raises `FrozenInstanceError`. `result.horizon_tcm_mps == 0.0` (M6b placeholder). `result.per_cycle_tcm_mps == (0.0,) * n_cycles_propagated`. `result.v3_status` in the 5-valued enum. | bool/exact |

### 4.2 Multi-entry regression set

The 5 regression entries in `M6B_REGRESSION_IDS` are the literature-anchored E-M cyclers most likely to close on real ephemeris with single-rev Lambert. Per-entry expected behaviour:

| Catalogue id | Expected outcome | Expected drift band | Source / citation |
|---|---|---|---|
| `aldrin-classic-em-k1-outbound` | `closes=True` | 50,000 – 200,000 km | Byrnes / Longuski / Aldrin 1993; Rogers 2012 Table 1. The free-return reference. |
| `aldrin-classic-em-k1-inbound` | `closes=True` | 50,000 – 200,000 km | Same lineage; inbound-leg counterpart. The two together exercise the dynamic-frame symmetry. |
| `mcconaghy-2006-em-k2` | `closes=True` | 100,000 – 200,000 km | McConaghy / Landau / Longuski 2006 JSR. Single-rev formulation per the published table. |
| `russell-ocampo-2.1.1+2-case2` | `closes=True` | 50,000 – 150,000 km | Russell & Ocampo 2004 — strict-net 2-syn cycler. |
| `russell-ocampo-2.5.1+0` | `closes=True` | 100,000 – 200,000 km | Russell & Ocampo 2004 — wider-net 2-syn. |

**Regression test** (`test_real_closure_regression_set`):

```python
@pytest.mark.parametrize("entry_id", M6B_REGRESSION_IDS)
def test_real_closure_regression_set(entry_id: str) -> None:
    entries = load_m6b_entries()
    entry = next(e for e in entries if e["id"] == entry_id)
    if entry_id in EXPECTED_SKIPS:
        pytest.skip(EXPECTED_SKIPS[entry_id])
    ephem = Ephemeris("astropy")
    result = verify_real_closure(
        entry, n_cycles=N_CYCLES_DEFAULT, ephem=ephem,
        signature_priority_date=_parse_date(entry["priority_date"]),
    )
    assert result.closes, (
        f"{entry_id}: max_drift_km={result.max_drift_km}, "
        f"v3_status={result.v3_status}"
    )
    assert result.max_drift_km < REAL_DRIFT_TOLERANCE_KM
```

### 4.3 The 200,000 km real-ephemeris drift tolerance — derivation

The M6b tolerance is **binding** at the milestone level — not test-tunable. Derivation:

1. **M6a's tolerance as anchor.** M6a's `DRIFT_TOLERANCE_KM = 50,000` is derived in M6a plan §4.3 as ~0.013° of geometric breathing per lap at Mars's mean radius (~1.5 AU). That tolerance assumed (a) idealised circular-coplanar geometry where breathing is near-zero, (b) 3 laps with the gate cycler having the M5 optimiser's `closure_residual_kms < 0.5` km/s budget integrated over a lap. M6b's setting violates (a): real-ephemeris breathing is the dominant signal, not the noise.

2. **Real eccentricity scaling.** Earth's eccentricity is 0.0167; Mars's is 0.093. Per lap (≈ 2.135 yr for k=1 Aldrin), the planets' real heliocentric positions drift from their circular-coplanar idealisation by `e * a`:
   - Earth: `0.0167 * 1.496e8 km = 2.5 × 10⁶ km` peak excursion
   - Mars: `0.093 * 2.279e8 km = 2.1 × 10⁷ km` peak excursion
   The spacecraft's drift in the rotating frame is a *fraction* of these because the cycler's flyby geometry partially self-corrects, but the fraction is not negligible: published medium-fidelity work (Russell 2004 Tables 4.7–4.9; McConaghy 2006 §3) reports cycle-to-cycle drift in the 10⁴ – 2 × 10⁵ km range for the Aldrin and 2-syn E-M cyclers before TCMs.

3. **Choice of 200,000 km.** Calibration to the upper end of the published medium-fidelity range:
   - 200,000 km / 1.5 AU ≈ 8.9 × 10⁻⁴ rad ≈ 0.051° of geometric breathing per cycle at Mars's mean radius.
   - This is 4× M6a's 50,000 km tolerance and 4× the per-lap breathing band — accounting for (a) M6b's 2 cycles vs M6a's 3 laps (so total integrated drift over the test horizon is comparable), and (b) the real-vs-idealised eccentricity amplification.

4. **N=2 cycle choice.** Trade-off rationale:
   - **n_cycles=1** is degenerate — closing once means matching the start state at the end. Doesn't measure drift.
   - **n_cycles=2** covers 4.3 yr (k=1) to 8.5 yr (k=2). Drift visible if present.
   - **n_cycles=3** covers 6.4 yr to 12.8 yr. Doubles the CI runtime; more drift visible but mostly diagnostic, not gate-relevant.
   - **n_cycles=5** is spec §12(a)'s horizon for TCM optimisation. M7's batch concern.
   M6b binds at n_cycles=2 to (a) keep CI runtime under 60 s per test, (b) measure drift over a meaningful horizon, (c) leave M7's 5-cycle batch horizon open without coupling.

5. **Rejection power.** The degenerate-closure basin (V∞ > ~11 km/s, spec §10) produces open trajectories that diverge by AU/cycle = 1.5 × 10⁸ km/cycle. 200,000 km tolerance rejects these by 3 orders of magnitude.

6. **M7 compatibility.** When M7's TCM-budget optimiser runs, its target is to drive post-TCM drift well below `REAL_DRIFT_TOLERANCE_KM` — typically to 10,000–20,000 km. The M6b tolerance is the **untreated** threshold (pre-TCM); the M7 tolerance is post-TCM. The two coexist as named constants once M7 lands.

The value is exposed as `REAL_DRIFT_TOLERANCE_KM: Final[float] = 200_000.0` at module scope. M7 may add a tighter `POST_TCM_DRIFT_TOLERANCE_KM: Final[float] = 20_000.0` constant for its own gate.

### 4.4 EXPECTED_SKIPS registry

```python
EXPECTED_SKIPS: Final[dict[str, str]] = {
    # Multi-revolution Lambert blocker (see §5 risk #2)
    "s1l1-2syn-em-cpom": (
        "multi-rev Lambert blocker: published S1L1 cycler has n_revs=1 "
        "on intermediate E-E and M-M legs; M1 Lambert solver supports "
        "single-rev only. M-future / stretch."
    ),
    # VEM is M8 scope
    "jones-2017-vem-triple-family": (
        "VEM 3-body real-closure is M8 scope (VEM campaign + viz)."
    ),
    "vem-emeeve-3syn": (
        "VEM 3-body real-closure is M8 scope."
    ),
    # Non-Sun primaries: out of find_real_windows scope
    # (auto-detected by loader filter, not registered here; the
    # filter rejects them at load time so they never reach
    # the regression test).
}
```

**Why this is its own subsection:** the registry IS the M6b architectural statement. "Here is the list of catalogue entries M6b deliberately does not verify, with reason." Every M7 or M-future contributor reads this list to know what's outstanding.

### 4.5 Catalogue-loader tests

| Test | Assertion |
|---|---|
| `test_loader_filters_v1_pass_circular_coplanar_ballistic_sun_only` | `load_m6b_entries()` returns only entries with `model_assumption in (None, "circular-coplanar")` AND `trajectory_regime in (None, "ballistic")` AND `primary in (None, "Sun")`. Count is between 180 and 215 (the exact number depends on M6a's V1-pass-set; bounded for stability). |
| `test_loader_excludes_cr3bp_entries` | No entry in the loader output has `model_assumption == "cr3bp"`. |
| `test_loader_excludes_analytic_ephemeris_entries` | No entry has `model_assumption == "analytic-ephemeris"`. |
| `test_loader_excludes_non_sun_primaries` | No entry has `primary in ("Earth", "Jupiter", "Saturn")`. |
| `test_m6b_regression_ids_all_in_loader` | Every id in `M6B_REGRESSION_IDS` is present in `load_m6b_entries()`. |
| `test_expected_skips_ids_documented` | Every id in `EXPECTED_SKIPS` has a non-empty reason string. |

### 4.6 Helper-level tests

| Test | Assertion |
|---|---|
| `test_construct_real_ephemeris_cycler_aldrin` | `construct_real_ephemeris_cycler(aldrin_entry, ephem, launch_window)` returns a `Cycler` whose `encounters[0].t == _dt_to_t_sec(launch_window.departure_date)` and whose `bodies == ("E", "M", "E")` (Aldrin's encoded sequence). The leg's V∞-mismatch < 1 km/s (single-rev Lambert closes on Aldrin's geometry at the matched epoch). |
| `test_construct_raises_on_multi_rev_leg` | When the catalogue entry's `legs[j].n_revs > 0`, `construct_real_ephemeris_cycler` raises `MultiRevLambertRequiredError(catalogue_id, leg_index=j)`. |
| `test_resolve_real_t_start_prefers_priority_window` | Given the Aldrin signature + `priority_date = 1985-10-28`, `_resolve_real_t_start` returns a `t_sec` value within ±5 years of the priority date. |
| `test_resolve_real_t_start_returns_none_when_no_window` | If `find_real_windows` returns an empty list, `_resolve_real_t_start` returns `None` (caller uses this to short-circuit). |
| `test_check_vinf_continuity_diagnostic` | For the constructed Aldrin cycler, `_check_vinf_continuity(cycler, ephem)` returns a tuple of per-interior-encounter `||V∞_in| - |V∞_out||` values. Aldrin (k=1, 2-encounter) has no interior encounter, so returns `()`. The S2L1 McConaghy entry has 1 interior encounter; returns length-1 tuple. |

### 4.7 Phase-match plumbing test

| Test | Assertion |
|---|---|
| `test_find_real_windows_for_aldrin_signature_within_priority_window` | `find_real_windows(aldrin_signature, Ephemeris("astropy"), (1980-01-01, 1995-01-01), n=3)` returns a non-empty list; at least one window's `departure_date` is within ±5 years of `priority_date = 1985-10-28`. |

### 4.8 Tolerance summary

| Layer | Quantity | Tolerance |
|---|---|---|
| Aldrin cycler real-ephemeris closure over 2 cycles (M6b binding gate) | km drift | < 200,000 (`REAL_DRIFT_TOLERANCE_KM`) |
| Open-trajectory rejection (perturbed Aldrin) | km drift | > 1,000,000 |
| `construct_real_ephemeris_cycler` Aldrin V∞ mismatch | km/s | < 1.0 |
| Composition assertion (verify_long_term_stability calls) | Call count | == 1 exact |
| Regression set per-entry drift | km | < 200,000 for the 5 entries; xfail for `s1l1-2syn-em-cpom` |
| Determinism across re-runs | float | bitwise (modulo ephemeris caching) |

### 4.9 Test runtime budget

Per gate test, on the dev CI runner (single core, astropy DE440 lookup ~10 ms per call):

| Test | Astropy calls | Kepler iters | Expected runtime |
|---|---|---|---|
| `test_aldrin_cycler_periodic_over_2_cycles_astropy` | ~600 (200 samples × 2 cycles × ~1.5 anchor lookups per sample) | ~600 | 6 – 10 s |
| `test_real_drift_rejects_open_trajectory` | Same | Same | 6 – 10 s |
| `test_real_closure_regression_set` (5 entries) | 5x | 5x | 30 – 50 s |
| `test_2syn_em_cpom_periodic_over_2_cycles_astropy` (xfail) | Fails fast on multi-rev | 0 | < 1 s |
| Total M6b real-closure suite | — | — | 50 – 80 s |

Within the CI budget. If the regression set grows past 10 entries, the parametrised tests should switch to `pytest.mark.slow` per the M0 CI policy.

---

## 5. Risks

| # | Risk | Likelihood | Impact | Mitigation in M6b |
|---|---|---|---|---|
| 1 | **Aldrin cycler drift exceeds 200,000 km on real ephemeris (binding gate fails).** Possible causes: (a) `find_real_windows` picks an epoch outside Aldrin's true favourable launch window so the geometry is mismatched from the start; (b) the dynamic frame's `bodies=("E","M")` anchor is wrong (should be Earth-only for Aldrin?); (c) M6a's `verify_long_term_stability` has an undiscovered bug for n_laps=2 (M6a tested at n_laps=3). | medium | **high — milestone-blocking** | First, the test reports `result.max_drift_km`, `result.v3_status`, `result.per_encounter_vinf_mismatch_kms`, and the `t_start` window picked. Mitigations in priority order: (a) loosen `find_real_windows` `mismatch_cap_kms` from 5.0 to 10.0 and check if a better window exists; (b) try `frame_bodies=("E",)` per the M6a hand-off ambiguity; (c) check `verify_long_term_stability(cycler, n_laps=2)` works on circular ephemeris first (no real eccentricity to mask bugs). **Do NOT widen `REAL_DRIFT_TOLERANCE_KM`** — that masks the underlying issue. If drift is 200k–500k, document and proceed with N=2; if drift is > 500k, escalate. |
| 2 | **Multi-rev Lambert blocker affects more catalogue entries than expected.** The catalogue's `legs[].n_revs` field may be sparsely populated; `M6B_REGRESSION_IDS` was chosen by inspecting which entries explicitly have `n_revs=0` on every leg. If a regression entry's published n_revs is actually > 0 (mis-extracted), the test fails with `MultiRevLambertRequiredError`. | medium | medium | The regression test catches this cleanly: `result.v3_status == "v3-skipped-multirev"` is a precise failure mode. Move the affected entry to `EXPECTED_SKIPS` with reason `"multi-rev Lambert blocker"`. Document a follow-up issue for M-future to implement multi-rev Lambert; the catalogue entries needing it become the validation set. |
| 3 | **`find_real_windows` returns no window for an entry's priority_date range.** The function rejects mismatches above `mismatch_cap_kms=5.0`; if a literature cycler's published signature differs from real geometry by > 5 km/s (e.g. McConaghy's published V∞ is for an idealised geometry that doesn't survive eccentricity), the function returns an empty list. | medium | medium | `_resolve_real_t_start` returns `None` in this case; `verify_real_closure` returns `v3_status="v3-no-real-window"`. The catalogue entry is documented as a regression-set member that hit the no-window case. M-future may relax `mismatch_cap_kms` per-entry; M6b's gate cap is 20.0 (per the M6a xfail test's example) which is generous. |
| 4 | **astropy DE440 ephemeris call dominates regression suite runtime.** Per §4.9, the 5-entry regression set is 30–50 s. If the catalogue's circular-coplanar entry count grows (the 213 entries in scope now) and someone parametrises the regression test over all 213, runtime blows past CI budget. | low (in M6b) | medium (in M7) | M6b deliberately limits the regression set to 5 hand-picked entries. The full-catalogue verification belongs to M7's batch runner (offline, not CI). Document the 5-entry vs full-catalogue distinction in `_catalogue_loader_m6b.py`'s module docstring. |
| 5 | **M5 binding gate brokenness (task #54) leaks into M6b via `find_cyclers`.** If a future M6b user calls `find_cyclers` to produce a cycler for `verify_real_closure`, they hit M5's broken state. | low | low | M6b's binding gate uses the catalogue path, not `find_cyclers`. The composition map in §3.6 explicitly documents that `find_cyclers` is not called by M6b. If a future user wires it, they hit task #54 directly — the M5 issue, not an M6b regression. |
| 6 | **`RealClosureResult` shape conflicts with M7's TCM-budget needs.** If M7 needs additional fields (e.g. per-encounter TCM, per-leg ΔV breakdown), the M6b-locked shape forces M7 to either reshape or add a parallel dataclass. | low | low | The frozen dataclass is intentionally minimal. M7 can extend by either (a) adding a new field with a default value, backward-compatible, or (b) adding a new dataclass `TCMReport` and returning both `RealClosureResult` and `TCMReport` from a future top-level function. M7's plan will decide. The locked-at-M6b shape captures only the V2-real gate's outputs. |
| 7 | **Dynamic-frame transform float-precision degrades over multi-cycle horizon.** M6a's frame uses `θ(t) = atan2(r_b0_y, r_b0_x)` read from the ephemeris; over `n_cycles * cycler_period ≈ 4.3 yr * 1.5 = 6.4 yr ≈ 2 × 10⁸ s`, float-precision in `atan2` is ~1e-16 rad — well below the propagator's noise floor. Risk is low. | very low | very low | Already mitigated by M6a's design choice (read θ, not integrate ω). No M6b action needed; documented for completeness. |
| 8 | **Catalogue YAML schema drift** between M6b authorship and execution. If the parallel schema-v3 backfill agent renames `model_assumption` or `trajectory_regime`, M6b's loader filter breaks. | low | medium | The M6b loader uses `entry.get("field", default)` on all v2 fields, treating missing fields as "circular-coplanar" / "ballistic" / "Sun" per the schema-v2 default rules. A field rename (not omission) would break; the loader test `test_loader_filters_v1_pass_circular_coplanar_ballistic_sun_only` catches this. The fix is to update the field name in the loader. |
| 9 | **The Aldrin entry's `priority_date` is in the past (1985)** — `find_real_windows` may not have meaningful planetary geometry that far back if the astropy backend's epoch precision degrades pre-J2000 (it doesn't, but worth flagging). | low | low | DE440 covers 1550-01-01 to 2650-01-22 per the astropy docs. 1985 is well inside. The test asserts a window within ±5 yr of 1985, which is firmly in DE440's range. |
| 10 | **Frozen dataclass + dict-input ambiguity.** `verify_real_closure(cycler, ...)` accepts `cycler: Cycler | dict`. If the dict has missing fields (`bodies`, `legs`, `priority_date`), construction fails with a less-informative error than the catalogue-loader path. | low | low | The dict path documents required fields in the docstring; missing fields raise `ValueError` via `phase_signature_from_catalogue_entry` (already implemented in M6 slice). M6b's `construct_real_ephemeris_cycler` follows the same error pattern. |
| 11 | **M6b's V3 placeholder fields confuse M7 consumers.** `result.horizon_tcm_mps == 0.0` could be mistaken by M7 for "this cycler needs no TCMs" rather than "M6b didn't compute TCMs." | medium | low | Documented on the dataclass: "M6b: 0.0; M7 populates." The M7 catalogue ingestion script must check `v3_status` is `"v3-real-closure-pass"` AND the producer is M7's TCM-budget runner before writing V3 fields into the catalogue. M7's plan will encode this rule. |
| 12 | **Concurrent agents editing `data/catalogue.yaml`.** If the Schema-v2 backfill agent or another agent is modifying entries while M6b's loader reads them, mid-write reads can return malformed YAML. | low | low | M6b's loader uses standard PyYAML which reads the whole file before parsing; partial reads raise `yaml.YAMLError` rather than silently producing garbage data. The loader's caller (test) treats parse errors as test setup failure. Sequential coordination is the parent agent's job. |
| 13 | **Two-cycle horizon is too short to expose slow drift.** A cycler that's stable over 2 cycles but unstable over 5 cycles passes M6b but would fail M7. | low (M6b) | medium (M7) | M7's batch runner uses `n_cycles=5`; the M6b gate is structurally short. The dataclass supports `n_cycles >= 2`, so M7 can extend the same regression set to 5 cycles by changing the parametrisation. Document the cycle-count escalation in the M6b → M7 hand-off. |

---

## 6. Dependency additions

**None.** M6b uses only:

- `numpy` — already in deps.
- `scipy` — already in deps (M5 dep).
- `astropy` — already in deps (M6 slice introduced).
- `pyyaml` — already in deps (catalogue ingest uses it; the M6b loader does too).
- in-house M0–M6a modules.

No new top-level dependencies. No edits to `pyproject.toml`. No `uv.lock` regeneration.

---

## 7. Order of work

The `todo.md` mirrors this with checkboxes.

1. **Re-read predecessor docs.** Re-read spec §12.1 (the bridge `phase_match` ships), §12.2 (representation framework — V1 idealised vs V2 real-ephemeris instance), §14 (V2 multi-cycle gate, V3 TCM gate), §16.1 (catalogue schema v2 fields M6b reads), §16.2 (canonical signature — M6b's outputs inherit the idealised parent's signature, do not generate their own). Re-read M6a's plan §3.3 (StabilityReport shape) and §4.3 (50k tolerance derivation), §5 risks, and the M6a hand-off note in `phases/m6a-idealized-closure-verification/todo.md`. Re-read M6 slice's `search/phase_match.py` end-to-end. Confirm `verify_long_term_stability` with `n_laps=2` works on circular ephemeris (M6a tested at n_laps=3; M6b's n_cycles=2 default needs sanity check).
2. **Catalogue scan / loader scaffold.** Author `tests/data/__init__.py` (empty) + `tests/data/_catalogue_loader_m6b.py` with `load_m6b_entries()` body and `M6B_REGRESSION_IDS` tuple. Author `tests/data/test_catalogue_loader_m6b.py` with the §4.5 filter-tests. Run `uv run pytest tests/data/` green before any verify code lands.
3. **Phase-match plumbing test.** Add `test_find_real_windows_for_aldrin_signature_within_priority_window` to `tests/search/test_phase_match.py`. Confirm green — this is the binding precondition for the M6b gate; if `find_real_windows` doesn't return a window for Aldrin's priority date, M6b cannot construct anything.
4. **Create `verify/real_closure.py` skeleton.** Module docstring with spec §12.1, §12.2, §14 V2-real references. Constants (`REAL_DRIFT_TOLERANCE_KM`, `N_CYCLES_DEFAULT`). Exceptions (`MultiRevLambertRequiredError`, `RealClosureConstructionError`). `RealClosureResult` frozen dataclass per §3.3 with all 11 fields. `EXPECTED_SKIPS` registry. All four public functions stubbed with full docstrings + `NotImplementedError` bodies. Update `verify/__init__.py` re-exports. Run `uv run mypy src/cyclerfinder/verify/real_closure.py` clean on the skeleton.
5. **Implement helpers in dependency order.** Each helper gets a paired unit test before the next:
   - `_resolve_real_t_start(signature, ephem, priority_date)` — `test_resolve_real_t_start_prefers_priority_window`, `test_resolve_real_t_start_returns_none_when_no_window`.
   - `construct_real_ephemeris_cycler(entry, ephem, launch_window)` — `test_construct_real_ephemeris_cycler_aldrin`, `test_construct_raises_on_multi_rev_leg`.
   - `_check_vinf_continuity(cycler, ephem)` — `test_check_vinf_continuity_diagnostic`.
6. **Implement `verify_real_closure`.** Compose the helpers per §3.1.2 pseudocode. Land:
   - `test_real_closure_uses_m6a_machinery` first (the composition assertion — fail-fast if M6b reimplemented propagation).
   - `test_aldrin_cycler_periodic_over_2_cycles_astropy` second (**M6b BINDING GATE**). If it fails, escalate per risk #1 (DO NOT widen `REAL_DRIFT_TOLERANCE_KM`).
   - `test_2syn_em_cpom_periodic_over_2_cycles_astropy` (xfail) third.
   - `test_real_drift_rejects_open_trajectory` fourth.
   - `test_real_closure_result_frozen_and_v3_fields_locked` fifth.
7. **Run the regression set.** Parametrise `test_real_closure_regression_set` over `M6B_REGRESSION_IDS` minus `EXPECTED_SKIPS`. If any entry fails:
   - If `v3_status == "v3-skipped-multirev"`: move to `EXPECTED_SKIPS`; document.
   - If `v3_status == "v3-no-real-window"`: investigate `find_real_windows` per risk #3; if a manual mismatch_cap relax works, parametrise per-entry; else move to `EXPECTED_SKIPS`.
   - If `v3_status == "v3-real-closure-fail"`: investigate per risk #1; if drift is 200k–500k, document as borderline; if > 500k, treat as gate failure and escalate.
8. **Edit `search/optimize.py::optimise_cell_ephemeris` error message.** Update the `NotImplementedError` text per §3.4. Confirm `tests/search/test_optimize.py` still passes (the error-message text should not be asserted; only the exception class).
9. **Run the full local quality gate:** `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`.
10. **Commit** as `m6b: real-ephemeris closure verification (verify_real_closure + Aldrin gate over 2 cycles)`. Push; confirm CI green.
11. **Update `docs/overview.md`.** §4 milestone table: M6b → completed; M7 → planned. Add note that real-closure verification ships independently of TCM-budget (the M6b/M7 boundary).
12. **Hand-off note** appended to `todo.md` under `## Hand-off to M7` covering: actual `max_drift_km` reproduced for the Aldrin gate, regression-set per-entry drifts, EXPECTED_SKIPS final list (with reasons), the locked `RealClosureResult` shape that M7 inherits, the M7 catalogue-writer contract for V2 / V3 fields.

The order is "catalogue loader → phase-match plumbing → real_closure scaffold → helpers → gate → regression" deliberately: each step exercises an existing module before the next ships, the M6b binding gate lands after the helper-level tests so any failure has a localisable cause, and the regression set runs last as the broad coverage.

---

## 8. Exit checklist (the gate, restated)

Before declaring M6b done:

- [ ] `uv run pytest tests/verify/test_real_closure.py` green; the M6b binding gate `test_aldrin_cycler_periodic_over_2_cycles_astropy` passes with `closes=True`, `max_drift_km < 200,000` km, `n_cycles_propagated == 2`, `v3_status == "v3-real-closure-pass"`.
- [ ] `uv run pytest tests/verify/test_real_closure.py::test_2syn_em_cpom_periodic_over_2_cycles_astropy` is **xfail** (the multi-rev blocker). If it xpasses (someone landed multi-rev Lambert), document and consider tightening; the xfail is non-strict so xpass doesn't break the gate.
- [ ] `uv run pytest tests/verify/test_real_closure.py::test_real_closure_regression_set` green for all 5 entries minus EXPECTED_SKIPS.
- [ ] `uv run pytest tests/data/test_catalogue_loader_m6b.py` green.
- [ ] `uv run pytest tests/search/test_phase_match.py` green (including the new Aldrin-window test).
- [ ] `uv run pytest` green overall (no regression of M0–M6a tests; M6a `test_propagate.py` unchanged).
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean under `strict = true` — including `RealClosureResult`, `Final[float]` on `REAL_DRIFT_TOLERANCE_KM`, NDArray annotations on construction helper outputs.
- [ ] CI green on the M6b commit.
- [ ] `docs/overview.md` updated: M6b status = `completed`; M7 row updated to reflect M6b's hand-off; note added in §4 explaining the M6b/M7 boundary (geometric closure ships at M6b; TCM budget is M7).
- [ ] `## Hand-off to M7` section appended to `phases/m6b-real-ephemeris-closure/todo.md` covering:
  - The exact `max_drift_km` reproduced for the Aldrin gate (vs the 200,000 km bound), and which cycle pair was worst.
  - The regression-set per-entry results: `closes` status, `max_drift_km`, `v3_status`, any escalations.
  - The final `EXPECTED_SKIPS` list with reasons — this is the input M7's "outstanding entries" follow-up consumes.
  - Whether `find_real_windows`'s `mismatch_cap_kms` defaults were sufficient for the regression set, or whether per-entry overrides were needed.
  - Per-test wall-clock runtime for the M6b binding gate and the regression set — informs M7's batch-runner compute budget.
  - The frame-bodies decision: was the M6a `_resolve_frame_bodies` policy (`("E", "M")`) right for the Aldrin gate, or did the gate require an override?
  - Whether `verify_long_term_stability` at `n_laps=2` showed any anomaly vs the M6a-tested `n_laps=3` (sanity check; M6a tested at 3, M6b uses 2).
  - The locked `RealClosureResult` shape (all 11 fields) so M7 knows it does NOT need to reshape the dataclass — only populate `horizon_tcm_mps`, `per_cycle_tcm_mps`, and write `v3_status` to the catalogue.
  - The contract for M7's catalogue writer: it MUST check `v3_status == "v3-real-closure-pass"` before writing the V2-real fields into the catalogue; it MUST NOT write any V3 field until M7's TCM-budget runner produces a new `RealClosureResult` with populated TCM fields.
  - The recommended M7 first steps: (1) read this hand-off + spec §12(a) + spec §14 V3, (2) wire `optimise_cell_ephemeris` against `verify_real_closure` as drift-feasibility check, (3) define TCM-budget computation (which encounters / what guidance / what cost), (4) extend `RealClosureResult` via additive field or define `TCMReport` parallel dataclass.

(Writing the M7 plan doc is the first task of M7, not an M6b exit criterion.)
