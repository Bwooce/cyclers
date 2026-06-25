# Low-thrust / DSM releg genome — IMPLEMENTATION PLAN (#449)

**Date:** 2026-06-25
**Design draft:** `docs/superpowers/specs/2026-06-25-449-lowthrust-dsm-releg-genome-design-draft.md`
**Status:** PLAN — bite-sized TDD tasks. No code written by this doc.

Each task is self-contained: a failing test first, the minimal implementation,
a verify step, and a single pathspec commit. Reuse over rebuild — the DSM leg
solver (`search/dsm_leg.py`, #307), the SF low-thrust leg solver
(`core/sims_flanagan.py` + `search/lowthrust.py`, #309), and the VILM cost model
(`search/vilm.py`) already exist and are tested; the new work is the **releg
SWAP + moon-tour driver**, not a new optimiser.

Conventions: `uv run pytest <path> -q` to run tests; `uv run ruff check . &&
uv run ruff format --check .` before every commit (`feedback_run_ruff_before_commit`);
pathspec commits only, NEVER `git add -A` (sibling #450 is committing to this same
tree); no `Co-Authored-By`. Prefix multi-step bash with `date -Iseconds`.

---

## Task 0 — (OPTIONAL, gate) confirm the swap-seam contract is dependency-injectable

**Why first:** the design rests on `_close_one_phasing` already taking an injected
`lambert` callable (line 463) and `_cycle_residual` calling `_lambert` directly.
Confirm before building so the swap is truly a one-call-site change.

- **Verify (read-only, no commit):** `grep -n 'lambert' src/cyclerfinder/search/discovery_campaign.py src/cyclerfinder/data/validation/v2_moontour.py`. Confirm discovery-side line ~494 uses the injected `lambert` param and validation-side line ~275 uses the module `_lambert`.
- **Output:** a one-line note in the PR/issue confirming both seams; if the validation side is NOT injectable, Task 4 adds a tiny injection param (still minimal).
- **No commit** (read-only confirmation).

---

## Task 1 — `RelegResult` + `Releg` protocol + `BallisticReleg` (regression-preserving)

**Goal:** the swap interface and the zero-ΔV backend that reproduces today's
ballistic leg exactly.

- **Failing test:** `tests/search/test_releg_solver.py::test_ballistic_releg_matches_lambert_path`
  - Build a fixed Jovian Io→Europa leg (planet-frame moon states from
    `core.satellites` + a `_moon_state`-style helper, a chosen tof, mu=Jupiter).
  - Assert `BallisticReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)` returns
    a `RelegResult` whose `vinf_out`/`vinf_in` equal the values computed by calling
    `core.lambert.lambert` directly and selecting the lowest-energy branch (the
    current `_close_one_phasing` logic, lines 494-501), and `dv_kms == 0.0`,
    `feasible is True`.
- **Minimal impl:** `src/cyclerfinder/search/releg_solver.py`
  - `@dataclass(frozen=True) RelegResult(vinf_out, vinf_in, dv_kms, feasible)`.
  - `class Releg(Protocol)` with `solve(...) -> RelegResult`.
  - `class BallisticReleg` wrapping `core.lambert.lambert` + lowest-energy
    branch pick (lift the exact selection logic from `_close_one_phasing`).
- **Verify:** `uv run pytest tests/search/test_releg_solver.py -q`; ruff.
- **Commit (pathspec):** `git add src/cyclerfinder/search/releg_solver.py tests/search/test_releg_solver.py && git commit -m "search: #449 Releg protocol + BallisticReleg (regression-preserving leg seam)"`

---

## Task 2 — `DsmReleg` backend (powered, V∞-retargeting) + VILM-floor golden

**Goal:** the primary powered backend; one DSM per leg, retargets arrival V∞,
reports delivered ΔV, golden-anchored to the VILM floor.

- **Failing test 1 (capability):** `test_dsm_releg_retargets_vinf`
  - On a leg where the ballistic arrival V∞ ≠ the requested arrival V∞, assert
    `DsmReleg().solve(..., vinf_target_in=<value>)` returns `feasible`, an arrival
    V∞ within tolerance of the target, and `dv_kms > 0`.
- **Failing test 2 (golden, sourced):** `test_dsm_releg_dv_geq_vilm_floor`
  - For a Ganymede→Europa leg, assert the delivered `dv_kms` is **≥**
    `search.vilm.vilm_dv_min("Ganymede", "Europa") − tol` and within the powered
    band. EXPECTED side cites Campagnola-Russell Part-1 (digest
    `2026-06-05-endgame-tisserand-mining.md`) — not a self-computed number.
- **Minimal impl:** add `class DsmReleg(Releg)` to `releg_solver.py`
  - Wrap `search.dsm_leg.dsm_leg(r0=r_a, v0=<departure v from a seed branch>,
    tof=tof_s, eta=<optimised>, target_r=r_b, mu=mu, max_revs=n_rev)`; optimise
    `eta` (and rev branch) to hit `vinf_target_in`; `dv_kms = dsm_dv_kms`.
- **Verify:** `uv run pytest tests/search/test_releg_solver.py -q`; ruff.
- **Commit (pathspec):** `... -m "search: #449 DsmReleg backend — V∞-retargeting one-DSM leg, VILM-floor golden"`

---

## Task 3 — `data/golden/campagnola_endgame_releg.yaml` + golden-load test

**Goal:** freeze the sourced EXPECTED values once, consumed by Tasks 2/5.

- **Failing test:** `tests/data/test_golden_campagnola_releg.py::test_golden_values_sourced`
  - Load the YAML; assert it carries Part-1 Table 1 (no-GA) + Table 2 (with-GA)
    ΔV_min targets for ≥2 moon pairs, the Europa endgame 154/147 m/s, and a
    `disjoint_contour_pairs` list (Ariel-Umbriel etc.) with `bridgeable: false`.
    Assert each value has a `source` field tracing to the paper/mining-note line.
- **Minimal impl:** write `data/golden/campagnola_endgame_releg.yaml` with values
  transcribed from `vilm.py`'s already-validated outputs' *sources* (NOT from
  re-running `vilm.py` — cite the paper directly per `feedback_golden_tests_sourced_only`).
- **Verify:** `uv run pytest tests/data/test_golden_campagnola_releg.py -q`; ruff.
- **Commit (pathspec):** `... -m "data: #449 sourced golden for Campagnola-Russell endgame relegs"`

---

## Task 4 — `releg_moontour.py` driver: close a cycle with a chosen backend + VILM-floor prefilter

**Goal:** loop the legs of a tour skeleton with a `Releg` backend, enforce
post-retarget V∞-continuity + closed-cycle wrap, sum per-cycle ΔV, skip
structurally-dead legs cheaply.

- **Failing test 1 (positive control):** `tests/search/test_releg_moontour.py::test_jovian_positive_control_closes_powered`
  - Io-Europa-Ganymede-Io skeleton (registry positive control, links at vinf=4):
    with `DsmReleg`, assert the cycle closes (continuity residual below ballistic
    gate post-retarget) with total `dv_kms` inside the powered dv-band.
- **Failing test 2 (prefilter / structural):** `test_uranus_disjoint_prefiltered_empty`
  - Ariel-Umbriel-Ariel skeleton: assert the VILM-floor prefilter marks it
    `unbridgeable` (via `vilm.min_vinf_for_vilm` / `vilm_dv_min` exceeding the
    band) and the driver returns an EMPTY verdict WITHOUT running the DSM solve.
- **Minimal impl:** `src/cyclerfinder/search/releg_moontour.py`
  - `close_powered_cycle(primary, sequence, leg_tofs_days, n_revs, releg: Releg,
    phasing, *, dv_band) -> PoweredCycleVerdict` (frozen dataclass:
    per-leg ΔV, total ΔV, continuity residual, feasible, prefilter_skipped).
  - Reuse `_moon_state` (factor a shared helper or import from
    `discovery_campaign`) + `core.satellites` for moon states; reuse the wrap
    continuity definition from `_close_one_phasing`.
- **Verify:** `uv run pytest tests/search/test_releg_moontour.py -q`; ruff.
- **Commit (pathspec):** `... -m "search: #449 releg_moontour driver — powered cycle close + VILM-floor prefilter"`

---

## Task 5 — wire the swap into the discovery seam (injection, regression-locked)

**Goal:** prove `BallisticReleg` injected at `_close_one_phasing` reproduces
today's result bit-for-bit, and `DsmReleg` is a drop-in.

- **Failing test:** `tests/search/test_releg_swap_seam.py::test_ballistic_releg_swap_matches_baseline`
  - Run `RepeatedMoonTarget._close_one_phasing` on a fixed skeleton with (a) the
    current `lambert` callable and (b) a thin adapter that routes the same call
    through `BallisticReleg`; assert identical `(feasible, worst, vinf, tofs)`.
- **Minimal impl:** add a `releg`-aware adapter (a callable matching the existing
  injected-`lambert` signature, or a tiny `leg_solver` param) so `BallisticReleg`
  is call-compatible. Do NOT change default behaviour — the default stays
  `lambert` so all existing tests pass untouched.
- **Verify:** `uv run pytest tests/search/test_releg_swap_seam.py tests/search/test_discovery_campaign*.py -q` (regression); ruff.
- **Commit (pathspec):** `... -m "search: #449 inject Releg at moontour discovery seam (BallisticReleg = baseline)"`

---

## Task 6 — releg-aware V2 moontour gate (validation seam) + powered dv-band classification

**Goal:** the V2 moontour gauntlet accepts a powered cycle and classifies it by
the powered dv-band, not the ballistic closure floor.

- **Failing test:** `tests/data/test_v2_moontour_powered.py::test_powered_cycle_passes_v2_in_band`
  - A powered Jovian positive-control candidate: assert a releg-aware
    `run_v2_moontour(..., releg=DsmReleg(), dv_band="powered")` passes ≥3 cycles,
    drift-bounded, with per-cycle ΔV inside `verify/dv_band_acceptance.py`'s
    powered window; and a too-expensive cycle is FAILED (honest negative).
- **Minimal impl:** thread an optional `releg`/`leg_solver` + `dv_band` param into
  `v2_moontour._cycle_residual` (default = today's `_lambert`, ballistic — all
  existing v2_moontour tests stay green) and couple the ΔV sum to
  `dv_band_acceptance`.
- **Verify:** `uv run pytest tests/data/test_v2_moontour*.py -q` (regression + new); ruff.
- **Commit (pathspec):** `... -m "data: #449 releg-aware V2 moontour gate + powered dv-band classification"`

---

## Task 7 — `LowThrustReleg` backend (optional second backend, gated behind DSM)

**Goal:** the Sims-Flanagan low-thrust leg as a swappable backend; bracket-golden
only (the DSM branch carries the tight golden).

- **Failing test:** `tests/search/test_releg_solver.py::test_lowthrust_releg_brackets_dsm`
  - On a feasible positive-control leg, assert `LowThrustReleg().solve(...)` is
    `feasible` and its `dv_kms` brackets the DSM/VILM-floor result within a
    documented tolerance (SF is a different transcription — bracket, not equal).
  - Mark with a generous timeout; if SF cannot converge in bounded time, the test
    `xfail`s with a recorded reason (ship DSM-only per design §7).
- **Minimal impl:** add `class LowThrustReleg(Releg)` wrapping
  `core.sims_flanagan.SimsFlanaganLeg` + `search.lowthrust.solve_leg_min_dv`;
  boundary states pinned to the moon encounter states retargeted to the needed V∞.
- **Verify:** `uv run pytest tests/search/test_releg_solver.py -q`; ruff.
- **Commit (pathspec):** `... -m "search: #449 LowThrustReleg backend (Sims-Flanagan, bracket-golden, gated behind DSM)"`

---

## Task 8 — ratchets + docs + registry re-stamp plumbing

**Goal:** keep the frozen-census ratchets green and document the capability.

- **Run ALL ratchets** (`feedback_catalogue_edits_run_all_ratchets`):
  `uv run pytest tests/data tests/search -q`. This plan adds NO catalogue.yaml row
  (no discovery run here — that is a separate campaign issue), so census ratchets
  should be untouched; confirm they stay green.
- **Empty-region re-stamp helper test:** `test_powered_empty_restamp_records_method`
  — assert the driver, on a powered-empty region (Uranus), produces a
  capability-subsumption re-stamp record (method+version+git_sha) suitable for
  appending to `empty_regions.jsonl` — WITHOUT writing it here (writeback is a
  campaign-issue action, not a plan action).
- **Docs:** add a short capability note to `data/OUTSTANDING.md` /
  the relevant capability index referencing the design draft
  (`feedback_update_docs_proactively`).
- **Verify:** full `uv run pytest tests/data tests/search -q`; ruff.
- **Commit (pathspec):** `... -m "search/docs: #449 powered-empty re-stamp record + capability note (no writeback)"`

---

## Out of scope for this plan (separate issues)

- **The discovery CAMPAIGN** (actually relegging the `repeated-moon-*-sweep`
  skeletons at scale and re-stamping the registry) — a follow-on run issue, not a
  build issue. This plan ships the *capability + golden*; the campaign spends it.
- **3D / inclined relegs** (the Amalthea-inclination half of the re-open keys) —
  a sibling capability (the registry names "3D/inclined relegs" separately from
  "low-thrust relegs"); out of scope for #449.
- **#450 DA/HOTM** — independent (design §5); no coupling.

---

## Definition of done

- DSM releg backend + driver build green; the Campagnola-Russell leveraging golden
  passes; the Jovian positive control closes inside the powered band; the Uranian
  disjoint case is honestly reported empty.
- `BallisticReleg` reproduces the baseline ballistic leg bit-for-bit (no
  regression in discovery_campaign / v2_moontour).
- All `tests/data` + `tests/search` ratchets green.
- The "blocked by #450" tag on #449 is dropped (design §5 finding).
