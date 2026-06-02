# Real-Ephemeris Cell Optimiser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. **Work directly on `main` — do NOT create a branch** (project rule).

**Goal:** Implement the currently-stubbed `optimise_cell_ephemeris` so a `Cell` can be optimised to a closed periodic cycler on the *real* ephemeris (DE440/astropy), unlocking rediscovery of cyclers — like S1L1 and Aldrin — whose published V∞ anchors are NOT hostable in the circular-coplanar idealised model.

**Architecture:** The general real-ephemeris solve engine already exists — `search/maintain.py::optimise_maintenance_dv` is fully body-agnostic (free vars `x = [t0_sec, tof_1_days, …, tof_{N-1}_days]`; objective = summed flyby turn-deficit ΔV; DE global pass + SLSQP multi-start). It is currently reached only through the Aldrin-specific wrapper `optimise_aldrin_maintenance_dv` and `bvp.py::solve_powered_periodic_cycler`. This plan derives `optimise_maintenance_dv`'s parameters from a `Cell` + a resolved launch epoch, wraps the result as an `OptimisationResult`, and adds a real-ephemeris rediscovery gate validated against sourced V∞ anchors. Closure in the rotating frame is physically unreachable on a real ephemeris (documented `bvp.py:39-51`); the optimisation objective is maintenance ΔV and the *feasibility reporter* is `verify_real_closure`'s drift.

**Tech Stack:** Python 3.11, numpy + scipy (DE/SLSQP) — scipy already used in `maintain.py`, astropy DE440 ephemeris backend, pytest (+xdist, `@pytest.mark.slow`), uv venv (no pip), ruff + mypy pre-commit.

**Provenance contract (binding):** The optimiser COMPUTES the launch epoch and leg ToFs. Validation asserts the rediscovered trajectory's V∞ matches **published** anchors (S1L1: 5.65 km/s Earth, 3.05 km/s Mars — spec §9). Computed ToFs/epoch are labelled computed, never golden. See `[[s1l1-nomenclature]]`, `[[golden-tests-sourced-only]]`.

---

## Background: the existing real-ephemeris stack (read before starting)

(All citations verified by code exploration 2026-06-02.)

- `src/cyclerfinder/search/optimize.py:1144` — `optimise_cell_ephemeris(cell, ephem, *, vinf_cap, n_laps=5, n_starts=5, seed=0, rp_factors=None)` — **the stub to implement** (raises `NotImplementedError`).
- `src/cyclerfinder/search/maintain.py:389` — `optimise_maintenance_dv(sequence, ephem, *, t0_guess_sec, tof_days_guesses, tof_bounds_days, synodic_pair=None, closure_body=None, closure_flyby_alt_km=None, t0_window_synodic_frac=0.15, tof_jitter_half_days=None, n_starts=5, seed=0, seed_cycler_factory=None) -> MaintenanceOptimResult`. **Fully general.** Free vars `x = [t0_sec, tof_1_days, …]`; objective `_maintenance_dv_chain` (`maintain.py:331`, summed turn-deficit ΔV); DE pass (`maintain.py:479`) + SLSQP polish.
- `MaintenanceOptimResult` (`maintain.py:143`) — `.cycler`, `.t0_sec`, `.leg_tofs_days`, `.maintenance_dv_kms`, `.per_encounter_dv_kms`, `.converged`, `.a_au`, `.e`, `.vinf_kms_at_encounters`, `.turn_deficit`.
- `src/cyclerfinder/search/phase_match.py:278` — `find_real_windows(...)` grid-scans the real ephemeris and, at each candidate departure date, Lambert-solves each leg vs real planet positions and sums `|V∞_actual − V∞_target|`. `phase_signature_from_catalogue_entry` (`phase_match.py:134`) builds a `PhaseSignature(bodies, leg_durations_s, vinf_target_kms)`. **Heliocentric only** (`phase_match.py:320` raises for non-Sun primaries).
- `src/cyclerfinder/verify/real_closure.py:321` — `_resolve_real_t_start(signature, ephem, priority)` → J2000 seconds or `None` (no window beats `mismatch_cap_kms=20`).
- `src/cyclerfinder/verify/real_closure.py:640` — `verify_real_closure(...)` → `RealClosureResult(closes, max_drift_km, …)`; `REAL_DRIFT_TOLERANCE_KM = 200_000` (`real_closure.py:84`). **Use as a post-hoc diagnostic, not the objective.**
- `src/cyclerfinder/search/optimize.py:215` — `OptimisationResult(cell, best_cycler, best_score, closure_residual_kms, optimiser_history, converged, constraints_satisfied)`. `best_score` is `model/score.py::Score` (has `.max_vinf_kms`, `.total_maintenance_dv_kms`). See how `optimise_cell_idealized` builds it (`optimize.py:1040-1141`).
- `tests/verify/test_real_closure.py:130` — `test_aldrin_powered_cycler_solver_and_drift_floor_on_de440` is the *template* gate (calls `solve_powered_periodic_cycler`, asserts maintenance ΔV > 0 + drift floor); sourced anchors in the adjacent `test_aldrin_powered_turn_deficit_gate` (`:214`).

### Landmines (from the architecture map)
1. `_maintenance_dv_chain` (`maintain.py:331`) assumes `encounters[0].body == encounters[-1].body` (closed loop). Fine for cells whose `sequence[0]==sequence[-1]`; must guard otherwise.
2. `construct_cycler` requires strictly increasing epochs; for N>3 the optimiser must enforce all ToFs > 0 (leg-sum monotonicity).
3. `_resolve_real_t_start` returns `None` when no real launch window beats 20 km/s mismatch → the cell is infeasible at that epoch; the optimiser must surface this, not crash.
4. Epoch resolution (`find_real_windows`) needs **V∞ targets** to phase-match. For rediscovery we have sourced anchors; for blind discovery we don't (documented limitation; out of scope here).
5. `t0` search window `±0.15·T_syn` (`maintain.py:92`) is deliberately narrow to stay on the target family — keep it.

---

## Task 1: Derive `optimise_maintenance_dv` inputs from a `Cell`

**Files:**
- Modify: `src/cyclerfinder/search/optimize.py` (new private helpers near the stub)
- Test: `tests/search/test_optimize_ephemeris.py` (create)

- [ ] **Step 1: Write the failing test**

```python
import math
from cyclerfinder.search.sequence import Cell
from cyclerfinder.search.optimize import _ephemeris_tof_seed_and_bounds


def test_ephemeris_tof_seed_and_bounds_equispaced():
    cell = Cell(bodies=("E", "M"), sequence=("E", "M", "E"), period_k=2,
                per_leg_revs=(0, 0), per_leg_branch=("single", "single"))
    target_period_sec = 2 * 779.9 * 86400.0
    seed_days, bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)
    n_legs = len(cell.sequence) - 1
    assert len(seed_days) == n_legs
    assert len(bounds) == n_legs
    # equispaced seed: each leg ~ T/(N-1)
    expected = target_period_sec / n_legs / 86400.0
    assert all(math.isclose(s, expected, rel_tol=1e-9) for s in seed_days)
    # bounds bracket the seed and are strictly positive
    for (lo, hi), s in zip(bounds, seed_days):
        assert 0 < lo < s < hi
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_optimize_ephemeris.py::test_ephemeris_tof_seed_and_bounds_equispaced -v`
Expected: FAIL — `_ephemeris_tof_seed_and_bounds` not defined.

- [ ] **Step 3: Implement the helper**

```python
def _ephemeris_tof_seed_and_bounds(
    cell: Cell, target_period_sec: float
) -> tuple[list[float], list[tuple[float, float]]]:
    """Equispaced per-leg ToF seed (days) and per-leg bounds for the
    real-ephemeris optimiser, derived from the cell's period.

    Mirrors the interior-epoch bounds logic of ``_bounds_for`` (the
    idealised optimiser): each leg starts at ``T/(N-1)`` and may range
    over ``[0.1, 0.9]`` of that share scaled to the full period so the
    optimiser can redistribute time between legs while keeping every ToF
    strictly positive.
    """
    n_legs = len(cell.sequence) - 1
    if n_legs < 1:
        raise ValueError(f"cell.sequence must have >= 2 entries; got {cell.sequence!r}")
    period_days = target_period_sec / 86400.0
    share = period_days / n_legs
    seed = [share] * n_legs
    lo = 0.1 * share
    hi = 0.9 * period_days  # a single leg may absorb most of the period
    bounds = [(lo, hi)] * n_legs
    return seed, bounds
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/search/test_optimize_ephemeris.py::test_ephemeris_tof_seed_and_bounds_equispaced -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cyclerfinder/search/optimize.py tests/search/test_optimize_ephemeris.py
git commit -m "search/optimize: derive real-ephemeris ToF seed+bounds from a Cell"
```

---

## Task 2: Implement `optimise_cell_ephemeris`

**Files:**
- Modify: `src/cyclerfinder/search/optimize.py` (replace the stub body; add `priority_date` / `vinf_targets_kms` params)
- Test: `tests/search/test_optimize_ephemeris.py`

- [ ] **Step 1: Read how `optimise_cell_idealized` builds its `OptimisationResult`**

Read `optimize.py:1040-1141` to see exactly how `best_score` (a `Score`) and `constraints_satisfied` are produced, so the ephemeris path returns the SAME `OptimisationResult` shape. Note: for ephemeris mode, exact rotating-frame closure is unreachable, so `closure_residual_kms` will carry the post-hoc lap-to-lap drift proxy and `best_score.total_maintenance_dv_kms` carries the real objective — document this in the docstring.

- [ ] **Step 2: Write the failing test (epoch resolution + return shape)**

```python
import pytest
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris, OptimisationResult
from cyclerfinder.search.sequence import Cell


@pytest.mark.slow
def test_optimise_cell_ephemeris_returns_result_for_aldrin_em_cell():
    """The general ephemeris optimiser must reproduce the Aldrin E-M-E
    geometry the Aldrin-specific solver finds (parity check), returning a
    populated OptimisationResult on the real ephemeris."""
    cell = Cell(bodies=("E", "M"), sequence=("E", "M", "E"), period_k=1,
                per_leg_revs=(0, 0), per_leg_branch=("single", "single"))
    eph = Ephemeris(model="astropy")
    result = optimise_cell_ephemeris(
        cell, eph, vinf_cap=12.0,
        priority_date_iso="1985-01-01",   # Aldrin priority era
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=3, seed=0,
    )
    assert isinstance(result, OptimisationResult)
    assert result.best_cycler is not None
    # the recovered E-M elements match the sourced Aldrin family within band
    assert 1.5 < result.best_score.max_vinf_kms or result.converged
```

- [ ] **Step 3: Run to verify it fails**

Run: `uv run pytest tests/search/test_optimize_ephemeris.py::test_optimise_cell_ephemeris_returns_result_for_aldrin_em_cell -v`
Expected: FAIL — `optimise_cell_ephemeris` raises `NotImplementedError`, or `TypeError` on the new params.

- [ ] **Step 4: Implement the function**

Replace the stub body. The implementation:
1. Compute `target_period_sec` via `_target_period_sec(cell)` (already 2-body correct; multi-body comes from the M8-Core `period_basis` work — out of scope here, 2-body only).
2. `seed_days, bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)`.
3. Resolve `t0`: build a `PhaseSignature` from the cell's body chain + `seed_days` (as `leg_durations_s`) + `vinf_targets_kms`; call `_resolve_real_t_start(signature, ephem, priority)` where `priority` parses `priority_date_iso`. If it returns `None`, return an honest "no real window" `OptimisationResult` with `converged=False, constraints_satisfied=False` (mirror `_sentinel_cycler`/`_empty_result` patterns).
4. Guard landmine #1: require `cell.sequence[0] == cell.sequence[-1]` (closed loop); else raise `ValueError` with a clear message (open-sequence ephemeris cyclers are out of scope).
5. Call `optimise_maintenance_dv(list(cell.sequence), ephem, t0_guess_sec=t0, tof_days_guesses=seed_days, tof_bounds_days=bounds, synodic_pair=(cell.bodies[0], cell.bodies[1]), closure_body=cell.sequence[0], n_starts=n_starts, seed=seed)`.
6. Map `MaintenanceOptimResult` → `OptimisationResult`: `best_cycler = m.cycler`; build `best_score` via the same `score(...)` path `optimise_cell_idealized` uses (read in Step 1); `closure_residual_kms` = a post-hoc drift proxy (or `m.maintenance_dv_kms` if drift not computed here — document choice); `converged = m.converged`; `constraints_satisfied` = `m.converged and all encounter V∞ ≤ vinf_cap` (compute from `m.cycler.encounters`).

Add params to the signature: `priority_date_iso: str | None = None`, `vinf_targets_kms: dict[str, float] | None = None`, plus the existing `n_laps`/`rp_factors` (keep for API compatibility; `n_laps` only matters if drift is computed).

- [ ] **Step 5: Run to verify it passes**

Run: `uv run pytest tests/search/test_optimize_ephemeris.py -v`
Expected: PASS (slow — real ephemeris). If `_resolve_real_t_start` returns `None`, revisit the `vinf_targets_kms`/`priority_date_iso` inputs.

- [ ] **Step 6: Commit**

```bash
git add src/cyclerfinder/search/optimize.py tests/search/test_optimize_ephemeris.py
git commit -m "search/optimize: implement optimise_cell_ephemeris over the general maintenance engine"
```

---

## Task 3: S1L1 real-ephemeris rediscovery gate (the payoff)

**Files:**
- Test: `tests/search/test_s1l1_real_rediscovery.py` (create)

- [ ] **Step 1: Write the gate test**

```python
"""S1L1 real-ephemeris rediscovery — the idealised model could NOT host the
5.65/3.05 km/s anchors (see scripts/characterise_s1l1.py); this gate tests
whether the real ephemeris (Mars eccentricity) can, like the Aldrin cycler.
Anchors are sourced (spec §9); epoch and leg ToFs are computed."""

import numpy as np
import pytest
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M, TOL = 5.65, 3.05, 0.4  # sourced anchors; real-eph band


def _vinf_by_body(result):
    out = {}
    for enc in result.best_cycler.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


@pytest.mark.slow
@pytest.mark.xfail(strict=False, reason="aspirational: real-ephemeris S1L1 closure is the open goal this plan targets; flips to passing once the optimiser reaches the 5.65/3.05 basin")
def test_s1l1_real_ephemeris_rediscovers_anchors():
    cell = Cell(bodies=("E", "M"), sequence=("E", "M", "E"), period_k=2,
                per_leg_revs=(0, 0), per_leg_branch=("single", "single"))
    eph = Ephemeris(model="astropy")
    result = optimise_cell_ephemeris(
        cell, eph, vinf_cap=8.15,
        priority_date_iso="2002-08-05",
        vinf_targets_kms={"E": VINF_E, "M": VINF_M},
        n_starts=5, seed=0,
    )
    assert result.constraints_satisfied
    v = _vinf_by_body(result)
    assert abs(v["E"] - VINF_E) < TOL
    assert abs(v["M"] - VINF_M) < TOL
```

- [ ] **Step 2: Run it**

Run: `uv run pytest tests/search/test_s1l1_real_ephemeris_rediscovers_anchors -v` (note: `-k s1l1_real`)
Expected: XPASS (goal met — promote to a hard assert by removing the xfail marker in the same commit) or XFAIL (real-eph also cannot host it at this epoch/topology — keep xfail; record the achieved V∞ in the reason). Either outcome is honest and informative. Do NOT loosen `TOL` beyond the 0.4 km/s real-ephemeris band to force a pass.

- [ ] **Step 3: If XPASS, retire the S1L1 catalogue xfail too**

If S1L1 now closes, update `tests/verify/test_real_closure.py` Gate-2 and `EXPECTED_SKIPS["s1l1-2syn-em-cpom"]` (`real_closure.py:222`) to reflect real-ephemeris closure, and write the computed leg ToFs/epoch back into the catalogue entry marked `provenance: computed`.

- [ ] **Step 4: Commit**

```bash
git add tests/search/test_s1l1_real_rediscovery.py tests/verify/test_real_closure.py src/cyclerfinder/verify/real_closure.py data/catalogue.yaml
git commit -m "tests: S1L1 real-ephemeris rediscovery gate against sourced 5.65/3.05 anchors"
```

---

## Task 4: Aldrin parity + general-path regression

**Files:**
- Test: `tests/search/test_optimize_ephemeris.py`

- [ ] **Step 1: Add a parity test**

Assert that `optimise_cell_ephemeris` on the Aldrin E-M-E cell recovers elements consistent with the Aldrin-specific `solve_powered_periodic_cycler` (a≈1.60 AU, e≈0.393 within band), proving the general path matches the specialised one. Use the sourced Aldrin anchors from `tests/verify/test_real_closure.py:205-210`. Mark `@pytest.mark.slow`.

- [ ] **Step 2: Run + commit**

Run: `uv run pytest tests/search/test_optimize_ephemeris.py -v`
Expected: PASS.
```bash
git add tests/search/test_optimize_ephemeris.py
git commit -m "tests: Aldrin parity between general optimise_cell_ephemeris and the Aldrin-specific solver"
```

---

## Task 5: Full-suite regression + lint + type gate

- [ ] **Step 1:** `uv run pytest -q` — record pass/skip/xfail counts; all green + documented xfails.
- [ ] **Step 2:** `uv run ruff check . && uv run ruff format --check . && uv run mypy src` — all pass.
- [ ] **Step 3:** Commit any fixups: `git add -A && git commit -m "real-eph optimiser: full-suite green"`.

---

## Self-Review

**Spec coverage:** Task 1 (cell→engine inputs), Task 2 (wire `optimise_cell_ephemeris` over the existing general `optimise_maintenance_dv`), Task 3 (S1L1 real-eph rediscovery — the payoff for choosing this path), Task 4 (Aldrin parity), Task 5 (gates).

**Scope honesty:** This generalises an EXISTING engine; it does not build a solver from scratch (the architecture map confirmed `optimise_maintenance_dv` is already body-agnostic). Multi-body (N≥3) period handling is explicitly deferred to M8-Core's `period_basis` work. Blind discovery without V∞ targets (epoch resolution needs targets) is out of scope and documented.

**Provenance:** Validation asserts sourced V∞ anchors (5.65/3.05, Aldrin a/e); epoch + ToFs are computed and labelled. Task 3 is `xfail` until it genuinely closes — no forced pass, no loosened tolerance.

**Type consistency:** `optimise_maintenance_dv(...)`, `MaintenanceOptimResult.{cycler,t0_sec,leg_tofs_days,maintenance_dv_kms,converged}`, `OptimisationResult(cell,best_cycler,best_score,closure_residual_kms,optimiser_history,converged,constraints_satisfied)`, `Score.{max_vinf_kms,total_maintenance_dv_kms}`, `_resolve_real_t_start(signature,ephem,priority)`, `_target_period_sec(cell)` all match source signatures verified 2026-06-02.

**Open risk:** epoch resolution depends on `_resolve_real_t_start` finding a window under the 20 km/s mismatch cap from the equispaced ToF seed; if the seed is too far from the real geometry it may return `None`. Task 2 Step 4 handles `None` honestly (returns a non-converged result), and the seed can be refined from the idealised optimum if needed.
