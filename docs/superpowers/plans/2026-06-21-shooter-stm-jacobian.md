# Analytic STM Jacobian for the n-body shooter — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `nbody.shooter.shoot` an STM (state-transition-matrix) Jacobian from per-leg REBOUND co-integrated variational propagation, replacing the `1+6n`-re-propagation finite-difference Jacobian, so the full multiple-shooting corrector becomes tractable on the multi-rev SnLm cyclers — then run the literal-Russell-parent verdict.

**Architecture:** Add `with_stm=True` to `RestrictedNBody.propagate` (REBOUND `add_variation`); assemble the block-bidiagonal multiple-shooting Jacobian from per-leg STMs as an opt-in `jacobian="stm"` in `shoot` (FD stays default + oracle); validate against FD + complex-step + analytic-Keplerian; run a detached, checkpointed, days-tolerant literal-parent batch.

**Tech Stack:** Python 3, REBOUND 5.0.0 (native variational particles), numpy, scipy, pytest, uv. `from __future__ import annotations`, full hints, ruff+mypy clean.

**Spec:** `docs/superpowers/specs/2026-06-21-shooter-stm-jacobian-design.md`

**Reused APIs (read before coding):**
- `nbody/propagator.py::RestrictedNBody.propagate(r0_km, v0_km_s, t0_sec, t1_sec, *, bodies=(), accuracy=1e-10, ephem=None, cache=None, max_steps, max_wall_sec) -> NBodyArc` (read `NBodyArc` dataclass ~line 41; perturbers ride on rails via `_install_rails_forces` custom `additional_forces`). REBOUND import is inside the method.
- `core/cr3bp.py` — `propagate(with_stm=True)` + `cr3bp_stm_eom` (the project's variational-STM precedent to mirror).
- `nbody/shooter.py::shoot` (the LM loop, `residual_of_x`, `_x_to_states`/`_states_to_x`, `_fd_jacobian`, `defect_residual`). Read `defect_residual` (~line 199) for the EXACT residual layout (the block Jacobian must match it row-for-row).
- `search/shooter_russell_seed.py::russell_shooting_seed`, `scripts/shooter_russell_batch.py` (Phase 3 reuses these).

**REBOUND variational note (load-bearing):** `sim.add_variation(order=1)` adds a variational particle set; seed its spacecraft particle with a unit basis perturbation, co-integrate, read the variation's final spacecraft state = one column of Φ. Six variations (one per basis vector) co-integrated in ONE `sim.integrate` give the full 6×6 Φ. **CRITICAL RISK:** the rails perturbers are custom `additional_forces`; REBOUND's built-in variational equations cover in-sim gravity but may NOT include the custom-force gradient. So Φ is exact for Sun-only (`bodies=()`) but may be wrong with perturbers — Task 2's parity test MUST cover the perturbed case and is the go/no-go gate.

---

## Task 1: STM propagation (`propagate(with_stm=True)`)

**Files:** Modify `src/cyclerfinder/nbody/propagator.py`; Test `tests/nbody/test_propagator_stm.py`

- [ ] **Step 1: Write the failing test (Sun-only Φ vs analytic Keplerian + FD)**
```python
from __future__ import annotations

import numpy as np

from cyclerfinder.nbody.propagator import RestrictedNBody


def _fd_stm(prop, r0, v0, t0, t1, *, h=1.0) -> np.ndarray:
    x0 = np.concatenate([r0, v0])
    base = prop.propagate(x0[:3], x0[3:], t0, t1)
    f0 = np.concatenate([base.r_final_km, base.v_final_km_s])
    jac = np.empty((6, 6))
    for j in range(6):
        xp = x0.copy(); xp[j] += h
        a = prop.propagate(xp[:3], xp[3:], t0, t1)
        fj = np.concatenate([a.r_final_km, a.v_final_km_s])
        jac[:, j] = (fj - f0) / h
    return jac


def test_sun_only_stm_matches_fd() -> None:
    prop = RestrictedNBody()
    # heliocentric ~1 AU circular-ish leg, Sun-only (bodies=())
    r0 = np.array([1.496e8, 0.0, 0.0]); v0 = np.array([0.0, 29.78, 0.0])
    t0, t1 = 0.0, 120.0 * 86400.0
    arc = prop.propagate(r0, v0, t0, t1, with_stm=True)
    assert arc.stm is not None and arc.stm.shape == (6, 6)
    fd = _fd_stm(prop, r0, v0, t0, t1, h=1.0)
    rel = np.linalg.norm(arc.stm - fd) / np.linalg.norm(fd)
    assert rel < 1e-3, f"Sun-only STM vs FD rel={rel}"


def test_with_stm_false_is_unchanged() -> None:
    prop = RestrictedNBody()
    r0 = np.array([1.496e8, 0.0, 0.0]); v0 = np.array([0.0, 29.78, 0.0])
    a = prop.propagate(r0, v0, 0.0, 1e7)
    assert getattr(a, "stm", None) is None  # default path: no STM
```
(Read `NBodyArc`'s actual final-state field names — the test uses `r_final_km`/`v_final_km_s`; ADAPT to the real names.)

- [ ] **Step 2: Run → FAIL** `uv run pytest tests/nbody/test_propagator_stm.py -v` (no `with_stm` / no `.stm`).
- [ ] **Step 3: Implement** — add `with_stm: bool = False` to `propagate`. Add an optional `stm: np.ndarray | None = None` field to `NBodyArc` (default None → existing callers unchanged). When `with_stm`: after building the sim + spacecraft particle, add 6 `sim.add_variation(order=1)` sets, seed each with one unit basis perturbation on the spacecraft particle (3 position, 3 velocity), `sim.integrate(t1)`, read each variation's final spacecraft `[x,y,z,vx,vy,vz]` as a column of the 6×6 `stm`; attach to the returned `NBodyArc`. Keep `with_stm=False` byte-identical (no variations added).
- [ ] **Step 4: Run → PASS** (both tests); ruff + mypy clean.
- [ ] **Step 5: Commit** `git add -A && git commit -m "nbody/#388: STM propagation via REBOUND co-integrated variational particles"`

---

## Task 2: STM parity gate — perturbed case + complex-step (GO/NO-GO)

**Files:** Test `tests/nbody/test_propagator_stm.py` (append); possibly Modify `propagator.py` if the perturbed gradient must be added

- [ ] **Step 1: Append the perturbed-parity + CSD tests**
```python
import pytest

from cyclerfinder.core.ephemeris import Ephemeris


@pytest.mark.slow
def test_perturbed_stm_matches_fd() -> None:
    prop = RestrictedNBody()
    ephem = Ephemeris("astropy")
    r0 = np.array([1.496e8, 0.0, 0.0]); v0 = np.array([0.0, 29.78, 0.0])
    t0, t1 = 0.0, 200.0 * 86400.0
    kw = dict(bodies=("E", "M"), ephem=ephem)
    arc = prop.propagate(r0, v0, t0, t1, with_stm=True, **kw)
    x0 = np.concatenate([r0, v0]); f0v = []
    base = prop.propagate(r0, v0, t0, t1, **kw)
    f0 = np.concatenate([base.r_final_km, base.v_final_km_s])
    fd = np.empty((6, 6))
    for j in range(6):
        xp = x0.copy(); xp[j] += 1.0
        a = prop.propagate(xp[:3], xp[3:], t0, t1, **kw)
        fd[:, j] = (np.concatenate([a.r_final_km, a.v_final_km_s]) - f0) / 1.0
    rel = np.linalg.norm(arc.stm - fd) / np.linalg.norm(fd)
    assert rel < 5e-3, f"PERTURBED STM vs FD rel={rel} — does REBOUND variation include the rails additional_forces gradient?"
```
- [ ] **Step 2: Run** `uv run pytest tests/nbody/test_propagator_stm.py -m slow -v`
- [ ] **Step 3: GO/NO-GO.** If the perturbed parity PASSES (rel < 5e-3), the REBOUND variation includes (or adequately approximates) the perturber gradient — proceed; commit. **If it FAILS**, REBOUND's variational equations omit the custom rails-force gradient — STOP and report BLOCKED with the actual rel error and the per-column discrepancy. Do NOT hack the tolerance. (Resolution, decided then: either add the analytic perturber-gradient term to the variational seeding, or accept FD for perturbed legs / Sun-only-STM as a preconditioner — a design decision for the user, not a silent choice.)
- [ ] **Step 4: (if GO) Commit** `git add -A && git commit -m "nbody/#388: STM perturbed-parity gate (REBOUND variation vs FD with perturbers)"`

---

## Task 3: STM-assembled multiple-shooting Jacobian (`jacobian="stm"`)

**Files:** Modify `src/cyclerfinder/nbody/shooter.py`; Test `tests/nbody/test_shooter_stm_jacobian.py`

- [ ] **Step 1: Write the failing parity test (STM Jacobian vs FD Jacobian)**
```python
from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.slow
def test_stm_jacobian_matches_fd_jacobian() -> None:
    # Build a small ShootingSeed (reuse an existing shooter test fixture / a 2-node
    # near-circular seed). Read tests/nbody/test_shooter_*.py for a ready fixture.
    from cyclerfinder.nbody import shooter
    seed, ephem, bodies = _small_shooting_fixture()  # define from existing tests
    x0 = shooter._states_to_x(seed.node_states)
    f0 = _residual(seed, x0, ephem, bodies)
    fd = shooter._fd_jacobian(lambda x: _residual(seed, x, ephem, bodies), x0, f0,
                              column_eval=shooter._serial_columns)
    stm_jac = shooter._stm_jacobian(seed, x0, ephem=ephem, bodies=bodies)
    rel = np.linalg.norm(stm_jac - fd) / np.linalg.norm(fd)
    assert rel < 5e-3, f"STM Jacobian vs FD rel={rel}"
```
(Define `_small_shooting_fixture` and `_residual` from an existing `tests/nbody/test_shooter_*.py` fixture — read those tests for a cheap 2–3 node seed.)
- [ ] **Step 2: Run → FAIL** (no `_stm_jacobian`).
- [ ] **Step 3: Implement `_stm_jacobian(seed, x, *, ephem, bodies, accuracy, max_wall_sec)`** in shooter.py: read `defect_residual`'s exact layout, then for each leg propagate with `with_stm=True` to get `Φ_i`; assemble the block-bidiagonal Jacobian — for each full-state continuity defect `c_i = leg_i(node_i) − node_{i+1}`: `∂c_i/∂node_i = Φ_i`, `∂c_i/∂node_{i+1} = −I_6`; fill boundary/periodicity rows to match `defect_residual` exactly. Add `jacobian: Literal["fd","stm"] = "fd"` to `shoot`; when `"stm"`, pass `jac=lambda x: _stm_jacobian(...)` to `least_squares`. `"fd"` path byte-identical.
- [ ] **Step 4: Run → PASS** parity; add a second test that `shoot(jacobian="stm")` reaches the same fixed point as `shoot(jacobian="fd")` on the cheap fixture (defect norms within tol) — `@pytest.mark.slow`. ruff + mypy clean.
- [ ] **Step 5: Commit** `git add -A && git commit -m "nbody/#388: block-bidiagonal STM Jacobian for shoot (jacobian='stm', FD stays oracle)"`

---

## Task 4: Literal-parent batch (detached, checkpointed) + results note

**Files:** Modify `scripts/shooter_russell_batch.py`; Create `docs/notes/2026-06-21-shooter-stm-batch-results.md`

- [ ] **Step 1: Wire `jacobian="stm"` + incremental checkpointing into the batch.** Modify `scripts/shooter_russell_batch.py`: pass `jacobian="stm"` to `shooter.shoot`; **write each (row, epoch) record to the JSONL runlog immediately after it completes** (open in append mode, flush per record) so a kill/restart loses nothing. Do NOT cap K for monitorability (per `feedback_long_runs_acceptable`) — run the best phase-error epochs over LaunchWindow 1..21 to completion. Add a `--resume` guard that skips (row,epoch) pairs already in the runlog.
- [ ] **Step 2: Launch detached + checkpointed.** `setsid nohup uv run python scripts/shooter_russell_batch.py > data/runs/shooter-stm-batch.log 2>&1 &` (survives agent reaping). Record the PID + runlog path. Do NOT block on it. Report progress by reading the JSONL runlog / `.log` tail, however many days it takes. (Per the spec, this is the days-tolerant run.)
- [ ] **Step 3: Write `docs/notes/2026-06-21-shooter-stm-batch-results.md`** as results arrive: per-row table (converged / defect / emerged n-body V∞ vs sourced anchor E&M / anchor-match / bend-feasible / winning epoch / wall), the STM-vs-FD speedup actually observed, and the verdict. If any row (esp. `mcconaghy-2006-em-k2` or V3 regression rows) closes + anchor-matches → PROPOSED V0→V1, HELD, no writeback. If none → the decisive characterized negative: the full multiple-shooting corrector, seeded from the literal Russell parent and now run to convergence, does not land the sourced family. No catalogue writeback either way.
- [ ] **Step 4: Commit** (script + note + whatever runlog exists) `git add scripts/shooter_russell_batch.py docs/notes/2026-06-21-shooter-stm-batch-results.md data/runs/ && git commit -m "search/#388: STM-Jacobian literal-parent batch (detached, checkpointed) + results note"`

---

## Final verification
- [ ] `uv run pytest tests/nbody/test_propagator_stm.py tests/nbody/test_shooter_stm_jacobian.py -q` (incl `-m slow`) — PASS (or Task 2 BLOCKED reported, halting before Phase 3).
- [ ] `uv run ruff check . && uv run ruff format --check .` — clean.
- [ ] Confirm NO `data/catalogue.yaml` / `validate.py` edit.
- [ ] Report: the STM parity results (Task 2 GO/NO-GO), the Jacobian speedup, and — once the detached batch produces rows — the literal-parent verdict.
