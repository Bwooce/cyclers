# Seed the full shooter from the Russell parent — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Feed the existing `nbody/shooter.py` full-state multiple-shooting corrector the constructed Russell ψ generic-return parent (via `correct._vinf_nodes` → `seed_from_conic`) and run `shoot()` on the real n-body model for the catalogue SnLm rows — the decisive seeding test of the #135/S1L1 family-selection wall.

**Architecture:** One small new module `search/shooter_russell_seed.py` (a seed adapter) + a batch script. Pure reuse of `narc_continuation`, `correct._vinf_nodes`, `nbody/shooter` (`seed_from_conic`, `shoot`), `cycler_assembly`.

**Tech Stack:** Python 3, numpy, scipy, pytest, uv. `from __future__ import annotations`, full hints, ruff+mypy clean.

**Spec:** `docs/superpowers/specs/2026-06-21-shooter-russell-seed-design.md`

**Reused APIs (read before coding):**
- `search/narc_continuation.py`: `russell_parent_to_ballistic_seed(model, cycler, row) -> NarcSeed` (`.sequence,.per_leg_revs,.per_leg_branch,.tof_seed_days,.period_sec,.vinf_anchor_e_kms,.vinf_anchor_m_kms`); `candidate_epochs(ephem, target_phase, *, launch_window_synodics, grid)`.
- `search/correct.py`: `_vinf_nodes(*, sequence, per_leg_revs, per_leg_branch, t0_sec, free_tof_days, slack_leg, period_days, ephem, mu_central=MU_SUN_KM3_S2) -> dict[str,np.ndarray]` (keys `b{i}_in`/`b{i}_out`). `DAY_S` is defined in correct.py.
- `nbody/shooter.py`: `seed_from_conic(*, sequence, vinf_nodes, t0_sec, tofs_days, slack_leg, period_days, ephem) -> ShootingSeed`; `shoot(...) -> ShootResult` (READ its exact signature + `ShootResult` fields: converged flag, defect/residual, per-node/per-encounter emerged V∞; use the actual names).
- `search/cycler_assembly.py`: `descriptor_to_phsi`, `assemble_cycler`; `search/generic_return.py`: `RussellModel`; `cyclerfinder.core.ephemeris.Ephemeris`; `cyclerfinder.data.catalog.load_catalog`.

---

## Task 1: Seed adapter `russell_shooting_seed`

**Files:** Create `src/cyclerfinder/search/shooter_russell_seed.py`; Test `tests/search/test_shooter_russell_seed.py`

- [ ] **Step 1: Write the failing test**
```python
from __future__ import annotations

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import russell_parent_to_ballistic_seed
from cyclerfinder.search.shooter_russell_seed import russell_shooting_seed


def test_adapter_builds_shooting_seed() -> None:
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    phsi = descriptor_to_phsi(row)
    assert phsi is not None
    cyc = assemble_cycler(m, phsi)
    assert cyc is not None
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    ephem = Ephemeris("astropy")
    sseed = russell_shooting_seed(seed, t0_sec=0.0, ephem=ephem)
    n = len(seed.sequence)
    assert len(sseed.node_states) == n          # one (6,) state per encounter
    assert all(s.shape == (6,) for s in sseed.node_states)
    assert len(sseed.tofs) == n - 1
    assert sseed.sequence == seed.sequence
```

- [ ] **Step 2: Run → FAIL** `uv run pytest tests/search/test_shooter_russell_seed.py -v`

- [ ] **Step 3: Implement** `src/cyclerfinder/search/shooter_russell_seed.py`:
```python
from __future__ import annotations

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody import shooter
from cyclerfinder.search import correct
from cyclerfinder.search.narc_continuation import NarcSeed


def russell_shooting_seed(seed: NarcSeed, t0_sec: float, ephem: Ephemeris) -> shooter.ShootingSeed:
    """Bridge a Russell parent NarcSeed -> a multiple-shooting ShootingSeed.

    Pins the longest leg as the slack leg (the ballistic_correct convention),
    builds the per-encounter v_inf VECTOR nodes via ``correct._vinf_nodes`` (one
    Lambert per leg), then maps them to a ShootingSeed via ``shooter.seed_from_conic``.
    Raises ValueError if the Lambert node solve is degenerate (caught by the batch).
    """
    tofs = list(seed.tof_seed_days)
    slack_leg = int(np.argmax(tofs))
    free_tofs = [t for i, t in enumerate(tofs) if i != slack_leg]
    period_days = seed.period_sec / correct.DAY_S
    nodes = correct._vinf_nodes(
        sequence=seed.sequence,
        per_leg_revs=seed.per_leg_revs,
        per_leg_branch=seed.per_leg_branch,
        t0_sec=t0_sec,
        free_tof_days=free_tofs,
        slack_leg=slack_leg,
        period_days=period_days,
        ephem=ephem,
    )
    return shooter.seed_from_conic(
        sequence=seed.sequence,
        vinf_nodes=nodes,
        t0_sec=t0_sec,
        tofs_days=tofs,
        slack_leg=slack_leg,
        period_days=period_days,
        ephem=ephem,
    )
```
(If `_vinf_nodes` expects the FULL tof list rather than free-only, or `seed_from_conic` wants different arg names, ADAPT to the real signatures you read — the intent is fixed: parent → v∞-vector nodes → ShootingSeed. Note: importing the private `correct._vinf_nodes` is acceptable intra-package reuse; if ruff flags it, add a targeted noqa as the sibling modules do.)

- [ ] **Step 4: Run → PASS**; `uv run ruff check src/cyclerfinder/search/shooter_russell_seed.py tests/search/test_shooter_russell_seed.py && uv run ruff format --check ...` clean; mypy clean.
- [ ] **Step 5: Commit** `git add -A && git commit -m "search/#388: seed adapter Russell parent -> shooter ShootingSeed"`

---

## Task 2: Batch + results note

**Files:** Create `scripts/shooter_russell_batch.py`; Create `docs/notes/2026-06-21-shooter-russell-results.md`; (optional smoke test in `tests/search/test_shooter_russell_seed.py`)

- [ ] **Step 1: Optional shoot smoke test** (APPEND to the test file) — confirms `shoot()` runs from the adapter's seed without error (converged or not):
```python
def test_shoot_runs_from_russell_seed() -> None:
    import pytest

    from cyclerfinder.nbody import shooter
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    seed = russell_parent_to_ballistic_seed(m, assemble_cycler(m, descriptor_to_phsi(row)), row)
    ephem = Ephemeris("astropy")
    sseed = russell_shooting_seed(seed, t0_sec=0.0, ephem=ephem)
    # READ shooter.shoot's real signature; call it minimally. Mark slow if needed.
    res = shooter.shoot(sseed, ephem=ephem)   # adapt kwargs to the real signature
    assert res is not None
```
Mark `@pytest.mark.slow` if `shoot()` is heavy. If the smoke is too slow/finicky for CI, drop it and rely on the batch — note that choice.

- [ ] **Step 2: Write `scripts/shooter_russell_batch.py`** — `main()`:
  - `m = RussellModel()`; `ephem = Ephemeris("astropy")`; `import time`.
  - For each `e` in `load_catalog().entries`: `phsi = descriptor_to_phsi(e.raw)`; skip None. `cyc = assemble_cycler(m, phsi)`; skip None. `try: seed = russell_parent_to_ballistic_seed(m, cyc, e.raw) except ValueError: continue`.
  - `epochs = candidate_epochs(ephem, 0.0, launch_window_synodics=range(1, 22), grid=100)[:K]` with `K = 3` (DOCUMENT the cap — full n-body shoot per epoch is expensive; only the best 3 phase-error epochs are run, logged as a coverage cap).
  - For each epoch: `try: sseed = russell_shooting_seed(seed, t0_sec=epoch, ephem=ephem); res = shooter.shoot(sseed, ephem=ephem, ...) except Exception: continue`. Keep the lowest-defect `ShootResult`.
  - Anchor compare: emerged per-encounter n-body V∞ (from `ShootResult`) min-distance to `seed.vinf_anchor_e_kms` / `_m_kms`; `anchor_match = bothE&M <= 0.5`. `promote = e.id=="mcconaghy-2006-em-k2" and converged and anchor_match and bend_feasible` (use the ShootResult's real converged/bend fields).
  - Print per-row line + append record. Write `data/runs/shooter-russell-{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}.jsonl`. Summary; "HELD — no writeback". `if __name__=="__main__": main()`.
  - Run: `timeout 3000 uv run python scripts/shooter_russell_batch.py 2>&1 | tail -40`. If a single row's shoot exceeds ~5 min or the whole thing the timeout, reduce K to 1 and/or `launch_window_synodics` to `range(1,6)`, document the reduction, re-run.

- [ ] **Step 3: Write `docs/notes/2026-06-21-shooter-russell-results.md`** factually: per-row table (converged / defect / emerged n-body V∞ vs sourced anchor E&M / anchor-match / bend-feasible / winning epoch), the K/launch-window cap used, and the verdict. If any row (esp. mcconaghy-2006-em-k2 or the V3 regression rows) closes + anchor-matches: PROPOSED V0→V1, HELD, no writeback. If none: state plainly — feeding the full shooter the *constructed Russell parent* still lands off-anchor → the seeding/family-selection wall is confirmed even from the literal parent (the decisive close of the #388/#135/S1L1 thread). Compare emerged V∞ to the narc-continuation low-energy basin.

- [ ] **Step 4: Commit** `git add scripts/shooter_russell_batch.py docs/notes/2026-06-21-shooter-russell-results.md data/runs/ src/ tests/ && git commit -m "search/#388: shooter-from-Russell-parent batch + results note (held, no writeback)"`

---

## Final verification
- [ ] `uv run pytest tests/search/test_shooter_russell_seed.py -q` — PASS.
- [ ] `uv run ruff check . && uv run ruff format --check .` — clean.
- [ ] Confirm NO `data/catalogue.yaml` / `validate.py` edit.
- [ ] Report the decisive verdict: does the literal Russell parent close on the n-body model (held V0→V1) or confirm the seeding wall.
