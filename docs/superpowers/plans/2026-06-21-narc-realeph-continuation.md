# N-arc Real-Ephemeris Continuation Lane Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Continue the validated Russell idealized generic-return parents to the real DE440 ephemeris via the Russell §5.4 homotopy + N-arc corrector, seeded by the new generator, and attempt held V0→V1 closure on the catalogue SnLm rows.

**Architecture:** One new module `search/narc_continuation.py` (bridge + phase-angle epoch + N-arc ramp driver) + a batch script, reusing `search/continuation.py::ramped_ephemeris`/`LADDER`, `search/correct.py::ballistic_correct`, and the new `search/generic_return.py` / `search/cycler_assembly.py`.

**Tech Stack:** Python 3, numpy, scipy, pytest, uv. `from __future__ import annotations`, full type hints, ruff+mypy clean (pre-commit enforces).

**Spec:** `docs/superpowers/specs/2026-06-21-narc-realeph-continuation-design.md`

**Key reused APIs (read before coding):**
- `search/correct.py::ballistic_correct(sequence, per_leg_revs, per_leg_branch, t0_seed_sec, tof_seed_days, period_sec, ephem, *, vinf_cap, residual_mode="vector", tol_kms=0.1) -> BallisticClosureResult` (has `.constraints_satisfied`, `.max_residual_kms`/residual, per-encounter V∞; read the dataclass at correct.py:39). `per_leg_branch` uses LAMBERT labels (`"single"/"low"/"high"`), via `_pick`.
- `search/continuation.py::ramped_ephemeris(lam_e, lam_i, lam_p=0.0) -> Ephemeris` and `LADDER`.
- `search/cycler_assembly.py::descriptor_to_phsi(row) -> PhsiSpec|None`, `assemble_cycler(model, phsi) -> Cycler|None`.
- `search/generic_return.py::RussellModel` (`.tu_days`=58.1324409), `Cycler` (fields: `p,h,s,i,generic_return,turn_angles,ar,tr,vinf_e,vinf_m`; `generic_return` has `n_revs`, `branch` in `"slow"/"fast"`).
- `search/dsm_descriptor_seed.py::seed_dsm_chain_from_descriptor(row)` — its `.per_leg_tof_days` is the SOURCED per-leg ToF list (days) we reuse as the bridge's ToF seed.

**Constants:** real Earth-Mars synodic = `1/(1/1.0 - 1/1.8808)` yr ≈ 2.1354 yr (real Mars sidereal 1.8808 yr); `YEAR_DAYS=365.25`, `DAY_S=86400`.

---

## Task 1: The bridge — Russell parent → N-arc real-eph seed

**Files:** Create `src/cyclerfinder/search/narc_continuation.py`; Test `tests/search/test_narc_continuation.py`

- [ ] **Step 1: Write the failing test**
```python
from __future__ import annotations

from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import NarcSeed, russell_parent_to_ballistic_seed


def test_bridge_builds_real_seed() -> None:
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    phsi = descriptor_to_phsi(row)
    assert phsi is not None
    cyc = assemble_cycler(m, phsi)
    assert cyc is not None
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    assert isinstance(seed, NarcSeed)
    assert seed.sequence == ("E", "E", "M", "M")
    assert len(seed.per_leg_revs) == len(seed.sequence) - 1
    assert len(seed.per_leg_branch) == len(seed.sequence) - 1
    assert len(seed.tof_seed_days) == len(seed.sequence) - 1
    # sourced per-leg ToFs (the 4.991gG2 published values), days
    assert abs(seed.tof_seed_days[0] - 533.70) < 1.0
    assert abs(seed.tof_seed_days[1] - 150.0) < 1.0
    assert abs(seed.tof_seed_days[2] - 1026.21) < 1.0
    # period_sec from REAL synodic * p (p=2 for this two-synodic row)
    real_syn_yr = 1.0 / (1.0 / 1.0 - 1.0 / 1.8808)
    assert abs(seed.period_sec - phsi.p * real_syn_yr * 365.25 * 86400.0) < 1.0e7
    assert seed.vinf_anchor_e_kms > 0 and seed.vinf_anchor_m_kms > 0
    # lambert branch labels, not slow/fast
    assert all(b in ("single", "low", "high") for b in seed.per_leg_branch)
```

- [ ] **Step 2: Run → FAIL** `uv run pytest tests/search/test_narc_continuation.py -v` (ImportError).

- [ ] **Step 3: Implement** `NarcSeed` (frozen dataclass: `sequence: tuple[str,...]`, `per_leg_revs: tuple[int,...]`, `per_leg_branch: tuple[str,...]`, `tof_seed_days: tuple[float,...]`, `period_sec: float`, `vinf_anchor_e_kms: float`, `vinf_anchor_m_kms: float`) and `russell_parent_to_ballistic_seed(model, cycler, row) -> NarcSeed`:
  - `sequence = tuple(row["sequence_canonical"].split("-"))`.
  - `tof_seed_days`: from `seed_dsm_chain_from_descriptor(row).per_leg_tof_days` (the sourced published ToFs). If that returns None, raise `ValueError` (caller skips).
  - `per_leg_revs` / `per_leg_branch` per leg `i` (`a,b = sequence[i], sequence[i+1]`):
    - transit leg (`a != b`): `(0, "single")`.
    - same-body resonant leg (`a == b`): revs = `cycler.generic_return.n_revs`; branch = map the Russell `"slow"/"fast"` to a lambert label — `"low" if cycler.generic_return.branch == "slow" else "high"` (provisional; the driver smoke + batch will confirm). For a 0-rev resonant leg use `"single"`.
  - `real_syn_yr = 1.0/(1.0/1.0 - 1.0/1.8808)`; `period_sec = phsi.p * real_syn_yr * 365.25 * 86400.0` — get `p` from `descriptor_to_phsi(row)` (call it, or accept `phsi` is recoverable from `cycler.p`; use `cycler.p`).
  - `vinf_anchor_e_kms`/`vinf_anchor_m_kms` from `row["vinf_kms_at_encounters"]` (E and M entries).

- [ ] **Step 4: Run → PASS**; `uv run ruff check src/cyclerfinder/search/narc_continuation.py tests/search/test_narc_continuation.py && uv run ruff format --check ...` clean.
- [ ] **Step 5: Commit** `git add -A && git commit -m "search/#388: bridge Russell parent -> N-arc real-eph seed"`

---

## Task 2: Phase-angle epoch derivation (Russell §5.3)

**Files:** `narc_continuation.py`; test same file

- [ ] **Step 1: Write the failing test**
```python
def test_candidate_epochs_match_phase() -> None:
    import numpy as np

    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.narc_continuation import candidate_epochs, parent_phase_angle

    ephem = Ephemeris("astropy")
    target_phase = 1.0  # rad, the parent's beginning E-M relative phase
    epochs = candidate_epochs(ephem, target_phase, launch_window_synodics=range(1, 6), grid=50)
    assert epochs  # at least one epoch found in the window
    for t0 in epochs:
        assert abs(parent_phase_angle(ephem, t0) - target_phase) < 0.15  # within ~9 deg


def test_parent_phase_angle_range() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.narc_continuation import parent_phase_angle

    ph = parent_phase_angle(Ephemeris("astropy"), 0.0)
    assert -3.1416 <= ph <= 3.1416
```

- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement**:
  - `parent_phase_angle(ephem, t0_sec) -> float`: the Earth-Mars relative phase at epoch — `r_E, _ = ephem.state("E", t0)`; `r_M, _ = ephem.state("M", t0)`; return the signed in-ecliptic angle from Earth's heliocentric direction to Mars's (`atan2` of the cross/dot of the xy projections), in `(-pi, pi]`.
  - `candidate_epochs(ephem, target_phase, *, launch_window_synodics=range(1, 22), grid=100) -> list[float]`: real synodic `T_syn = 2.1354 yr` in seconds; for each `w` in `launch_window_synodics`, scan `grid` equally-spaced epochs across `[w*T_syn, (w+1)*T_syn]` after J2000 (t=0), evaluate `parent_phase_angle`, and pick the epoch minimising `|wrap(parent_phase_angle - target_phase)|` in that window (one per window). Return the list sorted by that phase error (best-first). This is Russell's LaunchWindow 1..21 + 100-interval grid.

- [ ] **Step 4: Run → PASS**; ruff + mypy clean.
- [ ] **Step 5: Commit** `git add -A && git commit -m "search/#388: phase-angle epoch derivation (Russell SS5.3 LaunchWindow + 100-grid)"`

---

## Task 3: N-arc continuation driver (ramp + ballistic_correct)

**Files:** `narc_continuation.py`; test same file

- [ ] **Step 1: Write the failing test** (smoke: runs end-to-end at the circular rung)
```python
def test_narc_driver_circular_rung_smoke() -> None:
    import cyclerfinder.search.continuation as cont
    from cyclerfinder.data.catalog import load_catalog
    from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
    from cyclerfinder.search.generic_return import RussellModel
    from cyclerfinder.search.narc_continuation import (
        narc_continuation_correct,
        russell_parent_to_ballistic_seed,
    )

    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    cyc = assemble_cycler(m, descriptor_to_phsi(row))
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    # circular ramped backend (lam_e=lam_i=0) as the final ephemeris isolates plumbing
    circ = cont.ramped_ephemeris(0.0, 0.0)
    res = narc_continuation_correct(
        seed, ladder=(1,), final_ephemeris=circ, epochs=[0.0], vinf_cap=99.0
    )
    assert res is not None
    assert hasattr(res, "converged")
    assert hasattr(res, "max_residual_kms")
    assert hasattr(res, "emerged_vinf_kms")
```

- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `NarcContinuationResult` (frozen: `converged: bool`, `max_residual_kms: float`, `emerged_vinf_kms: tuple[float,...]`, `t0_sec: float`, `tof_days: tuple[float,...]`, `bend_feasible: bool`, `winning_epoch_sec: float`) and `narc_continuation_correct(seed, *, ladder=cont.LADDER, final_ephemeris=None, epochs=None, vinf_cap, residual_mode="vector", tol_kms=0.1) -> NarcContinuationResult`:
  - `final_ephemeris` default `Ephemeris("astropy")`; `epochs` default `None` → caller supplies (the batch computes them via `candidate_epochs`); if `None`, use `[0.0]`.
  - For each `t0` in `epochs`: walk the homotopy — build the ramp rung ephemerides from `ladder` (reuse `cont._ramp_schedule(nstep)` if available, else ramp `lam_e` 0→1 then `lam_i` 0→1 in `nstep` steps via `ramped_ephemeris`, final step = `final_ephemeris`). At each rung call `ballistic_correct(seed.sequence, seed.per_leg_revs, seed.per_leg_branch, t0_seed_sec=cur_t0, tof_seed_days=cur_tofs, period_sec=seed.period_sec, ephem=rung_ephem, vinf_cap=vinf_cap, residual_mode=residual_mode, tol_kms=tol_kms)`; re-seed `cur_t0, cur_tofs` from the converged result for the next rung (read `BallisticClosureResult` fields for the converged t0/tofs).
  - Keep the lowest-`max_residual_kms` converged final-ephemeris result across epochs. Set `emerged_vinf_kms` from the result's per-encounter V∞, `bend_feasible` from `.constraints_satisfied`, `winning_epoch_sec` from the best epoch.
  - Never raise on a corrector failure: a failing rung/epoch is skipped; if none converge, return a `converged=False` result with the best (lowest-residual) attempt's data.

- [ ] **Step 4: Run → PASS** (the smoke runs at the circular rung); ruff + mypy clean.
- [ ] **Step 5: Commit** `git add -A && git commit -m "search/#388: N-arc continuation driver (homotopy ramp + ballistic_correct per rung)"`

---

## Task 4: Batch + results note

**Files:** Create `scripts/narc_continuation_batch.py`; Create `docs/notes/2026-06-21-narc-continuation-results.md`

- [ ] **Step 1: Write `scripts/narc_continuation_batch.py`**:
  - For each `e` in `load_catalog().entries`: `phsi = descriptor_to_phsi(e.raw)`; skip if None. `cyc = assemble_cycler(m, phsi)`; skip if None. `seed = russell_parent_to_ballistic_seed(m, cyc, e.raw)`.
  - `target_phase = parent_phase_angle(...)` — actually the parent's REQUIRED phase comes from the idealized parent; for v1 use the assembled cycler's departure phase (document that this is the parent's beginning E-M relative phase; if not directly available, use `target_phase=0.0` and rely on the LaunchWindow scan — document the simplification).
  - `epochs = candidate_epochs(Ephemeris("astropy"), target_phase, launch_window_synodics=range(1,22), grid=100)`.
  - `res = narc_continuation_correct(seed, epochs=epochs, vinf_cap=VINF_CEILING_KMS)` (import the cap from `core.constants`).
  - Per-row line: `id [vlevel] converged={res.converged} res={res.max_residual_kms:.3e} vinf_emerged={...} vinf_anchor_E={seed.vinf_anchor_e_kms:.2f}/M={seed.vinf_anchor_m_kms:.2f} anchor_match={...} bend={res.bend_feasible}`. `anchor_match` = min |emerged − anchor| ≤ 0.5 for both E and M.
  - Flag `mcconaghy-2006-em-k2` PROPOSED V0→V1 (HELD) iff `res.converged and anchor_match and res.bend_feasible`.
  - Runlog `data/runs/narc-continuation-<time.strftime stamp>.jsonl`; summary counts; print "HELD — no writeback". `if __name__ == "__main__": main()`.
  - Run it: `timeout 3000 uv run python scripts/narc_continuation_batch.py` (real-eph homotopy over rows × epochs × rungs — may take many minutes; if it exceeds the timeout, reduce `launch_window_synodics` to `range(1,8)` and note it).

- [ ] **Step 2: Write the results note** `docs/notes/2026-06-21-narc-continuation-results.md`: per-row outcome (converged / emerged V∞ vs sourced anchor / bend-feasible / winning epoch), whether any row (esp. `mcconaghy-2006-em-k2`) closes on DE440 within 0.5 km/s, and the held verdict. State plainly if no row closes — the honest characterized terminal negative (true in-basin Russell parent reached DE440 but did not close). NO writeback.

- [ ] **Step 3: Commit** `git add scripts/narc_continuation_batch.py docs/notes/2026-06-21-narc-continuation-results.md data/runs/ && git commit -m "search/#388: N-arc continuation batch + results note (held, no writeback)"`

---

## Final verification
- [ ] `uv run pytest tests/search/test_narc_continuation.py -q` — PASS.
- [ ] `uv run ruff check . && uv run ruff format --check .` — clean.
- [ ] Confirm NO `data/catalogue.yaml` / `validate.py` edit (held-for-review gate).
- [ ] Report the DE440 closure verdict per row (esp. mcconaghy-2006-em-k2) and the held promotion outcome.
