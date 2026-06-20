# Multi-rev Descriptor-Seed Lane Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the DSM descriptor-seed closure lane represent multi-revolution resonant same-body legs by seeding each leg at its published arc ToF and enabling the already-built multi-rev Lambert branch selection (currently bypassed with `max_revs=0`).

**Architecture:** Three edits to `src/cyclerfinder/search/dsm_descriptor_seed.py` (sourced seed ToFs + rev cap + bracketed bounds) plus threading `max_revs` into the corrector call, then a batch re-run. No genome code changes — `dsm_leg.py` already implements multi-rev branch selection and is tested.

**Tech Stack:** Python 3, numpy, scipy, pytest, uv (run everything via `uv run`). Pre-commit runs ruff + mypy; keep `from __future__ import annotations` and full type hints.

**Spec:** `docs/superpowers/specs/2026-06-20-multirev-descriptor-seed-design.md`

**Sourced reference values** (computed from the catalogue + `core.constants.PLANETS`; `P_E = 365.26 d`, `P_M = 686.99 d`):

| row | sequence | seed ToFs (d) | max_revs |
|---|---|---|---|
| `russell-ch4-4.991gG2` | E-E-M-M | `[533.70, 150.0, 1026.21]` | 2 |
| `russell-ch4-9.353Gg2` | E-E-M-M | `[629.62, 85.0, 930.26]` | 2 |

Mapping rule: same-body legs (`a == b`) consume `free_return_arcs[*].tof_years` (non-null) in list order × 365.25; cross-body transit legs consume `invariants.transit_times_days` in order.

---

### Task 1: Sourced seed ToFs + per-leg rev cap

**Files:**
- Modify: `src/cyclerfinder/search/dsm_descriptor_seed.py` (the `DsmChainSeed` dataclass at lines 32-41; the seed-building block in `seed_dsm_chain_from_descriptor` at lines 110-171)
- Test: `tests/search/test_dsm_descriptor_seed.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/search/test_dsm_descriptor_seed.py`:

```python
def test_seed_tofs_are_sourced_not_slack() -> None:
    # The resonant same-body legs seed at the PUBLISHED arc ToF (free_return_arcs
    # tof_years x 365.25); the transit leg seeds at the sourced
    # invariants.transit_times_days value -- not the old slack heuristic.
    seed = dds.seed_dsm_chain_from_descriptor(_row("russell-ch4-4.991gG2"))
    assert seed is not None
    n_legs = len(seed.sequence) - 1
    tof = seed.x0[4 : 4 + n_legs]
    assert np.allclose(tof, [533.70, 150.0, 1026.21], atol=0.5)
    assert np.allclose(seed.per_leg_tof_days, [533.70, 150.0, 1026.21], atol=0.5)


def test_seed_max_revs_from_published_tof() -> None:
    # max_revs = max over legs of floor(arc_tof_days / body_period_days) + 1,
    # capped at Russell's 6-body-period generic-return ceiling. For an E-E-M-M
    # two-synodic row this is 2.
    for rid in ("russell-ch4-4.991gG2", "russell-ch4-9.353Gg2"):
        seed = dds.seed_dsm_chain_from_descriptor(_row(rid))
        assert seed is not None
        assert seed.max_revs == 2, rid
```

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py::test_seed_tofs_are_sourced_not_slack tests/search/test_dsm_descriptor_seed.py::test_seed_max_revs_from_published_tof -v`
Expected: FAIL — `AttributeError: 'DsmChainSeed' object has no attribute 'per_leg_tof_days'` (and `max_revs`).

- [ ] **Step 3: Add the two fields to `DsmChainSeed`**

In `src/cyclerfinder/search/dsm_descriptor_seed.py`, the dataclass currently ends:

```python
    transit_branch: str
    vinf_anchor_kms: float  # the row's sourced Russell-table V_inf cell
```

Append two fields:

```python
    transit_branch: str
    vinf_anchor_kms: float  # the row's sourced Russell-table V_inf cell
    per_leg_tof_days: tuple[float, ...] = ()  # SOURCED per-leg seed ToFs (audit)
    max_revs: int = 0  # global Lambert rev cap (Russell-sourced; see seed builder)
```

- [ ] **Step 4: Add a body-period helper + Russell cap constant**

First extend the constants import at the top of the file (the current import is
`from cyclerfinder.core.ephemeris import Ephemeris`) by adding above it:

```python
from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
```

Then, near the top of the module (after the `V1_TOLERANCE_KMS` / `YEAR_DAYS` /
`DAY_S` constants — `DAY_S = 86400.0` and `YEAR_DAYS = 365.25` already exist), add:

```python
# Russell 2004 SS2.1: the generic-return search caps ToF at 6 body periods; this is
# the global ceiling on how many revolutions the lane will enumerate.
RUSSELL_GENERIC_RETURN_BODY_PERIOD_CAP = 6


def _body_period_days(body: str) -> float:
    """Sidereal period of a planet (days), from its semi-major axis."""
    a_km = PLANETS[body].sma_au * AU_KM
    return float(2.0 * np.pi * np.sqrt(a_km**3 / MU_SUN_KM3_S2) / DAY_S)
```

- [ ] **Step 5: Replace the slack-heuristic seed-ToF block with sourced mapping**

In `seed_dsm_chain_from_descriptor`, the current block is:

```python
    g_tof_days = float(g_tof_yr * YEAR_DAYS)
    big_g_tof_days = float(big_g_tof_yr * YEAR_DAYS)
    transit_days = float(arc.tof_g_days)
    # Distribute remaining time across non-transit legs; at minimum 30 d.
    slack_days = max(30.0, big_g_tof_days - transit_days)
    tof_seed_days: list[float] = []
    for i in range(n_legs):
        body_a, body_b = sequence[i], sequence[i + 1]
        # Transit leg: E->M or M->E
        if {body_a, body_b} == {"E", "M"}:
            tof_seed_days.append(transit_days)
        else:
            tof_seed_days.append(max(slack_days, g_tof_days))
```

Replace it with the sourced mapping:

```python
    # SOURCED per-leg seed ToFs (spec 2026-06-20): a same-body resonant leg seeds at
    # its PUBLISHED arc ToF (free_return_arcs tof_years x 365.25, in list order); a
    # cross-body transit leg seeds at the row's sourced invariants.transit_times_days
    # (in order). Falls back to the computed arc transit (arc.tof_g_days) only if the
    # row has no transit_times_days entry left.
    arc_tofs_days = [
        float(a["tof_years"]) * YEAR_DAYS
        for a in (row.get("free_return_arcs") or [])
        if a.get("tof_years") is not None
    ]
    transit_days_list = [
        float(t) for t in ((row.get("invariants") or {}).get("transit_times_days") or [])
    ]
    big_g_tof_days = float(big_g_tof_yr * YEAR_DAYS)  # retained for the bounds cap below
    tof_seed_days = []
    per_leg_rev_cap = []
    arc_i = 0
    transit_i = 0
    for i in range(n_legs):
        body_a, body_b = sequence[i], sequence[i + 1]
        if body_a == body_b:
            # Same-body resonant return -> the next published arc ToF.
            if arc_i >= len(arc_tofs_days):
                return None  # more resonant legs than descriptor arcs -> cannot seed
            leg_tof = arc_tofs_days[arc_i]
            arc_i += 1
            rev_body = body_a
        else:
            # Cross-body transit -> the next sourced transit time (else computed arc).
            if transit_i < len(transit_days_list):
                leg_tof = transit_days_list[transit_i]
                transit_i += 1
            else:
                leg_tof = float(arc.tof_g_days)
            inner = body_a if PLANETS[body_a].sma_au <= PLANETS[body_b].sma_au else body_b
            rev_body = inner
        tof_seed_days.append(leg_tof)
        period_days = _body_period_days(rev_body)
        per_leg_rev_cap.append(int(np.floor(leg_tof / period_days)) + 1)
    max_revs = min(max(per_leg_rev_cap), RUSSELL_GENERIC_RETURN_BODY_PERIOD_CAP)
```

(`PLANETS` is already imported via the constants import added in Step 4.)

- [ ] **Step 6: Populate the two new seed fields**

At the `return DsmChainSeed(...)` near the end of the function, the current call ends:

```python
        transit_branch=str(arc.branch),
        vinf_anchor_kms=float(vinf_m),
    )
```

Change to:

```python
        transit_branch=str(arc.branch),
        vinf_anchor_kms=float(vinf_m),
        per_leg_tof_days=tuple(tof_seed_days),
        max_revs=int(max_revs),
    )
```

- [ ] **Step 7: Run the new + existing seed tests**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py -v`
Expected: PASS — the two new tests pass; `test_seed_built_for_reachable_descriptor_row` and `test_no_descriptor_row_returns_none` still pass (they assert shape/None, unaffected by ToF values).

- [ ] **Step 8: Commit**

```bash
git add src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git commit -m "search/#388: seed DSM legs at published arc ToFs + Russell rev cap"
```

---

### Task 2: Bracket resonant-leg ToF bounds at the published value

**Files:**
- Modify: `src/cyclerfinder/search/dsm_descriptor_seed.py` (the bounds block in `seed_dsm_chain_from_descriptor`, currently lines ~129-148)
- Test: `tests/search/test_dsm_descriptor_seed.py`

- [ ] **Step 1: Write the failing test**

```python
def test_resonant_leg_tof_bounds_bracket_published() -> None:
    # A same-body resonant leg's ToF box brackets its published seed ToF (0.7x..1.3x)
    # so the corrector cannot collapse it to the degenerate near-zero-ToF single-rev
    # solution; the transit leg keeps the sequence-keyed bound.
    seed = dds.seed_dsm_chain_from_descriptor(_row("russell-ch4-4.991gG2"))
    assert seed is not None
    n_legs = len(seed.sequence) - 1
    lower = seed.bounds.lower[4 : 4 + n_legs]
    upper = seed.bounds.upper[4 : 4 + n_legs]
    # leg 0 (E->E) is resonant: published 533.70 d -> [373.6, 693.8]
    assert lower[0] == __import__("pytest").approx(0.7 * 533.70, abs=1.0)
    assert upper[0] == __import__("pytest").approx(1.3 * 533.70, abs=1.0)
    # leg 2 (M->M) resonant: published 1026.21 d
    assert lower[2] == __import__("pytest").approx(0.7 * 1026.21, abs=1.0)
    assert upper[2] == __import__("pytest").approx(1.3 * 1026.21, abs=1.0)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py::test_resonant_leg_tof_bounds_bracket_published -v`
Expected: FAIL — the resonant-leg lower bound is the generic `30.0` (sequence-keyed), not `0.7 * published`.

- [ ] **Step 3: Bracket the resonant-leg bounds**

The current bounds block is:

```python
    tof_upper_cap = max(big_g_tof_days * 2.0, 2000.0)
    upper_capped = bounds_raw.upper.copy()
    for i in range(n_legs):
        ui = 4 + i  # tof slot i in the flat bounds vector
        if not np.isfinite(upper_capped[ui]):
            upper_capped[ui] = tof_upper_cap
    bounds = dsm_leg.DsmBounds(lower=bounds_raw.lower, upper=upper_capped)
```

Replace with (also bracketing the resonant-leg lower/upper around the sourced ToF):

```python
    tof_upper_cap = max(big_g_tof_days * 2.0, 2000.0)
    upper_capped = bounds_raw.upper.copy()
    lower_capped = bounds_raw.lower.copy()
    for i in range(n_legs):
        ui = 4 + i  # tof slot i in the flat bounds vector
        body_a, body_b = sequence[i], sequence[i + 1]
        if body_a == body_b:
            # Resonant same-body leg: bracket the PUBLISHED seed ToF (spec Change 3)
            # so the corrector stays in the multi-rev resonant basin and cannot reach
            # the degenerate near-zero-ToF single-rev region.
            published = tof_seed_days[i]
            lower_capped[ui] = 0.7 * published
            upper_capped[ui] = 1.3 * published
        elif not np.isfinite(upper_capped[ui]):
            upper_capped[ui] = tof_upper_cap
    bounds = dsm_leg.DsmBounds(lower=lower_capped, upper=upper_capped)
```

- [ ] **Step 4: Run the test**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py -v`
Expected: PASS — the new bracket test passes; all earlier seed tests still pass.

- [ ] **Step 5: Commit**

```bash
git add src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git commit -m "search/#388: bracket resonant-leg ToF bounds at the published value"
```

---

### Task 3: Thread max_revs into the closer + surface emerged revs

**Files:**
- Modify: `src/cyclerfinder/search/dsm_descriptor_seed.py` (`DsmClosureResult` dataclass lines 174-210; `close_row_dsm` body lines 247-285)
- Test: `tests/search/test_dsm_descriptor_seed.py`

- [ ] **Step 1: Write the failing test**

```python
def test_closure_result_carries_seed_max_revs_and_n_revs() -> None:
    # close_row_dsm threads the seed's Russell rev cap into the corrector and reports
    # the emerged per-leg revolution count for the runlog/audit.
    row = _row("russell-ch4-4.991gG2")
    ephem = Ephemeris("astropy")
    res = dds.close_row_dsm(row, ephem)
    assert res.max_revs_used == 2
    assert isinstance(res.n_revs_per_leg, tuple)
    assert len(res.n_revs_per_leg) == 3  # E-E-M-M -> 3 legs
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py::test_closure_result_carries_seed_max_revs_and_n_revs -v`
Expected: FAIL — `AttributeError: 'DsmClosureResult' object has no attribute 'max_revs_used'`.

- [ ] **Step 3: Add the two fields to `DsmClosureResult`**

The dataclass currently ends:

```python
    anchor_match: bool
    hyperbolic_impossible: bool
    seed: DsmChainSeed | None
```

Change to:

```python
    anchor_match: bool
    hyperbolic_impossible: bool
    seed: DsmChainSeed | None
    max_revs_used: int = 0  # the Russell rev cap passed to the corrector
    n_revs_per_leg: tuple[int, ...] = ()  # EMERGED per-leg back-arc rev count (audit)
```

- [ ] **Step 4: Thread max_revs through `close_row_dsm` + populate the fields**

In `close_row_dsm`, the `None`-seed early return currently is:

```python
    seed = seed_dsm_chain_from_descriptor(row)
    if seed is None:
        return DsmClosureResult(
            converged=False,
            max_residual_kms=float("nan"),
            dv_dsm_kms=(),
            vinf_per_encounter_kms=(),
            vinf_anchor_kms=float("nan"),
            anchor_match=False,
            hyperbolic_impossible=False,
            seed=None,
        )
```

Leave it unchanged (the new fields default to `0` / `()`).

The corrector call currently is:

```python
    res = dsm_leg.dsm_chain_correct(
        seed.x0,
        sequence=seed.sequence,
        ephem=ephem,
        bounds=seed.bounds,
        tol_kms=tol_kms,
        charge_flyby_continuity=True,
        gradient=gradient,
    )
```

Add the rev cap:

```python
    res = dsm_leg.dsm_chain_correct(
        seed.x0,
        sequence=seed.sequence,
        ephem=ephem,
        bounds=seed.bounds,
        tol_kms=tol_kms,
        max_revs=seed.max_revs,
        charge_flyby_continuity=True,
        gradient=gradient,
    )
```

The final `return DsmClosureResult(...)` currently ends:

```python
        anchor_match=best <= V1_TOLERANCE_KMS,
        hyperbolic_impossible=hyper,
        seed=seed,
    )
```

Change to:

```python
        anchor_match=best <= V1_TOLERANCE_KMS,
        hyperbolic_impossible=hyper,
        seed=seed,
        max_revs_used=int(seed.max_revs),
        n_revs_per_leg=tuple(int(n) for n in res.n_revs_per_leg),
    )
```

- [ ] **Step 5: Run the test**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py::test_closure_result_carries_seed_max_revs_and_n_revs -v`
Expected: PASS. (This test calls the real ephemeris; it does not require convergence, only the new fields.)

- [ ] **Step 6: Run the full seed/closer test module + ruff/mypy**

Run: `uv run pytest tests/search/test_dsm_descriptor_seed.py tests/search/test_dsm_leg.py -q && uv run ruff check src/cyclerfinder/search/dsm_descriptor_seed.py && uv run ruff format --check src/cyclerfinder/search/dsm_descriptor_seed.py`
Expected: all PASS / no lint errors.

- [ ] **Step 7: Commit**

```bash
git add src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git commit -m "search/#388: thread Russell rev cap into closer + surface emerged revs"
```

---

### Task 4: Batch re-run + results note

**Files:**
- Modify: `scripts/dsm_closure_batch.py` (the per-row print + the record dict, lines ~80-101)
- Create: `docs/notes/2026-06-20-dsm-closure-batch-multirev-results.md`

- [ ] **Step 1: Add emerged revs to the batch record + print line**

In `scripts/dsm_closure_batch.py`, the record dict currently includes:

```python
        rec = {
            "id": rid,
            "validation_level": vlevel,
            "converged": res.converged,
            "max_residual_kms": res.max_residual_kms,
            "vinf_anchor_kms": res.vinf_anchor_kms,
            "vinf_per_encounter_kms": list(res.vinf_per_encounter_kms),
            "anchor_match": res.anchor_match,
            "dsm_dv_kms": list(res.dv_dsm_kms),
            "hyperbolic_impossible": res.hyperbolic_impossible,
            "proposed_promotion": "V0->V1" if promote else None,
            "wall_s": round(dt, 1),
        }
```

Add the two new audit fields:

```python
            "hyperbolic_impossible": res.hyperbolic_impossible,
            "max_revs_used": res.max_revs_used,
            "n_revs_per_leg": list(res.n_revs_per_leg),
            "proposed_promotion": "V0->V1" if promote else None,
            "wall_s": round(dt, 1),
```

The per-row print currently is:

```python
        print(
            f"[{time.time() - t0:6.0f}s] {rid:24s} [{vlevel}] "
            f"conv={res.converged!s:5} res={res.max_residual_kms:.3e} "
            f"anchor={res.vinf_anchor_kms:.3f} match={res.anchor_match!s:5} "
            f"dsmdV={sum(res.dv_dsm_kms):.3f}{tag}",
            flush=True,
        )
```

Change the last format line to include revs:

```python
        print(
            f"[{time.time() - t0:6.0f}s] {rid:24s} [{vlevel}] "
            f"conv={res.converged!s:5} res={res.max_residual_kms:.3e} "
            f"anchor={res.vinf_anchor_kms:.3f} match={res.anchor_match!s:5} "
            f"dsmdV={sum(res.dv_dsm_kms):.3f} revs={list(res.n_revs_per_leg)}{tag}",
            flush=True,
        )
```

- [ ] **Step 2: Run the batch on real DE440**

Run: `timeout 1800 uv run python scripts/dsm_closure_batch.py`
Expected: completes; prints per-row `conv=`, `res=`, `revs=[...]` and a summary. The success signal is a validated **regression** row (NOT S1L1 — i.e. NOT `4.991gG2` / `mcconaghy-2006-em-k2`) converging with `res < 0.1` and `match=True` and `revs` showing >0 on its resonant legs. (Note `9.353Gg2` is V1; the cleanest non-S1L1 V3 regression row is `russell-ch4-8.049gGf2`.)

- [ ] **Step 3: Write the results note**

Create `docs/notes/2026-06-20-dsm-closure-batch-multirev-results.md` recording, factually from the run output: the per-row table (converged / residual / anchor-match / emerged revs / DSM dV), comparison against the 2026-06-20 single-rev baseline (0/9, residuals 20-77 km/s), whether the regression rows now reconverge (the golden cross-check), and the S1L1 promotion-target outcome — **held for review, no writeback** regardless. If no regression row reconverges, state that plainly: it means the obstruction is deeper than rev count (per the spec's error-handling gate — stop, do not loosen tolerance).

- [ ] **Step 4: Commit**

```bash
git add scripts/dsm_closure_batch.py docs/notes/2026-06-20-dsm-closure-batch-multirev-results.md data/runs/
git commit -m "search/#388: multi-rev closure batch re-run + results note"
```

---

## Final verification

- [ ] Run: `uv run pytest tests/search/ -q` — Expected: PASS (no regressions in the search suite).
- [ ] Run: `uv run ruff check . && uv run ruff format --check .` — Expected: clean.
- [ ] Confirm: no `data/catalogue.yaml` edit and no `validation_level` change anywhere (the honesty gate — any promotion is held for session review, not written here).
- [ ] Report the regression-row reconvergence verdict and the S1L1 outcome to the user.
