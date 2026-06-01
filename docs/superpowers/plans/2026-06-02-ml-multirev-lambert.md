# M-L · Multi-Rev Lambert Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `lambert(..., max_revs=N)` return the full multi-revolution solution set (the `n=0` single-rev solution plus, for each `n∈[1,N]`, the two `low`/`high` branches when they exist), then thread it through real-ephemeris construction so catalogue legs with `n_revs>0` can be built and verified instead of skipped.

**Architecture:** The universal-variable Lambert iteration already implemented for single-rev (`core/lambert.py`) extends to multi-rev by searching the universal variable `z` in per-revolution domains `z ∈ ((2πn)², (2π(n+1))²)`. Inside each `n≥1` domain the time-of-flight curve `t(z)` is convex with a single minimum `t_min(n)`; if the requested `tof > t_min(n)` two solutions exist (one on each side of the minimum → `low`/`high` branches), otherwise revolution `n` is infeasible and contributes nothing. We keep the existing single-rev fast path bit-for-bit unchanged for `max_revs=0`, add a multi-rev solver alongside it, and validate every multi-rev solution against `lamberthub`'s independent `izzo2015`/`gooding1990` multi-rev solvers (the same sourced-crosscheck discipline the single-rev gate already uses). Downstream, the `n_revs>0` rejection gates are removed and `n_revs`/`branch` are threaded from catalogue leg metadata into the solver.

**Tech stack:** Python 3.11, numpy (solver inner loop, numpy-only — no SciPy on the production path), `lamberthub` (dev-dependency crosscheck), pytest + xdist, uv-managed venv, ruff + mypy.

---

## Branch-selection contract (settled)

The multi-rev minimum-time point `z_min(n)` splits revolution `n` into two solutions:

- **`branch="low"`** — the solution with `z < z_min(n)` (smaller universal variable, lower transfer energy / the shorter-semimajor-axis member). Maps to `lamberthub`'s `low_path=True`.
- **`branch="high"`** — the solution with `z > z_min(n)` (larger universal variable). Maps to `lamberthub`'s `low_path=False`.

The `low`/`high` ↔ `low_path` mapping is **asserted empirically** in Task 4 against `lamberthub` rather than assumed; if the convention is inverted the fix is a one-line swap and the test is the arbiter.

**Catalogue default:** catalogue legs carry `n_revs` but **no `branch` field** (`catalog.py:583`). When a leg has `n_revs>0` and no explicit branch, construction requests **`branch="low"`** (matches `lamberthub`'s default and the lower-energy physical convention). An optional per-leg `branch` string in the catalogue YAML overrides this. This default is documented in `construct_real_ephemeris_cycler`'s docstring (Task 6).

---

## Solution-ordering & return contract (settled)

`lambert(r1, r2, tof, *, max_revs=N)` returns `list[LambertSolution]` ordered:

1. `LambertSolution(n_revs=0, branch="single", …)` — always first, always present on success (unchanged from today).
2. then, for `n = 1, 2, …, N` in ascending order, the feasible branches for that `n` in the order `low` then `high`.

A revolution `n` that is infeasible for the given `tof` (i.e. `tof <= t_min(n)`) contributes **zero** entries — it is silently skipped, never an error. So `lambert(..., max_revs=2)` returns between 1 and 5 solutions. `max_revs=0` returns exactly the length-1 single-rev list, byte-for-byte identical to today.

---

## File structure

| File | Responsibility | Change |
|---|---|---|
| `src/cyclerfinder/core/lambert.py` | UV Lambert solver | Add multi-rev solver path + min-time bracketing; refactor single-rev core into a shared `_solve_uv_branch` helper reused by both; extend `lambert()` to loop revolutions. Single-rev result unchanged. |
| `tests/core/test_lambert_multirev.py` | Multi-rev unit + crosscheck gate | **New file.** All multi-rev unit tests + the lamberthub multi-rev crosscheck gate. |
| `tests/core/test_lambert.py` | Existing single-rev tests | Update `test_lambert_max_revs_stub` (the stub assertion is now wrong) → rename to `test_lambert_max_revs_returns_multirev`. |
| `src/cyclerfinder/verify/real_closure.py` | Real-ephemeris construction | Delete the `MultiRevLambertRequiredError` raise (`:519-524`); thread `n_revs` (and optional `branch`) from leg metadata into `lambert()`; set `Leg.n_revs`/`branch` from the chosen solution. |
| `src/cyclerfinder/verify/crosscheck.py` | V1 lamberthub crosscheck | Allow `n_revs>0` legs (`:128-132`): call `izzo2015`/`gooding1990` with `M=leg.n_revs`, `low_path=(leg.branch!="high")`; stop skipping multi-rev legs in `crosscheck_cycler` (`:182-183`). |
| `tests/verify/test_real_closure.py` | Construction gate | Repurpose `test_construct_raises_on_multi_rev_leg` (`:506`) → `test_construct_builds_multi_rev_leg` (now builds, no longer raises). The `n_revs=1` → `v3-skipped-multirev` test (`:686`) becomes a *closes-or-construction-error* test. |
| `src/cyclerfinder/verify/real_closure.py` | Skip bookkeeping | Remove `MultiRevLambertRequiredError`-driven `v3-skipped-multirev` status path (`:746-751`); the docstring/`__all__` references. Decide retention of the exception class (kept as deprecated-unused or deleted — Task 8). |

**Out of scope for M-L** (explicitly): inclination/3D-frame work (M-3D), N-encounter loader reclassification (M-N), the `optimise_cell_ephemeris` stub (M-ED). M-L only makes the *primitive* and its single-leg construction path multi-rev-capable.

---

## Tasks

### Task 1: Multi-rev time-of-flight minimum finder

**Files:**
- Modify: `src/cyclerfinder/core/lambert.py`
- Test: `tests/core/test_lambert_multirev.py` (create)

The per-revolution time curve `t(z)` on `z ∈ ((2πn)², (2π(n+1))²)` is convex with a single interior minimum. We need `t_min(n)` and its location `z_min(n)` to (a) decide feasibility and (b) split the two branches. Implement a numpy-only golden-section / bounded minimiser on `t(z)` within the open interval.

- [ ] **Step 1: Write the failing test**

```python
# tests/core/test_lambert_multirev.py
"""Multi-revolution Lambert: unit tests + lamberthub crosscheck gate.

The correctness gate compares the in-house multi-rev solver against
lamberthub.izzo2015 / gooding1990 with matching M (revs) and low_path
(branch) — the same sourced-crosscheck discipline as the single-rev gate
in test_lambert.py. EXPECTED values trace to an independent third-party
solver, never to a value our own solver computed.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.lambert import _min_time_of_revolution

from .conftest import Leg


def test_min_time_of_revolution_is_a_lower_bound(leg_long: Leg) -> None:
    """t_min(n=1) is finite, positive, and below the n=1 domain's endpoints.

    The endpoints z=(2*pi*n)^2 and z=(2*pi*(n+1))^2 are time singularities
    (t -> +inf); the interior minimum must be strictly smaller than t at a
    sample just inside each endpoint.
    """
    r1_n = float(np.linalg.norm(leg_long.r1))
    r2_n = float(np.linalg.norm(leg_long.r2))
    cos_dnu = float(np.dot(leg_long.r1, leg_long.r2) / (r1_n * r2_n))
    from math import acos, sqrt

    dnu = acos(max(min(cos_dnu, 1.0), -1.0))
    sin_dnu = float(np.sin(dnu))
    a_coef = sin_dnu * sqrt(r1_n * r2_n / (1.0 - cos_dnu))

    z_min, t_min = _min_time_of_revolution(1, a_coef, r1_n, r2_n, MU_SUN_KM3_S2)
    assert (2.0 * np.pi) ** 2 < z_min < (4.0 * np.pi) ** 2
    assert t_min > 0.0
    assert np.isfinite(t_min)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/core/test_lambert_multirev.py::test_min_time_of_revolution_is_a_lower_bound -v`
Expected: FAIL with `ImportError: cannot import name '_min_time_of_revolution'`.

- [ ] **Step 3: Write minimal implementation**

Add to `src/cyclerfinder/core/lambert.py` (after `_dt_dz`):

```python
def _min_time_of_revolution(
    n: int, a_coef: float, r1_n: float, r2_n: float, mu: float
) -> tuple[float, float]:
    """Return ``(z_min, t_min)`` for revolution ``n >= 1``.

    On the open interval ``z in ((2*pi*n)**2, (2*pi*(n+1))**2)`` the
    universal-variable time-of-flight ``t(z)`` is convex with a single
    interior minimum. We locate it with a bounded golden-section search
    (numpy-only; no SciPy on the production path). ``t_min`` is the shortest
    time of flight achievable with exactly ``n`` full revolutions; a
    requested ``tof <= t_min`` means revolution ``n`` is infeasible.
    """
    lo = (2.0 * np.pi * n) ** 2
    hi = (2.0 * np.pi * (n + 1)) ** 2
    # Stay strictly inside the open interval: t -> +inf at both endpoints.
    span = hi - lo
    a = lo + 1.0e-6 * span
    b = hi - 1.0e-6 * span

    def _t(z_in: float) -> float:
        t_z, _ = _t_of_z(z_in, a_coef, r1_n, r2_n, mu)
        return t_z

    inv_phi = (sqrt(5.0) - 1.0) / 2.0  # 1/golden ratio
    c = b - inv_phi * (b - a)
    d = a + inv_phi * (b - a)
    fc = _t(c)
    fd = _t(d)
    for _ in range(200):
        if abs(b - a) < 1.0e-9 * (abs(a) + abs(b)):
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = _t(c)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = _t(d)
    z_min = 0.5 * (a + b)
    t_min, _ = _t_of_z(z_min, a_coef, r1_n, r2_n, mu)
    return z_min, t_min
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/core/test_lambert_multirev.py::test_min_time_of_revolution_is_a_lower_bound -v`
Expected: PASS.

- [ ] **Step 5: Commit** (only on explicit user request — otherwise hold)

```bash
uv run ruff check src/cyclerfinder/core/lambert.py tests/core/test_lambert_multirev.py
uv run ruff format --check src/cyclerfinder/core/lambert.py tests/core/test_lambert_multirev.py
git add src/cyclerfinder/core/lambert.py tests/core/test_lambert_multirev.py
git commit -m "core/lambert: add multi-rev time-of-flight minimum finder"
```

---

### Task 2: Extract a shared single-branch UV Newton solver

**Files:**
- Modify: `src/cyclerfinder/core/lambert.py`
- Test: `tests/core/test_lambert.py` (existing single-rev gates are the regression guard)

The existing `lambert()` body contains the bracket-bootstrap + Newton iteration that, given a bracket `[z_lo, z_hi]` known to contain the root and a sign-monotone `f(z)=t(z)-tof`, returns the converged `z`. Multi-rev needs to call this same Newton machinery twice per revolution (once on each side of `z_min`). Refactor the bracketed-Newton core into `_solve_uv_branch(z_lo, z_hi, a_coef, r1_n, r2_n, tof, mu)` returning the converged `z`, and have the single-rev path call it. **No behavioural change** — the single-rev gates must stay green bit-for-bit.

- [ ] **Step 1: Run the existing single-rev gate to capture the green baseline**

Run: `uv run pytest tests/core/test_lambert.py tests/core/test_lambert_kepler_consistency.py -v`
Expected: PASS (record this is the baseline the refactor must preserve).

- [ ] **Step 2: Extract the helper**

Add to `src/cyclerfinder/core/lambert.py` a function that takes a known sign-changing bracket and runs the existing safeguarded Newton loop (lines 333-370 of the current file), returning `z`:

```python
def _solve_uv_branch(
    z_lo: float,
    z_hi: float,
    a_coef: float,
    r1_n: float,
    r2_n: float,
    tof: float,
    mu: float,
) -> float:
    """Safeguarded Newton on f(z)=t(z)-tof within a known bracket [z_lo, z_hi].

    Assumes t(z) is monotone on the bracket (true on each side of a
    revolution's minimum, and across the whole single-rev range). Falls back
    to bisection when a Newton step leaves the bracket or drives y(z) < 0.
    Raises LambertConvergenceError on failure. This is the exact iteration
    the single-rev path used before multi-rev; extracted so the multi-rev
    path can reuse it per branch.
    """

    def _y_only(z_in: float) -> float:
        c = stumpff_c(z_in)
        s = stumpff_s(z_in)
        return r1_n + r2_n + a_coef * (z_in * s - 1.0) / sqrt(c)

    z = 0.5 * (z_lo + z_hi)
    residual: float = 0.0
    for _it in range(_NEWTON_MAX_ITER):
        try:
            t_z, y = _t_of_z(z, a_coef, r1_n, r2_n, mu)
        except ValueError:
            z = 0.5 * (z_lo + z_hi)
            continue
        residual = t_z - tof
        if abs(residual) / tof < _NEWTON_TOL_REL:
            return z
        if residual < 0.0:
            z_lo = z
        else:
            z_hi = z
        dt = _dt_dz(z, y, a_coef, mu)
        z_next = z - residual / dt if dt != 0.0 else z
        if not (z_lo < z_next < z_hi) or _y_only(z_next) < 0.0:
            z_next = 0.5 * (z_lo + z_hi)
        if abs(z_next - z) < _NEWTON_TOL_DZ:
            return z_next
        z = z_next
    raise LambertConvergenceError(z, residual)
```

Then replace the inline Newton loop in `lambert()` (current lines 329-370) with a call to `_solve_uv_branch(z_lo, z_hi, a_coef, r1_n, r2_n, tof, mu)` after the existing single-rev bracket bootstrap establishes `z_lo`/`z_hi`. Keep the existing bracket bootstrap (lines 271-327) and the final Lagrange-coefficient block (372-390) unchanged.

- [ ] **Step 3: Run the single-rev gate to verify no regression**

Run: `uv run pytest tests/core/test_lambert.py tests/core/test_lambert_kepler_consistency.py -v`
Expected: PASS (identical to Step 1 baseline).

- [ ] **Step 4: Commit** (hold unless user asked)

```bash
uv run ruff check src/cyclerfinder/core/lambert.py
uv run ruff format --check src/cyclerfinder/core/lambert.py
git add src/cyclerfinder/core/lambert.py
git commit -m "core/lambert: extract shared bracketed-Newton branch solver"
```

---

### Task 3: Multi-rev solver loop in `lambert()`

**Files:**
- Modify: `src/cyclerfinder/core/lambert.py`
- Test: `tests/core/test_lambert_multirev.py`

Extend `lambert()` to, after computing the single-rev solution, loop `n = 1..max_revs`: for each `n` compute `z_min(n), t_min(n)` (Task 1); if `tof <= t_min(n)` skip; else solve two branches with `_solve_uv_branch` on `[(2πn)²+eps, z_min]` (high-time side → `low` branch candidate) and `[z_min, (2π(n+1))²-eps]`, build a `LambertSolution` from each converged `z` via the shared Lagrange-coefficient computation, and append. Factor the Lagrange-coefficient → `(v1,v2)` computation into `_velocities_from_z(z, ...)` so single-rev and multi-rev share it.

- [ ] **Step 1: Write the failing test (count + structure, no self-computed values)**

```python
def test_lambert_max_revs_returns_n0_plus_branches(leg_long: Leg) -> None:
    """A 500-day Earth->Mars leg admits n=0 and (at least) n=1 low/high.

    Asserts STRUCTURE only — counts, n_revs labels, branch labels, shapes.
    Numerical correctness is the crosscheck gate (test below), whose
    EXPECTED side is lamberthub, not our own solver.
    """
    from cyclerfinder.core.lambert import lambert

    sols = lambert(leg_long.r1, leg_long.r2, leg_long.tof, max_revs=1)
    assert sols[0].n_revs == 0 and sols[0].branch == "single"
    n1 = [s for s in sols if s.n_revs == 1]
    branches = sorted(s.branch for s in n1)
    assert branches == ["high", "low"]  # both branches feasible for 500 d
    for s in sols:
        assert s.v1.shape == (3,) and s.v2.shape == (3,)
        assert s.v1.dtype == np.float64


def test_lambert_infeasible_revolution_is_skipped(leg_short: Leg) -> None:
    """A 50-day Earth->Earth arc is far too short for even n=1: only single-rev.

    t_min(1) for a 50-day short arc vastly exceeds 50 days, so revolution 1
    is infeasible and contributes zero solutions (skipped, not an error).
    """
    from cyclerfinder.core.lambert import lambert

    sols = lambert(leg_short.r1, leg_short.r2, leg_short.tof, max_revs=3)
    assert all(s.n_revs == 0 for s in sols)
    assert len(sols) == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/core/test_lambert_multirev.py -k "returns_n0_plus_branches or infeasible" -v`
Expected: FAIL (`lambert` still returns only the single-rev solution).

- [ ] **Step 3: Implement the multi-rev loop**

In `src/cyclerfinder/core/lambert.py`, factor velocity construction:

```python
def _velocities_from_z(
    z: float, a_coef: float, r1_arr: Vec3, r2_arr: Vec3, r1_n: float, r2_n: float, mu: float
) -> tuple[Vec3, Vec3]:
    """Lagrange-coefficient velocities at the converged universal variable z."""
    c = stumpff_c(z)
    s = stumpff_s(z)
    y = r1_n + r2_n + a_coef * (z * s - 1.0) / sqrt(c)
    f = 1.0 - y / r1_n
    g = a_coef * sqrt(y / mu)
    g_dot = 1.0 - y / r2_n
    v1 = (r2_arr - f * r1_arr) / g
    v2 = (g_dot * r2_arr - r1_arr) / g
    return v1.astype(np.float64, copy=False), v2.astype(np.float64, copy=False)
```

Replace the single-rev return block (current 372-390) to use `_velocities_from_z`, then append the multi-rev loop before `return`:

```python
    v1, v2 = _velocities_from_z(z, a_coef, r1_arr, r2_arr, r1_n, r2_n, mu)
    solutions = [LambertSolution(n_revs=0, branch="single", v1=v1, v2=v2)]

    eps = 1.0e-6
    for n in range(1, max_revs + 1):
        z_min, t_min = _min_time_of_revolution(n, a_coef, r1_n, r2_n, mu)
        if tof <= t_min:
            continue  # revolution n infeasible for this tof; skip silently
        lo_endpoint = (2.0 * np.pi * n) ** 2
        hi_endpoint = (2.0 * np.pi * (n + 1)) ** 2
        span = hi_endpoint - lo_endpoint
        # "low" branch: z in (endpoint_lo, z_min); "high": z in (z_min, endpoint_hi).
        # The empirical low/high <-> low_path mapping is asserted in the
        # crosscheck test; if inverted, swap these two labels.
        for branch_label, z_a, z_b in (
            ("low", lo_endpoint + eps * span, z_min),
            ("high", z_min, hi_endpoint - eps * span),
        ):
            try:
                z_sol = _solve_uv_branch(z_a, z_b, a_coef, r1_n, r2_n, tof, mu)
            except LambertConvergenceError:
                continue
            vb1, vb2 = _velocities_from_z(z_sol, a_coef, r1_arr, r2_arr, r1_n, r2_n, mu)
            solutions.append(
                LambertSolution(n_revs=n, branch=branch_label, v1=vb1, v2=vb2)
            )
    return solutions
```

Note: `_solve_uv_branch` assumes monotone `t(z)` on the bracket — true on each side of `z_min`. Each half-interval brackets exactly one root because `t(z_a)` and `t(z_b)` straddle `tof` (one endpoint is the `+inf` singularity side, the other is `t_min < tof`).

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/core/test_lambert_multirev.py -k "returns_n0_plus_branches or infeasible" -v`
Expected: PASS.

- [ ] **Step 5: Update the stale single-rev stub test**

In `tests/core/test_lambert.py`, replace `test_lambert_max_revs_stub` (lines 45-49):

```python
def test_lambert_max_revs_returns_multirev(leg_long: Leg) -> None:
    """``max_revs >= 1`` now returns multi-rev solutions when feasible."""
    sols = lambert(leg_long.r1, leg_long.r2, leg_long.tof, max_revs=1)
    assert sols[0].n_revs == 0
    assert any(s.n_revs == 1 for s in sols)
```

- [ ] **Step 6: Run the full core lambert suite**

Run: `uv run pytest tests/core/test_lambert.py tests/core/test_lambert_multirev.py -v`
Expected: PASS.

- [ ] **Step 7: Commit** (hold unless asked)

```bash
uv run ruff check src/cyclerfinder/core/lambert.py tests/core/
uv run ruff format --check src/cyclerfinder/core/lambert.py tests/core/
git add src/cyclerfinder/core/lambert.py tests/core/test_lambert.py tests/core/test_lambert_multirev.py
git commit -m "core/lambert: implement multi-rev solver loop (n>=1 low/high branches)"
```

---

### Task 4: Multi-rev crosscheck gate against lamberthub (the correctness gate)

**Files:**
- Test: `tests/core/test_lambert_multirev.py`

This is the binding correctness gate and the source-attested check: in-house multi-rev velocities must agree with `lamberthub`'s independent `izzo2015` and `gooding1990` multi-rev solvers to < 1e-3 m/s, for matching `M=n_revs` and `low_path` per branch. This also empirically nails the `low`/`high` ↔ `low_path` mapping from Task 3.

- [ ] **Step 1: Write the crosscheck gate**

```python
@pytest.mark.parametrize("branch,low_path", [("low", True), ("high", False)])
def test_multirev_crosscheck_against_lamberthub(
    leg_long: Leg, branch: str, low_path: bool
) -> None:
    """In-house n=1 branch agrees with lamberthub izzo+gooding < 1e-3 m/s.

    EXPECTED values come from lamberthub (independent third-party solvers),
    not from our own solver — satisfies the sourced-golden discipline.
    If this fails only by swapping branch<->low_path, the low/high labels in
    lambert() are inverted: swap them there, not here.
    """
    from lamberthub import gooding1990, izzo2015  # type: ignore[import-untyped]

    from cyclerfinder.core.lambert import lambert

    sols = lambert(leg_long.r1, leg_long.r2, leg_long.tof, max_revs=1)
    mine = next(s for s in sols if s.n_revs == 1 and s.branch == branch)

    r1 = np.asarray(leg_long.r1, dtype=np.float64)
    r2 = np.asarray(leg_long.r2, dtype=np.float64)
    v1_izzo, v2_izzo = izzo2015(
        MU_SUN_KM3_S2, r1, r2, leg_long.tof, M=1, prograde=True, low_path=low_path
    )
    v1_g, v2_g = gooding1990(
        MU_SUN_KM3_S2, r1, r2, leg_long.tof, M=1, prograde=True, low_path=low_path
    )
    worst_mps = 1000.0 * max(
        float(np.linalg.norm(mine.v1 - v1_izzo)),
        float(np.linalg.norm(mine.v2 - v2_izzo)),
        float(np.linalg.norm(mine.v1 - v1_g)),
        float(np.linalg.norm(mine.v2 - v2_g)),
    )
    assert worst_mps < 1.0e-3, (branch, worst_mps)
```

- [ ] **Step 2: Run the gate**

Run: `uv run pytest tests/core/test_lambert_multirev.py::test_multirev_crosscheck_against_lamberthub -v`
Expected: PASS for both branches. **If both fail**, the convergence tolerance or branch geometry is wrong — debug the solver, do not loosen the 1e-3 gate. **If they pass only when swapped** (`low`↔`high`), invert the two branch labels in `lambert()` (Task 3, Step 3) and re-run; the test asserts the mapping.

- [ ] **Step 3: Commit** (hold unless asked)

```bash
uv run ruff check tests/core/test_lambert_multirev.py
uv run ruff format --check tests/core/test_lambert_multirev.py
git add tests/core/test_lambert_multirev.py
git commit -m "core/lambert: multi-rev crosscheck gate vs lamberthub (sourced)"
```

---

### Task 5: Add `max_revs` to `lambert_crosscheck` helper

**Files:**
- Modify: `src/cyclerfinder/core/lambert.py` (`lambert_crosscheck`, lines 393-440)
- Test: `tests/core/test_lambert_multirev.py`

`lambert_crosscheck` currently hardwires single-rev (`M=0`, picks `[0]`). Generalise it to accept `n_revs`/`branch` so M-N and M-ED can reuse it for multi-rev legs.

- [ ] **Step 1: Write the failing test**

```python
def test_lambert_crosscheck_multirev(leg_long: Leg) -> None:
    """lambert_crosscheck(..., n_revs=1, branch='low') agrees < 1e-3 m/s."""
    from cyclerfinder.core.lambert import lambert_crosscheck

    res = lambert_crosscheck(
        leg_long.r1, leg_long.r2, leg_long.tof, n_revs=1, branch="low"
    )
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/core/test_lambert_multirev.py::test_lambert_crosscheck_multirev -v`
Expected: FAIL (`lambert_crosscheck` has no `n_revs`/`branch` kwargs → `TypeError`).

- [ ] **Step 3: Implement**

Update `lambert_crosscheck` signature to `(..., prograde=True, n_revs=0, branch="single")`. Select the in-house solution matching `(n_revs, branch)` from `lambert(..., max_revs=n_revs)`; call `izzo2015`/`gooding1990` with `M=n_revs`, `low_path=(branch != "high")`. Keep the default path (`n_revs=0`) identical to today.

- [ ] **Step 4: Run to verify pass + no single-rev regression**

Run: `uv run pytest tests/core/test_lambert.py tests/core/test_lambert_multirev.py -v`
Expected: PASS (the three existing single-rev crosscheck gates still green).

- [ ] **Step 5: Commit** (hold unless asked)

```bash
uv run ruff check src/cyclerfinder/core/lambert.py tests/core/test_lambert_multirev.py
uv run ruff format --check src/cyclerfinder/core/lambert.py tests/core/test_lambert_multirev.py
git add src/cyclerfinder/core/lambert.py tests/core/test_lambert_multirev.py
git commit -m "core/lambert: generalise lambert_crosscheck to multi-rev"
```

---

### Task 6: Thread multi-rev through real-ephemeris construction

**Files:**
- Modify: `src/cyclerfinder/verify/real_closure.py` (`construct_real_ephemeris_cycler`, `:516-594`)
- Test: `tests/verify/test_real_closure.py` (`test_construct_raises_on_multi_rev_leg` → repurpose)

Delete the `n_revs>0` raise; for each leg request `lambert(..., max_revs=n_revs)` and select the catalogue branch (default `"low"`, optional leg `branch` override); populate `Leg.n_revs`/`branch` from the chosen solution.

- [ ] **Step 1: Repurpose the existing raise-test into a build-test**

In `tests/verify/test_real_closure.py`, replace `test_construct_raises_on_multi_rev_leg` (`:506-524`):

```python
def test_construct_builds_multi_rev_leg(astropy_ephem: Ephemeris) -> None:
    """A leg with n_revs=1 now builds a multi-rev Leg instead of raising.

    Uses a long tof so revolution 1 is feasible (a 200-day E->M leg is too
    short for n=1; widen to ~780 d so t_min(1) < tof).
    """
    entry = {
        "id": "synthetic-multirev",
        "bodies": ["E", "M"],
        "legs": [{"from": "E", "to": "M", "tof_days": 780.0, "n_revs": 1}],
        "period": {"years": 2.135},
    }
    cyc = construct_real_ephemeris_cycler(entry, astropy_ephem, 0.0)
    assert cyc.legs[0].n_revs == 1
    assert cyc.legs[0].branch in ("low", "high")
```

(Verify the 780 d feasibility assumption when running; if `t_min(1)` still exceeds it for the real E→M geometry, bump `tof_days` until revolution 1 is feasible — the test documents a *feasible* multi-rev leg, the exact tof is not load-bearing.)

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/verify/test_real_closure.py::test_construct_builds_multi_rev_leg -v`
Expected: FAIL (`MultiRevLambertRequiredError` raised).

- [ ] **Step 3: Implement the threading**

In `construct_real_ephemeris_cycler`, replace the per-leg block (`:518-548`):

```python
    leg_vels: list[tuple[NDArray[np.float64], NDArray[np.float64], int, str]] = []
    for j, leg in enumerate(legs_meta):
        n_revs = int(leg.get("n_revs", 0) or 0)
        requested_branch = leg.get("branch") or ("single" if n_revs == 0 else "low")
        r1, _v1_planet = planet_states[j]
        r2, _v2_planet = planet_states[j + 1]
        tof = encounter_times[j + 1] - encounter_times[j]
        try:
            sols = lambert(r1, r2, tof, mu=MU_SUN_KM3_S2, max_revs=n_revs)
        except (LambertConvergenceError, LambertGeometryError) as exc:
            raise RealClosureConstructionError(
                catalogue_id=cat_id if cat_id is None else str(cat_id),
                leg_index=j,
                cause=exc,
            ) from exc
        chosen = next(
            (s for s in sols if s.n_revs == n_revs and s.branch == requested_branch), None
        )
        if chosen is None:
            raise RealClosureConstructionError(
                catalogue_id=cat_id if cat_id is None else str(cat_id),
                leg_index=j,
                cause=ValueError(
                    f"no Lambert solution n_revs={n_revs} branch={requested_branch!r}; "
                    f"available={[(s.n_revs, s.branch) for s in sols]}"
                ),
            )
        leg_vels.append(
            (
                np.asarray(chosen.v1, dtype=np.float64),
                np.asarray(chosen.v2, dtype=np.float64),
                chosen.n_revs,
                chosen.branch,
            )
        )
```

Update the Encounter-building unpack (`:556,560,564-565`) to ignore the two new tuple fields (`leg_vels[0][0]` etc. still index velocities). Update the Leg-building loop (`:580-594`) to read `n_revs`/`branch` from the tuple instead of hardcoding `0`/`"single"`. Remove the now-unused `MultiRevLambertRequiredError` import/raise. Update the docstring: remove the `MultiRevLambertRequiredError` Raises entry, document the branch default.

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/verify/test_real_closure.py::test_construct_builds_multi_rev_leg -v`
Expected: PASS.

- [ ] **Step 5: Commit** (hold unless asked)

```bash
uv run ruff check src/cyclerfinder/verify/real_closure.py tests/verify/test_real_closure.py
uv run ruff format --check src/cyclerfinder/verify/real_closure.py tests/verify/test_real_closure.py
git add src/cyclerfinder/verify/real_closure.py tests/verify/test_real_closure.py
git commit -m "verify/real_closure: build multi-rev legs instead of skipping"
```

---

### Task 7: Allow multi-rev legs through the V1 crosscheck

**Files:**
- Modify: `src/cyclerfinder/verify/crosscheck.py` (`crosscheck_leg` `:127-141`, `crosscheck_cycler` `:180-185`)
- Test: `tests/verify/test_crosscheck.py` (add a multi-rev leg case)

- [ ] **Step 1: Write the failing test** (locate the existing crosscheck test file first; add a case that builds a `Cycler` with one `n_revs=1` leg and asserts `crosscheck_leg` returns a `passed=True` result rather than raising). If `tests/verify/test_crosscheck.py` does not exist, find the crosscheck tests via `grep -rl crosscheck_leg tests/` and add the case there.

```python
def test_crosscheck_leg_multirev_passes(astropy_ephem) -> None:
    """A multi-rev leg crosschecks against lamberthub with matching M/low_path."""
    from cyclerfinder.verify.real_closure import construct_real_ephemeris_cycler
    from cyclerfinder.verify.crosscheck import crosscheck_leg

    entry = {
        "id": "synthetic-multirev",
        "bodies": ["E", "M"],
        "legs": [{"from": "E", "to": "M", "tof_days": 780.0, "n_revs": 1}],
        "period": {"years": 2.135},
    }
    cyc = construct_real_ephemeris_cycler(entry, astropy_ephem, 0.0)
    res = crosscheck_leg(cyc.legs[0], cyc, astropy_ephem, leg_index=0)
    assert res.passed, res.max_diff_mps
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest <crosscheck-test-file>::test_crosscheck_leg_multirev_passes -v`
Expected: FAIL (`NotImplementedError` raised at `crosscheck.py:128`).

- [ ] **Step 3: Implement**

In `crosscheck_leg`, remove the `n_revs>0` `NotImplementedError` (`:128-132`). Select the in-house solution matching the leg: `mine = next(s for s in lambert(r1, r2, tof_sec, mu=mu, prograde=True, max_revs=leg.n_revs) if s.n_revs == leg.n_revs and s.branch == leg.branch)`. Call `izzo2015`/`gooding1990` with `M=leg.n_revs, low_path=(leg.branch != "high")`. In `crosscheck_cycler`, delete the `if leg.n_revs > 0: continue` skip (`:182-183`) so multi-rev legs are crosschecked too; update the docstring.

- [ ] **Step 4: Run to verify pass + no single-rev regression**

Run: `uv run pytest <crosscheck-test-file> -v`
Expected: PASS (existing single-rev crosscheck tests still green).

- [ ] **Step 5: Commit** (hold unless asked)

```bash
uv run ruff check src/cyclerfinder/verify/crosscheck.py <crosscheck-test-file>
uv run ruff format --check src/cyclerfinder/verify/crosscheck.py <crosscheck-test-file>
git add src/cyclerfinder/verify/crosscheck.py <crosscheck-test-file>
git commit -m "verify/crosscheck: crosscheck multi-rev legs vs lamberthub"
```

---

### Task 8: Retire the `v3-skipped-multirev` gate & decide the exception's fate

**Files:**
- Modify: `src/cyclerfinder/verify/real_closure.py` (`verify_real_closure` `:746-751`, docstring `:659-660,204`, `__all__` `:824`, EXPECTED_SKIPS `:241-243`)
- Test: `tests/verify/test_real_closure.py` (`:686` `v3-skipped-multirev` test)

With construction now building multi-rev legs, `MultiRevLambertRequiredError` is never raised. Decide: **delete** the exception + its `except` branch, OR keep it as a guarded fallback. Recommendation: **delete** the `except MultiRevLambertRequiredError` branch in `verify_real_closure` and the class; a multi-rev leg now either constructs (and closure is measured normally) or raises `RealClosureConstructionError` like any other Lambert failure. Keep `MultiRevLambertRequiredError` exported as a deprecated alias only if other repos import it — grep first.

- [ ] **Step 1: Grep for external dependence on the symbol**

Run: `grep -rn "MultiRevLambertRequiredError\|v3-skipped-multirev" src/ tests/`
Expected: enumerate all references so none dangle after removal.

- [ ] **Step 2: Repurpose the skip test**

The `n_revs=1` → `v3-skipped-multirev` test (`:686`) becomes: a feasible multi-rev catalogue entry yields a `RealClosureResult` whose `v3_status` is no longer `"v3-skipped-multirev"` (it either closes or reports a real drift/construction outcome). Rewrite the assertion accordingly.

- [ ] **Step 3: Run to verify failure** (old test still asserts the skip status)

Run: `uv run pytest tests/verify/test_real_closure.py -k multirev -v`
Expected: FAIL on the old skip assertion.

- [ ] **Step 4: Implement removal**

Delete the `except MultiRevLambertRequiredError → v3-skipped-multirev` branch in `verify_real_closure` (`:746-751`); delete the class (`:118-131`) and its `__all__` entry; remove `s1l1-2syn-em-cpom` from `EXPECTED_SKIPS` if its skip reason was multi-rev (`:241-243`) — **but** if S1L1 now constructs yet fails to *close* (likely, since real closure needs the full 4-encounter geometry / M-N), re-add it to EXPECTED_SKIPS with an honest new reason ("constructs on DE440 but does not close as a 3-encounter slice; full S1L1 closure is M-N/M-ED") rather than leaving a failing test. Update module docstring references.

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/verify/test_real_closure.py -v`
Expected: PASS (the S1L1 XFAIL either flips to passing if it closes, or is consciously re-marked XFAIL with a non-multirev reason — see Task 9).

- [ ] **Step 6: Commit** (hold unless asked)

```bash
uv run ruff check src/cyclerfinder/verify/real_closure.py tests/verify/test_real_closure.py
uv run ruff format --check src/cyclerfinder/verify/real_closure.py tests/verify/test_real_closure.py
git add src/cyclerfinder/verify/real_closure.py tests/verify/test_real_closure.py
git commit -m "verify/real_closure: retire multi-rev skip gate"
```

---

### Task 9: Bring in the golden orbit — attempt S1L1 closure on DE440

**Files:**
- Modify: `tests/verify/test_real_closure.py` (`test_2syn_em_cpom_periodic_over_2_cycles_astropy` `:266-289`)

This is the payoff task: with multi-rev Lambert landed, attempt to actually construct + verify the published S1L1 2-syn E-M family (`s1l1-2syn-em-cpom`) on real ephemeris and either **promote it to a passing golden gate** or document honestly why it still cannot close (handing the residual to M-N/M-ED).

- [ ] **Step 1: Construct + verify S1L1 and observe the outcome**

Run the currently-XFAIL gate **without** the xfail marker locally (temporarily) to see what actually happens now:

Run: `uv run pytest tests/verify/test_real_closure.py::test_2syn_em_cpom_periodic_over_2_cycles_astropy -v -rX`
Observe: does it now construct and close (`result.closes`, `max_drift_km < REAL_DRIFT_TOLERANCE_KM`), or construct-but-not-close, or still error?

- [ ] **Step 2: Branch on the result**

- **If it closes:** remove the `@pytest.mark.xfail` marker — the gate becomes a passing golden orbit. Add a comment tagging it as a sourced golden anchor (priority date 2002-08-14, source AIAA 2002-4420 per the existing catalogue entry). This is the second golden orbit after Aldrin.
- **If it constructs but drift exceeds tolerance, or only 2 of the published 4 encounters are tabulated** (the known S1L1 data_gap — `s1l1-2syn-em-cpom` is the E-E-M-M / S1L1 4-encounter topology, not a 3-encounter slice): keep `xfail(strict=True)` but **rewrite the reason** to the real, post-multi-rev blocker (e.g. "multi-rev Lambert lands, S1L1 now constructs on DE440; full closure needs the 4-encounter E-E-M-M geometry + ephemeris convergence — M-N/M-ED"). Do NOT leave the stale "multi-rev blocker" reason. Flip `strict=False`→`strict=True` only if construction is now deterministic.

- [ ] **Step 3: Run the verify suite**

Run: `uv run pytest tests/verify/ -v`
Expected: green (S1L1 either passes as golden, or xfails with an honest updated reason).

- [ ] **Step 4: Commit** (hold unless asked)

```bash
uv run ruff check tests/verify/test_real_closure.py
uv run ruff format --check tests/verify/test_real_closure.py
git add tests/verify/test_real_closure.py
git commit -m "verify: attempt S1L1 closure on DE440 post multi-rev Lambert"
```

---

### Task 10: Full-suite regression + lint + type gate

**Files:** none (verification only)

- [ ] **Step 1: Full fast suite**

Run: `uv run pytest -q`
Expected: no new failures; the multi-rev XFAILs that were blocked purely on the solver are resolved (passing or honestly re-reasoned). The M5-optimiser XFAILs (`test_2syn_em_rediscovers_5_65_kms_earth` etc.) remain XFAIL — those are M-ED, not M-L.

- [ ] **Step 2: Lint + format + types**

Run: `uv run ruff check . && uv run ruff format --check . && uv run mypy src/`
Expected: clean.

- [ ] **Step 3: Final commit** (hold unless asked)

Only if the user has asked to commit. Otherwise report the green suite and stop.

---

## Self-review

- **Spec coverage:** M-L roadmap objective = "`lambert(..., max_revs=N)` returns full solution set; construction selects `(n_revs, branch)`; the `n_revs>0` rejection gates removed." → Tasks 1-4 (solver), 5 (crosscheck helper), 6 (construction selection — note `construct.py` general path already selects; the real-ephemeris path is Task 6), 7-8 (gate removal), 9 (golden orbit payoff). Covered.
- **Placeholder scan:** every code step has concrete code; the only deferred specifics are the empirically-verified branch-label mapping (Task 4 is the arbiter) and the exact feasible `tof_days` for the synthetic multi-rev leg (Task 6 Step 1 says to bump until feasible). Both are flagged as test-arbitrated, not hand-waved.
- **Type consistency:** `LambertSolution(n_revs:int, branch:str, v1:Vec3, v2:Vec3)` used identically across tasks; `_solve_uv_branch`/`_min_time_of_revolution`/`_velocities_from_z` signatures consistent between definition (Tasks 1-3) and call sites.
- **Discipline:** every golden/crosscheck EXPECTED value is `lamberthub` (independent solver) — never our own output. Commits are held pending explicit user request. ruff+mypy gate before each commit. No `--no-verify`, no AI attribution.

## Dependencies & handoff

M-L has no upstream deps. It unblocks **M-N** (the E-E-M-M intermediate loop legs need multi-rev to construct) and feeds **M-ED** (multi-rev crosscheck is reused by ephemeris-mode verification). After M-L: detail M-N (ready), and run the design pass for M-3D / M-ED.
