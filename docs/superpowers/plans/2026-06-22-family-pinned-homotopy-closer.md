# Family-pinned Penalty Homotopy Closer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Drive the full STM Earth-Mars cycler shooter toward the *published* V∞ family by adding an opt-in V∞-anchor penalty to the corrector, ramping it to zero over a homotopy ladder, and running a detached/checkpointed batch whose verdict (PROPOSED V0→V1 or stronger negative) is read only from the unpenalized solve.

**Architecture:** A V∞-anchor penalty is appended to the existing `defect_residual` and `_stm_jacobian` behind opt-in args (default path byte-identical). A small ladder/warm-start driver (`family_pinned_shoot`) sweeps the penalty weight `λv: W → 0`, warm-starting each rung from the previous. A sibling batch script runs it over the 4 descriptor rows × best-phase epochs, checkpointed, and records the `λv=0` verdict. HELD throughout — no catalogue/validate writeback.

**Tech Stack:** Python 3.11, numpy, scipy `least_squares(method="lm")`, REBOUND/IAS15 (via `nbody/propagator.py`), `uv` for env/test, `pytest` (`@pytest.mark.slow` for n-body solves), ruff + mypy pre-commit.

**Concurrency note (binding):** Another agent is editing `data/catalogue.yaml`, `data/catalogue.schema.json`, `src/cyclerfinder/data/`, and `tests/data/` + `tests/search/test_resonance_network.py` in the same working tree. This plan touches ONLY `src/cyclerfinder/nbody/shooter.py`, `src/cyclerfinder/search/family_pinned_shoot.py` (new), `tests/nbody/test_shooter_family_pinned.py` (new), `tests/search/test_family_pinned_shoot.py` (new), and `scripts/shooter_family_pinned_batch.py` (new). **Every commit MUST use explicit pathspecs (never `git add -A`/`git add .`).** Re-run `git status` before committing and stage only this plan's files.

---

## File Structure

- **Modify** `src/cyclerfinder/nbody/shooter.py` — add `vinf_anchors` / `vinf_weight` opt-in args to `defect_residual`, `_stm_jacobian`, and `shoot` (Tasks 1–3). The penalty is a basin-selector; the default (no anchors / zero weight) path is byte-identical.
- **Create** `src/cyclerfinder/search/family_pinned_shoot.py` — the ladder/warm-start homotopy driver + `FamilyPinnedResult` (Task 4).
- **Create** `scripts/shooter_family_pinned_batch.py` — detached/checkpointed batch + verdict (Task 5).
- **Create** `tests/nbody/test_shooter_family_pinned.py` — penalty residual + Jacobian + shoot tests (Tasks 1–3).
- **Create** `tests/search/test_family_pinned_shoot.py` — driver tests (Task 4).

The penalty math (used in Tasks 1 and 2):
- Per encounter node `i` with body `b = sequence[i]` that has an anchor `a = vinf_anchors[b]`:
  `v∞_i = node_states[i][3:] − v_planet(b, epochs[i])`, `|v∞_i| = ‖v∞_i‖`.
  Penalty residual row: `sqrt(vinf_weight) · (|v∞_i| − a)`.
  Penalty Jacobian row (wrt node `i`'s 6-state): `sqrt(vinf_weight) · [0,0,0, v̂∞_i]` where `v̂∞_i = v∞_i/|v∞_i|` (zero row if `|v∞_i| < 1e-12`).
- Penalty rows are appended AFTER the existing leg-defect + hinge + wrap blocks, in node order, only for nodes whose body is in `vinf_anchors`.

---

## Task 1: V∞-anchor penalty in `defect_residual`

**Files:**
- Modify: `src/cyclerfinder/nbody/shooter.py` (the `defect_residual` function)
- Test: `tests/nbody/test_shooter_family_pinned.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/nbody/test_shooter_family_pinned.py`:

```python
from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import (  # noqa: E402
    ShootingSeed,
    defect_residual,
)


def _two_body_seed() -> tuple[ShootingSeed, Ephemeris]:
    """A 3-node Sun-only seed sampled from one continuous two-body arc.

    Mirrors the fixture in tests/nbody/test_shooter_stm_jacobian.py: cheap,
    deterministic, bodies=() so the penalty (which depends only on node velocity
    vs planet velocity from the ephemeris) is exercised without perturber cost.
    """
    ephem = Ephemeris("circular")
    prop = RestrictedNBody("rebound")
    r0 = np.array([1.30e8, 0.0, 0.0])
    v0 = np.array([0.0, 26.0, 1.5])
    day = 86400.0
    t_nodes = [0.0, 120.0 * day, 260.0 * day]
    states: list[np.ndarray] = [np.concatenate([r0, v0])]
    cur_r, cur_v, cur_t = r0, v0, 0.0
    for t1 in t_nodes[1:]:
        arc = prop.propagate(cur_r, cur_v, t0_sec=cur_t, t1_sec=t1, bodies=(), accuracy=1e-11)
        states.append(np.concatenate([arc.r_km, arc.v_km_s]))
        cur_r, cur_v, cur_t = arc.r_km, arc.v_km_s, t1
    zero = np.zeros(3)
    seed = ShootingSeed(
        node_states=states,
        epochs=t_nodes,
        tofs=[120.0, 140.0],
        sequence=("E", "M", "E"),
        slack_leg=1,
        period_days=260.0,
        vinf_in=[zero, zero, zero],
        vinf_out=[zero, zero, zero],
    )
    return seed, ephem


def test_vinf_penalty_off_is_identical() -> None:
    seed, ephem = _two_body_seed()
    base = defect_residual(seed, ephem=ephem, bodies=())
    none_set = defect_residual(seed, ephem=ephem, bodies=(), vinf_anchors=None, vinf_weight=0.0)
    np.testing.assert_array_equal(base, none_set)
    # anchors present but zero weight -> still identical (no rows appended)
    zero_w = defect_residual(
        seed, ephem=ephem, bodies=(), vinf_anchors={"E": 5.0, "M": 6.0}, vinf_weight=0.0
    )
    np.testing.assert_array_equal(base, zero_w)


def test_vinf_penalty_rows_values() -> None:
    seed, ephem = _two_body_seed()
    w = 4.0
    anchors = {"E": 5.0, "M": 6.0}
    base = defect_residual(seed, ephem=ephem, bodies=())
    res = defect_residual(seed, ephem=ephem, bodies=(), vinf_anchors=anchors, vinf_weight=w)
    pen = res[len(base) :]
    # sequence is E, M, E -> all three nodes carry an anchor -> 3 penalty rows
    assert pen.shape == (3,)
    sw = float(np.sqrt(w))
    expected = []
    for i, body in enumerate(seed.sequence):
        _, v_pl = ephem.state(body, seed.epochs[i])
        mag = float(np.linalg.norm(np.asarray(seed.node_states[i][3:]) - np.asarray(v_pl)))
        expected.append(sw * (mag - anchors[body]))
    np.testing.assert_allclose(pen, np.asarray(expected), rtol=1e-12, atol=1e-9)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/nbody/test_shooter_family_pinned.py -q -o addopts="" -p no:cacheprovider`
Expected: FAIL — `defect_residual() got an unexpected keyword argument 'vinf_anchors'`.

- [ ] **Step 3: Add the penalty args + block to `defect_residual`**

In `src/cyclerfinder/nbody/shooter.py`, change the `defect_residual` signature from:

```python
def defect_residual(
    seed: ShootingSeed,
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
) -> NDArray[np.float64]:
```

to (add the two new keyword args at the end):

```python
def defect_residual(
    seed: ShootingSeed,
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
    vinf_anchors: Mapping[str, float] | None = None,
    vinf_weight: float = 0.0,
) -> NDArray[np.float64]:
```

Then, immediately BEFORE the final `return np.asarray(res, dtype=np.float64)`, insert the penalty block:

```python
    # 4. V∞-anchor penalty (family-pinned homotopy, #388). Opt-in: only when
    #    anchors are supplied AND the weight is positive. Each penalty row pulls a
    #    node's V∞ magnitude toward the SOURCED anchor for its body; the homotopy
    #    driver ramps vinf_weight to zero so the recorded V∞ emerges from the
    #    unpenalized (weight=0) solve. Appended after leg/hinge/wrap, in node order.
    if vinf_anchors and vinf_weight > 0.0:
        sw = float(np.sqrt(vinf_weight))
        for i, body in enumerate(seed.sequence):
            anchor = vinf_anchors.get(body)
            if anchor is None:
                continue
            s_i = seed.node_states[i]
            _, v_pl = ephem.state(body, seed.epochs[i])
            vinf_mag = float(
                np.linalg.norm(np.asarray(s_i[3:], dtype=np.float64) - np.asarray(v_pl, dtype=np.float64))
            )
            res.append(sw * (vinf_mag - float(anchor)))
```

Confirm `Mapping` is already imported at the top of the file (it is:
`from collections.abc import Callable, Mapping, Sequence`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/nbody/test_shooter_family_pinned.py -q -o addopts="" -p no:cacheprovider`
Expected: PASS (2 passed).

- [ ] **Step 5: Lint, typecheck, commit (pathspec-scoped)**

```bash
uv run ruff check src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
uv run ruff format src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
uv run mypy src/cyclerfinder/nbody/shooter.py
git status   # confirm only shooter.py + the new test are staged below
git add src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
git commit -m "nbody/#388: V_inf-anchor penalty rows in defect_residual (opt-in)"
```

---

## Task 2: V∞-anchor penalty rows in `_stm_jacobian`

**Files:**
- Modify: `src/cyclerfinder/nbody/shooter.py` (the `_stm_jacobian` function)
- Test: `tests/nbody/test_shooter_family_pinned.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/nbody/test_shooter_family_pinned.py` (add these imports to the existing import block at the top: `_fd_jacobian`, `_seed_with_states`, `_serial_columns`, `_states_to_x`, `_x_to_states`, `_stm_jacobian`):

```python
from cyclerfinder.nbody.shooter import (  # noqa: E402  (extend the existing import)
    ShootingSeed,
    _fd_jacobian,
    _seed_with_states,
    _serial_columns,
    _states_to_x,
    _stm_jacobian,
    _x_to_states,
    defect_residual,
)


@pytest.mark.slow
def test_stm_jacobian_with_penalty_matches_fd() -> None:
    """The augmented STM Jacobian (incl. penalty rows) matches the FD oracle."""
    seed, ephem = _two_body_seed()
    bodies: tuple[str, ...] = ()
    anchors = {"E": 5.0, "M": 6.0}
    w = 4.0
    x0 = _states_to_x(seed.node_states)

    def resid(x: np.ndarray) -> np.ndarray:
        trial = _seed_with_states(seed, _x_to_states(x, len(seed.sequence)))
        return defect_residual(
            trial, ephem=ephem, bodies=bodies, accuracy=1e-11, vinf_anchors=anchors, vinf_weight=w
        )

    f0 = resid(x0)
    fd = _fd_jacobian(resid, x0, f0, column_eval=_serial_columns)
    stm = _stm_jacobian(
        seed, x0, ephem=ephem, bodies=bodies, accuracy=1e-11, vinf_anchors=anchors, vinf_weight=w
    )
    assert stm.shape == fd.shape
    rel = np.linalg.norm(stm - fd) / np.linalg.norm(fd)
    print(f"penalty STM Jacobian vs FD rel = {rel:.3e}")
    assert rel < 5e-3, f"penalty STM Jacobian vs FD rel={rel}"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/nbody/test_shooter_family_pinned.py::test_stm_jacobian_with_penalty_matches_fd -q -o addopts="" -p no:cacheprovider`
Expected: FAIL — `_stm_jacobian() got an unexpected keyword argument 'vinf_anchors'`.

- [ ] **Step 3: Add the penalty args + rows to `_stm_jacobian`**

In `src/cyclerfinder/nbody/shooter.py`, change the `_stm_jacobian` signature from:

```python
def _stm_jacobian(
    seed: ShootingSeed,
    x: NDArray[np.float64],
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
) -> NDArray[np.float64]:
```

to:

```python
def _stm_jacobian(
    seed: ShootingSeed,
    x: NDArray[np.float64],
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
    vinf_anchors: Mapping[str, float] | None = None,
    vinf_weight: float = 0.0,
) -> NDArray[np.float64]:
```

Find the row-count setup (currently):

```python
    n_leg = (n - 1) * _STATE_DIM
    n_hinge = max(0, n - 2)
    n_rows = n_leg + n_hinge + _STATE_DIM
    n_cols = n * _STATE_DIM
    jac = np.zeros((n_rows, n_cols), dtype=np.float64)
    eye6 = np.eye(_STATE_DIM)
```

Replace it with a version that reserves penalty rows (matching `defect_residual`'s order — leg, hinge, wrap, then penalty):

```python
    n_leg = (n - 1) * _STATE_DIM
    n_hinge = max(0, n - 2)
    # Penalty node indices (those whose body has an anchor), in node order — must
    # match the residual's penalty-row order in defect_residual.
    pen_nodes: list[int] = []
    if vinf_anchors and vinf_weight > 0.0:
        pen_nodes = [i for i, b in enumerate(seed.sequence) if b in vinf_anchors]
    n_pen = len(pen_nodes)
    n_rows = n_leg + n_hinge + _STATE_DIM + n_pen
    n_cols = n * _STATE_DIM
    jac = np.zeros((n_rows, n_cols), dtype=np.float64)
    eye6 = np.eye(_STATE_DIM)
```

The existing leg / wrap fills are unchanged (they index by `n_leg`/`n_hinge` and the
wrap slice `slice(n_leg + n_hinge, n_leg + n_hinge + _STATE_DIM)` — if the current
code uses `n_rows` to locate the wrap slice, change the wrap slice to the explicit
`slice(n_leg + n_hinge, n_leg + n_hinge + _STATE_DIM)` so the trailing penalty rows
are not overwritten). Confirm the wrap fill reads:

```python
    wrap = slice(n_leg + n_hinge, n_leg + n_hinge + _STATE_DIM)
    jac[wrap, 0:_STATE_DIM] = -eye6
    jac[wrap, (n - 1) * _STATE_DIM : n * _STATE_DIM] = eye6
```

Then, immediately before `return jac`, add the penalty Jacobian rows:

```python
    # Penalty rows: d(sqrt(w)*|v∞_i|)/d(node_i velocity) = sqrt(w) * v̂∞_i, zero
    # on the position block and on every other node. Matches the residual's
    # penalty rows (defect_residual block 4). A near-zero-V∞ node -> zero row.
    if n_pen:
        sw = float(np.sqrt(vinf_weight))
        base = n_leg + n_hinge + _STATE_DIM
        for k, i in enumerate(pen_nodes):
            s_i = states[i]
            _, v_pl = ephem.state(seed.sequence[i], seed.epochs[i])
            dv = np.asarray(s_i[3:], dtype=np.float64) - np.asarray(v_pl, dtype=np.float64)
            mag = float(np.linalg.norm(dv))
            if mag > 1e-12:
                vhat = dv / mag
                jac[base + k, i * _STATE_DIM + 3 : i * _STATE_DIM + 6] = sw * vhat
```

(`states = _x_to_states(x, n)` already exists earlier in the function — reuse it.)

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/nbody/test_shooter_family_pinned.py::test_stm_jacobian_with_penalty_matches_fd -s -q -o addopts="" -p no:cacheprovider`
Expected: PASS; printed `penalty STM Jacobian vs FD rel` well below 5e-3.

- [ ] **Step 5: Lint, typecheck, commit (pathspec-scoped)**

```bash
uv run ruff check src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
uv run ruff format src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
uv run mypy src/cyclerfinder/nbody/shooter.py
git status
git add src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
git commit -m "nbody/#388: V_inf-anchor penalty rows in STM Jacobian (gravity-gradient analogue)"
```

---

## Task 3: Thread the penalty through `shoot`

**Files:**
- Modify: `src/cyclerfinder/nbody/shooter.py` (the `shoot` function)
- Test: `tests/nbody/test_shooter_family_pinned.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/nbody/test_shooter_family_pinned.py` (extend the shooter import with `shoot`):

```python
from cyclerfinder.nbody.shooter import shoot  # noqa: E402  (add to the import block)


def test_shoot_penalty_off_matches_plain() -> None:
    """vinf_weight=0 leaves shoot byte-for-byte unchanged."""
    seed, ephem = _two_body_seed()
    plain = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=20, jacobian="stm")
    off = shoot(
        seed,
        ephem=ephem,
        bodies=(),
        accuracy=1e-11,
        max_nfev=20,
        jacobian="stm",
        vinf_anchors={"E": 5.0, "M": 6.0},
        vinf_weight=0.0,
    )
    assert plain.defect_norm == off.defect_norm
    for a, b in zip(plain.corrected_states, off.corrected_states, strict=True):
        np.testing.assert_array_equal(a, b)


def test_shoot_penalty_rejects_parallel_fd() -> None:
    """A positive penalty with the parallel FD path is unsupported -> ValueError."""
    seed, ephem = _two_body_seed()
    with pytest.raises(ValueError, match="vinf penalty"):
        shoot(
            seed,
            ephem=ephem,
            bodies=(),
            accuracy=1e-11,
            max_nfev=5,
            n_jobs=4,
            vinf_anchors={"E": 5.0},
            vinf_weight=1.0,
        )


@pytest.mark.slow
def test_shoot_penalty_biases_vinf_toward_anchor() -> None:
    """A strong penalty pulls a node's corrected V∞ toward the anchor."""
    seed, ephem = _two_body_seed()
    # Natural (unpenalized) V∞ at node 0 from the fixture.
    base = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=40, jacobian="stm")
    base_e = base.vinf_per_encounter_kms[0]
    target = base_e + 5.0  # pull the Earth-node V∞ 5 km/s away from natural
    pinned = shoot(
        seed,
        ephem=ephem,
        bodies=(),
        accuracy=1e-11,
        max_nfev=40,
        jacobian="stm",
        vinf_anchors={"E": target},
        vinf_weight=50.0,
    )
    # The penalized solve's Earth-node V∞ moves toward the target vs the baseline.
    assert abs(pinned.vinf_per_encounter_kms[0] - target) < abs(base_e - target)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/nbody/test_shooter_family_pinned.py -q -o addopts="" -p no:cacheprovider -k "penalty_off or rejects_parallel"`
Expected: FAIL — `shoot() got an unexpected keyword argument 'vinf_anchors'`.

- [ ] **Step 3: Add the penalty args + wiring + guard to `shoot`**

In `src/cyclerfinder/nbody/shooter.py`, extend the `shoot` signature (add the two args after `progress`):

```python
    progress: Callable[[str, int, float, float], None] | None = None,
    vinf_anchors: Mapping[str, float] | None = None,
    vinf_weight: float = 0.0,
) -> ShootResult:
```

Immediately after `n = len(seed.sequence)` at the top of the body, add the guard:

```python
    if vinf_anchors and vinf_weight > 0.0 and n_jobs > 1:
        raise ValueError(
            "vinf penalty (family-pinned homotopy) is supported only with n_jobs=1 "
            "(the parallel FD worker path does not carry anchors); use jacobian='stm'."
        )
```

Update `residual_of_x` to pass the anchors through:

```python
    def residual_of_x(x: NDArray[np.float64]) -> NDArray[np.float64]:
        states = _x_to_states(x, n)
        trial = _seed_with_states(seed, states)
        return defect_residual(
            trial,
            ephem=ephem,
            bodies=bodies,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
            vinf_anchors=vinf_anchors,
            vinf_weight=vinf_weight,
        )
```

And update the STM Jacobian closure `jac_of_x_stm` to pass them:

```python
        def jac_of_x_stm(x: NDArray[np.float64]) -> NDArray[np.float64]:
            t_eval = _time.monotonic()
            jac = _stm_jacobian(
                seed,
                x,
                ephem=ephem,
                bodies=bodies,
                accuracy=accuracy,
                max_wall_sec=max_wall_sec,
                vinf_anchors=vinf_anchors,
                vinf_weight=vinf_weight,
            )
            if progress is not None:
                _j_count[0] += 1
                progress("J", _j_count[0], _last_norm[0], _time.monotonic() - t_eval)
            return jac
```

Note: the `least_squares` acceptance / `converged` logic uses only the leg-defect
block (`final_res[:n_leg_defects]`), which is unchanged by appended penalty rows —
leave it as is. The recorded `vinf_per_encounter_kms` is computed from
`corrected_states` (the emerged value), independent of the penalty — leave as is.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/nbody/test_shooter_family_pinned.py -q -o addopts="" -p no:cacheprovider`
Expected: PASS (all, incl. the slow bias test).

- [ ] **Step 5: Confirm the existing STM tests still pass (no regression)**

Run: `uv run pytest tests/nbody/test_shooter_stm_jacobian.py tests/nbody/test_propagator_stm.py -q -o addopts="" -p no:cacheprovider`
Expected: PASS (6 passed).

- [ ] **Step 6: Lint, typecheck, commit (pathspec-scoped)**

```bash
uv run ruff check src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
uv run ruff format src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
uv run mypy src/cyclerfinder/nbody/shooter.py
git status
git add src/cyclerfinder/nbody/shooter.py tests/nbody/test_shooter_family_pinned.py
git commit -m "nbody/#388: thread V_inf penalty through shoot (n_jobs=1 guard, default path unchanged)"
```

---

## Task 4: The ladder/warm-start homotopy driver

**Files:**
- Create: `src/cyclerfinder/search/family_pinned_shoot.py`
- Test: `tests/search/test_family_pinned_shoot.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/search/test_family_pinned_shoot.py`:

```python
from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import ShootingSeed  # noqa: E402
from cyclerfinder.search.family_pinned_shoot import (  # noqa: E402
    FamilyPinnedResult,
    family_pinned_shoot,
)


def _two_body_seed() -> tuple[ShootingSeed, Ephemeris]:
    ephem = Ephemeris("circular")
    prop = RestrictedNBody("rebound")
    r0 = np.array([1.30e8, 0.0, 0.0])
    v0 = np.array([0.0, 26.0, 1.5])
    day = 86400.0
    t_nodes = [0.0, 120.0 * day, 260.0 * day]
    states: list[np.ndarray] = [np.concatenate([r0, v0])]
    cur_r, cur_v, cur_t = r0, v0, 0.0
    for t1 in t_nodes[1:]:
        arc = prop.propagate(cur_r, cur_v, t0_sec=cur_t, t1_sec=t1, bodies=(), accuracy=1e-11)
        states.append(np.concatenate([arc.r_km, arc.v_km_s]))
        cur_r, cur_v, cur_t = arc.r_km, arc.v_km_s, t1
    zero = np.zeros(3)
    seed = ShootingSeed(
        node_states=states,
        epochs=t_nodes,
        tofs=[120.0, 140.0],
        sequence=("E", "M", "E"),
        slack_leg=1,
        period_days=260.0,
        vinf_in=[zero, zero, zero],
        vinf_out=[zero, zero, zero],
    )
    return seed, ephem


@pytest.mark.slow
def test_family_pinned_shoot_returns_result_and_trace() -> None:
    seed, ephem = _two_body_seed()
    res = family_pinned_shoot(
        seed,
        ephem=ephem,
        bodies=(),
        vinf_anchors={"E": 6.0, "M": 6.0},
        weight_ladder=(10.0, 1.0, 0.0),
        accuracy=1e-11,
        max_nfev=20,
    )
    assert isinstance(res, FamilyPinnedResult)
    # final rung is the unpenalized (weight=0) solve
    assert res.final_weight == 0.0
    assert res.final.sequence == seed.sequence
    # one trace entry per ladder rung
    assert [w for (w, _d, _v) in res.trace] == [10.0, 1.0, 0.0]
    # anchor_retention is finite
    assert np.isfinite(res.anchor_retention_kms)


@pytest.mark.slow
def test_family_pinned_shoot_retains_set_anchor() -> None:
    """If the anchor IS the seed's natural node V∞, the λv->0 solve retains it."""
    seed, ephem = _two_body_seed()
    # natural Earth-node V∞ of the seed (node 0, body E)
    _, v_pl = ephem.state("E", seed.epochs[0])
    natural_e = float(np.linalg.norm(np.asarray(seed.node_states[0][3:]) - np.asarray(v_pl)))
    res = family_pinned_shoot(
        seed,
        ephem=ephem,
        bodies=(),
        vinf_anchors={"E": natural_e},
        weight_ladder=(10.0, 1.0, 0.0),
        accuracy=1e-11,
        max_nfev=20,
    )
    # final (unpenalized) Earth-node V∞ stays near the natural/anchor value
    assert abs(res.final.vinf_per_encounter_kms[0] - natural_e) < 1.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/search/test_family_pinned_shoot.py -q -o addopts="" -p no:cacheprovider`
Expected: FAIL — `ModuleNotFoundError: cyclerfinder.search.family_pinned_shoot`.

- [ ] **Step 3: Implement the driver**

Create `src/cyclerfinder/search/family_pinned_shoot.py`:

```python
"""Family-pinned penalty homotopy closer (#388).

Drives the full STM multiple-shooting corrector (``nbody.shooter.shoot``) toward
the PUBLISHED V∞ family by ramping a V∞-anchor penalty weight from a calibrated
``W`` down to zero. Each rung warm-starts from the previous rung's corrected
states (the ``continuation.continuation_correct`` ladder pattern). The final
``weight == 0`` rung is the verdict solve: its emerged V∞ is the recorded value
(the penalty is a basin-selector ramped to zero, never the recorded number — the
golden-discipline requirement of ``feedback_golden_tests_sourced_only``).

See ``docs/superpowers/specs/2026-06-22-family-pinned-homotopy-closer-design.md``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from cyclerfinder.nbody.shooter import ShootResult, ShootingSeed, _seed_with_states, shoot

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris


@dataclass(frozen=True)
class FamilyPinnedResult:
    """Result of a family-pinned penalty homotopy run.

    ``final`` is the unpenalized (``weight == 0``) ShootResult — the verdict solve
    whose emerged V∞ is the recorded value. ``trace`` is one
    ``(weight, defect_norm, vinf_per_encounter_kms)`` tuple per ladder rung.
    ``anchor_retention_kms`` is the change in the best per-anchor V∞ residual from
    the FIRST penalized rung to the final ``weight == 0`` rung — small means V∞
    held near the published family as the penalty lifted; large means it snapped
    off-anchor (the stronger characterized negative).
    """

    final: ShootResult
    final_weight: float
    trace: list[tuple[float, float, list[float]]]
    anchor_retention_kms: float
    vinf_anchors: dict[str, float]


def _best_anchor_residual(
    vinf_per_encounter_kms: Sequence[float], vinf_anchors: Mapping[str, float]
) -> float:
    """Worst-over-anchors of the best-matching encounter V∞ (km/s)."""
    worst = 0.0
    for anchor in vinf_anchors.values():
        best = min((abs(v - float(anchor)) for v in vinf_per_encounter_kms), default=float("inf"))
        worst = max(worst, best)
    return worst


def family_pinned_shoot(
    seed: ShootingSeed,
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    vinf_anchors: Mapping[str, float],
    weight_ladder: Sequence[float] = (40.0, 10.0, 2.5, 0.5, 0.0),
    accuracy: float = 1e-9,
    max_nfev: int = 100,
    max_wall_sec: float = 30.0,
    progress: Callable[[str, int, float, float], None] | None = None,
) -> FamilyPinnedResult:
    """Ramp the V∞-anchor penalty from ``weight_ladder[0]`` down to (a final) 0.

    Each rung runs ``shoot(jacobian="stm", vinf_anchors=..., vinf_weight=w)``,
    warm-started from the previous rung's corrected states. ``weight_ladder`` MUST
    end at 0.0 (asserted) so the verdict solve is unpenalized. Returns a
    :class:`FamilyPinnedResult`.
    """
    ladder = [float(w) for w in weight_ladder]
    if not ladder or ladder[-1] != 0.0:
        raise ValueError("weight_ladder must be non-empty and end at 0.0 (the verdict rung)")
    anchors = {str(k): float(v) for k, v in vinf_anchors.items()}

    cur = seed
    trace: list[tuple[float, float, list[float]]] = []
    last: ShootResult | None = None
    first_penalized_resid: float | None = None

    for w in ladder:
        res = shoot(
            cur,
            ephem=ephem,
            bodies=bodies,
            accuracy=accuracy,
            max_nfev=max_nfev,
            max_wall_sec=max_wall_sec,
            jacobian="stm",
            vinf_anchors=anchors,
            vinf_weight=w,
            progress=progress,
        )
        vinf = list(res.vinf_per_encounter_kms)
        trace.append((w, res.defect_norm, vinf))
        if w > 0.0 and first_penalized_resid is None:
            first_penalized_resid = _best_anchor_residual(vinf, anchors)
        # warm-start the next rung from this rung's corrected states
        cur = _seed_with_states(cur, res.corrected_states)
        last = res

    assert last is not None
    final_resid = _best_anchor_residual(list(last.vinf_per_encounter_kms), anchors)
    retention = (
        abs(final_resid - first_penalized_resid)
        if first_penalized_resid is not None
        else float("nan")
    )
    return FamilyPinnedResult(
        final=last,
        final_weight=ladder[-1],
        trace=trace,
        anchor_retention_kms=retention,
        vinf_anchors=anchors,
    )


__all__ = ["FamilyPinnedResult", "family_pinned_shoot"]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/search/test_family_pinned_shoot.py -q -o addopts="" -p no:cacheprovider`
Expected: PASS (2 passed).

- [ ] **Step 5: Lint, typecheck, commit (pathspec-scoped)**

```bash
uv run ruff check src/cyclerfinder/search/family_pinned_shoot.py tests/search/test_family_pinned_shoot.py
uv run ruff format src/cyclerfinder/search/family_pinned_shoot.py tests/search/test_family_pinned_shoot.py
uv run mypy src/cyclerfinder/search/family_pinned_shoot.py
git status
git add src/cyclerfinder/search/family_pinned_shoot.py tests/search/test_family_pinned_shoot.py
git commit -m "search/#388: family-pinned penalty homotopy driver (ladder + warm-start)"
```

---

## Task 5: Detached/checkpointed batch + verdict

**Files:**
- Create: `scripts/shooter_family_pinned_batch.py`
- Test: (none — script; verified by a manual smoke command in Step 3)

- [ ] **Step 1: Implement the batch script**

Create `scripts/shooter_family_pinned_batch.py` (mirrors `scripts/shooter_russell_batch.py`'s resilience: per-(row,epoch) append+fsync checkpoint, `--resume`, heartbeat; but calls `family_pinned_shoot` and records the `λv=0` verdict + anchor-retention):

```python
"""#388 — family-pinned penalty homotopy closer batch (detached, checkpointed).

For each SnLm descriptor row x best-phase epoch, ramps the V∞-anchor penalty from
a calibrated W down to 0 (``search.family_pinned_shoot.family_pinned_shoot``) and
records the UNPENALIZED (λv=0) verdict: converged?, defect, emerged V∞ vs the
SOURCED anchors, anchor-match, anchor-retention, bend. HELD — no writeback. A row
is flagged PROPOSED V0->V1 only for ``mcconaghy-2006-em-k2`` AND only if the λv=0
solve converges, anchor-matches within 0.5 km/s of BOTH anchors, and is
bend-feasible — recorded, never applied.

The penalty weight ladder top ``W`` is a TUNING value (not sourced): calibrated so
the penalty bites without swamping continuity. See the design spec.

Launch detached:
    setsid nohup uv run python scripts/shooter_family_pinned_batch.py --resume \
        > data/runs/shooter-family-pinned.log 2>&1 &
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.family_pinned_shoot import family_pinned_shoot
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import (
    candidate_epochs,
    russell_parent_to_ballistic_seed,
)
from cyclerfinder.search.shooter_russell_seed import russell_shooting_seed

K: int = 3
LAUNCH_WINDOW_SYNODICS: range = range(1, 22)
EPOCH_GRID: int = 100

# Penalty weight ladder top (TUNING value, not sourced — see design spec §Open
# questions). Calibrated so the penalty residual is comparable to the converged
# continuity residual on row 9.353Gg2.
WEIGHT_LADDER: tuple[float, ...] = (40.0, 10.0, 2.5, 0.5, 0.0)
SHOOT_ACCURACY: float = 1e-9
SHOOT_MAX_NFEV: int = 100
LEG_WALL_BUDGET_SEC: float = 30.0
ANCHOR_MATCH_TOL_KMS: float = 0.5

RUNLOG = Path("data/runs/shooter-family-pinned.jsonl")


def _load_done(runlog: Path) -> set[tuple[str, int]]:
    done: set[tuple[str, int]] = set()
    if not runlog.exists():
        return done
    with runlog.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "id" in rec and "epoch_index" in rec:
                done.add((str(rec["id"]), int(rec["epoch_index"])))
    return done


def _append(runlog: Path, rec: dict[str, Any]) -> None:
    runlog.parent.mkdir(parents=True, exist_ok=True)
    with runlog.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", action="store_true", help="skip (row,epoch) pairs in runlog")
    args = parser.parse_args()
    sys.stdout.reconfigure(line_buffering=True)

    m = RussellModel()
    ephem = Ephemeris("astropy")
    done = _load_done(RUNLOG) if args.resume else set()
    if args.resume:
        print(f"[resume] {len(done)} (row,epoch) pairs already in {RUNLOG}")

    t_wall0 = time.monotonic()
    catalog = load_catalog()

    for e in catalog.entries:
        wall = time.monotonic() - t_wall0
        try:
            phsi = descriptor_to_phsi(e.raw)
            if phsi is None:
                continue
            cyc = assemble_cycler(m, phsi)
            if cyc is None:
                continue
            try:
                seed = russell_parent_to_ballistic_seed(m, cyc, e.raw)
            except ValueError:
                continue

            vlevel = str(e.raw.get("validation_level", "V0"))
            bodies = tuple(dict.fromkeys(seed.sequence))
            anchors = {"E": seed.vinf_anchor_e_kms, "M": seed.vinf_anchor_m_kms}
            epochs = candidate_epochs(
                ephem, 0.0, launch_window_synodics=LAUNCH_WINDOW_SYNODICS, grid=EPOCH_GRID
            )[:K]

            for epoch_index, t0 in enumerate(epochs):
                if (e.id, epoch_index) in done:
                    continue
                wall = time.monotonic() - t_wall0
                t_shoot0 = time.monotonic()

                def _hb(
                    kind: str,
                    count: int,
                    defect_norm: float,
                    elapsed: float,
                    _id: str = e.id,
                    _ep: int = epoch_index,
                ) -> None:
                    now = time.monotonic() - t_wall0
                    print(
                        f"[{now:.0f}s]   {_id:24s} ep{_ep} {kind}{count} "
                        f"defect={defect_norm:.3e} ({elapsed:.0f}s)"
                    )

                row_error: str | None = None
                res = None
                try:
                    sseed = russell_shooting_seed(seed, t0_sec=t0, ephem=ephem)
                    res = family_pinned_shoot(
                        sseed,
                        ephem=ephem,
                        bodies=bodies,
                        vinf_anchors=anchors,
                        weight_ladder=WEIGHT_LADDER,
                        accuracy=SHOOT_ACCURACY,
                        max_nfev=SHOOT_MAX_NFEV,
                        max_wall_sec=LEG_WALL_BUDGET_SEC,
                        progress=_hb,
                    )
                except Exception as exc:  # honest per-(row,epoch) record, never raised
                    row_error = f"{type(exc).__name__}: {exc}"

                shoot_wall = time.monotonic() - t_shoot0
                if res is None:
                    _append(
                        RUNLOG,
                        {
                            "id": e.id,
                            "validation_level": vlevel,
                            "sequence": list(seed.sequence),
                            "epoch_index": epoch_index,
                            "epoch_sec": t0,
                            "shot": False,
                            "error": row_error,
                            "shoot_wall_sec": shoot_wall,
                            "anchor_e_kms": seed.vinf_anchor_e_kms,
                            "anchor_m_kms": seed.vinf_anchor_m_kms,
                        },
                    )
                    print(f"[{wall:.0f}s] {e.id:24s} [{vlevel}] ep{epoch_index} NO SHOOT ({row_error})")
                    continue

                final = res.final
                vinf = list(final.vinf_per_encounter_kms)
                best_e = min((abs(v - seed.vinf_anchor_e_kms) for v in vinf), default=float("inf"))
                best_m = min((abs(v - seed.vinf_anchor_m_kms) for v in vinf), default=float("inf"))
                anchor_match = best_e <= ANCHOR_MATCH_TOL_KMS and best_m <= ANCHOR_MATCH_TOL_KMS
                converged = final.converged
                bend = final.bend_feasible
                promote = e.id == "mcconaghy-2006-em-k2" and converged and anchor_match and bend

                _append(
                    RUNLOG,
                    {
                        "id": e.id,
                        "validation_level": vlevel,
                        "sequence": list(seed.sequence),
                        "epoch_index": epoch_index,
                        "epoch_sec": t0,
                        "shot": True,
                        "final_weight": res.final_weight,
                        "converged": converged,
                        "defect_norm": final.defect_norm,
                        "seed_defect_norm": final.seed_defect_norm,
                        "vinf_per_encounter_kms": vinf,
                        "anchor_e_kms": seed.vinf_anchor_e_kms,
                        "anchor_m_kms": seed.vinf_anchor_m_kms,
                        "best_e_residual_kms": best_e,
                        "best_m_residual_kms": best_m,
                        "anchor_match": anchor_match,
                        "anchor_retention_kms": res.anchor_retention_kms,
                        "bend_feasible": bend,
                        "trace": res.trace,
                        "promote_proposed_held": promote,
                        "shoot_wall_sec": shoot_wall,
                        "error": row_error,
                    },
                )
                flag = " *** PROPOSED V0->V1 (HELD) ***" if promote else ""
                print(
                    f"[{wall:.0f}s] {e.id:24s} [{vlevel}] ep{epoch_index} "
                    f"conv={converged} defect={final.defect_norm:.3e} "
                    f"vinf={[round(v, 2) for v in vinf]} "
                    f"anchorE/M={best_e:.2f}/{best_m:.2f} retain={res.anchor_retention_kms:.2f} "
                    f"match={anchor_match} bend={bend} ({shoot_wall:.0f}s){flag}"
                )
        except Exception as exc:  # one bad row must not abort the batch
            _append(RUNLOG, {"id": e.id, "shot": False, "error": f"row-fatal {type(exc).__name__}: {exc}"})
            print(f"[{wall:.0f}s] {e.id:24s} ROW-FATAL {type(exc).__name__}: {exc}")

    print()
    print("=" * 72)
    print(f"Batch complete. Runlog: {RUNLOG} (K={K}, weight_ladder={WEIGHT_LADDER})")
    print("HELD — no writeback (data/catalogue.yaml and validate.py untouched).")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Verify the script parses and imports cleanly**

Run: `uv run python -c "import ast; ast.parse(open('scripts/shooter_family_pinned_batch.py').read()); print('parse OK')"`
Expected: `parse OK`.
Run: `uv run ruff check scripts/shooter_family_pinned_batch.py && uv run ruff format scripts/shooter_family_pinned_batch.py`
Expected: `All checks passed!` then formatting applied/clean.

- [ ] **Step 3: One-row time-boxed smoke (confirms the end-to-end family-pinned path returns a verdict, not an error)**

Run (time-boxed; uses a short ladder + small nfev so it returns in minutes):

```bash
timeout 600 uv run python -u -c "
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import candidate_epochs, russell_parent_to_ballistic_seed
from cyclerfinder.search.shooter_russell_seed import russell_shooting_seed
from cyclerfinder.search.family_pinned_shoot import family_pinned_shoot
m=RussellModel(); ephem=Ephemeris('astropy'); cat=load_catalog()
for e in cat.entries:
    p=descriptor_to_phsi(e.raw)
    if p is None: continue
    c=assemble_cycler(m,p)
    if c is None: continue
    try: seed=russell_parent_to_ballistic_seed(m,c,e.raw)
    except ValueError: continue
    bodies=tuple(dict.fromkeys(seed.sequence))
    t0=candidate_epochs(ephem,0.0,launch_window_synodics=range(1,3),grid=20)[0]
    ss=russell_shooting_seed(seed,t0_sec=t0,ephem=ephem)
    r=family_pinned_shoot(ss,ephem=ephem,bodies=bodies,
        vinf_anchors={'E':seed.vinf_anchor_e_kms,'M':seed.vinf_anchor_m_kms},
        weight_ladder=(10.0,0.0),accuracy=1e-9,max_nfev=3,max_wall_sec=30.0)
    print('SMOKE',e.id,'final_w',r.final_weight,'conv',r.final.converged,
          'vinf',[round(v,2) for v in r.final.vinf_per_encounter_kms],
          'retain',round(r.anchor_retention_kms,3))
    break
"
```
Expected: a `SMOKE ...` line with a finite emerged V∞ list (not an exception). If it exits 124 (timeout) note it and re-run with `weight_ladder=(10.0,0.0)` already minimal — the path is validated by reaching the first heartbeat; record that.

- [ ] **Step 4: Commit (pathspec-scoped)**

```bash
git status
git add scripts/shooter_family_pinned_batch.py
git commit -m "search/#388: family-pinned homotopy batch (detached, checkpointed, --resume, heartbeat)"
```

- [ ] **Step 5: Launch the detached run + report**

(Execution-time, after all tasks land and the smoke is green.)

```bash
mkdir -p data/runs
setsid nohup uv run python scripts/shooter_family_pinned_batch.py --resume \
    > data/runs/shooter-family-pinned.log 2>&1 &
```
Then monitor via the heartbeat lines / `data/runs/shooter-family-pinned.jsonl`, and
write `docs/notes/2026-06-22-family-pinned-results.md` (per-row table: λv=0
converged / defect / emerged V∞ vs anchor / anchor-match / **anchor-retention** /
bend; the verdict). PROPOSED V0→V1 only under the golden gate; HELD either way.

---

## Final verification
- [ ] `uv run pytest tests/nbody/test_shooter_family_pinned.py tests/search/test_family_pinned_shoot.py tests/nbody/test_shooter_stm_jacobian.py tests/nbody/test_propagator_stm.py -q` (incl `-m slow`) — PASS.
- [ ] `uv run ruff check <this plan's files>` and `uv run ruff format --check <this plan's files>` — clean (do NOT run repo-wide ruff; the other agent's untracked scratch carries unrelated errors).
- [ ] `uv run mypy src/cyclerfinder/nbody/shooter.py src/cyclerfinder/search/family_pinned_shoot.py` — clean.
- [ ] Confirm NO edit to `data/catalogue.yaml` / `src/cyclerfinder/data/validate.py` by this work (`git log --stat` for this plan's commits shows only the 5 files above).
- [ ] Confirm every commit used explicit pathspecs (no `git add -A`).
```
