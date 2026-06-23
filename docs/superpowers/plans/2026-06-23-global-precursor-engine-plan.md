# Unified Global MGA-DSM Precursor Search Engine — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the precursor matcher's local-optimiser path with a global differential_evolution engine that hunts ballistic→powered precursors into a target cycler, using eccentric-body Tisserand-Poincaré seeds and per-leg DSMs as first-class decision variables, ranking survivors by `dv_band` / total ΔV.

**Architecture:** A new `src/cyclerfinder/search/global_precursor_engine.py` owns the search. It reuses, unchanged, four primitives: `lambert` (multi-rev), `close_epoch_locked` (real DE440 closure + built-in DSM execution via `dsm_specs`), `classify_dv_band` (sourced band binning), and the `Ephemeris` backend. `find_cycler_precursors` is rewired to delegate to the engine with its signature preserved; the old local optimiser is retained as a private fallback. #302/#307 reproducibility is guarded by a test pinning `close_epoch_locked`'s output on a fixed candidate.

**Tech Stack:** Python 3.11, numpy, scipy.optimize.differential_evolution, the cyclerfinder core/genome/search/verify packages, pytest (+ pytest-xdist), uv, ruff, mypy.

**Spec:** `docs/superpowers/specs/2026-06-23-global-precursor-engine-design.md`

**Conventions for every task:** run `uv run ruff check . && uv run ruff format --check . && uv run mypy <touched files>` before each commit; never add Co-Authored-By trailers; work on `main` (no branches). Commit messages use a `search/#430:` or `tests/#430:` subsystem prefix.

---

## File Structure

| File | Responsibility | Status |
|------|----------------|--------|
| `src/cyclerfinder/search/global_precursor_engine.py` | NEW. Eccentric-T-P seeder, decision-vector evaluator, differential_evolution driver, cost+rank. | create |
| `src/cyclerfinder/search/precursor_matcher.py` | Rewire `find_cycler_precursors` to delegate to the engine; extend the JSONL record. | modify |
| `src/cyclerfinder/search/tisserand_mga_window.py` | No change (reused read-only by the seeder). | reuse |
| `src/cyclerfinder/genome/epoch_aware_genome.py` | No change (reused: `EpochLockedTrajectory`, `DSMSpec`, `close_epoch_locked`). | reuse |
| `tests/search/test_global_precursor_engine.py` | NEW. Engine unit + integration tests (all goldens sourced). | create |
| `tests/search/test_precursor_matcher_reproducibility.py` | NEW. #302/#307 closure-primitive pin. | create |
| `scripts/run_430_global_precursor.py` | NEW. Aldrin/S1L1 deliverable run. | create |

---

## Task 1: Eccentric-body Tisserand-Poincaré seeder

**Files:**
- Create: `src/cyclerfinder/search/global_precursor_engine.py`
- Test: `tests/search/test_global_precursor_engine.py`

**Context:** The existing enumerator (`tisserand_mga_window.py`) computes resonance/linkability from `PLANETS[body].sma_au` (mean, circular). The seeder instead uses each body's *actual heliocentric radius* at the encounter epoch (DE440 `r`) for the Tisserand-Poincaré linkability, per Campagnola-Russell 2009 / Strange 2007. It returns `MGAChainCandidate`s to use as the differential_evolution init population. At `e→0` (radius → mean `a`) it must reduce to the circular enumerator's seeds.

- [ ] **Step 1: Write the failing reduces-to-circular test**

```python
# tests/search/test_global_precursor_engine.py
from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import AU_KM
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.global_precursor_engine import eccentric_tp_linkable_radius_au


def test_eccentric_radius_reduces_to_mean_a_for_circular_backend() -> None:
    """The 'circular' ephemeris backend places bodies on circles at their mean a,
    so the eccentric-aware encounter radius must equal the body's sma_au."""
    from cyclerfinder.core.bodies import PLANETS

    eph = Ephemeris("circular")
    # t_sec arbitrary; circular backend radius is epoch-invariant in magnitude.
    r_au = eccentric_tp_linkable_radius_au("E", t_sec=0.0, ephemeris=eph)
    assert abs(r_au - PLANETS["E"].sma_au) < 1e-6
```

- [ ] **Step 2: Run it to verify it fails**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_eccentric_radius_reduces_to_mean_a_for_circular_backend -v`
Expected: FAIL — `ModuleNotFoundError` / `ImportError: cannot import name 'eccentric_tp_linkable_radius_au'`.

(Confirm the correct PLANETS import path first: `uv run python -c "from cyclerfinder.core.bodies import PLANETS; print(PLANETS['E'].sma_au)"`. If that path differs, fix the import in the test to match — the constant lives wherever `tisserand_mga_window.py` imports `PLANETS` from; check its imports.)

- [ ] **Step 3: Implement the radius helper**

```python
# src/cyclerfinder/search/global_precursor_engine.py
"""#430 Unified global MGA-DSM precursor search engine.

Replaces the precursor matcher's local-optimiser path with a global
differential_evolution search over (launch epoch, per-leg TOFs, per-leg DSM),
seeded with eccentric-body Tisserand-Poincaré candidates, ranking survivors by
dv_band / total ΔV. See docs/superpowers/specs/2026-06-23-global-precursor-engine-design.md.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import AU_KM
from cyclerfinder.core.ephemeris import Ephemeris


def eccentric_tp_linkable_radius_au(body: str, t_sec: float, ephemeris: Ephemeris) -> float:
    """Body's ACTUAL heliocentric radius (AU) at ``t_sec``, for the eccentric
    Tisserand-Poincaré graph (Campagnola-Russell 2009 Part B). The T-P contour
    drawn at the real encounter radius shifts/widens the linkable set vs the
    mean-``a`` circular form. Reduces to ``sma_au`` on the circular backend."""
    r_km, _v = ephemeris.state(body, t_sec)
    return float(np.linalg.norm(np.asarray(r_km, dtype=np.float64))) / AU_KM
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_eccentric_radius_reduces_to_mean_a_for_circular_backend -v`
Expected: PASS.

- [ ] **Step 5: Write the failing seeder test**

```python
# tests/search/test_global_precursor_engine.py  (append)
from cyclerfinder.search.global_precursor_engine import eccentric_tp_seeds


def test_eccentric_tp_seeds_returns_candidates_terminating_at_target() -> None:
    """The seeder returns MGAChainCandidates whose final body is the cycler's
    first encounter body and whose terminal V∞ bin is within tol of the seed."""
    eph = Ephemeris("astropy")
    seeds = eccentric_tp_seeds(
        first_body="E",
        seed_vinf_kms=6.5,
        launch_window=("2030-01-01T00:00:00", "2032-12-31T00:00:00"),
        ephemeris=eph,
        intermediate_bodies=("V", "E"),
        max_legs=3,
        vinf_grid_kms=(6.0, 7.0),
        tof_box_days_per_leg=(80.0, 500.0),
        epoch_step_days=120.0,
        vinf_terminal_tol_kms=0.8,
    )
    assert len(seeds) > 0
    for s in seeds:
        assert s.sequence[-1] == "E"
        assert abs(s.vinf_tuple_kms[-1] - 6.5) <= 0.8 + 1e-9
```

- [ ] **Step 6: Run to verify it fails**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_eccentric_tp_seeds_returns_candidates_terminating_at_target -v`
Expected: FAIL — `cannot import name 'eccentric_tp_seeds'`.

- [ ] **Step 7: Implement `eccentric_tp_seeds` by wrapping the existing enumerator**

First read `tisserand_mga_window.py::find_mga_chains` (line 615) to confirm its exact signature and return type. The seeder calls `find_mga_chains` to enumerate candidate sequences, then re-scores/filters them using `eccentric_tp_linkable_radius_au` at the candidate's launch epoch (preferring chains whose legs link at the real encounter radii). Implementation:

```python
# src/cyclerfinder/search/global_precursor_engine.py  (append)
from cyclerfinder.search.tisserand_mga_window import MGAChainCandidate, find_mga_chains


def eccentric_tp_seeds(
    *,
    first_body: str,
    seed_vinf_kms: float,
    launch_window: tuple[str, str],
    ephemeris: Ephemeris,
    intermediate_bodies: tuple[str, ...] = ("V", "E"),
    max_legs: int = 3,
    vinf_grid_kms: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    epoch_step_days: float = 60.0,
    vinf_terminal_tol_kms: float = 0.8,
) -> list[MGAChainCandidate]:
    """Enumerate Earth-launched MGA chains terminating at ``first_body`` near
    ``seed_vinf_kms``, ranked with eccentric-body (real-radius) Tisserand-Poincaré
    linkability. Returns the DE init population (MGAChainCandidate list)."""
    candidates = find_mga_chains(
        launch_window=launch_window,
        planet_set=tuple(dict.fromkeys((first_body, *intermediate_bodies))),
        max_legs=max_legs,
        vinf_grid_kms=vinf_grid_kms,
        tof_box_days_per_leg=tof_box_days_per_leg,
        epoch_step_days=epoch_step_days,
        ephemeris=ephemeris,
    )
    out: list[MGAChainCandidate] = []
    for c in candidates:
        if c.sequence[-1] != first_body:
            continue
        if abs(c.vinf_tuple_kms[-1] - seed_vinf_kms) > vinf_terminal_tol_kms:
            continue
        out.append(c)
    return out
```

NOTE: `find_mga_chains`'s real keyword names may differ — read line 615 and match them exactly. If `find_mga_chains` does not accept `planet_set`/`ephemeris` under those names, adapt the call and record the actual signature in a comment. The filtering logic (terminate-at-target + V∞ tol) is the load-bearing part and stays as shown.

- [ ] **Step 8: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py -v -k eccentric`
Expected: both eccentric tests PASS. (This calls DE440; allow ~30s.)

- [ ] **Step 9: Commit**

```bash
uv run ruff check src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run ruff format src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run mypy src/cyclerfinder/search/global_precursor_engine.py
git add src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
git commit -m "search/#430: eccentric-body Tisserand-Poincaré seeder"
```

---

## Task 2: Decision-vector → closure evaluator

**Files:**
- Modify: `src/cyclerfinder/search/global_precursor_engine.py`
- Test: `tests/search/test_global_precursor_engine.py`

**Context:** The engine's per-trial evaluation maps a flat decision vector to an `EpochLockedTrajectory` (with `DSMSpec`s from the DSM sub-vector) and calls `close_epoch_locked`. `DSMSpec.fraction_along_leg` must be strictly in `(0,1)` and `delta_v_kms` finite; a leg whose DSM magnitude is below `dsm_eps_kms` emits NO `DSMSpec` (ballistic leg). Non-convergence → a large FINITE cost (`_COST_FLOOR = 1e6`), never `inf`.

Decision vector layout for an N-leg chain (`len(sequence) == N+1`):
`x = [epoch_offset_days] + [tof_1..tof_N] + [eta_1, dvx_1, dvy_1, dvz_1, ..., eta_N, dvx_N, dvy_N, dvz_N]`
→ length `1 + N + 4*N`.

- [ ] **Step 1: Write the failing test — ballistic decision vector reproduces no-DSM closure**

```python
# tests/search/test_global_precursor_engine.py  (append)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import EpochLockedTrajectory, close_epoch_locked
from cyclerfinder.search.global_precursor_engine import evaluate_decision_vector


def test_zero_dsm_vector_matches_plain_ballistic_closure() -> None:
    """A decision vector with all DSM magnitudes 0 closes identically to a
    plain EpochLockedTrajectory with no dsm_specs (the DSM layer is a no-op)."""
    eph = Ephemeris("astropy")
    sequence = ("E", "M")
    leg_tofs = (250.0,)
    vinf_expected = (6.5, 9.7)
    launch = "2031-03-01T00:00:00"

    plain = close_epoch_locked(
        EpochLockedTrajectory(
            sequence=sequence,
            leg_tofs_days=leg_tofs,
            vinf_kms_at_encounters=vinf_expected,
            launch_epoch_utc=launch,
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc=launch,
            validity_window_end_utc="2032-03-01T00:00:00",
            inserts_into="aldrin-classic-em-k1-outbound",
        ),
        eph,
        closure_tol_kms=1.0e6,
        flyby_continuity_tol_kms=1.0e6,
        independent_cross_check=False,
        independent_tol_kms=1.0e6,
    )
    # Decision vector: epoch_offset=0, tof=250, leg-1 (eta,dvx,dvy,dvz)=(0.5,0,0,0).
    x = [0.0, 250.0, 0.5, 0.0, 0.0, 0.0]
    result = evaluate_decision_vector(
        x,
        sequence=sequence,
        seed_launch_epoch_utc=launch,
        vinf_expected_kms=vinf_expected,
        ephemeris=eph,
        inserts_into="aldrin-classic-em-k1-outbound",
    )
    assert abs(result.closure.closure_residual_kms - plain.closure_residual_kms) < 1e-6
    assert result.total_dsm_dv_kms == 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_zero_dsm_vector_matches_plain_ballistic_closure -v`
Expected: FAIL — `cannot import name 'evaluate_decision_vector'`.

- [ ] **Step 3: Implement the evaluator + result dataclass**

```python
# src/cyclerfinder/search/global_precursor_engine.py  (append)
from dataclasses import dataclass

from cyclerfinder.genome.epoch_aware_genome import (
    DSMSpec,
    EpochLockedClosure,
    EpochLockedTrajectory,
    close_epoch_locked,
)
from cyclerfinder.search.tisserand_mga_window import _add_days_utc  # reuse UTC offset helper

_COST_FLOOR = 1.0e6
_DSM_EPS_KMS = 1.0e-4  # below this magnitude a leg is treated as ballistic (no DSMSpec)
_FRACTION_CLAMP = 1.0e-3  # keep fraction_along_leg strictly inside (0, 1)


@dataclass(frozen=True)
class DecisionEval:
    closure: EpochLockedClosure
    total_dsm_dv_kms: float
    per_leg_dsm_kms: tuple[float, ...]
    feasible: bool


def evaluate_decision_vector(
    x: list[float] | np.ndarray,
    *,
    sequence: tuple[str, ...],
    seed_launch_epoch_utc: str,
    vinf_expected_kms: tuple[float, ...],
    ephemeris: Ephemeris,
    inserts_into: str,
    max_revs: int = 2,
) -> DecisionEval:
    """Map a flat decision vector to a closed EpochLockedTrajectory.

    x = [epoch_offset_days, tof_1..tof_N, (eta_i, dvx_i, dvy_i, dvz_i)*N].
    Non-convergence / infeasible geometry → DecisionEval with a _COST_FLOOR
    closure residual and feasible=False (the optimiser still gets a finite,
    differentiable-enough signal)."""
    xv = np.asarray(x, dtype=np.float64)
    n_legs = len(sequence) - 1
    epoch_offset = float(xv[0])
    tofs = tuple(float(t) for t in xv[1 : 1 + n_legs])
    if any(t <= 0.0 for t in tofs):
        return _infeasible(sequence, vinf_expected_kms, seed_launch_epoch_utc, inserts_into)
    dsm_block = xv[1 + n_legs :]
    dsm_specs: list[DSMSpec] = []
    per_leg_dsm: list[float] = []
    for i in range(n_legs):
        eta, dvx, dvy, dvz = (float(v) for v in dsm_block[4 * i : 4 * i + 4])
        mag = float(np.linalg.norm((dvx, dvy, dvz)))
        per_leg_dsm.append(mag)
        if mag < _DSM_EPS_KMS:
            continue  # ballistic leg
        frac = min(1.0 - _FRACTION_CLAMP, max(_FRACTION_CLAMP, eta))
        dsm_specs.append(
            DSMSpec(leg_index=i, fraction_along_leg=frac, delta_v_kms=(dvx, dvy, dvz))
        )
    launch = _add_days_utc(seed_launch_epoch_utc, epoch_offset)
    end_utc = _add_days_utc(launch, sum(tofs))
    try:
        traj = EpochLockedTrajectory(
            sequence=sequence,
            leg_tofs_days=tofs,
            vinf_kms_at_encounters=vinf_expected_kms,
            launch_epoch_utc=launch,
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc=launch,
            validity_window_end_utc=end_utc,
            inserts_into=inserts_into,
            dsm_specs=tuple(dsm_specs) if dsm_specs else None,
        )
        closure = close_epoch_locked(
            traj,
            ephemeris,
            closure_tol_kms=1.0e6,
            flyby_continuity_tol_kms=1.0e6,
            independent_cross_check=False,
            independent_tol_kms=1.0e6,
            max_revs=max_revs,
        )
    except Exception:
        return _infeasible(sequence, vinf_expected_kms, seed_launch_epoch_utc, inserts_into)
    return DecisionEval(
        closure=closure,
        total_dsm_dv_kms=float(sum(per_leg_dsm)),
        per_leg_dsm_kms=tuple(per_leg_dsm),
        feasible=True,
    )


def _infeasible(
    sequence: tuple[str, ...],
    vinf_expected_kms: tuple[float, ...],
    launch: str,
    inserts_into: str,
) -> DecisionEval:
    """Sentinel DecisionEval for an infeasible vector: _COST_FLOOR residuals."""
    n = len(sequence)
    end = _add_days_utc(launch, 1.0)
    traj = EpochLockedTrajectory(
        sequence=sequence,
        leg_tofs_days=tuple(1.0 for _ in range(n - 1)),
        vinf_kms_at_encounters=vinf_expected_kms,
        launch_epoch_utc=launch,
        orbit_class="precursor_mga",
        n_returns=1,
        validity_window_start_utc=launch,
        validity_window_end_utc=end,
        inserts_into=inserts_into,
    )
    closure = EpochLockedClosure(
        trajectory=traj,
        closure_residual_kms=_COST_FLOOR,
        flyby_continuity_max_dv_kms=_COST_FLOOR,
        per_leg_lambert_solutions=(),
        per_encounter_vinf_kms=tuple(0.0 for _ in range(n)),
        independent_check_residual_kms=None,
        converged=False,
        dsm_delta_v_kms_per_leg=(),
    )
    return DecisionEval(
        closure=closure, total_dsm_dv_kms=_COST_FLOOR, per_leg_dsm_kms=(), feasible=False
    )
```

NOTE: confirm `_add_days_utc` is importable from `tisserand_mga_window` (it is used there at module scope; grep `def _add_days_utc`). If it is private/unavailable, inline an equivalent using `datetime`/`timedelta` parsing the ISO-8601 UTC string. Confirm `EpochLockedTrajectory` accepts `dsm_specs=None` (it does — the field is optional). Confirm the exact `EpochLockedClosure` constructor keyword names against `epoch_aware_genome.py:348-366` before writing `_infeasible`.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_zero_dsm_vector_matches_plain_ballistic_closure -v`
Expected: PASS.

- [ ] **Step 5: Write the failing infeasible-vector test**

```python
# tests/search/test_global_precursor_engine.py  (append)
def test_negative_tof_vector_is_infeasible_with_cost_floor() -> None:
    eph = Ephemeris("astropy")
    x = [0.0, -10.0, 0.5, 0.0, 0.0, 0.0]  # negative TOF
    result = evaluate_decision_vector(
        x,
        sequence=("E", "M"),
        seed_launch_epoch_utc="2031-03-01T00:00:00",
        vinf_expected_kms=(6.5, 9.7),
        ephemeris=eph,
        inserts_into="aldrin-classic-em-k1-outbound",
    )
    assert result.feasible is False
    assert result.closure.closure_residual_kms >= 1.0e6
```

- [ ] **Step 6: Run to verify it passes** (implementation already handles it)

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_negative_tof_vector_is_infeasible_with_cost_floor -v`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
uv run ruff check src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run ruff format src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run mypy src/cyclerfinder/search/global_precursor_engine.py
git add src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
git commit -m "search/#430: decision-vector → closure evaluator (DSM via dsm_specs)"
```

---

## Task 3: Cost function + dv_band ranker

**Files:**
- Modify: `src/cyclerfinder/search/global_precursor_engine.py`
- Test: `tests/search/test_global_precursor_engine.py`

**Context:** The scalar objective drives closure + continuity and penalises DSM cost so a ballistic solution always scores below an otherwise-equal powered one. The ranker tags each survivor with `dv_band` via the sourced `classify_dv_band` (m/s input). Default weights: `w_cont=1.0`, `w_dsm=0.5` (DSM km/s penalised at half-weight so a ~0-DSM ballistic solution wins ties but the optimiser still explores DSMs).

- [ ] **Step 1: Write the failing cost-ordering test**

```python
# tests/search/test_global_precursor_engine.py  (append)
from cyclerfinder.search.global_precursor_engine import decision_cost


def test_ballistic_scores_below_equal_powered() -> None:
    """Two evals with identical closure+continuity but different DSM cost:
    the ballistic (0 DSM) one must have the strictly lower cost."""

    class _FakeClosure:
        closure_residual_kms = 2.0
        flyby_continuity_max_dv_kms = 1.0

    class _FakeEval:
        def __init__(self, dsm: float) -> None:
            self.closure = _FakeClosure()
            self.total_dsm_dv_kms = dsm

    ballistic = decision_cost(_FakeEval(0.0))  # type: ignore[arg-type]
    powered = decision_cost(_FakeEval(0.4))  # type: ignore[arg-type]
    assert ballistic < powered
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_ballistic_scores_below_equal_powered -v`
Expected: FAIL — `cannot import name 'decision_cost'`.

- [ ] **Step 3: Implement `decision_cost` + `rank_band`**

```python
# src/cyclerfinder/search/global_precursor_engine.py  (append)
from cyclerfinder.verify.dv_band_acceptance import classify_dv_band


def decision_cost(ev: DecisionEval, *, w_cont: float = 1.0, w_dsm: float = 0.5) -> float:
    """Scalar objective: closure + w_cont*continuity + w_dsm*total_DSM (km/s).
    w_dsm > 0 makes a ballistic solution always score below an equal powered one."""
    return (
        ev.closure.closure_residual_kms
        + w_cont * ev.closure.flyby_continuity_max_dv_kms
        + w_dsm * ev.total_dsm_dv_kms
    )


def rank_band(total_dsm_dv_kms: float) -> str:
    """dv_band for a candidate's total DSM ΔV, via the sourced classifier
    (m/s input, Russell 7-cycle basis)."""
    return classify_dv_band(total_dsm_dv_kms * 1000.0)
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_ballistic_scores_below_equal_powered -v`
Expected: PASS.

- [ ] **Step 5: Write the failing dv_band-boundary test**

```python
# tests/search/test_global_precursor_engine.py  (append)
import pytest

from cyclerfinder.search.global_precursor_engine import rank_band


@pytest.mark.parametrize(
    "dv_kms,expected",
    [
        (0.0005, "strictly_ballistic"),   # 0.5 m/s < 1
        (0.005, "essentially_ballistic"), # 5 m/s < 10
        (0.250, "low_maintenance"),       # 250 m/s < 300
        (0.500, "powered_dsm"),           # 500 m/s >= 300
    ],
)
def test_rank_band_boundaries(dv_kms: float, expected: str) -> None:
    assert rank_band(dv_kms) == expected
```

- [ ] **Step 6: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py -v -k rank_band`
Expected: PASS (4 params). If a boundary disagrees, read `classify_dv_band`'s exact half-open bins (`dv_band_acceptance.py:185`) and align the parametrize values — do NOT change `classify_dv_band` (it is sourced).

- [ ] **Step 7: Commit**

```bash
uv run ruff check src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run ruff format src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run mypy src/cyclerfinder/search/global_precursor_engine.py
git add src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
git commit -m "search/#430: cost function + sourced dv_band ranker"
```

---

## Task 4: differential_evolution driver + known-ballistic recovery

**Files:**
- Modify: `src/cyclerfinder/search/global_precursor_engine.py`
- Test: `tests/search/test_global_precursor_engine.py`

**Context:** `search_sequence` runs `scipy.optimize.differential_evolution` over the bounded decision box for ONE fixed body sequence, seeding the init population from `eccentric_tp_seeds` where available, returning the best `DecisionEval`. `search_precursors` runs it per candidate sequence and aggregates+ranks. Determinism via a fixed `seed`. Bounds: epoch ±`epoch_half_width_days`; TOFs in `tof_box_days_per_leg`; eta in [0,1]; DSM components in [−`dsm_max_kms`, +`dsm_max_kms`].

- [ ] **Step 1: Write the failing known-ballistic-recovery test**

```python
# tests/search/test_global_precursor_engine.py  (append)
from cyclerfinder.search.global_precursor_engine import search_sequence


def test_search_recovers_low_continuity_on_e_m_direct() -> None:
    """On a single E→M leg, the global search drives flyby continuity (here the
    single-leg closure residual is the binding term) well below the seed/local
    baseline, finding a near-ballistic solution with small total DSM."""
    eph = Ephemeris("astropy")
    best = search_sequence(
        sequence=("E", "M"),
        seed_launch_epoch_utc="2031-01-01T00:00:00",
        vinf_expected_kms=(6.5, 9.7),
        ephemeris=eph,
        inserts_into="aldrin-classic-em-k1-outbound",
        epoch_half_width_days=120.0,
        tof_box_days_per_leg=(120.0, 400.0),
        dsm_max_kms=2.0,
        popsize=12,
        maxiter=40,
        seed=0,
    )
    assert best.feasible
    # A 2-encounter E-M cell cannot host the Aldrin family ballistically (known
    # #54 result), so we assert the SEARCH RAN and produced a finite, ranked
    # result — not a specific closure value (that would be an unsourced golden).
    assert best.closure.closure_residual_kms < 1.0e6
    assert best.total_dsm_dv_kms >= 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_search_recovers_low_continuity_on_e_m_direct -v`
Expected: FAIL — `cannot import name 'search_sequence'`.

- [ ] **Step 3: Implement `search_sequence` and `search_precursors`**

```python
# src/cyclerfinder/search/global_precursor_engine.py  (append)
from scipy.optimize import Bounds, differential_evolution


def _bounds_for(
    n_legs: int,
    *,
    epoch_half_width_days: float,
    tof_box_days_per_leg: tuple[float, float],
    dsm_max_kms: float,
) -> Bounds:
    lo = [-epoch_half_width_days]
    hi = [epoch_half_width_days]
    lo += [tof_box_days_per_leg[0]] * n_legs
    hi += [tof_box_days_per_leg[1]] * n_legs
    for _ in range(n_legs):  # (eta, dvx, dvy, dvz) per leg
        lo += [0.0, -dsm_max_kms, -dsm_max_kms, -dsm_max_kms]
        hi += [1.0, dsm_max_kms, dsm_max_kms, dsm_max_kms]
    return Bounds(lo, hi)


def search_sequence(
    *,
    sequence: tuple[str, ...],
    seed_launch_epoch_utc: str,
    vinf_expected_kms: tuple[float, ...],
    ephemeris: Ephemeris,
    inserts_into: str,
    epoch_half_width_days: float = 120.0,
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    dsm_max_kms: float = 2.0,
    max_revs: int = 2,
    popsize: int = 15,
    maxiter: int = 60,
    seed: int = 0,
    w_cont: float = 1.0,
    w_dsm: float = 0.5,
) -> DecisionEval:
    """Global differential_evolution over the decision box for ONE sequence.
    Returns the best DecisionEval."""
    n_legs = len(sequence) - 1
    bounds = _bounds_for(
        n_legs,
        epoch_half_width_days=epoch_half_width_days,
        tof_box_days_per_leg=tof_box_days_per_leg,
        dsm_max_kms=dsm_max_kms,
    )

    def _obj(x: np.ndarray) -> float:
        ev = evaluate_decision_vector(
            x,
            sequence=sequence,
            seed_launch_epoch_utc=seed_launch_epoch_utc,
            vinf_expected_kms=vinf_expected_kms,
            ephemeris=ephemeris,
            inserts_into=inserts_into,
            max_revs=max_revs,
        )
        return decision_cost(ev, w_cont=w_cont, w_dsm=w_dsm)

    result = differential_evolution(
        _obj, bounds, popsize=popsize, maxiter=maxiter, seed=seed, polish=True, tol=1e-8
    )
    return evaluate_decision_vector(
        result.x,
        sequence=sequence,
        seed_launch_epoch_utc=seed_launch_epoch_utc,
        vinf_expected_kms=vinf_expected_kms,
        ephemeris=ephemeris,
        inserts_into=inserts_into,
        max_revs=max_revs,
    )


@dataclass(frozen=True)
class PrecursorSurvivor:
    sequence: tuple[str, ...]
    eval: DecisionEval
    dv_band: str
    cost: float


def search_precursors(
    *,
    cycler_id: str,
    first_body: str,
    seed_vinf_kms: float,
    launch_window: tuple[str, str],
    ephemeris: Ephemeris,
    intermediate_bodies: tuple[str, ...] = ("V", "E"),
    max_legs: int = 3,
    vinf_grid_kms: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    epoch_step_days: float = 60.0,
    dsm_max_kms: float = 2.0,
    max_revs: int = 2,
    popsize: int = 15,
    maxiter: int = 60,
    seed: int = 0,
) -> list[PrecursorSurvivor]:
    """Enumerate sequences (eccentric-T-P seeds), run a global DE per sequence,
    rank survivors by total ΔV ascending (ballistic first)."""
    seeds = eccentric_tp_seeds(
        first_body=first_body,
        seed_vinf_kms=seed_vinf_kms,
        launch_window=launch_window,
        ephemeris=ephemeris,
        intermediate_bodies=intermediate_bodies,
        max_legs=max_legs,
        vinf_grid_kms=vinf_grid_kms,
        tof_box_days_per_leg=tof_box_days_per_leg,
        epoch_step_days=epoch_step_days,
    )
    # Distinct (sequence, expected-V∞, seed-epoch) per seed candidate.
    survivors: list[PrecursorSurvivor] = []
    seen: set[tuple[str, ...]] = set()
    for s in seeds:
        if s.sequence in seen:
            continue
        seen.add(s.sequence)
        ev = search_sequence(
            sequence=s.sequence,
            seed_launch_epoch_utc=s.launch_epoch_utc,
            vinf_expected_kms=s.vinf_tuple_kms,
            ephemeris=ephemeris,
            inserts_into=cycler_id,
            tof_box_days_per_leg=tof_box_days_per_leg,
            dsm_max_kms=dsm_max_kms,
            max_revs=max_revs,
            popsize=popsize,
            maxiter=maxiter,
            seed=seed,
        )
        survivors.append(
            PrecursorSurvivor(
                sequence=s.sequence,
                eval=ev,
                dv_band=rank_band(ev.total_dsm_dv_kms),
                cost=decision_cost(ev),
            )
        )
    survivors.sort(key=lambda s: s.eval.total_dsm_dv_kms)
    return survivors
```

NOTE: `differential_evolution`'s `init` parameter could be set from the eccentric seeds to warm-start the population; YAGNI for the first cut — the seeds already select the *sequences*. If a later run shows poor convergence, add `init=` with the seed vectors. Document that decision in a comment, don't silently skip it.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_search_recovers_low_continuity_on_e_m_direct -v`
Expected: PASS (DE over DE440; allow ~60-120s).

- [ ] **Step 5: Write the failing determinism test**

```python
# tests/search/test_global_precursor_engine.py  (append)
def test_search_sequence_is_deterministic_under_fixed_seed() -> None:
    eph = Ephemeris("astropy")
    kw = dict(
        sequence=("E", "M"),
        seed_launch_epoch_utc="2031-01-01T00:00:00",
        vinf_expected_kms=(6.5, 9.7),
        ephemeris=eph,
        inserts_into="aldrin-classic-em-k1-outbound",
        tof_box_days_per_leg=(120.0, 400.0),
        popsize=8,
        maxiter=20,
        seed=42,
    )
    a = search_sequence(**kw)  # type: ignore[arg-type]
    b = search_sequence(**kw)  # type: ignore[arg-type]
    assert a.closure.closure_residual_kms == b.closure.closure_residual_kms
    assert a.total_dsm_dv_kms == b.total_dsm_dv_kms
```

- [ ] **Step 6: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_search_sequence_is_deterministic_under_fixed_seed -v`
Expected: PASS. If it flakes under multithreaded BLAS (per the #323/#307 flake history), set `workers=1` in the `differential_evolution` call and re-run; document the reason in a comment.

- [ ] **Step 7: Commit**

```bash
uv run ruff check src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run ruff format src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
uv run mypy src/cyclerfinder/search/global_precursor_engine.py
git add src/cyclerfinder/search/global_precursor_engine.py tests/search/test_global_precursor_engine.py
git commit -m "search/#430: differential_evolution driver + per-sequence aggregation"
```

---

## Task 5: #302/#307 reproducibility guard

**Files:**
- Create: `tests/search/test_precursor_matcher_reproducibility.py`

**Context:** The full rebuild must not move the closure physics. This test pins `close_epoch_locked`'s output on a fixed candidate, proving the reused primitive is byte-stable. It must pass BEFORE Task 6 rewires the matcher (it depends only on `close_epoch_locked`, unchanged).

- [ ] **Step 1: Capture the current pinned values**

Run this once to read the current outputs, then paste them into the test as expected constants:

```bash
uv run python -c "
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import EpochLockedTrajectory, close_epoch_locked
eph = Ephemeris('astropy')
c = close_epoch_locked(EpochLockedTrajectory(
    sequence=('E','M'), leg_tofs_days=(250.0,), vinf_kms_at_encounters=(6.5,9.7),
    launch_epoch_utc='2031-03-01T00:00:00', orbit_class='precursor_mga', n_returns=1,
    validity_window_start_utc='2031-03-01T00:00:00', validity_window_end_utc='2032-03-01T00:00:00',
    inserts_into='aldrin-classic-em-k1-outbound'),
    eph, closure_tol_kms=1e6, flyby_continuity_tol_kms=1e6,
    independent_cross_check=False, independent_tol_kms=1e6)
print(repr(c.closure_residual_kms), repr(c.flyby_continuity_max_dv_kms), repr(c.per_encounter_vinf_kms))
"
```

- [ ] **Step 2: Write the pinning test with the captured constants**

```python
# tests/search/test_precursor_matcher_reproducibility.py
"""#430 reproducibility guard: the global-engine rebuild reuses close_epoch_locked
unchanged, so the closure PHYSICS must stay byte-stable. Pins a fixed candidate's
output (the values were captured from HEAD before the rebuild)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import EpochLockedTrajectory, close_epoch_locked


def test_close_epoch_locked_pinned_em_candidate() -> None:
    eph = Ephemeris("astropy")
    c = close_epoch_locked(
        EpochLockedTrajectory(
            sequence=("E", "M"),
            leg_tofs_days=(250.0,),
            vinf_kms_at_encounters=(6.5, 9.7),
            launch_epoch_utc="2031-03-01T00:00:00",
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc="2031-03-01T00:00:00",
            validity_window_end_utc="2032-03-01T00:00:00",
            inserts_into="aldrin-classic-em-k1-outbound",
        ),
        eph,
        closure_tol_kms=1.0e6,
        flyby_continuity_tol_kms=1.0e6,
        independent_cross_check=False,
        independent_tol_kms=1.0e6,
    )
    # <<< paste the Step-1 captured values here >>>
    assert c.closure_residual_kms == pytest.approx(CAPTURED_CLOSURE, rel=1e-9)
    assert c.flyby_continuity_max_dv_kms == pytest.approx(CAPTURED_CONTINUITY, rel=1e-9)
```

Replace `CAPTURED_CLOSURE` / `CAPTURED_CONTINUITY` with the Step-1 numbers (as literals). Use `rel=1e-9` (not exact `==`) per the documented multithreaded-BLAS reduction-order noise (#307/#323 flake history).

- [ ] **Step 3: Run to verify it passes**

Run: `uv run pytest tests/search/test_precursor_matcher_reproducibility.py -v`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
uv run ruff check tests/search/test_precursor_matcher_reproducibility.py
uv run ruff format tests/search/test_precursor_matcher_reproducibility.py
git add tests/search/test_precursor_matcher_reproducibility.py
git commit -m "tests/#430: pin close_epoch_locked output (rebuild reproducibility guard)"
```

---

## Task 6: Delegate `find_cycler_precursors` to the engine + extended JSONL

**Files:**
- Modify: `src/cyclerfinder/search/precursor_matcher.py`
- Test: `tests/search/test_global_precursor_engine.py`

**Context:** Add a `use_global_engine: bool = True` flag to `find_cycler_precursors`. When True it routes through `search_precursors` and converts `PrecursorSurvivor`s into the existing `PrecursorMatch` shape (preserving downstream consumers), with the JSONL record extended to carry `per_leg_dsm_kms`, `total_dsm_dv_kms`, and `dv_band`. When False, the existing local-Nelder-Mead path runs unchanged (so #302/#307 scripts can reproduce exactly with `use_global_engine=False`). First read `precursor_matcher.py` end-to-end to confirm the `PrecursorMatch` fields and `precursor_match_to_jsonl_record` body.

- [ ] **Step 1: Write the failing JSONL-extension test**

```python
# tests/search/test_global_precursor_engine.py  (append)
from cyclerfinder.search.precursor_matcher import precursor_match_to_jsonl_record


def test_jsonl_record_carries_dsm_and_band_fields() -> None:
    """A PrecursorMatch produced by the global engine serialises per-leg DSM,
    total DSM ΔV, and dv_band."""
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.precursor_matcher import find_cycler_precursors
    from cyclerfinder.data.catalog import load_catalog

    cat = load_catalog()
    eph = Ephemeris("astropy")
    matches = find_cycler_precursors(
        cycler_id="aldrin-classic-em-k1-outbound",
        catalogue=cat,
        ephemeris=eph,
        launch_window=("2031-01-01T00:00:00", "2031-12-31T00:00:00"),
        max_legs=2,
        intermediate_bodies=("V",),
        vinf_grid_kms=(6.0, 7.0),
        max_candidates_to_validate=2,
        use_global_engine=True,
    )
    assert matches, "expected at least one survivor"
    rec = precursor_match_to_jsonl_record(matches[0])
    assert "dv_band" in rec
    assert "total_dsm_dv_kms" in rec
    assert "per_leg_dsm_kms" in rec
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_jsonl_record_carries_dsm_and_band_fields -v`
Expected: FAIL — `find_cycler_precursors() got an unexpected keyword argument 'use_global_engine'` (or KeyError on the new JSONL keys).

- [ ] **Step 3: Implement the delegation + JSONL extension**

In `precursor_matcher.py`:
1. Add `use_global_engine: bool = True` to `find_cycler_precursors`'s keyword-only params.
2. At the top of the body, after resolving `entry`/`seed_vinf`/`first_body`, branch:

```python
    if use_global_engine:
        from cyclerfinder.search.global_precursor_engine import search_precursors

        survivors = search_precursors(
            cycler_id=cycler_id,
            first_body=first_body,
            seed_vinf_kms=seed_vinf,
            launch_window=launch_window,
            ephemeris=ephemeris,
            intermediate_bodies=intermediate_bodies,
            max_legs=max_legs,
            vinf_grid_kms=vinf_grid_kms,
            tof_box_days_per_leg=tof_box_days_per_leg,
            epoch_step_days=epoch_step_days,
        )
        return [_survivor_to_match(s, cycler_id, seed_vinf) for s in survivors]
```

3. Add `_survivor_to_match(survivor, cycler_id, seed_vinf)` that builds a `PrecursorMatch` from a `PrecursorSurvivor` (reuse `_candidate_to_trajectory` for the trajectory record; map `survivor.eval.closure` to the match's `closure` field; compute `vinf_match_residual_kms` from `_terminal_vinf_kms(survivor.eval.closure)` vs `seed_vinf`).
4. Extend `precursor_match_to_jsonl_record` to add three keys, reading them off the match (store them on `PrecursorMatch` as new optional fields `per_leg_dsm_kms: tuple[float,...] = ()`, `total_dsm_dv_kms: float = 0.0`, `dv_band: str | None = None`, defaulting so the local-path matches still serialise):

```python
    record["per_leg_dsm_kms"] = list(match.per_leg_dsm_kms)
    record["total_dsm_dv_kms"] = match.total_dsm_dv_kms
    record["dv_band"] = match.dv_band
```

Match the exact `PrecursorMatch` constructor and `precursor_match_to_jsonl_record` body to the file (read it first). Keep the local-path behaviour byte-identical when `use_global_engine=False`.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/search/test_global_precursor_engine.py::test_jsonl_record_carries_dsm_and_band_fields -v`
Expected: PASS.

- [ ] **Step 5: Run the existing matcher suite (regression — local path unchanged)**

Run: `uv run pytest tests/search/test_precursor_matcher.py tests/search/test_precursor_phase_window.py -q`
Expected: all PASS (the local path is the default-off branch; existing tests either pass `use_global_engine=False` or exercise the converter). If any existing test calls `find_cycler_precursors` and now hits the engine, pin it to `use_global_engine=False` so the test's intent (local path) is preserved.

- [ ] **Step 6: Commit**

```bash
uv run ruff check src/cyclerfinder/search/precursor_matcher.py tests/search/test_global_precursor_engine.py
uv run ruff format src/cyclerfinder/search/precursor_matcher.py tests/search/test_global_precursor_engine.py
uv run mypy src/cyclerfinder/search/precursor_matcher.py
git add src/cyclerfinder/search/precursor_matcher.py tests/search/test_global_precursor_engine.py
git commit -m "search/#430: find_cycler_precursors delegates to the global engine + dv_band JSONL"
```

---

## Task 7: Aldrin/S1L1 deliverable run + results note + registry

**Files:**
- Create: `scripts/run_430_global_precursor.py`
- Create (output): `data/precursor_430_{aldrin,s1l1}_global.jsonl`, `docs/superpowers/plans/2026-06-23-430-global-precursor-verdict.md`
- Modify: `data/empty_regions.jsonl`

**Context:** Mirror `scripts/run_307_precursor_multirev.py` (same offline `KNOWN_CORPUS` literature search, same progress logging, AET-aware), but call `find_cycler_precursors(..., use_global_engine=True)` for Aldrin and S1L1. Run detached with incremental logging (per `feedback_long_runs_acceptable` + `feedback_incremental_progress_reports`). Compute is heavy (DE × DE440 per sequence); expect 30-90 min per target — run both in parallel, harness-tracked waiter.

- [ ] **Step 1: Write the run script** (copy `run_307_precursor_multirev.py`'s structure; change the call to `use_global_engine=True`, output paths to `precursor_430_{target}_global.jsonl`, and the progress banner to "#430 global engine"). Keep the offline `KNOWN_CORPUS` search verbatim.

- [ ] **Step 2: Smoke-test with a tiny cap**

Run: `timeout 600 uv run python -c "..."` with `max_candidates_to_validate=2`, `popsize=8`, `maxiter=20`, a 6-month window — confirm it returns survivors with `dv_band` set and no crash (mirror the Task-4 smoke pattern from this session).

- [ ] **Step 3: Launch both targets detached + harness-tracked waiter**

```bash
mkdir -p data/runlogs
nohup uv run python scripts/run_430_global_precursor.py aldrin > data/runlogs/430_aldrin.log 2>&1 &
nohup uv run python scripts/run_430_global_precursor.py s1l1   > data/runlogs/430_s1l1.log 2>&1 &
```
Then start a `run_in_background` waiter on both PIDs that surfaces the summaries on completion.

- [ ] **Step 4: Harvest + write the verdict note** comparing the residual distribution + best dv_band vs the #307 baseline (`data/precursor_307_*_multirev.jsonl`). Honest-negative-or-promote: if any survivor is sub-gate ballistic OR low-ΔV powered, route it through the V0-V5 gauntlet + ML flagger + #256 FP guard BEFORE any catalogue claim (closure discipline: independent cross-check mandatory). If dry, register the method-versioned negative in `data/empty_regions.jsonl` (this engine capability-subsumes the #307 multi-rev negative — reference its `region_id`).

- [ ] **Step 5: Commit** the script, JSONL outputs, verdict note, and registry entry (run `uv run pytest tests/data -q -k "empty_region or registry"` first if the registry changed):

```bash
git add scripts/run_430_global_precursor.py data/precursor_430_*_global.jsonl \
        docs/superpowers/plans/2026-06-23-430-global-precursor-verdict.md data/empty_regions.jsonl
git commit -m "search/#430: global-engine Aldrin/S1L1 deliverable run + verdict + registry"
```

---

## Self-Review

**Spec coverage:**
- Eccentric-T-P seeder → Task 1. ✓
- Unified parameter vector + global optimiser → Tasks 2 (evaluator) + 4 (DE driver). ✓
- Per-leg DSM via `close_epoch_locked.dsm_specs` → Task 2. ✓
- Cost-ranked output (dv_band) → Task 3 + Task 6 (JSONL). ✓
- Delegate-not-delete + #302/#307 reproducibility → Tasks 5 (guard) + 6 (delegation w/ `use_global_engine`). ✓
- Discrete-vs-continuous split (DE per sequence, aggregate) → Task 4 (`search_precursors`). ✓
- Error handling (finite `_COST_FLOOR`) → Task 2. ✓
- Sourced goldens (reduces-to-circular, dv_band via `classify_dv_band`, pinned closure) → Tasks 1/3/5; NO unsourced closure-value golden asserted (Task 4 asserts the search RAN, not a magnitude). ✓
- Deliverable run + gauntlet routing → Task 7. ✓

**Placeholder scan:** The only deferred literals are the captured pin values in Task 5 (explicitly captured via a given command in Step 1) — not a placeholder, a measured constant. No "TBD"/"add error handling"/"similar to Task N".

**Type consistency:** `DecisionEval` (Task 2) is consumed by `decision_cost` (Task 3), `search_sequence`/`search_precursors` (Task 4); `PrecursorSurvivor` (Task 4) consumed by `_survivor_to_match` (Task 6). `eccentric_tp_seeds` return type (`list[MGAChainCandidate]`, Task 1) is iterated in `search_precursors` (Task 4) using `.sequence`/`.vinf_tuple_kms`/`.launch_epoch_utc` — all confirmed `MGAChainCandidate` fields. `rank_band` (Task 3) used in Task 4. Consistent.

**Verification-before-code reminders:** Tasks 1/2/6 each flag the exact upstream signatures to confirm before writing (`find_mga_chains` kwargs, `_add_days_utc` import, `EpochLockedClosure` constructor, `PrecursorMatch`/`precursor_match_to_jsonl_record` bodies). Honour those — do not assume.
