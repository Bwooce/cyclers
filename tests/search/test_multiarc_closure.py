"""Multi-arc closure harness (#248) — wrapper, seed enumerator, smoke.

Tests the three pieces of ``cyclerfinder.search.multiarc_closure``:

1. :func:`safe_chain_residual` returns a finite PENALTY (not a crash) when a
   decision vector drives the ephemeris epoch search outside the kernel's coverage
   (the unguarded ``OutOfRangeError`` that aborted the inline probe).
2. :func:`resonant_return_tof_grid` / :func:`resonant_return_seeds` enumerate the
   DISCRETE small-N resonant returns (descriptor value + integer multiples), not a
   +-% scale of one value.
3. A fast smoke: :func:`close_multiarc_row` runs end-to-end on the closest E-E-M-M
   row with ``n_starts=1`` and returns a FINITE residual without crashing (it
   asserts the harness RUNS, NOT that it converges).
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search import dsm_leg, multiarc_closure

DAY_S = 86400.0


def _row(row_id: str) -> dict[str, Any]:
    return load_catalog().by_id[row_id].raw


def _body_period_days(body: str) -> float:
    a_km = PLANETS[body].sma_au * AU_KM
    return float(2.0 * np.pi * np.sqrt(a_km**3 / MU_SUN_KM3_S2) / DAY_S)


# ---------------------------------------------------------------------------
# 1. Epoch-range-safe wrapper: penalty, not crash.
# ---------------------------------------------------------------------------


def test_safe_residual_penalises_out_of_range() -> None:
    """An epoch outside DE440 (1549-2650) returns a penalty, never an exception.

    A departure epoch in ~year 3500 is well past DE440's upper bound, so the very
    first ``ephem.state`` raises ``jplephem`` ``OutOfRangeError``. The wrapper must
    catch it and return ``(PENALTY_RESIDUAL_KMS, False, None)`` — the optimiser is
    then free to explore wild epochs without aborting the campaign.
    """
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy")
    seq = ("E", "M")
    # t0 ~ year 3500 (TDB seconds since J2000): far outside DE440 coverage.
    t0_out = (3500.0 - 2000.0) * 365.25 * DAY_S
    x = dsm_leg.dsm_chain_decision_vector(
        t0_sec=t0_out,
        vinf_out0_kms=3.0,
        alpha0=0.0,
        beta0=0.0,
        tof_days_per_leg=(200.0,),
        eta_per_leg=(0.5,),
    )
    res_kms, converged, result = multiarc_closure.safe_chain_residual(
        x,
        sequence=seq,
        ephem=eph,
        charge_flyby_continuity=False,
    )
    assert np.isfinite(res_kms)
    assert res_kms == multiarc_closure.PENALTY_RESIDUAL_KMS
    assert converged is False
    assert result is None


def test_safe_residual_in_range_returns_canonical_metric() -> None:
    """An in-range seed returns the corrector's own canonical metric + result.

    Bounds keep the bounded least-squares epoch search inside DE440 (an unbounded
    solver wanders out of range and is correctly penalised — that is the other
    test); with a finite epoch box the corrector runs and the wrapper passes its
    canonical fields straight through.
    """
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy")
    seq = ("E", "M")
    t0 = 27.0 * 365.25 * DAY_S  # ~2027, well inside DE440
    x = dsm_leg.dsm_chain_decision_vector(
        t0_sec=t0,
        vinf_out0_kms=3.0,
        alpha0=0.0,
        beta0=0.0,
        tof_days_per_leg=(200.0,),
        eta_per_leg=(0.5,),
    )
    bounds = dsm_leg.sequence_keyed_bounds(
        sequence=seq,
        t0_window_sec=(t0 - 300.0 * DAY_S, t0 + 300.0 * DAY_S),
    )
    res_kms, converged, result = multiarc_closure.safe_chain_residual(
        x,
        sequence=seq,
        ephem=eph,
        bounds=bounds,
        charge_flyby_continuity=False,
    )
    assert result is not None
    # The triple is the CANONICAL metric: equals the corrector's own fields.
    assert res_kms == pytest.approx(float(result.max_residual_kms))
    assert converged == bool(result.converged)
    assert np.isfinite(res_kms)


# ---------------------------------------------------------------------------
# 2. Discrete resonant-return enumeration (descriptor value + integer multiples).
# ---------------------------------------------------------------------------


def test_resonant_grid_yields_descriptor_and_integer_multiples() -> None:
    """The grid is the body period x {1..N} PLUS the descriptor value, de-duped."""
    body = "M"
    descriptor = 537.3  # an arbitrary descriptor arc duration (days)
    n = 4
    grid = multiarc_closure.resonant_return_tof_grid(body, descriptor, n_resonant=n)

    period = _body_period_days(body)
    # Every integer multiple 1..N of the body period appears.
    for k in range(1, n + 1):
        assert any(abs(g - k * period) < 1e-3 for g in grid), (
            f"missing {k}x Mars period ({k * period:.1f} d) in {grid}"
        )
    # The descriptor value itself appears.
    assert any(abs(g - descriptor) < 1e-3 for g in grid)
    # Sorted, de-duplicated.
    assert list(grid) == sorted(set(grid))


def test_resonant_grid_dedupes_descriptor_equal_to_a_multiple() -> None:
    """When the descriptor equals a period multiple it is not double-counted."""
    body = "E"
    period = _body_period_days(body)
    grid = multiarc_closure.resonant_return_tof_grid(body, 2.0 * period, n_resonant=3)
    # 3 multiples; the descriptor (== 2x period) collapses onto the 2x entry.
    assert len(grid) == 3


def test_resonant_return_seeds_enumerate_multistart_grid() -> None:
    """The seed enumerator yields >1 distinct decision vector for an E-E-M-M row.

    Each resonant leg is crossed over its discrete resonant-return grid, so the
    Cartesian product is a genuine multi-start set (design point 4) — not the single
    charged seed that was the #244 root cause.
    """
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy")
    warnings.filterwarnings("ignore")
    row = _row("mcconaghy-2006-em-k2")
    t_center = 27.0 * 365.25 * DAY_S
    seeds = list(multiarc_closure.resonant_return_seeds(row, eph, t_center, n_resonant=3))
    assert len(seeds) >= 2, "expected a multi-start grid, got a single seed"
    # Distinct decision vectors (the resonant ToF actually varies across seeds — the
    # identical-across-seeds bug the design calls out).
    xs = [tuple(np.round(s.x0, 6)) for s in seeds]
    assert len(set(xs)) == len(xs)
    # Every seed has the same sequence and a transit ToF in a sane window.
    seq = seeds[0].sequence
    assert all(s.sequence == seq for s in seeds)
    assert all(50.0 < s.transit_tof_days < 1000.0 for s in seeds)


def test_resonant_return_seeds_empty_for_no_descriptor_row() -> None:
    """A row without a per-arc descriptor yields a clean empty iterator."""
    catalog = load_catalog()
    ocampo = next(e for e in catalog.entries if e.id.startswith("russell-ocampo"))
    eph = Ephemeris("circular")
    seeds = list(multiarc_closure.resonant_return_seeds(ocampo.raw, eph, 0.0))
    assert seeds == []


# ---------------------------------------------------------------------------
# 3. Fast smoke: the driver runs end-to-end (RUNS, not converges).
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_close_multiarc_row_smoke_runs() -> None:
    """``close_multiarc_row`` runs end-to-end on the closest row with n_starts=1.

    Asserts the harness RETURNS a finite residual without crashing — NOT that it
    converges (convergence is the coordinator's bounded campaign, deferred). This is
    the proof-of-life smoke for the #248 deliverable.
    """
    pytest.importorskip("astropy")
    warnings.filterwarnings("ignore")
    eph = Ephemeris("astropy")
    row = _row("mcconaghy-2006-em-k2")
    report = multiarc_closure.close_multiarc_row(row, eph, n_starts=1, gradient="lambert")
    assert report.row_id == "mcconaghy-2006-em-k2"
    assert report.n_starts_run == 1
    assert report.n_seeds_available >= 1
    assert np.isfinite(report.best_max_residual_kms)
    assert len(report.sequence) >= 4  # E-E-M-M is a 4-body, 3-leg multi-arc row
    # converged is whatever it is; we only assert it RAN and is a clean bool.
    assert isinstance(report.converged, bool)
