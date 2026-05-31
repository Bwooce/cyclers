"""Tests for :mod:`cyclerfinder.core.lambert`.

Includes the M1 gate: ``lambert_crosscheck(...)["max_diff_mps"] < 1e-3`` on
three distinct legs (Aldrin medium, short Earth-to-Earth, long Earth-to-Mars).

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §4.4, §4.1.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.lambert import (
    LambertGeometryError,
    LambertSolution,
    lambert,
    lambert_crosscheck,
)

from .conftest import Leg

# ---------------------------------------------------------------------------
# Standalone behaviour
# ---------------------------------------------------------------------------


def test_lambert_returns_list_singleton(leg_aldrin: Leg) -> None:
    """Single-rev request returns a length-1 list with the expected metadata."""
    sols = lambert(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof)
    assert isinstance(sols, list)
    assert len(sols) == 1
    sol = sols[0]
    assert sol.n_revs == 0
    assert sol.branch == "single"
    assert sol.v1.shape == (3,)
    assert sol.v2.shape == (3,)
    assert sol.v1.dtype == np.float64
    assert sol.v2.dtype == np.float64


def test_lambert_max_revs_stub(leg_aldrin: Leg) -> None:
    """``max_revs >= 1`` is a documented M4 stub: still returns just the single-rev."""
    sols = lambert(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof, max_revs=2)
    assert len(sols) == 1
    assert sols[0].n_revs == 0


def test_lambert_retrograde(leg_aldrin: Leg) -> None:
    """Retrograde request matches lamberthub with the same flag set."""
    from lamberthub import izzo2015  # type: ignore[import-untyped]

    sols = lambert(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof, prograde=False)
    sol = sols[0]
    v1_ref, v2_ref = izzo2015(
        MU_SUN_KM3_S2,
        np.asarray(leg_aldrin.r1, dtype=np.float64),
        np.asarray(leg_aldrin.r2, dtype=np.float64),
        leg_aldrin.tof,
        M=0,
        prograde=False,
    )
    diff = max(
        float(np.linalg.norm(sol.v1 - v1_ref)),
        float(np.linalg.norm(sol.v2 - v2_ref)),
    )
    assert diff * 1000.0 < 1.0e-3


def test_lambert_zero_tof_raises(leg_aldrin: Leg) -> None:
    """Non-positive ``tof`` is a :class:`ValueError`."""
    with pytest.raises(ValueError):
        lambert(leg_aldrin.r1, leg_aldrin.r2, 0.0)
    with pytest.raises(ValueError):
        lambert(leg_aldrin.r1, leg_aldrin.r2, -100.0)


def test_lambert_180_deg_raises() -> None:
    """A pure 180-degree transfer raises :class:`LambertGeometryError`."""
    r1 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    r2 = np.array([-1.5 * AU_KM, 0.0, 0.0], dtype=np.float64)
    with pytest.raises(LambertGeometryError):
        lambert(r1, r2, 200.0 * SECONDS_PER_DAY)


def test_lambert_zero_magnitude_raises() -> None:
    """Zero-magnitude endpoints are a :class:`ValueError`."""
    r1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    r2 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    with pytest.raises(ValueError):
        lambert(r1, r2, 100.0 * SECONDS_PER_DAY)


def test_lambert_solution_dataclass_frozen() -> None:
    """:class:`LambertSolution` is frozen — direct assignment is rejected."""
    sol = LambertSolution(
        n_revs=0,
        branch="single",
        v1=np.zeros(3, dtype=np.float64),
        v2=np.zeros(3, dtype=np.float64),
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        sol.n_revs = 1  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Gate: lamberthub cross-check on three legs
# ---------------------------------------------------------------------------


def test_aldrin_leg_cross_check(leg_aldrin: Leg) -> None:
    """Aldrin E->M ~146 d: agreement with izzo+gooding < 1e-3 m/s."""
    res = lambert_crosscheck(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof)
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]


def test_short_arc_cross_check(leg_short: Leg) -> None:
    """Earth->Earth short arc 50 d: agreement < 1e-3 m/s."""
    res = lambert_crosscheck(leg_short.r1, leg_short.r2, leg_short.tof)
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]


def test_long_arc_cross_check(leg_long: Leg) -> None:
    """Earth->Mars long arc 500 d: agreement < 1e-3 m/s."""
    res = lambert_crosscheck(leg_long.r1, leg_long.r2, leg_long.tof)
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]
