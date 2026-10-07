"""#968: the Jovian lane's periodicity wrap must be rotation-aware.

The wrap residual in :func:`cyclerfinder.nbody.jovian.jovian_defect_residual` and
:func:`cyclerfinder.nbody.jovian_ideal.subarc_defect_residual` compares the end node
and the start node relative to the home moon. A cycler repeats ROTATED by the home
moon's angular advance over the period, so the relative state at the end must be
compared with the start relative state rotated into the moon's new orbital frame.

Fixture (exact by construction, no published value needed): a massless "home moon"
on a circular orbit about Jupiter alone, and a spacecraft on a circular orbit in the
same plane with period 4/7 of the moon's. Over T = (4/3) P_moon the moon advances
480 deg (120 deg mod 360) and the spacecraft 840 deg (also 120 deg mod 360), so the
spacecraft state at t0 + T is exactly the t0 state rotated by the moon's advance:
the orbit is periodic in the moon's rotating frame, with a non-zero advance. A
coplanar case and an inclined (non-equatorial plane) case are both pinned.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from cyclerfinder.nbody.jovian import (
    MU_JUPITER_KM3_S2,
    JovianRailsCache,
    jovian_defect_residual,
)
from cyclerfinder.nbody.jovian_ideal import SubarcSeed, subarc_defect_residual
from cyclerfinder.nbody.jovian_stm import jovian_stm_jacobian, subarc_stm_jacobian
from cyclerfinder.nbody.shooter import ShootingSeed, _states_to_x

Vec = NDArray[np.float64]

A_MOON_KM = 1.0704e6
P_MOON_S = 2.0 * math.pi * math.sqrt(A_MOON_KM**3 / MU_JUPITER_KM3_S2)
A_SC_KM = A_MOON_KM * (4.0 / 7.0) ** (2.0 / 3.0)
T_S = (4.0 / 3.0) * P_MOON_S
T0_S = 1.0e8
HOME = "Ganymede"  # the name only; the moon is massless here (moons=() in the force)


def _plane(incl_rad: float) -> tuple[Vec, Vec]:
    """Orthonormal in-plane basis (e1, e2) of a plane inclined about the x axis."""
    e1 = np.array([1.0, 0.0, 0.0])
    e2 = np.array([0.0, math.cos(incl_rad), math.sin(incl_rad)])
    return e1, e2


class _CircularEphem:
    """One massless moon on a circular orbit in an inclined plane (duck-typed ephemeris)."""

    def __init__(self, incl_rad: float, phase_rad: float) -> None:
        self.e1, self.e2 = _plane(incl_rad)
        self.phase = phase_rad
        self.n = 2.0 * math.pi / P_MOON_S

    def state(self, moon: str, t_sec: float) -> tuple[Vec, Vec]:
        th = self.phase + self.n * (t_sec - T0_S)
        v = A_MOON_KM * self.n
        r = A_MOON_KM * (math.cos(th) * self.e1 + math.sin(th) * self.e2)
        vv = v * (-math.sin(th) * self.e1 + math.cos(th) * self.e2)
        return r, vv


def _sc_state(incl_rad: float, t_sec: float) -> Vec:
    e1, e2 = _plane(incl_rad)
    n = math.sqrt(MU_JUPITER_KM3_S2 / A_SC_KM**3)
    th = 0.3 + n * (t_sec - T0_S)
    r = A_SC_KM * (math.cos(th) * e1 + math.sin(th) * e2)
    v = A_SC_KM * n * (-math.sin(th) * e1 + math.cos(th) * e2)
    return np.concatenate([r, v])


def _case(incl_rad: float) -> tuple[_CircularEphem, ShootingSeed, JovianRailsCache]:
    eph = _CircularEphem(incl_rad, phase_rad=1.1)
    t1 = T0_S + T_S
    zero = np.zeros(3)
    seed = ShootingSeed(
        node_states=[_sc_state(incl_rad, T0_S), _sc_state(incl_rad, t1)],
        epochs=[T0_S, t1],
        tofs=[T_S / 86400.0],
        sequence=(HOME, HOME),
        slack_leg=0,
        period_days=T_S / 86400.0,
        vinf_in=[zero, zero],
        vinf_out=[zero, zero],
    )
    cache = JovianRailsCache((), eph, T0_S, t1)  # type: ignore[arg-type]
    return eph, seed, cache


def _frame(r: Vec, v: Vec) -> NDArray[np.float64]:
    rh = r / np.linalg.norm(r)
    h = np.cross(r, v)
    hh = h / np.linalg.norm(h)
    return np.column_stack([rh, np.cross(hh, rh), hh])


INCLS = [0.0, math.radians(25.0)]


@pytest.mark.parametrize("incl", INCLS)
def test_fixture_is_rotating_frame_periodic_with_nonzero_advance(incl: float) -> None:
    """Independent check of the fixture: translation wrap is large, rotated wrap is zero."""
    eph, seed, _ = _case(incl)
    r0m, v0m = eph.state(HOME, seed.epochs[0])
    r1m, v1m = eph.state(HOME, seed.epochs[1])
    adv = math.degrees(math.acos(float(np.dot(r0m, r1m)) / A_MOON_KM**2))
    assert adv == pytest.approx(120.0, abs=1e-9)
    s0, s1 = seed.node_states
    rel0 = np.concatenate([s0[:3] - r0m, s0[3:] - v0m])
    rel1 = np.concatenate([s1[:3] - r1m, s1[3:] - v1m])
    rot = _frame(r1m, v1m) @ _frame(r0m, v0m).T
    assert float(np.linalg.norm(rel1[:3] - rel0[:3])) > 1.0e5  # km: translation-only fails
    assert float(np.linalg.norm(rel1[:3] - rot @ rel0[:3])) < 1e-6
    assert float(np.linalg.norm(rel1[3:] - rot @ rel0[3:])) < 1e-12


@pytest.mark.parametrize("incl", INCLS)
def test_lane_wrap_is_zero_on_rotating_frame_periodic_orbit(incl: float) -> None:
    eph, seed, cache = _case(incl)
    res = jovian_defect_residual(seed, ephem=eph, cache=cache, moons=())  # type: ignore[arg-type]
    leg, wrap = res[:6], res[-6:]
    assert float(np.linalg.norm(leg)) < 1e-3  # REBOUND Kepler leg closes (sanity)
    assert float(np.linalg.norm(wrap[:3])) < 1e-6
    assert float(np.linalg.norm(wrap[3:])) < 1e-6  # velocity rows carry the 1e3 weight


@pytest.mark.parametrize("incl", INCLS)
def test_subarc_wrap_is_zero_on_rotating_frame_periodic_orbit(incl: float) -> None:
    eph, seed, cache = _case(incl)
    sub = SubarcSeed(
        node_states=list(seed.node_states),
        epochs=list(seed.epochs),
        encounter_idx=(0, 1),
        sequence=seed.sequence,
        vinf_in=list(seed.vinf_in),
        vinf_out=list(seed.vinf_out),
        n_subarcs=1,
    )
    res = subarc_defect_residual(
        sub,
        sub.node_states,
        ephem=eph,  # type: ignore[arg-type]
        cache=cache,
        moons=(),
    )
    assert float(np.linalg.norm(res[-6:])) < 1e-6


@pytest.mark.parametrize("incl", INCLS)
def test_stm_jacobian_wrap_rows_match_residual(incl: float) -> None:
    """The analytic wrap rows equal the exact (linear) derivative of the wrap residual."""
    eph, seed, _ = _case(incl)
    x = _states_to_x(seed.node_states)
    jac = jovian_stm_jacobian(seed, x, ephem=eph, moons=())
    sub = SubarcSeed(
        node_states=list(seed.node_states),
        epochs=list(seed.epochs),
        encounter_idx=(0, 1),
        sequence=seed.sequence,
        vinf_in=list(seed.vinf_in),
        vinf_out=list(seed.vinf_out),
        n_subarcs=1,
    )
    jac_sub = subarc_stm_jacobian(sub, x, ephem=eph, moons=())
    r0m, v0m = eph.state(HOME, seed.epochs[0])
    r1m, v1m = eph.state(HOME, seed.epochs[1])
    rot = _frame(r1m, v1m) @ _frame(r0m, v0m).T
    w = np.diag([1.0, 1.0, 1.0, 1e3, 1e3, 1e3])
    blk = np.zeros((6, 6))
    blk[:3, :3] = rot
    blk[3:, 3:] = rot
    for j in (jac, jac_sub):
        np.testing.assert_allclose(j[-6:, 0:6], -w @ blk, atol=1e-12)
        np.testing.assert_allclose(j[-6:, 6:12], w, atol=1e-12)
