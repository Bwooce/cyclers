"""Structural and metric tests for the M3 model dataclasses.

No physics integration here — hand-crafted minimal cyclers exercise the
dataclass invariants (frozenness, defaults) and the four metric methods
on :class:`Cycler`. The physics-end-to-end tests live in
``test_aldrin.py`` and ``test_construct_2syn_em.py``.

Plan: ``docs/phases/m3-model-construct/plan.md`` §4.2.
"""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from math import sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.frames import synodic_omega
from cyclerfinder.model import Cycler, Encounter, Leg

TOL_DV_ZERO_KMS: float = 1.0e-12
TOL_RADIAL_AU: float = 1.0e-6


def _zero3() -> np.ndarray:
    return np.zeros(3, dtype=np.float64)


def _v(x: float, y: float, z: float) -> np.ndarray:
    return np.array([x, y, z], dtype=np.float64)


def test_encounter_is_frozen() -> None:
    enc = Encounter(
        body="E",
        t=0.0,
        r=_zero3(),
        v_planet=_zero3(),
        vinf_in=_zero3(),
        vinf_out=_zero3(),
    )
    with pytest.raises(FrozenInstanceError):
        enc.body = "M"  # type: ignore[misc]


def test_leg_is_frozen() -> None:
    leg = Leg(
        from_body="E",
        to_body="M",
        t_depart=0.0,
        t_arrive=1.0,
        v_depart=_zero3(),
        v_arrive=_zero3(),
    )
    with pytest.raises(FrozenInstanceError):
        leg.t_depart = 99.0  # type: ignore[misc]


def test_leg_defaults() -> None:
    leg = Leg(
        from_body="E",
        to_body="M",
        t_depart=0.0,
        t_arrive=1.0,
        v_depart=_zero3(),
        v_arrive=_zero3(),
    )
    assert leg.n_revs == 0
    assert leg.branch == "single"


def test_cycler_is_frozen() -> None:
    cyc = Cycler(bodies=["E"], period=0.0, encounters=[], legs=[])
    with pytest.raises(FrozenInstanceError):
        cyc.period = 1.0  # type: ignore[misc]


def test_cycler_max_vinf_returns_largest_magnitude() -> None:
    encs = [
        Encounter("E", 0.0, _zero3(), _zero3(), _v(3.0, 0.0, 0.0), _v(3.0, 0.0, 0.0)),
        Encounter("M", 1.0, _zero3(), _zero3(), _v(0.0, 5.0, 0.0), _v(0.0, 5.0, 0.0)),
        Encounter("E", 2.0, _zero3(), _zero3(), _v(0.0, 0.0, 4.0), _v(0.0, 0.0, 4.0)),
    ]
    cyc = Cycler(bodies=["E", "M", "E"], period=2.0, encounters=encs, legs=[])
    assert cyc.max_vinf() == pytest.approx(5.0)


def test_cycler_maintenance_dv_zero_when_vinf_in_equals_vinf_out() -> None:
    encs = [
        Encounter("E", 0.0, _zero3(), _zero3(), _v(1.0, 2.0, 3.0), _v(1.0, 2.0, 3.0)),
        Encounter("M", 1.0, _zero3(), _zero3(), _v(0.5, -0.5, 0.1), _v(0.5, -0.5, 0.1)),
    ]
    cyc = Cycler(bodies=["E", "M"], period=1.0, encounters=encs, legs=[])
    assert cyc.maintenance_dv() < TOL_DV_ZERO_KMS


def test_cycler_maintenance_dv_sums_discontinuities() -> None:
    encs = [
        Encounter("E", 0.0, _zero3(), _zero3(), _v(1.0, 0.0, 0.0), _v(2.0, 0.0, 0.0)),
        Encounter("M", 1.0, _zero3(), _zero3(), _v(0.0, 3.0, 0.0), _v(0.0, 7.0, 0.0)),
    ]
    cyc = Cycler(bodies=["E", "M"], period=1.0, encounters=encs, legs=[])
    # |2-1| + |7-3| = 1 + 4 = 5.
    assert cyc.maintenance_dv() == pytest.approx(5.0)


def test_cycler_radial_span_hohmann() -> None:
    """A circular departure -> arrival Hohmann transfer has known peri/apo.

    Set r_dep = 1 AU, r_arr = 2 AU; the Hohmann transfer ellipse has
    semi-major axis 1.5 AU, perihelion 1 AU, aphelion 2 AU. The departure
    velocity magnitude is ``sqrt(mu * (2/r_dep - 1/a))``.
    """
    r_dep_km = 1.0 * AU_KM
    r_arr_km = 2.0 * AU_KM
    a_km = 0.5 * (r_dep_km + r_arr_km)
    v_dep_mag = sqrt(MU_SUN_KM3_S2 * (2.0 / r_dep_km - 1.0 / a_km))
    v_arr_mag = sqrt(MU_SUN_KM3_S2 * (2.0 / r_arr_km - 1.0 / a_km))

    enc_dep = Encounter(
        "E",
        0.0,
        r=_v(r_dep_km, 0.0, 0.0),
        v_planet=_zero3(),
        vinf_in=_zero3(),
        vinf_out=_zero3(),
    )
    enc_arr = Encounter(
        "M",
        1.0e7,
        r=_v(-r_arr_km, 0.0, 0.0),
        v_planet=_zero3(),
        vinf_in=_zero3(),
        vinf_out=_zero3(),
    )
    leg = Leg(
        from_body="E",
        to_body="M",
        t_depart=0.0,
        t_arrive=1.0e7,
        v_depart=_v(0.0, v_dep_mag, 0.0),
        v_arrive=_v(0.0, -v_arr_mag, 0.0),
    )
    cyc = Cycler(bodies=["E", "M"], period=1.0e7, encounters=[enc_dep, enc_arr], legs=[leg])
    peri_au, apo_au = cyc.radial_span()
    assert abs(peri_au - 1.0) < TOL_RADIAL_AU
    assert abs(apo_au - 2.0) < TOL_RADIAL_AU


def test_closure_residual_uses_default_earth_omega() -> None:
    """``omega=None`` must match ``omega=synodic_omega('E')`` numerically."""
    # Minimal but representative cycler: two encounters at separated epochs
    # with non-trivial r and v_depart so the rotating-frame transform does
    # real work.
    enc0 = Encounter(
        "E",
        0.0,
        r=_v(1.0 * AU_KM, 0.0, 0.0),
        v_planet=_v(0.0, 29.78, 0.0),
        vinf_in=_v(0.0, 5.0, 0.0),
        vinf_out=_v(0.0, 5.0, 0.0),
    )
    enc1 = Encounter(
        "M",
        146.0 * 86400.0,
        r=_v(0.5 * AU_KM, 1.3 * AU_KM, 0.0),
        v_planet=_v(-22.0, 10.0, 0.0),
        vinf_in=_v(0.0, 6.0, 0.0),
        vinf_out=_v(0.0, 6.0, 0.0),
    )
    leg = Leg(
        from_body="E",
        to_body="M",
        t_depart=0.0,
        t_arrive=146.0 * 86400.0,
        v_depart=_v(0.0, 31.0, 0.0),
        v_arrive=_v(-21.0, 12.0, 0.0),
    )
    cyc = Cycler(bodies=["E", "M"], period=146.0 * 86400.0, encounters=[enc0, enc1], legs=[leg])
    r_default = cyc.closure_residual(omega_rad_per_s=None)
    r_explicit = cyc.closure_residual(omega_rad_per_s=synodic_omega("E"))
    assert r_default == pytest.approx(r_explicit)


def test_closure_residual_is_finite_and_nonnegative() -> None:
    enc0 = Encounter(
        "E",
        0.0,
        r=_v(1.0 * AU_KM, 0.0, 0.0),
        v_planet=_zero3(),
        vinf_in=_zero3(),
        vinf_out=_zero3(),
    )
    enc1 = Encounter(
        "E",
        1.0e7,
        r=_v(1.0 * AU_KM, 0.0, 0.0),
        v_planet=_zero3(),
        vinf_in=_zero3(),
        vinf_out=_zero3(),
    )
    leg = Leg(
        from_body="E",
        to_body="E",
        t_depart=0.0,
        t_arrive=1.0e7,
        v_depart=_v(0.0, 30.0, 0.0),
        v_arrive=_v(0.0, 30.0, 0.0),
    )
    cyc = Cycler(bodies=["E", "E"], period=1.0e7, encounters=[enc0, enc1], legs=[leg])
    r = cyc.closure_residual()
    assert r >= 0.0
    import math

    assert math.isfinite(r)
