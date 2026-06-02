"""Tests for the inclined circular ephemeris backend (3D, STAGE 2).

The analytic ``circular`` backend gains orbital inclination and longitude of
ascending node so Venus (~3.39 deg) and Mars (~1.85 deg) carry real 3D states.
The ``inc_deg == 0.0`` path is exact-float short-circuited to the flat code, so
existing planar behaviour is byte-identical.

SOURCED orbital elements: Standish & Williams, "Approximate Positions of the
Planets", JPL Solar System Dynamics, Table 1 (J2000 elements, 1800-2050 AD).
https://ssd.jpl.nasa.gov/planets/approx_pos.html
"""

from __future__ import annotations

from dataclasses import replace
from math import cos, pi, radians, sin, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SECONDS_PER_DAY,
    PlanetData,
)
from cyclerfinder.core.ephemeris import (
    Ephemeris,
    _CircularBackend,
    _InclinedCircularBackend,
)

# SOURCED — Standish & Williams Table 1 (J2000).
_VENUS_INC_DEG = 3.39467605
_VENUS_LAN_DEG = 76.67984255
_MARS_INC_DEG = 1.84969142
_MARS_LAN_DEG = 49.55953891


def _period_s(planet: PlanetData) -> float:
    a_km = planet.sma_au * AU_KM
    return 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)


def _inclined_planet(code: str, inc_deg: float, lan_deg: float) -> PlanetData:
    """A copy of a PLANETS record with inclination/node set (test-local)."""
    return replace(PLANETS[code], inc_deg=inc_deg, lan_deg=lan_deg)


def test_inclined_circular_i0_byte_identical_to_flat() -> None:
    """inc_deg=0 inclined path == flat path by exact numpy equality. # INVARIANT."""
    flat = _CircularBackend()
    incl = _InclinedCircularBackend()
    for body in PLANETS:
        for t_days in (0.0, 73.0, 365.0):
            r_flat, v_flat = flat.state(body, t_days * SECONDS_PER_DAY)
            r_inc, v_inc = incl.state(body, t_days * SECONDS_PER_DAY)
            assert np.array_equal(r_flat, r_inc)
            assert np.array_equal(v_flat, v_inc)


def test_existing_planar_tests_unchanged() -> None:
    """i=0 bodies stay exactly in the ecliptic via the public backend. # INVARIANT."""
    eph = Ephemeris(model="circular")
    for body in PLANETS:
        for t_days in (0.0, 73.0, 365.0, 1000.0):
            r, v = eph.state(body, t_days * SECONDS_PER_DAY)
            assert r[2] == 0.0
            assert v[2] == 0.0


def test_venus_inclined_z_nonzero() -> None:
    """Venus (SOURCED inc/lan) has z != 0; |z| ~ sma*sin(inc) at quarter-period.

    SOURCED elements: Standish & Williams Table 1 (J2000).
    """
    venus = _inclined_planet("V", _VENUS_INC_DEG, _VENUS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"V": venus})
    t_quarter = _period_s(venus) / 4.0
    r, _ = backend.state("V", t_quarter)
    a_km = venus.sma_au * AU_KM
    z_max = a_km * sin(radians(_VENUS_INC_DEG))  # ~3.0e6 km
    assert abs(r[2]) > 1.0e5
    assert abs(r[2]) == pytest.approx(z_max, rel=1.0e-9)


def test_mars_inclined_z_nonzero() -> None:
    """Mars (SOURCED inc/lan) has z != 0. SOURCED Standish & Williams Table 1."""
    mars = _inclined_planet("M", _MARS_INC_DEG, _MARS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"M": mars})
    t_quarter = _period_s(mars) / 4.0
    r, _ = backend.state("M", t_quarter)
    assert abs(r[2]) > 1.0e5


def test_inclined_circular_speed_preserved() -> None:
    """|v| == sqrt(mu_sun / |r|): inclination preserves orbital energy. # INVARIANT."""
    venus = _inclined_planet("V", _VENUS_INC_DEG, _VENUS_LAN_DEG)
    mars = _inclined_planet("M", _MARS_INC_DEG, _MARS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"V": venus, "M": mars})
    for code in ("V", "M"):
        period = _period_s(backend._planets[code])
        for frac in (0.0, 0.13, 0.5, 0.77):
            r, v = backend.state(code, frac * period)
            r_n = float(np.linalg.norm(r))
            v_n = float(np.linalg.norm(v))
            assert v_n == pytest.approx(sqrt(MU_SUN_KM3_S2 / r_n), rel=1.0e-12)


def test_inclined_circular_angular_momentum_direction() -> None:
    """h = r x v parallel to orbit normal n_hat within 1e-12 rel. # INVARIANT.

    n_hat = (-sin(lan)*sin(inc), cos(lan)*sin(inc), cos(inc)). Catches a wrong
    rotation order/sign (which would flip h).
    """
    venus = _inclined_planet("V", _VENUS_INC_DEG, _VENUS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"V": venus})
    inc = radians(_VENUS_INC_DEG)
    lan = radians(_VENUS_LAN_DEG)
    n_hat = np.array([-sin(lan) * sin(inc), cos(lan) * sin(inc), cos(inc)], dtype=np.float64)
    period = _period_s(venus)
    for frac in (0.1, 0.35, 0.6, 0.9):
        r, v = backend.state("V", frac * period)
        h = np.cross(r, v)
        h_hat = h / np.linalg.norm(h)
        # Same direction (not anti-parallel): dot == +1.
        assert float(np.dot(h_hat, n_hat)) == pytest.approx(1.0, abs=1.0e-12)


def test_inclined_circular_period_closes() -> None:
    """Within 1 km of t=0 position after exactly one period. # INVARIANT."""
    venus = _inclined_planet("V", _VENUS_INC_DEG, _VENUS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"V": venus})
    period = _period_s(venus)
    r0, _ = backend.state("V", 0.0)
    r1, _ = backend.state("V", period)
    assert float(np.linalg.norm(r1 - r0)) < 1.0


def test_inclined_venus_node_in_ecliptic_plane() -> None:
    """At the ascending node (in-plane angle 0) the body lies in the ecliptic.

    SOURCED geometric consequence of the orbital-element definition: at t=0 the
    in-plane state lies along the node line (reference direction), so after the
    R_x rotation the z-component is exactly 0.
    """
    venus = _inclined_planet("V", _VENUS_INC_DEG, _VENUS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"V": venus})
    r, _ = backend.state("V", 0.0)
    assert abs(r[2]) < 1.0e-3


def test_ephemeris_model_property_still_circular() -> None:
    """No public interface change. # INVARIANT."""
    assert Ephemeris(model="circular").model == "circular"


def test_astropy_backend_unaffected() -> None:
    """Smoke: astropy backend still instantiates and returns a 3-tuple state."""
    ephem = Ephemeris(model="astropy")
    r, v = ephem.state("E", 0.0)
    assert r.shape == (3,)
    assert v.shape == (3,)
    # Earth ~ 1 AU from the Sun.
    assert float(np.linalg.norm(r)) == pytest.approx(AU_KM, rel=0.05)
