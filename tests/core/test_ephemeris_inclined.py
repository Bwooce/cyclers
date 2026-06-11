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
from numpy.typing import NDArray

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

    INDEPENDENTLY-DERIVED expected (Standish ascending-node convention):
    n_hat = (sin(lan)*sin(inc), -cos(lan)*sin(inc), cos(inc)).
    Physical check: polar orbit i=90 deg, lan=0 -> n_hat = (0, -1, 0), so at the
    node (r along +x) v = n_hat x r_hat = +z: the body crosses the ecliptic
    going NORTH — the definition of the ASCENDING node. The previous expected
    value here was derived with the SAME mirrored R_x(-inc) the backend used
    (a circular test), which hid a plane-mirror bug (normals 2*inc from DE440).
    """
    venus = _inclined_planet("V", _VENUS_INC_DEG, _VENUS_LAN_DEG)
    backend = _InclinedCircularBackend(planets={"V": venus})
    inc = radians(_VENUS_INC_DEG)
    lan = radians(_VENUS_LAN_DEG)
    n_hat = np.array([sin(lan) * sin(inc), -cos(lan) * sin(inc), cos(inc)], dtype=np.float64)
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


def _de440_orbit_normal(body: str, t_sec: float) -> NDArray[np.float64]:
    """Unit orbit normal h_hat = (r x v)/|r x v| from the DE440 backend."""
    r, v = Ephemeris(model="astropy").state(body, t_sec)
    h = np.cross(r, v)
    return np.asarray(h / np.linalg.norm(h), dtype=np.float64)


def _angle_deg(a: NDArray[np.float64], b: NDArray[np.float64]) -> float:
    return float(np.degrees(np.arccos(np.clip(float(np.dot(a, b)), -1.0, 1.0))))


@pytest.mark.parametrize(
    ("body", "inc_deg", "lan_deg"),
    [("V", _VENUS_INC_DEG, _VENUS_LAN_DEG), ("M", _MARS_INC_DEG, _MARS_LAN_DEG)],
)
def test_inclined_orbit_normal_anchored_to_de440(body: str, inc_deg: float, lan_deg: float) -> None:
    """ANCHOR: inclined-backend orbit normal matches DE440's to <= 0.01 deg.

    The DE440 normal is computed INDEPENDENTLY (h = r x v from the astropy
    backend) — no shared rotation code. Under the pre-fix mirrored R_x(-inc)
    the disagreement was ~2*inc (Venus 6.7893 deg, Mars 3.6995 deg), so this
    anchor pins the ascending/descending-node sign unambiguously.
    """
    planet = _inclined_planet(body, inc_deg, lan_deg)
    backend = _InclinedCircularBackend(planets={body: planet})
    period = _period_s(planet)
    r0, v0 = backend.state(body, 0.13 * period)
    h = np.cross(r0, v0)
    n_model = h / np.linalg.norm(h)
    # The DE440 normal direction is (very nearly) constant; one epoch suffices.
    n_de440 = _de440_orbit_normal(body, 100.0 * SECONDS_PER_DAY)
    assert _angle_deg(np.asarray(n_model, dtype=np.float64), n_de440) < 0.01


def test_continuation_ramp_lambda1_plane_matches_de440() -> None:
    """REGRESSION (2i tilt defect): ramped backend at lam=1 matches DE440's plane.

    At lam_i=1 the continuation's _RampedElementsBackend tilts Mars by its J2000
    mean inclination about its mean node, so its orbit normal must agree with
    DE440's to a small fraction of inc (1.8497 deg). Under the pre-fix mirrored
    R_x(-i) in _tilt the angle was 2i = 3.6995 deg — defeating the homotopy ramp
    (the "small" final mean-elements -> DE440 step was the full plane flip).
    Threshold 0.05 deg ~ inc/37: fails loudly at 2i, passes the fixed geometry.
    """
    from cyclerfinder.search.continuation import ramped_ephemeris

    ramp = ramped_ephemeris(1.0, 1.0, 1.0)
    r, v = ramp.state("M", 123.0 * SECONDS_PER_DAY)
    h = np.cross(r, v)
    n_ramp = h / np.linalg.norm(h)
    n_de440 = _de440_orbit_normal("M", 123.0 * SECONDS_PER_DAY)
    angle = _angle_deg(np.asarray(n_ramp, dtype=np.float64), n_de440)
    inc_deg = 1.84972647778  # Russell Table 5.4 J2000 mean inclination (lam=1 target)
    assert angle < 0.05
    assert angle < inc_deg / 30.0  # "much less than inc": the 2i mirror is 75x over
