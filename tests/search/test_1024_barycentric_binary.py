"""#1024: the barycentric circular-model option for binaries (production scales: km, km^3/s^2).

Expected values: the plu060 kernel's mean barycentric Charon radius (17,464 km, #998 note sec. 1,
JPL ephemeris, not this code), the registry GMs, and identities (zero total momentum, the
mass-fraction-0 limit equal to the default model).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.search.two_working_body import CircularSystem, Cycle, LambertLeg, cycle_flybys
from cyclerfinder.search.two_working_body_binary import (
    barycentric_binary,
    primary_min_distance_km,
    with_mass_fraction,
)
from cyclerfinder.search.two_working_body_enum import leg_extent

DAY = 86400.0
MU_SYS = 975.5  # Pluto + Charon system GM (registry PRIMARIES["Pluto"])
GM_CHARON = 106.1  # registry
P_CHARON = 6.38722 * DAY  # fitted to plu060 (#998 note sec. 1)
P_NIX = 24.85472 * DAY


def _kepler_a(period_s: float) -> float:
    return float((MU_SYS * (period_s / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0))


def _pluto() -> CircularSystem:
    bodies = {
        "Charon": (_kepler_a(P_CHARON), P_CHARON, 0.0),
        "Nix": (_kepler_a(P_NIX), P_NIX, 0.0),
    }
    return CircularSystem(MU_SYS, bodies, frozenset({"Nix"}))


def test_mass_fraction_zero_reproduces_default_states() -> None:
    base = _pluto()
    bary = barycentric_binary(base, "Charon", 0.0)
    for t in (0.0, 1.234e5, 3.3e6):
        for code in ("Charon", "Nix"):
            r0, v0 = base.state(code, t)
            r1, v1 = bary.state(code, t)
            assert np.allclose(r1, r0, rtol=1e-12, atol=0.0)
            assert np.allclose(v1, v0, rtol=1e-12, atol=1e-15)
        assert np.allclose(bary.primary_state(t)[0], 0.0)


def test_pluto_charon_barycentric_radius_and_speed() -> None:
    mu_c = GM_CHARON / MU_SYS
    bary = barycentric_binary(_pluto(), "Charon", mu_c)
    r, v = bary.state("Charon", 0.0)
    rn = float(np.linalg.norm(r))
    assert abs(rn - 17_464.0) < 1.0  # plu060 mean barycentric radius
    assert math.isclose(float(np.linalg.norm(v)), 2.0 * math.pi * rn / P_CHARON, rel_tol=1e-14)
    # not Keplerian about the system GM at that radius
    assert float(np.linalg.norm(v)) < 0.9 * math.sqrt(MU_SYS / rn)
    assert not bary.is_keplerian("Charon") and bary.is_keplerian("Nix")
    # the small moon is unchanged
    assert np.array_equal(bary.state("Nix", 7.0e5)[0], _pluto().state("Nix", 7.0e5)[0])


def test_primary_opposite_with_zero_total_momentum() -> None:
    gm_pluto = MU_SYS - GM_CHARON
    bary = barycentric_binary(_pluto(), "Charon", GM_CHARON / MU_SYS)
    for t in (0.0, 2.0e5, 1.1e6):
        rc, vc = bary.state("Charon", t)
        rp, vp = bary.primary_state(t)
        assert abs(float(np.linalg.norm(rp)) - 2131.5) < 1.0
        assert np.allclose(gm_pluto * rp + GM_CHARON * rc, 0.0, atol=1e-6)
        assert np.allclose(gm_pluto * vp + GM_CHARON * vc, 0.0, atol=1e-9)


def test_with_mass_fraction_ramp_endpoints() -> None:
    base = _pluto()
    a_rel = base.bodies["Charon"][0]
    mu_c = GM_CHARON / MU_SYS
    one = barycentric_binary(base, "Charon", mu_c)
    zero = with_mass_fraction(one, "Charon", a_rel, 0.0)
    half = with_mass_fraction(one, "Charon", a_rel, 0.5 * mu_c)
    assert zero.bodies["Charon"][0] == a_rel
    assert math.isclose(half.bodies["Charon"][0], (1.0 - 0.5 * mu_c) * a_rel, rel_tol=1e-15)
    assert half.kinematic == frozenset({"Charon"})


def test_primary_min_distance_equals_r_min_at_mass_fraction_zero() -> None:
    """With the primary at the origin, the closest approach to it is the cycle's r_min."""
    base = _pluto()
    bary = barycentric_binary(base, "Charon", 0.0)
    syn = 1.0 / abs(1.0 / P_CHARON - 1.0 / P_NIX)
    # a plain two-Lambert-leg cycle with dates far from any degenerate geometry
    legs = (LambertLeg("Charon", "Nix", 0, "single"), LambertLeg("Nix", "Charon", 0, "single"))
    cyc = Cycle(legs, 3 * syn)
    x = np.array([0.37 * DAY, 9.1 * DAY])
    fl = cycle_flybys(bary, cyc, x)
    if fl is None:
        pytest.skip("no flyby directions at these dates")
    r_min, _ = leg_extent(bary, cyc, x, flybys=fl)
    pa = primary_min_distance_km(bary, cyc, x, fl)
    assert abs(pa.min_km - r_min) < 1e-3 * r_min
