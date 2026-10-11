"""#1022d: consumers of the registry's moon circle (campaign constants, moon_circular) use the
stored mean motion for the four Pluto small moons and are bit-identical for every other moon.

Expected values for the four come from Brozovic and Jacobson (2024) Table 11 (n in deg/day);
the bitwise reference for the other moons is tests/core/data/registry_before_1022d.json, recorded
from the pre-change code (a regression reference, not a golden of correctness).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import RepeatedMoonTarget
from cyclerfinder.search.two_working_body import moon_circular

REPO = Path(__file__).resolve().parents[2]
TABLE_11_N = {"Styx": 17.855, "Nix": 14.484, "Kerberos": 11.191, "Hydra": 9.424}
MU = 975.5


def _before() -> dict[str, dict[str, float]]:
    ref = json.loads((REPO / "tests/core/data/registry_before_1022d.json").read_text())
    return dict(ref["moons"])


@pytest.mark.parametrize("primary", ["Jupiter", "Saturn", "Uranus", "Neptune", "Mars"])
def test_campaign_constants_unchanged_for_other_primaries(primary: str) -> None:
    ref = _before()
    consts = RepeatedMoonTarget(primary=primary)._moon_consts()
    assert consts
    for moon, (a, n) in consts.items():
        assert (a, n) == (ref[moon]["sma_km"], ref[moon]["circle_rate_rad_day"])


def test_campaign_constants_for_pluto_use_the_stored_rates() -> None:
    consts = RepeatedMoonTarget(primary="Pluto")._moon_consts()
    assert set(consts) == {"Charon", *TABLE_11_N}
    ref = _before()
    assert consts["Charon"] == (ref["Charon"]["sma_km"], ref["Charon"]["circle_rate_rad_day"])
    for moon, n_deg in TABLE_11_N.items():
        a, n = consts[moon]
        n_rad_s = math.radians(n_deg) / 86400.0
        assert n == math.radians(n_deg)
        assert a == pytest.approx((MU / n_rad_s**2) ** (1.0 / 3.0), rel=1e-12)


def test_moon_circular_unchanged_for_europa_and_sourced_for_nix() -> None:
    sysj = moon_circular(PRIMARIES["Jupiter"], ["Europa"])
    a = SATELLITES["Europa"].sma_km
    assert sysj.bodies["Europa"] == (a, 2.0 * math.pi * math.sqrt(a**3 / PRIMARIES["Jupiter"]), 0.0)
    sysp = moon_circular(MU, ["Nix"])
    a_nix, period, theta0 = sysp.bodies["Nix"]
    assert period == pytest.approx(360.0 / TABLE_11_N["Nix"] * 86400.0, rel=1e-15)
    assert a_nix == pytest.approx((MU * (period / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0), rel=1e-15)
    assert theta0 == 0.0
