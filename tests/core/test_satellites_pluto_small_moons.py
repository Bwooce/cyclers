"""#1037: the Pluto small-moon registry entries pin to their published source values.

The EXPECTED side below is transcribed from the JPL SSD satellite pages (ssd.jpl.nasa.gov/sats/
phys_par/ and /sats/elem/, ephemeris PLU060; page reference M. Brozovic and R. A. Jacobson 2024,
AJ 167:256), accessed 2026-10-11. It is not computed from the registry. Styx and Kerberos GM are
printed as upper limits ("< 0.0003", "< 0.0002"); the registry stores the limit.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

REPO = Path(__file__).resolve().parents[2]

# name: (GM km^3/s^2, GM is an upper limit, mean radius km, a km, printed period d, printed e)
JPL_PAGES_2026_10_11 = {
    "Styx": (0.0003, True, 5.2, 43200.0, 20.16, 0.025),
    "Nix": (0.0015, False, 18.0, 49300.0, 24.85, 0.015),
    "Kerberos": (0.0002, True, 6.0, 58300.0, 32.17, 0.010),
    "Hydra": (0.0020, False, 18.5, 65200.0, 38.20, 0.009),
}


@pytest.mark.parametrize("moon", sorted(JPL_PAGES_2026_10_11))
def test_registry_matches_jpl_pages(moon: str) -> None:
    gm, _is_limit, radius, sma, _period, _ecc = JPL_PAGES_2026_10_11[moon]
    sat = SATELLITES[moon]
    assert sat.primary == "Pluto"
    assert sat.mu_km3_s2 == gm
    assert sat.radius_eq_km == radius
    assert sat.sma_km == sma
    assert not sat.retrograde


@pytest.mark.parametrize("moon", sorted(JPL_PAGES_2026_10_11))
def test_registry_period_within_documented_offset_of_printed_period(moon: str) -> None:
    """The printed periods are PLU060 mean periods; the registry period follows from the printed
    mean a by Kepler III about the system GM. They differ by 1-4 percent (the printed a is not the
    Kepler radius of the printed period), which this test pins so that a change is noticed."""
    from math import pi, sqrt

    sat = SATELLITES[moon]
    kepler_d = 2.0 * pi * sqrt(sat.sma_km**3 / PRIMARIES["Pluto"]) / 86400.0
    printed = JPL_PAGES_2026_10_11[moon][4]
    assert 0.010 < kepler_d / printed - 1.0 < 0.045


def test_task_998_driver_uses_registry_values() -> None:
    spec = importlib.util.spec_from_file_location(
        "run_998_enumerate_1037", REPO / "scripts" / "run_998_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for moon in ("Styx", "Nix", "Kerberos", "Hydra"):
        gm, _lim, radius, _sma, _per, _e = JPL_PAGES_2026_10_11[moon]
        body = mod.pluto_system(moon).body(moon)
        assert body.mu_km3_s2 == gm
        assert body.radius_km == radius
        assert body.massless
    # the Charon constants the driver carries are the registry's and the JPL page's
    # (GM 106.1 +- 0.3, mean radius 606.0 +- 0.5 km)
    assert mod.CHARON_GM == SATELLITES["Charon"].mu_km3_s2 == 106.1
    assert mod.CHARON_RADIUS == SATELLITES["Charon"].radius_eq_km == 606.0
    assert SATELLITES["Charon"].safe_alt_km == mod.CHARON_FLOOR
