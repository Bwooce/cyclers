"""#1037: the Pluto small-moon registry entries pin to their published source values.

The EXPECTED side below is transcribed from the JPL SSD satellite pages (ssd.jpl.nasa.gov/sats/
phys_par/ and /sats/elem/, ephemeris PLU060; page reference M. Brozovic and R. A. Jacobson 2024,
AJ 167:256), accessed 2026-10-11. It is not computed from the registry. Styx and Kerberos GM are
printed as upper limits ("< 0.0003", "< 0.0002"); the registry stores the limit.
"""

from __future__ import annotations

import importlib.util
import json
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


#: Brozovic and Jacobson (2024), AJ 167:256, Table 11 (mean precessing ellipse): n in deg/day,
#: printed P in days, printed a in km. Transcribed from
#: data/sources/brozovic-jacobson-2024-tables.yaml (table 11); the expected side of the tests
#: below, not computed from the registry.
TABLE_11 = {
    "Styx": (17.855, 20.16, 42410.0),
    "Nix": (14.484, 24.85, 48689.0),
    "Kerberos": (11.191, 32.17, 57749.0),
    "Hydra": (9.424, 38.20, 64719.0),
}
MU_SYSTEM = 975.5  # registry PRIMARIES["Pluto"]


@pytest.mark.parametrize("moon", sorted(TABLE_11))
def test_stored_mean_motion_is_table_11(moon: str) -> None:
    n, printed_p, _a = TABLE_11[moon]
    sat = SATELLITES[moon]
    assert sat.mean_motion_sourced_deg_day == n
    assert sat.mean_motion_deg_day == n
    assert "Table 11" in sat.mean_motion_source
    # P = 360/n to the printed digits: half a unit of the last printed digit (0.005 d) plus the
    # rounding of n itself (0.0005 deg/day is 0.0009 d at Nix)
    assert abs(360.0 / n - printed_p) < 0.006


@pytest.mark.parametrize("moon", sorted(TABLE_11))
def test_kepler_radius_of_the_stored_rate(moon: str) -> None:
    from math import radians

    n, _p, table_a = TABLE_11[moon]
    sat = SATELLITES[moon]
    n_rad_s = radians(n) / 86400.0
    assert sat.kepler_sma_km == pytest.approx((MU_SYSTEM / n_rad_s**2) ** (1.0 / 3.0), rel=1e-12)
    # The paper's own Table 11 a is not the Kepler radius of n: 0.2-0.6 percent above it.
    assert 0.002 < table_a / sat.kepler_sma_km - 1.0 < 0.006


@pytest.mark.parametrize("moon", sorted(TABLE_11))
def test_kepler_period_from_the_printed_osculating_a_stays_longer(moon: str) -> None:
    """sma_km is Table 10's mean osculating a (page value). Its Kepler period is 1-4.5 percent
    longer than the printed (Table 11) period; this documents osculating against mean and
    pins sma_km."""
    from math import pi, sqrt

    sat = SATELLITES[moon]
    kepler_d = 2.0 * pi * sqrt(sat.sma_km**3 / PRIMARIES["Pluto"]) / 86400.0
    printed = JPL_PAGES_2026_10_11[moon][4]
    assert 0.010 < kepler_d / printed - 1.0 < 0.045


def test_orbit_circle_of_the_four_and_of_every_other_moon() -> None:
    import math

    from cyclerfinder.core.satellites import orbit_circle

    for moon, (n, _p, _a) in TABLE_11.items():
        radius, rate = orbit_circle(moon)
        assert radius == SATELLITES[moon].kepler_sma_km
        assert rate == math.radians(n)
        # v = n r agrees with sqrt(mu / r) for the Kepler radius
        v = rate / 86400.0 * radius
        assert v == pytest.approx(math.sqrt(MU_SYSTEM / radius), rel=1e-12)
    # every other moon: bit-identical to the registry before #1022d (fixture recorded from the
    # pre-change code, a regression reference and not a golden of correctness)
    ref = json.loads((REPO / "tests" / "core" / "data" / "registry_before_1022d.json").read_text())
    assert set(ref["moons"]) == set(SATELLITES)
    for moon, before in ref["moons"].items():
        if moon in TABLE_11:
            continue
        assert SATELLITES[moon].sma_km == before["sma_km"]
        assert SATELLITES[moon].mean_motion_deg_day == before["mean_motion_deg_day"]
        assert SATELLITES[moon].kepler_sma_km == before["sma_km"]
        assert orbit_circle(moon) == (before["sma_km"], before["circle_rate_rad_day"])


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


def test_998_driver_periods_agree_with_the_stored_rates() -> None:
    """The #998 driver keeps its plu060-fitted periods (more precise, the model's own input); the
    registry stores the sourced Table 11 rate. Both are kept; neither may drift from the other."""
    spec = importlib.util.spec_from_file_location(
        "run_998_enumerate_1022d", REPO / "scripts" / "run_998_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for moon, (n, _p, _a) in TABLE_11.items():
        assert mod.PERIODS_D[moon] == pytest.approx(360.0 / n, rel=5e-5)
