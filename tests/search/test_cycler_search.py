from __future__ import annotations

import math

from cyclerfinder.search.cycler_search import (
    aphelion_ratio,
    cycler_tof,
    generic_return_aphelion,
    max_earth_flyby_bend,
    turn_ratio,
)
from cyclerfinder.search.generic_return import RussellModel, returns_at_vinf


def test_cycler_tof_eq_3_1() -> None:
    m = RussellModel()
    tau = m.synodic_yr("E", "M")
    assert abs(cycler_tof(m, p=4, h=3, s=1) - (tau * 4 - 3 / 2) / 1) < 1e-9
    assert abs(cycler_tof(m, p=4, h=9, s=2) - (tau * 4 - 9 / 2) / 2) < 1e-9


def test_aphelion_ratio() -> None:
    assert abs(aphelion_ratio(1.64) - 1.64 / 1.52) < 1e-12


def test_turn_ratio() -> None:
    assert abs(turn_ratio(max_allowable=1.0, omega_max=0.5) - 2.0) < 1e-12


def test_max_earth_flyby_bend_positive() -> None:
    m = RussellModel()
    b = max_earth_flyby_bend(m, vinf=0.2)  # canonical AU/TU
    assert 0.0 < b < math.pi


def test_generic_return_aphelion_reconstructs() -> None:
    m = RussellModel()
    rs = returns_at_vinf(m, "E", 0.1838, dtheta_deg=1.0, max_revs_cap=6)
    assert rs
    aph = generic_return_aphelion(m, "E", rs[0])
    assert aph > 0.0  # AU
