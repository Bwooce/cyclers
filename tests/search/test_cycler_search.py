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


def test_search_small_returns_gated_cyclers() -> None:
    from cyclerfinder.search.cycler_search import search_cyclers
    from cyclerfinder.search.generic_return import RussellModel

    cs = search_cyclers(RussellModel(), p_max=2, ar_min=0.9, tr_min=0.85, max_revs_cap=6)
    assert isinstance(cs, list)
    for c in cs:
        assert c.ar > 0.9 and c.tr > 0.85
        assert c.p >= 1 and c.s >= 1


def test_generic_returns_at_tof_single_solve() -> None:
    from cyclerfinder.search.cycler_search import generic_returns_at_tof
    from cyclerfinder.search.generic_return import RussellModel

    rs = generic_returns_at_tof(RussellModel(), "E", 3.0, max_revs_cap=4)
    assert all(r.vinf > 0 for r in rs)
