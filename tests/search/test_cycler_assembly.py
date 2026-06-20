# ruff: noqa: N802, N806, RUF059 -- Russell's v_F/v_B/a_E notation (verbatim spec test)
from __future__ import annotations

import math

from cyclerfinder.search.cycler_assembly import (
    full_rev_circle_z,
    full_rev_feasible_vF,
    half_rev_components,
)
from cyclerfinder.search.generic_return import RussellModel


def test_full_rev_feasible_vF_resonant() -> None:
    m = RussellModel()
    # N=M=1 (1:1 resonance) at Earth -> v_F == v_B (Earth circular speed): a_F = a_B
    vF = full_rev_feasible_vF(m, "E", n=1, big_m=1)
    assert abs(vF - m.body_circular_speed("E")) < 1e-9


def test_full_rev_circle_z_eq_2_17() -> None:
    m = RussellModel()
    vB = m.body_circular_speed("E")
    vinf = 0.1838
    vF = vB
    z = full_rev_circle_z(vF, vinf, vB)
    assert abs(z - (vF**2 - vinf**2 - vB**2) / (2 * vB)) < 1e-12


def test_half_rev_components_eq_2_18_2_19() -> None:
    m = RussellModel()
    aE = m.sma_au("E")
    # r1 = r2 = aE (same-body half-rev), transfer sma a = aE -> v_Hr should be ~0 (a=a_min)
    vhr, vht = half_rev_components(m, "E", a=aE)
    expected_vht = math.sqrt(2 * m.mu_sun * aE / (aE**2 + aE * aE))
    assert abs(vht - expected_vht) < 1e-9
    assert vht > 0


def test_f_count_table_3_2() -> None:
    from cyclerfinder.search.cycler_assembly import f_count

    assert [f_count(h) for h in range(0, 9)] == [1, 2, 2, 2, 3, 4, 4, 4, 5]


def test_omega_minimax_branches() -> None:
    import numpy as np

    from cyclerfinder.search.cycler_assembly import omega_minimax
    from cyclerfinder.search.generic_return import RussellModel

    m = RussellModel()
    vinf = 0.1838
    vinf_in = np.array([vinf, 0.0, 0.0])  # incoming v∞ (canonical)
    w1 = omega_minimax(m, "E", vinf, vinf_in, f_j=1)
    w2 = omega_minimax(m, "E", vinf, vinf_in, f_j=2)
    w3 = omega_minimax(m, "E", vinf, vinf_in, f_j=3)
    for w in (w1, w2, w3):
        assert 0.0 < w < np.pi


def test_group_half_years_sums_to_h() -> None:
    import numpy as np

    from cyclerfinder.search.cycler_assembly import group_half_years
    from cyclerfinder.search.generic_return import RussellModel

    hs, omega_max = group_half_years(
        RussellModel(), "E", 0.1838, np.array([0.1838, 0.0, 0.0]), h=10, s=3
    )
    assert sum(hs) == 10
    assert len(hs) == 3
    assert omega_max > 0.0
