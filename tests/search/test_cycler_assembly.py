# ruff: noqa: N802, N806, RUF059 -- Russell's v_F/v_B/a_E notation (verbatim spec test)
from __future__ import annotations

import math

import pytest

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


@pytest.mark.slow
def test_golden_cycler_4_9_2_m1_grouping() -> None:
    """Russell Cycler 4.9.2.-1 worked example (dissertation p.81): h=9,s=2 groups
    as (9,0); first return needs f=6 re-initiating flybys, second needs f=1."""
    import numpy as np

    from cyclerfinder.search.cycler_assembly import f_count, group_half_years
    from cyclerfinder.search.generic_return import RussellModel

    m = RussellModel()
    # a representative incoming v∞ for the grouping decision (direction immaterial to the
    # {9,0} split, which is driven by the omega_c >= omega_minimax test piling onto group 1)
    vinf = 0.1838
    vinf_in = np.array([vinf * 0.3, vinf * 0.95, 0.0])
    hs, omega_max = group_half_years(m, "E", vinf, vinf_in, h=9, s=2)
    assert hs == (9, 0)
    assert f_count(9) == 6  # six re-initiating flybys (the [83,45,45,45,45,83] list)
    assert f_count(0) == 1  # single 24-deg flyby
    assert omega_max > 0.0


def test_descriptor_to_phsi_maps_known_row() -> None:
    from cyclerfinder.data.catalog import load_catalog
    from cyclerfinder.search.cycler_assembly import descriptor_to_phsi

    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    spec = descriptor_to_phsi(row)
    assert spec is not None
    assert spec.s >= 1
    assert spec.vinf_e_kms > 0.0
    assert spec.vinf_m_kms > 0.0
    assert spec.sequence == ("E", "E", "M", "M")


def test_descriptor_to_phsi_none_for_ocampo() -> None:
    from cyclerfinder.data.catalog import load_catalog
    from cyclerfinder.search.cycler_assembly import descriptor_to_phsi

    cat = load_catalog()
    ocampo = next(e for e in cat.entries if e.id.startswith("russell-ocampo"))
    assert descriptor_to_phsi(ocampo.raw) is None
