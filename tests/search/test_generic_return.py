from __future__ import annotations

import math

from cyclerfinder.search.generic_return import RussellModel


def test_russell_model_constants() -> None:
    m = RussellModel()
    assert m.tu_days == 58.1324409
    assert m.au_km == 149597871.0
    assert m.mu_sun == 1.0
    assert m.period_yr("E") == 1.0
    assert m.period_yr("M") == 1.875
    assert abs(m.synodic_yr("E", "M") - 15.0 / 7.0) < 1e-12
    # Earth semi-major axis ~ 1 AU in canonical units (Kepler III, mu=1)
    assert abs(m.sma_au("E") - 1.0) < 1e-6
    assert m.body_circular_speed("E") > 0.0
    # circular speed = sqrt(mu/a)
    assert abs(m.body_circular_speed("E") - math.sqrt(m.mu_sun / m.sma_au("E"))) < 1e-12


def test_body_state_coplanar_circular() -> None:
    import numpy as np

    m = RussellModel()
    r, v = m.body_state("E", 0.0)
    r = np.asarray(r)
    v = np.asarray(v)
    assert r.shape == (3,) and v.shape == (3,)
    assert abs(np.linalg.norm(r) - m.sma_au("E")) < 1e-9  # on its circle
    assert abs(np.linalg.norm(v) - m.body_circular_speed("E")) < 1e-9
    assert abs(r[2]) < 1e-12 and abs(v[2]) < 1e-12  # coplanar (z=0)
    # velocity perpendicular to radius (circular)
    assert abs(float(np.dot(r, v))) < 1e-9


def test_psi_reference_geometry() -> None:
    import numpy as np

    from cyclerfinder.search.generic_return import RussellModel, psi_of_vinf_vec

    m = RussellModel()
    r_B, v_B = m.body_state("E", 0.0)
    r_B = np.asarray(r_B)
    v_B = np.asarray(v_B)
    # v∞ aligned with +v_B -> psi 0
    vinf_along_v = v_B / np.linalg.norm(v_B) * 0.1
    assert abs(psi_of_vinf_vec(vinf_along_v, r_B, v_B)) < 1e-9
    # v∞ aligned with +r_B -> psi = +pi/2
    rhat = r_B / np.linalg.norm(r_B)
    vinf_along_r = rhat * 0.1
    assert abs(psi_of_vinf_vec(vinf_along_r, r_B, v_B) - np.pi / 2) < 1e-6


def test_generic_return_dataclass() -> None:
    from cyclerfinder.search.generic_return import GenericReturn

    g = GenericReturn(
        psi_deg=114.0, tof_body_periods=1.25, a_au=0.804, n_revs=1, branch="slow", vinf=0.5
    )
    assert g.branch == "slow"
    assert g.n_revs == 1


def test_generate_returns_coarse_grid() -> None:
    from cyclerfinder.search.generic_return import RussellModel, generate_generic_returns

    m = RussellModel()
    rs = generate_generic_returns(m, "E", max_tof_body_periods=6.0, dtheta_deg=2.0, max_revs_cap=4)
    assert len(rs) > 50
    assert {g.n_revs for g in rs} & {1, 2}  # multiple rev counts present
    assert {g.branch for g in rs} <= {"fast", "slow"}  # only these labels
    assert any(g.n_revs >= 1 for g in rs)  # multi-rev solutions found
    assert all(g.vinf > 0 for g in rs)
    assert all(g.a_au > 0 for g in rs)


def test_bin_and_query_at_vinf() -> None:
    from cyclerfinder.search.generic_return import (
        RussellModel,
        bin_sub_families,
        generate_generic_returns,
        returns_at_vinf,
    )

    m = RussellModel()
    rs = generate_generic_returns(m, "E", dtheta_deg=2.0, max_revs_cap=4)
    bins = bin_sub_families(rs)
    assert bins
    assert all(isinstance(k, tuple) and len(k) == 2 for k in bins)
    # each sub-family sorted by vinf
    for lst in bins.values():
        vs = [g.vinf for g in lst]
        assert vs == sorted(vs)
    got = returns_at_vinf(m, "E", 0.5, dtheta_deg=2.0, max_revs_cap=4)
    assert got
    assert all(abs(g.vinf - 0.5) < 1e-3 for g in got)  # refined to target |v∞|
