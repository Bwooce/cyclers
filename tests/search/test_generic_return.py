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
