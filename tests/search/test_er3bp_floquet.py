from __future__ import annotations

import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.er3bp_floquet import er3bp_monodromy, floquet_classify


def test_monodromy_is_6x6_real() -> None:
    sys = ER3BPSystem(mu=0.0121550, e=0.0001, primary_name="E", secondary_name="M")
    state0 = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
    m = er3bp_monodromy(state0, period_f=2.0 * np.pi, system=sys)
    assert m.shape == (6, 6)
    assert np.isfinite(m).all()


def test_floquet_classify_on_broucke_orbit() -> None:
    """Broucke Orbit-1 monodromy classifies to a definite tag; symplectic
    eigenvalues have magnitude-product ~1."""
    sys = ER3BPSystem(mu=0.0121550, e=0.0001, primary_name="E", secondary_name="M")
    state0 = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
    m = er3bp_monodromy(state0, period_f=2.0 * np.pi, system=sys)
    res = floquet_classify(m)
    assert res.stability_tag in {"stable", "unstable", "marginal"}
    assert np.isclose(float(np.prod(np.abs(res.eigenvalues))), 1.0, atol=1e-2)
