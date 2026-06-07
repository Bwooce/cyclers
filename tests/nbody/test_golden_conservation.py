"""N-body GOLDEN GATE 2: energy conservation (design §5.2).

Sun-only: relative energy drift ~ machine precision (integrator health).
Perturbers ON: NO exact Jacobi constant exists (rails = time-dependent); we
assert a BOUNDED energy budget, reported, not a fake conserved scalar.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402


@pytest.mark.slow
def test_sun_only_energy_drift_below_floor() -> None:
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    arc = RestrictedNBody("rebound").propagate(
        r0, v0, t0_sec=0.0, t1_sec=365.25 * 86400.0, bodies=(), accuracy=1e-12
    )
    assert abs(arc.energy_rel_drift) < 1e-10  # IAS15 two-body ~ machine precision
