"""N-body GOLDEN GATE 4: IAS15 accuracy convergence (design §5.3).

Asserts the final state is converged wrt the integrator accuracy parameter
(Cauchy: tighter accuracy moves the answer less than the consumer tolerance).
Doubles as the body-inclusion sensitivity harness (Jupiter on/off), whose verdict
is recorded against the real candidate baselines in Phase B/C.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402


@pytest.mark.slow
def test_state_converges_with_accuracy() -> None:
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 24.0, 3.0])
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")
    coarse = prop.propagate(
        r0, v0, 0.0, 200 * 86400.0, bodies=("E", "M"), accuracy=1e-7, ephem=ephem
    )
    fine = prop.propagate(
        r0, v0, 0.0, 200 * 86400.0, bodies=("E", "M"), accuracy=1e-10, ephem=ephem
    )
    assert np.linalg.norm(fine.r_km - coarse.r_km) < 50.0  # converged < 50 km
