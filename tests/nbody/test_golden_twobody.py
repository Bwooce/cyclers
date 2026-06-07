"""N-body GOLDEN GATE 1: Sun-only n-body == core/kepler.propagate (design §5.1).

SOURCE-FREE cross-implementation golden: REBOUND/IAS15 and core/kepler must both
reduce to the SAME analytic two-body solution. A disagreement is a §0 frame/time/
units slip in the glue, not a physics choice. Neither side is EXPECTED; the
analytic Kepler problem is.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core import kepler  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402


@pytest.mark.slow
def test_sun_only_matches_kepler_propagate() -> None:
    r0 = np.array([1.496e8, 0.0, 0.0])  # ~1 AU
    v0 = np.array([0.0, 29.78, 0.0])  # ~circular
    dt = 120.0 * 86400.0  # 120 days
    kep_r, kep_v = kepler.propagate(r0, v0, dt)
    arc = RestrictedNBody("rebound").propagate(
        r0, v0, t0_sec=0.0, t1_sec=dt, bodies=(), accuracy=1e-12
    )
    assert np.linalg.norm(arc.r_km - kep_r) < 1.0  # < 1 km
    assert np.linalg.norm(arc.v_km_s - kep_v) < 1e-5
