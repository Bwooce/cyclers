"""N-body Phase A: propagator interface + rails-force shape (plan Phase A)."""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.nbody.forces import rails_acceleration


def test_sun_only_acceleration_is_central_inverse_square() -> None:
    # No perturbers -> pure -mu_sun r_hat / r^2.
    r = np.array([1.5e8, 0.0, 0.0])  # ~1 AU on +x
    a = rails_acceleration(r, t_sec=0.0, bodies=(), ephem=None)
    expected = -MU_SUN_KM3_S2 / (1.5e8**2)
    assert a[0] == pytest.approx(expected, rel=1e-12)
    assert a[1] == 0.0 and a[2] == 0.0


def test_propagator_backend_selectable() -> None:
    pytest.importorskip("rebound")
    from cyclerfinder.nbody.propagator import RestrictedNBody

    prop = RestrictedNBody("rebound")
    assert prop.backend == "rebound"
    with pytest.raises(ValueError, match="backend"):
        RestrictedNBody("tudat")  # deferred slot: not wired yet (design Q1)


@pytest.mark.slow
def test_anchor_err_km_is_reserved_zero_contract() -> None:
    """RestrictedNBody seeds at the exact (r0, v0) and never anchors against an
    ephemeris, so anchor_err_km is a RESERVED slot that is always 0.0. Callers
    that need a closure error compute it themselves (rung.py does). Pin the
    contract so a future backend cannot silently start populating it unannounced.
    """
    pytest.importorskip("rebound")
    from cyclerfinder.nbody.propagator import RestrictedNBody

    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    arc = RestrictedNBody("rebound").propagate(
        r0, v0, t0_sec=0.0, t1_sec=10.0 * 86400.0, bodies=(), accuracy=1e-10
    )
    assert arc.anchor_err_km == 0.0
