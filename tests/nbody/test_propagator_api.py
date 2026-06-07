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
