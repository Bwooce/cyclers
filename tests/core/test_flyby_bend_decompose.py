"""M-3D Phase 3: flyby bend in/out-of-plane attribution (plan §3). DIAGNOSTIC
ONLY. Physics invariants (no sourced V_inf anchor): a purely in-plane bend has
zero out-of-plane component, and the total recovers bend_angle.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.flyby import bend_angle, bend_decompose

_NZ = np.array([0.0, 0.0, 1.0])  # ecliptic normal


def test_inplane_bend_has_zero_outofplane() -> None:
    vin = np.array([5.0, 0.0, 0.0])
    vout = np.array([5.0 * np.cos(0.3), 5.0 * np.sin(0.3), 0.0])
    d_in, d_out = bend_decompose(vin, vout, _NZ)
    assert d_in == pytest.approx(0.3, abs=1e-9)
    assert abs(d_out) < 1e-9


def test_pure_outofplane_bend_has_zero_inplane() -> None:
    vin = np.array([5.0, 0.0, 0.0])
    vout = np.array([5.0 * np.cos(0.2), 0.0, 5.0 * np.sin(0.2)])
    d_in, d_out = bend_decompose(vin, vout, _NZ)
    assert abs(d_in) < 1e-9
    assert abs(d_out) == pytest.approx(0.2, abs=1e-9)


def test_decomposition_components_consistent_with_total_bend() -> None:
    vin = np.array([4.0, 1.0, 0.5])
    vout = np.array([3.5, 1.8, -0.4])
    d_in, d_out = bend_decompose(vin, vout, _NZ)
    total = bend_angle(vin, vout)
    # In/out-of-plane are orthogonal contributions; quadrature recovers the
    # total to small-angle order (consistency invariant, not a sourced anchor).
    assert np.hypot(d_in, d_out) == pytest.approx(total, rel=0.15)


def test_zero_vector_raises() -> None:
    with pytest.raises(ValueError):
        bend_decompose(np.zeros(3), np.array([1.0, 0.0, 0.0]), _NZ)
