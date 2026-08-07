"""Unit tests for cyclerfinder.core.ccr4bp_titan_hyperion."""

from __future__ import annotations

import math

from cyclerfinder.core.ccr4bp_titan_hyperion import (
    L_KM,
    saturn_titan_hyperion_default,
    v_unit_km_s,
)


def test_saturn_titan_hyperion_default() -> None:
    sys = saturn_titan_hyperion_default()
    assert sys.mu > 0.0
    assert sys.mu_gan > 0.0
    assert sys.a_gan > 1.0
    assert sys.omega_gan < 0.0  # Synodic rate in rotating frame is negative for outer moon
    assert math.isclose(L_KM, 1221870.0, rel_tol=1e-5)
    v_unit = v_unit_km_s()
    assert 5.0 < v_unit < 10.0
