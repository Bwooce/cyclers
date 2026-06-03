"""Tests for :func:`cyclerfinder.core.kepler.coe_to_rv`.

These are geometry-contract tests: the EXPECTED values are COMPUTED from the
defining formulae of the perifocal->inertial conversion (vis-viva, the radius
equation, and the orthogonality of position and velocity for a circular orbit),
not taken from any external source. They assert that the implementation honours
its own mathematical contract.
"""

from __future__ import annotations

from math import sqrt

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.kepler import coe_to_rv


def test_coe_to_rv_circular() -> None:
    """Circular orbit (e=0) at 1 AU: |r|=AU, |v|=sqrt(mu/r), r.v=0.

    EXPECTED COMPUTED — geometry contract. For a circular orbit the radius
    equals ``a``, the speed equals the circular velocity ``sqrt(mu/r)``, and
    position is perpendicular to velocity.
    """
    r, v = coe_to_rv(a_km=AU_KM, e=0.0, true_anom_rad=0.0, mu=MU_SUN_KM3_S2)
    assert abs(float(np.linalg.norm(r)) - AU_KM) < 1.0
    assert abs(float(np.linalg.norm(v)) - sqrt(MU_SUN_KM3_S2 / AU_KM)) < 1e-6
    assert abs(float(np.dot(r, v))) < 1e-3


def test_coe_to_rv_perihelion() -> None:
    """Elliptic orbit at perihelion (nu=0): |r| = a(1-e).

    EXPECTED COMPUTED — geometry contract (the radius equation at nu=0).
    """
    a = 1.30 * AU_KM
    e = 0.257
    r, v = coe_to_rv(a_km=a, e=e, true_anom_rad=0.0, mu=MU_SUN_KM3_S2)
    assert abs(float(np.linalg.norm(r)) - a * (1.0 - e)) < 1.0
    # At perihelion the speed equals vis-viva: sqrt(mu (2/r - 1/a)).
    r_peri = a * (1.0 - e)
    v_peri = sqrt(MU_SUN_KM3_S2 * (2.0 / r_peri - 1.0 / a))
    assert abs(float(np.linalg.norm(v)) - v_peri) < 1e-6
    # Position and velocity are orthogonal at periapsis.
    assert abs(float(np.dot(r, v))) < 1e-3


def test_coe_to_rv_arg_peri_rotation() -> None:
    """Argument of periapsis rotates the state about +z without changing norms.

    EXPECTED COMPUTED — geometry contract (a rigid rotation preserves lengths).
    """
    a = 1.30 * AU_KM
    e = 0.257
    nu = 0.7
    r0, v0 = coe_to_rv(a_km=a, e=e, true_anom_rad=nu, arg_peri_rad=0.0)
    r1, v1 = coe_to_rv(a_km=a, e=e, true_anom_rad=nu, arg_peri_rad=1.234)
    assert abs(float(np.linalg.norm(r0)) - float(np.linalg.norm(r1))) < 1e-6
    assert abs(float(np.linalg.norm(v0)) - float(np.linalg.norm(v1))) < 1e-9
    assert r0[2] == 0.0
    assert r1[2] == 0.0
