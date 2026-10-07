"""#968/#1004: GM and radius overrides on ``jovian_stm.propagate_with_stm``.

The default path must be unchanged, a GM of zero must reduce to the Jupiter-only (Kepler)
propagation, and a scaled GM must change the arc (the knob is live).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.nbody.jovian_ideal import IdealJovianEphemeris
from cyclerfinder.nbody.jovian_stm import propagate_with_stm

A_G = SATELLITES["Ganymede"].sma_km
EPH = IdealJovianEphemeris({"Ganymede": A_G}, {"Ganymede": 0.0}, {"Ganymede": 0.0})
# Start 20,000 km behind Ganymede on its orbit, slightly faster: a close pass within a day.
R0 = np.array([A_G * np.cos(-0.02), A_G * np.sin(-0.02), 0.0])
V_G = np.sqrt(1.26686534e8 / A_G)
V0 = np.array([-V_G * np.sin(-0.02), V_G * np.cos(-0.02), 0.0]) * 1.01
T1 = 2.0 * 86400.0


def test_default_equals_registry_override() -> None:
    a = propagate_with_stm(R0, V0, 0.0, T1, ephem=EPH, moons=("Ganymede",))
    b = propagate_with_stm(
        R0,
        V0,
        0.0,
        T1,
        ephem=EPH,
        moons=("Ganymede",),
        mu_overrides={"Ganymede": SATELLITES["Ganymede"].mu_km3_s2},
        radius_overrides={"Ganymede": SATELLITES["Ganymede"].radius_eq_km},
    )
    for x, y in zip(a, b, strict=True):
        np.testing.assert_array_equal(x, y)


def test_zero_gm_is_jupiter_only() -> None:
    a = propagate_with_stm(
        R0, V0, 0.0, T1, ephem=EPH, moons=("Ganymede",), mu_overrides={"Ganymede": 0.0}
    )
    b = propagate_with_stm(R0, V0, 0.0, T1, ephem=EPH, moons=())
    np.testing.assert_allclose(a[0], b[0], atol=1e-6)
    np.testing.assert_allclose(a[1], b[1], atol=1e-11)
    np.testing.assert_allclose(a[2], b[2], atol=1e-8)


def test_scaled_gm_is_live_and_linear_when_weak() -> None:
    """Small GM scales perturb the Kepler arc linearly (the knob reaches the force)."""
    mu = SATELLITES["Ganymede"].mu_km3_s2
    kep = propagate_with_stm(R0, V0, 0.0, T1, ephem=EPH, moons=())
    d = []
    for s in (1e-5, 2e-5):
        arc = propagate_with_stm(
            R0, V0, 0.0, T1, ephem=EPH, moons=("Ganymede",), mu_overrides={"Ganymede": s * mu}
        )
        d.append(arc[0] - kep[0])
    assert float(np.linalg.norm(d[0])) > 1e-3  # km
    np.testing.assert_allclose(d[1], 2.0 * d[0], rtol=0.02, atol=1e-3)
