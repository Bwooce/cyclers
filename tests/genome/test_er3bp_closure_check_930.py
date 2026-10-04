"""#930: the corrector's independent full-period closure check must fail loudly.

Controls (expected values are physics or printed numbers, never values this code computed):

* mu = 0: the elliptic problem is Kepler motion about the primary at the origin (Gomez & Olle 1986,
  digest 2026-10-04 section on the mu = 0 limit). A test particle on an a = 1 Kepler ellipse (the
  primaries' mean motion is 1) closes in the pulsating frame after f = 2 pi, whatever the particle
  eccentricity. Started at the periapsis of both orbits it is symmetric (y = 0, x' = 0 at f = 0 and
  at f = pi), so the symmetric half-period corrector applies. With a half period that is not pi
  Newton cannot converge here (measured: no symmetric crossing exists for 1.05 pi), so the wrong
  period is exercised through a loose Newton tolerance, which accepts the unchanged state, and the
  independent check must then reject it. A state perturbed by 1e-4 and accepted at tol = 1e-3
  reproduces the Gomez & Olle 1991 situation (accepted by the crossing residual, 4e-5 and 2e-2
  away from closing).
* Peng & Xu 2015 (ASR 55:1015) M5N2 halo, mu = 0.0122, e = 0.0554, 15-digit row of Neelakantan &
  Ramanan 2022 Table 8, period 4 pi: a known printed elliptic orbit must be accepted.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_periodic import (
    ClosureError,
    ConvergenceError,
    correct_er3bp_periodic,
)

_E_PRIMARIES = 0.0554
_E_PARTICLE = 0.3


def _kepler_pulsating_state(e_p: float, e_part: float) -> np.ndarray:
    """Pulsating-frame state at f = 0 of an a = 1 Kepler orbit (mu = 0), both at periapsis.

    Primaries: r(0) = 1 - e_p, df/dt = sqrt(p) / r^2 with p = 1 - e_p^2 (total mass 1). Pulsating
    coordinates are the inertial position rotated by -f and divided by r(f); at f = 0 with r' = 0:
    x = X / r, x' = V_x / (r df/dt) + Y / r, y' = V_y / (r df/dt) - X / r.
    """
    r0 = 1.0 - e_p
    fdot = math.sqrt(1.0 - e_p * e_p) / (r0 * r0)
    xp = 1.0 - e_part
    vy = math.sqrt((1.0 + e_part) / (1.0 - e_part))
    return np.array([xp / r0, 0.0, 0.0, 0.0, vy / (r0 * fdot) - xp / r0, 0.0])


def test_closure_error_is_a_convergence_error() -> None:
    assert issubclass(ClosureError, ConvergenceError)


def test_mu0_kepler_right_period_passes() -> None:
    system = ER3BPSystem(0.0, _E_PRIMARIES, "P0", "P1")
    s0 = _kepler_pulsating_state(_E_PRIMARIES, _E_PARTICLE)
    orbit = correct_er3bp_periodic(system, s0, math.pi, tol=1e-9)
    assert orbit.independent_residual < 1e-6
    assert float(np.abs(orbit.state0 - s0).max()) < 1e-7
    assert orbit.period_f == pytest.approx(2.0 * math.pi)


def test_mu0_kepler_wrong_period_is_rejected() -> None:
    system = ER3BPSystem(0.0, _E_PRIMARIES, "P0", "P1")
    s0 = _kepler_pulsating_state(_E_PRIMARIES, _E_PARTICLE)
    with pytest.raises(ClosureError) as info:
        correct_er3bp_periodic(system, s0, 1.05 * math.pi, tol=1.0)
    assert info.value.independent_error > 1e-2
    assert info.value.orbit.independent_residual == info.value.independent_error


def test_mu0_kepler_small_state_error_accepted_by_loose_tolerance_is_rejected() -> None:
    system = ER3BPSystem(0.0, _E_PRIMARIES, "P0", "P1")
    s0 = _kepler_pulsating_state(_E_PRIMARIES, _E_PARTICLE)
    s_bad = s0 + np.array([1e-4, 0.0, 0.0, 0.0, 1e-4, 0.0])
    with pytest.raises(ClosureError) as info:
        correct_er3bp_periodic(system, s_bad, math.pi, tol=1e-3)
    assert info.value.independent_error > 1e-5


def test_mu0_wrong_period_is_returned_only_when_the_check_is_disabled() -> None:
    system = ER3BPSystem(0.0, _E_PRIMARIES, "P0", "P1")
    s0 = _kepler_pulsating_state(_E_PRIMARIES, _E_PARTICLE)
    orbit = correct_er3bp_periodic(system, s0, 1.05 * math.pi, tol=1.0, independent_tol=None)
    assert orbit.independent_residual > 1e-2


def test_printed_m5n2_halo_is_accepted() -> None:
    system = ER3BPSystem(0.0122, 0.0554, "Earth", "Moon")
    nr = np.array([0.851666641652152, 0.0, 0.183285539178136, 0.0, 0.25828972225268, 0.0])
    orbit = correct_er3bp_periodic(
        system,
        nr,
        2.0 * math.pi,
        free_vars=(0, 2, 4),
        residual_indices=(1, 3, 5),
        tol=1e-9,
    )
    assert float(np.abs(orbit.state0 - nr).max()) < 1e-7
    assert orbit.independent_residual < 1e-6  # measured 1.1e-8, default independent_tol 1e-5
