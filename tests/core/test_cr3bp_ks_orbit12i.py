"""#928 stage 4 / #929 (b): the closest passage of the #896 Font-Nunes-Simo orbit 12i, by plain
DOP853, the Sundman propagator and the KS propagator.

Source of the orbit: J. Font, A. Nunes and C. Simo, "A numerical study of the orbits of second
species of the planar circular RTBP", Celest. Mech. Dyn. Astron. 103, 143-162 (2009), Table 4,
p.155, Fig. 12i: mu = 1e-4, C_J = 2.8 (with the constant mu(1 - mu)), encounter circle radius
mu^(2/5), (phi, psi) as printed below. Frame and angle conventions as read in
``tests/core/test_cr3bp_font_nunes_simo_second_species.py`` (paper frame = project frame turned by
pi; phi the position angle on the circle about the small primary, psi the velocity angle; the
printed point is an exit from the disk). The papers used Levi-Civita regularisation.

The passage that ends at the printed point (entered about 0.104 time units earlier) comes within
8.72e-6 of the small primary; 8.7e-6 is the project's #896 measurement, not a printed number.
Measured on 2026-10-05 for that passage (entry state from the KS run at rtol 1e-13, then forward
to the exit; angles in rad against the printed exit point):

==========  =======  ==========  ==============  ======
integrator  rtol     |d angle|   |Jacobi error|  nfev
==========  =======  ==========  ==============  ======
KS          1e-12    3.4e-13     8.9e-16         600
KS          1e-13    1.2e-13     4.4e-16         696
plain       1e-12    2.7e-10     4.1e-10         3257
plain       1e-13    2.0e-10     3.0e-10         3473
plain       2.3e-14  2.1e-10     3.3e-10         3797
Sundman r2  1e-12    5.7e-11     1.5e-10         1229
Sundman r2  1e-13    3.0e-10     4.5e-10         1601
Sundman r2  2.3e-14  2.3e-10     3.1e-10         1865
==========  =======  ==========  ==============  ======

(plain and Sundman with atol = rtol, as in the #896 test; 2.3e-14 is the smallest rtol DOP853
accepts.) Backward from the printed point to the entry at rtol 1e-13 the Sundman run's Jacobi error
is 1.18e-9 and its angle error 8.3e-10 against KS: the #896 passage tolerance of 1e-9
(``_passage_tolerance``) holds for the angles but not for the Jacobi constant in that direction.
The KS run conserves C to rounding, returns to the printed point to 6e-14 after the round trip, and
uses a fifth of the plain run's evaluations.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom, jacobi_constant
from cyclerfinder.core.cr3bp_ks import (
    KSModel,
    MoonCentredCR3BP,
    integrate_ks,
    ks_state_from_physical,
    physical_from_ks_state,
    propagate_ks,
)
from cyclerfinder.core.cr3bp_regularized import sundman_rhs

FloatArray = NDArray[np.float64]

MU = 1.0e-4
C_J = 2.8
RADIUS: float = MU**0.4
X_M = 1.0 - MU
C_PROJECT = C_J - MU * (1.0 - MU)
PHI_12I = "3.27824276059703769776947746511"  # 2009 Table 4, p.155, as printed
PSI_12I = "3.21924695055804868387324783072"


def _state(phi: float, psi: float) -> FloatArray:
    """Project-frame state on the circle for printed (phi, psi) (paper frame turned by pi)."""
    pos, vel = phi + math.pi, psi + math.pi
    s = np.array([X_M + RADIUS * math.cos(pos), RADIUS * math.sin(pos), 0.0, 0.0, 0.0, 0.0])
    speed = math.sqrt(jacobi_constant(s, MU) - C_PROJECT)
    s[3], s[4] = speed * math.cos(vel), speed * math.sin(vel)
    return s


def _angles(s: FloatArray) -> FloatArray:
    return np.array([math.atan2(-s[1], -(s[0] - X_M)), math.atan2(-s[4], -s[3])])


def _ks_entry(rtol: float) -> tuple[FloatArray, float, float]:
    """KS backward from the printed exit to the entry crossing: (state, time, closest r)."""
    m = MoonCentredCR3BP(MU)
    y0 = ks_state_from_physical(m, _state(float(PHI_12I), float(PSI_12I)))

    def circle(s: float, y: FloatArray, model: KSModel) -> float:
        return float(y[:4] @ y[:4]) - RADIUS

    circle.terminal = True  # type: ignore[attr-defined]
    circle.direction = 1.0  # type: ignore[attr-defined]  # backward in s, r grows: - to +

    def extremum(s: float, y: FloatArray, model: KSModel) -> float:
        return float(y[:4] @ y[4:8])

    sol = integrate_ks(m, y0, (0.0, -1e3), rtol=rtol, atol=1e-16, events=[circle, extremum])
    assert sol.y_events is not None
    ye = sol.y_events[0][0]
    r_min = min(float(y[:4] @ y[:4]) for y in sol.y_events[1])
    return physical_from_ks_state(m, ye), float(ye[9]), r_min


def test_ks_reaches_the_closest_passage_of_12i() -> None:
    """Closest approach 8.7239e-6 (#896 measured 8.7e-6); C conserved to rounding (measured
    4.4e-16 to 8.9e-16); entry agrees between rtol 1e-12 and 1e-13 (measured 3e-15 in time);
    the round trip returns to the printed point (measured 5.6e-14)."""
    s0 = _state(float(PHI_12I), float(PSI_12I))
    entry, t_entry, r_min = _ks_entry(1e-13)
    assert 8.70e-6 < r_min < 8.75e-6
    assert -0.105 < t_entry < -0.103
    assert abs(math.hypot(entry[0] - X_M, entry[1]) - RADIUS) < 1e-15
    assert abs(jacobi_constant(entry, MU) - C_PROJECT) < 1e-14
    entry12, t12, _ = _ks_entry(1e-12)
    assert np.abs(_angles(entry12) - _angles(entry)).max() < 1e-11
    assert abs(t12 - t_entry) < 1e-11
    arc = propagate_ks(MoonCentredCR3BP(MU), entry, -t_entry, rtol=1e-13, atol=1e-16)
    assert np.abs(arc.state - s0).max() < 1e-12
    assert arc.r_min == pytest.approx(r_min, rel=1e-9)


def _exit_event(t: float, y: FloatArray, *_args: object) -> float:
    return float(math.hypot(y[0] - X_M, y[1]) - RADIUS)


_exit_event.terminal = True  # type: ignore[attr-defined]
_exit_event.direction = 1.0  # type: ignore[attr-defined]


@pytest.mark.parametrize("rtol", [1e-12, 1e-13, 2.3e-14])
def test_three_way_comparison_through_the_passage(rtol: float) -> None:
    """From the KS entry state forward to the exit: KS hits the printed exit to 1e-12 rad with
    C exact to 1e-14; plain DOP853 and the Sundman (r2) form land within the #896 passage
    tolerance 1e-9 rad of it, with Jacobi errors of 1.5e-10 to 4.5e-10 (bound 1e-9); KS takes
    under a third of the plain run's evaluations. Table of measurements in the module
    docstring."""
    s_exit = _state(float(PHI_12I), float(PSI_12I))
    entry, t_entry, _ = _ks_entry(1e-13)
    ks = propagate_ks(MoonCentredCR3BP(MU), entry, -t_entry, rtol=max(rtol, 1e-13), atol=1e-16)
    assert np.abs(_angles(ks.state) - _angles(s_exit)).max() < 1e-12
    assert abs(jacobi_constant(ks.state, MU) - C_PROJECT) < 1e-14

    plain = solve_ivp(
        cr3bp_eom,
        (0, 1),
        entry,
        method="DOP853",
        rtol=rtol,
        atol=rtol,
        args=(MU,),
        events=_exit_event,
    )
    assert plain.y_events is not None and plain.t_events is not None
    p_exit = plain.y_events[0][0]
    assert np.abs(_angles(p_exit) - _angles(s_exit)).max() < 1e-9
    assert abs(float(plain.t_events[0][0]) + t_entry) < 1e-9
    assert abs(jacobi_constant(p_exit, MU) - C_PROJECT) < 1e-9

    sund = solve_ivp(
        sundman_rhs,
        (0, 1e5),
        np.append(entry, 0.0),
        method="DOP853",
        rtol=rtol,
        atol=rtol,
        args=(MU, "r2"),
        events=_exit_event,
    )
    assert sund.y_events is not None
    s_end = sund.y_events[0][0]
    assert np.abs(_angles(s_end[:6]) - _angles(s_exit)).max() < 1e-9
    assert abs(s_end[6] + t_entry) < 1e-9
    assert abs(jacobi_constant(s_end[:6], MU) - C_PROJECT) < 1e-9

    assert 3 * ks.nfev < plain.nfev
