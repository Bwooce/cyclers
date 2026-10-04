"""More published controls from Neelakantan & Ramanan (2022): the two circular-problem halo orbits
and the multi-revolution "M4N2 Lyapunov" row.

Source: Neelakantan & Ramanan, "Two-impulse transfer to multi-revolution halo orbits in the
Earth-Moon elliptic restricted three body problem framework", Journal of Astrophysics and Astronomy
43:50 (2022), DOI 10.1007/s12036-022-09830-x. Filed in the private paper corpus as
neelakantan-ramanan-2022-two-impulse-transfer-multi-revolution-halo-orbits-earth-moon-ER3BP-
jaa-43-50-doi-10.1007-s12036-022-09830-x.pdf. Table 8 rows 1-4 are tested in
``test_er3bp_neelakantan_2022.py``; this file does not repeat them.

1. Circular-problem halo orbits (``core/cr3bp.py``). The paper prints two L1 northern halo states
   [x0, 0, z0, 0, ydot0, 0] with periods in days only:
   - p7, section 5.1, the e = 0 benchmark (Az = 15,000 km): [0.8229570002125, 0.0,
     0.04241855682667, 0.0, 0.152044631998602, 0.0], period 12.05848 days;
   - p8, section 5.2.1, the circular comparison orbit for M5N2 (Az = 47,924 km):
     [0.845272317414636, 0.0, 0.170480760011887, 0.0, 0.265160667222108, 0.0], period 11.459994
     days.
   The paper states mu = 0.0122 (p5) for its Earth-Moon work. With that mass ratio neither state
   crosses the x-z plane perpendicularly at its first return ((xdot, zdot) residual norms 0.031
   and 0.0037); that is recorded below as a strict expected failure. Solving for the mass ratio
   at which each printed state does cross perpendicularly gives 0.012277471000011 (benchmark) and
   0.01227747100007 (comparison orbit), independently and to eleven digits the same value. This
   value, mu = 0.012277471, is not printed anywhere in the paper; it is identified here from the
   two states. (It is the Earth-Moon mass ratio commonly used in older textbook examples; that
   attribution is ours, not the paper's.) With it both states are periodic halo orbits in the
   project's model. Control: with the project's Earth-Moon value 0.0121505843 neither closes.

   The printed periods in days cannot be checked directly because the paper prints no time unit.
   Their ratio is unit-free: printed 12.05848 / 11.459994 = 1.0522239; measured at mu = 0.012277471
   2.7533715 / 2.6183733 = 1.0515584, a relative difference of 6.3e-4, well beyond the seven or
   eight printed digits. The implied time units are 4.3795 and 4.3768 days; neither matches a
   standard Earth-Moon unit (4.3425 days from 384,400 km and G(M_E + M_M); 4.3483 days from the
   sidereal month). Recorded as a strict expected failure.

2. Table 8 (p14), "M4N2 Lyapunov": x0 = 0.804504659626012, z0 = 0.00000000000000020,
   ydot0 = 0.31826866733409, period 4 pi (digits re-read from the PDF). With mu = 0.0122 and
   e = 0.0554 it does not close from periapsis (closure 4.51; residual (y, xdot) at the half
   period (-2.28, -1.86)) or from apoapsis (closure 2.0; half-period residual (0.71, 0.067)). From
   periapsis the state leaves the L1 region within about half a revolution of the primaries.
   The sibling row M2N1, an orbit of the same region, has a monodromy multiplier of 2.56e5 per
   2 pi; a 15-digit state on such an orbit misses the half-period (2 pi) crossing by about 1e-10,
   so rounding of the printed digits cannot explain the miss. A Newton differential correction in
   (x0, ydot0) at the half-period symmetry (target y = xdot = 0 at f0 + 2 pi, ``er3bp_stm_eom``)
   started from the printed state converged only to orbits far away: (dx0, dydot0) =
   (+0.094, -0.266) from periapsis and (-0.486, +1.106) from apoapsis. The one nearby
   symmetric orbit found, at (dx0, dydot0) = (-1.8e-3, -1.07e-3) from periapsis, crosses
   perpendicularly at f = pi, so its period is 2 pi (N = 1) and it is not the row's orbit. The
   other mass ratios tried (0.012277471, 0.0121505843) and e = 0.0549 or 0 do not help either. No
   periodic orbit with the row's period was found near the printed state; the direct closure is
   kept as a strict expected failure.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.er3bp as er3bp
from cyclerfinder.core.cr3bp import cr3bp_eom

FloatArray = NDArray[np.float64]

_MU_PRINTED = 0.0122
_MU_IDENTIFIED = 0.012277471  # not printed in the paper; see module docstring
_MU_PROJECT_EM = 0.0121505843
_ECC = 0.0554

# (label, state, printed period in days)
_HALOS = [
    (
        "benchmark Az 15000 km (p7)",
        (0.8229570002125, 0.0, 0.04241855682667, 0.0, 0.152044631998602, 0.0),
        12.05848,
    ),
    (
        "M5N2 circular comparison (p8)",
        (0.845272317414636, 0.0, 0.170480760011887, 0.0, 0.265160667222108, 0.0),
        11.459994,
    ),
]


def _y_zero_event() -> Callable[[float, FloatArray], float]:
    def ev(t: float, y: FloatArray, *_args: object) -> float:
        return float(y[1])

    ev.terminal = True  # type: ignore[attr-defined]
    return ev


def _first_return(state0: FloatArray, mu: float) -> tuple[float, FloatArray]:
    """Time and state at the first x-z plane crossing after leaving the start."""
    kick = solve_ivp(
        cr3bp_eom, (0.0, 1e-3), state0, args=(mu,), method="DOP853", rtol=1e-13, atol=1e-13
    )
    sol = solve_ivp(
        cr3bp_eom,
        (1e-3, 10.0),
        kick.y[:, -1],
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        events=_y_zero_event(),
    )
    assert sol.t_events is not None and sol.y_events is not None
    assert sol.t_events[0].size > 0
    return float(sol.t_events[0][0]), np.asarray(sol.y_events[0][0], dtype=np.float64)


def _closure(state0: FloatArray, mu: float, period: float) -> float:
    sol = solve_ivp(
        cr3bp_eom, (0.0, period), state0, args=(mu,), method="DOP853", rtol=1e-13, atol=1e-13
    )
    return float(np.linalg.norm(sol.y[:, -1] - state0))


@pytest.fixture(scope="module")
def halo_half_periods() -> dict[str, tuple[float, FloatArray]]:
    return {label: _first_return(np.array(state), _MU_IDENTIFIED) for label, state, _ in _HALOS}


@pytest.mark.parametrize(("label", "state", "days"), _HALOS)
def test_printed_halo_is_periodic_at_mu_0_012277471(
    label: str,
    state: tuple[float, ...],
    days: float,
    halo_half_periods: dict[str, tuple[float, FloatArray]],
) -> None:
    """Measured: half-period (xdot, zdot) residuals 3.9e-12 and 3.3e-12; closures after twice the
    half period (2.7533715 and 2.6183733) 1.8e-10 and 2.7e-11."""
    state0 = np.array(state)
    half_t, half = halo_half_periods[label]
    assert math.hypot(half[3], half[5]) < 1e-9, label
    assert _closure(state0, _MU_IDENTIFIED, 2.0 * half_t) < 1e-8, label


@pytest.mark.xfail(
    strict=True,
    reason=(
        "#896: with the paper's stated mu = 0.0122 the printed circular-problem halo states are "
        "not periodic: half-period (xdot, zdot) residual norms 0.031 (benchmark) and 0.0037 "
        "(M5N2 comparison); they close at the unprinted mu = 0.012277471 instead"
    ),
)
@pytest.mark.parametrize(("label", "state", "days"), _HALOS)
def test_printed_halo_is_periodic_at_stated_mu(
    label: str, state: tuple[float, ...], days: float
) -> None:
    _, half = _first_return(np.array(state), _MU_PRINTED)
    assert math.hypot(half[3], half[5]) < 1e-9, label


@pytest.mark.parametrize(("label", "state", "days"), _HALOS)
def test_printed_halo_is_not_periodic_at_project_earth_moon_mu(
    label: str, state: tuple[float, ...], days: float
) -> None:
    """Control: half-period (xdot, zdot) residual norms 0.054 and 0.0061 at mu = 0.0121505843."""
    _, half = _first_return(np.array(state), _MU_PROJECT_EM)
    assert math.hypot(half[3], half[5]) > 1e-3, label


@pytest.mark.xfail(
    strict=True,
    reason=(
        "#896: printed period ratio 12.05848/11.459994 = 1.0522239 but the measured ratio at "
        "mu = 0.012277471 is 2.7533715/2.6183733 = 1.0515584 (relative 6.3e-4); implied time "
        "units 4.3795 and 4.3768 days disagree"
    ),
)
def test_printed_periods_in_days_have_the_measured_ratio(
    halo_half_periods: dict[str, tuple[float, FloatArray]],
) -> None:
    (label_a, _, days_a), (label_b, _, days_b) = _HALOS
    measured = halo_half_periods[label_a][0] / halo_half_periods[label_b][0]
    assert measured / (days_a / days_b) == pytest.approx(1.0, abs=1e-6)


_M4N2_LYAPUNOV = np.array([0.804504659626012, 0.0, 0.0, 0.0, 0.31826866733409, 0.0])


@pytest.mark.xfail(
    strict=True,
    reason=(
        "#896: Table 8 'M4N2 Lyapunov' does not close in the ER3BP (mu 0.0122, e 0.0554): "
        "closure 4.51 from periapsis, 2.0 from apoapsis; half-period (y, xdot) (-2.28, -1.86) and "
        "(0.71, 0.067); a 2.56e5-per-2pi multiplier would give ~1e-10, so not digit rounding; "
        "Newton at the half-period symmetry lands far away or on a 2pi orbit. #925: no 4pi "
        "symmetric orbit within 0.01 of the printed state (bounded solves never close; unbounded "
        "ones reach unrelated 4pi orbits at x0 0.769 and 0.694), while the sibling M2N1 row "
        "re-corrects in place and the M5N2 row and its printed multipliers reproduce "
        "(test_er3bp_peng_xu_2015.py): a slip or mislabel in the paper's row, not the model"
    ),
)
@pytest.mark.parametrize("f0", [0.0, math.pi], ids=["periapsis", "apoapsis"])
def test_table_8_m4n2_lyapunov_closes(f0: float) -> None:
    period = 4.0 * math.pi
    sol = solve_ivp(
        er3bp.er3bp_eom,
        (f0, f0 + period),
        _M4N2_LYAPUNOV,
        args=(_MU_PRINTED, _ECC),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    assert sol.success
    assert sol.sol is not None
    half = sol.sol(f0 + 0.5 * period)
    assert math.hypot(half[1], half[3]) < 1e-6
    assert float(np.linalg.norm(sol.y[:, -1] - _M4N2_LYAPUNOV)) < 1e-6
