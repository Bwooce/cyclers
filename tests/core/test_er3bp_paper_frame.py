"""ER3BP paper-frame (Antoniadou & Libert 2018 Eq. 1) EOM + corrector tests.

Promotes the #442 capability deliverables: the validated Eq. 1 EOM in the
paper's non-pulsating rotating frame, and the 4-vector full-period corrector
that tracks the CONNECTED 3/1 (pi,0) resonant family from e2=0 to e2=0.90.

Reference: Antoniadou, K.I. & Libert, A.-S. (2018), arXiv:1805.00288.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom
from cyclerfinder.core.er3bp_paper_frame import (
    correct_resonant_member,
    osculating_a1,
    paper_frame_eom,
)

MU = 0.001
THETA0 = math.pi
PERIOD = 2.0 * math.pi

# Published Antoniadou & Libert 2018 Fig. 3 DS-map value for the 3/1 resonant
# family's osculating semi-major-axis ratio a1/a2 (= a1 since a2 = 1). SOURCED:
# Fig. 3 caption / Fig. 11(e) header (a1 = 0.480674). This is the EXPECTED side
# of the golden test and must trace to the paper, never to our own code.
PUBLISHED_A1 = 0.4807

# Converged CONNECTED 3/1 (pi,0) family members (#442 verdict doc UPDATE (c)),
# IC = [x, y=0, vx=0, vy, theta0=pi], T=2pi. Used as seeds + regression anchors.
LADDER: dict[float, tuple[float, float]] = {
    0.00: (0.4793554969, 0.9624369118),
    0.30: (0.6909309271, 0.5087470978),
    0.55: (0.8217405562, 0.3046047985),
    0.85: (0.9457500071, 0.0255730480),
    0.90: (0.9341012757, -0.3487825823),
}


def test_paper_frame_reduces_to_cr3bp_at_e2_zero() -> None:
    """VALIDATION: at e2=0 the paper-frame EOM reproduces cr3bp_eom over a period.

    OUR-vs-OUR parity: integrate an arbitrary bounded planar IC under both the
    e2=0 paper-frame EOM and the planar slice of cr3bp_eom; the trajectories must
    agree to < 1e-9.
    """
    state4 = [0.85, 0.0, 0.1, 0.42]  # x, y, vx, vy
    t_final = 6.0  # ~ one synodic period

    s0_paper = np.array([state4[0], state4[1], state4[2], state4[3], 0.0])
    sol_p = solve_ivp(
        paper_frame_eom,
        (0.0, t_final),
        s0_paper,
        args=(MU, 0.0),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol_p.success

    s0_cr = np.array([state4[0], state4[1], 0.0, state4[2], state4[3], 0.0])
    sol_c = solve_ivp(
        cr3bp_eom,
        (0.0, t_final),
        s0_cr,
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol_c.success
    assert sol_p.sol is not None and sol_c.sol is not None

    tg = np.linspace(0.0, t_final, 400)
    yp = sol_p.sol(tg)[:4]  # x, y, vx, vy
    yc = sol_c.sol(tg)[[0, 1, 3, 4]]  # x, y, vx, vy (drop z, vz)
    assert float(np.max(np.abs(yp - yc))) < 1e-9


def test_connected_member_regression_anchor_e2_zero() -> None:
    """REGRESSION: the e2=0 connected member matches the #442 verdict-doc IC.

    Converge from the published seed and assert the IC + independent closure
    residual reproduce verdict-doc UPDATE (c): x=0.4793554969, vy=0.9624369118,
    independent residual < 1e-8.
    """
    x0, vy0 = LADDER[0.00]
    res = correct_resonant_member(x0, vy0, theta0=THETA0, mu=MU, e2=0.0)
    assert float(res["residual"]) < 1e-8
    assert abs(float(res["x"]) - 0.4793554969) < 1e-7
    assert abs(float(res["vy"]) - 0.9624369118) < 1e-7


@pytest.mark.parametrize("e2", list(LADDER.keys()))
def test_connected_family_osculating_a1_matches_published(e2: float) -> None:
    """GOLDEN: connected 3/1 (pi,0) family osc a1 matches the published 0.4807.

    Converge each ladder member in the paper frame and assert its osculating
    semi-major axis lands at the published Antoniadou & Libert 2018 Fig. 3
    DS-map value (a1 = 0.4807). The converged family sits at osc a1 ~ 0.4806
    throughout (#442 verdict), so a 2e-3 tolerance is comfortably met while
    still being a meaningful match to the sourced value.
    """
    x0, vy0 = LADDER[e2]
    res = correct_resonant_member(x0, vy0, theta0=THETA0, mu=MU, e2=e2)
    assert float(res["residual"]) < 1e-8  # member genuinely closes
    a1 = osculating_a1(float(res["x"]), float(res["vy"]), THETA0, MU, e2)
    assert abs(a1 - PUBLISHED_A1) < 2e-3, f"e2={e2}: osc a1={a1:.6f} vs published {PUBLISHED_A1}"
