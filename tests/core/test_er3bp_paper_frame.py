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
    correct_doubly_symmetric_member,
    correct_resonant_member,
    osculating_a1,
    osculating_e1,
    paper_frame_eom,
    primary_kepler,
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


# Published Antoniadou & Libert 2018 3/1 ISOLATED stable family (I_c)
# representative, Fig. 11(e) DS-map header (config (theta3,theta1)=(pi,0)):
# a1/a2 = 0.480674, e1 = 0.659951 (display-rounded at the nominal e2 ~ 0.90).
# SOURCED: arXiv:1805.00288 Fig. 11(e) header + Sec. 5.3. These are the EXPECTED
# side of the golden test and trace to the paper, never to our own code.
PUBLISHED_IC_A1 = 0.480674
PUBLISHED_IC_E1 = 0.659951


def test_isolated_ic_stable_member_reached_and_closed() -> None:
    """GOLDEN (#457): the 3/1 ISOLATED stable family (I_c) representative is
    reached + closed by the joint doubly-symmetric corrector.

    Seed the published (pi,0) config (P1 at pericentre, theta0=pi) at e2=0.91 and
    converge with :func:`correct_doubly_symmetric_member` (the joint residual
    ``[y(T/2),vx(T/2),y(T),vx(T)]`` that selects the true doubly-symmetric member
    where the plain 2-var full-period objective floors ~2e-6). Assert:
      * the member CLOSES (independent full-period residual < 1e-8);
      * it is genuinely DOUBLY-symmetric (perpendicular at t=0 and t=T/2);
      * its osculating (a1, e1) match the published I_c representative
        (a1=0.480674, e1=0.659951) -- the exact member sits at e2=0.91 vs the
        header's nominal e2~0.90, so a small per-member offset is expected and a
        4e-3 tolerance is a meaningful match to the sourced values;
      * it is STABLE (monodromy eigenvalues on the unit circle), matching the
        published stable I_c segment.
    """
    e2 = 0.91
    half_period = math.pi
    # Build the published (pi,0) P1-pericentre perpendicular-crossing seed.
    r_prim, _rdot, thetadot, _thetaddot = primary_kepler(THETA0, e2)
    x_star = -MU * r_prim
    r1 = PUBLISHED_IC_A1 * (1.0 - PUBLISHED_IC_E1)
    v1 = math.sqrt((1.0 - MU) * (2.0 / r1 - 1.0 / PUBLISHED_IC_A1))
    x_seed = x_star + r1
    vy_seed = v1 - thetadot * x_seed

    res = correct_doubly_symmetric_member(
        x_seed, vy_seed, theta0=THETA0, mu=MU, e2=e2, half_period=half_period
    )
    x, vy = float(res["x"]), float(res["vy"])

    # (1) member genuinely closes over the full period.
    assert float(res["residual"]) < 1e-8, f"independent residual {res['residual']}"

    # (2) doubly-symmetric: perpendicular crossing at t=0 (IC) and at t=T/2.
    sol = solve_ivp(
        paper_frame_eom,
        (0.0, half_period),
        np.array([x, 0.0, 0.0, vy, THETA0]),
        args=(MU, e2),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    assert sol.success
    y_half, vx_half = float(sol.y[1, -1]), float(sol.y[2, -1])
    assert abs(y_half) < 1e-6 and abs(vx_half) < 1e-6, (
        f"t=T/2 not perpendicular: y={y_half}, vx={vx_half}"
    )

    # (3) osculating (a1, e1) match the published representative.
    a1 = osculating_a1(x, vy, THETA0, MU, e2)
    e1 = osculating_e1(x, vy, THETA0, MU, e2)
    assert abs(a1 - PUBLISHED_IC_A1) < 4e-3, f"a1={a1:.6f} vs published {PUBLISHED_IC_A1}"
    assert abs(e1 - PUBLISHED_IC_E1) < 4e-3, f"e1={e1:.6f} vs published {PUBLISHED_IC_E1}"

    # (4) STABLE: full-period monodromy eigenvalues on the unit circle.
    period = 2.0 * half_period

    def _prop(s0: np.ndarray) -> np.ndarray:
        out = solve_ivp(
            paper_frame_eom,
            (0.0, period),
            s0,
            args=(MU, e2),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        )
        return np.asarray(out.y[:4, -1])

    base = np.array([x, 0.0, 0.0, vy, THETA0])
    stm = np.zeros((4, 4))
    h = 1e-7
    for j in range(4):
        sp = base.copy()
        sp[j] += h
        sm = base.copy()
        sm[j] -= h
        stm[:, j] = (_prop(sp) - _prop(sm)) / (2.0 * h)
    eig_mod_max = float(np.max(np.abs(np.linalg.eigvals(stm))))
    assert eig_mod_max < 1.01, f"|eig|max={eig_mod_max:.5f} (expected stable, on unit circle)"
