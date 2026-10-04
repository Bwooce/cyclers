"""Gates for the #889 PERTBP stroboscopic machinery (identities and printed values only).

Golden values come from the published papers, filed in the private paper corpus as
    kumar-anderson-delallave-2021-highorder-resonant-manifold-expansions-cnsns-arxiv-2109.14800.pdf
("CNSNS": Sec. 3 mass ratio, Table 1 orbits) and
    kumar-anderson-delallave-2025-gpu-connections-tori-perturbed-crtbp-siam-ads-24-219-arxiv-2109.14814.pdf
("SIAM": Sec. 8.3 rotation numbers and periods); everything else is a mathematical identity.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cyclerfinder.core.er3bp_paper_frame import paper_frame_eom
from cyclerfinder.search import pertbp_strob_889 as m

MU = m.MU_JE_CNSNS
E = m.ECC_EUROPA
RNG = np.random.default_rng(889)


@pytest.fixture(scope="module", autouse=True)
def _threads() -> None:
    m.set_threads(2)


def _sample_states(k: int) -> np.ndarray:
    """States on 3:4 / 5:6-like ellipses away from both primaries."""
    out = []
    for _ in range(k):
        r = RNG.uniform(1.1, 1.4)
        a = RNG.uniform(0, 2 * math.pi)
        vt = RNG.uniform(-0.35, -0.15)  # retrograde-ish in the rotating frame
        vr = RNG.uniform(-0.1, 0.1)
        x, y = r * math.cos(a), r * math.sin(a)
        vx = vr * math.cos(a) - vt * math.sin(a)
        vy = vr * math.sin(a) + vt * math.cos(a)
        out.append(m.velocity_to_momentum(0.0, np.array([x, y, vx, vy]), 0.0))
    return np.array(out)


def test_integrator_matches_scipy_dop853() -> None:
    s = _sample_states(1)[0]

    def f(t: float, z: np.ndarray) -> np.ndarray:
        return m.rhs(t, z, MU, E)

    y, _ = m.flow(s, 0.7, 0.7 + 4 * math.pi, MU, E)
    sol = solve_ivp(f, (0.7, 0.7 + 4 * math.pi), s, method="DOP853", rtol=1e-13, atol=1e-13)
    assert np.max(np.abs(y - sol.y[:, -1])) < 1e-11


def test_e0_conserves_jacobi() -> None:
    """At e = 0 the model is the PCRTBP: Jacobi constant conserved to 1e-10."""
    s = _sample_states(4)
    c0 = m.jacobi(s, MU)
    y, _ = m.flow(s, 0.0, 10 * 2 * math.pi, MU, 0.0)
    assert np.max(np.abs(m.jacobi(y, MU) - c0)) < 1e-10


def test_e0_rhs_is_pcrtbp_in_momentum_form() -> None:
    """At e = 0, n = 1, rho = 1: SIAM Eq. 3.5 reduces to Eq. 3.2 term by term."""
    s = _sample_states(1)[0]
    x, y, px, py = s
    r1 = math.hypot(x + MU, y)
    r2 = math.hypot(x - 1 + MU, y)
    expect = np.array(
        [
            px + y,
            py - x,
            py - (1 - MU) * (x + MU) / r1**3 - MU * (x - 1 + MU) / r2**3,
            -px - (1 - MU) * y / r1**3 - MU * y / r2**3,
        ]
    )
    assert np.max(np.abs(m.rhs(1.234, s, MU, 0.0) - expect)) < 1e-15


def test_matches_independent_lagrangian_er3bp_frame() -> None:
    """The Hamiltonian flow agrees with the independently derived non-pulsating
    Antoniadou-Libert Lagrangian EOM (``core.er3bp_paper_frame``), positions identical
    and velocities by ``xdot = px + n y``, ``ydot = py - n x``. ``e`` is exaggerated to
    0.2 so that a wrong ``n(t)`` or ``rho(t)`` cannot hide."""
    ee = 0.2
    s = _sample_states(1)[0]
    t0, t1 = 0.4, 0.4 + 3.0
    y_h, _ = m.flow(s, t0, t1, MU, ee)
    # true anomaly at t0 for the 5th coordinate of the Lagrangian model
    ecc = t0
    for _ in range(60):
        ecc -= (ecc - ee * math.sin(ecc) - t0) / (1 - ee * math.cos(ecc))
    f0 = 2 * math.atan2(
        math.sqrt(1 + ee) * math.sin(ecc / 2), math.sqrt(1 - ee) * math.cos(ecc / 2)
    )
    v0 = m.momentum_to_velocity(t0, s, ee)
    sol = solve_ivp(
        lambda t, z: paper_frame_eom(t, z, MU, ee),
        (t0, t1),
        np.array([*v0, f0]),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    v_l = sol.y[:4, -1]
    v_h = m.momentum_to_velocity(t1, y_h, ee)
    assert np.max(np.abs(v_l - v_h)) < 1e-9


def test_map_inverse() -> None:
    s = _sample_states(3)
    f1, _ = m.strob(s, 1, MU, E)
    back, _ = m.strob(f1, -1, MU, E)
    assert np.max(np.abs(back - s)) < 1e-11


def test_periodic_in_time() -> None:
    """The vector field is 2 pi periodic: flowing from 2 pi k equals flowing from 0."""
    s = _sample_states(2)
    a, _ = m.flow(s, 0.3, 2.3, MU, E)
    b, _ = m.flow(s, 0.3 + 6 * math.pi, 2.3 + 6 * math.pi, MU, E)
    assert np.max(np.abs(a - b)) < 1e-12


def test_stm_matches_central_differences() -> None:
    s = _sample_states(1)[0]
    _, phi = m.strob(s, 1, MU, E, with_stm=True)
    assert phi is not None
    h = 1e-6
    fd = np.empty((4, 4))
    for j in range(4):
        d = np.zeros(4)
        d[j] = h
        p, _ = m.strob(s + d, 1, MU, E)
        q, _ = m.strob(s - d, 1, MU, E)
        fd[:, j] = (p - q) / (2 * h)
    assert np.max(np.abs(fd - phi)) / np.max(np.abs(phi)) < 1e-7


def test_map_is_symplectic_in_canonical_coordinates() -> None:
    s = _sample_states(2)
    _, phi = m.strob(s, 2, MU, E, with_stm=True)
    assert phi is not None
    for p in phi:
        assert m.symplectic_defect(p) < 1e-10


def test_time_reversal_symmetry_at_periapse() -> None:
    """CMDA Sec. 5.3: with t = 0 at periapse, M F^-1(M z) = F(z), M = diag(1,-1,-1,1)."""
    s = _sample_states(2)
    f1, _ = m.strob(s, 1, MU, E)
    g, _ = m.strob(m.REVERSOR * s, -1, MU, E)
    assert np.max(np.abs(m.REVERSOR * g - f1)) < 1e-11


@pytest.mark.parametrize("name", ["3:4", "5:6"])
def test_cnsns_table1_orbits_reproduced(name: str) -> None:
    """CNSNS Table 1 at the CNSNS mass ratio: period, ydot, multipliers, C = 3.0024."""
    tab = m.TABLE1_3_4 if name == "3:4" else m.TABLE1_5_6
    po = m.po_at_fixed_x(MU, tab["x"], tab["ydot"], tab["T"])
    ydot = po.state[3] - po.state[0]
    assert abs(ydot - tab["ydot"]) < 1e-12
    assert abs(po.period - tab["T"]) < 1e-8
    # Table 1 prints 16 digits; the 5:6 multipliers agree to 1.1e-7 relative, the
    # 3:4 ones to 2e-9 (the printed 5:6 period itself agrees to 4e-9 only).
    assert abs(po.lam_u / tab["lam_u"] - 1) < 1e-6
    assert abs(po.lam_s / tab["lam_s"] - 1) < 1e-6
    assert abs(po.jacobi - m.TABLE1_JACOBI) < 1e-11


def test_table1_periods_give_the_printed_rotation_numbers() -> None:
    """SIAM Sec. 8.3 prints omega_u = 1.558039, omega_s = 1.030011 and the periods
    38.3281 (5:6) and 25.3376-25.3394 (3:4); omega = 4 pi^2 / T (CMDA Sec. 4.9)."""
    assert round(m.omega_from_period(m.TABLE1_3_4["T"]), 6) == m.PRINTED_OMEGA_U
    assert round(m.omega_from_period(m.TABLE1_5_6["T"]), 6) == m.PRINTED_OMEGA_S
    # SIAM's printed 5:6 period matches the Table 1 orbit, not 4 pi^2 / 1.030011 (= 38.32815):
    # the printed omegas are roundings of the Table 1 values.
    assert round(m.TABLE1_5_6["T"], 4) == 38.3281
    assert round(m.period_from_omega(m.PRINTED_OMEGA_S), 4) != 38.3281
    assert 25.3376 <= m.TABLE1_3_4["T"] <= 25.3394


def test_unperturbed_circle_is_invariant() -> None:
    """At e = 0 the periodic orbit sampled at T theta / 2 pi is an invariant circle
    of the map (identity, up to interpolation), and is reversible."""
    tab = m.TABLE1_3_4
    po = m.po_at_fixed_x(MU, tab["x"], tab["ydot"], tab["T"])
    c = m.seed_circle(po, 511)
    on, off = m.circle_residual(c)
    assert on < 1e-11 and off < 1e-10
    assert m.reversibility_defect(c) < 1e-10
