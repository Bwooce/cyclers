"""#928 stage 3: the KS transition matrix (variational equations, projection, fixed-time
correction) against independent references.

Recipe: Stiefel & Scheifele 1971 digest section 12 (the book prints no KS transition matrix):
10 x 10 variational equations of the (u, w, h, t) system, the lift Jacobian, the projection to
(x, xdot), and the fixed-time correction ``- zdot_f dt_f/dz0``. References:

* P = 0: ``core.kepler_stm.shepperd_stm`` (analytic, 3D) and the closed form of Deprit &
  Deprit-Bartholome 1968 (Bull. Astron. 3:315, Tables I and II, p.328; planar), coded below from
  the digest ``docs/notes/2026-10-04-digest-deprit-deprit-bartholome-1968-kepler-matrizants.md``
  section 3 with the printed slip in b44 corrected (inner factor e2 = Y r^2 + G x + 3 Q t, as in
  b41 to b43; digest section 3.2, checked there against finite differences).
* CR3BP passes at q = 1e-2, 1e-3, 1e-4 of the Moon: central finite differences of the KS flow
  and of the Cartesian flow (``cr3bp_eom``, DOP853 rtol 2.3e-14); symplecticity in the canonical
  rotating-frame variables (x, p), p = v + e_z x x; determinant 1; invariance under the fibre
  angle and fibre gauge of the lift.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom
from cyclerfinder.core.cr3bp_ks import KSModel, MoonCentredCR3BP, propagate_ks
from cyclerfinder.core.kepler_stm import shepperd_stm

FloatArray = NDArray[np.float64]

MU_EM = 0.0121505856
RTOL, ATOL = 1e-13, 1e-15


@pytest.mark.parametrize(
    ("k2", "r0", "v0", "dt"),
    [
        (1.0, [1.0, 0.0, 0.0], [0.0, 1.2, 0.1], 3.0),  # ellipse
        (1.0, [1.0, 0.2, 0.0], [0.3, 1.5, 0.2], 2.0),  # hyperbola
        (1.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 7.0),  # circle, more than one revolution
        (MU_EM, [0.05, 0.01, -0.02], [-0.1, 0.6, 0.2], 0.3),  # lunar k2, hyperbolic
        (1.0, [1.0, 0.0, 0.0], [0.0, math.sqrt(2.0) * (1 - 1e-9), 0.0], 3.0),  # near-parabolic
    ],
)
def test_kepler_stm_against_shepperd(
    k2: float, r0: list[float], v0: list[float], dt: float
) -> None:
    """P = 0: the projected KS matrix equals Shepperd's analytic STM. Measured 6.4e-16 to
    9.5e-14 relative to the largest entry (the largest, 9.5e-14, on the 7-unit circle)."""
    r0a, v0a = np.array(r0), np.array(v0)
    arc = propagate_ks(
        KSModel(k2), np.concatenate([r0a, v0a]), dt, rtol=RTOL, atol=ATOL, with_stm=True
    )
    assert arc.stm is not None
    _, _, phi = shepperd_stm(r0a, v0a, dt, mu=k2)
    assert np.abs(arc.stm - phi).max() < 1e-12 * max(1.0, float(np.abs(phi).max()))


def _deprit_a(s: FloatArray, t: float, mu: float) -> FloatArray:
    """Deprit & Deprit-Bartholome 1968 Table I, A(t); s = (x, y, X, Y)."""
    x, y, vx, vy = s
    r3 = math.hypot(x, y) ** 3
    return np.array(
        [
            [2 * x - 3 * vx * t, -y, -y * vy, -x * vy + 2 * y * vx],
            [2 * y - 3 * vy * t, x, -y * vx + 2 * x * vy, -x * vx],
            [-vx + 3 * mu * x * t / r3, -vy, -(vy**2) + mu * y * y / r3, vx * vy - mu * x * y / r3],
            [-vy + 3 * mu * y * t / r3, vx, vx * vy - mu * x * y / r3, -(vx**2) + mu * x * x / r3],
        ]
    )


def _deprit_b(s: FloatArray, t: float, mu: float) -> FloatArray:
    """Table II, B(t) = A(t)^-1, with b44's inner factor corrected to e2 (printed slip)."""
    x, y, vx, vy = s
    r = math.hypot(x, y)
    r3 = r**3
    h = 0.5 * (vx * vx + vy * vy) - mu / r
    g = x * vy - y * vx
    p = g * vy - mu * x / r
    q = -g * vx - mu * y / r
    e1 = vx * r * r - g * y + 3 * p * t
    e2 = vy * r * r + g * x + 3 * q * t
    k = mu / (2 * h * g * g)
    return np.array(
        [
            [-(mu / (2 * h)) * x / r3, -(mu / (2 * h)) * y / r3, -vx / (2 * h), -vy / (2 * h)],
            [
                (-vx + 3 * mu * x * t / r3) / g,
                (-vy + 3 * mu * y * t / r3) / g,
                -(2 * x - 3 * vx * t) / g,
                -(2 * y - 3 * vy * t) / g,
            ],
            [
                -k * (x / r3) * e1,
                -k * (y / r3) * e1 + 1 / g,
                -(vx / (2 * h * g * g)) * e1 + x * x / (g * g),
                -(vy / (2 * h * g * g)) * e1 + x * y / (g * g),
            ],
            [
                -k * (x / r3) * e2 - 1 / g,
                -k * (y / r3) * e2,
                -(vx / (2 * h * g * g)) * e2 + x * y / (g * g),
                -(vy / (2 * h * g * g)) * e2 + y * y / (g * g),
            ],
        ]
    )


@pytest.mark.parametrize(
    ("mu", "s0", "dt"),
    [
        (1.0, [1.0, 0.0, 0.0, 1.15], 2.5),  # prograde ellipse
        (1.0, [0.7, -0.4, 0.5, 0.9], 9.0),  # several revolutions
        (1.0, [1.0, 0.0, 0.0, -1.3], 2.0),  # retrograde (G < 0)
        (1.0, [1.0, 0.2, 0.3, 1.5], 2.0),  # hyperbola
        (MU_EM, [0.02, 0.0, 0.05, 1.1], 0.2),  # lunar hyperbolic pass, q ~ 0.02
    ],
)
def test_kepler_stm_against_deprit_closed_form(mu: float, s0: list[float], dt: float) -> None:
    """Planar P = 0: the KS matrix equals A(t) B(t0) of Deprit's Tables I and II. Measured
    2.0e-15 to 6.6e-14 relative to the largest entry (Shepperd against the same closed form:
    2.2e-15 to 5.0e-13)."""
    st = np.array([s0[0], s0[1], 0.0, s0[2], s0[3], 0.0])
    arc = propagate_ks(KSModel(mu), st, dt, rtol=RTOL, atol=ATOL, with_stm=True)
    assert arc.stm is not None
    idx = [0, 1, 3, 4]
    sf = arc.state[idx]
    closed = _deprit_a(sf, dt, mu) @ _deprit_b(np.array(s0), 0.0, mu)
    # A(t0) B(t0) = I, to rounding relative to |A||B| entry by entry. An absolute 1e-12 cannot
    # hold for the lunar pass: B has an entry of 2.15e4, whose storage alone carries u * 2.15e4
    # = 2.4e-12. The product adds at most gamma_4 |A||B| (gamma_n = n u / (1 - n u); Higham
    # 2002, Accuracy and Stability of Numerical Algorithms, 2nd ed., section 3.5, eq. 3.13);
    # the rest of the 16 u is margin for rounding inside the entry formulas, which has no
    # strict bound in |A||B| where terms cancel. Measured: at most 1.48 u |A||B| over the five
    # cases. The printed b44 slip (e1 in place of e2) gives 2e15 u |A||B| or more.
    a0, b0 = _deprit_a(np.array(s0), 0.0, mu), _deprit_b(np.array(s0), 0.0, mu)
    u = np.finfo(np.float64).eps / 2
    assert np.all(np.abs(a0 @ b0 - np.eye(4)) <= 16 * u * (np.abs(a0) @ np.abs(b0)))
    ks4 = arc.stm[np.ix_(idx, idx)]
    assert np.abs(ks4 - closed).max() < 2e-12 * max(1.0, float(np.abs(closed).max()))
    # planar data: no coupling into z
    assert np.abs(arc.stm[np.ix_([2, 5], idx)]).max() == 0.0


def _pass_start(q: float) -> FloatArray:
    """3D state about 0.07 from the Moon whose pass has pericentre about q (as in
    test_cr3bp_ks.py)."""
    m = MoonCentredCR3BP(MU_EM)
    vp = math.sqrt(0.05**2 + 2 * MU_EM / q)
    c7, s7, c4, s4 = math.cos(0.7), math.sin(0.7), math.cos(0.4), math.sin(0.4)
    peri = np.array([1 - MU_EM + q * c7, q * s7, 0.0, -vp * s7 * c4, vp * c7 * c4, vp * s4])
    return propagate_ks(m, peri, -0.08, rtol=RTOL, atol=ATOL).state


_T_PASS = 0.16
_J6 = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
_T_XP = np.block(
    [[np.eye(3), np.zeros((3, 3))], [np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0.0]]), np.eye(3)]]
)  # (x, v) -> (x, p), p = v + e_z x x


@pytest.mark.parametrize(
    ("q", "ks_fd", "cart_fd"), [(1e-2, 2e-9, 1e-8), (1e-3, 2e-9, 3e-7), (1e-4, 2e-9, 1e-5)]
)
def test_pass_stm_against_finite_differences(q: float, ks_fd: float, cart_fd: float) -> None:
    """Through a pass of pericentre q, relative to the largest entry (|Phi| about 20):

    * against central differences (step 1e-7) of the KS flow: measured 5.5e-10, 4.2e-10, 3.7e-10;
    * against central differences (step 1e-5) of the Cartesian flow at rtol 2.3e-14: measured
      1.9e-9 (step 1e-6), 4.1e-8, 1.25e-6. The Cartesian reference degrades as q falls (its own
      integration error, divided by the step); the bounds follow it.
    """
    m = MoonCentredCR3BP(MU_EM)
    s0 = _pass_start(q)
    arc = propagate_ks(m, s0, _T_PASS, rtol=RTOL, atol=ATOL, with_stm=True)
    assert arc.stm is not None
    scale = float(np.abs(arc.stm).max())
    fd = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = 1e-7
        plus = propagate_ks(m, s0 + e, _T_PASS, rtol=RTOL, atol=ATOL).state
        minus = propagate_ks(m, s0 - e, _T_PASS, rtol=RTOL, atol=ATOL).state
        fd[:, j] = (plus - minus) / 2e-7
    assert np.abs(arc.stm - fd).max() < ks_fd * scale
    step = 1e-6 if q == 1e-2 else 1e-5

    def cart(s: FloatArray) -> FloatArray:
        sol = solve_ivp(
            cr3bp_eom, (0, _T_PASS), s, method="DOP853", rtol=2.3e-14, atol=2.3e-17, args=(MU_EM,)
        )
        return np.asarray(sol.y[:, -1], dtype=np.float64)

    for j in range(6):
        e = np.zeros(6)
        e[j] = step
        fd[:, j] = (cart(s0 + e) - cart(s0 - e)) / (2 * step)
    assert np.abs(arc.stm - fd).max() < cart_fd * scale


@pytest.mark.parametrize("q", [1e-2, 1e-3, 1e-4, 1e-6])
def test_pass_stm_symplectic_and_gauge_invariant(q: float) -> None:
    """Phi_xp^T J Phi_xp = J in (x, p) (measured below 3.1e-16 of |Phi|^2), det Phi = 1
    (measured 1.2e-14), and the matrix does not depend on the fibre angle or on a fibre
    component added to the lift Jacobian (measured 6e-15 of |Phi|)."""
    m = MoonCentredCR3BP(MU_EM)
    s0 = _pass_start(q)
    arc = propagate_ks(m, s0, _T_PASS, rtol=RTOL, atol=ATOL, with_stm=True)
    assert arc.stm is not None
    scale = float(np.abs(arc.stm).max())
    phi_xp = _T_XP @ arc.stm @ np.linalg.inv(_T_XP)
    assert np.abs(phi_xp.T @ _J6 @ phi_xp - _J6).max() < 1e-13 * scale**2
    assert np.linalg.det(arc.stm) == pytest.approx(1.0, abs=1e-12)
    other = propagate_ks(
        m,
        s0,
        _T_PASS,
        rtol=RTOL,
        atol=ATOL,
        with_stm=True,
        fibre_angle=1.3,
        gauge=np.array([0.3, -1.0, 2.0, 0.5, 0.1, -0.7]),
    )
    assert other.stm is not None
    assert np.abs(other.stm - arc.stm).max() < 1e-12 * scale
    # without the gauge, the lift's fibre angle alone
    turned = propagate_ks(m, s0, _T_PASS, rtol=RTOL, atol=ATOL, with_stm=True, fibre_angle=-2.0)
    assert turned.stm is not None
    assert np.abs(turned.stm - arc.stm).max() < 1e-12 * scale


def test_backward_stm_is_the_inverse() -> None:
    """Phi(-T) at the end state times Phi(T) is the identity (measured 2.3e-13)."""
    m = MoonCentredCR3BP(MU_EM)
    s0 = _pass_start(1e-3)
    fwd = propagate_ks(m, s0, _T_PASS, rtol=RTOL, atol=ATOL, with_stm=True)
    back = propagate_ks(m, fwd.state, -_T_PASS, rtol=RTOL, atol=ATOL, with_stm=True)
    assert fwd.stm is not None and back.stm is not None
    assert np.abs(back.stm @ fwd.stm - np.eye(6)).max() < 1e-11


def test_stm_needs_an_autonomous_potential_model() -> None:
    class Forced(KSModel):
        def force(self, t: float, x: FloatArray) -> FloatArray:
            return np.array([1e-3, 0.0, 0.0])

    with pytest.raises(NotImplementedError):
        propagate_ks(Forced(1.0), np.array([1.0, 0, 0, 0, 1, 0]), 0.5, with_stm=True)
