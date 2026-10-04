"""#896/#893: the concentric circular restricted four-body model against a published torus.

Source: B. Kumar, R. L. Anderson, R. de la Llave and B. Gunter, "Computation and analysis of
Jupiter-Europa and Jupiter-Ganymede resonant orbits in the planar concentric circular restricted
4-body problem", AAS/AIAA Astrodynamics Specialist Conference, AAS 21-651 (2021),
arXiv:2109.14815; filed in the private paper corpus as
kumar-anderson-delallave-gunter-2021-europa-ganymede-resonant-orbits-ccr4bp-AAS-21-651-arxiv-2109.14815.pdf.

What the paper prints, and what is checked here:

* p. 2, Eqs. (1)-(2): the model. Europa (m2) and Ganymede (m3) move on circles about Jupiter at
  the Kepler rates ``Omega_i = sqrt(G(m1 + m_i) / r_1i^3)``; units ``r12 = G(m1 + m2) = 1``;
  ``mu = m2/(m1 + m2)``, ``mu3 = m3/(m1 + m2)``; Ganymede's synodic angle is
  ``theta3 = (Omega3 - 1) t + theta3_0``; the Ganymede terms are a direct term plus the indirect
  term ``-mu3 (cos theta3, sin theta3) / r13^2``. ``core/ccr4bp.py`` has the same terms with
  ``mu_gan = mu3``, ``omega_gan = Omega3 - 1`` and ``a_gan = r13``, so the paper's system is built
  here from the paper's own constants. One difference: the paper centres Ganymede's circle on
  Jupiter, ``(-mu + r13 cos theta3, r13 sin theta3)``, and ``core/ccr4bp.py`` centres it on the
  Jupiter-Europa barycentre. Solving the torus below with the paper's centring instead (scratch
  computation, not kept as a test) moves the closest approach by 0.015 km.
* p. 8, Table 1: GM of Jupiter, Europa and Ganymede and the two moons' periods. The two GM
  values in the "Europa" and "Ganymede" rows are interchanged (9.887e12 m^3/s^2 is Ganymede's,
  3.201e12 is Europa's); with the rows read the other way round the table reproduces both mass
  ratios the paper prints elsewhere to all 16 digits (``mu3`` on p. 8, ``mu_bar2`` on p. 12).
* p. 8: "Starting from a periodic orbit with Jacobi constant value 3.0041 and omega = 3.097849
  ... the final value of mu3 = 7.804102777055038e-5 ... the closest approach to Europa decreases
  from 22052 km to 18721 km". The Figure 1 title gives the same rotation number to 19 digits,
  3.097848962221668715. The rotation number fixes the period of the Jupiter-Europa 3:4 orbit
  through Eq. (9), ``omega = 2 pi Omega1 / |Omega3 - 1|``, so ``T = 2 pi Tp / omega`` with
  ``Tp = 2 pi / |Omega3 - 1|``; with ``Omega3`` fixed the paper's Eq. (1) makes ``r13`` change
  with ``mu3`` (p. 6).
* p. 9, Figure 2: the unstable multiplier of the stroboscopic map, read from the plot as about
  7.18 at ``mu3 = 0`` and about 6.915 at the end of the curve (``mu3`` near 7.8e-5). These are
  figure readings, so they are checked to 0.02.
* p. 12: in the Jupiter-Ganymede frame the 3:4 torus frequency tends to
  ``Omega1 = 0.503493`` and the rotation number to ``omega = 3.119948`` as Europa's mass goes to
  zero; both follow from Table 1's periods alone.

How the torus is computed. The 2-torus is parameterized by the orbit angle (frequency
``Omega1 = 2 pi / T``) and Ganymede's synodic angle. Its section at orbit angle zero is a curve
``K(phi)`` indexed by Ganymede's phase ``phi``; the flow over exactly ``T`` maps it to itself with
``phi`` advanced by ``(Omega3 - 1) T``. ``K`` is discretized at N phases (trigonometric
interpolation for the shift), the invariance equations plus one phase condition (mean ``y`` = 0)
are solved by damped Gauss-Newton, and ``mu3`` is raised from 0 to the paper's value in eight
steps with a secant predictor. For speed the Newton iterations integrate all N points together
with a vectorized copy of ``core/ccr4bp.py``'s right-hand side, which a test below checks against
``ccr4bp.ccr4bp_eom`` point by point; the converged curve is then checked for invariance with
``ccr4bp.propagate_ccr4bp`` itself, and the closest approach is measured by propagating with
``core/ccr4bp.py``. The paper's invariant circle is the stroboscopic section at Ganymede phase
zero (Ganymede on the +x axis, aligned with Europa); its point nearest Europa is found by flowing
the section curve points ``K(-(Omega3 - 1) t)`` for time ``t`` near zero.

Measured (2026-10-04): base orbit x0 = 1.03283365, C = 3.0041057, closest approach 22051.69 km,
stroboscopic multiplier 7.17997. Torus at the paper's mu3 with N = 127: invariance residual
3e-10 in the project model, closest approach 18721.37 km (N = 95: 18722.13; N = 191: 18721.26),
stroboscopic multiplier 6.9206. With Ganymede's synodic rate reversed (control) the same
computation gives 21774.8 km.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp

FloatArray = NDArray[np.float64]

# --- Table 1 (p. 8), digit for digit. The GM rows are interchanged in the paper. ---
GM_JUPITER = 1.2668653785779600e17
GM_PRINTED_EUROPA_ROW = 9.8869974284299492e12
GM_PRINTED_GANYMEDE_ROW = 3.2009998067205903e12
PERIOD_EUROPA_S = 3.0689648366400000e5
PERIOD_GANYMEDE_S = 6.1808096312640002e5

# --- Printed values ---
MU3_PRINTED = 7.804102777055038e-5  # p. 8
MU_BAR2_PRINTED = 2.5265115494603433e-5  # p. 12
OMEGA_ROT_PRINTED = 3.097848962221668715  # Figure 1 title (text: 3.097849)
JACOBI_PRINTED = 3.0041  # p. 8
CLOSEST_MU3_ZERO_KM = 22052.0  # p. 8
CLOSEST_MU3_PHYSICAL_KM = 18721.0  # p. 8
OMEGA1_LIMIT_PRINTED = 0.503493  # p. 12
OMEGA_LIMIT_PRINTED = 3.119948  # p. 12
MULTIPLIER_MU3_ZERO_FIG2 = 7.18  # p. 9, Figure 2, read from the plot
MULTIPLIER_MU3_PHYSICAL_FIG2 = 6.915  # p. 9, Figure 2, read from the plot

# --- Derived from the paper's constants (Europa = the smaller GM) ---
GM_EUROPA = GM_PRINTED_GANYMEDE_ROW
GM_GANYMEDE = GM_PRINTED_EUROPA_ROW
MU = GM_EUROPA / (GM_JUPITER + GM_EUROPA)
MU3 = GM_GANYMEDE / (GM_JUPITER + GM_EUROPA)
OMEGA3 = PERIOD_EUROPA_S / PERIOD_GANYMEDE_S  # Ganymede's inertial rate in Europa units
T_FORCING = 2.0 * math.pi / abs(OMEGA3 - 1.0)  # Tp, stroboscopic period
T_ORBIT = 2.0 * math.pi * T_FORCING / OMEGA_ROT_PRINTED  # Eq. (9) with Omega1 = 2 pi / T
LENGTH_UNIT_KM = (
    (GM_JUPITER + GM_EUROPA) ** (1.0 / 3.0) * (PERIOD_EUROPA_S / (2.0 * math.pi)) ** (2.0 / 3.0)
) / 1.0e3

N_CURVE = 127
N_MASS_STEPS = 8
_IDX4 = [0, 1, 3, 4]


def _r13(mu3: float) -> float:
    """Ganymede's radius from the paper's Eq. (1) with Omega3 held fixed (p. 6)."""
    return float(((1.0 - MU + mu3) / OMEGA3**2) ** (1.0 / 3.0))


def _system(mu3: float, theta0: float = 0.0, sense: float = 1.0) -> ccr4bp.CCR4BPSystem:
    return ccr4bp.CCR4BPSystem(
        mu=MU,
        mu_gan=mu3,
        a_gan=_r13(mu3),
        omega_gan=sense * (OMEGA3 - 1.0),
        theta_gan0=theta0,
    )


# ---------------------------------------------------------------------------
# Vectorized right-hand side (same terms as ccr4bp.ccr4bp_eom, planar) for the Newton solves.
# ---------------------------------------------------------------------------


def _batch_rhs(
    mu3: float, a_gan: float, omega_gan: float, phis: FloatArray, with_stm: bool
) -> Callable[[float, FloatArray], FloatArray]:
    n = phis.size

    def rhs(t: float, y: FloatArray) -> FloatArray:
        s = y[: 4 * n].reshape(n, 4)
        x, yy, vx, vy = s[:, 0], s[:, 1], s[:, 2], s[:, 3]
        th = phis + omega_gan * t
        cth, sth = np.cos(th), np.sin(th)
        gx, gy = a_gan * cth, a_gan * sth
        d1x, d2x, dgx, dgy = x + MU, x - 1.0 + MU, x - gx, yy - gy
        r1s, r2s, r3s = d1x * d1x + yy * yy, d2x * d2x + yy * yy, dgx * dgx + dgy * dgy
        r1c, r2c, r3c = r1s * np.sqrt(r1s), r2s * np.sqrt(r2s), r3s * np.sqrt(r3s)
        ax = (
            2.0 * vy
            + x
            - (1.0 - MU) * d1x / r1c
            - MU * d2x / r2c
            - mu3 * dgx / r3c
            - mu3 * cth / a_gan**2
        )
        ay = (
            -2.0 * vx
            + yy
            - (1.0 - MU) * yy / r1c
            - MU * yy / r2c
            - mu3 * dgy / r3c
            - mu3 * sth / a_gan**2
        )
        ds = np.stack([vx, vy, ax, ay], axis=1).ravel()
        if not with_stm:
            return ds
        r15, r25, r35 = r1s * r1c, r2s * r2c, r3s * r3c
        common = 1.0 - (1.0 - MU) / r1c - MU / r2c - mu3 / r3c
        uxx = common + 3 * (1 - MU) * d1x**2 / r15 + 3 * MU * d2x**2 / r25 + 3 * mu3 * dgx**2 / r35
        uyy = common + 3 * (1 - MU) * yy**2 / r15 + 3 * MU * yy**2 / r25 + 3 * mu3 * dgy**2 / r35
        uxy = 3 * (1 - MU) * d1x * yy / r15 + 3 * MU * d2x * yy / r25 + 3 * mu3 * dgx * dgy / r35
        a = np.zeros((n, 4, 4))
        a[:, 0, 2] = 1.0
        a[:, 1, 3] = 1.0
        a[:, 2, 0], a[:, 2, 1], a[:, 3, 0], a[:, 3, 1] = uxx, uxy, uxy, uyy
        a[:, 2, 3], a[:, 3, 2] = 2.0, -2.0
        dphi = a @ y[4 * n :].reshape(n, 4, 4)
        return np.concatenate([ds, dphi.ravel()])

    return rhs


def _batch_flow(
    u: FloatArray, mu3: float, omega_gan: float, phis: FloatArray, with_stm: bool
) -> tuple[FloatArray, FloatArray | None]:
    n = phis.size
    y0 = u.ravel()
    if with_stm:
        y0 = np.concatenate([y0, np.tile(np.eye(4), (n, 1, 1)).ravel()])
    rhs = _batch_rhs(mu3, _r13(mu3), omega_gan, phis, with_stm)
    sol = solve_ivp(rhs, (0.0, T_ORBIT), y0, method="DOP853", rtol=1e-12, atol=1e-12)
    yf = sol.y[:, -1]
    if with_stm:
        return yf[: 4 * n].reshape(n, 4), yf[4 * n :].reshape(n, 4, 4)
    return yf.reshape(n, 4), None


def _shift_matrix(n: int, shift: float) -> FloatArray:
    """Values at ``phi_j + shift`` from values at ``phi_j`` by trigonometric interpolation."""
    k = np.fft.fftfreq(n, 1.0 / n)
    f = np.fft.fft(np.eye(n), axis=0)
    return np.asarray(np.real(np.fft.ifft(np.exp(1j * k * shift)[:, None] * f, axis=0)))


def _interp(u: FloatArray, phi: float) -> FloatArray:
    n = u.shape[0]
    k = np.fft.fftfreq(n, 1.0 / n)
    c = np.fft.fft(u, axis=0)
    return np.asarray(np.real(np.sum(c * np.exp(1j * k * phi)[:, None], axis=0) / n))


def _solve_curve(u: FloatArray, mu3: float, omega_gan: float) -> tuple[FloatArray, float]:
    n = u.shape[0]
    phis = 2.0 * np.pi * np.arange(n) / n
    shift = _shift_matrix(n, omega_gan * T_ORBIT)
    res_max = math.inf
    for _ in range(20):
        uf, stm = _batch_flow(u, mu3, omega_gan, phis, True)
        assert stm is not None
        res = np.concatenate([(uf - shift @ u).ravel(), [float(np.mean(u[:, 1]))]])
        res_max = float(np.max(np.abs(res)))
        if res_max < 1e-8:
            break
        jac = np.zeros((4 * n + 1, 4 * n))
        for j in range(n):
            jac[4 * j : 4 * j + 4, 4 * j : 4 * j + 4] = stm[j]
        jac[: 4 * n] -= np.kron(shift, np.eye(4))
        jac[4 * n, 1::4] = 1.0 / n
        step = np.linalg.lstsq(jac, -res, rcond=None)[0].reshape(n, 4)
        lam = 1.0
        while True:
            trial = u + lam * step
            uft, _ = _batch_flow(trial, mu3, omega_gan, phis, False)
            r_trial = max(
                float(np.max(np.abs(uft - shift @ trial))), abs(float(np.mean(trial[:, 1])))
            )
            if r_trial < res_max or lam < 1e-3:
                break
            lam *= 0.5
        u = trial
    return u, res_max


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BaseOrbit:
    state: FloatArray  # (x, y, z, vx, vy, vz) at the perpendicular crossing nearest Europa
    jacobi: float
    closest_km: float
    multiplier_t: float  # unstable multiplier over one orbit period T


@dataclass(frozen=True)
class Torus:
    curve: FloatArray  # (N, 4) section curve K(phi), phi = 2 pi j / N
    mu3: float
    omega_gan: float
    newton_residual: float


def _half_period_crossing(x0: float, vy0: float) -> tuple[FloatArray, FloatArray]:
    y0 = np.concatenate([np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), np.eye(6).ravel()])
    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom,
        (0.0, 0.5 * T_ORBIT),
        y0,
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
    )
    return sol.y[:6, -1], sol.y[6:, -1].reshape(6, 6)


@pytest.fixture(scope="module")
def base_orbit() -> BaseOrbit:
    """The Jupiter-Europa 3:4 orbit with the period fixed by the printed rotation number.

    Seed: the printed closest approach on the +x axis beyond Europa, with the rotating-frame
    speed of a Kepler 3:4 ellipse at that perijove. Two-variable symmetric Newton:
    ``y = vx = 0`` at ``T/2``."""
    x0 = 1.0 - MU + CLOSEST_MU3_ZERO_KM / LENGTH_UNIT_KM
    sma = (4.0 / 3.0) ** (2.0 / 3.0)
    vy0 = math.sqrt(2.0 / x0 - 1.0 / sma) - x0
    for _ in range(20):
        sf, m = _half_period_crossing(x0, vy0)
        r = np.array([sf[1], sf[3]])
        jac = np.array([[m[1, 0], m[1, 4]], [m[3, 0], m[3, 4]]])
        d = np.linalg.solve(jac, -r)
        x0, vy0 = x0 + d[0], vy0 + d[1]
        if np.max(np.abs(d)) < 1e-13:
            break
    s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom,
        (0.0, T_ORBIT),
        np.concatenate([s0, np.eye(6).ravel()]),
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert np.max(np.abs(sol.y[:6, -1] - s0)) < 1e-8, "base orbit does not close"
    tt = np.linspace(0.0, T_ORBIT, 100001)
    dense = sol.sol
    assert dense is not None
    yy = dense(tt)
    closest = float(np.min(np.hypot(yy[0] - (1.0 - MU), yy[1]))) * LENGTH_UNIT_KM
    mono = sol.y[6:, -1].reshape(6, 6)[np.ix_(_IDX4, _IDX4)]
    lam = float(np.max(np.abs(np.linalg.eigvals(mono))))
    return BaseOrbit(s0, cr3bp.jacobi_constant(s0, MU), closest, lam)


def _continue_torus(base: BaseOrbit, n: int, sense: float) -> Torus:
    omega_gan = sense * (OMEGA3 - 1.0)
    u = np.tile(base.state[_IDX4], (n, 1))
    u_prev: FloatArray | None = None
    res = math.inf
    mu3 = 0.0
    for frac in np.linspace(0.0, 1.0, N_MASS_STEPS + 1)[1:]:
        mu3 = float(frac) * MU3
        guess = u if u_prev is None else 2.0 * u - u_prev
        u_new, res = _solve_curve(guess, mu3, omega_gan)
        u_prev, u = u, u_new
    return Torus(u, mu3, omega_gan, res)


@pytest.fixture(scope="module")
def torus(base_orbit: BaseOrbit) -> Torus:
    return _continue_torus(base_orbit, N_CURVE, 1.0)


def _section_closest_km(tor: Torus) -> float:
    """Closest approach to Europa of the stroboscopic circle at Ganymede phase zero, measured
    with ``core/ccr4bp.py``: the circle point at orbit angle ``Omega1 t`` is the section curve
    point at Ganymede phase ``-omega_gan t`` flowed for time ``t``."""

    def dist(t: float) -> float:
        phi = (-tor.omega_gan * t) % (2.0 * math.pi)
        u = _interp(tor.curve, phi)
        s6 = np.array([u[0], u[1], 0.0, u[2], u[3], 0.0])
        if t != 0.0:
            sys_t = ccr4bp.CCR4BPSystem(
                mu=MU, mu_gan=tor.mu3, a_gan=_r13(tor.mu3), omega_gan=tor.omega_gan, theta_gan0=phi
            )
            s6 = ccr4bp.propagate_ccr4bp(sys_t, s6, t).state_f
        return math.hypot(s6[0] - (1.0 - MU), s6[1])

    best = minimize_scalar(dist, bounds=(-0.05, 0.05), method="bounded", options={"xatol": 1e-7})
    return float(best.fun) * float(LENGTH_UNIT_KM)


def _strob_multiplier(tor: Torus, n_steps: int = 4000) -> float:
    """Unstable multiplier of the stroboscopic map: growth rate of the linear cocycle
    ``v -> DPhi_T(K(phi)) v``, ``phi -> phi + omega_gan T``, converted from per-``T`` to
    per-``Tp``. The monodromy blocks at the grid phases are trigonometrically interpolated."""
    n = tor.curve.shape[0]
    phis = 2.0 * np.pi * np.arange(n) / n
    _, stm = _batch_flow(tor.curve, tor.mu3, tor.omega_gan, phis, True)
    assert stm is not None
    coef = np.fft.fft(stm.reshape(n, 16), axis=0)
    k = np.fft.fftfreq(n, 1.0 / n)
    v = np.array([1.0, 0.0, 0.0, 0.0])
    log_sum, phi, burn = 0.0, 0.0, 200
    for i in range(n_steps):
        m = np.real(np.sum(coef * np.exp(1j * k * phi)[:, None], axis=0) / n).reshape(4, 4)
        v = m @ v
        nv = float(np.linalg.norm(v))
        v /= nv
        if i >= burn:
            log_sum += math.log(nv)
        phi = (phi + tor.omega_gan * T_ORBIT) % (2.0 * math.pi)
    lam_t = math.exp(log_sum / (n_steps - burn))
    return float(lam_t ** (T_FORCING / T_ORBIT))


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


def test_table1_reproduces_printed_mass_ratios_with_gm_rows_interchanged() -> None:
    """Table 1 gives mu3 (p. 8) and mu_bar2 (p. 12) to 16 digits once its two moon GM values
    are read in each other's row; read as printed it gives the two ratios swapped."""
    assert pytest.approx(MU3_PRINTED, rel=1e-15) == MU3
    assert pytest.approx(MU_BAR2_PRINTED, rel=1e-15) == GM_EUROPA / (GM_JUPITER + GM_GANYMEDE)
    as_printed = GM_PRINTED_GANYMEDE_ROW / (GM_JUPITER + GM_PRINTED_EUROPA_ROW)
    assert as_printed == pytest.approx(MU_BAR2_PRINTED, rel=1e-15)
    assert abs(as_printed / MU3_PRINTED - 1.0) > 0.6


def test_printed_frequencies_follow_from_table1_periods() -> None:
    """p. 12: Omega1 -> 0.503493 and omega -> 3.119948 in the Jupiter-Ganymede frame for the
    3:4 Jupiter-Europa torus as Europa's mass goes to zero. In that frame time is scaled by
    Ganymede's rate, Europa's rate is ``Omega2 = T_G/T_E``, the 3:4 orbit's period is 8 pi in
    Europa units, and Eq. (9) gives the rotation number."""
    omega2 = PERIOD_GANYMEDE_S / PERIOD_EUROPA_S
    omega1 = (2.0 * math.pi / (8.0 * math.pi)) * omega2
    assert omega1 == pytest.approx(OMEGA1_LIMIT_PRINTED, abs=5e-7)
    assert 2.0 * math.pi * omega1 / abs(omega2 - 1.0) == pytest.approx(
        OMEGA_LIMIT_PRINTED, abs=5e-7
    )


def test_model_rate_formula_is_the_papers_eq1() -> None:
    """``ccr4bp.two_body_synodic_rate`` is the paper's Eq. (1) in the paper's units: with
    ``r13`` from Eq. (1) at the paper's ``Omega3`` it returns ``Omega3 - 1``."""
    rate = ccr4bp.two_body_synodic_rate(MU, MU3, _r13(MU3))
    assert rate == pytest.approx(OMEGA3 - 1.0, abs=1e-15)


def test_project_default_constants_are_close_to_the_papers() -> None:
    """The registry-built default system against the paper's Table 1 system. Measured
    (2026-10-04): mass ratios agree to 5.4e-4 (Europa) and 0.9e-4 (Ganymede) relative; the
    default's Ganymede synodic rate is -0.5035527 against the paper's -0.5034688 from the
    printed periods, 1.7e-4 relative, because the default derives the rate from the registry
    semi-major axes by Kepler's law rather than from the observed periods. Over one 3:4 orbit
    (T = 25.3) that is a phase difference of 2.1e-3 rad. The bounds are a reconciliation
    record, not printed numbers."""
    d = ccr4bp.jupiter_europa_ganymede_default()
    assert d.mu == pytest.approx(MU, rel=1e-3)
    assert d.mu_gan == pytest.approx(MU3, rel=1e-3)
    assert d.omega_gan == pytest.approx(OMEGA3 - 1.0, rel=3e-4)


def test_vectorized_rhs_matches_project_model() -> None:
    """The batch right-hand side used inside the Newton solves is ``ccr4bp.ccr4bp_eom``."""
    rng = np.random.default_rng(896)
    phis = rng.uniform(0.0, 2.0 * math.pi, 5)
    u = np.column_stack(
        [
            rng.uniform(0.6, 1.4, 5),
            rng.uniform(-0.4, 0.4, 5),
            rng.uniform(-0.2, 0.2, 5),
            rng.uniform(-0.2, 0.2, 5),
        ]
    )
    t = 3.7
    batch = _batch_rhs(MU3, _r13(MU3), OMEGA3 - 1.0, phis, False)(t, u.ravel()).reshape(5, 4)
    for j in range(5):
        s6 = np.array([u[j, 0], u[j, 1], 0.0, u[j, 2], u[j, 3], 0.0])
        ref = ccr4bp.ccr4bp_eom(t, s6, _system(MU3, phis[j]))
        np.testing.assert_allclose(batch[j], ref[_IDX4], rtol=1e-13, atol=1e-15)


# ---------------------------------------------------------------------------
# The orbit at mu3 = 0 (Jupiter-Europa three-body problem)
# ---------------------------------------------------------------------------


def test_base_orbit_jacobi_constant_matches_print(base_orbit: BaseOrbit) -> None:
    """The 3:4 orbit with the period fixed by the printed rotation number has the printed
    Jacobi constant 3.0041 (measured 3.0041057)."""
    assert base_orbit.jacobi == pytest.approx(JACOBI_PRINTED, abs=5e-5)


def test_base_orbit_closest_approach_matches_print(base_orbit: BaseOrbit) -> None:
    """22052 km printed; measured 22051.69 km."""
    assert base_orbit.closest_km == pytest.approx(CLOSEST_MU3_ZERO_KM, abs=0.5)


def test_base_orbit_multiplier_matches_figure2_start(base_orbit: BaseOrbit) -> None:
    """At mu3 = 0 the stroboscopic map's multiplier is the orbit's multiplier raised to
    ``Tp/T``; Figure 2 starts at about 7.18 (measured 7.1800)."""
    strob = base_orbit.multiplier_t ** (T_FORCING / T_ORBIT)
    assert strob == pytest.approx(MULTIPLIER_MU3_ZERO_FIG2, abs=0.02)


# ---------------------------------------------------------------------------
# The torus at the paper's Ganymede mass
# ---------------------------------------------------------------------------


def test_torus_is_invariant_in_project_model(torus: Torus) -> None:
    """Every section-curve point, propagated for T with ``ccr4bp.propagate_ccr4bp``, lands on
    the curve at the shifted Ganymede phase (measured 3e-10)."""
    n = torus.curve.shape[0]
    assert torus.mu3 == MU3
    shift = _shift_matrix(n, torus.omega_gan * T_ORBIT)
    target = shift @ torus.curve
    worst = 0.0
    for j in range(n):
        u = torus.curve[j]
        s6 = np.array([u[0], u[1], 0.0, u[2], u[3], 0.0])
        sf = ccr4bp.propagate_ccr4bp(_system(MU3, 2.0 * math.pi * j / n), s6, T_ORBIT).state_f
        worst = max(worst, float(np.max(np.abs(sf[_IDX4] - target[j]))))
    assert worst < 1e-8


def test_torus_closest_approach_matches_print(torus: Torus) -> None:
    """The paper's 18721 km (p. 8); measured 18721.37 km with N = 127."""
    assert _section_closest_km(torus) == pytest.approx(CLOSEST_MU3_PHYSICAL_KM, abs=0.5)


def test_torus_multiplier_matches_figure2_end(torus: Torus) -> None:
    """Figure 2 ends at about 6.915; measured 6.9206."""
    assert _strob_multiplier(torus) == pytest.approx(MULTIPLIER_MU3_PHYSICAL_FIG2, abs=0.02)


def test_reversed_ganymede_rate_misses_the_print(base_orbit: BaseOrbit) -> None:
    """Control: the same computation with Ganymede's synodic rate reversed (the error class of
    #891) gives 21774.8 km, about 3000 km from the print. A coarser curve (N = 63, about 10 km
    truncation error) is enough to show it."""
    wrong = _continue_torus(base_orbit, 63, -1.0)
    assert abs(_section_closest_km(wrong) - CLOSEST_MU3_PHYSICAL_KM) > 1000.0
