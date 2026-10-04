"""#896: the periodic orbit that replaces L2 in the bicircular problem, against a printed table.

Rosales, J.J., Jorba, A., Jorba-Cusco, M. (2021), "Families of Halo-like invariant tori around L2
in the Earth-Moon Bicircular Problem", Celestial Mechanics and Dynamical Astronomy 133:16,
doi:10.1007/s10569-021-10012-0. Filed in the private paper corpus as
rosales-jorba-jorba-cusco-2021-families-halo-like-invariant-tori-l2-earth-moon-bicircular-problem-cmda-133-16-doi-10.1007-s10569-021-10012-0.pdf

What the paper prints and how it was obtained:

- Table 1 (page 4): mu = 0.012150581623433, m_s = 328900.55, omega_s = 0.925195985518289,
  a_s = 388.811143023351121. Eq. (1) (page 3): Earth at x = mu, Sun at (a_s cos theta,
  -a_s sin theta), theta = omega_s t, that is, clockwise in the rotating frame.
- Section 2 (pages 4 to 6), eq. (2): H^eps = H_RTBP + eps H_S, eps multiplying the Sun's mass;
  "a pseudo-arclength continuation scheme with respect to an artificial parameter", multiple
  shooting with four sections, period T = 2 pi / omega_s. Fig. 2 and the text on page 6:
  "Moving from L2 to the right [towards the Moon], although the parameter eps becomes negative
  (which has no physical sense), it decreases until it hits a turning point, and then increases
  to become positive and reach eps = 1, that is, the BCP." Fig. 3: the orbit "revolves around
  the L2 twice in one period".
- Table 2 (page 7): eigenvalues of the monodromy matrix of that orbit at eps = 1:
  776607.104649077169597, 1.660211640235458 and 0.865694004478591 - 0.500573561636870 i
  (the other three are the reciprocals); stability saddle x saddle x center (page 7).

This module's frame has the Earth at -mu, which is the paper's frame rotated by pi; the Sun is
then at angle pi at t = 0 (``theta_sun0 = pi``) and still moves clockwise. The test repeats the
paper's route: pseudo-arclength continuation from L2 in eps through negative Sun mass, multiple
shooting with eight segments on ``bcr4bp.bcr4bp_stm_eom``, the eps-sensitivity of each segment
integrated alongside. The multipliers are taken, as on page 6 of the paper, from the spectrum of
the cyclic block matrix of the segment transition matrices (its eigenvalues are the 8th roots of
the multipliers), not from the product monodromy matrix: plain eigenvalues of the product resolve
the unit-modulus pair only to about 1e-10 because the product has norm 1e6.

Measured (2026-10-04): lambda_1 = 776607.1046491 (paper agreement 1e-13 relative), lambda_2 =
1.6602116402356 (1e-13), pair 0.86569400447866 - 0.50057356163672 i (1e-13 absolute). Spread of
the same computation over 4, 8 and 16 segments and integration tolerances 1e-13 and 1e-14: at
most 4e-13 relative for lambda_1, 2e-13 for lambda_2, 3e-13 absolute for the pair. The orbit is
planar; lambda_2 belongs to the out-of-plane block (an out-of-plane saddle) and the unit-modulus
pair to the in-plane block, together with lambda_1. The continuation reaches eps = 1 at x(0) =
1.101104 in this frame (paper Fig. 2: about -1.10 in its frame) after a turning point near
eps = -0.065.

The Sun's phase at t = 0 is not tested by this orbit: changing theta_sun0 by delta is the same
system shifted in time by delta / omega_s, and the continuation starts at L2, which does not
depend on time, so every phase leads to the same orbit and the same multipliers. The Sun's
sense of motion is tested: with the Sun moving counter-clockwise the route towards the Moon
starts into positive eps (the paper says negative) and folds away from eps = 1, and the route on
the far side of L2 reaches eps = 1 at a different orbit (measured lambda_1 = 1.0577e6, unit
pair at argument 0.317 rad against 0.524).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import cyclerfinder.core.bcr4bp as bcr4bp

FloatArray = NDArray[np.float64]

# Table 1, page 4.
_MU = 0.012150581623433
_M_SUN = 328900.55
_OMEGA_SUN = 0.925195985518289
_A_SUN = 388.811143023351121

# Table 2, page 7.
_LAMBDA_1 = 776607.104649077169597
_LAMBDA_2 = 1.660211640235458
_LAMBDA_3 = complex(0.865694004478591, -0.500573561636870)

_PERIOD = 2.0 * math.pi / _OMEGA_SUN
_N_SEG = 8
_DT = _PERIOD / _N_SEG
_RTOL = 1e-13
_EYE = np.eye(6).ravel()


def _system(eps: float, sense: float) -> bcr4bp.BCR4BPSystem:
    """sense = 1: the Sun clockwise as in the paper; sense = -1: reversed."""
    return bcr4bp.BCR4BPSystem(
        mu=_MU,
        mu_sun=eps * _M_SUN,
        a_sun_nondim=_A_SUN,
        omega_sun_nondim=sense * _OMEGA_SUN,
        theta_sun0=math.pi,
    )


def _rhs(
    t: float,
    y: FloatArray,
    system: bcr4bp.BCR4BPSystem,
    zero: bcr4bp.BCR4BPSystem,
    unit: bcr4bp.BCR4BPSystem,
) -> FloatArray:
    """State, transition matrix and the derivative of the state with respect to eps."""
    s = y[:6]
    out = bcr4bp.bcr4bp_stm_eom(t, np.concatenate([s, _EYE]), system)
    jac = out[6:].reshape(6, 6)
    # The Sun's term is linear in its mass, so d f / d eps is the Sun's term at eps = 1.
    d_eps = bcr4bp.bcr4bp_eom(t, s, unit) - bcr4bp.bcr4bp_eom(t, s, zero)
    return np.concatenate([out[:6], (jac @ y[6:42].reshape(6, 6)).ravel(), jac @ y[42:] + d_eps])


def _shoot(z: FloatArray, sense: float) -> tuple[FloatArray, FloatArray, list[FloatArray]]:
    """Multiple-shooting residual and Jacobian (with the eps column) for z = (nodes, eps)."""
    nodes = z[:-1].reshape(_N_SEG, 6)
    system = _system(float(z[-1]), sense)
    zero, unit = _system(0.0, sense), _system(1.0, sense)
    residual = np.zeros(6 * _N_SEG)
    jac = np.zeros((6 * _N_SEG, 6 * _N_SEG + 1))
    stms: list[FloatArray] = []
    for i in range(_N_SEG):
        sol = solve_ivp(
            _rhs,
            (i * _DT, (i + 1) * _DT),
            np.concatenate([nodes[i], _EYE, np.zeros(6)]),
            args=(system, zero, unit),
            method="DOP853",
            rtol=_RTOL,
            atol=_RTOL,
        )
        end = sol.y[:, -1]
        j = (i + 1) % _N_SEG
        stms.append(end[6:42].reshape(6, 6))
        residual[6 * i : 6 * i + 6] = end[:6] - nodes[j]
        jac[6 * i : 6 * i + 6, 6 * i : 6 * i + 6] = stms[-1]
        jac[6 * i : 6 * i + 6, 6 * j : 6 * j + 6] -= np.eye(6)
        jac[6 * i : 6 * i + 6, -1] = end[42:]
    return residual, jac, stms


def _x_l2() -> float:
    def d_omega(x: float) -> float:
        return (
            x
            - (1.0 - _MU) * (x + _MU) / abs(x + _MU) ** 3
            - _MU * (x - 1.0 + _MU) / abs(x - 1.0 + _MU) ** 3
        )

    return float(brentq(d_omega, 1.05, 1.3, xtol=1e-15))


def _l2_start(sense: float, toward_moon: bool) -> tuple[FloatArray, FloatArray]:
    """L2 at eps = 0 and the unit tangent of the branch, oriented by the direction in x."""
    z = np.concatenate([np.tile([_x_l2(), 0.0, 0.0, 0.0, 0.0, 0.0], _N_SEG), [0.0]])
    _, jac, _ = _shoot(z, sense)
    tau: FloatArray = np.linalg.svd(jac)[2][-1]
    # In this frame the Moon is at x = 1 - mu < x(L2): towards the Moon is decreasing x.
    if (tau[0] < 0.0) != toward_moon:
        tau = -tau
    return z, tau


@dataclass
class _Branch:
    nodes: FloatArray
    stms: list[FloatArray]
    path: list[tuple[float, float]]  # (eps, x at t = 0) at each continuation step


def _continue_to_bcp(sense: float, toward_moon: bool) -> _Branch:
    """Pseudo-arclength continuation from L2 (eps = 0) to the first point with eps = 1."""
    z, tau = _l2_start(sense, toward_moon)
    ds = 0.02
    z_prev = z
    path = [(0.0, float(z[0]))]
    for _ in range(40):
        predicted = z + ds * tau
        zc = predicted.copy()
        for _ in range(10):
            residual, jac, _ = _shoot(zc, sense)
            if float(np.linalg.norm(residual)) < 1e-10:
                break
            full = np.vstack([jac, tau])
            zc = zc - np.linalg.solve(full, np.concatenate([residual, [tau @ (zc - predicted)]]))
        else:
            ds *= 0.5
            if ds < 1e-4:
                raise AssertionError(f"continuation stalled at eps = {z[-1]}")
            continue
        new_tau = np.linalg.solve(np.vstack([jac, tau]), np.eye(len(z))[-1])
        tau, z_prev, z = new_tau / np.linalg.norm(new_tau), z, zc
        path.append((float(z[-1]), float(z[0])))
        if z[-1] >= 1.0:
            break
        ds = min(ds * 1.6, 0.3)
    else:
        raise AssertionError(f"eps = 1 not reached in 40 steps; path {path}")
    # Secant to eps = 1, then Newton with eps held at 1.
    z = z_prev + (1.0 - z_prev[-1]) / (z[-1] - z_prev[-1]) * (z - z_prev)
    z[-1] = 1.0
    for _ in range(10):
        residual, jac, stms = _shoot(z, sense)
        if float(np.linalg.norm(residual)) < 1e-12:
            return _Branch(z[:-1].reshape(_N_SEG, 6), stms, path)
        z[:-1] -= np.linalg.solve(jac[:, :-1], residual)
    raise AssertionError("Newton at eps = 1 did not converge")


@dataclass
class _Multipliers:
    big: complex
    big_inverse: complex
    saddle: complex
    saddle_inverse: complex
    unit: complex  # the member with negative imaginary part, as printed
    counts: tuple[int, ...]


def _cyclic_multipliers(stms: list[FloatArray]) -> _Multipliers:
    """Multipliers from the cyclic block matrix of the segment maps (each one n times)."""
    n = len(stms)
    cyc = np.zeros((6 * n, 6 * n))
    for i, stm in enumerate(stms):
        j = (i + 1) % n
        cyc[6 * j : 6 * j + 6, 6 * i : 6 * i + 6] = stm
    lam = np.linalg.eigvals(cyc) ** n
    real = np.abs(lam.imag) < 1e-6 * np.abs(lam)
    groups = [
        lam[real & (np.abs(lam) > 10.0)],
        lam[real & (np.abs(lam) < 1e-3)],
        lam[real & (np.abs(lam) > 1.2) & (np.abs(lam) < 10.0)],
        lam[real & (np.abs(lam) > 0.1) & (np.abs(lam) < 0.9)],
        lam[~real & (np.abs(np.abs(lam) - 1.0) < 1e-3) & (lam.imag < 0.0)],
    ]
    big, big_inverse, saddle, saddle_inverse, unit = (complex(np.mean(g)) for g in groups)
    return _Multipliers(
        big, big_inverse, saddle, saddle_inverse, unit, tuple(len(g) for g in groups)
    )


@pytest.fixture(scope="module")
def l2_replacement() -> _Branch:
    return _continue_to_bcp(sense=1.0, toward_moon=True)


def test_route_from_l2_passes_through_negative_sun_mass(l2_replacement: _Branch) -> None:
    """Fig. 2 and page 6: towards the Moon eps first goes negative, turns, and reaches 1."""
    eps = [p[0] for p in l2_replacement.path]
    x0 = [p[1] for p in l2_replacement.path]
    assert eps[1] < 0.0
    assert min(eps) < 0.0
    first_positive = next(i for i, e in enumerate(eps) if e > 0.0)
    assert min(eps[:first_positive]) < 0.0
    assert all(b < a for a, b in itertools.pairwise(x0))
    # Graph read from Fig. 2: eps = 1 at x about -1.10 in the paper's frame.
    assert abs(float(l2_replacement.nodes[0, 0]) - 1.10) < 0.01


def test_orbit_loops_twice_around_l2(l2_replacement: _Branch) -> None:
    """Fig. 3 caption: the orbit revolves around L2 twice in one period of the Sun."""
    x_l2 = _x_l2()
    sol = solve_ivp(
        bcr4bp.bcr4bp_eom,
        (0.0, _PERIOD),
        l2_replacement.nodes[0],
        args=(_system(1.0, 1.0),),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, _PERIOD, 801),
    )
    angles = np.unwrap(np.arctan2(sol.y[1], sol.y[0] - x_l2))
    assert round(abs(angles[-1] - angles[0]) / (2.0 * math.pi)) == 2
    assert float(np.max(np.abs(sol.y[2]))) < 1e-12  # planar


def test_table_2_multipliers(l2_replacement: _Branch) -> None:
    """Table 2, page 7, to the 1e-12 level that the computation's own spread (4e-13) supports."""
    m = _cyclic_multipliers(l2_replacement.stms)
    assert m.counts == (_N_SEG,) * 5  # saddle x saddle x center, each multiplier n times
    assert math.isclose(m.big.real, _LAMBDA_1, rel_tol=2e-12)
    assert math.isclose(m.saddle.real, _LAMBDA_2, rel_tol=2e-12)
    assert abs(m.unit - _LAMBDA_3) < 1e-12
    # Symplectic structure: reciprocal pairs and a pair on the unit circle.
    assert math.isclose((m.big * m.big_inverse).real, 1.0, abs_tol=1e-9)
    assert math.isclose((m.saddle * m.saddle_inverse).real, 1.0, abs_tol=1e-12)
    assert abs(abs(m.unit) - 1.0) < 1e-12


def test_saddle_lambda_2_is_out_of_plane(l2_replacement: _Branch) -> None:
    monodromy = np.eye(6)
    for stm in l2_replacement.stms:
        monodromy = stm @ monodromy
    vertical = np.linalg.eigvals(monodromy[np.ix_([2, 5], [2, 5])]).real
    assert math.isclose(float(np.max(vertical)), _LAMBDA_2, rel_tol=1e-11)


def test_reversed_sun_sense_starts_into_positive_sun_mass() -> None:
    """Page 6: from L2 towards the Moon eps becomes negative. Reversed, it becomes positive."""
    assert _l2_start(1.0, toward_moon=True)[1][-1] < 0.0
    assert _l2_start(-1.0, toward_moon=True)[1][-1] > 0.0


def test_reversed_sun_sense_gives_a_different_orbit() -> None:
    """With the Sun counter-clockwise, the far-side route reaches eps = 1 at another orbit.

    Measured: lambda_1 = 1.0577e6, unit pair 0.95018 - 0.31170 i (argument 0.317 rad against
    0.524 printed), x(0) = 1.18329.
    """
    branch = _continue_to_bcp(sense=-1.0, toward_moon=False)
    m = _cyclic_multipliers(branch.stms)
    assert abs(m.big.real / _LAMBDA_1 - 1.0) > 0.1
    assert abs(abs(np.angle(m.unit)) - abs(np.angle(_LAMBDA_3))) > 0.1
