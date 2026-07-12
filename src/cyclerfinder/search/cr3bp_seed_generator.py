"""Analytic CR3BP planar periodic-orbit seed generator at arbitrary mu (#435).

Builds a small-amplitude *linear* seed for a planar periodic orbit and refines
it with the existing fixed-Jacobi symmetric corrector
(:func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`).
Three families:

* :func:`lyapunov_seed` -- planar Lyapunov orbit about a collinear libration
  point (L1/L2), from the standard collinear linearization (Szebehely 1967;
  Koon, Lo, Marsden & Ross, *Dynamical Systems, the Three-Body Problem and
  Space Mission Design*, 2011, Ch. 2.5-2.7).
* :func:`lyapunov_seed_3d` -- vertical-Lyapunov / tulip-branch 3D orbit from
  the *first-order linearized* out-of-plane mode (NOT the halo branch).
* :func:`richardson_halo_seed` -- genuine halo-branch (Class I / Class II)
  orbit from Richardson's (1980) closed-form *third-order* analytic
  construction (#580; see :func:`richardson_halo_ic` for the pure,
  iteration-free analytic map and :func:`richardson_halo_coefficients` for
  the underlying Appendix I coefficients).
* :func:`dro_seed` -- distant retrograde orbit about the secondary, from a
  retrograde circular two-body guess clamped to stay bound.

Frame convention (matches :mod:`cyclerfinder.core.cr3bp`): nondimensional
rotating frame, primary at ``(-mu, 0, 0)``, secondary at ``(1-mu, 0, 0)``.

This feeds the #435 high-eccentricity Sun-planet ER3BP discovery, where a
converged CR3BP planar orbit at arbitrary mu is the continuation anchor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from cyclerfinder.core.cr3bp import CR3BPSystem, jacobi_constant
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    correct_general_periodic_3d,
)
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x

# Amplitude ladder (smallest-first) for the Lyapunov linear seed. The corrector's
# convergence basin is non-monotonic in amplitude and μ-dependent: at some μ the
# nominal 1e-3 lands in a non-convergent gap while smaller amplitudes (deeper in
# the linear regime, genuine small-amplitude Lyapunov, period ~ 2π/ω) converge
# cleanly. We try the caller's amplitude first, then walk this ladder.
_LYAPUNOV_AMPLITUDE_LADDER = (1.0e-4, 5.0e-5, 1.0e-5, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)

# Vertical (z) amplitude ladder (smallest-first) for the 3D vertical-Lyapunov /
# halo seed. The 3D corrector's convergence basin is non-monotonic in the
# vertical amplitude and μ-dependent; we try the caller's amplitude_z first,
# then walk this ladder (smallest-first, deepest in the linear regime).
_LYAPUNOV_Z_AMPLITUDE_LADDER = (0.02, 0.01, 5.0e-3, 0.03, 0.05, 0.08)


def _collinear_c2(x_l: float, mu: float) -> float:
    """In-plane linear coefficient c2 at a collinear point.

    ``c2 = (1-mu)/|x_L + mu|^3 + mu/|x_L - (1-mu)|^3`` -- the standard
    collinear-linearization coefficient (Koon-Lo-Marsden-Ross 2011, Eq. 2.5.x):
    the sum of the two primaries' gravity-gradient terms at ``x_L``.
    """
    d1 = abs(x_l + mu)
    d2 = abs(x_l - (1.0 - mu))
    return (1.0 - mu) / d1**3 + mu / d2**3


def lyapunov_seed(
    system: CR3BPSystem,
    *,
    point: str = "L1",
    amplitude: float = 1e-3,
) -> tuple[NDArray[np.float64], float]:
    """Converged planar Lyapunov orbit about a collinear point, at arbitrary mu.

    Construct the textbook small-amplitude linear seed at the collinear point
    ``point`` ("L1" or "L2") and refine it with the fixed-Jacobi symmetric
    corrector. ``amplitude`` is the nondimensional in-plane x-amplitude ``Ax``
    of the linear seed (small, e.g. ``1e-3``).

    The collinear linearization gives the in-plane frequency

        omega = sqrt((2 - c2 + sqrt(9 c2^2 - 8 c2)) / 2)

    and a planar Lyapunov initial condition (x-axis perpendicular crossing)

        x0  = x_L - Ax,   y0=z0=vx0=vz0 = 0,
        vy0 = -Ax * omega * tau,   tau = -(omega^2 + 1 + 2 c2) / (2 omega).

    Returns
    -------
    (state0, period)
        The corrected 6-state (perpendicular x-axis crossing) and the FULL
        nondimensional period ``T = 2 * t_half``.

    Raises
    ------
    ValueError
        If ``point`` is not a collinear point with an in-plane oscillatory mode
        (radicand ``9 c2^2 - 8 c2 < 0``), or if the corrector fails to converge.
    """
    mu = system.mu
    x_l = lagrange_collinear_x(mu, point)
    c2 = _collinear_c2(x_l, mu)
    radicand = 9.0 * c2 * c2 - 8.0 * c2
    if radicand < 0.0:
        raise ValueError(
            f"lyapunov_seed: no in-plane oscillatory mode at {point} "
            f"(9 c2^2 - 8 c2 = {radicand:.3e} < 0, c2={c2:.6f}, mu={mu:.3e})"
        )
    omega = math.sqrt((2.0 - c2 + math.sqrt(radicand)) / 2.0)
    tau = -(omega * omega + 1.0 + 2.0 * c2) / (2.0 * omega)
    period_guess = 2.0 * math.pi / omega

    # Try the requested amplitude first, then the smallest-first ladder. The
    # corrector's basin is non-monotonic in (amplitude, μ); walking the ladder
    # makes the generator reliable across Sun-Mercury/Mars/Pluto/Earth-Moon.
    tried: list[float] = []
    last_residual = math.inf
    last_iter = 0
    for ax in (float(amplitude), *_LYAPUNOV_AMPLITUDE_LADDER):
        if any(abs(ax - t) <= 1e-15 for t in tried):
            continue
        tried.append(ax)
        x0 = x_l - ax
        vy0 = -ax * omega * tau
        seed_state = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0], dtype=np.float64)
        jacobi = jacobi_constant(seed_state, mu)
        orbit = correct_symmetric_fixed_jacobi(
            system,
            x0_guess=x0,
            jacobi=jacobi,
            period_guess=period_guess,
            ydot0_sign=math.copysign(1.0, vy0) if vy0 != 0.0 else 1.0,
            half_crossings=1,
        )
        if orbit.converged:
            state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
            return state0, orbit.period
        last_residual = orbit.crossing_residual
        last_iter = orbit.n_iter
    raise ValueError(
        f"lyapunov_seed: corrector did not converge at {point} for any amplitude "
        f"in {tried} (last residual={last_residual:.3e}, n_iter={last_iter}, mu={mu:.3e})"
    )


def lyapunov_seed_3d(
    system: CR3BPSystem,
    *,
    point: str = "L1",
    amplitude_z: float = 0.02,
    amplitude_x: float = 1e-3,
) -> tuple[NDArray[np.float64], float]:
    """Converged genuinely-3D vertical-Lyapunov / halo orbit at a collinear point.

    The existing :func:`lyapunov_seed` produces only PLANAR orbits. This builds
    the textbook small-amplitude *vertical* (out-of-plane) linear seed at the
    collinear point ``point`` ("L1" or "L2") and refines it with the full-3D
    symmetric (tulip / NRHO style) corrector.

    The out-of-plane linear mode at a collinear point decouples from the in-plane
    motion as ``z'' + c2 z = 0`` (Koon, Lo, Marsden & Ross 2011 §2.5), so the
    vertical frequency is ``omega_z = sqrt(c2)`` with
    ``c2 = _collinear_c2(x_L, mu)``. The seed IC is the perpendicular x-z-plane
    crossing

        x0  = x_L,   y0 = 0,   z0 = amplitude_z,
        vx0 = 0,     vy0 = <small in-plane>,   vz0 = 0,

    with the in-plane ``vy0`` taken from the planar-Lyapunov linear expression at
    ``amplitude_x`` (small; the corrector refines it). Refinement uses
    :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
    with the symmetric tulip free vars ``(z0, ydot0, T)`` and the
    perpendicular-half-period residual ``(y, xdot, zdot)`` at ``T/2``.

    Parameters
    ----------
    system :
        CR3BP system (only ``system.mu`` is used).
    point :
        Collinear point ("L1" or "L2") with an oscillatory vertical mode.
    amplitude_z :
        Nondimensional out-of-plane amplitude of the linear seed. If the
        corrector is fragile here, a small smallest-first ladder is walked.
    amplitude_x :
        Nondimensional in-plane x-amplitude used only to seed ``vy0``.

    Returns
    -------
    (state0, period)
        The corrected genuinely-3D 6-state (``|z0| > 0``) and the FULL
        nondimensional period.

    Raises
    ------
    ValueError
        If no amplitude in the ladder converges (message names point, mu, and
        the last residual).
    """
    mu = system.mu
    x_l = lagrange_collinear_x(mu, point)
    c2 = _collinear_c2(x_l, mu)
    if c2 <= 0.0:
        raise ValueError(
            f"lyapunov_seed_3d: no vertical oscillatory mode at {point} "
            f"(c2={c2:.6e} <= 0, mu={mu:.3e})"
        )
    omega_z = math.sqrt(c2)
    period_guess = 2.0 * math.pi / omega_z

    # In-plane vy0 seed from the planar-Lyapunov linear expression (the corrector
    # refines it). Only used if the in-plane radicand is non-negative; otherwise
    # start vy0 small.
    radicand = 9.0 * c2 * c2 - 8.0 * c2
    if radicand >= 0.0:
        omega_xy = math.sqrt((2.0 - c2 + math.sqrt(radicand)) / 2.0)
        tau = -(omega_xy * omega_xy + 1.0 + 2.0 * c2) / (2.0 * omega_xy)
        vy0_seed = -float(amplitude_x) * omega_xy * tau
    else:
        vy0_seed = float(amplitude_x)

    tried: list[float] = []
    last_residual = math.inf
    last_iter = 0
    for az in (float(amplitude_z), *_LYAPUNOV_Z_AMPLITUDE_LADDER):
        if any(abs(az - t) <= 1e-15 for t in tried):
            continue
        tried.append(az)
        seed_state = np.array([x_l, 0.0, az, 0.0, vy0_seed, 0.0], dtype=np.float64)
        orbit = correct_general_periodic_3d(
            system,
            state0_guess=seed_state,
            period_guess=period_guess,
            free_vars=FREE_VARS_SYMMETRIC_TULIP,
            residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
            is_half_period_residual=True,
        )
        if orbit.converged and not orbit.degenerate_planar:
            return orbit.state0, orbit.T_TU
        last_residual = orbit.independent_closure_residual
        last_iter = orbit.n_iter
    raise ValueError(
        f"lyapunov_seed_3d: corrector did not converge to a non-planar orbit at "
        f"{point} for any amplitude_z in {tried} "
        f"(last residual={last_residual:.3e}, n_iter={last_iter}, mu={mu:.3e})"
    )


def _richardson_gamma_l(mu: float, point: str) -> float:
    """Solve Richardson's collinear-point quintic for ``gamma_L`` (#580).

    ``gamma_L`` is the libration point's normalized distance to its nearer
    primary (Richardson 1980, p. 244: ``gamma_L = r/a1``, with ``r`` the
    L-point-to-nearer-primary distance and ``a1`` the primary separation).
    Richardson's Appendix I gives every coefficient as a closed-form function
    of ``gamma_L`` (via ``c2, c3, c4``), so ``gamma_L`` is solved directly
    from the standard collinear-point quintic rather than converted from
    ``x_L`` -- this also covers L3, which
    :func:`cyclerfinder.search.reachable_representatives.lagrange_collinear_x`
    does not (that helper only supports L1/L2). The quintics are the
    textbook Euler quintic for each collinear point (e.g. Szebehely 1967;
    same equations underlying the JPL/Fortran-Astrodynamics-Toolkit
    ``compute_libration_points`` routine).
    """
    if point == "L1":
        eq = lambda g: (  # noqa: E731
            g**5 - (3.0 - mu) * g**4 + (3.0 - 2.0 * mu) * g**3 - mu * g**2 + 2.0 * mu * g - mu
        )
        lo, hi = 1.0e-8, 1.0 - 1.0e-8
    elif point == "L2":
        eq = lambda g: (  # noqa: E731
            g**5 + (3.0 - mu) * g**4 + (3.0 - 2.0 * mu) * g**3 - mu * g**2 - 2.0 * mu * g - mu
        )
        lo, hi = 1.0e-8, 1.0 - 1.0e-8
    elif point == "L3":
        eq = lambda g: (  # noqa: E731
            g**5
            + (2.0 + mu) * g**4
            + (1.0 + 2.0 * mu) * g**3
            - (1.0 - mu) * g**2
            - 2.0 * (1.0 - mu) * g
            - (1.0 - mu)
        )
        lo, hi = 0.2, 1.5
    else:
        raise ValueError(f"point must be 'L1', 'L2', or 'L3', got {point!r}")
    return float(brentq(eq, lo, hi))


def _richardson_x_l(mu: float, point: str, gamma_l: float) -> float:
    """Collinear point x-coordinate (standard barycentric frame) from ``gamma_L``."""
    if point == "L1":
        return (1.0 - mu) - gamma_l
    if point == "L2":
        return (1.0 - mu) + gamma_l
    return -(mu + gamma_l)  # L3


def _richardson_c_n(n: int, mu: float, point: str, gamma_l: float) -> float:
    """Richardson's ``c_n`` Legendre-expansion coefficient (Eq. 8a-b).

    ``c_n`` depends only on ``mu`` and ``gamma_L`` and feeds every downstream
    Appendix I coefficient. Verified numerically (see #580 golden test)
    against Richardson (1980) Table I, p. 253, for all three collinear
    points -- including the sign convention for L2 (``c_n`` alternates sign
    with ``n``, so ``c3`` is NEGATIVE at L2; the digest's Table I
    transcription dropped that minus sign, corrected 2026-07-12 against a
    fresh direct read of the source PDF, see the #580 OUTSTANDING resolution
    note).
    """
    if point == "L1":
        return float(
            (mu + (-1) ** n * (1.0 - mu) * (gamma_l / (1.0 - gamma_l)) ** (n + 1)) / gamma_l**3
        )
    if point == "L2":
        return float(
            ((-1) ** n * (mu + (1.0 - mu) * (gamma_l / (1.0 + gamma_l)) ** (n + 1))) / gamma_l**3
        )
    if point == "L3":
        return (1.0 - mu + mu * (gamma_l / (1.0 + gamma_l)) ** (n + 1)) / gamma_l**3
    raise ValueError(f"point must be 'L1', 'L2', or 'L3', got {point!r}")


@dataclass(frozen=True)
class RichardsonHaloCoefficients:
    """Richardson's (1980) Appendix I closed-form halo-construction constants.

    All fields are pure functions of ``(mu, point)`` -- no amplitude or
    iteration involved. Field names and equation references match Richardson
    (1980), *Celestial Mechanics* 22, 241-253, Appendix I (p. 250-252)
    exactly; see ``docs/notes/2026-07-12-digest-richardson-1980-collinear-
    halo-analytic.md`` section 4 for the sourced Table I golden values these
    are validated against (#580).

    Attributes
    ----------
    point, mu, gamma_l :
        Inputs echoed back (``gamma_l`` per :func:`_richardson_gamma_l`).
    lam :
        ``lambda`` -- the in-plane linearized frequency, the positive root of
        the characteristic quartic ``lambda^4 + (c2-2) lambda^2 -
        (2 c2 + 1)(c2 - 1) = 0`` (Richardson's Eq. 4, restated in Appendix I).
    k :
        y/x in-plane amplitude ratio, ``k = 2 lambda / (lambda^2 + 1 - c2)``
        (Appendix I).
    delta :
        ``Delta = lambda^2 - c2`` (Eq. 12), the out-of-plane frequency
        correction that Eq. 18's amplitude constraint balances.
    c2, c3, c4 :
        Legendre-expansion coefficients (Eq. 8a-b), functions of
        ``mu, gamma_L`` alone.
    s1, s2 :
        Frequency-correction coefficients (Appendix I; feed
        ``omega = 1 + s1 Ax^2 + s2 Az^2``, Eq. 17).
    l1, l2 :
        Amplitude-constraint coefficients (Appendix I; feed Eq. 18,
        ``l1 Ax^2 + l2 Az^2 + Delta = 0``).
    a1, a2, d1, d2 :
        Intermediate Appendix I constants feeding ``s1, s2, l1, l2``.
    a21..a32, b21..b32, d21, d31, d32 :
        The Eq. 20a-c series coefficients (Appendix I closed-form
        recursions in ``c2, c3, c4, k, lambda``).
    """

    point: str
    mu: float
    gamma_l: float
    lam: float
    k: float
    delta: float
    c2: float
    c3: float
    c4: float
    s1: float
    s2: float
    l1: float
    l2: float
    a1: float
    a2: float
    d1: float
    d2: float
    a21: float
    a22: float
    a23: float
    a24: float
    a31: float
    a32: float
    b21: float
    b22: float
    b31: float
    b32: float
    d21: float
    d31: float
    d32: float


def richardson_halo_coefficients(mu: float, point: str) -> RichardsonHaloCoefficients:
    """Compute Richardson's (1980) Appendix I halo coefficients at ``(mu, point)``.

    Pure, iteration-free closed-form evaluation (only ``gamma_L`` requires a
    1D root-find, Eq. 8a-b's own quintic -- not a shooting/Newton corrector).
    ``point`` is ``"L1"``, ``"L2"``, or ``"L3"``. See
    :class:`RichardsonHaloCoefficients` for field definitions and the #580
    golden test (``tests/search/test_cr3bp_seed_generator.py``) for the
    Table I cross-check at Sun-Earth ``mu=3.04036e-6``.

    Raises
    ------
    ValueError
        If ``point`` is not ``"L1"``/``"L2"``/``"L3"``, or if the in-plane
        radicand ``9 c2^2 - 8 c2`` is negative (no oscillatory in-plane mode
        -- does not occur for any physical collinear point but guarded for
        safety).
    """
    gamma_l = _richardson_gamma_l(mu, point)
    c2 = _richardson_c_n(2, mu, point, gamma_l)
    c3 = _richardson_c_n(3, mu, point, gamma_l)
    c4 = _richardson_c_n(4, mu, point, gamma_l)

    radicand = 9.0 * c2 * c2 - 8.0 * c2
    if radicand < 0.0:
        raise ValueError(
            f"richardson_halo_coefficients: no in-plane oscillatory mode at {point} "
            f"(9 c2^2 - 8 c2 = {radicand:.3e} < 0, c2={c2:.6f}, mu={mu:.3e})"
        )
    lam = math.sqrt((2.0 - c2 + math.sqrt(radicand)) / 2.0)
    k = 2.0 * lam / (lam * lam + 1.0 - c2)
    delta = lam * lam - c2

    d1 = 16.0 * lam**4 + 4.0 * lam**2 * (c2 - 2.0) - 2.0 * c2**2 + c2 + 1.0
    d2 = 81.0 * lam**4 + 9.0 * lam**2 * (c2 - 2.0) - 2.0 * c2**2 + c2 + 1.0
    d3 = 2.0 * lam * (lam * (1.0 + k * k) - 2.0 * k)

    a21 = 3.0 * c3 * (k * k - 2.0) / 4.0 / (1.0 + 2.0 * c2)
    a23 = -3.0 * lam * c3 * (3.0 * k**3 * lam - 6.0 * k * (k - lam) + 4.0) / 4.0 / k / d1
    b21 = -3.0 * c3 * lam * (3.0 * lam * k - 4.0) / 2.0 / d1
    s1 = (
        1.5 * c3 * (2.0 * a21 * (k * k - 2.0) - a23 * (k * k + 2.0) - 2.0 * k * b21)
        - 0.375 * c4 * (3.0 * k**4 - 8.0 * k * k + 8.0)
    ) / d3

    a22 = 3.0 * c3 / 4.0 / (1.0 + 2.0 * c2)
    a24 = -3.0 * c3 * lam * (2.0 + 3.0 * lam * k) / 4.0 / k / d1
    b22 = 3.0 * lam * c3 / d1
    d21 = -c3 / 2.0 / lam**2
    s2 = (
        1.5 * c3 * (2.0 * a22 * (k * k - 2.0) + a24 * (k * k + 2.0) + 2.0 * k * b22 + 5.0 * d21)
        + 0.375 * c4 * (12.0 - k * k)
    ) / d3

    a1 = -1.5 * c3 * (2.0 * a21 + a23 + 5.0 * d21) - 0.375 * c4 * (12.0 - k * k)
    a2 = 1.5 * c3 * (a24 - 2.0 * a22) + 1.125 * c4
    l1 = 2.0 * s1 * lam * lam + a1
    l2 = 2.0 * s2 * lam * lam + a2

    a31 = (
        -9.0 * lam * (c3 * (k * a23 - b21) + k * c4 * (1.0 + 0.25 * k * k)) / d2
        + (9.0 * lam * lam + 1.0 - c2)
        * (3.0 * c3 * (2.0 * a23 - k * b21) + c4 * (2.0 + 3.0 * k * k))
        / 2.0
        / d2
    )
    a32 = (
        -9.0 * lam * (4.0 * c3 * (k * a24 - b22) + k * c4) / 4.0 / d2
        - 3.0 * (9.0 * lam * lam + 1.0 - c2) * (c3 * (k * b22 + d21 - 2.0 * a24) - c4) / 2.0 / d2
    )
    b31 = (
        3.0 * lam * (3.0 * c3 * (k * b21 - 2.0 * a23) - c4 * (2.0 + 3.0 * k * k))
        + (9.0 * lam * lam + 1.0 + 2.0 * c2)
        * (12.0 * c3 * (k * a23 - b21) + 3.0 * k * c4 * (4.0 + k * k))
        / 8.0
    ) / d2
    b32 = (
        3.0 * lam * (3.0 * c3 * (k * b22 + d21 - 2.0 * a24) - 3.0 * c4)
        + (9.0 * lam * lam + 1.0 + 2.0 * c2) * (12.0 * c3 * (k * a24 - b22) + 3.0 * c4 * k) / 8.0
    ) / d2
    d31 = 3.0 * (4.0 * c3 * a24 + c4) / 64.0 / lam**2
    d32 = 3.0 * (4.0 * c3 * (a23 - d21) + c4 * (4.0 + k * k)) / 64.0 / lam**2

    return RichardsonHaloCoefficients(
        point=point,
        mu=mu,
        gamma_l=gamma_l,
        lam=lam,
        k=k,
        delta=delta,
        c2=c2,
        c3=c3,
        c4=c4,
        s1=s1,
        s2=s2,
        l1=l1,
        l2=l2,
        a1=a1,
        a2=a2,
        d1=d1,
        d2=d2,
        a21=a21,
        a22=a22,
        a23=a23,
        a24=a24,
        a31=a31,
        a32=a32,
        b21=b21,
        b22=b22,
        b31=b31,
        b32=b32,
        d21=d21,
        d31=d31,
        d32=d32,
    )


def richardson_halo_ic(
    mu: float,
    point: str,
    amplitude_z: float,
    *,
    branch: str = "I",
    phi: float = 0.0,
) -> tuple[NDArray[np.float64], float]:
    """Richardson's (1980) third-order analytic halo IC -- no correction (#580).

    Evaluates Eq. 20a-c (and their ``tau1``-derivatives) at ``tau1 = phi``
    and converts from Richardson's libration-point-centered, ``gamma_L``-
    normalized units to the standard CR3BP barycentric nondimensional
    rotating frame (matches :mod:`cyclerfinder.core.cr3bp`: primary at
    ``(-mu, 0, 0)``, secondary at ``(1-mu, 0, 0)``). This is the genuinely
    *closed-form, iteration-free* map ``(mu, point, A_z) -> IC`` -- the only
    root-find involved is ``gamma_L``'s quintic (Eq. 8a-b), not a
    shooting/Newton periodicity correction. Callers wanting a CONVERGED
    periodic orbit should feed the result to
    :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
    directly, or use :func:`richardson_halo_seed`.

    Parameters
    ----------
    mu :
        CR3BP mass ratio.
    point :
        Collinear point, ``"L1"``, ``"L2"``, or ``"L3"``.
    amplitude_z :
        Nondimensional CR3BP out-of-plane amplitude (same convention as
        :func:`lyapunov_seed_3d`'s ``amplitude_z``, i.e. in the standard
        barycentric distance unit, NOT Richardson's own ``gamma_L``-scaled
        ``Az``). Internally ``Az_richardson = amplitude_z / gamma_L``
        (Richardson normalizes distance so the L-point-to-nearer-primary
        separation is 1; the appendix's ``A_z`` is a multiple of that unit).
        Sign selects the north/south mirror member.
    branch :
        ``"I"`` (n=1) or ``"II"`` (n=3) -- Richardson's two solution
        branches (Eq. 19/21); ``delta_n = 2 - n`` multiplies the z-series in
        Eq. 20c.
    phi :
        Phase angle (Eq. 20's ``tau1 = lambda*tau + phi``, evaluated at
        ``tau=0``). At ``phi=0`` every ``sin(n*phi)`` term vanishes
        identically, so ``y0 = xdot0 = zdot0 = 0`` EXACTLY -- the IC is
        already a perpendicular x-z-plane crossing, the exact form
        :data:`cyclerfinder.search.cr3bp_general_periodic_3d.FREE_VARS_SYMMETRIC_TULIP`
        assumes. This is why ``phi=0`` (the default, and what
        :func:`richardson_halo_seed` always uses) is the natural choice for
        feeding the existing tulip/NRHO-style corrector.

    Returns
    -------
    (state0, period_guess)
        The UNCORRECTED analytic 6-state and Richardson's own third-order
        period estimate ``2 pi / (lambda * omega)`` (nondim TU), where
        ``omega = 1 + s1 Ax^2 + s2 Az^2`` is the frequency correction
        (Eq. 17).

    Raises
    ------
    ValueError
        If ``branch`` is not ``"I"``/``"II"``, or if ``amplitude_z`` is
        infeasible: Eq. 18's amplitude constraint
        ``Ax = sqrt((-Delta - l2 Az^2) / l1)`` has no real solution.
    """
    if branch == "I":
        n_switch = 1
    elif branch == "II":
        n_switch = 3
    else:
        raise ValueError(f"branch must be 'I' or 'II', got {branch!r}")
    delta_n = float(2 - n_switch)

    c = richardson_halo_coefficients(mu, point)
    az = float(amplitude_z) / c.gamma_l
    az2 = az * az

    term = (-c.delta - c.l2 * az2) / c.l1
    if term < 0.0:
        raise ValueError(
            f"richardson_halo_ic: infeasible amplitude_z={amplitude_z:.3e} at {point} "
            f"(l1={c.l1:.3e}, l2={c.l2:.3e}, delta={c.delta:.3e}, "
            f"Ax^2=(-delta-l2*Az^2)/l1={term:.3e} < 0 => no real Ax; Eq. 18)"
        )
    ax = math.sqrt(term)
    ax2 = ax * ax
    ax3 = ax2 * ax
    omega = 1.0 + c.s1 * ax2 + c.s2 * az2

    tau1 = float(phi)
    ct1, st1 = math.cos(tau1), math.sin(tau1)
    c2t1, s2t1 = math.cos(2.0 * tau1), math.sin(2.0 * tau1)
    c3t1, s3t1 = math.cos(3.0 * tau1), math.sin(3.0 * tau1)

    # Eq. 20a-c, evaluated at tau1 (Richardson's own libration-point-centered,
    # gamma_L-normalized units).
    x = (
        c.a21 * ax2
        + c.a22 * az2
        - ax * ct1
        + (c.a23 * ax2 - c.a24 * az2) * c2t1
        + (c.a31 * ax3 - c.a32 * ax * az2) * c3t1
    )
    y = (
        c.k * ax * st1
        + (c.b21 * ax2 - c.b22 * az2) * s2t1
        + (c.b31 * ax3 - c.b32 * ax * az2) * s3t1
    )
    z = (
        delta_n * az * ct1
        + delta_n * c.d21 * ax * az * (c2t1 - 3.0)
        + delta_n * (c.d32 * az * ax2 - c.d31 * ax3) * c3t1
    )
    # d/dtau1 of the above (Richardson units; converted to nondim TU below).
    vx = (
        ax * st1
        - 2.0 * (c.a23 * ax2 - c.a24 * az2) * s2t1
        - 3.0 * (c.a31 * ax3 - c.a32 * ax * az2) * s3t1
    )
    vy = (
        c.k * ax * ct1
        + 2.0 * (c.b21 * ax2 - c.b22 * az2) * c2t1
        + 3.0 * (c.b31 * ax3 - c.b32 * ax * az2) * c3t1
    )
    vz = (
        -delta_n * az * st1
        - 2.0 * delta_n * c.d21 * ax * az * s2t1
        - 3.0 * delta_n * (c.d32 * az * ax2 - c.d31 * ax3) * s3t1
    )

    # Convert Richardson (gamma_L-normalized, L-point-centered) -> standard
    # CR3BP barycentric nondim frame: positions scale by gamma_L, velocities
    # by gamma_L*lambda*omega (chain rule, tau1 = lambda*omega*t + const in
    # nondim TU), and x is offset by the libration point's own x-coordinate.
    x_l = _richardson_x_l(mu, point, c.gamma_l)
    scale_v = c.gamma_l * c.lam * omega
    state0 = np.array(
        [
            x * c.gamma_l + x_l,
            y * c.gamma_l,
            z * c.gamma_l,
            vx * scale_v,
            vy * scale_v,
            vz * scale_v,
        ],
        dtype=np.float64,
    )
    period_guess = 2.0 * math.pi / (c.lam * omega)
    return state0, period_guess


def richardson_halo_seed(
    system: CR3BPSystem,
    *,
    point: str = "L1",
    amplitude_z: float = 0.02,
    branch: str = "I",
) -> tuple[NDArray[np.float64], float]:
    """Converged genuine halo-branch orbit via Richardson's (1980) analytic seed.

    Fills the halo-branch (Class I / Class II) seed-generator gap: unlike
    :func:`lyapunov_seed_3d` (first-order linearized vertical mode, tulip/
    vertical-Lyapunov branch), this seeds from Richardson's closed-form
    THIRD-order construction (:func:`richardson_halo_ic`, evaluated at
    ``phi=0`` so the IC is already an exact perpendicular x-z-plane crossing)
    and refines it with the same
    :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
    symmetric tulip corrector -- with ``x0`` FIXED at Richardson's own
    analytic x-offset (NOT at ``x_L``), which is a materially better guess
    for genuinely halo-shaped (large-amplitude) 3D orbits than the
    linearized seed's ``x0 = x_L``.

    Parameters
    ----------
    system :
        CR3BP system (only ``system.mu`` is used).
    point :
        Collinear point (``"L1"``, ``"L2"``, or ``"L3"``) with an oscillatory
        vertical mode.
    amplitude_z :
        Nondimensional out-of-plane amplitude fed to :func:`richardson_halo_ic`
        (same convention as :func:`lyapunov_seed_3d`).
    branch :
        ``"I"`` or ``"II"`` -- Richardson's two solution branches.

    Returns
    -------
    (state0, period)
        The corrected genuinely-3D 6-state (``|z0| > 0``) and the FULL
        nondimensional period.

    Raises
    ------
    ValueError
        If ``amplitude_z`` is infeasible for the analytic construction
        (see :func:`richardson_halo_ic`), or if the corrector fails to
        converge to a non-planar orbit.
    """
    mu = system.mu
    state0_guess, period_guess = richardson_halo_ic(mu, point, amplitude_z, branch=branch, phi=0.0)
    orbit = correct_general_periodic_3d(
        system,
        state0_guess=state0_guess,
        period_guess=period_guess,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    if orbit.converged and not orbit.degenerate_planar:
        return orbit.state0, orbit.T_TU
    raise ValueError(
        f"richardson_halo_seed: corrector did not converge to a non-planar orbit at "
        f"{point} branch={branch} for amplitude_z={amplitude_z:.3e} "
        f"(residual={orbit.independent_closure_residual:.3e}, n_iter={orbit.n_iter}, mu={mu:.3e})"
    )


def dro_seed(
    system: CR3BPSystem,
    *,
    amplitude: float = 5e-2,
) -> tuple[NDArray[np.float64], float]:
    """Converged distant retrograde orbit (DRO) about the secondary.

    Retrograde co-orbital seed near the secondary (at ``(1-mu, 0, 0)``): a
    perpendicular x-axis crossing exterior to the secondary at radius
    ``amplitude``, with a retrograde circular two-body guess
    ``v ~ sqrt(mu / amplitude)`` (clamped to stay bound), refined with the
    fixed-Jacobi symmetric corrector.

    Returns
    -------
    (state0, period)
        The corrected 6-state and FULL nondimensional period.

    Raises
    ------
    ValueError
        If the corrector fails to converge.
    """
    mu = system.mu
    a = float(amplitude)
    x0 = (1.0 - mu) + a
    # Retrograde circular two-body speed about the secondary, clamped so the
    # guess stays bound (escape ~ sqrt(2 mu / a)); retrograde => negative vy.
    v_circ = math.sqrt(mu / a)
    vy0 = -v_circ
    seed_state = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0], dtype=np.float64)
    jacobi = jacobi_constant(seed_state, mu)
    period_guess = 2.0 * math.pi * a / v_circ

    orbit = correct_symmetric_fixed_jacobi(
        system,
        x0_guess=x0,
        jacobi=jacobi,
        period_guess=period_guess,
        ydot0_sign=-1.0,
        half_crossings=1,
    )
    if not orbit.converged:
        raise ValueError(
            f"dro_seed: corrector did not converge "
            f"(residual={orbit.crossing_residual:.3e}, n_iter={orbit.n_iter})"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    return state0, orbit.period
