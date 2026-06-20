"""Russell (2004) full-rev / half-rev re-initiating-return geometry.

Phase B of the Russell generic-return generator (#388): the analytic
re-initiating returns that close a cycler back onto the flyby body. The
equations are transcribed verbatim from the design spec Appendix A.3
(Russell §2.7.1-2.7.2, Eqs 2.12-2.25).

Frame (Appendix A.3): origin at the tip of ``v_B``, z-axis along ``v_B``,
y-axis opposed to the body's angular-momentum vector; the half-rev primed
frame has x' along ``r_B``. ``gamma`` is the body flight-path angle, which is
0 for a circular body orbit (``sin gamma = 0``, ``cos gamma = 1``), so the
circular-coplanar case the :class:`RussellModel` provides simplifies the
half-rev intersection considerably.

All quantities are in Russell's canonical units (mu_sun = 1, AU, TU).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from cyclerfinder.search.generic_return import RussellModel

_OMEGA_EPS = 1e-12


def _acos_clamped(x: float) -> float:
    """``math.acos`` with its argument clamped to ``[-1, 1]``.

    The turn-angle equations (Russell Eqs 3.2-3.8) combine products of sines and
    cosines that can drift a hair outside ``[-1, 1]`` from floating-point error;
    clamping keeps :func:`math.acos` in-domain without changing the geometry.
    """
    return math.acos(max(-1.0, min(1.0, x)))


def _open_omega(omega: float) -> float:
    """Keep a turn angle inside the open interval ``(0, pi)``.

    The min-max turn ``omega`` is a supremum: the degenerate endpoints (a 0 or
    pi turn) are unreachable by a finite flyby. When the geometry lands exactly
    on a boundary -- e.g. ``omega_c = pi`` for an exactly-radial incoming v_inf,
    where ``phi_GR = 0`` -- nudge the value just inside the open interval so the
    achievable-turn semantics hold. This is a numerical guard on the endpoints
    only, not a tolerance loosening of the formula.
    """
    return min(math.pi - _OMEGA_EPS, max(_OMEGA_EPS, omega))


def full_rev_feasible_vF(  # noqa: N802 -- Russell's v_F notation (spec API)
    model: RussellModel, body: str, n: int, big_m: int
) -> float:
    """Feasible full-rev post-flyby speed ``v_F`` (Russell Eq 2.13).

    Eq 2.13: ``v_F = sqrt(2 mu / r - mu * (N/M)^(2/3) / a_B)``, derived from the
    resonant transfer semi-major axis ``a_F = a_B * (N/M)^(2/3)`` substituted
    into vis-viva (Eq 2.12). Here ``N`` (``n``) is the spacecraft revolution
    count and ``M`` (``big_m``) the body revolution count over the same span;
    for a circular body the body distance ``r`` equals its semi-major axis
    ``a_B``.

    Raises
    ------
    ValueError
        If the argument of the square root is negative, i.e. the ``(N, M)``
        resonance is infeasible for this body.
    """
    mu = model.mu_sun
    a_b = model.sma_au(body)
    r = a_b  # circular body: distance from primary equals semi-major axis
    arg = 2.0 * mu / r - mu * (n / big_m) ** (2.0 / 3.0) / a_b
    if arg < 0.0:
        raise ValueError(f"infeasible full-rev resonance (N={n}, M={big_m}) for body {body!r}")
    return math.sqrt(arg)


def full_rev_circle_z(vF: float, vinf: float, vB: float) -> float:  # noqa: N803
    """Full-rev circle z-height ``z_F`` (Russell Eq 2.17).

    Eq 2.17: ``z_F = (v_F^2 - v_inf^2 - v_B^2) / (2 v_B)``. The full-rev circle
    is the slice of the v_inf sphere (Eq 2.15) at ``z = z_F``.
    """
    return (vF * vF - vinf * vinf - vB * vB) / (2.0 * vB)


def half_rev_components(model: RussellModel, body: str, a: float) -> tuple[float, float]:
    """Half-rev radial / transverse outbound velocity (Russell Eqs 2.18/2.19).

    With ``r1 = r2 = a_B`` (same-body half-rev to the flyby body):

    - Eq 2.18 radial: ``v_Hr^2 = mu * [2/(r1+r2) - 1/a]``
    - Eq 2.19 transverse: ``v_Htheta^2 = 2 mu r2 / (r1^2 + r1 r2)``
      (independent of ``a`` / TOF).

    The positive square root of Eq 2.18 is taken; below ``a_min`` (=
    ``(r1+r2)/2`` in the symmetric case) the argument goes negative -- the
    transfer cannot reach -- and ``v_Hr`` is returned as ``0.0``.

    Returns
    -------
    tuple[float, float]
        ``(v_Hr, v_Htheta)``.
    """
    mu = model.mu_sun
    r1 = model.sma_au(body)
    r2 = r1
    arg_r = mu * (2.0 / (r1 + r2) - 1.0 / a)
    v_hr = math.sqrt(arg_r) if arg_r > 0.0 else 0.0
    v_ht = math.sqrt(2.0 * mu * r2 / (r1 * r1 + r1 * r2))
    return v_hr, v_ht


def half_rev_intersection(
    model: RussellModel, body: str, a: float, vinf: float
) -> tuple[np.ndarray, np.ndarray]:
    """The two half-rev post-flyby v_inf tips (Russell Eqs 2.23-2.25 + K).

    Intersection of the half-rev circle (Eqs 2.20/2.21) with the v_inf sphere
    (Eq 2.22), transcribed verbatim from Appendix A.3 with ``gamma`` = body
    flight-path angle (0 for a circular orbit -> ``sin gamma = 0``,
    ``cos gamma = 1``):

    - K = ``2 v_Hr v_B sin g + v_B^2 (1 - 2 sin^2 g) - v_Htheta^2 - v_Hr^2 + v_inf^2``
    - x (2.23) = ``[v_Hr - v_B sin g] cos g + K sin g / (2 v_B cos g)``
    - y (2.24) = ``+/- sqrt{ v_inf^2 - [v_Hr - v_B sin g]^2 - K^2/(2 v_B cos g)^2 }``
    - z (2.25) = ``[v_Hr - v_B sin g] sin g - K / (2 v_B)``

    Returns the two ``(x, y, z)`` tips (the ``+`` and ``-`` y roots).

    Raises
    ------
    ValueError
        If the Eq 2.24 radicand is negative -- the circle and sphere do not
        intersect (no half-rev return at this ``a``, ``v_inf``).
    """
    v_b = model.body_circular_speed(body)
    v_hr, v_ht = half_rev_components(model, body, a)

    gamma = 0.0  # circular body orbit: flight-path angle is zero
    sin_g = math.sin(gamma)
    cos_g = math.cos(gamma)

    k = (
        2.0 * v_hr * v_b * sin_g
        + v_b * v_b * (1.0 - 2.0 * sin_g * sin_g)
        - v_ht * v_ht
        - v_hr * v_hr
        + vinf * vinf
    )

    base = v_hr - v_b * sin_g
    k_term = k / (2.0 * v_b * cos_g)

    x = base * cos_g + k_term * sin_g
    z = base * sin_g - k / (2.0 * v_b)

    radicand = vinf * vinf - base * base - k_term * k_term
    if radicand < 0.0:
        raise ValueError(
            f"no half-rev intersection (a={a}, vinf={vinf}) for body {body!r}: "
            "Eq 2.24 radicand is negative"
        )
    y = math.sqrt(radicand)

    tip_plus = np.array([x, y, z], dtype=np.float64)
    tip_minus = np.array([x, -y, z], dtype=np.float64)
    return tip_plus, tip_minus


def f_count(h_j: int) -> int:
    """Number of flybys ``f_j`` for ``h_j`` half-years (Russell Table 3.2).

    Even ``h_j``: ``f_j = h_j // 2 + 1``. Odd ``h_j``: ``f_j = 2 * INT(h_j/4 + 1)``.
    Reproduces the tabulated row ``h_j = 0..8 -> [1, 2, 2, 2, 3, 4, 4, 4, 5]``.
    """
    if h_j % 2 == 0:
        return h_j // 2 + 1
    return 2 * int(h_j / 4 + 1)


def omega_minimax(
    model: RussellModel,
    body: str,
    vinf: float,
    vinf_in_vec: np.ndarray,
    f_j: int,
) -> float:
    """Turn-angle min-max ``omega_minimax`` (Russell Eqs 3.2-3.8 + decision rule).

    Computes the smallest achievable maximum per-flyby turn angle for a run of
    ``f_j`` flybys that re-points the incoming v_inf vector ``vinf_in_vec`` (in
    Russell's canonical frame) onto the outgoing free-return direction. The
    body-velocity unit direction ``v_hat_e`` is taken from
    :meth:`RussellModel.body_state` at ``theta = 0`` (= ``+y``).

    Equations (Appendix A.1):

    - ``phi_FR = -asin(v_inf / (2 v_e))``                              (3.3)
    - ``phi_GR = pi/2 - acos((vinf_in . v_hat_e) / v_inf)``           (3.4)
    - ``omega_MIN = acos(cos phi_FR cos phi_GR + sin phi_FR sin phi_GR)`` (3.2)
    - ``omega_a = acos(cos^2 phi_FR cos lambda_a + sin^2 phi_FR)``,
      ``lambda_a = pi / (f_j - 2)``                                   (3.5)
    - ``omega_b``: solve Eq 3.6 for ``lambda`` in ``(0, pi/2)``, then
      ``omega_b = acos(cos phi_GR cos phi_FR cos lambda + sin phi_GR sin phi_FR)`` (3.7)
    - ``omega_c = pi - 2 |phi_GR|``                                    (3.8)

    Decision rule: ``f_j = 1 -> omega_c``; ``f_j = 2 -> acos(sin phi_GR sin phi_FR)``;
    ``f_j > 2 -> omega_MIN if omega_MIN >= omega_a else omega_b``.
    """
    v_e = model.body_circular_speed(body)
    _, v_body = model.body_state(body, 0.0)
    v_hat_e = v_body / np.linalg.norm(v_body)

    phi_FR = -math.asin(max(-1.0, min(1.0, vinf / (2.0 * v_e))))  # noqa: N806
    cos_in = float(np.dot(vinf_in_vec, v_hat_e)) / vinf
    phi_GR = math.pi / 2.0 - _acos_clamped(cos_in)  # noqa: N806

    if f_j == 1:
        return _open_omega(math.pi - 2.0 * abs(phi_GR))  # omega_c (3.8)

    if f_j == 2:
        return _open_omega(_acos_clamped(math.sin(phi_GR) * math.sin(phi_FR)))

    # f_j > 2.
    omega_MIN = _acos_clamped(  # noqa: N806
        math.cos(phi_FR) * math.cos(phi_GR) + math.sin(phi_FR) * math.sin(phi_GR)
    )
    lambda_a = math.pi / (f_j - 2)
    omega_a = _acos_clamped(math.cos(phi_FR) ** 2 * math.cos(lambda_a) + math.sin(phi_FR) ** 2)

    if omega_MIN >= omega_a:
        return _open_omega(omega_MIN)

    # omega_b branch: solve Eq 3.6 for lambda in (0, pi/2), then Eq 3.7.
    def eq_3_6(lam: float) -> float:
        lambda_b = (math.pi - 2.0 * lam) / (f_j - 2)
        omega = _acos_clamped(
            math.cos(phi_GR) * math.cos(phi_FR) * math.cos(lam)
            + math.sin(phi_GR) * math.sin(phi_FR)
        )
        omega_a_lam = _acos_clamped(
            math.cos(phi_FR) ** 2 * math.cos(lambda_b) + math.sin(phi_FR) ** 2
        )
        return omega - omega_a_lam

    lo, hi = 1e-9, math.pi / 2.0 - 1e-9
    if eq_3_6(lo) * eq_3_6(hi) > 0.0:
        # No sign change in the bracket: fall back to omega_a (Step 3 fallback).
        return _open_omega(omega_a)
    lam_root = brentq(eq_3_6, lo, hi)
    return _open_omega(
        _acos_clamped(
            math.cos(phi_GR) * math.cos(phi_FR) * math.cos(lam_root)
            + math.sin(phi_GR) * math.sin(phi_FR)
        )
    )  # omega_b (3.7)


def group_half_years(
    model: RussellModel,
    body: str,
    vinf: float,
    vinf_in_vec: np.ndarray,
    h: int,
    s: int,
) -> tuple[tuple[int, ...], float]:
    """Group ``h`` half-years across ``s`` generic returns (Russell §3.6).

    Tries the equal split ``h_j = h // s``. At that ``h_j`` it compares the
    single-flyby cap ``omega_c`` (= :func:`omega_minimax` with ``f_j = 1``) to the
    ``omega_minimax`` of the ``f_count(h_j)`` flybys in a group:

    - if ``omega_c >= omega_minimax`` -> equal split with the remainder on the
      LAST group: ``[h // s] * (s - 1) + [h // s + h % s]``;
    - otherwise pile everything on the first group: ``[h] + [0] * (s - 1)``.

    Returns ``(tuple_of_h_j, omega_max)`` where ``omega_max`` is the maximum
    ``omega_minimax`` over the groups (= the first group's, which dominates).
    """
    base = h // s
    rem = h % s

    omega_c = omega_minimax(model, body, vinf, vinf_in_vec, f_j=1)
    omega_mm = omega_minimax(model, body, vinf, vinf_in_vec, f_j=f_count(base))

    groups = [base] * (s - 1) + [base + rem] if omega_c >= omega_mm else [h] + [0] * (s - 1)

    omega_max = max(
        omega_minimax(model, body, vinf, vinf_in_vec, f_j=f_count(h_j)) for h_j in groups
    )
    return tuple(groups), omega_max
