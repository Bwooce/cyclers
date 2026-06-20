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

from cyclerfinder.search.generic_return import RussellModel


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
