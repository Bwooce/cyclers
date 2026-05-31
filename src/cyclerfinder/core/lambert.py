"""Universal-variable Lambert two-point boundary-value solver.

Given two heliocentric position vectors ``r1, r2`` and a time of flight ``tof``,
finds the heliocentric velocity vectors ``v1`` (at ``r1``) and ``v2`` (at
``r2``) such that a Keplerian conic connects them in exactly ``tof`` seconds.
The single-revolution case is fully implemented in M1; the multi-revolution
branches return an empty list and are M4's responsibility (see plan §3.2.1).

Algorithm
---------
The implementation follows Vallado, *Fundamentals of Astrodynamics and
Applications*, 4th ed., Algorithm 5.2 (LAMBERTUNIV), with the Stumpff functions
shared via :mod:`cyclerfinder.core._stumpff`. The inner loop is numpy-only —
no SciPy or external solver — to keep production-path imports minimal.

Cross-check
-----------
:func:`lambert_crosscheck` re-solves a leg with ``lamberthub.izzo2015`` and
``lamberthub.gooding1990`` (dev dependency only) and reports the worst
per-component velocity disagreement. This is the M1 gate.

References
----------
* Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th ed.,
  Microcosm Press, 2013, Algorithm 5.2.
* Bate, R. R., Mueller, D. D., and White, J. E., *Fundamentals of
  Astrodynamics*, Dover, 1971, §5.5.

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §3.2.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, sqrt
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core._stumpff import stumpff_c, stumpff_s
from cyclerfinder.core.constants import MU_SUN_KM3_S2

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

_NEWTON_TOL_REL: float = 1.0e-12
_NEWTON_TOL_DZ: float = 1.0e-12
_NEWTON_MAX_ITER: int = 60


class LambertError(Exception):
    """Base class for Lambert solver errors."""


class LambertGeometryError(LambertError):
    """Inputs describe a degenerate geometry the single-rev UV form can't handle.

    The two classic cases are a transfer angle of exactly ``0`` (the two
    position vectors coincide on a ray from the focus) and exactly ``180`` deg
    (``r2 = -k * r1`` for some ``k > 0``); for the latter the in-plane direction
    of the transfer is ambiguous and only a multi-rev / out-of-plane treatment
    resolves it. Multi-rev support is the M4 deliverable.
    """


class LambertConvergenceError(LambertError):
    """Newton iteration on the universal variable z failed to converge.

    Attributes
    ----------
    z:
        Last iterate of z.
    residual:
        Last value of ``t(z) - tof`` in seconds.
    """

    def __init__(self, z: float, residual: float) -> None:
        super().__init__(
            f"Lambert universal-variable Newton failed to converge: "
            f"z={z:.6e}, residual={residual:.6e} s"
        )
        self.z = z
        self.residual = residual


@dataclass(frozen=True)
class LambertSolution:
    """A single Lambert solution.

    Attributes
    ----------
    n_revs:
        Revolution count. ``0`` for the single-revolution branch; positive
        integers for multi-rev (M4).
    branch:
        ``"single"`` when ``n_revs == 0``; ``"low"`` or ``"high"`` for the two
        multi-rev branches once those land. The label is informational and may
        be refined when multi-rev support is added; downstream code that needs
        to filter solutions should switch on ``n_revs`` and ``branch`` together.
    v1:
        Heliocentric velocity at ``r1``, shape ``(3,)`` float64, km/s.
    v2:
        Heliocentric velocity at ``r2``, shape ``(3,)`` float64, km/s.
    """

    n_revs: int
    branch: str
    v1: Vec3
    v2: Vec3


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _t_of_z(z: float, a_coef: float, r1_n: float, r2_n: float, mu: float) -> tuple[float, float]:
    """Return ``(t, y)`` given current ``z``. Used inside the Newton loop.

    Vallado eqs. (Algorithm 5.2 main loop): y(z) and t(z). The argument
    ``a_coef`` is Vallado's ``A``, renamed to satisfy lower-case naming
    conventions while keeping the algebra recognisable.
    """
    c = stumpff_c(z)
    s = stumpff_s(z)
    y = r1_n + r2_n + a_coef * (z * s - 1.0) / sqrt(c)
    t = ((y / c) ** 1.5) * s / sqrt(mu) + a_coef * sqrt(y / mu)
    return t, y


def _dt_dz(z: float, y: float, a_coef: float, mu: float) -> float:
    """Analytic derivative dt/dz used for Newton steps.

    Vallado eq. (7-15) form, valid away from ``z = 0``. Near ``z = 0`` we
    fall back to a Maclaurin truncation to avoid the ``1/z`` divergence.
    """
    c = stumpff_c(z)
    s = stumpff_s(z)
    if abs(z) > 1.0e-5:
        # Standard form.
        sqrt_y_over_c = sqrt(y / c)
        y_over_c_15 = float((y / c) ** 1.5)
        term1 = (
            y_over_c_15 * ((1.0 / (2.0 * z)) * (c - 1.5 * s / c) + 0.75 * (s * s) / c) / sqrt(mu)
        )
        term2 = (a_coef / 8.0) * (3.0 * s * sqrt_y_over_c / c + a_coef * sqrt(c / y)) / sqrt(mu)
        return term1 + term2
    # Near z=0 fallback (Vallado): expansion about z=0.
    # dt/dz ≈ (sqrt(2)/40) y0^1.5 + (A/8) * (sqrt(y0) + A * sqrt(1/(2*y0)))  / sqrt(mu)
    y_15 = float(y**1.5)
    return (
        (sqrt(2.0) / 40.0) * y_15 + (a_coef / 8.0) * (sqrt(y) + a_coef * sqrt(1.0 / (2.0 * y)))
    ) / sqrt(mu)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def lambert(
    r1: Vec3,
    r2: Vec3,
    tof: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
    max_revs: int = 0,
) -> list[LambertSolution]:
    """Solve the Lambert problem for the leg ``r1 -> r2`` in time ``tof``.

    Parameters
    ----------
    r1, r2:
        Position vectors of the endpoints in the same inertial frame, ``(3,)``
        float64, km. Must be non-zero and non-collinear (a 180-degree transfer
        is singular for the single-rev universal-variable form; see
        :class:`LambertGeometryError`).
    tof:
        Time of flight, seconds. Must be strictly positive.
    mu:
        Central-body gravitational parameter, km^3/s^2. Defaults to the
        heliocentric :data:`cyclerfinder.core.constants.MU_SUN_KM3_S2`.
    prograde:
        ``True`` selects the short-way prograde transfer (transfer angle
        measured counter-clockwise about ``+z``); ``False`` selects the
        retrograde transfer.
    max_revs:
        Maximum number of full revolutions to consider. **M1 returns at most
        one solution (``n_revs == 0``); multi-revolution branches land in
        M4.** The parameter is accepted now so the interface is stable
        through M4.

    Returns
    -------
    list[LambertSolution]
        In M1, a length-1 list ``[LambertSolution(n_revs=0, branch="single",
        v1=..., v2=...)]``. The empty list is never returned in M1 success
        paths; errors raise instead.

    Raises
    ------
    ValueError
        On non-positive ``tof`` or zero-magnitude endpoint vectors.
    LambertGeometryError
        On 0-degree or 180-degree transfer geometry.
    LambertConvergenceError
        If Newton iteration fails within :data:`_NEWTON_MAX_ITER` steps.
    """
    if tof <= 0.0:
        raise ValueError(f"tof must be positive, got {tof}")
    # max_revs is reserved for M4; it must be a non-negative int but is
    # otherwise ignored in M1 (documented above). Validate the sign only.
    if max_revs < 0:
        raise ValueError(f"max_revs must be non-negative, got {max_revs}")

    r1_arr = np.asarray(r1, dtype=np.float64)
    r2_arr = np.asarray(r2, dtype=np.float64)

    r1_n = float(np.linalg.norm(r1_arr))
    r2_n = float(np.linalg.norm(r2_arr))
    if r1_n == 0.0 or r2_n == 0.0:
        raise ValueError("r1, r2 must be non-zero position vectors")

    # Transfer angle delta_nu in [0, 2*pi).
    cos_dnu = float(np.dot(r1_arr, r2_arr) / (r1_n * r2_n))
    cos_dnu = max(min(cos_dnu, 1.0), -1.0)  # clamp for numerical safety
    dnu = acos(cos_dnu)

    cross_z = float(r1_arr[0] * r2_arr[1] - r1_arr[1] * r2_arr[0])

    # Short/long-way selection per Vallado §7.6:
    # prograde + cross_z > 0  -> dnu unchanged (short way)
    # prograde + cross_z < 0  -> dnu = 2*pi - dnu (long way)
    # retrograde reverses both.
    if prograde:
        if cross_z < 0.0:
            dnu = 2.0 * np.pi - dnu
    else:
        if cross_z > 0.0:
            dnu = 2.0 * np.pi - dnu

    # Singular geometries.
    if abs(dnu) < 1.0e-12 or abs(dnu - np.pi) < 1.0e-9 or abs(dnu - 2.0 * np.pi) < 1.0e-12:
        raise LambertGeometryError(
            "Lambert single-rev universal-variable form is singular for 0- or "
            "180-degree transfer angles; multi-rev support arrives in M4."
        )

    sin_dnu = float(np.sin(dnu))
    # `a_coef` is Vallado's ``A`` (renamed to keep lower-case naming).
    a_coef = sin_dnu * sqrt(r1_n * r2_n / (1.0 - cos_dnu))
    if a_coef == 0.0:
        # Defensive; the dnu checks above should have caught this.
        raise LambertGeometryError("Lambert geometry yields A = 0; transfer is singular.")

    # ------------------------------------------------------------------
    # Bracket z so that y(z) >= 0, then Newton-iterate.
    # ------------------------------------------------------------------

    # Single-rev upper bound: the universal-variable t(z) goes to infinity at
    # z = (2*pi)^2 ~ 39.48 (this is the boundary at which the transfer becomes
    # a full revolution). Stepping beyond it pushes Newton onto the wrong
    # branch entirely. We keep z strictly below the singularity.
    z_high_single_rev = 4.0 * np.pi * np.pi
    z = 0.0

    # Bisection bootstrap: find a z with y(z) > 0 and walk Newton from there.
    # If y(z=0) < 0 we widen z upward until positive.
    def _y_only(z_in: float) -> float:
        c = stumpff_c(z_in)
        s = stumpff_s(z_in)
        return r1_n + r2_n + a_coef * (z_in * s - 1.0) / sqrt(c)

    y0 = _y_only(z)
    if y0 < 0.0:
        # Walk up until y becomes positive.
        step = 0.1
        z_try = z + step
        while _y_only(z_try) < 0.0 and z_try < z_high_single_rev:
            z_try += step
        z = z_try

    # Establish a sign-changing bracket [z_lo, z_hi] for f(z) = t(z) - tof,
    # within the single-rev valid range. t(z) is monotone-increasing inside
    # the valid range, so a bracket exists whenever the problem is well-posed.
    def _f_of_z(z_in: float) -> float:
        t_z, _y_unused = _t_of_z(z_in, a_coef, r1_n, r2_n, mu)
        return t_z - tof

    z_lo = -50.0  # well into hyperbolic territory
    # Step z_lo down further if t(z_lo) is still > tof (extremely short ToF).
    while True:
        try:
            f_lo = _f_of_z(z_lo)
        except ValueError:
            z_lo *= 0.5
            continue
        if f_lo < 0.0:
            break
        z_lo *= 2.0
        if z_lo < -1.0e6:
            # Pathological geometry / tof; let Newton raise below.
            break

    # Approach z_high from below; t -> +inf as z -> z_high_single_rev.
    eps_high = 1.0e-6
    z_hi = z_high_single_rev - eps_high

    # Initial Newton seed: midpoint of bracket if available, else current z.
    if _y_only(z) < 0.0 or not (z_lo <= z <= z_hi):
        z = 0.5 * (z_lo + z_hi)

    residual: float = 0.0
    converged = False
    for _it in range(_NEWTON_MAX_ITER):
        try:
            t_z, y = _t_of_z(z, a_coef, r1_n, r2_n, mu)
        except ValueError:
            # y(z) negative or sqrt(C(z)) issue: bisect into bracket.
            z = 0.5 * (z_lo + z_hi)
            continue

        residual = t_z - tof
        if abs(residual) / tof < _NEWTON_TOL_REL:
            converged = True
            break

        # Tighten bracket using sign of residual (monotone t(z)).
        if residual < 0.0:
            z_lo = z
        else:
            z_hi = z

        dt = _dt_dz(z, y, a_coef, mu)
        z_next = z - residual / dt if dt != 0.0 else z

        # Reject the Newton step if it leaves the bracket or pushes y < 0;
        # fall back to bisection in that case. This is Vallado's standard
        # safeguard for the universal-variable Lambert iteration.
        if not (z_lo < z_next < z_hi) or _y_only(z_next) < 0.0:
            z_next = 0.5 * (z_lo + z_hi)

        if abs(z_next - z) < _NEWTON_TOL_DZ:
            z = z_next
            converged = True
            break
        z = z_next

    if not converged:
        raise LambertConvergenceError(z, residual)

    # Final Lagrange coefficients.
    c = stumpff_c(z)
    s = stumpff_s(z)
    y = r1_n + r2_n + a_coef * (z * s - 1.0) / sqrt(c)
    f = 1.0 - y / r1_n
    g = a_coef * sqrt(y / mu)
    g_dot = 1.0 - y / r2_n

    v1 = (r2_arr - f * r1_arr) / g
    v2 = (g_dot * r2_arr - r1_arr) / g

    return [
        LambertSolution(
            n_revs=0,
            branch="single",
            v1=v1.astype(np.float64, copy=False),
            v2=v2.astype(np.float64, copy=False),
        )
    ]


def lambert_crosscheck(
    r1: Vec3,
    r2: Vec3,
    tof: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
) -> dict[str, Any]:
    """Compare the in-house Lambert with ``lamberthub.izzo2015`` and ``gooding1990``.

    Used by the M1 gate test and (later) by M3's V0/V1 closure verification.
    ``lamberthub`` is a development dependency; importing it inside this
    function keeps it off the production import path.

    Returns
    -------
    dict
        Keys: ``"mine"`` (the in-house :class:`LambertSolution`), ``"izzo"``
        (a ``(v1, v2)`` tuple from ``lamberthub.izzo2015``), ``"gooding"``
        (likewise from ``lamberthub.gooding1990``), and ``"max_diff_mps"``
        (worst velocity component disagreement across the four pairings,
        expressed in metres per second).
    """
    # Local import: lamberthub is a dev dependency.
    from lamberthub import gooding1990, izzo2015  # type: ignore[import-untyped]

    mine = lambert(r1, r2, tof, mu=mu, prograde=prograde)[0]

    r1_arr = np.asarray(r1, dtype=np.float64)
    r2_arr = np.asarray(r2, dtype=np.float64)

    v1_izzo, v2_izzo = izzo2015(mu, r1_arr, r2_arr, tof, M=0, prograde=prograde)
    v1_g, v2_g = gooding1990(mu, r1_arr, r2_arr, tof, M=0, prograde=prograde)

    diffs_km_s = [
        float(np.linalg.norm(mine.v1 - v1_izzo)),
        float(np.linalg.norm(mine.v2 - v2_izzo)),
        float(np.linalg.norm(mine.v1 - v1_g)),
        float(np.linalg.norm(mine.v2 - v2_g)),
    ]
    max_diff_mps = 1000.0 * max(diffs_km_s)

    return {
        "mine": mine,
        "izzo": (v1_izzo, v2_izzo),
        "gooding": (v1_g, v2_g),
        "max_diff_mps": max_diff_mps,
    }
