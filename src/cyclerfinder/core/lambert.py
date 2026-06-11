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
from math import acos, copysign, expm1, log1p, sqrt
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core._stumpff import stumpff_c, stumpff_s
from cyclerfinder.core.constants import MU_SUN_KM3_S2

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

_ROOT_TOL_REL: float = 1.0e-12
_ROOT_TOL_DZ: float = 1.0e-12
_ROOT_MAX_ITER: int = 60
_BRACKET_MAX_WIDEN_ITERS: int = 100
"""Bound on the bracket-finder widen loop. Well-posed geometries find the
z_lo bracket in < 10 widen iterations; this cap exists to fail fast on
pathological multi-start grid perturbations rather than loop forever (the
optimiser's exception handler drops the start on LambertConvergenceError)."""


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


def _y_of_z(z: float, a_coef: float, r1_n: float, r2_n: float) -> float:
    """Vallado's ``y(z)``. Standalone so the Newton safeguard can reject a
    step that would drive ``y`` negative without raising."""
    c = stumpff_c(z)
    s = stumpff_s(z)
    return r1_n + r2_n + a_coef * (z * s - 1.0) / sqrt(c)


def _dt_dz(z: float, y: float, a_coef: float, mu: float) -> float:
    """Analytic derivative dt/dz used for single-rev Newton steps.

    Vallado eq. (7-15) form, valid away from ``z = 0``. Near ``z = 0`` we
    fall back to a Maclaurin truncation to avoid the ``1/z`` divergence. This
    derivative is only valid on the single-revolution ``z`` branch; the
    multi-rev branches use the derivative-free :func:`_solve_uv_branch`.
    """
    c = stumpff_c(z)
    s = stumpff_s(z)
    if abs(z) > 1.0e-5:
        # Standard form. Differentiating sqrt(mu)*t = (y/C)^1.5 * S + A*sqrt(y)
        # with the Stumpff identities dC/dz = (1 - z*S - 2*C)/(2z) and
        # dS/dz = (C - 3*S)/(2z) gives
        #   sqrt(mu) dt/dz = (y/C)^1.5 [ (1/(2z))(C - 3S/(2C)) + (3/4) S^2/C ]
        #                  + (A/8) [ 3 S sqrt(y) / C + A sqrt(C/y) ].
        # (A historical defect computed the first piece of term2 as
        # 3*S*sqrt(y)/C^1.5 — one spurious sqrt(C) — which was inconsistent
        # with the z->0 Maclaurin limit below and made Newton oscillate on
        # long-way transfers; see task #205.)
        y_over_c_15 = float((y / c) ** 1.5)
        term1 = (
            y_over_c_15 * ((1.0 / (2.0 * z)) * (c - 1.5 * s / c) + 0.75 * (s * s) / c) / sqrt(mu)
        )
        term2 = (a_coef / 8.0) * (3.0 * s * sqrt(y) / c + a_coef * sqrt(c / y)) / sqrt(mu)
        return term1 + term2
    # Near z=0 fallback (Vallado): expansion about z=0.
    y_15 = float(y**1.5)
    return (
        (sqrt(2.0) / 40.0) * y_15 + (a_coef / 8.0) * (sqrt(y) + a_coef * sqrt(1.0 / (2.0 * y)))
    ) / sqrt(mu)


def _solve_single_rev_newton(
    z_seed: float,
    z_lo: float,
    z_hi: float,
    a_coef: float,
    r1_n: float,
    r2_n: float,
    tof: float,
    mu: float,
) -> float:
    """Safeguarded Newton root of ``f(z)=t(z)-tof`` on the single-rev branch.

    Reproduces Vallado's standard universal-variable iteration: a Newton step
    using the analytic :func:`_dt_dz`, with bracket tightening on the residual
    sign and a bisection fallback whenever a step would leave the bracket or
    drive ``y`` negative. The analytic derivative is valid throughout the
    single-rev ``z`` domain, so Newton converges quadratically here; the
    multi-rev branches use the derivative-free :func:`_solve_uv_branch`
    instead.

    Raises :class:`LambertConvergenceError` if the iteration does not converge
    within :data:`_ROOT_MAX_ITER` steps.
    """
    z = z_seed
    if _y_of_z(z, a_coef, r1_n, r2_n) < 0.0 or not (z_lo <= z <= z_hi):
        z = 0.5 * (z_lo + z_hi)

    residual: float = 0.0
    for _it in range(_ROOT_MAX_ITER):
        try:
            t_z, y = _t_of_z(z, a_coef, r1_n, r2_n, mu)
        except ValueError:
            # y(z) negative or sqrt(C(z)) issue: bisect into bracket.
            z = 0.5 * (z_lo + z_hi)
            continue

        residual = t_z - tof
        if abs(residual) / tof < _ROOT_TOL_REL:
            return z

        # Tighten bracket using sign of residual (monotone t(z)).
        if residual < 0.0:
            z_lo = z
        else:
            z_hi = z

        dt = _dt_dz(z, y, a_coef, mu)
        z_next = z - residual / dt if dt != 0.0 else z

        # Reject the Newton step if it leaves the bracket or pushes y < 0;
        # fall back to bisection in that case (Vallado's standard safeguard).
        if not (z_lo < z_next < z_hi) or _y_of_z(z_next, a_coef, r1_n, r2_n) < 0.0:
            z_next = 0.5 * (z_lo + z_hi)

        if abs(z_next - z) < _ROOT_TOL_DZ:
            return z_next
        z = z_next

    raise LambertConvergenceError(z, residual)


def _find_single_rev_bracket(
    a_coef: float,
    r1_n: float,
    r2_n: float,
    tof: float,
    mu: float,
    z_high_single_rev: float,
) -> tuple[float, float, int]:
    """Establish a sign-changing bracket ``[z_lo, z_hi]`` for ``f(z)=t(z)-tof``.

    ``t(z)`` is monotone-increasing across the single-revolution valid domain,
    whose lower edge is the *floor* ``z_floor`` where ``y(z) -> 0`` (below it
    ``y < 0`` and the conic is undefined). At the floor ``t -> 0``, so
    ``f(z_floor) = -tof < 0`` for any positive ``tof``: the floor is always a
    valid lower bracket endpoint. The upper endpoint is just below
    ``z_high_single_rev`` where ``t -> +inf``.

    The previous implementation fixed ``z_lo = -50`` and widened by repeated
    doubling, halving toward ``z = 0`` whenever ``y(z_lo) < 0`` raised. For
    geometries whose floor sits near ``z = 0`` (small radius ratio / small
    transfer angle) that halve-toward-zero walk oscillated and exhausted
    :data:`_BRACKET_MAX_WIDEN_ITERS`, raising on short-but-feasible transfers
    (see ``tests/core/test_lambert.py`` task #56 cases). Anchoring ``z_lo`` at
    the floor via bisection brackets every well-posed geometry in
    ``O(log2(range / tol))`` steps and fails fast (clean raise) only when the
    geometry has no valid ``z = 0`` interior at all.

    Returns ``(z_lo, z_hi, widen_iters)`` where ``widen_iters`` is the number of
    floor-search iterations consumed (for diagnostics / regression on the walk
    cost).
    """
    eps_high = 1.0e-6
    z_hi = z_high_single_rev - eps_high

    # Stumpff C(z)/S(z) grow like cosh/sinh(sqrt(-z)) on the hyperbolic side and
    # overflow float64 near z ~ -6.5e5; stay well clear so y/t stay evaluable.
    z_floor_limit = -5.0e5

    def _y_only(z_in: float) -> float:
        c = stumpff_c(z_in)
        s = stumpff_s(z_in)
        return r1_n + r2_n + a_coef * (z_in * s - 1.0) / sqrt(c)

    def _f_of_z(z_in: float) -> float:
        t_z, _y_unused = _t_of_z(z_in, a_coef, r1_n, r2_n, mu)
        return t_z - tof

    # The interior anchor: z = 0 is valid (y > 0) for every well-posed single-
    # rev geometry. If it is not, the transfer has no valid universal-variable
    # interior here -- raise so the optimiser drops the start cleanly.
    if _y_only(0.0) <= 0.0:
        raise LambertConvergenceError(0.0, float("nan"))

    # t(z) is monotone-increasing across the valid domain and t(0) may be above
    # or below tof. If f(0) < 0 already, z = 0 is itself a lower bracket end.
    if _f_of_z(0.0) < 0.0:
        return 0.0, z_hi, 0

    widen_iters = 0

    # Expand a negative anchor geometrically. Stop early the moment we reach a
    # *valid* z with f(z) < 0 (that is a usable z_lo -- no need to find the
    # exact floor). Track the most-negative *invalid* (y < 0) point so we can
    # bisect to the floor if the whole expansion stayed valid-but-f>=0.
    z_anchor = -1.0
    z_valid = 0.0  # most-negative point known valid (y > 0); z = 0 to start
    z_invalid: float | None = None
    while z_anchor > z_floor_limit:
        widen_iters += 1
        if _y_only(z_anchor) <= 0.0:
            z_invalid = z_anchor
            break
        if _f_of_z(z_anchor) < 0.0:
            return z_anchor, z_hi, widen_iters
        z_valid = z_anchor  # valid, but f still >= 0; remember it
        z_anchor *= 2.0

    if z_invalid is None:
        # Reached the stumpff-safe limit with y still > 0 and f >= 0 the whole
        # way: t cannot be driven below tof in the safe domain -> infeasible.
        raise LambertConvergenceError(z_anchor, float("nan"))

    # Bisect [z_invalid (y<0), z_valid (y>0)] toward the floor, returning as soon
    # as a valid point with f < 0 is found (the floor itself has t -> 0 < tof).
    for _ in range(_BRACKET_MAX_WIDEN_ITERS):
        widen_iters += 1
        z_mid = 0.5 * (z_invalid + z_valid)
        if z_mid in (z_valid, z_invalid):
            break  # converged to floating-point resolution
        if _y_only(z_mid) > 0.0:
            z_valid = z_mid
            if _f_of_z(z_mid) < 0.0:
                return z_mid, z_hi, widen_iters
        else:
            z_invalid = z_mid

    return z_valid, z_hi, widen_iters


def _bracket_diagnostics(
    r1: Vec3,
    r2: Vec3,
    tof: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
) -> dict[str, float]:
    """Expose the single-rev bracket-finder's internals for tests.

    Returns ``{"z_lo", "z_hi", "widen_iters"}``. Not part of the public API;
    used by ``tests/core/test_lambert.py`` to assert the widen walk stays well
    under :data:`_BRACKET_MAX_WIDEN_ITERS` on geometries that previously
    approached the cap.
    """
    r1_arr = np.asarray(r1, dtype=np.float64)
    r2_arr = np.asarray(r2, dtype=np.float64)
    r1_n = float(np.linalg.norm(r1_arr))
    r2_n = float(np.linalg.norm(r2_arr))
    cos_dnu = max(min(float(np.dot(r1_arr, r2_arr) / (r1_n * r2_n)), 1.0), -1.0)
    dnu = acos(cos_dnu)
    cross_z = float(r1_arr[0] * r2_arr[1] - r1_arr[1] * r2_arr[0])
    if prograde:
        if cross_z < 0.0:
            dnu = 2.0 * np.pi - dnu
    else:
        if cross_z > 0.0:
            dnu = 2.0 * np.pi - dnu
    sin_dnu = float(np.sin(dnu))
    a_coef = sin_dnu * sqrt(r1_n * r2_n / (1.0 - cos_dnu))
    z_lo, z_hi, widen_iters = _find_single_rev_bracket(
        a_coef, r1_n, r2_n, tof, mu, 4.0 * np.pi * np.pi
    )
    return {"z_lo": z_lo, "z_hi": z_hi, "widen_iters": float(widen_iters)}


def _min_time_of_revolution(
    n: int, a_coef: float, r1_n: float, r2_n: float, mu: float
) -> tuple[float, float]:
    """Return ``(z_min, t_min)`` for revolution ``n >= 1``.

    On the open interval ``z in ((2*pi*n)**2, (2*pi*(n+1))**2)`` the
    universal-variable time-of-flight ``t(z)`` is convex with a single
    interior minimum. We locate it with a bounded golden-section search
    (numpy-only; no SciPy on the production path). ``t_min`` is the shortest
    time of flight achievable with exactly ``n`` full revolutions; a
    requested ``tof <= t_min`` means revolution ``n`` is infeasible.
    """
    lo = (2.0 * np.pi * n) ** 2
    hi = (2.0 * np.pi * (n + 1)) ** 2
    # Stay strictly inside the open interval: t -> +inf at both endpoints.
    span = hi - lo
    a = lo + 1.0e-6 * span
    b = hi - 1.0e-6 * span

    def _t(z_in: float) -> float:
        t_z, _ = _t_of_z(z_in, a_coef, r1_n, r2_n, mu)
        return t_z

    inv_phi = (sqrt(5.0) - 1.0) / 2.0  # 1/golden ratio
    c = b - inv_phi * (b - a)
    d = a + inv_phi * (b - a)
    fc = _t(c)
    fd = _t(d)
    for _ in range(200):
        if abs(b - a) < 1.0e-9 * (abs(a) + abs(b)):
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = _t(c)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = _t(d)
    z_min = 0.5 * (a + b)
    t_min, _ = _t_of_z(z_min, a_coef, r1_n, r2_n, mu)
    return z_min, t_min


def _solve_uv_branch(
    z_lo: float,
    z_hi: float,
    a_coef: float,
    r1_n: float,
    r2_n: float,
    tof: float,
    mu: float,
) -> float:
    """Derivative-free root of ``f(z)=t(z)-tof`` on a sign-changing bracket.

    Uses the Illinois method (false-position with the anti-stall
    modification): guaranteed to converge on any bracket where ``f(z_lo)``
    and ``f(z_hi)`` have opposite signs, with no dependence on the analytic
    ``dt/dz``. This matters because the closed-form derivative is only valid
    in the single-rev ``z`` domain; on the multi-rev branches (either side of
    a revolution's time minimum) it has the wrong magnitude and sign, so a
    Newton iteration crawls or diverges there. The bracket is sign-monotone
    on each branch (the whole single-rev range, and each side of the
    multi-rev minimum), so Illinois is both robust and fast here.

    The residual is iterated in log-compressed form,
    ``g = sign(t - tof) * log1p(|t - tof| / tof)``, which has exactly the same
    root and sign as the raw residual but compresses the huge values near the
    revolution-boundary time singularities (``t(z) ~ 1e29 s`` at a bracket
    endpoint maps to ``g ~ 48``). With the raw residual, false position against
    such an endpoint takes steps scaled by ``~tof / 1e29`` and the Illinois
    halving needs ``log2(1e29 / tof) ~ 70`` iterations just to deflate the
    stale endpoint — past :data:`_ROOT_MAX_ITER`, silently dropping a feasible
    high branch (task #205 defect B). Near the root ``g ~ (t - tof)/tof``, so
    the convergence test ``|g| < _ROOT_TOL_REL`` is the same relative-residual
    criterion as before: precision is unchanged, only the stall is removed.

    Raises :class:`LambertConvergenceError` if the supplied endpoints do not
    bracket a root or the iteration fails to converge.
    """

    def _g(z_in: float) -> float:
        t_z, _y = _t_of_z(z_in, a_coef, r1_n, r2_n, mu)
        resid = t_z - tof
        return copysign(log1p(abs(resid) / tof), resid)

    def _residual_seconds(g_val: float) -> float:
        """Invert the log compression for error reporting (seconds)."""
        return copysign(expm1(abs(g_val)) * tof, g_val)

    a, b = z_lo, z_hi
    try:
        fa = _g(a)
        fb = _g(b)
    except ValueError as exc:
        raise LambertConvergenceError(0.5 * (a + b), float("nan")) from exc

    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if (fa > 0.0) == (fb > 0.0):
        # Endpoints share a sign: no guaranteed root in this bracket.
        raise LambertConvergenceError(
            0.5 * (a + b), min(abs(_residual_seconds(fa)), abs(_residual_seconds(fb)))
        )

    c = 0.5 * (a + b)
    for _it in range(_ROOT_MAX_ITER):
        # False-position estimate; fall back to bisection if degenerate or
        # out of the current bracket.
        if fb != fa:
            c = b - fb * (b - a) / (fb - fa)
        lo, hi = (a, b) if a < b else (b, a)
        if not (lo < c < hi):
            c = 0.5 * (a + b)
        try:
            fc = _g(c)
        except ValueError:
            c = 0.5 * (a + b)
            fc = _g(c)

        # |g| ~ |t - tof| / tof near the root: same relative criterion as the
        # raw-residual form.
        if abs(fc) < _ROOT_TOL_REL or abs(b - a) < _ROOT_TOL_DZ:
            return c

        if (fc > 0.0) == (fb > 0.0):
            # Root lies between a and c: discard b, but down-weight the now
            # stale endpoint a (the Illinois modification) to avoid stalling.
            b, fb = c, fc
            fa *= 0.5
        else:
            # Root lies between c and b: a <- b, b <- c.
            a, fa = b, fb
            b, fb = c, fc

    raise LambertConvergenceError(c, _residual_seconds(fc))


def _velocities_from_z(
    z: float,
    a_coef: float,
    r1_arr: Vec3,
    r2_arr: Vec3,
    r1_n: float,
    r2_n: float,
    mu: float,
) -> tuple[Vec3, Vec3]:
    """Lagrange-coefficient velocities at the converged universal variable z."""
    c = stumpff_c(z)
    s = stumpff_s(z)
    y = r1_n + r2_n + a_coef * (z * s - 1.0) / sqrt(c)
    f = 1.0 - y / r1_n
    g = a_coef * sqrt(y / mu)
    g_dot = 1.0 - y / r2_n
    v1 = (r2_arr - f * r1_arr) / g
    v2 = (g_dot * r2_arr - r1_arr) / g
    return v1.astype(np.float64, copy=False), v2.astype(np.float64, copy=False)


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
        Maximum number of full revolutions to consider. ``0`` returns only
        the single-revolution (direct) transfer. For each ``n`` in
        ``[1, max_revs]`` two further solutions may be appended -- the ``low``
        and ``high`` branches -- when the time of flight admits revolution
        ``n``; infeasible revolutions are skipped silently.

    Returns
    -------
    list[LambertSolution]
        Always begins with the single-revolution solution
        ``LambertSolution(n_revs=0, branch="single", ...)``, followed by any
        feasible multi-revolution branches in ascending ``n``, ``low`` before
        ``high``. The empty list is never returned on success paths; errors
        raise instead.

    Raises
    ------
    ValueError
        On non-positive ``tof`` or zero-magnitude endpoint vectors.
    LambertGeometryError
        On 0-degree or 180-degree transfer geometry.
    LambertConvergenceError
        If the root solver fails to converge within :data:`_ROOT_MAX_ITER` steps.
    """
    if tof <= 0.0:
        raise ValueError(f"tof must be positive, got {tof}")
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
    # within the single-rev valid range. t(z) is monotone-increasing inside the
    # valid range, so a bracket exists whenever the problem is well-posed; the
    # floor-anchored finder locates z_lo via bisection rather than a fixed-start
    # linear widen walk (see _find_single_rev_bracket for the rationale and the
    # task #56 regression cases).
    z_lo, z_hi, _widen_iters = _find_single_rev_bracket(
        a_coef, r1_n, r2_n, tof, mu, z_high_single_rev
    )

    # Single-rev branch: safeguarded Newton inside the established bracket. The
    # analytic dt/dz is valid here, so Newton converges quadratically; only the
    # multi-rev branches below need the derivative-free Illinois solver.
    z = _solve_single_rev_newton(z, z_lo, z_hi, a_coef, r1_n, r2_n, tof, mu)
    v1, v2 = _velocities_from_z(z, a_coef, r1_arr, r2_arr, r1_n, r2_n, mu)
    solutions = [LambertSolution(n_revs=0, branch="single", v1=v1, v2=v2)]

    # Multi-rev branches. For each n in [1, max_revs] the time curve t(z) on
    # z in ((2*pi*n)^2, (2*pi*(n+1))^2) is convex with a single minimum
    # t_min(n); when tof > t_min(n) two solutions exist (one on each side of
    # the minimum -> low/high branches). tof <= t_min(n) means revolution n is
    # infeasible and contributes nothing (skipped silently, never an error).
    eps = 1.0e-6
    for n in range(1, max_revs + 1):
        z_min, t_min = _min_time_of_revolution(n, a_coef, r1_n, r2_n, mu)
        if tof <= t_min:
            continue
        lo_endpoint = (2.0 * np.pi * n) ** 2
        hi_endpoint = (2.0 * np.pi * (n + 1)) ** 2
        rev_span = hi_endpoint - lo_endpoint
        # "low" branch: z in (lo_endpoint, z_min); "high": z in (z_min, hi_endpoint).
        # The empirical low/high <-> lamberthub low_path mapping is asserted in
        # the multi-rev crosscheck test; if inverted, swap these two labels.
        for branch_label, z_a, z_b in (
            ("low", lo_endpoint + eps * rev_span, z_min),
            ("high", z_min, hi_endpoint - eps * rev_span),
        ):
            try:
                z_sol = _solve_uv_branch(z_a, z_b, a_coef, r1_n, r2_n, tof, mu)
            except LambertConvergenceError:
                continue
            vb1, vb2 = _velocities_from_z(z_sol, a_coef, r1_arr, r2_arr, r1_n, r2_n, mu)
            solutions.append(LambertSolution(n_revs=n, branch=branch_label, v1=vb1, v2=vb2))

    return solutions


def lambert_crosscheck(
    r1: Vec3,
    r2: Vec3,
    tof: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
    n_revs: int = 0,
    branch: str = "single",
) -> dict[str, Any]:
    """Compare the in-house Lambert with ``lamberthub.izzo2015`` and ``gooding1990``.

    Used by the M1 gate test and by M3's V0/V1 closure verification. With the
    default ``n_revs=0`` this crosschecks the single-revolution transfer; pass
    ``n_revs>=1`` with ``branch in {"low", "high"}`` to crosscheck a
    multi-revolution branch (``low`` maps to ``lamberthub``'s
    ``low_path=True``, ``high`` to ``low_path=False``).
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

    sols = lambert(r1, r2, tof, mu=mu, prograde=prograde, max_revs=n_revs)
    mine = next(s for s in sols if s.n_revs == n_revs and s.branch == branch)

    r1_arr = np.asarray(r1, dtype=np.float64)
    r2_arr = np.asarray(r2, dtype=np.float64)

    low_path = branch != "high"
    v1_izzo, v2_izzo = izzo2015(
        mu, r1_arr, r2_arr, tof, M=n_revs, prograde=prograde, low_path=low_path
    )
    v1_g, v2_g = gooding1990(
        mu, r1_arr, r2_arr, tof, M=n_revs, prograde=prograde, low_path=low_path
    )

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
