"""Transit vs. non-transit branch classification for libration-point manifolds (#547).

Establishes a from-first-principles, textbook-validated POSITIVE CONTROL for the
Conley-McGehee / Koon-Lo-Marsden-Ross "tube dynamics" picture in the planar CR3BP,
resolving the transit-vs-non-transit branch-classification blocker that task #534
flagged but never closed for the ``qp_tori``/``qp_torus_heteroclinic``
linking-number method family.

The classical picture (Conley 1968, McGehee 1969, KLMR 2006 Ch. 2-4; digested in
``docs/notes/2026-06-21-digest-klmr-2006-book.md``): at an energy just above the
L1 (or L2) threshold the bounded region near the libration point forms a
"bottleneck" (neck). The unstable manifold of the L1 Lyapunov orbit has TWO
branches (``+`` / ``-`` along the Floquet unstable eigenvector). One branch is a
genuine TRANSIT trajectory that threads through the neck into the secondary's
realm (crossing the ``x = 1 - mu`` surface of section at the secondary); the
OTHER branch is NON-TRANSIT -- it falls back toward the primary and never
reaches that surface. Which branch transits is decided here EMPIRICALLY (propagate
both, see which crosses), not guessed from an eigenvector-component sign -- exactly
the classification #534 said was the missing piece.

This module reuses the #314 planar-Floquet machinery in
``genome/heteroclinic_cycle.py`` (``LyapunovNode``, ``_planar_floquet_pair``,
``_seed_on_manifold``) plus the fixed-Jacobi Lyapunov corrector in
``search/cr3bp_periodic.py``. It adds only (a) a robust small-amplitude
libration-Lyapunov seed from the collinear point's linearization, and (b) the
transit classification itself.

Positive control (pinned in ``tests/genome/test_transit_manifold.py``): the
Earth-Moon L1 Lyapunov orbit at ``C ~ 3.1869`` (between ``C_L2 = 3.1722`` and
``C_L1 = 3.1883``, so the L1 neck is open and the L2 neck is closed -- the Moon
realm is bounded, the classic KLMR L1-gateway setting). Its ``+`` unstable branch
transits into the Moon realm (crosses ``x = 1 - mu``, approaches the Moon to
< 0.03 nondim); its ``-`` branch stays interior and never crosses. Pure
CR3BP + heteroclinic_cycle dependency only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.heteroclinic_cycle import (
    LyapunovNode,
    _planar_floquet_pair,
    _seed_on_manifold,
)
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x


def _collinear_linear_planar(mu: float, x_l: float) -> tuple[float, float, float]:
    """Return ``(nu, kappa, t_lin)`` for the planar linearization at a collinear point.

    ``c2`` is the classical second-order pseudo-potential coefficient on the
    x-axis at the collinear point ``x_l`` (KLMR 2006 Eq. 2.3.6):

        c2 = (1 - mu) / |x_l + mu|^3 + mu / |x_l - 1 + mu|^3.

    The planar linearized motion has one real saddle pair ``+/- lam`` and one
    imaginary center pair ``+/- i nu`` with (KLMR Eq. 2.3.10-2.3.12)

        nu^2 = ( c2 - 2 + sqrt(9 c2^2 - 8 c2) ) / 2,
        lam^2 = ( 2 - c2 + sqrt(9 c2^2 - 8 c2) ) / 2   [only used for context].

    The center pair sets the small-amplitude Lyapunov period ``t_lin = 2 pi / nu``.
    ``kappa = (nu^2 + 1 + 2 c2) / (2 nu)`` is the linear in-plane amplitude ratio
    that seeds ``ydot0`` from the x-amplitude (KLMR Eq. 2.3.17); it need only be
    good enough to start the differential corrector on the right family.
    """
    c2 = (1.0 - mu) / abs(x_l + mu) ** 3 + mu / abs(x_l - 1.0 + mu) ** 3
    disc = 9.0 * c2 * c2 - 8.0 * c2
    if disc <= 0.0:
        raise ValueError(f"non-oscillatory linearization at x_l={x_l}: 9c2^2-8c2={disc}")
    nu = math.sqrt((c2 - 2.0 + math.sqrt(disc)) / 2.0)
    kappa = (nu * nu + 1.0 + 2.0 * c2) / (2.0 * nu)
    t_lin = 2.0 * math.pi / nu
    return nu, kappa, t_lin


def libration_lyapunov(
    system: cr3bp.CR3BPSystem,
    point: str,
    amplitude: float,
    *,
    tol: float = 1e-11,
    x0_halfwidth: float = 0.15,
) -> LyapunovNode:
    """Correct a small planar Lyapunov orbit at collinear ``point`` (``L1``/``L2``).

    Seeds a small-amplitude perpendicular-crossing IC from the linearization
    (``x0 = x_l -/+ amplitude`` for ``L1``/``L2``; ``ydot0`` from the linear
    amplitude ratio) and corrects it at FIXED Jacobi with
    :func:`correct_symmetric_fixed_jacobi`, holding ``x0`` within ``x0_halfwidth``
    of the libration point so the strongly-unstable Newton step cannot slide onto
    a distant, unrelated family member. Returns a :class:`LyapunovNode` carrying
    the corrected IC, period, Jacobi constant, and the planar Floquet
    unstable/stable eigenvector pair.

    ``amplitude`` is the x-axis half-amplitude in nondimensional units; keep it
    small (~1e-3 to 1e-2) so the seed stays inside the corrector's basin. The
    resulting Jacobi constant (an OUTPUT of the chosen amplitude) is read from the
    returned node's ``jacobi`` field.
    """
    if point not in ("L1", "L2"):
        raise ValueError(f"point must be 'L1' or 'L2', got {point!r}")
    if amplitude <= 0.0:
        raise ValueError(f"amplitude must be positive, got {amplitude}")
    mu = system.mu
    x_l = lagrange_collinear_x(mu, point)
    nu, kappa, t_lin = _collinear_linear_planar(mu, x_l)
    # L1 orbit sits on the primary side (x0 < x_l), L2 beyond the secondary
    # (x0 > x_l). The perpendicular-crossing ydot0 sign follows the family.
    if point == "L1":
        x0 = x_l - amplitude
        ydot0_sign = 1.0
    else:
        x0 = x_l + amplitude
        ydot0_sign = -1.0
    ydot0_lin = math.copysign(kappa * nu * amplitude, ydot0_sign)
    jacobi = cr3bp.jacobi_constant(np.array([x0, 0.0, 0.0, 0.0, ydot0_lin, 0.0]), mu)
    orbit = correct_symmetric_fixed_jacobi(
        system,
        x0,
        jacobi,
        t_lin,
        ydot0_sign=ydot0_sign,
        tol=tol,
        x0_bounds=(x_l - x0_halfwidth, x_l + x0_halfwidth),
    )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    _lu, v_u, _ls, v_s = _planar_floquet_pair(system, state0, orbit.period)
    return LyapunovNode(
        label=point,
        state0=state0,
        period=orbit.period,
        jacobi=cr3bp.jacobi_constant(state0, mu),
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=orbit.converged,
    )


@dataclass(frozen=True)
class BranchTransit:
    """Diagnostics for one unstable-manifold branch propagated to a section.

    ``transits`` is ``True`` iff the branch reaches ``surface_x`` (crosses into
    the secondary's realm) within ``t_max``. ``n_crossings`` counts the crossings
    of ``x = surface_x``; ``first_crossing_time`` is the time of the first
    (``nan`` if none). ``min_secondary_distance`` is the closest approach to the
    secondary ``(1 - mu, 0, 0)`` over the propagation; ``x_min`` / ``x_max`` bound
    the excursion in x (a non-transit interior branch has ``x_max < surface_x``).
    """

    branch: int
    transits: bool
    n_crossings: int
    first_crossing_time: float
    min_secondary_distance: float
    x_min: float
    x_max: float


def classify_unstable_branch(
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode,
    branch: int,
    *,
    surface_x: float,
    t_max: float,
    eps: float = 1e-6,
    tau: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    max_step_frac: float = 1.0 / 400.0,
) -> BranchTransit:
    """Classify one unstable-manifold branch of ``node`` as transit or non-transit.

    Seeds ``state(tau) + branch * eps * v_u(tau)`` on the unstable manifold (via
    the #314 :func:`_seed_on_manifold`) and propagates FORWARD to ``t_max``,
    recording every crossing of ``x = surface_x``. A branch that reaches the
    surface is a transit trajectory (``transits=True``); one that never does is
    non-transit. ``branch`` must be ``+1`` or ``-1`` (the two ends of the unstable
    eigen-direction). Bounded integration -- never hangs, never fabricates a
    crossing.
    """
    if branch not in (+1, -1):
        raise ValueError(f"branch must be +1 or -1, got {branch!r}")
    mu = system.mu
    seed = _seed_on_manifold(
        system, node, tau=tau, direction="unstable", branch=branch, epsilon=eps
    )

    def _x_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[0] - surface_x)

    _x_event.terminal = False  # type: ignore[attr-defined]
    _x_event.direction = 0.0  # type: ignore[attr-defined]

    horizon = abs(float(t_max))
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, horizon),
        np.asarray(seed, float),
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_x_event,
        max_step=horizon * max_step_frac,
    )
    xs = sol.y[0]
    r2 = np.sqrt((sol.y[0] - (1.0 - mu)) ** 2 + sol.y[1] ** 2 + sol.y[2] ** 2)
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    n_cross = len(t_events)
    first_t = float(t_events[0]) if n_cross > 0 else float("nan")
    return BranchTransit(
        branch=branch,
        transits=n_cross > 0,
        n_crossings=n_cross,
        first_crossing_time=first_t,
        min_secondary_distance=float(r2.min()),
        x_min=float(xs.min()),
        x_max=float(xs.max()),
    )


__all__ = [
    "BranchTransit",
    "classify_unstable_branch",
    "libration_lyapunov",
]
