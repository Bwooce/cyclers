"""Stable/unstable manifold generation for QP-tori (#522, Owen & Baresi 2024).

``genome/qp_tori.py`` (#290) solves for the invariant CIRCLE ``K(theta)`` of a
quasi-periodic torus (GMOS) but does not compute the per-point linearized
stability needed to find the torus's own stable/unstable manifolds -- this
module adds that, following Owen & Baresi's description (Sec 2.3-2.4; see
``docs/notes/2026-07-03-digest-owen-baresi-2024-knot-theory-heteroclinic.md``):
"Each of these invariant circles has an associated Floquet matrix ... used to
find [the stable/unstable manifold directions] for points in the set."

At any torus point ``(theta_long, theta_trans)``, the "Floquet matrix" is the
STM of the CR3BP flow over one stroboscopic period ``t_strob`` starting at
that point -- a genuine, well-defined 6x6 linearization regardless of which
point on the torus it is evaluated at (this module evaluates it at every
grid point directly, rather than Owen & Baresi's simplification of only
computing it once per longitude slice and reusing along latitude -- more
direct, at the cost of one STM propagation per requested grid point).

The torus is generically of SADDLE type (Kumar/Anderson/de la Llave 2025:
"most unstable periodic orbits persist as whiskered tori") -- a real,
reciprocal (lambda, 1/lambda) eigenvalue pair transverse to the torus, with
|lambda_u| > 1 > |lambda_s| = 1/lambda_u, coexisting with the trivial
unit-eigenvalue pair from the CR3BP's own symplectic/energy structure
(``search/bifurcation_detector.py`` documents this reciprocal-pair structure
for periodic-orbit monodromies; the same symplectic argument applies to any
CR3BP STM, periodic or not, since it is the Jacobian of a Hamiltonian flow).
A torus point with NO real eigenvalue pair outside a tolerance of the unit
circle is reported as locally non-hyperbolic (``lam_u`` / ``lam_s`` ``None``)
rather than guessed at.

Pure CR3BP + qp_tori dependency only. Does not depend on
``search/linking_number.py`` or ``search/torus_map_contours.py`` (those stay
CR3BP-agnostic); a caller (a future ``scripts/run_*.py`` or test) wires this
module's manifold-endpoint grids into those generic primitives.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, evaluate_torus


@dataclass(frozen=True)
class LocalStability:
    """Local stable/unstable structure at one torus point.

    ``lam_u`` (``> 1``) / ``lam_s`` (``< 1``, ``= 1 / lam_u`` to within
    solver precision) are ``None`` if no real eigenvalue pair outside
    ``hyperbolicity_tol`` of the unit circle was found at this point (the
    point is locally non-hyperbolic under this single-period STM measure).
    """

    state: NDArray[np.float64]
    stm: NDArray[np.float64]
    lam_u: float | None
    vec_u: NDArray[np.float64] | None
    lam_s: float | None
    vec_s: NDArray[np.float64] | None


def torus_point_stm(
    torus: QPTorus,
    theta_long: float,
    theta_trans: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """State and STM-over-``t_strob`` at a torus point.

    The STM is the local linearization of the stroboscopic map AT this
    specific point (not a global monodromy of the torus as a whole -- see
    module docstring).
    """
    state = evaluate_torus(torus, theta_long, theta_trans)
    arc = cr3bp.propagate(torus.system, state, torus.t_strob, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    return state, arc.stm


def local_stability(
    state: NDArray[np.float64],
    stm: NDArray[np.float64],
    *,
    hyperbolicity_tol: float = 1e-3,
    prev_vec_u: NDArray[np.float64] | None = None,
    prev_vec_s: NDArray[np.float64] | None = None,
) -> LocalStability:
    """Extract the real stable/unstable eigenpair from a torus-point STM.

    Picks the largest-magnitude real eigenvalue ``> 1 + hyperbolicity_tol``
    as unstable and the smallest-magnitude real eigenvalue ``< 1 -
    hyperbolicity_tol`` as stable (the two are expected to be reciprocal by
    the CR3BP's symplectic structure, but are extracted independently
    rather than assuming it, since this is a per-point numerical STM, not
    an exact analytic monodromy).

    Eigenvector SIGN is fixed by continuity with ``prev_vec_u``/
    ``prev_vec_s`` when supplied (dot-product sign flip) -- callers walking
    a grid of nearby torus points should thread the previous point's
    eigenvectors through to keep a consistent, non-jumping bundle field
    (the same continuity discipline already used for periodic-orbit Floquet
    eigenvectors elsewhere in this codebase, e.g.
    ``genome/heteroclinic_cycle.py``'s ``_planar_floquet_pair``).
    """
    eigvals, eigvecs = np.linalg.eig(stm)
    real_mask = np.abs(np.imag(eigvals)) < 1e-8
    lam_u: float | None = None
    vec_u: NDArray[np.float64] | None = None
    lam_s: float | None = None
    vec_s: NDArray[np.float64] | None = None
    best_u = 1.0 + hyperbolicity_tol
    best_s = 1.0 - hyperbolicity_tol
    for i in range(eigvals.shape[0]):
        if not real_mask[i]:
            continue
        lam = float(np.real(eigvals[i]))
        abs_lam = abs(lam)
        if abs_lam > best_u:
            best_u = abs_lam
            lam_u = lam
            vec_u = np.real(eigvecs[:, i]).astype(np.float64)
        if 0.0 < abs_lam < best_s:
            best_s = abs_lam
            lam_s = lam
            vec_s = np.real(eigvecs[:, i]).astype(np.float64)

    if vec_u is not None:
        vec_u = vec_u / np.linalg.norm(vec_u)
        if prev_vec_u is not None and float(vec_u @ prev_vec_u) < 0.0:
            vec_u = -vec_u
    if vec_s is not None:
        vec_s = vec_s / np.linalg.norm(vec_s)
        if prev_vec_s is not None and float(vec_s @ prev_vec_s) < 0.0:
            vec_s = -vec_s

    return LocalStability(state=state, stm=stm, lam_u=lam_u, vec_u=vec_u, lam_s=lam_s, vec_s=vec_s)


@dataclass(frozen=True)
class ManifoldGrid:
    """Endpoint states of manifold trajectories propagated to a surface of section.

    ``origins`` is the ``(n_long, n_lat, 2)`` grid of ``(theta_long,
    theta_trans)`` values each trajectory was seeded from. ``endpoints`` is
    the matching ``(n_long, n_lat, 6)`` grid of states AT the surface of
    section (``nan``-filled where the trajectory never crossed within
    ``t_max``). ``hyperbolic`` is a ``(n_long, n_lat)`` boolean grid marking
    which origins had a valid local stable/unstable direction to perturb
    along at all.
    """

    origins: NDArray[np.float64]
    endpoints: NDArray[np.float64]
    hyperbolic: NDArray[np.bool_]


def _crossing_state(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    *,
    direction: float,
    surface_x: float,
    t_max: float,
    rtol: float,
    atol: float,
) -> NDArray[np.float64] | None:
    """Propagate ``state0`` (forward if ``direction>0`` else backward) until
    ``x == surface_x``, or return ``None`` if no crossing within ``|t_max|``.
    """
    from scipy.integrate import solve_ivp

    def _x_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[0] - surface_x)

    _x_event.terminal = True  # type: ignore[attr-defined]
    _x_event.direction = 0.0  # type: ignore[attr-defined]

    t_span = (0.0, math.copysign(t_max, direction))
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        state0,
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_x_event,
        max_step=abs(t_max) / 20.0,
    )
    if sol.t_events is None or len(sol.t_events[0]) == 0:
        return None
    return np.asarray(sol.y_events[0][0], dtype=np.float64)


def _compute_longitude_row(
    args: tuple[
        int,
        float,
        NDArray[np.float64],
        QPTorus,
        str,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
    ],
) -> tuple[int, NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    (
        i,
        theta_long,
        thetas_lat,
        torus,
        branch,
        eps,
        sign,
        surface_x,
        t_max,
        hyperbolicity_tol,
        rtol,
        atol,
    ) = args
    n_lat = len(thetas_lat)
    row_origins = np.zeros((n_lat, 2))
    row_endpoints = np.full((n_lat, 6), np.nan)
    row_hyperbolic = np.zeros(n_lat, dtype=bool)

    prev_vec = None
    for j, theta_trans in enumerate(thetas_lat):
        row_origins[j, 0] = theta_long
        row_origins[j, 1] = theta_trans
        state, stm = torus_point_stm(
            torus, float(theta_long), float(theta_trans), rtol=rtol, atol=atol
        )
        stab = local_stability(
            state,
            stm,
            hyperbolicity_tol=hyperbolicity_tol,
            prev_vec_u=prev_vec if branch == "unstable" else None,
            prev_vec_s=prev_vec if branch == "stable" else None,
        )
        vec = stab.vec_u if branch == "unstable" else stab.vec_s
        if vec is None:
            continue
        if prev_vec is None and vec[0] * sign < 0.0:
            vec = -vec
        row_hyperbolic[j] = True
        prev_vec = vec
        perturbed = state + eps * vec
        direction = 1.0 if branch == "unstable" else -1.0
        crossing = _crossing_state(
            torus.system,
            perturbed,
            direction=direction,
            surface_x=surface_x,
            t_max=t_max,
            rtol=rtol,
            atol=atol,
        )
        if crossing is not None:
            row_endpoints[j, :] = crossing

    return i, row_origins, row_endpoints, row_hyperbolic


def torus_manifold_grid(
    torus: QPTorus,
    *,
    n_long: int,
    n_lat: int,
    branch: str,
    eps: float = 1e-6,
    sign: float = 1.0,
    surface_x: float,
    t_max: float,
    hyperbolicity_tol: float = 1e-3,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> ManifoldGrid:
    """Generate the stable (``branch="stable"``) or unstable
    (``branch="unstable"``) manifold's surface-of-section endpoint grid.

    For each of the ``n_long x n_lat`` torus points, perturbs by
    ``sign * eps * vec_{s,u}`` and propagates BACKWARD (stable) or FORWARD
    (unstable) to ``x = surface_x``, recording the crossing state.
    Eigenvector continuity is threaded along the LATITUDE (``theta_trans``)
    direction within each longitude row (the direction Owen & Baresi's own
    per-longitude "invariant circle" grouping matches).
    """
    if branch not in ("stable", "unstable"):
        raise ValueError(f"branch must be 'stable' or 'unstable', got {branch!r}")
    thetas_long = 2.0 * math.pi * np.arange(n_long) / n_long
    thetas_lat = 2.0 * math.pi * np.arange(n_lat) / n_lat

    origins = np.zeros((n_long, n_lat, 2), dtype=np.float64)
    endpoints = np.full((n_long, n_lat, 6), np.nan, dtype=np.float64)
    hyperbolic = np.zeros((n_long, n_lat), dtype=np.bool_)

    if n_long >= 8:
        import concurrent.futures

        tasks = [
            (
                i,
                theta_long,
                thetas_lat,
                torus,
                branch,
                eps,
                sign,
                surface_x,
                t_max,
                hyperbolicity_tol,
                rtol,
                atol,
            )
            for i, theta_long in enumerate(thetas_long)
        ]
        max_workers = min(12, n_long)
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(_compute_longitude_row, tasks))

        for i, row_origins, row_endpoints, row_hyperbolic in results:
            origins[i] = row_origins
            endpoints[i] = row_endpoints
            hyperbolic[i] = row_hyperbolic
    else:
        for i, theta_long in enumerate(thetas_long):
            prev_vec: NDArray[np.float64] | None = None
            for j, theta_trans in enumerate(thetas_lat):
                origins[i, j, 0] = theta_long
                origins[i, j, 1] = theta_trans
                state, stm = torus_point_stm(
                    torus, float(theta_long), float(theta_trans), rtol=rtol, atol=atol
                )
                stab = local_stability(
                    state,
                    stm,
                    hyperbolicity_tol=hyperbolicity_tol,
                    prev_vec_u=prev_vec if branch == "unstable" else None,
                    prev_vec_s=prev_vec if branch == "stable" else None,
                )
                vec = stab.vec_u if branch == "unstable" else stab.vec_s
                if vec is None:
                    continue
                if prev_vec is None and vec[0] * sign < 0.0:
                    vec = -vec
                hyperbolic[i, j] = True
                prev_vec = vec
                perturbed = state + eps * vec
                direction = 1.0 if branch == "unstable" else -1.0
                crossing = _crossing_state(
                    torus.system,
                    perturbed,
                    direction=direction,
                    surface_x=surface_x,
                    t_max=t_max,
                    rtol=rtol,
                    atol=atol,
                )
                if crossing is not None:
                    endpoints[i, j, :] = crossing

    return ManifoldGrid(origins=origins, endpoints=endpoints, hyperbolic=hyperbolic)


__all__ = [
    "LocalStability",
    "ManifoldGrid",
    "local_stability",
    "torus_manifold_grid",
    "torus_point_stm",
]
