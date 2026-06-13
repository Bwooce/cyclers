"""Single-impulse reachable-set SPIKE in the CR3BP (Zhou-Armellin 2025).

A focused prototype of the *single-impulse reachable set* method of

    X. Zhou, R. Armellin, D. Qiao, X. Li, "Single-Impulse Reachable Set in
    Arbitrary Dynamics Using Polynomials," arXiv:2502.11280v1 [astro-ph.IM],
    16 Feb 2025.

See the mining note ``docs/notes/2026-06-07-zhou-2025-da-reachable-sets-mining.md``
for the full method transcription and the sourced benchmark numbers.

WHAT THIS SPIKE IS (and is NOT):

* It is the **sampling** realization of the method, NOT the differential-algebra
  (DA) / polynomial realization. Zhou build a high-order Taylor map
  ``x ~ T_x(alpha, beta)`` of the final state in the two impulse angles by
  integrating the CR3BP flow in DA arithmetic (DACE / DACEyPy), then extract the
  reachable-set (RS) *boundary* as the envelope of the projected curve family
  (Eqs. 8-45). That requires a DA-evaluable force model -- a dependency our pure
  float REBOUND/scipy stack does not host (mining note Sec. 3). This spike
  instead samples the impulse-angle domain directly and propagates each sample
  with the existing float propagator ``core.cr3bp.propagate``. The DA layer is the
  >84% CPU optimization (mining note Sec. 2.5); the *geometry it computes* -- the
  reachable footprint on a plane orthogonal to the nominal velocity -- is exactly
  what this spike computes by brute force. So this is the method's geometry
  without the method's speed.

* The impulse model is faithful: a single bounded impulse of **arbitrary
  direction** at a fixed epoch, magnitude on the max sphere ``||dv|| = dv_max``
  (Zhou Eq. 5: the RS boundary lives on the max-impulse sphere), parameterized by
  elevation ``alpha in [-pi/2, pi/2]`` and azimuth ``beta in [-pi, pi]`` in the
  local velocity frame (Eq. 4). This CHANGES the energy (Jacobi constant), which
  is the essential contrast with the Braik-Ross energy-preserving heading-change
  maneuver in ``reachable_network.py`` (which stays on one C_J manifold).

* The auxiliary-plane projection is faithful (Eq. 8): the RS is characterized on
  a plane orthogonal to the nominal velocity, with axes built from the angular
  momentum ``h = r x v``, the velocity ``v``, and ``h x v``. The 3-D reachable
  point cloud projects to a 2-D footprint ``(x_p, z_p)`` whose boundary is the RS
  envelope. This spike reports the footprint's convex-hull area and extent rather
  than the exact alpha-shape envelope (the envelope/root-find/local-poly layer of
  Eqs. 34-51 is the DA-specific part not reproduced here).

CROSS-CHECK (independent recompute, NOT a circular golden): the method invariant
the paper validates is *RS-boundary containment* -- the Monte-Carlo cloud of
random feasible impulses lies inside the RS boundary, with error index
``P = d_max^2 / S_RS`` well under 0.1% (mining note Sec. 4). Here the independent
truth set is a dense uniform-random MC cloud of impulses sampled over the full
ball ``||dv|| <= dv_max``; the spike's reachable footprint is the convex hull of a
coarse (alpha, beta) grid on the max sphere. Containment of the MC cloud by the
grid footprint (within a tolerance) is the cross-check. This needs no published
state vector and asserts no value our own code declared golden.

Pure: math / numpy / scipy + ``cyclerfinder.core.cr3bp``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import ConvexHull

import cyclerfinder.core.cr3bp as cr3bp

# ---------------------------------------------------------------------------
# Impulse parameterization on the max-magnitude sphere (Zhou Eqs. 4-6).
# ---------------------------------------------------------------------------


def impulse_vector(dv_mag: float, alpha: float, beta: float) -> NDArray[np.float64]:
    """Impulse ``dv`` in spherical angles in the *local velocity frame* (Eq. 4).

    The impulse is expressed in a right-handed frame whose first axis is the
    nominal velocity direction; ``alpha`` is the elevation out of the
    velocity-normal plane and ``beta`` the azimuth about the velocity axis. The
    returned vector is in that local frame and must be rotated into the rotating
    CR3BP frame by :func:`velocity_frame` before being added to the state
    velocity. With ``alpha = beta = 0`` the impulse is purely along-track.

    Zhou parameterize the boundary at ``||dv|| = dv_max`` (Eq. 5); a general
    ``dv_mag <= dv_max`` is allowed here for interior MC sampling.
    """
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    # alpha: elevation toward the (h x v) and h axes; beta: azimuth in the
    # plane orthogonal to v. Component along v is cos(alpha)*cos(beta) so that
    # (0, 0) is purely prograde.
    return dv_mag * np.array([ca * cb, ca * sb, sa], dtype=np.float64)


def velocity_frame(r: NDArray[np.float64], v: NDArray[np.float64]) -> NDArray[np.float64]:
    """Right-handed local frame columns ``[e_v, e_n, e_h]`` from state ``(r, v)``.

    Zhou's auxiliary-plane transform ``T(x_bar)`` (Eq. 8) is built from the
    angular momentum ``h = r x v``, the velocity ``v``, and ``h x v``. The same
    triad is used both to inject the impulse (the impulse local frame) and to
    project the reachable cloud onto the plane orthogonal to the nominal velocity.

    Columns:
      * ``e_v`` = v / |v|              (nominal velocity direction)
      * ``e_h`` = h / |h|             (orbit-normal, h = r x v)
      * ``e_n`` = e_h x e_v           (in-plane, completing the right-handed set)

    Raises
    ------
    ValueError
        If ``v`` or ``h = r x v`` is (numerically) zero, so no frame exists.
    """
    vnorm = float(np.linalg.norm(v))
    if vnorm < 1e-14:
        raise ValueError("velocity_frame: |v| ~ 0, no velocity frame")
    e_v = v / vnorm
    h = np.cross(r, v)
    hnorm = float(np.linalg.norm(h))
    if hnorm < 1e-14:
        raise ValueError("velocity_frame: |r x v| ~ 0, rectilinear state")
    e_h = h / hnorm
    e_n = np.cross(e_h, e_v)
    e_n /= float(np.linalg.norm(e_n))
    return np.column_stack([e_v, e_n, e_h])


def apply_impulse(
    state6: NDArray[np.float64], dv_local: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Add a local-frame impulse ``dv_local`` to the velocity of ``state6``.

    The impulse (expressed in the :func:`velocity_frame` of the state) is rotated
    into the rotating CR3BP frame and added to the velocity; position is
    unchanged (instantaneous impulse).
    """
    r = np.asarray(state6, float)[:3]
    v = np.asarray(state6, float)[3:6]
    frame = velocity_frame(r, v)
    dv_rot = frame @ np.asarray(dv_local, float)
    out = np.asarray(state6, float).copy()
    out[3:6] = v + dv_rot
    return out


# ---------------------------------------------------------------------------
# Reachable cloud: propagate single impulses over a horizon.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReachableCloud:
    """Propagated single-impulse final states and their auxiliary-plane footprint.

    ``final_states`` are the 6-D rotating-frame states at the horizon ``t_f`` for
    each sampled impulse. ``footprint`` is the 2-D projection ``(x_p, z_p)`` of the
    final *positions* onto the plane orthogonal to the nominal (un-maneuvered)
    final velocity, relative to the nominal final position (Zhou Eq. 8-11); the
    nominal final state is row 0 of ``final_states`` if ``include_nominal`` was set.
    """

    final_states: NDArray[np.float64]  # (N, 6)
    footprint: NDArray[np.float64]  # (N, 2)
    nominal_final: NDArray[np.float64]  # (6,)


def _project_footprint(
    final_states: NDArray[np.float64], nominal_final: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Project final positions onto the plane orthogonal to the nominal velocity.

    Zhou Eq. 8-11: the RS is studied on an auxiliary plane through the nominal
    final position, orthogonal to the nominal final velocity. The two in-plane
    axes are ``e_n`` and ``e_h`` of the nominal final state's velocity frame; the
    footprint coordinate of a maneuvered final position ``r`` is
    ``(x_p, z_p) = ((r - r_bar) . e_n, (r - r_bar) . e_h)``.
    """
    r_bar = nominal_final[:3]
    v_bar = nominal_final[3:6]
    frame = velocity_frame(r_bar, v_bar)
    e_n = frame[:, 1]
    e_h = frame[:, 2]
    dr = final_states[:, :3] - r_bar
    x_p = dr @ e_n
    z_p = dr @ e_h
    return np.column_stack([x_p, z_p])


def reachable_cloud(
    system: cr3bp.CR3BPSystem,
    seed_state: NDArray[np.float64],
    dv_mag: float,
    t_f: float,
    *,
    angle_samples: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None,
    n_alpha: int = 13,
    n_beta: int = 25,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> ReachableCloud:
    """Single-impulse reachable cloud on the max sphere over ``(alpha, beta)``.

    From ``seed_state``, apply an impulse of magnitude ``dv_mag`` (nondimensional)
    in every sampled ``(alpha, beta)`` direction on the sphere, propagate each for
    nondimensional time ``t_f`` with :func:`cyclerfinder.core.cr3bp.propagate`, and
    project the final positions onto the plane orthogonal to the *nominal* final
    velocity (no-impulse propagation), per Zhou Eqs. 4-11.

    If ``angle_samples`` is given it is a pair ``(alphas, betas)`` of 1-D arrays
    meshed into the grid; otherwise a uniform ``n_alpha x n_beta`` grid on
    ``[-pi/2, pi/2] x [-pi, pi]`` is used.

    Failed propagations (collision / integrator failure) are dropped.
    """
    if angle_samples is None:
        alphas = np.linspace(-math.pi / 2.0, math.pi / 2.0, n_alpha)
        betas = np.linspace(-math.pi, math.pi, n_beta)
    else:
        alphas, betas = angle_samples
    aa, bb = np.meshgrid(alphas, betas, indexing="ij")
    angles = np.column_stack([aa.ravel(), bb.ravel()])

    nominal_arc = cr3bp.propagate(system, seed_state, t_f, rtol=rtol, atol=atol)
    nominal_final = nominal_arc.state_f

    finals: list[NDArray[np.float64]] = []
    for alpha, beta in angles:
        dv_local = impulse_vector(dv_mag, float(alpha), float(beta))
        man0 = apply_impulse(seed_state, dv_local)
        try:
            arc = cr3bp.propagate(system, man0, t_f, rtol=rtol, atol=atol)
        except RuntimeError:
            continue
        finals.append(arc.state_f)
    final_states = np.asarray(finals, dtype=np.float64)
    footprint = _project_footprint(final_states, nominal_final)
    return ReachableCloud(
        final_states=final_states, footprint=footprint, nominal_final=nominal_final
    )


def monte_carlo_cloud(
    system: cr3bp.CR3BPSystem,
    seed_state: NDArray[np.float64],
    dv_max: float,
    t_f: float,
    nominal_final: NDArray[np.float64],
    *,
    n_samples: int = 400,
    on_sphere: bool = True,
    seed: int = 0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64]:
    """Independent MC truth cloud of single impulses, projected to the footprint.

    Draws ``n_samples`` impulses with uniformly random directions (uniform on the
    unit sphere) and magnitude either fixed at ``dv_max`` (``on_sphere=True``, the
    RS boundary) or uniform in ``[0, dv_max]`` (``on_sphere=False``, the full
    feasible ball, used as the containment truth set), propagates each, and
    projects onto the same nominal-velocity-orthogonal plane. This is the
    *independent recompute* the spike cross-checks against the structured-grid
    footprint -- it shares no anchor with the grid sampling.
    """
    rng = np.random.default_rng(seed)
    finals: list[NDArray[np.float64]] = []
    for _ in range(n_samples):
        u = rng.normal(size=3)
        u /= float(np.linalg.norm(u))
        mag = dv_max if on_sphere else dv_max * float(rng.uniform(0.0, 1.0))
        man0 = np.asarray(seed_state, float).copy()
        man0[3:6] = man0[3:6] + mag * u
        try:
            arc = cr3bp.propagate(system, man0, t_f, rtol=rtol, atol=atol)
        except RuntimeError:
            continue
        finals.append(arc.state_f)
    return _project_footprint(np.asarray(finals, dtype=np.float64), nominal_final)


# ---------------------------------------------------------------------------
# Footprint characterization + containment cross-check.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FootprintMetrics:
    """Reachable-footprint summary on the auxiliary plane.

    ``area`` is the convex-hull area ``S_RS`` of the footprint; ``extent`` is the
    bounding-box half-extents ``(dx_p, dz_p)``; ``centroid`` is the hull centroid.
    """

    area: float
    extent: tuple[float, float]
    centroid: tuple[float, float]
    n_points: int


def footprint_metrics(footprint: NDArray[np.float64]) -> FootprintMetrics:
    """Convex-hull area / extent / centroid of a 2-D footprint cloud."""
    pts = np.asarray(footprint, dtype=np.float64)
    n = pts.shape[0]
    if n < 3:
        return FootprintMetrics(area=0.0, extent=(0.0, 0.0), centroid=(0.0, 0.0), n_points=n)
    hull = ConvexHull(pts)
    cx, cz = float(pts[:, 0].mean()), float(pts[:, 1].mean())
    dx = float(pts[:, 0].max() - pts[:, 0].min()) / 2.0
    dz = float(pts[:, 1].max() - pts[:, 1].min()) / 2.0
    return FootprintMetrics(
        area=float(hull.volume),  # 2-D ConvexHull.volume is the area
        extent=(dx, dz),
        centroid=(cx, cz),
        n_points=n,
    )


def _point_in_hull(point: NDArray[np.float64], hull: ConvexHull, tol: float) -> bool:
    """True iff ``point`` is inside (or within ``tol`` of) the hull's facets.

    Uses the hull's facet half-space equations ``A x + b <= 0``; a point is inside
    when every facet residual is ``<= tol`` (tol allows a small slack so MC points
    on the boundary are not spuriously rejected).
    """
    a = hull.equations[:, :-1]
    b = hull.equations[:, -1]
    return bool(np.all(a @ point + b <= tol))


@dataclass(frozen=True)
class ContainmentResult:
    """Result of the MC-cloud-in-grid-footprint containment cross-check.

    ``contained_fraction`` is the fraction of MC points inside the grid footprint
    hull (within ``tol``); ``max_outside_dist`` is the largest signed facet
    violation ``d_max`` over MC points; ``error_index`` is Zhou's
    ``P = d_max^2 / S_RS`` (Eq. 55) -- the paper's reported accuracy gauge.
    """

    contained_fraction: float
    max_outside_dist: float
    error_index: float
    grid_area: float
    n_mc: int


def containment_crosscheck(
    grid_footprint: NDArray[np.float64],
    mc_footprint: NDArray[np.float64],
    *,
    tol: float = 0.0,
) -> ContainmentResult:
    """Cross-check: is the MC cloud contained by the grid footprint envelope?

    The grid footprint (structured ``(alpha, beta)`` sampling on the max sphere)
    defines a convex-hull boundary; the MC cloud (independent random impulses) is
    the truth set. Reports the contained fraction, the max outside distance
    ``d_max``, and Zhou's error index ``P = d_max^2 / S_RS`` (Eq. 55). This is the
    method invariant the paper validates (mining note Sec. 4): the MC cloud lies
    inside the RS boundary with ``P`` well under 0.1%.

    NOTE: the convex hull is an OVER-approximation of the true (possibly
    non-convex) RS envelope, so containment of an MC cloud by the grid hull is a
    *necessary* check on the grid's coverage, not a proof the grid traces the
    exact boundary -- consistent with the spike's "geometry without the DA
    envelope" scope.
    """
    grid_pts = np.asarray(grid_footprint, dtype=np.float64)
    mc_pts = np.asarray(mc_footprint, dtype=np.float64)
    hull = ConvexHull(grid_pts)
    a = hull.equations[:, :-1]
    b = hull.equations[:, -1]
    n_in = 0
    d_max = 0.0
    for p in mc_pts:
        resid = a @ p + b
        worst = float(resid.max())
        if worst <= tol:
            n_in += 1
        else:
            d_max = max(d_max, worst)
    area = float(hull.volume)
    error_index = (d_max * d_max) / area if area > 0.0 else math.inf
    frac = n_in / mc_pts.shape[0] if mc_pts.shape[0] else 0.0
    return ContainmentResult(
        contained_fraction=frac,
        max_outside_dist=d_max,
        error_index=error_index,
        grid_area=area,
        n_mc=mc_pts.shape[0],
    )


__all__ = [
    "ContainmentResult",
    "FootprintMetrics",
    "ReachableCloud",
    "apply_impulse",
    "containment_crosscheck",
    "footprint_metrics",
    "impulse_vector",
    "monte_carlo_cloud",
    "reachable_cloud",
    "velocity_frame",
]
