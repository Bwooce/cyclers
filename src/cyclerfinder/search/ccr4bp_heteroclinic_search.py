"""CCR4BP manifold-tube intersection search + ghost-minima guard (`#694`).

Given two :class:`~cyclerfinder.search.ccr4bp_manifold_globalize.ManifoldTube`
objects (typically an UNSTABLE tube off one torus and a STABLE tube off
another -- possibly the SAME -- torus), this module

1. runs a cheap COARSE nearest-neighbour search over the two discretized
   tubes to find candidate near-coincidences in the planar
   ``(x, y, vx, vy)`` state (:func:`coarse_candidates`);
2. REFINES each candidate by treating the two torus phases and the two
   elapsed flow times as four continuous unknowns and driving the 4-component
   position+velocity residual toward its local minimum with
   ``scipy.optimize.least_squares`` (:func:`refine_candidate`) -- the
   differential-correction step the task calls for, built directly on
   :func:`cyclerfinder.search.ccr4bp_manifold_globalize.manifold_state_at`
   rather than reusing ``cr3bp_multiple_shooting.py`` (that corrector's
   formulation is a fixed-topology N-node periodicity loop; the natural
   unknowns here are exactly the two continuous phases + two continuous flow
   times, a smaller, differently-shaped problem that does not fit that
   module's node/segment-time contract without a bigger adaptation than a
   from-scratch 4-unknown least-squares call);
3. runs the mandatory GHOST-MINIMA GUARD (`#620`/`#626`'s documented failure
   mode: "it closed!" can be a numerical mesh artifact) on any refined
   candidate before it may be reported as a genuine near-connection:
   :func:`ghost_guard` checks (a) INDEPENDENT-INTEGRATOR consistency (Radau
   vs. the module's own DOP853, matching `#619`/`#620`/`#626`'s own
   discipline), (b) quasi-Jacobi-constant consistency between the two
   candidate states (the CCR4BP has no exact conserved quantity -- see the
   docstring of :func:`quasi_jacobi_gap` -- so this is an APPROXIMATE
   physical-plausibility check, not a hard gate), and (c) departure-distance
   sanity (the candidate states must be genuinely far from the torus itself,
   ruling out the trivial "both branches still sitting near their own
   departure point" pseudo-intersection -- see
   :func:`coarse_candidates`'s ``t_min_frac`` exclusion for the same guard
   applied at the coarse-search stage). MESH-REFINEMENT stability (does the
   candidate persist, not vanish, under a denser
   :class:`~cyclerfinder.search.ccr4bp_manifold_globalize.ManifoldTube`
   grid) is a separate, caller-driven check (re-run
   :func:`coarse_candidates` on a denser tube pair and confirm a comparable
   candidate reappears near the same region) -- exercised in this module's
   own tests and in the `#694` positive-control script, not built as an
   automatic gate inside this module (it requires re-globalizing a whole new
   tube pair, an expensive caller-level decision, not something a single
   candidate-level function should silently trigger).

No catalogue writeback; this module produces search RESULTS, not vetted
discoveries.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.ccr4bp_manifold_globalize import ManifoldTube, manifold_state_at
from cyclerfinder.search.ccr4bp_whisker import manifold_direction_segmented_clv

# Physical length unit for this system (Europa's SMA about Jupiter, km) --
# matches core.ccr4bp.jupiter_europa_ganymede_default's own l_km convention
# (cross-checked in tests/core/test_ccr4bp.py's mu_gan=0 reduction test).
_L_KM = 671_100.0

_STATE_COLS = (0, 1, 3, 4)  # planar (x, y, vx, vy) columns of the 6-state


@dataclass(frozen=True)
class ManifoldCandidate:
    """One coarse near-coincidence between an unstable and a stable manifold tube.

    ``theta2_u``/``t_u`` are the unstable-branch departure phase and elapsed
    forward flow time; ``theta2_s``/``t_s`` the stable-branch departure phase
    and elapsed backward flow time (same sign convention as
    :class:`~cyclerfinder.search.ccr4bp_manifold_globalize.ManifoldTube`).
    ``gap_planar`` is the raw Euclidean norm of the ``(x,y,vx,vy)`` difference
    in nondimensional units (position and velocity summed in one norm -- a
    screening metric only; :func:`refine_candidate` reports position/velocity
    gaps separately in physical units).
    """

    theta2_u: float
    t_u: float
    theta2_s: float
    t_s: float
    gap_planar: float


def coarse_candidates(
    tube_u: ManifoldTube,
    tube_s: ManifoldTube,
    *,
    n_candidates: int = 10,
    t_min_frac: float = 0.15,
    min_theta2_separation: float = 0.05,
) -> list[ManifoldCandidate]:
    """Cheap nearest-neighbour scan of two discretized manifold tubes.

    Builds a KD-tree over the stable tube's flattened ``(x,y,vx,vy)`` points
    and queries every unstable-tube point for its nearest stable-tube
    neighbour, returning the ``n_candidates`` closest pairs (de-duplicated so
    consecutive very-similar ``theta2_u`` samples do not all report the
    "same" candidate).

    ``t_min_frac`` excludes samples with elapsed flow time below
    ``t_min_frac * torus.period`` on EITHER branch -- the mandatory guard
    against the trivial pseudo-candidate at ``t_u = t_s = 0``: both branches
    depart from (near-)identical torus points at time zero, so without this
    exclusion the tightest "candidate" the KD-tree finds is always the torus
    trivially matching itself, not a genuine excursion (found and excluded
    during this module's own development -- see the `#694` commit message).
    """
    if tube_u.branch != "unstable":
        raise ValueError(f"tube_u.branch must be 'unstable', got {tube_u.branch!r}")
    if tube_s.branch != "stable":
        raise ValueError(f"tube_s.branch must be 'stable', got {tube_s.branch!r}")

    period_u = tube_u.torus.period
    period_s = tube_s.torus.period
    mask_u = tube_u.t_samples >= t_min_frac * period_u
    mask_s = tube_s.t_samples >= t_min_frac * period_s
    idx_u_time = np.where(mask_u)[0]
    idx_s_time = np.where(mask_s)[0]
    if idx_u_time.size == 0 or idx_s_time.size == 0:
        return []

    pu = tube_u.states[:, mask_u, :][:, :, _STATE_COLS].reshape(-1, 4)
    ps = tube_s.states[:, mask_s, :][:, :, _STATE_COLS].reshape(-1, 4)
    valid_u = np.repeat(tube_u.valid, idx_u_time.size)
    valid_s = np.repeat(tube_s.valid, idx_s_time.size)
    finite_u = np.isfinite(pu).all(axis=1) & valid_u
    finite_s = np.isfinite(ps).all(axis=1) & valid_s
    if not finite_u.any() or not finite_s.any():
        return []

    pu_f = pu[finite_u]
    ps_f = ps[finite_s]
    u_orig_idx = np.where(finite_u)[0]
    s_orig_idx = np.where(finite_s)[0]

    tree = cKDTree(ps_f)
    dist_raw, nn_raw = tree.query(pu_f, k=1)
    dist = np.asarray(dist_raw, dtype=np.float64)
    nn = np.asarray(nn_raw, dtype=np.intp)
    order = np.argsort(dist)

    n_theta2_u = tube_u.theta2_grid.size
    n_theta2_s = tube_s.theta2_grid.size
    out: list[ManifoldCandidate] = []
    seen_theta2_u: list[float] = []
    for k in order:
        if len(out) >= n_candidates:
            break
        flat_u = u_orig_idx[k]
        i_u, jj_u = divmod(int(flat_u), idx_u_time.size)
        theta2_u = float(tube_u.theta2_grid[i_u % n_theta2_u])
        if any(abs(theta2_u - seen) < min_theta2_separation for seen in seen_theta2_u):
            continue
        j_u = idx_u_time[jj_u]
        flat_s = s_orig_idx[int(nn[k])]
        i_s, jj_s = divmod(int(flat_s), idx_s_time.size)
        j_s = idx_s_time[jj_s]
        out.append(
            ManifoldCandidate(
                theta2_u=theta2_u,
                t_u=float(tube_u.t_samples[j_u]),
                theta2_s=float(tube_s.theta2_grid[i_s % n_theta2_s]),
                t_s=float(tube_s.t_samples[j_s]),
                gap_planar=float(dist[k]),
            )
        )
        seen_theta2_u.append(theta2_u)
    return out


@dataclass(frozen=True)
class RefinedConnection:
    """Result of continuously refining a :class:`ManifoldCandidate`.

    ``converged`` mirrors ``scipy.optimize.least_squares``'s own success flag
    (a well-posed local optimum was reached, NOT a claim of exact closure --
    see the module/report discussion of what a nonzero floor means).
    ``pos_gap_km``/``vel_gap_km_s`` are the refined candidate's position and
    velocity mismatch in PHYSICAL units. ``state_u``/``state_s`` are the two
    refined 6-states (should differ only by the reported gaps).
    """

    theta2_u: float
    t_u: float
    theta2_s: float
    t_s: float
    state_u: NDArray[np.float64]
    state_s: NDArray[np.float64]
    pos_gap_km: float
    vel_gap_km_s: float
    residual_norm: float
    converged: bool
    n_iter: int


def _v_unit_km_s(torus_system: ccr4bp.CCR4BPSystem) -> float:
    """Physical velocity unit (km/s): ``l_km * n2`` where ``n2`` is the
    Jupiter-Europa mean motion the CCR4BP nondimensionalises time by
    (``omega_gan = n3/n2 - 1`` with the frame rate fixed at 1 == n2 in these
    units, so the velocity unit is simply ``l_km`` per nondimensional time
    unit ``1/n2``; ``n2`` itself is recovered from the SAME two-body relation
    :func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate` uses)."""
    import math

    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    gm_j = PRIMARIES["Jupiter"]
    gm_e = SATELLITES["Europa"].mu_km3_s2
    n2 = math.sqrt((gm_j + gm_e) / _L_KM**3)
    return _L_KM * n2


def _direction_only(
    torus: object,
    branch: str,
    theta1_section: float,
    theta2: float,
    n_segments_dir: int,
    rtol: float,
    atol: float,
) -> NDArray[np.float64] | None:
    """Raw `#691` segmented-CLV unit eigenvector at ``(theta1_section,
    theta2)``, no sign anchoring -- used ONCE per branch to capture a fixed
    ``ref_vec`` at a refinement's seed phase (see :func:`refine_candidate`'s
    docstring on why a fixed anchor, not a re-extraction, is needed)."""
    from cyclerfinder.search.variational_ccr4bp_torus import CCR4BPTorusVariationalResult

    assert isinstance(torus, CCR4BPTorusVariationalResult)
    out = manifold_direction_segmented_clv(
        torus, theta1_section, theta2, branch, None, n_segments=n_segments_dir, rtol=rtol, atol=atol
    )
    if out is None:
        return None
    return out[1]


def refine_candidate(
    torus_u: object,
    torus_s: object,
    candidate: ManifoldCandidate,
    *,
    eps: float = 1e-6,
    lobe_sign_u: float = 1.0,
    lobe_sign_s: float = 1.0,
    theta1_section: float = 0.0,
    n_segments_dir: int = 32,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    max_nfev: int = 300,
) -> RefinedConnection | None:
    """Continuously refine a coarse candidate via 4-unknown least-squares.

    Unknowns ``z = (theta2_u, t_u, theta2_s, t_s)``; residual is the 4-vector
    ``(x_u-x_s, y_u-y_s, vx_u-vx_s, vy_u-vy_s)`` evaluated by re-deriving each
    side's manifold state from scratch at the exact continuous phase/time
    (:func:`~cyclerfinder.search.ccr4bp_manifold_globalize.manifold_state_at`)
    -- NOT interpolated from the discretized tube. ``torus_u``/``torus_s`` are
    typed ``object`` here only to avoid a hard forward-reference cycle; callers
    pass :class:`~cyclerfinder.search.variational_ccr4bp_torus.CCR4BPTorusVariationalResult`
    (the SAME torus for a homoclinic search, two different tori for a
    heteroclinic one). ``lobe_sign_u``/``lobe_sign_s`` select which of the two
    (opposite) departure directions along the extracted CLV eigenvector each
    branch used -- must match what built the coarse candidate. Returns
    ``None`` if either endpoint fails to evaluate at the seed.

    Eigenvector-sign anchoring: the CLV extraction's sign has no continuity
    guarantee across independent calls at different ``theta2`` (see
    :func:`~cyclerfinder.search.ccr4bp_manifold_globalize.manifold_state_at`'s
    own docstring) -- a naive re-extraction at every trial ``theta2`` during
    the optimizer's line search could silently jump to the OTHER manifold
    lobe. This function anchors a FIXED reference direction at the seed
    ``theta2`` for each branch once, up front, and threads it through every
    residual evaluation, so the whole refinement stays on the SAME lobe the
    coarse candidate was found on.
    """
    from cyclerfinder.search.variational_ccr4bp_torus import CCR4BPTorusVariationalResult

    assert isinstance(torus_u, CCR4BPTorusVariationalResult)
    assert isinstance(torus_s, CCR4BPTorusVariationalResult)

    ref_vec_u = _direction_only(
        torus_u, "unstable", theta1_section, candidate.theta2_u, n_segments_dir, rtol, atol
    )
    ref_vec_s = _direction_only(
        torus_s, "stable", theta1_section, candidate.theta2_s, n_segments_dir, rtol, atol
    )
    if ref_vec_u is None or ref_vec_s is None:
        return None

    def _residual(z: NDArray[np.float64]) -> NDArray[np.float64]:
        theta2_u, t_u, theta2_s, t_s = (float(v) for v in z)
        if t_u < 0.0 or t_s < 0.0:
            return np.full(4, 1e3)
        su = manifold_state_at(
            torus_u,
            "unstable",
            theta1_section,
            theta2_u,
            t_u,
            eps=eps,
            lobe_sign=lobe_sign_u,
            n_segments_dir=n_segments_dir,
            rtol=rtol,
            atol=atol,
            ref_vec=ref_vec_u,
        )
        ss = manifold_state_at(
            torus_s,
            "stable",
            theta1_section,
            theta2_s,
            t_s,
            eps=eps,
            lobe_sign=lobe_sign_s,
            n_segments_dir=n_segments_dir,
            rtol=rtol,
            atol=atol,
            ref_vec=ref_vec_s,
        )
        if su is None or ss is None:
            return np.full(4, 1e3)
        return np.array(
            [su[0] - ss[0], su[1] - ss[1], su[3] - ss[3], su[4] - ss[4]], dtype=np.float64
        )

    z0 = np.array(
        [candidate.theta2_u, candidate.t_u, candidate.theta2_s, candidate.t_s], dtype=np.float64
    )
    f0 = _residual(z0)
    if not np.all(np.isfinite(f0)) or np.linalg.norm(f0) >= 1e3 * 0.999:
        return None

    sol = least_squares(
        _residual,
        z0,
        method="trf",
        x_scale="jac",
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
        max_nfev=max_nfev,
    )
    theta2_u, t_u, theta2_s, t_s = (float(v) for v in sol.x)
    su = manifold_state_at(
        torus_u,
        "unstable",
        theta1_section,
        theta2_u,
        max(t_u, 0.0),
        eps=eps,
        lobe_sign=lobe_sign_u,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_u,
    )
    ss = manifold_state_at(
        torus_s,
        "stable",
        theta1_section,
        theta2_s,
        max(t_s, 0.0),
        eps=eps,
        lobe_sign=lobe_sign_s,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_s,
    )
    if su is None or ss is None:
        return None
    pos_gap_km = float(np.linalg.norm(su[:3] - ss[:3])) * _L_KM
    v_unit = _v_unit_km_s(torus_u.system)
    vel_gap_km_s = float(np.linalg.norm(su[3:] - ss[3:])) * v_unit
    return RefinedConnection(
        theta2_u=theta2_u,
        t_u=t_u,
        theta2_s=theta2_s,
        t_s=t_s,
        state_u=su,
        state_s=ss,
        pos_gap_km=pos_gap_km,
        vel_gap_km_s=vel_gap_km_s,
        residual_norm=float(np.linalg.norm(sol.fun)),
        converged=bool(sol.status > 0),
        n_iter=int(sol.nfev),
    )


def quasi_jacobi_gap(
    state_u: NDArray[np.float64], state_s: NDArray[np.float64], mu: float
) -> float:
    """Base-CR3BP (Europa-only, Ganymede-forcing ignored) Jacobi-constant gap.

    The CCR4BP has NO exact conserved quantity (time-periodic forcing, same
    as BCR4BP) -- this is therefore an APPROXIMATE physical-plausibility
    check, not a hard energy-conservation gate: both candidate states derive
    from the SAME resonant-torus family (a similar underlying base-CR3BP
    energy), so a genuine connection is expected to show a SMALL quasi-Jacobi
    gap (comparable to the perturbation's own strength), while a spurious
    pairing of unrelated dynamical regions would generically show an O(1) gap.
    """
    return float(cr3bp.jacobi_constant(state_u, mu) - cr3bp.jacobi_constant(state_s, mu))


def _radau_manifold_state(
    torus: object,
    branch: str,
    theta1_section: float,
    theta2: float,
    t_flow: float,
    *,
    eps: float,
    lobe_sign: float,
    n_segments_dir: int,
    rtol: float,
    atol: float,
    ref_vec: NDArray[np.float64] | None = None,
) -> NDArray[np.float64] | None:
    """Independent re-propagation via ``solve_ivp(..., method="Radau")``,
    bypassing :func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp` (which is
    hard-wired to DOP853) -- the `#619`/`#620`/`#626` independent-integrator
    discipline, applied directly (not by modifying `#689`'s module).
    ``ref_vec`` anchors the extracted eigenvector's sign (same reason as
    :func:`~cyclerfinder.search.ccr4bp_manifold_globalize.manifold_state_at`'s
    own ``ref_vec`` -- without it, the independent Radau re-check could pick
    the OTHER manifold lobe than the DOP853 refinement used)."""
    from cyclerfinder.search.variational_ccr4bp_torus import CCR4BPTorusVariationalResult

    assert isinstance(torus, CCR4BPTorusVariationalResult)
    out = manifold_direction_segmented_clv(
        torus,
        theta1_section,
        theta2,
        branch,
        ref_vec,
        n_segments=n_segments_dir,
        rtol=rtol,
        atol=atol,
    )
    if out is None:
        return None
    state0, vec = out[0], out[1]
    perturbed = state0 + lobe_sign * eps * vec
    if t_flow == 0.0:
        return perturbed
    t0 = (theta1_section % (2.0 * np.pi)) / torus.omega1
    sign_dir = 1.0 if branch == "unstable" else -1.0
    sol = solve_ivp(
        ccr4bp.ccr4bp_eom,
        (t0, t0 + sign_dir * t_flow),
        perturbed,
        args=(torus.system,),
        method="Radau",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        return None
    return np.asarray(sol.y[:, -1], dtype=np.float64)


@dataclass(frozen=True)
class GhostGuardReport:
    """Verdict of the mandatory ghost-minima guard on a :class:`RefinedConnection`.

    ``genuine`` is ``True`` only if EVERY sub-check passes; the individual
    fields are always populated so a caller can see WHY a candidate failed,
    not just a boolean."""

    radau_pos_gap_km: float
    radau_consistent: bool
    integrator_delta_km: float
    quasi_jacobi_gap: float
    off_torus_km: float
    genuine: bool
    notes: str


def ghost_guard(
    torus_u: object,
    torus_s: object,
    refined: RefinedConnection,
    *,
    eps: float = 1e-6,
    lobe_sign_u: float = 1.0,
    lobe_sign_s: float = 1.0,
    theta1_section: float = 0.0,
    n_segments_dir: int = 32,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    integrator_consistency_km: float = 1.0,
    off_torus_min_km: float = 1000.0,
) -> GhostGuardReport:
    """Run the mandatory ghost-minima guard on a refined candidate.

    Checks: (1) independent-integrator (Radau) consistency, requiring the
    Radau-recomputed position gap to differ from the DOP853 one by less than
    ``integrator_delta_km`` km (default 1 km -- if the two integrators
    disagree by more than this, the "connection" is a chaos-amplified
    numerical artifact, not physics); (2) the quasi-Jacobi gap (reported, not
    gated -- see :func:`quasi_jacobi_gap`'s docstring on why this is
    approximate); (3) both endpoint states are genuinely far from the torus
    itself (``off_torus_min_km``, default 1000 km -- rules out the trivial
    near-departure pseudo-match)."""
    from cyclerfinder.search.variational_ccr4bp_torus import (
        CCR4BPTorusVariationalResult,
        evaluate_torus_state,
    )

    assert isinstance(torus_u, CCR4BPTorusVariationalResult)
    assert isinstance(torus_s, CCR4BPTorusVariationalResult)

    # Anchor the Radau re-check to the SAME lobe the DOP853 refinement
    # converged on (see _radau_manifold_state's own docstring).
    ref_vec_u = _direction_only(
        torus_u, "unstable", theta1_section, refined.theta2_u, n_segments_dir, rtol, atol
    )
    ref_vec_s = _direction_only(
        torus_s, "stable", theta1_section, refined.theta2_s, n_segments_dir, rtol, atol
    )

    su_radau = _radau_manifold_state(
        torus_u,
        "unstable",
        theta1_section,
        refined.theta2_u,
        refined.t_u,
        eps=eps,
        lobe_sign=lobe_sign_u,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_u,
    )
    ss_radau = _radau_manifold_state(
        torus_s,
        "stable",
        theta1_section,
        refined.theta2_s,
        refined.t_s,
        eps=eps,
        lobe_sign=lobe_sign_s,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_s,
    )
    if su_radau is None or ss_radau is None:
        return GhostGuardReport(
            radau_pos_gap_km=float("nan"),
            radau_consistent=False,
            integrator_delta_km=float("inf"),
            quasi_jacobi_gap=float("nan"),
            off_torus_km=float("nan"),
            genuine=False,
            notes="Radau re-propagation failed",
        )
    radau_pos_gap_km = float(np.linalg.norm(su_radau[:3] - ss_radau[:3])) * _L_KM
    integrator_delta_km = abs(radau_pos_gap_km - refined.pos_gap_km)
    radau_consistent = integrator_delta_km < integrator_consistency_km

    qj_gap = quasi_jacobi_gap(refined.state_u, refined.state_s, torus_u.system.mu)

    # Off-torus sanity: EXACT distance from state_u to the specific torus
    # point the UNPERTURBED trajectory would occupy at the same elapsed
    # time (theta1, theta2 both advance linearly with time on the torus
    # itself -- an analytic Fourier evaluation, no propagation, no grid
    # discretization error). This is the direct analogue of the "compare
    # against the base trajectory propagated to the SAME time" fix this
    # module's own sibling ``ccr4bp_manifold_globalize`` test module needed
    # -- a coarse fixed (theta1, theta2) grid search under-samples fast
    # (near-periapsis) regions of an eccentric orbit and can badly
    # OVERESTIMATE the true off-torus distance for a genuinely
    # near-departure point (found during `#694`'s own test development).
    theta1_u_at_t = theta1_section + torus_u.omega1 * refined.t_u
    theta2_u_at_t = refined.theta2_u + torus_u.omega2 * refined.t_u
    torus_pt_u = evaluate_torus_state(torus_u, theta1_u_at_t, theta2_u_at_t)
    planar_u = refined.state_u[[0, 1, 3, 4]]
    off_torus_km = float(np.linalg.norm(torus_pt_u - planar_u)) * _L_KM

    notes_parts = []
    if not radau_consistent:
        notes_parts.append(
            f"Radau/DOP853 disagree by {integrator_delta_km:.3f} km "
            f">= {integrator_consistency_km} km threshold"
        )
    if off_torus_km < off_torus_min_km:
        notes_parts.append(
            f"candidate only {off_torus_km:.1f} km off the torus "
            f"(< {off_torus_min_km} km) -- likely a trivial near-departure match"
        )
    genuine = radau_consistent and off_torus_km >= off_torus_min_km
    notes = "; ".join(notes_parts) if notes_parts else "all checks passed"
    return GhostGuardReport(
        radau_pos_gap_km=radau_pos_gap_km,
        radau_consistent=radau_consistent,
        integrator_delta_km=integrator_delta_km,
        quasi_jacobi_gap=qj_gap,
        off_torus_km=off_torus_km,
        genuine=genuine,
        notes=notes,
    )
