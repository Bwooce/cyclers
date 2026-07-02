"""#405 cross-system (Sun-Earth <-> Earth-Moon) heteroclinic-cycle framework.

Task 1: the inter-system frame bridge. Maps a 6-state between the SE-rotating CR3BP
frame, a common Earth-centered inertial frame (km, km/s), and the EM-rotating CR3BP
frame, parameterized by ``theta`` (the relative phase between the Sun-Earth line and
the Earth-Moon line). Correctness is gated by a round-trip identity AND a physical
Moon-position test (see tests/genome/test_cross_system_cycle.py).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode
from cyclerfinder.search.cr3bp_general_periodic_3d import correct_general_periodic_3d


def se_earth_system() -> cr3bp.CR3BPSystem:
    """Build the Sun-Earth CR3BP system (Earth is the SE secondary at 1-mu_SE).

    Thin alias for the shared registry, which now serves Sun-primary heliocentric
    pairs from ``core.constants.PLANETS`` + ``MU_SUN_KM3_S2`` (see
    ``cr3bp.cr3bp_system``, #409). Kept as a named convenience for the #405/#411
    cross-system call sites; the registry is the single source of truth.
    """
    return cr3bp.cr3bp_system("Sun", "Earth")


def em_moon_system() -> cr3bp.CR3BPSystem:
    """Build the Earth-Moon CR3BP system via the shared registry (Moon is registered)."""
    return cr3bp.cr3bp_system("Earth", "Moon")


def _rot_z(angle: float) -> NDArray[np.float64]:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


@dataclass(frozen=True)
class FrameBridge:
    """Transforms 6-states between SE-rot, Earth-centered inertial (km, km/s), EM-rot.

    SE-rot and EM-rot are nondimensional CR3BP rotating frames; the inertial frame is
    Earth-centered, dimensional. ``theta`` is the inertial angle of the EM-rot x-axis
    minus that of the SE-rot x-axis (the SE-rot x-axis is taken at inertial angle 0).
    """

    se: cr3bp.CR3BPSystem
    em: cr3bp.CR3BPSystem

    def se_rot_to_inertial(self, s: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        r = np.asarray(s[:3], float).copy()
        v = np.asarray(s[3:], float).copy()
        r = r - np.array([1.0 - self.se.mu, 0.0, 0.0])  # Earth-centered, SE-rot, nondim
        lkm, vunit = self.se.l_km, self.se.l_km / self.se.t_s
        pos_km = r * lkm
        vel_rot = v * vunit
        omega = np.array([0.0, 0.0, 1.0 / self.se.t_s])
        vel_in = vel_rot + np.cross(omega, pos_km)
        # SE-rot x-axis is taken at inertial angle 0, so no theta rotation here.
        return np.concatenate([pos_km, vel_in])

    def inertial_to_se_rot(self, x: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        pos_km = np.asarray(x[:3], float).copy()
        vel_in = np.asarray(x[3:], float).copy()
        omega = np.array([0.0, 0.0, 1.0 / self.se.t_s])
        vel_rot = vel_in - np.cross(omega, pos_km)
        lkm, vunit = self.se.l_km, self.se.l_km / self.se.t_s
        r = pos_km / lkm + np.array([1.0 - self.se.mu, 0.0, 0.0])
        v = vel_rot / vunit
        return np.concatenate([r, v])

    def em_rot_to_inertial(self, s: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        r = np.asarray(s[:3], float).copy()
        v = np.asarray(s[3:], float).copy()
        r = r - np.array([-self.em.mu, 0.0, 0.0])  # Earth-centered, EM-rot, nondim
        lkm, vunit = self.em.l_km, self.em.l_km / self.em.t_s
        pos_emrot = r * lkm
        vel_rot = v * vunit
        omega = np.array([0.0, 0.0, 1.0 / self.em.t_s])
        vel_in_emrot = vel_rot + np.cross(omega, pos_emrot)
        rot = _rot_z(theta)  # EM-rot x-axis leads SE-rot by theta in inertial frame
        return np.concatenate([rot @ pos_emrot, rot @ vel_in_emrot])

    def inertial_to_em_rot(self, x: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        rot_inv = _rot_z(-theta)
        pos_emrot = rot_inv @ np.asarray(x[:3], float)
        vel_in_emrot = rot_inv @ np.asarray(x[3:], float)
        omega = np.array([0.0, 0.0, 1.0 / self.em.t_s])
        vel_rot = vel_in_emrot - np.cross(omega, pos_emrot)
        lkm, vunit = self.em.l_km, self.em.l_km / self.em.t_s
        r = pos_emrot / lkm + np.array([-self.em.mu, 0.0, 0.0])
        v = vel_rot / vunit
        return np.concatenate([r, v])


# --- Task 3: cross-system connection corrector ------------------------------
#
# A connection between an Earth-Moon orbit and a Sun-Earth orbit cannot be matched
# on a rotating-frame Poincaré section: the two frames rotate at different rates, so
# "{y=0}" means different things in each. The match must happen in the COMMON
# Earth-centered inertial frame. We define the patch section as an inertial plane
# {x_inertial = X0}; both manifolds are propagated in their own rotating frame,
# each sample is transformed to inertial via the bridge, and the FIRST crossing of
# the section is detected (with linear refinement). The residual is the inertial
# POSITION gap at the section (km); the inertial VELOCITY gap there is the patch ΔV
# that a spacecraft must apply to transfer from the EM unstable manifold onto the SE
# stable manifold. Free variables: the two manifold phases (tau_u, tau_s) and the
# inter-system phase theta (the relative orientation of the EM line vs the SE line).
#
# Section default X0 = 1.5e6 km: the SE-L1/L2 orbits sit at ~1.6 Mkm from Earth and
# the EM-L2 unstable manifold reaches out past ~1.7 Mkm in inertial X, so both
# manifolds cross this plane (the #316 sunward-neck probe distance).

_PATCH_X0_KM_DEFAULT = 1.5e6


@dataclass(frozen=True)
class CrossConnection:
    """An EM-orbit -> SE-orbit connection matched on an inertial patch section.

    ``residual`` is the inertial POSITION gap (km) between the EM unstable manifold
    and the SE stable manifold at the section; ``patch_dv_kms`` is the inertial
    VELOCITY gap there (km/s, the impulsive patch maneuver). ``patch_state_inertial``
    is the 6-state (km, km/s) at the converged EM-side crossing.
    """

    label_from: str
    label_to: str
    c_em: float
    c_se: float
    theta: float
    tau_u: float
    tau_s: float
    patch_state_inertial: NDArray[np.float64]
    patch_dv_kms: float
    residual: float  # inertial position gap, km
    converged: bool
    n_iter: int
    notes: str = ""
    transit_time: float = float("nan")  # seconds, |t_u| + |t_s| at the converged crossing


@dataclass(frozen=True)
class Periodic3DNode:
    """A 3D periodic orbit serving as a cycle node.

    ``state0`` is the full 6-vector IC ``(x0, y0, z0, vx0, vy0, vz0)``;
    ``unstable_eigvec`` / ``stable_eigvec`` are the full 6-vectors of the
    Floquet saddle pair, used to seed the manifolds.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    unstable_eigvec: NDArray[np.float64]  # 6-vector
    stable_eigvec: NDArray[np.float64]  # 6-vector
    converged: bool


def _seed_on_manifold_3d(
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode | Periodic3DNode,
    *,
    tau: float,
    direction: str,
    branch: int,
    epsilon: float,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64]:
    """Manifold seed at phase ``tau`` along a planar or 3D ``node``."""
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    if branch not in (+1, -1):
        raise ValueError(f"branch must be +1 or -1; got {branch!r}")
    tau = float(tau) % float(node.period)
    if tau <= 0.0:
        state_tau = np.asarray(node.state0, float)
        phi = np.eye(6)
    else:
        arc = cr3bp.propagate(system, node.state0, tau, with_stm=True, rtol=rtol, atol=atol)
        assert arc.stm is not None
        state_tau = arc.state_f
        phi = arc.stm
    v = node.unstable_eigvec if direction == "unstable" else node.stable_eigvec
    v6 = v if len(v) == 6 else np.array([v[0], v[1], 0.0, v[2], v[3], 0.0], dtype=np.float64)
    v_tau = phi @ v6
    n = float(np.linalg.norm(v_tau))
    if n > 0.0:
        v_tau = v_tau / n
    return (state_tau + float(branch) * float(epsilon) * v_tau).astype(np.float64)


def _manifold_crossing_timed(
    bridge: FrameBridge,
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode | Periodic3DNode,
    *,
    side: str,
    tau: float,
    direction: str,
    branch: int,
    theta: float,
    x0_km: float,
    epsilon: float,
    max_time: float,
    n_samples: int = 1200,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    method: str = "DOP853",
) -> tuple[NDArray[np.float64], float] | None:
    """Core crossing detector: ``(inertial_state, t_cross)`` or ``None``.

    ``t_cross`` is the SIGNED integration time at the crossing (positive for unstable
    forward, negative for stable backward); its magnitude is the manifold flight time
    to/from the patch, needed for the #411 time-consistent cycle theta closure. See
    :func:`_manifold_inertial_at_section` for the geometry.
    """
    if side not in {"em", "se"}:
        raise ValueError(f"side must be 'em' or 'se'; got {side!r}")
    try:
        seed = _seed_on_manifold_3d(
            system, node, tau=tau, direction=direction, branch=branch, epsilon=epsilon
        )
    except (RuntimeError, ValueError):
        return None  # orbit not valid at this energy/phase (e.g. singular STM)
    horizon = abs(float(max_time))
    t_span = (0.0, horizon) if direction == "unstable" else (0.0, -horizon)
    t_eval = np.linspace(t_span[0], t_span[1], n_samples)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        np.asarray(seed, float),
        args=(system.mu,),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
    )
    if not sol.success:
        return None

    def _to_inertial(state: NDArray[np.float64]) -> NDArray[np.float64]:
        if side == "em":
            return bridge.em_rot_to_inertial(state, theta=theta)
        return bridge.se_rot_to_inertial(state, theta=theta)

    prev_in: NDArray[np.float64] | None = None
    prev_g = 0.0
    prev_t = 0.0
    for k in range(sol.y.shape[1]):
        inert = _to_inertial(sol.y[:, k])
        g = float(inert[0] - x0_km)
        if prev_in is not None and prev_g != g and (prev_g <= 0.0 <= g or g <= 0.0 <= prev_g):
            frac = prev_g / (prev_g - g)  # linear interpolation of the crossing
            crossing = (prev_in + frac * (inert - prev_in)).astype(np.float64)
            t_cross = prev_t + frac * (float(sol.t[k]) - prev_t)
            return crossing, float(t_cross)
        prev_in, prev_g, prev_t = inert, g, float(sol.t[k])
    return None


def _manifold_inertial_at_section(
    bridge: FrameBridge,
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode | Periodic3DNode,
    *,
    side: str,
    tau: float,
    direction: str,
    branch: int,
    theta: float,
    x0_km: float,
    epsilon: float,
    max_time: float,
    n_samples: int = 1200,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    method: str = "DOP853",
) -> NDArray[np.float64] | None:
    """Inertial 6-state where ``node``'s manifold first crosses {x_inertial = x0_km}.

    Seeds the manifold in-system, propagates forward (unstable) or backward (stable)
    over a bounded horizon, samples densely, transforms each sample to the common
    inertial frame (EM side via ``em_rot_to_inertial``, SE side via
    ``se_rot_to_inertial``), and returns the first sign-changing crossing of the
    section plane (linearly refined). Returns ``None`` if the manifold never reaches
    the plane within ``max_time`` (bounded — never hangs, never fabricates a crossing).
    Thin state-only wrapper over :func:`_manifold_crossing_timed`.
    """
    res = _manifold_crossing_timed(
        bridge,
        system,
        node,
        side=side,
        tau=tau,
        direction=direction,
        branch=branch,
        theta=theta,
        x0_km=x0_km,
        epsilon=epsilon,
        max_time=max_time,
        n_samples=n_samples,
        rtol=rtol,
        atol=atol,
        method=method,
    )
    return None if res is None else res[0]


def _cross_residual(
    bridge: FrameBridge,
    orbit_from: LyapunovNode | Periodic3DNode,
    orbit_to: LyapunovNode | Periodic3DNode,
    *,
    from_side: str,
    to_side: str,
    tau_u: float,
    tau_s: float,
    theta: float,
    branch_u: int,
    branch_s: int,
    x0_km: float,
    epsilon: float,
    max_time_u: float,
    max_time_s: float,
) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None, NDArray[np.float64] | None]:
    """Return ``(pos_gap_km(3), patch_state_inertial(6), vel_gap_kms(3))`` or Nones.

    ``pos_gap`` / ``vel_gap`` are the inertial position/velocity differences of the
    ``orbit_from`` UNSTABLE manifold and the ``orbit_to`` STABLE manifold at the
    section; ``patch_state_inertial`` is the from-side crossing 6-state. Each manifold
    is propagated in ITS OWN system (``from_side``/``to_side`` in {"em","se"}), so the
    matcher is direction-agnostic (EM->SE forward AND SE->EM return both correct).
    """
    sys_of = {"em": bridge.em, "se": bridge.se}
    p_u = _manifold_inertial_at_section(
        bridge,
        sys_of[from_side],
        orbit_from,
        side=from_side,
        tau=tau_u,
        direction="unstable",
        branch=branch_u,
        theta=theta,
        x0_km=x0_km,
        epsilon=epsilon,
        max_time=max_time_u,
    )
    p_s = _manifold_inertial_at_section(
        bridge,
        sys_of[to_side],
        orbit_to,
        side=to_side,
        tau=tau_s,
        direction="stable",
        branch=branch_s,
        theta=theta,
        x0_km=x0_km,
        epsilon=epsilon,
        max_time=max_time_s,
    )
    if p_u is None or p_s is None:
        return None, None, None
    pos_gap = (p_u[:3] - p_s[:3]).astype(np.float64)
    vel_gap = (p_u[3:] - p_s[3:]).astype(np.float64)
    return pos_gap, p_u.astype(np.float64), vel_gap


def _scan_cross_starts(
    resid: Callable[
        [float, float, float],
        tuple[NDArray[np.float64] | None, NDArray[np.float64] | None, NDArray[np.float64] | None],
    ],
    period_u: float,
    period_s: float,
    *,
    n_theta: int,
    n_tau: int,
    pos_tol_km: float = 2.0e5,
) -> tuple[float, float, float]:
    """Coarse 3-D grid search that prioritizes low velocity gap (patch dV).

    Among cells that co-reach the section with a position gap below pos_tol_km,
    selects the one that minimizes the velocity gap. If no cells are below
    pos_tol_km, falls back to the one with the smallest position gap.
    """
    thetas = [2.0 * math.pi * (i + 0.5) / n_theta for i in range(n_theta)]
    us = [period_u * (i + 0.5) / n_tau for i in range(n_tau)]
    ss = [period_s * (j + 0.5) / n_tau for j in range(n_tau)]

    best_tu, best_ts, best_th = 0.5 * period_u, 0.5 * period_s, math.pi
    candidates = []

    for th in thetas:
        for tu in us:
            for ts in ss:
                pos_gap, _, vel_gap = resid(tu, ts, th)
                if pos_gap is None or vel_gap is None:
                    continue
                rn_pos = float(np.linalg.norm(pos_gap))
                rn_dv = float(np.linalg.norm(vel_gap))
                candidates.append((rn_pos, rn_dv, tu, ts, th))

    if not candidates:
        return best_tu, best_ts, best_th

    good_pos = [c for c in candidates if c[0] < pos_tol_km]
    if good_pos:
        good_pos.sort(key=lambda x: x[1])
        _, _, best_tu, best_ts, best_th = good_pos[0]
    else:
        candidates.sort(key=lambda x: x[0])
        _, _, best_tu, best_ts, best_th = candidates[0]

    return best_tu, best_ts, best_th


def correct_cross_connection(
    bridge: FrameBridge,
    orbit_from: LyapunovNode | Periodic3DNode,
    orbit_to: LyapunovNode | Periodic3DNode,
    *,
    label_from: str,
    label_to: str,
    k: int = 1,
    epsilon: float = 1e-6,
    branch_u: int = +1,
    branch_s: int = -1,
    x0_km: float = _PATCH_X0_KM_DEFAULT,
    tol_km: float = 1e2,
    max_iter: int = 40,
    fd_step: float = 1e-6,
    max_time_factor: float = 8.0,
    scan_n: int = 12,
    scan_n_tau: int = 4,
) -> CrossConnection:
    """Match Wu(orbit_from, EM) to Ws(orbit_to, SE) on inertial plane {x=x0_km}.

    Residual = inertial POSITION gap (km); patch ΔV = inertial VELOCITY gap (km/s).
    Free vars (tau_u, tau_s, theta). A coarse 3-D scan (theta densest, ``scan_n``
    cells) seeds a 3x3 FD-Jacobian Newton with backtracking line-search. theta is
    wrapped to [0, 2π), the phases to their periods. NEVER raises for "no connection"
    — returns ``converged=False`` with a diagnostic note and always sets
    ``patch_dv_kms`` (inf if no crossing pair was ever found).

    The EM manifold must travel ~1.5 Mkm sunward to reach the SE neighbourhood; if it
    cannot reach the section at this energy, no connection exists there (a real
    Phase-A non-closure, surfaced as ``converged=False``).
    """
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1; got {max_iter}")
    # Direction-aware: each manifold is propagated in its OWN system. Infer the side of
    # each node from its label ("EM-..." -> em, "SE-..." -> se) so the matcher is correct
    # for BOTH the EM->SE forward leg AND the SE->EM return leg (the prior code hardwired
    # from=EM/to=SE, silently mis-propagating the return — #411).
    from_side = "em" if orbit_from.label.upper().startswith("EM") else "se"
    to_side = "em" if orbit_to.label.upper().startswith("EM") else "se"
    max_time_u = max_time_factor * orbit_from.period
    max_time_s = max_time_factor * orbit_to.period
    # Record c_em / c_se by SIDE (not by from/to position) so the result is consistent
    # in both directions.
    c_em = orbit_from.jacobi if from_side == "em" else orbit_to.jacobi
    c_se = orbit_from.jacobi if from_side == "se" else orbit_to.jacobi

    def _resid(
        tu: float, ts: float, th: float
    ) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None, NDArray[np.float64] | None]:
        return _cross_residual(
            bridge,
            orbit_from,
            orbit_to,
            from_side=from_side,
            to_side=to_side,
            tau_u=tu,
            tau_s=ts,
            theta=th,
            branch_u=branch_u,
            branch_s=branch_s,
            x0_km=x0_km,
            epsilon=epsilon,
            max_time_u=max_time_u,
            max_time_s=max_time_s,
        )

    tau_u, tau_s, theta = _scan_cross_starts(
        _resid, orbit_from.period, orbit_to.period, n_theta=scan_n, n_tau=scan_n_tau
    )

    def _fail(notes: str) -> CrossConnection:
        return CrossConnection(
            label_from=label_from,
            label_to=label_to,
            c_em=c_em,
            c_se=c_se,
            theta=theta,
            tau_u=tau_u,
            tau_s=tau_s,
            patch_state_inertial=np.zeros(6),
            patch_dv_kms=float("inf"),
            residual=float("inf"),
            converged=False,
            n_iter=0,
            notes=notes,
        )

    pos_gap, patch_state, vel_gap = _resid(tau_u, tau_s, theta)
    if pos_gap is None:
        return _fail("manifolds never co-reach the inertial patch section")

    period_u, period_s = orbit_from.period, orbit_to.period
    two_pi = 2.0 * math.pi
    n_iter = 0
    for n_iter in range(1, max_iter + 1):  # noqa: B007  (n_iter reported in result)
        rn = float(np.linalg.norm(pos_gap))
        if rn < tol_km:
            break
        jac = np.zeros((3, 3), dtype=np.float64)
        steps = [(fd_step, 0.0, 0.0), (0.0, fd_step, 0.0), (0.0, 0.0, fd_step)]
        ok = True
        for j, (du, ds, dth) in enumerate(steps):
            gp, _, _ = _resid(tau_u + du, tau_s + ds, theta + dth)
            if gp is None:
                ok = False
                break
            jac[:, j] = (gp - pos_gap) / fd_step
        if not ok:
            break  # FD probe left a section branch; stop and report best-so-far
        try:
            step = np.linalg.solve(jac, -pos_gap)
        except np.linalg.LinAlgError:
            step, *_ = np.linalg.lstsq(jac, -pos_gap, rcond=None)
        alpha = 1.0
        improved = False
        for _ in range(25):
            tu_t = (tau_u + alpha * float(step[0])) % period_u
            ts_t = (tau_s + alpha * float(step[1])) % period_s
            th_t = (theta + alpha * float(step[2])) % two_pi
            gp_t, ps_t, vg_t = _resid(tu_t, ts_t, th_t)
            if gp_t is not None and float(np.linalg.norm(gp_t)) < rn:
                tau_u, tau_s, theta = tu_t, ts_t, th_t
                pos_gap, patch_state, vel_gap = gp_t, ps_t, vg_t
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break

    final_rn = float(np.linalg.norm(pos_gap)) if pos_gap is not None else float("inf")
    patch_dv = float(np.linalg.norm(vel_gap)) if vel_gap is not None else float("inf")
    converged = pos_gap is not None and final_rn < tol_km

    # Transit time at the converged crossing: |t_u| (unstable flight to patch) + |t_s|
    # (stable flight from patch), in seconds. Needed for the #411 time-consistent cycle
    # theta closure. Computed only on convergence (one extra propagation pair); nan otherwise.
    transit_time = float("nan")
    if converged:
        sys_of = {"em": bridge.em, "se": bridge.se}
        cu = _manifold_crossing_timed(
            bridge,
            sys_of[from_side],
            orbit_from,
            side=from_side,
            tau=tau_u,
            direction="unstable",
            branch=branch_u,
            theta=theta,
            x0_km=x0_km,
            epsilon=epsilon,
            max_time=max_time_u,
        )
        cs = _manifold_crossing_timed(
            bridge,
            sys_of[to_side],
            orbit_to,
            side=to_side,
            tau=tau_s,
            direction="stable",
            branch=branch_s,
            theta=theta,
            x0_km=x0_km,
            epsilon=epsilon,
            max_time=max_time_s,
        )
        if cu is not None and cs is not None:
            # Crossing times are each in the manifold's OWN nondim time; scale each by its
            # system's t_s to seconds (cu is from_side, cs is to_side -- possibly different
            # systems with different t_s).
            transit_time = abs(cu[1]) * sys_of[from_side].t_s + abs(cs[1]) * sys_of[to_side].t_s

    return CrossConnection(
        label_from=label_from,
        label_to=label_to,
        c_em=c_em,
        c_se=c_se,
        theta=theta,
        tau_u=tau_u,
        tau_s=tau_s,
        patch_state_inertial=patch_state if patch_state is not None else np.zeros(6),
        patch_dv_kms=patch_dv,
        residual=final_rn,
        converged=converged,
        n_iter=n_iter,
        notes="" if converged else "did not reach tol_km",
        transit_time=transit_time,
    )


# --- Task 4: bounded closure search + CrossCycle ----------------------------
#
# A cross-system CYCLE chains a forward EM->SE connection with a SE->EM RETURN
# connection so the trajectory comes back to the starting EM orbit. Spatial closure
# (both inertial position legs match) is necessary but not sufficient: the relative
# phase ``theta`` between the SE line and the EM line must ALSO return commensurately
# (periodic-up-to-rotation). theta_closure_residual measures the phase mismatch
# |theta_return - theta_forward| reduced to [0, pi]; closure requires it small.
#
# This simultaneous spatial+phase closure is genuinely hard. The #316 survey notes
# the natural cross-system closure is the ~19yr Metonic (235:19) commensurability, far
# outside any single-revolution bounded grid. A non-closing grid is therefore an
# honest CLEAN NEGATIVE (an acceptable Phase-A result), surfaced as ``closed=False``
# with a diagnostic note — NEVER fabricated closure or loosened thresholds.

# Per-libration Lyapunov seeds (x0 guess, ydot0 sign). EM/SE families need the -1
# branch (the +1 sign collapses onto a different, shorter-period member). L2 seeds are
# the Task-3 working values; L1 seeds are the sunward-side equivalents.
# L2 seeds are the Task-3 verified-converging values (EM-L2 x0=1.18 -> x0_corr~1.182
# T~3.42; SE-L2 x0=1.009 -> x0_corr~1.0108 T~3.07). The nearby x0=1.155 / 1.0101 wander
# off the branch (Newton lands on a different-period member or fails). L1 seeds are the
# sunward-side equivalents.
_EM_SEEDS: dict[str, tuple[float, float]] = {
    "EM-L1": (0.85, -1.0),
    "EM-L2": (1.18, -1.0),
}
_SE_SEEDS: dict[str, tuple[float, float]] = {
    "SE-L1": (0.9893, -1.0),
    "SE-L2": (1.009, -1.0),
}
# Period guesses by family (nondim); only used as Newton seeds.
_EM_PERIOD_GUESS = 3.4
_SE_PERIOD_GUESS = 3.06


@dataclass(frozen=True)
class CrossCycle:
    """A candidate cross-system cycle: EM->SE forward leg chained to a SE->EM return.

    ``closed`` is True only when BOTH legs converged AND ``theta_closure_residual`` is
    below the commensurability tolerance (periodic-up-to-rotation). ``max_leg_residual``
    is the larger of the two legs' inertial position gaps (km). ``independent_residual``
    is reserved for an independent cross-check (Task 5; nan here). When not closed,
    ``notes`` records what was observed (which leg failed / best theta residual).
    """

    connections: list[CrossConnection]
    c_em: float
    c_se: float
    libration_pair: tuple[str, str]
    theta_closure_residual: float
    closed: bool
    max_leg_residual: float
    independent_residual: float
    notes: str = ""


def _theta_gap(theta_a: float, theta_b: float) -> float:
    """|theta_a - theta_b| reduced to [0, pi] (relative-phase commensurability gap)."""
    d = (theta_a - theta_b) % (2.0 * math.pi)
    return float(min(d, 2.0 * math.pi - d))


def _wrap_pi(x: float) -> float:
    """Wrap an angle (rad) to [-pi, pi) (exactly +pi maps to -pi; same point)."""
    return float((x + math.pi) % (2.0 * math.pi) - math.pi)


def theta_commensurability(
    gap_rad: float,
    dtheta_em_rad: float,
    dtheta_se_rad: float,
    *,
    n_max: int = 60,
    tol_rad: float = 5.0e-2,
) -> tuple[int, int, float, bool]:
    """Smallest-residual non-negative (n_em, n_se) nulling the cross-cycle phase gap.

    Necessary-condition feasibility check for closing a cross-system cycle by
    multi-revolution phasing (#411): adding ``n_em`` revolutions on the EM orbit and
    ``n_se`` on the SE orbit advances the relative phase by ``n_em*dtheta_em +
    n_se*dtheta_se``; closure needs ``gap + n_em*dtheta_em + n_se*dtheta_se ≡ 0
    (mod 2π)``. Grid-searches ``0 <= n_em, n_se <= n_max`` and returns
    ``(n_em, n_se, residual_rad, feasible)`` for the minimal mod-2π residual, with
    ``feasible = residual < tol_rad``. This is NECESSARY-not-sufficient: the
    multi-rev connections themselves must still re-converge (sufficiency). It asserts
    NO orbit — only whether a commensurate phase closure exists at tractable counts.
    """
    two_pi = 2.0 * math.pi
    best_n_em, best_n_se, best_res = 0, 0, float("inf")
    for n_em in range(n_max + 1):
        partial = gap_rad + n_em * dtheta_em_rad
        for n_se in range(n_max + 1):
            val = (partial + n_se * dtheta_se_rad) % two_pi
            res = min(val, two_pi - val)
            if res < best_res:
                best_n_em, best_n_se, best_res = n_em, n_se, res
    return best_n_em, best_n_se, float(best_res), bool(best_res < tol_rad)


def _real_unit_6d(v: NDArray[np.complex128]) -> NDArray[np.float64]:
    vr = np.real(v)
    n = float(np.linalg.norm(vr))
    if n < 1e-14:
        vr = np.real(v) + np.imag(v)
        n = float(np.linalg.norm(vr))
    if n < 1e-14:
        raise ValueError("eigenvector is numerically zero")
    return (vr / n).astype(np.float64)


def _3d_floquet_pair(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, NDArray[np.float64], float, NDArray[np.float64]]:
    """Return (|lam_u|, v_u, |lam_s|, v_s) for the 3D monodromy saddle pair.

    Integrates the 6x6 STM over one period, computes eigenvalues, and returns the
    largest-magnitude eigenvalue (> 1) and smallest-magnitude eigenvalue (< 1)
    with their real-normalised 6D right eigenvectors.
    """
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    eigs, vecs = np.linalg.eig(arc.stm)
    mags = np.abs(eigs)

    # Filter out complex conjugate pairs that are close to the unit circle
    hyperbolic_indices = []
    for i, mag in enumerate(mags):
        if abs(mag - 1.0) > 1e-3:
            hyperbolic_indices.append(i)

    if len(hyperbolic_indices) >= 2:
        idx_u = max(hyperbolic_indices, key=lambda i: mags[i])
        idx_s = min(hyperbolic_indices, key=lambda i: mags[i])
    else:
        idx_u = int(np.argmax(mags))
        idx_s = int(np.argmin(mags))

    lam_u = eigs[idx_u]
    lam_s = eigs[idx_s]

    v_u = _real_unit_6d(vecs[:, idx_u])
    v_s = _real_unit_6d(vecs[:, idx_s])

    return float(np.abs(lam_u)), v_u, float(np.abs(lam_s)), v_s


def _build_em_node(em: cr3bp.CR3BPSystem, label: str, c_em: float) -> LyapunovNode:
    x0, sign = _EM_SEEDS.get(label, (1.155, -1.0))
    return LyapunovNode.from_libration(
        em, x0_guess=x0, jacobi=c_em, period_guess=_EM_PERIOD_GUESS, label=label, ydot0_sign=sign
    )


def _build_se_node(se: cr3bp.CR3BPSystem, label: str, c_se: float) -> LyapunovNode:
    x0, sign = _SE_SEEDS.get(label, (1.0101, -1.0))
    return LyapunovNode.from_libration(
        se, x0_guess=x0, jacobi=c_se, period_guess=_SE_PERIOD_GUESS, label=label, ydot0_sign=sign
    )


def _build_em_node_3d(em: cr3bp.CR3BPSystem, label: str, c_em: float, z0: float) -> Periodic3DNode:
    planar = _build_em_node(em, label, c_em)
    if abs(z0) < 1e-9:
        _, v_u, _, v_s = _3d_floquet_pair(em, planar.state0, planar.period)
        return Periodic3DNode(
            label=label,
            state0=planar.state0,
            period=planar.period,
            jacobi=planar.jacobi,
            unstable_eigvec=v_u,
            stable_eigvec=v_s,
            converged=planar.converged,
        )

    state0_guess = np.array(
        [planar.state0[0], 0.0, z0, 0.0, planar.state0[4], 0.0], dtype=np.float64
    )
    orbit = correct_general_periodic_3d(
        em,
        state0_guess,
        planar.period,
        free_vars=(0, 4, 6),
        residual_indices=(1, 3, 5),
        is_half_period_residual=True,
    )
    jac = cr3bp.jacobi_constant(orbit.state0, em.mu)
    _, v_u, _, v_s = _3d_floquet_pair(em, orbit.state0, orbit.T_TU)

    v_u_planar = np.array(
        [
            planar.unstable_eigvec[0],
            planar.unstable_eigvec[1],
            0.0,
            planar.unstable_eigvec[2],
            planar.unstable_eigvec[3],
            0.0,
        ],
        dtype=np.float64,
    )
    if np.dot(v_u, v_u_planar) < 0.0:
        v_u = -v_u
    v_s_planar = np.array(
        [
            planar.stable_eigvec[0],
            planar.stable_eigvec[1],
            0.0,
            planar.stable_eigvec[2],
            planar.stable_eigvec[3],
            0.0,
        ],
        dtype=np.float64,
    )
    if np.dot(v_s, v_s_planar) < 0.0:
        v_s = -v_s

    return Periodic3DNode(
        label=label,
        state0=orbit.state0,
        period=orbit.T_TU,
        jacobi=jac,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=orbit.converged,
    )


def _build_se_node_3d(se: cr3bp.CR3BPSystem, label: str, c_se: float, z0: float) -> Periodic3DNode:
    planar = _build_se_node(se, label, c_se)
    if abs(z0) < 1e-9:
        _, v_u, _, v_s = _3d_floquet_pair(se, planar.state0, planar.period)
        return Periodic3DNode(
            label=label,
            state0=planar.state0,
            period=planar.period,
            jacobi=planar.jacobi,
            unstable_eigvec=v_u,
            stable_eigvec=v_s,
            converged=planar.converged,
        )

    state0_guess = np.array(
        [planar.state0[0], 0.0, z0, 0.0, planar.state0[4], 0.0], dtype=np.float64
    )
    orbit = correct_general_periodic_3d(
        se,
        state0_guess,
        planar.period,
        free_vars=(0, 4, 6),
        residual_indices=(1, 3, 5),
        is_half_period_residual=True,
    )
    jac = cr3bp.jacobi_constant(orbit.state0, se.mu)
    _, v_u, _, v_s = _3d_floquet_pair(se, orbit.state0, orbit.T_TU)

    v_u_planar = np.array(
        [
            planar.unstable_eigvec[0],
            planar.unstable_eigvec[1],
            0.0,
            planar.unstable_eigvec[2],
            planar.unstable_eigvec[3],
            0.0,
        ],
        dtype=np.float64,
    )
    if np.dot(v_u, v_u_planar) < 0.0:
        v_u = -v_u
    v_s_planar = np.array(
        [
            planar.stable_eigvec[0],
            planar.stable_eigvec[1],
            0.0,
            planar.stable_eigvec[2],
            planar.stable_eigvec[3],
            0.0,
        ],
        dtype=np.float64,
    )
    if np.dot(v_s, v_s_planar) < 0.0:
        v_s = -v_s

    return Periodic3DNode(
        label=label,
        state0=orbit.state0,
        period=orbit.T_TU,
        jacobi=jac,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=orbit.converged,
    )


def crosscheck_cross_cycle(
    bridge: FrameBridge,
    cycle: CrossCycle,
    *,
    method: str = "Radau",
    rtol: float = 1e-11,
    atol: float = 1e-11,
    epsilon: float = 1e-6,
    max_time_factor: float = 8.0,
) -> CrossCycle:
    """Re-derive each converged leg's patch-state with an independent integrator.

    For each connection in ``cycle.connections`` that ``converged``, re-runs
    ``_manifold_inertial_at_section`` at the stored ``(theta, tau_u, tau_s, branch)``
    using ``method`` (default Radau; the corrector used DOP853) and compares the
    re-derived inertial patch POSITION to the stored ``conn.patch_state_inertial[:3]``.
    The maximum position disagreement (km) across all re-derived legs is
    ``independent_residual``.  A leg that fails to re-derive contributes ``inf`` (a
    real failure, surfaced — never silently dropped).  If NO leg converged (the
    Task-4 clean-negative case) there is nothing to re-derive and the result is
    ``inf`` (finite, non-nan — the test requires not-nan).

    Returns a new ``CrossCycle`` (dataclasses.replace) with ``independent_residual``
    filled; all other fields are unchanged.  Mirrors ``heteroclinic_cycle.crosscheck_cycle``.
    """
    import dataclasses

    em, se = bridge.em, bridge.se
    worst: float = 0.0
    any_converged = False

    for conn in cycle.connections:
        if not conn.converged:
            continue
        any_converged = True
        # Determine which system and direction this leg belongs to.
        # Forward leg: label_from is an EM label (EM-L*) → unstable from EM.
        # Return leg: label_from is a SE label (SE-L*) → unstable from SE.
        if conn.label_from.startswith("EM"):
            system = em
            side = "em"
            label = conn.label_from
            try:
                node = _build_em_node(em, label, conn.c_em)
            except Exception:
                worst = float("inf")
                continue
            if not node.converged:
                worst = float("inf")
                continue
            max_time = max_time_factor * node.period
            recheck = _manifold_inertial_at_section(
                bridge,
                system,
                node,
                side=side,
                tau=conn.tau_u,
                direction="unstable",
                branch=+1,  # forward leg uses branch_u=+1 (corrector default)
                theta=conn.theta,
                x0_km=float(conn.patch_state_inertial[0])
                if np.any(conn.patch_state_inertial[:3] != 0.0)
                else _PATCH_X0_KM_DEFAULT,
                epsilon=epsilon,
                max_time=max_time,
                method=method,
                rtol=rtol,
                atol=atol,
            )
        else:
            # Return leg: SE unstable.
            system = se
            side = "se"
            label = conn.label_from
            try:
                node = _build_se_node(se, label, conn.c_se)
            except Exception:
                worst = float("inf")
                continue
            if not node.converged:
                worst = float("inf")
                continue
            max_time = max_time_factor * node.period
            recheck = _manifold_inertial_at_section(
                bridge,
                system,
                node,
                side=side,
                tau=conn.tau_u,
                direction="unstable",
                branch=+1,
                theta=conn.theta,
                x0_km=float(conn.patch_state_inertial[0])
                if np.any(conn.patch_state_inertial[:3] != 0.0)
                else _PATCH_X0_KM_DEFAULT,
                epsilon=epsilon,
                max_time=max_time,
                method=method,
                rtol=rtol,
                atol=atol,
            )

        if recheck is None:
            worst = float("inf")
            continue
        pos_gap = float(np.linalg.norm(recheck[:3] - conn.patch_state_inertial[:3]))
        worst = max(worst, pos_gap)

    if not any_converged:
        worst = float("inf")

    return dataclasses.replace(cycle, independent_residual=worst)


def crosscheck_cross_cycle_3d(
    bridge: FrameBridge,
    cycle: CrossCycleClosure,
    z_em: float,
    z_se: float,
    *,
    method: str = "Radau",
    rtol: float = 1e-11,
    atol: float = 1e-11,
    epsilon: float = 1e-6,
    max_time_factor: float = 8.0,
) -> float:
    """Re-derive each converged 3D leg's patch-state with an independent integrator.

    Returns the maximum position disagreement (km) across all re-derived legs.
    """
    em, se = bridge.em, bridge.se
    worst: float = 0.0
    any_converged = False

    for conn in [cycle.forward, cycle.ret]:
        if not conn.converged:
            continue
        any_converged = True
        if conn.label_from.startswith("EM"):
            system = em
            side = "em"
            label = conn.label_from
            try:
                node = _build_em_node_3d(em, label, conn.c_em, z_em)
            except Exception:
                worst = float("inf")
                continue
            if not node.converged:
                worst = float("inf")
                continue
            max_time = max_time_factor * node.period
            recheck = _manifold_inertial_at_section(
                bridge,
                system,
                node,
                side=side,
                tau=conn.tau_u,
                direction="unstable",
                branch=+1,
                theta=conn.theta,
                x0_km=float(conn.patch_state_inertial[0])
                if np.any(conn.patch_state_inertial[:3] != 0.0)
                else _PATCH_X0_KM_DEFAULT,
                epsilon=epsilon,
                max_time=max_time,
                method=method,
                rtol=rtol,
                atol=atol,
            )
        else:
            system = se
            side = "se"
            label = conn.label_from
            try:
                node = _build_se_node_3d(se, label, conn.c_se, z_se)
            except Exception:
                worst = float("inf")
                continue
            if not node.converged:
                worst = float("inf")
                continue
            max_time = max_time_factor * node.period
            recheck = _manifold_inertial_at_section(
                bridge,
                system,
                node,
                side=side,
                tau=conn.tau_u,
                direction="unstable",
                branch=+1,
                theta=conn.theta,
                x0_km=float(conn.patch_state_inertial[0])
                if np.any(conn.patch_state_inertial[:3] != 0.0)
                else _PATCH_X0_KM_DEFAULT,
                epsilon=epsilon,
                max_time=max_time,
                method=method,
                rtol=rtol,
                atol=atol,
            )

        if recheck is None:
            worst = float("inf")
            continue
        pos_gap = float(np.linalg.norm(recheck[:3] - conn.patch_state_inertial[:3]))
        worst = max(worst, pos_gap)

    if not any_converged:
        worst = float("inf")

    return worst


def search_cross_cycle(
    bridge: FrameBridge,
    *,
    c_em_grid: tuple[float, ...],
    c_se_grid: tuple[float, ...],
    libration_pairs: tuple[tuple[str, str], ...],
    max_attempts: int = 2,
    theta_tol: float = 1e-2,
    tol_km: float = 1e2,
    return_scan_n: int = 4,
    return_scan_n_tau: int = 2,
    return_max_time_factor: float = 4.0,
    **conn_kwargs: object,
) -> list[CrossCycle]:
    """Bounded energy x libration grid search for a closed cross-system cycle.

    For each ``(c_em, c_se, (em_lib, se_lib))`` grid point: build the two Lyapunov
    nodes, correct the forward EM->SE connection and the (harder) SE->EM return
    connection, then test spatial AND phase closure. The return leg may need a
    different stable/unstable branch than the forward leg (like #314's L2->L1 return);
    up to ``max_attempts`` branch/k variations are tried and the best (least position
    gap) is kept.

    BOUNDED COST: the return leg (SE unstable -> EM stable) integrates an SE-system
    manifold over many SE periods (each SE time unit ~58 days) and is ~25x costlier
    than the forward leg, so it is run on a coarser scan and a shorter horizon
    (``return_scan_n`` / ``return_scan_n_tau`` / ``return_max_time_factor``). This is
    a deliberately bounded Phase-A search, not an exhaustive one: a coarser grid that
    finds no crossing is still an honest non-closure (the manifolds genuinely never
    co-reach the section at this horizon), not an artifact to hide.

    Returns one ``CrossCycle`` per grid point. NEVER raises for non-closure — a grid
    that does not close is an honest CLEAN NEGATIVE (closed=False with a diagnostic
    note). Thresholds are NEVER loosened to manufacture a closed=True.
    """
    em, se = bridge.em, bridge.se
    # Return-leg branch/k variations to try (forward defaults are branch_u=+1, branch_s=-1).
    return_variants: list[dict[str, int]] = [
        {"branch_u": +1, "branch_s": -1, "k": 1},
        {"branch_u": -1, "branch_s": +1, "k": 1},
        {"branch_u": +1, "branch_s": +1, "k": 1},
        {"branch_u": -1, "branch_s": -1, "k": 1},
    ][: max(1, int(max_attempts))]

    results: list[CrossCycle] = []
    for c_em in c_em_grid:
        for c_se in c_se_grid:
            for em_lib, se_lib in libration_pairs:
                pair = (em_lib, se_lib)
                try:
                    em_node = _build_em_node(em, em_lib, c_em)
                    se_node = _build_se_node(se, se_lib, c_se)
                except Exception as exc:  # report, never abort the grid
                    results.append(
                        CrossCycle(
                            connections=[],
                            c_em=c_em,
                            c_se=c_se,
                            libration_pair=pair,
                            theta_closure_residual=float("nan"),
                            closed=False,
                            max_leg_residual=float("inf"),
                            independent_residual=float("nan"),
                            notes=f"node construction failed: {exc}",
                        )
                    )
                    continue
                if not (em_node.converged and se_node.converged):
                    results.append(
                        CrossCycle(
                            connections=[],
                            c_em=c_em,
                            c_se=c_se,
                            libration_pair=pair,
                            theta_closure_residual=float("nan"),
                            closed=False,
                            max_leg_residual=float("inf"),
                            independent_residual=float("nan"),
                            notes=(
                                f"Lyapunov node did not converge "
                                f"(EM={em_node.converged}, SE={se_node.converged})"
                            ),
                        )
                    )
                    continue

                # Forward leg: EM unstable -> SE stable (Task-3 working direction).
                fwd = correct_cross_connection(
                    bridge,
                    em_node,
                    se_node,
                    label_from=em_lib,
                    label_to=se_lib,
                    tol_km=tol_km,
                    **conn_kwargs,  # type: ignore[arg-type]
                )

                # Return leg: SE unstable -> EM stable (harder). Try branch/k variants,
                # keep the one with the smallest inertial position gap.
                ret: CrossConnection | None = None
                for variant in return_variants:
                    cand = correct_cross_connection(
                        bridge,
                        se_node,
                        em_node,
                        label_from=se_lib,
                        label_to=em_lib,
                        tol_km=tol_km,
                        branch_u=variant["branch_u"],
                        branch_s=variant["branch_s"],
                        k=variant["k"],
                        scan_n=return_scan_n,
                        scan_n_tau=return_scan_n_tau,
                        max_time_factor=return_max_time_factor,
                        **conn_kwargs,  # type: ignore[arg-type]
                    )
                    if ret is None or cand.residual < ret.residual:
                        ret = cand
                    if cand.converged:
                        break
                assert ret is not None  # max_attempts >= 1 guarantees one candidate

                connections = [fwd, ret]
                theta_res = _theta_gap(ret.theta, fwd.theta)
                max_leg = max(fwd.residual, ret.residual)
                legs_ok = fwd.converged and ret.converged
                phase_ok = theta_res < theta_tol
                closed = legs_ok and phase_ok

                if closed:
                    notes = ""
                elif not fwd.converged:
                    notes = f"forward leg did not converge: residual={fwd.residual:.3e} km"
                elif not ret.converged:
                    notes = f"return leg did not converge: residual={ret.residual:.3e} km"
                else:
                    notes = f"theta not commensurate: residual={theta_res:.3e} rad"

                results.append(
                    CrossCycle(
                        connections=connections,
                        c_em=c_em,
                        c_se=c_se,
                        libration_pair=pair,
                        theta_closure_residual=theta_res,
                        closed=closed,
                        max_leg_residual=max_leg,
                        independent_residual=float("nan"),
                        notes=notes,
                    )
                )
    return results


# --- Task 5 (#411): time-consistent single-revolution cross-cycle corrector ------
#
# search_cross_cycle's theta_closure_residual is the FROZEN-orientation proxy
# |theta_ret - theta_fwd|; it ignores the relative phase advancing at omega_rel during
# the legs' transits and the orbit dwells. The true periodic-up-to-rotation cycle must
# satisfy TWO phase-consistency conditions (one per handoff), accounting for that time:
#
#   R1 = wrap[ theta_ret - theta_fwd - omega_rel*(t_fwd + n_se*T_se) ]   (fwd -> ret)
#   R2 = wrap[ theta_fwd - theta_ret - omega_rel*(t_ret + n_em*T_em) ]   (ret -> fwd)
#
# theta_fwd, theta_ret are OUTPUTS of each leg's spatial corrector (it solves theta to
# make the manifolds meet), so (R1, R2) are two scalar functions of the continuous
# amplitude knobs (c_em, c_se). Two equations, two unknowns -> isolated solutions. The
# #411 amplitude-knob analysis (docs/notes/2026-06-20-411-amplitude-theta-closure.md)
# shows each Delta-theta(C) sweeps the full circle within the manifold shadow budget, so
# a single-revolution (n_em=n_se=1) solution generically EXISTS. This corrector finds it
# by a damped finite-difference Newton over (c_em, c_se), holding the return-leg branch
# fixed (selected once at the seed) so the residual is a smooth function of the knobs.
#
# SUFFICIENCY remains the open scientific question: theta may close where the patch dV is
# high. The result carries total_patch_dv_kms so a high-dV closure is reported honestly,
# not hidden. Independent Radau verification is via crosscheck on the leg connections.


@dataclass(frozen=True)
class CrossCycleClosure:
    """Result of the #411 time-consistent single-rev cross-system cycle corrector.

    ``closed`` is True only when BOTH legs converged spatially AND the phase-consistency
    residual norm ``theta_residual_norm = hypot(R1, R2)`` is below ``theta_tol_rad``
    (periodic-up-to-rotation). ``total_patch_dv_kms`` is the sum of both legs' patch
    maneuvers -- the cost of the cycle, reported even when high so a theta-closed but
    expensive cycle is surfaced as a clean (non-)result, never hidden. ``cycle_time_s``
    is t_fwd + n_se*T_se + t_ret + n_em*T_em.
    """

    c_em: float
    c_se: float
    n_em: int
    n_se: int
    libration_pair: tuple[str, str]
    forward: CrossConnection
    ret: CrossConnection
    r1_rad: float
    r2_rad: float
    theta_residual_norm: float
    max_leg_residual_km: float
    total_patch_dv_kms: float
    cycle_time_s: float
    closed: bool
    n_iter: int
    notes: str = ""


_RETURN_VARIANTS: tuple[dict[str, int], ...] = (
    {"branch_u": +1, "branch_s": -1, "k": 1},
    {"branch_u": -1, "branch_s": +1, "k": 1},
    {"branch_u": +1, "branch_s": +1, "k": 1},
    {"branch_u": -1, "branch_s": -1, "k": 1},
)


def correct_cross_cycle(
    bridge: FrameBridge,
    *,
    em_lib: str,
    se_lib: str,
    c_em0: float,
    c_se0: float,
    n_em: int = 1,
    n_se: int = 1,
    theta_tol_rad: float = 1e-2,
    tol_km: float = 1e2,
    max_iter: int = 10,
    solver: str = "newton",
    c_em_bounds: tuple[float, float] | None = None,
    c_se_bounds: tuple[float, float] | None = None,
    fd_em: float = 2e-4,
    fd_se: float = 5e-7,
    max_attempts: int = 4,
    return_scan_n: int = 4,
    return_scan_n_tau: int = 2,
    return_max_time_factor: float = 4.0,
    on_iter: Callable[[int, float, float, float, float], None] | None = None,
    **conn_kwargs: object,
) -> CrossCycleClosure:
    """Solve the time-consistent single-rev cross-system cycle over (c_em, c_se).

    Damped finite-difference Newton on the two phase-consistency residuals (R1, R2)
    [see module Task-5 note], with the return-leg branch fixed after seed selection so
    the residual is smooth in the amplitude knobs. ``fd_se`` is ~400x smaller than
    ``fd_em``: the SE per-rev phase advance is hypersensitive to c_se (Delta-theta_se mod
    2pi sweeps the full circle over ~1e-4 in C, since omega_rel*T_se ~ 38 rad before the
    mod), so a larger step would wrap the phase and corrupt the finite-difference
    derivative. NEVER raises for non-closure -- returns
    ``closed=False`` with diagnostics; thresholds are NEVER loosened to manufacture
    closure, and ``total_patch_dv_kms`` is always reported.
    """
    em, se = bridge.em, bridge.se
    omega_rel = 1.0 / em.t_s - 1.0 / se.t_s
    variants = list(_RETURN_VARIANTS[: max(1, int(max_attempts))])

    def _solve(
        c_em: float, c_se: float, variant: dict[str, int] | None
    ) -> tuple[CrossConnection, CrossConnection, float, float, dict[str, int]] | None:
        """Build nodes, solve both legs; return (fwd, ret, T_em_s, T_se_s, variant)."""
        try:
            em_node = _build_em_node(em, em_lib, c_em)
            se_node = _build_se_node(se, se_lib, c_se)
        except (RuntimeError, ValueError):
            return None
        if not (em_node.converged and se_node.converged):
            return None
        fwd = correct_cross_connection(
            bridge,
            em_node,
            se_node,
            label_from=em_lib,
            label_to=se_lib,
            tol_km=tol_km,
            **conn_kwargs,  # type: ignore[arg-type]
        )
        chosen = variant
        ret_kwargs = dict(conn_kwargs)
        ret_kwargs.pop("scan_n", None)
        ret_kwargs.pop("scan_n_tau", None)
        ret_kwargs.pop("max_time_factor", None)
        if chosen is None:
            best: CrossConnection | None = None
            for v in variants:
                cand = correct_cross_connection(
                    bridge,
                    se_node,
                    em_node,
                    label_from=se_lib,
                    label_to=em_lib,
                    tol_km=tol_km,
                    branch_u=v["branch_u"],
                    branch_s=v["branch_s"],
                    k=v["k"],
                    scan_n=return_scan_n,
                    scan_n_tau=return_scan_n_tau,
                    max_time_factor=return_max_time_factor,
                    **ret_kwargs,  # type: ignore[arg-type]
                )
                if best is None or cand.residual < best.residual:
                    best, chosen = cand, v
                if cand.converged:
                    break
            assert best is not None and chosen is not None
            ret = best
        else:
            ret = correct_cross_connection(
                bridge,
                se_node,
                em_node,
                label_from=se_lib,
                label_to=em_lib,
                tol_km=tol_km,
                branch_u=chosen["branch_u"],
                branch_s=chosen["branch_s"],
                k=chosen["k"],
                scan_n=return_scan_n,
                scan_n_tau=return_scan_n_tau,
                max_time_factor=return_max_time_factor,
                **ret_kwargs,  # type: ignore[arg-type]
            )
        return fwd, ret, em_node.period * em.t_s, se_node.period * se.t_s, chosen

    def _resid(
        fwd: CrossConnection, ret: CrossConnection, t_em: float, t_se: float
    ) -> tuple[float, float]:
        r1 = _wrap_pi(ret.theta - fwd.theta - omega_rel * (fwd.transit_time + n_se * t_se))
        r2 = _wrap_pi(fwd.theta - ret.theta - omega_rel * (ret.transit_time + n_em * t_em))
        return r1, r2

    def _result(
        c_em: float,
        c_se: float,
        fwd: CrossConnection,
        ret: CrossConnection,
        t_em: float,
        t_se: float,
        n_iter: int,
        notes: str,
    ) -> CrossCycleClosure:
        legs_ok = fwd.converged and ret.converged
        if legs_ok:
            r1, r2 = _resid(fwd, ret, t_em, t_se)
            res_norm = math.hypot(r1, r2)
            cycle_time = fwd.transit_time + n_se * t_se + ret.transit_time + n_em * t_em
        else:
            r1 = r2 = res_norm = cycle_time = float("nan")
        closed = legs_ok and res_norm < theta_tol_rad
        return CrossCycleClosure(
            c_em=c_em,
            c_se=c_se,
            n_em=n_em,
            n_se=n_se,
            libration_pair=(em_lib, se_lib),
            forward=fwd,
            ret=ret,
            r1_rad=r1,
            r2_rad=r2,
            theta_residual_norm=res_norm,
            max_leg_residual_km=max(fwd.residual, ret.residual),
            total_patch_dv_kms=fwd.patch_dv_kms + ret.patch_dv_kms,
            cycle_time_s=cycle_time,
            closed=closed,
            n_iter=n_iter,
            notes=notes,
        )

    seed = _solve(c_em0, c_se0, None)
    if seed is None:
        raise ValueError(f"seed node construction failed at c_em={c_em0}, c_se={c_se0}")
    fwd, ret, t_em, t_se, variant = seed
    if not (fwd.converged and ret.converged) and solver != "feasibility_ls":
        # feasibility_ls relaxes the seed gate and drives leg closure jointly.
        return _result(c_em0, c_se0, fwd, ret, t_em, t_se, 0, "seed legs did not converge")

    if solver == "feasibility_ls":
        # #496 TRUE feasibility-first corrector (Braik-Ross 2025 blueprint, Phase 1 + 2
        # jointly): fold leg POSITION GAPS into the bounded least_squares residual so the
        # optimizer drives leg closure AND phase closure simultaneously over (c_em, c_se),
        # starting from an INFEASIBLE seed. The seed gate is relaxed — we do NOT require
        # both legs to converge before the solver runs.
        #
        # Residual (4-component):
        #   [fwd.residual_km / _LEG_SCALE,          <- forward leg spatial gap (O(1) at seed)
        #    ret.residual_km / _LEG_SCALE,           <- return  leg spatial gap (O(1) at seed)
        #    R1,                                     <- fwd->ret phase residual (rad)
        #    R2]                                     <- ret->fwd phase residual (rad)
        # R1/R2 are set to ±pi when a leg does not converge (max wrapped residual), so trf
        # steers toward leg closure first, then phase closure. c_se is clamped BELOW the
        # Canalias bifurcation (3.000863625) so the solver cannot wander off the SE-L2
        # family. FD steps are larger than bounded_ls (leg gaps are less sensitive to c_se
        # than phase, but the FD step must still be commensurate with the manifold response).
        _leg_scale = 3.0e5  # km: normalises ~2e5 km infeasible gap to O(1)
        lo_fs = [
            c_em_bounds[0] if c_em_bounds is not None else c_em0 - 0.05,
            c_se_bounds[0] if c_se_bounds is not None else max(3.0, c_se0 - 1.0e-3),
        ]
        hi_fs = [
            c_em_bounds[1] if c_em_bounds is not None else c_em0 + 0.05,
            c_se_bounds[1] if c_se_bounds is not None else c_se0 + 1.0e-3,
        ]
        if lo_fs[0] >= hi_fs[0] or lo_fs[1] >= hi_fs[1]:
            return _result(
                c_em0, c_se0, fwd, ret, t_em, t_se, 0, "feasibility_ls: degenerate bounds"
            )
        # Larger FD steps to capture leg-gap gradients (leg position is less sensitive to
        # c_se than phase, so fd_se_feas >> fd_se is safe here).
        _fd_em_feas = max(fd_em, 1.0e-3)
        _fd_se_feas = max(fd_se, 1.0e-4)
        nfev_fs = [0]

        def _feas_resid(x: NDArray[np.float64]) -> NDArray[np.float64]:
            nfev_fs[0] += 1
            c_em_v, c_se_v = float(x[0]), float(x[1])
            cand = _solve(c_em_v, c_se_v, variant)
            if cand is None:
                # Node construction failed — hard infeasible, steer back.
                return np.array([math.pi, math.pi, math.pi, math.pi], dtype=np.float64)
            fwd_c, ret_c, t_em_c, t_se_c, _ = cand
            # Leg position gaps (normalised to O(1) at the infeasible seed).
            f_gap = min(fwd_c.residual, _leg_scale) / _leg_scale
            r_gap = min(ret_c.residual, _leg_scale) / _leg_scale
            # Phase residuals: valid only when both legs converge spatially.
            if fwd_c.converged and ret_c.converged:
                rr1, rr2 = _resid(fwd_c, ret_c, t_em_c, t_se_c)
            else:
                rr1 = 0.0 if fwd_c.converged else math.pi
                rr2 = 0.0 if ret_c.converged else math.pi
            if on_iter is not None:
                on_iter(nfev_fs[0], c_em_v, c_se_v, rr1, rr2)
            return np.array([f_gap, r_gap, rr1, rr2], dtype=np.float64)

        sol_fs = least_squares(
            _feas_resid,
            x0=np.array([c_em0, c_se0], dtype=np.float64),
            bounds=(lo_fs, hi_fs),
            method="trf",
            diff_step=[_fd_em_feas / abs(c_em0), _fd_se_feas / abs(c_se0)],
            max_nfev=max(8, max_iter * 4),
            xtol=1e-12,
        )
        ce_fs, cs_fs = float(sol_fs.x[0]), float(sol_fs.x[1])
        final_fs = _solve(ce_fs, cs_fs, variant)
        if final_fs is None or not (final_fs[0].converged and final_fs[1].converged):
            # Could not find a leg-closing seed — clean negative.
            best_notes = (
                f"feasibility_ls: no leg-closing (c_em,c_se) found in bounds "
                f"c_em∈[{lo_fs[0]:.4f},{hi_fs[0]:.4f}] "
                f"c_se∈[{lo_fs[1]:.7f},{hi_fs[1]:.7f}]; reporting seed"
            )
            return _result(c_em0, c_se0, fwd, ret, t_em, t_se, nfev_fs[0], best_notes)
        return _result(
            ce_fs, cs_fs, final_fs[0], final_fs[1], final_fs[2], final_fs[3], nfev_fs[0], ""
        )

    if solver == "bounded_ls":
        # #496 stall fix (Braik-Ross 2025 / Ross-Scheeres blueprint: bound + feasibility-
        # first). The unbounded damped Newton below STALLS because its c_se step falls off
        # the SE Lyapunov family near the Canalias bifurcation -- no improving step keeps
        # both legs converged, so it breaks at ~0.59 rad (#411). A bounded trust-region
        # least_squares on (R1, R2) cannot leave the validated family band; infeasible
        # probes return the max wrapped residual so trf steers back. c_se band defaults
        # TIGHT (off-bifurcation); widen via c_se_bounds only within the stable family.
        lo = [
            c_em_bounds[0] if c_em_bounds is not None else c_em0 - 0.05,
            c_se_bounds[0] if c_se_bounds is not None else c_se0 - 1.0e-3,
        ]
        hi = [
            c_em_bounds[1] if c_em_bounds is not None else c_em0 + 0.05,
            c_se_bounds[1] if c_se_bounds is not None else c_se0 + 1.0e-3,
        ]
        if lo[0] >= hi[0] or lo[1] >= hi[1]:
            return _result(c_em0, c_se0, fwd, ret, t_em, t_se, 0, "bounded_ls: degenerate bounds")
        nfev = [0]

        def _ls_resid(x: NDArray[np.float64]) -> NDArray[np.float64]:
            nfev[0] += 1
            cand = _solve(float(x[0]), float(x[1]), variant)
            if cand is None or not (cand[0].converged and cand[1].converged):
                return np.array([math.pi, math.pi], dtype=np.float64)  # infeasible -> steer back
            rr1, rr2 = _resid(cand[0], cand[1], cand[2], cand[3])
            if on_iter is not None:
                on_iter(nfev[0], float(x[0]), float(x[1]), rr1, rr2)
            return np.array([rr1, rr2], dtype=np.float64)

        sol = least_squares(
            _ls_resid,
            x0=np.array([c_em0, c_se0], dtype=np.float64),
            bounds=(lo, hi),
            method="trf",
            diff_step=[fd_em / abs(c_em0), fd_se / abs(c_se0)],
            max_nfev=max(8, max_iter * 4),
            xtol=1e-12,
        )
        ce, cs = float(sol.x[0]), float(sol.x[1])
        final = _solve(ce, cs, variant)
        if final is None or not (final[0].converged and final[1].converged):
            return _result(
                c_em0,
                c_se0,
                fwd,
                ret,
                t_em,
                t_se,
                nfev[0],
                "bounded_ls: final eval off-family; reporting seed",
            )
        return _result(ce, cs, final[0], final[1], final[2], final[3], nfev[0], "")

    c = np.array([c_em0, c_se0], dtype=np.float64)
    fd = np.array([fd_em, fd_se], dtype=np.float64)
    r1, r2 = _resid(fwd, ret, t_em, t_se)
    f = np.array([r1, r2], dtype=np.float64)
    if on_iter is not None:
        on_iter(0, float(c[0]), float(c[1]), float(f[0]), float(f[1]))
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        if math.hypot(*f) < theta_tol_rad:
            break
        jac = np.zeros((2, 2), dtype=np.float64)
        ok = True
        for j in range(2):
            cp = c.copy()
            cp[j] += fd[j]
            probe = _solve(float(cp[0]), float(cp[1]), variant)
            if probe is None or not (probe[0].converged and probe[1].converged):
                ok = False
                break
            fp1, fp2 = _resid(probe[0], probe[1], probe[2], probe[3])
            jac[:, j] = (np.array([fp1, fp2]) - f) / fd[j]
        if not ok:
            break
        try:
            step = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            break
        alpha = 1.0
        improved = False
        for _ in range(20):
            ct = c + alpha * step
            cand = _solve(float(ct[0]), float(ct[1]), variant)
            if cand is not None and cand[0].converged and cand[1].converged:
                cf1, cf2 = _resid(cand[0], cand[1], cand[2], cand[3])
                if math.hypot(cf1, cf2) < math.hypot(*f):
                    c = ct
                    fwd, ret, t_em, t_se = cand[0], cand[1], cand[2], cand[3]
                    f = np.array([cf1, cf2])
                    improved = True
                    break
            alpha *= 0.5
        if on_iter is not None:
            on_iter(n_iter, float(c[0]), float(c[1]), float(f[0]), float(f[1]))
        if not improved:
            break

    res_norm = math.hypot(*f)
    notes = "" if res_norm < theta_tol_rad else f"theta residual {res_norm:.3e} rad > tol"
    return _result(float(c[0]), float(c[1]), fwd, ret, t_em, t_se, n_iter, notes)


def correct_cross_cycle_3d(
    bridge: FrameBridge,
    *,
    em_lib: str,
    se_lib: str,
    c_em0: float,
    c_se0: float,
    z_em: float,
    z_se: float,
    n_em: int = 1,
    n_se: int = 1,
    theta_tol_rad: float = 1e-2,
    tol_km: float = 1e2,
    max_iter: int = 10,
    solver: str = "bounded_ls",
    c_em_bounds: tuple[float, float] | None = None,
    c_se_bounds: tuple[float, float] | None = None,
    fd_em: float = 2e-4,
    fd_se: float = 5e-7,
    max_attempts: int = 4,
    return_scan_n: int = 4,
    return_scan_n_tau: int = 2,
    return_max_time_factor: float = 4.0,
    on_iter: Callable[[int, float, float, float, float], None] | None = None,
    **conn_kwargs: object,
) -> CrossCycleClosure:
    """Solve the time-consistent single-rev cross-system cycle in 3D.

    Solves over (c_em, c_se), holding z_em, z_se fixed. Damped finite-difference Newton
    or least-squares on the two phase-consistency residuals (R1, R2), building 3D periodic
    nodes (Halo/NRHO orbits) with the specified out-of-plane amplitudes.
    """
    em, se = bridge.em, bridge.se
    omega_rel = 1.0 / em.t_s - 1.0 / se.t_s
    variants = list(_RETURN_VARIANTS[: max(1, int(max_attempts))])

    def _solve(
        c_em: float, c_se: float, variant: dict[str, int] | None
    ) -> tuple[CrossConnection, CrossConnection, float, float, dict[str, int]] | None:
        """Build nodes, solve both legs; return (fwd, ret, T_em_s, T_se_s, variant)."""
        try:
            em_node = _build_em_node_3d(em, em_lib, c_em, z_em)
            se_node = _build_se_node_3d(se, se_lib, c_se, z_se)
        except (RuntimeError, ValueError):
            return None
        if not (em_node.converged and se_node.converged):
            return None
        fwd = correct_cross_connection(
            bridge,
            em_node,
            se_node,
            label_from=em_lib,
            label_to=se_lib,
            tol_km=tol_km,
            **conn_kwargs,  # type: ignore[arg-type]
        )
        chosen = variant
        ret_kwargs = dict(conn_kwargs)
        ret_kwargs.pop("scan_n", None)
        ret_kwargs.pop("scan_n_tau", None)
        ret_kwargs.pop("max_time_factor", None)
        if chosen is None:
            best: CrossConnection | None = None
            for v in variants:
                cand = correct_cross_connection(
                    bridge,
                    se_node,
                    em_node,
                    label_from=se_lib,
                    label_to=em_lib,
                    tol_km=tol_km,
                    branch_u=v["branch_u"],
                    branch_s=v["branch_s"],
                    k=v["k"],
                    scan_n=return_scan_n,
                    scan_n_tau=return_scan_n_tau,
                    max_time_factor=return_max_time_factor,
                    **ret_kwargs,  # type: ignore[arg-type]
                )
                if best is None or cand.residual < best.residual:
                    best, chosen = cand, v
                if cand.converged:
                    break
            assert best is not None and chosen is not None
            ret = best
        else:
            ret = correct_cross_connection(
                bridge,
                se_node,
                em_node,
                label_from=se_lib,
                label_to=em_lib,
                tol_km=tol_km,
                branch_u=chosen["branch_u"],
                branch_s=chosen["branch_s"],
                k=chosen["k"],
                scan_n=return_scan_n,
                scan_n_tau=return_scan_n_tau,
                max_time_factor=return_max_time_factor,
                **ret_kwargs,  # type: ignore[arg-type]
            )
        return fwd, ret, em_node.period * em.t_s, se_node.period * se.t_s, chosen

    def _resid(
        fwd: CrossConnection, ret: CrossConnection, t_em: float, t_se: float
    ) -> tuple[float, float]:
        r1 = _wrap_pi(ret.theta - fwd.theta - omega_rel * (fwd.transit_time + n_se * t_se))
        r2 = _wrap_pi(fwd.theta - ret.theta - omega_rel * (ret.transit_time + n_em * t_em))
        return r1, r2

    def _result(
        c_em: float,
        c_se: float,
        fwd: CrossConnection,
        ret: CrossConnection,
        t_em: float,
        t_se: float,
        n_iter: int,
        notes: str,
    ) -> CrossCycleClosure:
        legs_ok = fwd.converged and ret.converged
        if legs_ok:
            r1, r2 = _resid(fwd, ret, t_em, t_se)
            res_norm = math.hypot(r1, r2)
            cycle_time = fwd.transit_time + n_se * t_se + ret.transit_time + n_em * t_em
        else:
            r1 = r2 = res_norm = cycle_time = float("nan")
        closed = legs_ok and res_norm < theta_tol_rad
        return CrossCycleClosure(
            c_em=c_em,
            c_se=c_se,
            n_em=n_em,
            n_se=n_se,
            libration_pair=(em_lib, se_lib),
            forward=fwd,
            ret=ret,
            r1_rad=r1,
            r2_rad=r2,
            theta_residual_norm=res_norm,
            max_leg_residual_km=max(fwd.residual, ret.residual),
            total_patch_dv_kms=fwd.patch_dv_kms + ret.patch_dv_kms,
            cycle_time_s=cycle_time,
            closed=closed,
            n_iter=n_iter,
            notes=notes,
        )

    seed = _solve(c_em0, c_se0, None)
    if seed is None:
        raise ValueError(f"seed node construction failed at c_em={c_em0}, c_se={c_se0}")
    fwd, ret, t_em, t_se, variant = seed
    if not (fwd.converged and ret.converged) and solver != "feasibility_ls":
        return _result(c_em0, c_se0, fwd, ret, t_em, t_se, 0, "seed legs did not converge")

    if solver == "feasibility_ls":
        _leg_scale = 3.0e5
        lo_fs = [
            c_em_bounds[0] if c_em_bounds is not None else c_em0 - 0.05,
            c_se_bounds[0] if c_se_bounds is not None else max(3.0, c_se0 - 1.0e-3),
        ]
        hi_fs = [
            c_em_bounds[1] if c_em_bounds is not None else c_em0 + 0.05,
            c_se_bounds[1] if c_se_bounds is not None else c_se0 + 1.0e-3,
        ]
        if lo_fs[0] >= hi_fs[0] or lo_fs[1] >= hi_fs[1]:
            return _result(
                c_em0, c_se0, fwd, ret, t_em, t_se, 0, "feasibility_ls: degenerate bounds"
            )
        _fd_em_feas = max(fd_em, 1.0e-3)
        _fd_se_feas = max(fd_se, 1.0e-4)
        nfev_fs = [0]

        def _feas_resid(x: NDArray[np.float64]) -> NDArray[np.float64]:
            nfev_fs[0] += 1
            c_em_v, c_se_v = float(x[0]), float(x[1])
            cand = _solve(c_em_v, c_se_v, variant)
            if cand is None:
                return np.array([math.pi, math.pi, math.pi, math.pi], dtype=np.float64)
            fwd_c, ret_c, t_em_c, t_se_c, _ = cand
            f_gap = min(fwd_c.residual, _leg_scale) / _leg_scale
            r_gap = min(ret_c.residual, _leg_scale) / _leg_scale
            if fwd_c.converged and ret_c.converged:
                rr1, rr2 = _resid(fwd_c, ret_c, t_em_c, t_se_c)
            else:
                rr1 = 0.0 if fwd_c.converged else math.pi
                rr2 = 0.0 if ret_c.converged else math.pi
            if on_iter is not None:
                on_iter(nfev_fs[0], c_em_v, c_se_v, rr1, rr2)
            return np.array([f_gap, r_gap, rr1, rr2], dtype=np.float64)

        sol_fs = least_squares(
            _feas_resid,
            x0=np.array([c_em0, c_se0], dtype=np.float64),
            bounds=(lo_fs, hi_fs),
            method="trf",
            diff_step=[_fd_em_feas / abs(c_em0), _fd_se_feas / abs(c_se0)],
            max_nfev=max(8, max_iter * 4),
            xtol=1e-12,
        )
        ce_fs, cs_fs = float(sol_fs.x[0]), float(sol_fs.x[1])
        final_fs = _solve(ce_fs, cs_fs, variant)
        if final_fs is None or not (final_fs[0].converged and final_fs[1].converged):
            best_notes = (
                f"feasibility_ls: no leg-closing (c_em,c_se) found in bounds "
                f"c_em∈[{lo_fs[0]:.4f},{hi_fs[0]:.4f}] "
                f"c_se∈[{lo_fs[1]:.7f},{hi_fs[1]:.7f}]; reporting seed"
            )
            return _result(c_em0, c_se0, fwd, ret, t_em, t_se, nfev_fs[0], best_notes)
        return _result(
            ce_fs, cs_fs, final_fs[0], final_fs[1], final_fs[2], final_fs[3], nfev_fs[0], ""
        )

    if solver == "bounded_ls":
        lo = [
            c_em_bounds[0] if c_em_bounds is not None else c_em0 - 0.05,
            c_se_bounds[0] if c_se_bounds is not None else c_se0 - 1.0e-3,
        ]
        hi = [
            c_em_bounds[1] if c_em_bounds is not None else c_em0 + 0.05,
            c_se_bounds[1] if c_se_bounds is not None else c_se0 + 1.0e-3,
        ]
        if lo[0] >= hi[0] or lo[1] >= hi[1]:
            return _result(c_em0, c_se0, fwd, ret, t_em, t_se, 0, "bounded_ls: degenerate bounds")
        nfev = [0]

        def _ls_resid(x: NDArray[np.float64]) -> NDArray[np.float64]:
            nfev[0] += 1
            cand = _solve(float(x[0]), float(x[1]), variant)
            if cand is None or not (cand[0].converged and cand[1].converged):
                return np.array([math.pi, math.pi], dtype=np.float64)
            rr1, rr2 = _resid(cand[0], cand[1], cand[2], cand[3])
            if on_iter is not None:
                on_iter(nfev[0], float(x[0]), float(x[1]), rr1, rr2)
            return np.array([rr1, rr2], dtype=np.float64)

        sol = least_squares(
            _ls_resid,
            x0=np.array([c_em0, c_se0], dtype=np.float64),
            bounds=(lo, hi),
            method="trf",
            diff_step=[fd_em / abs(c_em0), fd_se / abs(c_se0)],
            max_nfev=max(8, max_iter * 4),
            xtol=1e-12,
        )
        ce, cs = float(sol.x[0]), float(sol.x[1])
        final = _solve(ce, cs, variant)
        if final is None or not (final[0].converged and final[1].converged):
            return _result(
                c_em0,
                c_se0,
                fwd,
                ret,
                t_em,
                t_se,
                nfev[0],
                "bounded_ls: final eval off-family; reporting seed",
            )
        return _result(ce, cs, final[0], final[1], final[2], final[3], nfev[0], "")

    c = np.array([c_em0, c_se0], dtype=np.float64)
    fd = np.array([fd_em, fd_se], dtype=np.float64)
    r1, r2 = _resid(fwd, ret, t_em, t_se)
    f = np.array([r1, r2], dtype=np.float64)
    if on_iter is not None:
        on_iter(0, float(c[0]), float(c[1]), float(f[0]), float(f[1]))
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        if math.hypot(*f) < theta_tol_rad:
            break
        jac = np.zeros((2, 2), dtype=np.float64)
        ok = True
        for j in range(2):
            cp = c.copy()
            cp[j] += fd[j]
            probe = _solve(float(cp[0]), float(cp[1]), variant)
            if probe is None or not (probe[0].converged and probe[1].converged):
                ok = False
                break
            fp1, fp2 = _resid(probe[0], probe[1], probe[2], probe[3])
            jac[:, j] = (np.array([fp1, fp2]) - f) / fd[j]
        if not ok:
            break
        try:
            step = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            break
        alpha = 1.0
        improved = False
        for _ in range(20):
            ct = c + alpha * step
            cand = _solve(float(ct[0]), float(ct[1]), variant)
            if cand is not None and cand[0].converged and cand[1].converged:
                cf1, cf2 = _resid(cand[0], cand[1], cand[2], cand[3])
                if math.hypot(cf1, cf2) < math.hypot(*f):
                    c = ct
                    fwd, ret, t_em, t_se = cand[0], cand[1], cand[2], cand[3]
                    f = np.array([cf1, cf2])
                    improved = True
                    break
            alpha *= 0.5
        if on_iter is not None:
            on_iter(n_iter, float(c[0]), float(c[1]), float(f[0]), float(f[1]))
        if not improved:
            break

    res_norm = math.hypot(*f)
    notes = "" if res_norm < theta_tol_rad else f"theta residual {res_norm:.3e} rad > tol"
    return _result(float(c[0]), float(c[1]), fwd, ret, t_em, t_se, n_iter, notes)
