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

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode, _seed_on_manifold


def se_earth_system() -> cr3bp.CR3BPSystem:
    """Build the Sun-Earth CR3BP system (Earth is the SE secondary at 1-mu_SE).

    The shared registry (``core.satellites``) carries Earth only as a *primary* and
    has no Sun entry, so ``cr3bp.cr3bp_system("Sun", "Earth")`` cannot construct the
    SE pair. We assemble it here from the same sourced constants used everywhere else:
    ``MU_SUN_KM3_S2`` (IAU 2015 nominal solar GM), Earth's planet-only GM, and Earth's
    heliocentric SMA. Convention matches ``cr3bp.cr3bp_system``: pair GM = G(m1+m2),
    ``l_km`` = secondary SMA about primary, ``t_s`` = sqrt(l_km^3 / G(m1+m2)).
    """
    earth = PLANETS["E"]
    gm_pair = MU_SUN_KM3_S2 + earth.mu_km3_s2
    mu = earth.mu_km3_s2 / gm_pair
    l_km = earth.sma_au * AU_KM
    t_s = math.sqrt(l_km**3 / gm_pair)
    return cr3bp.CR3BPSystem(mu=mu, primary="Sun", secondary="Earth", l_km=l_km, t_s=t_s)


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


def _manifold_inertial_at_section(
    bridge: FrameBridge,
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode,
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
    """
    if side not in {"em", "se"}:
        raise ValueError(f"side must be 'em' or 'se'; got {side!r}")
    seed = _seed_on_manifold(
        system, node, tau=tau, direction=direction, branch=branch, epsilon=epsilon
    )
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
    for k in range(sol.y.shape[1]):
        inert = _to_inertial(sol.y[:, k])
        g = float(inert[0] - x0_km)
        if prev_in is not None and prev_g != g and (prev_g <= 0.0 <= g or g <= 0.0 <= prev_g):
            frac = prev_g / (prev_g - g)  # linear interpolation of the crossing
            return (prev_in + frac * (inert - prev_in)).astype(np.float64)
        prev_in, prev_g = inert, g
    return None


def _cross_residual(
    bridge: FrameBridge,
    orbit_from: LyapunovNode,
    orbit_to: LyapunovNode,
    em: cr3bp.CR3BPSystem,
    se: cr3bp.CR3BPSystem,
    *,
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

    ``pos_gap`` and ``vel_gap`` are the inertial position/velocity differences of the
    EM unstable and SE stable manifolds at the section; ``patch_state_inertial`` is
    the EM-side crossing 6-state.
    """
    p_u = _manifold_inertial_at_section(
        bridge,
        em,
        orbit_from,
        side="em",
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
        se,
        orbit_to,
        side="se",
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
) -> tuple[float, float, float]:
    """Coarse 3-D grid over (tau_u, tau_s, theta) -> least inertial position-gap cell.

    theta is the most sensitive variable (it rotates one whole frame relative to the
    other) so it is gridded most densely. Falls back to the geometric centres if every
    cell misses the section.
    """
    thetas = [2.0 * math.pi * (i + 0.5) / n_theta for i in range(n_theta)]
    us = [period_u * (i + 0.5) / n_tau for i in range(n_tau)]
    ss = [period_s * (j + 0.5) / n_tau for j in range(n_tau)]
    best = float("inf")
    best_tu, best_ts, best_th = 0.5 * period_u, 0.5 * period_s, math.pi
    for th in thetas:
        for tu in us:
            for ts in ss:
                pos_gap, _, _ = resid(tu, ts, th)
                if pos_gap is None:
                    continue
                rn = float(np.linalg.norm(pos_gap))
                if rn < best:
                    best, best_tu, best_ts, best_th = rn, tu, ts, th
    return best_tu, best_ts, best_th


def correct_cross_connection(
    bridge: FrameBridge,
    orbit_from: LyapunovNode,
    orbit_to: LyapunovNode,
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
    em, se = bridge.em, bridge.se
    max_time_u = max_time_factor * orbit_from.period
    max_time_s = max_time_factor * orbit_to.period
    c_em = orbit_from.jacobi
    c_se = orbit_to.jacobi

    def _resid(
        tu: float, ts: float, th: float
    ) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None, NDArray[np.float64] | None]:
        return _cross_residual(
            bridge,
            orbit_from,
            orbit_to,
            em,
            se,
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
    )
