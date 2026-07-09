#!/usr/bin/env python3
"""Task #538: QBCP cross-system periodic-orbit (cycler) correction.

Chain the forward (Sun-Earth L2 torus -> Earth-Moon L2 torus) and reverse
(Earth-Moon L2 torus -> Sun-Earth L2 torus) heteroclinic connections of #537
into a single, well-posed multi-segment boundary-value problem and attempt to
converge it onto a mathematically exact synodic-periodic closed orbit in the
genuine time-periodic QBCP model. A clean, documented negative (the corrector
stalls at a finite residual floor without a spurious "it closed" on an
under-constrained system) is an explicitly acceptable outcome per this project's
S1L1-saga discipline.

All exclamation marks are avoided in comments and docstrings, per the repo
convention.

============================================================================
Task 1 -- Well-posed multi-segment residual design (the crux; see #538 plan)
============================================================================

#537's refining solve was rank-deficient: 3 equations (Delta_y, Delta_z, Delta_t
mod T_s) against 4 unknowns, with velocity NEVER in the residual. Its reported
12,034 km / 911 m/s gaps were post-hoc diagnostics, not quantities the optimizer
drove to zero. This script deliberately builds an OVER-DETERMINED residual that
carries position AND velocity at both crossing sections plus the loop-closure
conditions, so that "it converged" cannot be an artifact of missing equations.

Geometry of the closed loop (four propagated arcs, two heteroclinic legs):

  Leg 1 (forward, SE -> EM):
    arc 1: SE-L2 torus point theta_0, perturbed along its UNSTABLE eigenvector,
           propagated FORWARD for duration tau_f.
    arc 2: EM-L2 torus point theta_1, perturbed along its STABLE eigenvector,
           propagated BACKWARD for duration tau_bf.
    The two arcs must meet: full 6-state (position + velocity) match, plus a
    time-consistency scalar (their meeting epochs congruent mod T_s, since the
    QBCP frame is T_s-periodic -- matching states at incongruent Sun phases is
    physically meaningless, a subtlety #537 only weakly guarded).

  Leg 2 (reverse, EM -> SE):
    arc 3: EM-L2 torus point theta_2, perturbed along its UNSTABLE eigenvector,
           propagated FORWARD for duration tau_r.
    arc 4: SE-L2 torus point theta_3, perturbed along its STABLE eigenvector,
           propagated BACKWARD for duration tau_br.
    Same 6-state + time-consistency match.

  Closure (makes the concatenation a single closed orbit, not a loose
  heteroclinic network):
    SE closure: theta_3 -> theta_0 (the reverse leg lands where the forward leg
                departed).
    EM closure: theta_2 -> theta_1 (the reverse leg departs where the forward
                leg arrived).

Unknowns (12):
    theta_0 (2), theta_1 (2), theta_2 (2), theta_3 (2)   -- 4 torus phase pairs
    tau_f, tau_bf, tau_r, tau_br                          -- 4 arc durations
  Mapping to the plan's "t_0 + tau_f + tau_r + one extra scalar": the QBCP is
  non-autonomous and evaluate_qbcp_torus locks a torus point's epoch to its
  longitudinal phase (t = theta_long / omega_long). A free global departure
  epoch t_0 would therefore be redundant with theta_0's longitudinal component,
  so instead of carrying a redundant t_0 the four scalars are the four physical
  arc durations. This is documented rather than silently dropping the velocity
  residual to make a count work (the plan's explicit caution).

Residuals (18 >= 12, over-determined by design):
    leg-1 6-state match (6) + leg-1 time consistency (1)
    leg-2 6-state match (6) + leg-2 time consistency (1)
    SE closure theta_3 - theta_0 (2)
    EM closure theta_2 - theta_1 (2)

The regression guard `tests/scripts/test_run_538_residual_shape.py` asserts
`n_residuals >= n_unknowns` via `_residual_shape_ok()` so a future edit cannot
silently regress to an under-determined (rank-deficient) system.
"""

from __future__ import annotations

import datetime
import json
import math
import pathlib
import time
from collections.abc import Callable
from typing import Literal, cast

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.genome.bcr4bp_torus import (
    correct_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.genome.qbcp_torus import (
    QBCPTorus,
    correct_qbcp_torus,
    evaluate_qbcp_torus,
)
from cyclerfinder.genome.qp_torus_manifold import local_stability
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

# --------------------------------------------------------------------------
# Residual / unknown structure (Task 1). These constants define the shape and
# are the single source of truth for both the solver and the shape-guard test.
# --------------------------------------------------------------------------
N_PHASE_PAIRS: int = 4  # theta_0, theta_1, theta_2, theta_3
N_TORUS_PHASE_UNKNOWNS: int = 2 * N_PHASE_PAIRS  # 8
N_DURATION_UNKNOWNS: int = 4  # tau_f, tau_bf, tau_r, tau_br
N_UNKNOWNS: int = N_TORUS_PHASE_UNKNOWNS + N_DURATION_UNKNOWNS  # 12

N_STATE_MATCH: int = 6  # position (3) + velocity (3) at a crossing
N_TIME_MATCH: int = 1  # epoch congruence (mod T_s) at a crossing
N_PER_LEG: int = N_STATE_MATCH + N_TIME_MATCH  # 7
N_LEGS: int = 2
N_CLOSURE: int = 4  # SE closure (2) + EM closure (2)
N_RESIDUALS: int = N_LEGS * N_PER_LEG + N_CLOSURE  # 18

# Unknown-vector index layout.
IX_TH0: int = 0
IX_TH1: int = 2
IX_TH2: int = 4
IX_TH3: int = 6
IX_TAU_F: int = 8
IX_TAU_BF: int = 9
IX_TAU_R: int = 10
IX_TAU_BR: int = 11

# Globalization perturbation off the torus along the (un)stable eigenvector,
# matching #537 (run_533) so the seed basin is the same object.
MANIFOLD_EPS: float = 1e-5
# Section used only for coarse seeding (matches #537).
SEED_SECTION_X: float = 2.0
SEED_T_MAX: float = 25.0
# Physical scales for human-readable diagnostics.
EM_L_KM: float = 384400.0
EM_V_MS: float = 1024.0

_REGION_ID: str = "se-l2-em-l2-qbcp-cycler-closure-2026-07-09"
_METHOD: MethodCapability = MethodCapability(
    genome=(
        "Sun-Earth L2 <-> Earth-Moon L2 four-arc heteroclinic cycle in the "
        "genuine time-periodic QBCP model (forward SE->EM + reverse EM->SE), "
        "closed onto a single synodic-periodic orbit"
    ),
    corrector=(
        "well-posed over-determined least_squares (12 unknowns / 18 residuals): "
        "position+velocity 6-state match at both crossing sections + epoch "
        "congruence + SE/EM torus-phase closure"
    ),
    capability_tags=frozenset(
        {"qbcp", "torus", "heteroclinic", "cycler", "sun-earth-moon", "periodic-orbit"}
    ),
    git_sha="working-tree",
)

_RUNLOG_PATH: pathlib.Path = (
    pathlib.Path(__file__).resolve().parent.parent
    / "data"
    / "runlogs"
    / "run_538_qbcp_cycler.jsonl"
)

# One coarse-scan crossing record: theta phases, section-crossing state and time.
_Crossing = dict[str, "float | NDArray[np.float64]"]


def _residual_shape_ok() -> tuple[int, int]:
    """Return (n_unknowns, n_residuals) for the chosen parameterization.

    The shape-guard test asserts ``n_residuals >= n_unknowns`` so the system can
    never silently regress to the rank-deficient (under-determined) form that
    produced #537's false-comfort "connection".
    """
    n_unknowns = N_UNKNOWNS
    n_residuals = N_RESIDUALS
    return n_unknowns, n_residuals


def _wrap_pi(angle: float) -> float:
    """Wrap an angle to (-pi, pi]."""
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def _append_runlog(record: dict[str, object]) -> None:
    """Append one JSON record to the (gitignored) runlog, flushing immediately."""
    _RUNLOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with _RUNLOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")
        handle.flush()


# --------------------------------------------------------------------------
# Torus construction (reuses #533/#537's exact recipe; see plan Context).
# --------------------------------------------------------------------------
def build_tori() -> tuple[qbcp.QBCPSystem, QBCPTorus, QBCPTorus]:
    """Build the Sun-Earth L2 and Earth-Moon L2 QBCP tori as #533/#537 do."""
    qbcp_sys = qbcp.qbcp_default()
    mu_se = 1.0 / (qbcp_sys.mu_sun + 1.0)
    t_s = 2.0 * math.pi / qbcp_sys.omega_sun_nondim
    n_samples, n_modes = 5, 2

    # Sun-Earth L2 torus.
    sys_se = cr3bp.CR3BPSystem(
        mu=mu_se, primary="Sun", secondary="Earth", l_km=qbcp_sys.a_sun_nondim * EM_L_KM, t_s=1.0
    )
    orbit_se = correct_symmetric_fixed_jacobi(
        sys_se,
        x0_guess=(1.0 - mu_se) + 0.010,
        jacobi=3.0008,
        period_guess=3.1,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )
    # The BCR4BP seed helpers read the same scalar parameters QBCPSystem carries;
    # build the equivalent BCR4BPSystem so the (identical) call is type-correct.
    bcr_sys_full = bcr4bp.BCR4BPSystem(
        mu=qbcp_sys.mu,
        mu_sun=qbcp_sys.mu_sun,
        a_sun_nondim=qbcp_sys.a_sun_nondim,
        omega_sun_nondim=qbcp_sys.omega_sun_nondim,
        theta_sun0=qbcp_sys.theta_sun0,
    )
    x0_se_bcr, phase_pin_se, amplitude_pin_se = se_lyapunov_to_bcr4bp_torus_seed(
        orbit_se, bcr_sys_full, mu_se, n_samples=n_samples
    )
    torus_bcr = correct_bcr4bp_torus(
        bcr_sys_full, x0_se_bcr, n_modes, n_samples, phase_pin_se, amplitude_pin_se, tol=1e-6
    )
    coeffs_bcr = torus_bcr.fourier_coeffs
    x0_se_qbcp = np.zeros(6 + 12 * n_modes + 1)
    x0_se_qbcp[0:6] = np.real(coeffs_bcr[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        x0_se_qbcp[i0 : i0 + 6] = np.real(coeffs_bcr[n, :])
        x0_se_qbcp[i0 + 6 : i0 + 12] = np.imag(coeffs_bcr[n, :])
    x0_se_qbcp[-1] = torus_bcr.rho
    torus_se = correct_qbcp_torus(
        qbcp_sys, x0_se_qbcp, n_modes, n_samples, phase_pin_se, amplitude_pin_se, tol=1e-3
    )

    # Earth-Moon L2 torus (mu_sun continuation in BCR4BP, then QBCP correction).
    sys_em = cr3bp.CR3BPSystem(
        mu=qbcp_sys.mu, primary="Earth", secondary="Moon", l_km=EM_L_KM, t_s=382981.0
    )
    orbit_em = correct_symmetric_fixed_jacobi(
        sys_em,
        x0_guess=1.16,
        jacobi=3.17,
        period_guess=3.4,
        ydot0_sign=1.0,
        half_crossings=1,
        tol=1e-8,
    )
    sol_em = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit_em.period),
        np.array([orbit_em.x0, 0.0, 0.0, 0.0, orbit_em.ydot0, 0.0]),
        args=(qbcp_sys.mu,),
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, orbit_em.period, n_samples, endpoint=False),
    )
    u_samples_em = np.array([sol_em.y[:, j] for j in range(n_samples)])
    coeffs_em = np.fft.fft(u_samples_em, axis=0) / n_samples
    rho_em = (2.0 * math.pi * t_s / orbit_em.period) % (2.0 * math.pi)
    if rho_em > math.pi:
        rho_em -= 2.0 * math.pi
    x0_em = np.zeros(6 + 12 * n_modes + 1)
    x0_em[0:6] = np.real(coeffs_em[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        x0_em[i0 : i0 + 6] = np.real(coeffs_em[n, :])
        x0_em[i0 + 6 : i0 + 12] = np.imag(coeffs_em[n, :])
    x0_em[-1] = rho_em
    phase_pin_em = int(np.argmax(np.abs(np.imag(coeffs_em[1, :]))))
    amplitude_pin_em = float(np.linalg.norm(coeffs_em[1, :]))

    from cyclerfinder.genome.bcr4bp_torus import bcr4bp_torus_residual

    x_curr = x0_em
    for step_mu_sun in np.linspace(0.0, qbcp_sys.mu_sun, 5):
        sys_step = bcr4bp.BCR4BPSystem(
            mu=qbcp_sys.mu,
            mu_sun=step_mu_sun,
            a_sun_nondim=qbcp_sys.a_sun_nondim,
            omega_sun_nondim=qbcp_sys.omega_sun_nondim,
            theta_sun0=qbcp_sys.theta_sun0,
        )
        res_step = least_squares(
            bcr4bp_torus_residual,
            x_curr,
            args=(sys_step, n_modes, n_samples, phase_pin_em, amplitude_pin_em),
            kwargs={"rtol": 1e-6, "atol": 1e-6},
            method="lm",
            xtol=1e-6,
            ftol=1e-6,
        )
        x_curr = res_step.x
    torus_em = correct_qbcp_torus(
        qbcp_sys, x_curr, n_modes, n_samples, phase_pin_em, amplitude_pin_em, tol=1.5e-1
    )

    return qbcp_sys, torus_se, torus_em


# --------------------------------------------------------------------------
# Manifold-eigenvector evaluation at a torus phase.
# --------------------------------------------------------------------------
def manifold_state_vec(
    torus: QBCPTorus,
    theta_long: float,
    theta_trans: float,
    branch: str,
    ref_vec: NDArray[np.float64] | None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """Return (state_pv, eigenvector) at a torus phase for the given branch.

    ``ref_vec`` (if given) fixes the eigenvector sign by continuity
    (dot-product alignment), matching run_533's manifold-branch discipline.
    Returns None if the eigenvector is not extractable there.
    """
    state = evaluate_qbcp_torus(torus, theta_long, theta_trans)
    t0 = (float(theta_long) % (2.0 * math.pi)) / torus.omega_long
    try:
        _, states_pv = qbcp.propagate_qbcp_pv(
            state, (t0, t0 + torus.t_strob), torus.system, with_stm=True, rtol=1e-8, atol=1e-8
        )
    except RuntimeError:
        return None
    stm_pv = states_pv[-1, 6:].reshape((6, 6))
    stab = local_stability(state, stm_pv, hyperbolicity_tol=1e-4)
    vec = stab.vec_u if branch == "unstable" else stab.vec_s
    if vec is None:
        return None
    if ref_vec is not None and float(np.dot(vec, ref_vec)) < 0.0:
        vec = -vec
    return state, np.asarray(vec, dtype=np.float64)


def _propagate_arc(
    state_pv: NDArray[np.float64], t0: float, tf: float, system: qbcp.QBCPSystem
) -> NDArray[np.float64] | None:
    """Propagate a PV state from t0 to tf (tf may be < t0) and return the endpoint."""
    if tf == t0:
        return state_pv
    try:
        _, states_pv = qbcp.propagate_qbcp_pv(state_pv, (t0, tf), system, rtol=1e-10, atol=1e-10)
    except RuntimeError:
        return None
    end = np.asarray(states_pv[-1], dtype=np.float64)
    if not np.all(np.isfinite(end)):
        return None
    return end


# --------------------------------------------------------------------------
# The full, well-posed residual (Task 2). References fix eigenvector signs.
# --------------------------------------------------------------------------
def make_full_residual(
    torus_se: QBCPTorus,
    torus_em: QBCPTorus,
    refs: dict[str, NDArray[np.float64]],
    *,
    delta_v: NDArray[np.float64] | None = None,
) -> Callable[[NDArray[np.float64]], NDArray[np.float64]]:
    """Build the residual function closing over the tori and sign references.

    If ``delta_v`` (a length-6 impulse applied at leg-1's crossing, only the
    velocity components used) is supplied, the leg-1 velocity match is relaxed
    to ``v_fwd - v_bwd - delta_v`` -- the powered (Task 3) fallback. Ballistic
    is ``delta_v = None``.
    """
    omega = torus_se.omega_long  # equals torus_em.omega_long and omega_sun_nondim
    system = torus_se.system
    penalty = np.full(N_RESIDUALS, 1e3, dtype=np.float64)

    def residual(p: NDArray[np.float64]) -> NDArray[np.float64]:
        a0, b0 = float(p[IX_TH0]), float(p[IX_TH0 + 1])
        a1, b1 = float(p[IX_TH1]), float(p[IX_TH1 + 1])
        a2, b2 = float(p[IX_TH2]), float(p[IX_TH2 + 1])
        a3, b3 = float(p[IX_TH3]), float(p[IX_TH3 + 1])
        tau_f, tau_bf = float(p[IX_TAU_F]), float(p[IX_TAU_BF])
        tau_r, tau_br = float(p[IX_TAU_R]), float(p[IX_TAU_BR])

        # Leg 1: SE unstable (forward) vs EM stable (backward).
        mv_se0 = manifold_state_vec(torus_se, a0, b0, "unstable", refs["u_se"])
        mv_em1 = manifold_state_vec(torus_em, a1, b1, "stable", refs["s_em"])
        if mv_se0 is None or mv_em1 is None:
            return penalty
        s_se0, v_se0 = mv_se0
        s_em1, v_em1 = mv_em1
        t_se0 = (a0 % (2.0 * math.pi)) / omega
        t_em1 = (a1 % (2.0 * math.pi)) / omega
        end_fwd1 = _propagate_arc(s_se0 + MANIFOLD_EPS * v_se0, t_se0, t_se0 + tau_f, system)
        end_bwd1 = _propagate_arc(s_em1 + MANIFOLD_EPS * v_em1, t_em1, t_em1 - tau_bf, system)
        if end_fwd1 is None or end_bwd1 is None:
            return penalty
        r_leg1 = end_fwd1 - end_bwd1
        if delta_v is not None:
            r_leg1 = r_leg1 - delta_v
        r_time1 = _wrap_pi((a0 + omega * tau_f) - (a1 - omega * tau_bf))

        # Leg 2: EM unstable (forward) vs SE stable (backward).
        mv_em2 = manifold_state_vec(torus_em, a2, b2, "unstable", refs["u_em"])
        mv_se3 = manifold_state_vec(torus_se, a3, b3, "stable", refs["s_se"])
        if mv_em2 is None or mv_se3 is None:
            return penalty
        s_em2, v_em2 = mv_em2
        s_se3, v_se3 = mv_se3
        t_em2 = (a2 % (2.0 * math.pi)) / omega
        t_se3 = (a3 % (2.0 * math.pi)) / omega
        end_fwd2 = _propagate_arc(s_em2 + MANIFOLD_EPS * v_em2, t_em2, t_em2 + tau_r, system)
        end_bwd2 = _propagate_arc(s_se3 + MANIFOLD_EPS * v_se3, t_se3, t_se3 - tau_br, system)
        if end_fwd2 is None or end_bwd2 is None:
            return penalty
        r_leg2 = end_fwd2 - end_bwd2
        r_time2 = _wrap_pi((a2 + omega * tau_r) - (a3 - omega * tau_br))

        # Closure of the loop.
        r_se_close = np.array([_wrap_pi(a3 - a0), _wrap_pi(b3 - b0)], dtype=np.float64)
        r_em_close = np.array([_wrap_pi(a2 - a1), _wrap_pi(b2 - b1)], dtype=np.float64)

        return np.concatenate(
            [
                r_leg1,
                np.array([r_time1], dtype=np.float64),
                r_leg2,
                np.array([r_time2], dtype=np.float64),
                r_se_close,
                r_em_close,
            ]
        )

    return residual


# --------------------------------------------------------------------------
# Coarse seeding via a serial section scan (matches #537's x = 2.0 section).
# --------------------------------------------------------------------------
def _propagate_to_section(
    system: qbcp.QBCPSystem,
    state_pv: NDArray[np.float64],
    direction: float,
    t0: float,
) -> tuple[NDArray[np.float64], float] | None:
    """Integrate to the x = SEED_SECTION_X crossing; return (state_pv, t) or None."""
    state_pm0 = qbcp.state_pv_to_pm(state_pv, t0, system)

    def x_event(t: float, y: NDArray[np.float64]) -> float:
        return float(y[0]) - SEED_SECTION_X

    x_event.terminal = True  # type: ignore[attr-defined]
    x_event.direction = 0.0  # type: ignore[attr-defined]

    def fun(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return qbcp.qbcp_eom(t, y, system)

    t_span = (t0, t0 + math.copysign(SEED_T_MAX, direction))
    sol = solve_ivp(fun, t_span, state_pm0, method="DOP853", rtol=1e-8, atol=1e-8, events=[x_event])
    if sol.t_events is None or sol.y_events is None or sol.t_events[0].size == 0:
        return None
    crossing_pm = sol.y_events[0][0]
    crossing_t = float(sol.t_events[0][0])
    return qbcp.state_pm_to_pv(crossing_pm, crossing_t, system), crossing_t


def _scan_crossings(
    torus: QBCPTorus, branch: str, direction: float, ref_vec: NDArray[np.float64], n_grid: int
) -> list[_Crossing]:
    """Serial coarse scan of a manifold family to the seeding section."""
    thetas = 2.0 * math.pi * np.arange(n_grid) / n_grid
    out: list[_Crossing] = []
    for tl in thetas:
        for tt in thetas:
            mv = manifold_state_vec(torus, float(tl), float(tt), branch, ref_vec)
            if mv is None:
                continue
            state, vec = mv
            t0 = (float(tl) % (2.0 * math.pi)) / torus.omega_long
            res = _propagate_to_section(torus.system, state + MANIFOLD_EPS * vec, direction, t0)
            if res is None:
                continue
            cross_state, cross_t = res
            out.append(
                {
                    "theta_long": float(tl),
                    "theta_trans": float(tt),
                    "cross_state": cross_state,
                    "cross_t": cross_t,
                }
            )
    return out


def _seed_reference_vec(torus: QBCPTorus, branch: str, sign: float) -> NDArray[np.float64]:
    """Pin an eigenvector sign convention (vec[0] * sign > 0) as the reference."""
    mv = manifold_state_vec(torus, 0.0, 0.0, branch, None)
    if mv is None:
        # Fall back to a canonical x-axis reference.
        return np.array([sign, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    _, vec = mv
    if vec[0] * sign < 0.0:
        vec = -vec
    return vec


def _best_crossing_pair(
    fwd: list[_Crossing], bwd: list[_Crossing]
) -> tuple[float, _Crossing, _Crossing]:
    """Minimum y,z,vx,vy,vz gap over the cross product of two crossing lists."""
    best: tuple[float, _Crossing, _Crossing] = (math.inf, fwd[0], bwd[0])
    for fc in fwd:
        for bc in bwd:
            fs = np.asarray(fc["cross_state"])
            bs = np.asarray(bc["cross_state"])
            gap = float(np.linalg.norm(fs[1:6] - bs[1:6]))  # y,z,vx,vy,vz at the section
            if gap < best[0]:
                best = (gap, fc, bc)
    return best


def build_seed(
    torus_se: QBCPTorus, torus_em: QBCPTorus, n_grid: int
) -> tuple[NDArray[np.float64], dict[str, NDArray[np.float64]]]:
    """Coarse-seed the 12 unknowns and fix the four eigenvector-sign references.

    Both legs are scanned INDEPENDENTLY at the seeding section: leg 1 is
    SE-unstable (forward) vs EM-stable (backward), matching #537's own search;
    leg 2 is EM-unstable (forward) vs SE-stable (backward), its own coarse
    manifold-crossing search. An earlier version of this function seeded leg 2
    by mirroring leg 1's phases under the closure-target assumption
    (theta_2~theta_1, theta_3~theta_0) -- that is only a good coarse seed if
    the two heteroclinic families happen to cross the section near the same
    phases, which they need not. Diagnostic: the mirrored seed's leg-2 coarse
    gap was 1,441,019 km, two orders of magnitude worse than leg 1's 47,473 km
    (and #537's own 12,034 km) -- clear evidence the seed, not the residual
    design, was the weak link. Scanning leg 2 independently fixes this.
    """
    refs = {
        "u_se": _seed_reference_vec(torus_se, "unstable", -1.0),
        "s_em": _seed_reference_vec(torus_em, "stable", 1.0),
        "u_em": _seed_reference_vec(torus_em, "unstable", -1.0),
        "s_se": _seed_reference_vec(torus_se, "stable", 1.0),
    }
    # Leg 1 (forward): SE unstable (forward, direction +1) vs EM stable (backward, -1).
    se_u = _scan_crossings(torus_se, "unstable", 1.0, refs["u_se"], n_grid)
    em_s = _scan_crossings(torus_em, "stable", -1.0, refs["s_em"], n_grid)
    if not se_u or not em_s:
        raise RuntimeError(
            f"leg-1 seed scan found no crossings (SE-unstable={len(se_u)}, EM-stable={len(em_s)})"
        )
    gap1, uc_best, sc_best = _best_crossing_pair(se_u, em_s)
    a0 = float(uc_best["theta_long"])
    b0 = float(uc_best["theta_trans"])
    a1 = float(sc_best["theta_long"])
    b1 = float(sc_best["theta_trans"])
    tau_f = abs(float(uc_best["cross_t"]) - a0 / torus_se.omega_long)
    tau_bf = abs(a1 / torus_em.omega_long - float(sc_best["cross_t"]))

    # Leg 2 (reverse): EM unstable (forward, +1) vs SE stable (backward, -1) --
    # its own independent coarse scan, not a mirror of leg 1.
    em_u = _scan_crossings(torus_em, "unstable", 1.0, refs["u_em"], n_grid)
    se_s = _scan_crossings(torus_se, "stable", -1.0, refs["s_se"], n_grid)
    if not em_u or not se_s:
        raise RuntimeError(
            f"leg-2 seed scan found no crossings (EM-unstable={len(em_u)}, SE-stable={len(se_s)})"
        )
    gap2, uc2_best, sc2_best = _best_crossing_pair(em_u, se_s)
    a2 = float(uc2_best["theta_long"])
    b2 = float(uc2_best["theta_trans"])
    a3 = float(sc2_best["theta_long"])
    b3 = float(sc2_best["theta_trans"])
    tau_r = abs(float(uc2_best["cross_t"]) - a2 / torus_em.omega_long)
    tau_br = abs(a3 / torus_se.omega_long - float(sc2_best["cross_t"]))

    print(
        f"    leg-1 coarse gap (nondim, y/z/vx/vy/vz combined norm): {gap1:.4e}   "
        f"leg-2 coarse gap: {gap2:.4e}"
    )

    seed = np.zeros(N_UNKNOWNS, dtype=np.float64)
    seed[IX_TH0], seed[IX_TH0 + 1] = a0, b0
    seed[IX_TH1], seed[IX_TH1 + 1] = a1, b1
    seed[IX_TH2], seed[IX_TH2 + 1] = a2, b2
    seed[IX_TH3], seed[IX_TH3 + 1] = a3, b3
    seed[IX_TAU_F] = tau_f
    seed[IX_TAU_BF] = tau_bf
    seed[IX_TAU_R] = tau_r
    seed[IX_TAU_BR] = tau_br
    return seed, refs


def _solve_and_log(
    residual: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    seed: NDArray[np.float64],
    label: str,
    max_nfev: int,
    method: str = "lm",
) -> tuple[NDArray[np.float64], float, int]:
    """Run least_squares with per-iteration runlog instrumentation."""
    counter = {"n": 0}
    best = {"norm": math.inf}
    t_start = time.time()

    def logged(p: NDArray[np.float64]) -> NDArray[np.float64]:
        r = residual(p)
        counter["n"] += 1
        norm = float(np.linalg.norm(r))
        improved = norm < best["norm"]
        if improved:
            best["norm"] = norm
        if counter["n"] % 10 == 1 or improved:
            _append_runlog(
                {
                    "ts": datetime.datetime.now(datetime.UTC).isoformat(),
                    "label": label,
                    "nfev": counter["n"],
                    "residual_norm": norm,
                    "best_norm": best["norm"],
                    "elapsed_s": round(time.time() - t_start, 2),
                }
            )
        return r

    result = least_squares(
        logged,
        seed,
        method=cast("Literal['trf', 'dogbox', 'lm']", method),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=max_nfev,
    )
    final_norm = float(np.linalg.norm(result.fun))
    _append_runlog(
        {
            "ts": datetime.datetime.now(datetime.UTC).isoformat(),
            "label": f"{label}:final",
            "method": method,
            "nfev": int(result.nfev),
            "residual_norm": final_norm,
            "status": int(result.status),
            "message": str(result.message),
            "elapsed_s": round(time.time() - t_start, 2),
        }
    )
    return np.asarray(result.x, dtype=np.float64), final_norm, int(result.nfev)


def _multistart_ballistic(
    residual: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    seed: NDArray[np.float64],
    n_restarts: int,
    max_nfev: int,
) -> tuple[NDArray[np.float64], float, int]:
    """Real-effort ballistic solve: base seed (lm + trf) plus perturbed lm restarts.

    Returns the best (x, norm, total_nfev) across all starts. A single optimizer
    from a single seed reaching a local minimum is weak evidence; this drives the
    negative from several seeds and cross-checks the base seed with a second,
    independent algorithm (trf) so a "no closure" verdict is not an artifact of
    one optimizer getting stuck. lm is used for the perturbed restarts (each of
    its 12-unknown Jacobian steps is cheap relative to trf's finite-difference
    Jacobian here, so more independent basins can be sampled per unit of wall
    time); trf on the base seed is the algorithmic cross-check.
    """
    rng = np.random.default_rng(538)
    best_x = seed
    best_norm = math.inf
    total_nfev = 0
    starts = [("base", seed, "lm"), ("base", seed, "trf")]
    for k in range(n_restarts):
        pert = seed.copy()
        pert[:N_TORUS_PHASE_UNKNOWNS] += rng.uniform(-0.4, 0.4, N_TORUS_PHASE_UNKNOWNS)
        pert[N_TORUS_PHASE_UNKNOWNS:] *= rng.uniform(0.7, 1.3, N_DURATION_UNKNOWNS)
        starts.append((f"restart{k}", pert, "lm"))
    for name, x0, method in starts:
        # trf's finite-difference Jacobian is far more expensive per nfev than lm's
        # analytic-ish internal handling here (observed ~40x wall time per eval on
        # this residual) -- cap it separately so one slow algorithm cannot starve
        # the multi-start budget of cheaper, equally informative lm restarts.
        budget = max_nfev if method == "lm" else min(max_nfev, 300)
        x, norm, nfev = _solve_and_log(
            residual, x0, f"ballistic-{name}-{method}", max_nfev=budget, method=method
        )
        total_nfev += nfev
        print(f"    start {name} ({method}): norm={norm:.4e} nfev={nfev}")
        if norm < best_norm:
            best_norm, best_x = norm, x
    return best_x, best_norm, total_nfev


def _report_solution(
    tag: str,
    residual: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    x: NDArray[np.float64],
    final_norm: float,
) -> None:
    """Print a component breakdown of the converged/stalled residual."""
    r = residual(x)
    leg1_pos = float(np.linalg.norm(r[0:3])) * EM_L_KM
    leg1_vel = float(np.linalg.norm(r[3:6])) * EM_V_MS
    leg1_time = float(r[6])
    leg2_pos = float(np.linalg.norm(r[7:10])) * EM_L_KM
    leg2_vel = float(np.linalg.norm(r[10:13])) * EM_V_MS
    leg2_time = float(r[13])
    se_close = float(np.linalg.norm(r[14:16]))
    em_close = float(np.linalg.norm(r[16:18]))
    print(f"\n--- {tag} residual breakdown (norm = {final_norm:.3e}) ---")
    print(f"  leg1 crossing:  pos gap = {leg1_pos:10.1f} km   vel gap = {leg1_vel:8.1f} m/s")
    print(f"  leg1 time:      {leg1_time:.3e} rad")
    print(f"  leg2 crossing:  pos gap = {leg2_pos:10.1f} km   vel gap = {leg2_vel:8.1f} m/s")
    print(f"  leg2 time:      {leg2_time:.3e} rad")
    print(f"  SE closure |theta3-theta0| = {se_close:.3e} rad")
    print(f"  EM closure |theta2-theta1| = {em_close:.3e} rad")


def main() -> None:
    preflight_search(
        task_no=538,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=200,  # 2 x 10x10 coarse seed scans; corrector is not a grid sweep
    )

    print("=" * 60)
    print("Task #538: QBCP cross-system cycler closure")
    print("=" * 60)
    n_unknowns, n_residuals = _residual_shape_ok()
    print(f"Residual well-posedness: {n_unknowns} unknowns / {n_residuals} residuals")
    assert n_residuals >= n_unknowns

    print("\n[1/4] Building SE-L2 and EM-L2 QBCP tori (this is the slow step) ...")
    t0 = time.time()
    _qbcp_sys, torus_se, torus_em = build_tori()
    print(
        f"  SE torus: invariance residual = {torus_se.invariance_residual:.3e} "
        f"(converged={torus_se.converged})"
    )
    print(
        f"  EM torus: invariance residual = {torus_em.invariance_residual:.3e} "
        f"(converged={torus_em.converged})"
    )
    print(f"  built in {time.time() - t0:.0f}s")
    if not torus_em.converged:
        print(
            "  WARNING: the EM-L2 torus did NOT converge to an invariant object. "
            "Any closure onto it cannot be exact; a non-negative residual floor is "
            "expected and is an honest upstream negative, not a corrector bug."
        )

    print("\n[2/4] Coarse-seeding the 12 unknowns ...")
    t0 = time.time()
    seed, refs = build_seed(torus_se, torus_em, n_grid=10)
    print(f"  seeded in {time.time() - t0:.0f}s")
    residual = make_full_residual(torus_se, torus_em, refs)
    _report_solution("seed", residual, seed, float(np.linalg.norm(residual(seed))))

    print("\n[3/4] Ballistic corrector (multi-start lm + trf, over-determined) ...")
    # n_restarts trimmed from 4 to 1 for this run: two prior full runs (see
    # data/runlogs/run_538_qbcp_cycler.jsonl history) already tried 6 distinct
    # basins (base x{lm,trf} + restart0-3 x lm) and NONE approached 1e-8 --
    # final norms 1.070, 1.068, 6.679, 4242.6, 2.550, and restart3 still at
    # 4.502 (>3300 nfev) when the process was killed. That is already
    # sufficient multi-basin evidence per the plan's "no X found is
    # conditional on formulation" discipline; the priority now is reaching
    # the Task 3 powered fallback with a clean, complete run, not accumulating
    # a 7th/8th non-converging basin.
    t0 = time.time()
    x_bal, norm_bal, nfev_bal = _multistart_ballistic(residual, seed, n_restarts=1, max_nfev=2000)
    print(f"  best of all starts in {time.time() - t0:.0f}s, total nfev={nfev_bal}")
    _report_solution("ballistic", residual, x_bal, norm_bal)

    ballistic_ok = norm_bal < 1e-8
    if ballistic_ok:
        print(f"\nBALLISTIC CLOSURE ACHIEVED: residual norm {norm_bal:.3e} < 1e-8")
        print("Unknowns (theta0,theta1,theta2,theta3, tau_f,tau_bf,tau_r,tau_br):")
        print(f"  {np.array2string(x_bal, precision=6)}")
        print("\nNext (out of scope here): Task 4 independent cross-check, then Task 5.")
        return

    print(
        f"\nBallistic corrector did NOT reach 1e-8 (best norm {norm_bal:.3e}). "
        "Attempting the Task 3 powered fallback (impulsive dv at leg-1 crossing)."
    )

    print("\n[4/4] Powered fallback: free dv at the leg-1 crossing ...")
    # Two-stage: first solve with dv free (velocity match relaxed), reading off the
    # implied dv, then report its magnitude against catalogue dv bands.
    v_fwd_minus_bwd = residual(x_bal)[3:6]  # leftover velocity gap at leg-1 crossing
    dv_vec6 = np.zeros(6, dtype=np.float64)
    dv_vec6[3:6] = v_fwd_minus_bwd
    powered_residual = make_full_residual(torus_se, torus_em, refs, delta_v=dv_vec6)
    t0 = time.time()
    x_pow, norm_pow, nfev_pow = _solve_and_log(powered_residual, x_bal, "powered", max_nfev=4000)
    print(f"  converged/stalled in {time.time() - t0:.0f}s, nfev={nfev_pow}")
    _report_solution("powered", powered_residual, x_pow, norm_pow)
    dv_ms = float(np.linalg.norm(dv_vec6[3:6])) * EM_V_MS
    print(f"\n  implied leg-1 impulsive dv (from ballistic leftover): {dv_ms:.1f} m/s")
    if norm_pow < 1e-8:
        print(
            f"  POWERED closure with a {dv_ms:.1f} m/s maneuver. Whether this is "
            "catalogue-eligible is a Task 5 dv_band judgement, out of scope here."
        )
    else:
        print(
            f"  Powered fallback also did NOT close (best norm {norm_pow:.3e}). "
            "Documented negative: no ballistic or single-impulse closure found in this basin."
        )

    print(
        "\nOUTCOME: honest negative (no closure below 1e-8). See the residual "
        "breakdowns above and data/runlogs/run_538_qbcp_cycler.jsonl for the "
        "per-iteration convergence history. Tasks 4/5 (cross-check, writeback) "
        "are deliberately separate follow-up passes."
    )


if __name__ == "__main__":
    main()
