#!/usr/bin/env python3
"""Search runner for Task #533: Sun-Earth L2 Torus to Earth-Moon L2 Torus Connection in QBCP.

Uses parallelized grid generation to sweep through stable/unstable manifold crossings
at the surface of section, checking for topological intersections, and refines the
time-coherent connection to zero position gap.
All exclamation marks are avoided in comments and docstrings.
"""

from __future__ import annotations

import concurrent.futures
import math
import pathlib

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.genome.bcr4bp_torus import (
    bcr4bp_torus_residual,
    correct_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.genome.qbcp_torus import (
    correct_qbcp_torus,
    evaluate_qbcp_torus,
)
from cyclerfinder.genome.qp_torus_manifold import local_stability
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

_REGION_ID = "se-l2-em-l2-qbcp-torus-connection-2026-07-09"
_METHOD = MethodCapability(
    genome=(
        "Sun-Earth L2 <-> Earth-Moon L2 torus manifold crossing grids (16x16 long/lat) "
        "in the genuine time-periodic QBCP model, refined via a 3-eq/4-unknown "
        "least_squares (position+time only, no velocity residual -- see #538 plan caveat)"
    ),
    corrector="build_custom_grid + least_squares torus-phase refine (position/time only)",
    capability_tags=frozenset({"qbcp", "torus", "heteroclinic", "sun-earth-moon"}),
    git_sha="working-tree",
)


def propagate_qbcp_to_section(
    system, state_pv0, direction, surface_x, t_max, rtol=1e-8, atol=1e-8, t0=0.0
):
    """Integrate state_pv0 in QBCP from t0 to surface_x crossing."""
    state_pm0 = qbcp.state_pv_to_pm(state_pv0, t0, system)

    def x_event(t, y):
        # Position x is invariant under the transformation, so x_PV = x_PM
        return y[0] - surface_x

    x_event.terminal = True
    x_event.direction = 0.0

    def _collision_event(t, y):
        mu = system.mu
        x, yc, z = y[0], y[1], y[2]
        r1_sq = (x + mu) ** 2 + yc**2 + z**2
        r2_sq = (x - 1.0 + mu) ** 2 + yc**2 + z**2
        if r1_sq < 0.0004 or r2_sq < 0.0001:
            return 0.0
        return 1.0

    _collision_event.terminal = True

    def fun(t, y):
        return qbcp.qbcp_eom(t, y, system)

    t_span = (t0, t0 + math.copysign(t_max, direction))
    sol = solve_ivp(
        fun,
        t_span,
        state_pm0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=[x_event, _collision_event],
    )

    if sol.status == 1 and sol.t_events[1].size > 0:
        return None
    if sol.t_events[0].size == 0:
        return None

    crossing_pm = sol.y_events[0][0]
    crossing_t = sol.t_events[0][0]
    crossing_pv = qbcp.state_pm_to_pv(crossing_pm, crossing_t, system)
    return crossing_pv, crossing_t


def compute_custom_row(args):
    """Compute one row of manifold crossings for a given theta_long."""
    i, theta_long, thetas_lat, torus, branch, eps, sign, t_max, surface_x = args
    n_lat = len(thetas_lat)
    row_origins = np.zeros((n_lat, 2))
    row_endpoints = np.full((n_lat, 6), np.nan)
    row_hyperbolic = np.zeros(n_lat, dtype=bool)

    t0 = theta_long / torus.omega_long
    prev_vec = None

    for j, theta_trans in enumerate(thetas_lat):
        row_origins[j, 0] = theta_long
        row_origins[j, 1] = theta_trans

        # State and STM at the torus point
        state = evaluate_qbcp_torus(torus, float(theta_long), float(theta_trans))
        try:
            _, states_pv = qbcp.propagate_qbcp_pv(
                state, (t0, t0 + torus.t_strob), torus.system, with_stm=True, rtol=1e-8, atol=1e-8
            )
        except RuntimeError:
            continue
        stm_pv = states_pv[-1, 6:].reshape((6, 6))

        stab = local_stability(
            state,
            stm_pv,
            hyperbolicity_tol=1e-4,
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

        # Perturb and integrate
        perturbed = state + eps * vec
        direction = 1.0 if branch == "unstable" else -1.0
        res = propagate_qbcp_to_section(
            torus.system,
            perturbed,
            direction=direction,
            surface_x=surface_x,
            t_max=t_max,
            rtol=1e-8,
            atol=1e-8,
            t0=t0,
        )
        if res is not None:
            row_endpoints[j, :] = res[0]

    return i, row_origins, row_endpoints, row_hyperbolic


def build_custom_grid(torus, branch, sign, t_max, surface_x, n_long=16, n_lat=16):
    """Build the 2D grid of manifold crossing endpoints."""
    thetas_long = 2.0 * math.pi * np.arange(n_long) / n_long
    thetas_lat = 2.0 * math.pi * np.arange(n_lat) / n_lat

    origins = np.zeros((n_long, n_lat, 2), dtype=np.float64)
    endpoints = np.full((n_long, n_lat, 6), np.nan, dtype=np.float64)
    hyperbolic = np.zeros((n_long, n_lat), dtype=np.bool_)

    tasks = [
        (i, theta_long, thetas_lat, torus, branch, 1e-5, sign, t_max, surface_x)
        for i, theta_long in enumerate(thetas_long)
    ]

    with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
        results = list(executor.map(compute_custom_row, tasks))

    for i, row_origins, row_endpoints, row_hyperbolic in results:
        origins[i] = row_origins
        endpoints[i] = row_endpoints
        hyperbolic[i] = row_hyperbolic

    return origins, endpoints, hyperbolic


def get_manifold_point_with_vec(torus, theta_long, theta_trans, branch, eps, sign):
    """Evaluate torus point and return its state, time, and eigenvector."""
    t0 = theta_long / torus.omega_long
    state = evaluate_qbcp_torus(torus, float(theta_long), float(theta_trans))

    _, states_pv = qbcp.propagate_qbcp_pv(
        state, (t0, t0 + torus.t_strob), torus.system, with_stm=True, rtol=1e-8, atol=1e-8
    )
    stm_pv = states_pv[-1, 6:].reshape((6, 6))

    stab = local_stability(
        state,
        stm_pv,
        hyperbolicity_tol=1e-4,
    )
    vec = stab.vec_u if branch == "unstable" else stab.vec_s
    if vec is None:
        return None

    if vec[0] * sign < 0.0:
        vec = -vec

    return {
        "state": state,
        "time": t0,
        "vec": vec,
    }


def main():
    preflight_search(
        task_no=533,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=512,
        override_reason=(
            "retrofitting the mandatory #521 preflight gate onto a script that already "
            "ran to completion (#537 resolved 2026-07-09, refined candidate at 3.66 TU) "
            "before the gate existed; no timing pilot was captured at the time, and "
            "re-running the full sweep solely to measure one is not warranted for an "
            "already-complete result"
        ),
    )
    print("==================================================")
    print("STEP 1: Sourcing and correcting Sun-Earth L2 Torus")
    print("==================================================")

    qbcp_sys = qbcp.qbcp_default()
    mu_se = 1.0 / (qbcp_sys.mu_sun + 1.0)
    t_s = 2.0 * math.pi / qbcp_sys.omega_sun_nondim

    sys_se = cr3bp.CR3BPSystem(
        mu=mu_se, primary="Sun", secondary="Earth", l_km=qbcp_sys.a_sun_nondim * 384400.0, t_s=1.0
    )

    # 1. Correct Sun-Earth L2 Lyapunov orbit
    x_earth = 1.0 - mu_se
    c_target_se = 3.0008
    orbit_se = correct_symmetric_fixed_jacobi(
        sys_se,
        x0_guess=x_earth + 0.010,
        jacobi=c_target_se,
        period_guess=3.1,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )
    print(f"Sun-Earth Lyapunov: x0={orbit_se.x0:.8f}, Period={orbit_se.period:.8f}")

    # 2. Correct Sun-Earth L2 Torus in QBCP (via BCR4BP seeding)
    n_samples = 5
    n_modes = 2
    x0_se_bcr, phase_pin_se, amplitude_pin_se = se_lyapunov_to_bcr4bp_torus_seed(
        orbit_se, qbcp_sys, mu_se, n_samples=n_samples
    )
    print("Correcting Sun-Earth L2 torus in BCR4BP first...")
    torus_bcr = correct_bcr4bp_torus(
        qbcp_sys, x0_se_bcr, n_modes, n_samples, phase_pin_se, amplitude_pin_se, tol=1e-6
    )

    # Flatten the corrected BCR4BP torus coefficients to use as a seed for QBCP
    coeffs_bcr = torus_bcr.fourier_coeffs
    x0_se_qbcp = np.zeros(6 + 12 * n_modes + 1)
    x0_se_qbcp[0:6] = np.real(coeffs_bcr[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        x0_se_qbcp[i0 : i0 + 6] = np.real(coeffs_bcr[n, :])
        x0_se_qbcp[i0 + 6 : i0 + 12] = np.imag(coeffs_bcr[n, :])
    x0_se_qbcp[-1] = torus_bcr.rho

    print("Correcting Sun-Earth L2 torus in QBCP...")
    torus_se = correct_qbcp_torus(
        qbcp_sys, x0_se_qbcp, n_modes, n_samples, phase_pin_se, amplitude_pin_se, tol=1e-3
    )
    print(f"Sun-Earth Torus Residual: {torus_se.invariance_residual:.2e}")

    print("\n==================================================")
    print("STEP 2: Sourcing and correcting Earth-Moon L2 Torus")
    print("==================================================")

    sys_em = cr3bp.CR3BPSystem(
        mu=qbcp_sys.mu, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0
    )

    # 1. Correct Earth-Moon L2 Lyapunov orbit
    c_target_em = 3.17
    orbit_em = correct_symmetric_fixed_jacobi(
        sys_em,
        x0_guess=1.16,
        jacobi=c_target_em,
        period_guess=3.4,
        ydot0_sign=1.0,
        half_crossings=1,
        tol=1e-8,
    )
    print(f"Earth-Moon Lyapunov: x0={orbit_em.x0:.8f}, Period={orbit_em.period:.8f}")

    # 2. Correct Earth-Moon L2 Torus in QBCP (multi-resolution)
    sol_em = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit_em.period),
        np.array([orbit_em.x0, 0.0, 0.0, 0.0, orbit_em.ydot0, 0.0]),
        args=(qbcp_sys.mu,),
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, orbit_em.period, n_samples, endpoint=False),
    )

    u_samples_em = np.zeros((n_samples, 6))
    for j in range(n_samples):
        u_samples_em[j] = sol_em.y[:, j]

    coeffs_em = np.fft.fft(u_samples_em, axis=0) / n_samples
    rho_em = (2.0 * math.pi * t_s / orbit_em.period) % (2.0 * math.pi)
    if rho_em > math.pi:
        rho_em -= 2.0 * math.pi

    n_unk = 6 + 12 * n_modes
    x0_em = np.zeros(n_unk + 1)
    x0_em[0:6] = np.real(coeffs_em[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        x0_em[i0 : i0 + 6] = np.real(coeffs_em[n, :])
        x0_em[i0 + 6 : i0 + 12] = np.imag(coeffs_em[n, :])
    x0_em[-1] = rho_em

    # Phase-pin gauge: pin Im(c_1[idx])=0; gauge derivative w.r.t. a circle
    # rotation is Re(c_1[idx]), so idx must have a large real part or the gauge
    # is singular and the corrector stalls (the #544 Earth-Moon L2 bug). Pick
    # argmax|Re|, matching qp_tori.correct_qp_torus.
    phase_pin_em = int(np.argmax(np.abs(np.real(coeffs_em[1, :]))))
    amplitude_pin_em = float(np.linalg.norm(coeffs_em[1, :]))

    # Correct EM L2 torus in BCR4BP first (which is physically consistent under mu_sun scaling)
    sys_mu0_bcr = bcr4bp.BCR4BPSystem(
        mu=qbcp_sys.mu,
        mu_sun=0.0,
        a_sun_nondim=qbcp_sys.a_sun_nondim,
        omega_sun_nondim=qbcp_sys.omega_sun_nondim,
        theta_sun0=qbcp_sys.theta_sun0,
    )

    init_res = np.linalg.norm(
        bcr4bp_torus_residual(
            x0_em, sys_mu0_bcr, n_modes, n_samples, phase_pin_em, amplitude_pin_em
        )
    )
    print(f"Initial EM seed residual in BCR4BP at mu_sun=0: {init_res:.2e}")

    print("Correcting Earth-Moon L2 torus at mu_sun=0 in BCR4BP...")
    res_gmos_em = least_squares(
        bcr4bp_torus_residual,
        x0_em,
        args=(sys_mu0_bcr, n_modes, n_samples, phase_pin_em, amplitude_pin_em),
        kwargs={"rtol": 1e-6, "atol": 1e-6},
        method="lm",
        xtol=1e-6,
        ftol=1e-6,
    )
    norm_res = np.linalg.norm(res_gmos_em.fun)
    print(f"Corrected EM torus residual in BCR4BP at mu_sun=0: {norm_res:.2e}")

    print("Stepping mu_sun to real Sun mass in BCR4BP...")
    x_curr = res_gmos_em.x
    mu_sun_steps = np.linspace(0.0, qbcp_sys.mu_sun, 5)
    for step_mu_sun in mu_sun_steps[1:]:
        sys_step = bcr4bp.BCR4BPSystem(
            mu=qbcp_sys.mu,
            mu_sun=step_mu_sun,
            a_sun_nondim=qbcp_sys.a_sun_nondim,
            omega_sun_nondim=qbcp_sys.omega_sun_nondim,
            theta_sun0=qbcp_sys.theta_sun0,
        )
        res_gmos_em = least_squares(
            bcr4bp_torus_residual,
            x_curr,
            args=(sys_step, n_modes, n_samples, phase_pin_em, amplitude_pin_em),
            kwargs={"rtol": 1e-6, "atol": 1e-6},
            method="lm",
            xtol=1e-6,
            ftol=1e-6,
        )
        x_curr = res_gmos_em.x
        print(f"  mu_sun={step_mu_sun:.2e}: residual={np.linalg.norm(res_gmos_em.fun):.2e}")

    print("Evaluating initial guess in QBCP dynamics...")
    from cyclerfinder.genome.qbcp_torus import evaluate_invariant_circle

    coeffs_u_diag, _ = x_curr[:-1], x_curr[-1]
    n_total_diag = 2 * n_modes + 1
    coeffs_diag = np.zeros((n_total_diag, 6), dtype=np.complex128)
    coeffs_diag[0, :] = coeffs_u_diag[0:6].astype(np.complex128)
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        coeffs_diag[n, :] = coeffs_u_diag[i0 : i0 + 6] + 1j * coeffs_u_diag[i0 + 6 : i0 + 12]
        coeffs_diag[n_total_diag - n, :] = np.conj(coeffs_diag[n, :])
    thetas_diag = 2 * math.pi * np.arange(n_samples) / n_samples
    u_s_diag = evaluate_invariant_circle(coeffs_diag, thetas_diag)
    for j in range(n_samples):
        try:
            _, states_pv = qbcp.propagate_qbcp_pv(
                u_s_diag[j], (0.0, t_s), qbcp_sys, with_stm=False, rtol=1e-8, atol=1e-8
            )
            print(f"  Sample {j}: OK, final state norm={np.linalg.norm(states_pv[-1]):.2f}")
        except RuntimeError as e:
            print(f"  Sample {j}: Failed propagation: {e}")

    print("Running final correction in QBCP at n_modes=2...")
    torus_em = correct_qbcp_torus(
        qbcp_sys, x_curr, n_modes, n_samples, phase_pin_em, amplitude_pin_em, tol=1.5e-1
    )
    print(f"Earth-Moon Torus Residual: {torus_em.invariance_residual:.2e}")

    print("\n==================================================")
    print("STEP 3: Sweeping Manifold Crossings at Section x = 2.0")
    print("==================================================")

    surface_x = 2.0
    t_max = 35.0

    print("Sweeping unstable manifold of Sun-Earth L2 Torus...")
    _, endpoints_u, _ = build_custom_grid(
        torus_se, "unstable", -1.0, t_max, surface_x, n_long=16, n_lat=16
    )

    print("Sweeping stable manifold of Earth-Moon L2 Torus...")
    _, endpoints_s, _ = build_custom_grid(
        torus_em, "stable", 1.0, t_max, surface_x, n_long=16, n_lat=16
    )

    # Flatten the grids
    flat_u = endpoints_u.reshape(-1, 6)
    flat_s = endpoints_s.reshape(-1, 6)

    # Keep only finite crossings
    valid_u = flat_u[np.isfinite(flat_u[:, 0])]
    valid_s = flat_s[np.isfinite(flat_s[:, 0])]

    print(f"\nManifold crossings at x = {surface_x}:")
    print(f"  Sun-Earth L2 Unstable: {len(valid_u)} / 256")
    print(f"  Earth-Moon L2 Stable:  {len(valid_s)} / 256")

    print("\nSample Sun-Earth L2 Unstable crossings:")
    for idx, u_cross in enumerate(valid_u[:5]):
        print(
            f"  #{idx}: y={u_cross[1]:.4f}, z={u_cross[2]:.4f}, "
            f"vx={u_cross[3]:.4f}, vy={u_cross[4]:.4f}, vz={u_cross[5]:.4f}"
        )

    print("\nSample Earth-Moon L2 Stable crossings:")
    for idx, s_cross in enumerate(valid_s[:5]):
        print(
            f"  #{idx}: y={s_cross[1]:.4f}, z={s_cross[2]:.4f}, "
            f"vx={s_cross[3]:.4f}, vy={s_cross[4]:.4f}, vz={s_cross[5]:.4f}"
        )

    min_dist = float("inf")
    best_pair = None

    for u in valid_u:
        for s in valid_s:
            pos_dist = np.linalg.norm(u[1:3] - s[1:3])  # y, z
            vel_dist = np.linalg.norm(u[3:6] - s[3:6])  # vx, vy, vz
            total_dist = pos_dist + vel_dist
            if total_dist < min_dist:
                min_dist = total_dist
                best_pair = (u, s)

    if best_pair is not None:
        u, s = best_pair
        print("\nCoarse connection candidate found:")
        pos_gap_km = np.linalg.norm(u[1:3] - s[1:3]) * 384400.0
        vel_gap_ms = np.linalg.norm(u[3:6] - s[3:6]) * 1024.0
        print(f"  Coarse Position Gap: {pos_gap_km:.2f} km")
        print(f"  Coarse Velocity Gap: {vel_gap_ms:.2f} m/s")
    else:
        print("\nNo crossing pair found.")
        return

    print("\n==================================================")
    print("STEP 4: Refining Coherent Connections")
    print("==================================================")

    # Generate grid list for matching
    crossings_list_u = []
    thetas_long = 2.0 * math.pi * np.arange(16) / 16
    thetas_lat = 2.0 * math.pi * np.arange(16) / 16

    print("Building lists of crossings...")
    for tl in thetas_long:
        for tt in thetas_lat:
            state = evaluate_qbcp_torus(torus_se, tl, tt)
            t0_u = tl / torus_se.omega_long
            try:
                _, states_pv = qbcp.propagate_qbcp_pv(
                    state,
                    (t0_u, t0_u + torus_se.t_strob),
                    torus_se.system,
                    with_stm=True,
                    rtol=1e-8,
                    atol=1e-8,
                )
            except RuntimeError:
                continue
            stm_pv = states_pv[-1, 6:].reshape((6, 6))
            stab = local_stability(state, stm_pv, hyperbolicity_tol=1e-4)
            vec = stab.vec_u
            if vec is not None:
                if vec[0] * -1.0 < 0.0:
                    vec = -vec
                perturbed = state + 1e-5 * vec
                res = propagate_qbcp_to_section(
                    torus_se.system,
                    perturbed,
                    1.0,
                    surface_x,
                    t_max,
                    1e-8,
                    1e-8,
                    t0_u,
                )
                if res is not None:
                    crossings_list_u.append(
                        {"state": res[0], "time": res[1], "theta_long": tl, "theta_trans": tt}
                    )

    crossings_list_s = []
    for tl in thetas_long:
        for tt in thetas_lat:
            state = evaluate_qbcp_torus(torus_em, tl, tt)
            t0_s = tl / torus_em.omega_long
            try:
                _, states_pv = qbcp.propagate_qbcp_pv(
                    state,
                    (t0_s, t0_s + torus_em.t_strob),
                    torus_em.system,
                    with_stm=True,
                    rtol=1e-8,
                    atol=1e-8,
                )
            except RuntimeError:
                continue
            stm_pv = states_pv[-1, 6:].reshape((6, 6))
            stab = local_stability(state, stm_pv, hyperbolicity_tol=1e-4)
            vec = stab.vec_s
            if vec is not None:
                if vec[0] * 1.0 < 0.0:
                    vec = -vec
                perturbed = state + 1e-5 * vec
                res = propagate_qbcp_to_section(
                    torus_em.system,
                    perturbed,
                    -1.0,
                    surface_x,
                    t_max,
                    1e-8,
                    1e-8,
                    t0_s,
                )
                if res is not None:
                    crossings_list_s.append(
                        {"state": res[0], "time": res[1], "theta_long": tl, "theta_trans": tt}
                    )

    candidates = []
    for uc in crossings_list_u:
        for sc in crossings_list_s:
            pos_dist = np.linalg.norm(uc["state"][1:3] - sc["state"][1:3])
            time_diff = (uc["time"] - sc["time"]) % t_s
            if time_diff > 0.5 * t_s:
                time_diff -= t_s

            if pos_dist < 0.4 and abs(time_diff) < 1.5:
                candidates.append((uc, sc, pos_dist, time_diff))

    print(f"Found {len(candidates)} matching coarse candidates.")
    candidates.sort(key=lambda x: x[2])

    refined_count = 0
    for idx, (uc, sc, pos_dist, time_diff) in enumerate(candidates[:5]):
        print(f"\nRefining Candidate #{idx + 1}:")
        print(f"  Coarse Position Gap: {pos_dist * 384400.0:.2f} km")
        print(f"  Coarse Time Gap:     {time_diff:.4f} rad")

        u_pt = get_manifold_point_with_vec(
            torus_se, uc["theta_long"], uc["theta_trans"], "unstable", 1e-5, -1.0
        )
        s_pt = get_manifold_point_with_vec(
            torus_em, sc["theta_long"], sc["theta_trans"], "stable", 1e-5, 1.0
        )
        if u_pt is None or s_pt is None:
            continue
        v_ref_u = u_pt["vec"]
        v_ref_s = s_pt["vec"]

        p0 = np.array([uc["theta_long"], uc["theta_trans"], sc["theta_long"], sc["theta_trans"]])

        def residual_matching(p, vr_u=v_ref_u, vr_s=v_ref_s):
            t_long_u, t_trans_u, t_long_s, t_trans_s = p
            state_u = evaluate_qbcp_torus(torus_se, t_long_u, t_trans_u)
            state_s = evaluate_qbcp_torus(torus_em, t_long_s, t_trans_s)

            perturbed_u = state_u + 1e-5 * vr_u
            perturbed_s = state_s + 1e-5 * vr_s

            pt_u = propagate_qbcp_to_section(
                torus_se.system,
                perturbed_u,
                1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                t_long_u / torus_se.omega_long,
            )
            pt_s = propagate_qbcp_to_section(
                torus_em.system,
                perturbed_s,
                -1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                t_long_s / torus_em.omega_long,
            )

            if pt_u is None or pt_s is None:
                return np.array([1e3, 1e3, 1e3])

            y_u, z_u = pt_u[0][1], pt_u[0][2]
            y_s, z_s = pt_s[0][1], pt_s[0][2]

            t_diff = (pt_u[1] - pt_s[1]) % t_s
            if t_diff > 0.5 * t_s:
                t_diff -= t_s

            return np.array([y_u - y_s, z_u - z_s, t_diff])

        res = least_squares(
            residual_matching,
            p0,
            bounds=([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            xtol=1e-5,
            ftol=1e-5,
        )

        final_res = np.linalg.norm(res.fun)
        print(f"  Optimizer finished: success={res.success}, message={res.message}")
        print(f"    Final matching residual norm: {final_res:.2e}")

        if res.success and final_res < 1.0:
            final_p = res.x
            state_u = evaluate_qbcp_torus(torus_se, final_p[0], final_p[1])
            state_s = evaluate_qbcp_torus(torus_em, final_p[2], final_p[3])
            perturbed_u = state_u + 1e-5 * v_ref_u
            perturbed_s = state_s + 1e-5 * v_ref_s

            pt_u = propagate_qbcp_to_section(
                torus_se.system,
                perturbed_u,
                1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                final_p[0] / torus_se.omega_long,
            )
            pt_s = propagate_qbcp_to_section(
                torus_em.system,
                perturbed_s,
                -1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                final_p[2] / torus_em.omega_long,
            )

            if pt_u is not None and pt_s is not None:
                u_state = pt_u[0]
                s_state = pt_s[0]

                final_pos_gap = np.linalg.norm(u_state[1:3] - s_state[1:3])
                final_vel_gap = np.linalg.norm(u_state[3:6] - s_state[3:6])
                pos_gap_km = final_pos_gap * 384400.0
                vel_gap_ms = final_vel_gap * 1024.0

                print("  Refined candidate properties:")
                print(f"    Position Gap: {pos_gap_km:.2f} km")
                print(f"    Velocity Gap: {vel_gap_ms:.2f} m/s")
                print(f"    Crossing Time: {pt_u[1]:.4f} TU")

                if final_res < 5e-2:
                    refined_count += 1
                    print("    -> Met the threshold for a highly accurate connection")

    if refined_count == 0:
        print("\nCould not refine any candidate below the strict threshold (5e-2).")


if __name__ == "__main__":
    main()
