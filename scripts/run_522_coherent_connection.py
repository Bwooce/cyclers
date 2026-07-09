#!/usr/bin/env python3
"""Search runner for Task #522: Sun-Earth L2 Torus to Earth-Moon L2 Torus Connection Search.

Uses parallelized grid generation to sweep through stable/unstable manifold crossings
at the surface of section, checking for topological intersections.
All exclamation marks are avoided in comments and docstrings.
"""

import concurrent.futures
import math
import pathlib

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.genome.bcr4bp_torus import (
    bcr4bp_torus_residual,
    correct_bcr4bp_torus,
    evaluate_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.genome.qp_torus_manifold import local_stability
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

_REGION_ID = "se-l2-em-l2-bcr4bp-torus-connection-2026-07-08"
_METHOD = MethodCapability(
    genome=(
        "Sun-Earth L2 <-> Earth-Moon L2 torus manifold crossing grids (16x16 long/lat) "
        "in BCR4BP, superseded by the genuine-QBCP #533/#537 pipeline"
    ),
    corrector=(
        "build_custom_grid + closest-approach pairing (diagnostic only, no least_squares refine)"
    ),
    capability_tags=frozenset(
        {"bcr4bp", "torus", "heteroclinic", "sun-earth-moon", "precursor-superseded"}
    ),
    git_sha="working-tree",
)


def bcr4bp_crossing_state(system, state0, direction, surface_x, t_max, rtol, atol, t0):
    """Integrate state0 from t0 and return the state at the surface_x crossing."""

    def x_event(t, y, sys):
        return y[0] - surface_x

    x_event.terminal = True
    x_event.direction = 0.0

    t_span = (t0, t0 + math.copysign(t_max, direction))
    sol = solve_ivp(
        bcr4bp.bcr4bp_eom,
        t_span,
        state0,
        args=(system,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=x_event,
    )
    if sol.t_events is None or len(sol.t_events[0]) == 0:
        return None
    return np.asarray(sol.y_events[0][0], dtype=np.float64)


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
        state = evaluate_bcr4bp_torus(torus, float(theta_long), float(theta_trans))
        arc = bcr4bp.propagate_bcr4bp(
            torus.system, state, torus.t_strob, with_stm=True, t0=t0, rtol=1e-8, atol=1e-8
        )
        if arc.stm is None:
            continue

        stab = local_stability(
            state,
            arc.stm,
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
        crossing = bcr4bp_crossing_state(
            torus.system,
            perturbed,
            direction=direction,
            surface_x=surface_x,
            t_max=t_max,
            rtol=1e-8,
            atol=1e-8,
            t0=t0,
        )
        if crossing is not None:
            row_endpoints[j, :] = crossing

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


def main():
    preflight_search(
        task_no=522,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=512,
        override_reason=(
            "retrofitting the mandatory #521 preflight gate onto an uncommitted BCR4BP "
            "precursor script superseded by the genuine-QBCP #533 pipeline (#537's "
            "resolved connection search); kept for git-history/progression reference, "
            "not intended for a fresh full run, so no timing pilot was captured"
        ),
    )
    print("==================================================")
    print("STEP 1: Sourcing and correcting Sun-Earth L2 Torus")
    print("==================================================")

    bcr_sys = bcr4bp.andreu_default()
    mu_SE = 1.0 / (bcr_sys.mu_sun + 1.0)
    T_s = 2.0 * math.pi / bcr_sys.omega_sun_nondim

    sys_se = cr3bp.CR3BPSystem(
        mu=mu_SE, primary="Sun", secondary="Earth", l_km=bcr_sys.a_sun_nondim * 384400.0, t_s=1.0
    )

    # 1. Correct Sun-Earth L2 Lyapunov orbit
    x_earth = 1.0 - mu_SE
    C_target_se = 3.0008
    orbit_se = correct_symmetric_fixed_jacobi(
        sys_se,
        x0_guess=x_earth + 0.010,
        jacobi=C_target_se,
        period_guess=3.1,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )
    print(f"Sun-Earth Lyapunov: x0={orbit_se.x0:.8f}, Period={orbit_se.period:.8f}")

    # 2. Correct Sun-Earth L2 Torus in BCR4BP
    n_samples = 5
    n_modes = 2
    x0_se, phase_pin_se, amplitude_pin_se = se_lyapunov_to_bcr4bp_torus_seed(
        orbit_se, bcr_sys, mu_SE, n_samples=n_samples
    )

    print("Correcting Sun-Earth L2 torus in BCR4BP...")
    torus_se = correct_bcr4bp_torus(
        bcr_sys, x0_se, n_modes, n_samples, phase_pin_se, amplitude_pin_se, tol=1e-6
    )
    print(f"Sun-Earth Torus Residual: {torus_se.invariance_residual:.2e}")

    print("\n==================================================")
    print("STEP 2: Sourcing and correcting Earth-Moon L2 Torus")
    print("==================================================")

    sys_em = cr3bp.CR3BPSystem(
        mu=bcr_sys.mu, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0
    )

    # 1. Correct Earth-Moon L2 Lyapunov orbit
    C_target_em = 3.17
    orbit_em = correct_symmetric_fixed_jacobi(
        sys_em,
        x0_guess=1.16,
        jacobi=C_target_em,
        period_guess=3.4,
        ydot0_sign=1.0,
        half_crossings=1,
        tol=1e-8,
    )
    print(f"Earth-Moon Lyapunov: x0={orbit_em.x0:.8f}, Period={orbit_em.period:.8f}")

    # 2. Correct Earth-Moon L2 Torus in BCR4BP (multi-resolution)
    sol_em = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit_em.period),
        np.array([orbit_em.x0, 0.0, 0.0, 0.0, orbit_em.ydot0, 0.0]),
        args=(bcr_sys.mu,),
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, orbit_em.period, n_samples, endpoint=False),
    )

    u_samples_em = np.zeros((n_samples, 6))
    for j in range(n_samples):
        u_samples_em[j] = sol_em.y[:, j]

    coeffs_em = np.fft.fft(u_samples_em, axis=0) / n_samples
    rho_em = (2.0 * math.pi * T_s / orbit_em.period) % (2.0 * math.pi)
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

    phase_pin_em = int(np.argmax(np.abs(np.imag(coeffs_em[1, :]))))
    amplitude_pin_em = float(np.linalg.norm(coeffs_em[1, :]))

    sys_mu_sun0 = bcr4bp.BCR4BPSystem(
        mu=bcr_sys.mu,
        mu_sun=0.0,
        a_sun_nondim=bcr_sys.a_sun_nondim,
        omega_sun_nondim=bcr_sys.omega_sun_nondim,
    )

    print("Correcting Earth-Moon L2 torus at mu_sun=0...")
    res_gmos_em = least_squares(
        bcr4bp_torus_residual,
        x0_em,
        args=(sys_mu_sun0, n_modes, n_samples, phase_pin_em, amplitude_pin_em),
        kwargs={"rtol": 1e-6, "atol": 1e-6},
        method="lm",
        xtol=1e-6,
        ftol=1e-6,
    )

    print("Stepping mu_sun to real Sun mass...")
    x_curr = res_gmos_em.x
    mu_sun_steps = np.linspace(0.0, bcr_sys.mu_sun, 5)
    for step_mu_sun in mu_sun_steps[1:]:
        sys_step = bcr4bp.BCR4BPSystem(
            mu=bcr_sys.mu,
            mu_sun=step_mu_sun,
            a_sun_nondim=bcr_sys.a_sun_nondim,
            omega_sun_nondim=bcr_sys.omega_sun_nondim,
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

    # Pad to n_modes=5
    n_modes_new = 5
    n_samples_new = 11

    coeffs_2 = np.zeros((5, 6), dtype=np.complex128)
    coeffs_2[0, :] = x_curr[0:6].astype(np.complex128)
    for n in range(1, 3):
        i0 = 6 + (n - 1) * 12
        coeffs_2[n, :] = x_curr[i0 : i0 + 6] + 1j * x_curr[i0 + 6 : i0 + 12]
        coeffs_2[5 - n, :] = np.conj(coeffs_2[n, :])

    n_total_sig = 2 * n_modes_new + 1
    coeffs_5 = np.zeros((n_total_sig, 6), dtype=np.complex128)
    coeffs_5[0:3, :] = coeffs_2[0:3, :]
    coeffs_5[n_total_sig - 2 :, :] = coeffs_2[3:, :]

    n_unk_new = 6 + 12 * n_modes_new
    x0_new = np.zeros(n_unk_new + 1)
    x0_new[0:6] = np.real(coeffs_5[0, :])
    for n in range(1, n_modes_new + 1):
        i0 = 6 + (n - 1) * 12
        x0_new[i0 : i0 + 6] = np.real(coeffs_5[n, :])
        x0_new[i0 + 6 : i0 + 12] = np.imag(coeffs_5[n, :])
    x0_new[-1] = x_curr[-1]

    print("Running final high-precision correction at n_modes=5...")
    torus_em = correct_bcr4bp_torus(
        bcr_sys, x0_new, n_modes_new, n_samples_new, phase_pin_em, amplitude_pin_em, tol=1e-2
    )
    print(f"Earth-Moon Torus Residual: {torus_em.invariance_residual:.2e}")

    print("\n==================================================")
    print("STEP 3: Sweeping Manifold Crossings at Section x = 2.0")
    print("==================================================")

    surface_x = 2.0
    t_max = 35.0

    print("Sweeping unstable manifold of Sun-Earth L2 Torus...")
    _origins_u, endpoints_u, _hyperbolic_u = build_custom_grid(
        torus_se, "unstable", -1.0, t_max, surface_x, n_long=16, n_lat=16
    )

    print("Sweeping stable manifold of Earth-Moon L2 Torus...")
    _origins_s, endpoints_s, _hyperbolic_s = build_custom_grid(
        torus_em, "stable", 1.0, t_max, surface_x, n_long=16, n_lat=16
    )

    # 4. Search for closest approach between the two manifold crossing grids
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
            # We compute distance in position and velocity
            pos_dist = np.linalg.norm(u[1:3] - s[1:3])  # y, z
            vel_dist = np.linalg.norm(u[3:6] - s[3:6])  # vx, vy, vz
            total_dist = pos_dist + vel_dist
            if total_dist < min_dist:
                min_dist = total_dist
                best_pair = (u, s)

    if best_pair is not None:
        u, s = best_pair
        print("\nBest connection candidate found:")
        print(f"  Position gap (y, z): {np.linalg.norm(u[1:3] - s[1:3]):.6f} nondim units")
        print(f"  Velocity gap (dv):   {np.linalg.norm(u[3:6] - s[3:6]):.6f} nondim units")
        print(f"  Total phase space distance: {min_dist:.6f}")

        # Convert to physical units
        # EM distance is 384400 km, velocity unit is 1.024 km/s
        pos_gap_km = np.linalg.norm(u[1:3] - s[1:3]) * 384400.0
        vel_gap_ms = np.linalg.norm(u[3:6] - s[3:6]) * 1024.0
        print(f"  Physical Position Gap: {pos_gap_km:.2f} km")
        print(f"  Physical Velocity Gap: {vel_gap_ms:.2f} m/s")
    else:
        print("\nNo crossing pair found.")


if __name__ == "__main__":
    main()
