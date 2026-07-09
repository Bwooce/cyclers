#!/usr/bin/env python3
"""Search runner to find exact time-coherent connections between Sun-Earth L2 and Earth-Moon L2.

Solves the matching equations:
  y_U - y_S = 0
  z_U - z_S = 0
  (t_U - t_S) mod T_s = 0
using a multi-variable root finder, then computes the delta-V.
All exclamation marks are avoided in comments and docstrings.
"""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_torus import (
    bcr4bp_torus_residual,
    correct_bcr4bp_torus,
    evaluate_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.genome.qp_torus_manifold import local_stability
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi


def propagate_to_section(system, state0, direction, surface_x, t_max, rtol, atol, t0):
    """Integrate state0 from t0 and return (crossing_state, crossing_time) or None."""

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
    return np.asarray(sol.y_events[0][0], dtype=np.float64), float(sol.t_events[0][0])


def get_manifold_point(
    torus, theta_long, theta_trans, branch, eps, sign, direction, surface_x, t_max
):
    """Evaluate torus point, compute eigenvector, perturb, and propagate to section."""
    t0 = theta_long / torus.omega_long
    state = evaluate_bcr4bp_torus(torus, float(theta_long), float(theta_trans))

    # Propagate to get STM
    arc = bcr4bp.propagate_bcr4bp(
        torus.system, state, torus.t_strob, with_stm=True, t0=t0, rtol=1e-8, atol=1e-8
    )
    if arc.stm is None:
        return None

    stab = local_stability(
        state,
        arc.stm,
        hyperbolicity_tol=1e-4,
    )
    vec = stab.vec_u if branch == "unstable" else stab.vec_s
    if vec is None:
        return None

    if vec[0] * sign < 0.0:
        vec = -vec

    perturbed = state + eps * vec
    res = propagate_to_section(torus.system, perturbed, direction, surface_x, t_max, 1e-8, 1e-8, t0)
    if res is None:
        return None
    crossing_state, crossing_time = res
    return {
        "state": crossing_state,
        "time": crossing_time,
        "theta_long": theta_long,
        "theta_trans": theta_trans,
    }


def get_manifold_point_with_vec(torus, theta_long, theta_trans, branch, eps, sign):
    """Evaluate torus point and return its state, time, and eigenvector."""
    t0 = theta_long / torus.omega_long
    state = evaluate_bcr4bp_torus(torus, float(theta_long), float(theta_trans))

    arc = bcr4bp.propagate_bcr4bp(
        torus.system, state, torus.t_strob, with_stm=True, t0=t0, rtol=1e-8, atol=1e-8
    )
    if arc.stm is None:
        return None

    stab = local_stability(
        state,
        arc.stm,
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

    # 2. Correct Sun-Earth L2 Torus in BCR4BP
    n_samples = 5
    n_modes = 2
    x0_se, phase_pin_se, amplitude_pin_se = se_lyapunov_to_bcr4bp_torus_seed(
        orbit_se, bcr_sys, mu_SE, n_samples=n_samples
    )
    torus_se = correct_bcr4bp_torus(
        bcr_sys, x0_se, n_modes, n_samples, phase_pin_se, amplitude_pin_se, tol=1e-6
    )

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

    res_gmos_em = least_squares(
        bcr4bp_torus_residual,
        x0_em,
        args=(sys_mu_sun0, n_modes, n_samples, phase_pin_em, amplitude_pin_em),
        kwargs={"rtol": 1e-6, "atol": 1e-6},
        method="lm",
        xtol=1e-6,
        ftol=1e-6,
    )

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

    torus_em = correct_bcr4bp_torus(
        bcr_sys, x0_new, n_modes_new, n_samples_new, phase_pin_em, amplitude_pin_em, tol=1e-2
    )

    print("\n==================================================")
    print("STEP 3: Finding Coherent Connections")
    print("==================================================")

    surface_x = 2.0
    t_max = 35.0

    # Generate coarse grids
    n_long = 16
    n_lat = 16
    thetas_long = 2.0 * math.pi * np.arange(n_long) / n_long
    thetas_lat = 2.0 * math.pi * np.arange(n_lat) / n_lat

    print("Building grids...")
    crossings_u = []
    for tl in thetas_long:
        for tt in thetas_lat:
            pt = get_manifold_point(torus_se, tl, tt, "unstable", 1e-5, -1.0, 1.0, surface_x, t_max)
            if pt is not None:
                crossings_u.append(pt)

    crossings_s = []
    for tl in thetas_long:
        for tt in thetas_lat:
            pt = get_manifold_point(torus_em, tl, tt, "stable", 1e-5, 1.0, -1.0, surface_x, t_max)
            if pt is not None:
                crossings_s.append(pt)

    print(f"Generated {len(crossings_u)} unstable and {len(crossings_s)} stable crossings.")

    # Match candidates
    candidates = []
    for u in crossings_u:
        for s in crossings_s:
            pos_dist = np.linalg.norm(u["state"][1:3] - s["state"][1:3])
            time_diff = (u["time"] - s["time"]) % T_s
            if time_diff > 0.5 * T_s:
                time_diff -= T_s

            # If position gap < 0.4 and time gap < 1.5
            if pos_dist < 0.4 and abs(time_diff) < 1.5:
                candidates.append((u, s, pos_dist, time_diff))

    print(f"Found {len(candidates)} matching coarse candidates.")

    # Sort candidates by position gap
    candidates.sort(key=lambda x: x[2])

    # Run root finder on the best candidates
    refined_count = 0
    for idx, (u, s, pos_dist, time_diff) in enumerate(candidates[:5]):
        print(f"\nRefining Candidate #{idx + 1}:")
        print(f"  Coarse Position Gap: {pos_dist * 384400.0:.2f} km")
        print(f"  Coarse Time Gap:     {time_diff:.4f} rad")

        # Get reference eigenvectors at guess points
        u_pt = get_manifold_point_with_vec(
            torus_se, u["theta_long"], u["theta_trans"], "unstable", 1e-5, -1.0
        )
        s_pt = get_manifold_point_with_vec(
            torus_em, s["theta_long"], s["theta_trans"], "stable", 1e-5, 1.0
        )
        if u_pt is None or s_pt is None:
            continue
        v_ref_u = u_pt["vec"]
        v_ref_s = s_pt["vec"]

        # Initial guess for the 4 angles:
        # [theta_long_u, theta_trans_u, theta_long_s, theta_trans_s]
        p0 = np.array([u["theta_long"], u["theta_trans"], s["theta_long"], s["theta_trans"]])

        def residual_matching(p, vr_u=v_ref_u, vr_s=v_ref_s):
            t_long_u, t_trans_u, t_long_s, t_trans_s = p

            # Evaluate torus state
            state_u = evaluate_bcr4bp_torus(torus_se, t_long_u, t_trans_u)
            state_s = evaluate_bcr4bp_torus(torus_em, t_long_s, t_trans_s)

            # Perturb along frozen reference vectors
            perturbed_u = state_u + 1e-5 * vr_u
            perturbed_s = state_s + 1e-5 * vr_s

            # Propagate to section
            pt_u = propagate_to_section(
                torus_se.system,
                perturbed_u,
                1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                t_long_u / torus_se.omega_long,
            )
            pt_s = propagate_to_section(
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

            t_diff = (pt_u[1] - pt_s[1]) % T_s
            if t_diff > 0.5 * T_s:
                t_diff -= T_s

            return np.array([y_u - y_s, z_u - z_s, t_diff])

        res = least_squares(
            residual_matching,
            p0,
            bounds=([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            xtol=1e-5,
            ftol=1e-5,
        )

        initial_res = np.linalg.norm(residual_matching(p0))
        final_res = np.linalg.norm(res.fun)
        print(f"  Optimizer finished: success={res.success}, message={res.message}")
        print(f"    Initial residual norm: {initial_res:.2e}")
        print(f"    Final residual norm:   {final_res:.2e}")

        if res.success and final_res < 1e-2:
            refined_count += 1
            final_p = res.x

            state_u = evaluate_bcr4bp_torus(torus_se, final_p[0], final_p[1])
            state_s = evaluate_bcr4bp_torus(torus_em, final_p[2], final_p[3])
            perturbed_u = state_u + 1e-5 * v_ref_u
            perturbed_s = state_s + 1e-5 * v_ref_s

            pt_u = propagate_to_section(
                torus_se.system,
                perturbed_u,
                1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                final_p[0] / torus_se.omega_long,
            )
            pt_s = propagate_to_section(
                torus_em.system,
                perturbed_s,
                -1.0,
                surface_x,
                t_max,
                1e-8,
                1e-8,
                final_p[2] / torus_em.omega_long,
            )

            u_state = pt_u[0]
            s_state = pt_s[0]

            final_pos_gap = np.linalg.norm(u_state[1:3] - s_state[1:3])
            final_vel_gap = np.linalg.norm(u_state[3:6] - s_state[3:6])

            print("  SUCCESSFULLY REFINED:")
            print(f"    Position Gap: {final_pos_gap * 384400.0:.2f} km")
            print(f"    Velocity Gap: {final_vel_gap * 1024.0:.2f} m/s")
            print(f"    Crossing Time: {pt_u[1]:.4f} TU")

    if refined_count == 0:
        print("\nCould not refine any candidate to zero position gap.")


if __name__ == "__main__":
    main()
