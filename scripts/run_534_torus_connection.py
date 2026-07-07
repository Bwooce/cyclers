#!/usr/bin/env python3
"""Search runner for Task #534: Earth-Moon L1/L2 Quasi-Halo Torus Connection Search.

Uses parallelized grid generation to sweep through stable/unstable manifold crossings
at the surface of section, checking for topological intersections via linking numbers.
"""

import concurrent.futures
import math

import numpy as np
import scipy.integrate

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import correct_qp_torus
from cyclerfinder.genome.qp_torus_heteroclinic import build_manifold_grids, scan_linking_number
from cyclerfinder.genome.qp_torus_manifold import (
    ManifoldGrid,
    local_stability,
    torus_point_stm,
)
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.nrho_continuation import correct_symmetric_nrho


def ydot_crossing_state(
    system, state0, direction, surface_x, t_max, rtol, atol, target_ydot_sign=-1.0
):
    """Custom crossing function to collect the first crossing matching a specific ydot sign."""

    def x_event(t, y, mu):
        return y[0] - surface_x

    x_event.terminal = False
    x_event.direction = 0.0

    t_span = (0.0, math.copysign(t_max, direction))
    sol = scipy.integrate.solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        state0,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=x_event,
        max_step=abs(t_max) / 20.0,
    )
    if sol.t_events is None or len(sol.t_events[0]) == 0:
        return None
    for y_event in sol.y_events[0]:
        if y_event[4] * target_ydot_sign > 0.0:
            return np.asarray(y_event, dtype=np.float64)
    return None


def compute_custom_row(args):
    i, theta_long, thetas_lat, torus, branch, eps, sign, t_max, surface_x = args
    n_lat = len(thetas_lat)
    row_origins = np.zeros((n_lat, 2))
    row_endpoints = np.full((n_lat, 6), np.nan)
    row_hyperbolic = np.zeros(n_lat, dtype=bool)

    prev_vec = None
    for j, theta_trans in enumerate(thetas_lat):
        row_origins[j, 0] = theta_long
        row_origins[j, 1] = theta_trans
        state, stm = torus_point_stm(torus, float(theta_long), float(theta_trans))
        stab = local_stability(
            state,
            stm,
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
        perturbed = state + eps * vec
        direction = 1.0 if branch == "unstable" else -1.0
        crossing = ydot_crossing_state(
            torus.system,
            perturbed,
            direction=direction,
            surface_x=surface_x,
            t_max=t_max,
            rtol=1e-12,
            atol=1e-12,
            target_ydot_sign=-1.0,
        )
        if crossing is not None:
            row_endpoints[j, :] = crossing

    return i, row_origins, row_endpoints, row_hyperbolic


def build_custom_grid(torus, branch, sign, t_max, surface_x):
    n_long, n_lat = 16, 16
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

    return ManifoldGrid(origins=origins, endpoints=endpoints, hyperbolic=hyperbolic)


def main():
    sysm = cr3bp.CR3BPSystem(
        mu=0.012153643, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0
    )

    print("--------------------------------------------------")
    print("STEP 1: Sourcing and correcting L1/L2 NRHO orbits...")
    print("--------------------------------------------------")
    print("Correcting L1 NRHO...")
    res_l1 = correct_symmetric_nrho(sysm, 0.836314, -0.145689, 0.231349, 2.751255)
    mono_l1 = monodromy(
        sysm, np.array([res_l1.x0, 0.0, res_l1.z0, 0.0, res_l1.ydot0, 0.0]), res_l1.T_TU
    )
    eigs_l1 = floquet_multipliers(mono_l1)
    cands_l1 = [e for e in eigs_l1 if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    torus_l1 = correct_qp_torus(
        sysm,
        np.array([res_l1.x0, 0.0, res_l1.z0, 0.0, res_l1.ydot0, 0.0]),
        res_l1.T_TU,
        (cands_l1[0], cands_l1[1]),
        k=5,
        n_trans=8,
        initial_torus_amplitude=5e-4,
        tol=1e-5,
    )
    print(f"L1 Torus Converged: Inv. Residual = {torus_l1.invariance_residual:.2e}")

    print("Correcting L2 NRHO...")
    res_l2 = correct_symmetric_nrho(sysm, 1.023731, 0.183250, -0.106950, 1.533637)
    mono_l2 = monodromy(
        sysm, np.array([res_l2.x0, 0.0, res_l2.z0, 0.0, res_l2.ydot0, 0.0]), res_l2.T_TU
    )
    eigs_l2 = floquet_multipliers(mono_l2)
    cands_l2 = [e for e in eigs_l2 if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    torus_l2 = correct_qp_torus(
        sysm,
        np.array([res_l2.x0, 0.0, res_l2.z0, 0.0, res_l2.ydot0, 0.0]),
        res_l2.T_TU,
        (cands_l2[0], cands_l2[1]),
        k=7,
        n_trans=8,
        initial_torus_amplitude=5e-4,
        tol=1e-5,
    )
    print(f"L2 Torus Converged: Inv. Residual = {torus_l2.invariance_residual:.2e}")

    surface_x = 1.0 - sysm.mu

    print("\n--------------------------------------------------")
    print("STEP 2: Running Standard First-Crossing Search...")
    print("--------------------------------------------------")
    # This calls the library sequential build_manifold_grids, which now uses our
    # parallelized torus_manifold_grid backend internally.
    stable_grid, unstable_grid = build_manifold_grids(
        torus_l2,
        torus_l1,
        n_long=16,
        n_lat=16,
        eps=1e-5,
        surface_x=surface_x,
        t_max=15.0,
        stable_sign=-1.0,
        unstable_sign=1.0,
    )

    n_cross_u = np.sum(np.isfinite(unstable_grid.endpoints[:, :, 0]))
    n_cross_s = np.sum(np.isfinite(stable_grid.endpoints[:, :, 0]))
    print(f"Crossings found (first-crossing): unstable={n_cross_u}/256, stable={n_cross_s}/256")

    z_s = stable_grid.endpoints[:, :, 2]
    z_u = unstable_grid.endpoints[:, :, 2]

    # Check overlap in z
    z_s_finite = z_s[np.isfinite(z_s)]
    z_u_finite = z_u[np.isfinite(z_u)]
    if len(z_s_finite) > 0 and len(z_u_finite) > 0:
        overlap_min = max(z_s_finite.min(), z_u_finite.min())
        overlap_max = min(z_s_finite.max(), z_u_finite.max())
        print(
            f"z ranges: unstable [{z_u_finite.min():.5f}, {z_u_finite.max():.5f}], "
            f"stable [{z_s_finite.min():.5f}, {z_s_finite.max():.5f}]"
        )
        print(f"Overlap z range: {overlap_min:.5f} to {overlap_max:.5f}")

        if overlap_min < overlap_max:
            d_values = np.linspace(overlap_min, overlap_max, 50)
            result = scan_linking_number(
                stable_grid,
                unstable_grid,
                scanning_component="z",
                curve_components=("y", "ydot", "zdot"),
                d_values=d_values,
            )
            print("Linking numbers:", result.linking_numbers.tolist())
            print("Sign change locations:", result.sign_change_locations())
        else:
            print("No z-range overlap in first-crossing search.")
    else:
        print("No finite crossings found in first-crossing search.")

    print("\n--------------------------------------------------")
    print("STEP 3: Running Advanced ydot < 0 Multi-Crossing Search...")
    print("--------------------------------------------------")
    # For this search, we use a custom builder that forces ydot < 0 to bring
    # L1 and L2 into phase alignment

    print("Integrating L1 unstable manifold (t_max = 15.0)...")
    unstable_grid_adv = build_custom_grid(torus_l1, "unstable", 1.0, 15.0, surface_x)
    print("Integrating L2 stable manifold (t_max = 25.0)...")
    stable_grid_adv = build_custom_grid(torus_l2, "stable", -1.0, 25.0, surface_x)

    n_cross_u_adv = np.sum(np.isfinite(unstable_grid_adv.endpoints[:, :, 0]))
    n_cross_s_adv = np.sum(np.isfinite(stable_grid_adv.endpoints[:, :, 0]))
    print(
        f"Crossings found (ydot < 0 filter): "
        f"unstable={n_cross_u_adv}/256, stable={n_cross_s_adv}/256"
    )

    z_s_adv = stable_grid_adv.endpoints[:, :, 2]
    z_u_adv = unstable_grid_adv.endpoints[:, :, 2]
    z_s_adv_finite = z_s_adv[np.isfinite(z_s_adv)]
    z_u_adv_finite = z_u_adv[np.isfinite(z_u_adv)]

    if len(z_s_adv_finite) > 0 and len(z_u_adv_finite) > 0:
        overlap_min = max(z_s_adv_finite.min(), z_u_adv_finite.min())
        overlap_max = min(z_s_adv_finite.max(), z_u_adv_finite.max())
        print(
            f"z ranges: unstable [{z_u_adv_finite.min():.5f}, {z_u_adv_finite.max():.5f}], "
            f"stable [{z_s_adv_finite.min():.5f}, {z_s_adv_finite.max():.5f}]"
        )
        print(f"Overlap z range: {overlap_min:.5f} to {overlap_max:.5f}")

        if overlap_min < overlap_max:
            d_values = np.linspace(overlap_min, overlap_max, 50)
            result = scan_linking_number(
                stable_grid,
                unstable_grid,
                scanning_component="z",
                curve_components=("y", "ydot", "zdot"),
                d_values=d_values,
            )
            print("Linking numbers:", result.linking_numbers.tolist())
            print("Sign change locations:", result.sign_change_locations())
        else:
            print("No z-range overlap in first-crossing search.")
    else:
        print("No finite crossings found in first-crossing search.")

    print("\n--------------------------------------------------")
    print("STEP 3: Running Advanced ydot < 0 Multi-Crossing Search...")
    print("--------------------------------------------------")
    # For this search, we use a custom builder that forces ydot < 0 to bring
    # L1 and L2 into phase alignment

    print("Integrating L1 unstable manifold (t_max = 15.0)...")
    unstable_grid_adv = build_custom_grid(torus_l1, "unstable", 1.0, 15.0, surface_x)
    print("Integrating L2 stable manifold (t_max = 25.0)...")
    stable_grid_adv = build_custom_grid(torus_l2, "stable", -1.0, 25.0, surface_x)

    n_cross_u_adv = np.sum(np.isfinite(unstable_grid_adv.endpoints[:, :, 0]))
    n_cross_s_adv = np.sum(np.isfinite(stable_grid_adv.endpoints[:, :, 0]))
    print(
        f"Crossings found (ydot < 0 filter): "
        f"unstable={n_cross_u_adv}/256, stable={n_cross_s_adv}/256"
    )

    z_s_adv = stable_grid_adv.endpoints[:, :, 2]
    z_u_adv = unstable_grid_adv.endpoints[:, :, 2]
    z_s_adv_finite = z_s_adv[np.isfinite(z_s_adv)]
    z_u_adv_finite = z_u_adv[np.isfinite(z_u_adv)]

    if len(z_s_adv_finite) > 0 and len(z_u_adv_finite) > 0:
        overlap_min = max(z_s_adv_finite.min(), z_u_adv_finite.min())
        overlap_max = min(z_s_adv_finite.max(), z_u_adv_finite.max())
        print(
            f"z ranges: unstable [{z_u_adv_finite.min():.5f}, {z_u_adv_finite.max():.5f}], "
            f"stable [{z_s_adv_finite.min():.5f}, {z_s_adv_finite.max():.5f}]"
        )
        print(f"Overlap z range: {overlap_min:.5f} to {overlap_max:.5f}")

        if overlap_min < overlap_max:
            d_values = np.linspace(overlap_min, overlap_max, 50)
            result = scan_linking_number(
                stable_grid_adv,
                unstable_grid_adv,
                scanning_component="z",
                curve_components=("y", "ydot", "zdot"),
                d_values=d_values,
            )
            print("Linking numbers:", result.linking_numbers.tolist())
            print("Sign change locations (connections):", result.sign_change_locations())
        else:
            print("No z-range overlap in ydot < 0 search.")
    else:
        print("No finite crossings found in ydot < 0 search.")


if __name__ == "__main__":
    main()
