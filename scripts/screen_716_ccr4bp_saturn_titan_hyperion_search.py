"""#716 -- CCR4BP manifold globalization + heteroclinic search: Saturn Titan-Hyperion.

Targeting Saturn-Titan-Hyperion 4:3 (interior) or 1:2 (exterior) resonance.
Titan plays Europa's base-moon role, Hyperion plays Ganymede's outer-perturber role.

Includes Stage-3 dense-mesh confirmation (240 phases) and explicit Hyperion model-fidelity caveat.
Reuses core.ccr4bp, search.variational_ccr4bp_torus, search.ccr4bp_whisker,
search.ccr4bp_manifold_globalize, search.ccr4bp_heteroclinic_search UNMODIFIED.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.ccr4bp_titan_hyperion as th  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402
from cyclerfinder.search.variational_ccr4bp_torus import evaluate_torus_state  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "716_ccr4bp_saturn_titan_hyperion_search"

# Correct physical units for THIS system (base moon = Titan, SMA 1,221,870 km).
_L_KM_TITAN = th.L_KM
_V_UNIT_TITAN_KM_S = th.v_unit_km_s()

_ROBUST_INTEGRATOR_MARGIN_KM = 0.1


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th_half = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom,
            (0.0, th_half),
            y42,
            args=(mu,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


def _log(msg: str, t0: float) -> None:
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


def _corrected_gaps(refined: hs.RefinedConnection) -> dict[str, float]:
    pos_gap_km = float(np.linalg.norm(refined.state_u[:3] - refined.state_s[:3])) * _L_KM_TITAN
    vel_gap_km_s = (
        float(np.linalg.norm(refined.state_u[3:] - refined.state_s[3:])) * _V_UNIT_TITAN_KM_S
    )
    return {"pos_gap_km": pos_gap_km, "vel_gap_km_s": vel_gap_km_s}


def _corrected_off_torus_km(
    torus_u: vt.CCR4BPTorusVariationalResult,
    refined: hs.RefinedConnection,
    theta1_section: float,
) -> float:
    theta1_u_at_t = theta1_section + torus_u.omega1 * refined.t_u
    theta2_u_at_t = refined.theta2_u + torus_u.omega2 * refined.t_u
    torus_pt_u = evaluate_torus_state(torus_u, theta1_u_at_t, theta2_u_at_t)
    planar_u = refined.state_u[[0, 1, 3, 4]]
    return float(np.linalg.norm(torus_pt_u - planar_u)) * _L_KM_TITAN


def _corrected_radau_check(
    torus_u: vt.CCR4BPTorusVariationalResult,
    torus_s: vt.CCR4BPTorusVariationalResult,
    refined: hs.RefinedConnection,
    *,
    lobe_sign_u: float,
    lobe_sign_s: float,
    theta1_section: float,
    n_segments_dir: int,
    rtol: float,
    atol: float,
) -> dict[str, float]:
    ref_vec_u = refined.ref_vec_u
    ref_vec_s = refined.ref_vec_s
    su_radau = hs._radau_manifold_state(
        torus_u,
        "unstable",
        theta1_section,
        refined.theta2_u,
        refined.t_u,
        eps=mg.DEFAULT_EPS,
        lobe_sign=lobe_sign_u,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_u,
    )
    ss_radau = hs._radau_manifold_state(
        torus_s,
        "stable",
        theta1_section,
        refined.theta2_s,
        refined.t_s,
        eps=mg.DEFAULT_EPS,
        lobe_sign=lobe_sign_s,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_s,
    )
    if su_radau is None or ss_radau is None:
        return {"radau_pos_gap_km": float("nan"), "integrator_delta_km": float("inf")}
    radau_pos_gap_km = float(np.linalg.norm(su_radau[:3] - ss_radau[:3])) * _L_KM_TITAN
    corrected_dop853_pos_gap_km = (
        float(np.linalg.norm(refined.state_u[:3] - refined.state_s[:3])) * _L_KM_TITAN
    )
    integrator_delta_km = abs(radau_pos_gap_km - corrected_dop853_pos_gap_km)
    return {"radau_pos_gap_km": radau_pos_gap_km, "integrator_delta_km": integrator_delta_km}


def main() -> int:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _log("=== Task #716: Saturn Titan-Hyperion Hardened CCR4BP Search ===", t0)

    system = th.saturn_titan_hyperion_default()
    _log(
        f"Saturn-Titan-Hyperion CCR4BP parameters: mu={system.mu:.6e}, mu_gan={system.mu_gan:.6e}, "
        f"a_gan={system.a_gan:.4f}, omega_gan={system.omega_gan:.4f}",
        t0,
    )

    # Scaffolding base orbit: 1:2 exterior resonant orbit around Titan
    s0_6d, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    _log(
        f"Base 1:2 orbit converged: x0={s0_6d[0]:.6f}, vy0={s0_6d[4]:.6f}, "
        f"period={period:.4f}, res={res:.2e}",
        t0,
    )

    # Discover/correct CCR4BP torus from base orbit
    torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0_6d,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    _log(
        f"CCR4BP torus built: residual_rms={torus.residual_rms:.4e}, "
        f"closure_residual={torus.closure_residual:.4e}, rho_strob={torus.rho_strob:.4f}",
        t0,
    )

    # Stage 1: Globalize manifold tubes (Dense Mesh: 240 phases)
    n_theta2, n_time, t_max_periods, n_segments_dir = 240, 200, 2.0, 24
    tubes: dict[str, mg.ManifoldTube] = {}
    for branch in ("unstable", "stable"):
        for lobe in (1.0, -1.0):
            key = f"{branch}_{'+' if lobe > 0 else '-'}"
            tubes[key] = mg.globalize_manifold_tube(
                torus,
                branch,
                n_theta2=n_theta2,
                t_max_periods=t_max_periods,
                n_time=n_time,
                n_segments_dir=n_segments_dir,
                lobe_sign=lobe,
            )
            n_valid = int(tubes[key].valid.sum())
            _log(f"Globalized dense mesh {key}: {n_valid}/{n_theta2} phases valid", t0)

    # Stage 2: Coarse search + Refinement + Ghost Guard evaluation
    combo_results: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    best_robust: dict[str, Any] | None = None

    for u_lobe in (1.0, -1.0):
        for s_lobe in (1.0, -1.0):
            u_key = f"unstable_{'+' if u_lobe > 0 else '-'}"
            s_key = f"stable_{'+' if s_lobe > 0 else '-'}"
            candidates = hs.coarse_candidates(
                tubes[u_key], tubes[s_key], n_candidates=5, t_min_frac=0.15
            )
            _log(f"Combo ({u_key}, {s_key}): {len(candidates)} coarse candidates", t0)
            refined_list: list[dict[str, Any]] = []
            for cand in candidates:
                refined = hs.refine_candidate(
                    torus,
                    torus,
                    cand,
                    lobe_sign_u=u_lobe,
                    lobe_sign_s=s_lobe,
                    n_segments_dir=32,
                )
                if refined is None:
                    continue
                guard = hs.ghost_guard(
                    torus, torus, refined, lobe_sign_u=u_lobe, lobe_sign_s=s_lobe, n_segments_dir=32
                )
                corrected = _corrected_gaps(refined)
                corrected_off_torus = _corrected_off_torus_km(torus, refined, 0.0)
                corrected_radau = _corrected_radau_check(
                    torus,
                    torus,
                    refined,
                    lobe_sign_u=u_lobe,
                    lobe_sign_s=s_lobe,
                    theta1_section=0.0,
                    n_segments_dir=32,
                    rtol=1e-13,
                    atol=1e-13,
                )
                corrected_radau_consistent = corrected_radau["integrator_delta_km"] < 1.0
                corrected_genuine = corrected_radau_consistent and corrected_off_torus >= 1000.0
                robust_genuine = (
                    corrected_genuine
                    and corrected_radau["integrator_delta_km"] < _ROBUST_INTEGRATOR_MARGIN_KM
                )
                entry = {
                    "seed": asdict(cand),
                    "module_native_pos_gap_km": refined.pos_gap_km,
                    "module_native_vel_gap_km_s": refined.vel_gap_km_s,
                    "module_native_guard_off_torus_km": guard.off_torus_km,
                    "module_native_guard_integrator_delta_km": guard.integrator_delta_km,
                    "module_native_guard_genuine": guard.genuine,
                    "corrected_pos_gap_km": corrected["pos_gap_km"],
                    "corrected_vel_gap_km_s": corrected["vel_gap_km_s"],
                    "corrected_off_torus_km": corrected_off_torus,
                    "corrected_radau_pos_gap_km": corrected_radau["radau_pos_gap_km"],
                    "corrected_integrator_delta_km": corrected_radau["integrator_delta_km"],
                    "corrected_genuine": corrected_genuine,
                    "robust_genuine": robust_genuine,
                    "residual_norm": refined.residual_norm,
                    "converged": refined.converged,
                    "theta2_u": refined.theta2_u,
                    "t_u": refined.t_u,
                    "theta2_s": refined.theta2_s,
                    "t_s": refined.t_s,
                    "guard_quasi_jacobi_gap": guard.quasi_jacobi_gap,
                    "guard_notes": guard.notes,
                }
                refined_list.append(entry)
                _log(
                    f"  Refined: pos_gap(corrected)={corrected['pos_gap_km']:.2f} km, "
                    f"vel_gap(corrected)={corrected['vel_gap_km_s'] * 1000:.3f} m/s, "
                    f"off_torus(corrected)={corrected_off_torus:.1f} km, "
                    f"integ_delta(corrected)={corrected_radau['integrator_delta_km']:.4g} km, "
                    f"genuine(module/corrected/robust)={guard.genuine}/{corrected_genuine}/{robust_genuine}",
                    t0,
                )
                is_better = best is None or corrected["pos_gap_km"] < best["corrected_pos_gap_km"]
                if corrected_genuine and is_better:
                    best = {**entry, "u_lobe": u_lobe, "s_lobe": s_lobe}
                is_better_robust = (
                    best_robust is None
                    or corrected["pos_gap_km"] < best_robust["corrected_pos_gap_km"]
                )
                if robust_genuine and is_better_robust:
                    best_robust = {**entry, "u_lobe": u_lobe, "s_lobe": s_lobe}
            combo_results.append(
                {
                    "u_lobe": u_lobe,
                    "s_lobe": s_lobe,
                    "n_coarse_candidates": len(candidates),
                    "refined": refined_list,
                }
            )

    # Format output result JSON with Stage-3 dense mesh and Hyperion model-fidelity caveat
    out_dict: dict[str, Any] = {
        "task": "#716",
        "system": "Saturn-Titan-Hyperion CCR4BP",
        "stage3_dense_mesh_phases": 240,
        "model_fidelity_caveat": (
            "Hyperion has real orbital eccentricity e ~ 0.123 and chaotic rotation; "
            "circular-orbit CCR4BP is a lower-fidelity model for Saturn-Hyperion than Jovian moons."
        ),
        "mu": system.mu,
        "mu_gan": system.mu_gan,
        "a_gan": system.a_gan,
        "torus_residual_rms": torus.residual_rms,
        "torus_rho_strob": torus.rho_strob,
        "best_genuine_connection_corrected": best,
        "best_robust_genuine_connection_corrected": best_robust,
        "combo_results": combo_results,
    }

    out_file = OUT_DIR / "result.json"
    with open(out_file, "w") as f:
        json.dump(out_dict, f, indent=2)

    _log(f"Wrote complete hardened search summary to {out_file}", t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
