"""#696 -- CCR4BP manifold globalization + heteroclinic search: Jupiter-Io-Ganymede.

Genuine discovery attempt (not a capability-proof exercise -- that was #694's
JEG positive control, `scripts/screen_694_ccr4bp_heteroclinic_search.py`,
whose structure this script mirrors). #693's second-ranked novelty-cleared
CCR4BP candidate: Io-Ganymede, `mu_gan=7.80e-5` essentially IDENTICAL forcing
strength to the already-validated JEG system's own perturber term, 4.06
period ratio (a real 4:2:1 Laplace-chain consequence), e/Delta-i both under
#693's tractability bar, not found in 2 targeted literature checks.

Every module used (core.ccr4bp, search.variational_ccr4bp_torus,
search.ccr4bp_whisker, search.ccr4bp_manifold_globalize,
search.ccr4bp_heteroclinic_search) is #689-#694's code, reused UNMODIFIED.
The only new code is the system constructor (core.ccr4bp_io_ganymede,
#696) and this driver.

Base resonant orbit choice (documented in full in
tests/search/test_ccr4bp_torus_io_ganymede.py's module docstring): NOT the
naive "moon-ratio-matching" spacecraft:Io = 1:4 exterior orbit (a~2.5198,
which sits within 0.018 Io-SMA units of Ganymede's own orbital radius
a_gan~2.5377 -- propagating it under physical Ganymede forcing throws a
genuine near-collision RuntimeError, confirmed and regression-tested). This
script uses the OTHER natural "4:1"-labeled orbit, spacecraft:Io = 4:1
(interior, a~0.3968), which clears Ganymede by ~2.14 Io-SMA units and whose
CCR4BP torus converges CLEANER than JEG's own Europa 3:4 torus.

IMPORTANT physical-unit caveat (a genuine #694-pipeline generality finding,
not fixed here per this task's "reuse unmodified" scope): search.
ccr4bp_heteroclinic_search hardcodes a module-level `_L_KM = 671_100.0` (=
EUROPA's SMA) and `_v_unit_km_s` hardcodes GM_Europa -- both baked in for the
JEG positive control, NOT parameterized per-system. For Io-Ganymede (base
moon = Io, SMA 421,800 km) every km-denominated field that module returns
(`RefinedConnection.pos_gap_km`/`vel_gap_km_s`, `GhostGuardReport.
radau_pos_gap_km`/`off_torus_km`) is scaled by the WRONG physical unit --
about 0.628x too large (Europa's SMA / Io's SMA). This is large enough to
matter for the `off_torus_min_km=1000` ghost-guard gate, so this script
independently recomputes every physical-unit quantity using Io's OWN SMA/GM
(correct nondimensional-to-km conversion) from the raw nondimensional states
each RefinedConnection/GhostGuardReport already carries, and reports BOTH the
module-native (Europa-scaled) and corrected (Io-scaled) verdicts side by
side -- flagging any case where they disagree on "genuine".

No catalogue writeback -- a discovery-attempt script, not a vetted result.

Run:  uv run python scripts/screen_696_ccr4bp_io_ganymede_search.py
Outputs -> data/found/696_ccr4bp_io_ganymede_search/result.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.core.ccr4bp_io_ganymede import jupiter_io_ganymede_default  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402
from cyclerfinder.search.variational_ccr4bp_torus import evaluate_torus_state  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "696_ccr4bp_io_ganymede_search"

# Correct physical units for THIS system (base moon = Io, not Europa).
_L_KM_IO = SATELLITES["Io"].sma_km  # 421,800 km
_GM_J = PRIMARIES["Jupiter"]
_GM_IO = SATELLITES["Io"].mu_km3_s2
_N2_IO = math.sqrt((_GM_J + _GM_IO) / _L_KM_IO**3)  # Io's own mean motion, rad/s
_V_UNIT_IO_KM_S = _L_KM_IO * _N2_IO


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694's own scripts/tests."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
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


def _corrected_gaps(
    torus_u: vt.CCR4BPTorusVariationalResult,
    torus_s: vt.CCR4BPTorusVariationalResult,
    refined: hs.RefinedConnection,
) -> dict[str, float]:
    """Recompute pos/vel gap in km using Io's OWN physical units (not the
    module's hardcoded Europa constants) directly from the raw nondimensional
    states the module already returns -- no #694 module code touched."""
    pos_gap_km = float(np.linalg.norm(refined.state_u[:3] - refined.state_s[:3])) * _L_KM_IO
    vel_gap_km_s = (
        float(np.linalg.norm(refined.state_u[3:] - refined.state_s[3:])) * _V_UNIT_IO_KM_S
    )
    return {"pos_gap_km": pos_gap_km, "vel_gap_km_s": vel_gap_km_s}


def _corrected_off_torus_km(
    torus_u: vt.CCR4BPTorusVariationalResult,
    refined: hs.RefinedConnection,
    theta1_section: float,
) -> float:
    """Replicate ghost_guard's off-torus distance formula with Io's own km
    scale (see ghost_guard's own docstring for the formula this mirrors)."""
    theta1_u_at_t = theta1_section + torus_u.omega1 * refined.t_u
    theta2_u_at_t = refined.theta2_u + torus_u.omega2 * refined.t_u
    torus_pt_u = evaluate_torus_state(torus_u, theta1_u_at_t, theta2_u_at_t)
    planar_u = refined.state_u[[0, 1, 3, 4]]
    return float(np.linalg.norm(torus_pt_u - planar_u)) * _L_KM_IO


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
    """Re-run #694's own private Radau re-propagation helper (the same
    independent-integrator check ghost_guard performs internally) and rescale
    the resulting position gap with Io's own km unit instead of the module's
    hardcoded Europa constant. Calls hs._radau_manifold_state directly --
    the same private-helper access pattern this project's OWN test module
    (tests/search/test_ccr4bp_heteroclinic_search.py) already uses for
    hs._L_KM -- not a modification of #694's module."""
    ref_vec_u = hs._direction_only(
        torus_u, "unstable", theta1_section, refined.theta2_u, n_segments_dir, rtol, atol
    )
    ref_vec_s = hs._direction_only(
        torus_s, "stable", theta1_section, refined.theta2_s, n_segments_dir, rtol, atol
    )
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
    radau_pos_gap_km = float(np.linalg.norm(su_radau[:3] - ss_radau[:3])) * _L_KM_IO
    corrected_dop853_pos_gap_km = (
        float(np.linalg.norm(refined.state_u[:3] - refined.state_s[:3])) * _L_KM_IO
    )
    integrator_delta_km = abs(radau_pos_gap_km - corrected_dop853_pos_gap_km)
    return {"radau_pos_gap_km": radau_pos_gap_km, "integrator_delta_km": integrator_delta_km}


def main() -> dict[str, Any]:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    system = jupiter_io_ganymede_default()
    _log(
        f"system: mu={system.mu:.6e} mu_gan={system.mu_gan:.6e} a_gan={system.a_gan:.4f} "
        f"omega_gan={system.omega_gan:.4f}",
        t0,
    )
    s0, period, res = _resonant_symmetric_orbit(system.mu, 4, 1)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    _log(f"base 4:1 (interior) resonant orbit converged, perp residual {res:.2e}", t0)

    torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    _log(
        f"physical-mass 4:1 CCR4BP torus built: residual_rms={torus.residual_rms:.3e}, "
        f"closure_residual={torus.closure_residual:.3e}, rho_strob={torus.rho_strob:.4f}, "
        f"period={torus.period:.4f} TU",
        t0,
    )

    # ------------------------------------------------------------------
    # Stage 1: globalize both branches, both lobes.
    # ------------------------------------------------------------------
    n_theta2, n_time, t_max_periods, n_segments_dir = 60, 150, 2.0, 24
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
            _log(f"globalized {key}: {n_valid}/{n_theta2} phases valid", t0)

    # ------------------------------------------------------------------
    # Stage 2: coarse search + refine + ghost-guard, all 4 lobe-pair combos.
    # ------------------------------------------------------------------
    combo_results: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    for u_lobe in (1.0, -1.0):
        for s_lobe in (1.0, -1.0):
            u_key = f"unstable_{'+' if u_lobe > 0 else '-'}"
            s_key = f"stable_{'+' if s_lobe > 0 else '-'}"
            candidates = hs.coarse_candidates(
                tubes[u_key], tubes[s_key], n_candidates=5, t_min_frac=0.15
            )
            _log(f"combo ({u_key}, {s_key}): {len(candidates)} coarse candidates", t0)
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
                corrected = _corrected_gaps(torus, torus, refined)
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
                entry = {
                    "seed": asdict(cand),
                    # Module-native fields (Europa-SMA-scaled -- see module docstring caveat).
                    "module_native_pos_gap_km": refined.pos_gap_km,
                    "module_native_vel_gap_km_s": refined.vel_gap_km_s,
                    "module_native_guard_off_torus_km": guard.off_torus_km,
                    "module_native_guard_integrator_delta_km": guard.integrator_delta_km,
                    "module_native_guard_genuine": guard.genuine,
                    # Independently corrected (Io-SMA-scaled) physical values.
                    "corrected_pos_gap_km": corrected["pos_gap_km"],
                    "corrected_vel_gap_km_s": corrected["vel_gap_km_s"],
                    "corrected_off_torus_km": corrected_off_torus,
                    "corrected_radau_pos_gap_km": corrected_radau["radau_pos_gap_km"],
                    "corrected_integrator_delta_km": corrected_radau["integrator_delta_km"],
                    "corrected_genuine": corrected_genuine,
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
                    f"  refined: pos_gap(corrected)={corrected['pos_gap_km']:.2f} km, "
                    f"vel_gap(corrected)={corrected['vel_gap_km_s'] * 1000:.3f} m/s, "
                    f"off_torus(corrected)={corrected_off_torus:.1f} km, "
                    f"genuine(module/corrected)={guard.genuine}/{corrected_genuine}",
                    t0,
                )
                is_better = best is None or corrected["pos_gap_km"] < best["corrected_pos_gap_km"]
                if corrected_genuine and is_better:
                    best = {**entry, "u_lobe": u_lobe, "s_lobe": s_lobe}
            combo_results.append(
                {
                    "u_lobe": u_lobe,
                    "s_lobe": s_lobe,
                    "n_coarse_candidates": len(candidates),
                    "refined": refined_list,
                }
            )

    if best is not None:
        _log(
            f"best genuine (corrected) connection: pos_gap={best['corrected_pos_gap_km']:.3f} km",
            t0,
        )
    else:
        _log("NO genuine connection found across any lobe combo (corrected verdict)", t0)

    # ------------------------------------------------------------------
    # Stage 3: mesh-refinement stability check on the best candidate.
    # ------------------------------------------------------------------
    mesh_check: dict[str, Any] | None = None
    if best is not None:
        n_theta2_dense, n_time_dense = 120, 300
        tube_u_dense = mg.globalize_manifold_tube(
            torus,
            "unstable",
            n_theta2=n_theta2_dense,
            t_max_periods=t_max_periods,
            n_time=n_time_dense,
            n_segments_dir=n_segments_dir,
            lobe_sign=best["u_lobe"],
        )
        tube_s_dense = mg.globalize_manifold_tube(
            torus,
            "stable",
            n_theta2=n_theta2_dense,
            t_max_periods=t_max_periods,
            n_time=n_time_dense,
            n_segments_dir=n_segments_dir,
            lobe_sign=best["s_lobe"],
        )
        dense_candidates = hs.coarse_candidates(
            tube_u_dense, tube_s_dense, n_candidates=5, t_min_frac=0.15
        )
        dense_best_gap = min((c.gap_planar for c in dense_candidates), default=float("inf"))
        reref = hs.refine_candidate(
            torus,
            torus,
            hs.ManifoldCandidate(
                theta2_u=best["theta2_u"],
                t_u=best["t_u"],
                theta2_s=best["theta2_s"],
                t_s=best["t_s"],
                gap_planar=0.0,
            ),
            lobe_sign_u=best["u_lobe"],
            lobe_sign_s=best["s_lobe"],
            n_segments_dir=48,
        )
        reref_corrected_pos_gap_km = (
            float(np.linalg.norm(reref.state_u[:3] - reref.state_s[:3])) * _L_KM_IO
            if reref is not None
            else None
        )
        mesh_check = {
            "n_theta2_dense": n_theta2_dense,
            "n_time_dense": n_time_dense,
            "dense_grid_coarse_best_gap_planar": dense_best_gap,
            "rerefine_at_n_segments_dir_48_pos_gap_km_corrected": reref_corrected_pos_gap_km,
            "reref_matches_original": (
                abs(reref_corrected_pos_gap_km - best["corrected_pos_gap_km"]) < 5.0
                if reref_corrected_pos_gap_km is not None
                else False
            ),
        }
        _log(
            f"mesh-refinement check: dense-grid re-refine pos_gap(corrected)="
            f"{reref_corrected_pos_gap_km}",
            t0,
        )

    result = {
        "task": "#696",
        "system": "Jupiter-Io-Ganymede CCR4BP, Io 4:1 (interior) resonant torus (physical mass)",
        "physical_unit_caveat": (
            "search.ccr4bp_heteroclinic_search hardcodes Europa's SMA/GM (_L_KM=671100 km) for "
            "all km conversions -- a real, unfixed (per this task's reuse-unmodified scope) "
            "generality gap for a non-JEG base moon. module_native_* fields use that (WRONG, "
            "~1.59x too large) scale; corrected_* fields independently recompute the same "
            "quantities using Io's own SMA (421800 km) / GM. The 'genuine' verdict used for "
            "'best' selection is the CORRECTED one."
        ),
        "l_km_io": _L_KM_IO,
        "v_unit_io_km_s": _V_UNIT_IO_KM_S,
        "torus_residual_rms": torus.residual_rms,
        "torus_closure_residual": torus.closure_residual,
        "torus_rho_strob": torus.rho_strob,
        "torus_period_tu": torus.period,
        "globalization_params": {
            "n_theta2": n_theta2,
            "n_time": n_time,
            "t_max_periods": t_max_periods,
            "n_segments_dir": n_segments_dir,
            "eps": mg.DEFAULT_EPS,
        },
        "combo_results": combo_results,
        "best_genuine_connection_corrected": best,
        "mesh_refinement_check": mesh_check,
        "base_orbit_note": (
            "Base resonant orbit is spacecraft:Io=4:1 (interior, a~0.3968), NOT the "
            "moon-ratio-matching spacecraft:Io=1:4 exterior orbit (a~2.5198), which sits "
            "within 0.018 Io-SMA units of Ganymede's own orbital radius and throws a genuine "
            "near-collision RuntimeError under physical Ganymede forcing -- see "
            "tests/search/test_ccr4bp_torus_io_ganymede.py::"
            "test_1to4_resonant_orbit_collides_with_ganymede."
        ),
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result


if __name__ == "__main__":
    main()
