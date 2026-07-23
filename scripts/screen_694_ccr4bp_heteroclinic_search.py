"""#694 -- CCR4BP manifold globalization + mesh-intersection heteroclinic
search: the JEG (Jupiter-Europa-Ganymede) positive-control run.

Stage B's final sub-task per #686's plan. Runs the full pipeline built in
this task (search/ccr4bp_manifold_globalize.py + search/ccr4bp_heteroclinic_
search.py) against the ALREADY fully-validated Jupiter-Europa 3:4 CCR4BP
torus (#690/#691's own object) before this machinery is trusted on any novel
system (Io-Europa, Io-Ganymede -- explicitly NOT run here, per #694's scope).

Positive control scope: per #694's dispatch, the SINGLE-TORUS
(homoclinic-style: the torus's own unstable manifold meeting its own stable
manifold) search is used as the first, simpler positive control, with the
genuine two-torus (Ganymede-4:3 <-> Europa-3:4) heteroclinic cycle #686's own
shortlist item ultimately wants left as a documented future extension (see
the bottom of this file's output for what that would additionally require).

Both manifold LOBES (lobe_sign = +1 and -1, the two sides of each branch's
extracted eigendirection) are searched, giving 4 lobe-pair combinations. For
each, the coarse KD-tree search + continuous least-squares refinement +
mandatory ghost-minima guard (independent Radau integrator, quasi-Jacobi
consistency, exact off-torus divergence) run per
search/ccr4bp_heteroclinic_search.py. A denser re-globalization is run at the
end as the mesh-refinement stability check on the best candidate found.

No catalogue writeback -- capability-proof + discovery-attempt, not a vetted
result. Screen-grade naming (not "run_*.py") because this is an exploratory
capability run in the #688/#693 sense, not a preflight_search-gated
registry sweep.

Run:  uv run python scripts/screen_694_ccr4bp_heteroclinic_search.py
Outputs -> data/found/694_ccr4bp_heteroclinic_search/result.json
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

import cyclerfinder.core.ccr4bp as ccr4bp  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "694_ccr4bp_heteroclinic_search"


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691's own test modules
    (base-CR3BP resonant-orbit seed; no production code under test)."""
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


def main() -> dict[str, Any]:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    system = ccr4bp.jupiter_europa_ganymede_default()
    s0, period, res = _resonant_symmetric_orbit(system.mu, 3, 4)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    _log(f"base 3:4 resonant orbit converged, perp residual {res:.2e}", t0)

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
        f"physical-mass 3:4 CCR4BP torus built: residual_rms={torus.residual_rms:.3e}, "
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
                entry = {
                    "seed": asdict(cand),
                    "refined_pos_gap_km": refined.pos_gap_km,
                    "refined_vel_gap_km_s": refined.vel_gap_km_s,
                    "residual_norm": refined.residual_norm,
                    "converged": refined.converged,
                    "theta2_u": refined.theta2_u,
                    "t_u": refined.t_u,
                    "theta2_s": refined.theta2_s,
                    "t_s": refined.t_s,
                    "guard_radau_pos_gap_km": guard.radau_pos_gap_km,
                    "guard_integrator_delta_km": guard.integrator_delta_km,
                    "guard_quasi_jacobi_gap": guard.quasi_jacobi_gap,
                    "guard_off_torus_km": guard.off_torus_km,
                    "guard_genuine": guard.genuine,
                    "guard_notes": guard.notes,
                }
                refined_list.append(entry)
                _log(
                    f"  refined: pos_gap={refined.pos_gap_km:.2f} km, "
                    f"vel_gap={refined.vel_gap_km_s * 1000:.3f} m/s, "
                    f"genuine={guard.genuine}",
                    t0,
                )
                is_better = best is None or refined.pos_gap_km < best["refined_pos_gap_km"]
                if guard.genuine and is_better:
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
        _log(f"best genuine connection: pos_gap={best['refined_pos_gap_km']:.3f} km", t0)
    else:
        _log("NO genuine connection found across any lobe combo", t0)

    # ------------------------------------------------------------------
    # Stage 3: mesh-refinement stability check on the best candidate.
    # ------------------------------------------------------------------
    mesh_check: dict[str, Any] | None = None
    if best is not None:
        n_theta2_dense, n_time_dense = 120, 300
        u_key = f"unstable_{'+' if best['u_lobe'] > 0 else '-'}"
        s_key = f"stable_{'+' if best['s_lobe'] > 0 else '-'}"
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
        # Refine the SAME seed the sparse grid found (independent of tube
        # resolution -- manifold_state_at is a continuous primitive) as the
        # primary stability check, plus report the dense-grid coarse floor.
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
        mesh_check = {
            "n_theta2_dense": n_theta2_dense,
            "n_time_dense": n_time_dense,
            "dense_grid_coarse_best_gap_planar": dense_best_gap,
            "rerefine_at_n_segments_dir_48_pos_gap_km": (reref.pos_gap_km if reref else None),
            "reref_matches_original": (
                abs(reref.pos_gap_km - best["refined_pos_gap_km"]) < 5.0 if reref else False
            ),
        }
        reref_gap = mesh_check["rerefine_at_n_segments_dir_48_pos_gap_km"]
        _log(f"mesh-refinement check: dense-grid re-refine pos_gap={reref_gap}", t0)

    result = {
        "task": "#694",
        "system": "Jupiter-Europa-Ganymede CCR4BP, Europa 3:4 resonant torus (physical mass)",
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
        "best_genuine_connection": best,
        "mesh_refinement_check": mesh_check,
        "two_torus_extension_note": (
            "This positive control is the SINGLE-TORUS homoclinic case "
            "(Europa 3:4 unstable manifold meeting its OWN stable manifold), "
            "per #694's own dispatch (acceptable simpler first control). A "
            "genuine two-torus heteroclinic cycle (#686's own shortlist item: "
            "Ganymede-4:3 <-> Europa-3:4) would additionally require: (1) "
            "building a second CCR4BP torus for the Jupiter-Ganymede 4:3 "
            "resonant orbit under EUROPA's forcing (a different base system, "
            "mu/mu_gan roles swapped -- #689's CCR4BPSystem is already "
            "system-agnostic so this is a parameter change, not new code); "
            "(2) globalizing that torus's manifolds with this SAME module; "
            "(3) running coarse_candidates/refine_candidate/ghost_guard with "
            "torus_u != torus_s (already supported by this module's own "
            "signature -- refine_candidate and ghost_guard both take "
            "independent torus_u/torus_s arguments); (4) closing the CYCLE "
            "itself (Ganymede torus -> Europa torus -> back to Ganymede "
            "torus) requires phasing commensurability across BOTH legs, an "
            "additional constraint beyond a single connection. No code "
            "changes needed beyond calling this module twice with different "
            "torus objects; not run here per #694's explicit JEG-only scope."
        ),
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result


if __name__ == "__main__":
    main()
