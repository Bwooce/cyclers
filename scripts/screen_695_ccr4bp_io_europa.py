"""#695 -- CCR4BP manifold globalization + mesh-intersection heteroclinic
search: Jupiter Io-Europa, a GENUINE discovery attempt (not a positive
control -- that was `#694`'s own JEG run).

`#693`'s cross-solar-system CCR4BP moon-pair screening pass ranked Jupiter
Io-Europa as the TOP novelty-cleared candidate for the now-proven
`#689`-`#694` pipeline. This script runs that full pipeline -- system
construction (`#695`'s own ``core.ccr4bp_io_europa``), base resonant orbit,
`#690`'s torus corrector, `#691`'s segmented-CLV whisker extraction, `#694`'s
manifold globalization + intersection search + mandatory ghost-guard -- on
Io-Europa, reusing every `#689`-`#694` module UNMODIFIED.

Base resonant orbit: the spacecraft:Io = 1:2 symmetric resonance (period
ratio 2:1, exterior reading; the interior 2:1 reading does not converge --
see ``tests/search/test_ccr4bp_torus_io_europa.py``'s own documented
negative).

Physical-unit correction (this task's own finding)
----------------------------------------------------
`#694`'s ``search/ccr4bp_heteroclinic_search.py`` hardcodes its physical-unit
conversion (``_L_KM = 671_100.0``, Europa's own SMA, and ``_v_unit_km_s``'s
internal GM lookup) to the JEG system specifically -- despite every other
function in that module being genuinely system-agnostic (``torus_u``/
``torus_s`` are free parameters throughout). Applied UNMODIFIED to Io-Europa
(whose correct length unit is Io's SMA, 421,800 km -- a factor
421800/671100 = 0.6285 smaller), every "_km"-suffixed field
``refine_candidate``/``ghost_guard`` return (``pos_gap_km``,
``vel_gap_km_s``, ``radau_pos_gap_km``, ``integrator_delta_km``,
``off_torus_km``) is computed in the WRONG physical scale for this system,
and the ghost-guard's own absolute-km GATES (``integrator_consistency_km=1.0
km``, ``off_torus_min_km=1000.0 km``) are then compared against those
wrongly-scaled numbers. This is NOT a numerical bug in `#694`'s math -- the
underlying nondimensional residual (``residual_norm``, from
``refine_candidate``) is unaffected, since it never touches ``_L_KM`` -- it
is a hardcoded PHYSICAL CONSTANT baked into two helper functions that this
task's own dispatch requires be reused unmodified. This script therefore
calls `#694`'s functions exactly as published (getting their raw, JEG-scaled
"_km" fields for reproducibility/transparency), then INDEPENDENTLY
back-corrects every physical-unit field using this task's own correct
constants (``core.ccr4bp_io_europa.L_KM`` / ``.v_unit_km_s()``) via pure
post-hoc arithmetic (no private-function calls, no module edits): since
every wrong "_km" field is linear in the wrong constant
(``value_wrong_km = nondim_value * wrong_constant``), the correct value is
recovered by ``value_correct_km = value_wrong_km * (correct_constant /
wrong_constant)``. The corrected fields (prefixed ``corrected_``) are what
this script's own genuine/non-genuine verdict is based on -- NOT the raw
`#694`-native fields (which are kept in the output purely for audit/
reproducibility against `#694`'s own JEG result.json format).

Run:  uv run python scripts/screen_695_ccr4bp_io_europa.py
Outputs -> data/found/695_ccr4bp_io_europa_search/result.json
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

import cyclerfinder.core.ccr4bp_io_europa as ioeu  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "695_ccr4bp_io_europa_search"

# #694's own hardcoded (Europa-SMA-based) physical-unit constants -- needed
# ONLY to back out the nondimensional value from its "_km"-suffixed outputs
# (see module docstring). Duplicated here as bare constants (not imported
# privates) so this script does not depend on that module's internals.
_JEG_L_KM = 671_100.0
_JEG_GM_JUPITER = 1.26686534e8
_JEG_GM_EUROPA = 3202.739


def _jeg_v_unit_km_s() -> float:
    import math

    n2 = math.sqrt((_JEG_GM_JUPITER + _JEG_GM_EUROPA) / _JEG_L_KM**3)
    return _JEG_L_KM * n2


_LENGTH_SCALE = ioeu.L_KM / _JEG_L_KM
_VELOCITY_SCALE = ioeu.v_unit_km_s() / _jeg_v_unit_km_s()


def _corrected_km(wrong_km: float) -> float:
    return wrong_km * _LENGTH_SCALE


def _corrected_vel(wrong_km_s: float) -> float:
    return wrong_km_s * _VELOCITY_SCALE


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694's own test/driver
    copies (no production code under test)."""
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

    system = ioeu.jupiter_io_europa_default()
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    _log(f"base spacecraft:Io=1:2 resonant orbit converged, perp residual {res:.2e}", t0)

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
        f"physical-mass 1:2 CCR4BP torus built: residual_rms={torus.residual_rms:.3e}, "
        f"closure_residual={torus.closure_residual:.3e}, rho_strob={torus.rho_strob:.4f}, "
        f"period={torus.period:.4f} TU",
        t0,
    )

    # ------------------------------------------------------------------
    # Stage 1: globalize both branches, both lobes.
    #
    # t_max_periods=2.0 matches #694's own JEG default. This task's own
    # exploratory probe (scratch analysis, not committed as production code)
    # found off-torus divergence already reaches ~0.2-0.9 nondim units by
    # t=2 periods at this torus's strongest phase band (a narrow,
    # SHARPLY unstable region near theta2~2*pi with |lam_u| up to ~425 --
    # much less UNIFORMLY unstable than JEG's torus, which had |lam_u|
    # 9-1700 fairly evenly across phase, but comparably strong at its own
    # peak), so the same horizon used for JEG is a reasonable and
    # comparable starting point here too, not an under-resourced guess.
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
                tubes[u_key], tubes[s_key], n_candidates=8, t_min_frac=0.15
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
                corrected_off_torus_km = _corrected_km(guard.off_torus_km)
                corrected_integrator_delta_km = _corrected_km(guard.integrator_delta_km)
                corrected_genuine = (
                    corrected_integrator_delta_km < 1.0 and corrected_off_torus_km >= 1000.0
                )
                entry: dict[str, Any] = {
                    "seed": asdict(cand),
                    "refined_pos_gap_km_JEGSCALE_RAW": refined.pos_gap_km,
                    "refined_vel_gap_km_s_JEGSCALE_RAW": refined.vel_gap_km_s,
                    "corrected_pos_gap_km": _corrected_km(refined.pos_gap_km),
                    "corrected_vel_gap_km_s": _corrected_vel(refined.vel_gap_km_s),
                    "residual_norm": refined.residual_norm,  # nondim, unaffected by unit bug
                    "converged": refined.converged,
                    "theta2_u": refined.theta2_u,
                    "t_u": refined.t_u,
                    "theta2_s": refined.theta2_s,
                    "t_s": refined.t_s,
                    "guard_radau_pos_gap_km_JEGSCALE_RAW": guard.radau_pos_gap_km,
                    "guard_integrator_delta_km_JEGSCALE_RAW": guard.integrator_delta_km,
                    "corrected_radau_pos_gap_km": _corrected_km(guard.radau_pos_gap_km),
                    "corrected_integrator_delta_km": corrected_integrator_delta_km,
                    "guard_quasi_jacobi_gap": guard.quasi_jacobi_gap,  # nondim
                    "guard_off_torus_km_JEGSCALE_RAW": guard.off_torus_km,
                    "corrected_off_torus_km": corrected_off_torus_km,
                    "guard_genuine_JEGSCALE_RAW": guard.genuine,
                    "corrected_genuine": corrected_genuine,
                    "guard_notes": guard.notes,
                }
                refined_list.append(entry)
                _log(
                    f"  refined: corrected_pos_gap={entry['corrected_pos_gap_km']:.2f} km, "
                    f"corrected_vel_gap={entry['corrected_vel_gap_km_s'] * 1000:.3f} m/s, "
                    f"corrected_genuine={corrected_genuine}",
                    t0,
                )
                is_better = (
                    best is None or entry["corrected_pos_gap_km"] < best["corrected_pos_gap_km"]
                )
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
            f"best CORRECTED genuine connection: pos_gap={best['corrected_pos_gap_km']:.3f} km",
            t0,
        )
    else:
        _log("NO genuine connection found across any lobe combo (corrected units)", t0)

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
        reref_corrected_km = _corrected_km(reref.pos_gap_km) if reref else None
        mesh_check = {
            "n_theta2_dense": n_theta2_dense,
            "n_time_dense": n_time_dense,
            "dense_grid_coarse_best_gap_planar": dense_best_gap,
            "rerefine_at_n_segments_dir_48_corrected_pos_gap_km": reref_corrected_km,
            "reref_matches_original": (
                abs(reref_corrected_km - best["corrected_pos_gap_km"]) < 5.0
                if reref_corrected_km is not None
                else False
            ),
        }
        _log(f"mesh-refinement check: dense-grid re-refine pos_gap={reref_corrected_km}", t0)

    result = {
        "task": "#695",
        "system": "Jupiter-Io-Europa CCR4BP, Io spacecraft:Io=1:2 resonant torus (physical mass)",
        "unit_correction_note": (
            "#694's search/ccr4bp_heteroclinic_search.py hardcodes its km-unit "
            "conversion to the JEG system (Europa SMA = 671100 km); fields "
            "suffixed _JEGSCALE_RAW are #694's own unmodified output (wrong "
            "physical scale for THIS system, kept for audit); fields prefixed "
            "corrected_ are this script's own post-hoc rescaling to the "
            "correct Io-SMA-based unit (core.ccr4bp_io_europa.L_KM = "
            f"{ioeu.L_KM} km). Length scale factor applied: {_LENGTH_SCALE:.6f}, "
            f"velocity scale factor: {_VELOCITY_SCALE:.6f}. The nondimensional "
            "residual_norm and quasi_jacobi_gap fields are UNAFFECTED by this "
            "bug (never touch the hardcoded km constant) and are reported "
            "as-is from #694's own code."
        ),
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
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result


if __name__ == "__main__":
    main()
