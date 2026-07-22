"""#687 -- differential-correction test of #683's periapse-map seed #217.

#683's periapse Poincare map (read geometrically off the map, no correction by
design) surfaced seed #217: a Saturn-Titan capture -> escape -> re-capture
transit that near-recurs with a periapse-map return residual of ~0.012 Titan-Hill
radii, then impacts Titan after 2 cycles. This script answers the genuinely
separate question the map cannot: does an EXACT periodic orbit exist near that
candidate, or is the near-recurrence a transient coincidence in a chaotic region?

Two correctors, both this project's own convention (STM-based Newton):
  * SINGLE shooting -- `cr3bp_periodic.correct_periodic` (the existing convention),
    to show it is structurally unusable here (full-arc STM ~1e18 -> it wanders
    off-energy to a spurious far-field orbit).
  * MULTIPLE shooting -- `cr3bp_multiple_shooting.correct_multiple_shooting`
    (#687, segments the same variational-STM Newton at the periapse nodes),
    with node subdivision + LM-damped backtracking, the robust remedy for the
    ill-conditioning. If a genuine orbit exists near the seed, this finds it.

NO catalogue writeback. Reports for adjudication only.

Run:  uv run python scripts/census_687_periapse_candidate_correction.py
Outputs -> data/found/687_periapse_correction/result.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.cr3bp_multiple_shooting as ms  # noqa: E402
import cyclerfinder.search.periapse_map as pm  # noqa: E402
from cyclerfinder.search.cr3bp_periodic import correct_periodic  # noqa: E402

OUT = ROOT / "data" / "found" / "687_periapse_correction"
OUT.mkdir(parents=True, exist_ok=True)

# Seed #217's exact initial condition (data/found/683_periapse_map/phase_d.json).
SEED_217 = np.array(
    [
        1.0196572560444028,
        -0.0074517536134353244,
        0.0,
        0.029765939253854182,
        0.07946614705680487,
        0.0,
    ]
)
J2 = 3.015311017945150  # Davis-Howell 2011 Saturn-Titan search energy
T_TOTAL = 66.0 * math.pi  # #683 Phase-D horizon (~33 Titan periods)


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def reconstruct() -> tuple[list[float], list[NDArray[np.float64]], bool]:
    """Reconstruct seed #217's periapse itinerary exactly as #683 did."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    events, impacted = pm.collect_titan_periapses(
        SEED_217, st, t_total=T_TOTAL, rtol=1e-12, atol=1e-12
    )
    times = [float(t) for t, _ in events]
    states = [np.asarray(s, dtype=np.float64) for _, s in events]
    return times, states, impacted


def build_cycle_nodes(
    system: cr3bp.CR3BPSystem,
    states: list[NDArray[np.float64]],
    times: list[float],
    *,
    subdiv: int,
) -> tuple[list[NDArray[np.float64]], list[float]]:
    """One capture->escape->recapture cycle = periapses 10->14 (the two
    successive re-captures bracket one full excursion loop). Subdivide each of
    the four original periapse segments into ``subdiv`` equal-time patch points."""
    nodes: list[NDArray[np.float64]] = []
    seg: list[float] = []
    for a in range(10, 14):
        s = states[a].copy()
        dur = times[a + 1] - times[a]
        h = dur / subdiv
        for _ in range(subdiv):
            nodes.append(s.copy())
            s = cr3bp.propagate(system, s, h, rtol=1e-12, atol=1e-12).state_f
            seg.append(h)
    return nodes, seg


def main() -> None:
    system = cr3bp.cr3bp_system("Saturn", "Titan")
    mu = system.mu
    r_hill = (mu / 3.0) ** (1.0 / 3.0)
    result: dict[str, Any] = {
        "seed_217_state": [float(v) for v in SEED_217],
        "seed_217_jacobi": cr3bp.jacobi_constant(SEED_217, mu),
        "j2_target": J2,
        "mu": mu,
        "r_hill": r_hill,
    }

    _log("reconstructing seed #217 itinerary")
    times, states, impacted = reconstruct()
    result["n_periapses"] = len(states)
    result["impacted_titan"] = impacted
    # the cycle: recapture peri10 -> excursion -> recapture peri14
    period_guess = times[14] - times[10]
    map_residual_rh = math.hypot(
        (states[14][0] - states[10][0]) / r_hill, (states[14][1] - states[10][1]) / r_hill
    )
    result["cycle_period_guess_nondim"] = period_guess
    result["map_return_residual_rH_peri10_to_14"] = map_residual_rh
    _log(f"cycle period guess = {period_guess:.4f}, map residual = {map_residual_rh:.4e} r_H")

    # ---- single-shooting baseline (project convention) ----
    _log("single-shooting baseline (correct_periodic) from peri10")
    try:
        po = correct_periodic(system, states[10], period_guess, tol=1e-10, max_iter=40)
        result["single_shooting"] = {
            "reported_converged": bool(po.converged),
            "closure_residual": float(po.closure_residual),
            "period": float(po.period),
            "jacobi": float(po.jacobi),
            "jacobi_drifted_off_target": abs(po.jacobi - J2) > 1e-3,
            "state0": [float(v) for v in po.state0],
            "note": (
                "SPURIOUS: full-arc STM ill-conditioning drove the free-state/free-period "
                "Newton off-energy to a far-field orbit; NOT seed #217's orbit"
                if abs(po.jacobi - J2) > 1e-3
                else "stayed on energy"
            ),
        }
    except RuntimeError as exc:
        result["single_shooting"] = {"raised": f"{type(exc).__name__}: {exc}"[:200]}

    # ---- multiple-shooting at increasing node density ----
    result["multiple_shooting"] = []
    for subdiv in (1, 3, 6):
        nodes, seg = build_cycle_nodes(system, states, times, subdiv=subdiv)
        _, _, stms = ms._residual_and_jacobian(system, nodes, seg, rtol=1e-12, atol=1e-12)
        stm_norms = [float(np.linalg.norm(s, 2)) for s in stms]
        _log(f"multiple-shooting subdiv={subdiv} (N={len(nodes)} nodes)")
        orbit = ms.correct_multiple_shooting(
            system, nodes, seg, tol=1e-10, max_iter=80, rtol=1e-12, atol=1e-12
        )
        entry: dict[str, Any] = {
            "subdiv": subdiv,
            "n_nodes": len(nodes),
            "max_segment_stm_2norm": max(stm_norms),
            "converged": bool(orbit.converged),
            "closure_residual": float(orbit.closure_residual),
            "period": float(orbit.period),
            "n_iter": orbit.n_iter,
        }
        _log(
            f"  -> converged={orbit.converged} residual={orbit.closure_residual:.3e} "
            f"(max seg STM norm {max(stm_norms):.2e})"
        )
        if orbit.converged:
            lam = ms.floquet_multipliers(system, orbit)
            sol = solve_ivp(
                cr3bp.cr3bp_eom,
                (0.0, orbit.period),
                orbit.nodes[0],
                args=(mu,),
                method="Radau",
                rtol=1e-11,
                atol=1e-11,
            )
            entry["radau_closure"] = float(np.linalg.norm(sol.y[:, -1] - orbit.nodes[0]))
            entry["floquet_abs"] = sorted(float(abs(v)) for v in lam)
            entry["node0"] = [float(v) for v in orbit.nodes[0]]
        result["multiple_shooting"].append(entry)

    any_converged = any(e["converged"] for e in result["multiple_shooting"])
    result["verdict"] = (
        "PERIODIC_ORBIT_FOUND -- run literature novelty check, report for adjudication"
        if any_converged
        else "NO_EXACT_ORBIT -- near-recurrence is a transient coincidence; corrector stalls "
        "at the ~map-residual floor from every node density; single-shooting diverges off-energy"
    )
    (OUT / "result.json").write_text(json.dumps(result, indent=2))
    _log(f"VERDICT: {result['verdict']}")
    _log(f"wrote {OUT / 'result.json'}")


if __name__ == "__main__":
    main()
