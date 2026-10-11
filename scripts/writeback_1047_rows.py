"""#1047: catalogue rows for the asymmetric Earth-Moon families A (5:2 prograde), B (5:3
retrograde) and D (7:3 prograde) of #1000 (owner decision 2026-10-11, "Rows for A, B, D only"):
our_status known-class-member (#875 (iii)), validation V1.

For each family the representative is the #1000 grid member that seeded its continuation (line 1
of data/1000_complement/continuation/<family>_p.jsonl). Its perigee state is re-closed with the
project's full-state single-shooting corrector (cr3bp_periodic.correct_periodic, the corrector of
the Casoliva Table 3 rows), and the V1 evidence follows the Casoliva rows' convention:
- closure residual of the corrected orbit;
- an independent Radau re-propagation: the Jacobi drift over one period (crosscheck_periodic), and
  the arc-by-arc agreement with DOP853 over 4p equal arcs (each Radau arc starts from the DOP853
  state); the whole-period Radau closure is recorded too, but on an orbit with |lambda| ~ 1e2-1e3 it
  measures the instability's amplification of integrator error, not the orbit (the reason #1000
  uses multiple shooting);
- the stability index from the monodromy, both from one full-period STM and from the product of
  p per-arc STMs (the #1000 note: single shooting over a lambda ~ 1e2-1e3 orbit can drift);
- a dense DOP853 re-propagation (run_1000_complement.classify): perigee, periselene, maximum Moon
  distance, windings, symmetry (perpendicular x-axis crossings);
- inertial p:q and sense from N = w_E + T / (2 pi);
- the mirror image (x, -y, -vx, vy) is a distinct orbit (asymmetric).
Writes data/1047_rows/writeback_v1.json. The catalogue edit is made in the catalogue window.

Usage: uv run python scripts/writeback_1047_rows.py
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core import cr3bp
from cyclerfinder.search import cr3bp_periodic as cp
from cyclerfinder.search.earth_moon_resonant_families import (
    earth_moon_system,
    planar_stability_index,
)

REPO = Path(__file__).resolve().parents[1]
CONT = REPO / "data" / "1000_complement" / "continuation"
OUT = REPO / "data" / "1047_rows" / "writeback_v1.json"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


sys.path.insert(0, str(REPO / "scripts"))
BASE = _load("run_1000_complement", REPO / "scripts" / "run_1000_complement.py")

#: family -> (row id, continuation file stem, p, q)
FAMILIES: dict[str, tuple[str, str, int, int]] = {
    "A": ("earth-moon-5-2-prograde-asymmetric-cycler-1000a-2026", "5_2_pro_wE3_wM0", 5, 2),
    "B": ("earth-moon-5-3-retrograde-asymmetric-cycler-1000b-2026", "5_3_ret_wE-8_wM-2", 5, 3),
    "D": ("earth-moon-7-3-prograde-asymmetric-cycler-1000d-2026", "7_3_pro_wE4_wM-1", 7, 3),
}


def arc_stm_product(system: cr3bp.CR3BPSystem, s6: np.ndarray, period: float, n: int) -> float:
    """Planar k = trace(M4) - 2 from the product of n equal-duration per-arc STMs."""
    m = np.eye(6)
    s = s6.copy()
    for _ in range(n):
        arc = cr3bp.propagate(system, s, period / n, with_stm=True, rtol=1e-12, atol=1e-12)
        assert arc.stm is not None
        m = np.asarray(arc.stm) @ m
        s = np.asarray(arc.state_f)
    m4 = m[np.ix_([0, 1, 3, 4], [0, 1, 3, 4])]
    return float(np.trace(m4) - 2.0)


def radau_arcs(system: cr3bp.CR3BPSystem, s6: np.ndarray, period: float, n: int) -> float:
    """Largest end-state difference between Radau and DOP853 over n equal arcs, each Radau arc
    started from the DOP853 state at the arc's start."""
    s = s6.copy()
    worst = 0.0
    for _ in range(n):
        d = cr3bp.propagate(system, s, period / n, rtol=1e-12, atol=1e-12)
        r = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, period / n),
            s,
            args=(system.mu,),
            method="Radau",
            rtol=1e-11,
            atol=1e-11,
        )
        worst = max(worst, float(np.linalg.norm(r.y[:, -1] - np.asarray(d.state_f))))
        s = np.asarray(d.state_f)
    return worst


def mirror_distance(s4: np.ndarray, period: float) -> float:
    """Minimum state-space distance from the mirror image of the perigee state to the orbit
    (dense DOP853 samples over one period); ~0 for a symmetric orbit."""
    mirror = np.array([s4[0], -s4[1], -s4[2], s4[3]])
    sol = solve_ivp(
        BASE.eom,
        (0.0, period),
        s4,
        args=(BASE.MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    tg = np.linspace(0.0, period, 200001)
    zz = sol.sol(tg)
    return float(np.min(np.linalg.norm(zz.T - mirror[None, :], axis=1)))


def evidence(fam: str) -> dict[str, Any]:
    rid, stem, p, q = FAMILIES[fam]
    seed = json.loads((CONT / f"{stem}_p.jsonl").read_text().splitlines()[0])
    assert seed["k"] == 0
    system = earth_moon_system()
    assert abs(system.mu - BASE.MU) < 1e-15
    s4 = np.array(seed["s_perigee"])
    s6 = np.array([s4[0], s4[1], 0.0, s4[2], s4[3], 0.0])
    orb = cp.correct_periodic(system, s6, seed["T"], tol=1e-10, max_iter=50)
    radau_ok, radau_dj = cp.crosscheck_periodic(system, orb)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orb.period),
        orb.state0,
        args=(system.mu,),
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    radau_closure = float(np.linalg.norm(sol.y[:, -1] - orb.state0))
    radau_arc_max = radau_arcs(system, orb.state0, orb.period, 4 * p)
    plain = cr3bp.propagate(system, orb.state0, orb.period, rtol=1e-12, atol=1e-12)
    closure_plain = float(np.linalg.norm(np.asarray(plain.state_f) - orb.state0))
    st = planar_stability_index(system, orb.state0, orb.period)
    k_arcs = arc_stm_product(system, orb.state0, orb.period, 4 * p)
    c4 = np.array([orb.state0[0], orb.state0[1], orb.state0[3], orb.state0[4]])
    cl = BASE.classify(c4, orb.period)
    n_in = cl["wind_E"] + orb.period / (2 * math.pi)
    return {
        "family": fam,
        "id": rid,
        "seed": {
            "file": f"data/1000_complement/continuation/{stem}_p.jsonl (line 1, k = 0)",
            "C": seed["C"],
            "T": seed["T"],
            "s_perigee": seed["s_perigee"],
            "closure_multiple_shooting_seed": seed["closure"],
        },
        "mass_ratio": system.mu,
        "lunit_km": system.l_km,
        "tunit_s": system.t_s,
        "state_nd": [float(x) for x in orb.state0],
        "period_nd": orb.period,
        "jacobi_constant": orb.jacobi,
        "converged": orb.converged,
        "closure_residual": orb.closure_residual,
        "closure_plain_dop853": closure_plain,
        "radau_ok_crosscheck_periodic": radau_ok,
        "radau_dj": radau_dj,
        "radau_closure_full_period": radau_closure,
        "radau_arc_max_diff": radau_arc_max,
        "radau_arcs": 4 * p,
        "k_par_full_period_stm": st.k_par,
        "k_perp": st.k_perp,
        "k_signed": st.k_signed,
        "k_par_arc_stm_product": k_arcs,
        "k_eig_agree": st.agree,
        "perigee_km": cl["perigee_km"],
        "perigee_alt_km": cl["perigee_alt_km"],
        "apogee_km": cl["apogee_km"],
        "periselene_km": cl["periselene_km"],
        "periselene_alt_km": cl["periselene_alt_km"],
        "max_moon_km": cl["max_moon_km"],
        "wind_E": cl["wind_E"],
        "wind_M": cl["wind_M"],
        "symmetric": cl["symmetric"],
        "perp_crossings": cl["perp_crossings"],
        "cycler_class_candidate": cl["cycler_class_candidate"],
        "period_days": orb.period * system.t_s / 86400.0,
        "inertial_N": n_in,
        "pq": [abs(round(n_in)), round(orb.period / (2 * math.pi))],
        "pq_expected": [p, q],
        "sense": "prograde" if n_in > 0 else "retrograde",
        "mirror_min_distance_nd": mirror_distance(c4, orb.period),
        "floors_km": {"earth": BASE.FLOOR_E, "moon": BASE.FLOOR_M},
        "geo_km": BASE.GEO_KM,
        "hill_km": BASE.HILL_KM,
    }


def main() -> None:
    out = {fam: evidence(fam) for fam in FAMILIES}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    for fam, e in out.items():
        print(
            fam,
            e["converged"],
            f"res {e['closure_residual']:.2e}",
            f"radau dj {e['radau_dj']:.1e} cl {e['radau_closure_full_period']:.1e}",
            f"arcs {e['radau_arc_max_diff']:.1e}",
            f"C {e['jacobi_constant']:.10f} T {e['period_nd']:.10f}",
            f"k {e['k_signed']:.3f} arcs {e['k_par_arc_stm_product']:.3f} perp {e['k_perp']:.3f}",
            f"pe {e['perigee_alt_km']:.0f} ps {e['periselene_alt_km']:.0f}",
            f"wE {e['wind_E']:.3f} wM {e['wind_M']:.3f} sym {e['symmetric']}",
            f"pq {e['pq']} {e['sense']} mirror {e['mirror_min_distance_nd']:.3e}",
            flush=True,
        )


if __name__ == "__main__":
    main()
