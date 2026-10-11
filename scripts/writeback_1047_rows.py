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
Writes data/1047_rows/writeback_v1.json; with --rows-out, also the draft catalogue rows (YAML text,
from the evidence file only). The catalogue edit itself is made in the catalogue window.

Usage: uv run python scripts/writeback_1047_rows.py [--rows-out <path>] [--rows-only]
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]
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


#: per family: (description, the #1048 full computed family (C range, ends, candidate range,
#: perigee/periselene altitudes, T, b), docs/notes/2026-10-08-1048-1000-coverage-gaps.md sec. 3)
FAMILY_TEXT: dict[str, tuple[str, str]] = {
    "A": (
        "5:2 prograde (inertial), winding about the Earth 3, about the Moon 0, asymmetric",
        "full computed C range 2.4022-3.0968 (low end: Earth floor; high end: loss of "
        "convergence, unexplained, 0.003 below Liang 2017's 5:2 C); cycler-class candidates "
        "over C 2.4034-2.9423 "
        "(perigee 225-35,752 km alt, periselene 11,637-35,997 km alt, T 12.696-13.093, b -275 to "
        "-103 from arc STMs)",
    ),
    "B": (
        "5:3 retrograde (inertial), winding about the Earth -8, about the Moon -2, asymmetric",
        "full computed C range 0.0159-0.9418 (low end: loss of convergence, unexplained; high end: "
        "Earth floor); cycler-class candidates over C 0.4203-0.9403 (perigee 200-35,530 km alt, "
        "periselene 13,265-22,444 km alt, T 18.877-19.241, b +97 to +149)",
    ),
    "D": (
        "7:3 prograde (inertial), winding about the Earth 4, about the Moon -1, asymmetric",
        "full computed C range 2.4136-2.5601 (low end: Earth floor; high end: loss of convergence, "
        "unexplained, the most fold-like stop); cycler-class candidates over C 2.4192-2.5601 "
        "(perigee 267-6,373 km alt, periselene 10,716-14,064 km alt, T 18.608-18.822, b +110 to "
        "+799)",
    ),
}
NAMES = {
    "A": "#1000 family A",
    "B": "#1000 family B",
    "D": "#1000 family D",
}
CASOLIVA_2010 = {
    "authors": [
        "Casoliva, J.",
        "Mondelo, J. M.",
        "Villac, B. F.",
        "Mease, K. D.",
        "Barrabes, E.",
        "Olle, M.",
    ],
    "year": 2010,
    "title": "Two Classes of Cycler Trajectories in the Earth-Moon System",
    "venue": "Journal of Guidance, Control, and Dynamics 33(5), pp. 1623-1640",
    "doi": "10.2514/1.46856",
    "url": "https://arc.aiaa.org/doi/10.2514/1.46856",
}


def row(e: dict[str, Any]) -> dict[str, Any]:
    fam = e["family"]
    we, wm = round(e["wind_E"]) + 0, round(e["wind_M"]) + 0
    p, q = e["pq"]
    desc, extent = FAMILY_TEXT[fam]
    sense = e["sense"]
    c, t, k = e["jacobi_constant"], e["period_nd"], e["k_signed"]
    name = (
        f"{NAMES[fam]}: ({p},{q}) {sense} asymmetric resonant Earth-Moon cycler "
        f"(CR3BP, C={c:.4f}, periselene {e['periselene_alt_km']:,.0f} km alt; known-class member)"
    )
    notes = (
        f"KNOWN-CLASS MEMBER (spec sec. 16.4; #875 (iii); owner decision 2026-10-11, 'Rows for A, "
        f"B, D only'). {NAMES[fam]}: {desc}. A computed member of the published class of near-"
        f"Keplerian p-q resonant Earth-Moon cyclers (Casoliva et al. 2010 Class 1, which includes "
        f"the asymmetric 7-3b/7-3c pair); NOT a literal copy of any published orbit and NOT "
        f"novel.\n\n"
        f"FAMILY (#1000, #1048, #1050): {extent}. This row is the #1000 grid member that seeded "
        f"the continuation ({e['seed']['file']}), re-closed with correct_periodic.\n\n"
        f"LITERATURE (#1047, docs/notes/2026-10-11-1047-earth-moon-families-literal-table-and-v3-"
        f"gate.md): literal table (every published Earth-Moon p:q orbit held as numbers: Casoliva "
        f"2010 Table 3, Vaquero 2013, Newton 1959, Casoliva 2008 seeds, Liang 2017, Liang 2020, "
        f"Hoelker-Winston 1968, Genova-Aldrin 2015) has NO member of this family. The v3 "
        f"literature gate returns 'published' at class level (Earth primary, Moon, topology "
        f"resonant); the gate does not compare p:q, sense or symmetry. Restrepo & Russell 2018: "
        f"excluded by construction (their database holds planar axisymmetric orbits only; this "
        f"orbit is asymmetric). Franz & Russell 2022: excluded by construction (their database "
        f"drops orbits that go more than about 350,000 km from the Moon; this orbit reaches "
        f"{e['max_moon_km']:,.0f} km).\n\n"
        f"VERIFICATION (V1; scripts/writeback_1047_rows.py, data/1047_rows/writeback_v1.json, "
        f"tests/search/test_1047_rows_writeback_v1.py): full-state single-shooting Newton "
        f"closure {e['closure_residual']:.1e}; independent Radau re-propagation keeps the "
        f"Jacobi constant to {e['radau_dj']:.1e} over one period and agrees with DOP853 arc by "
        f"arc ({e['radau_arcs']} arcs) to {e['radau_arc_max_diff']:.1e}. The orbit is strongly "
        f"unstable (k = {k:.1f}), so a whole-period closure depends on the integrator's step "
        f"sequence: {e['closure_plain_dop853']:.1e} for a plain DOP853 run and "
        f"{e['radau_closure_full_period']:.1e} for Radau (error amplification, not a closure "
        f"failure; the reason #1000 uses multiple shooting). Dense DOP853 re-propagation: perigee "
        f"{e['perigee_alt_km']:,.0f} km alt, periselene {e['periselene_alt_km']:,.0f} km alt "
        f"(inside the Hill radius {e['hill_km']:,.0f} km: a lunar encounter every period, "
        f"orbit_class cycler), winding about the Earth {we} and about the Moon "
        f"{wm}, no perpendicular x-axis crossing (asymmetric). NOT V2 (an unstable "
        f"orbit cannot satisfy V2-ballistic bounded drift).\n\n"
        f"MIRROR: the mirror image (x, -y, -vx, vy) of this orbit is a different orbit of the "
        f"same C and period (minimum state distance {e['mirror_min_distance_nd']:.3f} nondim); "
        f"it is not a separate row (owner decision: one row per family).\n\n"
        f"MODEL: planar CR3BP at this project's registry mu, l* and t*. No ephemeris claim.\n"
    )
    return {
        "id": e["id"],
        "name": name,
        "source": "this-project",
        "our_status": "known-class-member",
        "trajectory_regime": "ballistic",
        "model_assumption": "cr3bp",
        "cycler_class": "non-keplerian",
        "orbit_class": "cycler",
        "epoch_locked": False,
        "n_returns": "infinite",
        "validation_level": "V1",
        "delta_v_kms": None,
        "v_infinity_leveraging_dv_kms": None,
        "fleet_size": None,
        "flyby_mechanics": None,
        "primary": "Earth",
        "bodies": ["E", "Moon"],
        "sequence_canonical": "E-Moon",
        "sense": "n/a",
        "period": {
            "pair": "E-Moon",
            "k": 1,
            "years": None,
            "note": (
                f"T = {t!r} TU nondim (about 2 pi x {q} = {2 * math.pi * q:.4f}: the p-q "
                f"resonance); {e['period_days']:.6f} d via this project's registry t* "
                f"({e['tunit_s']:.4f} s). years: null (CR3BP, rotating-frame periodic orbit)."
            ),
        },
        "vinf_kms_at_encounters": [
            {
                "body": "E",
                "vinf_kms": None,
                "note": (
                    "CR3BP periodic orbit: conserved quantity is the Jacobi constant, not V_inf."
                ),
            },
            {"body": "Moon", "vinf_kms": None, "note": "Same -- Jacobi-constant model."},
        ],
        "orbit_elements": {
            "a_au": None,
            "e": None,
            "note": (
                "Planar CR3BP rotating-frame periodic orbit; Keplerian elements inapplicable. "
                "CR3BP identity tuple below."
            ),
            "cr3bp": {
                "family": (
                    f"#1000 family {fam}: ({p},{q}) {sense} asymmetric resonant Earth-Moon "
                    f"cycler (winding E {we}, M {wm}); member of the "
                    f"Casoliva 2010 Class 1 p-q resonant class"
                ),
                "mass_ratio": e["mass_ratio"],
                "libration_point": None,
                "jacobi_constant": e["jacobi_constant"],
                "period_nd": e["period_nd"],
                "stability_index": k,
                "state_nd": e["state_nd"],
                "lunit_km": e["lunit_km"],
                "tunit_s": e["tunit_s"],
                "perigee_km": e["perigee_km"],
                "periselene_km": e["periselene_km"],
                "apogee_km": e["apogee_km"],
                "max_moon_distance_km": e["max_moon_km"],
                "energy_regime": None,
                "tof_days_bounds": [e["period_days"], e["period_days"]],
            },
        },
        "source_ephemeris": None,
        "data_gaps": [
            {
                "path": "orbit_elements.cr3bp.state_nd",
                "kind": "derive",
                "note": (
                    "No published IC for this orbit or its family (#1047 literal table). "
                    "state_nd DERIVED: the #1000 multiple-shooting grid member that seeded the "
                    "family's continuation, re-closed with "
                    "cyclerfinder.search.cr3bp_periodic.correct_periodic at the registry mu. "
                    "Never a golden. This is why our_status is known-class-member, not "
                    "known-reproduction."
                ),
                "source_hint": (
                    "scripts/run_1000_shooting.py; scripts/run_1000_continue.py; "
                    "scripts/writeback_1047_rows.py; data/1047_rows/writeback_v1.json."
                ),
                "todo_ref": "#1047",
            },
            {
                "path": "orbit_elements.cr3bp.stability_index",
                "kind": "derive",
                "note": (
                    "Casoliva 2010 Eq. 6-8 convention, as on the casoliva-7-3b/7-3c rows: k = "
                    "max(|k_par|, |k_perp|) with its sign, from the full-period monodromy "
                    f"(earth_moon_resonant_families.planar_stability_index; k_par "
                    f"{e['k_par_full_period_stm']:.3f}, k_perp {e['k_perp']:.3f}). The product "
                    f"of {e['radau_arcs']} per-arc STMs gives k_par "
                    f"{e['k_par_arc_stm_product']:.3f}. Not comparable number-for-number with "
                    "the Barden half-period nu rows."
                ),
                "source_hint": (
                    "src/cyclerfinder/search/earth_moon_resonant_families.py::"
                    "planar_stability_index; data/1047_rows/writeback_v1.json."
                ),
                "todo_ref": "#1047",
            },
        ],
        "legs": [],
        "corroborating_sources": [
            {
                **copy.deepcopy(CASOLIVA_2010),
                "note": (
                    "The published CLASS source (not first_published): Class 1, near-Keplerian "
                    "p-q resonant Earth-Moon cyclers, including the asymmetric 7-3b/7-3c pair. "
                    f"Table 3 tabulates p-q = 1-2, 2-1, 3-2 and 7-3 only; it has no ({p},{q}) "
                    f"{sense} orbit."
                ),
            }
        ],
        "priority_date": "2026-10-08",
        "notes": notes,
    }


def write_rows(ev: dict[str, Any], path: Path) -> None:
    rows = [row(ev[fam]) for fam in FAMILIES]
    text = yaml.safe_dump(rows, sort_keys=False, width=100, allow_unicode=False)
    path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows-out", type=Path, default=None)
    ap.add_argument("--rows-only", action="store_true")
    args = ap.parse_args()
    if args.rows_only:
        out = json.loads(OUT.read_text())
    else:
        out = {fam: evidence(fam) for fam in FAMILIES}
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(out, indent=1) + "\n")
    if args.rows_out is not None:
        write_rows(out, args.rows_out)
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
