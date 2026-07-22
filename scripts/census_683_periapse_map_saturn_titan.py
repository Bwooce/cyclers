"""#683 -- periapse Poincaré-map cartography at Saturn-Titan.

Builds the one classical seedless discovery map this codebase genuinely lacked
(periapsis Poincaré maps, Davis-Howell lineage; see
``src/cyclerfinder/search/periapse_map.py``) and:

  1. CALIBRATION (Phase B): four hand-checkable correctness checks of the
     periapse machinery before any structure is trusted.
  2. POSITIVE CONTROL (Phase C): reproduce the Davis & Howell 2011 (Acta
     Astronautica 69, 1038-1049) Saturn-Titan and Sun-Saturn initial-condition
     periapse maps (its Figs. 3/4) -- the escape-L1 / escape-L2 / impact /
     captured lobe structure at the paper's own published Jacobi values. The
     paper is freely available from the authors' Purdue site and digested at
     ``docs/notes/2026-07-22-683-digest-davis-howell-2011-periapse-maps.md``.
  3. DISCOVERY (Phase D): search the captured region at Saturn-Titan (J=J2,
     both gateways open) for a REPEATING capture -> escape -> re-capture
     itinerary, read geometrically off the map's own lobe structure -- a
     genuinely different lens on the same transport question `#664`/`#685`'s
     set-oriented (statistical) pipeline probed at Sun-Earth.

NO catalogue writeback. A repeating-itinerary candidate (if any) is reported
for adjudication only; a clean negative is the expected, legitimate outcome
per the task's own moderate-low novelty ceiling.

Run (resumable; each phase skips if its output file already exists):
    uv run python scripts/census_683_periapse_map_saturn_titan.py
    uv run python scripts/census_683_periapse_map_saturn_titan.py --assemble

Outputs -> data/found/683_periapse_map/
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.periapse_map as pm  # noqa: E402

OUT = ROOT / "data" / "found" / "683_periapse_map"
OUT.mkdir(parents=True, exist_ok=True)

# Davis & Howell 2011 published search energies (their exact quoted values).
J1_SUN_SATURN = 3.0173046596239  # paper Fig. 3a/4a (J < J_L2, both necks open)
J2_SATURN_TITAN = 3.015311017945150  # paper Fig. 3b/4b/20 (J < J_L1, both necks open)


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------
# Phase A: L-point / regime cross-check
# --------------------------------------------------------------------------
def phase_a() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for tag, prim, sec, jval in [
        ("sun_saturn", "Sun", "Saturn", J1_SUN_SATURN),
        ("saturn_titan", "Saturn", "Titan", J2_SATURN_TITAN),
    ]:
        st = pm.build_periapse_map_system(prim, sec)
        out[tag] = {
            "mu": st.mu,
            "r_hill": st.r_hill,
            "x_l1": st.x_l1,
            "x_l2": st.x_l2,
            "c_l1": st.c_l1,
            "c_l2": st.c_l2,
            "j_search": jval,
            "j_below_c_l2": jval < st.c_l2,
            "j_below_c_l1": jval < st.c_l1,
            "gap_below_c_l2": st.c_l2 - jval,
        }
    return out


# --------------------------------------------------------------------------
# Phase B: calibration
# --------------------------------------------------------------------------
def phase_b() -> dict[str, Any]:
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    rng = np.random.default_rng(683)
    res: dict[str, Any] = {}

    # C1: periapsis detection is EXACT (rdot==0 to machine precision), not
    # sampling-based, and rddot>0 (a genuine minimum). Check on many next-peri.
    max_rdot = 0.0
    min_rddot = math.inf
    checked = 0
    while checked < 40:
        xp = rng.uniform(-0.6, 0.6)
        yp = rng.uniform(-0.6, 0.6)
        s0 = pm.construct_periapse_state(
            (1 - mu) + xp * st.r_hill, yp * st.r_hill, J2_SATURN_TITAN, mu
        )
        if s0 is None or not pm.is_periapsis(s0, mu):
            continue
        nxt = pm.next_periapse(s0, st)
        if nxt is None:
            continue
        max_rdot = max(max_rdot, abs(pm.secondary_radial_rate(nxt, mu)))
        min_rddot = min(min_rddot, pm.secondary_radial_accel(nxt, mu))
        checked += 1
    res["c1_max_abs_rdot_at_detected_periapsis"] = max_rdot
    res["c1_min_rddot_at_detected_periapsis"] = min_rddot
    res["c1_n_checked"] = checked

    # C2: Jacobi constant conserved across the periapse-to-periapse map.
    worst_dj = 0.0
    checked = 0
    while checked < 40:
        xp = rng.uniform(-0.6, 0.6)
        yp = rng.uniform(-0.6, 0.6)
        s0 = pm.construct_periapse_state(
            (1 - mu) + xp * st.r_hill, yp * st.r_hill, J2_SATURN_TITAN, mu
        )
        if s0 is None or not pm.is_periapsis(s0, mu):
            continue
        nxt = pm.next_periapse(s0, st)
        if nxt is None:
            continue
        worst_dj = max(worst_dj, abs(cr3bp.jacobi_constant(nxt, mu) - J2_SATURN_TITAN))
        checked += 1
    res["c2_worst_abs_dj_across_map"] = worst_dj

    # C3: round-trip parametrisation is exact -- a constructed periapse state
    # reproduces the target Jacobi constant AND is exactly on the section.
    worst_j = 0.0
    worst_rdot = 0.0
    for _ in range(200):
        xp = rng.uniform(-1.0, 1.0)
        yp = rng.uniform(-1.0, 1.0)
        s0 = pm.construct_periapse_state(
            (1 - mu) + xp * st.r_hill, yp * st.r_hill, J2_SATURN_TITAN, mu
        )
        if s0 is None:
            continue
        worst_j = max(worst_j, abs(cr3bp.jacobi_constant(s0, mu) - J2_SATURN_TITAN))
        worst_rdot = max(worst_rdot, abs(pm.secondary_radial_rate(s0, mu)))
    res["c3_worst_construct_j_err"] = worst_j
    res["c3_worst_construct_rdot"] = worst_rdot

    # C4: fate classifier is geometrically sane on hand-built cases.
    #  (a) a state fired straight at Titan impacts;
    #  (b) a state placed just inside the L2 neck with outward radial energy
    #      escapes through L2 (not L1); mirror for L1.
    # (a) periapse aimed to graze Titan: tiny periapse radius -> IMPACT.
    s_impact = pm.construct_periapse_state(
        (1 - mu) + 0.5 * st.secondary_radius_nd, 0.0, J2_SATURN_TITAN, mu
    )
    fate_impact = (
        pm.classify_fate(s_impact, st, max_revs=6)[0].name if s_impact is not None else "NONE"
    )
    res["c4a_grazing_periapsis_fate"] = fate_impact
    # (b) directional sanity: escape-L1 periapses sit on the L1 side, escape-L2
    # on the L2 side (centroid check over a coarse grid).
    esc_l1_x: list[float] = []
    esc_l2_x: list[float] = []
    for xp in np.linspace(-1.15, 1.15, 40):
        for yp in np.linspace(-0.9, 0.9, 34):
            s0 = pm.construct_periapse_state(
                (1 - mu) + xp * st.r_hill, yp * st.r_hill, J2_SATURN_TITAN, mu
            )
            if s0 is None or not pm.is_periapsis(s0, mu):
                continue
            f = pm.classify_fate(s0, st, max_revs=1)[0]
            if f is pm.PeriapseFate.ESCAPE_L1:
                esc_l1_x.append(xp)
            elif f is pm.PeriapseFate.ESCAPE_L2:
                esc_l2_x.append(xp)
    res["c4b_escape_l1_centroid_xp"] = float(np.mean(esc_l1_x)) if esc_l1_x else None
    res["c4b_escape_l2_centroid_xp"] = float(np.mean(esc_l2_x)) if esc_l2_x else None
    res["c4b_n_escape_l1"] = len(esc_l1_x)
    res["c4b_n_escape_l2"] = len(esc_l2_x)
    return res


# --------------------------------------------------------------------------
# Phase C: positive-control initial-condition maps
# --------------------------------------------------------------------------
def _ic_map(
    st: pm.PeriapseMapSystem,
    jacobi: float,
    *,
    max_revs: int,
    nx: int,
    ny: int,
    xr: tuple[float, float] = (-1.25, 1.25),
    yr: tuple[float, float] = (-1.05, 1.05),
) -> dict[str, Any]:
    mu = st.mu
    xs = np.linspace(xr[0], xr[1], nx)
    ys = np.linspace(yr[0], yr[1], ny)
    fate_grid = np.full((ny, nx), -1, dtype=np.int8)  # -1 outside contour
    code = {
        pm.PeriapseFate.IMPACT: 0,
        pm.PeriapseFate.ESCAPE_L1: 1,
        pm.PeriapseFate.ESCAPE_L2: 2,
        pm.PeriapseFate.CAPTURED: 3,
    }
    captured_states: list[list[float]] = []
    t0 = time.time()
    for j, yp in enumerate(ys):
        for i, xp in enumerate(xs):
            s0 = pm.construct_periapse_state((1 - mu) + xp * st.r_hill, yp * st.r_hill, jacobi, mu)
            if s0 is None or not pm.is_periapsis(s0, mu):
                continue
            fate, _ = pm.classify_fate(s0, st, max_revs=max_revs)
            fate_grid[j, i] = code[fate]
            if fate is pm.PeriapseFate.CAPTURED:
                captured_states.append([float(v) for v in s0])
        if j % 20 == 0:
            _log(f"    ic_map row {j}/{ny} ({time.time() - t0:.0f}s)")
    counts = {
        name: int(np.sum(fate_grid == c))
        for name, c in [("impact", 0), ("escape_l1", 1), ("escape_l2", 2), ("captured", 3)]
    }
    n_in = int(np.sum(fate_grid >= 0))

    # Centroids of the escape sets on the map (directional lobe check).
    def _centroid(cval: int) -> tuple[float | None, float | None]:
        idx = np.argwhere(fate_grid == cval)
        if idx.size == 0:
            return None, None
        xp_c = float(np.mean([xs[i] for _, i in idx]))
        yp_c = float(np.mean([ys[j] for j, _ in idx]))
        return xp_c, yp_c

    l1c = _centroid(1)
    l2c = _centroid(2)
    return {
        "jacobi": jacobi,
        "max_revs": max_revs,
        "nx": nx,
        "ny": ny,
        "xr": list(xr),
        "yr": list(yr),
        "n_in_contour": n_in,
        "counts": counts,
        "fractions": {k: (v / n_in if n_in else 0.0) for k, v in counts.items()},
        "escape_l1_centroid_xp_yp": list(l1c),
        "escape_l2_centroid_xp_yp": list(l2c),
        "captured_states": captured_states,
        "fate_grid": fate_grid.tolist(),
        "wall_s": time.time() - t0,
    }


def phase_c() -> dict[str, Any]:
    out: dict[str, Any] = {}
    st_st = pm.build_periapse_map_system("Saturn", "Titan")
    st_ss = pm.build_periapse_map_system("Sun", "Saturn")
    _log("  Saturn-Titan J2, 1 rev")
    out["saturn_titan_1rev"] = _ic_map(st_st, J2_SATURN_TITAN, max_revs=1, nx=160, ny=140)
    _log("  Saturn-Titan J2, 6 rev")
    out["saturn_titan_6rev"] = _ic_map(st_st, J2_SATURN_TITAN, max_revs=6, nx=160, ny=140)
    _log("  Sun-Saturn J1, 1 rev (cross-check)")
    out["sun_saturn_1rev"] = _ic_map(st_ss, J1_SUN_SATURN, max_revs=1, nx=140, ny=120)
    return out


# --------------------------------------------------------------------------
# Phase D: capture -> escape -> re-capture repeating-itinerary search
# --------------------------------------------------------------------------
_PHASE_D_T_TOTAL = 66.0 * math.pi  # ~33 Titan periods (Davis-Howell 2011 Fig.-18 horizon)
_PHASE_D_MAX_SEEDS = 900
_PHASE_D_CHUNK = 300  # seeds processed per invocation (keeps each call < 10 min)


def _phase_d_seeds(captured_states: list[list[float]]) -> list[list[float]]:
    seeds = list(captured_states)
    if len(seeds) > _PHASE_D_MAX_SEEDS:
        step = len(seeds) / float(_PHASE_D_MAX_SEEDS)
        seeds = [seeds[int(k * step)] for k in range(_PHASE_D_MAX_SEEDS)]
    return seeds


def _analyse_seed(
    s0: NDArray[np.float64], st: pm.PeriapseMapSystem, near_rp: float, excursion_rp: float
) -> tuple[bool, bool, float]:
    """Analyse one captured seed's long-term (fixed-time-horizon) itinerary.

    Returns ``(had_excursion, recaptured_after_excursion,
    min_excursion_return_residual_rH)``.

    Each periapsis is a closest-approach-to-Titan. A periapse with
    ``r_p <= near_rp`` is a genuine near-Titan (captured) periapse; a periapse
    with ``r_p >= excursion_rp`` (beyond the Hill radius) is an EXCURSION -- the
    trajectory's closest approach that revolution was outside Titan's Hill
    sphere, i.e. it is orbiting Saturn, not Titan. A capture -> escape ->
    re-capture itinerary is: near-Titan periapse(s), then >=1 excursion
    periapse, then near-Titan periapse(s) again. The residual is the closest
    periapse-map (Hill-radius-unit) return between a pre-excursion near-Titan
    periapse and a post-excursion near-Titan periapse -- a near-zero residual
    would be a repeating temporary-capture cycler at Titan.
    """
    events, _impacted = pm.collect_titan_periapses(s0, st, t_total=_PHASE_D_T_TOTAL)
    if len(events) < 3:
        return False, False, math.inf
    coords: list[tuple[float, float, float]] = []  # (x_p, y_p, r_p/rH)
    for _t, state in events:
        x_p, y_p, r_p, _ = pm.periapse_map_coords(state, st)
        coords.append((x_p, y_p, r_p / st.r_hill))
    had_excursion = any(rp >= excursion_rp for _, _, rp in coords)
    if not had_excursion:
        return False, False, math.inf
    near_before: tuple[float, float] | None = None
    seen_excursion = False
    recaptured = False
    seed_best = math.inf
    for x_p, y_p, rp in coords:
        if rp >= excursion_rp:
            seen_excursion = True
            continue
        if rp <= near_rp:
            if seen_excursion and near_before is not None:
                recaptured = True
                d = math.hypot(x_p - near_before[0], y_p - near_before[1])
                seed_best = min(seed_best, d)
                seen_excursion = False  # reset; look for the next excursion cycle
            near_before = (x_p, y_p)
    return had_excursion, recaptured, seed_best


def phase_d_chunk(captured_states: list[list[float]]) -> bool:
    """Process one chunk of captured seeds; checkpoint to phase_d_state.json.
    Returns True when ALL seeds are done (and writes phase_d.json)."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    near_rp = 0.6  # near-Titan periapse ceiling (Hill-radius units)
    excursion_rp = 1.0  # excursion floor: closest approach beyond the Hill radius
    seeds = _phase_d_seeds(captured_states)
    n_seeds = len(seeds)

    state_path = OUT / "phase_d_state.json"
    if state_path.exists():
        state = json.loads(state_path.read_text())
    else:
        state = {
            "cursor": 0,
            "n_excursion": 0,
            "n_escape_recapture": 0,
            "residuals": [],
            "best_residual": None,
            "best_record": None,
        }
    start = int(state["cursor"])
    end = min(start + _PHASE_D_CHUNK, n_seeds)
    t0 = time.time()
    for k in range(start, end):
        s0 = np.array(seeds[k], dtype=np.float64)
        had_excursion, recaptured, seed_best = _analyse_seed(s0, st, near_rp, excursion_rp)
        if had_excursion:
            state["n_excursion"] += 1
        if recaptured:
            state["n_escape_recapture"] += 1
            state["residuals"].append(seed_best)
            if state["best_residual"] is None or seed_best < state["best_residual"]:
                state["best_residual"] = seed_best
                state["best_record"] = {
                    "seed_index": k,
                    "seed_state": [float(v) for v in s0],
                    "min_excursion_return_residual_rH": seed_best,
                }
        if (k - start) % 50 == 0:
            _log(f"    phase_d seed {k}/{n_seeds} ({time.time() - t0:.0f}s)")
    state["cursor"] = end
    state_path.write_text(json.dumps(state))
    _log(f"  phase_d chunk done: {end}/{n_seeds}")
    if end < n_seeds:
        return False
    residuals = state["residuals"]
    _write(
        "phase_d.json",
        {
            "n_captured_seeds_probed": n_seeds,
            "t_total_nondim": _PHASE_D_T_TOTAL,
            "t_total_titan_periods": _PHASE_D_T_TOTAL / (2.0 * math.pi),
            "near_titan_rp_over_rH": near_rp,
            "excursion_rp_over_rH": excursion_rp,
            "n_excursion_trajectories": state["n_excursion"],
            "n_escape_recapture_trajectories": state["n_escape_recapture"],
            "min_excursion_return_residual_rH": state["best_residual"],
            "median_excursion_return_residual_rH": (
                float(np.median(residuals)) if residuals else None
            ),
            "best_record": state["best_record"],
        },
    )
    return True


# --------------------------------------------------------------------------
def _write(name: str, obj: Any) -> None:
    (OUT / name).write_text(json.dumps(obj, indent=2))


def _read(name: str) -> Any:
    return json.loads((OUT / name).read_text())


def assemble() -> None:
    a = _read("phase_a.json")
    b = _read("phase_b.json")
    c = _read("phase_c.json")
    d = _read("phase_d.json")
    print("\n===== #683 PERIAPSE-MAP SATURN-TITAN SUMMARY =====\n")
    print("Phase A -- L-point / regime cross-check:")
    for tag, v in a.items():
        print(
            f"  {tag}: mu={v['mu']:.6e} C_L1={v['c_l1']:.9f} C_L2={v['c_l2']:.9f} "
            f"J={v['j_search']:.9f} (J<C_L2={v['j_below_c_l2']}, gap={v['gap_below_c_l2']:.2e})"
        )
    print("\nPhase B -- calibration:")
    for k, v in b.items():
        print(f"  {k}: {v}")
    print("\nPhase C -- positive-control IC maps (fate fractions):")
    for tag, m in c.items():
        fr = {k: round(x, 3) for k, x in m["fractions"].items()}
        print(
            f"  {tag}: n_in={m['n_in_contour']} fractions={fr} "
            f"escL1_centroid_xp={m['escape_l1_centroid_xp_yp'][0]} "
            f"escL2_centroid_xp={m['escape_l2_centroid_xp_yp'][0]}"
        )
    print("\nPhase D -- capture->escape->re-capture repeating-itinerary search:")
    for k, v in d.items():
        if k == "best_record":
            continue
        print(f"  {k}: {v}")
    if d.get("best_record"):
        print(f"  best_record: {d['best_record']}")
    print("\n=================================================\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemble", action="store_true")
    args = ap.parse_args()
    if args.assemble:
        assemble()
        return

    if not (OUT / "phase_a.json").exists():
        _log("Phase A: L-point / regime cross-check")
        _write("phase_a.json", phase_a())
    if not (OUT / "phase_b.json").exists():
        _log("Phase B: calibration")
        _write("phase_b.json", phase_b())
    if not (OUT / "phase_c.json").exists():
        _log("Phase C: positive-control IC maps")
        _write("phase_c.json", phase_c())
    if not (OUT / "phase_d.json").exists():
        _log("Phase D: capture->escape->re-capture search (chunked)")
        c = _read("phase_c.json")
        captured = c["saturn_titan_6rev"]["captured_states"]
        done = False
        while not done:
            done = phase_d_chunk(captured)
    _log("ALL PHASES DONE")
    assemble()


if __name__ == "__main__":
    main()
