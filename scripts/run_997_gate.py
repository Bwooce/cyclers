"""#997 known-class gate: compare the continued families (data/997_lineage/*.jsonl) with the
catalogue's Earth-Moon CR3BP rows, and summarise the families and cycler-class candidates.

For every planar Earth-Moon catalogue row with a state and a period, the row's state is
propagated in the project CR3BP at the row's own mass ratio. The script records C, T, the
perigee and periselene over one period, and the perpendicular x-axis crossings (y = 0,
|xdot| < 1e-6), which shows whether the row is x-axis symmetric. Each row is then matched
against every #997 member by C, T and the perpendicular-crossing state. A match means
|dC| < 2e-3, |dT|/T < 2e-3 and a perpendicular crossing within 2e-3 in (x, ydot) of a
member's start or T/2 crossing.

Output: data/997_lineage/gate.json and data/997_lineage/summary.json.
"""

from __future__ import annotations

import glob
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom, jacobi_constant
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "997_lineage"
L_KM = 384400.0
R_E, R_M = 6378.137, 1737.4

# Segment table: which checkpoint files make up which family (each file is one direction from
# one seed; overlapping segments are the same orbits). Pre-registered family ids (note sec. 0.2).
FAMILY_OF = {
    "F1-0": "F1 Schwaniger retrograde cislunar",
    "F2-1": "F2 Newton 1/2 retrograde",
    "F2-2": "F2 Newton 1/2 retrograde",
    "F2-3": "F2 Newton 1/2 retrograde",
    "F2-4": "F2 Newton 1/2 retrograde",
    "F2-5": "F2 Newton 1/2 retrograde",
    "F3-6": "F3 Newton 1/2 direct",
    "F3-7": "F3 Newton 1/2 direct",
    "F3-8": "F3 Newton 1/2 direct",
    "F3-9": "F3 Newton 1/2 direct",
    "F3-10": "F3 Newton 1/2 direct",
    "F3-11": "F3 Newton 1/2 direct",
    "F4-2_5-d-rho0.06": "F4 2/5 direct (converged onto the doubled F3 orbit)",
    "F4-3_7-r-rho0.06": "F4 3/7 retrograde",
    "F4-2_3-r-rho0.02": "F4 2/3 retrograde",
    "F4-2_3-r-rho0.06": "F4 2/3 retrograde",
    "F4-2_3-r-rho0.2": "F4 2/3 retrograde",
    "F4-3_4-r-rho0.1": "F4 3/4 retrograde",
    "F5-fig90-n4": "F5 (mu homotopy left the Hoelker-Winston orbit; not an F5 member)",
}


def perpendicular_crossings(state: np.ndarray, period: float, mu: float) -> list[list[float]]:
    def ev(t: float, y: np.ndarray, m: float) -> float:
        return float(y[1])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, period),
        state,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=ev,
    )
    out = []
    for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True):
        if abs(y[3]) < 1e-6 and abs(y[2]) < 1e-9 and abs(y[5]) < 1e-9:
            out.append([float(t), float(y[0]), float(y[4])])
    if abs(state[1]) < 1e-9 and abs(state[3]) < 1e-6:
        out.insert(0, [0.0, float(state[0]), float(state[4])])
    return out


def extrema(state: np.ndarray, period: float, mu: float) -> dict[str, float]:
    def d1(t: float, y: np.ndarray, m: float) -> float:
        return float((y[0] + m) * y[3] + y[1] * y[4] + y[2] * y[5])

    def d2(t: float, y: np.ndarray, m: float) -> float:
        return float((y[0] - 1 + m) * y[3] + y[1] * y[4] + y[2] * y[5])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, period),
        state,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=[d1, d2],
    )
    r1 = [math.dist(y[:3], (-mu, 0, 0)) for y in sol.y_events[0]]
    r2 = [math.dist(y[:3], (1 - mu, 0, 0)) for y in sol.y_events[1]]
    return {
        "perigee_alt_km": min(r1) * L_KM - R_E if r1 else float("nan"),
        "periselene_alt_km": min(r2) * L_KM - R_M if r2 else float("nan"),
        "closure": float(np.linalg.norm(sol.y[:, -1] - state)),
    }


def max_moon_distance_km(m: dict[str, Any]) -> float:
    """Largest distance from the Moon over one period (planar member from its start state)."""
    mu = m["mu"]
    s0 = np.array([m["x0"], 0.0, 0.0, 0.0, m["ydot0"], 0.0])

    def d2(t: float, y: np.ndarray, mm: float) -> float:
        return float((y[0] - 1 + mm) * y[3] + y[1] * y[4])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, m["T"]),
        s0,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=d2,
    )
    return max(math.hypot(y[0] - 1 + mu, y[1]) for y in sol.y_events[0]) * L_KM


def load_members() -> list[dict[str, Any]]:
    rows = []
    for f in sorted(glob.glob(str(DATA / "F*_[pm].jsonl"))):
        seed = Path(f).stem.rsplit("_", 1)[0]
        for line in Path(f).read_text().splitlines():
            m = json.loads(line)
            m["seed"] = seed
            m["file"] = Path(f).name
            m["family"] = FAMILY_OF.get(seed, seed)
            rows.append(m)
    return rows


def main() -> None:
    preflight_search(
        task_no=997,
        region_id="em-both-primary-lineage-known-class-gate",
        method=MethodCapability(
            genome="catalogue Earth-Moon rows vs #997 continued families",
            corrector="none (propagation and matching only)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    members = load_members()
    chains: dict[str, list[dict[str, Any]]] = {}
    for m in members:
        chains.setdefault(m["file"], []).append(m)
    cat = yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())
    gate = []
    for r in cat:
        c = ((r.get("orbit_elements") or {}).get("cr3bp")) or {}
        st, per, mu = c.get("state_nd"), c.get("period_nd"), c.get("mass_ratio")
        if r.get("primary") != "Earth" or not st or not per or not mu:
            continue
        s = np.array(st, dtype=float)
        if abs(s[2]) > 1e-9 or abs(s[5]) > 1e-9:
            continue  # spatial rows: outside this planar lineage
        jc = jacobi_constant(s, mu)
        perp = perpendicular_crossings(s, per, mu)
        ext = extrema(s, per, mu)
        hits = []
        for fname, chain in chains.items():
            for a, b in itertools.pairwise(chain):
                if (a["C"] - jc) * (b["C"] - jc) > 0 or a["C"] == b["C"]:
                    continue
                w = (jc - a["C"]) / (b["C"] - a["C"])
                ip = {k: a[k] + w * (b[k] - a[k]) for k in ("x0", "ydot0", "x_half", "T")}
                ip.update(
                    {
                        k: a[k] + w * (b[k] - a[k])
                        for k in ("perigee_alt_km", "periselene_alt_km", "b_h")
                    }
                )
                close = any(
                    abs(p[1] - ip["x0"]) < 1e-3 and abs(p[2] - ip["ydot0"]) < 1e-3 for p in perp
                ) or any(abs(p[1] - ip["x_half"]) < 1e-3 for p in perp)
                hits.append(
                    {
                        "family": a["family"],
                        "file": fname,
                        "between_k": [a["k"], b["k"]],
                        "dT": ip["T"] - per,
                        "interp": ip,
                        "crossing_match": bool(close and abs(ip["T"] - per) / per < 1e-3),
                    }
                )
        best = sorted(hits, key=lambda h: (not h["crossing_match"], abs(h["dT"])))[:3]
        gate.append(
            {
                "id": r["id"],
                "orbit_class": r.get("orbit_class"),
                "C": jc,
                "T": per,
                "x_axis_symmetric": bool(perp),
                "perpendicular_crossings": perp,
                **ext,
                "matches": best,
            }
        )
        tag = "MATCH " + best[0]["family"] if best and best[0]["crossing_match"] else "-"
        print(
            f"{r['id']:45s} C={jc:8.4f} T={per:8.4f} sym={bool(perp)!s:5s} "
            f"pe={ext['perigee_alt_km']:9.0f} ps={ext['periselene_alt_km']:8.0f} {tag}"
        )
    (DATA / "gate.json").write_text(json.dumps(gate, indent=1) + "\n")

    # family summary and cycler-class candidates
    fams: dict[str, dict[str, Any]] = {}
    for m in members:
        f = fams.setdefault(
            m["family"],
            {"n_members": 0, "C": [], "T": [], "stops": [], "candidates": []},
        )
        f["n_members"] += 1
        f["C"].append(m["C"])
        f["T"].append(m["T"])
        if m.get("stop"):
            f["stops"].append({"seed": m["seed"], "stop": m["stop"], "C": m["C"]})
        if m["cycler_class_candidate"]:
            m["r2max_km"] = max_moon_distance_km(m)
            f["candidates"].append(
                {
                    k: m[k]
                    for k in (
                        "seed",
                        "k",
                        "C",
                        "T",
                        "T_days",
                        "perigee_alt_km",
                        "periselene_alt_km",
                        "b_h",
                        "b_v",
                        "x0",
                        "ydot0",
                        "r2max_km",
                    )
                }
            )
    summary = {}
    for name, f in fams.items():
        cand = f["candidates"]
        summary[name] = {
            "n_members": f["n_members"],
            "C_range": [min(f["C"]), max(f["C"])],
            "T_range": [min(f["T"]), max(f["T"])],
            "stops": f["stops"],
            "n_candidates": len(cand),
            "candidate_C_range": [min(c["C"] for c in cand), max(c["C"] for c in cand)]
            if cand
            else None,
            "candidates": cand,
            # Franz-Russell 2022 discard orbits ever > ~350,000 km from the Moon (note sec. 0.6)
            "candidate_min_of_max_moon_distance_km": min(c["r2max_km"] for c in cand)
            if cand
            else None,
            "candidates_excluded_by_franz_russell": all(c["r2max_km"] > 350000.0 for c in cand)
            if cand
            else None,
        }
        print(f"{name}: {f['n_members']} members, {len(cand)} candidates")
    (DATA / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
