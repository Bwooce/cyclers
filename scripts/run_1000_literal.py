"""#1000: explicit literal collision table of the asymmetric families A (5:2 prograde),
B (5:3 retrograde) and D (7:3 prograde) against every published Earth-Moon p:q orbit held as
numbers: Casoliva 2010 Table 3 (all 16 rows, from the #780 module transcription), Vaquero 2013
(the 6 catalogue rows), Newton 1959 Table 1 (11 rows, #997 seeds), Casoliva 2008 Table 2
seeds (mu = 1e-6; digest), Liang, Xu & Xu 2017 (two orbits; digest) and Liang et al. 2020
(one 2:1 member; digest). Hoelker-Winston 1968 and Genova-Aldrin 2015 print no orbit with these
p:q values (no n* of 5/2, 5/3 or 7/3; Genova-Aldrin give no C or state) and enter as one row each.

For an orbit with a state, the inertial revolutions about the Earth per period are
N = w_E + T / (2 pi) (w_E = winding number about the Earth in the rotating frame), so
p:q = |N| : round(T / (2 pi)) and the sense is sign(N).

Match rule (note sec. 5.1b): same p:q and sense, |dC| < 0.01 to the family's computed C range,
|dT|/T < 1e-3, perigee and periselene altitudes within 10%, Earth-Moon mu.
Output: data/1000_complement/literal_check.json and literal_check.md.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.earth_moon_resonant_families import TABLE3_ROWS, table3_seed_state

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_complement as base

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1000_complement"
FAMILIES = {
    "A": {
        "pq": (5, 2),
        "sense": 1,
        "C": (2.4109, 2.5916),
        "T": (12.6975, 12.7303),
        "pe": (495, 8748),
        "ps": (11772, 16348),
    },
    "B": {
        "pq": (5, 3),
        "sense": -1,
        "C": (0.6703, 0.9418),
        "T": (18.86, 18.91),
        "pe": (200, 14548),
        "ps": (17876, 22444),
    },
    "D": {
        "pq": (7, 3),
        "sense": 1,
        "C": (2.4136, 2.5601),
        "T": (18.65, 18.82),
        "pe": (267, 6373),
        "ps": (10716, 14064),
    },
}


def topo(s4: np.ndarray, period: float) -> dict[str, Any]:
    tg = np.linspace(0.0, period, 40001)
    sol = solve_ivp(
        base.eom,
        (0.0, period),
        s4,
        args=(base.MU,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=tg,
    )
    x, y = sol.y[0], sol.y[1]
    w_e = (np.unwrap(np.arctan2(y, x + base.MU))[-1] - math.atan2(y[0], x[0] + base.MU)) / (
        2 * math.pi
    )
    r1 = np.hypot(x + base.MU, y) * base.L_KM
    r2 = np.hypot(x - 1 + base.MU, y) * base.L_KM
    q = round(period / (2 * math.pi))
    n_in = w_e + period / (2 * math.pi)
    return {
        "wE": float(w_e),
        "N_inertial": float(n_in),
        "pq": [abs(round(n_in)), q],
        "sense": 1 if n_in > 0 else -1,
        "perigee_alt_km": float(r1.min() - base.R_E),
        "periselene_alt_km": float(r2.min() - base.R_M),
    }


def verdict(row: dict[str, Any]) -> dict[str, str]:
    out = {}
    for name, f in FAMILIES.items():
        pq = tuple(row.get("pq") or ())
        if not pq:
            out[name] = row.get("why_no_pq", "no p:q")
            continue
        if pq != f["pq"]:
            out[name] = f"p:q {pq[0]}:{pq[1]} differs"
            continue
        if row.get("sense") is not None and row["sense"] != f["sense"]:
            out[name] = "same p:q, opposite sense"
            continue
        if row.get("mu_note"):
            out[name] = f"same p:q; class-level only ({row['mu_note']})"
            continue
        c = row.get("C")
        same = (
            "same p:q and sense" if row.get("sense") is not None else "same p:q (sense not printed)"
        )
        if c is None or not (f["C"][0] - 0.01 <= c <= f["C"][1] + 0.01):
            out[name] = f"{same}; C {c} outside the family's computed range"
            continue
        out[name] = "CANDIDATE MATCH: check T, perigee, periselene"
    return out


def main() -> None:
    preflight_search(
        task_no=1000,
        region_id="em-complement-literal-collision-table",
        method=MethodCapability(
            genome="published Earth-Moon p:q orbits vs #1000 families A, B, D",
            corrector="none (propagation of published states)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=40,
    )
    rows: list[dict[str, Any]] = []
    for r in TABLE3_ROWS:
        s6 = table3_seed_state(r)
        s4 = np.array([s6[0], s6[1], s6[3], s6[4]])
        t = topo(s4, r.period)
        rows.append(
            {
                "source": "Casoliva et al. 2010 Table 3",
                "orbit": r.designation,
                "printed_pq": f"{r.p}-{r.q}",
                "C": r.c_j,
                "T": r.period,
                **t,
            }
        )
    cat = {x["id"]: x for x in yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())}
    for rid in [k for k in cat if k.startswith("vaquero-")]:
        c = cat[rid]["orbit_elements"]["cr3bp"]
        s6 = np.array(c["state_nd"])
        t = topo(np.array([s6[0], s6[1], s6[3], s6[4]]), c["period_nd"])
        rows.append(
            {
                "source": "Vaquero 2013 (catalogue row)",
                "orbit": rid,
                "C": c["jacobi_constant"],
                "T": c["period_nd"],
                **t,
            }
        )
    seeds = json.loads((REPO / "data" / "997_lineage" / "seeds.json").read_text())
    for s in seeds.get("F2", []) + seeds.get("F3", []):
        t = topo(np.array([s["x0"], 0.0, 0.0, s["ydot0"]]), 2 * s["t_half"])
        rows.append(
            {
                "source": "Newton 1959 Table 1 (registry mu, #997 seed)",
                "orbit": s["id"],
                "C": base.omega_eff(s["x0"], 0.0) - s["ydot0"] ** 2,
                "T": 2 * s["t_half"],
                **t,
            }
        )
    for d, pq, cj, tt in [
        ("12a", (1, 2), None, None),
        ("21a", (2, 1), None, None),
        ("23a/b", (2, 3), None, None),
        ("32a/b", (3, 2), None, None),
        ("52a", (5, 2), 1.0461882704974470, 12.5651492405106922),
        ("54a/b", (5, 4), None, None),
        ("73a", (7, 3), 0.8957501590757784, 18.8492803402344329),
    ]:
        rows.append(
            {
                "source": "Casoliva et al. 2008 Table 2 (seeds)",
                "orbit": d,
                "pq": list(pq),
                "sense": None,
                "C": cj,
                "T": tt,
                "mu_note": "mu = 1e-6 seed, not Earth-Moon",
            }
        )
    rows.append(
        {
            "source": "Liang, Xu & Xu 2017",
            "orbit": "5:2 PLPO",
            "pq": [5, 2],
            "sense": None,
            "C": 3.0996,
            "T": None,
            "note": "apoapsis about 0.82 L, no lunar pass",
        }
    )
    rows.append(
        {
            "source": "Liang, Xu & Xu 2017",
            "orbit": "7:3 PLPO",
            "pq": [7, 3],
            "sense": None,
            "C": 3.1858,
            "T": None,
            "note": "apoapsis about 0.75 L, no lunar pass",
        }
    )
    rows.append(
        {
            "source": "Liang, Xu, Peng & Xu 2020",
            "orbit": "2:1 resonant cycler orbit",
            "pq": [2, 1],
            "sense": None,
            "C": 2.0934,
            "T": None,
        }
    )
    rows.append(
        {
            "source": "Hoelker & Winston 1968",
            "orbit": "all 39 captions",
            "pq": None,
            "why_no_pq": "no n* of 5/2, 5/3 or 7/3 printed; mu = 1/80",
        }
    )
    rows.append(
        {
            "source": "Genova & Aldrin 2015 (AAS 15-794)",
            "orbit": "5 cyclers",
            "pq": None,
            "why_no_pq": "2:1, 3:1 and the Arenstorf 4-leaf; no C or state",
        }
    )
    for r in rows:
        r["verdict"] = verdict(r)
    OUT.joinpath("literal_check.json").write_text(json.dumps(rows, indent=1) + "\n")
    lines = [
        "| source | orbit | p:q (sense) | C | T | perigee alt (km) | periselene alt (km) | "
        "vs A (5:2 pro) | vs B (5:3 ret) | vs D (7:3 pro) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        pq = r.get("pq")
        sense = {1: "pro", -1: "ret", None: "?"}[r.get("sense")]
        pqs = f"{pq[0]}:{pq[1]} ({sense})" if pq else "-"
        fmt = lambda v, f: f.format(v) if isinstance(v, int | float) else "-"  # noqa: E731
        lines.append(
            f"| {r['source']} | {r['orbit']} | {pqs} | {fmt(r.get('C'), '{:.4f}')} | "
            f"{fmt(r.get('T'), '{:.3f}')} | {fmt(r.get('perigee_alt_km'), '{:.0f}')} | "
            f"{fmt(r.get('periselene_alt_km'), '{:.0f}')} | {r['verdict']['A']} | "
            f"{r['verdict']['B']} | {r['verdict']['D']} |"
        )
    OUT.joinpath("literal_check.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
