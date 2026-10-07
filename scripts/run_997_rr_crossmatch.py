"""#997: cross-match the cycler-class candidates against the Restrepo & Russell (2018) Earth-Moon
planar axisymmetric periodic-orbit database. Positive controls: the catalogue rows found on the
#997 families (casoliva-7-3a, casoliva-2-1b, vaquero-21-*).

Database: Restrepo, R.L. & Russell, R.P., Celest. Mech. Dyn. Astr. (2018) 130:49,
doi 10.1007/s10569-018-9844-6 (Apache 2.0). Read in place from a local copy, outside the repo
(default path below; --db to override). Columns per the database's columnNames_Details.txt:
col 3 = J - 3, col 6 = mu, cols 14-19 = state at the first perpendicular crossing in the
MOON-centred frame (primary at x = -1), col 20 = period, cols 30-31 = second perpendicular
crossing (x, vy) in the local-search files only. J is checked to equal the project's C
(no mu(1 - mu) term) after the shift X = x + 1 - mu.

Match rule (stated before reading the results): an RR record matches a target orbit if
|J - C| < 5e-3, |T_RR - T| / T < 5e-3, and one of the RR record's perpendicular crossings
(x, vy) lies within 5e-3 (in both x and vy, barycentric) of one of the target's perpendicular
crossings. The tolerance absorbs the mass-ratio difference: RR global searches use
mu = 0.0121437 (GM 398600/4900), the local searches 0.0121506; the registry is 0.0121505844.
Output: data/997_lineage/rr_crossmatch.json.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "997_lineage"
DEFAULT_DB = (
    Path.home()
    / "dev/references/restrepo-russell-2018/POs_database_v2019-12-28"
    / "POs_database_all-2017-12-14/Earth_Sys/Moon"
)
TOL_C, TOL_T, TOL_X = 5e-3, 5e-3, 5e-3


def load_rr(db: Path) -> tuple[np.ndarray, list[str]]:
    rows, names = [], []
    for f in sorted(glob.glob(str(db / "*" / "gridPOrun.out"))):
        a = np.loadtxt(f, comments="%")
        if a.shape[1] < 31:  # global files carry no second crossing
            a = np.hstack([a, np.full((a.shape[0], 31 - a.shape[1]), np.nan)])
        rows.append(a)
        names += [Path(f).parent.name] * a.shape[0]
    return np.vstack(rows), names


def perp_crossings(x0: float, yd0: float, period: float, mu: float) -> list[tuple[float, float]]:
    """All perpendicular x-axis crossings (x, ydot) of a planar symmetric orbit over one period."""

    def ev(t: float, y: np.ndarray, m: float) -> float:
        return float(y[1])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, period),
        np.array([x0, 0.0, 0.0, 0.0, yd0, 0.0]),
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=ev,
    )
    out = [(x0, yd0)]
    for y in sol.y_events[0]:
        if abs(y[3]) < 1e-6:
            out.append((float(y[0]), float(y[4])))
    return out


def match(
    rr: np.ndarray, names: list[str], c: float, t: float, crossings: list[tuple[float, float]]
) -> dict[str, Any]:
    j = rr[:, 2] + 3.0
    mu_rr = rr[:, 5]
    sel = np.where((np.abs(j - c) < TOL_C) & (np.abs(rr[:, 19] - t) / t < TOL_T))[0]
    hits = []
    for i in sel:
        rx = [(rr[i, 13] + 1 - mu_rr[i], rr[i, 17])]
        if not math.isnan(rr[i, 29]):
            rx.append((rr[i, 29] + 1 - mu_rr[i], rr[i, 30]))
        d = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a in rx for b in crossings)
        if d < TOL_X:
            hits.append(
                {
                    "file": names[i],
                    "row": int(i),
                    "J": float(j[i]),
                    "T": float(rr[i, 19]),
                    "mu": float(mu_rr[i]),
                    "dmax_crossing": float(d),
                    "b_h": float(rr[i, 9]),
                    "b_v": float(rr[i, 10]),
                    "Ncross": float(rr[i, 25]),
                    "Mrot": float(rr[i, 27]),
                    "Nrot": float(rr[i, 28]),
                }
            )
    hits.sort(key=lambda h: h["dmax_crossing"])
    return {"n_prefilter": len(sel), "n_match": len(hits), "best": hits[:3]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    args = ap.parse_args()
    preflight_search(
        task_no=997,
        region_id="em-both-primary-lineage-rr2018-crossmatch",
        method=MethodCapability(
            genome="#997 candidates vs Restrepo-Russell 2018 Earth-Moon database",
            corrector="none (lookup)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    rr, names = load_rr(args.db)
    j_all = rr[:, 2] + 3.0
    print(f"RR Earth-Moon: {rr.shape[0]} records, J range {j_all.min():.4f} .. {j_all.max():.4f}")
    out: dict[str, Any] = {
        "db_records": int(rr.shape[0]),
        "db_J_range": [float(j_all.min()), float(j_all.max())],
        "db_mu_values": sorted({float(m) for m in rr[:, 5]}),
        "tolerances": {"C": TOL_C, "T_rel": TOL_T, "crossing": TOL_X},
        "controls": {},
        "candidates": {},
    }
    gate = json.loads((DATA / "gate.json").read_text())
    for g in gate:
        if not (g["matches"] and g["matches"][0]["crossing_match"]):
            continue
        cr = [(p[1], p[2]) for p in g["perpendicular_crossings"]]
        res = match(rr, names, g["C"], g["T"], cr)
        res["C"], res["T"] = g["C"], g["T"]
        res["inside_db_J_range"] = bool(j_all.min() <= g["C"] <= j_all.max())
        out["controls"][g["id"]] = res
        print(f"control {g['id']:40s} C={g['C']:.4f} -> {res['n_match']} RR matches")
    summary = json.loads((DATA / "summary.json").read_text())
    for fam, v in summary.items():
        rows = []
        for cnd in v["candidates"]:
            cr = perp_crossings(cnd["x0"], cnd["ydot0"], cnd["T"], 0.01215058439469525)
            res = match(rr, names, cnd["C"], cnd["T"], cr)
            rows.append({"seed": cnd["seed"], "k": cnd["k"], "C": cnd["C"], "T": cnd["T"], **res})
        if rows:
            n_in = sum(1 for r in rows if j_all.min() <= r["C"] <= j_all.max())
            n_hit = sum(1 for r in rows if r["n_match"] > 0)
            out["candidates"][fam] = {
                "n": len(rows),
                "n_inside_db_J_range": n_in,
                "n_matched": n_hit,
                "members": rows,
            }
            print(f"{fam}: {len(rows)} candidates, {n_in} inside RR J range, {n_hit} matched")
    # Family level: every computed member (not only candidates) inside RR's J range. RR samples a
    # family at isolated grid points, so a family is "in RR" if some of its members match.
    fam_level: dict[str, Any] = {}
    for f in sorted(glob.glob(str(DATA / "F*_[pm].jsonl"))):
        hits_c = []
        n_in = 0
        for line in Path(f).read_text().splitlines():
            m = json.loads(line)
            if not (j_all.min() <= m["C"] <= j_all.max()):
                continue
            n_in += 1
            cr = perp_crossings(m["x0"], m["ydot0"], m["T"], m["mu"])
            if match(rr, names, m["C"], m["T"], cr)["n_match"] > 0:
                hits_c.append(round(m["C"], 4))
        fam_level[Path(f).name] = {"n_inside_J_range": n_in, "matched_C": hits_c}
        print(f"{Path(f).name}: {n_in} members inside RR J range, {len(hits_c)} matched")
    out["family_level"] = fam_level
    (DATA / "rr_crossmatch.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
