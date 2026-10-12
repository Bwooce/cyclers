"""#1055 amendment 14 check (ideal model): the V_inf DIRECTION of each continuous-gravity branch
against the chain tool's patched-conic member, per node, in and out.

Amendment 10 (V_inf magnitude, r_p / sigma, turn, epoch) cannot see a move along a crank of a
resonant leg when the neighbouring Lambert leg arrives almost along the moon's velocity: the turn
then hardly changes with the crank. This script reports, for every ideal-model member with a
#1025 / #1055 sigma file, the angle between the branch's V_inf in / out vectors (the two-body
asymptotes of the node state at that sigma, moon GM scaled by sigma) and the chain tool's (the same
asymptotes of the seed at sigma = 1), and for the resonant legs the pump and crank about the moon's
velocity (crank 0 = radially outward in the moons' plane, 90 = towards +z).

    uv run python scripts/check_1055_vinf_direction.py            # every data/1025_sigma member
    uv run python scripts/check_1055_vinf_direction.py gc4-19     # one member

Sources per member: data/1025_sigma/<m>.json (its sigma-0.02 point; the identity stage does not
store its z) and, for gc4-19, data/1055_gc4-19/limit.json (0.02, 0.01, 0.005, 0.0025). An O(sigma)
offset is about 0.03-0.2 deg at sigma 0.02; a crank-sibling stays at degrees as sigma falls.
Read-only; writes data/1055_vinf_direction.json.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.satellites import SATELLITES

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B = _load("run_1025_sigma_batch", ROOT / "scripts" / "run_1025_sigma_batch.py")
RB = _load("run_968_rungb", ROOT / "scripts" / "run_968_rungb.py")


def _ang(a: Any, b: Any) -> float:
    c = float(a @ b) / (float(np.linalg.norm(a)) * float(np.linalg.norm(b)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def vinf_vectors(moons: list[str], z: Any, sigma: float) -> list[tuple[Any, Any]]:
    out = []
    for k, moon in enumerate(moons):
        r, v = np.asarray(z[7 * k : 7 * k + 3]), np.asarray(z[7 * k + 3 : 7 * k + 6])
        mu = sigma * SATELLITES[moon].mu_km3_s2
        out.append((RB.asymptote(mu, r, v, "in"), RB.asymptote(mu, r, v, "out")))
    return out


def pump_crank(p: Any, moon: str, t: float, vec: Any) -> tuple[float, float]:
    _, vm = p.circ.state(moon, t)
    vh = np.asarray(vm) / float(np.linalg.norm(vm))
    zh = np.array([0.0, 0.0, 1.0])
    rh = np.cross(vh, zh)
    pump = _ang(vec, vh)
    crank = math.degrees(math.atan2(float(vec @ zh), float(vec @ rh)))
    return pump, crank


def member(name: str) -> dict[str, Any]:
    p = B.build(name)
    moons = list(p.moons)
    seed = np.asarray(p.seed)
    chain = vinf_vectors(moons, seed, 1.0)
    pts: list[tuple[float, Any]] = []
    src = ROOT / "data" / "1025_sigma" / f"{name}.json"
    if src.exists():
        qs = json.loads(src.read_text())["points"]
        if qs:  # gc4-19's #1025 file is empty (no start)
            pts.append((float(qs[0]["sigma"]), np.asarray(qs[0]["z"])))
    lim = ROOT / "data" / f"1055_{name}" / "limit.json"
    if lim.exists():  # the #1055 corrected branch going down (gc4-19)
        pts = [
            (float(q["sigma"]), np.asarray(q["z"])) for q in json.loads(lim.read_text())["points"]
        ]
    rows = []
    for sigma, z in pts:
        vv = vinf_vectors(moons, z, sigma)
        nodes = []
        for k, moon in enumerate(moons):
            t = float(z[7 * k + 6])
            nodes.append(
                {
                    "moon": moon,
                    "angle_in_deg": _ang(vv[k][0], chain[k][0]),
                    "angle_out_deg": _ang(vv[k][1], chain[k][1]),
                    "pump_crank_in": pump_crank(p, moon, t, vv[k][0]),
                    "pump_crank_out": pump_crank(p, moon, t, vv[k][1]),
                }
            )
        rows.append({"sigma": sigma, "nodes": nodes})
    chain_pc = [
        {
            "moon": moon,
            "pump_crank_in": pump_crank(p, moon, float(seed[7 * k + 6]), chain[k][0]),
            "pump_crank_out": pump_crank(p, moon, float(seed[7 * k + 6]), chain[k][1]),
        }
        for k, moon in enumerate(moons)
    ]
    out: dict[str, Any] = {"member": name, "chain": chain_pc, "points": rows}
    if len(rows) >= 3:  # Richardson limit of the angles from the two smallest sigmas
        a, b = rows[-1]["nodes"], rows[-2]["nodes"]
        out["limit_angles_deg"] = [
            (2 * x["angle_in_deg"] - y["angle_in_deg"], 2 * x["angle_out_deg"] - y["angle_out_deg"])
            for x, y in zip(a, b, strict=True)
        ]
    out["max_angle_at_first_point_deg"] = (
        max(max(n["angle_in_deg"], n["angle_out_deg"]) for n in rows[0]["nodes"]) if rows else None
    )
    return out


def _mx(a: list[float] | None) -> float | None:
    return round(max(a), 4) if a else None


def down(name: str) -> dict[str, Any]:
    """Sec. 14.2: the #1025 solver going down (0.01, 0.005, 0.0025) from the committed sigma-0.02
    point, the angles at each converged point, the Richardson limit, ratios and the verdict."""
    p = B.build(name)
    moons = list(p.moons)
    chain = vinf_vectors(moons, np.asarray(p.seed), 1.0)
    q = json.loads((ROOT / "data" / "1025_sigma" / f"{name}.json").read_text())["points"][0]
    z, s_prev = np.asarray(q["z"]), float(q["sigma"])

    def angles(zz: Any, s: float) -> list[float]:
        vv = vinf_vectors(moons, zz, s)
        return [_ang(vv[k][j], chain[k][j]) for k in range(len(moons)) for j in (0, 1)]

    pts = [{"sigma": s_prev, "converged": True, "angles": angles(z, s_prev)}]
    for sg in (0.01, 0.005, 0.0025):
        p.set_sigma(sg)
        zn, info = B.solve(p, B.scale_offsets(p, z, sg / s_prev), 40, 60, False)
        ok = bool(B.converged(info))
        pts.append(
            {
                "sigma": sg,
                "converged": ok,
                "angles": angles(zn, sg) if ok else None,
                "z": zn.tolist() if ok else None,
            }
        )
        if not ok:
            break
        z, s_prev = zn, sg
    good = [x for x in pts if x["converged"]]
    out: dict[str, Any] = {"member": name, "moons": moons, "points": pts}
    if len(good) < 3:
        out["verdict"] = "UNDECIDED"
        return out
    a, b = good[-1]["angles"], good[-2]["angles"]
    lim = [2 * x - y for x, y in zip(a, b, strict=True)]
    counted = [i for i, x in enumerate(good[0]["angles"]) if x >= 0.02]
    ratios = {
        i: [good[n]["angles"][i] / good[n - 1]["angles"][i] for n in range(1, len(good))]
        for i in counted
    }
    lim_ok = all(x < 0.1 for x in lim)
    band_ok = all(0.35 <= r <= 0.65 for rs in ratios.values() for r in rs)
    if lim_ok and band_ok:
        verdict = "PASS"
    elif lim_ok:
        verdict = "IRREGULAR"
    elif any(lim[i] >= 0.1 and ratios.get(i, [1.0])[-1] > 0.65 for i in range(len(lim))):
        verdict = "CRANK-SIBLING"
    else:
        verdict = "FAIL (limit >= 0.1 deg, last ratio <= 0.65)"
    out.update(
        limit_deg=lim,
        ratios={str(k): v for k, v in ratios.items()},
        sigmas_for_limit=[good[-2]["sigma"], good[-1]["sigma"]],
        verdict=verdict,
    )
    return out


def main() -> None:
    if sys.argv[1:2] == ["--down"]:
        odir = ROOT / "data" / "1055_vinf_limit"
        odir.mkdir(exist_ok=True)
        for name in sys.argv[2:]:
            r = down(name)
            (odir / f"{name}.json").write_text(json.dumps(r, indent=1))
            print(
                f"{name}: {r['verdict']} limit max {max(r.get('limit_deg', [float('nan')])):.4f} "
                f"max angle per sigma "
                f"{[(x['sigma'], _mx(x['angles'])) for x in r['points']]}",
                flush=True,
            )
        return
    names = sys.argv[1:] or sorted(
        f.stem
        for f in (ROOT / "data" / "1025_sigma").glob("*.json")
        if f.stem not in ("analysis", "summary")
    )
    res = []
    for name in names:
        r = member(name)
        res.append(r)
        print(
            f"{name}: max V_inf-direction angle at the first point"
            f" = {r['max_angle_at_first_point_deg']}"
        )
        for row in r["points"]:
            print(
                f"  sigma {row['sigma']:.4g}: "
                + " | ".join(
                    f"{n['moon'][0]} {n['angle_in_deg']:.3f}/{n['angle_out_deg']:.3f}"
                    for n in row["nodes"]
                )
            )
        if "limit_angles_deg" in r:
            print(
                "  Richardson limit:",
                [(round(a, 3), round(b, 3)) for a, b in r["limit_angles_deg"]],
            )
    if len(names) > 1:
        (ROOT / "data" / "1055_vinf_direction.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
