"""#942 positive control: Hollister & Menning 1970 Table 3 with the two-working-body corrector.

PRE-REGISTERED (written 2026-10-05 before any comparison was run):

Model: fixed mean-element inclined-elliptic heliocentric ephemeris
(``MeanElementSystem``); flyby floor for the turn: H&M's own acceptance is
"beyond 1.1 planet radii" (p.1193); Rmin is "minimum distance to planet centre"
(Table 3 footnote), compared with ``rmin_radii``.

Free directions: each full-revolution V-infinity is chosen to minimise the
largest turn of the flybys that bracket it (minimax). This is OUR hypothesis
for the convention of H&M's ref. 16 (not held); it is tested, not tuned.

Symmetric-return branch: chosen by the junction-residual sum at the printed
dates (``pick_branches``), never by theta or Rmin.

Stage F (forward, printed dates): per encounter, |V_inf| arrival and departure
vs printed V_r. PASS if >= 90 % of encounters are within 0.006 EMOS.

Stage C (corrected): ``correct_dates`` from (a) the printed dates, (b) the
printed dates plus uniform noise of +/-7 d (seed 942), (c) for orbit 1, the
Table 2 orbit-1H approximate dates. PASS per orbit if it converges
(max |residual| < 1e-9 km/s) and >= 90 % of encounters are within: date 3 d,
V_r 0.005 EMOS, theta 3.0 deg, Rmin 10 %. Overall PASS: all 15 orbits pass
from seed (a); seeds (b)/(c) report the fraction that return to the same
solution (dates within 1 d of the seed-(a) solution).

AMENDMENT 2026-10-05 (after stage-1 result 0/15; stage 1 stays reported):
``--periodic`` runs the same test on the model H&M state they used, "assuming
exact periodicity of the solar system" over 16 yr (p.1194): Earth period
5844/16 d, Venus 5844/26 d, mean longitudes exact at the anchor JD 2443363
(the middle of orbit 1's span, fixed before the rerun). Tolerances unchanged.
Diagnosis behind the amendment: stage 1 dates slip progressively (Venus
full-revolution steps of 224.70 d against the printed 225 d).

AMENDMENT 2 2026-10-05 (after the periodic rerun; earlier results stay reported):
(a) the minimax leaves the inner flyby of a block with two free directions
undetermined (orbits 4 and 8 missed only there); ties are now broken by raising
the smallest ratio with the largest held fixed (``_balance``), which cannot
change a gate verdict; (b) date comparisons use the print-error fixes of
``hollister_menning_1970._PRINT_ERROR_FIXES``; (c) a printed row whose own
(V_r, theta, Rmin) triple disagrees with the flyby formula by > 10 % is
"source-inconsistent"; the post-hoc fraction excludes those rows.

Usage: ``uv run python scripts/run_942_hm_control.py --out <dir> [--periodic]``
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.hollister_menning_1970 import (
    Row,
    build_cycle,
    load_table3,
    pick_branches,
)
from cyclerfinder.search.two_working_body import (
    EMOS_KMS,
    MeanElementSystem,
    correct_dates,
    cycle_flybys,
    date_residual,
    encounter_self_consistency,
    eval_lambert_legs,
    rmin_radii,
)

DAY = SECONDS_PER_DAY
F_TOL_EMOS = 0.006
C_TOL = {"date": 3.0, "vr": 0.005, "theta": 3.0, "rmin_rel": 0.10}
PASS_FRAC = 0.90
#: Table 2, orbit 1H (inclined-elliptic approximation, p.1195).
TABLE2_1H = [
    441,
    806,
    971,
    1196,
    1421,
    1592,
    1957,
    2125,
    2350,
    2575,
    2797,
    3163,
    3316,
    3541,
    3765,
    3935,
    4300,
    4471,
    4696,
    4921,
    5077,
    5442,
    5664,
    5889,
    6114,
    6285,
]


def log(msg: str) -> None:
    print(f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}", flush=True)


def match_rows(rows: list[Row], flybys: list, period_d: float) -> list[tuple[Row, object]]:
    """Pair each printed row (except the closing repeat) with the flyby nearest in date."""
    out = []
    for r in rows[:-1]:
        best = None
        for f in flybys:
            d = f.t_s / DAY
            dd = min(abs(d - r.date), abs(d - period_d - r.date), abs(d + period_d - r.date))
            if f.body == r.planet and (best is None or dd < best[0]):
                best = (dd, f)
        out.append((r, best))
    return out


def compare(system, rows, flybys, period_d):
    res = []
    for r, best in match_rows(rows, flybys, period_d):
        if best is None:
            res.append({"planet": r.planet, "date": r.date, "missing": True})
            continue
        dd, f = best
        body = system.body(f.body)
        th = math.radians(f.turn_deg)
        res.append(
            {
                "planet": r.planet,
                "date": r.date,
                "date_ours": f.t_s / DAY,
                "ddate": dd,
                "vr": r.vr_emos,
                "vr_ours": f.vinf_kms / EMOS_KMS,
                "theta": r.theta_deg,
                "theta_ours": f.turn_deg,
                "rmin": r.rmin_radii,
                "rmin_ours": rmin_radii(body, f.vinf_kms, th),
            }
        )
    return res


def encounter_pass(e: dict) -> bool:
    if e.get("missing"):
        return False
    return (
        e["ddate"] <= C_TOL["date"]
        and abs(e["vr_ours"] - e["vr"]) <= C_TOL["vr"]
        and abs(e["theta_ours"] - e["theta"]) <= C_TOL["theta"]
        and abs(e["rmin_ours"] / e["rmin"] - 1.0) <= C_TOL["rmin_rel"]
    )


def source_consistent(system, e: dict) -> bool:
    """Is the PRINTED (V_r, theta, Rmin) triple consistent with the flyby formula
    ``r_p = mu/v^2 (1/sin(theta/2) - 1)`` and registry constants, within 10 %?
    Independent of our solver: it tests the source row against itself."""
    body = system.body(e["planet"])
    v = e["vr"] * EMOS_KMS
    rp = rmin_radii(body, v, math.radians(e["theta"]))
    return bool(abs(rp / e["rmin"] - 1.0) <= 0.10)


def forward(system, orbit, rows, branches):
    cyc, x = build_cycle(system, orbit, rows, branches)
    legs = eval_lambert_legs(system, cyc, x)
    assert legs is not None
    r = date_residual(system, cyc, x)
    # arrival/departure magnitudes vs printed V_r of the encounter at that date
    by_date = {round(rw.date): rw for rw in rows}
    errs = []
    for ev in legs:
        for t, v in ((ev.t_dep, ev.vinf_dep), (ev.t_arr, ev.vinf_arr)):
            d = t / DAY
            rw = by_date.get(round(d)) or by_date.get(round(d - 5844.0))
            if rw is None:
                continue
            errs.append(abs(float(np.linalg.norm(v)) / EMOS_KMS - rw.vr_emos))
    return {"max_junction_residual_kms": float(np.max(np.abs(r))), "vr_errs": errs}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--orbits", type=str, default="1-15")
    ap.add_argument("--periodic", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    lo, hi = (int(v) for v in args.orbits.split("-"))
    preflight_search(
        task_no=942,
        region_id="hollister-menning-1970-table3-positive-control",
        method=MethodCapability(
            genome="Hollister-Menning Earth-Venus periodic swing-by chains (Table 3 structures)",
            corrector="two_working_body.correct_dates (H&M date residual, least squares)",
            capability_tags=frozenset({"ballistic", "patched-conic", "3d", "inclined-elliptic"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=hi - lo + 1,
    )
    system = (
        MeanElementSystem(periods_days={"E": 5844.0 / 16, "V": 5844.0 / 26}, anchor_jd=2443363.0)
        if args.periodic
        else MeanElementSystem()
    )
    table = load_table3()
    rng = np.random.default_rng(942)
    summary = {}
    t_start = time.time()
    orbits = list(range(lo, hi + 1))
    for n_done, orbit in enumerate(orbits):
        rows = table[orbit]
        branches, _ = pick_branches(system, orbit, rows)
        fw = forward(system, orbit, rows, branches)
        f_frac = float(np.mean(np.array(fw["vr_errs"]) <= F_TOL_EMOS))
        cyc, x0 = build_cycle(system, orbit, rows, branches)
        sol = correct_dates(system, cyc, x0)
        # Compared at the least-squares point even when it is not an exact zero
        # (reported separately; the pass needs ``converged``).
        fl = cycle_flybys(system, cyc, sol.x)
        cmp_a = compare(system, rows, fl, 5844.0) if fl else []
        for e in cmp_a:
            e["source_consistent"] = source_consistent(system, e)
        frac_a = float(np.mean([encounter_pass(e) for e in cmp_a])) if cmp_a else 0.0
        cons = [e for e in cmp_a if e.get("source_consistent")]
        frac_cons = float(np.mean([encounter_pass(e) for e in cons])) if cons else 0.0
        hm_sum_abs_emos = float(np.sum(np.abs(sol.residual))) / EMOS_KMS
        miss = encounter_self_consistency(system, cyc, sol.x)
        # seed (b): perturbed
        x_b = x0 + rng.uniform(-7.0, 7.0, size=x0.size) * DAY
        sol_b = correct_dates(system, cyc, x_b)
        same_b = bool(
            sol_b.converged and sol.converged and np.max(np.abs(sol_b.x - sol.x)) / DAY < 1.0
        )
        rec = {
            "orbit": orbit,
            "sy_branches": list(branches),
            "forward": {
                "max_junction_residual_kms": fw["max_junction_residual_kms"],
                "vr_frac_within": f_frac,
                "vr_err_max_emos": float(max(fw["vr_errs"])),
            },
            "forward_pass": f_frac >= PASS_FRAC,
            "converged": sol.converged,
            "max_residual_kms": sol.max_abs_residual_kms,
            "max_date_shift_d": float(np.max(np.abs(sol.x - x0)) / DAY),
            "self_consistency_miss_km": miss,
            "frac_encounters_within_tol": frac_a,
            "control_pass": bool(sol.converged and frac_a >= PASS_FRAC),
            "hm_sum_abs_residual_emos": hm_sum_abs_emos,
            "hm_tolerance_met": hm_sum_abs_emos <= 0.005,
            "n_source_inconsistent_rows": len(cmp_a) - len(cons),
            "posthoc_frac_within_tol_source_consistent_rows": frac_cons,
            "perturbed_seed_converged": sol_b.converged,
            "perturbed_seed_same_solution": same_b,
            "encounters": cmp_a,
        }
        if orbit == 1:
            rows_c = [
                Row(r.planet, float(d), r.vr_emos, r.theta_deg, r.rmin_radii)
                for r, d in zip(rows, TABLE2_1H, strict=True)
            ]
            cyc_c, x_c = build_cycle(system, 1, rows_c, branches)
            sol_c = correct_dates(system, cyc_c, x_c)
            rec["table2_seed_converged"] = sol_c.converged
            rec["table2_seed_same_solution"] = bool(
                sol_c.converged and sol.converged and np.max(np.abs(sol_c.x - sol.x)) / DAY < 1.0
            )
        summary[orbit] = rec
        (args.out / f"orbit_{orbit:02d}.json").write_text(json.dumps(rec, indent=1, default=float))
        el = time.time() - t_start
        eta = el / (n_done + 1) * (len(orbits) - n_done - 1)
        log(
            f"orbit {orbit}: fwd_frac={f_frac:.2f} conv={sol.converged} "
            f"res={sol.max_abs_residual_kms:.1e} "
            f"shift={rec['max_date_shift_d']:.2f}d match={frac_a:.2f} pass={rec['control_pass']} "
            f"pert_same={same_b} miss={miss:.1e}km hm_sum={hm_sum_abs_emos:.4f}EMOS "
            f"posthoc={frac_cons:.2f}({len(cmp_a) - len(cons)} incons) "
            f"[{n_done + 1}/{len(orbits)} eta {eta:.0f}s]"
        )
    (args.out / "summary.json").write_text(
        json.dumps(
            {k: {kk: vv for kk, vv in v.items() if kk != "encounters"} for k, v in summary.items()},
            indent=1,
            default=float,
        )
    )
    n_pass = sum(v["control_pass"] for v in summary.values())
    log(f"DONE control_pass {n_pass}/{len(summary)}")


if __name__ == "__main__":
    main()
