"""#1007 EXPLORATORY pre-check of the DRAFT same-phase drift metric (sec. 8.4), awaiting the owner.

Not a spec or catalogue change. For each chain: D'_i = max_j |t_(i+L)j - t_ij - L P_mean| and
W'_i = max_j | |V_inf,dep|_(i+L)j - |V_inf,dep|_ij |, i = 0 .. n-L-1; band D' <= 5 % of the chain's
shortest Lambert leg, W' <= 10 % of its smallest |V_inf,dep|; no growth (second-half max <= 2 x
first-half max, floors 0.5 d / 0.05 km/s). Junction-only continuity is reported beside it.

Chains: gc-1 and gc-2 (existing 10-cycle jup365 chains, L = 4); H&M 1970 orbit 1 from its printed
dates over two H&M periods (10 cycles, Standish, direct LM; L = 5); the broken-date negative (gc-1
epoch 0, 2nd Lambert start + 0.05 d).

Usage: uv run python scripts/v2_samephase_1007.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.search.hollister_menning_1970 import load_table3
from cyclerfinder.search.two_working_body import Cycle, LambertLeg, eval_lambert_legs

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
OUT = REPO / "data" / "1007_v2" / "samephase_exploratory.json"
LABEL = "exploratory, under the draft metric, awaiting the owner"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


V2 = _load("v2_chain_1007", REPO / "scripts" / "v2_chain_1007.py")
CHAIN = V2.CHAIN


def legs_data(
    cell: str, key: str, n: int, mode: str, x0: float, y: list[float], real: Any
) -> tuple[np.ndarray, list[float], list[float], float]:
    """(Lambert start dates in days, |V_inf,dep| per Lambert leg, flight times in days,
    junction-only continuity in km/s)."""
    circ, a, b = CHAIN.ENUM.cell_system(cell)
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    legs_chain = one.legs * n
    nl = sum(isinstance(lg, LambertLeg) for lg in one.legs)
    yy = np.asarray(y, dtype=float)
    x = np.concatenate([[x0], yy[: n * nl - 1]])
    if mode == "ramp":
        cyc = Cycle(legs_chain, float(yy[n * nl - 1]) * DAY)
        lam = eval_lambert_legs(real, cyc, x * DAY)
        assert lam is not None
        res = CHAIN.date_chain_residual(real, legs_chain, x0 * DAY, yy)
        return (
            x,
            [float(np.linalg.norm(lg.vinf_dep)) for lg in lam],
            [(lg.t_arr - lg.t_dep) / DAY for lg in lam],
            float(np.max(np.abs(res))),
        )
    sysm = CHAIN.Blend(circ, real, np.eye(3), 0.0, 0.0, 1.0)
    ev = CHAIN.chain_eval(sysm, legs_chain, x0 * DAY, yy)
    assert ev is not None
    mags, q = [], 0
    for blk in ev.block_flybys:
        q += 3 * (len(blk) - 1)
        mags.append(abs(float(ev.residual[q])))
        q += 1
    vdep, tof = [], []
    for frm, t0, v0, _to, t1 in ev.segments[: n * nl]:
        _, w = real.state(frm, t0)
        vdep.append(float(np.linalg.norm(np.asarray(v0) - w)))
        tof.append((t1 - t0) / DAY)
    return x, vdep, tof, max(mags)


def same_phase(label: str, x: np.ndarray, vdep: list[float], tof: list[float], cont: float,
               n: int, nl: int, period_d: float, big_l: int, start: int = 0,
               quiet: bool = False) -> dict[str, Any]:  # fmt: skip
    """``start`` = 0: the draft as written (sec. 8.4 option (a)); 1: option (b), the open
    chain's start cycle skipped."""
    p_mean = period_d / n
    m = n - big_l
    d = [
        max(abs(x[(i + big_l) * nl + j] - x[i * nl + j] - big_l * p_mean) for j in range(nl))
        for i in range(start, m)
    ]
    w = [
        max(abs(vdep[(i + big_l) * nl + j] - vdep[i * nl + j]) for j in range(nl))
        for i in range(start, m)
    ]
    band_d, band_w = 0.05 * min(tof), 0.10 * min(vdep)
    h = max(1, len(d) // 2)

    def growth(v: list[float], floor: float) -> tuple[bool, float]:
        first, second = max(v[:h]), max(v[h:]) if v[h:] else 0.0
        r = second / first if first > 0 else (float("inf") if second > 0 else 0.0)
        return second <= floor or second <= 2.0 * first, r

    gd_ok, gd = growth(d, 0.5)
    gw_ok, gw = growth(w, 0.05)
    ok = max(d) <= band_d and max(w) <= band_w and gd_ok and gw_ok
    out = {
        "label": label, "status": LABEL, "L": big_l, "cycles": n, "comparisons": len(d),
        "start_cycle": start,
        "junction_continuity_kms": cont, "Dprime_d": d, "Wprime_kms": w,
        "band_D_d": band_d, "band_W_kms": band_w, "growth_D": gd, "growth_W": gw,
        "drift_pass": bool(ok), "continuity_pass": bool(cont < 1e-6),
    }  # fmt: skip
    if not quiet:
        print(
            f"{label}: L={big_l} n={n} cont={cont:.1e} D'max={max(d):.3f}/{band_d:.3f} d "
            f"W'max={max(w):.4f}/{band_w:.3f} km/s growth D {gd:.2f} W {gw:.2f} -> drift "
            f"{'PASS' if ok else 'FAIL'}, continuity {'PASS' if cont < 1e-6 else 'FAIL'} [{LABEL}]",
            flush=True,
        )
    return out


def ev_rows() -> dict[str, Any]:
    """Sec. 8.6: the ev rows' new 10-cycle Standish chains, options (a) and (b)."""
    real = CHAIN.real_ephemeris("ev", "mean")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    res: dict[str, Any] = {"status": LABEL}
    for row, sub in (("ev-A", "ev10_evA"), ("ev-B", "ev10_evB"), ("ev-C", "ev10_evC")):
        key = wb[row]["key"]
        circ, a, b = CHAIN.ENUM.cell_system("ev")
        _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
        nl = sum(isinstance(lg, LambertLeg) for lg in one.legs)
        out = []
        for e in json.loads((REPO / "data" / "1007_v2" / sub / "realeph_chain.json").read_text()):
            conv = [st for st in e["steps"] if st.get("converged")]
            reached = bool(conv) and conv[-1]["lambda"] >= 1.0
            rec: dict[str, Any] = {
                "epoch_jd": e["epoch_jd"],
                "reached_lambda1": reached,
                "last_converged_lambda": e.get("last_converged_lambda"),
                "gate": conv[-1].get("status") if reached else None,
                "worst_ratio": conv[-1].get("worst_ratio") if reached else None,
            }
            if reached:
                y = e["final_y_dates"]
                x, v, t, c = legs_data("ev", key, 10, "ramp", e["x0_days"], y, real)
                for opt, st in (("a", 0), ("b", 1)):
                    rec[opt] = same_phase(f"{row} JD {e['epoch_jd']:.1f} ({opt})", x, v, t, c, 10,
                                          nl, y[10 * nl - 1], 5, st, quiet=True)  # fmt: skip
                    rec[opt]["v2_draft_pass"] = bool(
                        rec[opt]["drift_pass"]
                        and rec[opt]["continuity_pass"]
                        and rec["gate"] == "pass"
                    )
            print(
                f"{row} JD {e['epoch_jd']:.1f}: lambda1={reached} gate={rec['gate']} "
                f"worst={rec['worst_ratio']} "
                + (
                    " ".join(
                        f"({o}) D'={max(rec[o]['Dprime_d']):.2f}/{rec[o]['band_D_d']:.2f} "
                        f"W'={max(rec[o]['Wprime_kms']):.3f}/{rec[o]['band_W_kms']:.3f} "
                        f"g={rec[o]['growth_D']:.2f}/{rec[o]['growth_W']:.2f} "
                        f"{'PASS' if rec[o]['v2_draft_pass'] else 'FAIL'}"
                        for o in ("a", "b")
                    )
                    if reached
                    else "no 10-cycle chain"
                )
                + f" [{LABEL}]",
                flush=True,
            )
            out.append(rec)
        res[row] = out
        for o in ("a", "b"):
            k = sum(1 for r in out if r.get(o, {}).get("v2_draft_pass"))
            print(f"{row}: option ({o}) {k}/5 [{LABEL}]", flush=True)
    return res


def _final(label: str, j: dict[str, Any], sp: dict[str, Any]) -> dict[str, Any]:
    """Sec. 8.7: criteria 1-3 from judge (amended rules), criterion 4 = same-phase option (b)."""
    ok = bool(
        j.get("c1_continuity") and j.get("c2_gate") and j.get("c3_refly") and sp["drift_pass"]
    )
    rec = {
        "label": label,
        "c1_continuity": j.get("c1_continuity"),
        "junction_kms": j.get("junction_vinf_max_kms", j.get("continuity_max")),
        "c2_gate": j.get("c2_gate"),
        "gate": j.get("gate"),
        "c3_refly": j.get("c3_refly"),
        "refly_max_km": j.get("refly_max_km"),
        "c4_same_phase": sp,
        "pass": ok,
    }
    print(
        f"{label}: c1={j.get('c1_continuity')} c2={j.get('c2_gate')} c3={j.get('c3_refly')} "
        f"c4={sp['drift_pass']} (D'={max(sp['Dprime_d']):.3f}/{sp['band_D_d']:.3f} d "
        f"W'={max(sp['Wprime_kms']):.4f}/{sp['band_W_kms']:.3f} km/s "
        f"g={sp['growth_D']:.2f}/{sp['growth_W']:.2f}) -> {'PASS' if ok else 'FAIL'}",
        flush=True,
    )
    return rec


def rejudge_gc() -> dict[str, Any]:
    """Sec. 8.7: the final-wording re-judge of gc-1 and gc-2, controls first."""
    V2.RULES = "amended"
    out: dict[str, Any] = {"rules": "final wording (owner 2026-10-11): same-phase option (b)"}
    # positive control: H&M orbit 1 from the printed dates, 10 cycles
    real_ev = CHAIN.real_ephemeris("ev", "mean")
    key = "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"
    t3 = load_table3()[1]
    one_period = []
    for k in range(5):
        one_period += [t3[5 * k + 1].date, t3[5 * k + 4].date]
    starts = one_period + [v + 5844.0 for v in one_period]
    x0 = starts[0]
    y_seed = np.array([*starts[1:], 2 * 5844.0])
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    sol = least_squares(
        lambda yy: CHAIN.date_chain_residual(real_ev, one.legs * 10, x0 * DAY, yy),
        y_seed, method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200 * len(y_seed),
    )  # fmt: skip
    j = V2.judge("ev", key, [474.59719334719347, 1100.6183800415804], 10, "ramp", x0,
                 list(sol.x), real_ev, "H&M orbit 1")  # fmt: skip
    x, v, t, c = legs_data("ev", key, 10, "ramp", x0, list(sol.x), real_ev)
    sp = same_phase("H&M", x, v, t, c, 10, 2, sol.x[-1], 5, 1, quiet=True)
    out["positive"] = _final("CONTROL+ H&M orbit 1 (10 cycles)", j, sp)
    # negative control
    real_gc = CHAIN.real_ephemeris("gc", "spice")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    g1 = wb["gc-1"]
    e = json.loads((REPO / "data" / "943_gc1_realeph" / "e0" / "realeph_chain.json").read_text())[0]
    y_bad = list(e["final_y"])
    y_bad[0] += 0.05
    j = V2.judge("gc", g1["key"], g1["x_days"], 10, "blend", e["x0_days"], y_bad, real_gc, "neg")
    x, v, t, c = legs_data("gc", g1["key"], 10, "blend", e["x0_days"], y_bad, real_gc)
    sp = same_phase("neg", x, v, t, c, 10, 3, y_bad[29], 4, 1, quiet=True)
    out["negative"] = _final("CONTROL- gc-1 e0 broken", j, sp)
    pos_ok, neg_ok = out["positive"]["pass"], not out["negative"]["pass"]
    out["controls_ok"] = {"positive": pos_ok, "negative": neg_ok}
    if not (pos_ok and neg_ok):
        print("STOP: a control did not behave; rows not judged", flush=True)
        return out
    rows = []
    for i in range(5):
        e = json.loads(
            (REPO / "data" / "943_gc1_realeph" / f"e{i}" / "realeph_chain.json").read_text()
        )[0]
        lab = f"gc-1 JD {e['epoch_jd']:.1f}"
        j = V2.judge("gc", g1["key"], g1["x_days"], 10, "blend", e["x0_days"], e["final_y"],
                     real_gc, lab)  # fmt: skip
        x, v, t, c = legs_data("gc", g1["key"], 10, "blend", e["x0_days"], e["final_y"], real_gc)
        sp = same_phase(lab, x, v, t, c, 10, 3, e["final_y"][29], 4, 1, quiet=True)
        rows.append(_final(lab, j, sp))
    out["gc-1"] = rows
    g2 = wb["gc-2"]
    rows = []
    for e in json.loads((REPO / "data" / "942_943_ladder_realeph_chains.json").read_text())["gc2"]:
        xs = e["final_x_days"]
        y = list(xs[1:])
        lab = f"gc-2 JD {e['epoch_jd']:.1f}"
        j = V2.judge("gc", g2["key"], g2["x_days"], 10, "blend", xs[0], y, real_gc, lab)
        x, v, t, c = legs_data("gc", g2["key"], 10, "blend", xs[0], y, real_gc)
        sp = same_phase(lab, x, v, t, c, 10, 4, y[39], 4, 1, quiet=True)
        rows.append(_final(lab, j, sp))
    out["gc-2"] = rows
    for name in ("gc-1", "gc-2"):
        k = sum(r["pass"] for r in out[name])
        print(f"{name}: V2 (final wording) {'PASS' if k >= 3 else 'FAIL'} ({k}/5)", flush=True)
    return out


def _slope(v: list[float]) -> float:
    return float(np.polyfit(np.arange(len(v)), np.asarray(v), 1)[0])


def explore_c() -> dict[str, Any]:
    """Sec. 8.9 (EXPLORATORY data, not a verdict): option (C) evidence. H&M orbit 1 from its printed
    dates over 15 cycles (48 yr); gc-1 / gc-2 slopes on their existing 10-cycle chains."""
    V2.RULES = "amended"
    out: dict[str, Any] = {"status": "exploratory data only, not a verdict (sec. 8.9)"}
    real_ev = CHAIN.real_ephemeris("ev", "mean")
    key = "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"
    t3 = load_table3()[1]
    one_period = []
    for k in range(5):
        one_period += [t3[5 * k + 1].date, t3[5 * k + 4].date]
    starts = [v + q * 5844.0 for q in range(3) for v in one_period]
    x0 = starts[0]
    y_seed = np.array([*starts[1:], 3 * 5844.0])
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    sol = least_squares(
        lambda yy: CHAIN.date_chain_residual(real_ev, one.legs * 15, x0 * DAY, yy),
        y_seed, method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200 * len(y_seed),
    )  # fmt: skip
    j = V2.judge("ev", key, [474.59719334719347, 1100.6183800415804], 15, "ramp", x0,
                 list(sol.x), real_ev, "H&M orbit 1, 15 cycles")  # fmt: skip
    x, v, t, c = legs_data("ev", key, 15, "ramp", x0, list(sol.x), real_ev)
    sp = same_phase("H&M 15", x, v, t, c, 15, 2, sol.x[-1], 5, 1, quiet=True)
    hm = {
        "max_shift_from_print_d": float(np.max(np.abs(sol.x[:-1] - y_seed[:-1]))),
        "chain_period_d": float(sol.x[-1]),
        "continuity_junction_kms": c,
        "gate": j.get("gate"),
        "refly_max_km": j.get("refly_max_km"),
        "c1": j.get("c1_continuity"), "c2": j.get("c2_gate"), "c3": j.get("c3_refly"),
        "Dprime_d": sp["Dprime_d"], "Wprime_kms": sp["Wprime_kms"],
        "band_D_d": sp["band_D_d"], "band_W_kms": sp["band_W_kms"],
        "slope_D_d_per_cycle": _slope(sp["Dprime_d"]),
        "slope_W_kms_per_cycle": _slope(sp["Wprime_kms"]),
    }  # fmt: skip
    out["hm_orbit1_15"] = hm
    print(
        f"H&M 15 cycles: shift from print {hm['max_shift_from_print_d']:.2f} d, period "
        f"{hm['chain_period_d']:.1f} d (3 x 5844 = 17532), c1={hm['c1']} c2={hm['c2']} "
        f"({hm['gate']}) c3={hm['c3']} (re-fly {hm['refly_max_km']:.2e} km)",
        flush=True,
    )
    print("  D' =", [round(q, 3) for q in hm["Dprime_d"]], "band", round(hm["band_D_d"], 3),
          "slope", round(hm["slope_D_d_per_cycle"], 4), flush=True)  # fmt: skip
    print("  W' =", [round(q, 3) for q in hm["Wprime_kms"]], "band", round(hm["band_W_kms"], 3),
          "slope", round(hm["slope_W_kms_per_cycle"], 4), flush=True)  # fmt: skip
    real_gc = CHAIN.real_ephemeris("gc", "spice")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    for row, nl in (("gc-1", 3), ("gc-2", 4)):
        recs = []
        if row == "gc-1":
            chains = []
            for i in range(5):
                e = json.loads(
                    (REPO / "data" / "943_gc1_realeph" / f"e{i}" / "realeph_chain.json").read_text()
                )[0]
                chains.append((e["epoch_jd"], e["x0_days"], e["final_y"]))
        else:
            chains = [
                (e["epoch_jd"], e["final_x_days"][0], list(e["final_x_days"][1:]))
                for e in json.loads(
                    (REPO / "data" / "942_943_ladder_realeph_chains.json").read_text()
                )["gc2"]
            ]
        for jd, x0g, y in chains:
            x, v, t, c = legs_data("gc", wb[row]["key"], 10, "blend", x0g, y, real_gc)
            spg = same_phase(row, x, v, t, c, 10, nl, y[10 * nl - 1], 4, 1, quiet=True)
            rec = {
                "epoch_jd": jd, "Dprime_d": spg["Dprime_d"], "Wprime_kms": spg["Wprime_kms"],
                "band_D_d": spg["band_D_d"], "band_W_kms": spg["band_W_kms"],
                "slope_D_d_per_cycle": _slope(spg["Dprime_d"]),
                "slope_W_kms_per_cycle": _slope(spg["Wprime_kms"]),
            }  # fmt: skip
            recs.append(rec)
            print(
                f"{row} JD {jd:.1f}: D' max {max(rec['Dprime_d']):.3f}/{rec['band_D_d']:.3f} slope "
                f"{rec['slope_D_d_per_cycle']:+.4f}; W' max {max(rec['Wprime_kms']):.4f}/"
                f"{rec['band_W_kms']:.3f} slope {rec['slope_W_kms_per_cycle']:+.5f}",
                flush=True,
            )
        out[row] = recs
    return out


FINAL = "final wording (owner 2026-10-11): band only, ends skipped, slope as diagnostic"


def final_c4(x: np.ndarray, vdep: list[float], tof: list[float], n: int, nl: int,
             period_d: float, big_l: int) -> dict[str, Any]:  # fmt: skip
    """Criterion 4 under the final wording (sec. 8.4b): same-phase D' and W' over
    i = 1 .. n-L-2 (the first and last cycle excluded), band only; slopes are diagnostics."""
    p_mean = period_d / n
    idx = list(range(1, n - big_l - 1))
    assert n >= 2 * big_l and idx, (n, big_l)
    d = [max(abs(x[(i + big_l) * nl + j] - x[i * nl + j] - big_l * p_mean) for j in range(nl))
         for i in idx]  # fmt: skip
    w = [max(abs(vdep[(i + big_l) * nl + j] - vdep[i * nl + j]) for j in range(nl)) for i in idx]
    band_d, band_w = 0.05 * min(tof), 0.10 * min(vdep)
    return {
        "comparisons": idx, "Dprime_d": d, "Wprime_kms": w, "band_D_d": band_d,
        "band_W_kms": band_w,
        "slope_D_d_per_cycle": _slope(d) if len(d) > 1 else 0.0,
        "slope_W_kms_per_cycle": _slope(w) if len(w) > 1 else 0.0,
        "pass": bool(max(d) <= band_d and max(w) <= band_w),
    }  # fmt: skip


def final_record(label: str, cell: str, key: str, x_tpl: list[float], n: int, mode: str,
                 x0: float, y: list[float], real: Any, nl: int, big_l: int,
                 extra: dict[str, Any] | None = None) -> dict[str, Any]:  # fmt: skip
    V2.RULES = "amended"
    j = V2.judge(cell, key, x_tpl, n, mode, x0, y, real, label)
    if j.get("error"):
        rec = {"label": label, "pass": False, "error": j["error"], "rules": FINAL}
        print(f"{label}: {j['error']} -> FAIL", flush=True)
        return rec
    x, v, t, c = legs_data(cell, key, n, mode, x0, y, real)
    c4 = final_c4(x, v, t, n, nl, y[n * nl - 1], big_l)
    rec = {
        "label": label, "rules": FINAL,
        "c1_continuity": j["c1_continuity"], "junction_kms": c,
        "c2_gate": j["c2_gate"], "gate": j["gate"],
        "c3_refly": j["c3_refly"], "refly_max_km": j["refly_max_km"],
        "c4": c4,
        "pass": bool(j["c1_continuity"] and j["c2_gate"] and j["c3_refly"] and c4["pass"]),
    } | (extra or {})  # fmt: skip
    print(
        f"{label}: c1={int(j['c1_continuity'])} ({c:.1e}) c2={int(j['c2_gate'])} "
        f"({j['gate'].get('worst_ratio', float('nan')):.3f}) c3={int(j['c3_refly'])} "
        f"({j['refly_max_km']:.1e} km) c4={int(c4['pass'])} D'={max(c4['Dprime_d']):.3f}/"
        f"{c4['band_D_d']:.3f} W'={max(c4['Wprime_kms']):.4f}/{c4['band_W_kms']:.3f} "
        f"[slopes {c4['slope_D_d_per_cycle']:+.4f} d, {c4['slope_W_kms_per_cycle']:+.5f} km/s per "
        f"cycle] -> {'PASS' if rec['pass'] else 'FAIL'}",
        flush=True,
    )
    return rec


def _hm_chain(n: int) -> tuple[float, list[float], dict[str, Any]]:
    real_ev = CHAIN.real_ephemeris("ev", "mean")
    key = "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"
    t3 = load_table3()[1]
    one_period = []
    for k in range(5):
        one_period += [t3[5 * k + 1].date, t3[5 * k + 4].date]
    reps = n // 5
    starts = [v + q * 5844.0 for q in range(reps) for v in one_period]
    x0 = starts[0]
    y_seed = np.array([*starts[1:], reps * 5844.0])
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    sol = least_squares(
        lambda yy: CHAIN.date_chain_residual(real_ev, one.legs * n, x0 * DAY, yy),
        y_seed, method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200 * len(y_seed),
    )  # fmt: skip
    return (
        x0,
        list(sol.x),
        {"max_shift_from_print_d": float(np.max(np.abs(sol.x[:-1] - y_seed[:-1])))},
    )


def final_controls() -> dict[str, Any]:
    out: dict[str, Any] = {"rules": FINAL}
    real_ev = CHAIN.real_ephemeris("ev", "mean")
    x0, y, extra = _hm_chain(10)
    out["positive"] = final_record(
        "CONTROL+ H&M orbit 1 (10 cycles)", "ev", "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s",
        [474.59719334719347, 1100.6183800415804], 10, "ramp", x0, y, real_ev, 2, 5, extra,
    )  # fmt: skip
    real_gc = CHAIN.real_ephemeris("gc", "spice")
    g1 = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())["gc-1"]
    e = json.loads((REPO / "data" / "943_gc1_realeph" / "e0" / "realeph_chain.json").read_text())[0]
    y_bad = list(e["final_y"])
    y_bad[0] += 0.05
    out["negative"] = final_record("CONTROL- gc-1 e0 broken", "gc", g1["key"], g1["x_days"], 10,
                                   "blend", e["x0_days"], y_bad, real_gc, 3, 4)  # fmt: skip
    ok = out["positive"]["pass"] and not out["negative"]["pass"]
    out["controls_ok"] = ok
    print(f"CONTROLS {'behave' if ok else 'DO NOT behave: STOP'}", flush=True)
    return out


def final_gc() -> dict[str, Any]:
    out: dict[str, Any] = {"rules": FINAL}
    real_gc = CHAIN.real_ephemeris("gc", "spice")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    g1, g2 = wb["gc-1"], wb["gc-2"]
    rows = []
    for i in range(5):
        e = json.loads(
            (REPO / "data" / "943_gc1_realeph" / f"e{i}" / "realeph_chain.json").read_text()
        )[0]
        rows.append(final_record(f"gc-1 JD {e['epoch_jd']:.1f}", "gc", g1["key"], g1["x_days"], 10,
                                 "blend", e["x0_days"], e["final_y"], real_gc, 3, 4))  # fmt: skip
    out["gc-1"] = rows
    rows = []
    for e in json.loads((REPO / "data" / "942_943_ladder_realeph_chains.json").read_text())["gc2"]:
        xs = e["final_x_days"]
        rows.append(final_record(f"gc-2 JD {e['epoch_jd']:.1f}", "gc", g2["key"], g2["x_days"], 10,
                                 "blend", xs[0], list(xs[1:]), real_gc, 4, 4))  # fmt: skip
    out["gc-2"] = rows
    for name in ("gc-1", "gc-2"):
        k = sum(r["pass"] for r in out[name])
        print(f"{name}: V2 (final wording) {'PASS' if k >= 3 else 'FAIL'} ({k}/5)", flush=True)
    return out


def final_ev() -> dict[str, Any]:
    out: dict[str, Any] = {"rules": FINAL}
    real = CHAIN.real_ephemeris("ev", "mean")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    for row, sub in (("ev-A", "ev10_evA"), ("ev-B", "ev10_evB"), ("ev-C", "ev10_evC")):
        key, xt = wb[row]["key"], wb[row]["x_days"]
        circ, a, b = CHAIN.ENUM.cell_system("ev")
        _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
        nl = sum(isinstance(lg, LambertLeg) for lg in one.legs)
        rows = []
        for e in json.loads((REPO / "data" / "1007_v2" / sub / "realeph_chain.json").read_text()):
            conv = [st for st in e["steps"] if st.get("converged")]
            if not (conv and conv[-1]["lambda"] >= 1.0):
                lab = f"{row} JD {e['epoch_jd']:.1f}"
                rows.append(
                    {"label": lab, "pass": False, "error": "no 10-cycle chain", "rules": FINAL}
                )
                print(f"{row} JD {e['epoch_jd']:.1f}: no 10-cycle chain -> FAIL", flush=True)
                continue
            rows.append(final_record(f"{row} JD {e['epoch_jd']:.1f}", "ev", key, xt, 10, "ramp",
                                     e["x0_days"], e["final_y_dates"], real, nl, 5))  # fmt: skip
        out[row] = rows
        k = sum(r["pass"] for r in rows)
        print(f"{row}: V2 (final wording) {'PASS' if k >= 3 else 'FAIL'} ({k}/5)", flush=True)
    return out


def main() -> None:
    for flag, fn, name in (
        ("--final-controls", final_controls, "final_controls"),
        ("--final-gc", final_gc, "final_gc"),
        ("--final-ev", final_ev, "final_ev"),
    ):
        if flag in sys.argv:
            r = fn()
            path = REPO / "data" / "1007_v2" / f"{name}.json"
            path.write_text(json.dumps(r, indent=1, default=float))
            print("wrote", path.relative_to(REPO), flush=True)
            return
    if "--explore-c" in sys.argv:
        r = explore_c()
        path = REPO / "data" / "1007_v2" / "explore_c_15cycles.json"
        path.write_text(json.dumps(r, indent=1, default=float))
        print("wrote", path.relative_to(REPO), flush=True)
        return
    if "--rejudge-gc" in sys.argv:
        r = rejudge_gc()
        path = REPO / "data" / "1007_v2" / "rejudge_gc_final.json"
        path.write_text(json.dumps(r, indent=1, default=float))
        print("wrote", path.relative_to(REPO), flush=True)
        return
    if "--ev-rows" in sys.argv:
        r = ev_rows()
        path = REPO / "data" / "1007_v2" / "samephase_ev_rows_exploratory.json"
        path.write_text(json.dumps(r, indent=1, default=float))
        print("wrote", path.relative_to(REPO), flush=True)
        return
    res: dict[str, Any] = {"status": LABEL}
    real_gc = CHAIN.real_ephemeris("gc", "spice")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    rows = []
    g1 = wb["gc-1"]
    for i in range(5):
        e = json.loads(
            (REPO / "data" / "943_gc1_realeph" / f"e{i}" / "realeph_chain.json").read_text()
        )[0]
        x, v, t, c = legs_data("gc", g1["key"], 10, "blend", e["x0_days"], e["final_y"], real_gc)
        rows.append(same_phase(f"gc-1 JD {e['epoch_jd']:.1f}", x, v, t, c, 10, 3,
                               e["final_y"][29], 4))  # fmt: skip
    res["gc-1"] = rows
    g2 = wb["gc-2"]
    rows = []
    ladder = json.loads((REPO / "data" / "942_943_ladder_realeph_chains.json").read_text())["gc2"]
    for e in ladder:
        xs = e["final_x_days"]
        y = list(xs[1:])
        x, v, t, c = legs_data("gc", g2["key"], 10, "blend", xs[0], y, real_gc)
        rows.append(same_phase(f"gc-2 JD {e['epoch_jd']:.1f}", x, v, t, c, 10, 4, y[39], 4))
    res["gc-2"] = rows
    # broken-date negative
    e = json.loads((REPO / "data" / "943_gc1_realeph" / "e0" / "realeph_chain.json").read_text())[0]
    y_bad = list(e["final_y"])
    y_bad[0] += 0.05
    x, v, t, c = legs_data("gc", g1["key"], 10, "blend", e["x0_days"], y_bad, real_gc)
    res["negative"] = same_phase("NEGATIVE gc-1 e0 broken", x, v, t, c, 10, 3, y_bad[29], 4)
    # H&M orbit 1, printed dates, two periods
    real_ev = CHAIN.real_ephemeris("ev", "mean")
    key = "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"
    t3 = load_table3()[1]
    one_period = []
    for k in range(5):
        one_period += [t3[5 * k + 1].date, t3[5 * k + 4].date]
    starts = one_period + [s + 5844.0 for s in one_period]
    x0 = starts[0]
    y_seed = np.array([*starts[1:], 2 * 5844.0])
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    legs = one.legs * 10
    sol = least_squares(
        lambda yy: CHAIN.date_chain_residual(real_ev, legs, x0 * DAY, yy),
        y_seed, method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200 * len(y_seed),
    )  # fmt: skip
    x, v, t, c = legs_data("ev", key, 10, "ramp", x0, list(sol.x), real_ev)
    hm = same_phase("POSITIVE H&M orbit 1 (printed, 10 cycles)", x, v, t, c, 10, 2, sol.x[-1], 5)
    hm["max_shift_from_print_d"] = float(np.max(np.abs(sol.x[:-1] - y_seed[:-1])))
    hm["chain_period_d"] = float(sol.x[-1])
    print(f"H&M: max shift from print {hm['max_shift_from_print_d']:.3f} d, chain period "
          f"{hm['chain_period_d']:.3f} d (2 x 5844 = 11688)", flush=True)  # fmt: skip
    res["hm_orbit1"] = hm
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print("wrote", OUT.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
