"""#885 -- bridge from a published Uranus tour end-state into a catalogued quasi-cycler.

Staged; every stage writes JSONL under ``data/found/885_uranus_bridge/`` and prints
timestamped progress. Run each stage in the foreground; none takes more than a few
minutes on 4 workers.

    uv run python scripts/screen_885_uranus_tour_to_cycler_bridge.py --stage turns
    uv run python scripts/screen_885_uranus_tour_to_cycler_bridge.py --stage control
    uv run python scripts/screen_885_uranus_tour_to_cycler_bridge.py --stage ideal

Stages
------
turns    Demanded turn at every encounter of the six catalogued Uranian (1,1)
         quasi-cyclers, against the maximum ballistic bend (``core.flyby.max_bend``)
         at 50 km, in the rows' own construction (circular-coplanar moons,
         ``phase0`` = 30 deg, sorted-moon ``rel_offset`` convention of
         ``data/validation/v2_moontour.py``). Includes a flyable control (a
         same-conic continuation, whose demanded turn must be zero).
control  The pre-registered positive control of the Landau-Davis-Karimi tour step
         against three equatorial legs of their Figure 6 (100 km altitude).
ideal    Circular-coplanar Tisserand points of the published stretches and the
         cyclers, the published-V-infinity consistency check, and the phase-free
         flyby-count ladder from the published end-states to the cycler points.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.search import uranus_bridge_885 as ub

OUT_DIR = Path("data/found/885_uranus_bridge")
GAUNTLET_566 = Path("data/gauntlet_566_five_representatives.jsonl")
SCAN_567 = Path("data/scan_567_epoch_robustness.jsonl")
PHASE0_DEG = 30.0  # #566 meta "phase0_deg": 29.999999999999996


def log(msg: str) -> None:
    print(f"[{dt.datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)


def write_jsonl(name: str, rows: list[dict[str, Any]]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    with path.open("w") as fh:
        for r in rows:
            fh.write(json.dumps(r, default=float) + "\n")
    log(f"wrote {path} ({len(rows)} rows)")
    return path


def git_sha() -> str:
    import subprocess

    return subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=False
    ).stdout.strip()


# ---------------------------------------------------------------------------
# Stage: turns
# ---------------------------------------------------------------------------

CATALOGUE_IDS = {
    "#312 Umbriel-Oberon-Umbriel": "umbriel-oberon-1-1-uranian-quasi-cycler-2026",
    "Titania-Oberon-Titania": "titania-oberon-1-1-uranian-quasi-cycler-2026",
    "Ariel-Umbriel-Ariel": "ariel-umbriel-1-1-uranian-quasi-cycler-2026",
    "Ariel-Titania-Ariel": "ariel-titania-1-1-uranian-quasi-cycler-2026",
    "Ariel-Oberon-Ariel": "ariel-oberon-1-1-uranian-quasi-cycler-2026",
    "Umbriel-Titania-Umbriel": "umbriel-titania-1-1-uranian-quasi-cycler-2026",
}


def _row_phases(seq: list[str], rel_offset_deg: float) -> ub.MoonPhases:
    """The v2_moontour convention: sorted distinct moons, first at phase0, second +rel."""
    s = sorted(set(seq))
    th = {m: 0.0 for m in ub.TOUR_MOONS}
    th[s[0]] = math.radians(PHASE0_DEG)
    th[s[1]] = math.radians(PHASE0_DEG + rel_offset_deg)
    return ub.MoonPhases(theta0=th)


def _turn_rows() -> list[dict[str, Any]]:
    cands = json.loads(SCAN_567.open().readline())["candidates"]
    rows: list[dict[str, Any]] = []
    for c in cands:
        label = next(k for k in CATALOGUE_IDS if c["label"].startswith(k))
        seq = c["sequence"]
        a, b = seq[0], seq[1]
        ph = _row_phases(seq, c["rel_offset_deg"])
        leg_s = c["tof_days"][0] * ub.DAY_S
        nrev = int(c["n_revs"][0])
        cg = ub.CyclerGeometry(a, b, leg_s, 0.0, nrev)
        # Encounters A(0) -> B(leg) -> A(2 leg) -> B(3 leg): legs solved in the
        # row's own way (Lambert between circular moon positions, its n_rev).
        out0, in1 = cg.leg_vinf(ph, 0.0, a, b)
        out1, in2 = cg.leg_vinf(ph, leg_s, b, a)
        out2, _ = cg.leg_vinf(ph, 2 * leg_s, a, b)
        encs = []
        for moon, vin, vout, t in ((b, in1, out1, leg_s), (a, in2, out2, 2 * leg_s)):
            turn = abs(ub.signed_angle(vin, vout))
            vmag = float(np.linalg.norm(vin))
            bmax50 = ub.max_bend_rad(moon, vmag, 50.0)
            bmax0 = ub.max_bend_rad(moon, vmag, 0.0)
            encs.append(
                {
                    "moon": moon,
                    "t_days": t / ub.DAY_S,
                    "vinf_in_kms": vmag,
                    "vinf_out_kms": float(np.linalg.norm(vout)),
                    "demanded_turn_deg": math.degrees(turn),
                    "max_bend_50km_deg": math.degrees(bmax50),
                    "max_bend_surface_deg": math.degrees(bmax0),
                    "demanded_over_max_50km": turn / bmax50,
                    "impulse_beyond_bend_kms": ub.entry_correction_kms(vin, vout, bmax50),
                }
            )
        rows.append(
            {
                "kind": "turn_check",
                "catalogue_id": CATALOGUE_IDS[label],
                "candidate_id": c["candidate_id"],
                "sequence": seq,
                "rel_offset_deg": c["rel_offset_deg"],
                "phase0_deg": PHASE0_DEG,
                "n_revs": c["n_revs"],
                "tof_days": c["tof_days"],
                "source_vinf_kms": c["vinf_kms"],
                "reproduced_vinf_kms": [
                    float(np.linalg.norm(out0)),
                    float(np.linalg.norm(in1)),
                    float(np.linalg.norm(in2)),
                ],
                "encounters": encs,
                "convention": (
                    "circular-coplanar moons, sorted distinct moons: theta[sorted[0]] = phase0, "
                    "theta[sorted[1]] = phase0 + rel_offset (v2_moontour.run_v2_moontour); "
                    "each leg the row's Lambert at its n_rev, min |v1 - v_moon| branch"
                ),
            }
        )
    return rows


def _turn_control() -> dict[str, Any]:
    """A flyable case: continue ONE conic through two moons; demanded turn must be 0.

    The Ellison et al. Oberon (3.448) / Ariel (4.325) orbit is placed so that it
    meets Oberon (inbound), then Ariel (inbound), then the Ariel radius again
    (outbound); Oberon and Ariel are phased to sit at their crossings. Lambert legs
    between those points must recover the same conic, so the demanded turn at the
    Ariel node is zero to solver precision.
    """
    a, e = ub.orbit_from_vinf_pair("Oberon", 3.448, "Ariel", 4.325)
    el = ub.Elements(
        a=a, e=e, omega=0.0, h=math.sqrt(ub.MU_URANUS * a * (1 - e * e)), mean_anomaly=0.0
    )

    def crossing(moon: str, outbound: bool) -> float:
        r = ub.moon_sma(moon)
        p = a * (1 - e * e)
        nu = math.acos((p / r - 1.0) / e)
        nu = nu if outbound else 2 * math.pi - nu
        ea = 2 * math.atan2(
            math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2)
        )
        return (ea - e * math.sin(ea)) % (2 * math.pi)

    n = el.n
    # Oberon inbound (before periapsis), then Ariel inbound, then Ariel outbound.
    m_o = crossing("Oberon", False) - 2 * math.pi
    m_a1 = crossing("Ariel", False) - 2 * math.pi
    m_a2 = crossing("Ariel", True)
    t_o, t_a1, t_a2 = 0.0, (m_a1 - m_o) / n, (m_a2 - m_o) / n
    r_o, _ = ub.state_at_mean_anomaly(el, m_o)
    r_a1, _ = ub.state_at_mean_anomaly(el, m_a1)
    r_a2, _ = ub.state_at_mean_anomaly(el, m_a2)
    th = {m: 0.0 for m in ub.TOUR_MOONS}
    th["Oberon"] = math.atan2(r_o[1], r_o[0])
    # Ariel is phased to sit at r_a1 at t_a1. The second leg ends at the conic's
    # own next Ariel-radius point r_a2 (no moon need be there): the check is only
    # that Lambert legs on either side of a node recover one conic (turn = 0).
    th["Ariel"] = math.atan2(r_a1[1], r_a1[0]) - ub.moon_mean_motion("Ariel") * t_a1
    ph = ub.MoonPhases(theta0=th)
    from cyclerfinder.core.lambert import lambert

    def lam(r1: np.ndarray, r2: np.ndarray, tof: float) -> tuple[np.ndarray, np.ndarray]:
        s = lambert(np.r_[r1, 0.0], np.r_[r2, 0.0], tof, mu=ub.MU_URANUS, prograde=True)[0]
        return np.asarray(s.v1[:2]), np.asarray(s.v2[:2])

    ro, _ = ph.state("Oberon", t_o)
    ra1, vma1 = ph.state("Ariel", t_a1)
    _, v_arr = lam(ro, ra1, t_a1 - t_o)
    v_dep, _ = lam(ra1, r_a2, t_a2 - t_a1)
    vin, vout = v_arr - vma1, v_dep - vma1
    turn = abs(ub.signed_angle(vin, vout))
    return {
        "kind": "turn_control_same_conic",
        "description": (
            "Ellison Oberon 3.448 / Ariel 4.325 conic continued through Ariel to its next "
            "Ariel-radius crossing; Lambert legs between the positions; demanded turn at the "
            "Ariel node must be 0"
        ),
        "vinf_in_kms": float(np.linalg.norm(vin)),
        "vinf_out_kms": float(np.linalg.norm(vout)),
        "demanded_turn_deg": math.degrees(turn),
        "position_check_km": float(np.linalg.norm(ra1 - r_a1)),
        "passes": bool(math.degrees(turn) < 1e-3),
    }


def stage_turns() -> None:
    log("stage turns: start")
    rows = _turn_rows()
    for r in rows:
        e1, e2 = r["encounters"]
        log(
            f"{r['catalogue_id']}: reproduced V-inf "
            f"{[f'{v:.4f}' for v in r['reproduced_vinf_kms']]} "
            f"(source {[f'{v:.4f}' for v in r['source_vinf_kms']]}); "
            f"turn at {e1['moon']} {e1['demanded_turn_deg']:.1f} deg "
            f"vs max {e1['max_bend_50km_deg']:.2f}; "
            f"at {e2['moon']} {e2['demanded_turn_deg']:.1f} vs {e2['max_bend_50km_deg']:.2f}"
        )
    ctl = _turn_control()
    log(
        f"turn control (same conic): demanded turn {ctl['demanded_turn_deg']:.2e} deg, "
        f"passes={ctl['passes']}"
    )
    meta = {
        "_meta": True,
        "task": "#885 stage turns",
        "git_sha": git_sha(),
        "phase0_deg": PHASE0_DEG,
        "sources": [str(SCAN_567), str(GAUNTLET_566)],
    }
    write_jsonl("turn_check_six_rows.jsonl", [meta, ctl, *rows])


# ---------------------------------------------------------------------------
# Stage: control
# ---------------------------------------------------------------------------

KERNELS: list[str] = []


def _kernels() -> list[str]:
    from cyclerfinder.data.validation.v4_uranus_strict import (
        DEFAULT_LSK_PATH,
        DEFAULT_PCK_PATH,
        DEFAULT_URA_PATH,
    )

    return [str(DEFAULT_LSK_PATH), str(DEFAULT_PCK_PATH), str(DEFAULT_URA_PATH)]


CONTROL_LEGS = [
    # label, dep, vinf_in, arr, vinf_out, tof_days, departure date (published; hour assumed),
    # published dv m/s
    ("C1", "Ariel", 3.990, "Ariel", 3.897, 26.0, "2050-03-04", 15.0),
    ("C2", "Umbriel", 4.095, "Umbriel", 4.070, 23.0, "2050-04-23", 5.0),
    ("C3", "Ariel", 3.897, "Umbriel", 4.095, 24.0, "2050-03-30", 8.0),
]


def _control_job(
    args: tuple[str, str, float, str, float, float, str, float, int, float, float, int],
) -> dict[str, Any]:
    label, dep, vin, arr, vout, tof_d, date, pub, hour, dal_deg, dt_d, n_bend = args
    ph = ub.phases_from_ura111(f"{date}T{hour:02d}:00:00", _kernels())
    r0, vm0 = ph.state(dep, 0.0)
    that = vm0 / np.linalg.norm(vm0)
    rhat = r0 / np.linalg.norm(r0)
    dmax = ub.max_bend_rad(dep, vin, 100.0)
    ts = np.arange(tof_d - 1.0, tof_d + 1.0 + 1e-9, dt_d)
    best: dict[str, Any] | None = None
    n_cross = 0
    for sgn in (1.0, -1.0):
        for al_deg in np.arange(0.0, 180.0 + 1e-9, dal_deg):
            al = math.radians(al_deg)
            vinf_in = vin * (math.cos(al) * that + sgn * math.sin(al) * rhat)
            for bend in np.linspace(-dmax, dmax, n_bend):
                vo = ub.rotate(vinf_in, bend)
                prev: dict[tuple[int, str], ub.Leg] = {}
                for t in ts:
                    legs = ub.landau_leg(
                        dep_moon=dep,
                        t_dep=0.0,
                        vinf_out=vo,
                        arr_moon=arr,
                        t_arr=float(t) * ub.DAY_S,
                        phases=ph,
                        rp_floor_km=ub.RP_FLOOR_KM,
                    )
                    cur = {(lg.n_revs, lg.branch): lg for lg in legs}
                    for key, lg in cur.items():
                        p = prev.get(key)
                        if p is None:
                            continue
                        f1 = float(np.linalg.norm(p.vinf_in)) - vout
                        f2 = float(np.linalg.norm(lg.vinf_in)) - vout
                        if f1 * f2 > 0.0 or f1 == f2:
                            continue
                        w = f1 / (f1 - f2)
                        dv = p.dv_kms + (lg.dv_kms - p.dv_kms) * w
                        n_cross += 1
                        if best is None or dv < best["dv_ms"] / 1000.0:
                            best = {
                                "dv_ms": dv * 1000.0,
                                "pump_in_deg": al_deg * sgn,
                                "bend_deg": math.degrees(bend),
                                "bend_max_deg": math.degrees(dmax),
                                "tof_days": float(t) - dt_d + w * dt_d,
                                "n_revs": key[0],
                                "branch": key[1],
                                "burn_radius_km": float(np.linalg.norm(lg.r_burn)),
                            }
                    prev = cur
    return {
        "kind": "control_leg",
        "label": label,
        "published": {
            "dep": dep,
            "vinf_in": vin,
            "arr": arr,
            "vinf_out": vout,
            "tof_days": tof_d,
            "date": date,
            "dv_ms": pub,
        },
        "departure_hour_utc_assumed": hour,
        "grid": {
            "pump_step_deg": dal_deg,
            "tof_step_days": dt_d,
            "n_bend_samples": n_bend,
            "tof_window_days": [tof_d - 1, tof_d + 1],
            "alt_km": 100.0,
            "rp_floor_km": ub.RP_FLOOR_KM,
        },
        "n_solutions_meeting_vinf_out": n_cross,
        "best": best,
    }


def stage_control(workers: int) -> None:
    log("stage control: start")
    jobs = []
    for label, dep, vin, arr, vout, tof, date, pub in CONTROL_LEGS:
        hours = [12] if label != "C3" else [0, 6, 12, 18]
        for h in hours:
            jobs.append((label, dep, vin, arr, vout, tof, date, pub, h, 0.5, 0.02, 9))
    t0 = time.time()
    rows = []
    with Pool(workers) as pool:
        for r in pool.imap_unordered(_control_job, jobs):
            b = r["best"]
            best_txt = "none" if b is None else f"{b['dv_ms']:.2f} m/s"
            log(
                f"{r['label']} hour {r['departure_hour_utc_assumed']:02d}: best {best_txt} "
                f"vs published {r['published']['dv_ms']} "
                f"({time.time() - t0:.0f} s)"
            )
            rows.append(r)
    # Pre-registered verdict (note Section 2). C1, C2 use the 12:00 UTC row; C3 the
    # minimum over departure hours (the hour scan was added after the coarse run;
    # the 12:00 value is also reported).
    verdict_rows = []
    for label, *_rest, pub in CONTROL_LEGS:
        cands = [r for r in rows if r["label"] == label and r["best"] is not None]
        best12 = next(
            (r["best"]["dv_ms"] for r in cands if r["departure_hour_utc_assumed"] == 12), None
        )
        bmin = min((r["best"]["dv_ms"] for r in cands), default=None)
        verdict_rows.append(
            {
                "label": label,
                "published_ms": pub,
                "ours_min_ms": bmin,
                "ours_12h_ms": best12,
                "upper_ok": bmin is not None and bmin <= pub + 3.0,
                "lower_ok": bmin is not None and bmin >= 0.5 * pub - 3.0,
            }
        )
    upper_all = all(v["upper_ok"] for v in verdict_rows)
    lower_two = sum(v["lower_ok"] for v in verdict_rows) >= 2
    summary = {
        "kind": "control_verdict",
        "per_leg": verdict_rows,
        "pass": upper_all and lower_two,
        "criterion": (
            "pre-registered: all legs ours <= pub + 3 m/s AND >= 2 legs ours >= 0.5 pub - 3 m/s"
        ),
    }
    log(f"control verdict: {'PASS' if summary['pass'] else 'FAIL'} {verdict_rows}")
    meta = {
        "_meta": True,
        "task": "#885 stage control",
        "git_sha": git_sha(),
        "source": "Landau, Davis & Karimi 2023 (AAS 23-460) Figure 6, events 14-20",
    }
    write_jsonl(
        "control_landau2023_fig6.jsonl",
        [meta, summary, *sorted(rows, key=lambda r: (r["label"], r["departure_hour_utc_assumed"]))],
    )


# ---------------------------------------------------------------------------
# Stage: ideal
# ---------------------------------------------------------------------------

PUBLISHED_STRETCHES = {
    # name: (moon1, vinf1, moon2 (the LAST flyby of the stretch), vinf2,
    #        published V-inf at other moons)
    "ellison2025_table8_oberon_ariel": (
        "Oberon",
        3.448,
        "Ariel",
        4.325,
        {"Umbriel": [4.392, 4.396], "Titania": [3.880, 3.922], "Miranda": [3.311, 3.370]},
    ),
    "landau2025_table3_umbriel_oberon": (
        "Umbriel",
        4.1,
        "Oberon",
        2.6,
        {"Ariel": [3.9, 4.3], "Titania": [3.7, 4.2], "Miranda": [2.7, 3.7]},
    ),
}
CYCLERS = {
    "ariel-oberon-1-1-uranian-quasi-cycler-2026": (
        "Ariel",
        1.520866047614147,
        "Oberon",
        1.8285940380726622,
    ),
    "umbriel-oberon-1-1-uranian-quasi-cycler-2026": (
        "Umbriel",
        0.9199258810725036,
        "Oberon",
        0.9604309791298091,
    ),
}


def _orbit_row(name: str, m1: str, v1: float, m2: str, v2: float) -> dict[str, Any]:
    a, e = ub.orbit_from_vinf_pair(m1, v1, m2, v2)
    allm = ("Miranda", *ub.TOUR_MOONS)
    return {
        "kind": "tisserand_point",
        "name": name,
        "defined_by": {m1: v1, m2: v2},
        "a_km": a,
        "e": e,
        "rp_km": a * (1 - e),
        "ra_km": a * (1 + e),
        "period_days": 2 * math.pi * math.sqrt(a**3 / ub.MU_URANUS) / ub.DAY_S,
        "specific_energy_km2s2": -ub.MU_URANUS / (2 * a),
        "h_km2s": math.sqrt(ub.MU_URANUS * a * (1 - e * e)),
        "vinf_at_moons_kms": {m: ub.vinf_at(m, a, e) for m in allm},
        "pump_deg_at_moons": {
            m: (math.degrees(p[1]) if (p := ub.arrival_pump(m, a, e)) is not None else None)
            for m in ub.TOUR_MOONS
        },
        "max_bend_50km_deg": {
            m: (
                math.degrees(ub.max_bend_rad(m, v, 50.0))
                if not math.isnan(v := ub.vinf_at(m, a, e))
                else None
            )
            for m in ub.TOUR_MOONS
        },
    }


def _ladder_job(args: tuple[str, str, str, float, tuple[str, ...], float, int]) -> dict[str, Any]:
    name, pub, cyc, alt, moons, pump_bucket_deg, max_front = args
    m1, v1, m2, v2, _ = PUBLISHED_STRETCHES[pub]
    a, e = ub.orbit_from_vinf_pair(m1, v1, m2, v2)
    st = ub.arrival_pump(m2, a, e)
    assert st is not None
    c1, w1, c2, w2 = CYCLERS[cyc]
    ac, ec = ub.orbit_from_vinf_pair(c1, w1, c2, w2)
    targets = {}
    for cm in (c1, c2):
        tp = ub.arrival_pump(cm, ac, ec)
        assert tp is not None
        targets[cm] = tp
    t0 = time.time()
    res = ub.tisserand_ladder(
        (m2, st[0], st[1]),
        targets,
        moons=moons,
        alt_km=alt,
        pump_bucket_rad=math.radians(pump_bucket_deg),
        max_front=max_front,
    )
    path = []
    for s in res.path:
        orb = ub.orbit_from_pump(s[0], s[1], s[2])
        path.append(
            {
                "moon": s[0],
                "vinf_kms": s[1],
                "pump_deg": math.degrees(s[2]),
                "max_bend_deg": math.degrees(ub.max_bend_rad(s[0], s[1], alt)),
            }
        )
        del orb
    return {
        "kind": "ladder",
        "name": name,
        "from": pub,
        "to": cyc,
        "alt_km": alt,
        "moons": list(moons),
        "rp_floor_km": ub.RP_FLOOR_KM,
        "pump_bucket_deg": pump_bucket_deg,
        "n_flybys": res.n_flybys,
        "n_states": res.n_states,
        "wall_s": time.time() - t0,
        "path": path,
    }


def stage_ideal(workers: int) -> None:
    log("stage ideal: start")
    rows: list[dict[str, Any]] = []
    for name, (m1, v1, m2, v2, other) in PUBLISHED_STRETCHES.items():
        r = _orbit_row(name, m1, v1, m2, v2)
        r["published_vinf_other_moons_kms"] = other
        rows.append(r)
        log(
            f"{name}: a={r['a_km']:.0f} e={r['e']:.4f} rp={r['rp_km']:.0f} "
            f"P={r['period_days']:.2f} d "
            f"V-inf {{{', '.join(f'{k}: {v:.3f}' for k, v in r['vinf_at_moons_kms'].items())}}}"
        )
    for name, (m1, v1, m2, v2) in CYCLERS.items():
        r = _orbit_row(name, m1, v1, m2, v2)
        rows.append(r)
        log(
            f"{name}: a={r['a_km']:.0f} e={r['e']:.4f} rp={r['rp_km']:.0f} "
            f"P={r['period_days']:.2f} d"
        )
    jobs = [
        (
            "AO_4moons_50km",
            "ellison2025_table8_oberon_ariel",
            "ariel-oberon-1-1-uranian-quasi-cycler-2026",
            50.0,
            ub.TOUR_MOONS,
            0.05,
            600_000,
        ),
        (
            "AO_4moons_100km",
            "ellison2025_table8_oberon_ariel",
            "ariel-oberon-1-1-uranian-quasi-cycler-2026",
            100.0,
            ub.TOUR_MOONS,
            0.05,
            600_000,
        ),
        (
            "AO_ArielOberon_50km",
            "ellison2025_table8_oberon_ariel",
            "ariel-oberon-1-1-uranian-quasi-cycler-2026",
            50.0,
            ("Ariel", "Oberon"),
            0.05,
            600_000,
        ),
        (
            "UO_4moons_50km",
            "landau2025_table3_umbriel_oberon",
            "umbriel-oberon-1-1-uranian-quasi-cycler-2026",
            50.0,
            ub.TOUR_MOONS,
            0.1,
            1_000_000,
        ),
    ]
    t0 = time.time()
    with Pool(workers) as pool:
        for r in pool.imap_unordered(_ladder_job, jobs):
            log(
                f"ladder {r['name']}: flybys={r['n_flybys']} states={r['n_states']} "
                f"({r['wall_s']:.0f} s; {time.time() - t0:.0f} s elapsed)"
            )
            rows.append(r)
    meta = {
        "_meta": True,
        "task": "#885 stage ideal",
        "git_sha": git_sha(),
        "model": (
            "circular-coplanar patched conics about Uranus; "
            "registry moon constants (core/satellites.py)"
        ),
    }
    write_jsonl("ideal_model.jsonl", [meta, *rows])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--stage", required=True, choices=["turns", "control", "ideal"])
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args(argv)
    workers = min(4, max(1, args.workers))
    if args.stage == "turns":
        stage_turns()
    elif args.stage == "control":
        stage_control(workers)
    else:
        stage_ideal(workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
