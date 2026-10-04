"""#895 staged driver: ballistic Titania-Oberon arcs in a URA111 force model.

Pre-registration: ``docs/notes/2026-10-04-895-titania-oberon-realeph.md`` section 1.

This is a targeted correction around one known seed (the `#890` orbit) at five pre-registered
epochs, not a sweep of a region, so it does not call ``preflight_search`` (the ratchet in
``tests/scripts/test_scripts_call_preflight.py`` covers ``scripts/run_*.py``; the `#890` driver
``screen_890_*.py`` is the same kind of script).

Stages (each resumable, each invocation bounded by ``--budget`` seconds of wall time):

* ``control``  force-model positive control P1 and the table interpolation check;
* ``p2``       the `#890` model orbit rebuilt in this frame (P2);
* ``arc``      homotopy for one epoch and cycle count (checkpointed after every step);
* ``verify``   criteria (c) to (f) on a converged arc;
* ``sunoff``   the Sun's effect (one Newton solve with the Sun off, from a converged arc).

Run with ``uv run python scripts/screen_895_titania_oberon_realeph.py <stage> [options]``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

import cyclerfinder.search.titania_oberon_realeph_895 as m

OUT = Path(__file__).resolve().parents[1] / "data" / "found" / "895_titania_oberon_realeph"
CACHE = Path(os.environ.get("CF895_CACHE", "/tmp/cf895_cache"))
REFINE_890 = OUT.parent / "890_titania_oberon_candidate" / "refine.json"

#: Pre-registered target dates (section 1.4), 00:00 TDB.
EPOCHS: dict[str, str] = {
    "E1": "2030-01-12 00:00:00 TDB",
    "E2": "2031-06-13 00:00:00 TDB",
    "E3": "2035-07-01 00:00:00 TDB",
    "E4": "2040-07-01 00:00:00 TDB",
    "E5": "2045-07-01 00:00:00 TDB",
}
FIT_MARGIN_D = 30.0
TABLE_MARGIN_D = 60.0
MAX_CYCLES = 12
NOMINAL_CYCLE_D = 123.2


def log(msg: str) -> None:
    print(f"[{datetime.now(UTC).strftime('%H:%M:%S')}Z] {msg}", flush=True)


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def dump(name: str, obj: dict[str, Any]) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    obj = {"_meta": {"task": "#895", "git_sha": git_sha(), "stage": name}, **obj}
    path = OUT / name
    path.write_text(json.dumps(obj, indent=1, default=float) + "\n")
    log(f"wrote {path}")
    return path


# ------------------------------------------------------------------------------------------- #
# Epoch set-up: table, circles, conjunction
# ------------------------------------------------------------------------------------------- #


def epoch_setup(tag: str) -> tuple[m.EphemerisTable, m.Circles, float]:
    t_ref = m.et_of(EPOCHS[tag])
    t0 = -TABLE_MARGIN_D * m.DAY_S
    t1 = (MAX_CYCLES * NOMINAL_CYCLE_D + TABLE_MARGIN_D + 20.0) * m.DAY_S
    table = m.build_table(t_ref, t0, t1, cache_dir=CACHE)
    circles = m.fit_circles(
        table, -FIT_MARGIN_D * m.DAY_S, (MAX_CYCLES * NOMINAL_CYCLE_D + FIT_MARGIN_D) * m.DAY_S
    )
    t_c = circles.conjunction_near(0.0)
    return table, circles, t_c


def stage_control(args: argparse.Namespace) -> None:
    rows: list[dict[str, Any]] = []
    interp: list[dict[str, Any]] = []
    epochs_info: dict[str, Any] = {}
    variants: dict[str, dict[str, Any]] = {
        "full": {},
        "no_sun": {"sun": False},
        "no_j4": {"j4": 0.0},
        "french_j2": {"j2": m.J2_FRENCH_2024},
        "ura111_pole": {"pole": m.pole_vector(*m.POLE_URA111_RA_DEC_DEG)},
        "no_self_mass": {"self_mass": False},
        "no_zonal": {"j2": 0.0, "j4": 0.0},
        "former_repo_j2": {"j2": m.J2_FORMER_REPO, "j4": 0.0},
    }
    sp = m.load_kernels()
    for tag in EPOCHS:
        table, circles, t_c = epoch_setup(tag)
        ang = math.degrees(math.acos(float(np.clip(-np.dot(circles.ez, m.POLE_IAU), -1, 1))))
        epochs_info[tag] = {
            "t_ref_et": table.t_ref_et,
            "t_c_s": t_c,
            "t_c_tdb": m.tdb_string(table.t_ref_et + t_c),
            "circles": circles.to_dict(),
            "angle_ez_from_minus_iau_pole_deg": ang,
        }
        log(f"{tag}: conjunction {epochs_info[tag]['t_c_tdb']}, ez vs -pole {ang:.4f} deg")
        # interpolation check at mid-knots over the first 40 days
        for k, naif in enumerate(m.NAIF_IDS):
            worst = 0.0
            model = m.full_model(table)
            for i in range(200):
                t = table.t0 + (i * 17 + 0.5) * table.h
                ref = np.asarray(sp.spkgeo(naif, table.t_ref_et + t, "J2000", 799)[0])
                p, _ = model.body_states(t)
                worst = max(worst, float(np.linalg.norm(p[k] - ref[:3])))
            interp.append(
                {"epoch": tag, "body": m.BODY_NAMES[k], "max_mid_knot_err_m": worst * 1e3}
            )
        for name, kw in variants.items():
            if args.quick and name not in ("full", "no_zonal", "former_repo_j2"):
                continue
            for k in m.MOONS:
                err = m.moon_positive_control(table, k, t_c, (30.0, 123.0), **kw)
                rows.append(
                    {"epoch": tag, "variant": name, "moon": m.BODY_NAMES[k], "err30_km": err[0],
                     "err123_km": err[1]}
                )  # fmt: skip
            log(
                f"{tag} {name}: "
                + ", ".join(
                    f"{r['moon']} {r['err30_km']:.2f}/{r['err123_km']:.1f}" for r in rows[-5:]
                )
            )
    full = [r for r in rows if r["variant"] == "full"]
    p1_pass = all(r["err30_km"] <= 10.0 and r["err123_km"] <= 40.0 for r in full)
    disc = {}
    for v in ("no_zonal", "former_repo_j2"):
        sel = [r for r in rows if r["variant"] == v and r["moon"] in ("Miranda", "Ariel")]
        disc[v] = all(r["err30_km"] > 10.0 for r in sel)
    dump(
        "control.json",
        {
            "criterion": "P1: every moon, every epoch: <= 10 km at 30 d and <= 40 km at 123 d "
            "(variant 'full'); discrimination: Miranda and Ariel > 10 km at 30 d for "
            "'no_zonal' and 'former_repo_j2' at every epoch",
            "P1_pass": p1_pass,
            "discrimination_pass": disc,
            "epochs": epochs_info,
            "interpolation": interp,
            "rows": rows,
        },
    )
    log(f"P1 pass {p1_pass}; discrimination {disc}")


def stage_p2(args: argparse.Namespace) -> None:
    c = m.registry_890()
    s890 = np.asarray(json.loads(REFINE_890.read_text())["refined_state"])
    y0 = m.state_890_to_inertial(c, s890)
    t_cyc = c.cycle_s
    out: dict[str, Any] = {"cycle_d": t_cyc / m.DAY_S}
    for label, model in (
        ("model_890", m.model_890(c)),
        ("oberon_about_uranus", m.model_890(c, oberon_about_barycentre=False)),
    ):
        res: dict[str, Any] = {}
        for rtol in (1e-12, 1e-13):
            half = m.propagate(model, 0.0, t_cyc / 2, y0, rtol=rtol, atol_r=1e-10, atol_v=1e-16)
            rot = m.inertial_to_rot_890(c, t_cyc / 2, half.state)
            full = m.propagate(
                model, 0.0, t_cyc, y0, rtol=rtol, atol_r=1e-10, atol_v=1e-16, record=True,
                max_steps=40000, allow_incomplete=True,
            )  # fmt: skip
            r: dict[str, Any] = {
                "half_cycle_y_m": rot[1] * c.a_t * 1e3,
                "half_cycle_vx_mm_s": rot[2] * c.a_t * c.n_t * 1e6,
                "completed_cycle": full.complete,
            }
            if full.complete:
                end = m.inertial_to_rot_890(c, t_cyc, full.state)
                r["return_m"] = float(np.linalg.norm(end[:2] - s890[:2]) * c.a_t * 1e3)
                r["return_cm_s"] = float(np.linalg.norm(end[2:] - s890[2:]) * c.a_t * c.n_t * 1e5)
            assert full.rec_t is not None and full.rec_y is not None
            tr = m.hermite_traj(full.rec_t, full.rec_y, model)
            for k in (m.I_TITANIA, m.I_OBERON):

                def bst(
                    t: float, k: int = k, model: m.ForceModel = model
                ) -> tuple[np.ndarray, np.ndarray]:
                    pp, vv = model.body_states(t)
                    return pp[k], vv[k]

                enc = m.find_minima(tr, bst, 0.0, float(full.rec_t[-1]), body=k, max_dist=30000.0)
                r[m.BODY_NAMES[k]] = [{"t_d": e.t / m.DAY_S, **e.osculating()} for e in enc]
            res[f"rtol_{rtol:g}"] = r
        out[label] = res
        for sol_name, meth, rt in (
            ("scipy_DOP853_1e-13", "DOP853", 1e-13),
            ("scipy_Radau_1e-12", "Radau", 1e-12),
            ("scipy_LSODA_1e-13", "LSODA", 1e-13),
        ):
            if label != "model_890":
                continue
            sol = m.propagate_scipy(model, 0.0, t_cyc, y0, method=meth, rtol=rt)
            end = m.inertial_to_rot_890(c, t_cyc, sol.y[:, -1])
            res[sol_name] = {
                "return_m": float(np.linalg.norm(end[:2] - s890[:2]) * c.a_t * 1e3),
                "return_cm_s": float(np.linalg.norm(end[2:] - s890[2:]) * c.a_t * c.n_t * 1e5),
            }
        log(f"{label}: {json.dumps(res, default=float)[:600]}")
    dump("p2.json", out)


# ------------------------------------------------------------------------------------------- #
# Arc: continuation with checkpoints
# ------------------------------------------------------------------------------------------- #


def arc_paths(tag: str, n: int, variant: str = "") -> tuple[Path, Path]:
    d = OUT / "arcs"
    d.mkdir(parents=True, exist_ok=True)
    stem = f"{tag}_N{n}{variant}"
    return d / f"{stem}.json", d / f"{stem}.npz"


def save_ckpt(
    path_json: Path, path_npz: Path, info: dict[str, Any], st: m.ContinuationState
) -> None:
    arrays: dict[str, Any] = {"times": st.times, "x_cur": st.x_cur}
    if st.x_prev is not None:
        arrays["x_prev"] = st.x_prev
    np.savez(path_npz, **arrays)
    info = dict(info)
    info["state"] = {
        "param": st.param,
        "p_cur": st.p_cur,
        "p_prev": st.p_prev,
        "step": st.step,
        "attempts": st.attempts,
        "status": st.status,
        "signature": st.signature,
    }
    info.setdefault("logs", {})[st.param] = st.log
    info["_meta"] = {
        "task": "#895",
        "git_sha": git_sha(),
        "saved_utc": datetime.now(UTC).isoformat(),
    }
    path_json.write_text(json.dumps(info, indent=1, default=float) + "\n")


def load_ckpt(path_json: Path, path_npz: Path) -> tuple[dict[str, Any], m.ContinuationState]:
    info = json.loads(path_json.read_text())
    z = np.load(path_npz)
    s = info["state"]
    st = m.ContinuationState(
        param=s["param"],
        times=z["times"],
        x_cur=z["x_cur"],
        p_cur=s["p_cur"],
        x_prev=z.get("x_prev", None),
        p_prev=s["p_prev"],
        step=s["step"],
        attempts=s["attempts"],
        status=s["status"],
        signature=[tuple(x) for x in s["signature"]],
        log=info.get("logs", {}).get(s["param"], []),
    )
    return info, st


def stage_arc(args: argparse.Namespace) -> None:
    tag, n = args.epoch, args.cycles
    t_start = time.time()
    pj, pz = arc_paths(tag, n, args.variant)
    table, circles, t_c = epoch_setup(tag)
    times = m.arc_times(circles, t_c, n)
    expected = m.expected_signature(n)
    c890 = m.registry_890()
    model_kw: dict[str, Any] = {}

    def fl_of(model: m.ForceModel, x: np.ndarray) -> list[m.FlybyInfo]:
        return m.arc_flybys(model, times, x, circles.ez)

    def sig_of(model: m.ForceModel, x: np.ndarray) -> list[tuple[str, int, int]]:
        return m.branch_signature(fl_of(model, x))

    def alts(model: m.ForceModel, x: np.ndarray) -> list[tuple[str, float, float]]:
        return [(m.BODY_NAMES[f.body], round(f.t, 1), round(f.alt_km, 1)) for f in fl_of(model, x)]

    if pj.exists():
        info, st = load_ckpt(pj, pz)
        log(f"resumed {pj.name}: {st.param}={st.p_cur:.6g} step {st.step:.4g} status {st.status}")
    else:
        info = {
            "epoch": tag,
            "cycles": n,
            "t_ref_et": table.t_ref_et,
            "t_c_s": t_c,
            "t_c_tdb": m.tdb_string(table.t_ref_et + t_c),
            "circles": circles.to_dict(),
            "n_nodes": len(times),
            "expected_signature": expected,
        }
        if "_D1" in args.variant:
            _, pzp = arc_paths(tag, 1, "_D1periodic")
            zp = np.load(pzp.with_name(pzp.stem + "_final.npz"))
            tt, x1 = m.tile_periodic(zp["times"], zp["x"], zp["rot"], -2, len(times) - 3)
            info["route"] = "D1: tiled periodic orbit of the circular model (note section 3)"
            info["tile_time_mismatch_s"] = float(np.max(np.abs(tt - times)))
            r1 = m.newton(m.homotopy_model(table, circles, 0.0), times, x1)
            info["lam0_from_periodic"] = {
                "converged": r1.converged, "reason": r1.reason, "history": r1.history,
                "flybys": alts(m.homotopy_model(table, circles, 0.0), r1.x) if r1.converged else [],
            }  # fmt: skip
            sig1 = sig_of(m.homotopy_model(table, circles, 0.0), r1.x) if r1.converged else []
            log(f"D1 lam=0 from the tiled periodic orbit: converged {r1.converged}, "
                f"{len(r1.history) - 1} it, signature ok {sig1 == expected}")  # fmt: skip
            status = "running" if (r1.converged and sig1 == expected) else "failed"
            st = m.ContinuationState("lam", times, r1.x, 0.0, signature=sig1, status=status)
            save_ckpt(pj, pz, info, st)
            return stage_arc_loop(args, info, st, pj, pz, t_start, table, circles, t_c, times)
        orb = m.orbit_890(json.loads(REFINE_890.read_text())["refined_state"], c890)
        x0 = m.seed_from_890(orb, circles, t_c, times)
        res = m.newton(m.homotopy_model(table, circles, 0.0, **model_kw), times, x0)
        info["lam0_direct"] = {
            "converged": res.converged,
            "reason": res.reason,
            "history": res.history,
        }
        log(f"lam=0 from the scaled #890 seed: converged {res.converged} ({res.reason})")
        sig = sig_of(m.homotopy_model(table, circles, 0.0), res.x) if res.converged else []
        if res.converged and sig == expected:
            st = m.ContinuationState("lam", times, res.x, 0.0, signature=sig)
            info["route"] = "direct"
        else:
            info["route"] = "fallback: constants homotopy sigma at lam = 0 (pre-registration 1.3)"
            c0 = m.circles_890(c890, circles, t_c)
            xs = m.seed_from_890(orb, c0, t_c, times)
            m0 = m.constants_blend_model(table, c890, c0, circles, 0.0)
            r0 = m.newton(m0, times, xs)
            info["sigma0"] = {"converged": r0.converged, "reason": r0.reason, "history": r0.history}
            sig0 = sig_of(m0, r0.x) if r0.converged else []
            info["sigma0"]["signature"] = sig0
            info["sigma0"]["flybys"] = alts(m0, r0.x)
            log(f"sigma=0: converged {r0.converged}, signature ok {sig0 == expected}")
            if not r0.converged or sig0 != expected:
                st = m.ContinuationState("sigma", times, r0.x, 0.0, status="failed", signature=sig0)
                save_ckpt(pj, pz, info, st)
                log("FAILED at sigma = 0")
                return
            st = m.ContinuationState("sigma", times, r0.x, 0.0, signature=sig0)
        save_ckpt(pj, pz, info, st)

    stage_arc_loop(args, info, st, pj, pz, t_start, table, circles, t_c, times)


def stage_arc_loop(
    args: argparse.Namespace,
    info: dict[str, Any],
    st: m.ContinuationState,
    pj: Path,
    pz: Path,
    t_start: float,
    table: m.EphemerisTable,
    circles: m.Circles,
    t_c: float,
    times: np.ndarray,
) -> None:
    c890 = m.registry_890()
    model_kw: dict[str, Any] = {}

    def fl_of(model: m.ForceModel, x: np.ndarray) -> list[m.FlybyInfo]:
        return m.arc_flybys(model, times, x, circles.ez)

    c0 = m.circles_890(c890, circles, t_c)
    while time.time() - t_start < args.budget:
        if st.status == "done" and st.param == "sigma":
            info["sigma_done"] = {"attempts": st.attempts}
            info.setdefault("logs", {})["sigma"] = st.log
            st = m.ContinuationState("lam", times, st.x_cur, 0.0, signature=st.signature)
            log("sigma continuation done; starting lam")
        if st.status != "running":
            break
        if st.param == "sigma":
            st_model = lambda p: m.constants_blend_model(table, c890, c0, circles, p)  # noqa: E731
        else:
            st_model = lambda p: m.homotopy_model(table, circles, p, **model_kw)  # noqa: E731
        t0 = time.time()
        m.continuation_step(st, st_model, fl_of)
        e = st.log[-1]
        log(
            f"{st.param}={e['value']:.6g} step {e['step']:.4g}: conv {e['converged']} "
            f"it {e['iterations']} acc {e['accepted']} res {e['history'][-1]['max_r_km']:.3g} km "
            f"alts {[a[2] for a in e.get('flybys', [])]} "
            f"{'BRANCH ' + str(e.get('signature')) if e.get('branch_event') else ''}"
            f"({time.time() - t0:.1f} s)"
        )
        save_ckpt(pj, pz, info, st)
    if st.status == "done" and st.param == "lam":
        np.savez(pz.with_name(pz.stem + "_final.npz"), times=st.times, x=st.x_cur)
        log(f"DONE: {pj.name} converged at lam = 1")
    elif st.status == "failed":
        log(f"FAILED: {pj.name} at {st.param} = {st.p_cur:.6g}")
    else:
        log(f"paused: {st.param} = {st.p_cur:.6g}, step {st.step:.4g}")


# ------------------------------------------------------------------------------------------- #
# Deviation D1: periodic orbit of the circular model, continued from #890 in sigma
# ------------------------------------------------------------------------------------------- #


def stage_d1periodic(args: argparse.Namespace) -> None:
    tag = args.epoch
    t_start = time.time()
    pj, pz = arc_paths(tag, 1, "_D1periodic")
    table, circles, t_c = epoch_setup(tag)
    c890 = m.registry_890()
    c0 = m.circles_890(c890, circles, t_c)

    def model_of(p: float) -> m.ForceModel:
        return m.constants_blend_model(table, c890, c0, circles, p)

    def layout_of(p: float) -> tuple[np.ndarray, np.ndarray]:
        return m.periodic_layout(c0, circles, p, t_c)

    if pj.exists():
        info, st = load_ckpt(pj, pz)
        log(f"resumed {pj.name}: sigma={st.p_cur:.6g} step {st.step:.4g} status {st.status}")
    else:
        info = {"epoch": tag, "route": "D1", "t_c_s": t_c, "circles": circles.to_dict()}
        orb = m.orbit_890(json.loads(REFINE_890.read_text())["refined_state"], c890)
        times, rot = layout_of(0.0)
        x0 = m.seed_from_890(orb, c0, t_c, times[:-1])
        m0 = model_of(0.0)
        ev0 = m.periodic_eval(m0, times, x0, rot)
        wrap = float(np.linalg.norm(ev0.jumps[-1, :3]))
        info["sigma0_seed_jumps"] = {
            "max_r_km": ev0.max_r,
            "max_v_kms": ev0.max_v,
            "wrap_r_km": wrap,
        }
        r0 = m.periodic_newton(m0, times, x0, rot)
        fl0 = m.periodic_flybys(m0, times, r0.x, rot, circles.ez)
        info["sigma0"] = {
            "converged": r0.converged, "iterations": len(r0.history) - 1, "history": r0.history,
            "flybys": [(m.BODY_NAMES[f.body], f.t, f.alt_km) for f in fl0],
        }  # fmt: skip
        sig0 = m.branch_signature(fl0)
        ok = r0.converged and len(r0.history) - 1 <= 2 and sig0 == m.expected_signature(1)
        info["sigma0"]["control_pass"] = ok
        log(
            f"sigma=0 control: seed jumps {ev0.max_r:.3g} km "
            f"(wrap {info['sigma0_seed_jumps']['wrap_r_km']:.3g}), "
            f"Newton {len(r0.history) - 1} it, "
            f"conv {r0.converged}, sig ok {sig0 == m.expected_signature(1)}"
        )
        st = m.ContinuationState("sigma", times, r0.x, 0.0, signature=sig0,
                                 status="running" if ok else "failed")  # fmt: skip
        save_ckpt(pj, pz, info, st)
    while time.time() - t_start < args.budget and st.status == "running":
        t0 = time.time()
        m.periodic_continuation_step(st, model_of, layout_of, circles.ez)
        e = st.log[-1]
        log(
            f"sigma={e['value']:.6g} step {e['step']:.4g}: conv {e['converged']} "
            f"it {e['iterations']} acc {e['accepted']} "
            f"alts {[a[2] for a in e.get('flybys', [])]} ({time.time() - t0:.1f} s)"
        )
        save_ckpt(pj, pz, info, st)
    if st.status == "done":
        times, rot = layout_of(1.0)
        fl = m.periodic_flybys(model_of(1.0), times, st.x_cur, rot, circles.ez)
        info["sigma1_flybys"] = [f.to_dict() for f in fl]
        save_ckpt(pj, pz, info, st)
        np.savez(pz.with_name(pz.stem + "_final.npz"), times=times, x=st.x_cur, rot=rot)
        log(
            "DONE periodic at sigma = 1: "
            + ", ".join(f"{m.BODY_NAMES[f.body]} {f.alt_km:.1f}" for f in fl)
        )
    else:
        log(f"{st.status}: sigma = {st.p_cur:.6g}")


# ------------------------------------------------------------------------------------------- #
# Verification: criteria (c) to (f) of the pre-registration
# ------------------------------------------------------------------------------------------- #

VERIFIERS = {
    "DOP853_1e-13": ("DOP853", 1e-13),
    "LSODA_1e-13": ("LSODA", 1e-13),
    "LSODA_1e-12": ("LSODA", 1e-12),
}
A_MEAN_KM = {0: 129900.0, 1: 190900.0, 2: 266000.0, 3: 436300.0, 4: 583500.0}


def load_final(tag: str, n: int, variant: str = "") -> tuple[np.ndarray, np.ndarray]:
    _, pz = arc_paths(tag, n, variant)
    z = np.load(pz.with_name(pz.stem + "_final.npz"))
    return z["times"], z["x"]


def flyby_nodes(n: int) -> list[tuple[int, int]]:
    """(node index of the nominal flyby, moon) for the 2N + 1 flybys."""
    out = []
    for j in range(2 * n + 1):
        out.append((2 + 12 * j, m.I_TITANIA if j % 2 == 0 else m.I_OBERON))
    return out


def stage_verify(args: argparse.Namespace) -> None:
    tag, n = args.epoch, args.cycles
    table, circles, t_c = epoch_setup(tag)
    times, x = load_final(tag, n, args.variant)
    model = m.homotopy_model(table, circles, 1.0)
    k_seg = len(times) - 1
    t_ref = table.t_ref_et
    out: dict[str, Any] = {"epoch": tag, "cycles": n, "t_c_tdb": m.tdb_string(t_ref + t_c)}
    # solver STMs and solver junctions
    ev = m.shoot_eval(model, times, x, stm=True)
    assert ev.stms is not None
    out["solver_junctions"] = {"max_r_m": ev.max_r * 1e3, "max_v_mm_s": ev.max_v * 1e6}
    # (c) junctions under each verifier
    jumps: dict[str, np.ndarray] = {}
    sols: dict[int, Any] = {}
    for name, (meth, rtol) in VERIFIERS.items():
        jj = np.empty((k_seg, 6))
        t0 = time.time()
        for i in range(k_seg):
            sol = m.propagate_scipy(
                model,
                times[i],
                times[i + 1],
                x[i],
                method=meth,
                rtol=rtol,
                dense=name.startswith("DOP"),
            )
            jj[i] = sol.y[:, -1] - x[i + 1]
            if name.startswith("DOP"):
                sols[i] = sol.sol
        jumps[name] = jj
        log(
            f"(c) {name}: max {np.linalg.norm(jj[:, :3], axis=1).max() * 1e3:.3f} m, "
            f"{np.linalg.norm(jj[:, 3:], axis=1).max() * 1e6:.4f} mm/s ({time.time() - t0:.0f} s)"
        )
    out["c_junctions"] = {
        name: {
            "max_r_m": float(np.linalg.norm(jj[:, :3], axis=1).max() * 1e3),
            "max_v_mm_s": float(np.linalg.norm(jj[:, 3:], axis=1).max() * 1e6),
            "per_junction_r_m": (np.linalg.norm(jj[:, :3], axis=1) * 1e3).tolist(),
        }
        for name, jj in jumps.items()
    }
    lsoda_self = np.linalg.norm(jumps["LSODA_1e-13"][:, :3] - jumps["LSODA_1e-12"][:, :3], axis=1)
    out["c_lsoda_self_convergence_max_m"] = float(lsoda_self.max() * 1e3)
    c_pass = (
        all(
            out["c_junctions"][nm]["max_r_m"] < 1.0 and out["c_junctions"][nm]["max_v_mm_s"] < 1.0
            for nm in ("DOP853_1e-13", "LSODA_1e-13")
        )
        and out["c_lsoda_self_convergence_max_m"] <= 0.3
    )
    out["c_pass"] = c_pass

    # (d) every flyby integrated through
    def stm_between(a: int, b: int) -> np.ndarray:
        assert ev.stms is not None
        phi = np.eye(6)
        for i in range(a, b):
            phi = ev.stms[i] @ phi
        return phi

    d_rows = []
    for idx, moon in flyby_nodes(n):
        a, b = max(0, idx - 6), min(k_seg, idx + 6)
        row: dict[str, Any] = {"node": idx, "moon": m.BODY_NAMES[moon], "from": a, "to": b,
                               "span_d": (times[b] - times[a]) / m.DAY_S}  # fmt: skip
        ends = {}
        for name in ("DOP853_1e-13", "LSODA_1e-13"):
            meth, rtol = VERIFIERS[name]
            sol = m.propagate_scipy(model, times[a], times[b], x[a], method=meth, rtol=rtol)
            land = sol.y[:, -1] - x[b]
            ends[name] = sol.y[:, -1]
            pred = np.zeros(6)
            for kk in range(a + 1, b + 1):
                pred += stm_between(kk, b) @ jumps[name][kk - 1]
            row[name] = {
                "landing_r_km": float(np.linalg.norm(land[:3])),
                "landing_v_cm_s": float(np.linalg.norm(land[3:]) * 1e5),
                "pred_r_km": float(np.linalg.norm(pred[:3])),
                "pred_v_cm_s": float(np.linalg.norm(pred[3:]) * 1e5),
                "landing_minus_pred_r_km": float(np.linalg.norm(land[:3] - pred[:3])),
            }
        diff = float(np.linalg.norm(ends["DOP853_1e-13"][:3] - ends["LSODA_1e-13"][:3]))
        row["integrators_end_diff_km"] = diff
        ok = True
        posed = True
        for name in ("DOP853_1e-13", "LSODA_1e-13"):
            r = row[name]
            posed &= r["pred_r_km"] <= 0.3 and r["pred_v_cm_s"] <= 0.3
            ok &= r["landing_r_km"] <= 1.0 and r["landing_v_cm_s"] <= 1.0
            ok &= r["landing_r_km"] <= r["pred_r_km"] + 3 * diff + 0.010
        row["well_posed"] = bool(posed)
        row["pass"] = bool(ok and posed)
        d_rows.append(row)
        log(
            f"(d) node {idx} {row['moon']}: DOP {row['DOP853_1e-13']['landing_r_km']:.3f} km "
            f"(pred {row['DOP853_1e-13']['pred_r_km']:.3f}), "
            f"LSODA {row['LSODA_1e-13']['landing_r_km']:.3f} km, "
            f"diff {diff:.3f} km, pass {row['pass']}"
        )
    out["d_flyby_through"] = d_rows
    out["d_pass"] = all(r["pass"] for r in d_rows)

    # (e), (f) encounters on the verified DOP853 trajectory, moons from SPICE
    def traj(t: float) -> np.ndarray:
        i = int(np.searchsorted(times, t, side="right") - 1)
        i = min(max(i, 0), k_seg - 1)
        return np.asarray(sols[i](t))

    allmin: dict[str, list[dict[str, Any]]] = {}
    flybys = []
    for k in m.MOONS:

        def bstate(t: float, k: int = k) -> tuple[np.ndarray, np.ndarray]:
            s6 = m.spice_state(m.NAIF_IDS[k], t_ref + t)
            return s6[:3], s6[3:]

        encs = m.find_minima(traj, bstate, float(times[0]), float(times[-1]), dt=1800.0, body=k)
        rh = m.hill_radius_km(k, A_MEAN_KM[k])
        rows: list[dict[str, Any]] = []
        for e in encs:
            o = e.osculating()
            pm, _ = bstate(e.t)
            rows.append({
                "t_s": e.t, "tdb": m.tdb_string(t_ref + e.t), "dist_km": e.dist_km,
                "hill_radii": e.dist_km / rh, "side": 1 if float(np.dot(e.rel_pos, pm)) > 0 else -1,
                "sense": 1 if float(np.dot(np.cross(e.rel_pos, e.rel_vel), circles.ez)) > 0 else -1,
                "z_ref_km": float(np.dot(traj(e.t)[:3], circles.ez)),
                "moon_z_ref_km": float(np.dot(pm, circles.ez)),
                **o,
            })  # fmt: skip
        allmin[m.BODY_NAMES[k]] = rows
    for idx, moon in flyby_nodes(n):
        tn = times[idx]
        cands = [r for r in allmin[m.BODY_NAMES[moon]] if abs(r["t_s"] - tn) < 3 * m.DAY_S]
        cands.sort(key=lambda r: r["dist_km"])
        f = cands[0] if cands else None
        flybys.append({"node": idx, "moon": m.BODY_NAMES[moon], **(f or {"missing": True})})
    e_pass = all(
        (not f.get("missing")) and f["altitude_km"] >= 50.0 and f["energy_km2s2"] > 0
        for f in flybys
    )
    targeted = {(f["moon"], round(f["t_s"], 3)) for f in flybys if not f.get("missing")}
    others = []
    closest = {}
    for name, rows in allmin.items():
        if rows:
            closest[name] = min(r["dist_km"] for r in rows)
        for r in rows:
            if (name, round(r["t_s"], 3)) in targeted:
                continue
            if r["hill_radii"] <= 2.0:
                others.append({"moon": name, **r})

    # Uranus: minima of |r|
    def ustate(t: float) -> tuple[np.ndarray, np.ndarray]:
        return np.zeros(3), np.zeros(3)

    umin = m.find_minima(traj, ustate, float(times[0]), float(times[-1]), dt=3600.0, body=-1)
    ts = np.linspace(times[0], times[-1], 20000)
    rr = np.array([np.linalg.norm(traj(t)[:3]) for t in ts])
    closest["Uranus"] = float(min([e.dist_km for e in umin] + [float(rr.min())]))
    out["uranus_distance_range_km"] = [float(rr.min()), float(rr.max())]
    f_pass = (
        all(r["altitude_km"] >= 50.0 for rows in allmin.values() for r in rows)
        and closest["Uranus"] >= m.R_REF_KM + 1000.0
    )
    out["e_flybys"] = flybys
    out["e_pass"] = e_pass
    out["f_other_within_2_hill"] = others
    out["f_closest_km"] = closest
    out["f_pass"] = f_pass
    out["all_minima"] = allmin
    out["pass_c_to_f"] = bool(c_pass and out["d_pass"] and e_pass and f_pass)
    for f in flybys:
        if not f.get("missing"):
            log(
                f"(e) {f['moon']:8s} {f['tdb']} alt {f['altitude_km']:8.1f} km "
                f"v {f['speed_rel_kms']:.4f} "
                f"vinf {f['vinf_kms']:.4f} e {f['ecc']:.3f} side {f['side']} z {f['z_ref_km']:.0f}"
            )
    log(f"(f) others within 2 Hill radii: {len(others)}; closest {json.dumps(closest)}")
    log(f"c {c_pass} d {out['d_pass']} e {e_pass} f {f_pass}")
    dump(f"verify_{tag}_N{n}{args.variant}.json", out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("stage", choices=["control", "p2", "arc", "verify", "d1periodic"])
    ap.add_argument("--epoch", default="E1", choices=sorted(EPOCHS))
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--variant", default="")
    ap.add_argument("--budget", type=float, default=420.0)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if "_D2" in args.variant:
        m.NEWTON_TOL.update({"r": 5e-5, "v": 5e-8})
        log("deviation D2: Newton tolerance 5 cm and 5e-8 km/s")
    t = time.time()
    stages = {
        "control": stage_control,
        "p2": stage_p2,
        "arc": stage_arc,
        "verify": stage_verify,
        "d1periodic": stage_d1periodic,
    }
    stages[args.stage](args)
    log(f"done in {time.time() - t:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
