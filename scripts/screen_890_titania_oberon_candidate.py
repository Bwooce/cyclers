# ruff: noqa: E501
"""#890: is the Titania-Oberon-Titania closure a periodic orbit of the four-body model?

Pre-registration: docs/notes/2026-10-04-890-titania-oberon-candidate.md section 1.

Stages (each writes data/found/890_titania_oberon_candidate/<stage>.json[l]):

    uv run python scripts/screen_890_titania_oberon_candidate.py --stage rebuild
    uv run python scripts/screen_890_titania_oberon_candidate.py --stage decisive
    uv run python scripts/screen_890_titania_oberon_candidate.py --stage direct
    uv run python scripts/screen_890_titania_oberon_candidate.py --stage continue [--max-steps N]
    uv run python scripts/screen_890_titania_oberon_candidate.py --stage verify
    uv run python scripts/screen_890_titania_oberon_candidate.py --stage sensitivity
    uv run python scripts/screen_890_titania_oberon_candidate.py --stage variants

``continue`` checkpoints every converged lam to continuation.jsonl and resumes
from the last line, so it can be run repeatedly under a wall-clock budget.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.search import two_moon_periodic_890 as tm
from cyclerfinder.verify.turn_gate import Encounter, demanded_turn_gate
from cyclerfinder.verify.turn_gate_closures import symmetric_closure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "found" / "890_titania_oberon_candidate"
STORED = ROOT / "data" / "found" / "888_turn_gate" / "candidates.jsonl"
N_INTERIOR = 11
LAM_LADDER = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
T0 = time.time()


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{stamp} +{time.time() - T0:7.1f}s] {msg}", flush=True)


def git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:  # pragma: no cover
        return "unknown"


def jsonable(x: Any) -> Any:
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        f = float(x)
        return f if math.isfinite(f) else str(f)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def write_json(name: str, obj: dict[str, Any]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    obj = {"_meta": {"task": "#890", "git_sha": git_sha(), "stage": name}, **obj}
    (OUT / name).write_text(json.dumps(jsonable(obj), indent=1) + "\n")
    log(f"wrote {OUT / name}")


def stored_record() -> dict[str, Any]:
    for line in STORED.read_text().splitlines():
        g = json.loads(line)
        if g.get("kind") == "candidate" and g["sequence"] == ["Titania", "Oberon", "Titania"]:
            return g
    raise RuntimeError("stored candidate not found")


def closure() -> Any:
    g = stored_record()
    c = symmetric_closure(
        "Uranus",
        "Titania",
        "Oberon",
        tof_days=g["tof_days_per_leg"],
        n_rev=(g["n_rev"][0], g["n_rev"][1]),
        rel_offset_deg=g["rel_offset_deg"],
        branches=tuple(g["branches"]),  # type: ignore[arg-type]
    )
    assert c is not None
    return c


def conic_elements(r: np.ndarray, v: np.ndarray, gm: float) -> dict[str, float]:
    rn = float(np.linalg.norm(r))
    energy = 0.5 * float(v @ v) - gm / rn
    a = -gm / (2.0 * energy)
    h = float(r[0] * v[1] - r[1] * v[0])
    evec = (float(v @ v) - gm / rn) * r / gm - float(r @ v) * v / gm
    e = float(np.linalg.norm(evec))
    return {
        "a_km": a,
        "e": e,
        "rp_km": a * (1 - e),
        "ra_km": a * (1 + e),
        "period_days": 2 * math.pi * math.sqrt(a**3 / gm) / 86400.0,
        "argp_deg": math.degrees(math.atan2(evec[1], evec[0])),
        "h": h,
    }


# ---------------------------------------------------------------------------
# rebuild
# ---------------------------------------------------------------------------


def stage_rebuild() -> None:
    g = stored_record()
    c = closure()
    geom = tm.closure_geometry(c)
    rep = demanded_turn_gate(c.encounters)
    el = conic_elements(geom.r_dep, geom.v_dep, geom.gm_sys)
    vinf = list(c.stored_convention_vinf)
    diff_vinf = max(abs(a - b) for a, b in zip(vinf, g["vinf_kms"], strict=True))
    turns = {e.body: e.demanded_turn_deg for e in rep.encounters}
    req_alt = {e.body: e.required_alt_km for e in rep.encounters}
    log(f"vinf {vinf} (stored diff {diff_vinf:.2e}); conic {el}")
    log(f"turns {turns}; required altitudes {req_alt}")
    # frame checks
    m0 = tm.TwoMoonModel(lam=0.0)
    k1 = tm.kepler_propagate(geom.r_dep, geom.v_dep, 86400.0, geom.gm_sys)
    t1 = m0.tau(1.0)
    s = tm.rot_from_inertial(m0, t1, k1[:2], k1[2:])
    r_rt, v_rt = tm.inertial_from_rot(m0, t1, s)
    t_end = m0.tau(geom.tof_days)
    a = tm.propagate(m0, s, t1, t_end)
    r, v = tm.inertial_from_rot(m0, t_end, a.state)
    k = tm.kepler_propagate(k1[:2], k1[2:], (geom.tof_days - 1.0) * 86400.0, geom.gm_sys)
    ro, _ = m0.moon_inertial("Oberon", t_end)
    m1 = tm.TwoMoonModel(lam=1.0)
    th_t = m1.system.omega_gan * 2.5 * m1.forcing_period
    frame = {
        "round_trip_pos_km": float(np.max(np.abs(r_rt - k1[:2]))),
        "round_trip_vel_kms": float(np.max(np.abs(v_rt - k1[2:]))),
        "lam0_leg_pos_err_rel": float(np.max(np.abs(r - k[:2]))) / m0.length_km,
        "lam0_leg_vel_err_rel": float(np.max(np.abs(v - k[2:]))) / m0.vel_unit_kms,
        "lam0_leg_end_to_oberon_km": float(np.linalg.norm(r - ro)),
        "lam0_half_period_days": m0.days(2.5 * m0.forcing_period),
        "lam1_forcing_period_days": m1.days(m1.forcing_period),
        "lam1_cycle_days": m1.days(5 * m1.forcing_period),
        "lam1_oberon_angle_at_T_deg": math.degrees(th_t % (2 * math.pi)),
    }
    log(f"frame checks {frame}")
    # leg moon distances (local minima) along the conic, both moons
    write_json(
        "rebuild.json",
        {
            "stored_vinf_kms": g["vinf_kms"],
            "rebuilt_vinf_kms": vinf,
            "max_vinf_diff_kms": diff_vinf,
            "residual_kms": c.residual_kms,
            "conic_leg0": el,
            "gate_project_floor": rep.as_dict(),
            "local_vinf_kms": {
                "titania_in": geom.base_in_local,
                "titania_out": geom.base_out_local,
                "oberon_in": geom.pert_in_local,
                "oberon_out": geom.pert_out_local,
            },
            "frame_checks": frame,
        },
    )


# ---------------------------------------------------------------------------
# helpers on trajectories
# ---------------------------------------------------------------------------


def full_cycle(model: tm.TwoMoonModel, s0: np.ndarray, tau0: float, method: str = "DOP853") -> Any:
    rtol = 1e-13 if method == "DOP853" else 1e-12
    return tm.propagate(
        model, s0, tau0, tau0 + 5 * model.forcing_period, method=method, rtol=rtol, atol=1e-13,
        dense=True,
    )  # fmt: skip


def encounter_table(model: tm.TwoMoonModel, sol: Any) -> list[dict[str, Any]]:
    out = []
    for body in ("Titania", "Oberon"):
        within = 2.0 * tm.hill_radius_km("Uranus", body)
        out += tm.find_encounters(model, sol, bodies=[body], within_km=within, n_samples=60000)
    out.sort(key=lambda e: e["tau"])
    return out


def closure_error(
    model: tm.TwoMoonModel, s_start: np.ndarray, s_end: np.ndarray
) -> dict[str, float]:
    d = s_end - s_start
    return {
        "pos_km": float(np.linalg.norm(d[:2])) * model.length_km,
        "vel_kms": float(np.linalg.norm(d[2:])) * model.vel_unit_kms,
    }


def fmt_enc(e: dict[str, Any]) -> str:
    return (
        f"{e['body']:8s} t {e['t_days']:8.3f} d  d {e['dist_km']:10.1f} km "
        f"({e['hill_radii']:.2f} R_H)  v_rel {e['rel_speed_kms']:.4f}  e_osc {e['osc_ecc']:.3f}"
    )


# ---------------------------------------------------------------------------
# decisive
# ---------------------------------------------------------------------------


def stage_decisive() -> None:
    c = closure()
    geom = tm.closure_geometry(c)
    m = tm.TwoMoonModel(lam=1.0)
    res: dict[str, Any] = {}
    # (a) literal stored conic from the leg-0 apoapsis nearest mid-leg
    gm = geom.gm_sys
    el = conic_elements(geom.r_dep, geom.v_dep, gm)
    ts = np.linspace(0, geom.tof_days * 86400.0, 20001)
    rr = []
    for t in ts[::50]:
        st = tm.kepler_propagate(geom.r_dep, geom.v_dep, float(t), gm)
        rr.append(float(np.linalg.norm(st[:2])))
    rr_a = np.array(rr)
    tt = ts[::50]
    peaks = [
        i for i in range(1, len(rr_a) - 1) if rr_a[i] >= rr_a[i - 1] and rr_a[i] >= rr_a[i + 1]
    ]
    mid = geom.tof_days * 86400.0 / 2
    i_apo = min(peaks, key=lambda i: abs(tt[i] - mid))
    # refine apoapsis time by golden search on r
    lo, hi = tt[i_apo - 1], tt[i_apo + 1]
    for _ in range(80):
        m1_, m2_ = lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618
        r1 = np.linalg.norm(tm.kepler_propagate(geom.r_dep, geom.v_dep, m1_, gm)[:2])
        r2 = np.linalg.norm(tm.kepler_propagate(geom.r_dep, geom.v_dep, m2_, gm)[:2])
        if r1 > r2:
            hi = m2_
        else:
            lo = m1_
    t_apo = 0.5 * (lo + hi)
    st = tm.kepler_propagate(geom.r_dep, geom.v_dep, t_apo, gm)
    tau_apo = m.tau(t_apo / 86400.0)
    rp, vp = tm.planet_offset_inertial(m, tau_apo)
    s_apo = tm.rot_from_inertial(m, tau_apo, st[:2] + rp, st[2:] + vp)
    log(f"(a) apoapsis at t = {t_apo / 86400:.3f} d, r = {np.linalg.norm(st[:2]):.1f} km")
    for label, s0, tau0 in (("a_literal_conic_from_apoapsis", s_apo, tau_apo),):
        sol = full_cycle(m, s0, tau0).sol
        enc = encounter_table(m, sol)
        for e in enc:
            log("   " + fmt_enc(e))
        ce = closure_error(m, s0, sol.y[:4, -1])
        log(f"   closure after one cycle: {ce}")
        res[label] = {"tau0": tau0, "t0_days": m.days(tau0), "state0": s0, "encounters": enc,
                      "closure": ce, "conic": el}  # fmt: skip
    # (b) hyperbola periapsis at Titania (honest patched-conic realisation)
    s_b, hyp = tm.flyby_node_state(m, "Titania", 0.0, geom.base_in_local, geom.base_out_local)
    log(f"(b) Titania hyperbola: rp {hyp.rp_km:.1f} km, turn {math.degrees(hyp.turn_rad):.2f} deg, "
        f"state {s_b}")  # fmt: skip
    s_b0 = np.array([s_b[0], 0.0, 0.0, s_b[3]])
    sol = full_cycle(m, s_b0, 0.0).sol
    enc = encounter_table(m, sol)
    for e in enc:
        log("   " + fmt_enc(e))
    ce = closure_error(m, s_b0, sol.y[:4, -1])
    log(f"   closure after one cycle: {ce}")
    # where is the trajectory at tau = T relative to the intended Oberon periapsis?
    half = 2.5 * m.forcing_period
    s_t, _hyp_o = tm.flyby_node_state(m, "Oberon", half, geom.pert_in_local, geom.pert_out_local)
    at_t = sol.sol(half)[:4]
    res["b_titania_periapsis"] = {
        "state0": s_b0,
        "dropped_y_vx": [s_b[1], s_b[2]],
        "hyperbola": {"rp_km": hyp.rp_km, "turn_deg": math.degrees(hyp.turn_rad)},
        "encounters": enc,
        "closure": ce,
        "miss_of_intended_oberon_periapsis_at_T": closure_error(m, s_t, at_t),
    }
    miss = res["b_titania_periapsis"]["miss_of_intended_oberon_periapsis_at_T"]
    log(f"   miss at T of intended Oberon periapsis: {miss}")
    write_json("decisive.json", res)


# ---------------------------------------------------------------------------
# correction
# ---------------------------------------------------------------------------


def run_newton(model: tm.TwoMoonModel, sh: tm.SymmetricShooter, z0: np.ndarray, **kw: Any) -> Any:
    return sh.newton(z0, log=log, **kw)


def flyby_summary(model: tm.TwoMoonModel, z: np.ndarray, sh: tm.SymmetricShooter) -> dict[str, Any]:
    s0, _mids, s_end = sh.unpack(z)
    out: dict[str, Any] = {}
    for label, body, s, tau in (("titania", "Titania", s0, 0.0), ("oberon", "Oberon", s_end, sh.T)):
        dr, dv = tm.relative_inertial(model, body, tau, s)
        gm = model.gm_base if body == "Titania" else model.gm_pert
        osc = tm.osculating_flyby(gm, dr, dv)
        rb, _ = model.moon_inertial(body, tau)
        radial_sign = float(np.sign(dr @ rb))
        out[label] = {
            "rp_km": float(np.linalg.norm(dr)),
            "rp_over_lam_km": float(np.linalg.norm(dr)) / model.lam,
            "periapsis_side": "outside" if radial_sign > 0 else "inside",
            "vp_kms": float(np.linalg.norm(dv)),
            "turn_sense": float(np.sign(dr[0] * dv[1] - dr[1] * dv[0])),
            **{f"osc_{k}": v for k, v in osc.items()},
        }
    return out


def stage_direct() -> None:
    c = closure()
    geom = tm.closure_geometry(c)
    m = tm.TwoMoonModel(lam=1.0)
    sh, z0, info = tm.symmetric_guess(m, geom, N_INTERIOR)
    log(f"direct lam=1 guess info {info}")
    r = run_newton(m, sh, z0, max_iter=30)
    log(f"direct: {r.message}, final {r.final_residual:.3e}")
    obj: dict[str, Any] = {"guess_info": info, "converged": r.converged, "message": r.message,
                           "residual_history": r.residual_history, "z": r.z}  # fmt: skip
    if r.converged:
        obj["flybys"] = flyby_summary(m, r.z, sh)
        log(f"   flybys {obj['flybys']}")
    write_json("direct.json", obj)


def load_continuation() -> list[dict[str, Any]]:
    p = OUT / "continuation.jsonl"
    if not p.exists():
        return []
    return [json.loads(x) for x in p.read_text().splitlines() if x.strip()]


def append_continuation(row: dict[str, Any]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "continuation.jsonl").open("a") as fh:
        fh.write(json.dumps(jsonable(row)) + "\n")
        fh.flush()


def stage_continue(max_steps: int, budget_s: float) -> None:
    c = closure()
    geom = tm.closure_geometry(c)
    rows = load_continuation()
    if not rows:
        lam = LAM_LADDER[0]
        m = tm.TwoMoonModel(lam=lam)
        sh, z0, info = tm.symmetric_guess(m, geom, N_INTERIOR)
        log(f"lam {lam}: guess {info}")
        r = run_newton(m, sh, z0, max_iter=30, tol=1e-10)
        row = {"kind": "step", "lam": lam, "converged": r.converged, "message": r.message,
               "residual_history": r.residual_history, "final_residual": r.final_residual,
               "jac_det_sign": r.jac_logdet_sign, "z": r.z, "guess_info": info}  # fmt: skip
        if r.converged:
            row["flybys"] = flyby_summary(m, r.z, sh)
            log(f"   flybys {row['flybys']}")
        append_continuation(row)
        if not r.converged:
            log("positive control FAILED at the first lam; stopping")
            return
        rows = load_continuation()
    last = [x for x in rows if x["kind"] == "step" and x["converged"]][-1]
    lam_prev = float(last["lam"])
    z_prev = np.array(last["z"])
    sign_prev = float(last["jac_det_sign"])
    targets = [x for x in LAM_LADDER if x > lam_prev + 1e-15]
    step_count = 0
    lam_target = targets[0] if targets else None
    while lam_target is not None and step_count < max_steps and time.time() - T0 < budget_s:
        m = tm.TwoMoonModel(lam=lam_target)
        sh = tm.SymmetricShooter(m, 2.5 * m.forcing_period, N_INTERIOR)
        m_prev = tm.TwoMoonModel(lam=lam_prev)
        zg = tm.rescale_flyby_nodes(z_prev, m_prev, m, N_INTERIOR)
        log(f"lam {lam_target}: from lam {lam_prev} (flyby offsets scaled by lam)")
        r = run_newton(m, sh, zg, max_iter=20, tol=1e-10)
        row = {"kind": "step", "lam": lam_target, "converged": r.converged, "message": r.message,
               "residual_history": r.residual_history, "final_residual": r.final_residual,
               "jac_det_sign": r.jac_logdet_sign, "z": r.z}  # fmt: skip
        step_count += 1
        if r.converged:
            row["flybys"] = flyby_summary(m, r.z, sh)
            fb = row["flybys"]
            ti, ob = fb["titania"], fb["oberon"]
            log(
                f"   converged; Titania rp {ti['rp_km']:.2f} km ({ti['periapsis_side']}, "
                f"turn {ti['osc_turn_deg']:.2f}), Oberon rp {ob['rp_km']:.2f} km "
                f"({ob['periapsis_side']}, turn {ob['osc_turn_deg']:.2f}); "
                f"det sign {r.jac_logdet_sign:+.0f}"
            )
            if r.jac_logdet_sign != sign_prev:
                log("   det J changed sign: a fold or a bifurcation lies between the two lam")
                row["det_sign_change"] = True
            append_continuation(row)
            z_prev, lam_prev, sign_prev = r.z, lam_target, r.jac_logdet_sign
            targets = [x for x in LAM_LADDER if x > lam_prev + 1e-15]
            lam_target = targets[0] if targets else None
        else:
            append_continuation({**row, "kind": "failed_step"})
            dl = lam_target - lam_prev
            if dl / lam_prev < 1e-4:
                log(f"   step collapsed at lam {lam_prev}: stopping (fold or failure)")
                append_continuation({"kind": "stop", "lam": lam_prev, "reason": "step collapse"})
                return
            lam_target = lam_prev + 0.5 * dl
    log(f"continuation paused at lam {lam_prev} after {step_count} steps")


# ---------------------------------------------------------------------------
# verify (pre-registered criteria V1 to V6, note section 1.6)
# ---------------------------------------------------------------------------

KM = 1.0
CM_S = 1e-5  # km/s


def converged_solution(lam: float = 1.0) -> tuple[tm.TwoMoonModel, tm.SymmetricShooter, np.ndarray]:
    rows = [r for r in load_continuation() if r["kind"] == "step" and r["converged"]]
    row = [r for r in rows if abs(float(r["lam"]) - lam) < 1e-12][-1]
    m = tm.TwoMoonModel(lam=lam)
    sh = tm.SymmetricShooter(m, 2.5 * m.forcing_period, N_INTERIOR)
    return m, sh, np.array(row["z"])


def apoapsis_times(model: tm.TwoMoonModel, sol: Any, n: int = 20000) -> list[float]:
    ts = np.linspace(sol.t[0], sol.t[-1], n)
    ys = sol.sol(ts)
    sy = model.system
    r = np.hypot(ys[0] + sy.mu, ys[1])
    idx = np.nonzero((r[1:-1] > r[:-2]) & (r[1:-1] > r[2:]))[0] + 1
    return [float(ts[i]) for i in idx]


def stage_verify() -> None:
    m, sh, z = converged_solution(1.0)
    res, _ = sh.evaluate(z, with_jac=False)
    v1 = float(np.max(np.abs(res)))
    log(f"V1 multiple-shooting residual {v1:.3e}")
    s0, _mids, _s_end = sh.unpack(z)
    period = 5 * m.forcing_period
    out: dict[str, Any] = {"z": z, "V1_residual": v1, "V1_pass": v1 <= 1e-10}
    # monodromy from one STM integration over the full cycle
    arc = tm.propagate(m, s0, 0.0, period, with_stm=True)
    assert arc.stm is not None
    eig = np.linalg.eigvals(arc.stm)
    lam_max = float(np.max(np.abs(eig)))
    log(f"monodromy eigenvalues {eig}; |Lambda|max {lam_max:.4e}; det {np.linalg.det(arc.stm):.6f}")
    out["monodromy_eigenvalues"] = [[float(e.real), float(e.imag)] for e in eig]
    out["monodromy_max_abs"] = lam_max
    out["monodromy_det"] = float(np.linalg.det(arc.stm))
    runs: dict[str, Any] = {}
    sols: dict[str, Any] = {}
    for method in ("DOP853", "Radau"):
        t1 = time.time()
        sol = full_cycle(m, s0, 0.0, method=method).sol
        sols[method] = sol
        end = sol.y[:4, -1]
        half = sol.sol(sh.T)[:4]
        ce = closure_error(m, s0, end)
        half_sym = {
            "y_km": abs(float(half[1])) * m.length_km,
            "vx_kms": abs(float(half[2])) * m.vel_unit_kms,
        }
        full_nd = float(np.max(np.abs(end - s0)))
        direct = ce["pos_km"] <= 1 * KM and ce["vel_kms"] <= 1 * CM_S
        fallback = (
            lam_max > 1e6
            and half_sym["y_km"] <= 1 * KM
            and half_sym["vx_kms"] <= 1 * CM_S
            and full_nd <= 10 * lam_max * (1e-10 + 1e-12)
        )
        runs[method] = {
            "closure": ce,
            "closure_nondim_max": full_nd,
            "half_cycle_symmetry": half_sym,
            "pass_direct": direct,
            "pass_fallback": fallback,
            "pass": bool(direct or fallback),
            "nfev": int(sol.nfev),
            "wall_s": time.time() - t1,
        }
        log(f"{method}: closure {ce}, half-cycle {half_sym}, pass {direct or fallback}")
    # DOP853 vs Radau at tau = T (sol.sol on dense output)
    a, b = sols["DOP853"].sol(sh.T)[:4], sols["Radau"].sol(sh.T)[:4]
    agree = closure_error(m, a, b)
    runs["DOP853_vs_Radau_at_T"] = agree
    v2 = runs["DOP853"]["pass"]
    v3 = runs["Radau"]["pass"] and agree["pos_km"] <= 1 * KM and agree["vel_kms"] <= 1 * CM_S
    out["V2_closure"] = runs["DOP853"]
    out["V3_closure_radau"] = runs["Radau"]
    out["V3_agreement_at_T"] = agree
    out["V2_pass"], out["V3_pass"] = bool(v2), bool(v3)
    log(f"V2 {v2}; V3 {v3} (agreement at T {agree})")
    # encounters
    sol = sols["DOP853"]
    enc = encounter_table(m, sol)
    targeted_taus = {"Titania": [0.0, period], "Oberon": [sh.T]}
    tol_t = 0.05
    targeted, other = [], []
    for e in enc:
        is_t = any(abs(e["tau"] - t) < tol_t for t in targeted_taus[e["body"]])
        (targeted if is_t else other).append(e)
        log(("   T " if is_t else "   * ") + fmt_enc(e))
    out["encounters_targeted"] = targeted
    out["encounters_other_within_2RH"] = other
    # V4
    v4 = []
    for e in targeted:
        soi = tm.laplace_soi_km("Uranus", e["body"])
        ok = e["alt_km"] >= 50.0 and e["dist_km"] <= soi and e["osc_ecc"] > 1.0
        v4.append({"body": e["body"], "t_days": e["t_days"], "alt_km": e["alt_km"],
                   "dist_km": e["dist_km"], "soi_km": soi, "osc_ecc": e["osc_ecc"],
                   "osc_vinf_kms": e["osc_vinf_kms"], "osc_turn_deg": e["osc_turn_deg"],
                   "pass": ok})  # fmt: skip
    out["V4"] = v4
    out["V4_pass"] = bool(all(x["pass"] for x in v4) and len(v4) == 3)
    # V5: gate with SOI-crossing V-infinity vectors; and osculating asymptotes (not independent)
    encs_soi, encs_osc, soi_rows = [], [], []
    for e in targeted:
        body = e["body"]
        if e["tau"] < 0.5 or e["tau"] > period - 0.5:
            # the Titania flyby straddles the cycle boundary: integrate a window around tau = 0
            w = tm.propagate(m, s0, 0.0, 1.5, dense=True).sol
            wb = tm.propagate(m, s0, 0.0, -1.5, dense=True).sol
            cross_out = tm.soi_crossing_velocities(
                m, w, body, 0.0, tm.laplace_soi_km("Uranus", body)
            )
            cross_in = tm.soi_crossing_velocities(
                m, wb, body, 0.0, tm.laplace_soi_km("Uranus", body)
            )
            cin, cout = cross_in["in"], cross_out["out"]
            if e["tau"] > 0.5:
                continue  # the same flyby, counted once
        else:
            cr = tm.soi_crossing_velocities(
                m, sol, body, e["tau"], tm.laplace_soi_km("Uranus", body)
            )
            cin, cout = cr["in"], cr["out"]
        if cin is None or cout is None:
            soi_rows.append({"body": body, "error": "no SOI crossing found"})
            continue
        encs_soi.append(Encounter.for_body(body, cin["dv_kms"], cout["dv_kms"], label="soi"))
        tp = e["tau"] if e["tau"] < period - 0.5 else 0.0
        sp = s0 if tp == 0.0 else sol.sol(tp)[:4]
        dr, dv = tm.relative_inertial(m, body, tp, sp)
        gm = m.gm_base if body == "Titania" else m.gm_pert
        a_in, a_out = tm.osculating_asymptotes(gm, dr, dv)
        encs_osc.append(Encounter.for_body(body, a_in, a_out, label="osculating"))
        soi_rows.append({"body": body, "soi_in_t_days": cin["t_days"], "soi_out_t_days": cout["t_days"],
                         "soi_in_vinf_kms": float(np.linalg.norm(cin["dv_kms"])),
                         "soi_out_vinf_kms": float(np.linalg.norm(cout["dv_kms"])),
                         "soi_crossing_duration_days": cout["t_days"] - cin["t_days"]})  # fmt: skip
    rep_soi = demanded_turn_gate(encs_soi)
    rep_osc = demanded_turn_gate(encs_osc)
    out["V5_soi_rows"] = soi_rows
    out["V5_gate_soi"] = rep_soi.as_dict()
    out["V5_gate_osculating_not_independent"] = rep_osc.as_dict()
    out["V5_pass"] = bool(rep_soi.turn_feasible and len(encs_soi) == 2)
    for et in rep_soi.encounters:
        log(f"V5 SOI gate {et.body}: in {et.vinf_in_kms:.4f} out {et.vinf_out_kms:.4f} demanded "
            f"{et.demanded_turn_deg:.2f} avail {et.available_bend_deg:.2f} ratio {et.ratio:.3f}")  # fmt: skip
    # V6
    r_pl = np.hypot(sol.y[0] + m.system.mu, sol.y[1]) * m.length_km
    out["min_dist_uranus_km"] = float(r_pl.min())
    out["max_dist_uranus_km"] = float(r_pl.max())
    out["V6_pass"] = len(other) == 0
    # identity: apoapses between flybys, turn sense, sides
    apos = apoapsis_times(m, sol)
    leg0 = [t for t in apos if 0 < t < sh.T]
    leg1 = [t for t in apos if sh.T < t < period]
    out["apoapses_per_leg"] = [len(leg0), len(leg1)]
    out["flybys"] = flyby_summary(m, z, sh)
    log(f"apoapses per leg {out['apoapses_per_leg']}; flybys {out['flybys']}")
    write_json("verify.json", out)


# ---------------------------------------------------------------------------
# sensitivity (V7, reported)
# ---------------------------------------------------------------------------


def next_periapsis(model: tm.TwoMoonModel, s: np.ndarray, tau0: float, body: str,
                   tau_guess: float) -> tuple[float, float]:  # fmt: skip
    sol = tm.propagate(model, s, tau0, tau_guess + 0.6, dense=True).sol
    from scipy.optimize import minimize_scalar

    def d(t: float) -> float:
        return float(np.linalg.norm(tm.relative_inertial(model, body, t, sol.sol(t)[:4])[0]))

    ts = np.linspace(tau_guess - 0.5, tau_guess + 0.5, 4001)
    dd = [d(float(t)) for t in ts]
    i = int(np.argmin(dd))
    r = minimize_scalar(d, bounds=(ts[max(i - 1, 0)], ts[min(i + 1, 4000)]), method="bounded",
                        options={"xatol": 1e-12})  # fmt: skip
    return float(r.fun), float(r.x)


def stage_sensitivity() -> None:
    m, sh, z = converged_solution(1.0)
    s0, _mids, _ = sh.unpack(z)
    period = 5 * m.forcing_period
    sol = tm.propagate(m, s0, 0.0, period, dense=True).sol
    apos = apoapsis_times(m, sol)
    rows = []
    for body, tau_f in (("Oberon", sh.T), ("Titania", period)):
        ta = max(t for t in apos if t < tau_f - 0.5)
        s_a = sol.sol(ta)[:4]
        base_rp, base_t = next_periapsis(m, s_a, ta, body, tau_f)
        r_in, v_in = tm.inertial_from_rot(m, ta, s_a)
        rh = r_in / np.linalg.norm(r_in)
        th = np.array([-rh[1], rh[0]])
        log(
            f"{body}: apoapsis {m.days(tau_f - ta):.3f} d before the flyby; nominal rp {base_rp:.2f} km"
        )
        for kind, vec, mag in (("pos_radial_km", rh, 1.0), ("pos_along_km", th, 1.0),
                               ("vel_radial_ms", rh, 1e-3), ("vel_along_ms", th, 1e-3)):  # fmt: skip
            for scale in (1.0, 0.1):
                if kind.startswith("pos"):
                    sp = tm.rot_from_inertial(m, ta, r_in + scale * mag * vec, v_in)
                else:
                    sp = tm.rot_from_inertial(m, ta, r_in, v_in + scale * mag * vec)
                rp, tp = next_periapsis(m, sp, ta, body, tau_f)
                row = {"flyby": body, "perturbation": kind, "size": scale,
                       "apoapsis_lead_days": m.days(tau_f - ta),
                       "d_rp_km": rp - base_rp, "d_t_s": (tp - base_t) * m.time_unit_s,
                       "d_rp_km_per_unit": (rp - base_rp) / scale}  # fmt: skip
                rows.append(row)
                log(f"   {kind} x{scale}: d rp {rp - base_rp:+.3f} km, d t {row['d_t_s']:+.1f} s")
    # growth per half cycle and full cycle from the monodromy
    arc = tm.propagate(m, s0, 0.0, period, with_stm=True)
    assert arc.stm is not None
    write_json("sensitivity.json", {"rows": rows,
                                    "monodromy_singular_values": np.linalg.svd(arc.stm)[1]})  # fmt: skip


# ---------------------------------------------------------------------------
# variants from the stored #888 enumeration (step 6)
# ---------------------------------------------------------------------------


def stage_variants() -> None:
    p = ROOT / "data" / "found" / "888_turn_gate" / "extended_uranus.jsonl"
    rows = []
    for line in p.read_text().splitlines():
        g = json.loads(line)
        if g.get("kind") != "closure":
            continue
        if {g["anchor"], g["flyby"]} != {"Titania", "Oberon"}:
            continue
        rows.append({k: g.get(k) for k in ("anchor", "flyby", "tof_days", "n_rev", "branches",
                     "rel_offset_deg", "vinf_kms", "worst_ratio_project_floor",
                     "pass_project_floor", "n_commensurate_int")})  # fmt: skip
    rows.sort(key=lambda r: r["worst_ratio_project_floor"] or 1e9)
    for r in rows[:15]:
        log(f"   {r['anchor']}-{r['flyby']} n {r['n_commensurate_int']} tof {r['tof_days']:.2f} "
            f"n_rev {r['n_rev']} br {r['branches']} off {r['rel_offset_deg']} "
            f"worst ratio {r['worst_ratio_project_floor']:.3f}")  # fmt: skip
    write_json("variants.json", {"n_titania_oberon_closures": len(rows), "closures": rows})


# ---------------------------------------------------------------------------
# refine (POST HOC, added after verify): one full-cycle Newton step with the
# monodromy, then the V2/V3 closure tests again
# ---------------------------------------------------------------------------


def refined_start() -> tuple[tm.TwoMoonModel, np.ndarray, np.ndarray, list[dict[str, float]]]:
    m, sh, z = converged_solution(1.0)
    s0, _mids, _ = sh.unpack(z)
    period = 5 * m.forcing_period
    x = s0.copy()
    hist = []
    for _ in range(2):
        a = tm.propagate(m, x, 0.0, period, with_stm=True)
        assert a.stm is not None
        d = a.state - x
        hist.append(closure_error(m, x, a.state))
        x = x + np.linalg.solve(np.eye(4) - a.stm, d)
    return m, s0, x, hist


def stage_refine() -> None:
    m, s0, x, hist = refined_start()
    out: dict[str, Any] = {
        "post_hoc": True,
        "newton_closure_history": hist,
        "shift_from_shooting_solution": closure_error(m, s0, x),
        "refined_state": x,
    }
    log(f"refine: closure history {hist}; shift {out['shift_from_shooting_solution']}")
    ends = {}
    for method in ("DOP853", "Radau"):
        sol = full_cycle(m, x, 0.0, method=method).sol
        ce = closure_error(m, x, sol.y[:4, -1])
        ends[method] = sol
        out[method] = {"closure": ce, "pass": ce["pos_km"] <= 1 * KM and ce["vel_kms"] <= 1 * CM_S}
        log(f"refine {method}: closure {ce}")
    t_half = 2.5 * m.forcing_period
    out["DOP853_vs_Radau_at_T"] = closure_error(
        m, ends["DOP853"].sol(t_half)[:4], ends["Radau"].sol(t_half)[:4]
    )
    out["DOP853_vs_Radau_at_end"] = closure_error(
        m, ends["DOP853"].y[:4, -1], ends["Radau"].y[:4, -1]
    )
    log(
        f"refine: DOP853 vs Radau at T {out['DOP853_vs_Radau_at_T']}, end {out['DOP853_vs_Radau_at_end']}"
    )
    write_json("refine.json", out)


# ---------------------------------------------------------------------------
# realeph: first real-ephemeris look (no correction)
# ---------------------------------------------------------------------------

URA_KERNELS = (
    Path.home() / "GMAT" / "R2022a" / "data" / "time" / "SPICELeapSecondKernel.tls",
    Path.home() / "GMAT" / "R2022a" / "data" / "planetary_ephem" / "spk" / "uranian" / "ura111.bsp",
)
URANUS_J2 = 3.34343e-3  # Jacobson 2014 (as in data/validation/v4_uranus.py)
URANUS_R_EQ = 25559.0
MOONS5 = ("Miranda", "Ariel", "Umbriel", "Titania", "Oberon")


def stage_realeph(n_epochs: int) -> None:
    import spiceypy as spice

    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    for k in URA_KERNELS:
        spice.furnsh(str(k))

    def st(moon: str, et: float) -> np.ndarray:
        s, _ = spice.spkezr(moon.upper(), et, "J2000", "NONE", "URANUS")
        return np.asarray(s, dtype=np.float64)

    m, _s0, x, _ = refined_start()
    gm_moons = {mo: SATELLITES[mo].mu_km3_s2 for mo in MOONS5}
    gm_u = PRIMARIES["Uranus"] - sum(gm_moons.values())
    # conjunction epochs: Titania and Oberon at the same longitude in Titania's plane
    et0 = float(spice.str2et("2030-01-01T00:00:00"))

    def rel_lon(et: float) -> float:
        t, o = st("Titania", et), st("Oberon", et)
        h = np.cross(t[:3], t[3:])
        h /= np.linalg.norm(h)
        xh = t[:3] / np.linalg.norm(t[:3])
        yh = np.cross(h, xh)
        return math.atan2(float(o[:3] @ yh), float(o[:3] @ xh))

    epochs = []
    et = et0
    prev = rel_lon(et)
    while len(epochs) < n_epochs:
        et_n = et + 3600.0
        cur = rel_lon(et_n)
        if prev > 0.0 >= cur and abs(cur - prev) < 1.0:
            a, b = et, et_n
            for _ in range(50):
                mid = 0.5 * (a + b)
                if rel_lon(mid) > 0.0:
                    a = mid
                else:
                    b = mid
            epochs.append(0.5 * (a + b))
            et_n += 20 * 86400.0  # skip ahead, then keep scanning
            cur = rel_lon(et_n)
        et, prev = et_n, cur
    # the moon-relative flyby state at tau = 0
    dr, dv = tm.relative_inertial(m, "Titania", 0.0, x)
    period_s = 5 * m.forcing_period * m.time_unit_s
    rows = []
    for et_c in epochs:
        t = st("Titania", et_c)
        h = np.cross(t[:3], t[3:])
        zh = h / np.linalg.norm(h)
        xh = t[:3] / np.linalg.norm(t[:3])
        yh = np.cross(zh, xh)
        rot = np.column_stack([xh, yh, zh])
        r0 = t[:3] + rot @ np.array([dr[0], dr[1], 0.0])
        v0 = t[3:] + rot @ np.array([dv[0], dv[1], 0.0])

        def rhs(ts: float, y: np.ndarray, et_c: float = et_c, zh: np.ndarray = zh) -> np.ndarray:
            r = y[:3]
            rn = float(np.linalg.norm(r))
            a = -gm_u * r / rn**3
            zc = float(r @ zh)
            c = -1.5 * gm_u * URANUS_J2 * URANUS_R_EQ**2 / rn**5
            a = a + c * ((1 - 5 * zc**2 / rn**2) * r + (2 * zc) * zh)
            # J2 about the axis zh: a = c [ (1 - 5 z^2/r^2) r + 2 z zh ]
            for mo in MOONS5:
                rm = st(mo, et_c + ts)[:3]
                d = r - rm
                a = a - gm_moons[mo] * (
                    d / float(np.linalg.norm(d)) ** 3 + rm / float(np.linalg.norm(rm)) ** 3
                )
            return np.concatenate([y[3:], a])

        t1 = time.time()
        from scipy.integrate import solve_ivp

        sol = solve_ivp(rhs, (0.0, period_s), np.concatenate([r0, v0]), method="DOP853",
                        rtol=1e-11, atol=1e-6, dense_output=True)  # fmt: skip
        ts = np.linspace(0.0, period_s, 40001)
        ys = sol.sol(ts)
        res: dict[str, Any] = {"epoch_utc": spice.et2utc(et_c, "ISOC", 0), "wall_s": None}
        for mo in ("Titania", "Oberon"):
            dd = np.array(
                [np.linalg.norm(ys[:3, i] - st(mo, et_c + tt)[:3]) for i, tt in enumerate(ts)]
            )
            hill = tm.hill_radius_km("Uranus", mo)
            idx = np.nonzero((dd[1:-1] < dd[:-2]) & (dd[1:-1] < dd[2:]))[0] + 1
            mins = [
                {"t_days": float(ts[i] / 86400), "dist_km": float(dd[i])}
                for i in idx
                if dd[i] < 3 * hill
            ]
            res[mo] = {"start_dist_km": float(dd[0]), "minima_within_3RH": mins,
                       "closest_overall_km": float(dd[1:].min())}  # fmt: skip
        tit_end = st("Titania", et_c + period_s)
        res["end_dist_to_titania_km"] = float(np.linalg.norm(ys[:3, -1] - tit_end[:3]))
        res["wall_s"] = time.time() - t1
        log(f"epoch {res['epoch_utc']}: Titania minima {res['Titania']['minima_within_3RH']}; "
            f"Oberon minima {res['Oberon']['minima_within_3RH']}; end-to-Titania "
            f"{res['end_dist_to_titania_km']:.0f} km [{res['wall_s']:.0f} s]")  # fmt: skip
        rows.append(res)
    write_json("realeph.json", {"model": "Uranus point mass (GM_sys minus the five moons) + J2 about "
                                "Titania's orbit normal + Miranda/Ariel/Umbriel/Titania/Oberon point "
                                "masses from URA111 (indirect terms included), no correction",
                                "epochs": rows})  # fmt: skip


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["rebuild", "decisive", "direct", "continue", "verify", "sensitivity",
                             "variants", "refine", "realeph"])  # fmt: skip
    ap.add_argument("--max-steps", type=int, default=100)
    ap.add_argument("--budget-s", type=float, default=420.0)
    ap.add_argument("--n-epochs", type=int, default=3)
    args = ap.parse_args(argv)
    log(f"#890 stage {args.stage} (git {git_sha()})")
    if args.stage == "rebuild":
        stage_rebuild()
    elif args.stage == "decisive":
        stage_decisive()
    elif args.stage == "direct":
        stage_direct()
    elif args.stage == "continue":
        stage_continue(args.max_steps, args.budget_s)
    elif args.stage == "verify":
        stage_verify()
    elif args.stage == "sensitivity":
        stage_sensitivity()
    elif args.stage == "variants":
        stage_variants()
    elif args.stage == "refine":
        stage_refine()
    elif args.stage == "realeph":
        stage_realeph(args.n_epochs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
