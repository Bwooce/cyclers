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
from cyclerfinder.verify.turn_gate import demanded_turn_gate
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
# main
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["rebuild", "decisive", "direct", "continue", "verify", "sensitivity",
                             "variants"])  # fmt: skip
    ap.add_argument("--max-steps", type=int, default=100)
    ap.add_argument("--budget-s", type=float, default=420.0)
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
    else:
        raise SystemExit(f"stage {args.stage} not implemented yet")
    return 0


if __name__ == "__main__":
    sys.exit(main())
