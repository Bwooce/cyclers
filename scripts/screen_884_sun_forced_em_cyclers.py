"""#884 -- which catalogued Earth-Moon cycler families survive the Sun (BCR4BP).

Staged, resumable driver (each stage writes one JSON per task under
``data/found/884_sun_forced_em_cyclers/<stage>/``; re-running skips finished tasks):

  --stage controls   positive control against Brown et al. 2025: planar L1 Lyapunov at
                     T* = Tg/2 (a = 2) and planar L2 Lyapunov at T* = Tg (a = 1).
  --stage families   walk every catalogued Earth-Moon family (symmetric seeds) in the CR3BP
                     and correct its members at the low-order commensurate periods.
  --stage melnikov   Melnikov scan per commensurate member (zeros, amplitude, checks).
  --stage continue   continue each Melnikov zero in eps = mu_sun / mu_sun_phys to 1.
  --stage summary    assemble summary.json.

Model: :mod:`cyclerfinder.core.bcr4bp` constants (Andreu / Rosales-Jorba), homotopy on the
whole solar term, Sun synodic rate fixed, forcing period Tg = 2 pi / omega_S. See
``docs/notes/2026-10-04-884-sun-forced-em-cyclers.md``.

Run: ``uv run python scripts/screen_884_sun_forced_em_cyclers.py --stage <name> [--workers 4]``
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.search import sun_forced_periodic_884 as sf  # noqa: E402

OUT = ROOT / "data" / "found" / "884_sun_forced_em_cyclers"
EM_MU_ROWS = (0.0121505, 0.0121506)  # catalogue rows at the Earth-Moon mass ratio


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {msg}", flush=True)


def ratio_of(period: float, tg: float, max_den: int = 6) -> Fraction:
    return Fraction(period / tg).limit_denominator(max_den)


# ---------------------------------------------------------------------------
# Seeds.
# ---------------------------------------------------------------------------


def _perpendicular_crossing(mu: float, state: list[float], period: float) -> list[float] | None:
    """A symmetric start (x, 0, z, 0, vy, 0) on the orbit, or None if none is found."""
    from scipy.integrate import solve_ivp

    from cyclerfinder.core import cr3bp

    x = np.asarray(state, dtype=float)
    if abs(x[1]) < 1e-9 and abs(x[3]) < 1e-6 and abs(x[5]) < 1e-6:
        return [float(x[0]), float(x[2]), float(x[4])]

    def ev(t: float, y: Any) -> float:
        return float(y[1])

    sol = solve_ivp(
        lambda t, y: cr3bp.cr3bp_eom(t, y, mu),
        (0, period),
        x,
        events=ev,
        rtol=1e-12,
        atol=1e-12,
        method="DOP853",
    )
    for y in sol.y_events[0]:
        if abs(y[3]) < 1e-6 and abs(y[5]) < 1e-6:
            return [float(y[0]), float(y[2]), float(y[4])]
    return None


def catalogue_seeds() -> list[dict[str, Any]]:
    rows = yaml.safe_load((ROOT / "data" / "catalogue.yaml").read_text())
    seeds: list[dict[str, Any]] = []
    for r in rows:
        c = ((r.get("orbit_elements") or {}).get("cr3bp")) or {}
        if c.get("period_nd") is None or not c.get("state_nd"):
            continue
        mu = float(c["mass_ratio"])
        if not (EM_MU_ROWS[0] < mu < EM_MU_ROWS[1]):
            continue
        sym = _perpendicular_crossing(mu, c["state_nd"], float(c["period_nd"]))
        seeds.append(
            {
                "row": r["id"],
                "period": float(c["period_nd"]),
                "jacobi": c.get("jacobi_constant"),
                "seed": sym,
                "orbit_class": r.get("orbit_class"),
            }
        )
    return seeds


# One seed per catalogued family (the walk covers the rest of the family).
FAMILY_SEEDS: dict[str, str] = {
    "C11 (Braik-Ross C11a)": "braik-ross-c11a-cycler-2026",
    "C11 (Braik-Ross C11b)": "braik-ross-c11b-cycler-2026",
    "C11 (Ross-RT 1:1)": "ross-rt-em-cycler-11-2025",
    "C21 planar (Ross-RT 2:1)": "ross-rt-em-cycler-21-2025",
    "C21 spatial (#438)": "em-cycler-21-3d-spatial-2026",
    "C21 3D corridor (#682)": "braik-ross-c21-3d-corridor-05-2026",
    "C31 (Ross-RT 3:1)": "ross-rt-em-cycler-31-2025",
    "C32 (Ross-RT 3:2)": "ross-rt-em-cycler-32-2025",
    "C32 (Braik-Ross)": "braik-ross-c32-cycler-2026",
    "C33 (Ross-RT 3:3)": "ross-rt-em-cycler-33-2025",
    "C32 3D corridor (#682)": "braik-ross-c32-3d-corridor-05-2026",
    "L1 Lyapunov 3D corridor (#682)": "lyapunov3d-l1-corridor-04-2026",
    "R21-S planar resonant (Braik-Ross)": "braik-ross-planar-r21-s-corridor-2026",
    "R31-S planar resonant (Braik-Ross)": "braik-ross-planar-r31-s-corridor-2026",
    "R52-S planar resonant (Braik-Ross)": "braik-ross-planar-r52-s-corridor-2026",
    "Casoliva 1:2 (c)": "casoliva-1-2c-em-resonant-po-2010",
    "Casoliva 1:2 (d)": "casoliva-1-2d-em-resonant-po-2010",
    "Casoliva 2:1 (a)": "casoliva-2-1a-em-resonant-po-2010",
    "Casoliva 2:1 (b)": "casoliva-2-1b-em-resonant-po-2010",
    "Casoliva 3:2 (c)": "casoliva-3-2c-em-resonant-po-2010",
    "Casoliva 7:3 (a)": "casoliva-7-3a-em-cycler-2010",
    "Casoliva 7:3 (b)": "casoliva-7-3b-em-cycler-2010",
    "Casoliva 7:3 (c)": "casoliva-7-3c-em-cycler-2010",
    "Vaquero 2:1": "vaquero-21-c246-em-cycler-2013",
    "Vaquero 3:1": "vaquero-31-c254-em-cycler-2013",
}


def targets(tg: float) -> list[float]:
    """Half-integer multiples of Tg (a = 1, 2) plus the coordinator's a = 3 and a = 6 cases."""
    out = [k * tg / 2 for k in range(1, 13)]
    out += [8 * tg / 3, 5 * tg / 6]
    return out


# ---------------------------------------------------------------------------
# Stage workers.
# ---------------------------------------------------------------------------


def _write(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=float))
    tmp.replace(path)


def _cr3bp_info(model: sf.ForcedModel, member: sf.SymmetricOrbit) -> dict[str, Any]:
    arc = sf.propagate(model, 0.0, 0.0, member.state, 0.0, member.period, variational=True)
    assert arc.stm is not None
    eig = np.linalg.eigvals(arc.stm)
    prob = sf.ShootingProblem(model, member.period, 0.0, 1)
    sol = sf.propagate(model, 0.0, 0.0, member.state, 0.0, member.period, dense=True)
    ts = np.linspace(0, member.period, int(400 * member.period))
    st = sol.sol(ts)
    mu = model.mu
    dm = np.sqrt((st[0] - 1 + mu) ** 2 + st[1] ** 2 + st[2] ** 2).min() * sf.EM_LENGTH_KM
    del prob
    return {
        "max_abs_floquet_cr3bp": float(np.max(np.abs(eig))),
        "periselene_km_cr3bp": float(dm),
        "jacobi": sf.jacobi(member.state, mu),
    }


def family_task(name: str, seed: dict[str, Any]) -> dict[str, Any]:
    model = sf.default_model()
    tg = model.tg
    tl = targets(tg)
    s = seed["seed"]
    rec: dict[str, Any] = {"family": name, "row": seed["row"], "seed_period": seed["period"]}
    if s is None:
        rec["status"] = "no perpendicular crossing (asymmetric orbit); skipped"
        return rec
    members: list[dict[str, Any]] = []
    walks = []
    for direction in (+1.0, -1.0):
        found, info = sf.walk_family(
            model,
            np.array(s),
            seed["period"],
            tl,
            direction,
            max_steps=250,
            t_window=(max(0.5, seed["period"] - 10), seed["period"] + 10),
        )
        walks.append({"direction": direction, **{k: v for k, v in info.items()}})
        for m in found:
            fr = ratio_of(m.period, tg)
            entry = {
                "period": m.period,
                "ratio_T_over_Tg": f"{fr.numerator}/{fr.denominator}",
                "a": fr.denominator,
                "n": fr.numerator,
                "state": m.state.tolist(),
                "residual": m.residual,
                **_cr3bp_info(model, m),
            }
            if not any(
                np.max(np.abs(np.array(e["state"]) - m.state)) < 1e-6
                and abs(e["period"] - m.period) < 1e-9
                for e in members
            ):
                members.append(entry)
    rec.update({"status": "ok", "walks": walks, "members": members})
    return rec


def melnikov_task(member: dict[str, Any]) -> dict[str, Any]:
    model = sf.default_model()
    a, n = member["a"], member["n"]
    p = n * model.tg
    x0 = np.array(member["state"])
    th, mel, zeros = sf.find_zeros(model, x0, p, a)
    span = 2 * math.pi / a
    probe = np.array([0.3 * span, 0.7 * span])
    shifted = sf.melnikov_scan(model, x0, p, probe + span)
    base = sf.melnikov_scan(model, x0, p, probe)
    var = sf.melnikov_variational(model, x0, p, float(probe[0]))
    amp = float(np.max(np.abs(mel)))
    # Derivative sign at each zero (alternation <-> elliptic / hyperbolic pairing).
    slopes = []
    for z in zeros:
        h = 1e-4
        f = sf.melnikov_scan(model, x0, p, np.array([z - h, z + h]))
        slopes.append(float((f[1] - f[0]) / (2 * h)))
    return {
        **member,
        "forced_period": p,
        "melnikov_amplitude": amp,
        "melnikov_zeros": zeros,
        "melnikov_slopes": slopes,
        "periodicity_check": float(np.max(np.abs(shifted - base))),
        "variational_vs_quadrature": float(abs(var - base[0])),
        "theta_grid": th.tolist()[::4],
        "mel_grid": mel.tolist()[::4],
    }


def continue_task(member: dict[str, Any], theta0: float, max_steps: int) -> dict[str, Any]:
    model = sf.default_model()
    p = member["forced_period"]
    nseg = max(2, math.ceil(p / 1.5))
    prob = sf.ShootingProblem(model, p, theta0, nseg)
    xs = prob.nodes_from_orbit(np.array(member["state"]))
    t0 = time.time()
    br = sf.continue_in_eps(prob, xs, max_steps=max_steps)
    rec: dict[str, Any] = {
        "family": member["family"],
        "row": member["row"],
        "ratio": member["ratio_T_over_Tg"],
        "a": member["a"],
        "period_cr3bp": member["period"],
        "forced_period": p,
        "theta0": theta0,
        "n_seg": nseg,
        "stop_reason": br.stop_reason,
        "eps_path": br.eps,
        "max_abs_floquet_path": br.max_abs_eig,
        "folds": br.folds,
        "eps_max": max(br.eps) if br.eps else None,
        "wall_s": time.time() - t0,
    }
    if br.stop_reason == "reached_target":
        xs1 = br.nodes[-1]
        rec["diagnostics_eps1"] = sf.orbit_diagnostics(prob, xs1, 1.0)
        rec["nodes_eps1"] = xs1.tolist()
        rec["shift_from_cr3bp_max"] = float(np.max(np.abs(xs1 - xs)))
    elif br.nodes:
        rec["nodes_last"] = br.nodes[-1].tolist()
    return rec


def control_task(name: str, seed: list[float], seed_period: float, target: float, a: int) -> Any:
    model = sf.default_model()
    member, info = sf.continue_to_period(model, np.array(seed), seed_period, target)
    if member is None:
        return {"control": name, "status": "member not reached", "info": info}
    p = target * a
    _, mel, zeros = sf.find_zeros(model, member.state, p, a)
    out: dict[str, Any] = {
        "control": name,
        "a": a,
        "period_cr3bp": member.period,
        "forced_period": p,
        "state": member.state.tolist(),
        "melnikov_zeros": zeros,
        "melnikov_amplitude": float(np.max(np.abs(mel))),
        "branches": [],
    }
    nseg = max(2, math.ceil(p / 1.5))
    trial = [*zeros, 0.37 * 2 * math.pi / a, 0.81 * 2 * math.pi / a]
    for th0 in trial:
        prob = sf.ShootingProblem(model, p, th0, nseg)
        xs = prob.nodes_from_orbit(member.state)
        br = sf.continue_in_eps(prob, xs, max_steps=200)
        b: dict[str, Any] = {
            "theta0": th0,
            "is_melnikov_zero": th0 in zeros,
            "stop_reason": br.stop_reason,
            "eps_max": max(br.eps) if br.eps else None,
            "folds": br.folds,
        }
        if br.stop_reason == "reached_target":
            b["diagnostics_eps1"] = sf.orbit_diagnostics(prob, br.nodes[-1], 1.0)
            # Sun phase at the orbit's x-axis crossing nearest t = 0 (symmetry test).
            b["start_state_eps1"] = br.nodes[-1][0].tolist()
        out["branches"].append(b)
    return out


# ---------------------------------------------------------------------------
# Stage drivers.
# ---------------------------------------------------------------------------


def run_pool(tasks: list[tuple[Path, Any, tuple[Any, ...]]], workers: int, budget_s: float) -> None:
    todo = [t for t in tasks if not t[0].exists()]
    log(f"{len(tasks)} tasks, {len(todo)} to do, workers={workers}")
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {}
        it = iter(todo)
        for _ in range(workers):
            nxt = next(it, None)
            if nxt is not None:
                futs[ex.submit(nxt[1], *nxt[2])] = nxt
        while futs:
            for fut in as_completed(list(futs)):
                path, _, args = futs.pop(fut)
                try:
                    res = fut.result()
                    _write(path, res)
                    done += 1
                    log(f"done {path.name} ({done}/{len(todo)}, {time.time() - t0:.0f} s)")
                except Exception as exc:
                    _write(path, {"error": f"{type(exc).__name__}: {exc}", "args": str(args)[:300]})
                    log(f"ERROR {path.name}: {exc}")
                if time.time() - t0 < budget_s:
                    nxt = next(it, None)
                    if nxt is not None:
                        futs[ex.submit(nxt[1], *nxt[2])] = nxt
                break
    remaining = [t for t in tasks if not t[0].exists()]
    log(f"stage chunk finished; {len(remaining)} tasks remain")


def stage_controls(workers: int, budget: float) -> None:
    model = sf.default_model()
    tg = model.tg
    # Seeds: Braik-Ross Table 2 LL1 / LL2 representatives (data/golden/...family_ics.yaml).
    tasks = [
        (
            OUT / "controls" / "L1_lyapunov_a2.json",
            control_task,
            (
                "L1 planar Lyapunov, T* = Tg/2",
                [0.8115256290557147, 0.0, 0.2561843220006502],
                2.946253150022597,
                tg / 2,
                2,
            ),
        ),
        (
            OUT / "controls" / "L2_lyapunov_a1.json",
            control_task,
            (
                "L2 planar Lyapunov, T* = Tg",
                [1.100554841329441, 0.0, 0.2702199577829618],
                3.476412471522109,
                tg,
                1,
            ),
        ),
    ]
    run_pool(tasks, workers, budget)


def stage_families(workers: int, budget: float) -> None:
    seeds = {s["row"]: s for s in catalogue_seeds()}
    tasks = []
    for name, row in FAMILY_SEEDS.items():
        tasks.append((OUT / "families" / f"{row}.json", family_task, (name, seeds[row])))
    run_pool(tasks, workers, budget)


def _members() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for f in sorted((OUT / "families").glob("*.json")):
        rec = json.loads(f.read_text())
        for m in rec.get("members", []):
            out.append({"family": rec["family"], "row": rec["row"], **m})
    # Deduplicate across seeds of the same family (same period, same state).
    uniq: list[dict[str, Any]] = []
    for m in out:
        if not any(
            abs(u["period"] - m["period"]) < 1e-9
            and np.max(np.abs(np.array(u["state"]) - np.array(m["state"]))) < 1e-6
            for u in uniq
        ):
            uniq.append(m)
    return uniq


def _mkey(m: dict[str, Any]) -> str:
    return f"{m['row']}__{m['n']}-{m['a']}__x{m['state'][0]:+.5f}"


def stage_melnikov(workers: int, budget: float) -> None:
    tasks = [(OUT / "melnikov" / f"{_mkey(m)}.json", melnikov_task, (m,)) for m in _members()]
    run_pool(tasks, workers, budget)


def stage_continue(workers: int, budget: float, max_a: int, max_steps: int) -> None:
    tasks = []
    for f in sorted((OUT / "melnikov").glob("*.json")):
        m = json.loads(f.read_text())
        if "error" in m or m["a"] > max_a:
            continue
        for i, z in enumerate(m["melnikov_zeros"]):
            path = OUT / "continue" / f"{f.stem}__z{i}.json"
            tasks.append((path, continue_task, (m, z, max_steps)))
    run_pool(tasks, workers, budget)


def stage_summary() -> None:
    rows = []
    for f in sorted((OUT / "continue").glob("*.json")):
        r = json.loads(f.read_text())
        if "error" in r:
            rows.append({"task": f.stem, "error": r["error"]})
            continue
        d = r.get("diagnostics_eps1") or {}
        rows.append(
            {
                "task": f.stem,
                "family": r["family"],
                "ratio": r["ratio"],
                "a": r["a"],
                "period_cr3bp": r["period_cr3bp"],
                "theta0": r["theta0"],
                "stop_reason": r["stop_reason"],
                "eps_max": r["eps_max"],
                "folds": r["folds"],
                "closure_dop853": d.get("closure_dop853"),
                "closure_radau": d.get("closure_radau"),
                "max_abs_floquet": d.get("max_abs_floquet"),
                "stability": d.get("stability"),
                "periselene_km": d.get("periselene_km"),
                "cycler_class": d.get("cycler_class"),
            }
        )
    fams = []
    for f in sorted((OUT / "families").glob("*.json")):
        r = json.loads(f.read_text())
        fams.append(
            {
                "family": r.get("family"),
                "row": r.get("row"),
                "status": r.get("status"),
                "period_range_walked": [
                    min((w.get("T_min", math.inf) for w in r.get("walks", [])), default=None),
                    max((w.get("T_max", -math.inf) for w in r.get("walks", [])), default=None),
                ],
                "walk_stops": [w.get("reason") for w in r.get("walks", [])],
                "members": [
                    {
                        k: m[k]
                        for k in (
                            "ratio_T_over_Tg",
                            "period",
                            "jacobi",
                            "periselene_km_cr3bp",
                            "max_abs_floquet_cr3bp",
                        )
                    }
                    for m in r.get("members", [])
                ],
            }
        )
    mels = []
    for f in sorted((OUT / "melnikov").glob("*.json")):
        m = json.loads(f.read_text())
        if "error" in m:
            continue
        mels.append(
            {
                k: m[k]
                for k in (
                    "family",
                    "ratio_T_over_Tg",
                    "a",
                    "period",
                    "melnikov_amplitude",
                    "melnikov_zeros",
                    "melnikov_slopes",
                    "periodicity_check",
                    "variational_vs_quadrature",
                )
            }
        )
    _write(OUT / "summary.json", {"families": fams, "melnikov": mels, "continuations": rows})
    log(f"summary: {len(fams)} families, {len(mels)} Melnikov scans, {len(rows)} branches")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage",
        required=True,
        choices=["controls", "families", "melnikov", "continue", "summary"],
    )
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument(
        "--budget-s",
        type=float,
        default=420.0,
        help="stop submitting new tasks after this many seconds",
    )
    ap.add_argument("--max-a", type=int, default=3)
    ap.add_argument("--max-steps", type=int, default=150)
    args = ap.parse_args()
    log(f"#884 stage {args.stage}; Tg = {sf.default_model().tg:.9f} TU")
    w = min(args.workers, 4)
    if args.stage == "controls":
        stage_controls(w, args.budget_s)
    elif args.stage == "families":
        stage_families(w, args.budget_s)
    elif args.stage == "melnikov":
        stage_melnikov(w, args.budget_s)
    elif args.stage == "continue":
        stage_continue(w, args.budget_s, args.max_a, args.max_steps)
    else:
        stage_summary()


if __name__ == "__main__":
    main()
