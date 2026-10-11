"""#1055: degenerate-family seed correction (pre-registration:
docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md).

    uv run python scripts/run_1055_seedfix.py --target gc-1|gc4-19 --stage fix|start|neg|cont
    uv run python scripts/run_1055_seedfix.py --target ganeur316|gancal1|gc1e2 \
        --stage fix|start|neg|cont|down|ias15

- fix: detection (sigma 0.02 vs 0.005) and the family root at sigma0 = 0.02; writes fix.json.
- start: plain damped Newton on the full system at sigma0 from z* (the corrected seed).
- neg: the negative control: Newton from the family moved off the root (alpha* +- 0.5 L).
- cont: the object's own continuation from the start (ideal: the #1025 driver; jup365: the #1044
  open chain, stages sigma, down, ias15).
Data in data/1055_<target>/; live log in .../live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.degenerate_family import (
    bordered_solve,
    column_scales,
    detect_degenerate,
    family_frame,
    family_root_1d,
    family_root_nd,
    scaled_svd,
)

ROOT = Path(__file__).resolve().parents[1]
SIGMA0 = 0.02
TOL = 3e-3  # complement rows, weighted (amendment 5: above the sigma-0.02 noise floor)
G_FLOOR = 1e-3  # amendment 5: a root is |g| below the start's noise-floor rule (weighted)
IDEAL = ("gc-1", "gc4-19")
FLOOR_KM = {"Ganymede": 100.0, "Callisto": 200.0, "Europa": 100.0}  # project floors
REAL = {"ganeur316": "968_rungb", "gancal1": "1044_gancal1", "gc1e2": "1044_gc1"}


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class Target:
    """A uniform face over the two formulations: rj(z) at a sigma, scale, describe, newton."""

    def __init__(self, name: str, out: Path) -> None:
        self.name, self.out = name, out
        self.ckpt_tag = "start"
        self.max_wall = 430.0
        out.mkdir(parents=True, exist_ok=True)
        if name in IDEAL:
            self.b = _load("run_1025_sigma_batch", ROOT / "scripts" / "run_1025_sigma_batch.py")
            self.b.OUT = out
            self.p = self.b.build(name)
            self.seed = self.p.seed
            self.moons = list(self.p.moons)
            self.chain_vinf = list(self.p.pc_vinf)
        else:
            rb = _load("run_968_rungb", ROOT / "scripts" / "run_968_rungb.py")
            rb.OUT = out
            rb.END_MODE = "direction_seedmag"
            if name != "ganeur316":  # #1044 settings; GanEur#316 keeps the #1039 (b') defaults
                rb.PER_CYCLE = 4
                rb.FORCE = ("Ganymede", "Callisto")
                rb.SCALED = ("Ganymede", "Callisto")
                rb.SEED_CAP_KM = 1.0e12
            seed = out / "seed_chain.json"
            if not seed.exists():
                shutil.copy(ROOT / "data" / REAL[name] / "seed_chain.json", seed)
            self.rb = rb
            self.c = rb.build(1)
            self.seed = self.c.seed
            self.moons = list(self.c.moons)
            self.chain_vinf = list(self.c.chain_vinf)

    def log(self, msg: str) -> None:
        (self.out / "live").mkdir(parents=True, exist_ok=True)
        line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {self.name}: {msg}"
        print(line, flush=True)
        with (self.out / "live" / "runlog.txt").open("a") as fh:
            fh.write(line + "\n")

    def set_sigma(self, s: float) -> None:
        (self.p if hasattr(self, "p") else self.c).set_sigma(s)

    def scale(self, z: Any, f: float) -> Any:
        if hasattr(self, "p"):
            return self.b.scale_offsets(self.p, z, f)
        return self.rb.scale_offsets(self.c, z, f)

    def rj(self, z: Any) -> tuple[Any, Any]:
        if hasattr(self, "p"):
            r, j, _ = self.b.residual_and_jac(self.p, z)
        else:
            r, j, _ = self.rb.residual_and_jac(self.c, z)
        if not (np.all(np.isfinite(r)) and np.all(np.isfinite(j))):
            raise RuntimeError("non-finite residual")
        return r, j

    def newton(self, z: Any, iters: int) -> tuple[Any, dict[str, Any], bool]:
        if hasattr(self, "p"):
            zz, info = self.b.solve(self.p, z, iters, 0, False)
            return zz, info, bool(self.b.converged(info))
        # amendment 3: the start may be noise-floor-limited, as rungb's own start (#1039, #1044).
        # One rungb Newton iteration at a time, so that the rule is checked after each, and the
        # iterate is checkpointed (a call stopped by the time limit resumes from it).
        ck = self.out / f"newton_{self.ckpt_tag}.npy"
        zz = np.load(ck) if ck.exists() else np.asarray(z, dtype=np.float64)
        info: dict[str, Any] = {"dr": [np.inf], "dv": [np.inf]}
        ok = False
        for _ in range(iters):
            zz, info = self.rb.newton(self.c, zz, 1, f"{self.name} newton")
            np.save(ck, zz)
            ends = max(info.get("end_asymptote_kms", [0.0]))
            noise = max(info["dr"]) < 1e-2 and max(info["dv"]) < 1e-6
            ok = (bool(self.rb.converged(info)) or noise) and ends < 1e-6
            if ok:
                break
        return zz, info, ok

    def vinf(self, z: Any, s: float) -> list[float]:
        if hasattr(self, "p"):
            return [n["vinf_kms"] for n in self.b.describe(self.p, z, s)]
        return [n["vinf_kms"] for n in self.rb.describe(self.c, z, s)["nodes"]]


def stage_fix(t: Target) -> None:
    """Detection, then the family root. Resumable (amendment 6): the detection and the 1-d family
    state live in fix_state.json; a call stops after --max-wall seconds and the next continues."""
    sp = t.out / "fix_state.json"
    fs: dict[str, Any] = json.loads(sp.read_text()) if sp.exists() else {}
    t.set_sigma(SIGMA0)
    z0 = t.scale(t.seed, SIGMA0)
    if "q" not in fs:
        t.set_sigma(SIGMA0 / 4)
        _, j_lo = t.rj(t.scale(t.seed, SIGMA0 / 4))
        t.set_sigma(SIGMA0)
        _, j_hi = t.rj(z0)
        q = detect_degenerate(j_hi, j_lo, 4.0)
        s_hi = scaled_svd(j_hi, column_scales(j_hi))[1]
        s_lo = scaled_svd(j_lo, column_scales(j_lo))[1]
        t.log(f"detect: q {q}; sv(0.02) {s_hi[-4:].tolist()} sv(0.005) {s_lo[-4:].tolist()}")
        fs.update(q=q, sv_hi=s_hi[-4:].tolist(), sv_lo=s_lo[-4:].tolist(), family={})
        sp.write_text(json.dumps(fs))
    q = int(fs["q"])
    rec: dict[str, Any] = {"q": q, "sv_hi": fs["sv_hi"], "sv_lo": fs["sv_lo"]}
    t0 = time.monotonic()
    if q == 0:
        rec.update(status="no degenerate direction (identity)", z_star=z0.tolist(), alpha_star=0.0)
    elif q == 1:

        def save(state: dict[str, Any]) -> None:
            fs["family"] = state
            sp.write_text(json.dumps(fs))

        def domain_ok(z: Any) -> bool:
            """Amendment 7: every flyby node at or above the sigma-scaled surface plus floor."""
            for k, moon in enumerate(t.moons):
                floor = SATELLITES[moon].radius_eq_km + FLOOR_KM[moon]
                if float(np.linalg.norm(z[7 * k : 7 * k + 3])) < SIGMA0 * floor:
                    return False
            return True

        res = family_root_1d(
            t.rj,
            z0,
            TOL,
            log=t.log,
            g_floor=G_FLOOR,
            state=fs["family"],
            save=save,
            deadline=t0 + t.max_wall,
            domain_ok=domain_ok,
        )
        if res.status == "paused":
            t.log(f"fix: paused in phase {fs['family'].get('phase')}; run the stage again")
            return
        rec.update(status=res.status, alpha_star=res.alpha_star, newton_length=res.newton_length)
        rec["scan"] = res.scan
        rec["z_star"] = None if res.z_star is None else res.z_star.tolist()
    else:
        res = family_root_nd(t.rj, z0, q, TOL, log=t.log, g_floor=G_FLOOR)
        rec.update(status=res.status, alpha_star=res.alpha_star, newton_length=res.newton_length)
        rec["scan"] = res.scan
        rec["z_star"] = None if res.z_star is None else res.z_star.tolist()
    d = column_scales(t.rj(z0)[1])
    rec["wall_s"] = time.monotonic() - t0
    if rec["z_star"] is not None:
        zs = np.asarray(rec["z_star"])
        rec["vinf_star"] = t.vinf(zs, SIGMA0)
        rec["vinf_seed"] = t.vinf(z0, SIGMA0)
        rec["move_scaled"] = float(np.linalg.norm(d * (zs - z0)))
    (t.out / "fix.json").write_text(json.dumps(rec))
    t.log(f"fix: {rec['status']} alpha* {rec.get('alpha_star')} wall {rec['wall_s']:.0f} s")


def stage_start(t: Target) -> None:
    if not (t.out / "fix.json").exists():
        t.log("start: no fix.json (the fix stage has not finished)")
        raise SystemExit(1)
    rec = json.loads((t.out / "fix.json").read_text())
    if rec["z_star"] is None:
        t.log("start: no corrected seed")
        return
    t.set_sigma(SIGMA0)
    z, info, ok = t.newton(np.asarray(rec["z_star"]), 40)
    out = {"converged": ok, "dr": max(info["dr"]), "dv": max(info["dv"])}
    if ok:
        out["z"] = z.tolist()
        out["vinf"] = t.vinf(z, SIGMA0)
    (t.out / "start.json").write_text(json.dumps(out))
    t.log(f"start: converged={ok} dr {out['dr']:.2e} dv {out['dv']:.2e} vinf {out.get('vinf')}")


def stage_neg(t: Target) -> None:
    """Negative control: the family moved off the root by +- 0.5 L along each of the q null
    directions of the base frame (arclength rows about the corrected seed); Newton (40) must FAIL
    to start from every such point."""
    if not (t.out / "fix.json").exists():
        t.log("neg: no fix.json (the fix stage has not finished)")
        raise SystemExit(1)
    rec = json.loads((t.out / "fix.json").read_text())
    q = int(rec["q"])
    if q == 0 or rec["z_star"] is None:
        t.log("neg: needs q >= 1 and a root")
        return
    t.set_sigma(SIGMA0)
    z0 = t.scale(t.seed, SIGMA0)
    fr = family_frame(t.rj, z0, q, TOL)  # the same base and vectors as the fix stage
    d, u_n, v_n = fr.d, fr.u_n, fr.v_n
    zs = np.asarray(rec["z_star"])
    big_l = float(rec["newton_length"])
    out = {}
    for i in range(q):
        for sgn, lab in ((1.0, "plus"), (-1.0, "minus")):
            key = f"v{i}_{lab}"
            a = np.zeros(q)
            a[i] = sgn * 0.5 * big_l
            br = bordered_solve(t.rj, zs + (v_n @ a) / d, d, u_n, v_n, zs, a, TOL)
            row: dict[str, Any] = {"offset": a.tolist(), "family_point": bool(br.converged)}
            if br.converged:
                t.ckpt_tag = f"neg_{key}"
                _, info, ok = t.newton(br.z, 40)
                row.update(started=ok, dr=max(info["dr"]), dv=max(info["dv"]))
            out[key] = row
            t.log(f"neg {key}: {row}")
    (t.out / "neg.json").write_text(json.dumps(out))


def stage_neg2(t: Target) -> None:
    """Amendment 8 (q = 1): the negative control in the LOCAL null frame at the root (as amendment
    7 traces the family), oriented along the scan's direction of travel, on BOTH sides: beyond the
    root (0.5 L, then 0.25 L if not constructible) and back towards the base (-0.5 L, then -0.25 L).
    The first constructed point on each side is started with Newton (40) and must FAIL."""
    rec = json.loads((t.out / "fix.json").read_text())
    assert int(rec["q"]) == 1 and rec["z_star"] is not None
    t.set_sigma(SIGMA0)
    zs = np.asarray(rec["z_star"])
    big_l = float(rec["newton_length"])
    _, js = t.rj(zs)
    d = column_scales(js)
    u, _, vt = scaled_svd(js, d)
    u_n, v_n = u[:, -1:], vt[-1:].T
    fz = json.loads((t.out / "fix_state.json").read_text())["family"]
    vk = np.asarray(fz.get("seg_v", fz["vk"]), dtype=np.float64)
    sgn_scan = float(fz.get("seg_sgn", -1.0))
    # orient v along the scan's direction of travel: positive offsets lie beyond the root
    if float(v_n[:, 0] @ vk) * sgn_scan < 0.0:
        v_n = -v_n
    out: dict[str, Any] = {"frame": "local at the root, + = beyond the root", "tried": []}
    for side, fracs in (("beyond", (0.5, 0.25)), ("back", (-0.5, -0.25))):
        for frac in fracs:
            a = np.array([frac * big_l])
            br = bordered_solve(t.rj, zs + (v_n @ a) / d, d, u_n, v_n, zs, a, TOL)
            row: dict[str, Any] = {
                "side": side,
                "frac_L": frac,
                "offset": a.tolist(),
                "family_point": bool(br.converged),
            }
            if br.converged:
                t.ckpt_tag = f"neg2_{frac}"
                _, info, ok = t.newton(br.z, 40)
                row.update(started=ok, dr=max(info["dr"]), dv=max(info["dv"]))
            out["tried"].append(row)
            t.log(f"neg2 {side} {frac}: {row}")
            (t.out / "neg2.json").write_text(json.dumps(out))
            if br.converged:
                break
    (t.out / "neg2.json").write_text(json.dumps(out))


def stage_limit(t: Target) -> None:
    """Ideal model: the branch's sigma -> 0 member. From the corrected start, converge going down
    (0.01, 0.005, 0.0025, 0.00125) at the floors and store every point's node epochs, V_inf and z;
    the limit is the Richardson value 2 X(s) - X(2 s) from the two smallest sigmas (limit.json)."""
    st = json.loads((t.out / "start.json").read_text())
    z, s_prev = np.asarray(st["z"]), SIGMA0
    rows = [{"sigma": SIGMA0, "z": z.tolist()}]
    for sg in (0.01, 0.005, 0.0025, 0.00125):
        t.set_sigma(sg)
        zn, info = t.b.solve(t.p, t.scale(z, sg / s_prev), 40, 60, False)
        if not t.b.converged(info):
            t.log(f"limit: no convergence at sigma {sg}")
            break
        z, s_prev = zn, sg
        rows.append({"sigma": sg, "z": z.tolist()})
    for row in rows:
        zz = np.asarray(row["z"])
        t.set_sigma(row["sigma"])
        row["epochs_days"] = [float(zz[7 * k + 6]) / 86400.0 for k in range(len(t.moons))]
        row["vinf"] = t.vinf(zz, row["sigma"])
    out: dict[str, Any] = {"moons": t.moons, "points": rows}
    if len(rows) >= 3:
        a, b = rows[-1], rows[-2]
        out["limit_epochs_days"] = [
            2 * x - y for x, y in zip(a["epochs_days"], b["epochs_days"], strict=True)
        ]
        out["limit_vinf"] = [2 * x - y for x, y in zip(a["vinf"], b["vinf"], strict=True)]
    (t.out / "limit.json").write_text(json.dumps(out))
    t.log(f"limit: epochs {out.get('limit_epochs_days')} vinf {out.get('limit_vinf')}")


def stage_cont(t: Target) -> None:
    st = json.loads((t.out / "start.json").read_text())
    assert st["converged"]
    z = np.asarray(st["z"])
    if hasattr(t, "p"):
        b = t.b
        path = t.out / f"{t.name}.json"
        if not path.exists():
            t.set_sigma(SIGMA0)
            _, _, info = b.residual_and_jac(t.p, z, False)
            pt = b._point(t.p, z, SIGMA0, info)
            path.write_text(json.dumps({"member": t.name, "points": [pt], "stage": "identity"}))
        b.run_member(t.name, time.monotonic() + t.max_wall)
        return
    rb = t.rb
    path = t.out / "sigma_n1.json"
    if not path.exists():
        t.set_sigma(SIGMA0)
        path.write_text(
            json.dumps(
                {"points": [{"sigma": SIGMA0, "z": z.tolist(), **rb.describe(t.c, z, SIGMA0)}]}
            )
        )
    rb.stage_sigma(1)


def main() -> None:
    preflight_search(
        task_no=1055,
        region_id="degenerate-family-seed-correction",
        method=MethodCapability(
            genome="full-rev Jovian chains: ideal gc-1, gc4-19; jup365 GanEur#316, GanCal#1, gc-1",
            corrector="Lyapunov-Schmidt family root + forward-backward shooting, DOP853 + STM",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=[*IDEAL, *REAL], required=True)
    ap.add_argument(
        "--stage",
        choices=["fix", "start", "neg", "cont", "down", "ias15", "limit", "neg2"],
        required=True,
    )
    ap.add_argument("--max-wall", type=float, default=430.0, help="ideal cont: seconds per call")
    args = ap.parse_args()
    t = Target(args.target, ROOT / "data" / f"1055_{args.target}")
    t.max_wall = args.max_wall
    if args.stage == "down":
        t.rb.stage_down(1)
        return
    if args.stage == "ias15":
        t.rb.stage_ias15(1)
        return
    {
        "fix": stage_fix,
        "start": stage_start,
        "neg": stage_neg,
        "cont": stage_cont,
        "limit": stage_limit,
        "neg2": stage_neg2,
    }[args.stage](t)


if __name__ == "__main__":
    main()
