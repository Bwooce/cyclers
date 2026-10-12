"""#1054 amendment (2026-10-12): rounding sensitivity of the free-C multiple-shooting landing point.

Why: after the macOS 27.0.1 update (2026-10-12 09:38 AEDT) two converged seeds in
tests/scripts/test_1054_seed_design.py no longer reproduce their recorded C and T bit-for-bit
(relative 6e-8 and 1.6e-9), with the same code. Note: docs/notes/2026-10-11-1054-seed-design-fix.md
sec. 10.1.

    uv run python scripts/probe_1054_float_sensitivity.py bound     # derived drift bound vs record
    uv run python scripts/probe_1054_float_sensitivity.py ulp       # 1-ulp ensembles (C spread)
    uv run python scripts/probe_1054_float_sensitivity.py control   # wrong-row positive control
    uv run python scripts/probe_1054_float_sensitivity.py stable    # status/it/n_halved under 1 ulp

The drift bound (`drift_bound`), used by the test:

* With C free, the shooting system (4n continuity rows + 1 perigee-phase row, 5n unknowns) has a
  family of solutions, so the stopping rule (residual < CONT_TOL) does not constrain the landing
  point along the family: two runs can both stop and still differ there.
* Each Newton step is the minimum-norm least-squares solution of J dx = -f. An error df in the
  evaluated residual moves the step by J+ df, at most |df|_2 / sigma_min(J) (smallest singular value
  that lstsq keeps). Near convergence the Jacobi-integral redundancy makes sigma_min collapse
  (about 4e-8 at the last step of the two seeds), and the component of the error along the family
  is never removed by later steps.
* What differs between two environments is the residual's ROUNDING, not its truncation error:
  both run the same algorithm and (essentially) the same integrator steps, so the truncation error
  (what rtol = atol = 1e-12 controls) is common to both and cancels. The per-component error eps
  is therefore the arc end point's rounding noise, measured at the converged nodes as
  max |arc(s + d) - arc(s) - Phi d| over 1-ulp perturbations d (`arc_noise`);
  |df|_2 <= sqrt(4n) eps.
  (The truncation-scale choice eps_arc = atol + rtol max|y| gives a bound 30-70 times larger that
  accepts a wrong-perigee-row regression; `control` prints both.)
* Summed over the steps taken: |dx| <= sum_k sqrt(4n) eps / sigma_min,k. Then
  |dC| <= |grad C(s0)|_2 (|dx| + eps) (C is evaluated after propagation to the first perigee, where
  it is conserved to the same noise), and |dT| <= sqrt(n) |dx| when k_min = p (T is the sum of the
  n durations).

The test holds these bounds as fixed per-seed constants (printed by `bound`), never recomputed from
the solve under test: a regression that lowered sigma_min would otherwise widen its own tolerance.
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sh = importlib.import_module("run_1000_shooting")

RTOL = ATOL = 1e-12  # run_1000_shooting.arc (DOP853)
EPS = float(np.finfo(float).eps)

# the converged seeds whose C and T the test compares: (label, p, q, sense, rp, omega, record)
GAPS_FILE = REPO / "data" / "1000_complement" / "shooting" / "gaps.jsonl"
ARM_E_FILE = REPO / "data" / "1054_seed_design" / "arm_E.jsonl"


def _rec(path: Path, keys: tuple[str, ...], seed: tuple[int, ...]) -> dict[str, Any]:
    for line in path.read_text().splitlines():
        r = json.loads(line)
        if tuple(r[k] for k in keys) == seed:
            return dict(r)
    raise KeyError(seed)


def cases() -> list[tuple[str, int, int, int, float, float, dict[str, Any]]]:
    g = _rec(GAPS_FILE, ("p", "q", "sense", "ir", "iw"), (3, 1, 1, 7, 0))
    e = _rec(ARM_E_FILE, ("p", "q", "sense", "irn", "k"), (7, 3, -1, 3, 11))
    return [
        ("gaps (3,1,1,7,0) legacy", 3, 1, 1, sh.RP[7], sh.OMEGA[0], g),
        ("arm E (7,3,-1,3,11) f1", 7, 3, -1, sh.RP_F1[3], sh.OMEGA[11], e),
    ]


@contextmanager
def capture_jacobians(store: list[np.ndarray]) -> Iterator[None]:
    """Record every matrix `shoot` passes to np.linalg.lstsq (one per Newton step taken)."""
    orig = np.linalg.lstsq

    def wrapped(a: Any, b: Any, rcond: Any = None) -> Any:
        store.append(np.array(a, copy=True))
        return orig(a, b, rcond=rcond)

    np.linalg.lstsq = wrapped
    try:
        yield
    finally:
        np.linalg.lstsq = orig


def sigma_min_kept(a: np.ndarray) -> float:
    """Smallest singular value np.linalg.lstsq keeps (rcond=None: cutoff eps max(M, N) s_max)."""
    s = np.linalg.svd(a, compute_uv=False)
    return float(s[s > EPS * max(a.shape) * s[0]].min())


def drift_bound(p: int, q: int, sense: int, rp: float, om: float) -> dict[str, Any]:
    """Solve the seed (free C, no F3) and return the record plus the derived drift bounds."""
    jacs: list[np.ndarray] = []
    with capture_jacobians(jacs):
        rec = sh.solve_values({}, p, q, sense, rp, om, None)
    if rec["status"] != "converged":
        return {"rec": rec}
    n = p
    nodes = np.array(rec["nodes"])
    sig = [sigma_min_kept(a) for a in jacs]
    _, g = sh.jac_c(np.array(rec["s0"]))
    gn = float(np.linalg.norm(g))
    out: dict[str, Any] = {"rec": rec, "sigma_min": sig}
    for name, eps in (
        ("rtol", ATOL + RTOL * float(np.max(np.abs(nodes)))),
        ("noise", arc_noise(rec)),
    ):
        dx = sum(math.sqrt(4 * n) * eps / s for s in sig)
        out[name] = {"eps": eps, "dx": dx, "tol_C": gn * (dx + eps), "tol_T": math.sqrt(n) * dx}
    return out


def arc_noise(rec: dict[str, Any], n_samples: int = 10) -> float:
    """Max rounding noise of the arc end points at the converged nodes:
    |arc(s + d) - arc(s) - Phi d| over 1-ulp random-sign perturbations d (seeded)."""
    rng = np.random.default_rng(1054)
    worst = 0.0
    for x, tau in zip(rec["nodes"], rec["taus"], strict=True):
        s = np.array(x)
        res0 = sh.arc(s, tau)
        assert res0 is not None
        sf, phi = res0
        for _ in range(n_samples):
            d = s * EPS * rng.choice([-1.0, 1.0], size=4)
            res1 = sh.arc(s + d, tau)
            assert res1 is not None
            worst = max(worst, float(np.max(np.abs(res1[0] - sf - phi @ d))))
    return worst


def cmd_bound() -> None:
    for label, p, q, sense, rp, om, r in cases():
        b = drift_bound(p, q, sense, rp, om)
        rec = b["rec"]
        dc, dt = abs(rec["C"] - r["C"]), abs(rec["T"] - r["T"])
        sig = ", ".join(f"{s:.2e}" for s in b["sigma_min"])
        print(f"{label}: status {rec['status']} it {rec['it']} k_min {rec['k_min']}")
        print(f"  sigma_min per step: {sig}")
        for name in ("rtol", "noise"):
            x = b[name]
            print(f"  [{name}] eps {x['eps']:.2e}  |dx| bound {x['dx']:.2e}")
            for v, d in (("C", dc), ("T", dt)):
                tol = x[f"tol_{v}"]
                print(f"    |d{v}| record vs now {d:.3e}  bound {tol:.3e}  ratio {d / tol:.3f}")


def _perturbed_solve(job: tuple[int, int, int, float, float, int, str]) -> dict[str, Any]:
    """Solve with 1-ulp random sign perturbations of the skeleton nodes ("seed") or of every Newton
    step ("step"); k is the random seed."""
    p, q, sense, rp, om, k, kind = job
    rng = np.random.default_rng(k)
    orig_nodes: Callable[..., Any] = sh.skeleton_nodes
    orig_lstsq = np.linalg.lstsq

    def nodes_p(*a: Any) -> Any:
        nodes, tp = orig_nodes(*a)
        return [x * (1 + EPS * rng.choice([-1.0, 1.0], size=4)) for x in nodes], tp

    def lstsq_p(a: Any, b: Any, rcond: Any = None) -> Any:
        x = orig_lstsq(a, b, rcond=rcond)
        return (x[0] * (1 + EPS * rng.choice([-1.0, 1.0], size=x[0].shape)), *x[1:])

    if kind == "seed":
        sh.skeleton_nodes = nodes_p  # type: ignore[attr-defined]
    else:
        np.linalg.lstsq = lstsq_p
    try:
        out: dict[str, Any] = sh.solve_values({}, p, q, sense, rp, om, None)
        return out
    finally:
        sh.skeleton_nodes = orig_nodes  # type: ignore[attr-defined]
        np.linalg.lstsq = orig_lstsq


def cmd_ulp(n_runs: int) -> None:
    with Pool(3) as pool:
        for label, p, q, sense, rp, om, r in cases():
            base = sh.solve_values({}, p, q, sense, rp, om, None)
            for kind in ("seed", "step"):
                jobs = [(p, q, sense, rp, om, k, kind) for k in range(n_runs)]
                outs = pool.map(_perturbed_solve, jobs)
                cs = np.array([o["C"] for o in outs if o["status"] == "converged"])
                its = sorted({o["it"] for o in outs})
                print(
                    f"{label} 1-ulp {kind}: {len(cs)}/{n_runs} converged, it {its}, "
                    f"max|C - C_now| {np.max(np.abs(cs - base['C'])):.3e}, std {cs.std():.3e}, "
                    f"|C_rec - C_now| {abs(r['C'] - base['C']):.3e}",
                    flush=True,
                )


def cmd_control(tols: list[float] | None) -> None:
    """A real regression (wrong perigee row or omega) must fall outside the C tolerance: the test's
    constants when given (--tol-c, one per case), else the noise bound; the rtol bound is shown."""
    for i_case, (label, p, q, sense, rp, om, r) in enumerate(cases()):
        b = drift_bound(p, q, sense, rp, om)
        tol = tols[i_case] if tols else b["noise"]["tol_C"]
        print(f"{label}: rtol-scale tol_C {b['rtol']['tol_C']:.3e} (shown, not used)")
        alts = {
            "legacy RP rows": sh.RP,
            "F1 rows": sh.RP_F1,
        }
        ir = (sh.RP if "legacy" in label else sh.RP_F1).index(rp)
        iw = sh.OMEGA.index(om)
        print(f"{label}: tol_C {tol:.3e}, record C {r['C']!r}")
        trials: list[tuple[str, float, float]] = []
        for name, rows in alts.items():
            for j in (ir - 1, ir, ir + 1):
                if 0 <= j < len(rows) and rows[j] != rp:
                    trials.append((f"{name}[{j}], omega[{iw}]", rows[j], om))
        for jw in (iw - 1, iw + 1):
            trials.append(
                (f"same row, omega[{jw % len(sh.OMEGA)}]", rp, sh.OMEGA[jw % len(sh.OMEGA)])
            )
        for name, rp2, om2 in trials:
            out = sh.solve_values({}, p, q, sense, rp2, om2, None)
            if out["status"] == "converged":
                dc = abs(out["C"] - r["C"])
                verdict = "REJECTED" if dc > tol else "ACCEPTED (control fails)"
                if dc <= b["rtol"]["tol_C"]:
                    verdict += "; the rtol-scale bound would accept it"
                print(f"  {name}: converged it {out['it']}, |dC| {dc:.3e} -> {verdict}")
            else:
                print(f"  {name}: {out['status']} it {out['it']} -> rejected by status")


def _stable_job(job: tuple[str, int, int, int, float, float, float | None, int]) -> str:
    """k = -1: the unperturbed reference; k >= 0: 1-ulp perturbed skeleton nodes."""
    name, p, q, sense, rp, om, drop, k = job
    rng = np.random.default_rng(max(k, 0))
    orig_nodes = sh.skeleton_nodes

    def nodes_p(*a: Any) -> Any:
        nodes, tp = orig_nodes(*a)
        return [x * (1 + EPS * rng.choice([-1.0, 1.0], size=4)) for x in nodes], tp

    if k >= 0:
        sh.skeleton_nodes = nodes_p  # type: ignore[attr-defined]
    try:
        out = sh.solve_values({}, p, q, sense, rp, om, None, drop)
    finally:
        sh.skeleton_nodes = orig_nodes  # type: ignore[attr-defined]
    return json.dumps([name, k, [out.get(x) for x in ("status", "it", "n_halved")]])


def cmd_stable(n_runs: int) -> None:
    """Do status, iteration count and halving count of every test seed survive 1-ulp seed noise?"""
    sys.path.insert(0, str(REPO / "tests" / "scripts"))
    t = importlib.import_module("test_1054_seed_design")
    seeds: list[tuple[str, int, int, int, float, float, float | None]] = []
    for p, q, sense, ir, iw in t.REPRO_SEEDS + t.GAPS_SEEDS:
        seeds.append(
            (f"legacy {(p, q, sense, ir, iw)}", p, q, sense, sh.RP[ir], sh.OMEGA[iw], None)
        )
    for p, q, sense, ir, iw in t.F3_SEEDS:
        seeds.append((f"F3 {(p, q, sense, ir, iw)}", p, q, sense, sh.RP[ir], sh.OMEGA[iw], 0.2))
    for p, q, sense, ir, iw in t.ARM_E_SEEDS:
        seeds.append((f"f1 {(p, q, sense, ir, iw)}", p, q, sense, sh.RP_F1[ir], sh.OMEGA[iw], None))
    jobs = [(*sd, k) for sd in seeds for k in range(-1, n_runs)]
    ref: dict[str, Any] = {}
    flips = 0
    with Pool(3) as pool:
        for i, line in enumerate(pool.imap(_stable_job, jobs), 1):
            name, k, key = json.loads(line)
            if k < 0:
                ref[name] = key
            elif key != ref[name]:
                flips += 1
                print(f"FLIP {name} k={k}: {ref[name]} -> {key}", flush=True)
            if i % 25 == 0:
                print(f"  {i}/{len(jobs)} solves, {flips} flips", flush=True)
    print(f"stable: {len(seeds)} seeds x {n_runs} perturbed solves, {flips} flips")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["bound", "ulp", "control", "stable"])
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--tol-c", type=float, nargs="*", default=None)
    args = ap.parse_args()
    if args.cmd == "bound":
        cmd_bound()
    elif args.cmd == "ulp":
        cmd_ulp(args.n)
    elif args.cmd == "control":
        cmd_control(args.tol_c)
    else:
        cmd_stable(args.n)


if __name__ == "__main__":
    main()
