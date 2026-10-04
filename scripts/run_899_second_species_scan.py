"""#899 step 2: the Casoliva et al. 2010 seed-and-continuation scan (a REPRODUCTION run).

For each resonance p-q, a grid of Jacobi constants inside the Eq. 15 interval and both Eq. 16
directions: second-species seed at mu = 1e-6 (Eqs. 14-17), corrected at fixed C_J (optionally
also members of each 2008 Table 2 seed's family at mu = 1e-6, sampled in C_J, and two-arc
chains of returning collision arcs joined by matched lunar flybys on a C_J grid); the
three-step continuation to the paper's mu (mass at fixed C_J first, by default); at that mu, a
walk in C_J both ways and the correction of every crossing of T = 2 pi q; each crossing is
compared with every printed 2010 Table 3 row of the same p-q by the same-orbit test (all y = 0
crossings and their mirror images). The 2008 Table 2 seeds are run first as named seeds.

Output (incremental, flushed): ``<out>/runlog.jsonl`` (one record per seed: stage times,
continuation legs with their end states and stability transitions, crossings, row matches),
``<out>/members.jsonl`` (every continuation member: mu, C, T, periselene, perigee, k_par,
k_perp, Floquet regime, nu_eff and the validity numbers), and ``<out>/progress.txt``. A rerun
with the same ``--out`` skips seeds already in the runlog.

    uv run python scripts/run_899_second_species_scan.py --pq 7-3 2-1 1-2 3-2 --n-c 100 \\
        --out data/runlogs/899_scan
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import time
from typing import Any

import numpy as np

from cyclerfinder.core.floquet_classes import find_family_transitions
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search import earth_moon_resonant_families as emrf
from cyclerfinder.search import second_species_continuation as ssc

_REGION_ID = "earth-moon-casoliva-2010-second-species-mass-continuation-reproduction"
_METHOD = MethodCapability(
    genome=(
        "planar CR3BP p-q resonant second-species seeds at mu = 1e-6 (Casoliva 2008 Eqs. "
        "14-17 and Table 2), continued in the mass ratio to the paper's Earth-Moon value"
    ),
    corrector=(
        "second_species_continuation: multiple shooting with the planar Levi-Civita "
        "propagator, continue_mu / continue_jacobi / three_step / resonant_crossings"
    ),
    capability_tags=frozenset(
        {"cr3bp", "planar", "earth-moon", "second-species", "mu-continuation"}
    ),
    git_sha="working-tree",
)


def _diag_dict(d: ssc.Diagnostics) -> dict[str, Any]:
    return {
        "mu": d.mu,
        "C": d.jacobi,
        "T": d.period,
        "periselene": d.periselene,
        "perigee": d.perigee,
        "k_par": d.k_par,
        "k_perp": d.k_perp,
        "floquet_regime": d.floquet_regime,
        "nu_eff": d.nu_eff,
        "validity_speed": d.validity_speed,
        "validity_rp": d.validity_rp,
        "impact_moon": d.impact_moon,
        "impact_earth": d.impact_earth,
    }


def _transitions(leg: ssc.LegResult, param: str) -> list[dict[str, Any]]:
    if len(leg.members) < 2:
        return []
    mons = [ssc.reduced_monodromy(m.diag.k_par, m.diag.k_perp) for m in leg.members]
    if param == "mu":
        params = [math.log(m.diag.mu) for m in leg.members]
    else:
        params = [m.diag.jacobi for m in leg.members]
    out = []
    for tr in find_family_transitions(mons, params):
        out.append(
            {
                "kind": tr.kind,  # on the normalised scale: k = +-1 is Casoliva's +-2
                "param": "log_mu" if param == "mu" else "C",
                "lo": tr.lo,
                "hi": tr.hi,
                "step": tr.hi - tr.lo,
                "region_before": tr.region_before,
                "region_after": tr.region_after,
            }
        )
    return out


class Writer:
    def __init__(self, out: pathlib.Path) -> None:
        out.mkdir(parents=True, exist_ok=True)
        self.runlog = (out / "runlog.jsonl").open("a")
        self.members = (out / "members.jsonl").open("a")
        self.progress = out / "progress.txt"
        self.seed_id = ""

    def member(self, m: ssc.Member) -> None:
        rec = {"seed": self.seed_id, "leg": m.leg, "step": m.step, **_diag_dict(m.diag)}
        self.members.write(json.dumps(rec) + "\n")
        self.members.flush()

    def record(self, rec: dict[str, Any]) -> None:
        self.runlog.write(json.dumps(rec) + "\n")
        self.runlog.flush()


def _done(out: pathlib.Path) -> set[str]:
    path = out / "runlog.jsonl"
    if not path.exists():
        return set()
    done = set()
    for line in path.read_text().splitlines():
        try:
            done.add(json.loads(line)["seed"])
        except (ValueError, KeyError):
            continue
    return done


def run_seed(
    seed_id: str,
    p: int,
    q: int,
    make_orbit: Any,
    args: argparse.Namespace,
    w: Writer,
    all_hits: list[tuple[str, ssc.MSOrbit]],
) -> dict[str, Any]:
    w.seed_id = seed_id
    rec: dict[str, Any] = {"seed": seed_id, "p": p, "q": q, "t_start": time.time()}
    t0 = time.time()
    try:
        orbit = make_orbit()
    except ssc._STEP_FAILURES as exc:
        rec.update(stage="seed_failed", error=str(exc)[:200], t_seed=time.time() - t0)
        return rec
    rec["seed_orbit"] = _diag_dict(ssc.diagnose(orbit))
    rec["t_seed"] = time.time() - t0
    t1 = time.time()
    path = ssc.three_step(
        orbit,
        args.mu_target,
        first_fix=args.first_fix,
        max_switches=args.max_switches,
        log=w.member,
    )
    rec["t_continuation"] = time.time() - t1
    rec["legs"] = [
        {
            "leg": name,
            "reason": leg.reason,
            "n": len(leg.members),
            "rejected_newton": leg.rejected_newton,
            "rejected_jump": leg.rejected_jump,
            "end": _diag_dict(leg.members[-1].diag) if leg.members else None,
            "transitions": _transitions(leg, "C" if name == "C@mu" else "mu"),
        }
        for name, leg in path.legs
    ]
    if not path.reached or path.final is None:
        rec["stage"] = "continuation_short"
        return rec
    t2 = time.time()
    hits, walks = ssc.resonant_crossings(path.final, q, max_steps=args.walk_steps, log=w.member)
    rec["t_walk"] = time.time() - t2
    rec["walks"] = {
        str(k): {
            "reason": v.reason,
            "n": len(v.members),
            "C_range": [
                min(m.diag.jacobi for m in v.members),
                max(m.diag.jacobi for m in v.members),
            ]
            if v.members
            else None,
            "transitions": _transitions(v, "C"),
        }
        for k, v in walks.items()
    }
    rows = [r for r in emrf.TABLE3_ROWS if (r.p, r.q) == (p, q)]
    rec["hits"] = []
    for h in hits:
        dup = next((sid for sid, o in all_hits if ssc.same_orbit(o, h.orbit)), None)
        matches = {}
        for row in rows:
            dist, mirrored, _ = ssc.orbit_distance(
                h.orbit, np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
            )
            matches[row.designation] = {
                "dist": dist,
                "mirrored": mirrored,
                "dC": h.diag.jacobi - row.c_j,
                "k_printed": row.k,
                "k_casoliva": ssc.casoliva_k(row.designation, h.diag),
            }
        rec["hits"].append(
            {
                "orbit": _diag_dict(h.diag),
                "state0": h.orbit.nodes[0].tolist(),
                "duplicate_of": dup,
                "matches": matches,
            }
        )
        if dup is None:
            all_hits.append((seed_id, h.orbit))
    rec["stage"] = "done"
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pq", nargs="+", default=["7-3", "2-1", "1-2", "3-2"])
    ap.add_argument("--n-c", type=int, default=100, help="C_J grid points per p-q")
    ap.add_argument("--branches", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--first-fix", choices=["jacobi", "period"], default="jacobi")
    ap.add_argument("--max-switches", type=int, default=6)
    ap.add_argument("--walk-steps", type=int, default=300)
    ap.add_argument("--mu-target", type=float, default=ssc.CASOLIVA_MU_2010)
    ap.add_argument("--no-table2", action="store_true", help="skip the 2008 Table 2 seeds")
    ap.add_argument(
        "--table2-walk-dc",
        type=float,
        default=0.0,
        help="also seed from members of each Table 2 seed's family at mu = 1e-6, every dC",
    )
    ap.add_argument("--table2-walk-steps", type=int, default=200)
    ap.add_argument(
        "--n-chain", type=int, default=60, help="C_J grid points for two-arc chain seeds"
    )
    ap.add_argument("--max-minutes", type=float, default=math.inf)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("data/runlogs/899_scan"))
    args = ap.parse_args()

    pqs = [tuple(int(v) for v in s.split("-")) for s in args.pq]
    jobs: list[tuple[str, int, int, Any]] = []
    if not args.no_table2:
        for s in ssc.TABLE2_SEEDS:
            if (s.p, s.q) in pqs:

                def mk(s: ssc.Table2Seed = s) -> ssc.MSOrbit:
                    start = ssc.MSOrbit(ssc.SEED_MU, np.array([s.state_project()]), s.period)
                    return ssc.correct(ssc.renode(start, 12), fix_jacobi=s.c_j)

                jobs.append((f"table2:{s.designation}", s.p, s.q, mk))
    if args.table2_walk_dc > 0.0:
        # members of each 2008 seed's family at mu = 1e-6, sampled every table2_walk_dc in C
        for s in ssc.TABLE2_SEEDS:
            if (s.p, s.q) not in pqs:
                continue
            start = ssc.MSOrbit(ssc.SEED_MU, np.array([s.state_project()]), s.period)
            base = ssc.correct(ssc.renode(start, 12), fix_jacobi=s.c_j)
            for direction in (1.0, -1.0):
                last_c = base.jacobi
                leg = ssc.continue_jacobi(
                    base, direction, max_steps=args.table2_walk_steps, ds_max=0.05
                )
                for k, m in enumerate(leg.members):
                    if abs(m.diag.jacobi - last_c) >= args.table2_walk_dc:
                        last_c = m.diag.jacobi
                        sid = f"walk:{s.designation}:{direction:+.0f}:{k}:C={m.diag.jacobi:.6f}"
                        jobs.append((sid, s.p, s.q, (lambda o=m.orbit: o)))
                print(
                    f"family walk of {s.designation} at mu = 1e-6, direction {direction:+.0f}: "
                    f"{leg.reason}, {len(leg.members)} members",
                    flush=True,
                )
    for p, q in pqs:
        if args.n_c <= 0:
            continue
        lo, hi = ssc.jacobi_interval(p, q)
        grid = lo + (hi - lo) * (np.arange(args.n_c) + 0.5) / args.n_c
        for c in grid:
            for br in args.branches:

                def mk2(p: int = p, q: int = q, c: float = float(c), br: int = br) -> ssc.MSOrbit:
                    return ssc.corrected_seed(p, q, c, br)

                jobs.append((f"grid:{p}-{q}:C={c:.6f}:b{br}", p, q, mk2))
    for p, q in pqs:
        if args.n_chain <= 0:
            continue
        lo, hi = ssc.jacobi_interval(p, q)
        grid = lo + (hi - lo) * (np.arange(args.n_chain) + 0.5) / args.n_chain
        for arcs in ssc.two_arc_chains(p, q):
            label = "+".join(f"({a.i},{a.j},{a.sign:+d})" for a in arcs)
            for c in grid:
                try:
                    for a in arcs:
                        ssc.arc_relative_velocity(a, float(c))
                except ValueError:
                    continue  # an arc of the chain does not reach the Moon at this C

                def mk3(
                    arcs: tuple[ssc.ReturningArc, ...] = arcs, c: float = float(c)
                ) -> ssc.MSOrbit:
                    return ssc.corrected_chain(list(arcs), c)

                jobs.append((f"chain:{p}-{q}:{label}:C={c:.6f}", p, q, mk3))

    preflight_search(
        task_no=899,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=len(jobs),
        override_reason=(
            "reproduction of published orbits (Casoliva et al. 2010 Table 3) by the published "
            "seed-and-continuation method, a validation of the #899 pipeline, not a discovery "
            "sweep; seeds per resonance are a fixed C_J grid"
        ),
    )

    w = Writer(args.out)
    done = _done(args.out)
    all_hits: list[tuple[str, ssc.MSOrbit]] = []
    todo = [j for j in jobs if j[0] not in done]
    t_start = time.time()
    print(f"{len(jobs)} seeds, {len(done)} already done, {len(todo)} to run", flush=True)
    for i, (sid, p, q, mk) in enumerate(todo):
        if (time.time() - t_start) / 60.0 > args.max_minutes:
            print("max-minutes reached; rerun with the same --out to resume", flush=True)
            break
        rec = run_seed(sid, p, q, mk, args, w, all_hits)
        rec["t_total"] = time.time() - rec["t_start"]
        w.record(rec)
        elapsed = time.time() - t_start
        eta = elapsed / (i + 1) * (len(todo) - i - 1)
        best = []
        for h in rec.get("hits", []):
            closest = min(h["matches"].items(), key=lambda kv: kv[1]["dist"], default=None)
            if closest is not None:
                best.append(f"{closest[0]}:{closest[1]['dist']:.1e}")
        line = (
            f"{time.strftime('%Y-%m-%dT%H:%M:%S')} [{i + 1}/{len(todo)}] {sid} "
            f"{rec.get('stage')} {rec['t_total']:.1f}s hits={best} "
            f"elapsed={elapsed / 60:.1f}min eta={eta / 60:.1f}min"
        )
        print(line, flush=True)
        w.progress.write_text(line + "\n")


if __name__ == "__main__":
    main()
