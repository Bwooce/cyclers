"""#1023 EGGIE in the consistent ideal Galilean model (pre-registration:
docs/notes/2026-10-08-1023-jovian-void-rerun.md sec. 3).

    uv run python scripts/run_1023_eggie.py --stage pc      # step 1: patched-conic roots

Steps 2-3 (continuous gravity) were not run: no patched-conic root near Table 4 (note sec. 3.6).

Checkpoints in data/1023_eggie/; live log data/1023_eggie/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2
from cyclerfinder.search.resonant_conic import (
    ideal_moon_smas,
    ideal_t_syn,
    ideal_t_syn_consistent,
)
from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    LambertLeg,
    correct_dates,
    cycle_flybys,
    encounter_self_consistency,
    gate_cycle,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1023_eggie"
DAY = 86400.0
IO, EUR, GAN = "Io", "Europa", "Ganymede"
MOONS = (IO, EUR, GAN)
NODE_MOONS = (EUR, GAN, GAN, IO)
TABLE4 = {EUR: 9.12, GAN: 7.07, IO: 8.38}
ALT_FLOOR_KM = 25.0
Arr = NDArray[np.float64]


def _log(msg: str) -> None:
    (OUT / "live").mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with (OUT / "live" / "runlog.txt").open("a") as fh:
        fh.write(line + "\n")


def system(variant: str = "consistent", laplace_deg: float = 180.0) -> tuple[CircularSystem, float]:
    """The ideal model and the cycle period T for a variant (note secs. 3.1 and 6)."""
    if variant == "tsyn705":  # T_syn = 7.05 d as printed; smas rebuilt for it (sec. 6 (i))
        t_syn = 7.05 * DAY
        d = math.radians(5.2)
        k = {IO: 8.0, EUR: 4.0, GAN: 2.0}
        smas = {
            m: (MU_JUPITER_KM3_S2 / ((k[m] * math.pi + d) / t_syn) ** 2) ** (1.0 / 3.0)
            for m in MOONS
        }
    elif variant == "coded":  # the OLD coded model (sec. 6 (iii))
        t_syn = ideal_t_syn()
        smas = ideal_moon_smas()
    else:
        t_syn = ideal_t_syn_consistent()
        smas = ideal_moon_smas()
    # Io = Europa = 0; Ganymede from the Laplace angle lambda_I - 3 lambda_E + 2 lambda_G.
    th0 = {IO: 0.0, EUR: 0.0, GAN: math.radians(laplace_deg) / 2.0}
    bodies = {
        m: (smas[m], 2.0 * math.pi * math.sqrt(smas[m] ** 3 / MU_JUPITER_KM3_S2), th0[m])
        for m in MOONS
    }
    return CircularSystem(MU_JUPITER_KM3_S2, bodies), 4.0 * t_syn


def stage_pc(variant: str = "consistent") -> None:
    angles = [float(a) for a in range(0, 360, 15)] if variant == "laplace" else [180.0]
    out = []
    for ang in angles:
        circ, t = system("consistent" if variant == "laplace" else variant, ang)
        roots = search(circ, t)
        for q in roots:
            q["laplace_deg"] = ang
        out.extend(roots)
        near = roots[0]["dist_table4"] if roots else float("nan")
        _log(f"pc {variant} laplace={ang}: {len(roots)} roots, nearest {near:.3f} km/s")
    out.sort(key=lambda q: q["dist_table4"])
    name = "pc_roots.json" if variant == "consistent" else f"pc_roots_{variant}.json"
    (OUT / name).write_text(json.dumps(out, indent=1))
    gp = [q for q in out if q["gate25"] == "pass"]
    best = out[0]["dist_table4"] if out else float("nan")
    _log(
        f"pc {variant}: {len(out)} roots; nearest {best:.3f} "
        f"km/s ({json.dumps(out[0]['vinf']) if out else '-'}); gate-passing {len(gp)}, nearest "
        f"{gp[0]['dist_table4'] if gp else float('nan'):.3f}"
    )


def search(circ: CircularSystem, t: float) -> list[dict]:  # type: ignore[type-arg]
    seeds_days = np.array([0.0, 1.59, 10.19, 17.53])
    revs = [(0, "single"), (1, "low"), (1, "high"), (2, "low"), (2, "high")]
    roots = []
    for r1, r2, r3, r4 in itertools.product(revs[:1], revs[1:3], revs, revs):
        cyc = Cycle(
            (
                LambertLeg(EUR, GAN, *r1),
                LambertLeg(GAN, GAN, *r2),
                LambertLeg(GAN, IO, *r3),
                LambertLeg(IO, EUR, *r4),
            ),
            t,
        )
        for shift in np.linspace(0.0, t / 4.0, 12, endpoint=False):
            sol = correct_dates(circ, cyc, seeds_days * DAY + shift, tol_kms=1e-9)
            if not sol.converged:
                continue
            fl = cycle_flybys(circ, cyc, sol.x)
            if fl is None:
                continue
            miss = encounter_self_consistency(circ, cyc, sol.x)
            if miss >= 1.0:
                continue
            vinf = {}
            for f in fl:
                vinf.setdefault(f.body, f.vinf_kms)
            dist = max(abs(vinf[m] - TABLE4[m]) for m in MOONS)
            gate25 = gate_cycle(circ, fl, alt_floor_km=ALT_FLOOR_KM)
            gate = gate_cycle(circ, fl)
            key = f"{r1}{r2}{r3}{r4}"
            row = {
                "legs": key,
                "x_days": (sol.x / DAY).tolist(),
                "vinf": vinf,
                "dist_table4": dist,
                "turns_deg": [f.turn_deg for f in fl],
                "flyby_bodies": [f.body for f in fl],
                "gate25": gate25.status,
                "gate_project": gate.status,
                "miss_km": miss,
            }
            if not any(
                abs(row["vinf"][IO] - q["vinf"][IO]) < 1e-6 and q["legs"] == key for q in roots
            ):
                roots.append(row)
    roots.sort(key=lambda q: q["dist_table4"])
    return roots


def main() -> None:
    preflight_search(
        task_no=1023,
        region_id="eggie-consistent-ideal-model",
        method=MethodCapability(
            genome="EGGIE (one structure), consistent ideal Galilean model",
            corrector="correct_dates (patched conic) then forward-backward shooting, DOP853 + STM",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["pc"], required=True)
    ap.add_argument(
        "--variant", choices=["consistent", "tsyn705", "laplace", "coded"], default="consistent"
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "pc":
        stage_pc(args.variant)


if __name__ == "__main__":
    main()
