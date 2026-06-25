"""#466 energy-moving QP-GMOS continuation campaign — the decisive descent.

Walks the #290 smoke QP-torus DOWN the #296 parent (1,1) energy ladder from
C_J~3.12785 (parent step 112) toward the #320 SILVER Bracket-2 region at
C_J~3.03196 (parent step 8), re-converging the QP-torus at each parent member.
Writes incremental JSONL (one member per row) + a timestamped progress line so the
multi-hour walk is monitorable + resumable. Report-only: NO catalogue writeback; any
candidate quasi-cycler is flagged for human gauntlet review, NEVER self-admitted.

    uv run python scripts/run_466_energy_walk.py --stride 1 --target-step 8
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
from typing import Any

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.genome.qp_tori_energy_walk import (
    EnergyWalkMember,
    walk_energy,
)

ROOT = Path(__file__).resolve().parents[1]

EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=EM_MU, primary="earth", secondary="moon", l_km=EM_L_KM, t_s=EM_T_S)


def _build_smoke_seed() -> QPTorus:
    sub = ROOT / "data" / "family_296_3d_subfamilies_299.jsonl"
    par = ROOT / "data" / "family_296_3d_em_11.jsonl"
    bracket: dict[str, Any] | None = None
    with sub.open() as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            br = d.get("bracket")
            if br and br.get("classification") == "neimark_sacker":
                bracket = br
                break
    assert bracket is not None
    parent: dict[str, Any] | None = None
    with par.open() as fh:
        for line in fh:
            o = json.loads(line)
            if o.get("type") == "header":
                continue
            if o.get("step_index") == int(bracket["step_a"]):
                parent = o
                break
    assert parent is not None
    system = _em_system()
    return correct_qp_torus(
        system,
        np.asarray(parent["state_nd"], dtype=np.float64),
        float(parent["T_TU"]),
        (
            complex(bracket["eig_a_re"], bracket["eig_a_im"]),
            complex(bracket["eig_b_re"], bracket["eig_b_im"]),
        ),
        k=int(bracket["k"]),
        n_long=16,
        n_trans=2,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="466_campaign_smoke",
    )


def _member_row(m: EnergyWalkMember) -> dict[str, Any]:
    return {
        "jacobi": m.jacobi,
        "parent_step": m.extras.get("parent_step"),
        "rho": m.rho,
        "freq_ratio": m.freq_ratio,
        "is_irrational": m.is_practically_irrational,
        "residual_norm": m.residual_norm,
        "independent_residual": m.extras.get("independent_residual"),
        "near_resonance": (
            None
            if m.near_resonance is None
            else {"p": m.near_resonance.p, "q": m.near_resonance.q, "d": m.near_resonance.distance}
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1, help="parent step_index stride")
    ap.add_argument("--target-step", type=int, default=8, help="parent step to descend to")
    ap.add_argument("--max-steps", type=int, default=200)
    ap.add_argument("--out", type=str, default="data/family_466_energy_walk.jsonl")
    args = ap.parse_args()

    t0 = dt.datetime.now()
    print(
        f"[{t0.isoformat(timespec='seconds')}] building smoke seed (parent step 112)...", flush=True
    )
    torus = _build_smoke_seed()
    print(
        f"[{dt.datetime.now().isoformat(timespec='seconds')}] seed built in "
        f"{(dt.datetime.now() - t0).total_seconds():.0f}s; descending toward step "
        f"{args.target_step} (stride {args.stride})",
        flush=True,
    )

    out = ROOT / args.out
    fh = out.open("w")

    def on_step(m: EnergyWalkMember) -> None:
        fh.write(json.dumps(_member_row(m)) + "\n")
        fh.flush()
        ts = dt.datetime.now().isoformat(timespec="seconds")
        indep = float(m.extras.get("independent_residual", float("nan")))
        print(
            f"[{ts}] step={int(m.extras.get('parent_step', -1))} C_J={m.jacobi:.6f} "
            f"rho={m.rho:.6f} ratio={m.freq_ratio:.6f} irr={m.is_practically_irrational} "
            f"res={m.residual_norm:.2e} indep={indep:.2e}",
            flush=True,
        )

    fam = walk_energy(
        torus,
        direction="down",
        step_stride=args.stride,
        max_steps=args.max_steps,
        target_step=args.target_step,
        seed_step=112,
        on_step=on_step,
    )
    fh.close()
    print(
        f"DONE: {len(fam.members)} members, terminated={fam.terminated_reason}, "
        f"stop_C_J={fam.stop_cj:.6f}, "
        f"{len(fam.resonance_crossings)} resonance crossings",
        flush=True,
    )
    cj = [m.jacobi for m in fam.members]
    print(
        f"C_J span: {min(cj):.6f} .. {max(cj):.6f} (delta {max(cj) - min(cj):.6f}); "
        f"reached SILVER region (C_J<=3.0320)? {min(cj) <= 3.0320}",
        flush=True,
    )


if __name__ == "__main__":
    main()
