"""#333 QP-GMOS family continuation campaign.

Seeds off the #290 smoke torus (the first accepted #299 Neimark-Sacker bracket,
k=4, off the #296 Braik-Ross (1,1) Earth-Moon family), walks both directions in
Jacobi energy, writes incremental JSONL (one member per row) via the on_step
callback so a multi-hour walk is monitorable + resumable, and prints a per-step
timestamped progress line. Report-only: NO catalogue writeback.

    uv run python scripts/run_333_qp_family.py --seed smoke --ds 5e-3 --max-steps 200
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
from cyclerfinder.genome.qp_tori_arclength import (
    QPTorusFamilyMember,
    continue_qp_family_arclength,
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
        notes="333_campaign_smoke",
    )


def _member_row(m: QPTorusFamilyMember) -> dict[str, Any]:
    return {
        "jacobi": m.jacobi,
        "arclength_s": m.arclength_s,
        "rho": m.rho,
        "freq_ratio": m.freq_ratio,
        "is_irrational": m.is_practically_irrational,
        "fold_index": m.fold_index,
        "residual_norm": m.residual_norm,
        "near_resonance": (
            None
            if m.near_resonance is None
            else {
                "p": m.near_resonance.p,
                "q": m.near_resonance.q,
                "d": m.near_resonance.distance,
            }
        ),
        "independent_residual": m.extras.get("independent_residual"),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", choices=["smoke"], default="smoke")
    ap.add_argument("--ds", type=float, default=5e-3)
    ap.add_argument("--max-steps", type=int, default=200)
    args = ap.parse_args()

    torus = _build_smoke_seed()
    out = ROOT / "data" / f"family_333_qp_{args.seed}.jsonl"
    fh = out.open("w")

    def on_step(m: QPTorusFamilyMember) -> None:
        fh.write(json.dumps(_member_row(m)) + "\n")
        fh.flush()
        ts = dt.datetime.now().isoformat(timespec="seconds")
        print(
            f"[{ts}] member C_J={m.jacobi:.6f} rho={m.rho:.6f} "
            f"ratio={m.freq_ratio:.6f} irr={m.is_practically_irrational} "
            f"res={m.residual_norm:.2e}",
            flush=True,
        )

    fam = continue_qp_family_arclength(
        torus,
        ds=args.ds,
        max_steps=args.max_steps,
        direction="both",
        corrector_tol=1e-8,
        on_step=on_step,
    )
    fh.close()
    print(
        f"DONE: {len(fam.members)} members, {len(fam.folds)} folds, "
        f"{len(fam.resonance_crossings)} resonance crossings, "
        f"terminated={fam.terminated_reason}",
        flush=True,
    )


if __name__ == "__main__":
    main()
