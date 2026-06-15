"""Generate ``data/family_296_3d_em_11.jsonl`` — the 3D Earth-Moon (1,1) family
reproduced by the Phase 2 family tracer (#296).

Walks pseudo-arclength from the #287 spike's converged 3D seed forward and
backward at step 0.01 (4-D arclength), recording every accepted member with
full state vector + period + Jacobi C + monodromy summary + Floquet
multipliers. NOT a catalogue writeback — a family REGISTRY (Phase 3 will
adjudicate catalogue admission via literature_check + V0-V5).

The seed is the spike's converged member at z0_guess = 0.05 from the planar
Braik-Ross C11a IC (data/spike_287.jsonl). The spike's family extent in
``(x0, z0, T, C)`` is the cross-check; this tracer should reach further.

Usage::

  uv run python scripts/family_296_3d_em_11_generator.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_3d_family_tracer import (
    Family3DMember,
    continue_general_3d_family,
)

# Spike seed (#287 converged 3D member; data/spike_287.jsonl).
SPIKE_X0 = -0.8116406668238195
SPIKE_Z0 = -0.2408102083477011
SPIKE_YDOT0 = -0.10629710963669947
SPIKE_T = 10.204301970414399

EM_MU = 1.2150584270572e-2  # Braik-Ross 2026 Table 1
EM_L_KM = 384400.0
EM_T_S = 375699.8

OUT_PATH = Path("/home/bruce/dev/cyclers/data/family_296_3d_em_11.jsonl")


def _member_to_dict(member: Family3DMember, *, system: cr3bp.CR3BPSystem) -> dict:
    """Serialize a Family3DMember to a JSON-ready dict.

    Includes the full state vector, period (TU + days), Jacobi C, monodromy
    summary stats (det + condition number; the full 36-component matrix is
    kept available but flattened), the 6 Floquet multipliers, and the
    stability tag.
    """
    orb = member.orbit
    rec = {
        "step_index": int(member.step_index),
        "arc_length": float(member.arc_length),
        "stability_tag": member.stability_tag,
        # Core orbit.
        "state_nd": orb.state0.tolist(),
        "T_TU": float(orb.T_TU),
        "T_days": float(orb.T_TU * (system.t_s / 86400.0)),
        "jacobi_constant": float(orb.jacobi),
        "degenerate_planar": bool(orb.degenerate_planar),
        # Closure diagnostics.
        "corrector_residual": float(orb.corrector_residual),
        "independent_closure_residual": float(orb.independent_closure_residual),
        "n_iter": int(orb.n_iter),
    }
    if member.monodromy is not None:
        rec["monodromy"] = member.monodromy.tolist()
        rec["monodromy_det"] = float(np.linalg.det(member.monodromy))
    if member.floquet is not None:
        rec["floquet_real"] = [float(z.real) for z in member.floquet]
        rec["floquet_imag"] = [float(z.imag) for z in member.floquet]
        rec["floquet_abs"] = [float(abs(z)) for z in member.floquet]
    return rec


def main() -> None:
    """Run the pseudo-arclength walk and write the family registry."""
    system = cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )
    seed = np.array(
        [SPIKE_X0, 0.0, SPIKE_Z0, 0.0, SPIKE_YDOT0, 0.0],
        dtype=np.float64,
    )

    print(f"[{time.strftime('%H:%M:%S')}] Starting pseudo-arclength walk")
    print(
        f"  seed:  x0={SPIKE_X0:+.5f}  z0={SPIKE_Z0:+.5f}  "
        f"ydot0={SPIKE_YDOT0:+.5f}  T={SPIKE_T:.4f} TU"
    )
    print("  step:  0.01 (arclength in (x0, z0, ydot0, T))")
    print("  n_steps_max: 200 per direction")

    n_step_reports = 0

    def progress(member: Family3DMember) -> None:
        nonlocal n_step_reports
        n_step_reports += 1
        if n_step_reports % 10 == 0:
            ts = time.strftime("%H:%M:%S")
            orb = member.orbit
            print(
                f"  [{ts}] step={member.step_index:+4d}  "
                f"x0={orb.state0[0]:+.5f}  z0={orb.state0[2]:+.5f}  "
                f"T={orb.T_TU:.4f}  C={orb.jacobi:.4f}  "
                f"tag={member.stability_tag}  "
                f"closure={orb.independent_closure_residual:.2e}"
            )

    fam = continue_general_3d_family(
        system,
        seed,
        SPIKE_T,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=200,
        direction="both",
        corrector_tol=1e-10,
        closure_tol=1e-6,
        fold_detection=True,
        monodromy_eval=True,
        on_step=progress,
    )

    print(f"\n[{time.strftime('%H:%M:%S')}] Walk complete.")
    print(
        f"  members: {len(fam.members)} "
        f"(forward {fam.n_steps_forward}, backward {fam.n_steps_backward})"
    )
    print(f"  forward terminated: {fam.forward_termination}")
    print(f"  backward terminated: {fam.backward_termination}")
    print(f"  folds detected: {len(fam.folds)}")

    closures = [m.orbit.independent_closure_residual for m in fam.members]
    print(f"  closure max: {max(closures):.2e}, median: {np.median(closures):.2e}")
    x0s = np.array([m.orbit.state0[0] for m in fam.members])
    z0s = np.array([m.orbit.state0[2] for m in fam.members])
    ts = np.array([m.orbit.T_TU for m in fam.members])
    cs = np.array([m.orbit.jacobi for m in fam.members])
    print(f"  x0 extent: {x0s.min():+.5f} .. {x0s.max():+.5f}")
    print(f"  z0 extent: {z0s.min():+.5f} .. {z0s.max():+.5f}")
    print(f"  T extent:  {ts.min():.4f} .. {ts.max():.4f} TU")
    print(f"  C extent:  {cs.min():.4f} .. {cs.max():.4f}")

    # Stability tag census.
    from collections import Counter

    tag_counts = Counter(m.stability_tag for m in fam.members)
    print(f"  stability tags: {dict(tag_counts)}")

    # Write family JSONL: header + each member.
    print(f"\n[{time.strftime('%H:%M:%S')}] Writing {len(fam.members)} members to {OUT_PATH}")
    with OUT_PATH.open("w") as f:
        # Header record (single line; identifies the walk).
        header = {
            "type": "header",
            "issue": 296,
            "phase": "phase2_3d_family_tracer",
            "system": {
                "primary": "earth",
                "secondary": "moon",
                "mu": EM_MU,
                "l_km": EM_L_KM,
                "t_s": EM_T_S,
            },
            "seed": {"x0": SPIKE_X0, "z0": SPIKE_Z0, "ydot0": SPIKE_YDOT0, "T_TU": SPIKE_T},
            "continuation": fam.continuation_mode,
            "step": fam.step,
            "n_members": len(fam.members),
            "n_folds": len(fam.folds),
            "forward_termination": fam.forward_termination,
            "backward_termination": fam.backward_termination,
            "x0_extent": [float(x0s.min()), float(x0s.max())],
            "z0_extent": [float(z0s.min()), float(z0s.max())],
            "T_TU_extent": [float(ts.min()), float(ts.max())],
            "C_extent": [float(cs.min()), float(cs.max())],
            "closure_max": float(max(closures)),
            "closure_median": float(np.median(closures)),
            "stability_tags": dict(tag_counts),
            "metadata": fam.metadata,
            "discipline": (
                "Likely-rediscovery of Antoniadou-Voyatzis 2018 spatial CR3BP "
                "(KNOWN_CORPUS commit 568d8a4); Phase 3 will run literature_check."
            ),
        }
        f.write(json.dumps(header) + "\n")
        for m in fam.members:
            rec = _member_to_dict(m, system=system)
            f.write(json.dumps(rec) + "\n")

    # Brief fold summary.
    if fam.folds:
        print(f"\nFolds detected ({len(fam.folds)}):")
        for fold in fam.folds:
            print(
                f"  step={fold.step_index:+4d}  param={fold.natural_param}  "
                f"tangent: {fold.tangent_before:+.4f} -> {fold.tangent_after:+.4f}"
            )

    print(f"\n[{time.strftime('%H:%M:%S')}] Done.")


if __name__ == "__main__":
    main()
