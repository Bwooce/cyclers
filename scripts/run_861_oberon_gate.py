"""#861: Uranus-Oberon fold-turning seeding-fix gate driver.

`#860` diagnosed why `#859`'s Stage A smoke test failed on Uranus-Oberon (the
naive `two_body_resonant_seed(x0_sign=-1)` + `continue_family` pipeline lands
on a stable, near-unit-circle family at every published ratio, never the
paper's own labeled unstable saddle family) and proposed two fixes: a
conjugate-apse seed (`jovian_resonant_families.two_body_conjugate_apse_seed`)
and pseudo-arclength fold-turning continuation in place of natural-parameter
continuation (`resonant_atlas_stage_a_prime.fold_turn_family`, wired in `#861`
part 1). This script runs BOTH fixes together on Uranus-Oberon ONLY, at its
six AAS 24-288 (`#728` digest) published resonant ratios, and reports the
gate: does the fold-turned family recover an unstable segment whose C-range
overlaps the paper's own printed range and whose topology (period, winding,
close approach) matches the label, for >=4 of 6 families?

Sourced constants (Anderson & Kumar 2024, AAS 24-288, `#728` digest, verbatim):

    mu = 3.54326e-5 (paper's own stated Uranus-Oberon mass ratio, matching
    this project's own DE440-registry value -- GM_Oberon/(GM_Uranus+GM_Oberon)
    = 205.3/(5.7945564e6+205.3) -- to <0.02% relative, so no "paper vs
    registry" divergence to flag here, unlike the Jupiter-Europa/
    Neptune-Triton sibling modules).

    Jacobi-constant ranges per family (p.7-8, "Unstable Oberon Resonant
    Periodic Orbits"):
        3:4 in [2.9916, 3.0261]   4:3 in [2.9836, 3.0279]
        4:5 in [2.9914, 3.0157]   5:4 in [2.9902, 3.0165]
        5:6 in [2.9921, 3.0104]   6:5 in [2.9949, 3.0109]

Foreground only, chunked (this project's own standing lesson,
`feedback_subagent_background_is_fatal`) -- run repeatedly (same output file)
to resume; already-computed (p, q, seed_kind) cells are skipped.

Output: ``data/found/861_resonant_seeding_oberon_gate/results.jsonl`` (one
JSON line per (p, q, seed_kind) cell, the checkpoint itself).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.resonant_atlas_stage_a_prime as sap

#: AAS 24-288 p.169-equivalent (`#728` digest Sec. 1, "L1/L2 libration-point
#: Jacobi constants at this mu = 3.54326e-5"), the paper's own stated mass
#: ratio -- verbatim.
OBERON_MU = 3.54326e-5

#: AAS 24-288 p.7-8 (`#728` digest, verbatim), the six published families'
#: own printed Jacobi-constant ranges.
OBERON_PUBLISHED_C_RANGES: dict[tuple[int, int], tuple[float, float]] = {
    (3, 4): (2.9916, 3.0261),
    (4, 5): (2.9914, 3.0157),
    (5, 6): (2.9921, 3.0104),
    (4, 3): (2.9836, 3.0279),
    (5, 4): (2.9902, 3.0165),
    (6, 5): (2.9949, 3.0109),
}

RATIOS: list[tuple[int, int]] = list(OBERON_PUBLISHED_C_RANGES)
SEED_KINDS: list[sap.SeedKind] = ["opposition", "conjugate_apse"]


def _oberon_system() -> cr3bp.CR3BPSystem:
    base = cr3bp.cr3bp_system("Uranus", "Oberon")
    return cr3bp.CR3BPSystem(
        mu=OBERON_MU, primary="Uranus", secondary="Oberon", l_km=base.l_km, t_s=base.t_s
    )


def _load_done(path: Path) -> dict[tuple[int, int, str], dict[str, Any]]:
    done: dict[tuple[int, int, str], dict[str, Any]] = {}
    if not path.exists():
        return done
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            done[(row["p"], row["q"], row["seed_kind"])] = row
    return done


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="data/found/861_resonant_seeding_oberon_gate/results.jsonl")
    ap.add_argument("--max-cells", type=int, default=1000, help="cap cells run THIS invocation")
    ap.add_argument("--max-steps", type=int, default=60)
    ap.add_argument("--ds0", type=float, default=0.012)
    ap.add_argument("--ds-max", type=float, default=0.04)
    ap.add_argument("--c-span", type=float, default=0.12)
    ap.add_argument("--record-every", type=int, default=2)
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    done = _load_done(out_path)
    system = _oberon_system()

    n_run = 0
    for p, q in RATIOS:
        for seed_kind in SEED_KINDS:
            key = (p, q, seed_kind)
            if key in done:
                continue
            if n_run >= args.max_cells:
                print(f"[861] hit --max-cells={args.max_cells}, stopping this invocation")
                return
            t0 = time.time()
            print(f"[861] running {p}:{q} {seed_kind} ...", flush=True)
            result = sap.fold_turn_family(
                system,
                p,
                q,
                seed_kind,
                system_key="uranus-oberon",
                max_steps=args.max_steps,
                ds0=args.ds0,
                ds_max=args.ds_max,
                c_span=args.c_span,
                record_every=args.record_every,
            )
            elapsed = time.time() - t0
            c_lo, c_hi = OBERON_PUBLISHED_C_RANGES[(p, q)]
            members_payload = [sap.fold_turn_member_to_dict(m) for m in result.members]
            unstable_in_range = [
                m for m in result.members if not m.stable and c_lo <= m.jacobi <= c_hi
            ]
            topology_matches = [
                m for m in unstable_in_range if m.period_matches_q and m.winding_matches_p
            ]
            row = {
                "p": p,
                "q": q,
                "seed_kind": seed_kind,
                "elapsed_s": elapsed,
                "seed_converged": result.seed_converged,
                "seed_jacobi": result.seed_jacobi,
                "half_crossings": result.half_crossings,
                "ydot0_sign": result.ydot0_sign,
                "n_members": len(result.members),
                "stop_reason_up": result.stop_reason_up,
                "stop_reason_down": result.stop_reason_down,
                "c_min": result.c_min,
                "c_max": result.c_max,
                "published_c_lo": c_lo,
                "published_c_hi": c_hi,
                "n_unstable": result.n_unstable,
                "n_unstable_in_published_range": len(unstable_in_range),
                "n_topology_matches": len(topology_matches),
                "max_abs_lambda": result.max_abs_lambda,
                "best_topology_match": (
                    sap.fold_turn_member_to_dict(topology_matches[0]) if topology_matches else None
                ),
                "members": members_payload,
            }
            with out_path.open("a") as f:
                f.write(json.dumps(row) + "\n")
            n_run += 1
            print(
                f"[861]   done in {elapsed:.1f}s: n_members={row['n_members']} "
                f"n_unstable={row['n_unstable']} n_in_range={row['n_unstable_in_published_range']} "
                f"n_topo_match={row['n_topology_matches']} max|lambda|={row['max_abs_lambda']:.3g} "
                f"C=[{row['c_min']:.4f},{row['c_max']:.4f}] stop=({row['stop_reason_up']},"
                f"{row['stop_reason_down']})",
                flush=True,
            )

    print(f"[861] all {len(RATIOS) * len(SEED_KINDS)} cells done -> {out_path}")


if __name__ == "__main__":
    main()
