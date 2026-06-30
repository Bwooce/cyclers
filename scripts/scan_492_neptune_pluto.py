"""#492 — novel moon-cycler discovery sweep: Neptune (Triton-flyby) + Pluto (Charon-flyby).

The genuinely-novel moon-cycler frontier (vs #491's V0-known ingestion). #489 established the
binding constraint: the flyby (repeated) body must be MASSIVE enough to turn the trajectory.
The UNSWEPT viable systems (Uranus was done in #285/#312) with a massive flyby body + NO
published cycler literature are Neptune (Triton, mu=1428) and Pluto (Charon, mu=106). This
runs the #285 prioritized repeated-moon scan on both, then we apply the #489 physical-sanity
filter (massive flyby) + a live lit-novelty check on any survivor.

NO catalogue writeback — SILVERs feed the gauntlet + human review (as #312 did for Uranus).
"""

from __future__ import annotations

from pathlib import Path

from cyclerfinder.search.saturn_uranus_campaign import run_prioritized_scan

ROOT = Path(__file__).resolve().parent.parent
N_REV_GRID = (0, 1, 2, 3, 4, 5)
SEQ_LENGTHS = (3,)  # two-body repeated cyclers X-Y-X


def main() -> int:
    systems = [
        ("Neptune", ["Triton", "Proteus"], ROOT / "data" / "scan_492_neptune.jsonl"),
        ("Pluto", ["Charon", "Nix", "Hydra"], ROOT / "data" / "scan_492_pluto.jsonl"),
    ]
    for primary, moons, out in systems:
        print(f"\n=== #492 scan: {primary} {moons} ===", flush=True)
        summary = run_prioritized_scan(
            primary=primary,
            moons=moons,
            seq_lengths=SEQ_LENGTHS,
            n_rev_grid=N_REV_GRID,
            out_path=out,
            git_sha="492-neptune-pluto",
        )
        d = summary.as_dict()["summary"] if "summary" in summary.as_dict() else summary.as_dict()
        print(
            f"  enumerated={d.get('enumerated')} closed={d.get('closed')} "
            f"SILVER={d.get('verdict_counts', {}).get('SILVER', d.get('silver_pre_guards'))} "
            f"near_miss={d.get('near_miss_count')}",
            flush=True,
        )
        top = d.get("top5_near_misses", [])
        for nm in top[:5]:
            print(
                f"    {'-'.join(nm['sequence'])} n_rev={nm['n_rev']} "
                f"res={nm['residual_kms']:.4f} maxVinf={nm['max_vinf_kms']:.2f}",
                flush=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
