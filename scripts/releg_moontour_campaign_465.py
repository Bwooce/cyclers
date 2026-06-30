"""#465 at-scale powered moon-tour discovery campaign.

The explicit follow-on to the multi-rev leveraging releg verdict
(docs/superpowers/plans/2026-06-26-465-multirev-leveraging-verdict.md): sweep moon-tour
skeletons x ToF-scale x phasing across the Galilean + Saturnian systems with the
MultiRevLeveragingReleg backend, recording which cycles CLOSE IN-BAND (total ΔV below the
3.5 km/s/cycle powered ceiling). Respects the reachability ceiling (V∞ ≲ ~2·V_M) honestly —
infeasible skeletons stay empty (prefilter or chain-stall), never a fabricated bridge.

Output: data/releg_moontour_campaign_465.jsonl (one row per skeleton x scale x phasing) +
a console summary of distinct in-band closures. Distinct in-band hits are then lit-checked
(separate step) to separate novel candidates from V0-known reproductions of the published
Campagnola/Strange endgame tours. No catalogue row is self-admitted.
"""

from __future__ import annotations

import json
import math
import time
from itertools import combinations

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.releg_moontour import close_powered_cycle
from cyclerfinder.search.releg_solver import MultiRevLeveragingReleg

OUT = "data/releg_moontour_campaign_465.jsonl"
IN_BAND_KMS = 3.5  # powered ceiling (#465 gate)

# Moons per system, inner->outer (adjacent pairs share Tisserand contours; the chain
# walks V∞ WITHIN a contour, so only contiguous-ring tours are bridgeable).
SYSTEMS: dict[str, list[str]] = {
    "Jupiter": ["Io", "Europa", "Ganymede", "Callisto"],
    "Saturn": ["Enceladus", "Tethys", "Dione", "Rhea", "Titan"],
}
TOF_SCALES = (1.0, 1.2, 1.5, 2.0, 2.5)
PHASE_STEPS = (0.5, 1.0, 1.5, 2.0)


def _period_days(m: str) -> float:
    s = SATELLITES[m]
    return 2.0 * math.pi * math.sqrt(s.sma_km**3 / PRIMARIES[s.primary]) / 86400.0


def geomean_tofs(sequence: tuple[str, ...], scale: float) -> tuple[float, ...]:
    return tuple(
        scale * math.sqrt(_period_days(sequence[k]) * _period_days(sequence[k + 1]))
        for k in range(len(sequence) - 1)
    )


def skeletons(moons: list[str]) -> list[tuple[str, ...]]:
    """Closed cyclic tours over 3- and 4-moon CONTIGUOUS sub-rings (adjacent pairs)."""
    out: list[tuple[str, ...]] = []
    idx = {m: i for i, m in enumerate(moons)}
    # 3-moon contiguous tours A-B-C-A and 4-moon A-B-C-D-A over adjacent rings
    for size in (3, 4):
        for combo in combinations(moons, size):
            # contiguous in the inner->outer ring?
            ii = sorted(idx[m] for m in combo)
            if ii[-1] - ii[0] != size - 1:
                continue
            ring = [moons[i] for i in ii]
            out.append((*ring, ring[0]))  # close the cycle
    return out


def main() -> int:
    releg = MultiRevLeveragingReleg()
    rows: list[dict] = []
    t0 = time.monotonic()
    n = 0
    for primary, moons in SYSTEMS.items():
        for seq in skeletons(moons):
            for scale in TOF_SCALES:
                tofs = geomean_tofs(seq, scale)
                for step in PHASE_STEPS:
                    phasing = {m: i * step for i, m in enumerate(dict.fromkeys(seq))}
                    n += 1
                    try:
                        v = close_powered_cycle(
                            primary=primary,
                            sequence=seq,
                            leg_tofs_days=tofs,
                            n_revs=tuple(0 for _ in range(len(seq) - 1)),
                            releg=releg,
                            phasing=phasing,
                            dv_band="powered_dsm",
                        )
                    except Exception as exc:
                        rows.append(
                            {
                                "primary": primary,
                                "sequence": list(seq),
                                "scale": scale,
                                "phase_step": step,
                                "error": repr(exc),
                            }
                        )
                        continue
                    rows.append(
                        {
                            "primary": primary,
                            "sequence": list(seq),
                            "scale": scale,
                            "phase_step": step,
                            "feasible": v.feasible,
                            "prefilter_skipped": v.prefilter_skipped,
                            "total_dv_kms": v.total_dv_kms
                            if math.isfinite(v.total_dv_kms)
                            else None,
                            "continuity_kms": v.continuity_residual_kms
                            if math.isfinite(v.continuity_residual_kms)
                            else None,
                            "in_band": bool(v.feasible and v.total_dv_kms < IN_BAND_KMS),
                            "dv_band": v.dv_band,
                        }
                    )
    with open(OUT, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    dt = time.monotonic() - t0
    print(f"[campaign] {n} skeleton x scale x phasing combos in {dt:.1f}s -> {OUT}")

    # Distinct in-band closures (best ΔV per skeleton).
    inband = [r for r in rows if r.get("in_band")]
    best: dict[tuple, dict] = {}
    for r in inband:
        key = (r["primary"], tuple(r["sequence"]))
        if key not in best or r["total_dv_kms"] < best[key]["total_dv_kms"]:
            best[key] = r
    print(
        f"\n=== IN-BAND powered moon-tour closures: {len(inband)} runs, "
        f"{len(best)} distinct skeletons ==="
    )
    for (primary, seq), r in sorted(best.items(), key=lambda kv: kv[1]["total_dv_kms"]):
        print(
            f"  {primary:8} {'-'.join(seq):40} ΔV={r['total_dv_kms']:.3f} km/s "
            f"(scale {r['scale']}, phase {r['phase_step']}, band {r['dv_band']})"
        )
    oob = [r for r in rows if r.get("feasible") and not r.get("in_band")]
    empties = [r for r in rows if r.get("prefilter_skipped")]
    print(f"\nfeasible-but-oob runs: {len(oob)}; prefilter-empty runs: {len(empties)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
