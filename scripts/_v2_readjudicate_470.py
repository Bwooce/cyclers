"""#473 STEP 0 — re-adjudicate the 10 #470 tours by DRIFT-SERIES SHAPE.

#470 ran the official run_v2_moontour on the #468 in-band tours and got 0/10,
but it ONLY compared max_drift to the strict 50,000 km floor — it never tested
whether a tour's drift is BOUNDED-OSCILLATING (a valid quasi_cycler) vs
MONOTONICALLY DIVERGENT (a genuine reject). This is the corrected adjudication:
re-run each #470 winning config through the OFFICIAL run_v2_moontour at
n_cycles=10 and classify EACH by drift-series shape.

A tour that bounds-and-returns is a quasi_cycler #470 WRONGLY rejected; one that
diverges monotonically stays a genuine reject.

Output: out/readjudicate_470_runlog.jsonl (one flushed line per tour) +
out/readjudicate_470_summary.json
"""

from __future__ import annotations

import importlib.util as ilu
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

from cyclerfinder.data.validation.v2_moontour import run_v2_moontour
from cyclerfinder.search.releg_moontour import close_powered_cycle
from cyclerfinder.search.releg_solver import MultiRevLeveragingReleg

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _drift_shape_473 import BOUNDED_SHAPES, classify_drift_shape  # noqa: E402

_spec = ilu.spec_from_file_location(
    "campaign_468", str(REPO / "scripts" / "campaign_468_multirev_tour.py")
)
C = ilu.module_from_spec(_spec)
sys.modules["campaign_468"] = C
_spec.loader.exec_module(C)

PROPOSALS = REPO / "data" / "admission_proposals_468.jsonl"
RUNLOG = REPO / "out" / "readjudicate_470_runlog.jsonl"
SUMMARY = REPO / "out" / "readjudicate_470_summary.json"
N_CYCLES = 10


def _now() -> str:
    return datetime.now(UTC).astimezone().isoformat(timespec="seconds")


def winning_config(sk):
    """Re-run the #468 campaign sweep; return (tofs, off, verdict) of cheapest in-band."""
    n_legs = len(sk.sequence) - 1
    best = None
    best_tofs = None
    best_off = None
    for scale in C.TOF_SCALES:
        for off in C.PHASE_SEEDS:
            tofs = C._geomean_tofs(sk.sequence, scale=scale)
            v = close_powered_cycle(
                primary=sk.system,
                sequence=sk.sequence,
                leg_tofs_days=tofs,
                n_revs=tuple([0] * n_legs),
                releg=MultiRevLeveragingReleg(),
                phasing=C._phasing(sk.sequence, off),
                dv_band="powered_dsm",
            )
            if v.prefilter_skipped:
                return None, None, None
            if (
                v.feasible
                and v.total_dv_kms < C.IN_BAND_KMS
                and (best is None or v.total_dv_kms < best.total_dv_kms)
            ):
                best = v
                best_tofs = tofs
                best_off = off
    return best_tofs, best_off, best


def main() -> None:
    target_seqs = {
        tuple(json.loads(line)["sequence"]) for line in PROPOSALS.read_text().splitlines()
    }
    RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    runlog = RUNLOG.open("w", encoding="utf-8")
    rows = []
    t0 = time.time()
    item_id = 0
    n_target = len(target_seqs)

    for sk in C.SKELETONS:
        if sk.sequence not in target_seqs:
            continue
        item_id += 1
        tofs, off, wv = winning_config(sk)
        if tofs is None:
            row = {"sequence": list(sk.sequence), "label": sk.label, "status": "no-inband-config"}
            rows.append(row)
            runlog.write(json.dumps(row, ensure_ascii=True) + "\n")
            runlog.flush()
            continue
        n_legs = len(sk.sequence) - 1
        # The #468 campaign phasing is offset*index per distinct moon (radians).
        # run_v2_moontour uses (phase0_deg, rel_offset_deg). For the 2-distinct-moon
        # tours the campaign off maps to rel_offset = off (rad) on the 2nd moon, so
        # we re-express. For 3-distinct-moon tours the campaign spreads moons by
        # off*index; run_v2_moontour's convention spreads the 3rd moon by
        # rel_off + 2pi*(j-1)/n. The phasings differ for 3-moon tours, so we drive
        # BOTH the official gauntlet (its own convention) AND report.
        # We pass phase0_deg=0 and rel_offset_deg from the campaign's first nonzero
        # spacing (off radians -> degrees) so the official run is phased like #470's
        # winning seed as closely as the official convention allows.
        import math

        rel_off_deg = math.degrees(off) if off else 0.0
        vinf = (wv.target_vinf_kms,) * len(sk.sequence)
        verdict = run_v2_moontour(
            candidate_id="-".join(sk.sequence),
            sequence=sk.sequence,
            vinf_tuple_kms=vinf,
            leg_tofs_days=tofs,
            rel_offset_deg=rel_off_deg,
            system=None,
            n_cycles=N_CYCLES,
            n_revs=tuple([0] * n_legs),
            phase0_deg=0.0,
            releg=MultiRevLeveragingReleg(),
            dv_band="powered_dsm",
        )
        drifts = [float(c.rendezvous_drift_kms) for c in verdict.per_cycle]
        cls = classify_drift_shape(drifts[1:])
        flipped = (
            verdict.n_cycles_completed >= N_CYCLES
            and cls["shape"] in BOUNDED_SHAPES
            and not verdict.passes_v2  # strict pass would be a STRICT cycler, not a flip
        )
        strict_cycler = verdict.passes_v2
        row = {
            "sequence": list(sk.sequence),
            "label": sk.label,
            "leg_tofs_days": [round(t, 3) for t in tofs],
            "campaign_off_rad": off,
            "target_vinf_kms": round(wv.target_vinf_kms, 3),
            "n_cycles_completed": verdict.n_cycles_completed,
            "passes_v2_strict": verdict.passes_v2,
            "max_drift_kms": round(verdict.max_drift_kms, 1),
            "per_cycle_drift_kms": [round(d, 1) for d in drifts],
            "shape": cls["shape"],
            "shape_detail": cls,
            "flips_to_bounded_quasi": flipped,
            "is_strict_cycler": strict_cycler,
            "verdict_470": "reject (max_drift > 50k floor)",
        }
        rows.append(row)
        runlog.write(
            json.dumps(
                {
                    "item_id": item_id,
                    "sub_step": sk.label,
                    "result": (
                        "strict_cycler"
                        if strict_cycler
                        else ("flip_bounded" if flipped else "divergent_reject")
                    ),
                    "residual": round(verdict.max_drift_kms, 1),
                    "ts": _now(),
                    "k_of_N": f"{item_id}/{n_target}",
                    "elapsed_s": round(time.time() - t0, 1),
                    "shape": cls["shape"],
                    "per_cycle_drift_kms": [round(d, 1) for d in drifts],
                },
                ensure_ascii=True,
            )
            + "\n"
        )
        runlog.flush()
        tag = (
            "FLIP->bounded"
            if flipped
            else ("STRICT-CYCLER" if strict_cycler else "divergent-reject")
        )
        print(
            f"[{_now()}] {sk.label}: shape={cls['shape']} "
            f"max_drift={verdict.max_drift_kms:,.0f} km {tag}"
        )
    runlog.close()

    n_flip = sum(1 for r in rows if r.get("flips_to_bounded_quasi"))
    n_strict = sum(1 for r in rows if r.get("is_strict_cycler"))
    summary = {
        "run_date": "2026-06-26",
        "n_cycles": N_CYCLES,
        "n_tours": len([r for r in rows if "shape" in r]),
        "n_flips_to_bounded_quasi": n_flip,
        "n_strict_cyclers": n_strict,
        "rows": rows,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\n[{_now()}] FLIPS to bounded quasi_cycler: {n_flip}; strict cyclers: {n_strict}")
    print(f"summary -> {SUMMARY}")


if __name__ == "__main__":
    main()
