"""#388 — family-pinned penalty homotopy closer batch (detached, checkpointed).

For each SnLm descriptor row x best-phase epoch, ramps the V∞-anchor penalty from
a calibrated W down to 0 (``search.family_pinned_shoot.family_pinned_shoot``) and
records the UNPENALIZED (λv=0) verdict: converged?, defect, emerged V∞ vs the
SOURCED anchors, anchor-match, anchor-retention, bend. HELD — no writeback. A row
is flagged PROPOSED V0->V1 only for ``mcconaghy-2006-em-k2`` AND only if the λv=0
solve converges, anchor-matches within 0.5 km/s of BOTH anchors, and is
bend-feasible — recorded, never applied.

The penalty weight ladder top ``W`` is a TUNING value (not sourced): calibrated so
the penalty bites without swamping continuity. See the design spec.

Launch detached:
    setsid nohup uv run python scripts/shooter_family_pinned_batch.py --resume \
        > data/runs/shooter-family-pinned.log 2>&1 &
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.family_pinned_shoot import family_pinned_shoot
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import (
    candidate_epochs,
    russell_parent_to_ballistic_seed,
)
from cyclerfinder.search.shooter_russell_seed import russell_shooting_seed

K: int = 3
LAUNCH_WINDOW_SYNODICS: range = range(1, 22)
EPOCH_GRID: int = 100

# Penalty weight ladder top (TUNING value, not sourced — see design spec §Open
# questions). Calibrated so the penalty residual is comparable to the converged
# continuity residual on row 9.353Gg2.
WEIGHT_LADDER: tuple[float, ...] = (40.0, 10.0, 2.5, 0.5, 0.0)
SHOOT_ACCURACY: float = 1e-9
SHOOT_MAX_NFEV: int = 100
LEG_WALL_BUDGET_SEC: float = 30.0
ANCHOR_MATCH_TOL_KMS: float = 0.5

RUNLOG = Path("data/runs/shooter-family-pinned.jsonl")


def _load_done(runlog: Path) -> set[tuple[str, int]]:
    done: set[tuple[str, int]] = set()
    if not runlog.exists():
        return done
    with runlog.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "id" in rec and "epoch_index" in rec:
                done.add((str(rec["id"]), int(rec["epoch_index"])))
    return done


def _append(runlog: Path, rec: dict[str, Any]) -> None:
    runlog.parent.mkdir(parents=True, exist_ok=True)
    with runlog.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", action="store_true", help="skip (row,epoch) pairs in runlog")
    args = parser.parse_args()
    sys.stdout.reconfigure(line_buffering=True)

    m = RussellModel()
    ephem = Ephemeris("astropy")
    done = _load_done(RUNLOG) if args.resume else set()
    if args.resume:
        print(f"[resume] {len(done)} (row,epoch) pairs already in {RUNLOG}")

    t_wall0 = time.monotonic()
    catalog = load_catalog()

    for e in catalog.entries:
        wall = time.monotonic() - t_wall0
        try:
            phsi = descriptor_to_phsi(e.raw)
            if phsi is None:
                continue
            cyc = assemble_cycler(m, phsi)
            if cyc is None:
                continue
            try:
                seed = russell_parent_to_ballistic_seed(m, cyc, e.raw)
            except ValueError:
                continue

            vlevel = str(e.raw.get("validation_level", "V0"))
            bodies = tuple(dict.fromkeys(seed.sequence))
            anchors = {"E": seed.vinf_anchor_e_kms, "M": seed.vinf_anchor_m_kms}
            epochs = candidate_epochs(
                ephem, 0.0, launch_window_synodics=LAUNCH_WINDOW_SYNODICS, grid=EPOCH_GRID
            )[:K]

            for epoch_index, t0 in enumerate(epochs):
                if (e.id, epoch_index) in done:
                    continue
                wall = time.monotonic() - t_wall0
                t_shoot0 = time.monotonic()

                def _hb(
                    kind: str,
                    count: int,
                    defect_norm: float,
                    elapsed: float,
                    _id: str = e.id,
                    _ep: int = epoch_index,
                ) -> None:
                    now = time.monotonic() - t_wall0
                    print(
                        f"[{now:.0f}s]   {_id:24s} ep{_ep} {kind}{count} "
                        f"defect={defect_norm:.3e} ({elapsed:.0f}s)"
                    )

                row_error: str | None = None
                res = None
                try:
                    sseed = russell_shooting_seed(seed, t0_sec=t0, ephem=ephem)
                    res = family_pinned_shoot(
                        sseed,
                        ephem=ephem,
                        bodies=bodies,
                        vinf_anchors=anchors,
                        weight_ladder=WEIGHT_LADDER,
                        accuracy=SHOOT_ACCURACY,
                        max_nfev=SHOOT_MAX_NFEV,
                        max_wall_sec=LEG_WALL_BUDGET_SEC,
                        progress=_hb,
                    )
                except Exception as exc:  # honest per-(row,epoch) record, never raised
                    row_error = f"{type(exc).__name__}: {exc}"

                shoot_wall = time.monotonic() - t_shoot0
                if res is None:
                    _append(
                        RUNLOG,
                        {
                            "id": e.id,
                            "validation_level": vlevel,
                            "sequence": list(seed.sequence),
                            "epoch_index": epoch_index,
                            "epoch_sec": t0,
                            "shot": False,
                            "error": row_error,
                            "shoot_wall_sec": shoot_wall,
                            "anchor_e_kms": seed.vinf_anchor_e_kms,
                            "anchor_m_kms": seed.vinf_anchor_m_kms,
                        },
                    )
                    print(
                        f"[{wall:.0f}s] {e.id:24s} [{vlevel}] ep{epoch_index} "
                        f"NO SHOOT ({row_error})"
                    )
                    continue

                final = res.final
                vinf = list(final.vinf_per_encounter_kms)
                best_e = min((abs(v - seed.vinf_anchor_e_kms) for v in vinf), default=float("inf"))
                best_m = min((abs(v - seed.vinf_anchor_m_kms) for v in vinf), default=float("inf"))
                anchor_match = best_e <= ANCHOR_MATCH_TOL_KMS and best_m <= ANCHOR_MATCH_TOL_KMS
                converged = final.converged
                bend = final.bend_feasible
                promote = e.id == "mcconaghy-2006-em-k2" and converged and anchor_match and bend

                _append(
                    RUNLOG,
                    {
                        "id": e.id,
                        "validation_level": vlevel,
                        "sequence": list(seed.sequence),
                        "epoch_index": epoch_index,
                        "epoch_sec": t0,
                        "shot": True,
                        "final_weight": res.final_weight,
                        "converged": converged,
                        "defect_norm": final.defect_norm,
                        "seed_defect_norm": final.seed_defect_norm,
                        "vinf_per_encounter_kms": vinf,
                        "anchor_e_kms": seed.vinf_anchor_e_kms,
                        "anchor_m_kms": seed.vinf_anchor_m_kms,
                        "best_e_residual_kms": best_e,
                        "best_m_residual_kms": best_m,
                        "anchor_match": anchor_match,
                        "anchor_retention_kms": res.anchor_retention_kms,
                        "bend_feasible": bend,
                        "trace": res.trace,
                        "promote_proposed_held": promote,
                        "shoot_wall_sec": shoot_wall,
                        "error": row_error,
                    },
                )
                flag = " *** PROPOSED V0->V1 (HELD) ***" if promote else ""
                print(
                    f"[{wall:.0f}s] {e.id:24s} [{vlevel}] ep{epoch_index} "
                    f"conv={converged} defect={final.defect_norm:.3e} "
                    f"vinf={[round(v, 2) for v in vinf]} "
                    f"anchorE/M={best_e:.2f}/{best_m:.2f} retain={res.anchor_retention_kms:.2f} "
                    f"match={anchor_match} bend={bend} ({shoot_wall:.0f}s){flag}"
                )
        except Exception as exc:  # one bad row must not abort the batch
            _append(
                RUNLOG,
                {"id": e.id, "shot": False, "error": f"row-fatal {type(exc).__name__}: {exc}"},
            )
            print(f"[{wall:.0f}s] {e.id:24s} ROW-FATAL {type(exc).__name__}: {exc}")

    print()
    print("=" * 72)
    print(f"Batch complete. Runlog: {RUNLOG} (K={K}, weight_ladder={WEIGHT_LADDER})")
    print("HELD — no writeback (data/catalogue.yaml and validate.py untouched).")


if __name__ == "__main__":
    main()
