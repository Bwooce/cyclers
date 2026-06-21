"""#388 — the decisive seeding test: full n-body shooter from the Russell parent.

Feeds the CONSTRUCTED, golden-validated Russell ψ generic-return parent (the
seed bridge ``search/shooter_russell_seed.russell_shooting_seed``) into the
existing full-state multiple-shooting corrector
(``nbody/shooter.shoot`` — the Jones/SNOPT analogue, full restricted-n-body
propagation) on the real (DE440 via astropy) ephemeris, for every
descriptor-bearing catalogue row.

This is the strongest remaining shot at breaking the #135 / S1L1
seeding/family-selection wall: it is the one lever never tried — seeding the
full shooter from the literal constructed parent rather than from the blind /
near-miss-survey basin. See ``docs/notes/2026-06-21-narc-continuation-results.md``
(post-run analysis) for the rationale.

**STM Jacobian (2026-06-21).** The earlier FD run was compute-infeasible: the
LM solver's internal finite-difference Jacobian is (6*n_nodes+1) full-state
residuals per iteration (~1000 s/iteration serially), which forced a K=1 /
max_nfev=2 coverage cap that could only see the seed basin, never convergence.
This batch uses the analytic block-bidiagonal STM Jacobian
(``shooter.shoot(jacobian="stm")``) — ONE co-integrated variational propagation
per leg instead of the 6*n_nodes+1 FD re-propagations — so the corrector can run
to convergence. The coverage cap is therefore removed: the methodologically
-correct K best-phase-error epochs over LaunchWindow 1..21 are shot to
completion, however long that takes (per ``feedback_long_runs_acceptable`` —
monitorability is not the measure of the best path).

**Detached + checkpointed.** Each (row, epoch) record is written to the JSONL
runlog immediately on completion (append + flush), so a kill/restart loses
nothing. ``--resume`` skips (row, epoch) pairs already present in the runlog.
Launch detached so the run survives agent reaping:

    setsid nohup uv run python scripts/shooter_russell_batch.py --resume \
        > data/runs/shooter-stm-batch.log 2>&1 &

Per the standing discipline: HELD — no writeback. This script never touches
``data/catalogue.yaml`` or ``validate.py``; it only writes a JSONL runlog and
prints a summary. A row is flagged ``PROPOSED V0->V1 (HELD)`` only for
``mcconaghy-2006-em-k2`` AND only if it converges, anchor-matches, and is
bend-feasible — and even then it is recorded, never applied.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.nbody import shooter
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import (
    candidate_epochs,
    russell_parent_to_ballistic_seed,
)
from cyclerfinder.search.shooter_russell_seed import russell_shooting_seed

# --- Epoch scan (methodologically correct; NO monitorability cap) -------------
# The best-K phase-error epochs over a LaunchWindow 1..21 synodic scan are shot to
# completion. With the STM Jacobian each shoot is tractable, so the earlier
# FD-era K=1 / range(1, 6) coverage cap is removed (see module docstring).
K: int = 3
LAUNCH_WINDOW_SYNODICS: range = range(1, 22)
EPOCH_GRID: int = 100

# --- shoot() budget -----------------------------------------------------------
# accuracy mirrors tests/nbody/test_shooter_jones_gate.py. With the analytic STM
# Jacobian the per-iteration cost is one residual (~40 s) + one variational
# propagation per leg, not the (6*n_nodes+1) FD residuals — so the corrector can
# run to a real fixed point. max_nfev is generous; the LM solve stops on its own
# ftol/xtol when it converges or stalls.
SHOOT_ACCURACY: float = 1e-9
SHOOT_MAX_NFEV: int = 100
LEG_WALL_BUDGET_SEC: float = 30.0
SHOOT_JACOBIAN: str = "stm"

# The anchor-match tolerance: an emerged per-encounter V∞ must land within this
# of the SOURCED E and M anchors for the row to count as anchor-matched.
ANCHOR_MATCH_TOL_KMS: float = 0.5

RUNLOG = Path("data/runs/shooter-stm-batch.jsonl")


def _load_done(runlog: Path) -> set[tuple[str, int]]:
    """(id, epoch_index) pairs already recorded in the runlog (for --resume)."""
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
    """Append one record and flush+fsync so a kill loses nothing."""
    runlog.parent.mkdir(parents=True, exist_ok=True)
    with runlog.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        import os

        os.fsync(fh.fileno())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="skip (row, epoch) pairs already present in the runlog",
    )
    args = parser.parse_args()

    # Line-buffer stdout so per-row progress flushes promptly when piped (each
    # shoot() is minutes long; buffered output would hide all progress).
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

            epochs = candidate_epochs(
                ephem,
                0.0,
                launch_window_synodics=LAUNCH_WINDOW_SYNODICS,
                grid=EPOCH_GRID,
            )[:K]

            for epoch_index, t0 in enumerate(epochs):
                if (e.id, epoch_index) in done:
                    continue
                wall = time.monotonic() - t_wall0
                t_shoot0 = time.monotonic()
                row_error: str | None = None
                res: shooter.ShootResult | None = None
                try:
                    sseed = russell_shooting_seed(seed, t0_sec=t0, ephem=ephem)
                    res = shooter.shoot(
                        sseed,
                        ephem=ephem,
                        bodies=bodies,
                        accuracy=SHOOT_ACCURACY,
                        max_nfev=SHOOT_MAX_NFEV,
                        max_wall_sec=LEG_WALL_BUDGET_SEC,
                        jacobian=SHOOT_JACOBIAN,  # type: ignore[arg-type]
                    )
                except Exception as exc:  # honest per-(row,epoch) record, never raised
                    row_error = f"{type(exc).__name__}: {exc}"

                shoot_wall = time.monotonic() - t_shoot0
                if res is None:
                    rec = {
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
                    }
                    _append(RUNLOG, rec)
                    print(
                        f"[{wall:.0f}s] {e.id:24s} [{vlevel}] ep{epoch_index} "
                        f"NO SHOOT ({shoot_wall:.0f}s, error={row_error})"
                    )
                    continue

                vinf = list(res.vinf_per_encounter_kms)
                best_e = min((abs(v - seed.vinf_anchor_e_kms) for v in vinf), default=float("inf"))
                best_m = min((abs(v - seed.vinf_anchor_m_kms) for v in vinf), default=float("inf"))
                anchor_match = best_e <= ANCHOR_MATCH_TOL_KMS and best_m <= ANCHOR_MATCH_TOL_KMS
                converged = res.converged
                bend = res.bend_feasible
                promote = e.id == "mcconaghy-2006-em-k2" and converged and anchor_match and bend

                rec = {
                    "id": e.id,
                    "validation_level": vlevel,
                    "sequence": list(seed.sequence),
                    "epoch_index": epoch_index,
                    "epoch_sec": t0,
                    "shot": True,
                    "converged": converged,
                    "defect_norm": res.defect_norm,
                    "seed_defect_norm": res.seed_defect_norm,
                    "vinf_per_encounter_kms": vinf,
                    "correction_dv_kms": res.correction_dv_kms,
                    "bend_feasible": bend,
                    "n_iterations": res.n_iterations,
                    "anchor_e_kms": seed.vinf_anchor_e_kms,
                    "anchor_m_kms": seed.vinf_anchor_m_kms,
                    "best_e_residual_kms": best_e,
                    "best_m_residual_kms": best_m,
                    "anchor_match": anchor_match,
                    "promote_proposed_held": promote,
                    "jacobian": SHOOT_JACOBIAN,
                    "shoot_wall_sec": shoot_wall,
                    "error": row_error,
                }
                _append(RUNLOG, rec)
                flag = " *** PROPOSED V0->V1 (HELD) ***" if promote else ""
                print(
                    f"[{wall:.0f}s] {e.id:24s} [{vlevel}] ep{epoch_index} "
                    f"conv={converged} defect={res.defect_norm:.3e} "
                    f"vinf={[round(v, 2) for v in vinf]} "
                    f"anchorE/M={best_e:.2f}/{best_m:.2f} "
                    f"match={anchor_match} bend={bend} "
                    f"({shoot_wall:.0f}s, nfev={res.n_iterations}){flag}"
                )
        except Exception as exc:  # one bad row must not abort the batch
            _append(
                RUNLOG,
                {"id": e.id, "shot": False, "error": f"row-fatal {type(exc).__name__}: {exc}"},
            )
            print(f"[{wall:.0f}s] {e.id:24s} ROW-FATAL {type(exc).__name__}: {exc}")

    print()
    print("=" * 72)
    print(
        f"Batch complete. Runlog: {RUNLOG} "
        f"(K={K}, launch_window_synodics={LAUNCH_WINDOW_SYNODICS}, grid={EPOCH_GRID}, "
        f"max_nfev={SHOOT_MAX_NFEV}, jacobian={SHOOT_JACOBIAN})"
    )
    print("HELD — no writeback (data/catalogue.yaml and validate.py untouched).")


if __name__ == "__main__":
    main()
