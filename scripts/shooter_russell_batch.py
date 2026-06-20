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

Coverage caps (DOCUMENTED — full n-body ``shoot()`` per epoch is expensive):

- ``K`` epochs per row: only the best-K phase-error candidate epochs (from the
  Russell §5.3 LaunchWindow scan) are actually shot. This is a logged coverage
  cap, not a physics choice; the lowest-defect ShootResult across the K epochs is
  retained per row.
- ``LAUNCH_WINDOW_SYNODICS`` / ``EPOCH_GRID`` size the candidate-epoch scan.

CAP REDUCTION (2026-06-21, measured): a single full-state ``defect_residual``
on these E-E-M-M rows (ToFs up to ~1026 days) costs ~40 s wall *regardless* of
``max_wall_sec`` (the individual legs each complete under budget, so the guard
never trips; it is the sum of 3 legs that dominates). The LM solver's internal
finite-difference Jacobian is (6*n_nodes + 1) ≈ 25 such residuals per iteration,
i.e. ~1000 s/iteration serially — intractable. The original plan's ``K=3`` /
``range(1, 22)`` was therefore reduced to ``K=1`` / ``range(1, 6)`` and the
Jacobian is parallelised across all cores (``N_JOBS``). This keeps each row to a
small number of LM iterations within the batch wall budget. The reduction is a
COVERAGE cap (fewer epochs probed), not a physics relaxation; it is recorded in
``docs/notes/2026-06-21-shooter-russell-results.md``.

Per the standing discipline: HELD — no writeback. This script never touches
``data/catalogue.yaml`` or ``validate.py``; it only writes a JSONL runlog and
prints a summary. A row is flagged ``PROPOSED V0->V1 (HELD)`` only for
``mcconaghy-2006-em-k2`` AND only if it converges, anchor-matches, and is
bend-feasible — and even then it is recorded, never applied.
"""

from __future__ import annotations

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

# --- Coverage caps (documented; see module docstring) -------------------------
# Full n-body shoot per epoch is expensive. Only the best-K phase-error epochs
# are run. If a single row exceeds the wall budget, K and the launch window are
# reduced (and that reduction is documented inline + in the results note).
# Reduced from the plan's K=3 / range(1, 22) after the measured ~40 s/residual
# cost (see module docstring "CAP REDUCTION"); K=1 shoots only the single
# best-phase epoch per row, range(1, 6) scans 5 LaunchWindows to find it.
K: int = 1
LAUNCH_WINDOW_SYNODICS: range = range(1, 6)
EPOCH_GRID: int = 100

# --- shoot() budget -----------------------------------------------------------
# accuracy mirrors tests/nbody/test_shooter_jones_gate.py. max_nfev is kept small
# (the convergence/family-selection signal, not a fine optimum, is the #388
# deliverable) and the Jacobian is parallelised across cores via N_JOBS.
SHOOT_ACCURACY: float = 1e-9
SHOOT_MAX_NFEV: int = 12
LEG_WALL_BUDGET_SEC: float = 5.0
N_JOBS: int = 16

# The anchor-match tolerance: an emerged per-encounter V∞ must land within this
# of the SOURCED E and M anchors for the row to count as anchor-matched.
ANCHOR_MATCH_TOL_KMS: float = 0.5


def main() -> None:
    # Line-buffer stdout so per-row progress flushes promptly when piped (each
    # shoot() is minutes long; buffered output would hide all progress).
    sys.stdout.reconfigure(line_buffering=True)

    m = RussellModel()
    ephem = Ephemeris("astropy")
    records: list[dict[str, Any]] = []

    t_wall0 = time.monotonic()
    catalog = load_catalog()

    descriptor_rows = 0
    converged_rows = 0
    anchor_matched_rows = 0
    promotions = 0

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

            descriptor_rows += 1
            vlevel = str(e.raw.get("validation_level", "V0"))
            bodies = tuple(dict.fromkeys(seed.sequence))

            epochs = candidate_epochs(
                ephem,
                0.0,
                launch_window_synodics=LAUNCH_WINDOW_SYNODICS,
                grid=EPOCH_GRID,
            )[:K]

            best: shooter.ShootResult | None = None
            best_epoch: float | None = None
            row_error: str | None = None
            for t0 in epochs:
                try:
                    sseed = russell_shooting_seed(seed, t0_sec=t0, ephem=ephem)
                    res = shooter.shoot(
                        sseed,
                        ephem=ephem,
                        bodies=bodies,
                        accuracy=SHOOT_ACCURACY,
                        max_nfev=SHOOT_MAX_NFEV,
                        max_wall_sec=LEG_WALL_BUDGET_SEC,
                        n_jobs=N_JOBS,
                    )
                except Exception as exc:  # honest per-epoch record, never raised
                    row_error = f"{type(exc).__name__}: {exc}"
                    continue
                if best is None or res.defect_norm < best.defect_norm:
                    best = res
                    best_epoch = t0

            if best is None:
                rec = {
                    "id": e.id,
                    "validation_level": vlevel,
                    "sequence": list(seed.sequence),
                    "shot": False,
                    "error": row_error,
                    "anchor_e_kms": seed.vinf_anchor_e_kms,
                    "anchor_m_kms": seed.vinf_anchor_m_kms,
                }
                records.append(rec)
                print(f"[{wall:.0f}s] {e.id:24s} [{vlevel}] NO SHOOT (error={row_error})")
                continue

            vinf = list(best.vinf_per_encounter_kms)
            best_e = min(abs(v - seed.vinf_anchor_e_kms) for v in vinf) if vinf else float("inf")
            best_m = min(abs(v - seed.vinf_anchor_m_kms) for v in vinf) if vinf else float("inf")
            anchor_match = best_e <= ANCHOR_MATCH_TOL_KMS and best_m <= ANCHOR_MATCH_TOL_KMS
            converged = best.converged
            bend = best.bend_feasible
            promote = e.id == "mcconaghy-2006-em-k2" and converged and anchor_match and bend

            if converged:
                converged_rows += 1
            if anchor_match:
                anchor_matched_rows += 1
            if promote:
                promotions += 1

            rec = {
                "id": e.id,
                "validation_level": vlevel,
                "sequence": list(seed.sequence),
                "shot": True,
                "converged": converged,
                "defect_norm": best.defect_norm,
                "seed_defect_norm": best.seed_defect_norm,
                "vinf_per_encounter_kms": vinf,
                "correction_dv_kms": best.correction_dv_kms,
                "bend_feasible": bend,
                "n_iterations": best.n_iterations,
                "best_epoch_sec": best_epoch,
                "anchor_e_kms": seed.vinf_anchor_e_kms,
                "anchor_m_kms": seed.vinf_anchor_m_kms,
                "best_e_residual_kms": best_e,
                "best_m_residual_kms": best_m,
                "anchor_match": anchor_match,
                "promote_proposed_held": promote,
                "error": row_error,
            }
            records.append(rec)
            flag = " *** PROPOSED V0->V1 (HELD) ***" if promote else ""
            print(
                f"[{wall:.0f}s] {e.id:24s} [{vlevel}] conv={converged} "
                f"defect={best.defect_norm:.3e} "
                f"vinf={[round(v, 2) for v in vinf]} "
                f"anchorE/M={best_e:.2f}/{best_m:.2f} "
                f"match={anchor_match} bend={bend}{flag}"
            )
        except Exception as exc:  # one bad row must not abort the batch
            records.append(
                {"id": e.id, "shot": False, "error": f"row-fatal {type(exc).__name__}: {exc}"}
            )
            print(f"[{wall:.0f}s] {e.id:24s} ROW-FATAL {type(exc).__name__}: {exc}")

    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    out = Path("data/runs") / f"shooter-russell-{stamp}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")

    print()
    print("=" * 72)
    print(
        f"Batch summary: descriptor rows {descriptor_rows}, converged "
        f"{converged_rows}, anchor-matched {anchor_matched_rows}, "
        f"proposed-promotions {promotions} (held)."
    )
    print(f"Runlog: {out}")
    print(
        f"K={K}, launch_window_synodics={LAUNCH_WINDOW_SYNODICS}, grid={EPOCH_GRID}, "
        f"max_nfev={SHOOT_MAX_NFEV}, n_jobs={N_JOBS}"
    )
    print("HELD — no writeback (data/catalogue.yaml and validate.py untouched).")


if __name__ == "__main__":
    main()
