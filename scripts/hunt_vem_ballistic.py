"""Dense parallel epoch x branch hunt for the Jones VEM family (task #110, Part 2).

Drives the parallel scan engine (search/scan.py) over a full 12.8-yr repeat
period at fine epoch resolution x rev/branch topologies, for both Jones VEM
member rows, and reports every distinct closed ballistic family found (V_inf
range, bend feasibility). This is the "density lever" the prototype's main()
loop used to reach the ~6.4 Mars S1L1 family, generalised to the VEM cells.

NON-GOLDEN: every V_inf printed is OUR computation. The sourced Jones multiset
is read only for the near-anchor proximity tag (a CROSS-CHECK, not fitted).

Run: uv run python -W ignore scripts/hunt_vem_ballistic.py [n_epochs] [ephem_model] [residual_mode]

``ephem_model`` (default ``astropy``) is passed straight through to
``scan_parallel`` / ``Ephemeris(model=...)``; the corrector and scan engine are
ephemeris-agnostic, so ``inclined-circular`` (M-3D, sourced J2000 inc/Ω) re-runs
the identical grid on the inclined backend (task #120, the 3D-hypothesis test).

``residual_mode`` (default ``magnitude``) selects the corrector residual:
``magnitude`` is the #110/#120 baseline; ``vector`` is the Jones-method
bend-feasibility-aware residual (task #122 Phase 1, the falsifier run).
"""

from __future__ import annotations

import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.search.scan import ScanResult, build_epoch_branch_grid, scan_parallel

CATALOGUE_PATH = Path(__file__).resolve().parent.parent / "data" / "catalogue.yaml"

J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
MEMBERS = ("jones-2017-vem-emevve-outbound", "jones-2017-vem-meevem-inbound")


def _row(entry_id: str) -> dict[str, Any]:
    for r in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if r["id"] == entry_id:
            return cast("dict[str, Any]", r)
    raise AssertionError(entry_id)


def _t_sec(iso: str) -> float:
    dt = datetime.fromisoformat(iso).replace(tzinfo=UTC)
    return (dt - J2000).total_seconds()


def _topologies(n_legs: int) -> list[tuple[tuple[int, ...], tuple[str, ...]]]:
    """A small rev/branch grid: all-direct, plus a multi-rev variant on each
    interior leg (low and high branch). The corrector refines from each."""
    base_revs = (0,) * n_legs
    base_branch = ("single",) * n_legs
    topos: list[tuple[tuple[int, ...], tuple[str, ...]]] = [(base_revs, base_branch)]
    for leg in range(n_legs):
        for branch in ("low", "high"):
            revs = list(base_revs)
            br = list(base_branch)
            revs[leg] = 1
            br[leg] = branch
            topos.append((tuple(revs), tuple(br)))
    return topos


def hunt(
    entry_id: str,
    n_epochs: int,
    ephem_model: str = "astropy",
    residual_mode: str = "magnitude",
) -> list[ScanResult]:
    row = _row(entry_id)
    seq = tuple(row["sequence_canonical"].split("-"))
    n_legs = len(seq) - 1
    period_years = float(row["period"]["years"])
    period_sec = period_years * DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY

    # First-cycle sourced segment ToFs (one cycle = n_legs legs), the anchor seed.
    segs = row["trajectory"]["segments"]
    full_seed = [float(segs[i]["tof_days"]) for i in range(n_legs)]
    slack_leg = int(np.argmax(full_seed))
    free_seed = [t for i, t in enumerate(full_seed) if i != slack_leg]

    # Sourced per-body V_inf cross-check anchors.
    sourced = sorted(float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"])

    # Epoch grid: full 12.8-yr period at fine resolution centred on epoch_tzero.
    t0_base = _t_sec(row["trajectory"]["epoch_tzero"])
    period_days = period_years * DAYS_PER_JULIAN_YEAR
    step = period_days / n_epochs
    t0_seeds = [t0_base + i * step * SECONDS_PER_DAY for i in range(n_epochs)]

    grid = build_epoch_branch_grid(
        sequence=seq,
        period_sec=period_sec,
        vinf_cap=8.0,
        t0_seeds_sec=t0_seeds,
        branch_topologies=_topologies(n_legs),
        tof_seed_days=free_seed,
        slack_leg=slack_leg,
        residual_mode=residual_mode,
    )
    print(
        f"\n== {entry_id}  seq={'-'.join(seq)}  legs={n_legs}  "
        f"period={period_years} yr  grid={len(grid)} points "
        f"({len(_topologies(n_legs))} topos x {n_epochs} epochs)  "
        f"ephem={ephem_model}  residual={residual_mode}"
    )
    t0 = time.perf_counter()
    results = scan_parallel(grid, ephem_model=ephem_model, max_workers=16)
    dt = time.perf_counter() - t0
    closed = [r for r in results if r.closed]
    print(f"   scan {dt:.1f}s  closed={len(closed)}/{len(grid)}")

    # Distinct families: round per-encounter max-V_inf signature.
    families: dict[tuple[float, float, bool], ScanResult] = {}
    for r in closed:
        vmax = round(max(r.result.vinf_per_encounter_kms), 1)
        vmin = round(min(r.result.vinf_per_encounter_kms), 1)
        key = (vmin, vmax, r.result.bend_feasible)
        if key not in families or r.max_residual_kms < families[key].max_residual_kms:
            families[key] = r
    print(f"   distinct closed families (Vinf_min, Vinf_max, bend_feasible): {len(families)}")
    for key in sorted(families):
        r = families[key]
        vinfs = [round(v, 2) for v in r.result.vinf_per_encounter_kms]
        bf = "BEND-OK" if key[2] else "powered"
        flag = " <-- bend-feasible < 10" if (key[2] and key[1] < 10.0) else ""
        print(
            f"     Vinf[{key[0]:.1f}..{key[1]:.1f}] {bf} "
            f"res={r.max_residual_kms:.3f} t0={key} vinfs={vinfs}{flag}"
        )

    # Best (lowest max-Vinf) closed solution and proximity to sourced.
    if closed:
        best = min(closed, key=lambda r: max(r.result.vinf_per_encounter_kms))
        bv = sorted(best.result.vinf_per_encounter_kms)
        print(
            f"   BEST (lowest max-Vinf) closed: max-Vinf={max(bv):.2f} "
            f"res={best.max_residual_kms:.3f} bend_feasible={best.result.bend_feasible}"
        )
        print(f"     our per-encounter Vinf (sorted): {[round(v, 2) for v in bv]}")
        print(f"     sourced multiset (n={len(sourced)}): {[round(v, 2) for v in sourced]}")
    return closed


def main() -> None:
    n_epochs = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    ephem_model = sys.argv[2] if len(sys.argv) > 2 else "astropy"
    residual_mode = sys.argv[3] if len(sys.argv) > 3 else "magnitude"
    for entry_id in MEMBERS:
        hunt(entry_id, n_epochs, ephem_model, residual_mode)


if __name__ == "__main__":
    main()
