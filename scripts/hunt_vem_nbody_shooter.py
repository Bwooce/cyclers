"""Phase 3a near-miss survey + n-body multiple-shooting hunt for the Jones VEM
family (task #133 / #136 Phase C).

Two stages, per the #135 verdict (seeding/basin, not solver) and the AAS 17-577
method deep-dive (broad-search near-Hohmann seeding -> SNOPT n-body correction):

  1. NEAR-MISS SURVEY (:func:`cyclerfinder.nbody.shooter.near_miss_survey`): scan
     an epoch grid over one 12.8-yr repeat period x rev/branch topologies, run the
     conic ``ballistic_correct`` at each point, and collect the lowest-V_inf chains
     within the relaxed near-miss tolerance (0.5 km/s). These are the shooting
     seeds (NOT a blind scan).
  2. N-BODY SHOOT (:func:`cyclerfinder.nbody.shooter.shoot`): multiple-shooting
     differential correction in restricted n-body from the best near-miss seed,
     driving the full-state defects toward ballistic (the SNOPT analogue).

NON-GOLDEN: every V_inf printed is OUR computation. The sourced Jones multiset is
read only for the proximity tag (a cross-check, never fitted). The headline gate
(tests/nbody/test_shooter_jones_gate.py) is the golden-disciplined comparison.

Run: uv run python -W ignore scripts/hunt_vem_nbody_shooter.py [n_epochs] [n_seeds]
"""

from __future__ import annotations

import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.shooter import (
    near_miss_survey,
    shoot,
    shooting_seed_from_near_miss,
)

CATALOGUE_PATH = Path(__file__).resolve().parent.parent / "data" / "catalogue.yaml"
J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
MEMBERS = ("jones-2017-vem-emevve-outbound", "jones-2017-vem-meevem-inbound")


def _row(entry_id: str) -> dict[str, Any]:
    for r in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if r["id"] == entry_id:
            return cast("dict[str, Any]", r)
    raise AssertionError(entry_id)


def _t_sec(iso: str) -> float:
    return (datetime.fromisoformat(iso).replace(tzinfo=UTC) - J2000).total_seconds()


def _sourced(entry_id: str) -> list[float]:
    row = _row(entry_id)
    return sorted(float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"])


def hunt(entry_id: str, n_epochs: int, n_seeds: int) -> None:
    row = _row(entry_id)
    seq = tuple(row["sequence_canonical"].split("-"))
    n_legs = len(seq) - 1
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY
    segs = row["trajectory"]["segments"]
    tof_seed = [float(segs[i]["tof_days"]) for i in range(n_legs)]
    t0_base = _t_sec(row["trajectory"]["epoch_tzero"])
    ephem = Ephemeris("astropy")

    print(f"\n== {entry_id}  seq={'-'.join(seq)}  period={row['period']['years']} yr")
    t0 = time.perf_counter()
    seeds = near_miss_survey(
        sequence=seq,
        period_sec=period_sec,
        t0_base_sec=t0_base,
        tof_seed_days=tof_seed,
        ephem=ephem,
        n_epochs=n_epochs,
        vinf_cap=8.0,
        near_miss_tol_kms=0.5,
    )
    dt_survey = time.perf_counter() - t0
    print(f"   near-miss survey: {len(seeds)} seeds within 0.5 km/s  ({dt_survey:.1f}s)")
    sourced = _sourced(entry_id)
    print(f"   sourced multiset (n={len(sourced)}): {[round(v, 2) for v in sourced]}")
    if not seeds:
        print("   NO near-miss seed found -> nothing to shoot (recorded finding).")
        return
    for s in seeds[:5]:
        print(
            f"     seed max-Vinf={s.max_vinf_kms:.2f} res={s.max_residual_kms:.3f} "
            f"bend={s.bend_feasible} vinfs={[round(v, 2) for v in s.vinf_per_encounter_kms]}"
        )

    for k, nm in enumerate(seeds[:n_seeds]):
        seed = shooting_seed_from_near_miss(nm, seq, period_sec, ephem)
        t0 = time.perf_counter()
        res = shoot(seed, ephem=ephem, bodies=("V", "E", "M"), accuracy=1e-9, max_nfev=80)
        print(
            f"   SHOOT seed#{k} (seed max-Vinf {nm.max_vinf_kms:.2f}): "
            f"converged={res.converged} defect {res.seed_defect_norm:.2e}->{res.defect_norm:.2e} "
            f"bend_feasible={res.bend_feasible} corrΔV={res.correction_dv_kms:.3f} "
            f"({time.perf_counter() - t0:.1f}s)"
        )
        corrected = sorted(round(v, 2) for v in res.vinf_per_encounter_kms)
        print(f"     corrected Vinf (sorted): {corrected}")


def main() -> None:
    n_epochs = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    for entry_id in MEMBERS:
        hunt(entry_id, n_epochs, n_seeds)


if __name__ == "__main__":
    main()
