"""N-body Phase C HEADLINE: Jones VEM ballistic rediscovery in n-body (plan Phase C).

GOLDEN: EXPECTED = catalogue SOURCED vinf_kms_at_encounters (AAS 17-577). The
shooter output is the side under test. Fold-corrected sorted-multiset compare,
VEM_VINF_TOL_KMS = 0.5. xfail until the near-miss-seeded shoot reaches the Jones
basin (honesty boundary + #135 verdict: convergence NOT assumed).
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

import pytest
import yaml  # type: ignore[import-untyped]

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.constants import (  # noqa: E402
    DAYS_PER_JULIAN_YEAR,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.shooter import (  # noqa: E402
    ShootResult,
    near_miss_survey,
    shoot,
    shooting_seed_from_near_miss,
)

CATALOGUE_PATH = Path(__file__).resolve().parents[2] / "data" / "catalogue.yaml"
VEM_VINF_TOL_KMS = 0.5
_VEM_MEMBERS = ("jones-2017-vem-emevve-outbound", "jones-2017-vem-meevem-inbound")
_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def _row(entry_id: str) -> dict[str, Any]:
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row["id"] == entry_id:
            return cast("dict[str, Any]", row)
    raise AssertionError(f"row {entry_id!r} not found")


def _sourced_multiset_one_period(entry_id: str) -> list[float]:
    """First-cycle sourced V_inf multiset (fold rule: the corrector solves ONE
    period -> len(sequence_canonical) encounters; the sourced table spans two)."""
    row = _row(entry_id)
    seq = tuple(row["sequence_canonical"].split("-"))
    n = len(seq)
    encs = row["vinf_kms_at_encounters"][:n]
    assert tuple(e["body"] for e in encs) == seq
    return sorted(float(e["vinf_kms"]) for e in encs)


def _t_sec(iso: str) -> float:
    return (datetime.fromisoformat(iso).replace(tzinfo=UTC) - _J2000).total_seconds()


def _shoot_for(entry_id: str) -> ShootResult:
    """Near-miss survey -> best low-V_inf seed -> n-body multiple-shooting solve.

    Bounded (small epoch grid, single topology) so the headline gate stays a
    tractable slow test; the dense survey is the offline driver
    (scripts/hunt_vem_nbody_shooter.py).
    """
    row = _row(entry_id)
    seq = tuple(row["sequence_canonical"].split("-"))
    n_legs = len(seq) - 1
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY
    segs = row["trajectory"]["segments"]
    tof_seed = [float(segs[i]["tof_days"]) for i in range(n_legs)]
    t0_base = _t_sec(row["trajectory"]["epoch_tzero"])
    ephem = Ephemeris("astropy")

    seeds = near_miss_survey(
        sequence=seq,
        period_sec=period_sec,
        t0_base_sec=t0_base,
        tof_seed_days=tof_seed,
        ephem=ephem,
        n_epochs=8,
        vinf_cap=8.0,
        near_miss_tol_kms=0.5,
    )
    if not seeds:
        pytest.skip(f"{entry_id}: no near-miss seed within tolerance (recorded finding)")
    seed = shooting_seed_from_near_miss(seeds[0], seq, period_sec, ephem)
    return shoot(seed, ephem=ephem, bodies=("V", "E", "M"), accuracy=1e-9, max_nfev=60)


@pytest.mark.slow
@pytest.mark.xfail(
    reason="n-body Jones shooter: rediscovery to the sourced AAS 17-577 multiset "
    "within 0.5 km/s. Open per the honesty boundary + #135 verdict (seeding/basin "
    "problem): the near-miss survey surfaces only the high-V_inf Lambert-chain "
    "basin (#110/#120/#122 floor ~11-18 km/s, 0 bend-feasible), and the n-body "
    "shoot from those seeds does not reach the Jones 2.5-3.9 km/s family. Flip "
    "ONLY when >=1 member converges within VEM_VINF_TOL_KMS; see plan Phase C "
    "Task C.4 + C.5 STOP/report branch.",
    strict=False,
)
@pytest.mark.parametrize("entry_id", _VEM_MEMBERS)
def test_jones_vem_nbody_rediscovers_sourced_multiset(entry_id: str) -> None:
    expected = _sourced_multiset_one_period(entry_id)
    result = _shoot_for(entry_id)
    assert result.converged
    got = sorted(result.vinf_per_encounter_kms)
    assert len(got) == len(expected)
    for g, x in zip(got, expected, strict=True):
        assert g == pytest.approx(x, abs=VEM_VINF_TOL_KMS)
