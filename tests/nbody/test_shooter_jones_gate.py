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
# BOUNDED GATE: both Jones VEM rows are 6-node / 5-leg itineraries (E-M-E-V-V-E
# and M-E-E-V-E-M). A full LM shoot is a 36-DOF finite-difference Jacobian over
# 5-leg REBOUND/IAS15 propagations; running BOTH members unbounded is what overran
# 60+ CPU-min (each leg could burn the propagator's 90 s wall budget on the
# high-V∞ basin's step-collapsing seeds). The headline question needs ONE member's
# worth of bounded signal, not two — so the gate shoots the single lowest-V∞-floor
# member and caps wall via max_nfev + a tight per-leg max_wall_sec.
_VEM_MEMBERS = ("jones-2017-vem-meevem-inbound",)
# BUDGET ARITHMETIC (binding, #133 Phase C): this member is 6 nodes / 5 legs, so
# the multiple-shooting state vector is 36-DOF. scipy LM's finite-difference
# Jacobian costs one residual eval per parameter, i.e. ONE Jacobian sweep =
# 36 residual_of_x calls = 36 * 5 = 180 leg propagations. The FD-Jacobian cost
# model is n_LM_steps * (n_params + 1) sweeps * n_legs * t_leg. To buy at least
# ONE genuine LM step the solver must afford a full 36-column Jacobian, so
# max_nfev must be >= ~38 (1 seed eval + 36 FD columns + a trial step); fewer and
# the solve is a no-op (Jacobian never completes, no step taken).
#   total residual calls ~= 1 (seed_res) + max_nfev (inside lm) + 1 (final_res)
#   worst-case wall      ~= total_calls * n_legs * max_wall_sec
# With max_nfev=38, max_wall_sec=3 s: (1+38+1) * 5 * 3 = 600 s = 10 min modelled
# worst case (every leg burning its full budget), which lands inside the 15-min
# wall cap with slack for the between-chunks budget-check granularity. The
# per-leg wall budget is the propagator's own divergence budget: a high-V∞ basin
# seed that collapses IAS15 steps short-circuits on this instead of grinding the
# 90 s default. Run on its own: `-m slow -o addopts="--timeout=900"`.
# NOTE: the 24-min figure a predecessor recorded here was from a run whose stdout
# was lost to a session kill — it was UNVERIFIED and is not relied on.
_LEG_WALL_BUDGET_SEC = 3.0
_SHOOT_MAX_NFEV = 38
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

    Bounded on every axis so the headline gate stays a tractable slow test: small
    epoch grid (n_epochs=8), single topology, single member, capped solver budget
    (max_nfev) AND a tight per-leg propagation wall budget (max_wall_sec) so a
    high-V∞ basin seed that collapses IAS15 steps short-circuits instead of
    grinding. The dense survey is the offline driver
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
    return shoot(
        seed,
        ephem=ephem,
        bodies=("V", "E", "M"),
        accuracy=1e-9,
        max_nfev=_SHOOT_MAX_NFEV,
        max_wall_sec=_LEG_WALL_BUDGET_SEC,
    )


@pytest.mark.slow
@pytest.mark.xfail(
    reason="n-body Jones shooter: rediscovery to the sourced AAS 17-577 multiset "
    "within 0.5 km/s. Open per the honesty boundary + #135 verdict (seeding/basin "
    "problem): for these 6-node VEM rows the near-miss survey surfaces only the "
    "high-V_inf Lambert-chain basin (max-V_inf ~29-33 km/s, 0 bend-feasible vs the "
    "sourced Jones 2.5-5.2 km/s family), and the bounded n-body shoot from those "
    "seeds does not reach it. Flip ONLY when the member converges within "
    "VEM_VINF_TOL_KMS; see plan Phase C Task C.4 + C.5 STOP/report branch.",
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
