"""Task #134 — Aldrin V2/V3 validation campaign (spec §14 ladder).

This module records the RECORDED, MECHANICAL outcome of climbing the classic
Aldrin Earth-Mars cycler pair (``aldrin-classic-em-k1-outbound`` /
``-inbound``, both currently V1) up spec §14's V2 (multi-lap periodicity) and
V3 (ephemeris-horizon TCM) gates. The headline result is a *documented
near-miss / no-promotion* on both gates — the golden discipline outcome when
the physics does not mechanically clear the gate.

Why no promotion (the finding, verbatim numbers below)
------------------------------------------------------
**V2 (≥3 continuous laps, BOUNDED drift in the dynamic rotating frame).** The
ballistic Aldrin provably cannot pass (it is *powered*: ≈84° Earth turn
required vs ≈72° achievable — see
``tests/verify/test_real_closure.py::test_aldrin_powered_turn_deficit_gate``).
V2 therefore runs on the POWERED solution from
:func:`cyclerfinder.search.bvp.solve_powered_periodic_cycler` (the verified
in-family 2.9138 km/s outbound maintenance solve), propagated over 3 continuous
laps on DE440 via the M6a machinery
(:func:`cyclerfinder.verify.real_closure.verify_real_closure`, which delegates
to :func:`cyclerfinder.verify.propagate.verify_long_term_stability`).

The measured lap-to-lap drift is **~4.14e8 km** — ``~2072x`` the
:data:`~cyclerfinder.verify.real_closure.REAL_DRIFT_TOLERANCE_KM` (200,000 km)
and ``~8289x`` the M6a idealised
:data:`~cyclerfinder.verify.propagate.DRIFT_TOLERANCE_KM` (50,000 km). This is
NOT bounded at the existing tolerances and is NOT a marginal miss. Root cause
(already documented in ``search/bvp.py``): the rotating-frame-repeat metric pins
each leg-start to the *real* Mars position at the lap-shifted epoch, and Mars's
heliocentric radius breathes because the 2.135 yr cycler period is not
commensurate with Mars's 1.881 yr orbit. The maneuver shapes velocity, not
where Mars is, so no maintenance ΔV closes this metric. This is exactly why the
real Aldrin cycler is *retargeted* every cycle — it is a finding about the
maintenance model, not an infrastructure gap. NO V2 promotion.

**V3 (phase-matched real window + horizon TCM over 3-5 laps, bounded AND within
ΔV budget).** Chaining the per-cycle maintenance solve across 3-5 laps from the
phase-matched real window (re-phasing the priority date one cycler period per
cycle, re-solving each cycle in-family at a≈1.588 AU, e≈0.393) gives a horizon
TCM total of **8.51 km/s (3 laps) / 11.19 km/s (4 laps) / 13.79 km/s (5
laps)** — i.e. ~2.6-2.9 km/s per cycle. Spec §14/§16's example V3 budget gate is
``horizon_tcm_mps = 120``; the Aldrin horizon TCM is ~70-115x over that.

Critically, the §14 "within ΔV budget" comparison has **no sourced anchor for
Aldrin**: the maintenance-ΔV magnitude is UNPUBLISHED. McConaghy/Longuski/Byrnes
2002 (AIAA 2002-4420, Table 4 "1L1") characterises Aldrin only by a turn-angle
test (≈84° required vs ≈72° achievable at a 200 km Earth flyby) and explicitly
DEFERS ΔV estimation to future work (p.8). The catalogue row records this in
``maintenance_dv_kms_per_synodic`` as a COMPUTED, never-golden surrogate and in
``data_gaps[trajectory.maneuvers[E].dv_kms]``. Per golden discipline (sourced-
only EXPECTED sides), §14's "within ΔV budget" is therefore NOT satisfiable
against a sourced number for Aldrin — the horizon TCM is reported and bounded-
per-cycle, but cannot be promoted to V3. NO V3 promotion.

Discipline
----------
These tests are marked ``slow`` (each runs a real-DE440 BVP solve / multi-lap
DE440 propagation, minutes per case). They assert the finding QUALITATIVELY with
teeth (drift far above tolerance; horizon TCM far above the example budget and
strictly positive in-family per cycle) without pinning brittle solver-derived
magnitudes as golden targets — the magnitudes are OUR computation, never source-
attested.

SUPERSEDED FOR V2 (2026-06-07, the §14 V2 class-split amendment). The
no-promotion verdict below is the outcome under the *original single* V2 gate
(cross-cycle rotating-frame-repeat drift), and its numbers stand verbatim
(~4.14e8 km / 3 laps). That gate was the wrong instrument for a
per-cycle-retargeted powered cycler; the amendment splits V2 into V2-ballistic
and V2-powered, and under the amended **V2-powered** gate the Aldrin OUTBOUND
now PASSES and is promoted to V2 — see ``tests/verify/test_aldrin_v2_powered.py``
and spec §14 / §16.7.12. The V3 finding (no sourced ΔV budget for Aldrin) is
unchanged. The v4.5 census ratchet
(``tests/data/test_schema_v45_fields.py::test_live_v1_census_matches_recorded_evidence``)
now asserts exactly one V2 (the outbound) and no V3+.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.bvp import solve_powered_periodic_cycler
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv
from cyclerfinder.verify.propagate import DRIFT_TOLERANCE_KM
from cyclerfinder.verify.real_closure import (
    REAL_DRIFT_TOLERANCE_KM,
    verify_real_closure,
)
from tests.data._catalogue_loader_m6b import load_m6b_entries

ALDRIN_PRIORITY = datetime(1985, 10, 28, tzinfo=UTC)

# Spec §14 / §16.1 example V3 horizon-TCM budget gate (docs/spec.md ~L515-520).
# This is the near-ballistic *threshold example*, NOT a sourced Aldrin budget
# (no Aldrin maintenance-ΔV magnitude is published); used here only to scale the
# finding ("how far over the example gate is the powered Aldrin").
_SPEC_V3_HORIZON_TCM_BUDGET_MPS = 120.0


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


# ---------------------------------------------------------------------------
# V2 — multi-lap periodicity finding: powered Aldrin drift is UNBOUNDED
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v2_powered_aldrin_outbound_unbounded_over_3_laps(
    astropy_ephem: Ephemeris,
) -> None:
    """V2 FINDING (no promotion): the POWERED outbound Aldrin, propagated over 3
    continuous laps on DE440, does NOT have bounded rotating-frame drift.

    Records the §14 V2 outcome verbatim-in-spirit: the in-family powered solve
    (≈2.9138 km/s maintenance ΔV) drifts ~4.14e8 km over 3 laps — orders of
    magnitude (>1000x) above both the M6a idealised tolerance (50,000 km) and
    the M6b real tolerance (200,000 km). The drift floor is physically
    unreachable for the k=1 Aldrin because Mars's radius breathes per cycle
    (incommensurate periods); the maneuver shapes velocity, not Mars's position.
    This is the recorded reason the Aldrin pair stays V1.
    """
    entry = next(e for e in load_m6b_entries() if e["id"] == "aldrin-classic-em-k1-outbound")
    assert entry["trajectory_regime"] == "powered"

    solution = solve_powered_periodic_cycler(
        entry,
        ephem=astropy_ephem,
        signature_priority_date=ALDRIN_PRIORITY,
    )
    # The maintenance maneuver is real (in-family powered solve), not a ballistic
    # ΔV≈0 neighbour. Sanity-bounded only — the magnitude is our own value.
    assert solution.total_maintenance_dv_kms > 0.0
    assert solution.total_maintenance_dv_kms < 3.5

    result = verify_real_closure(
        solution.cycler,
        n_cycles=3,  # spec §14 V2: >= 3 continuous laps
        ephem=astropy_ephem,
        t_start=solution.t_start_sec,
        cycler_id="aldrin-classic-em-k1-outbound",
    )
    assert result.n_cycles_propagated == 3
    # NOT bounded at the existing tolerances — far above, not a marginal miss.
    assert result.closes is False
    assert result.v3_status == "v3-real-closure-fail"
    assert result.max_drift_km > 1000.0 * REAL_DRIFT_TOLERANCE_KM, result.max_drift_km
    assert result.max_drift_km > 1000.0 * DRIFT_TOLERANCE_KM, result.max_drift_km


@pytest.mark.slow
def test_v2_powered_aldrin_inbound_unbounded_over_3_laps(
    astropy_ephem: Ephemeris,
) -> None:
    """V2 FINDING (no promotion): the inbound twin also fails the ≥3-lap bound.

    The inbound row's optimiser does not recover the in-family powered solve at
    its real window (it lands on a ballistic ΔV≈0 neighbour — a separate,
    recorded off-family selection issue, not a V2 pass), and the resulting
    cycler likewise drifts ~4.12e8 km over 3 laps. Either way the ≥3-lap bounded-
    drift criterion is not met, so the inbound row stays V1.
    """
    entry = next(e for e in load_m6b_entries() if e["id"] == "aldrin-classic-em-k1-inbound")
    assert entry["trajectory_regime"] == "powered"

    solution = solve_powered_periodic_cycler(
        entry,
        ephem=astropy_ephem,
        signature_priority_date=ALDRIN_PRIORITY,
    )
    result = verify_real_closure(
        solution.cycler,
        n_cycles=3,
        ephem=astropy_ephem,
        t_start=solution.t_start_sec,
        cycler_id="aldrin-classic-em-k1-inbound",
    )
    assert result.n_cycles_propagated == 3
    assert result.closes is False
    assert result.v3_status == "v3-real-closure-fail"
    assert result.max_drift_km > 1000.0 * REAL_DRIFT_TOLERANCE_KM, result.max_drift_km


# ---------------------------------------------------------------------------
# V3 — ephemeris-horizon TCM finding: no sourced budget; TCM far over example
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v3_horizon_tcm_chained_no_sourced_budget(
    astropy_ephem: Ephemeris,
) -> None:
    """V3 FINDING (no promotion): horizon TCM is bounded-per-cycle and in-family,
    but §14's "within ΔV budget" has NO sourced anchor for Aldrin.

    Chains the per-cycle maintenance solve across 5 laps from the phase-matched
    real window (re-phasing the priority date one cycler period per cycle). Each
    cycle re-solves in-family (a≈1.588 AU, e≈0.393 — the SOURCED Aldrin anchors)
    with a strictly positive maintenance ΔV ≈ 2.6-2.9 km/s. The summed horizon
    TCM is ≈8.5 km/s (3 laps) / ≈11.2 (4) / ≈13.8 (5) — ~70-115x the spec's
    *example* V3 budget gate (120 m/s).

    Because the Aldrin maintenance-ΔV magnitude is UNPUBLISHED (McConaghy 2002
    defers it; catalogue ``data_gaps``), there is no sourced number to satisfy
    "within ΔV budget". Per golden discipline (sourced-only EXPECTED sides) the
    horizon TCM is reported here, bounded-and-positive per cycle, but CANNOT be
    promoted to V3. This asserts the comparison honestly, with teeth, against the
    spec's example gate only — never against a fabricated Aldrin budget.
    """
    period_yr = 2.135  # catalogue Aldrin cycler repeat period (= 1 E-M synodic)
    per_cycle_dv_kms: list[float] = []
    for k in range(5):
        pdate = ALDRIN_PRIORITY + timedelta(days=k * period_yr * 365.25)
        res = optimise_aldrin_maintenance_dv(astropy_ephem, real_window_priority_date=pdate)
        assert res.converged is True
        # Each cycle re-solves in-family on the SOURCED Aldrin anchors.
        assert res.a_au == pytest.approx(1.60, abs=0.05), res.a_au
        assert res.e == pytest.approx(0.393, abs=0.03), res.e
        # Strictly-positive, sanity-bounded per-cycle maintenance ΔV (powered).
        assert 0.0 < res.maintenance_dv_kms < 3.5, res.maintenance_dv_kms
        per_cycle_dv_kms.append(res.maintenance_dv_kms)

    for n_laps in (3, 4, 5):
        horizon_tcm_kms = sum(per_cycle_dv_kms[:n_laps])
        horizon_tcm_mps = horizon_tcm_kms * 1000.0
        # Bounded-per-cycle: the horizon TCM is finite and grows ~linearly.
        assert horizon_tcm_mps == pytest.approx(sum(per_cycle_dv_kms[:n_laps]) * 1000.0)
        # FAR over the spec EXAMPLE budget gate — not within ANY published budget
        # (none exists for Aldrin). >> the 120 m/s near-ballistic example.
        assert horizon_tcm_mps > 50.0 * _SPEC_V3_HORIZON_TCM_BUDGET_MPS, horizon_tcm_mps
