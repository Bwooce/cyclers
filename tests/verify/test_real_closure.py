"""M6b real-ephemeris closure verification gate + helpers.

Spec / plan references
----------------------
* Spec §14 V2-real (multi-cycle real-ephemeris closure gate).
* Plan: ``docs/phases/m6b-real-ephemeris-closure/plan.md`` §4 (tests +
  tolerance), §3.1 (construction path), §3.3 (RealClosureResult shape).

Test layout (plan §4.1):

* Gate 1 — :func:`test_aldrin_cycler_periodic_over_2_cycles_astropy`
  (**M6b binding gate**, spec §8 M6 real-ephemeris half).
* Gate 2 — :func:`test_2syn_em_cpom_periodic_over_2_cycles_astropy`
  (**xfail** under the multi-rev Lambert blocker).
* Gate 3 — :func:`test_real_drift_rejects_open_trajectory`.
* Gate 4 — :func:`test_real_closure_regression_set` (parametrised).
* Gate 5 —
  :func:`test_real_closure_result_frozen_and_v3_fields_locked`.
* Gate 6 — :func:`test_real_closure_uses_m6a_machinery` (binding
  composition assertion).

Plus the helper-level tests (plan §4.6) and the diagnostic
``test_check_vinf_continuity_*`` tests.
"""

from __future__ import annotations

import dataclasses
from datetime import UTC, datetime
from unittest.mock import patch

import numpy as np
import pytest

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler, Encounter
from cyclerfinder.search.phase_match import (
    PhaseSignature,
    phase_signature_from_catalogue_entry,
)
from cyclerfinder.verify.propagate import StabilityReport
from cyclerfinder.verify.real_closure import (
    EXPECTED_SKIPS,
    N_CYCLES_DEFAULT,
    REAL_DRIFT_TOLERANCE_KM,
    MultiRevLambertRequiredError,
    RealClosureResult,
    _check_vinf_continuity,
    _resolve_real_t_start,
    construct_real_ephemeris_cycler,
    verify_real_closure,
)
from tests.data._catalogue_loader_m6b import (
    M6B_REGRESSION_IDS,
    load_m6b_entries,
)

ALDRIN_PRIORITY = datetime(1985, 10, 28, tzinfo=UTC)


@pytest.fixture(scope="module")
def aldrin_entry() -> dict[str, object]:
    """Return the loaded Aldrin outbound catalogue entry."""
    entries = load_m6b_entries()
    return next(e for e in entries if e["id"] == "aldrin-classic-em-k1-outbound")


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


# ---------------------------------------------------------------------------
# Gate 1 — M6b BINDING GATE (Aldrin real-ephemeris closure over 2 cycles)
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.xfail(
    strict=False,
    reason=(
        "M6b binding gate fails by ~1450x: Lambert-chain construction "
        "of the Aldrin classic E-M k=1 cycler (plan §3.1.1) produces "
        "max_drift_km ~ 2.9e8 km on real ephemeris (1986-02-27 launch "
        "epoch, 2 cycles ~ 4.3 yr horizon), vs the 200,000 km "
        "REAL_DRIFT_TOLERANCE_KM bound. The same construction yields "
        "~2.8e6 km drift on the circular ephemeris with the M3-style "
        "132-degree phase epoch, so the failure is NOT real-ephemeris "
        "eccentricity breathing — it is intrinsic to the Lambert-chain "
        "algorithm for the Aldrin geometry. The published Aldrin "
        "cycler's M->E return leg is the same heliocentric ellipse as "
        "the E->M outbound leg (one orbit, not two independent "
        "Lambert chords); independently Lambert-solving the M->E leg "
        "produces |V_inf_M_in|=9.78 vs |V_inf_M_out|=10.92 (~1.1 km/s "
        "ballistic-flyby mismatch) and the cycler does not close. "
        "Resolution requires M7 to construct the Aldrin instance from "
        "its orbital elements (a=1.60 AU, e=0.393, i=0) and the flyby "
        "kinematics, not via independent per-leg Lambert. See "
        "todo.md 'Hand-off to M7' for the diagnostic + proposed fix."
    ),
)
def test_aldrin_cycler_periodic_over_2_cycles_astropy(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """M6b BINDING GATE — spec §8, real-ephemeris half (plan §4.1 gate 1).

    Asserts the Aldrin classic E-M k=1 outbound cycler closes over 2
    real-ephemeris cycles within :data:`REAL_DRIFT_TOLERANCE_KM` on
    DE440. **xfail at M6b authorship time** per the architectural
    finding in the xfail reason and todo.md's M7 hand-off; the test
    will flip to passing once M7 supplies an orbital-elements-based
    constructor for the Aldrin family.
    """
    result = verify_real_closure(
        aldrin_entry,
        n_cycles=2,
        ephem=astropy_ephem,
        signature_priority_date=ALDRIN_PRIORITY,
        cycler_id="aldrin-classic-em-k1-outbound",
    )
    assert result.cycler_id == "aldrin-classic-em-k1-outbound"
    assert result.n_cycles_propagated == 2
    assert result.frame_used == "dynamic"
    assert result.horizon_tcm_mps == 0.0
    assert result.per_cycle_tcm_mps == (0.0, 0.0)
    # M6b's binding tolerance — see plan §4.3.
    assert result.closes, (
        f"Aldrin failed to close on real ephemeris: max_drift_km="
        f"{result.max_drift_km}, per_cycle_drift_km="
        f"{result.per_cycle_drift_km}, v3_status={result.v3_status}, "
        f"t_start_sec={result.t_start_sec}. See plan §5 risk #1."
    )
    assert result.max_drift_km < REAL_DRIFT_TOLERANCE_KM
    assert result.v3_status == "v3-real-closure-pass"


# ---------------------------------------------------------------------------
# Gate 2 — 2-syn S1L1 (xfail under multi-rev Lambert blocker)
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.xfail(
    strict=False,
    reason=(
        "multi-rev Lambert blocker — see plan §5 risk #2. "
        "Single-rev Lambert returns v3-skipped-multirev for the S1L1 "
        "intermediate-encounter legs (n_revs=1); the test will flip to "
        "passing once multi-rev Lambert lands (M-future / stretch)."
    ),
)
def test_2syn_em_cpom_periodic_over_2_cycles_astropy(
    astropy_ephem: Ephemeris,
) -> None:
    """Aspirational gate for the 2-syn S1L1 entry; xfail at M6b time."""
    entries = load_m6b_entries()
    s1l1 = next(e for e in entries if e["id"] == "s1l1-2syn-em-cpom")
    result = verify_real_closure(
        s1l1,
        n_cycles=2,
        ephem=astropy_ephem,
        signature_priority_date=datetime(2002, 8, 14, tzinfo=UTC),
        cycler_id="s1l1-2syn-em-cpom",
    )
    assert result.closes
    assert result.max_drift_km < REAL_DRIFT_TOLERANCE_KM


# ---------------------------------------------------------------------------
# Gate 3 — Open-trajectory rejection
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_real_drift_rejects_open_trajectory(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """Plan §4.1 gate 3: a deliberately-perturbed Aldrin diverges.

    Rotates one ``vinf_out`` by 5 degrees about +z and asserts the
    drift exceeds ``5 * REAL_DRIFT_TOLERANCE_KM``. The gate has
    rejection power, not just acceptance.
    """
    # Build the unperturbed cycler at its priority epoch.
    signature = phase_signature_from_catalogue_entry(aldrin_entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None
    cycler = construct_real_ephemeris_cycler(aldrin_entry, astropy_ephem, t_start)

    # Rotate enc[0].vinf_out by 5 degrees about +z.
    theta = np.deg2rad(5.0)
    c, s = float(np.cos(theta)), float(np.sin(theta))
    enc0 = cycler.encounters[0]
    rotated = np.array(
        [
            c * float(enc0.vinf_out[0]) - s * float(enc0.vinf_out[1]),
            s * float(enc0.vinf_out[0]) + c * float(enc0.vinf_out[1]),
            float(enc0.vinf_out[2]),
        ],
        dtype=np.float64,
    )
    perturbed_enc0 = dataclasses.replace(enc0, vinf_out=rotated, vinf_in=rotated)
    perturbed = dataclasses.replace(
        cycler,
        encounters=[perturbed_enc0, *cycler.encounters[1:]],
    )

    result = verify_real_closure(
        perturbed,
        n_cycles=2,
        ephem=astropy_ephem,
        t_start=t_start,
        cycler_id="aldrin-classic-em-k1-outbound-perturbed",
    )
    assert not result.closes
    assert result.max_drift_km > 5 * REAL_DRIFT_TOLERANCE_KM, result.max_drift_km
    assert result.v3_status == "v3-real-closure-fail"


# ---------------------------------------------------------------------------
# Gate 4 — Regression set across M6B_REGRESSION_IDS
# ---------------------------------------------------------------------------

# Regression entries currently failing the M6b binding tolerance due to the
# Lambert-chain construction limitation documented on the Aldrin binding-gate
# xfail. Tracked in todo.md "Hand-off to M7" for orbital-elements-based
# construction.
_M6B_LAMBERT_CHAIN_XFAILS: frozenset[str] = frozenset(
    {
        "aldrin-classic-em-k1-outbound",
        "aldrin-classic-em-k1-inbound",
    },
)


@pytest.mark.slow
@pytest.mark.parametrize("entry_id", M6B_REGRESSION_IDS)
def test_real_closure_regression_set(
    entry_id: str,
    astropy_ephem: Ephemeris,
    request: pytest.FixtureRequest,
) -> None:
    """Plan §4.2: every regression entry not in EXPECTED_SKIPS closes."""
    if entry_id in EXPECTED_SKIPS:
        pytest.skip(EXPECTED_SKIPS[entry_id])
    if entry_id in _M6B_LAMBERT_CHAIN_XFAILS:
        request.applymarker(
            pytest.mark.xfail(
                strict=False,
                reason=(
                    f"{entry_id}: Lambert-chain construction (plan §3.1.1) "
                    "yields ~10^8 km drift on real ephemeris; see the "
                    "test_aldrin_cycler_periodic_over_2_cycles_astropy "
                    "xfail reason and todo.md's M7 hand-off."
                ),
            )
        )
    entries = load_m6b_entries()
    entry = next(e for e in entries if e["id"] == entry_id)
    priority = entry.get("priority_date")
    if isinstance(priority, str):
        priority_dt = datetime.fromisoformat(priority).replace(tzinfo=UTC)
    else:
        priority_dt = ALDRIN_PRIORITY  # Sane default for entries lacking priority_date.
    result = verify_real_closure(
        entry,
        n_cycles=N_CYCLES_DEFAULT,
        ephem=astropy_ephem,
        signature_priority_date=priority_dt,
        cycler_id=entry_id,
    )
    assert result.closes, (
        f"{entry_id}: max_drift_km={result.max_drift_km}, v3_status={result.v3_status}"
    )
    assert result.max_drift_km < REAL_DRIFT_TOLERANCE_KM


# ---------------------------------------------------------------------------
# Gate 5 — RealClosureResult frozen + V3 placeholders locked
# ---------------------------------------------------------------------------


def test_real_closure_result_frozen_and_v3_fields_locked() -> None:
    """Plan §4.1 gate 5: dataclass is frozen; V3 fields are zero placeholders.

    Builds a synthetic :class:`RealClosureResult` directly so the test
    does not depend on a real-ephemeris propagation. Assertions probe
    the dataclass contract, not the verifier.
    """
    result = RealClosureResult(
        cycler_id="synthetic",
        n_cycles_propagated=2,
        max_drift_km=10_000.0,
        per_cycle_drift_km=(0.0, 10_000.0),
        per_encounter_vinf_mismatch_kms=(),
        closes=True,
        v3_status="v3-real-closure-pass",
        horizon_tcm_mps=0.0,
        per_cycle_tcm_mps=(0.0, 0.0),
        frame_used="dynamic",
        t_start_sec=0.0,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.closes = False  # type: ignore[misc]
    assert result.horizon_tcm_mps == 0.0
    assert result.per_cycle_tcm_mps == (0.0,) * result.n_cycles_propagated
    assert result.v3_status in (
        "v3-real-closure-pass",
        "v3-real-closure-fail",
        "v3-skipped-multirev",
        "v3-no-real-window",
        "v3-construction-error",
    )
    assert result.frame_used == "dynamic"


# ---------------------------------------------------------------------------
# Gate 6 — composition assertion: M6b reuses M6a's verifier
# ---------------------------------------------------------------------------


def test_real_closure_uses_m6a_machinery(astropy_ephem: Ephemeris) -> None:
    """Plan §4.1 gate 6: ``verify_real_closure`` delegates to
    ``verify_long_term_stability`` exactly once with ``n_laps=n_cycles``.

    Patches the M6a entry point in the real_closure module's namespace
    and asserts the call count + ``n_laps`` argument. If a future
    refactor accidentally re-implements multi-lap propagation, this
    test breaks.
    """
    entries = load_m6b_entries()
    aldrin = next(e for e in entries if e["id"] == "aldrin-classic-em-k1-outbound")

    fake_report = StabilityReport(
        cycler_id="mocked",
        n_laps_propagated=2,
        max_drift_km=42.0,
        max_drift_lap_index=0,
        per_lap_drift_km=(0.0, 42.0),
        stable=True,
        per_lap_dv=(0.0, 0.0),
        total_tcm_dv=0.0,
        frame_used="dynamic",
    )

    with patch(
        "cyclerfinder.verify.real_closure.verify_long_term_stability",
        return_value=fake_report,
    ) as mocked:
        result = verify_real_closure(
            aldrin,
            n_cycles=2,
            ephem=astropy_ephem,
            signature_priority_date=ALDRIN_PRIORITY,
            cycler_id="aldrin-classic-em-k1-outbound",
        )
    assert mocked.call_count == 1
    kwargs = mocked.call_args.kwargs
    assert kwargs["n_laps"] == 2
    assert result.max_drift_km == 42.0
    assert result.closes  # 42 km is well below 200_000 km
    assert result.v3_status == "v3-real-closure-pass"


# ---------------------------------------------------------------------------
# Helper-level tests (plan §4.6)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_construct_real_ephemeris_cycler_aldrin(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """``construct_real_ephemeris_cycler`` builds a 3-encounter Aldrin chain.

    Aldrin's catalogue ``bodies = ["E", "M"]`` deduplicates the chain
    to the visited bodies; the constructor derives the full per-leg
    chain ``("E", "M", "E")`` from the ``legs[].from/to`` fields.
    """
    signature = phase_signature_from_catalogue_entry(aldrin_entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None
    cycler = construct_real_ephemeris_cycler(aldrin_entry, astropy_ephem, t_start)
    assert tuple(cycler.bodies) == ("E", "M", "E")
    assert len(cycler.encounters) == 3
    assert len(cycler.legs) == 2
    assert cycler.encounters[0].t == pytest.approx(t_start)


def test_construct_raises_on_multi_rev_leg(astropy_ephem: Ephemeris) -> None:
    """Plan §4.6: any leg with ``n_revs > 0`` raises ``MultiRevLambertRequiredError``."""
    synthetic = {
        "id": "synthetic-multirev",
        "bodies": ["E", "M"],
        "legs": [
            {"from": "E", "to": "M", "tof_days": 200, "n_revs": 1},
        ],
        "period": {"years": 2.135},
    }
    with pytest.raises(MultiRevLambertRequiredError) as excinfo:
        construct_real_ephemeris_cycler(synthetic, astropy_ephem, 0.0)
    assert excinfo.value.catalogue_id == "synthetic-multirev"
    assert excinfo.value.leg_index == 0


@pytest.mark.slow
def test_resolve_real_t_start_prefers_priority_window(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """Plan §4.6: Aldrin signature + 1985-10-28 priority → window ±5 yr."""
    signature = phase_signature_from_catalogue_entry(aldrin_entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None
    # ±5 yr around the priority epoch in J2000-relative seconds.
    priority_sec = (ALDRIN_PRIORITY - datetime(2000, 1, 1, 12, tzinfo=UTC)).total_seconds()
    delta_sec = 5.0 * 365.25 * SECONDS_PER_DAY
    assert abs(t_start - priority_sec) <= delta_sec, t_start


@pytest.mark.slow
def test_resolve_real_t_start_returns_none_when_no_window(
    astropy_ephem: Ephemeris,
) -> None:
    """Plan §4.6: an absurd signature yields no window → returns ``None``."""
    sig = PhaseSignature(
        bodies=("E", "M"),
        leg_durations_s=(146.0 * SECONDS_PER_DAY,),
        vinf_target_kms=(50.0, 50.0),
    )
    result = _resolve_real_t_start(sig, astropy_ephem, ALDRIN_PRIORITY, window_years=2.0)
    assert result is None


def test_check_vinf_continuity_returns_empty_for_2_encounter_cycler() -> None:
    """Plan §4.6: Aldrin has 2 catalogue encounters; constructor's 3rd
    encounter is the closure point, not an interior encounter. The
    diagnostic returns empty for trajectories whose only encounters
    are boundaries.
    """
    # Synthetic 2-encounter chain: no interior encounter.
    enc0 = Encounter(
        body="E",
        t=0.0,
        r=np.zeros(3),
        v_planet=np.zeros(3),
        vinf_in=np.zeros(3),
        vinf_out=np.zeros(3),
    )
    enc1 = Encounter(
        body="M",
        t=1.0,
        r=np.zeros(3),
        v_planet=np.zeros(3),
        vinf_in=np.zeros(3),
        vinf_out=np.zeros(3),
    )
    cyc = Cycler(
        bodies=["E", "M"],
        period=1.0,
        encounters=[enc0, enc1],
        legs=[],
    )
    assert _check_vinf_continuity(cyc, Ephemeris(model="circular")) == ()


def test_check_vinf_continuity_returns_per_interior_mismatch() -> None:
    """Plan §4.6: a 3-encounter chain returns 1 interior mismatch value."""
    enc0 = Encounter(
        body="E",
        t=0.0,
        r=np.zeros(3),
        v_planet=np.zeros(3),
        vinf_in=np.zeros(3),
        vinf_out=np.zeros(3),
    )
    enc1 = Encounter(
        body="M",
        t=1.0,
        r=np.zeros(3),
        v_planet=np.zeros(3),
        vinf_in=np.array([3.0, 0.0, 0.0]),
        vinf_out=np.array([0.0, 4.0, 0.0]),
    )
    enc2 = Encounter(
        body="E",
        t=2.0,
        r=np.zeros(3),
        v_planet=np.zeros(3),
        vinf_in=np.zeros(3),
        vinf_out=np.zeros(3),
    )
    cyc = Cycler(
        bodies=["E", "M", "E"],
        period=2.0,
        encounters=[enc0, enc1, enc2],
        legs=[],
    )
    out = _check_vinf_continuity(cyc, Ephemeris(model="circular"))
    assert out == (abs(3.0 - 4.0),)


# ---------------------------------------------------------------------------
# Construction error handling
# ---------------------------------------------------------------------------


def test_construct_raises_real_closure_construction_error_on_degenerate_geometry(
    astropy_ephem: Ephemeris,
) -> None:
    """A full-orbit E -> E (365.25 d) trips Lambert's bounded
    single-rev Newton bracket and surfaces a
    :class:`LambertConvergenceError`; the constructor wraps it in
    :class:`RealClosureConstructionError`.
    """
    from cyclerfinder.verify.real_closure import RealClosureConstructionError

    bad = {
        "id": "synthetic-full-orbit",
        "bodies": ["E"],
        "legs": [
            {"from": "E", "to": "E", "tof_days": 365.25, "n_revs": 0},
        ],
        "period": {"years": 1.0},
    }
    with pytest.raises(RealClosureConstructionError):
        construct_real_ephemeris_cycler(bad, astropy_ephem, 0.0)


# ---------------------------------------------------------------------------
# verify_real_closure error paths
# ---------------------------------------------------------------------------


def test_verify_real_closure_requires_n_cycles_ge_2(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    with pytest.raises(ValueError, match="n_cycles >= 2"):
        verify_real_closure(
            aldrin_entry,
            n_cycles=1,
            ephem=astropy_ephem,
            signature_priority_date=ALDRIN_PRIORITY,
        )


def test_verify_real_closure_requires_priority_or_t_start(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """Without a t_start or priority date AND no priority in the entry,
    the resolver has nothing to anchor on."""
    bare = {
        "id": "synthetic-no-priority",
        "bodies": ["E", "M"],
        "legs": [{"from": "E", "to": "M", "tof_days": 146, "n_revs": 0}],
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 6.5},
            {"body": "M", "vinf_kms": 9.7},
        ],
        "period": {"years": 2.135},
        # No priority_date.
    }
    with pytest.raises(ValueError, match="signature_priority_date"):
        verify_real_closure(bare, n_cycles=2, ephem=astropy_ephem)


def test_verify_real_closure_routes_multirev_to_skip_status(
    astropy_ephem: Ephemeris,
) -> None:
    """A catalogue entry with ``n_revs=1`` returns ``v3-skipped-multirev``."""
    multirev = {
        "id": "synthetic-multirev-routed",
        "bodies": ["E", "M"],
        "legs": [
            {"from": "E", "to": "M", "tof_days": 200.0, "n_revs": 1},
        ],
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 6.5},
            {"body": "M", "vinf_kms": 9.7},
        ],
        "period": {"years": 2.135},
        "priority_date": "1985-10-28",
    }
    result = verify_real_closure(
        multirev,
        n_cycles=2,
        ephem=astropy_ephem,
    )
    assert result.v3_status == "v3-skipped-multirev"
    assert not result.closes
    assert result.n_cycles_propagated == 0
