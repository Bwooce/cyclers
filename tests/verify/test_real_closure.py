"""M6b real-ephemeris closure verification gate + helpers.

Spec / plan references
----------------------
* Spec §14 V2-real (multi-cycle real-ephemeris closure gate).
* Plan: ``docs/phases/m6b-real-ephemeris-closure/plan.md`` §4 (tests +
  tolerance), §3.1 (construction path), §3.3 (RealClosureResult shape).

Test layout (plan §4.1):

* Gate 1 — :func:`test_aldrin_cycler_periodic_over_2_cycles_astropy`
  (**aspirational / still-open** DE440 Lambert-chain closure; xfail —
  the binding Aldrin criterion moved to Gate 1b per task #72).
* Gate 1b — :func:`test_aldrin_powered_turn_deficit_gate` (**binding
  Aldrin gate**: sourced anchors + sourced 84°/72° turn deficit).
* Gate 2 — :func:`test_2syn_em_cpom_periodic_over_2_cycles_astropy`
  (**xfail** because the S1L1 entry's direct M->E return leg is undefined
  by construction — S/L are Earth-to-Earth resonant intervals, not
  Earth<->Mars transits; the topology must be re-modelled, not back-filled.
  Multi-rev Lambert is not the blocker).
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
from datetime import UTC, datetime, timedelta
from unittest.mock import patch

import numpy as np
import pytest

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler, Encounter
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv
from cyclerfinder.search.phase_match import (
    PhaseSignature,
    find_candidate_windows,
    leg_duration_seeds,
    phase_signature_from_catalogue_entry,
)
from cyclerfinder.verify.propagate import StabilityReport
from cyclerfinder.verify.real_closure import (
    EXPECTED_SKIPS,
    N_CYCLES_DEFAULT,
    REAL_DRIFT_TOLERANCE_KM,
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


def test_aldrin_ballistic_closure_fails_because_powered(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """NEGATIVE gate: a POWERED cycler must NOT close ballistically.

    The Aldrin Earth flyby needs an ~84° turn but can ballistically deliver
    only ~72° (the deficit asserted in
    :func:`test_aldrin_powered_turn_deficit_gate`). A purely ballistic
    Lambert-chain construction therefore cannot maintain the cycler, and on
    real DE440 ephemeris it drifts by orders of magnitude (~10^8 km over two
    cycles) — not a marginal tolerance miss. This is correct physics, not a
    deficiency: it is the reason the orbit is classified ``powered``.

    The POSITIVE counterpart — the powered cycler built WITH its maintenance
    maneuver by the flyby-constrained periodic solver — is
    :func:`test_aldrin_powered_cycler_solver_and_drift_floor_on_de440`.

    Root cause (diagnostic, 2026-06-01): the published Aldrin M->E return leg
    is the same heliocentric ellipse as the E->M outbound leg, but independent
    per-leg Lambert solves give |V_inf_M_in|=9.78 vs |V_inf_M_out|=10.92 km/s,
    and a free ballistic ellipse (P=2.024 yr) cannot stay in the 2.135 yr
    synodic resonance. Only the gravity assists (preserving |V_inf|, rotating
    it) plus the maintenance maneuver keep it periodic.
    """
    # Premise: this entry IS powered (84° required > 72° achievable at Earth).
    assert aldrin_entry["trajectory_regime"] == "powered"

    result = verify_real_closure(
        aldrin_entry,
        n_cycles=2,
        ephem=astropy_ephem,
        signature_priority_date=ALDRIN_PRIORITY,
        cycler_id="aldrin-classic-em-k1-outbound",
    )
    # Consequence: the ballistic construction does NOT close, and misses by
    # orders of magnitude (> 5x tolerance) — not a near-pass.
    assert result.closes is False
    assert result.max_drift_km > 5 * REAL_DRIFT_TOLERANCE_KM, (
        f"expected gross ballistic non-closure for the powered Aldrin cycler, "
        f"got max_drift_km={result.max_drift_km}"
    )
    assert result.v3_status == "v3-real-closure-fail"


def test_aldrin_powered_cycler_solver_and_drift_floor_on_de440(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """POSITIVE: the flyby-constrained periodic solver produces a genuine
    *powered* Aldrin cycler on real DE440, AND the rotating-frame drift floor
    documents why that cycler needs a per-cycle maintenance maneuver.

    This asserts the physically correct outcome, so it is green (no xfail):

    1. ``solve_powered_periodic_cycler`` yields a cycler with a strictly
       positive maintenance ΔV — the sourced Earth turn deficit (≈84° required
       vs ≈72° achievable at a 200 km flyby). A ballistic, ΔV=0 chain would be
       the *wrong* answer for this orbit (see
       :func:`test_aldrin_ballistic_closure_fails_because_powered`).
    2. Propagated over 2 cycles, the cycler does NOT meet the idealised
       rotating-frame repeat bound (``REAL_DRIFT_TOLERANCE_KM`` = 200,000 km),
       and misses it by orders of magnitude. That bound is a circular-coplanar
       idealisation eccentric Mars cannot satisfy: the propagator pins each
       leg-start to the *real* Mars position at the lap-shifted epoch, and
       Mars's heliocentric radius breathes ≈0.117 AU (≈1.75e7 km) per 2.135 yr
       cycle because the cycler period is not commensurate with Mars's 1.881 yr
       orbit. The maneuver shapes velocity, not where Mars is, so no maintenance
       ΔV closes this metric — which is exactly why a real Aldrin cycler is
       *retargeted* each cycle rather than being geometrically periodic.

    The faithful, *reachable* reproduction target — the SOURCED orbital anchors
    (a ≈ 1.60 AU, e ≈ 0.393) and the SOURCED Earth turn deficit — is asserted
    with teeth by the Phase-C gate below. The ΔV magnitude is unpublished
    (McConaghy 2002 defers it), so it is sanity-bounded here, never matched
    against a sourced number.
    """
    from cyclerfinder.search.bvp import solve_powered_periodic_cycler

    solution = solve_powered_periodic_cycler(
        aldrin_entry,
        ephem=astropy_ephem,
        signature_priority_date=ALDRIN_PRIORITY,
    )
    # The maintenance maneuver must be real (powered cycler), not a ballistic
    # ΔV≈0 neighbour. Sanity-bounded only — the magnitude is our own value.
    assert solution.total_maintenance_dv_kms > 0.0

    result = verify_real_closure(
        solution.cycler,
        n_cycles=2,
        ephem=astropy_ephem,
        t_start=solution.t_start_sec,
        cycler_id="aldrin-classic-em-k1-outbound",
    )
    # Idealised rotating-frame closure is physically unreachable for k=1 Aldrin;
    # the drift floor sits far above tolerance (measured ~4.14e8 km ~ 2072x at
    # n=2/n=3 — #134). Assert the regime qualitatively (>50x) so the gate has
    # teeth without pinning a brittle, solver-derived magnitude.
    assert result.closes is False
    assert result.v3_status == "v3-real-closure-fail"
    assert result.max_drift_km > 50.0 * REAL_DRIFT_TOLERANCE_KM


# ---------------------------------------------------------------------------
# Gate 1b — BINDING Aldrin reproduction gate (task #72, Phase C)
#
# Re-specs the milestone's binding Aldrin criterion away from Lambert-chain
# DE440 drift-closure (architecturally unreachable without the flyby-
# constrained BVP solver — see the xfail above, which stays an open goal,
# not a fabricated pass) and onto what the powered cycler CAN faithfully
# reproduce: the SOURCED orbital anchors and the SOURCED per-orbit turn
# deficit (McConaghy 84 deg required vs 72 deg achievable at a 200 km Earth
# flyby). The ΔV magnitude is unpublished, so it is sanity-bounded only —
# never matched against a sourced number. The gate has teeth: the honest
# bands (±2 deg on the turns, ~1-2 % on a/e) would fail if the optimiser
# slid off the Aldrin family onto a neighbouring ballistic cycler.
# ---------------------------------------------------------------------------

# Published, source-attested anchors (the only legitimate golden targets).
_ALDRIN_PUB_A_AU = 1.60
_ALDRIN_PUB_E = 0.393
_ALDRIN_PUB_VINF_E = 6.5
_ALDRIN_PUB_VINF_M = 9.7
_ALDRIN_PUB_EM_TOF_DAYS = 146.0
_ALDRIN_PUB_TURN_REQ_DEG = 84.0  # McConaghy 2002 Table 4 / dissertation, 1L1
_ALDRIN_PUB_TURN_MAX_DEG = 72.0  # achievable at a 200 km Earth flyby


def test_aldrin_powered_turn_deficit_gate(aldrin_entry: dict[str, object]) -> None:
    """BINDING Aldrin gate (task #72): sourced anchors + sourced turn deficit.

    Runs the periodic maintenance optimiser on the fast circular backend and
    asserts the recovered orbit reproduces the published Aldrin anchors AND the
    published Earth (geocentric) turn deficit (≈84° required vs ≈72° achievable
    → powered). The computed maintenance ΔV is sanity-bounded only. Also checks
    catalogue↔code consistency for the Aldrin outbound entry so the recorded
    regime / turn / altitude metadata cannot silently diverge from the model.
    """
    result = optimise_aldrin_maintenance_dv(Ephemeris("circular"), n_starts=5, seed=0)
    assert result.converged is True

    # --- Sourced orbital anchors (golden).
    assert result.a_au == pytest.approx(_ALDRIN_PUB_A_AU, abs=0.05)
    assert result.e == pytest.approx(_ALDRIN_PUB_E, abs=0.03)
    assert result.leg_tofs_days[0] == pytest.approx(_ALDRIN_PUB_EM_TOF_DAYS, abs=15.0)
    vinf = dict(result.vinf_kms_at_encounters)
    assert vinf["E"] == pytest.approx(_ALDRIN_PUB_VINF_E, abs=0.4)
    assert vinf["M"] == pytest.approx(_ALDRIN_PUB_VINF_M, abs=0.5)

    # --- Sourced per-orbit turn deficit (golden): the powered evidence.
    td = result.turn_deficit
    assert td is not None
    assert td.body == "E"
    assert td.turn_required_deg == pytest.approx(_ALDRIN_PUB_TURN_REQ_DEG, abs=2.0)
    assert td.turn_max_deg == pytest.approx(_ALDRIN_PUB_TURN_MAX_DEG, abs=2.0)
    assert td.deficit_deg > 0.0
    assert td.ballistically_feasible is False

    # --- Computed maintenance ΔV: positive (powered) and sanity-bounded only.
    assert result.maintenance_dv_kms > 0.0
    assert result.maintenance_dv_kms < 3.0

    # --- Catalogue ↔ code consistency (teeth on the recorded metadata).
    assert aldrin_entry["trajectory_regime"] == "powered"
    trajectory = aldrin_entry["trajectory"]
    assert isinstance(trajectory, dict)
    earth_man = next(m for m in trajectory["maneuvers"] if m["body"] == "E")
    assert earth_man["type"] == "flyby-powered"
    assert earth_man["periapsis_alt_km"] == 200
    assert earth_man["turning_angle_deg"] == pytest.approx(_ALDRIN_PUB_TURN_REQ_DEG, abs=2.0)
    assert earth_man["max_turning_angle_deg"] == pytest.approx(_ALDRIN_PUB_TURN_MAX_DEG, abs=2.0)
    # Recorded ΔV is COMPUTED (never golden): only require it present + sane.
    recorded_dv = aldrin_entry["maintenance_dv_kms_per_synodic"]
    assert isinstance(recorded_dv, (int, float))
    assert 0.0 < recorded_dv < 3.0


# ---------------------------------------------------------------------------
# Gate 2 — 2-syn S1L1 (xfail under multi-rev Lambert blocker)
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=False,
    reason=(
        "wrong S1L1 topology, not missing data — multi-rev Lambert is not "
        "the blocker (the solver lands). The catalogue entry posits a direct "
        "Mars->Earth return leg, but no such leg is defined by construction: "
        "in the McConaghy/Longuski/Byrnes nomenclature S and L denote "
        "consecutive Earth-to-Earth resonant intervals, not Earth<->Mars "
        "transits, and Mars is a flyby target of opportunity on the outbound "
        "arc. A single S1L1 vehicle gives one-way fast transit; the crewed "
        "return uses a mirrored conjugate L1S1 cycler. construct_real_"
        "ephemeris_cycler raises on the null return-leg tof_days. Forcing "
        "closure would require fabricating a ToF for a non-existent leg, "
        "which the golden-test discipline forbids. Flips to passing once the "
        "entry is re-modelled as outbound E->M plus the S1/L1 Earth-to-Earth "
        "resonant intervals. "
        "NOTE (2026-06-04): the 5.65/3.05 V∞ anchor is unverified-provenance "
        "(catalogue data_gap vinf_kms_at_encounters, s1l1-2syn-em-cpom): "
        "traces only to spec.md §9; unconfirmed in Patel 2019 / McConaghy 2006 "
        "/ Sanchez Net 2022 — see docs/notes/s1l1-target-topology-mining.md."
    ),
)
def test_2syn_em_cpom_periodic_over_2_cycles_astropy(
    astropy_ephem: Ephemeris,
) -> None:
    """Aspirational gate for the 2-syn S1L1 entry; xfail because the entry's
    direct M->E return leg is undefined by construction (S/L are Earth-to-
    Earth resonant intervals), not because of a missing ToF extraction.
    Additionally: the 5.65/3.05 km/s V∞ anchor is provenance-flagged
    (unverified-provenance data_gap) per docs/notes/s1l1-target-topology-mining.md
    — traces only to spec.md §9, unconfirmed in mined primary sources."""
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


@pytest.mark.parametrize("entry_id", M6B_REGRESSION_IDS)
def test_real_closure_regression_set(
    entry_id: str,
    astropy_ephem: Ephemeris,
) -> None:
    """Plan §4.2: every BALLISTIC regression entry not in EXPECTED_SKIPS closes.

    Powered entries are skipped: ballistic closure is not their target (by
    definition a powered cycler cannot close ballistically — see
    :func:`test_aldrin_ballistic_closure_fails_because_powered` for the
    negative result and
    :func:`test_aldrin_powered_cycler_solver_and_drift_floor_on_de440`
    for the powered-cycler construction + drift floor).
    """
    if entry_id in EXPECTED_SKIPS:
        pytest.skip(EXPECTED_SKIPS[entry_id])
    entries = load_m6b_entries()
    entry = next(e for e in entries if e["id"] == entry_id)
    if entry.get("trajectory_regime") == "powered":
        pytest.skip(
            f"{entry_id}: powered cycler — ballistic closure not expected; "
            "covered by test_aldrin_ballistic_closure_fails_because_powered "
            "(negative) + test_aldrin_powered_cycler_solver_and_drift_floor_on_de440 "
            "(positive: the solver builds the powered cycler with a real "
            "maintenance ΔV, and the idealised rotating-frame drift bound is "
            "physically unreachable for k=1 Aldrin — asserted as the drift "
            "floor)."
        )
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


def test_construct_builds_multi_rev_leg(astropy_ephem: Ephemeris) -> None:
    """A leg with ``n_revs=1`` now builds a multi-rev Leg instead of raising.

    Uses a 780 d Earth->Mars leg: t_min(1) for that geometry is ~630 d, so
    revolution 1 is feasible. The catalogue carries no ``branch`` field, so
    construction defaults to the ``low`` branch.
    """
    synthetic = {
        "id": "synthetic-multirev",
        "bodies": ["E", "M"],
        "legs": [
            {"from": "E", "to": "M", "tof_days": 780.0, "n_revs": 1},
        ],
        "period": {"years": 2.135},
    }
    cyc = construct_real_ephemeris_cycler(synthetic, astropy_ephem, 0.0)
    assert cyc.legs[0].n_revs == 1
    assert cyc.legs[0].branch == "low"


def test_resolve_real_t_start_picks_low_mismatch_window(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> None:
    """Plan §4.6 (STAGE 3): the robust resolver ranks candidate windows by
    V_inf mismatch, NOT calendar proximity (the proximity tie-break was the
    degenerate-basin bug).

    Two honest checks, both `# COMPUTED` (no sourced V_inf is re-asserted):

    1. The resolver returns the single lowest-mismatch window across the full
       ±10 yr range — i.e. its epoch equals the lowest-mismatch candidate, not
       the one nearest the priority date.
    2. That window's summed |V_inf| mismatch beats the 20 km/s resolver cap and
       lands inside the search range.
    """
    signature = phase_signature_from_catalogue_entry(aldrin_entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None

    # Independently reconstruct the candidate pool the resolver ranks.
    delta = timedelta(days=10.0 * 365.25)
    seeds = leg_duration_seeds(
        bodies=signature.bodies,
        primary_leg_durations_s=signature.leg_durations_s,
        vinf_target_kms=signature.vinf_target_kms,
        period_s=sum(signature.leg_durations_s),
    )
    windows = find_candidate_windows(
        seeds,
        astropy_ephem,
        (ALDRIN_PRIORITY - delta, ALDRIN_PRIORITY + delta),
        n=5 * len(seeds),
        mismatch_cap_kms=20.0,
    )
    assert windows, "expected at least one window within the resolver range"
    # COMPUTED: the resolver picks the lowest-mismatch candidate (windows[0]).
    best = windows[0]
    resolved_sec = (best.departure_date - datetime(2000, 1, 1, 12, tzinfo=UTC)).total_seconds()
    assert t_start == pytest.approx(resolved_sec)
    # COMPUTED: that window beats the resolver cap and is inside the range.
    assert best.mismatch_kms < 20.0
    assert abs(t_start - resolved_sec) <= delta.total_seconds()


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


def test_verify_real_closure_constructs_feasible_multirev(
    astropy_ephem: Ephemeris,
) -> None:
    """A feasible ``n_revs=1`` entry is built and measured, never skip-routed.

    With multi-rev Lambert landed there is no ``v3-skipped-multirev`` status;
    the 780 d Earth->Mars leg clears t_min(1) (~630 d) so construction
    succeeds and the result carries a real closure outcome.
    """
    multirev = {
        "id": "synthetic-multirev-routed",
        "bodies": ["E", "M"],
        "legs": [
            {"from": "E", "to": "M", "tof_days": 780.0, "n_revs": 1},
        ],
        "period": {"years": 2.135},
    }
    result = verify_real_closure(
        multirev,
        n_cycles=2,
        ephem=astropy_ephem,
        t_start=0.0,
    )
    assert result.v3_status != "v3-skipped-multirev"
    assert result.v3_status in (
        "v3-real-closure-pass",
        "v3-real-closure-fail",
        "v3-construction-error",
    )
