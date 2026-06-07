"""Tests for :mod:`cyclerfinder.search.maintain` — the Aldrin maintenance-ΔV
periodic optimiser.

These tests run on the fast ``Ephemeris("circular")`` backend so CI stays
quick; an ``astropy`` real-DE440 cross-check would be marked ``@pytest.mark.slow``
if added.

Golden assertions (items 2, 5) target ONLY the published, source-attested
Aldrin anchors (Rogers et al. 2012 Table 1; Russell 2004 Table 3.4; McConaghy /
Longuski / Byrnes 2002 AIAA 2002-4420 Table 4 row "1L1"):

* a = 1.60 AU, e = 0.393
* V∞_Earth = 6.5 km/s, V∞_Mars = 9.7-9.75 km/s
* Earth→Mars ToF = 146 days
* Earth flyby turn: 84° required vs 72° achievable (the "powered" deficit)

The *turn angles* are source-traceable, so the test asserts the computed Earth
turn against McConaghy's published 84° / 72° (item 5). The computed maintenance
ΔV in km/s, by contrast, has NO published counterpart (McConaghy 2002 defers
it), so it is REPORTED / sanity-bounded only (item 3) — never an exact match
target. No assertion anywhere relies on a ballistic ΔV == 0.

Tolerances are honest engineering bands (~1-2 % on a/e, a few tenths km/s on
V∞, ±15 d on ToF, ±2° on the turn angles), NOT widened to force a pass.
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.maintain import (
    MaintenanceOptimResult,
    _build_chain,
    _default_t0_guess,
    _maintenance_dv_chain,
    idealized_flyby_turn_deficit,
    optimise_aldrin_maintenance_dv,
    optimise_maintenance_dv,
)

# Published, source-attested anchors — the ONLY legitimate assertion targets.
_PUB_A_AU = 1.60
_PUB_E = 0.393
_PUB_VINF_E_KMS = 6.5
_PUB_VINF_M_KMS = 9.7
_PUB_EM_TOF_DAYS = 146.0
# McConaghy/Longuski/Byrnes 2002 Table 4 "1L1": the Earth flyby needs an 84°
# turn but can ballistically deliver only 72° — the reason Aldrin is powered.
_PUB_EARTH_TURN_REQ_DEG = 84.0
_PUB_EARTH_TURN_MAX_DEG = 72.0

# Honest tolerances (per task spec).
_TOL_A = 0.05
_TOL_E = 0.03
_TOL_VINF_E = 0.4
_TOL_VINF_M = 0.5
_TOL_TOF = 15.0
_TOL_TURN_DEG = 2.0


@pytest.fixture(scope="module")
def aldrin_result() -> MaintenanceOptimResult:
    """Optimise the Aldrin maintenance ΔV once on the circular backend."""
    return optimise_aldrin_maintenance_dv(Ephemeris("circular"), n_starts=5, seed=0)


def _vinf(result: MaintenanceOptimResult, body: str) -> float:
    """First ``|V∞|`` (km/s) for ``body`` across the encounter list."""
    for code, vinf in result.vinf_kms_at_encounters:
        if code == body:
            return vinf
    raise AssertionError(f"no encounter for body {body!r}")


# ---------------------------------------------------------------------------
# Item 1: the optimiser converges.
# ---------------------------------------------------------------------------


def test_optimiser_converges(aldrin_result: MaintenanceOptimResult) -> None:
    assert aldrin_result.converged is True


def test_result_is_frozen() -> None:
    """The result dataclass is immutable (frozen)."""
    from dataclasses import FrozenInstanceError

    r = optimise_aldrin_maintenance_dv(Ephemeris("circular"), n_starts=2, seed=0)
    with pytest.raises(FrozenInstanceError):
        r.a_au = 0.0  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Item 2: GOLDEN assertions against the published anchors only.
# ---------------------------------------------------------------------------


def test_semi_major_axis_matches_published(aldrin_result: MaintenanceOptimResult) -> None:
    assert aldrin_result.a_au == pytest.approx(_PUB_A_AU, abs=_TOL_A), (
        f"computed a={aldrin_result.a_au:.4f} AU vs published {_PUB_A_AU} AU"
    )


def test_eccentricity_matches_published(aldrin_result: MaintenanceOptimResult) -> None:
    assert aldrin_result.e == pytest.approx(_PUB_E, abs=_TOL_E), (
        f"computed e={aldrin_result.e:.4f} vs published {_PUB_E}"
    )


def test_vinf_earth_matches_published(aldrin_result: MaintenanceOptimResult) -> None:
    vinf_e = _vinf(aldrin_result, "E")
    assert vinf_e == pytest.approx(_PUB_VINF_E_KMS, abs=_TOL_VINF_E), (
        f"computed V∞_Earth={vinf_e:.3f} km/s vs published {_PUB_VINF_E_KMS} km/s"
    )


def test_vinf_mars_matches_published(aldrin_result: MaintenanceOptimResult) -> None:
    vinf_m = _vinf(aldrin_result, "M")
    assert vinf_m == pytest.approx(_PUB_VINF_M_KMS, abs=_TOL_VINF_M), (
        f"computed V∞_Mars={vinf_m:.3f} km/s vs published {_PUB_VINF_M_KMS} km/s"
    )


def test_em_tof_matches_published(aldrin_result: MaintenanceOptimResult) -> None:
    em_tof = aldrin_result.leg_tofs_days[0]
    assert em_tof == pytest.approx(_PUB_EM_TOF_DAYS, abs=_TOL_TOF), (
        f"computed E→M ToF={em_tof:.2f} d vs published {_PUB_EM_TOF_DAYS} d"
    )


# ---------------------------------------------------------------------------
# Item 3: REPORT-ONLY sanity bounds on the computed maintenance ΔV (km/s).
# This value is OUR computation (no published counterpart exists — McConaghy
# 2002 defers it), so it is sanity-bounded only and is NEVER matched against a
# sourced target. It must be strictly positive: Aldrin is powered.
# ---------------------------------------------------------------------------


def test_maintenance_dv_is_sane(aldrin_result: MaintenanceOptimResult) -> None:
    dv = aldrin_result.maintenance_dv_kms
    import math

    assert math.isfinite(dv), f"maintenance ΔV not finite: {dv}"
    # Strictly positive — Aldrin's Earth flyby cannot close ballistically. A
    # zero here would mean the optimiser slid onto a neighbouring ballistic
    # family, which is exactly the artifact this redesign removes.
    assert dv > 0.0, f"powered Aldrin must have positive maintenance ΔV, got {dv}"
    # Sanity bound only: an Aldrin maintenance ΔV above a Hohmann-like ~3 km/s
    # would indicate nonsense, not the cycler. NOT a match against a published
    # value.
    assert dv < 3.0, f"maintenance ΔV implausibly large: {dv} km/s"

    # The per-encounter breakdown must sum to the reported total.
    breakdown_sum = sum(v for _b, v in aldrin_result.per_encounter_dv_kms)
    assert breakdown_sum == pytest.approx(dv, abs=1.0e-9)


# ---------------------------------------------------------------------------
# Item 5: GOLDEN turn-angle assertion. The Earth flyby turn is a geometric
# consequence of the sourced (a, e), so 84° required / 72° achievable IS a
# source-traceable target (McConaghy 2002 Table 4 "1L1"). This replaces any
# reliance on a ballistic ΔV == 0 match.
# ---------------------------------------------------------------------------


def test_earth_turn_matches_published(aldrin_result: MaintenanceOptimResult) -> None:
    td = aldrin_result.turn_deficit
    assert td is not None, "expected a turn deficit at the Earth return flyby"
    assert td.body == "E"
    assert td.turn_required_deg == pytest.approx(_PUB_EARTH_TURN_REQ_DEG, abs=_TOL_TURN_DEG), (
        f"computed Earth turn required={td.turn_required_deg:.2f}° vs published "
        f"{_PUB_EARTH_TURN_REQ_DEG}°"
    )
    assert td.turn_max_deg == pytest.approx(_PUB_EARTH_TURN_MAX_DEG, abs=_TOL_TURN_DEG), (
        f"computed Earth turn max={td.turn_max_deg:.2f}° vs published {_PUB_EARTH_TURN_MAX_DEG}°"
    )
    # The deficit (required > achievable) is the reason Aldrin is powered.
    assert td.deficit_deg > 0.0
    assert td.ballistically_feasible is False


def test_idealized_turn_deficit_from_published_anchors() -> None:
    """Fed the published anchors directly (no optimiser), the geometric turn
    deficit reproduces McConaghy's 84° / 72° — confirming the number is a
    property of the sourced orbit, not of our search."""
    # The sourced 200 km Earth flyby altitude (McConaghy's dissertation bases
    # the achievable geocentric turn on it) — not the conservative 300 km
    # default — so the achievable turn reproduces the published ≈72°.
    td = idealized_flyby_turn_deficit(_PUB_A_AU, _PUB_E, "E", flyby_alt_km=200.0)
    assert td is not None
    assert td.turn_required_deg == pytest.approx(_PUB_EARTH_TURN_REQ_DEG, abs=_TOL_TURN_DEG)
    assert td.turn_max_deg == pytest.approx(_PUB_EARTH_TURN_MAX_DEG, abs=_TOL_TURN_DEG)
    assert td.vinf_kms == pytest.approx(_PUB_VINF_E_KMS, abs=_TOL_VINF_E)
    assert td.deficit_deg > 0.0
    assert td.ballistically_feasible is False


def test_turn_deficit_none_when_orbit_misses_body() -> None:
    """An orbit that never reaches Mars' radius yields no Mars turn deficit."""
    # a=1.0, e=0.1 → apoapsis 1.1 AU, well inside Mars' 1.52 AU orbit.
    assert idealized_flyby_turn_deficit(1.0, 0.1, "M") is None


# ---------------------------------------------------------------------------
# Item 4: a deliberately detuned ToF guess (far from 146 d), with the Aldrin
# launch phase held, still converges back toward the published anchors —
# proving the optimiser finds a real minimum rather than accepting anything.
# ---------------------------------------------------------------------------


def test_detuned_tof_guess_converges_back_to_anchors() -> None:
    t0 = _default_t0_guess(_PUB_EM_TOF_DAYS)  # Aldrin launch phase
    detuned = optimise_aldrin_maintenance_dv(
        Ephemeris("circular"),
        t0_guess_sec=t0,
        em_tof_days_guess=210.0,  # far from the 146 d anchor
        me_tof_days_guess=520.0,
        n_starts=5,
        seed=7,
    )
    assert detuned.converged is True
    # It must recover the published anchors within the same honest bands,
    # not park at the detuned starting guess (210 d).
    assert detuned.leg_tofs_days[0] == pytest.approx(_PUB_EM_TOF_DAYS, abs=_TOL_TOF), (
        f"detuned start did not return to the anchor: ToF={detuned.leg_tofs_days[0]:.2f} d"
    )
    assert detuned.a_au == pytest.approx(_PUB_A_AU, abs=_TOL_A)
    assert _vinf(detuned, "E") == pytest.approx(_PUB_VINF_E_KMS, abs=_TOL_VINF_E)
    assert _vinf(detuned, "M") == pytest.approx(_PUB_VINF_M_KMS, abs=_TOL_VINF_M)


# ---------------------------------------------------------------------------
# Structural generality (#75): the optimiser is body-/length-agnostic. These
# tests assert NO sourced number — only that the Aldrin wrapper is a faithful
# pass-through of the generalised core, and that the machinery accepts a
# sequence longer than the hardcoded E-M-E (forward map + objective handle an
# arbitrary closed chain). The Aldrin physics is validated by the tests above.
# ---------------------------------------------------------------------------


def test_aldrin_wrapper_matches_generalized_core() -> None:
    """``optimise_aldrin_maintenance_dv`` is a bit-for-bit pass-through of
    ``optimise_maintenance_dv`` with the Aldrin parameters pinned."""
    t0 = _default_t0_guess(_PUB_EM_TOF_DAYS)
    wrapped = optimise_aldrin_maintenance_dv(
        Ephemeris("circular"), t0_guess_sec=t0, n_starts=5, seed=0
    )
    from cyclerfinder.search.construct import build_aldrin_seed
    from cyclerfinder.search.maintain import (
        _ALDRIN_EARTH_FLYBY_ALT_KM,
        _TOF1_BOUNDS_DAYS,
        _TOF2_BOUNDS_DAYS,
    )

    generalized = optimise_maintenance_dv(
        ["E", "M", "E"],
        Ephemeris("circular"),
        t0_guess_sec=t0,
        tof_days_guesses=(146.0, 634.0),
        tof_bounds_days=(_TOF1_BOUNDS_DAYS, _TOF2_BOUNDS_DAYS),
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        tof_jitter_half_days=(20.0, 60.0),
        n_starts=5,
        seed=0,
        seed_cycler_factory=lambda: build_aldrin_seed(Ephemeris("circular")),
    )
    assert generalized.t0_sec == wrapped.t0_sec
    assert generalized.leg_tofs_days == wrapped.leg_tofs_days
    assert generalized.a_au == wrapped.a_au
    assert generalized.e == wrapped.e
    assert generalized.maintenance_dv_kms == wrapped.maintenance_dv_kms
    assert generalized.converged == wrapped.converged


def test_chain_forward_map_handles_four_encounters() -> None:
    """The forward map / ΔV accounting accept a sequence longer than E-M-E.

    Single-rev Lambert builds each leg, so a closed four-encounter chain
    constructs and yields a finite, non-negative maintenance ΔV. No sourced
    quantity is asserted — this only proves the machinery is length-agnostic.
    """
    import numpy as np

    sequence = ["E", "M", "V", "E"]
    # t0 = 0, three positive leg ToFs (days); arbitrary but non-degenerate.
    x = np.array([0.0, 200.0, 300.0, 250.0], dtype=np.float64)
    cycler = _build_chain(x, sequence, Ephemeris("circular"))
    assert cycler is not None, "single-rev Lambert should build a 4-encounter chain"
    assert [enc.body for enc in cycler.encounters] == sequence
    assert len(cycler.legs) == 3
    # Encounter epochs are the running cumulative sum of the leg ToFs.
    assert cycler.encounters[1].t == pytest.approx(200.0 * SECONDS_PER_DAY)
    assert cycler.encounters[2].t == pytest.approx(500.0 * SECONDS_PER_DAY)
    dv = _maintenance_dv_chain(cycler)
    import math

    assert math.isfinite(dv) and dv >= 0.0


def test_generalized_optimiser_runs_on_four_encounter_sequence() -> None:
    """End-to-end: the optimiser drives an arbitrary closed four-encounter
    sequence and returns a well-shaped, finite result (length-agnostic)."""
    import math

    import numpy as np

    sequence = ["E", "M", "V", "E"]
    seed = _build_chain(
        np.array([0.0, 200.0, 300.0, 250.0], dtype=np.float64),
        sequence,
        Ephemeris("circular"),
    )
    assert seed is not None
    result = optimise_maintenance_dv(
        sequence,
        Ephemeris("circular"),
        t0_guess_sec=0.0,
        tof_days_guesses=(200.0, 300.0, 250.0),
        tof_bounds_days=((120.0, 300.0), (180.0, 420.0), (150.0, 360.0)),
        n_starts=2,
        seed=0,
        seed_cycler_factory=lambda: seed,
    )
    assert len(result.leg_tofs_days) == 3
    assert math.isfinite(result.maintenance_dv_kms)
    assert result.maintenance_dv_kms >= 0.0
    assert [b for b, _ in result.vinf_kms_at_encounters] == sequence


# ---------------------------------------------------------------------------
# STAGE 1 — multi-rev / multi-encounter parameter plumbing
# ---------------------------------------------------------------------------


def test_build_chain_passes_revs_and_branch() -> None:
    """STAGE 1: ``_build_chain`` accepts explicit per-leg revs / branch and
    still returns a built ``Cycler`` without raising.

    Provenance: ``# COMPUTED`` — circular ephemeris; no velocity / V_inf golden
    is asserted, only that the multi-rev plumbing constructs a chain.
    """
    import numpy as np

    sequence = ["E", "M", "E"]
    # Direct outbound, then a long return leg with room for a 1-rev solution.
    x = np.array([0.0, 200.0, 900.0], dtype=np.float64)
    cycler = _build_chain(
        x,
        sequence,
        Ephemeris("circular"),
        per_leg_revs=(0, 1),
        per_leg_branch=("single", "low"),
    )
    assert cycler is not None
    assert [enc.body for enc in cycler.encounters] == sequence
    assert len(cycler.legs) == 2
    # The requested return-leg revolution count is honoured by construct_cycler.
    assert cycler.legs[1].n_revs == 1


def test_build_chain_none_defaults_match_single_rev() -> None:
    """STAGE 1: omitting per_leg_revs/per_leg_branch is byte-identical to the
    previous all-single-rev behaviour (Aldrin bit-reproducibility guard).

    Provenance: ``# COMPUTED`` — equality between two computed cyclers.
    """
    import numpy as np

    sequence = ["E", "M", "E"]
    x = np.array([0.0, 146.0, 634.0], dtype=np.float64)
    default = _build_chain(x, sequence, Ephemeris("circular"))
    explicit = _build_chain(
        x,
        sequence,
        Ephemeris("circular"),
        per_leg_revs=None,
        per_leg_branch=None,
    )
    assert default is not None and explicit is not None
    for ld, le in zip(default.legs, explicit.legs, strict=True):
        assert np.array_equal(ld.v_depart, le.v_depart)
        assert np.array_equal(ld.v_arrive, le.v_arrive)
        assert ld.n_revs == le.n_revs


def test_optimise_maintenance_dv_multirev_param_threaded() -> None:
    """STAGE 1: ``optimise_maintenance_dv`` threads explicit per-leg revs /
    branch through the whole solve and converges to a finite maintenance ΔV.

    Provenance: ``# COMPUTED`` — circular ephemeris; never asserts a sourced
    V_inf target, only convergence and finiteness of the computed cost.
    """
    import math

    import numpy as np

    sequence = ["E", "M", "E", "E"]
    # The middle E->E leg carries the 1-rev branch; a 1-rev Lambert leg needs a
    # ToF above its physical minimum (~700 d on the circular E-orbit), so seed
    # it there rather than below the floor.
    seed = _build_chain(
        np.array([0.0, 200.0, 760.0, 800.0], dtype=np.float64),
        sequence,
        Ephemeris("circular"),
        per_leg_revs=(0, 1, 0),
        per_leg_branch=("single", "low", "single"),
    )
    assert seed is not None
    result = optimise_maintenance_dv(
        sequence,
        Ephemeris("circular"),
        t0_guess_sec=0.0,
        tof_days_guesses=(200.0, 760.0, 800.0),
        tof_bounds_days=((150.0, 300.0), (700.0, 1100.0), (600.0, 1000.0)),
        per_leg_revs=(0, 1, 0),
        per_leg_branch=("single", "low", "single"),
        n_starts=2,
        seed=0,
        seed_cycler_factory=lambda: seed,
    )
    assert result.converged is True
    assert math.isfinite(result.maintenance_dv_kms)
    assert result.maintenance_dv_kms >= 0.0
    # The 1-rev return leg survives the solve.
    assert result.cycler.legs[1].n_revs == 1


# ---------------------------------------------------------------------------
# SLOW DE440 gate (finding #114): the real-ephemeris Aldrin solve must land
# IN-FAMILY when seeded from the sourced-anchor real launch window, not on the
# off-family degenerate basin (a ≈ 0.95 AU, e ≈ 0.99, V∞ ≈ 38 km/s, ΔV ≈ 55
# km/s) that the bare circular phase seed slides onto on real DE440.
#
# The seed is the V∞-anchor real launch-window resolver (the same one M6b's
# real-closure / BVP path uses), passed via ``real_window_priority_date``. This
# is the value that goes public on cyclers.space, so it is gated with teeth:
# a/e/V∞ near the SOURCED anchors AND the plausibility bar (ΔV < 3.0). The ΔV
# magnitude itself is OUR computation (unpublished) so it is bounded, not
# matched. Cross-checked: the independent ``solve_powered_periodic_cycler`` BVP
# path returns the same ΔV (2.9138 km/s) on this window.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_aldrin_de440_seeded_solve_is_in_family() -> None:
    from datetime import UTC, datetime

    result = optimise_aldrin_maintenance_dv(
        Ephemeris("astropy"),
        n_starts=5,
        seed=0,
        real_window_priority_date=datetime(1985, 10, 28, tzinfo=UTC),
    )

    assert result.converged is True, "DE440 seeded solve did not converge"

    # In-family on the SOURCED orbital anchors (same honest bands as the
    # circular surface above).
    assert result.a_au == pytest.approx(_PUB_A_AU, abs=_TOL_A), (
        f"DE440 a={result.a_au:.4f} AU off-family vs published {_PUB_A_AU} AU"
    )
    assert result.e == pytest.approx(_PUB_E, abs=_TOL_E), (
        f"DE440 e={result.e:.4f} off-family vs published {_PUB_E}"
    )
    vinf_e = _vinf(result, "E")
    vinf_m = _vinf(result, "M")
    assert vinf_e == pytest.approx(_PUB_VINF_E_KMS, abs=_TOL_VINF_E), (
        f"DE440 V∞_Earth={vinf_e:.3f} km/s off-family vs published {_PUB_VINF_E_KMS} km/s"
    )
    assert vinf_m == pytest.approx(_PUB_VINF_M_KMS, abs=_TOL_VINF_M), (
        f"DE440 V∞_Mars={vinf_m:.3f} km/s off-family vs published {_PUB_VINF_M_KMS} km/s"
    )

    # Plausibility bar (finding #114): the off-family basin reported ΔV ≈ 55
    # km/s. The in-family value must be strictly positive (Aldrin is powered)
    # and below the Hohmann-like sanity ceiling. NOT a match against a published
    # number — the ΔV magnitude is unpublished.
    dv = result.maintenance_dv_kms
    assert dv > 0.0, f"in-family Aldrin must have positive maintenance ΔV, got {dv}"
    assert dv < 3.0, f"DE440 maintenance ΔV implausibly large (off-family?): {dv} km/s"


# ---------------------------------------------------------------------------
# Oberth periapsis cost model (task #151) — opt-in, default unchanged.
#
# The default model ("asymptote") MUST leave the reported ΔV bit-identical;
# "oberth_periapsis" re-costs the SAME recovered cycler (anchors unchanged) with
# the Oberth credit. The ΔV magnitudes are OUR computation (DIAGNOSTIC), so the
# circular tests assert ORDERING / IDENTITY, not a sourced target.
# ---------------------------------------------------------------------------


def test_turn_deficit_default_model_is_asymptote() -> None:
    # Default flyby_cost_model reproduces the asymptote-rotation ΔV exactly.
    default = idealized_flyby_turn_deficit(_PUB_A_AU, _PUB_E, "E", flyby_alt_km=200.0)
    asym = idealized_flyby_turn_deficit(
        _PUB_A_AU, _PUB_E, "E", flyby_alt_km=200.0, flyby_cost_model="asymptote"
    )
    assert default is not None and asym is not None
    assert default.dv_kms == asym.dv_kms


def test_turn_deficit_oberth_cheaper_than_asymptote() -> None:
    asym = idealized_flyby_turn_deficit(
        _PUB_A_AU, _PUB_E, "E", flyby_alt_km=200.0, flyby_cost_model="asymptote"
    )
    oberth = idealized_flyby_turn_deficit(
        _PUB_A_AU, _PUB_E, "E", flyby_alt_km=200.0, flyby_cost_model="oberth_periapsis"
    )
    assert asym is not None and oberth is not None
    # Same geometry -> identical turn angles; only the charged ΔV differs.
    assert oberth.turn_required_deg == asym.turn_required_deg
    assert oberth.turn_max_deg == asym.turn_max_deg
    # Aldrin sits in the deep-well regime: Oberth credit is real.
    assert 0.0 < oberth.dv_kms < asym.dv_kms


def test_aldrin_wrapper_default_unchanged_oberth_cheaper() -> None:
    # The wrapper's default must reproduce the historic result; the opt-in model
    # only changes the reported ΔV, never the recovered anchors.
    base = optimise_aldrin_maintenance_dv(Ephemeris("circular"), n_starts=5, seed=0)
    oberth = optimise_aldrin_maintenance_dv(
        Ephemeris("circular"), n_starts=5, seed=0, flyby_cost_model="oberth_periapsis"
    )
    assert oberth.a_au == base.a_au
    assert oberth.e == base.e
    assert oberth.t0_sec == base.t0_sec
    assert oberth.leg_tofs_days == base.leg_tofs_days
    assert 0.0 < oberth.maintenance_dv_kms < base.maintenance_dv_kms


@pytest.mark.slow
def test_aldrin_de440_recost_under_both_models() -> None:
    """Re-cost the in-family DE440 Aldrin schedule under both models.

    Regression anchor: the asymptote model reproduces ≈2.9138 km/s (the value
    that cross-checks the BVP path). The Oberth value is OUR computation
    (DIAGNOSTIC/PROVISIONAL) — bounded and ordered, not matched.
    """
    from datetime import UTC, datetime

    d = datetime(1985, 10, 28, tzinfo=UTC)
    old = optimise_aldrin_maintenance_dv(
        Ephemeris("astropy"), n_starts=5, seed=0, real_window_priority_date=d
    )
    new = optimise_aldrin_maintenance_dv(
        Ephemeris("astropy"),
        n_starts=5,
        seed=0,
        real_window_priority_date=d,
        flyby_cost_model="oberth_periapsis",
    )
    # Anchors identical across models (objective unchanged).
    assert new.a_au == old.a_au
    assert new.e == old.e
    assert new.t0_sec == old.t0_sec
    # Asymptote regression anchor.
    assert old.maintenance_dv_kms == pytest.approx(2.9138, abs=1e-3)
    # Oberth re-cost: strictly cheaper, still strictly positive (powered cycler).
    assert 0.0 < new.maintenance_dv_kms < old.maintenance_dv_kms
    assert new.maintenance_dv_kms == pytest.approx(1.9336, abs=1e-3)
    # Both pass the publication plausibility bar.
    from cyclerfinder.verify.plausibility import QuantityKind, check_publishable

    assert check_publishable(QuantityKind.MAINTENANCE_DV_KMS, new.maintenance_dv_kms).ok
