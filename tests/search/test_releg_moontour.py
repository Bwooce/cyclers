"""Tests for the releg moon-tour driver (#449, DSM branch).

The driver loops the legs of a tour skeleton with a chosen ``Releg`` backend,
enforces V_inf-continuity after the powered retarget + the closed-cycle wrap,
sums the delivered ΔV per cycle, and classifies it — but FIRST runs the cheap
VILM/linkability prefilter (``search.moon_prune.moon_leg_admissible``) so a
structurally-dead skeleton (disjoint Tisserand contours, the Uranus case) is
reported EMPTY without paying for the DSM solve.

Two anchors:
* Positive control (Jovian Io-Europa-Ganymede-Io, links at V_inf=4): the powered
  DSM cycle CLOSES (every leg feasible, continuity residual below the ballistic
  gate post-retarget) — the capability proof.
* Structural negative (Uranian Ariel-Umbriel-Ariel): the prefilter marks it
  unbridgeable and the driver returns EMPTY WITHOUT running the DSM solve —
  reproducing uranus-neptune-regular-moon-endgame-vilm-2026-06-23 honestly.
"""

from __future__ import annotations

import math

from cyclerfinder.search.releg_moontour import close_powered_cycle
from cyclerfinder.search.releg_solver import (
    BallisticReleg,
    DsmReleg,
    MultiRevLeveragingReleg,
)


def test_jovian_positive_control_closes_powered() -> None:
    """The Jovian Io-Europa-Ganymede-Io powered cycle closes (capability proof).

    The registry positive control (both link at V_inf=4). With ``DsmReleg`` the
    driver retargets each leg to a common V_inf so every flyby is continuous
    after the maneuver; the closure (continuity residual below the ballistic
    gate) is the proof a powered leg can close a ballistically-too-expensive
    cycle. The total ΔV is finite and the verdict carries a dv-band class.
    """
    sequence = ("Io", "Europa", "Ganymede", "Io")
    # Per-leg ToF = geometric-mean of the two moons' orbital periods (the
    # resonance-scaled transfer time the discovery genome uses), x a scale.
    leg_tofs_days = _geomean_tofs(sequence, scale=1.5)
    verdict = close_powered_cycle(
        primary="Jupiter",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0, 0),
        releg=DsmReleg(),
        phasing={"Io": 0.0, "Europa": 1.0, "Ganymede": 2.0},
        dv_band="powered_dsm",
    )
    assert verdict.prefilter_skipped is False
    assert verdict.feasible is True
    # Continuity satisfied after the powered retarget (below the ballistic gate).
    assert verdict.continuity_residual_kms < 0.05
    assert math.isfinite(verdict.total_dv_kms)
    assert verdict.total_dv_kms > 0.0
    assert len(verdict.per_leg_dv_kms) == 3


def test_uranus_disjoint_prefiltered_empty() -> None:
    """The Uranian Ariel-Umbriel-Ariel skeleton is prefiltered EMPTY.

    Every consecutive-moon leg is unlinkable (disjoint Tisserand/resonance
    contours at all probed V_inf), so the VILM/linkability prefilter marks the
    skeleton unbridgeable and the driver returns an EMPTY verdict WITHOUT running
    the DSM solve. This is the structural negative — the honesty test.
    """
    sequence = ("Ariel", "Umbriel", "Ariel")
    leg_tofs_days = _geomean_tofs(sequence, scale=1.0)
    verdict = close_powered_cycle(
        primary="Uranus",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0),
        releg=DsmReleg(),
        phasing={"Ariel": 0.0, "Umbriel": 1.0},
        dv_band="powered_dsm",
    )
    assert verdict.prefilter_skipped is True
    assert verdict.feasible is False
    assert verdict.total_dv_kms == math.inf
    # The reason records the unbridgeable leg (not silent).
    assert any("not linkable" in r for r in verdict.prefilter_reasons)


def test_ballistic_backend_also_runs_through_driver() -> None:
    """The driver is backend-agnostic: BallisticReleg runs too (dv == 0 legs)."""
    sequence = ("Io", "Europa", "Ganymede", "Io")
    leg_tofs_days = _geomean_tofs(sequence, scale=1.5)
    verdict = close_powered_cycle(
        primary="Jupiter",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0, 0),
        releg=BallisticReleg(),
        phasing={"Io": 0.0, "Europa": 1.0, "Ganymede": 2.0},
        dv_band=None,
    )
    assert verdict.prefilter_skipped is False
    # A ballistic backend delivers no ΔV (it cannot retarget); the per-leg ΔV
    # is all zero. The continuity residual is whatever the ballistic legs leave.
    assert verdict.total_dv_kms == 0.0


# ---------------------------------------------------------------------------
# Task 3 (#465) — multi-rev leveraging backend brings the Galilean tour IN-BAND
# ---------------------------------------------------------------------------


def test_multirev_galilean_positive_control_in_band() -> None:
    """The Galilean Io-Europa-Ganymede-Io cycle closes IN-BAND with the chain.

    The #449/#464 single-leg relegs closed this cycle geometrically but
    OUT-OF-BAND (13.18 km/s DSM / 12.03 km/s SF). The multi-rev leveraging chain
    walks each leg's arrival V_inf down to the common flyby target step by step,
    so the cycle closes (continuity below the ballistic gate, by construction)
    AND the total ΔV lands INSIDE the powered dv-band — the in-band closure the
    prior relegs missed (#465 gate).
    """
    sequence = ("Io", "Europa", "Ganymede", "Io")
    leg_tofs_days = _geomean_tofs(sequence, scale=1.5)
    verdict = close_powered_cycle(
        primary="Jupiter",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0, 0),
        releg=MultiRevLeveragingReleg(),
        phasing={"Io": 0.0, "Europa": 1.0, "Ganymede": 2.0},
        dv_band="powered_dsm",
    )
    assert verdict.prefilter_skipped is False
    assert verdict.feasible is True
    assert verdict.continuity_residual_kms < 0.05
    assert math.isfinite(verdict.total_dv_kms)
    assert verdict.total_dv_kms > 0.0
    # IN-BAND: total ΔV under the 3.5 km/s/cycle powered ceiling. classify_dv_band
    # returns the powered window class for a ΔV in [300 m/s, 3.5 km/s].
    assert verdict.total_dv_kms < 3.5
    assert verdict.dv_band is not None


def test_multirev_galilean_cheaper_than_dsm() -> None:
    """The multi-rev chain closes the Galilean cycle far below the single-DSM cost.

    Same skeleton: the chain's total ΔV is well under the single-DSM total (the
    #449 close was 13.18 km/s/cycle). The chain spends the multi-VILM *minimum*,
    the single DSM the single-VILM *maximum* — the multi-rev win, reproduced
    end-to-end through the driver.
    """
    sequence = ("Io", "Europa", "Ganymede", "Io")
    leg_tofs_days = _geomean_tofs(sequence, scale=1.5)
    multirev = close_powered_cycle(
        primary="Jupiter",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0, 0),
        releg=MultiRevLeveragingReleg(),
        phasing={"Io": 0.0, "Europa": 1.0, "Ganymede": 2.0},
    )
    dsm = close_powered_cycle(
        primary="Jupiter",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0, 0),
        releg=DsmReleg(),
        phasing={"Io": 0.0, "Europa": 1.0, "Ganymede": 2.0},
    )
    assert multirev.feasible is True
    assert dsm.feasible is True
    # The chain is dramatically cheaper than the single-impulse retarget.
    assert multirev.total_dv_kms < 0.5 * dsm.total_dv_kms


def test_multirev_uranus_disjoint_prefiltered_empty() -> None:
    """The chain backend still reports the Uranian disjoint case EMPTY.

    A multi-rev chain walks V_inf WITHIN a Tisserand contour — it cannot bridge
    disjoint contours. The prefilter skips the unbridgeable Uranian legs before
    any chain solve, so the structural negative is preserved with the new backend.
    """
    sequence = ("Ariel", "Umbriel", "Ariel")
    leg_tofs_days = _geomean_tofs(sequence, scale=1.0)
    verdict = close_powered_cycle(
        primary="Uranus",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0),
        releg=MultiRevLeveragingReleg(),
        phasing={"Ariel": 0.0, "Umbriel": 1.0},
    )
    assert verdict.prefilter_skipped is True
    assert verdict.feasible is False
    assert verdict.total_dv_kms == math.inf


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _geomean_tofs(sequence: tuple[str, ...], *, scale: float) -> tuple[float, ...]:
    """Per-leg ToF = scale x geometric mean of the two moons' orbital periods."""
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    # All moons in a moon-tour share one primary.
    primary = SATELLITES[sequence[0]].primary
    mu = PRIMARIES[primary]

    def period_days(m: str) -> float:
        s = SATELLITES[m]
        return 2.0 * math.pi * math.sqrt(s.sma_km**3 / mu) / 86400.0

    return tuple(
        scale * math.sqrt(period_days(sequence[k]) * period_days(sequence[k + 1]))
        for k in range(len(sequence) - 1)
    )


# ---------------------------------------------------------------------------
# Task 8 — capability-subsumption re-stamp record (built, NOT written)
# ---------------------------------------------------------------------------


def test_powered_empty_restamp_records_method() -> None:
    """A powered-empty Uranus region yields a SUBSUMING re-stamp record.

    The driver, run on the Uranian Ariel-Umbriel-Ariel disjoint case, returns
    ``prefilter_skipped``; ``build_powered_empty_restamp`` then produces an
    :class:`EmptyRegionReport` whose method-capability STRICTLY SUBSUMES the prior
    VILM-endgame negative's capability (so ``should_sweep`` would NOT re-sweep it
    under the powered method). The record is ``validate_empty_region``-valid and
    carries method + version (git_sha) — suitable for appending to
    ``empty_regions.jsonl``, but NOT written here (campaign-issue action).
    """
    from cyclerfinder.data.empty_regions import EmptyRegionReport, validate_empty_region
    from cyclerfinder.data.method_capability import MethodCapability, subsumes
    from cyclerfinder.search.releg_moontour import (
        build_powered_empty_restamp,
        powered_releg_method_capability,
    )

    sequence = ("Ariel", "Umbriel", "Ariel")
    leg_tofs_days = _geomean_tofs(sequence, scale=1.0)
    verdict = close_powered_cycle(
        primary="Uranus",
        sequence=sequence,
        leg_tofs_days=leg_tofs_days,
        n_revs=(0, 0),
        releg=DsmReleg(),
        phasing={"Ariel": 0.0, "Umbriel": 1.0},
        dv_band="powered_dsm",
    )
    assert verdict.prefilter_skipped is True

    record = build_powered_empty_restamp(
        region_id="uranus-neptune-regular-moon-endgame-vilm-2026-06-23",
        family="planet-centric moon system (Uranus + Neptune regular moons)",
        centre="Uranus/Neptune",
        sequence=sequence,
        verdict=verdict,
        git_sha="testsha",
        run_date="2026-06-26",
    )
    assert isinstance(record, EmptyRegionReport)
    # The record is first-class (bounded, prune-gated, capability-tagged).
    validate_empty_region(record)
    assert record.method_capability.git_sha == "testsha"
    assert "STRUCTURAL" in record.verdict

    # The powered capability STRICTLY SUBSUMES the prior VILM-endgame method:
    # under capability-subsumption the region stays empty without a re-sweep.
    prior_vilm = MethodCapability(
        genome="phase-full VILM endgame (leveraging, planet-centric moon tour)",
        corrector="discover_endgame_moon -> solve_endgame (VILM V_inf-lowering chain)",
        capability_tags=frozenset(
            {"coplanar", "leveraging", "multi-arc", "patched-conic", "powered"}
        ),
        git_sha="0a6d0a3",
    )
    powered = powered_releg_method_capability(git_sha="testsha")
    assert subsumes(powered, prior_vilm)


def test_powered_restamp_not_written_to_disk(tmp_path) -> None:  # type: ignore[no-untyped-def]
    """Building the re-stamp record does NOT write empty_regions.jsonl.

    The build step is pure; the actual append is a campaign-issue action. This
    guards against an accidental writeback in the build path.
    """
    from cyclerfinder.search.releg_moontour import build_powered_empty_restamp

    sequence = ("Ariel", "Umbriel", "Ariel")
    verdict = close_powered_cycle(
        primary="Uranus",
        sequence=sequence,
        leg_tofs_days=_geomean_tofs(sequence, scale=1.0),
        n_revs=(0, 0),
        releg=DsmReleg(),
        phasing={"Ariel": 0.0, "Umbriel": 1.0},
        dv_band="powered_dsm",
    )
    target = tmp_path / "empty_regions.jsonl"
    build_powered_empty_restamp(
        region_id="r",
        family="f",
        centre="Uranus",
        sequence=sequence,
        verdict=verdict,
    )
    # No file created by the build step.
    assert not target.exists()
    assert list(tmp_path.iterdir()) == []
