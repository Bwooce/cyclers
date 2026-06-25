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
from cyclerfinder.search.releg_solver import BallisticReleg, DsmReleg


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
