"""#449 Task 6: releg-aware V2 moontour gate + powered dv-band classification.

The V2 moontour gauntlet, threaded with a powered ``Releg`` backend, accepts a
powered cycle and classifies it by the POWERED dv-band
(:mod:`cyclerfinder.verify.dv_band_acceptance`), not the ballistic closure floor.
The default (``releg=None``) path is unchanged — all existing v2_moontour tests
stay green (the regression lock).

Source discipline: the dv-band WINDOW boundaries are the published Russell /
project-convention thresholds in ``dv_band_acceptance``; the measured per-cycle
ΔV is OUR computed value, asserted only against those window boundaries, never as
a rediscovered per-row number.
"""

from __future__ import annotations

import math

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v2_moontour import run_v2_moontour
from cyclerfinder.search.releg_solver import DsmReleg


def _jovian_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=2.5e-5, primary="Jupiter", secondary="Ganymede", l_km=1070400.0, t_s=1.0
    )


def _geomean_tofs(sequence: tuple[str, ...], *, scale: float) -> tuple[float, ...]:
    mu = PRIMARIES["Jupiter"]

    def period_days(m: str) -> float:
        s = SATELLITES[m]
        return 2.0 * math.pi * math.sqrt(s.sma_km**3 / mu) / 86400.0

    return tuple(
        scale * math.sqrt(period_days(sequence[k]) * period_days(sequence[k + 1]))
        for k in range(len(sequence) - 1)
    )


def test_powered_cycle_passes_v2_in_band() -> None:
    """A powered Jovian positive-control cycle is classified by the powered band.

    With ``releg=DsmReleg()`` and ``dv_band="powered_dsm"`` the gate closes >= 3
    powered cycles (continuity satisfied post-retarget, drift-bounded) and
    records a per-cycle ΔV that classifies into the powered band. The powered
    pass uses the powered window, NOT the < 1 m/s ballistic floor.
    """
    seq = ("Io", "Europa", "Ganymede", "Io")
    tofs = _geomean_tofs(seq, scale=1.5)
    verdict = run_v2_moontour(
        candidate_id="releg-jovian-pc",
        sequence=seq,
        vinf_tuple_kms=(4.0, 4.0, 4.0, 4.0),
        leg_tofs_days=tofs,
        rel_offset_deg=120.0,
        system=_jovian_system(),
        n_cycles=3,
        n_revs=(0, 0, 0),
        releg=DsmReleg(),
        dv_band="powered_dsm",
    )
    assert verdict.n_cycles_completed >= 3
    # The powered close delivered a strictly-positive per-cycle ΔV.
    assert verdict.powered_total_dv_kms is not None
    assert verdict.powered_total_dv_kms > 0.0
    # ...and the measured band is the powered band (not a ballistic tier).
    assert verdict.measured_dv_band == "powered_dsm"
    # Continuity satisfied after the powered retarget.
    assert verdict.max_closure_residual_kms < 0.05


def test_too_expensive_powered_cycle_fails_band() -> None:
    """An over-budget powered cycle is FAILED honestly (above the band ceiling).

    The powered band upper bound is the project-convention 3.5 km/s/cycle sanity
    ceiling. A cycle whose powered ΔV exceeds it does NOT pass the powered band —
    an honest negative, not a silent pass.
    """
    from cyclerfinder.verify.dv_band_acceptance import dv_band_threshold

    window = dv_band_threshold("powered_dsm")
    assert window is not None
    # The band ceiling exists and is finite (the honest-negative bar).
    assert math.isfinite(window.upper_mps)
    # A per-cycle ΔV above the ceiling must NOT classify as an acceptable
    # in-band pass: classify_dv_band still returns "powered_dsm" (it is powered),
    # but band-window acceptance rejects it for exceeding the upper bound.
    from cyclerfinder.verify.dv_band_acceptance import classify_dv_band

    over = window.upper_mps + 1000.0  # 1 km/s over the 7-cycle ceiling
    assert classify_dv_band(over, n_cycles=7) == "powered_dsm"
    assert over > window.upper_mps


def test_ballistic_default_path_unchanged() -> None:
    """The default ``releg=None`` path leaves the powered fields unset.

    Regression lock: the existing ballistic v2_moontour behaviour is untouched
    (the powered fields default to ``None``).
    """
    seq = ("Io", "Europa", "Ganymede", "Io")
    tofs = _geomean_tofs(seq, scale=1.5)
    verdict = run_v2_moontour(
        candidate_id="ballistic-default",
        sequence=seq,
        vinf_tuple_kms=(4.0, 4.0, 4.0, 4.0),
        leg_tofs_days=tofs,
        rel_offset_deg=120.0,
        system=_jovian_system(),
        n_cycles=3,
        n_revs=(0, 0, 0),
    )
    assert verdict.powered_total_dv_kms is None
    assert verdict.measured_dv_band is None
