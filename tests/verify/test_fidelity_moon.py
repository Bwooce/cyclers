"""Tier-1 Phase 7: Axis-B fidelity persistence tracks Jovicentric V_inf across
rungs (plan Phase 7 Task 7.0). Same _moves_toward_band logic, about the primary."""

from __future__ import annotations

from cyclerfinder.verify.fidelity import persistence_for_primary


def test_jovicentric_vinf_persistence_runs_about_jupiter() -> None:
    report = persistence_for_primary(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        primary="Jupiter",
        tof_seed_days=(4.0, 3.4),
        period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        t0_seed_sec=0.6 * 86400.0,
        slack_leg=2,
    )
    assert report is not None
    # The tracked scalar is V_inf about JUPITER, not the Sun.
    assert report.quantity in ("vinf", "vinf_primary")
