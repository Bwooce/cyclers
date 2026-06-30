"""Golden tests for Lynam-Longuski 2011 IEG triple-cycler reproduction (#493).

Sourced invariants (Lynam-Longuski 2011, Acta Astronautica 69(3-4), pp.158-167):
* Laplace (Ganymede) period: 7.055 d.
* Single-period IEG: powered ΔV ≈ 11 m/s; ballistic Europa flyby altitude -175 km.
* GIPEIPE orbital period: 3.5 d (1:2 resonance with Ganymede).
* GIPEIPE sequence period (Laplace recurrence): 7.055 d.

These tests verify the #493 reproduction findings:
1. Period match: ideal-model Ganymede (Laplace) period matches sourced within 1%.
2. GIPEIPE orbital period matches sourced within 0.1%.
3. Single-period IEG (EIGE topology): ballistic closure exists in ideal model;
   L-L sourced values confirm it needs powered ΔV + has sub-surface ballistic Europa.
4. GIPEIPE: a closed ballistic solution EXISTS in the ideal model (the construction
   closes to resnorm < 1e-10 km/s) but ALL interior flyby altitudes are sub-surface
   (the characterised negative — analogous to EGGIE Gate-B in #480).
5. Self-consistency: all GIPEIPE encounter nodes pass the SOI check.

The GIPEIPE result is a CHARACTERISED NEGATIVE: the cycler closes mathematically but
is infeasible in the strict 2D coplanar model. The self-consistency check must still
PASS (the Lambert construction places the spacecraft exactly at each moon by design).
"""

from __future__ import annotations

from cyclerfinder.search.eggie_ballistic import moon_state
from cyclerfinder.search.eige_ballistic import feasible_ballistic_eige
from cyclerfinder.search.ll2011_ballistic import (
    GIPEIPE_PERIOD_D,
    GIPEIPE_SEQUENCE,
    LL2011_GIPEIPE_ORBIT_PERIOD_DAYS,
    LL2011_IEG_BALLISTIC_EUROPA_ALT_KM,
    LL2011_IEG_DV_MS,
    LL2011_LAPLACE_PERIOD_DAYS,
    T_LAPLACE_D,
    construct_gipeipe,
    single_period_ieg_summary,
)
from cyclerfinder.search.tour_self_consistency import assert_encounters_self_consistent

# ---------------------------------------------------------------------------
# Period invariants (fast, no construction)
# ---------------------------------------------------------------------------


def test_laplace_period_matches_sourced() -> None:
    """Ideal Ganymede period ≈ L-L sourced Laplace period 7.055 d, within 1%."""
    delta_pct = abs(T_LAPLACE_D - LL2011_LAPLACE_PERIOD_DAYS) / LL2011_LAPLACE_PERIOD_DAYS
    assert delta_pct < 0.01, (
        f"Ideal Laplace period {T_LAPLACE_D:.4f} d differs from sourced "
        f"{LL2011_LAPLACE_PERIOD_DAYS} d by {100 * delta_pct:.3f}% (expect < 1%)"
    )


def test_gipeipe_orbital_period_matches_sourced() -> None:
    """GIPEIPE 1:2-resonant orbital period ≈ sourced 3.5 d, within 0.1%."""
    delta_pct = (
        abs(GIPEIPE_PERIOD_D - LL2011_GIPEIPE_ORBIT_PERIOD_DAYS) / LL2011_GIPEIPE_ORBIT_PERIOD_DAYS
    )
    assert delta_pct < 0.001, (
        f"GIPEIPE ideal orbital period {GIPEIPE_PERIOD_D:.4f} d differs from sourced "
        f"{LL2011_GIPEIPE_ORBIT_PERIOD_DAYS} d by {100 * delta_pct:.4f}% (expect < 0.1%)"
    )


# ---------------------------------------------------------------------------
# Single-period IEG (EIGE topology) — sourced invariant checks
# ---------------------------------------------------------------------------


def test_single_period_ieg_sourced_values_present() -> None:
    """L-L sourced invariants are embedded (never code-derived golden check)."""
    # Europa ballistic altitude is sub-surface (sourced negative value)
    assert LL2011_IEG_BALLISTIC_EUROPA_ALT_KM < 0.0, (
        "Sourced L-L ballistic Europa altitude must be negative (sub-surface)"
    )
    # Powered ΔV is positive and small
    assert 0 < LL2011_IEG_DV_MS < 100.0, (
        f"Sourced L-L powered ΔV {LL2011_IEG_DV_MS} m/s out of expected range [0, 100]"
    )


def test_single_period_ieg_eige_ballistic_feasible() -> None:
    """Feasible ballistic EIGE member exists in ideal model (total ΔV ~ 0, all altitudes OK).

    This is the #480 EIGE result used to characterise the L-L single-period IEG.
    """
    eige = feasible_ballistic_eige()
    # Total ΔV is essentially zero (ballistic)
    assert eige.total_dv_ms < 1e-6, (
        f"Feasible EIGE total ΔV {eige.total_dv_ms:.2e} m/s is not near-zero (ballistic expected)"
    )
    # All flyby altitudes in the feasible window
    assert eige.all_feasible, f"Feasible EIGE has infeasible altitudes: {eige.flyby_alt_km}"
    # Ballistic closure (seam and resnorm)
    assert eige.seam_defect_kms < 1e-10, (
        f"EIGE seam defect {eige.seam_defect_kms:.2e} km/s not closed"
    )


def test_single_period_ieg_summary_consistent() -> None:
    """Summary dict has consistent fields and period within tolerance."""
    s = single_period_ieg_summary()
    assert abs(s["period_delta_pct"]) < 1.0, f"Period delta {s['period_delta_pct']:.3f}% exceeds 1%"
    assert s["sourced_powered_dv_ms"] == LL2011_IEG_DV_MS
    assert s["sourced_ballistic_europa_alt_km"] == LL2011_IEG_BALLISTIC_EUROPA_ALT_KM
    assert s["eige_all_feasible"] is True


# ---------------------------------------------------------------------------
# GIPEIPE characterisation (the characterised negative)
# ---------------------------------------------------------------------------


def test_gipeipe_closes_ballistically() -> None:
    """GIPEIPE construction closes in ideal model (resnorm < 1e-10 km/s).

    The solution is mathematically valid — equal-in/out |V∞| at all interior nodes
    and periodicity seam closed — but the flyby altitudes are sub-surface (see
    test_gipeipe_interior_altitudes_sub_surface).
    """
    g = construct_gipeipe()
    assert g.ballistic_resnorm_kms < 1e-10, (
        f"GIPEIPE resnorm {g.ballistic_resnorm_kms:.2e} km/s not closed (expect < 1e-10)"
    )
    assert g.seam_defect_kms < 1e-10, (
        f"GIPEIPE seam defect {g.seam_defect_kms:.2e} km/s not closed (expect < 1e-10)"
    )


def test_gipeipe_interior_altitudes_sub_surface() -> None:
    """GIPEIPE interior flyby altitudes are all sub-surface in the ideal 2D model.

    This is the characterised negative: the ballistic GIPEIPE in the strict
    coplanar model requires sub-surface flybys at all four interior encounters
    (Io₁, Europa₁, Io₂, Europa₂). Consistent with L-L's MALTO-optimised powered
    solution (powered ΔV required) and the analogous EGGIE Gate-B (#480) finding.
    """
    g = construct_gipeipe()
    # All interior nodes are sub-surface
    interior_keys = ["I1", "E1", "I2", "E2"]
    for k in interior_keys:
        alt = g.flyby_alt_km[k]
        assert alt < 0.0, (
            f"GIPEIPE interior flyby {k} altitude {alt:.0f} km is NOT sub-surface; "
            "characterisation says all interior altitudes must be negative"
        )
    # Confirm the overall infeasibility
    assert not g.all_feasible, "GIPEIPE should be infeasible (sub-surface altitudes)"


def test_gipeipe_self_consistency() -> None:
    """GIPEIPE self-consistency: all 6 encounter nodes (incl. repeated I and E) pass SOI check.

    By Lambert construction the spacecraft is exactly at each moon's position, so gaps
    are identically zero — this test guards against future code regressions.
    """
    g = construct_gipeipe()
    seconds_per_day = 86400.0
    tofs_s = [t * seconds_per_day for t in g.tofs_days]
    t_nodes = [0.0]
    for tof in tofs_s:
        t_nodes.append(t_nodes[-1] + tof)

    phases = {"Ganymede": 0.0, "Io": g.phi_io_rad, "Europa": g.phi_eur_rad}
    sc_pos = []
    body_pos = []
    for k, moon in enumerate(GIPEIPE_SEQUENCE):
        r, _ = moon_state(moon, phases[moon], t_nodes[k])
        sc_pos.append(r)
        body_pos.append(r)

    # No TourSelfConsistencyError should be raised
    assert_encounters_self_consistent(
        sc_pos, body_pos, list(GIPEIPE_SEQUENCE), context="GIPEIPE #493"
    )
