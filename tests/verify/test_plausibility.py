"""Publication-layer plausibility predicates (task #127, deliverable 2).

The 55.32 km/s near-miss is the regression fixture: feeding the off-family
maintenance-ΔV value MUST be refused; the in-family 2.9138 MUST pass. V_inf
predicates layer the physics ceiling, the data invariant, and the sourced
distribution.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.verify.plausibility import (
    MAINTENANCE_DV_CONVENTION_KMS,
    VINF_DATA_INVARIANT_KMS,
    PlausibilityVerdict,
    QuantityKind,
    check_publishable,
    sourced_vinf_max_kms,
)

# --- Maintenance-ΔV: the 55.32 regression fixture --------------------------


def test_off_family_5532_maintenance_dv_is_refused() -> None:
    """REGRESSION (#127): the 55.32 km/s off-family solve must be refused."""
    v = check_publishable(QuantityKind.MAINTENANCE_DV_KMS, 55.32)
    assert v.ok is False
    assert "55.32" in v.reason or "55.3" in v.reason
    assert "engineering bar" in v.reason


def test_in_family_29138_maintenance_dv_passes() -> None:
    """The published in-family value MUST stay publishable (nothing regresses)."""
    v = check_publishable(QuantityKind.MAINTENANCE_DV_KMS, 2.9138)
    assert v.ok is True


def test_maintenance_dv_at_bar_passes_just_above_refused() -> None:
    assert check_publishable(QuantityKind.MAINTENANCE_DV_KMS, MAINTENANCE_DV_CONVENTION_KMS).ok
    assert not check_publishable(
        QuantityKind.MAINTENANCE_DV_KMS, MAINTENANCE_DV_CONVENTION_KMS + 1e-6
    ).ok


def test_maintenance_dv_negative_and_nonfinite_refused() -> None:
    assert not check_publishable(QuantityKind.MAINTENANCE_DV_KMS, -0.5).ok
    assert not check_publishable(QuantityKind.MAINTENANCE_DV_KMS, math.inf).ok
    assert not check_publishable(QuantityKind.MAINTENANCE_DV_KMS, math.nan).ok


# --- V_inf: physics ceiling / data invariant / sourced distribution --------


def test_sourced_max_is_the_live_catalogue_value() -> None:
    # The live catalogue currently tops out at 20.3 km/s (Russell-Ocampo).
    assert sourced_vinf_max_kms() == pytest.approx(20.3, abs=0.01)


def test_russell_max_vinf_passes() -> None:
    v = check_publishable(QuantityKind.VINF_KMS, 20.3, {"body": "E"})
    assert v.ok is True


def test_over_ceiling_vinf_at_earth_refused_as_physics() -> None:
    v = check_publishable(QuantityKind.VINF_KMS, 80.0, {"body": "E"})
    assert v.ok is False
    assert "PHYSICALLY IMPOSSIBLE" in v.reason


def test_per_body_physics_ceiling_distinguishes_bodies() -> None:
    # Below the data invariant (50), the per-body physics ceiling is the binding
    # bar: 20 km/s is below the Uranus ceiling (16.4? no -> above). Use Neptune
    # (ceiling ~13.1) vs Earth (71.9): 25 km/s is physically impossible at
    # Neptune but fine (physics) at Earth.
    assert check_publishable(QuantityKind.VINF_KMS, 25.0, {"body": "E"}).ok is True
    v_n = check_publishable(QuantityKind.VINF_KMS, 25.0, {"body": "N"})
    assert v_n.ok is False
    assert "PHYSICALLY IMPOSSIBLE" in v_n.reason


def test_unit_error_vinf_refused_by_data_invariant() -> None:
    # 49 km/s is sub-Earth-ceiling (71.9) but >= the 50 data invariant after
    # rounding? use a value between the distribution bar and the data invariant.
    v = check_publishable(QuantityKind.VINF_KMS, VINF_DATA_INVARIANT_KMS, {"body": "V"})
    assert v.ok is False
    assert "data-layer invariant" in v.reason


def test_far_outlier_vinf_refused_by_distribution() -> None:
    # 40 km/s: sub-Earth-ceiling, sub-data-invariant, but way above sourced+headroom.
    v = check_publishable(QuantityKind.VINF_KMS, 40.0, {"body": "V"})
    assert v.ok is False
    assert "sourced-distribution bar" in v.reason


def test_within_headroom_vinf_passes_with_note() -> None:
    # 22 km/s: above sourced max 20.3 but within +5 headroom.
    v = check_publishable(QuantityKind.VINF_KMS, 22.0, {"body": "E"})
    assert v.ok is True
    assert "headroom" in v.reason


def test_no_body_uses_most_permissive_ceiling() -> None:
    # Without a body, the physics ceiling is Mercury's (~115.6); a value below it
    # but above distribution is refused by the distribution bar, not physics.
    v = check_publishable(QuantityKind.VINF_KMS, 100.0, {})
    assert v.ok is False
    # 100 < 115.6 Mercury ceiling -> not the physics branch; data invariant trips.
    assert "data-layer invariant" in v.reason


def test_unknown_quantity_kind_raises() -> None:
    with pytest.raises(ValueError, match="unknown quantity_kind"):
        check_publishable("frobnicate", 1.0)


def test_verdict_is_frozen_dataclass() -> None:
    v = PlausibilityVerdict(ok=True, reason="x")
    with pytest.raises((AttributeError, TypeError)):
        v.ok = False  # type: ignore[misc]
