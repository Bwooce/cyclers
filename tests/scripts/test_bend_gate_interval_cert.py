"""#625 -- tests for the generalized bend-gate interval-certification helper.

Covers the shared primitives in ``scripts/_bend_gate_interval_cert.py``
(factored out of ``scripts/certify_610_proteus_bend_interval.py``, which has
its own pre-existing test file covering the Proteus-specific
reproduction/positive-control checks -- not duplicated here). This file
focuses on the NEW pieces: (1) the primitives still behave correctly when
imported from their new shared location, (2) the structured
:func:`certify_bend_gate_over_box` wrapper is correct and its documented
"sup is realized at the box's lower corner" monotonicity property actually
holds, (3) a POSITIVE CONTROL proving the certifier can say "NOT certified"
(not a detector that always fires "safe" -- see
[[feedback_verify_gauntlet_with_positive_control]]), and (4) a genuine,
hand-verified real-world case where certification correctly FAILS (task
#625's own Sylvia finding -- Romulus's bend exceeds the gate at its real
survivor V_inf minimum), proving the certifier reports negative findings
honestly rather than only ever confirming what's expected.
"""

from __future__ import annotations

import math

import pytest

mp = pytest.importorskip("mpmath", reason="mpmath is an optional 'interval' extra (task #610/#625)")

import scripts._bend_gate_interval_cert as cert  # noqa: E402
from cyclerfinder.core.flyby import max_bend  # noqa: E402
from cyclerfinder.core.satellites import SATELLITES  # noqa: E402
from cyclerfinder.search.physical_sanity import DEFAULT_MIN_USEFUL_BEND_DEG  # noqa: E402


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 50
    mp.iv.dps = 50


@pytest.mark.parametrize("t", [0.001, 0.1, 0.5, 0.9, 0.999])
def test_rigorous_arcsin_matches_math_asin(t: float) -> None:
    iv = mp.iv
    enclosure = cert.rigorous_arcsin(iv, iv.mpf(t))
    ref = math.asin(t)
    assert abs(float(enclosure.a) - ref) < 1e-12
    assert abs(float(enclosure.b) - ref) < 1e-12


def test_bend_interval_encloses_point_reference() -> None:
    """Degenerate (point) box must contain the trusted non-interval reference."""
    iv = mp.iv
    mu, rp, vinf = 2.58342, 308.0, 1.823942965718726
    ref_deg = math.degrees(max_bend(mu, rp, vinf))
    enclosure = cert.bend_deg_interval(iv, mu, rp, rp, vinf, vinf)
    assert abs(float(enclosure.a) - ref_deg) < 1e-9
    assert abs(float(enclosure.b) - ref_deg) < 1e-9


def test_certify_bend_gate_over_box_hand_verified_corner() -> None:
    """Hand-verifiable case: since max_bend is monotonically decreasing in
    both r_p and V_inf (module docstring's own claim), sup(bend) over a box
    must equal the plain-math (non-interval) bend at the box's LOWER corner
    exactly -- computed independently here via math.asin, not copied from
    the interval code under test."""
    iv = mp.iv
    proteus = SATELLITES["Proteus"]
    rp_safe = proteus.radius_eq_km + proteus.safe_alt_km
    vinf_lo, vinf_hi = 1.823943, 13.117027

    result = cert.certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=proteus.mu_km3_s2,
        rp_lo_km=rp_safe,
        rp_hi_km=10.0 * rp_safe,
        vinf_lo_kms=vinf_lo,
        vinf_hi_kms=vinf_hi,
        gate_deg=DEFAULT_MIN_USEFUL_BEND_DEG,
        label="hand-verify",
    )

    # Independent hand computation at the lower corner (rp_lo, vinf_lo).
    hand_bend_rad = 2.0 * math.asin(proteus.mu_km3_s2 / (proteus.mu_km3_s2 + rp_safe * vinf_lo**2))
    hand_bend_deg = math.degrees(hand_bend_rad)

    assert abs(result["sup_bend_deg"] - hand_bend_deg) < 1e-6
    assert result["certified"] is True
    assert result["box"]["vinf_kms"] == [vinf_lo, vinf_hi]
    assert result["gate_deg"] == DEFAULT_MIN_USEFUL_BEND_DEG
    assert result["label"] == "hand-verify"
    assert "Does NOT certify" in result["scope_note"]


def test_sup_independent_of_upper_bounds() -> None:
    """Per the module docstring's monotonicity claim: widening the box's
    UPPER bounds must not change sup(bend) at all -- sup is realized at the
    lower corner regardless."""
    iv = mp.iv
    proteus = SATELLITES["Proteus"]
    rp_safe = proteus.radius_eq_km + proteus.safe_alt_km
    tight = cert.certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=proteus.mu_km3_s2,
        rp_lo_km=rp_safe,
        rp_hi_km=rp_safe,
        vinf_lo_kms=2.0,
        vinf_hi_kms=2.0,
        gate_deg=5.0,
    )
    wide = cert.certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=proteus.mu_km3_s2,
        rp_lo_km=rp_safe,
        rp_hi_km=1000.0 * rp_safe,
        vinf_lo_kms=2.0,
        vinf_hi_kms=500.0,
        gate_deg=5.0,
    )
    assert abs(tight["sup_bend_deg"] - wide["sup_bend_deg"]) < 1e-9


def test_positive_control_umbriel_is_not_certified() -> None:
    """POSITIVE CONTROL (mandatory per project discipline): a box known to
    exceed the gate must report certified=False, or the certifier would be
    vacuous."""
    iv = mp.iv
    umbriel = SATELLITES["Umbriel"]
    rp = umbriel.radius_eq_km + umbriel.safe_alt_km
    result = cert.certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=umbriel.mu_km3_s2,
        rp_lo_km=rp,
        rp_hi_km=rp,
        vinf_lo_kms=0.8,
        vinf_hi_kms=1.0,
        gate_deg=DEFAULT_MIN_USEFUL_BEND_DEG,
    )
    assert result["certified"] is False
    assert result["sup_bend_deg"] > DEFAULT_MIN_USEFUL_BEND_DEG


def test_genuine_negative_finding_sylvia_romulus_not_certified() -> None:
    """Task #625's own genuine (not forced-through) finding: Sylvia's Romulus
    moon, at its REAL residual-sub-gate-survivor minimum V_inf (reproduced
    in scripts/certify_625_bend_gate_registry.py), bends slightly ABOVE the
    5 deg gate -- so this specific body/box does NOT certify. Hand-verified
    against a plain-math computation, confirming the interval code reports
    this honestly rather than only ever confirming certification."""
    iv = mp.iv
    romulus = SATELLITES["Romulus"]
    rp_safe = romulus.radius_eq_km + romulus.safe_alt_km
    vinf_min = 0.007831273792325793  # reproduced real Sylvia Romulus->Remus survivor minimum

    hand_bend_deg = math.degrees(
        2.0 * math.asin(romulus.mu_km3_s2 / (romulus.mu_km3_s2 + rp_safe * vinf_min**2))
    )
    result = cert.certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=romulus.mu_km3_s2,
        rp_lo_km=rp_safe,
        rp_hi_km=rp_safe,
        vinf_lo_kms=vinf_min,
        vinf_hi_kms=1.0,
        gate_deg=DEFAULT_MIN_USEFUL_BEND_DEG,
    )
    assert abs(result["sup_bend_deg"] - hand_bend_deg) < 1e-6
    assert hand_bend_deg > DEFAULT_MIN_USEFUL_BEND_DEG
    assert result["certified"] is False


def test_have_mpmath_flag_true_when_mpmath_importable() -> None:
    assert cert.HAVE_MPMATH is True
