"""#610 -- tests for the Proteus bend interval-arithmetic certification POC.

Covers three things per this project's standing discipline (see
[[feedback_verify_gauntlet_with_positive_control]]): (1) the hand-built
interval arcsin is actually correct (cross-checked against ``math.asin``),
(2) the interval enclosure of ``max_bend`` matches the non-interval
:func:`cyclerfinder.core.flyby.max_bend` at point inputs, and (3) a POSITIVE
CONTROL -- a box known to exceed the 5 deg gate (Umbriel's real #312 SILVER
V_inf regime, which the module docstring of
``cyclerfinder.search.physical_sanity`` cites as ~14-16 deg) -- to prove the
certifier can also say "NOT certified" and isn't a detector that always fires
"safe". All tests skip cleanly if ``mpmath`` is not installed (same pattern
as the ``validation`` extras group; see pyproject.toml's ``interval`` extra).
"""

from __future__ import annotations

import math

import pytest

mp = pytest.importorskip("mpmath", reason="mpmath is an optional 'interval' extra (task #610 POC)")

import scripts.certify_610_proteus_bend_interval as certify610  # noqa: E402
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
    enclosure = certify610.rigorous_arcsin(iv, iv.mpf(t))
    ref = math.asin(t)
    # NOTE: ``math.asin`` (libm, float64) is itself only correctly-ROUNDED to
    # ~1 double ULP, not exact -- it is not a rigorous enclosure. The
    # mpmath.iv enclosure at 50 dps is far tighter than that, so the two can
    # legitimately differ by up to ~1 ULP of float64 (~1e-16 relative); a
    # strict containment assertion against the double reference would be
    # checking the wrong thing. Use a tolerance a couple of orders above
    # float64 ULP instead -- still >1e10 tighter than the 5 deg gate this
    # module cares about.
    assert abs(float(enclosure.a) - ref) < 1e-12
    assert abs(float(enclosure.b) - ref) < 1e-12
    assert float(enclosure.b) - float(enclosure.a) < 1e-12


@pytest.mark.parametrize(
    ("mu", "rp", "vinf"),
    [
        (2.58342, 308.0, 1.823942965718726),
        (2.58342, 308.0, 13.117027313506517),
        (85.1, 634.7, 0.92),
        (85.1, 634.7, 2.27),
    ],
)
def test_bend_interval_encloses_point_reference(mu: float, rp: float, vinf: float) -> None:
    """The interval enclosure at a degenerate (point) box must contain the
    non-interval reference computed by the project's own trusted
    :func:`cyclerfinder.core.flyby.max_bend` -- cross-checking the new
    interval code against the already-tested point implementation, not just
    itself."""
    iv = mp.iv
    ref_deg = math.degrees(max_bend(mu, rp, vinf))
    enclosure = certify610.bend_deg_interval(iv, mu, rp, rp, vinf, vinf)
    # ``max_bend`` (float64) is only correctly-rounded to ~1 double ULP, not
    # an exact/rigorous value -- see the tolerance note in
    # test_rigorous_arcsin_matches_math_asin. Use an absolute tolerance a
    # couple of orders above float64 ULP at this magnitude, not strict
    # containment against the double reference.
    assert abs(float(enclosure.a) - ref_deg) < 1e-9
    assert abs(float(enclosure.b) - ref_deg) < 1e-9
    assert float(enclosure.b) - float(enclosure.a) < 1e-9, (
        "point-box enclosure should be microscopically tight"
    )


def test_proteus_box_certifies_below_gate() -> None:
    """The exact data-grounded box (#599's real Proteus V_inf survivor range)
    must certify sup(bend) < the #324 gate -- the core positive result of
    this POC."""
    iv = mp.iv
    proteus = SATELLITES["Proteus"]
    rp_safe = proteus.radius_eq_km + proteus.safe_alt_km
    vinf_min, vinf_max, n_total, n_subgate = certify610._proteus_subgate_vinf_range()
    assert n_total == 1024
    assert n_subgate == 104
    enclosure = certify610.bend_deg_interval(
        iv, proteus.mu_km3_s2, rp_safe, rp_safe, vinf_min, vinf_max
    )
    assert float(enclosure.b) < DEFAULT_MIN_USEFUL_BEND_DEG


def test_widened_box_still_certifies_below_gate() -> None:
    """The widened, conservative box (10x periapsis range, V_inf down to a
    rigorously-derived margin below the true survivor minimum) must ALSO
    certify -- demonstrating the negative isn't a lucky escape confined to
    the exact grid #599 happened to evaluate."""
    iv = mp.iv
    proteus = SATELLITES["Proteus"]
    rp_safe = proteus.radius_eq_km + proteus.safe_alt_km
    enclosure = certify610.bend_deg_interval(
        iv, proteus.mu_km3_s2, rp_safe, 10.0 * rp_safe, 0.45, 20.0
    )
    assert float(enclosure.b) < DEFAULT_MIN_USEFUL_BEND_DEG


def test_positive_control_umbriel_is_not_certified_empty() -> None:
    """POSITIVE CONTROL (mandatory per project discipline): a box known to
    exceed the gate must NOT be reported as certified-empty, or the
    certifier would be vacuous (always answering "safe"). Umbriel's real
    #312 SILVER V_inf regime (~0.8-1.0 km/s) delivers ~14-16 deg of bend --
    comfortably over the 5 deg floor -- per
    ``cyclerfinder.search.physical_sanity``'s own module docstring."""
    iv = mp.iv
    umbriel = SATELLITES["Umbriel"]
    rp = umbriel.radius_eq_km + umbriel.safe_alt_km
    enclosure = certify610.bend_deg_interval(iv, umbriel.mu_km3_s2, rp, rp, 0.8, 1.0)
    sup = float(enclosure.b)
    assert sup > DEFAULT_MIN_USEFUL_BEND_DEG, (
        f"positive control failed: expected Umbriel's V_inf~0.8-1.0 km/s box to exceed "
        f"the {DEFAULT_MIN_USEFUL_BEND_DEG} deg gate, got sup={sup:.4f} -- the certifier "
        "may be vacuously always reporting 'safe'."
    )


def test_main_runs_and_certifies() -> None:
    assert certify610.main() == 0
