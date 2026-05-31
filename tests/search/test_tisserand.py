"""Tests for :mod:`cyclerfinder.search.tisserand`.

The Tisserand module supplies the M4 pruning gate (``linkable``); these
tests cover the round-trip identity, contour-shape sanity, the Aldrin-
neighbourhood linkable physics, and robustness under degenerate inputs.

Coplanar (``i = 0``) throughout — see module docstring.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.search.tisserand import (
    linkable,
    linkable_region,
    plot_tisserand,
    tisserand_to_vinf,
    vinf_contour,
    vinf_to_tisserand,
)

# ---------------------------------------------------------------------------
# Conversion round-trip
# ---------------------------------------------------------------------------


def test_vinf_to_tisserand_inverse() -> None:
    """``tisserand_to_vinf o vinf_to_tisserand`` is identity on (0, ~25] km/s."""
    for vinf in (0.5, 3.0, 5.0, 7.0, 15.0):
        t_p = vinf_to_tisserand("E", vinf)
        back = tisserand_to_vinf("E", t_p)
        assert back == pytest.approx(vinf, abs=1.0e-9), (
            f"round-trip failed at vinf={vinf}: -> T={t_p:.6f} -> {back}"
        )


def test_tisserand_to_vinf_above_3_returns_zero() -> None:
    """``T_p >= 3`` -> no real ``V_inf``; convention returns 0."""
    assert tisserand_to_vinf("E", 3.0) == 0.0
    assert tisserand_to_vinf("E", 4.5) == 0.0


def test_vinf_to_tisserand_at_zero_is_three() -> None:
    """At ``V_inf = 0`` the spacecraft co-orbits at the planet, so ``T = 3``."""
    assert vinf_to_tisserand("E", 0.0) == pytest.approx(3.0, abs=1.0e-12)


# ---------------------------------------------------------------------------
# Contour shape sanity
# ---------------------------------------------------------------------------


def test_contour_returns_empty_below_threshold() -> None:
    """A vanishingly small ``V_inf`` gives a contour that collapses to a point at
    ``(a=a_p, e=0)`` -- outside the default ``a_range`` only if filtered, but
    typically inside. We instead test the strictly impossible: at ``V_inf = 0``
    the cubic has the double-root ``u=1`` (not a sign change), so the brentq
    bracketing finds no root and the contour comes back empty. This is the
    documented "no real contour" regime.
    """
    a_au, e = vinf_contour("E", 0.0)
    # At V_inf = 0, T_p = 3, the cubic has a double root at u=1 (no sign
    # change), so the bracketing solver yields no roots -- empty contour.
    assert a_au.size == 0
    assert e.size == 0


def test_contour_non_empty_at_5kms() -> None:
    """A moderate ``V_inf`` produces a populated contour at Earth."""
    a_au, e = vinf_contour("E", 5.0, n_points=200)
    assert a_au.size >= 50, f"contour returned only {a_au.size} points"
    assert e.size == a_au.size


def test_contour_arrays_equal_length() -> None:
    """The two returned arrays always agree in length."""
    for body in ("V", "E", "M"):
        for vinf in (2.0, 5.0, 8.0):
            a_au, e = vinf_contour(body, vinf)
            assert a_au.shape == e.shape


def test_contour_e_within_unit_interval() -> None:
    """Every returned ``e`` is in ``[0, 1)``."""
    _, e = vinf_contour("E", 5.0)
    assert (e >= 0.0).all()
    assert (e < 1.0).all()


def test_contour_passes_near_planet_a_at_low_vinf() -> None:
    """At small V_inf, the contour's a-values cluster around ``a_p`` (Earth = 1 AU)."""
    a_au, _ = vinf_contour("E", 2.0)
    assert a_au.size > 0
    assert np.min(a_au) < 1.0
    assert np.max(a_au) > 1.0


def test_contour_rejects_negative_vinf() -> None:
    """Negative ``V_inf`` is nonsensical."""
    with pytest.raises(ValueError, match="non-negative"):
        vinf_contour("E", -1.0)


def test_contour_rejects_bad_a_range() -> None:
    """``a_min`` must be positive and less than ``a_max``."""
    with pytest.raises(ValueError, match="a_range_au"):
        vinf_contour("E", 5.0, a_range_au=(2.0, 1.0))


# ---------------------------------------------------------------------------
# linkable predicate -- spec §13.3 physics
# ---------------------------------------------------------------------------


def test_linkable_em_at_5_5_kms_true() -> None:
    """Aldrin neighbourhood: E and M ARE linkable at ``V_inf = 5.5 km/s``."""
    assert linkable("E", "M", 5.5) is True


def test_linkable_em_at_0_5_kms_false() -> None:
    """Far below the link threshold: E and M are NOT linkable at 0.5 km/s."""
    assert linkable("E", "M", 0.5) is False


def test_linkable_em_threshold_brackets() -> None:
    """There exists a V_inf below which E-M are not linkable, above which they are.

    Below 1.0 km/s the contours don't span; well above (5.5 km/s) they do.
    This is the bracketing the M4 enumerator will rely on.
    """
    assert linkable("E", "M", 0.5) is False
    assert linkable("E", "M", 5.5) is True


def test_linkable_symmetric() -> None:
    """Linkability is symmetric in the body pair."""
    pairs = [("E", "M"), ("V", "E"), ("V", "M")]
    for a, b in pairs:
        for vinf in (3.0, 5.0, 7.0):
            assert linkable(a, b, vinf) == linkable(b, a, vinf), (
                f"asymmetric at ({a},{b}, V_inf={vinf})"
            )


def test_linkable_region_em_non_empty() -> None:
    """At ``V_inf_cap = 12 km/s`` there are E-M linkable samples."""
    region = linkable_region("E", "M", 12.0, n_vinf=50)
    assert len(region) > 0


def test_linkable_region_empty_below_threshold() -> None:
    """At a tiny cap, the linkable region for E-M is empty."""
    region = linkable_region("E", "M", 0.1, n_vinf=10)
    assert region == []


@pytest.mark.parametrize("vinf", [0.0, 1.0e-6, 1.0e6])
def test_linkable_never_raises(vinf: float) -> None:
    """Pathological V_inf inputs return False without raising."""
    assert linkable("E", "M", vinf) is False


def test_linkable_region_negative_cap() -> None:
    """A non-positive cap returns the empty list."""
    assert linkable_region("E", "M", 0.0) == []
    assert linkable_region("E", "M", -1.0) == []


# ---------------------------------------------------------------------------
# Optional matplotlib plotting helper (viz extra)
# ---------------------------------------------------------------------------


def test_plot_runs_when_matplotlib_present() -> None:
    """If matplotlib is installed, ``plot_tisserand`` returns an Axes without raising."""
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")  # headless backend for CI

    ax = plot_tisserand(["E", "M"], [3.0, 5.0])
    from matplotlib.axes import Axes

    assert isinstance(ax, Axes)
