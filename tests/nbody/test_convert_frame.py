"""N-body Phase A: ICRS-equatorial <-> J2000-ecliptic reuses ephemeris.py's obliquity.

NON-GOLDEN cross-implementation check: both sides are OUR rotation; this asserts
the harness uses the SAME obliquity constant as core/ephemeris (design §0: a wrong
obliquity reads as a fake out-of-plane V_inf component).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.nbody.convert import ecliptic_to_icrs_eq, icrs_eq_to_ecliptic


def test_roundtrip_identity() -> None:
    v = np.array([1.0, 2.0, 3.0])
    back = ecliptic_to_icrs_eq(icrs_eq_to_ecliptic(v))
    assert np.allclose(back, v, atol=1e-12)


def test_uses_ephemeris_obliquity_constant() -> None:
    from cyclerfinder.core.ephemeris import _J2000_OBLIQUITY_RAD

    # +z equatorial maps to (0, +sin eps, +cos eps) in ecliptic about +x.
    z_eq = np.array([0.0, 0.0, 1.0])
    got = icrs_eq_to_ecliptic(z_eq)
    eps = _J2000_OBLIQUITY_RAD
    assert np.allclose(got, [0.0, np.sin(eps), np.cos(eps)], atol=1e-12)


def test_matches_ephemeris_backend_rotation() -> None:
    """The conversion is the same R_x(-eps) the astropy backend applies (design §0)."""
    from math import cos, sin

    from cyclerfinder.core.ephemeris import _J2000_OBLIQUITY_RAD

    eps = _J2000_OBLIQUITY_RAD
    expected = np.array([[1.0, 0.0, 0.0], [0.0, cos(eps), sin(eps)], [0.0, -sin(eps), cos(eps)]])
    v = np.array([4.0, -5.0, 6.0])
    assert np.allclose(icrs_eq_to_ecliptic(v), expected @ v, atol=1e-12)
