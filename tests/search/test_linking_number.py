"""Tests for the linking-number primitive (#522, Owen & Baresi 2024).

POSITIVE CONTROLS: the standard closed-form Hopf-link pair of unit circles
(textbook two-ring chain-link construction -- EXPECTED linking number +-1 is
a topological identity, not a value this module computed) and closed-form
unlinked/nested-coplanar pairs (EXPECTED linking number 0).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.search.linking_number import linking_number


def _circle(
    n: int,
    *,
    center: tuple[float, float, float],
    plane: str,
    radius: float = 1.0,
    reverse: bool = False,
) -> NDArray[np.float64]:
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    if reverse:
        t = -t
    a = radius * np.cos(t)
    b = radius * np.sin(t)
    zero = np.zeros_like(t)
    if plane == "xy":
        pts = np.stack([a, b, zero], axis=1)
    elif plane == "xz":
        pts = np.stack([a, zero, b], axis=1)
    elif plane == "yz":
        pts = np.stack([zero, a, b], axis=1)
    else:
        raise ValueError(plane)
    return pts + np.asarray(center, dtype=np.float64)


def test_hopf_link_has_linking_number_one() -> None:
    """Standard textbook chain-link pair: unit circle in the xy-plane at the
    origin, unit circle in the xz-plane centered at (0.3,0,0) -- threads
    through the first circle's disc at (-0.7,0,0) (deliberately OFF the
    disc's centroid at the origin, to avoid the fan triangulation's shared
    apex vertex -- piercing exactly through that degenerate point registers
    against multiple adjacent triangles and over-counts). Closed-form
    linking number +-1.
    """
    c1 = _circle(400, center=(0.0, 0.0, 0.0), plane="xy")
    c2 = _circle(400, center=(0.3, 0.0, 0.0), plane="xz")
    assert abs(linking_number(c1, c2)) == 1


def test_hopf_link_reversed_orientation_flips_sign() -> None:
    c1 = _circle(400, center=(0.0, 0.0, 0.0), plane="xy")
    c2 = _circle(400, center=(0.3, 0.0, 0.0), plane="xz")
    c2_rev = _circle(400, center=(0.3, 0.0, 0.0), plane="xz", reverse=True)
    l_fwd = linking_number(c1, c2)
    l_rev = linking_number(c1, c2_rev)
    assert l_fwd == -l_rev
    assert abs(l_fwd) == 1


def test_far_apart_circles_are_unlinked() -> None:
    c1 = _circle(200, center=(0.0, 0.0, 0.0), plane="xy")
    c2 = _circle(200, center=(50.0, 0.0, 0.0), plane="xz")
    assert linking_number(c1, c2) == 0


def test_nested_coplanar_circles_are_unlinked() -> None:
    """Two circles in the SAME plane, one nested inside the other -- close
    in appearance but topologically unlinked (no threading through the
    disc's normal direction).
    """
    c1 = _circle(300, center=(0.0, 0.0, 0.0), plane="xy", radius=2.0)
    c2 = _circle(300, center=(0.0, 0.0, 0.0), plane="xy", radius=0.5)
    assert linking_number(c1, c2) == 0


def test_side_by_side_coplanar_circles_are_unlinked() -> None:
    """Two circles in the same plane, non-overlapping -- also unlinked."""
    c1 = _circle(300, center=(0.0, 0.0, 0.0), plane="xy", radius=1.0)
    c2 = _circle(300, center=(3.0, 0.0, 0.0), plane="xy", radius=1.0)
    assert linking_number(c1, c2) == 0
