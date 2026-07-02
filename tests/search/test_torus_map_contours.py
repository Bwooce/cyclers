"""Tests for marching-squares level-curve extraction (#522).

POSITIVE CONTROLS: a closed-form circular level set (f(x,y)=x^2+y^2=1, a
plain non-periodic grid) and a closed-form periodic-grid case (a field
depending on one coordinate only, whose level crossings must form closed
loops by wrapping through the OTHER coordinate's periodicity) -- both
EXPECTED shapes are geometric identities, not values this module computed.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.search.torus_map_contours import marching_squares_contours


def test_circle_level_set_on_plain_grid() -> None:
    """f(x,y) = x^2 + y^2 on [-2,2]^2, level=1 -- the unit circle."""
    n = 121
    xs = np.linspace(-2.0, 2.0, n)
    ys = np.linspace(-2.0, 2.0, n)
    field = np.array([[x * x + y * y for y in ys] for x in xs])

    contours = marching_squares_contours(field, 1.0, periodic=(False, False))
    assert len(contours.polylines) == 1
    assert contours.closed[0]

    poly = contours.polylines[0]
    # Map grid-index coordinates back to physical (x, y).
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    x_phys = xs[0] + poly[:, 0] * dx
    y_phys = ys[0] + poly[:, 1] * dy
    radii = np.sqrt(x_phys**2 + y_phys**2)
    assert np.max(np.abs(radii - 1.0)) < 2.0 * max(dx, dy)

    # Arc length should be close to 2*pi (the unit circle's circumference).
    diffs = np.diff(np.stack([x_phys, y_phys], axis=1), axis=0)
    arc_length = float(np.sum(np.linalg.norm(diffs, axis=1)))
    assert abs(arc_length - 2.0 * np.pi) < 0.05


def test_periodic_grid_wraps_into_closed_loops() -> None:
    """field(i, j) = cos(2*pi*i/ni), constant in j -- level=0 crossings are
    two vertical-in-i lines, each closing into a loop via j-periodicity.
    """
    ni, nj = 60, 40
    i_idx = np.arange(ni)
    field = np.tile(np.cos(2.0 * np.pi * i_idx / ni)[:, None], (1, nj))

    contours = marching_squares_contours(field, 0.0, periodic=(True, True))
    assert len(contours.polylines) == 2
    assert all(contours.closed)
    for poly in contours.polylines:
        # Each loop spans the full j-range (closes by wrapping through j).
        assert np.max(poly[:, 1]) - np.min(poly[:, 1]) > 0.9 * nj
        # And stays at an i near the two theoretical crossings (ni/4, 3ni/4).
        i_mean = float(np.mean(poly[:, 0]))
        assert min(abs(i_mean - ni / 4.0), abs(i_mean - 3.0 * ni / 4.0)) < 1.0


def test_no_crossing_returns_empty() -> None:
    field = np.ones((20, 20)) * 5.0
    contours = marching_squares_contours(field, 100.0, periodic=(True, True))
    assert contours.polylines == []
