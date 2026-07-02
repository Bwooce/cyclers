"""Level-curve (contour) extraction on a regular 2D grid (#522, marching squares).

Owen & Baresi 2024's method (see ``docs/notes/2026-07-03-digest-owen-baresi-
2024-knot-theory-heteroclinic.md``) needs, for a chosen scanning state
component and a value ``D``, the set of torus angles ``(theta_long,
theta_trans)`` at which the torus map equals ``D`` -- a level curve on a
regular grid. This is the standard marching-squares algorithm (Lorensen &
Cline 1987), applied here to a grid that may be PERIODIC in either axis
(torus angles wrap at ``2*pi``), which a torus map genuinely is.

Pure: numpy only, no CR3BP/torus-specific code -- operates on a plain 2D
array. Callers map the returned grid-index-space polylines to actual
``(theta_long, theta_trans)`` angles via their own grid spacing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

# Standard marching-squares case table: for each of the 16 corner-sign
# combinations (bit 0 = bottom-left, bit 1 = bottom-right, bit 2 = top-right,
# bit 3 = top-left; bit set means corner value > level), the list of
# edge-pairs to connect. Edges: 0=bottom, 1=right, 2=top, 3=left.
# Cases 5 and 10 are the ambiguous saddle cases, resolved via the cell
# center value at call time (see ``_cell_segments``).
_EDGE_PAIRS: dict[int, list[tuple[int, int]]] = {
    0: [],
    1: [(3, 0)],
    2: [(0, 1)],
    3: [(3, 1)],
    4: [(1, 2)],
    5: [(3, 0), (1, 2)],  # ambiguous -- resolved by caller
    6: [(0, 2)],
    7: [(3, 2)],
    8: [(2, 3)],
    9: [(0, 2)],  # same pair as case 6, opposite orientation not tracked
    10: [(3, 0), (1, 2)],  # ambiguous -- resolved by caller (mirror of 5)
    11: [(1, 2)],
    12: [(3, 1)],
    13: [(0, 1)],
    14: [(3, 0)],
    15: [],
}


@dataclass(frozen=True)
class ContourSet:
    """Extracted level curves, as a list of open or closed polylines.

    Each polyline is an ``(n, 2)`` array of ``(i, j)`` FRACTIONAL grid-index
    coordinates (not yet scaled to physical angles). ``closed[k]`` is
    ``True`` iff polyline ``k``'s first and last points coincide (within
    floating-point tolerance) -- the case a genuine closed torus-map level
    curve should produce on a periodic grid.
    """

    polylines: list[NDArray[np.float64]]
    closed: list[bool]


def _edge_point(
    field: NDArray[np.float64], level: float, i0: int, j0: int, i1: int, j1: int
) -> tuple[float, float]:
    """Linear-interpolated crossing point on the grid edge ``(i0,j0)-(i1,j1)``."""
    v0 = field[i0, j0]
    v1 = field[i1, j1]
    denom = v1 - v0
    t = 0.5 if abs(denom) < 1e-300 else (level - v0) / denom
    t = min(max(t, 0.0), 1.0)
    return (i0 + t * (i1 - i0), j0 + t * (j1 - j0))


def _corner_indices(
    i: int, j: int, ni: int, nj: int, periodic: tuple[bool, bool]
) -> tuple[int, int, int, int]:
    """Return ``(i, i+1_wrapped, j, j+1_wrapped)`` for cell ``(i, j)``."""
    i1 = (i + 1) % ni if periodic[0] else i + 1
    j1 = (j + 1) % nj if periodic[1] else j + 1
    return i, i1, j, j1


def marching_squares_contours(
    field: NDArray[np.float64],
    level: float,
    *,
    periodic: tuple[bool, bool] = (True, True),
) -> ContourSet:
    """Extract level-``level`` contours of a 2D scalar ``field`` (shape ``(ni, nj)``).

    ``periodic`` controls whether the grid wraps at each axis's far edge
    (the natural case for torus angles). Returns raw line-segment endpoints
    stitched into polylines by shared-point matching; segments sharing a
    grid-edge crossing point are stitched EXACTLY (the crossing point for a
    shared edge is computed identically regardless of which adjacent cell
    computes it, so no distance-tolerance matching is needed).
    """
    ni, nj = field.shape
    ni_cells = ni if periodic[0] else ni - 1
    nj_cells = nj if periodic[1] else nj - 1

    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []

    # Corner-index -> local-edge-index -> point cache, keyed by (edge_axis,
    # i, j) so shared edges between adjacent cells resolve to the SAME
    # floating-point point object.
    edge_cache: dict[tuple[str, int, int], tuple[float, float]] = {}

    def _bottom(i: int, j0: int, j1: int) -> tuple[float, float]:
        key = ("h", i, j0)
        if key not in edge_cache:
            edge_cache[key] = _edge_point(field, level, i, j0, i, j1)
        return edge_cache[key]

    def _top(i1: int, j0: int, j1: int) -> tuple[float, float]:
        key = ("h", i1, j0)
        if key not in edge_cache:
            edge_cache[key] = _edge_point(field, level, i1, j0, i1, j1)
        return edge_cache[key]

    def _left(i0: int, i1: int, j: int) -> tuple[float, float]:
        key = ("v", i0, j)
        if key not in edge_cache:
            edge_cache[key] = _edge_point(field, level, i0, j, i1, j)
        return edge_cache[key]

    def _right(i0: int, i1: int, j1: int) -> tuple[float, float]:
        key = ("v", i0, j1)
        if key not in edge_cache:
            edge_cache[key] = _edge_point(field, level, i0, j1, i1, j1)
        return edge_cache[key]

    for i in range(ni_cells):
        for j in range(nj_cells):
            i0, i1, j0, j1 = _corner_indices(i, j, ni, nj, periodic)
            v_bl = field[i0, j0]
            v_br = field[i0, j1]
            v_tr = field[i1, j1]
            v_tl = field[i1, j0]
            case = (
                (1 if v_bl > level else 0)
                | (2 if v_br > level else 0)
                | (4 if v_tr > level else 0)
                | (8 if v_tl > level else 0)
            )
            if case in (0, 15):
                continue
            pairs = _EDGE_PAIRS[case]
            if case in (5, 10):
                center = 0.25 * (v_bl + v_br + v_tr + v_tl)
                # Standard asymptotic-decider tie-break: if the center lies
                # on the same side as the bottom-left corner, use the
                # "non-crossing" diagonal resolution; otherwise the other.
                same_side_as_bl = (center > level) == (v_bl > level)
                pairs = [(3, 0), (1, 2)] if same_side_as_bl else [(3, 2), (1, 0)]

            edge_points = {
                0: _bottom(i0, j0, j1),
                1: _right(i0, i1, j1),
                2: _top(i1, j0, j1),
                3: _left(i0, i1, j0),
            }
            for a, b in pairs:
                segments.append((edge_points[a], edge_points[b]))

    return _stitch_segments(segments)


def _stitch_segments(
    segments: list[tuple[tuple[float, float], tuple[float, float]]],
) -> ContourSet:
    """Link raw line segments into polylines by matching shared endpoints.

    Drops zero-length segments first: these arise when a grid VERTEX value
    exactly equals the level (a measure-zero coincidence, but common on
    synthetic test grids -- e.g. an axis-aligned point landing exactly on a
    circular level set), which the ``>`` corner classification treats as an
    edge case and can emit spurious self-loop segments for.
    """
    segments = [(a, b) for a, b in segments if a != b]

    adjacency: dict[tuple[float, float], list[tuple[float, float]]] = {}
    for a, b in segments:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    visited_edges: set[frozenset[tuple[float, float]]] = set()
    polylines: list[NDArray[np.float64]] = []
    closed_flags: list[bool] = []

    for start in list(adjacency.keys()):
        for neighbor in list(adjacency[start]):
            edge_key = frozenset((start, neighbor))
            if edge_key in visited_edges:
                continue
            # Walk a chain starting from this edge.
            chain = [start, neighbor]
            visited_edges.add(edge_key)
            current = neighbor
            while True:
                nxt = None
                for cand in adjacency.get(current, []):
                    ek = frozenset((current, cand))
                    if ek not in visited_edges:
                        nxt = cand
                        break
                if nxt is None:
                    break
                visited_edges.add(frozenset((current, nxt)))
                chain.append(nxt)
                current = nxt
                if current == chain[0]:
                    break
            arr = np.array(chain, dtype=np.float64)
            is_closed = bool(np.allclose(arr[0], arr[-1]))
            polylines.append(arr)
            closed_flags.append(is_closed)

    return ContourSet(polylines=polylines, closed=closed_flags)


__all__ = ["ContourSet", "marching_squares_contours"]
