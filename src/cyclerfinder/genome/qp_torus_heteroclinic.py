"""Linking-number heteroclinic-connection screen for two QP-tori (#522).

Wires together the pieces built for #522 (Owen & Baresi 2024's method — see
``docs/notes/2026-07-03-digest-owen-baresi-2024-knot-theory-heteroclinic.md``):

  1. ``genome/qp_torus_manifold.py``: per-point torus stability + manifold
     endpoint grids at a surface of section (CR3BP-specific).
  2. ``search/torus_map_contours.py``: marching-squares level-curve
     extraction on the manifold-endpoint grid, treated as a torus map
     (generic).
  3. ``search/linking_number.py``: linking number of the resulting reduced
     closed curves (generic).

Scanning a state-component value ``D`` and tracking where the linking
number CHANGES flags candidate heteroclinic-connection locations — a cheap
topological screen, not a converged trajectory (Owen & Baresi's own method
explicitly requires a differential-correction step afterward, which this
module does not implement; #522 Phase 1 is the screen only).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.genome.qp_tori import QPTorus
from cyclerfinder.genome.qp_torus_manifold import ManifoldGrid, torus_manifold_grid
from cyclerfinder.search.linking_number import linking_number
from cyclerfinder.search.torus_map_contours import marching_squares_contours

# State-vector component names, matching Owen & Baresi's own (x,y,z,xdot,
# ydot,zdot) convention (Sec 2.2, their Eq 2).
_STATE_NAMES = ("x", "y", "z", "xdot", "ydot", "zdot")


def _component_index(name: str) -> int:
    if name not in _STATE_NAMES:
        raise ValueError(f"unknown state component {name!r}; must be one of {_STATE_NAMES}")
    return _STATE_NAMES.index(name)


def _reduced_closed_curve(
    grid: ManifoldGrid,
    contour_grid_idx: NDArray[np.float64],
    component_names: tuple[str, str, str],
) -> NDArray[np.float64] | None:
    """Interpolate 3 state components at fractional grid-index points along
    a level-curve polyline, using bilinear interpolation on ``grid.endpoints``.

    Returns ``None`` if any sample point falls in a NaN (non-crossing)
    region of the grid -- the curve is not usable.
    """
    n_long, n_lat = grid.hyperbolic.shape
    idx = [_component_index(c) for c in component_names]
    out = np.zeros((contour_grid_idx.shape[0], 3), dtype=np.float64)
    for row, (fi, fj) in enumerate(contour_grid_idx):
        i0 = int(np.floor(fi)) % n_long
        i1 = (i0 + 1) % n_long
        j0 = int(np.floor(fj)) % n_lat
        j1 = (j0 + 1) % n_lat
        ti = fi - np.floor(fi)
        tj = fj - np.floor(fj)
        for k, comp in enumerate(idx):
            v00 = grid.endpoints[i0, j0, comp]
            v01 = grid.endpoints[i0, j1, comp]
            v10 = grid.endpoints[i1, j0, comp]
            v11 = grid.endpoints[i1, j1, comp]
            if not (
                np.isfinite(v00) and np.isfinite(v01) and np.isfinite(v10) and np.isfinite(v11)
            ):
                return None
            v0 = v00 * (1 - tj) + v01 * tj
            v1 = v10 * (1 - tj) + v11 * tj
            out[row, k] = v0 * (1 - ti) + v1 * ti
    return out


@dataclass(frozen=True)
class LinkingScanResult:
    """Linking-number evolution over a scan of the ``D`` scanning variable."""

    d_values: NDArray[np.float64]
    linking_numbers: NDArray[np.int64]

    def sign_change_locations(self) -> list[float]:
        """Midpoint ``D`` values where the linking number changed -- initial
        guesses for heteroclinic-connection locations (Owen & Baresi Sec 3.2:
        "an average of the D values before and after the change ... taken to
        be an initial guess").
        """
        changes = []
        for i in range(1, len(self.linking_numbers)):
            if self.linking_numbers[i] != self.linking_numbers[i - 1]:
                changes.append(0.5 * (float(self.d_values[i - 1]) + float(self.d_values[i])))
        return changes


def scan_linking_number(
    stable_grid: ManifoldGrid,
    unstable_grid: ManifoldGrid,
    *,
    scanning_component: str,
    curve_components: tuple[str, str, str],
    d_values: NDArray[np.float64],
) -> LinkingScanResult:
    """Sweep ``scanning_component`` over ``d_values``, computing the linking
    number of the stable/unstable reduced closed curves at each value.

    A ``D`` value where either manifold's torus map has no valid level
    curve (e.g. the level lies outside the map's range, or the region is
    NaN-filled from non-crossing trajectories) contributes linking number
    ``0`` at that ``D`` -- consistent with "no connection detectable there
    from this data", not asserted as a genuine absence.
    """
    d_idx = _component_index(scanning_component)
    n_long_s, n_lat_s = stable_grid.hyperbolic.shape
    n_long_u, n_lat_u = unstable_grid.hyperbolic.shape

    field_s = stable_grid.endpoints[:, :, d_idx]
    field_u = unstable_grid.endpoints[:, :, d_idx]

    linking_numbers = np.zeros(len(d_values), dtype=np.int64)
    for k, d in enumerate(d_values):
        curve_s = _first_closed_curve(
            field_s, float(d), n_long_s, n_lat_s, stable_grid, curve_components
        )
        curve_u = _first_closed_curve(
            field_u, float(d), n_long_u, n_lat_u, unstable_grid, curve_components
        )
        if curve_s is None or curve_u is None:
            linking_numbers[k] = 0
            continue
        linking_numbers[k] = linking_number(curve_s, curve_u)

    return LinkingScanResult(
        d_values=np.asarray(d_values, dtype=np.float64), linking_numbers=linking_numbers
    )


def _first_closed_curve(
    field: NDArray[np.float64],
    level: float,
    n_long: int,
    n_lat: int,
    grid: ManifoldGrid,
    curve_components: tuple[str, str, str],
) -> NDArray[np.float64] | None:
    if not np.any(np.isfinite(field)):
        return None
    field_filled = np.where(np.isfinite(field), field, np.nan)
    if np.all(np.isnan(field_filled)):
        return None
    finite_vals = field_filled[np.isfinite(field_filled)]
    if level < finite_vals.min() or level > finite_vals.max():
        return None
    # marching_squares_contours does not itself understand NaN; substitute
    # a value far outside the level range so NaN cells never falsely cross.
    safe_field = np.where(np.isfinite(field_filled), field_filled, level + 1e12)
    contours = marching_squares_contours(safe_field, level, periodic=(True, True))
    for poly, is_closed in zip(contours.polylines, contours.closed, strict=True):
        if not is_closed or poly.shape[0] < 4:
            continue
        curve = _reduced_closed_curve(grid, poly, curve_components)
        if curve is not None:
            return curve
    return None


def build_manifold_grids(
    stable_torus: QPTorus,
    unstable_torus: QPTorus,
    *,
    n_long: int,
    n_lat: int,
    eps: float,
    surface_x: float,
    t_max: float,
    stable_sign: float = 1.0,
    unstable_sign: float = 1.0,
) -> tuple[ManifoldGrid, ManifoldGrid]:
    """Convenience wrapper: stable manifold of ``stable_torus``, unstable
    manifold of ``unstable_torus``, both to the same surface of section.
    """
    stable_grid = torus_manifold_grid(
        stable_torus,
        n_long=n_long,
        n_lat=n_lat,
        branch="stable",
        eps=eps,
        sign=stable_sign,
        surface_x=surface_x,
        t_max=t_max,
    )
    unstable_grid = torus_manifold_grid(
        unstable_torus,
        n_long=n_long,
        n_lat=n_lat,
        branch="unstable",
        eps=eps,
        sign=unstable_sign,
        surface_x=surface_x,
        t_max=t_max,
    )
    return stable_grid, unstable_grid


__all__ = [
    "LinkingScanResult",
    "build_manifold_grids",
    "scan_linking_number",
]
