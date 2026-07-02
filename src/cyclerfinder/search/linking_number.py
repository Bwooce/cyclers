"""Linking number of two closed curves in 3D (#522, Owen & Baresi 2024).

The linking number of two disjoint closed curves is a topological invariant
counting how many times they wind around each other; it can only change if
the curves pass through one another (Owen & Baresi 2024, Astrodynamics
8:577-595, Sec 2.1 -- see ``docs/notes/2026-07-03-digest-owen-baresi-2024-
knot-theory-heteroclinic.md``). Their Sec 3.1 detects it by fan-triangulating
a surface bounded by one curve (from its centroid) and counting SIGNED
crossings of the other curve's line segments through that surface.

This module reproduces that algorithm's STRUCTURE (fan triangulation +
segment/triangle crossing count) but uses the standard, well-established
Moller-Trumbore ray-triangle intersection test rather than attempting to
reproduce the paper's own inside-triangle formula verbatim -- the scanned
PDF's three-dot-product description does not match the standard barycentric
same-side test (which needs cross products dotted with the face normal, not
plain edge-to-point dot products) closely enough to trust a byte-for-byte
transcription; Moller-Trumbore is provably correct and is itself in the same
computer-graphics collision-detection family Owen & Baresi cite (ref 1,
Amanatides & Choi, "Ray tracing triangular meshes"). Validated here against
closed-form linked/unlinked curve pairs, not against the paper's own numeric
results (which are for CR3BP tori, not raw linking-number arithmetic).

Pure: numpy only, no CR3BP-specific code -- callers supply arbitrary closed
polylines in R^3.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _segment_triangle_intersection_t(
    seg_start: NDArray[np.float64],
    seg_dir: NDArray[np.float64],
    v0: NDArray[np.float64],
    v1: NDArray[np.float64],
    v2: NDArray[np.float64],
    *,
    eps: float = 1e-12,
) -> float | None:
    """Moller-Trumbore intersection of the segment ``[seg_start, seg_start +
    seg_dir]`` with the triangle ``(v0, v1, v2)``.

    Returns the segment parameter ``t in [0, 1]`` at the crossing, or
    ``None`` if the segment (not the infinite line) misses the triangle's
    interior, or is parallel to its plane.
    """
    edge1 = v1 - v0
    edge2 = v2 - v0
    h = np.cross(seg_dir, edge2)
    a = float(edge1 @ h)
    if abs(a) < eps:
        return None
    f = 1.0 / a
    s = seg_start - v0
    u = f * float(s @ h)
    if u < 0.0 or u > 1.0:
        return None
    q = np.cross(s, edge1)
    v = f * float(seg_dir @ q)
    if v < 0.0 or u + v > 1.0:
        return None
    t = f * float(edge2 @ q)
    if 0.0 <= t <= 1.0:
        return t
    return None


def fan_triangulate(
    curve: NDArray[np.float64],
) -> list[tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]]:
    """Fan-triangulate a closed curve ``(n, 3)`` from its centroid.

    Returns one triangle ``(p_j, p_{j+1}, centroid)`` per curve edge
    (``n`` triangles for ``n`` points, wrapping the last edge back to the
    first point) -- a surface homeomorphic to a disc bounded by ``curve``
    (Owen & Baresi Sec 3.1: "the surface constructed from the triangles
    would therefore be homeomorphic to a disc", which is all the linking-
    number argument needs, regardless of how distorted the fan looks).
    """
    n = curve.shape[0]
    centroid = curve.mean(axis=0)
    return [(curve[j], curve[(j + 1) % n], centroid) for j in range(n)]


def linking_number(
    curve_a: NDArray[np.float64],
    curve_b: NDArray[np.float64],
    *,
    eps: float = 1e-12,
) -> int:
    """Linking number of two closed polyline curves ``(n_a, 3)``, ``(n_b, 3)``.

    Fan-triangulates ``curve_a`` from its centroid, then counts SIGNED
    crossings of every ``curve_b`` line segment through every triangle (sign
    = whether the segment crosses in the direction of, or against, the
    triangle's own normal ``(v1-v0) x (v2-v0)``). The result is symmetric
    under swapping which curve is triangulated (a standard topological fact);
    this implementation always triangulates ``curve_a`` for a well-defined,
    deterministic call convention.

    Degenerate inputs (either curve touching the other's triangulated
    surface exactly, i.e. ``a`` very close to ``0`` in the underlying
    Moller-Trumbore solve) are treated as non-intersections at that
    triangle/segment pair via ``eps`` -- consistent with Owen & Baresi's own
    use of the linking number as a robust screen (Sec 3.1: "This method...
    has proven robust during our use, despite what form the linking curves
    take").

    Known sharp edge case (found while validating this implementation): if
    ``curve_b`` passes exactly through ``curve_a``'s CENTROID (the shared
    apex vertex of every fan triangle), the crossing can register against
    several adjacent triangles at once and over-count. Real manifold-derived
    curves essentially never hit this measure-zero coincidence exactly, but
    a caller feeding synthetic/symmetric test curves should avoid centering
    one curve's piercing point exactly on the other's centroid.
    """
    triangles = fan_triangulate(curve_a)
    n_b = curve_b.shape[0]
    total = 0
    for k in range(n_b):
        q_k = curve_b[k]
        seg_dir = curve_b[(k + 1) % n_b] - q_k
        for v0, v1, v2 in triangles:
            t = _segment_triangle_intersection_t(q_k, seg_dir, v0, v1, v2, eps=eps)
            if t is None:
                continue
            normal = np.cross(v1 - v0, v2 - v0)
            total += 1 if float(seg_dir @ normal) > 0.0 else -1
    return total


__all__ = ["fan_triangulate", "linking_number"]
