"""Kustaanheimo-Stiefel (KS) map: the four-dimensional regularising coordinates (#928).

Convention: Kustaanheimo & Stiefel 1965 (J. reine angew. Math. 218:204, eqs. (5)-(17), (28),
(37)) and Stiefel & Scheifele 1971 (Linear and Regular Celestial Mechanics, eqs. (9,27)-(9,34),
(9,65), (9,69)-(9,71)); digests ``docs/notes/2026-10-05-digest-kustaanheimo-stiefel-1965-ks-
regularization.md`` (section 7) and ``docs/notes/2026-10-04-digest-stiefel-scheifele-1971-
linear-regular-celestial-mechanics.md``. All five held sources print the same matrix

    L(u) = [[u1, -u2, -u3,  u4],
            [u2,  u1, -u4, -u3],
            [u3,  u4,  u1,  u2],
            [u4, -u3,  u2, -u1]]

(Peters 1968's A* and Aarseth 1971's L^T are L3^T; Aarseth & Zare 1974's A_1 is 2 L3^T, where
L3 is the first three rows). With a physical 3-vector padded by a fourth component 0:

* position ``x = L(u) u`` (fourth component identically 0), ``|x| = u.u``;
* differential ``dx = 2 L3(u) du`` when ``du`` obeys the bilinear relation;
* regularising time ``dt = r ds`` (r = u.u) and ``u' = du/ds = (1/2) L3(u)^T xdot``;
* physical velocity ``xdot = 2 L3(u) u' / r``;
* bilinear relation (KS (7), S&S (9,34)) ``l(u, u') = u4 u1' - u3 u2' + u2 u3' - u1 u4' = 0``;
* generalised force ``Q = 2 L3(u)^T P`` (KS (28)).

The map is two-to-one per point of a circle: ``u`` and its rotations along the fibre (KS (14))
give the same ``x``. Planar data (x3 = 0) lift to ``u3 = u4 = 0`` on either inverse branch with
the default fibre angle, which is the Levi-Civita map (S&S p.35).

Pure numpy; no dynamics here (see :mod:`cyclerfinder.core.cr3bp_ks`).
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

# E[k] = dL/du_k (L is linear in u): L(u) = sum_k u_k E[k].
_E = np.zeros((4, 4, 4), dtype=np.float64)
# row 0: ( u1, -u2, -u3,  u4)
_E[0, 0, 0], _E[1, 0, 1], _E[2, 0, 2], _E[3, 0, 3] = 1.0, -1.0, -1.0, 1.0
# row 1: ( u2,  u1, -u4, -u3)
_E[1, 1, 0], _E[0, 1, 1], _E[3, 1, 2], _E[2, 1, 3] = 1.0, 1.0, -1.0, -1.0
# row 2: ( u3,  u4,  u1,  u2)
_E[2, 2, 0], _E[3, 2, 1], _E[0, 2, 2], _E[1, 2, 3] = 1.0, 1.0, 1.0, 1.0
# row 3: ( u4, -u3,  u2, -u1)
_E[3, 3, 0], _E[2, 3, 1], _E[1, 3, 2], _E[0, 3, 3] = 1.0, -1.0, 1.0, -1.0
_E.setflags(write=False)
KS_BASIS: FloatArray = _E
"""``KS_BASIS[k] = dL/du_k``; ``ks_matrix(u) = sum_k u[k] * KS_BASIS[k]``."""


def ks_matrix(u: ArrayLike) -> FloatArray:
    """The 4 x 4 KS matrix L(u) (KS 1965 eq. (5) matrix A; S&S (9,27))."""
    u1, u2, u3, u4 = (float(c) for c in np.asarray(u, dtype=np.float64))
    return np.array(
        [
            [u1, -u2, -u3, u4],
            [u2, u1, -u4, -u3],
            [u3, u4, u1, u2],
            [u4, -u3, u2, -u1],
        ],
        dtype=np.float64,
    )


def ks_matrix3(u: ArrayLike) -> FloatArray:
    """The 3 x 4 block L3(u) (first three rows of L): ``dx = 2 L3 du``."""
    return ks_matrix(u)[:3]


def ks_position(u: ArrayLike) -> FloatArray:
    """Physical position ``x = L3(u) u`` (KS (6); S&S (9,29))."""
    u1, u2, u3, u4 = (float(c) for c in np.asarray(u, dtype=np.float64))
    return np.array(
        [
            u1 * u1 - u2 * u2 - u3 * u3 + u4 * u4,
            2.0 * (u1 * u2 - u3 * u4),
            2.0 * (u1 * u3 + u2 * u4),
        ],
        dtype=np.float64,
    )


def ks_velocity(u: ArrayLike, up: ArrayLike) -> FloatArray:
    """Physical velocity ``xdot = 2 L3(u) u' / |u|^2`` (S&S (9,65)); u' = du/ds, dt = r ds.

    Undefined at ``u = 0`` (collision); raises ``ZeroDivisionError`` there.
    """
    u_arr = np.asarray(u, dtype=np.float64)
    r = float(u_arr @ u_arr)
    if r == 0.0:
        raise ZeroDivisionError("ks_velocity: the physical velocity is undefined at u = 0")
    return 2.0 * (ks_matrix3(u_arr) @ np.asarray(up, dtype=np.float64)) / r


def ks_velocity_lift(u: ArrayLike, xdot: ArrayLike) -> FloatArray:
    """Regularised velocity ``u' = (1/2) L3(u)^T xdot`` (KS (37); S&S (9,71)).

    The result satisfies the bilinear relation with ``u`` exactly (up to rounding).
    """
    return 0.5 * (ks_matrix3(u).T @ np.asarray(xdot, dtype=np.float64))


def ks_bilinear(u: ArrayLike, up: ArrayLike) -> float:
    """``l(u, u') = u4 u1' - u3 u2' + u2 u3' - u1 u4'`` (KS (7), S&S (9,34)); zero on motions.

    It equals the fourth component of ``L(u) u'``.
    """
    u1, u2, u3, u4 = (float(c) for c in np.asarray(u, dtype=np.float64))
    w1, w2, w3, w4 = (float(c) for c in np.asarray(up, dtype=np.float64))
    return u4 * w1 - u3 * w2 + u2 * w3 - u1 * w4


def ks_fibre_rotate(u: ArrayLike, phi: float) -> FloatArray:
    """Rotate ``u`` by ``phi`` along its fibre (KS 1965 eq. (14)); the image ``x`` is unchanged.

    v1 = u1 cos - u4 sin, v2 = u2 cos + u3 sin, v3 = -u2 sin + u3 cos, v4 = u1 sin + u4 cos.
    The map is linear, so the same rotation applied to ``u'`` keeps the pair on a motion.
    """
    u1, u2, u3, u4 = (float(c) for c in np.asarray(u, dtype=np.float64))
    c, s = math.cos(phi), math.sin(phi)
    return np.array(
        [u1 * c - u4 * s, u2 * c + u3 * s, -u2 * s + u3 * c, u1 * s + u4 * c], dtype=np.float64
    )


def ks_fibre_tangent(u: ArrayLike) -> FloatArray:
    """Tangent of the fibre circle through ``u``: ``(-u4, u3, -u2, u1)`` (KS (15))."""
    u1, u2, u3, u4 = (float(c) for c in np.asarray(u, dtype=np.float64))
    return np.array([-u4, u3, -u2, u1], dtype=np.float64)


def ks_lift(x: ArrayLike, phi: float = 0.0, branch: str = "auto") -> FloatArray:
    """A KS preimage ``u`` of the physical position ``x`` (KS (16); S&S (9,69), (9,70)).

    ``branch``: ``"pos"`` uses u1^2 + u4^2 = (r + x1)/2 (S&S (9,69), stable for x1 >= 0),
    ``"neg"`` uses u2^2 + u3^2 = (r - x1)/2 ((9,70), stable for x1 < 0), ``"auto"`` picks by the
    sign of x1. ``phi`` is the free angle of the first line of each branch: on ``"pos"``
    (u1, u4) = sqrt((r + x1)/2) (cos phi, sin phi), on ``"neg"`` (u2, u3) = sqrt((r - x1)/2)
    (cos phi, sin phi). ``phi = 0`` keeps planar data planar (u3 = u4 = 0). ``x = 0`` lifts to
    ``u = 0``.
    """
    x1, x2, x3 = (float(c) for c in np.asarray(x, dtype=np.float64))
    r = math.sqrt(x1 * x1 + x2 * x2 + x3 * x3)
    if r == 0.0:
        return np.zeros(4, dtype=np.float64)
    if branch == "auto":
        branch = "pos" if x1 >= 0.0 else "neg"
    c, s = math.cos(phi), math.sin(phi)
    if branch == "pos":
        rho = math.sqrt(0.5 * (r + x1))
        u1, u4 = rho * c, rho * s
        d = r + x1
        u2 = (x2 * u1 + x3 * u4) / d
        u3 = (x3 * u1 - x2 * u4) / d
    elif branch == "neg":
        rho = math.sqrt(0.5 * (r - x1))
        u2, u3 = rho * c, rho * s
        d = r - x1
        u1 = (x2 * u2 + x3 * u3) / d
        u4 = (x3 * u2 - x2 * u3) / d
    else:
        raise ValueError(f"ks_lift: branch must be 'auto', 'pos' or 'neg', got {branch!r}")
    return np.array([u1, u2, u3, u4], dtype=np.float64)


def ks_generalised_force(u: ArrayLike, force: ArrayLike) -> FloatArray:
    """Generalised force ``Q = 2 L3(u)^T P`` of a physical force ``P`` (KS 1965 eq. (28))."""
    return 2.0 * (ks_matrix3(u).T @ np.asarray(force, dtype=np.float64))
