"""#928 stage 1: the KS map against every held printed source, before any dynamics.

Test plan: Kustaanheimo-Stiefel 1965 digest
(``docs/notes/2026-10-05-digest-kustaanheimo-stiefel-1965-ks-regularization.md``) section 7. The
printed matrices are written out below entry by entry as each source prints them (not built from
the code), with the documented scale factor:

* KS 1965 eq. (5) (J. reine angew. Math. 218, p.205): the 4 x 4 matrix A = L(u).
* Stiefel & Scheifele 1971 eq. (9,27) (p.24): the same L(u).
* Peters 1968 A* (4 x 3): equals L3^T.
* Aarseth 1971 eq. (34) L^T (4 x 3): equals L3^T.
* Aarseth & Zare 1974 eq. (51) A_1 (4 x 3): equals 2 L3^T.

Numerical golden: Stiefel & Scheifele 1971 section 15, Example 2 (p.71, read on the page image
by the digest): the inverse KS map with u4 = 0 of the Example 1 data. Printed to six digits from
1971 hand arithmetic; the digest reproduced it to about 2e-6, hence the 5e-6 tolerance.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from cyclerfinder.core.ks import (
    KS_BASIS,
    ks_bilinear,
    ks_fibre_rotate,
    ks_fibre_tangent,
    ks_generalised_force,
    ks_lift,
    ks_matrix,
    ks_matrix3,
    ks_position,
    ks_velocity,
    ks_velocity_lift,
)

FloatArray = NDArray[np.float64]


def _ks1965_a(u: FloatArray) -> FloatArray:
    """KS 1965 eq. (5) / S&S (9,27), as printed."""
    u1, u2, u3, u4 = u
    return np.array([[u1, -u2, -u3, u4], [u2, u1, -u4, -u3], [u3, u4, u1, u2], [u4, -u3, u2, -u1]])


def _peters_a_star(u: FloatArray) -> FloatArray:
    """Peters 1968 A* (4 x 3), as printed (digest of Peters, section 1.2)."""
    u1, u2, u3, u4 = u
    return np.array([[u1, u2, u3], [-u2, u1, u4], [-u3, -u4, u1], [u4, -u3, u2]])


def _aarseth_lt(u: FloatArray) -> FloatArray:
    """Aarseth 1971 eq. (34) L^T (4 x 3), as printed."""
    u1, u2, u3, u4 = u
    return np.array([[u1, u2, u3], [-u2, u1, u4], [-u3, -u4, u1], [u4, -u3, u2]])


def _aarseth_zare_a1(q: FloatArray) -> FloatArray:
    """Aarseth & Zare 1974 eq. (51) A_1, as printed (with its factor 2)."""
    q1, q2, q3, q4 = q
    return 2.0 * np.array([[q1, q2, q3], [-q2, q1, q4], [-q3, -q4, q1], [q4, -q3, q2]])


def _draws(n: int = 8) -> list[tuple[FloatArray, FloatArray, FloatArray]]:
    rng = np.random.default_rng(1965)
    out: list[tuple[FloatArray, FloatArray, FloatArray]] = []
    while len(out) < n:
        u = rng.normal(size=4)
        if np.linalg.norm(u) < 1e-3:
            continue
        out.append((u, rng.normal(size=3), rng.normal(size=3)))
    return out


@pytest.mark.parametrize("k", range(8))
def test_map_matches_every_printed_source(k: int) -> None:
    u, xdot, force = _draws()[k]
    r = float(u @ u)
    x = ks_position(u)
    m3 = ks_matrix3(u)
    # 1. polynomial form (KS (6), S&S (9,29)) and the matrix form agree; |x| = u.u.
    u1, u2, u3, u4 = u
    poly = np.array(
        [u1**2 - u2**2 - u3**2 + u4**2, 2 * (u1 * u2 - u3 * u4), 2 * (u1 * u3 + u2 * u4)]
    )
    assert np.allclose(x, poly, rtol=0, atol=1e-14 * r)
    assert np.allclose(m3 @ u, x, rtol=0, atol=1e-14 * r)
    assert abs(float(ks_matrix(u)[3] @ u)) < 1e-14 * r  # fourth component identically zero
    assert abs(float(np.linalg.norm(x)) - r) < 1e-14 * r
    # 2. the code's matrix is every source's matrix up to its documented factor.
    assert np.array_equal(ks_matrix(u), _ks1965_a(u))
    assert np.array_equal(m3.T, _peters_a_star(u))
    assert np.array_equal(m3.T, _aarseth_lt(u))
    assert np.array_equal(2.0 * m3.T, _aarseth_zare_a1(u))
    assert np.array_equal(sum(u[j] * KS_BASIS[j] for j in range(4)), ks_matrix(u))
    assert np.allclose(m3 @ m3.T, r * np.eye(3), rtol=0, atol=1e-14 * r)  # (10), (11)
    assert np.allclose(ks_matrix(u).T @ ks_matrix(u), r * np.eye(4), rtol=0, atol=1e-14 * r)
    # dx = 2 L3 du (KS (8)); x is quadratic, so the central difference is exact up to rounding.
    du = ks_velocity_lift(u, force) * 1e-3
    assert abs(ks_bilinear(u, du)) < 1e-15 * r * 1e-3 * 10
    exact = ks_position(u + du) - ks_position(u - du)
    assert np.allclose(exact, 4.0 * m3 @ du, rtol=0, atol=4e-15 * r)
    # and the inverse du = L3^T dx / (2 r) (KS (17)) recovers du, which obeys the control.
    assert np.allclose(m3.T @ (2.0 * m3 @ du) / (2.0 * r), du, rtol=0, atol=1e-15 * r)
    # 3. velocity: u' = (1/2) L3^T xdot satisfies the control and round-trips.
    up = ks_velocity_lift(u, xdot)
    scale = math.sqrt(r) * float(np.linalg.norm(up))
    assert abs(ks_bilinear(u, up)) <= 1e-15 * scale * 4
    assert abs(ks_bilinear(u, up) - float(ks_matrix(u)[3] @ up)) < 1e-15 * scale * 4
    assert np.allclose(ks_velocity(u, up), xdot, rtol=0, atol=1e-14 * np.linalg.norm(xdot))
    # KS (37) written out as printed.
    printed37 = 0.5 * np.array(
        [
            u1 * xdot[0] + u2 * xdot[1] + u3 * xdot[2],
            -u2 * xdot[0] + u1 * xdot[1] + u4 * xdot[2],
            -u3 * xdot[0] - u4 * xdot[1] + u1 * xdot[2],
            u4 * xdot[0] - u3 * xdot[1] + u2 * xdot[2],
        ]
    )
    assert np.allclose(up, printed37, rtol=0, atol=1e-15 * scale * 4)
    # 5. forces: Q = 2 L3^T P (KS (28)) = Peters P_u = 2 A* P; work invariance Q.u' = 2 P.(L3 u').
    q = ks_generalised_force(u, force)
    assert np.allclose(q, 2.0 * _peters_a_star(u) @ force, rtol=0, atol=1e-14 * r)
    assert np.allclose(q, _aarseth_zare_a1(u) @ force, rtol=0, atol=1e-14 * r)
    assert float(q @ up) == pytest.approx(2.0 * float(force @ (m3 @ up)), rel=1e-13, abs=1e-14)
    # the book's u-equation force (r/2) L3^T P equals r Q / 4 (KS (46) with W' = 0).
    assert np.allclose(0.5 * r * m3.T @ force, 0.25 * r * q, rtol=0, atol=1e-14 * r * r)


@pytest.mark.parametrize("k", range(8))
def test_inverse_branches_and_fibre(k: int) -> None:
    u, xdot, _ = _draws()[k]
    x = ks_position(u)
    r = float(u @ u)
    for branch in ("pos", "neg"):
        for phi in (0.0, 0.7, -2.9):
            v = ks_lift(x, phi=phi, branch=branch)
            assert np.allclose(ks_position(v), x, rtol=0, atol=1e-14 * r), branch
            assert float(v @ v) == pytest.approx(r, rel=1e-14)
    # the original u lies on the fibre of its own image: some rotation of the lift reaches it.
    v0 = ks_lift(x)
    tangent = ks_fibre_tangent(v0)
    phi = math.atan2(float(tangent @ u), float(v0 @ u))
    assert np.allclose(ks_fibre_rotate(v0, phi), u, rtol=0, atol=1e-13 * math.sqrt(r))
    # fibre rotation leaves x unchanged; the velocity lift is covariant along the fibre,
    # (1/2) L3(R u)^T xdot = R ((1/2) L3(u)^T xdot), and the rotated pair obeys the control.
    for ang in (0.3, 1.9, -4.0):
        ur = ks_fibre_rotate(u, ang)
        assert np.allclose(ks_position(ur), x, rtol=0, atol=1e-14 * r)
        upr = ks_velocity_lift(ur, xdot)
        assert np.allclose(upr, ks_fibre_rotate(ks_velocity_lift(u, xdot), ang), atol=1e-14 * r)
        assert abs(ks_bilinear(ur, upr)) < 1e-14 * r
        assert np.allclose(ks_velocity(ur, upr), xdot, atol=1e-13 * np.linalg.norm(xdot))


@pytest.mark.parametrize("x1", [0.8, -0.8, 0.0, -1e-9])
def test_planar_data_lift_to_levi_civita(x1: float) -> None:
    """Planar x lifts to u3 = u4 = 0 exactly on both branches (S&S p.35: Levi-Civita)."""
    x = np.array([x1, 0.37, 0.0])
    xdot = np.array([0.2, -1.1, 0.0])
    for branch in ("pos", "neg"):
        u = ks_lift(x, branch=branch)
        up = ks_velocity_lift(u, xdot)
        assert u[2] == 0.0 and u[3] == 0.0
        assert up[2] == 0.0 and up[3] == 0.0
        # Levi-Civita: x1 + i x2 = (u1 + i u2)^2.
        z = complex(u[0], u[1]) ** 2
        assert z.real == pytest.approx(x1, abs=1e-15) and z.imag == pytest.approx(0.37, rel=1e-15)


def test_stiefel_scheifele_section15_example2() -> None:
    """S&S 1971 section 15 Example 2 (p.71): inverse KS with u4 = 0 of the Example 1 state.

    Printed: u = (0.865440, 0.595975, 0.087008, 0), u' = (-0.282049, 0.413169, 0.145277,
    0.058505), r = (u,u) = 1.111743, (u,u') = 0.014782, (u',u') = 0.274788. Six-digit 1971 hand
    arithmetic: tolerance 5e-6 (digest section 5).
    """
    x = np.array([0.38623, 1.03156, 0.15060])
    xdot = np.array([-0.90484, 0.33171, 0.24476])
    u = ks_lift(x, phi=0.0, branch="pos")
    up = ks_velocity_lift(u, xdot)
    assert np.allclose(u, [0.865440, 0.595975, 0.087008, 0.0], rtol=0, atol=5e-6)
    assert np.allclose(up, [-0.282049, 0.413169, 0.145277, 0.058505], rtol=0, atol=5e-6)
    assert float(u @ u) == pytest.approx(1.111743, abs=5e-6)
    assert float(u @ up) == pytest.approx(0.014782, abs=5e-6)
    assert float(up @ up) == pytest.approx(0.274788, abs=5e-6)
    # the printed bilinear value -0.12e-6 is input rounding; in double precision it is 0.
    assert abs(ks_bilinear(u, up)) < 1e-16


def test_velocity_undefined_at_collision() -> None:
    with pytest.raises(ZeroDivisionError):
        ks_velocity(np.zeros(4), np.array([1.0, 0, 0, 0]))
    assert np.array_equal(ks_lift(np.zeros(3)), np.zeros(4))
    with pytest.raises(ValueError):
        ks_lift(np.ones(3), branch="other")
