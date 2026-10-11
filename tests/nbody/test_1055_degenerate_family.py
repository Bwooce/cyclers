"""#1055: degenerate-family seed correction, on analytic toys at production scales (km near
1e6, km/s near 3, weights as in the Jovian lane: velocity rows x 1e3). Expected values are the
toys' analytic roots, set by construction (not computed by the code under test)."""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.nbody.degenerate_family import (
    bordered_solve,
    column_scales,
    detect_degenerate,
    family_root_1d,
    family_root_nd,
    scaled_svd,
)

W_ROOT = 4.0e5  # km: the solvability root along the family, O(1) away from the seed (w = 0)
U_ROOT = -2.0e5


def toy1(sigma: float, degenerate: bool = True):  # type: ignore[no-untyped-def]
    """z = (x km, vy km/s, w km). Rows: x fixed up to a curvature in w; vy fixed up to a slope in w;
    the family row is sigma-proportional (degenerate) or not."""
    f = sigma if degenerate else 1.0

    def rj(z):  # type: ignore[no-untyped-def]
        x, vy, w = z
        r = np.array(
            [
                x - 1.0e6 - 1e-9 * w**2,
                1e3 * (vy - 3.0 - 1e-6 * w),
                f * (1e3 * (w - W_ROOT) / 1e5 + 1e-2 * (vy - 3.0 - 1e-6 * w)),
            ]
        )
        j = np.array(
            [
                [1.0, 0.0, -2e-9 * w],
                [0.0, 1e3, -1e-3],
                [0.0, f * 1e-2, f * (1e3 / 1e5 - 1e-8)],
            ]
        )
        return r, j

    return rj


def toy2(sigma: float):  # type: ignore[no-untyped-def]
    """Two degenerate directions (w, u km), coupled."""

    def rj(z):  # type: ignore[no-untyped-def]
        x, vy, w, u = z
        r = np.array(
            [
                x - 1.0e6 - 1e-9 * w**2 - 5e-10 * u**2,
                1e3 * (vy - 3.0 - 1e-6 * w + 2e-6 * u),
                sigma * 1e3 * ((w - W_ROOT) + 0.2 * (u - U_ROOT)) / 1e5,
                sigma * 1e3 * ((u - U_ROOT) - 0.1 * (w - W_ROOT)) / 1e5,
            ]
        )
        j = np.array(
            [
                [1.0, 0.0, -2e-9 * w, -1e-9 * u],
                [0.0, 1e3, -1e-3, 2e-3],
                [0.0, 0.0, sigma * 1e-2, sigma * 0.2e-2],
                [0.0, 0.0, -sigma * 0.1e-2, sigma * 1e-2],
            ]
        )
        return r, j

    return rj


def test_detect_degenerate_one() -> None:
    z = np.array([1.0e6, 3.0, 0.0])
    assert detect_degenerate(toy1(0.02)(z)[1], toy1(0.005)(z)[1], 4.0) == 1


def test_detect_degenerate_two() -> None:
    z = np.array([1.0e6, 3.0, 0.0, 0.0])
    assert detect_degenerate(toy2(0.02)(z)[1], toy2(0.005)(z)[1], 4.0) == 2


def test_detect_well_conditioned_is_zero() -> None:
    z = np.array([1.0e6, 3.0, 0.0])
    j_hi = toy1(0.02, degenerate=False)(z)[1]
    j_lo = toy1(0.005, degenerate=False)(z)[1]
    assert detect_degenerate(j_hi, j_lo, 4.0) == 0


def test_bordered_solve_holds_arclength_row() -> None:
    rj = toy1(0.02)
    z0 = np.array([1.0e6, 3.0, 0.0])
    _, j0 = rj(z0)
    d = column_scales(j0)
    u, _, vt = scaled_svd(j0, d)
    alpha = np.array([0.37])
    br = bordered_solve(rj, z0 + vt[-1] * 0.37 / d, d, u[:, -1:], vt[-1:].T, z0, alpha, 1e-6)
    assert br.converged
    assert abs(float(vt[-1] @ (d * (br.z - z0))) - 0.37) < 1e-12
    assert br.complement_norm < 1e-6


def test_family_root_1d_finds_the_solvability_root() -> None:
    rj = toy1(0.02)
    z0 = np.array([1.0e6, 3.0, 0.0])
    res = family_root_1d(rj, z0, tol=1e-6)
    assert res.status == "root", res.status
    assert res.z_star is not None
    w = res.z_star[2]
    # the secant stops at |g| < 1e-3 |g(0)|: w within 1e-3 of the seed's distance from the root
    assert abs(w - W_ROOT) < 1e-3 * W_ROOT
    # hand-off: plain Newton from the corrected seed reaches the exact root
    z = res.z_star.copy()
    for _ in range(8):
        r, j = rj(z)
        z = z + np.linalg.solve(j, -r)
    assert abs(z[2] - W_ROOT) < 1e-9 * W_ROOT
    assert abs(z[0] - (1.0e6 + 1e-9 * W_ROOT**2)) < 1e-6


def test_family_root_nd_two_directions() -> None:
    rj = toy2(0.02)
    z0 = np.array([1.0e6, 3.0, 0.0, 0.0])
    res = family_root_nd(rj, z0, 2, tol=1e-6)
    assert res.status == "root", res.status
    assert res.z_star is not None
    assert res.z_star[2] == pytest.approx(W_ROOT, rel=2e-3)
    assert res.z_star[3] == pytest.approx(U_ROOT, rel=2e-3)


def test_well_conditioned_seed_needs_no_move() -> None:
    """A seed already at the root: g(0) = 0 within rounding, the root is the base."""
    rj = toy1(0.02)
    z0 = np.array([1.0e6 + 1e-9 * W_ROOT**2, 3.0 + 1e-6 * W_ROOT, W_ROOT])
    res = family_root_1d(rj, z0, tol=1e-6)
    assert res.status == "root at base"
    assert res.z_star is not None
    assert abs(res.z_star[2] - W_ROOT) < 1e-6 * W_ROOT


def test_family_frame_converges_an_off_family_seed() -> None:
    """A seed off the family (vy 10 m/s high): the frame's base has its complement converged (the
    residual lies along the left null vector) and is linearised there (amendment 4)."""
    from cyclerfinder.nbody.degenerate_family import family_frame

    rj = toy1(0.02)
    fr = family_frame(rj, np.array([1.0e6, 3.01, 0.0]), 1, 1e-6)
    assert fr.converged
    assert fr.base.complement_norm < 1e-6
    r, _ = rj(fr.base.z)
    u = fr.u_n[:, 0]
    assert float(np.max(np.abs(r - u * float(u @ r)))) < 1e-6  # the residual lies along u
