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


def test_family_root_1d_resumes_from_its_state() -> None:
    """Calls with an already-past deadline still progress one step each; the saved state carries
    the run to the same root as an uninterrupted run (amendment 6)."""
    import json

    rj = toy1(0.02)
    z0 = np.array([1.0e6, 3.0, 0.0])
    full = family_root_1d(rj, z0, tol=1e-6)
    state: dict = {}  # type: ignore[type-arg]
    res = None
    for _ in range(200):
        state = json.loads(json.dumps(state))  # through JSON, as the driver stores it
        res = family_root_1d(rj, z0, tol=1e-6, state=state, deadline=0.0)
        if res.status != "paused":
            break
    assert res is not None and res.status == "root"
    assert res.z_star is not None and full.z_star is not None
    assert abs(res.z_star[2] - full.z_star[2]) < 1e-6 * W_ROOT


R_CIRCLE = 1.0e5  # km
THETA_ROOT = 2.0  # rad: 115 deg round the circle from the seed at theta = 0


def toy_circle(sigma: float):  # type: ignore[no-untyped-def]
    """A strongly curved family: the circle |(x, y)| = R (km), with a sigma-proportional row whose
    root is at theta = THETA_ROOT. A straight (chord) scan from theta = 0 leaves the circle; the
    pseudo-arclength scan must follow it round (amendment 7)."""

    def rj(z):  # type: ignore[no-untyped-def]
        x, y = z
        rr = float(np.hypot(x, y))
        th = float(np.arctan2(y, x))
        r = np.array([rr - R_CIRCLE, sigma * 1e3 * (th - THETA_ROOT)])
        j = np.array([[x / rr, y / rr], [sigma * 1e3 * (-y / rr**2), sigma * 1e3 * (x / rr**2)]])
        return r, j

    return rj


def test_family_root_1d_follows_a_curved_family() -> None:
    rj = toy_circle(0.02)
    # in scaled units the arc is dominated by x (scale 1): the range is not the point here
    res = family_root_1d(rj, np.array([R_CIRCLE, 0.0]), tol=1e-6, max_lengths=1e5)
    assert res.status == "root", res.status
    assert res.z_star is not None
    th = float(np.arctan2(res.z_star[1], res.z_star[0]))
    assert abs(th - THETA_ROOT) < 1e-3 * THETA_ROOT
    assert abs(float(np.hypot(*res.z_star)) - R_CIRCLE) < 1e-3


def test_family_root_1d_stops_at_a_domain_bound() -> None:
    """A domain bound (here theta < 1 rad) ends the scan in that direction: no root, no crash."""
    rj = toy_circle(0.02)
    res = family_root_1d(
        rj,
        np.array([R_CIRCLE, 0.0]),
        tol=1e-6,
        domain_ok=lambda z: float(np.arctan2(z[1], z[0])) < 1.0,
    )
    assert res.status == "no sign change of g within the scanned range"


def test_weak_direction_move_at_jovian_scale() -> None:
    """A near-degenerate direction with an ephemeris-like fixed stiffness (eps) and a forcing
    proportional to sigma: the implied move along it is c sigma / eps, so it shrinks in proportion
    to sigma (the gc-1 e2 reading, sec. 13.1). Node offsets at km scale (7 unknowns per node)."""
    from cyclerfinder.nbody.degenerate_family import weak_direction_move

    eps, c, n = 2.7e-7, 1.0e-3, 7
    q, _ = np.linalg.qr(np.random.default_rng(1055).normal(size=(n, n)))
    j = q @ np.diag([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, eps]) @ q.T  # columns of norm ~1: no rescaling

    def system(sigma: float) -> tuple[np.ndarray, np.ndarray]:
        z = np.zeros(n)
        z[0:3] = sigma * np.array([2.0e4, 1.0e4, 5.0e3])  # node offset, km (proportional to sigma)
        r = c * sigma * q[:, -1]  # O(sigma) forcing along the weak left vector
        return r, z

    moves, ratios = [], []
    for sigma in (0.02, 0.005):
        r, z = system(sigma)
        wm = weak_direction_move(r, j, z, n_nodes=1)
        assert wm.sv_min < 1e-5 < wm.sv_next  # the weak direction survives column scaling
        moves.append(float(np.linalg.norm(wm.dz)))
        ratios.append(wm.move_over_offset[0])
    # the move is c sigma / eps up to the column scaling, so it is proportional to sigma
    assert abs(moves[0] / moves[1] - 4.0) < 1e-6
    assert abs(moves[0] / (c * 0.02 / eps) - 1.0) < 0.2  # c sigma / eps up to the column scales
    assert abs(ratios[0] / ratios[1] - 1.0) < 1e-6  # offsets scale with sigma too


def test_forced_scan_finds_the_root_of_a_near_degenerate_family() -> None:
    """Amendment 13: a weak direction of FIXED stiffness eps (it does not scale with sigma, so the
    detection gives q = 0) with an O(sigma) forcing. Scanning it by declaration (the 1-d family
    root) finds the member where the projected residual vanishes, w* = -c sigma / eps (km), set by
    construction."""
    eps, c, sigma = 2.7e-7, 1.0e-3, 0.02
    w_star = -c * sigma / eps  # -74.07 km

    def rj(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x, vy, w = z
        r = np.array(
            [
                x - 1.0e6 - 1e-9 * w**2,
                1e3 * (vy - 3.0 - 1e-6 * w),
                eps * w + c * sigma + 1e-2 * (vy - 3.0 - 1e-6 * w),
            ]
        )
        j = np.array(
            [
                [1.0, 0.0, -2e-9 * w],
                [0.0, 1e3, -1e-3],
                [0.0, 1e-2, eps - 1e-8],
            ]
        )
        return r, j

    z0 = np.array([1.0e6, 3.0, 0.0])
    j_hi = rj(z0)[1]
    assert detect_degenerate(j_hi, j_hi, 4.0) == 0  # not sigma-proportional: no detection
    res = family_root_1d(rj, z0, tol=1e-9, max_lengths=5.0)
    assert res.status == "root", res.status
    assert res.z_star is not None
    assert abs(res.z_star[2] - w_star) < 1e-3 * abs(w_star)
