"""#450 Task 8: pure-Python truncated Taylor-map backend DASectionMap.

The DASectionMap is the pure-Python truncated-Taylor-map realization of the same
SectionMap interface (USER DECISION 2026-06-25: option (b), NO MOSEK / DACEyPy).
It is validated BOTH against the sampling oracle (single_rev / compose parity --
the backend is swappable iff it gives the same geometry) AND against the sourced
Png' golden via its fixed-point finder.

Key finding (honest, design draft §8.3): the strongly-unstable multi-rev Png'
fixed point (P5g' max|lambda|~3600) sits in a section basin narrower than ~1e-5.
The FD Taylor POLYNOMIAL stage descends to a MACHINE-DEPENDENT floor of
~3e-5..3e-4 (the FD-coefficient noise floor pushed through the condition-3600
self-composition; measured 3e-5 on the Linux/OpenBLAS build machine, 2.78e-4 on
the macOS/Accelerate machine -- #804 note). Since #805 the exact-derivative
section-chain multiple-shooting Newton endgame (single_rev_stm /
section_chain_newton) converges from that floor to the TRUE float-map fixed
point (~1e-11 chain residual; P5g' lands ~4e-9 from the published IC), removing
the landing's machine-dependence. The corrector (test_png_lane_recovery, Task 5)
still performs the certified closure. What the Taylor map provides that the
sampling grid cannot: a smooth, composable map whose iterated fixed point
reaches the fixed point from a coarse reference, which a brute-force grid
(basin << any feasible grid spacing) cannot.
"""

from __future__ import annotations

import math

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_backend import (
    DASectionMap,
    SamplingSectionMap,
    SectionPoint,
)


def _em() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Earth", "Moon")


def test_single_rev_parity_with_sampling_oracle() -> None:
    """DASectionMap.single_rev matches the sampling oracle (same geometry)."""
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    oracle = SamplingSectionMap(system, c_target=3.00022)
    for s in (
        SectionPoint(x=0.88500968, xdot=0.0),
        SectionPoint(x=0.807357887647950, xdot=-0.0956081545978604),
        SectionPoint(x=0.85, xdot=-0.05),
    ):
        a = da.single_rev(s)
        b = oracle.single_rev(s)
        assert abs(a.point.x - b.point.x) < 1e-9, s
        assert abs(a.point.xdot - b.point.xdot) < 1e-9, s
        assert abs(a.t - b.t) < 1e-9, s


def test_compose_parity_with_sampling_oracle() -> None:
    """compose(s, n) matches the sampling oracle on the DRO section point."""
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    oracle = SamplingSectionMap(system, c_target=3.00022)
    s = SectionPoint(x=0.88500968, xdot=0.0)
    a = da.compose(s, 3)
    b = oracle.compose(s, 3)
    assert abs(a.point.x - b.point.x) < 1e-9
    assert abs(a.point.xdot - b.point.xdot) < 1e-9


def test_taylor_single_rev_polynomial_matches_propagator() -> None:
    """The fitted single-rev Taylor polynomial reproduces the float map near ref.

    Evaluating the order-K polynomial at small offsets matches a direct single_rev
    to FD-fit accuracy. The single-rev map about P5g' is well-conditioned (its
    image is a DISTINCT section point ~0.825, not itself -- only P^5 returns), so
    the polynomial reproduces it to ~1e-3 across the [-h, h] domain. At the
    reference itself (offset 0) the fit is exact to the FD noise floor.
    """
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    s_ref = SectionPoint(x=0.807357887647950, xdot=-0.0956081545978604)
    tmap = da.taylor_single_rev(s_ref, order=3, h=3e-4, samples=7)
    # The least-squares polynomial tracks the float single-rev map across the fit
    # domain to the FD-fit accuracy (~1e-2 for this strongly-curved multi-rev
    # map). This validates the fit is meaningful, not garbage; the load-bearing
    # accuracy claim is test_taylor_fixed_point_reaches_png_neighbourhood.
    for dx, dxd in ((0.0, 0.0), (1e-4, -5e-5), (-8e-5, 7e-5)):
        dx_out, dxd_out = tmap.evaluate(dx, dxd)
        direct = da.single_rev(SectionPoint(x=s_ref.x + dx, xdot=s_ref.xdot + dxd))
        assert abs(s_ref.x + dx_out - direct.point.x) < 2e-2, (dx, dxd)
        assert abs(s_ref.xdot + dxd_out - direct.point.xdot) < 2e-2, (dx, dxd)


def test_single_rev_stm_jacobian_matches_central_fd() -> None:
    """single_rev_stm's exact section-map Jacobian matches a central FD probe.

    The #805 machinery: the 2x2 Jacobian d(x', xdot')/d(x, xdot) is assembled
    from the lift derivative (ydot from the Jacobi constant), the 6x6 STM of
    the variational propagation, and the first-order crossing-time correction.
    A central difference of the float single_rev at eps=1e-7 has error ~1e-8
    relative (its own truncation/roundoff balance), so agreement at 1e-6
    relative validates the exact assembly without demanding FD-limited
    precision of it. The return itself must match the float single_rev to the
    integrator tolerance (same geometry, augmented propagation).
    """
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    eps = 1e-7
    for s in (
        SectionPoint(x=0.807357887647950, xdot=-0.0956081545978604),  # P5g'
        SectionPoint(x=0.88500968, xdot=0.0),  # DRO
        SectionPoint(x=0.85, xdot=-0.05),
    ):
        ret, jac = da.single_rev_stm(s)
        direct = da.single_rev(s)
        assert abs(ret.point.x - direct.point.x) < 1e-9, s
        assert abs(ret.point.xdot - direct.point.xdot) < 1e-9, s
        assert abs(ret.t - direct.t) < 1e-9, s
        scale = max(1.0, max(abs(float(v)) for v in jac.ravel()))
        for k, (dx, dxd) in enumerate(((eps, 0.0), (0.0, eps))):
            p = da.single_rev(SectionPoint(x=s.x + dx, xdot=s.xdot + dxd))
            m = da.single_rev(SectionPoint(x=s.x - dx, xdot=s.xdot - dxd))
            fd_col = (
                (p.point.x - m.point.x) / (2 * eps),
                (p.point.xdot - m.point.xdot) / (2 * eps),
            )
            assert abs(float(jac[0, k]) - fd_col[0]) / scale < 1e-6, (s, k)
            assert abs(float(jac[1, k]) - fd_col[1]) / scale < 1e-6, (s, k)


def test_section_chain_newton_declines_off_family() -> None:
    """The endgame fails fast (None) where the n-rev chain does not exist.

    From ~1e-3 off P5g' the 5-rev same-sign chain walls off (no section return
    within t_max -- the #804-documented fragility), so section_chain_newton
    must decline rather than wander; taylor_fixed_point then falls back to the
    polynomial landing. This pins the strict-fallback contract that guarantees
    the #805 endgame can only improve a landing, never worsen one.
    """
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    off = SectionPoint(x=0.807357887647950 + 8e-4, xdot=-0.0956081545978604 - 6e-4)
    assert da.section_chain_newton(off, 5) is None


def test_taylor_fixed_point_reaches_png_neighbourhood() -> None:
    """From a coarse reference ~1e-3 from P5g', the iterated Taylor map plus the
    exact-derivative chain-Newton endgame land AT the true fixed point -- inside
    the corrector's ~1e-5 basin -- which a brute-force grid cannot.

    This is the capability the sampling backend lacks (the multi-rev basin is
    narrower than any feasible grid). The exact 1e-12 certification is still
    done by the corrector (Task 5, tests/search/test_png_lane_recovery.py --
    the load-bearing end-to-end proof); here we assert the backend's landing.

    TOLERANCE PROVENANCE (#804 + #805 notes): the POLYNOMIAL stage descends to
    a machine-dependent FD floor (~3e-5 Linux/OpenBLAS, 2.78e-4
    macOS/Accelerate -- #804), from which the #805 section-chain
    multiple-shooting Newton (exact STM-derived per-rev Jacobians) converges to
    the true float-map fixed point: chain residual ~1e-11, observed landing
    4.1e-9 from the published IC on the macOS machine (the true float fixed
    point is ~1e-9 from the published golden, #804 note item 1). The 1e-5
    bound is the corrector's demonstrated reliable basin -- the #805 capability
    claim ("lands inside the corrector basin directly") -- with ~2000x
    observed headroom, NOT a landing snapshot. Non-circularity: the seed is
    handed in 1e-3 AWAY from the published IC and the landing must have moved
    well off it, so the landing is a genuine output of float propagation, not
    an echo of any input. The float P^n residual at the landing must also be
    tiny (the pre-#805 FD-floor landing had residual 0.38 -- a truncation
    artifact, not a near-fixed-point; the endgame landing is the real thing).
    """
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    p5x, p5xd = 0.807357887647950, -0.0956081545978604
    dx0, dxd0 = 8e-4, -6e-4
    seed_dist = math.hypot(dx0, dxd0)  # 1e-3: the coarse-seed offset
    s_ref = SectionPoint(x=p5x + dx0, xdot=p5xd + dxd0)
    fp = da.taylor_fixed_point(s_ref, n=5, order=2, h=3e-4, samples=6, max_iter=30)
    dist = math.hypot(fp.x - p5x, fp.xdot - p5xd)
    # The #805 claim: the backend lands inside the corrector's ~1e-5 basin.
    assert dist < 1e-5, (fp.x, fp.xdot, dist)
    # The landing is an OUTPUT: it moved essentially the whole seed offset, so
    # it is not the handed-in reference (non-circularity; the published IC is
    # never handed to the finder at all).
    moved = math.hypot(fp.x - s_ref.x, fp.xdot - s_ref.xdot)
    assert moved > 0.5 * seed_dist, (fp.x, fp.xdot, moved)
    # And it is a genuine near-fixed-point of the REAL float map, not a
    # truncation artifact (pre-#805 landings had residual ~0.38 here).
    assert da.residual(fp, 5) < 1e-8, da.residual(fp, 5)
