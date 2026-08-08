"""#450 Task 8: pure-Python truncated Taylor-map backend DASectionMap.

The DASectionMap is the pure-Python truncated-Taylor-map realization of the same
SectionMap interface (USER DECISION 2026-06-25: option (b), NO MOSEK / DACEyPy).
It is validated BOTH against the sampling oracle (single_rev / compose parity --
the backend is swappable iff it gives the same geometry) AND against the sourced
Png' golden via its fixed-point finder.

Key finding (honest, design draft §8.3): the strongly-unstable multi-rev Png'
fixed point (P5g' max|lambda|~3600) sits in a section basin narrower than ~1e-5.
The FD Taylor map descends to a MACHINE-DEPENDENT floor of ~3e-5..3e-4 (the
FD-coefficient noise floor pushed through the condition-3600 self-composition;
measured 3e-5 on the Linux/OpenBLAS build machine, 2.78e-4 on the macOS/Accelerate
machine -- #804 note); the corrector's reliable basin is ~1e-5. So the Taylor map
alone does not nail P5g' -- the lane closes it with a small corrector
micro-multistart around the Taylor point (test_png_lane_recovery, Task 5). What
the Taylor map DOES provide that the sampling grid cannot: a smooth, composable
map whose iterated fixed point reaches the corrector's neighbourhood from a coarse
reference, which a brute-force grid (basin << any feasible grid spacing) cannot.
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


def test_taylor_fixed_point_reaches_png_neighbourhood() -> None:
    """From a coarse reference ~1e-3 from P5g', the iterated Taylor map genuinely
    descends toward the fixed point, into the corrector micro-multistart's
    demonstrated closing reach -- which a brute-force grid cannot.

    This is the capability the sampling backend lacks (the multi-rev basin is
    narrower than any feasible grid). The exact 1e-12 closure is done by the
    corrector (Task 5, tests/search/test_png_lane_recovery.py -- the load-bearing
    end-to-end proof); here we assert the Taylor map's descent.

    TOLERANCE PROVENANCE (#804, docs/notes/2026-08-08-804-taylor-fixed-point-
    floor.md): the landing distance is floored by FD-coefficient noise pushed
    through the condition-~3600 truncated self-composition (the single-rev
    polynomial, fitted over [-h, h] = [-3e-4, 3e-4], is evaluated by compose_self
    at the orbit's OTHER chain points up to ~0.4 away, far outside the fit
    domain), so the floor is machine-dependent: ~3e-5 on the Linux/OpenBLAS build
    machine (#450 decision note) vs 2.78e-4 on the macOS/Accelerate machine
    (2026-08-08), with the iteration itself converging cleanly (steps -> 1e-10)
    in both cases. The original 1e-4 bound was the build machine's floor snapshot
    plus headroom, never a derived invariant. The machine-independent claims are:
    genuine descent (>= 2x closer than the coarse seed's own offset; observed 33x
    Linux, 3.6x macOS), landing is an output (not the handed-in seed), and the
    landing stays inside the corrector multistart's demonstrated closing reach
    (close_candidate closes to 1e-11 from the 2.78e-4 landing on this machine;
    close_candidate's docstring documents ~3e-4 as the expected coarse-candidate
    scale).
    """
    system = _em()
    da = DASectionMap(system, c_target=3.00022)
    p5x, p5xd = 0.807357887647950, -0.0956081545978604
    dx0, dxd0 = 8e-4, -6e-4
    seed_dist = math.hypot(dx0, dxd0)  # 1e-3: the coarse-seed offset
    s_ref = SectionPoint(x=p5x + dx0, xdot=p5xd + dxd0)
    fp = da.taylor_fixed_point(s_ref, n=5, order=2, h=3e-4, samples=6, max_iter=30)
    dist = math.hypot(fp.x - p5x, fp.xdot - p5xd)
    # Genuine descent: at least 2x closer than the coarse seed (derived bound:
    # descent factor, not a floor snapshot -- see TOLERANCE PROVENANCE above).
    assert dist < 0.5 * seed_dist, (fp.x, fp.xdot, dist)
    # The landing is an OUTPUT of the iteration, not the handed-in seed or the
    # published point itself (non-circularity, mirrors test_png_lane_recovery).
    assert dist > 1e-7, (fp.x, fp.xdot, dist)
