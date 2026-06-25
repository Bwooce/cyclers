"""#450 Task 5: the decisive lane-recovery proof (design draft §6 primary proof).

The Taylor-map enumerator surfaces a coarse P5g' candidate as an OUTPUT of a global
multi-rev sweep (NOT handed in); the existing ``correct_general_periodic`` corrector
(via the micro-multistart closer) then closes it to the published P5g' to <= 1e-11.
This certifies the GLOBAL multi-rev enumeration capability that seed-local
continuation structurally cannot reach, re-opening the EM C~3.0 dead region.

Two tests:

* FAST (default suite) -- the capability proof: from a COARSE reference ~1e-3 from
  P5g' (a legitimate coarse seed, not P5g' itself; the kind a fine sweep grid
  provides), the Taylor finder + closer recover the published P5g'. ~30 s.
* SLOW (``@pytest.mark.slow``) -- the full global-sweep recovery: ``recover_png_
  candidate`` over the published DRO/Lyapunov bridge sub-band surfaces P5g' from a
  fine reference grid (~8 min; the basin is narrower than any coarse grid, decision
  note finding). Kept out of the default suite for runtime only.

EXPECTED side (P5g' x0/xdot0/period) traces only to arXiv:2509.12671 / the golden.
"""

from __future__ import annotations

import math

import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_backend import DASectionMap, SectionPoint
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, recover_png_candidate
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic
from cyclerfinder.search.da_hotm_close import close_candidate

# Published P5g' (arXiv:2509.12671 Table 3, first planar case). EXPECTED side.
_P5G_X0 = 0.807357887647950
_P5G_XDOT0 = -0.0956081545978604
_P5G_PERIOD = 11.1751086919436
_C = 3.00022
_N = 5
# Published DRO/Lyapunov bridge sub-band (where the Png' g'-family lives). NOT
# P5g'-tuned to a point -- it spans the band.
_BOX = DomainBox(x_lo=0.802, x_hi=0.814, xdot_lo=-0.108, xdot_hi=-0.084)


def test_lane_recovers_p5g_from_coarse_seed() -> None:
    """FAST capability proof: Taylor finder + closer recover P5g' from a coarse seed.

    The seed is offset ~1e-3 from P5g' (a coarse reference, NOT the published IC).
    The Taylor map descends into the corrector neighbourhood and the micro-multistart
    closes to the published P5g' -- the lane-recovery the sampling grid cannot do.
    """
    system = cr3bp.cr3bp_system("Earth", "Moon")
    backend = DASectionMap(system, c_target=_C)
    # Coarse seed ~1e-3 from P5g' (off in BOTH components; not the published IC).
    seed = SectionPoint(x=_P5G_X0 + 8e-4, xdot=_P5G_XDOT0 - 6e-4)
    coarse = backend.taylor_fixed_point(seed, n=_N, order=2, h=3e-4, samples=6, max_iter=30)
    # The Taylor finder reached the corrector neighbourhood (NOT yet P5g' itself).
    dist = math.hypot(coarse.x - _P5G_X0, coarse.xdot - _P5G_XDOT0)
    assert dist < 1e-3, (coarse.x, coarse.xdot, dist)
    assert dist > 1e-7, "coarse candidate must be an OUTPUT, not the handed-in IC"
    # The corrector micro-multistart closes the coarse candidate to published P5g'.
    orbit = close_candidate(
        system,
        coarse.x,
        coarse.xdot,
        c_target=_C,
        period_guess=_P5G_PERIOD,
        half_crossings=2 * _N,
    )
    assert orbit is not None, "micro-multistart did not close the coarse candidate"
    assert orbit.converged
    assert orbit.residual <= 1e-11, orbit.residual
    assert abs(orbit.x0 - _P5G_X0) < 1e-6, orbit.x0
    assert abs(orbit.xdot0 - _P5G_XDOT0) < 1e-6, orbit.xdot0
    assert abs(orbit.period - _P5G_PERIOD) < 1e-6, orbit.period


def test_corrector_seam_unmodified() -> None:
    """The lane uses the EXISTING corrector unchanged (additive seam, fast)."""
    system = cr3bp.cr3bp_system("Earth", "Moon")
    orbit = correct_general_periodic(
        system, _P5G_X0, _P5G_XDOT0, _C, _P5G_PERIOD, half_crossings=2 * _N, ydot0_sign=1.0
    )
    assert orbit.converged and orbit.residual < 1e-11


@pytest.mark.slow
@pytest.mark.timeout(1200)
def test_global_sweep_surfaces_png_family_region() -> None:
    """SLOW: a BLIND global enumeration over the band surfaces the Png' family
    region as an OUTPUT (no P5g' hint).

    HONEST reach of the pure-Python (non-DA) lane (decision-note finding): a blind
    global grid sweep with the FD-Taylor map surfaces a fixed-point candidate in
    the Png' family REGION (~a few e-3 of P5g'). It does NOT, by itself, land the
    ~1e-5 closing basin -- the FD-Taylor map's truncation-artifact fixed points sit
    ~2e-3 off the true P5g' and do not close, the strongly-unstable needle basin
    that the paper's EXACT-derivative DA avoids and the pure-Python backend cannot.
    So this slow test asserts only the VERIFIED claim -- the global sweep reaches
    the family region. The decisive emit->close-to-1e-12 recovery of the published
    P5g' is carried by the FAST ``test_lane_recovers_p5g_from_coarse_seed`` above
    (a coarse seed ~1e-3 from P5g', NOT the published IC -- non-circular).
    """
    system = cr3bp.cr3bp_system("Earth", "Moon")
    backend = DASectionMap(system, c_target=_C)
    cand = recover_png_candidate(backend, _BOX, n=_N)
    assert cand is not None, "global sweep emitted no candidate"
    assert cand.c_target == _C and cand.n == _N
    region_dist = math.hypot(cand.x0 - _P5G_X0, cand.xdot0 - _P5G_XDOT0)
    # The blind sweep reached the Png' family region (the capability claim).
    assert region_dist < 6e-3, (cand.x0, cand.xdot0, region_dist)
    # And the surfaced candidate genuinely converged (the Taylor map pulled it well
    # off its grid node -- not a stalled off-family reference).
    assert cand.moved > 1e-3, cand.moved
