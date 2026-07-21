"""#677 -- tests for FINER SUB-BOXING (subdivide-before-propagate) of a composed map.

Stage 11 of the Wilczak-Zgliczynski proof-machinery build (``#668``-``#676``).  #676
pinned the precise cause of #675's composed-covering collapse: it is NOT the inter-return
box-hull *decorrelation* that #676's QR reframing removes (reframing changed neither N=2 nor
N=3 on Oterma), it is intra-leg Jacobian *variation* over a leg IC box that has already
inflated to ~35x the original h-set by the 2nd return, over which the section-map Jacobian
[DP] genuinely varies ~3.7-4.0x.  The mean-value image applies that whole (loose) [DP] once
to the offset, so its width tracks the box's OWN inflation, not the true image extent.

The fix this diagnosis calls for (the standard CAPD-style W-Z technique): never let a single
propagated box get near 35x.  Subdivide the h-set into a modest N x N grid of small
sub-boxes BEFORE propagating, run each sub-box through its OWN mean-value chain (each
sub-box's own local [DP] stays near-constant across its own small extent, so its mean-value
image stays tight relative to ITS size), and take the UNION of the sub-box images as the
return's enclosure.

Before trusting the mechanism on the real regularized CR3BP map (``scripts/
certify_677_wz_oterma_subdivided_covering.py``), the MECHANISM is validated here on a
composed map whose true image is known by dense sampling and whose Jacobian genuinely
varies over the box (a strong box->Jacobian coupling, mirroring the real flyby map):

  * ``union_interval_boxes`` SOUNDNESS -- the union hull encloses every sub-box, and a
    concrete disjoint example returns the exact componentwise hull;
  * the LOAD-BEARING mechanism control -- for the nonlinear map ``P(w) = M w + c[w0 w1, 0]``
    with a strong stretch (lambda=3) and coupling (c=300), naive SINGLE-box mean-value
    propagation inflates to ~3x (N=2) / ~5x (N=3) the true image extent (mirroring the
    Oterma collapse), while subdividing into N x 1 sub-boxes and unioning tracks the true
    extent to <=1.5x -- checked NUMERICALLY against a dense-sampled ground truth, both for
    tightness and for SOUNDNESS (the union encloses every sampled true-image point);
  * the sub-box-count TRADEOFF -- finer subdivision monotonically tightens the union
    (nu=2 -> nu=3 -> nu=5 each strictly better at N=3), the real feasibility-vs-cost knob.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#677)",
)

import numpy as np  # noqa: E402

import scripts._validated_taylor_integrator as vti  # noqa: E402


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _w(x: Any) -> float:
    return float(x.delta.b)


# --------------------------------------------------------------------------- #
# 1. Soundness of the union primitive.                                         #
# --------------------------------------------------------------------------- #
def test_union_interval_boxes_encloses_and_is_hull() -> None:
    """The union hull encloses every input box and equals the exact componentwise hull."""
    iv = mp.iv
    b1 = [iv.mpf([-1.0, 0.5]), iv.mpf([2.0, 3.0])]
    b2 = [iv.mpf([0.0, 2.0]), iv.mpf([-1.0, 2.5])]
    u = vti.union_interval_boxes(iv, [b1, b2])
    # encloses every input box
    for b in (b1, b2):
        for i in range(2):
            assert bool(u[i].a <= b[i].a) and bool(u[i].b >= b[i].b)
    # exact componentwise hull
    assert abs(float(u[0].a) - (-1.0)) < 1e-30 and abs(float(u[0].b) - 2.0) < 1e-30
    assert abs(float(u[1].a) - (-1.0)) < 1e-30 and abs(float(u[1].b) - 3.0) < 1e-30
    # a single box is returned unchanged
    single = vti.union_interval_boxes(iv, [b1])
    for i in range(2):
        assert bool(single[i].a <= b1[i].a) and bool(single[i].b >= b1[i].b)


# --------------------------------------------------------------------------- #
# The nonlinear stretching map with a genuine box->Jacobian coupling.          #
#   P(w) = M w + c*[w0 w1, 0],  M = R(theta) diag(lambda, 1/lambda).           #
# DP(w) = [[m00 + c w1, m01 + c w0],[m10, m11]] varies strongly across a box.   #
# --------------------------------------------------------------------------- #
_THETA, _LAM, _C = 0.7, 3.0, 300.0


def _m_pt() -> list[list[Any]]:
    r = np.array([[math.cos(_THETA), -math.sin(_THETA)], [math.sin(_THETA), math.cos(_THETA)]])
    m = r @ np.array([[_LAM, 0.0], [0.0, 1.0 / _LAM]])
    return [[_pt(float(m[i][j])) for j in range(2)] for i in range(2)]


def _p_point(mc: np.ndarray, w: list[float]) -> list[float]:
    return [mc[0, 0] * w[0] + mc[0, 1] * w[1] + _C * w[0] * w[1], mc[1, 0] * w[0] + mc[1, 1] * w[1]]


def _dp_over(m: list[list[Any]], box: list[Any]) -> list[list[Any]]:
    return [[m[0][0] + _pt(_C) * box[1], m[0][1] + _pt(_C) * box[0]], [m[1][0], m[1][1]]]


def _subbox_chain(
    iv: Any,
    m: list[list[Any]],
    mc: np.ndarray,
    a_lo: float,
    a_hi: float,
    ru: float,
    rs: float,
    n: int,
) -> list[Any]:
    """Independent single-box mean-value chain of the sub-box [a_lo,a_hi]*ru x [-1,1]*rs
    through ``n`` returns; returns its mean-value image box at return ``n`` (the #675 chain,
    seeded from a sub-box instead of the whole h-set)."""
    ac = 0.5 * (a_lo + a_hi)
    cenf = [ac * ru, 0.0]
    box = [iv.mpf([a_lo * ru, a_hi * ru]), iv.mpf([-rs, rs])]
    box_cen = list(cenf)
    dp_boxes: list[Any] = []
    for _ in range(n):
        dp = _dp_over(m, box)
        dp_boxes.append(dp)
        new_cenf = _p_point(mc, cenf)
        off_box = [box[0] - _pt(box_cen[0]), box[1] - _pt(box_cen[1])]
        box = vti.section_map_meanvalue_image(iv, [_pt(new_cenf[0]), _pt(new_cenf[1])], dp, off_box)
        box_cen = list(new_cenf)
        cenf = new_cenf
    comp = vti.compose_section_jacobians(iv, dp_boxes)
    off = [iv.mpf([(a_lo - ac) * ru, (a_hi - ac) * ru]), iv.mpf([-rs, rs])]
    return vti.section_map_meanvalue_image(iv, [_pt(cenf[0]), _pt(cenf[1])], comp, off)


def _true_box(mc: np.ndarray, ru: float, rs: float, n: int) -> list[list[float]]:
    """Dense-sampled true image box of P^n over the h-set [-1,1]ru x [-1,1]rs."""
    pts = []
    for ia in np.linspace(-1, 1, 81):
        for ib in np.linspace(-1, 1, 5):
            w = [ia * ru, ib * rs]
            for _ in range(n):
                w = _p_point(mc, w)
            pts.append(w)
    arr = np.array(pts)
    return [[float(arr[:, i].min()), float(arr[:, i].max())] for i in range(2)]


def _union_nu(
    iv: Any, m: list[list[Any]], mc: np.ndarray, ru: float, rs: float, n: int, nu: int
) -> list[Any]:
    edges = np.linspace(-1.0, 1.0, nu + 1)
    cells = [_subbox_chain(iv, m, mc, edges[i], edges[i + 1], ru, rs, n) for i in range(nu)]
    return vti.union_interval_boxes(iv, cells)


def test_subdivision_recovers_tightness_vs_singlebox() -> None:
    """LOAD-BEARING: single-box mean-value propagation inflates strongly over 2-3 returns
    (mirroring #676's Oterma collapse) while subdividing + unioning tracks the dense-sampled
    true image extent -- and the union rigorously ENCLOSES the true image (soundness)."""
    iv = mp.iv
    m = _m_pt()
    mc = np.array([[float(m[i][j].mid) for j in range(2)] for i in range(2)])
    ru, rs = 5e-3, 1e-5
    for n, single_lo, union_hi in ((2, 2.5, 1.5), (3, 4.0, 1.5)):
        true_b = _true_box(mc, ru, rs, n)
        true_hw = max(0.5 * (true_b[i][1] - true_b[i][0]) for i in range(2))
        single = _subbox_chain(iv, m, mc, -1.0, 1.0, ru, rs, n)  # nu=1 == whole box
        single_hw = max(0.5 * _w(single[i]) for i in range(2))
        union5 = _union_nu(iv, m, mc, ru, rs, n, 5)
        union_hw = max(0.5 * _w(union5[i]) for i in range(2))
        # single-box inflates well past the true extent (the collapse mechanism)
        assert single_hw / true_hw > single_lo, f"single-box should inflate at N={n}"
        # subdivision (nu=5) tracks the true extent
        assert union_hw / true_hw < union_hi, f"nu=5 union should track true at N={n}"
        assert union_hw < single_hw, f"subdivision must beat single-box at N={n}"
        # SOUNDNESS: the union hull encloses the whole dense-sampled true image
        for i in range(2):
            assert bool(union5[i].a <= _pt(true_b[i][0]).a) and bool(
                union5[i].b >= _pt(true_b[i][1]).b
            )


def test_finer_subdivision_monotonically_tightens() -> None:
    """The feasibility-vs-cost knob: at N=3 the union tightens strictly as nu grows
    (nu=2 -> nu=3 -> nu=5), each remaining a sound enclosure -- more sub-boxes (more
    propagations) buy a tighter enclosure."""
    iv = mp.iv
    m = _m_pt()
    mc = np.array([[float(m[i][j].mid) for j in range(2)] for i in range(2)])
    ru, rs, n = 5e-3, 1e-5, 3
    hw = {}
    for nu in (2, 3, 5):
        u = _union_nu(iv, m, mc, ru, rs, n, nu)
        hw[nu] = max(0.5 * _w(u[i]) for i in range(2))
    assert hw[2] > hw[3] > hw[5], f"finer subdivision should tighten monotonically: {hw}"
    single = _subbox_chain(iv, m, mc, -1.0, 1.0, ru, rs, n)
    single_hw = max(0.5 * _w(single[i]) for i in range(2))
    assert hw[2] < single_hw, "even the coarsest subdivision beats the single box"
