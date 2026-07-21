"""#676 -- tests for the INTER-RETURN QR REFRAMING (discrete analogue of #669's QR).

Stage 10 of the Wilczak-Zgliczynski proof-machinery build (``#668``-``#675``).  #675's
composed covering COLLAPSED (ratio 0.284 -> 0.027 at N=2) because the composition carried
the reachable set from return to return as an axis-aligned interval BOX: the leg-(k+1) IC
box is built component-wise (``section_map_meanvalue_image``), which is the axis-aligned
hull of a thin, sheared true image, and ``rigorous_section_map`` then re-inflates from that
hull -- the SAME wrapping effect #669's within-flow QR reframing was built to defeat,
reappearing at the DISCRETE return-to-return boundary.

Before trusting the reframing on the real regularized CR3BP map, the MECHANISM is validated
here on maps whose composed image is known in closed form:

  * ``qr_rebasing_factor`` / ``reframe_image_qr`` / ``reframe_jacobian_qr`` SOUNDNESS -- the
    reframed enclosure rigorously contains the true image / composed Jacobian;
  * the LOAD-BEARING mechanism control -- for a linear stretch+rotation map (the textbook
    wrapping example, #669's own positive control was a rigid rotation of a wide box) whose
    N-fold image is an EXACTLY-known parallelepiped, an un-reframed box-hull composition
    wraps measurably over 3 returns while the reframed composition tracks the true extent;
  * the covering-relevant mechanism -- for a nonlinear map with a genuine box->Jacobian
    coupling (DP widens with the IC box, exactly the real flyby map's behaviour), image
    reframing keeps the leg IC box tighter -> tighter DP_k -> the composed covering ratio at
    N=3 is BETTER than the un-reframed box-hull composition's;
  * the DESIGN DECISION -- reframing on the single-leg image (``reframe_image_qr``) is what
    robustly helps; separately QR-reframing the ACCUMULATED composed Jacobian
    (``reframe_jacobian_qr``) adds its own rebasing wrapping that, at the low return counts
    reachable here, does NOT help on top of image reframing -- so the driver reframes the
    image and keeps the composed Jacobian as the raw chain-rule product.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#676)",
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


def _rot_stretch(theta: float, lam: float) -> list[list[Any]]:
    """Point interval matrix  M = R(theta) @ diag(lam, 1/lam)  (rotate a hyperbolic stretch)."""
    r = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    m = r @ np.array([[lam, 0.0], [0.0, 1.0 / lam]])
    return [[_pt(float(m[i][j])) for j in range(2)] for i in range(2)]


# --------------------------------------------------------------------------- #
# 1. Soundness of the reframing primitives.                                    #
# --------------------------------------------------------------------------- #
def test_qr_rebasing_factor_orthogonal_and_encloses() -> None:
    """A is orthogonal (A^T A = I) and A @ K rigorously encloses M."""
    iv = mp.iv
    m = [[_pt(3.0), _pt(1.4)], [_pt(-2.1), _pt(0.8)]]
    reb = vti.qr_rebasing_factor(iv, m)
    assert reb is not None
    a, kfac = reb
    # A orthogonal: A^T A = I (A is the numpy float QR factor, so to float-QR precision)
    ata = vti._imatmul(iv, [[a[j][i] for j in range(2)] for i in range(2)], a)
    for i in range(2):
        for j in range(2):
            target = 1.0 if i == j else 0.0
            assert abs(float(ata[i][j].mid) - target) < 1e-12
    # A @ K rigorously encloses M
    ak = vti._imatmul(iv, a, kfac)
    for i in range(2):
        for j in range(2):
            assert bool(ak[i][j].a <= m[i][j].a) and bool(ak[i][j].b >= m[i][j].b)


def test_reframe_image_qr_encloses_true_image() -> None:
    """The reframed parallelepiped hull contains the true image of any realisation in the box."""
    iv = mp.iv
    m = _rot_stretch(0.5, 1.6)
    frame = [[_pt(1.0), _pt(0.0)], [_pt(0.0), _pt(1.0)]]
    rbox = [iv.mpf([-1e-3, 1e-3]), iv.mpf([-1e-5, 1e-5])]
    center = [_pt(0.0), _pt(0.0)]
    res = vti.reframe_image_qr(iv, m, frame, rbox, center)
    assert res is not None
    c_k, a_k, r_k = res
    box = vti.enclosure_box(iv, c_k, a_k, r_k)
    # concrete realisation: the box corner (+ru, +rs) mapped exactly by the point M
    corner = [_pt(1e-3), _pt(1e-5)]
    true_img = vti._imatvec(iv, m, corner)
    for i in range(2):
        assert bool(box[i].a <= true_img[i].a) and bool(box[i].b >= true_img[i].b)


def test_reframe_jacobian_qr_encloses_composed() -> None:
    """The reframed accumulated Jacobian A_new @ (K @ comp_rest) encloses DP_2 @ DP_1."""
    iv = mp.iv
    dp1 = [[_pt(2.0), _pt(0.3)], [_pt(-0.4), _pt(0.6)]]
    dp2 = [[_pt(1.7), _pt(-0.5)], [_pt(0.2), _pt(0.9)]]
    reb = vti.qr_rebasing_factor(iv, dp1)
    assert reb is not None
    comp_frame, comp_rest = reb
    res = vti.reframe_jacobian_qr(iv, dp2, comp_frame, comp_rest)
    assert res is not None
    a_new, rest_new = res
    reframed = vti._imatmul(iv, a_new, rest_new)
    raw = vti.compose_section_jacobians(iv, [dp1, dp2])  # DP_2 @ DP_1
    for i in range(2):
        for j in range(2):
            assert bool(reframed[i][j].a <= raw[i][j].a) and bool(reframed[i][j].b >= raw[i][j].b)


# --------------------------------------------------------------------------- #
# 2. Load-bearing mechanism: linear wrapping with an EXACTLY known composed image. #
# --------------------------------------------------------------------------- #
def _true_hull(iv: Any, m: list[list[Any]], box0: list[Any], k: int) -> list[Any]:
    """Tightest box hull of the true image M^k @ box0 (single matmul of the exact M^k)."""
    m_np = np.array([[float(m[i][j].mid) for j in range(2)] for i in range(2)])
    mk = np.linalg.matrix_power(m_np, k)
    mk_iv = [[_pt(float(mk[i][j])) for j in range(2)] for i in range(2)]
    return vti._imatvec(iv, mk_iv, box0)


def _box_prop(iv: Any, m: list[list[Any]], box0: list[Any], k: int) -> list[Any]:
    """Un-reframed box propagation: box-hull of M @ (...) at every step (wraps)."""
    box = list(box0)
    for _ in range(k):
        box = vti._imatvec(iv, m, box)
    return box


def _reframed(iv: Any, m: list[list[Any]], box0: list[Any], k: int) -> list[Any]:
    """Reframed propagation via reframe_image_qr (carry a point frame + tight box)."""
    center = [_pt(0.0), _pt(0.0)]
    frame = [[_pt(1.0), _pt(0.0)], [_pt(0.0), _pt(1.0)]]
    rbox = list(box0)
    for _ in range(k):
        center = vti._imatvec(iv, m, center)
        res = vti.reframe_image_qr(iv, m, frame, rbox, center)
        assert res is not None
        center, frame, rbox = res
    return vti.enclosure_box(iv, center, frame, rbox)


def test_reframing_tracks_true_extent_while_boxprop_wraps() -> None:
    """Linear stretch+rotation: reframed hull tracks the true extent through 3 returns while
    the un-reframed box-hull composition wraps measurably -- the mechanism, in closed form."""
    iv = mp.iv
    m = _rot_stretch(0.5, 1.6)
    box0 = [iv.mpf([-1e-3, 1e-3]), iv.mpf([-1e-5, 1e-5])]
    for k in (1, 2, 3):
        true_w = max(_w(c) for c in _true_hull(iv, m, box0, k))
        box_w = max(_w(c) for c in _box_prop(iv, m, box0, k))
        rf_w = max(_w(c) for c in _reframed(iv, m, box0, k))
        # reframed always soundly encloses the true extent
        assert rf_w >= true_w * (1.0 - 1e-9)
        # reframed stays essentially AT the true extent; box-prop does not
        assert rf_w <= 1.02 * true_w, f"reframed should track true at k={k}"
        if k == 1:
            assert abs(box_w - true_w) < 1e-18  # one step: no wrapping yet, identical
        if k == 3:
            assert box_w > 1.3 * true_w, "box-prop should wrap measurably by k=3"
            assert rf_w < box_w, "reframing should beat box-prop by k=3"


# --------------------------------------------------------------------------- #
# 3. Covering-relevant mechanism: nonlinear box->Jacobian coupling (the real behaviour). #
# --------------------------------------------------------------------------- #
def _covering_chain(
    iv: Any,
    theta: float,
    lam: float,
    c: float,
    ru: float,
    rs: float,
    n_returns: int,
    *,
    reframe: bool,
) -> list[float]:
    """Composed-covering separation/width ratio per return for the nonlinear map
    ``P(w) = M w + c*[w0 w1, 0]`` (fixed point at 0, DP widens with the IC box).

    ``reframe=False`` reproduces #675's box-hull box chain + raw composed Jacobian;
    ``reframe=True`` uses ``reframe_image_qr`` for the box chain (composed Jacobian stays the
    raw chain-rule product -- the driver's design choice).  Returns the ratio at each return.
    """
    m = _rot_stretch(theta, lam)

    def p_center(w: list[Any]) -> list[Any]:
        mv = vti._imatvec(iv, m, w)
        return [mv[0] + _pt(c) * w[0] * w[1], mv[1]]

    def dp_over(box: list[Any]) -> list[list[Any]]:
        return [[m[0][0] + _pt(c) * box[1], m[0][1] + _pt(c) * box[0]], [m[1][0], m[1][1]]]

    def foff(name: str) -> list[Any]:
        a = {"box": iv.mpf([-1, 1]), "left": iv.mpf([-1, -1]), "right": iv.mpf([1, 1])}[name]
        return [a * _pt(ru), iv.mpf([-1, 1]) * _pt(rs)]

    c0 = [_pt(0.0), _pt(0.0)]
    box_ic = [c0[0] + foff("box")[0], c0[1] + foff("box")[1]]
    center = list(c0)
    frame = [[_pt(1.0), _pt(0.0)], [_pt(0.0), _pt(1.0)]]
    rbox = list(foff("box"))
    box_center = list(c0)
    dp_boxes: list[Any] = []
    ratios: list[float] = []
    for _ in range(n_returns):
        dp_k = dp_over(box_ic)
        dp_boxes.append(dp_k)
        wk_hat = p_center(center)
        center = wk_hat
        if reframe:
            res = vti.reframe_image_qr(iv, dp_k, frame, rbox, wk_hat)
            assert res is not None
            _, frame, rbox = res
            box_ic = vti.enclosure_box(iv, wk_hat, frame, rbox)
        else:
            off = [box_ic[i] - box_center[i] for i in range(2)]
            box_ic = vti.section_map_meanvalue_image(iv, wk_hat, dp_k, off)
            box_center = wk_hat
        jac = vti.compose_section_jacobians(iv, dp_boxes)
        left = vti.section_map_meanvalue_image(iv, wk_hat, jac, foff("left"))
        right = vti.section_map_meanvalue_image(iv, wk_hat, jac, foff("right"))
        sep = abs(float(left[0].mid) - float(right[0].mid))
        edgew = max(_w(left[0]), _w(right[0]))
        ratios.append(sep / edgew if edgew > 0 else float("inf"))
    return ratios


def test_image_reframing_improves_composed_covering_ratio_at_three_returns() -> None:
    """Nonlinear box->Jacobian coupling: image reframing keeps the leg IC box tighter, so the
    composed covering ratio at N=3 exceeds the un-reframed box-hull chain's -- the covering-
    relevant statement of the mechanism (identical at N=1, benefit appears once returns
    compound)."""
    iv = mp.iv
    theta, lam, c, ru, rs = 1.1, 1.5, 20.0, 1e-4, 1e-6
    box_ratios = _covering_chain(iv, theta, lam, c, ru, rs, 3, reframe=False)
    rf_ratios = _covering_chain(iv, theta, lam, c, ru, rs, 3, reframe=True)
    assert abs(box_ratios[0] - rf_ratios[0]) < 1e-9 * box_ratios[0]  # N=1 identical
    assert rf_ratios[2] > box_ratios[2], "image reframing should improve the N=3 covering ratio"
    assert rf_ratios[2] > 1.2 * box_ratios[2], "the N=3 improvement should be substantial"


def test_jacobian_reframing_does_not_help_on_top_of_image_reframing() -> None:
    """DESIGN DECISION control: additionally QR-reframing the ACCUMULATED composed Jacobian
    (reframe_jacobian_qr) adds its own rebasing wrapping that, at the low return counts
    reachable here, does not beat the raw chain-rule product on top of image reframing -- the
    reason the driver reframes the IMAGE (single-leg DP_k) but keeps the composed Jacobian raw."""
    iv = mp.iv
    theta, lam, c, ru, rs = 1.1, 1.5, 20.0, 1e-4, 1e-6
    m = _rot_stretch(theta, lam)

    def p_center(w: list[Any]) -> list[Any]:
        mv = vti._imatvec(iv, m, w)
        return [mv[0] + _pt(c) * w[0] * w[1], mv[1]]

    def dp_over(box: list[Any]) -> list[list[Any]]:
        return [[m[0][0] + _pt(c) * box[1], m[0][1] + _pt(c) * box[0]], [m[1][0], m[1][1]]]

    def foff(name: str) -> list[Any]:
        a = {"box": iv.mpf([-1, 1]), "left": iv.mpf([-1, -1]), "right": iv.mpf([1, 1])}[name]
        return [a * _pt(ru), iv.mpf([-1, 1]) * _pt(rs)]

    def ratio_at_n3(reframe_jac: bool) -> float:
        c0 = [_pt(0.0), _pt(0.0)]
        box_ic = [foff("box")[0], foff("box")[1]]
        center = list(c0)
        frame = [[_pt(1.0), _pt(0.0)], [_pt(0.0), _pt(1.0)]]
        rbox = list(foff("box"))
        dp_boxes: list[Any] = []
        comp_frame: list[list[Any]] | None = None
        comp_rest: list[list[Any]] | None = None
        r = 0.0
        for _ in range(3):
            dp_k = dp_over(box_ic)
            dp_boxes.append(dp_k)
            wk_hat = p_center(center)
            center = wk_hat
            res = vti.reframe_image_qr(iv, dp_k, frame, rbox, wk_hat)
            assert res is not None
            _, frame, rbox = res
            box_ic = vti.enclosure_box(iv, wk_hat, frame, rbox)
            if reframe_jac:
                if comp_frame is None or comp_rest is None:
                    reb = vti.qr_rebasing_factor(iv, dp_k)
                else:
                    reb = vti.reframe_jacobian_qr(iv, dp_k, comp_frame, comp_rest)
                assert reb is not None
                comp_frame, comp_rest = reb
                jac = vti._imatmul(iv, comp_frame, comp_rest)
            else:
                jac = vti.compose_section_jacobians(iv, dp_boxes)
            left = vti.section_map_meanvalue_image(iv, wk_hat, jac, foff("left"))
            right = vti.section_map_meanvalue_image(iv, wk_hat, jac, foff("right"))
            sep = abs(float(left[0].mid) - float(right[0].mid))
            edgew = max(_w(left[0]), _w(right[0]))
            r = sep / edgew if edgew > 0 else float("inf")
        return r

    raw_jac = ratio_at_n3(reframe_jac=False)
    reframed_jac = ratio_at_n3(reframe_jac=True)
    # image reframing + RAW composed Jacobian is at least as good as also reframing the Jacobian
    assert raw_jac >= reframed_jac, "reframing the accumulated Jacobian should not help at low N"
