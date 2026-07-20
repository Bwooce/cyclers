"""#674 -- tests for the correlation-preserving (mean-value) section-map image.

Stage 8 of the Wilczak-Zgliczynski proof-machinery build (``#668``-``#673``).  #673's
covering attempt hit a REPRESENTATION wall: propagating each face of an h-set
independently through the section map and box-hulling each image re-introduces the
wrapping effect, so the axis-aligned hull is ~9x too wide and masks the genuine
hyperbolic edge separation.  ``section_map_meanvalue_image`` replaces that with a
rigorous mean-value enclosure  ``P(face) subset P(w_hat) + [DP] (w0_face - w_hat)``
built from ONE tight linearization, which preserves the image's thin/sheared shape.

The load-bearing validation here:

  * exactness on a purely affine "section map" (mean value is exact to first order and
    the Jacobian is constant, so the enclosure is TIGHT and equals the true image), and
  * the mechanism itself: a strongly-expanding affine map on a thin h-set yields
    mean-value edge images whose hyperbolic SEPARATION exceeds their WIDTH, so the
    #673 covering checker certifies the covering -- whereas artificially inflating the
    Jacobian enclosure width (mimicking the box-hull wrapping #673 hit) makes the edges
    overlap and the SAME checker correctly reports non-covering.  So the helper's
    tightness is what flips the verdict, exactly as on the real regularized CR3BP map.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#674)",
)

import scripts._validated_taylor_integrator as vti  # noqa: E402


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _w(x: Any) -> float:
    return float(x.delta.b)


def test_meanvalue_is_exact_for_thin_inputs() -> None:
    """center + J@offset with thin (point) inputs reproduces the arithmetic exactly."""
    iv = mp.iv
    center = [_pt(2.0), _pt(-1.0)]
    jac = [[_pt(3.0), _pt(0.5)], [_pt(-2.0), _pt(4.0)]]
    offset = [_pt(0.5), _pt(-0.25)]  # dyadic -> binary-exact arithmetic
    img = vti.section_map_meanvalue_image(iv, center, jac, offset)
    # row0 = 2 + 3*0.5 + 0.5*(-0.25) = 3.375 ; row1 = -1 + (-2)*0.5 + 4*(-0.25) = -3.0
    assert bool(img[0].a <= _pt(3.375)) and bool(img[0].b >= _pt(3.375))
    assert bool(img[1].a <= _pt(-3.0)) and bool(img[1].b >= _pt(-3.0))
    assert _w(img[0]) < 1e-30 and _w(img[1]) < 1e-30  # tight (no over-approximation)


def test_meanvalue_is_inclusion_monotone_in_offset() -> None:
    """A wider offset box yields a superset image (soundness of the enclosure)."""
    iv = mp.iv
    center = [_pt(0.0), _pt(0.0)]
    jac = [[_pt(3.0), _pt(0.0)], [_pt(0.0), _pt(0.25)]]
    narrow = vti.section_map_meanvalue_image(iv, center, jac, [iv.mpf([-1e-3, 1e-3]), _pt(0.0)])
    wide = vti.section_map_meanvalue_image(iv, center, jac, [iv.mpf([-1e-2, 1e-2]), _pt(0.0)])
    assert bool(wide[0].a <= narrow[0].a) and bool(wide[0].b >= narrow[0].b)


# --------------------------------------------------------------------------- #
# Mechanism: tight mean-value image -> covering certifies; inflated -> fails.   #
# --------------------------------------------------------------------------- #
def _edge_images(iv: Any, jac: list[list[Any]], ru: float, rs: float) -> tuple[Any, Any, Any]:
    """Mean-value images (whole box, left edge a=-1, right edge a=+1) of a thin h-set.

    N is axis-aligned here (unstable axis = coord 0, stable axis = coord 1) with
    half-widths (ru, rs); center image at the origin.  offset of a face = a*ru*e0 +
    b*rs*e1 with a fixed on an edge, a,b in [-1,1] on the whole box.
    """
    center = [_pt(0.0), _pt(0.0)]

    def img(aval: Any) -> tuple[Any, Any]:
        off = [aval * iv.mpf(ru), iv.mpf([-1, 1]) * iv.mpf(rs)]
        out = vti.section_map_meanvalue_image(iv, center, jac, off)
        return out[0], out[1]

    box = img(iv.mpf([-1, 1]))
    left = img(iv.mpf([-1, -1]))
    right = img(iv.mpf([1, 1]))
    return left, right, box


def test_tight_meanvalue_image_certifies_covering() -> None:
    """Strong-stretch affine map, thin h-set: mean-value edges separate -> COVERS.

    jac stretches the unstable axis (coord 0) by 3 and contracts the stable axis by
    1/3.  With M placed at the image center (unstable half-width inside the edge gap),
    the two edges land strictly on opposite sides -> #673 checker certifies covering.
    """
    iv = mp.iv
    jac = [[_pt(3.0), _pt(0.0)], [_pt(0.0), _pt(1.0 / 3.0)]]
    ru, rs = 1e-6, 1e-8
    left, right, box = _edge_images(iv, jac, ru, rs)
    # unstable coord = coord 0; edge separation 2*ru*3 = 6e-6, per-edge width 2*rs*3 ~ tiny
    sep = abs(float(left[0].mid) - float(right[0].mid))
    edgew = max(_w(left[0]), _w(right[0]))
    assert sep > edgew  # disjoint edges (the property the box-hull destroyed in #673)
    hu = 0.4 * (sep - edgew)
    hs = 3.0 * (_w(box[1]) / 2.0 + rs)
    res = vti.covering_relation_2d(
        iv,
        m_center=[_pt(0.0), _pt(0.0)],
        m_uvec=[iv.mpf([hu, hu]), _pt(0.0)],
        m_svec=[_pt(0.0), iv.mpf([hs, hs])],
        image_left=[left[0], left[1]],
        image_right=[right[0], right[1]],
        image_box=[box[0], box[1]],
    )
    assert res["covers"] is True
    assert res["stable_containment"] is True and res["unstable_exit"] is True


def test_inflated_jacobian_width_breaks_covering() -> None:
    """Inflating the Jacobian enclosure width (mimicking box-hull wrapping) -> NO cover.

    Same map/h-set, but the unstable-row Jacobian entries carry a wide uncertainty
    (as the wrapped whole-box propagation did in #673): each edge image inflates past
    the genuine separation, the edges OVERLAP, and the SAME checker correctly reports
    non-covering.  This isolates the helper's tightness as the load-bearing ingredient.
    """
    iv = mp.iv
    infl = iv.mpf([-5.0, 5.0])  # huge width on the unstable row
    jac = [[iv.mpf(3.0) + infl, iv.mpf(0.0) + infl], [_pt(0.0), _pt(1.0 / 3.0)]]
    ru, rs = 1e-6, 1e-8
    left, right, box = _edge_images(iv, jac, ru, rs)
    sep = abs(float(left[0].mid) - float(right[0].mid))
    edgew = max(_w(left[0]), _w(right[0]))
    assert edgew > sep  # inflated edges overlap -> cannot show opposite-side exit
    hu = max(1e-9, 0.4 * abs(sep - edgew))
    hs = 3.0 * (_w(box[1]) / 2.0 + rs)
    res = vti.covering_relation_2d(
        iv,
        m_center=[_pt(0.0), _pt(0.0)],
        m_uvec=[iv.mpf([hu, hu]), _pt(0.0)],
        m_svec=[_pt(0.0), iv.mpf([hs, hs])],
        image_left=[left[0], left[1]],
        image_right=[right[0], right[1]],
        image_box=[box[0], box[1]],
    )
    assert res["covers"] is False
    assert res["unstable_exit"] is False
