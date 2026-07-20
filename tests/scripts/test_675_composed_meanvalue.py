"""#675 -- tests for the COMPOSED (multi-return) mean-value section-map image.

Stage 9 of the Wilczak-Zgliczynski proof-machinery build (``#668``-``#674``).  #674
certified a SINGLE-return covering only at a near-vacuous h-set size (ru=1e-8); at the
more meaningful ru=1e-6 it fails (ratio 0.28).  The genuine fix is to COMPOSE several
section returns: a hyperbolic map's stretch compounds multiplicatively (``P^N`` stretch
``~ lambda^N``), and the composed correlation-preserving image is built by the chain
rule from the per-leg Jacobian enclosures with #674's one-linearization discipline.

Before trusting any composed image on the real regularized CR3BP map, the composition
MACHINERY itself is validated here on maps whose N-fold stretch is known in closed form:

  * exactness -- an affine map with constant diagonal Jacobian ``diag(lambda_u, lambda_s)``
    composed N times has composed Jacobian ``diag(lambda_u^N, lambda_s^N)`` EXACTLY, and
    the composed mean-value image of a thin h-set separates by exactly ``2 ru lambda_u^N``
    with ~zero width (mean value is exact for affine maps, chain rule is exact for
    constant Jacobians);
  * soundness -- when each leg Jacobian carries an interval width (mimicking the genuine
    ``[DP]`` variation of the real flyby map), the composed enclosure rigorously CONTAINS
    the true composed image of any realization inside the leg intervals;
  * the mechanism, both directions -- if the per-leg enclosures stay TIGHT the composed
    edge separation still exceeds the composed width out to N=3 (covering certifies at
    every N), whereas if the legs are LOOSE the over-approximation compounds faster than
    the signal and the composed covering degrades below ratio 1 (covering fails) -- the
    precise ``why the real W-Z proof needs the enclosures kept tight`` fact.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#675)",
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


def _diag(a: float, d: float) -> list[list[Any]]:
    return [[_pt(a), _pt(0.0)], [_pt(0.0), _pt(d)]]


# --------------------------------------------------------------------------- #
# 1. compose_section_jacobians -- chain-rule product, exact for constant legs. #
# --------------------------------------------------------------------------- #
def test_compose_single_leg_is_identity() -> None:
    """One leg -> the composed Jacobian equals that leg's Jacobian unchanged."""
    iv = mp.iv
    dp = [[_pt(3.0), _pt(0.5)], [_pt(-2.0), _pt(4.0)]]
    comp = vti.compose_section_jacobians(iv, [dp])
    for i in range(2):
        for j in range(2):
            assert bool(comp[i][j].a <= dp[i][j].a) and bool(comp[i][j].b >= dp[i][j].b)
            assert _w(comp[i][j]) < 1e-30


def test_compose_affine_exact_stretch() -> None:
    """N constant diagonal legs -> composed Jacobian is diag(lambda_u^N, lambda_s^N) exact."""
    iv = mp.iv
    lam_u, lam_s = 3.0, 0.5  # dyadic -> binary-exact, so the enclosure is tight (width 0)
    for n in (1, 2, 3):
        legs = [_diag(lam_u, lam_s) for _ in range(n)]
        comp = vti.compose_section_jacobians(iv, legs)
        assert bool(comp[0][0].a <= _pt(lam_u**n)) and bool(comp[0][0].b >= _pt(lam_u**n))
        assert bool(comp[1][1].a <= _pt(lam_s**n)) and bool(comp[1][1].b >= _pt(lam_s**n))
        assert _w(comp[0][0]) < 1e-25 and _w(comp[1][1]) < 1e-25  # exact, no over-approximation
        assert bool(comp[0][1].b < 1e-30) and bool(comp[0][1].a > -1e-30)  # off-diagonal stays 0


def test_compose_order_is_last_leg_outermost() -> None:
    """Composition applies legs in flow order: DP_2 @ DP_1 (not DP_1 @ DP_2)."""
    iv = mp.iv
    # non-commuting: dp1 shears u<-s, dp2 shears s<-u; the two orders differ.
    dp1 = [[_pt(1.0), _pt(1.0)], [_pt(0.0), _pt(1.0)]]
    dp2 = [[_pt(1.0), _pt(0.0)], [_pt(1.0), _pt(1.0)]]
    comp = vti.compose_section_jacobians(iv, [dp1, dp2])  # == dp2 @ dp1
    # dp2 @ dp1 = [[1,1],[1,2]]
    assert abs(float(comp[0][0].mid) - 1.0) < 1e-30
    assert abs(float(comp[0][1].mid) - 1.0) < 1e-30
    assert abs(float(comp[1][0].mid) - 1.0) < 1e-30
    assert abs(float(comp[1][1].mid) - 2.0) < 1e-30


def test_compose_empty_raises() -> None:
    with pytest.raises(ValueError):
        vti.compose_section_jacobians(mp.iv, [])


# --------------------------------------------------------------------------- #
# 2. Composed mean-value image -- exact stretch, and soundness under widths.   #
# --------------------------------------------------------------------------- #
def _composed_edges(
    iv: Any, legs: list[list[list[Any]]], ru: float, rs: float
) -> tuple[Any, Any, Any]:
    """Composed mean-value images (whole box, left edge, right edge) of a thin h-set.

    N is axis-aligned (unstable axis = coord 0, stable = coord 1), half-widths (ru, rs),
    centre image at the origin.  Uses the SAME primitives as the real driver: the
    chain-rule composed Jacobian fed to ``section_map_meanvalue_image``.
    """
    center = [_pt(0.0), _pt(0.0)]
    jac = vti.compose_section_jacobians(iv, legs)

    def img(aval: Any) -> tuple[Any, Any]:
        off = [aval * iv.mpf(ru), iv.mpf([-1, 1]) * iv.mpf(rs)]
        out = vti.section_map_meanvalue_image(iv, center, jac, off)
        return out[0], out[1]

    return img(iv.mpf([-1, -1])), img(iv.mpf([1, 1])), img(iv.mpf([-1, 1]))


def test_composed_image_affine_exact_separation() -> None:
    """Affine diagonal map: composed edge separation is exactly 2 ru lambda_u^N, ~zero width."""
    iv = mp.iv
    lam_u, lam_s = 3.0, 0.5
    ru, rs = 1e-6, 1e-8
    for n in (1, 2, 3):
        legs = [_diag(lam_u, lam_s) for _ in range(n)]
        left, right, _ = _composed_edges(iv, legs, ru, rs)
        sep = abs(float(left[0].mid) - float(right[0].mid))
        assert abs(sep - 2.0 * ru * lam_u**n) < 1e-18  # exact geometric stretch
        assert _w(left[0]) < 1e-18 and _w(right[0]) < 1e-18  # mean value exact for affine


def test_composed_enclosure_is_sound_under_leg_widths() -> None:
    """Interval-width legs: the composed image contains the true image of any realisation."""
    iv = mp.iv
    eps = 0.3
    lam_u, lam_s = 3.0, 0.5
    ru, rs = 1e-6, 1e-8
    legs = [
        [[iv.mpf([lam_u - eps, lam_u + eps]), _pt(0.0)], [_pt(0.0), _pt(lam_s)]] for _ in range(3)
    ]
    left, right, _ = _composed_edges(iv, legs, ru, rs)
    # a concrete realisation (every leg at lam_u) -> true composed left-edge coord0
    true_left0 = -(lam_u**3) * ru
    true_right0 = (lam_u**3) * ru
    assert bool(left[0].a <= _pt(true_left0)) and bool(left[0].b >= _pt(true_left0))
    assert bool(right[0].a <= _pt(true_right0)) and bool(right[0].b >= _pt(true_right0))


# --------------------------------------------------------------------------- #
# 3. Mechanism (both directions): tight legs keep covering; loose legs lose it. #
# --------------------------------------------------------------------------- #
def _covers_at_n(iv: Any, legs: list[list[list[Any]]], ru: float, rs: float) -> dict[str, Any]:
    left, right, box = _composed_edges(iv, legs, ru, rs)
    sep = abs(float(left[0].mid) - float(right[0].mid))
    edgew = max(_w(left[0]), _w(right[0]))
    ratio = sep / edgew if edgew > 0 else float("inf")
    hu = 0.4 * (sep - edgew) if sep > edgew else max(1e-12, 0.25 * sep)
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
    res["ratio"] = ratio
    return res


def test_tight_legs_keep_covering_through_three_returns() -> None:
    """Small per-leg Jacobian width: composed separation still beats width out to N=3."""
    iv = mp.iv
    eps = 0.01  # tight legs (intermediate enclosures kept tight, as #669's QR reframing does)
    lam_u, lam_s = 3.0, 0.5
    ru, rs = 1e-6, 1e-8
    for n in (1, 2, 3):
        legs = [
            [[iv.mpf([lam_u - eps, lam_u + eps]), _pt(0.0)], [_pt(0.0), _pt(lam_s)]]
            for _ in range(n)
        ]
        res = _covers_at_n(iv, legs, ru, rs)
        assert res["ratio"] > 1.0, f"tight legs should keep ratio>1 at N={n}"
        assert res["covers"] is True, f"tight legs should still cover at N={n}"


def _loose_legs(iv: Any, eps: float, lam_u: float, lam_s: float, n: int) -> list[list[list[Any]]]:
    return [
        [[iv.mpf([lam_u - eps, lam_u + eps]), _pt(0.0)], [_pt(0.0), _pt(lam_s)]] for _ in range(n)
    ]


def test_loose_legs_break_covering_and_overapprox_compounds() -> None:
    """Wide per-leg Jacobian (uncertain-sign stretch) -> covering cannot certify, and the
    absolute over-approximation compounds with each composed return.

    ``eps > lambda_u`` makes each leg's unstable-stretch enclosure straddle zero -- the
    per-leg enclosure is NOT tight -- so the mean-value edge width exceeds the separation
    and the SAME checker correctly reports non-covering at every return.  Crucially the
    edge over-approximation width GROWS with N (the wrapping compounds), the honest reason
    a real W-Z chain must keep every intermediate h-set enclosure tight (the QR reframing).
    """
    iv = mp.iv
    eps, lam_u, lam_s = 4.0, 3.0, 0.5  # [lam_u-eps, lam_u+eps] = [-1, 7] straddles 0
    ru, rs = 1e-6, 1e-8
    widths = []
    for n in (1, 2, 3):
        res = _covers_at_n(iv, _loose_legs(iv, eps, lam_u, lam_s, n), ru, rs)
        assert res["ratio"] < 1.0, f"loose legs should have ratio<1 at N={n}"
        assert res["covers"] is False, f"loose legs should not cover at N={n}"
        left, right, _ = _composed_edges(iv, _loose_legs(iv, eps, lam_u, lam_s, n), ru, rs)
        widths.append(max(_w(left[0]), _w(right[0])))
    assert widths[2] > widths[1] > widths[0]  # over-approximation compounds with returns


def test_leg_tightness_is_load_bearing_at_three_returns() -> None:
    """Same lambda and N=3: TIGHT legs cover, LOOSE legs do not -- isolating leg-enclosure
    tightness as the load-bearing ingredient of the composed covering (parallel to #674)."""
    iv = mp.iv
    lam_u, lam_s = 3.0, 0.5
    ru, rs = 1e-6, 1e-8
    tight = [
        [[iv.mpf([lam_u - 0.01, lam_u + 0.01]), _pt(0.0)], [_pt(0.0), _pt(lam_s)]] for _ in range(3)
    ]
    loose = _loose_legs(iv, 4.0, lam_u, lam_s, 3)
    assert _covers_at_n(iv, tight, ru, rs)["covers"] is True
    assert _covers_at_n(iv, loose, ru, rs)["covers"] is False
