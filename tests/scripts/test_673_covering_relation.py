"""#673 -- tests for the rigorous 2D COVERING RELATION checker (the proof step).

Stage 7 of the Wilczak-Zgliczynski proof-machinery build (``#668``-``#672``).  The
W-Z existence proof rests on *covering relations* ``N ==f==> M`` between *h-sets*
verified through the Poincare section map.  ``covering_relation_2d`` is the checker
for the planar (u=s=1) case: given rigorous enclosures of the images of ``N``'s two
unstable edges and of the whole set ``N``, it decides -- in ``M``'s own coordinates
-- whether

  (S) the whole image is pinched strictly inside ``M`` in the STABLE direction, AND
  (U) the two unstable edges map strictly beyond OPPOSITE unstable ends of ``M``.

``(S) and (U)`` is the ZG "correctly aligned" sufficient condition for a covering
relation (a topological/degree statement), NOT an overlap test.

The load-bearing validation here (mandated by the dispatch): the checker gives the
CORRECT answer on synthetic maps where covering is OBVIOUSLY true or OBVIOUSLY false
by hand -- so it is neither trivially "always yes" nor "always no".  Positive control:
a hyperbolic diagonal map genuinely covers.  Negative controls: (a) a map that
EXPANDS the stable direction (cannot be pinched -> must fail (S)); (b) a map that does
NOT stretch the unstable direction past ``M`` (must fail (U)); (c) a map that shifts
the image entirely off ``M`` to one side (unstable edges on the SAME side -> fail (U)).
A rotated-frame positive control exercises the rigorous frame-inverse path.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#673)",
)

import scripts._validated_taylor_integrator as vti  # noqa: E402


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _matvec(iv: Any, a: list[list[float]], v: list[Any]) -> list[Any]:
    return [iv.mpf(a[i][0]) * v[0] + iv.mpf(a[i][1]) * v[1] for i in range(2)]


def _affine_images(
    iv: Any,
    amat: list[list[float]],
    shift: list[float],
    center: list[Any],
    uvec: list[Any],
    svec: list[Any],
) -> tuple[list[Any], list[Any], list[Any]]:
    """EXACT image enclosures of N's left edge, right edge, whole box under an affine map.

    N = { center + a*uvec + b*svec : a,b in [-1,1] };  f(p) = amat @ p + shift.  Because
    the map is affine, f(N) = f(center) + a*(amat@uvec) + b*(amat@svec); each parameter
    appears once, so this interval evaluation is TIGHT (no decorrelation), exactly the
    kind of enclosure the rigorous section map delivers for the real case.  Left edge is
    {a=-1}, right {a=+1}, box the full a,b in [-1,1].
    """
    fc = [_matvec(iv, amat, center)[i] + iv.mpf(shift[i]) for i in range(2)]
    fu = _matvec(iv, amat, uvec)
    fs = _matvec(iv, amat, svec)
    ab = iv.mpf([-1, 1])
    m1, p1 = iv.mpf([-1, -1]), iv.mpf([1, 1])

    def img(a: Any) -> list[Any]:
        return [fc[i] + a * fu[i] + ab * fs[i] for i in range(2)]

    return img(m1), img(p1), img(ab)


# Unit-box h-set used as both N (implicitly, via image construction) and M.
_UNIT = {
    "center": [_pt(0.0), _pt(0.0)],
    "uvec": [_pt(1.0), _pt(0.0)],
    "svec": [_pt(0.0), _pt(1.0)],
}


def _diag(lam: float, mu: float, shift: float = 0.0) -> tuple[list[list[float]], list[float]]:
    return [[lam, 0.0], [0.0, mu]], [shift, 0.0]


def _check(amat: list[list[float]], shift: list[float], n_hset: dict[str, Any]) -> dict[str, Any]:
    iv = mp.iv
    left, right, box = _affine_images(
        iv, amat, shift, n_hset["center"], n_hset["uvec"], n_hset["svec"]
    )
    return vti.covering_relation_2d(
        iv,
        m_center=_UNIT["center"],
        m_uvec=_UNIT["uvec"],
        m_svec=_UNIT["svec"],
        image_left=left,
        image_right=right,
        image_box=box,
    )


# --------------------------------------------------------------------------- #
# POSITIVE control: hyperbolic diagonal map genuinely covers.                   #
# --------------------------------------------------------------------------- #
def test_positive_control_hyperbolic_covers() -> None:
    """P(x,y)=(3x, y/3), N=M=unit box: stretches unstable x3, pinches stable /3 -> COVERS."""
    res = _check(*_diag(3.0, 1.0 / 3.0), _UNIT)
    assert res["covers"] is True
    assert res["stable_containment"] is True
    assert res["unstable_exit"] is True
    assert res["orientation"] == "preserving"
    # stable image is [-1/3, 1/3] strictly inside (-1,1); unstable edges at -+3
    assert res["stable_coord_image"][0] > -0.4 and res["stable_coord_image"][1] < 0.4
    assert res["unstable_coord_left"][1] < -1 and res["unstable_coord_right"][0] > 1


def test_positive_control_orientation_reversing_covers() -> None:
    """P(x,y)=(-3x, y/3): unstable edges swap sides -> still covers (reversing)."""
    res = _check(*_diag(-3.0, 1.0 / 3.0), _UNIT)
    assert res["covers"] is True
    assert res["orientation"] == "reversing"


def test_positive_control_rotated_frame_covers_via_local() -> None:
    """Hyperbolic map with N=M rotated 45deg -> covers, checked in M's OWN coordinates.

    N=M is the unit box rotated 45deg (unstable dir (1,1)/sqrt2, stable (-1,1)/sqrt2).
    The map is  A = B diag(3, 1/3) B^{-1}  (B=[uvec|svec] orthogonal) -- diag(3,1/3) in
    the LOCAL frame.  So in M-local coordinates the images are exact and tight:
    left edge (a=-1) -> unstable coord -3; right edge -> +3; whole box stable coord
    (1/3)*[-1,1] = [-1/3, 1/3].  This is the correlation-preserving path
    (``covering_relation_2d_local``); the ambient box-hull wrapper would spuriously
    inflate the stable coord under the rotation (documented caveat).
    """
    iv = mp.iv
    au_left = iv.mpf([-3, -3])
    au_right = iv.mpf([3, 3])
    as_box = iv.mpf(1) / iv.mpf(3) * iv.mpf([-1, 1])
    res = vti.covering_relation_2d_local(iv, au_left=au_left, au_right=au_right, as_box=as_box)
    assert res["covers"] is True
    assert res["stable_containment"] is True and res["unstable_exit"] is True
    assert res["orientation"] == "preserving"


def test_ambient_boxhull_wrapper_decorrelates_under_rotation() -> None:
    """DOCUMENTED caveat: the ambient box-hull wrapper can spuriously FAIL for rotated M.

    Same genuinely-covering rotated map as above, but fed through the ambient box-hull
    wrapper: the box hull of the sheared image decorrelates and inflates the stable
    coord, so the wrapper reports covers=False.  This guards the caveat -- a False from
    the wrapper is NOT proof of non-covering (a True always would be sound).
    """
    iv = mp.iv
    cf = float(1.0 / mp.sqrt(2))
    uvec = [_pt(cf), _pt(cf)]
    svec = [_pt(-cf), _pt(cf)]
    c2 = cf * cf
    d = 3.0 * c2 + (1.0 / 3.0) * c2
    o = 3.0 * c2 - (1.0 / 3.0) * c2
    amat = [[d, o], [o, d]]
    center = [_pt(0.0), _pt(0.0)]
    left, right, box = _affine_images(iv, amat, [0.0, 0.0], center, uvec, svec)
    res = vti.covering_relation_2d(
        iv,
        m_center=center,
        m_uvec=uvec,
        m_svec=svec,
        image_left=left,
        image_right=right,
        image_box=box,
    )
    assert res["covers"] is False  # box-hull decorrelation artifact, not true non-covering
    assert res["stable_containment"] is False


# --------------------------------------------------------------------------- #
# NEGATIVE controls: covering is OBVIOUSLY false by hand.                        #
# --------------------------------------------------------------------------- #
def test_negative_control_stable_expands_fails() -> None:
    """P(x,y)=(3x, 3y): stable direction EXPANDS -> image not pinched -> NOT cover."""
    res = _check(*_diag(3.0, 3.0), _UNIT)
    assert res["covers"] is False
    assert res["stable_containment"] is False  # the failing condition
    assert res["unstable_exit"] is True  # unstable stretch is fine; only (S) fails


def test_negative_control_unstable_too_weak_fails() -> None:
    """P(x,y)=(0.5x, y/3): unstable does NOT reach past M's ends -> NOT cover."""
    res = _check(*_diag(0.5, 1.0 / 3.0), _UNIT)
    assert res["covers"] is False
    assert res["unstable_exit"] is False
    assert res["stable_containment"] is True


def test_negative_control_shifted_off_target_fails() -> None:
    """P(x,y)=(3x+10, y/3): image shoved entirely to +x side -> edges same side -> NOT cover."""
    res = _check(*_diag(3.0, 1.0 / 3.0, shift=10.0), _UNIT)
    assert res["covers"] is False
    assert res["unstable_exit"] is False  # both edges land at +x (same side), no opposite exit
    # both unstable-coord edges are > +1 (7 and 13), i.e. same side
    assert res["unstable_coord_left"][0] > 1 and res["unstable_coord_right"][0] > 1


def test_overlap_alone_does_not_imply_covering() -> None:
    """A map whose image OVERLAPS M heavily but does not stretch across still fails.

    P(x,y)=(0.2x, 0.2y): image is a small box [-0.2,0.2]^2 sitting well inside M --
    huge overlap, but it does NOT cover (no unstable stretch).  Guards against the
    "overlap => covers" fallacy the dispatch explicitly warns about.
    """
    res = _check(*_diag(0.2, 0.2), _UNIT)
    assert res["covers"] is False
    assert res["unstable_exit"] is False
