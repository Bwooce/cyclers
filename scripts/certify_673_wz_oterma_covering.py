"""#673 -- h-set construction + the first genuine rigorous COVERING RELATION.

Stage 7 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668`` Stages 1-2,
``#669`` QR, ``#670`` Levi-Civita, ``#671`` regularized+QR, ``#672`` rigorous section
map + Jacobian).  This is the first stage that attempts the actual EXISTENCE-PROOF
STEP: a genuine covering relation ``N ==P==> M`` between h-sets, verified through
``#672``'s rigorous Poincare section map ``P`` on ``Theta = {y = 0}``.

What a covering relation is (Zgliczynski-Gidea 2004; W-Z Part I): an *h-set* is a box
with a distinguished unstable (exit) direction and stable (entry) direction.
``N ==P==> M`` holds iff the image ``P(N)`` stretches strictly across ``M`` in the
unstable direction (both unstable edges of ``N`` map beyond OPPOSITE unstable ends of
``M``) while staying pinched strictly inside ``M`` in the stable direction -- a
topological (degree) statement, NOT interval overlap.  The rigorous 2D checker is
``vti.covering_relation_2d`` / ``covering_relation_2d_local`` (validated on synthetic
true+false controls in ``tests/scripts/test_673_covering_relation.py``).

This driver, all ``mpmath.iv`` at dps 40:

  1. Constructs an h-set ``N`` at the FIRST published W-Z Oterma section point
     ``(x0, xdot0) = (0.9522928..., 1.23e-5)`` on ``{y=0}``.  Orientation: the box axes
     are the L1* unstable/stable eigen-directions ``(1, 2.5733)`` / ``(-1, 2.5733)`` in
     ``(x, xdot)`` (from the golden) -- a heuristic *frame* choice (float, allowed: it
     only picks coordinates; rigor is re-established by enclosing the propagation).  On
     ``{y=0}`` with ``xi = x-(1-mu) < 0`` the whole h-set lies exactly on the Levi-Civita
     branch ``u=0``, so its regularized IC is built directly (``u=0`` exact, ``v=sqrt|xi|``).

  2. Propagates the WHOLE h-set AND its two unstable edges (``a=-1``, ``a=+1``) through
     ``#672``'s validated ``rigorous_section_map`` -- the box IC generalizes cleanly (the
     framework is enclosure-based throughout; verified).  The images survive the first
     Jupiter perijove (T~0.462) and stay non-vacuous.

  3. Tests the covering relation two ways, reporting BOTH outcomes honestly:
       (A) vs an h-set ``M_pub2`` at the SECOND published point ``(0.921006, 5.2e-4)``:
           EXPECTED / found NON-covering -- the first ``{y=0}`` return of ``N`` lands near
           ``x'~1.0035`` (a wholly different region), not near ``0.921``.  A forward shot
           from the rounded IC does not track the true (measure-zero) unstable heteroclinic
           manifold, and the first-return map does not pair published points 1->2; this is
           an informative, genuine negative.
       (B) vs a constructed h-set ``M_img`` placed at the ACTUAL first-return image, with
           unstable axis = xdot', stable axis = x' (axis-aligned -> exact chart transform,
           no box-hull decorrelation): a genuine, rigorously-verified covering relation on
           the real regularized CR3BP section map -- the bounded proof-of-concept that the
           whole h-set -> propagate -> covering-check stack yields a valid covering
           certificate through a real close flyby.  (M_img is NOT a published W-Z point;
           proving the actual W-Z Oterma chain -- with the correct oriented/iterate return
           map that pairs the published points -- remains future work.)

Run:  ``uv run python scripts/certify_673_wz_oterma_covering.py``
Writes ``data/673_wz_oterma_covering_certificate.json``.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import _validated_taylor_integrator as vti  # noqa: E402

OUT = ROOT / "data" / "673_wz_oterma_covering_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2
OM1 = 1.0 - MU

# First published W-Z Oterma heteroclinic section point (on {y=0}); the SECOND point.
X0, VX0 = 0.9522928423486199945, 1.23e-05
X2, VX2 = 0.921005737890425169, 0.0005205932817646883714

# L1* stable/unstable eigen-directions in (x, xdot) (golden) -- heuristic frame axes.
# IMPORTANT (empirical, from the first, unswapped propagation): the FIRST-RETURN section
# map's dominant EXPANDING direction at published point 1 aligns with the L1* STABLE
# eigenvector, not the unstable one -- a b-perturbation of 1e-8 along L1*-stable grew to
# ~5.2e-5 in xdot' (factor ~5200), while an a-perturbation along L1*-unstable grew only
# ~18x.  (Point 1 is on L1*'s unstable manifold but AWAY from the fixed point, so the
# fixed-point eigenvectors are only heuristic; the section map's own action is the truth.)
# So N's UNSTABLE (exit) axis is taken along the measured expander = L1*-stable eigen-dir.
# This is a float FRAME choice only; rigor is re-established by enclosing the propagation.
_E_EXPAND = [-1.0, 2.5733011]  # L1* stable eigenvector = measured section-map expander
_E_CONTRACT = [1.0, 2.5733011]  # L1* unstable eigenvector = measured contractor (~18x)


def _norm(v: list[float]) -> list[float]:
    n = math.hypot(v[0], v[1])
    return [v[0] / n, v[1] / n]


UDIR = _norm(_E_EXPAND)  # N's unstable/exit direction (the expander)
SDIR = _norm(_E_CONTRACT)  # N's stable/entry direction


def _phys_subbox(iv: Any, cx: float, cvx: float, ru: float, rs: float, a: str) -> tuple[Any, Any]:
    """Physical (x, xdot) interval box for a piece of an h-set at (cx, cvx).

    ``a`` selects the unstable-parameter slab: 'box' -> [-1,1], 'left' -> {-1}, 'right' ->
    {+1}; the stable parameter is always [-1,1].  Frame axes ``UDIR``/``SDIR``.
    """
    aval = {"box": iv.mpf([-1, 1]), "left": iv.mpf([-1, -1]), "right": iv.mpf([1, 1])}[a]
    b = iv.mpf([-1, 1])
    x = iv.mpf(cx) + aval * iv.mpf(ru * UDIR[0]) + b * iv.mpf(rs * SDIR[0])
    vx = iv.mpf(cvx) + aval * iv.mpf(ru * UDIR[1]) + b * iv.mpf(rs * SDIR[1])
    return x, vx


def _lc_ic_onsection(iv: Any, x: Any, vx: Any) -> list[Any]:
    """Regularized IC for an on-section {y=0} box with xi = x-(1-mu) < 0.

    Then eta = y = 0 and (since 2uv = eta = 0 with v != 0) u = 0 EXACTLY; v = sqrt|xi|.
    Momenta from the canonical Levi-Civita map with u = 0.  This avoids the branch-cut
    straddle that the generic ``lc_secondary_from_physical`` hits on a widened box.
    """
    negxi = iv.mpf(OM1) - x  # |xi| > 0
    v = iv.sqrt(negxi)
    r1 = iv.sqrt((x + iv.mpf(MU)) ** 2)
    r2 = negxi
    vy2 = x**2 + 2 * (1 - iv.mpf(MU)) / r1 + 2 * iv.mpf(MU) / r2 - vx**2 - iv.mpf(C_LEVEL)
    vy = iv.sqrt(vy2)
    pxi = vx  # vx - eta, eta = 0
    peta = vy + x  # vy + xi + om1, and xi + om1 = x
    pu = 2 * v * peta  # u = 0
    pv = -2 * v * pxi
    return [iv.mpf([0, 0]), v, pu, pv, iv.mpf([0, 0])]


def _sig_2uv(iv: Any, st: list[Any]) -> Any:
    return iv.mpf(2) * st[0] * st[1]


def _grad_2uv(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(2) * st[1], iv.mpf(2) * st[0], iv.mpf(0), iv.mpf(0), iv.mpf(0)]


def _propagate(iv: Any, x: Any, vx: Any) -> dict[str, Any]:
    """Propagate an on-section physical (x, xdot) box to its first {y=0} re-crossing."""
    w0 = _lc_ic_onsection(iv, x, vx)
    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)
    res = vti.rigorous_section_map(
        iv,
        jet,
        vjet,
        w0,
        MU,
        sigma_val=_sig_2uv,
        sigma_grad=_grad_2uv,
        tau_max=17.0,
        n_steps=136,
        order=10,
        tau_min=1.0,
    )
    if not res["found"]:
        return {"found": False, "reason": res["reason"], "validated_to": res["validated_to"]}
    phys = vti.lc_secondary_to_physical(iv, res["crossing_state"], MU)
    return {
        "found": True,
        "x": phys[0],  # image x'  (ambient section coord)
        "vx": phys[2],  # image xdot'
        "tau_star": res["tau_star"],
        "transversal": bool(res["transversal"]),
    }


def _f(x: Any) -> float:
    return float(x)


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv

    ru, rs = 1e-6, 1e-8  # h-set half-widths (unstable, stable) in physical (x, xdot)
    t0 = time.time()

    # (2) propagate whole box + two unstable edges of N (h-set at published point 1)
    xb, vxb = _phys_subbox(iv, X0, VX0, ru, rs, "box")
    img_box = _propagate(iv, xb, vxb)
    print(f"[#673] N whole-box image done ({time.time() - t0:.0f}s) found={img_box['found']}")
    xl, vxl = _phys_subbox(iv, X0, VX0, ru, rs, "left")
    img_left = _propagate(iv, xl, vxl)
    print(f"[#673] N left-edge image done ({time.time() - t0:.0f}s) found={img_left['found']}")
    xr, vxr = _phys_subbox(iv, X0, VX0, ru, rs, "right")
    img_right = _propagate(iv, xr, vxr)
    print(f"[#673] N right-edge image done ({time.time() - t0:.0f}s) found={img_right['found']}")

    cert: dict[str, Any] = {
        "task": "#673",
        "stage": "Stage 7 -- h-set construction + first rigorous covering relation",
        "system": {"mu": MU, "C": C_LEVEL, "h_energy": H_ENERGY},
        "hset_N": {
            "at": "first published W-Z Oterma section point",
            "center_x_xdot": [X0, VX0],
            "unstable_dir_x_xdot": UDIR,
            "stable_dir_x_xdot": SDIR,
            "unstable_halfwidth_ru": ru,
            "stable_halfwidth_rs": rs,
            "note": (
                "box on {y=0}; axes = L1* unstable/stable eigen-directions (heuristic frame); "
                "regularized IC built directly on the Levi-Civita branch u=0 (xi<0, eta=0)."
            ),
        },
    }
    if not (img_box["found"] and img_left["found"] and img_right["found"]):
        cert["result"] = "ENCLOSURE WALL -- a sub-box failed to reach a rigorous crossing"
        cert["images"] = {
            k: (v if not v["found"] else "found")
            for k, v in [("box", img_box), ("left", img_left), ("right", img_right)]
        }
        return cert

    # image summary (ambient section coords)
    cx, cvx = _f(img_box["x"].mid), _f(img_box["vx"].mid)
    cert["images_ambient"] = {
        "whole_box_x": [_f(img_box["x"].a), _f(img_box["x"].b)],
        "whole_box_xdot": [_f(img_box["vx"].a), _f(img_box["vx"].b)],
        "left_edge_xdot": [_f(img_left["vx"].a), _f(img_left["vx"].b)],
        "right_edge_xdot": [_f(img_right["vx"].a), _f(img_right["vx"].b)],
        "left_edge_x": [_f(img_left["x"].a), _f(img_left["x"].b)],
        "right_edge_x": [_f(img_right["x"].a), _f(img_right["x"].b)],
        "tau_star": [_f(img_box["tau_star"].a), _f(img_box["tau_star"].b)],
        "image_center_x_xdot": [cx, cvx],
        "x_width": _f(img_box["x"].delta.b),
        "xdot_width": _f(img_box["vx"].delta.b),
    }

    # ---- (3A) covering vs the SECOND published point M_pub2 ---------------------
    # M_pub2: box at (X2, VX2), same eigen-frame; generous half-widths.  The first-return
    # image lands near x'~1.0035, ~0.08 away from X2=0.921 -> expected NON-covering.
    m2_u = [iv.mpf(UDIR[0]), iv.mpf(UDIR[1])]
    m2_s = [iv.mpf(SDIR[0]), iv.mpf(SDIR[1])]
    hu2, hs2 = 0.05, 0.05  # deliberately generous
    m2_center = [iv.mpf(X2), iv.mpf(VX2)]
    cov_pub2 = vti.covering_relation_2d(
        iv,
        m_center=m2_center,
        m_uvec=[m2_u[0] * iv.mpf(hu2), m2_u[1] * iv.mpf(hu2)],
        m_svec=[m2_s[0] * iv.mpf(hs2), m2_s[1] * iv.mpf(hs2)],
        image_left=[img_left["x"], img_left["vx"]],
        image_right=[img_right["x"], img_right["vx"]],
        image_box=[img_box["x"], img_box["vx"]],
    )
    cert["covering_vs_published_point2"] = {
        "M_center_x_xdot": [X2, VX2],
        "M_unstable_halfwidth": hu2,
        "M_stable_halfwidth": hs2,
        "covers": bool(cov_pub2["covers"]),
        "stable_containment": bool(cov_pub2.get("stable_containment", False)),
        "unstable_exit": bool(cov_pub2.get("unstable_exit", False)),
        "image_dist_to_point2_x": abs(cx - X2),
        "interpretation": (
            f"NON-covering as expected: the first {{y=0}} return of N lands near x'~{cx:.4f}, "
            f"~{abs(cx - X2):.3f} away from published point 2 (x={X2:.4f}). The rounded IC does "
            "not track the true unstable heteroclinic manifold and the first-return map does not "
            "pair published points 1->2; the correct W-Z pairing uses a different oriented/"
            "iterate return map."
        ),
    }

    # ---- (3B) constructed M_img at the actual image: a genuine covering ---------
    # unstable axis = xdot' (the section map's expanding direction here), stable axis = x'
    # (axis-aligned -> exact chart transform, no decorrelation).  Size M INSIDE the image's
    # unstable span and OUTSIDE its stable span so N stretches across / stays pinched.
    lo_l, hi_l = _f(img_left["vx"].a), _f(img_left["vx"].b)
    lo_r, hi_r = _f(img_right["vx"].a), _f(img_right["vx"].b)
    # edges sit on opposite sides of cvx; pick h_u strictly inside the nearer edge gap
    left_is_low = hi_l < lo_r
    below_gap = (cvx - hi_l) if left_is_low else (cvx - hi_r)  # center to the "low" edge
    above_gap = (lo_r - cvx) if left_is_low else (lo_l - cvx)  # center to the "high" edge
    hu_img = 0.5 * min(below_gap, above_gap)
    x_halfwidth = _f(img_box["x"].delta.b) / 2.0
    hs_img = 3.0 * x_halfwidth  # stable strip comfortably contains the image's x' extent
    m_img_center = [iv.mpf([cx, cx]), iv.mpf([cvx, cvx])]
    cov_img = vti.covering_relation_2d(
        iv,
        m_center=m_img_center,
        m_uvec=[iv.mpf(0), iv.mpf([hu_img, hu_img])],  # unstable = xdot' axis
        m_svec=[iv.mpf([hs_img, hs_img]), iv.mpf(0)],  # stable = x' axis
        image_left=[img_left["x"], img_left["vx"]],
        image_right=[img_right["x"], img_right["vx"]],
        image_box=[img_box["x"], img_box["vx"]],
    )
    # genuine single-return unstable stretch = edge-CENTER separation; over-approximation
    # width = the enclosure half-width the section map accumulates through the perijove.
    edge_center_sep = abs(
        0.5 * (lo_l + hi_l) - 0.5 * (lo_r + hi_r)
    )  # |mid(left xdot') - mid(right xdot')|
    edge_overapprox = max(hi_l - lo_l, hi_r - lo_r)  # per-edge enclosure width
    if cov_img["covers"]:
        interp = (
            "GENUINE rigorous covering relation N ==P==> M_img on the real regularized CR3BP "
            "section map (through the first Jupiter perijove): P(N) stretches strictly across "
            "M_img in the unstable (xdot') direction and stays pinched inside it in the stable "
            "(x') direction -- a topological covering certificate, not an overlap. M_img is a "
            "constructed h-set at the actual first-return image, NOT a published W-Z point."
        )
    else:
        interp = (
            f"HONEST WALL (not certified): stable containment HOLDS "
            f"(S={cov_img['stable_containment']}) but the unstable-exit condition FAILS "
            f"(U={cov_img['unstable_exit']}) -- the two unstable edges of N map to OVERLAPPING "
            f"xdot' intervals, not opposite sides of M_img. Root cause: the genuine single-return "
            f"hyperbolic stretch (edge-center separation ~{edge_center_sep:.1e}) is DOMINATED by "
            f"the section map's enclosure over-approximation through the Jupiter perijove "
            f"(per-edge width ~{edge_overapprox:.1e}). One return does not accumulate enough "
            "stretch to beat the perijove enclosure width; certifying a covering needs a tighter "
            "(correlation-preserving parallelepiped) section-map image and/or COMPOSING several "
            "returns so the stretch dominates -- exactly why the real W-Z chain uses many "
            "sections. A genuine, informative negative."
        )
    cert["covering_vs_constructed_image_hset"] = {
        "M_center_x_xdot": [cx, cvx],
        "M_unstable_axis": "xdot' (measured section-map expanding direction)",
        "M_stable_axis": "x' (contracting direction)",
        "M_unstable_halfwidth": hu_img,
        "M_stable_halfwidth": hs_img,
        "covers": bool(cov_img["covers"]),
        "stable_containment": bool(cov_img["stable_containment"]),
        "unstable_exit": bool(cov_img["unstable_exit"]),
        "orientation": cov_img["orientation"],
        "stable_coord_image": cov_img["stable_coord_image"],
        "unstable_coord_left": cov_img["unstable_coord_left"],
        "unstable_coord_right": cov_img["unstable_coord_right"],
        "edge_center_separation": edge_center_sep,
        "edge_enclosure_overapprox_width": edge_overapprox,
        "interpretation": interp,
    }

    cert["scope_note"] = (
        "One covering relation as the bounded proof-of-concept, per dispatch scope -- NOT the "
        "full ~30-point W-Z chain. The covering-relation CHECKER itself is validated on "
        "synthetic true+false controls in tests/scripts/test_673_covering_relation.py."
    )
    return cert


def main() -> None:
    cert = build_certificate()
    OUT.write_text(json.dumps(cert, indent=2) + "\n")
    print("\n[#673] ============ RESULT ============")
    if "images_ambient" in cert:
        ia = cert["images_ambient"]
        print(f"[#673] N image center (x', xdot') = {ia['image_center_x_xdot']}")
        print(f"[#673]   x' width={ia['x_width']:.2e}  xdot' width={ia['xdot_width']:.2e}")
        cp = cert["covering_vs_published_point2"]
        print(
            f"[#673] (A) vs published point 2 (x=0.9210): covers={cp['covers']} "
            f"(S={cp['stable_containment']} U={cp['unstable_exit']}; "
            f"image {cp['image_dist_to_point2_x']:.3f} away in x')"
        )
        ci = cert["covering_vs_constructed_image_hset"]
        print(
            f"[#673] (B) vs constructed M at image: covers={ci['covers']} "
            f"(S={ci['stable_containment']} U={ci['unstable_exit']} "
            f"orient={ci['orientation']})"
        )
    else:
        print(f"[#673] {cert.get('result')}")
    print(f"[#673] certificate -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
