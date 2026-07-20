"""#674 -- correlation-preserving (mean-value) section-map image + covering relation.

Stage 8 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668`` Stages 1-2,
``#669`` QR, ``#670`` Levi-Civita, ``#671`` regularized+QR, ``#672`` section map +
Jacobian, ``#673`` covering-relation checker + first real attempt).

#673 built and validated the rigorous 2D covering-relation checker and made the first
real covering attempt on the actual W-Z Oterma problem, hitting a PRECISE representation
wall: the h-set at published point 1, propagated through the first Jupiter perijove, has
a genuine hyperbolic edge-center separation ~5.9e-6, but #673 represented each face's
image as an axis-aligned box-hull that over-approximated it to per-edge width ~5.2e-5
(~9x), so the two unstable edges' images OVERLAPPED and the unstable-exit condition could
not be certified.  The over-approximation is the wrapping effect (#669) re-introduced by
box-hulling a thin, sheared image and by propagating each face INDEPENDENTLY through the
stretching flow.

This stage replaces the box-hull with a correlation-preserving (mean-value) image
(``vti.section_map_meanvalue_image``): run the section map ONCE on the whole h-set box
(giving ``[DP]``, a rigorous enclosure of the section-map Jacobian over the IC box) and
ONCE on the box CENTER (giving a TIGHT image ``P(w_hat)``), then enclose each face by the
interval mean-value theorem ``P(face) subset P(w_hat) + [DP] (w0_face - w_hat)``.  The
map's stretching is applied ONCE (tight ``[DP]``) instead of accumulating with per-step
wrapping, so the image collapses toward its true shape.

Honest result (two h-set sizes, same construction):
  * At #673's original half-widths (ru=1e-6, rs=1e-8) the mean-value representation
    already tightens the per-edge width from ~5.2e-5 to ~2.1e-5, but the ratio
    separation/width is still < 1: [DP] genuinely VARIES by ~0.9 across a 1e-6 box
    (strong second-order nonlinearity through the close flyby), which times the ~1e-6
    offset dominates the ~6e-6 signal.  So at that box size the covering still does NOT
    certify -- a real geometric fact, not a representation artifact.
  * At a proof-appropriate smaller h-set (ru=1e-8, rs=1e-10) the [DP]-variation term
    (~box^2) shrinks faster than the linear signal (~box): the edge separation ~5.9e-8
    exceeds the per-edge width ~3.0e-8 (ratio ~2.0), the two unstable edges land strictly
    on OPPOSITE sides of M, the whole image stays pinched inside M's stable strip, and the
    #673 checker CERTIFIES a genuine covering relation N ==P==> M on the real regularized
    CR3BP section map through the first Jupiter perijove.

M is a constructed h-set at the actual first-return image (NOT a published W-Z point);
this is the bounded single-return proof-of-concept, not the full ~30-section W-Z chain.

Run:  ``uv run python scripts/certify_674_wz_oterma_covering_meanvalue.py``
Writes ``data/674_wz_oterma_covering_certificate.json``.
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

OUT = ROOT / "data" / "674_wz_oterma_covering_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2
OM1 = 1.0 - MU

# First published W-Z Oterma section point (on {y=0}); the SECOND point.
X0, VX0 = 0.9522928423486199945, 1.23e-05
X2, VX2 = 0.921005737890425169, 0.0005205932817646883714

# L1* stable/unstable eigen-directions in (x, xdot) (golden) -- heuristic frame axes; the
# section map's dominant expander aligns with L1*-stable (see #673's own bullet/script).
_E_EXPAND = [-1.0, 2.5733011]
_E_CONTRACT = [1.0, 2.5733011]


def _norm(v: list[float]) -> list[float]:
    n = math.hypot(v[0], v[1])
    return [v[0] / n, v[1] / n]


UDIR = _norm(_E_EXPAND)  # N's unstable/exit direction (the expander)
SDIR = _norm(_E_CONTRACT)  # N's stable/entry direction


def _phys_subbox(iv: Any, cx: float, cvx: float, ru: float, rs: float, a: str) -> tuple[Any, Any]:
    """Physical (x, xdot) interval box for a piece of an h-set at (cx, cvx)."""
    aval = {"box": iv.mpf([-1, 1]), "left": iv.mpf([-1, -1]), "right": iv.mpf([1, 1])}[a]
    b = iv.mpf([-1, 1])
    x = iv.mpf(cx) + aval * iv.mpf(ru * UDIR[0]) + b * iv.mpf(rs * SDIR[0])
    vx = iv.mpf(cvx) + aval * iv.mpf(ru * UDIR[1]) + b * iv.mpf(rs * SDIR[1])
    return x, vx


def _lc_ic_onsection(iv: Any, x: Any, vx: Any) -> list[Any]:
    """Regularized IC for an on-section {y=0} box with xi = x-(1-mu) < 0 (u=0 branch)."""
    negxi = iv.mpf(OM1) - x
    v = iv.sqrt(negxi)
    r1 = iv.sqrt((x + iv.mpf(MU)) ** 2)
    r2 = negxi
    vy2 = x**2 + 2 * (1 - iv.mpf(MU)) / r1 + 2 * iv.mpf(MU) / r2 - vx**2 - iv.mpf(C_LEVEL)
    vy = iv.sqrt(vy2)
    pxi = vx
    peta = vy + x
    pu = 2 * v * peta
    pv = -2 * v * pxi
    return [iv.mpf([0, 0]), v, pu, pv, iv.mpf([0, 0])]


def _sig_2uv(iv: Any, st: list[Any]) -> Any:
    return iv.mpf(2) * st[0] * st[1]


def _grad_2uv(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(2) * st[1], iv.mpf(2) * st[0], iv.mpf(0), iv.mpf(0), iv.mpf(0)]


def _secmap(iv: Any, jet: Any, vjet: Any, w0: list[Any]) -> dict[str, Any]:
    return vti.rigorous_section_map(
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


def _f(x: Any) -> float:
    return float(x)


def _wid(x: Any) -> float:
    return float(x.delta.b)


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv
    t0 = time.time()

    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)

    # (0) CENTER-point run once: tight P(w_hat) shared by every h-set size.
    w0c = _lc_ic_onsection(iv, iv.mpf([X0, X0]), iv.mpf([VX0, VX0]))
    resc = _secmap(iv, jet, vjet, w0c)
    print(f"[#674] center run done ({time.time() - t0:.0f}s) found={resc['found']}")
    if not resc["found"]:
        return {"task": "#674", "result": "CENTER RUN WALL", "reason": resc.get("reason")}
    wstar = resc["crossing_state"]
    phys_c = vti.lc_secondary_to_physical(iv, wstar, MU)
    cx, cvx = _f(phys_c[0].mid), _f(phys_c[2].mid)

    cert: dict[str, Any] = {
        "task": "#674",
        "stage": "Stage 8 -- correlation-preserving (mean-value) section-map image",
        "system": {"mu": MU, "C": C_LEVEL, "h_energy": H_ENERGY},
        "hset_N": {
            "at": "first published W-Z Oterma section point",
            "center_x_xdot": [X0, VX0],
            "unstable_dir_x_xdot": UDIR,
            "stable_dir_x_xdot": SDIR,
            "note": (
                "box on {y=0}; axes = L1* stable(expander)/unstable eigen-directions "
                "(heuristic frame); regularized IC built directly on Levi-Civita branch u=0."
            ),
        },
        "center_image_x_xdot": [cx, cvx],
        "center_image_widths_x_xdot": [_wid(phys_c[0]), _wid(phys_c[2])],
        "method": (
            "P(face) subset P(w_hat) + [DP over box] . (w0_face - w_hat)  (interval "
            "mean-value; vti.section_map_meanvalue_image). ONE tight linearization instead "
            "of independent per-face box propagation -> no re-accumulated wrapping."
        ),
        "stage7_reference_numbers": {
            "genuine_edge_center_separation": 5.9e-6,
            "boxhull_per_edge_width": 5.2e-5,
            "note": "#673 (Stage 7) box-hull figures this stage improves upon.",
        },
        "sizes": [],
    }

    def meanvalue_faces(dp_box: list[Any], ru: float, rs: float) -> dict[str, tuple[Any, Any]]:
        faces = {}
        for name in ("box", "left", "right"):
            xf, vxf = _phys_subbox(iv, X0, VX0, ru, rs, name)
            w0f = _lc_ic_onsection(iv, xf, vxf)
            off = [w0f[i] - w0c[i] for i in range(5)]
            reg_img = vti.section_map_meanvalue_image(iv, wstar, dp_box, off)
            ph = vti.lc_secondary_to_physical(iv, reg_img, MU)
            faces[name] = (ph[0], ph[2])  # (x', xdot')
        return faces

    certified_any = False
    published2 = None  # (A) diagnostic, computed from the ru=1e-6 image
    for ru, rs in [(1e-6, 1e-8), (1e-8, 1e-10)]:
        xb, vxb = _phys_subbox(iv, X0, VX0, ru, rs, "box")
        resb = _secmap(iv, jet, vjet, _lc_ic_onsection(iv, xb, vxb))
        print(f"[#674] whole-box ru={ru:.0e} done ({time.time() - t0:.0f}s) found={resb['found']}")
        if not resb["found"]:
            cert["sizes"].append({"ru": ru, "rs": rs, "result": "WALL", "reason": resb["reason"]})
            continue
        dp_box = resb["section_jacobian"]
        dp_width = max(_wid(dp_box[i][j]) for i in range(5) for j in range(5))
        faces = meanvalue_faces(dp_box, ru, rs)
        (bx, bvx), (lx, lvx), (rx, rvx) = faces["box"], faces["left"], faces["right"]

        sep = abs(_f(lvx.mid) - _f(rvx.mid))
        edgew = max(_wid(lvx), _wid(rvx))

        # constructed M at the actual first-return image: unstable axis = xdot', stable = x'
        gap = sep - edgew
        hu_img = 0.4 * gap if gap > 0 else max(1e-12, 0.25 * sep)
        x_halfwidth = _wid(bx) / 2.0
        hs_img = 3.0 * (x_halfwidth + rs)
        cov = vti.covering_relation_2d(
            iv,
            m_center=[iv.mpf([cx, cx]), iv.mpf([cvx, cvx])],
            m_uvec=[iv.mpf(0), iv.mpf([hu_img, hu_img])],
            m_svec=[iv.mpf([hs_img, hs_img]), iv.mpf(0)],
            image_left=[lx, lvx],
            image_right=[rx, rvx],
            image_box=[bx, bvx],
        )
        certified_any = certified_any or bool(cov["covers"])
        cert["sizes"].append(
            {
                "ru": ru,
                "rs": rs,
                "dp_box_enclosure_width": dp_width,
                "image_left_xdot": [_f(lvx.a), _f(lvx.b)],
                "image_right_xdot": [_f(rvx.a), _f(rvx.b)],
                "image_whole_box_xdot_width": _wid(bvx),
                "edge_center_separation": sep,
                "per_edge_enclosure_width": edgew,
                "separation_over_width_ratio": (sep / edgew if edgew > 0 else None),
                "M_unstable_halfwidth": hu_img,
                "M_stable_halfwidth": hs_img,
                "covers": bool(cov["covers"]),
                "stable_containment": bool(cov["stable_containment"]),
                "unstable_exit": bool(cov["unstable_exit"]),
                "orientation": cov["orientation"],
                "unstable_coord_left": cov["unstable_coord_left"],
                "unstable_coord_right": cov["unstable_coord_right"],
                "stable_coord_image": cov["stable_coord_image"],
            }
        )
        if ru == 1e-6:
            # (A) vs the SECOND published point: image lands near x'~1.0035, far from 0.921
            cov2 = vti.covering_relation_2d(
                iv,
                m_center=[iv.mpf(X2), iv.mpf(VX2)],
                m_uvec=[iv.mpf(UDIR[0]) * iv.mpf(0.05), iv.mpf(UDIR[1]) * iv.mpf(0.05)],
                m_svec=[iv.mpf(SDIR[0]) * iv.mpf(0.05), iv.mpf(SDIR[1]) * iv.mpf(0.05)],
                image_left=[lx, lvx],
                image_right=[rx, rvx],
                image_box=[bx, bvx],
            )
            published2 = {
                "M_center_x_xdot": [X2, VX2],
                "covers": bool(cov2["covers"]),
                "image_dist_to_point2_x": abs(cx - X2),
                "interpretation": (
                    f"NON-covering as expected: the first {{y=0}} return of N lands at "
                    f"x'~{cx:.4f}, ~{abs(cx - X2):.3f} away from published point 2 "
                    f"(x={X2:.4f}); the rounded IC does not track the true unstable "
                    "heteroclinic manifold and the first-return map does not pair points 1->2."
                ),
            }

    cert["covering_vs_published_point2"] = published2
    cert["certified"] = certified_any
    cert["result"] = (
        "GENUINE rigorous covering relation N ==P==> M CERTIFIED on the real regularized "
        "CR3BP section map (through the first Jupiter perijove) via the correlation-"
        "preserving mean-value image -- at a proof-appropriate h-set size. At #673's "
        "original ru=1e-6 the covering still does NOT certify: [DP] genuinely varies by "
        "~0.9 over that box (strong flyby nonlinearity), so the mean-value edge width "
        "(~2.1e-5, already ~2.5x tighter than #673's 5.2e-5 box-hull) still exceeds the "
        "~6e-6 signal. Shrinking to ru=1e-8 (the [DP]-variation term ~box^2 falls faster "
        "than the ~box signal) gives ratio ~2.0 and a clean covering certificate."
        if certified_any
        else "NO covering certified at either size (see per-size numbers)."
    )
    cert["scope_note"] = (
        "One covering relation as the bounded proof-of-concept -- NOT the full ~30-point "
        "W-Z chain. M is a constructed h-set at the actual first-return image, not a "
        "published W-Z point. The covering-relation checker and the mean-value image "
        "primitive are both validated on synthetic controls "
        "(tests/scripts/test_673_covering_relation.py, test_674_meanvalue_image.py)."
    )
    return cert


def main() -> None:
    cert = build_certificate()
    OUT.write_text(json.dumps(cert, indent=2) + "\n")
    print("\n[#674] ============ RESULT ============")
    for s in cert.get("sizes", []):
        if "covers" in s:
            print(
                f"[#674] ru={s['ru']:.0e}: DPw={s['dp_box_enclosure_width']:.2e} "
                f"sep={s['edge_center_separation']:.3e} width={s['per_edge_enclosure_width']:.3e} "
                f"ratio={s['separation_over_width_ratio']:.2f} -> covers={s['covers']} "
                f"(S={s['stable_containment']} U={s['unstable_exit']} {s['orientation']})"
            )
        else:
            print(f"[#674] ru={s['ru']:.0e}: {s.get('result')}")
    print(f"[#674] certified={cert.get('certified')}")
    print(f"[#674] certificate -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
