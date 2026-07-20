"""#675 -- composed (multi-return) covering relation at a fixed, proof-meaningful h-set.

Stage 9 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668``-``#674``).

#674 certified a SINGLE-return covering relation on the real regularized CR3BP section
map, but only at a near-vacuous h-set size (ru=1e-8); at the more meaningful ru=1e-6 the
single-return covering fails (separation/width ratio 0.28) because the section-map
Jacobian ``[DP]`` genuinely varies ~0.9 across a 1e-6 box.  The genuine W-Z fix is to
COMPOSE several consecutive section returns at the FIXED ru=1e-6 size: a hyperbolic map's
stretch compounds multiplicatively (``P^N`` stretch ``~ lambda^N``), so if the per-leg
enclosures stay tight the compounded signal can dominate the compounded over-approximation
without shrinking the h-set.

Composition (all rigorous, via ``vti``):
  * CENTRE chain: run the section map on the box CENTRE to the 1st return (tight
    ``w1_hat = P(w_hat)``), then from ``w1_hat`` to the 2nd (``w2_hat``), then the 3rd --
    each ``wk_hat = P^k(w_hat)`` a tight enclosure.
  * BOX chain: run the section map on the whole h-set box ``W0`` (giving ``DP_1`` over
    ``W0`` and, via the mean-value image, an enclosure ``B1`` of ``P(W0)``); run it again
    on ``B1`` (giving ``DP_2`` over ``B1 superset P(W0)`` and an enclosure ``B2`` of
    ``P^2(W0)``); once more for ``DP_3``.
  * Composed image of a face of the ORIGINAL h-set: by the chain rule + interval
    mean-value theorem,
        ``P^N(face) subset P^N(w_hat) + (DP_N ... DP_1) . (w0_face - w_hat)``,
    with ``DP_N ... DP_1 = vti.compose_section_jacobians([DP_1, ..., DP_N])`` a rigorous
    enclosure of ``D(P^N)`` over ``W0`` (each factor ``DP_k`` encloses ``DP`` over an
    enclosure of ``P^{k-1}(W0)``).

Reports, at the FIXED ru=1e-6 size, how the genuine edge separation (signal) and the
mean-value edge width (over-approximation) each scale across N=1,2,3, and whether the
covering certifies at any N -- WITHOUT shrinking the h-set.  The composition MACHINERY
itself is validated on closed-form-stretch synthetic controls in
``tests/scripts/test_675_composed_meanvalue.py``.

Run:  ``uv run python scripts/certify_675_wz_oterma_composed_covering.py``
Writes ``data/675_wz_oterma_composed_covering_certificate.json`` (checkpointed per leg).
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

OUT = ROOT / "data" / "675_wz_oterma_composed_covering_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2
OM1 = 1.0 - MU

# First published W-Z Oterma section point (on {y=0}).
X0, VX0 = 0.9522928423486199945, 1.23e-05

# FIXED, proof-meaningful h-set size (#674's original; NOT shrunk to force a result).
RU, RS = 1e-6, 1e-8
N_RETURNS = 3

# L1* stable/unstable eigen-directions in (x, xdot) (heuristic frame; see #673/#674).
_E_EXPAND = [-1.0, 2.5733011]
_E_CONTRACT = [1.0, 2.5733011]


def _norm(v: list[float]) -> list[float]:
    n = math.hypot(v[0], v[1])
    return [v[0] / n, v[1] / n]


UDIR = _norm(_E_EXPAND)
SDIR = _norm(_E_CONTRACT)


def _phys_subbox(iv: Any, cx: float, cvx: float, ru: float, rs: float, a: str) -> tuple[Any, Any]:
    aval = {"box": iv.mpf([-1, 1]), "left": iv.mpf([-1, -1]), "right": iv.mpf([1, 1])}[a]
    b = iv.mpf([-1, 1])
    x = iv.mpf(cx) + aval * iv.mpf(ru * UDIR[0]) + b * iv.mpf(rs * SDIR[0])
    vx = iv.mpf(cvx) + aval * iv.mpf(ru * UDIR[1]) + b * iv.mpf(rs * SDIR[1])
    return x, vx


def _lc_ic_onsection(iv: Any, x: Any, vx: Any) -> list[Any]:
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
        tau_max=20.0,
        n_steps=160,
        order=10,
        tau_min=1.0,
    )


def _f(x: Any) -> float:
    return float(x)


def _wid(x: Any) -> float:
    return float(x.delta.b)


def _checkpoint(cert: dict[str, Any]) -> None:
    OUT.write_text(json.dumps(cert, indent=2) + "\n")


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv
    t0 = time.time()

    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)

    # Regularized IC of the ORIGINAL box centre and its three faces (fixed offsets).
    w0c = _lc_ic_onsection(iv, iv.mpf([X0, X0]), iv.mpf([VX0, VX0]))
    face_offsets: dict[str, list[Any]] = {}
    for name in ("box", "left", "right"):
        xf, vxf = _phys_subbox(iv, X0, VX0, RU, RS, name)
        w0f = _lc_ic_onsection(iv, xf, vxf)
        face_offsets[name] = [w0f[i] - w0c[i] for i in range(5)]

    cert: dict[str, Any] = {
        "task": "#675",
        "stage": "Stage 9 -- composed (multi-return) covering at fixed ru=1e-6",
        "system": {"mu": MU, "C": C_LEVEL, "h_energy": H_ENERGY},
        "hset_N": {
            "at": "first published W-Z Oterma section point",
            "center_x_xdot": [X0, VX0],
            "ru": RU,
            "rs": RS,
            "unstable_dir_x_xdot": UDIR,
            "stable_dir_x_xdot": SDIR,
            "note": "FIXED size, NOT shrunk to force a result (the whole point of Stage 9).",
        },
        "method": (
            "Composed mean-value image: P^N(face) subset P^N(w_hat) + "
            "(DP_N...DP_1).(w0_face - w_hat), DP_N...DP_1 = "
            "vti.compose_section_jacobians([DP_1..DP_N]) (chain rule, each DP_k over an "
            "enclosure of P^{k-1}(W0)); centre chain gives the tight P^N(w_hat)."
        ),
        "composition_control": "tests/scripts/test_675_composed_meanvalue.py",
        "returns": [],
    }
    _checkpoint(cert)

    # Chains.  centre_state = P^k(w_hat) tight; box_ic = enclosure of P^k(W0);
    # box_center = the tight centre about which box_ic is a mean-value image.
    centre_state = w0c
    box_ic = _lc_ic_onsection(iv, *_phys_subbox(iv, X0, VX0, RU, RS, "box"))
    box_center = w0c
    dp_boxes: list[Any] = []

    for k in range(1, N_RETURNS + 1):
        # --- CENTRE leg: P^{k} of the centre ---
        rc = _secmap(iv, jet, vjet, centre_state)
        print(f"[#675] return {k}: centre run done ({time.time() - t0:.0f}s) found={rc['found']}")
        if not rc["found"]:
            cert["returns"].append({"n": k, "result": "CENTRE WALL", "reason": rc.get("reason")})
            _checkpoint(cert)
            break
        centre_state = rc["crossing_state"]
        wk_hat = centre_state

        # --- BOX leg: DP_k over box_ic, and B_k = enclosure of P^k(W0) ---
        rb = _secmap(iv, jet, vjet, box_ic)
        print(f"[#675] return {k}: box run done ({time.time() - t0:.0f}s) found={rb['found']}")
        if not rb["found"]:
            cert["returns"].append({"n": k, "result": "BOX WALL", "reason": rb.get("reason")})
            _checkpoint(cert)
            break
        dp_k = rb["section_jacobian"]
        dp_boxes.append(dp_k)
        # B_k = mean-value image of box_ic under P about box_center: P(box_ic) subset
        #       P(box_center) + DP_k . (box_ic - box_center) = wk_hat + DP_k.(box_ic-box_center)
        box_offset_k = [box_ic[i] - box_center[i] for i in range(5)]
        b_next = vti.section_map_meanvalue_image(iv, wk_hat, dp_k, box_offset_k)
        box_ic = b_next
        box_center = wk_hat

        # --- composed image of the ORIGINAL faces at return k ---
        jac_comp = vti.compose_section_jacobians(iv, dp_boxes)  # DP_k ... DP_1 over W0
        dp_comp_w = max(_wid(jac_comp[i][j]) for i in range(5) for j in range(5))
        faces_phys: dict[str, tuple[Any, Any]] = {}
        for name in ("box", "left", "right"):
            reg_img = vti.section_map_meanvalue_image(iv, wk_hat, jac_comp, face_offsets[name])
            ph = vti.lc_secondary_to_physical(iv, reg_img, MU)
            faces_phys[name] = (ph[0], ph[2])  # (x', xdot')
        (bx, bvx), (lx, lvx), (rx, rvx) = (
            faces_phys["box"],
            faces_phys["left"],
            faces_phys["right"],
        )
        phys_c = vti.lc_secondary_to_physical(iv, wk_hat, MU)
        cx, cvx = _f(phys_c[0].mid), _f(phys_c[2].mid)

        sep = abs(_f(lvx.mid) - _f(rvx.mid))
        edgew = max(_wid(lvx), _wid(rvx))

        # constructed M at the actual k-th-return image (same construction as #674):
        # unstable axis = xdot', stable axis = x'.
        gap = sep - edgew
        hu_img = 0.4 * gap if gap > 0 else max(1e-12, 0.25 * sep)
        x_halfwidth = _wid(bx) / 2.0
        hs_img = 3.0 * (x_halfwidth + RS)
        cov = vti.covering_relation_2d(
            iv,
            m_center=[iv.mpf([cx, cx]), iv.mpf([cvx, cvx])],
            m_uvec=[iv.mpf(0), iv.mpf([hu_img, hu_img])],
            m_svec=[iv.mpf([hs_img, hs_img]), iv.mpf(0)],
            image_left=[lx, lvx],
            image_right=[rx, rvx],
            image_box=[bx, bvx],
        )
        cert["returns"].append(
            {
                "n": k,
                "tau_star_this_leg": _f(rc["tau_star"].mid),
                "center_image_x_xdot": [cx, cvx],
                "center_image_widths_x_xdot": [_wid(phys_c[0]), _wid(phys_c[2])],
                "composed_jacobian_enclosure_width": dp_comp_w,
                "leg_jacobian_enclosure_width": max(
                    _wid(dp_k[i][j]) for i in range(5) for j in range(5)
                ),
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
        _checkpoint(cert)
        print(
            f"[#675] return {k}: sep={sep:.3e} width={edgew:.3e} "
            f"ratio={(sep / edgew if edgew > 0 else float('nan')):.3f} covers={cov['covers']}"
        )

    certified = [r["n"] for r in cert["returns"] if r.get("covers")]
    cert["certified_returns"] = certified
    cert["certified"] = bool(certified)
    ratios = [
        (r["n"], r.get("separation_over_width_ratio"))
        for r in cert["returns"]
        if "separation_over_width_ratio" in r
    ]
    cert["ratio_scaling"] = ratios
    cert["result"] = _summarize(cert, ratios, certified)
    cert["scope_note"] = (
        "Composed 2-3 returns at FIXED ru=1e-6 (NOT shrunk). M is a constructed h-set at "
        "the actual k-th-return image, not a published W-Z point. Bounded proof-of-concept "
        "for the compounding mechanism, NOT the full ~30-section W-Z chain. Composition "
        "machinery validated on closed-form controls (test_675_composed_meanvalue.py)."
    )
    _checkpoint(cert)
    return cert


def _summarize(
    cert: dict[str, Any], ratios: list[tuple[int, float | None]], certified: list[int]
) -> str:
    walls = [r for r in cert["returns"] if "result" in r]
    parts = [f"ru={RU:g} (fixed). ratio/return: " + ", ".join(f"N={n}:{r:.3f}" for n, r in ratios)]
    if walls:
        parts.append(f"chain stopped: {walls[0]['result']} at N={walls[0]['n']}")
    if certified:
        parts.append(f"COVERING CERTIFIED at composed return(s) {certified} without shrinking ru.")
    else:
        parts.append(
            "NO composed covering certified at ru=1e-6 within the attempted returns: the "
            "over-approximation compounds at least as fast as the hyperbolic signal (see "
            "ratio scaling)."
        )
    return " | ".join(parts)


def main() -> None:
    cert = build_certificate()
    print("\n[#675] ============ RESULT ============")
    for r in cert.get("returns", []):
        if "covers" in r:
            print(
                f"[#675] N={r['n']}: DPcomp_w={r['composed_jacobian_enclosure_width']:.2e} "
                f"sep={r['edge_center_separation']:.3e} width={r['per_edge_enclosure_width']:.3e} "
                f"ratio={r['separation_over_width_ratio']:.3f} -> covers={r['covers']} "
                f"(S={r['stable_containment']} U={r['unstable_exit']} {r['orientation']})"
            )
        else:
            print(f"[#675] N={r['n']}: {r.get('result')} -- {r.get('reason')}")
    print(f"[#675] certified_returns={cert.get('certified_returns')}")
    print(f"[#675] {cert.get('result')}")
    print(f"[#675] certificate -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
