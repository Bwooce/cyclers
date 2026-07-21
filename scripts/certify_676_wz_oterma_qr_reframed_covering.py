"""#676 -- inter-return QR-reframed composed covering at the fixed ru=1e-6 h-set.

Stage 10 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668``-``#675``).

#675's composed covering COLLAPSED (ratio 0.284 -> 0.027 at N=2) because the composition
carries the reachable set from return to return as an axis-aligned interval BOX: the
leg-(k+1) IC box is built component-wise (``section_map_meanvalue_image``), the axis-aligned
hull of a thin, sheared true image, and ``rigorous_section_map`` re-inflates from that hull
-- the SAME wrapping effect #669's within-flow QR reframing defeats, reappearing at the
DISCRETE return-to-return boundary.

This applies the discrete analogue of #669's Lohner-QR reframing at the return boundary
(``vti.reframe_image_qr``): the leg's reachable set is carried as a parallelepiped
``center + A @ [r]`` with ``A`` a POINT orthogonal frame re-chosen each return from a QR of
the single-leg section Jacobian ``[DP_k]``, so ``[r]`` stays tight in coordinates aligned to
the current stretch and the box hull fed to the next leg's ``rigorous_section_map`` does not
accumulate the per-return wrapping.  The composed Jacobian (for the covering ratio) is kept
as the raw chain-rule product ``compose_section_jacobians`` -- separately QR-reframing the
ACCUMULATED Jacobian was shown (test_676) to add its own rebasing wrapping without net
benefit at these low return counts, so only the IMAGE is reframed.

To make N=3 reachable (the regime where the reframing benefit compounds), the section-search
window is widened to tau_max=40 with n_steps=320 -- the SAME step size h=0.125 as #675, so
the N=1 and N=2 crossings are found at the IDENTICAL fictitious times (16.57, 7.97) and are
directly comparable; the 3rd return simply occurs at tau~24 (past #675's tau_max=20 window,
which is exactly why #675 walled at N=3).  ru=1e-6 is NOT shrunk.

The driver runs BOTH the un-reframed box chain (reproducing #675) and the reframed chain in
the same run and reports the covering ratios side by side, so the effect of reframing is a
controlled comparison.  The reframing MECHANISM is validated on closed-form-ground-truth
synthetic controls in ``tests/scripts/test_676_qr_reframed_covering.py``.

Run:  ``uv run python scripts/certify_676_wz_oterma_qr_reframed_covering.py``
Writes ``data/676_wz_oterma_qr_reframed_covering_certificate.json`` (checkpointed per leg).
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

OUT = ROOT / "data" / "676_wz_oterma_qr_reframed_covering_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2
OM1 = 1.0 - MU

# First published W-Z Oterma section point (on {y=0}) -- identical to #675.
X0, VX0 = 0.9522928423486199945, 1.23e-05

# FIXED, proof-meaningful h-set size (#674/#675's original; NOT shrunk to force a result).
RU, RS = 1e-6, 1e-8
N_RETURNS = 3

# Wider section-search window than #675 (tau_max 20 -> 40) at the SAME step size h=0.125
# (n_steps 160 -> 320), so N=1/N=2 crossings are identical and N=3 (at tau~24) is reachable.
TAU_MAX, N_STEPS = 40.0, 320

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
        tau_max=TAU_MAX,
        n_steps=N_STEPS,
        order=10,
        tau_min=1.0,
    )


def _f(x: Any) -> float:
    return float(x)


def _wid(x: Any) -> float:
    return float(x.delta.b)


def _checkpoint(cert: dict[str, Any]) -> None:
    OUT.write_text(json.dumps(cert, indent=2) + "\n")


def _ratio_from_jac(
    iv: Any, wk_hat: list[Any], jac_comp: Any, face_offsets: dict[str, list[Any]]
) -> dict[str, Any]:
    """Composed covering ratio + verdict from a composed Jacobian (SAME as #675's #674 build)."""
    faces_phys: dict[str, tuple[Any, Any]] = {}
    for name in ("box", "left", "right"):
        reg_img = vti.section_map_meanvalue_image(iv, wk_hat, jac_comp, face_offsets[name])
        ph = vti.lc_secondary_to_physical(iv, reg_img, MU)
        faces_phys[name] = (ph[0], ph[2])
    (bx, bvx), (lx, lvx), (rx, rvx) = faces_phys["box"], faces_phys["left"], faces_phys["right"]
    phys_c = vti.lc_secondary_to_physical(iv, wk_hat, MU)
    cx, cvx = _f(phys_c[0].mid), _f(phys_c[2].mid)
    sep = abs(_f(lvx.mid) - _f(rvx.mid))
    edgew = max(_wid(lvx), _wid(rvx))
    gap = sep - edgew
    hu_img = 0.4 * gap if gap > 0 else max(1e-12, 0.25 * sep)
    hs_img = 3.0 * (_wid(bx) / 2.0 + RS)
    cov = vti.covering_relation_2d(
        iv,
        m_center=[iv.mpf([cx, cx]), iv.mpf([cvx, cvx])],
        m_uvec=[iv.mpf(0), iv.mpf([hu_img, hu_img])],
        m_svec=[iv.mpf([hs_img, hs_img]), iv.mpf(0)],
        image_left=[lx, lvx],
        image_right=[rx, rvx],
        image_box=[bx, bvx],
    )
    return {
        "center_image_x_xdot": [cx, cvx],
        "composed_jacobian_enclosure_width": max(
            _wid(jac_comp[i][j]) for i in range(5) for j in range(5)
        ),
        "edge_center_separation": sep,
        "per_edge_enclosure_width": edgew,
        "separation_over_width_ratio": (sep / edgew if edgew > 0 else None),
        "covers": bool(cov["covers"]),
        "stable_containment": bool(cov["stable_containment"]),
        "unstable_exit": bool(cov["unstable_exit"]),
        "orientation": cov["orientation"],
    }


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv
    t0 = time.time()

    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)

    w0c = _lc_ic_onsection(iv, iv.mpf([X0, X0]), iv.mpf([VX0, VX0]))
    face_offsets: dict[str, list[Any]] = {}
    for name in ("box", "left", "right"):
        xf, vxf = _phys_subbox(iv, X0, VX0, RU, RS, name)
        w0f = _lc_ic_onsection(iv, xf, vxf)
        face_offsets[name] = [w0f[i] - w0c[i] for i in range(5)]

    box_ic0 = _lc_ic_onsection(iv, *_phys_subbox(iv, X0, VX0, RU, RS, "box"))

    cert: dict[str, Any] = {
        "task": "#676",
        "stage": "Stage 10 -- inter-return QR-reframed composed covering at fixed ru=1e-6",
        "system": {"mu": MU, "C": C_LEVEL, "h_energy": H_ENERGY},
        "hset_N": {
            "at": "first published W-Z Oterma section point",
            "center_x_xdot": [X0, VX0],
            "ru": RU,
            "rs": RS,
            "unstable_dir_x_xdot": UDIR,
            "stable_dir_x_xdot": SDIR,
            "note": "FIXED size, NOT shrunk (identical to #674/#675).",
        },
        "search_window": {
            "tau_max": TAU_MAX,
            "n_steps": N_STEPS,
            "step_h": TAU_MAX / N_STEPS,
            "note": (
                "Wider than #675 (20/160) at the SAME h=0.125, so N=1/N=2 crossings are "
                "identical and N=3 (at tau~24, past #675's window) is reachable."
            ),
        },
        "method": (
            "Un-reframed BOX chain (reproduces #675: leg IC box = axis-aligned mean-value "
            "hull) vs REFRAMED chain (vti.reframe_image_qr: leg IC set carried as a "
            "parallelepiped center+A@[r], A a point orthogonal frame from QR of the single-leg "
            "[DP_k], box hull fed to the next leg). Composed Jacobian is the raw chain-rule "
            "product in both. Both run in one pass for a controlled comparison."
        ),
        "reframing_control": "tests/scripts/test_676_qr_reframed_covering.py",
        "returns": [],
    }
    _checkpoint(cert)

    centre_state = w0c
    # baseline (un-reframed) box chain state
    base_box_ic = box_ic0
    base_box_center = w0c
    base_dp_boxes: list[Any] = []
    base_wall = None
    # reframed box chain state (parallelepiped center + frame @ [rbox])
    rf_center = [vti._mid(w0c[i]) for i in range(5)]
    rf_frame: list[list[Any]] = [
        [iv.mpf(1) if i == j else iv.mpf(0) for j in range(5)] for i in range(5)
    ]
    rf_rbox = [box_ic0[i] - rf_center[i] for i in range(5)]
    rf_dp_boxes: list[Any] = []
    rf_wall = None

    for k in range(1, N_RETURNS + 1):
        rc = _secmap(iv, jet, vjet, centre_state)
        print(f"[#676] return {k}: centre done ({time.time() - t0:.0f}s) found={rc['found']}")
        if not rc["found"]:
            cert["returns"].append({"n": k, "result": "CENTRE WALL", "reason": rc.get("reason")})
            _checkpoint(cert)
            break
        centre_state = rc["crossing_state"]
        wk_hat = centre_state

        row: dict[str, Any] = {"n": k, "tau_star_this_leg": _f(rc["tau_star"].mid)}

        # --- baseline (un-reframed) leg ---
        if base_wall is None:
            rb = _secmap(iv, jet, vjet, base_box_ic)
            print(f"[#676] return {k}: baseline box done ({time.time() - t0:.0f}s) f={rb['found']}")
            if not rb["found"]:
                base_wall = rb.get("reason")
                row["baseline"] = {"result": "BOX WALL", "reason": base_wall}
            else:
                dp_k = rb["section_jacobian"]
                base_dp_boxes.append(dp_k)
                box_off = [base_box_ic[i] - base_box_center[i] for i in range(5)]
                base_box_ic = vti.section_map_meanvalue_image(iv, wk_hat, dp_k, box_off)
                base_box_center = wk_hat
                jac_comp = vti.compose_section_jacobians(iv, base_dp_boxes)
                res = _ratio_from_jac(iv, wk_hat, jac_comp, face_offsets)
                res["leg_jacobian_enclosure_width"] = max(
                    _wid(dp_k[i][j]) for i in range(5) for j in range(5)
                )
                row["baseline"] = res
        else:
            row["baseline"] = {"result": "BOX WALL (upstream)", "reason": base_wall}

        # --- reframed leg ---
        if rf_wall is None:
            rf_box_ic = vti.enclosure_box(iv, rf_center, rf_frame, rf_rbox)
            rr = _secmap(iv, jet, vjet, rf_box_ic)
            print(f"[#676] return {k}: reframed box done ({time.time() - t0:.0f}s) f={rr['found']}")
            if not rr["found"]:
                rf_wall = rr.get("reason")
                row["reframed"] = {"result": "BOX WALL", "reason": rf_wall}
            else:
                dp_k = rr["section_jacobian"]
                rf_dp_boxes.append(dp_k)
                reb = vti.reframe_image_qr(iv, dp_k, rf_frame, rf_rbox, wk_hat)
                if reb is None:
                    rf_wall = "reframe_image_qr: frame not rigorously invertible"
                    row["reframed"] = {"result": "REFRAME FAIL", "reason": rf_wall}
                else:
                    rf_center, rf_frame, rf_rbox = reb
                    jac_comp = vti.compose_section_jacobians(iv, rf_dp_boxes)
                    res = _ratio_from_jac(iv, wk_hat, jac_comp, face_offsets)
                    res["leg_jacobian_enclosure_width"] = max(
                        _wid(dp_k[i][j]) for i in range(5) for j in range(5)
                    )
                    res["reframed_box_ic_max_halfwidth"] = max(_wid(c) for c in rf_box_ic) / 2.0
                    row["reframed"] = res
        else:
            row["reframed"] = {"result": "BOX WALL (upstream)", "reason": rf_wall}

        cert["returns"].append(row)
        _checkpoint(cert)
        b = row["baseline"]
        r = row["reframed"]
        br = b.get("separation_over_width_ratio")
        rr_ratio = r.get("separation_over_width_ratio")
        print(
            f"[#676] return {k}: baseline ratio="
            f"{br:.4f} covers={b.get('covers')} | reframed ratio="
            f"{(rr_ratio if rr_ratio is not None else float('nan')):.4f} covers={r.get('covers')}"
            if br is not None
            else f"[#676] return {k}: baseline={b.get('result')} reframed={r.get('result')}"
        )

    _finalize(cert)
    _checkpoint(cert)
    return cert


def _finalize(cert: dict[str, Any]) -> None:
    def ratios(kind: str) -> list[tuple[int, float | None]]:
        out = []
        for row in cert["returns"]:
            if "n" in row and kind in row and "separation_over_width_ratio" in row[kind]:
                out.append((row["n"], row[kind]["separation_over_width_ratio"]))
        return out

    base = ratios("baseline")
    ref = ratios("reframed")
    cert["baseline_ratio_scaling"] = base
    cert["reframed_ratio_scaling"] = ref
    cert["baseline_certified"] = [
        row["n"] for row in cert["returns"] if row.get("baseline", {}).get("covers")
    ]
    cert["reframed_certified"] = [
        row["n"] for row in cert["returns"] if row.get("reframed", {}).get("covers")
    ]

    # verdict on whether reframing arrested #675's collapse
    verdict = "no multi-return comparison available"
    b2 = {n: r for n, r in base}.get(2)
    r2 = {n: r for n, r in ref}.get(2)
    b3 = {n: r for n, r in base}.get(3)
    r3 = {n: r for n, r in ref}.get(3)
    lines = []
    if b2 is not None and r2 is not None:
        rel = (r2 - b2) / b2 if b2 else 0.0
        lines.append(
            f"N=2 baseline={b2:.4f} reframed={r2:.4f} ({rel:+.1%}); "
            + ("reframing helps" if r2 > b2 * 1.02 else "no material change")
        )
    if b3 is not None and r3 is not None:
        rel = (r3 - b3) / b3 if b3 else 0.0
        lines.append(
            f"N=3 baseline={b3:.4f} reframed={r3:.4f} ({rel:+.1%}); "
            + ("reframing helps" if r3 > b3 * 1.02 else "no material change")
        )
    elif b3 is not None or r3 is not None:
        lines.append(f"N=3 baseline={b3} reframed={r3} (one chain walled at N=3)")
    if lines:
        verdict = " | ".join(lines)
    cert["verdict"] = verdict
    cert["result"] = (
        f"ru={RU:g} (fixed, NOT shrunk). "
        f"baseline ratio/return: {', '.join(f'N={n}:{r:.3f}' for n, r in base if r is not None)}. "
        f"reframed ratio/return: {', '.join(f'N={n}:{r:.3f}' for n, r in ref if r is not None)}. "
        f"{verdict}."
    )
    cert["scope_note"] = (
        "Composed 2-3 returns at FIXED ru=1e-6 (NOT shrunk). Inter-return QR reframing of the "
        "image (vti.reframe_image_qr); composed Jacobian kept as the raw chain-rule product. "
        "Wider tau window (40/320) at the SAME h=0.125 as #675 so N=1/N=2 are identical and N=3 "
        "is reachable. M is a constructed h-set at the k-th-return image, not a published W-Z "
        "point. Bounded proof-of-concept, NOT the full ~30-section W-Z chain. Reframing "
        "mechanism validated on closed-form controls (test_676_qr_reframed_covering.py)."
    )


def main() -> None:
    cert = build_certificate()
    print("\n[#676] ============ RESULT ============")
    for row in cert.get("returns", []):
        if "result" in row and "n" in row and "baseline" not in row:
            print(f"[#676] N={row['n']}: {row['result']} -- {row.get('reason')}")
            continue
        b = row.get("baseline", {})
        r = row.get("reframed", {})

        def fmt(d: dict[str, Any]) -> str:
            if "separation_over_width_ratio" in d:
                return (
                    f"ratio={d['separation_over_width_ratio']:.4f} covers={d['covers']} "
                    f"(jacW={d['composed_jacobian_enclosure_width']:.2e} "
                    f"legW={d.get('leg_jacobian_enclosure_width', float('nan')):.2e})"
                )
            return f"{d.get('result')}"

        print(f"[#676] N={row['n']}: baseline {fmt(b)}")
        print(f"[#676]        reframed {fmt(r)}")
    print(
        f"[#676] baseline_certified={cert.get('baseline_certified')} "
        f"reframed_certified={cert.get('reframed_certified')}"
    )
    print(f"[#676] {cert.get('result')}")
    print(f"[#676] certificate -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
