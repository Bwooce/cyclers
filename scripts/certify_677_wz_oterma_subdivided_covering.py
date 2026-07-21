"""#677 -- finer SUB-BOXING (subdivide-before-propagate) composed covering at fixed ru=1e-6.

Stage 11 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668``-``#676``).

#676 pinned the cause of #675's composed-covering collapse: NOT inter-return box-hull
decorrelation (QR reframing removed that mechanism cleanly on synthetic controls but did
not help on Oterma) -- it is intra-leg Jacobian VARIATION over a leg IC box that has already
inflated to ~35x the original h-set by the 2nd return, over which the section-map Jacobian
[DP] genuinely varies ~3.7-4.0x.  The single-box mean-value image applies that whole loose
[DP] once to the offset, so its width tracks the box's OWN inflation, not the true image.

This applies the standard CAPD-style W-Z remedy: never let a single propagated box get near
35x.  Subdivide the h-set's unstable extent into ``nu`` small sub-boxes BEFORE propagating,
run EACH sub-box through its OWN #675/#676 mean-value chain (each sub-box's own local [DP]
stays near-constant across its small extent, so its mean-value image stays tight relative to
ITS size), and take the UNION of the sub-box images (:func:`vti.union_interval_boxes`) as the
return's enclosure.  Each cell's whole chain runs in memory in one process; only per-return
FLOAT results are checkpointed, so the run is fully resumable at CELL granularity (never
backgrounded).

The covering ratio per return is computed EXACTLY as #675/#676 (``separation_over_width_ratio``
= edge-center separation / max edge enclosure width) but with the edge images taken from the
extreme cells' outer faces and the whole-set image from the union -- directly comparable to
#676's baseline (N=1 0.2838, N=2 0.027486, N=3 WALL).  ``nu=1`` reproduces #676's baseline
exactly (internal consistency check) and supplies the target h-set center M.  Ground-truth
edge points (a=+-1, propagated as pure points) give the true physical image extent to compare
the union enclosure against (#676's single box was ~35x the ORIGINAL h-set; is the union near
the TRUE extent?).  ru=1e-6 is NOT shrunk.

The subdivision MECHANISM is validated on closed-form / dense-sampled ground truth in
``tests/scripts/test_677_subdivision_covering.py`` (built and passing BEFORE this run).

Run (resumable; re-invoke until "ALL JOBS DONE", then --assemble):
    uv run python scripts/certify_677_wz_oterma_subdivided_covering.py           # one chunk
    uv run python scripts/certify_677_wz_oterma_subdivided_covering.py --assemble # final cert
State  -> data/677_wz_oterma_subdivided_state.json     (per-cell float results)
Cert   -> data/677_wz_oterma_subdivided_covering_certificate.json
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import _validated_taylor_integrator as vti  # noqa: E402

STATE = ROOT / "data" / "677_wz_oterma_subdivided_state.json"
CERT = ROOT / "data" / "677_wz_oterma_subdivided_covering_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2
OM1 = 1.0 - MU

X0, VX0 = 0.9522928423486199945, 1.23e-05  # first published W-Z Oterma section point
RU, RS = 1e-6, 1e-8  # FIXED, NOT shrunk
N_RETURNS = 3
TAU_MAX, N_STEPS = 40.0, 320  # SAME h=0.125 as #675/#676; tau_max 40 so N=3 (tau~24) reachable

RESOLUTIONS = [1, 2, 3, 5]  # nu (unstable-direction sub-boxes); ns=1 (stable extent already 1e-8)
BUDGET_S = float(os.environ.get("C677_BUDGET_S", "470"))  # per-invocation wall budget (s)

_E_EXPAND = [-1.0, 2.5733011]
_E_CONTRACT = [1.0, 2.5733011]


def _norm(v: list[float]) -> list[float]:
    n = math.hypot(v[0], v[1])
    return [v[0] / n, v[1] / n]


UDIR = _norm(_E_EXPAND)
SDIR = _norm(_E_CONTRACT)


def _phys(iv: Any, aval: Any, bval: Any) -> tuple[Any, Any]:
    """Physical (x, xdot) of h-set parameter point/interval (a, b) (a unstable, b stable)."""
    x = iv.mpf(X0) + aval * iv.mpf(RU * UDIR[0]) + bval * iv.mpf(RS * SDIR[0])
    vx = iv.mpf(VX0) + aval * iv.mpf(RU * UDIR[1]) + bval * iv.mpf(RS * SDIR[1])
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


def _iv2(iv: Any, lohi: list[float]) -> Any:
    return iv.mpf([lohi[0], lohi[1]])


def _cell_ic(
    iv: Any, a_lo: float, a_hi: float
) -> tuple[list[Any], list[Any], dict[str, list[Any]]]:
    """Center 5D state, box-IC 5D state, and 5D face offsets for the sub-box a in [a_lo,a_hi]."""
    a_c = 0.5 * (a_lo + a_hi)
    w0c = _lc_ic_onsection(iv, *_phys(iv, iv.mpf([a_c, a_c]), iv.mpf([0, 0])))
    box_ic = _lc_ic_onsection(iv, *_phys(iv, iv.mpf([a_lo, a_hi]), iv.mpf([-1, 1])))
    faces = {
        "box": (iv.mpf([a_lo, a_hi]), iv.mpf([-1, 1])),
        "left": (iv.mpf([a_lo, a_lo]), iv.mpf([-1, 1])),
        "right": (iv.mpf([a_hi, a_hi]), iv.mpf([-1, 1])),
    }
    offs: dict[str, list[Any]] = {}
    for name, (av, bv) in faces.items():
        wf = _lc_ic_onsection(iv, *_phys(iv, av, bv))
        offs[name] = [wf[i] - w0c[i] for i in range(5)]
    return w0c, box_ic, offs


def run_cell_chain(
    iv: Any, jet: Any, vjet: Any, a_lo: float, a_hi: float, *, point_only: bool
) -> list[dict[str, Any]]:
    """One sub-box's independent mean-value chain (centre point chain + box chain), all returns.

    Returns a per-return list of float results.  ``point_only`` runs ONLY the centre point
    chain (for ground-truth edge points a=+-1): image center per return, no box/DP/faces.
    """
    w0c, box_ic, offs = _cell_ic(iv, a_lo, a_hi)
    centre_state = w0c
    box_center = w0c
    dp_boxes: list[Any] = []
    box_wall: str | None = None
    out: list[dict[str, Any]] = []

    for k in range(1, N_RETURNS + 1):
        rc = _secmap(iv, jet, vjet, centre_state)
        if not rc["found"]:
            out.append({"n": k, "centre_wall": rc.get("reason")})
            break
        centre_state = rc["crossing_state"]
        wk_hat = centre_state
        phys_c = vti.lc_secondary_to_physical(iv, wk_hat, MU)
        row: dict[str, Any] = {
            "n": k,
            "tau_star": _f(rc["tau_star"].mid),
            "center_image_xvx": [_f(phys_c[0].mid), _f(phys_c[2].mid)],
        }
        if point_only:
            out.append(row)
            continue
        if box_wall is None:
            rb = _secmap(iv, jet, vjet, box_ic)
            if not rb["found"]:
                box_wall = rb.get("reason")
                row["box_wall"] = box_wall
            else:
                dp_k = rb["section_jacobian"]
                dp_boxes.append(dp_k)
                box_off = [box_ic[i] - box_center[i] for i in range(5)]
                box_ic = vti.section_map_meanvalue_image(iv, wk_hat, dp_k, box_off)
                box_center = wk_hat
                jac_comp = vti.compose_section_jacobians(iv, dp_boxes)
                imgs: dict[str, list[list[float]]] = {}
                for name in ("box", "left", "right"):
                    reg_img = vti.section_map_meanvalue_image(iv, wk_hat, jac_comp, offs[name])
                    ph = vti.lc_secondary_to_physical(iv, reg_img, MU)
                    imgs[name] = [[_f(ph[0].a), _f(ph[0].b)], [_f(ph[2].a), _f(ph[2].b)]]
                row["images_xvx"] = imgs  # {face: [[x_lo,x_hi],[vx_lo,vx_hi]]}
                row["box_ic_max_halfwidth"] = max(_wid(c) for c in box_ic) / 2.0
                row["composed_jac_width"] = max(
                    _wid(jac_comp[i][j]) for i in range(5) for j in range(5)
                )
                row["leg_jac_width"] = max(_wid(dp_k[i][j]) for i in range(5) for j in range(5))
        else:
            row["box_wall"] = f"upstream: {box_wall}"
        out.append(row)
    return out


def _job_list() -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []

    def cells_for(nu: int) -> list[dict[str, Any]]:
        edges = [(-1.0 + 2.0 * i / nu) for i in range(nu + 1)]
        return [
            {
                "id": f"nu{nu}_c{i}",
                "kind": "cell",
                "nu": nu,
                "ci": i,
                "a_lo": edges[i],
                "a_hi": edges[i + 1],
            }
            for i in range(nu)
        ]

    # coarse grids first, then cheap ground-truth edge points (so a complete decisive
    # nu1/nu2/nu3 result with a true-extent reference lands before the expensive fine grid)
    for nu in RESOLUTIONS:
        if nu <= 3:
            jobs.extend(cells_for(nu))
    jobs.append({"id": "gt_left", "kind": "point", "a_lo": -1.0, "a_hi": -1.0})
    jobs.append({"id": "gt_right", "kind": "point", "a_lo": 1.0, "a_hi": 1.0})
    for nu in RESOLUTIONS:
        if nu > 3:
            jobs.extend(cells_for(nu))
    return jobs


def _load_state() -> dict[str, Any]:
    if STATE.exists():
        return json.loads(STATE.read_text())  # type: ignore[no-any-return]
    return {"done": {}}


def run_chunk() -> None:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv
    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)

    state = _load_state()
    done = state["done"]
    jobs = _job_list()
    pending = [j for j in jobs if j["id"] not in done]
    print(f"[#677] {len(done)}/{len(jobs)} jobs done; {len(pending)} pending.")
    t0 = time.time()
    for job in pending:
        if time.time() - t0 > BUDGET_S:
            print(f"[#677] budget hit; exit to be re-invoked ({len(done)}/{len(jobs)} done).")
            break
        ts = time.time()
        res = run_cell_chain(
            iv, jet, vjet, job["a_lo"], job["a_hi"], point_only=(job["kind"] == "point")
        )
        done[job["id"]] = {"job": job, "returns": res}
        STATE.write_text(json.dumps(state, indent=2) + "\n")
        last = res[-1] if res else {}
        print(
            f"[#677] job {job['id']} done in {time.time() - ts:.0f}s "
            f"(reached n={last.get('n')}, {len(done)}/{len(jobs)} total)."
        )
    if len(done) == len(jobs):
        print("[#677] ALL JOBS DONE -- run with --assemble to build the certificate.")


def _cover_from_union(
    iv: Any,
    cx: float,
    cvx: float,
    left: dict[str, Any],
    right: dict[str, Any],
    box: dict[str, Any],
) -> dict[str, Any]:
    """#676's exact covering-ratio computation, on unioned physical image boxes."""
    lx, lvx = _iv2(iv, left["x"]), _iv2(iv, left["vx"])
    rx, rvx = _iv2(iv, right["x"]), _iv2(iv, right["vx"])
    bx, bvx = _iv2(iv, box["x"]), _iv2(iv, box["vx"])
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
        "edge_center_separation": sep,
        "per_edge_enclosure_width": edgew,
        "separation_over_width_ratio": (sep / edgew if edgew > 0 else None),
        "stable_box_width": _wid(bx),
        "covers": bool(cov["covers"]),
        "stable_containment": bool(cov["stable_containment"]),
        "unstable_exit": bool(cov["unstable_exit"]),
    }


def assemble() -> None:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv
    state = _load_state()
    done = state["done"]

    # index cells by (nu, ci)
    by_nu: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for rec in done.values():
        job = rec["job"]
        if job["kind"] != "cell":
            continue
        by_nu.setdefault(job["nu"], {})[job["ci"]] = rec["returns"]

    # global center chain = nu=1 cell (a in [-1,1], center a=0) -- the #676 whole-h-set center
    center_chain = by_nu[1][0]
    center_by_n = {r["n"]: r.get("center_image_xvx") for r in center_chain}

    cert: dict[str, Any] = {
        "task": "#677",
        "stage": "Stage 11 -- finer sub-boxing (subdivide-before-propagate) composed covering",
        "system": {"mu": MU, "C": C_LEVEL, "h_energy": H_ENERGY},
        "hset_N": {
            "center_x_xdot": [X0, VX0],
            "ru": RU,
            "rs": RS,
            "note": "FIXED size, NOT shrunk (identical to #674/#675/#676).",
        },
        "search_window": {"tau_max": TAU_MAX, "n_steps": N_STEPS, "step_h": TAU_MAX / N_STEPS},
        "method": (
            "Subdivide the h-set unstable extent into nu sub-boxes; each runs its OWN #675/#676 "
            "mean-value chain (centre point chain + box chain); union the sub-box images "
            "(vti.union_interval_boxes) per return. Covering ratio computed as #676 "
            "(sep/edgew) with edges from extreme cells' outer faces, whole-set from the union. "
            "M (target h-set) built from the nu=1 whole-h-set center image, identical to #676."
        ),
        "mechanism_control": "tests/scripts/test_677_subdivision_covering.py",
        "resolutions": {},
    }

    # ground-truth true physical extent per return (edge points a=+-1 + all cell centers)
    gt: dict[int, dict[str, Any]] = {}
    for n in range(1, N_RETURNS + 1):
        pts: list[list[float]] = []
        for key in ("gt_left", "gt_right"):
            if key in done:
                for r in done[key]["returns"]:
                    if r["n"] == n and "center_image_xvx" in r:
                        pts.append(r["center_image_xvx"])
        for cells in by_nu.values():
            for rows in cells.values():
                for r in rows:
                    if r["n"] == n and "center_image_xvx" in r:
                        pts.append(r["center_image_xvx"])
        if pts:
            xs = [p[0] for p in pts]
            vs = [p[1] for p in pts]
            gt[n] = {
                "x_extent": [min(xs), max(xs)],
                "vx_extent": [min(vs), max(vs)],
                "x_halfwidth": 0.5 * (max(xs) - min(xs)),
                "vx_halfwidth": 0.5 * (max(vs) - min(vs)),
                "n_points": len(pts),
            }
    cert["ground_truth_true_extent"] = gt

    for nu, cells in sorted(by_nu.items()):
        res_rows: list[dict[str, Any]] = []
        for n in range(1, N_RETURNS + 1):
            # collect per-cell images at return n; abort this return if any cell walled/missing
            cell_imgs: list[dict[str, Any]] = []
            walled = None
            for ci in range(nu):
                rows = cells.get(ci, [])
                rn = next((r for r in rows if r["n"] == n), None)
                if rn is None or "images_xvx" not in rn:
                    walled = (
                        (rn or {}).get("box_wall") or (rn or {}).get("centre_wall") or "missing"
                    )
                    break
                cell_imgs.append(
                    {
                        "ci": ci,
                        "img": rn["images_xvx"],
                        "box_hw": rn.get("box_ic_max_halfwidth"),
                        "jacw": rn.get("composed_jac_width"),
                    }
                )
            if walled is not None:
                res_rows.append({"n": n, "result": "WALL", "reason": walled})
                continue
            # union of all cells' "box" faces; extreme cells' outer faces for the edges
            box_boxes = []
            for c in cell_imgs:
                bx = c["img"]["box"]
                box_boxes.append([_iv2(iv, bx[0]), _iv2(iv, bx[1])])
            union = vti.union_interval_boxes(iv, box_boxes)
            box_u = {"x": [_f(union[0].a), _f(union[0].b)], "vx": [_f(union[1].a), _f(union[1].b)]}
            left_img = cell_imgs[0]["img"]["left"]  # leftmost cell, a=-1
            right_img = cell_imgs[-1]["img"]["right"]  # rightmost cell, a=+1
            left = {"x": left_img[0], "vx": left_img[1]}
            right = {"x": right_img[0], "vx": right_img[1]}
            cxvx = center_by_n.get(n)
            if cxvx is None:
                res_rows.append({"n": n, "result": "WALL", "reason": "global center missing"})
                continue
            cov = _cover_from_union(iv, cxvx[0], cxvx[1], left, right, box_u)
            union_hw = max(0.5 * (union[0].b - union[0].a), 0.5 * (union[1].b - union[1].a))
            row: dict[str, Any] = {
                "n": n,
                **cov,
                "union_image_max_halfwidth": _f(union_hw),
                "max_cell_box_ic_halfwidth": max(c["box_hw"] for c in cell_imgs),
                "max_composed_jac_width": max(c["jacw"] for c in cell_imgs),
            }
            if n in gt:
                true_hw = max(gt[n]["x_halfwidth"], gt[n]["vx_halfwidth"])
                row["union_over_true_extent"] = (_f(union_hw) / true_hw) if true_hw > 0 else None
            res_rows.append(row)
        cert["resolutions"][f"nu{nu}"] = res_rows

    # concise verdict
    def ratio(nu: int, n: int) -> float | None:
        for r in cert["resolutions"].get(f"nu{nu}", []):
            if r["n"] == n:
                v = r.get("separation_over_width_ratio")
                return float(v) if v is not None else None
        return None

    def over_true(nu: int, n: int) -> float | None:
        for r in cert["resolutions"].get(f"nu{nu}", []):
            if r["n"] == n:
                v = r.get("union_over_true_extent")
                return float(v) if v is not None else None
        return None

    lines = []
    for n in (2, 3):
        parts = [f"nu{nu}:{ratio(nu, n):.4f}" for nu in RESOLUTIONS if ratio(nu, n) is not None]
        ot = [
            f"nu{nu}:{over_true(nu, n):.2f}x" for nu in RESOLUTIONS if over_true(nu, n) is not None
        ]
        if parts:
            lines.append(f"N={n} ratio {{{', '.join(parts)}}} | union/true {{{', '.join(ot)}}}")
    cert["verdict"] = " || ".join(lines) if lines else "insufficient data"
    certified = {
        f"nu{nu}": [r["n"] for r in cert["resolutions"].get(f"nu{nu}", []) if r.get("covers")]
        for nu in RESOLUTIONS
    }
    cert["certified_returns"] = certified
    CERT.write_text(json.dumps(cert, indent=2) + "\n")
    print("[#677] ===== ASSEMBLED =====")
    print(f"[#677] {cert['verdict']}")
    print(f"[#677] certified_returns={certified}")
    print(f"[#677] certificate -> {CERT.relative_to(ROOT)}")


def main() -> None:
    if "--assemble" in sys.argv:
        assemble()
    else:
        run_chunk()


if __name__ == "__main__":
    main()
