"""#678 -- HIGHER-RESOLUTION subdivision covering (nu=8, nu=10) at fixed ru=1e-6.

Stage 12 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668``-``#677``).

#677 subdivided the h-set into ``nu`` sub-boxes before propagating, and found the N=2
composed covering ratio improving fairly consistently ~2-3x per subdivision step
(nu1 0.0275 -> nu2 0.0947 -> nu3 0.2020 -> nu5 0.4961) but never crossing the
certification threshold 1.0.  #677's own bullet flagged, EXPLICITLY as untested
speculation, that nu~8-10 MIGHT cross 1.0.  This stage tests that directly by extending
the exact #677 sub-box mean-value chain to nu=8 and nu=10, and reports whether the N=2
trend holds, saturates below 1.0, or reverses -- numerically, whatever it is.

It also re-tests N=3 at the higher resolutions.  #678's diagnosis of #677's nu5 N=3 wall
(instrumented replay, see ``tests/scripts/test_678_crossing_isolation_refine.py`` and the
OUTSTANDING.md writeup) showed the wall is NOT a within-step precision issue but a
multi-step crossing-time SMEAR: the N=2 box (hw~5e-5) inflates ~157x over the long tau~24
arc, spreading the crossing across ~5 integration steps, so no single step shows a strict
endpoint sign change.  The established remedy for box inflation in this stage-series is
finer subdivision; nu8/nu10 give a tighter N=2 box, so this run also answers empirically
whether a tighter box brings the N=3 crossing back within a single-step bracket (with the
new two-sided bracket refinement in ``isolate_section_crossing``) so N=3 can be evaluated
instead of walling.

Everything is identical to #677 except the resolutions and a finer, resumable
checkpoint: per (cell, RETURN) rather than per cell, plus in-process multiprocessing
across the (independent) sub-boxes -- so each synchronous invocation does a bounded chunk
and the run is never backgrounded.  ru=1e-6 is NOT shrunk.

Run (resumable; re-invoke until "ALL CELLS DONE", then --assemble):
    uv run python scripts/certify_678_wz_oterma_higherres_covering.py            # one chunk
    uv run python scripts/certify_678_wz_oterma_higherres_covering.py --assemble # final cert
State -> data/678_wz_oterma_higherres_state.json
Cert  -> data/678_wz_oterma_higherres_covering_certificate.json
"""

from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import _validated_taylor_integrator as vti  # noqa: E402
import certify_677_wz_oterma_subdivided_covering as base  # noqa: E402

STATE = ROOT / "data" / "678_wz_oterma_higherres_state.json"
CERT = ROOT / "data" / "678_wz_oterma_higherres_covering_certificate.json"
S677 = ROOT / "data" / "677_wz_oterma_subdivided_state.json"
C677 = ROOT / "data" / "677_wz_oterma_subdivided_covering_certificate.json"

RESOLUTIONS = [8, 10]  # higher nu than #677's {1,2,3,5}
N_RETURNS = 3
N_WORKERS = int(os.environ.get("C678_WORKERS", "4"))  # = M3 performance cores (pure-Python mpmath)
BUDGET_LAUNCH_S = float(os.environ.get("C678_BUDGET_S", "280"))  # stop launching new waves after
MAX_WAVE_S_GUESS = 260  # a single return's worst-case wall time (informs the Bash timeout)


def _iv_setup() -> Any:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    return mp.iv


def _new_cell(nu: int, ci: int, a_lo: float, a_hi: float) -> dict[str, Any]:
    return {
        "id": f"nu{nu}_c{ci}",
        "nu": nu,
        "ci": ci,
        "a_lo": a_lo,
        "a_hi": a_hi,
        "k": 1,  # next return to compute
        "centre_state": None,  # serialized 5-vec (set at first return)
        "box_ic": None,
        "box_center": None,
        "dp_boxes": [],  # list of serialized 5x5
        "centre_wall": None,
        "box_wall": None,
        "rows": [],  # per-return float results (as #677)
        "done": False,
    }


def _cells_for(nu: int) -> list[dict[str, Any]]:
    edges = [(-1.0 + 2.0 * i / nu) for i in range(nu + 1)]
    return [_new_cell(nu, i, edges[i], edges[i + 1]) for i in range(nu)]


def advance_one_return(cell: dict[str, Any]) -> dict[str, Any]:
    """Advance one sub-box cell by exactly ONE return (center chain + box chain).

    Operates on / returns a JSON-serializable cell dict; interval state
    (centre_state, box_ic, box_center, dp_boxes) is carried as pickled mpmath
    intervals inside the (json-stored) 'blob' produced by run_chunk's (de)serialize.
    Here it is passed already deserialized into live intervals under keys with the
    '_iv' suffix.  Mutates a copy and returns it.
    """
    iv = _iv_setup()
    jet = vti.make_cr3bp_lc_secondary_jet(base.H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(base.H_ENERGY)
    k = cell["k"]

    w0c, box_ic0, offs = base._cell_ic(iv, cell["a_lo"], cell["a_hi"])
    if k == 1:
        centre_state = w0c
        box_ic = box_ic0
        box_center = w0c
        dp_boxes: list[Any] = []
    else:
        centre_state = cell["_centre_iv"]
        box_ic = cell["_box_ic_iv"]
        box_center = cell["_box_center_iv"]
        dp_boxes = cell["_dp_boxes_iv"]

    # --- center point chain ---
    rc = base._secmap(iv, jet, vjet, centre_state)
    if not rc["found"]:
        cell["centre_wall"] = rc.get("reason")
        cell["rows"].append({"n": k, "centre_wall": rc.get("reason")})
        cell["done"] = True
        return cell
    centre_state = rc["crossing_state"]
    wk_hat = centre_state
    phys_c = vti.lc_secondary_to_physical(iv, wk_hat, base.MU)
    row: dict[str, Any] = {
        "n": k,
        "tau_star": base._f(rc["tau_star"].mid),
        "center_image_xvx": [base._f(phys_c[0].mid), base._f(phys_c[2].mid)],
    }

    # --- box chain (unless already walled upstream) ---
    if cell["box_wall"] is None:
        rb = base._secmap(iv, jet, vjet, box_ic)
        if not rb["found"]:
            cell["box_wall"] = rb.get("reason")
            row["box_wall"] = cell["box_wall"]
        else:
            dp_k = rb["section_jacobian"]
            dp_boxes = [*dp_boxes, dp_k]
            box_off = [box_ic[i] - box_center[i] for i in range(5)]
            box_ic = vti.section_map_meanvalue_image(iv, wk_hat, dp_k, box_off)
            box_center = wk_hat
            jac_comp = vti.compose_section_jacobians(iv, dp_boxes)
            imgs: dict[str, list[list[float]]] = {}
            for name in ("box", "left", "right"):
                reg_img = vti.section_map_meanvalue_image(iv, wk_hat, jac_comp, offs[name])
                ph = vti.lc_secondary_to_physical(iv, reg_img, base.MU)
                imgs[name] = [
                    [base._f(ph[0].a), base._f(ph[0].b)],
                    [base._f(ph[2].a), base._f(ph[2].b)],
                ]
            row["images_xvx"] = imgs
            row["box_ic_max_halfwidth"] = max(base._wid(c) for c in box_ic) / 2.0
            row["composed_jac_width"] = max(
                base._wid(jac_comp[i][j]) for i in range(5) for j in range(5)
            )
            row["leg_jac_width"] = max(base._wid(dp_k[i][j]) for i in range(5) for j in range(5))
    else:
        row["box_wall"] = f"upstream: {cell['box_wall']}"

    cell["rows"].append(row)
    # persist live interval state for the next return
    cell["_centre_iv"] = centre_state
    cell["_box_ic_iv"] = box_ic
    cell["_box_center_iv"] = box_center
    cell["_dp_boxes_iv"] = dp_boxes
    cell["k"] = k + 1
    # stop the cell once the box has walled (center-only continuation is not needed:
    # ground-truth extent points are reused from #677) or all returns are done
    if cell["box_wall"] is not None or cell["k"] > N_RETURNS:
        cell["done"] = True
    return cell


# ------- (de)serialization of live interval state to/from the JSON state ------- #
def _ser(iv: Any, cell: dict[str, Any]) -> dict[str, Any]:
    import pickle

    out = {k: v for k, v in cell.items() if not k.startswith("_")}
    if cell.get("_centre_iv") is not None:
        blob = {
            "centre": cell["_centre_iv"],
            "box_ic": cell["_box_ic_iv"],
            "box_center": cell["_box_center_iv"],
            "dp_boxes": cell["_dp_boxes_iv"],
        }
        out["_blob_hex"] = pickle.dumps(blob).hex()
    return out


def _des(iv: Any, cell: dict[str, Any]) -> dict[str, Any]:
    import pickle

    c = dict(cell)
    if "_blob_hex" in c:
        blob = pickle.loads(bytes.fromhex(c.pop("_blob_hex")))
        c["_centre_iv"] = blob["centre"]
        c["_box_ic_iv"] = blob["box_ic"]
        c["_box_center_iv"] = blob["box_center"]
        c["_dp_boxes_iv"] = blob["dp_boxes"]
    return c


def _load_state() -> dict[str, Any]:
    if STATE.exists():
        return json.loads(STATE.read_text())  # type: ignore[no-any-return]
    cells: dict[str, Any] = {}
    for nu in RESOLUTIONS:
        for c in _cells_for(nu):
            cells[c["id"]] = c
    return {"cells": cells}


def run_chunk() -> None:
    iv = _iv_setup()
    state = _load_state()
    cells = state["cells"]
    ndone = sum(1 for c in cells.values() if c["done"])
    print(f"[#678] {ndone}/{len(cells)} cells done; workers={N_WORKERS}", flush=True)
    ctx = get_context("spawn")
    t0 = time.time()
    with ctx.Pool(N_WORKERS) as pool:
        while True:
            pending = [c for c in cells.values() if not c["done"]]
            if not pending:
                print("[#678] ALL CELLS DONE -- run with --assemble to build the certificate.")
                break
            if time.time() - t0 > BUDGET_LAUNCH_S:
                print(f"[#678] launch budget hit; exit to be re-invoked ({ndone} cells done).")
                break
            wave = pending[:N_WORKERS]
            live = [_des(iv, c) for c in wave]
            ts = time.time()
            results = pool.map(advance_one_return, live)
            for r in results:
                cells[r["id"]] = _ser(iv, r)
            state["cells"] = cells
            STATE.write_text(json.dumps(state, indent=2) + "\n")
            ndone = sum(1 for c in cells.values() if c["done"])
            info = ", ".join(
                f"{r['id']}@n{r['rows'][-1]['n']}"
                f"{'(wall)' if r.get('box_wall') or r.get('centre_wall') else ''}"
                for r in results
            )
            print(
                f"[#678] wave {len(wave)} cells +1 return in {time.time() - ts:.0f}s "
                f"({ndone}/{len(cells)} done): {info}",
                flush=True,
            )


def _reuse_677_nu_le5(iv: Any) -> dict[str, Any]:
    """Pull #677's already-computed nu1/2/3/5 per-return rows + ground truth verbatim."""
    st = json.loads(S677.read_text())
    by_nu: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for rec in st["done"].values():
        job = rec["job"]
        if job.get("kind") != "cell":
            continue
        by_nu.setdefault(job["nu"], {})[job["ci"]] = rec["returns"]
    return by_nu


def assemble() -> None:
    iv = _iv_setup()
    state = _load_state()
    cells = state["cells"]

    by_nu: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for c in cells.values():
        by_nu.setdefault(c["nu"], {})[c["ci"]] = c["rows"]
    # also fold in #677's nu<=5 rows for a single combined verdict/trend line
    by_nu677 = _reuse_677_nu_le5(iv)
    for nu, cellrows in by_nu677.items():
        by_nu.setdefault(nu, cellrows)

    # ground-truth true extent: reuse #677's certificate (identical map, same edge pts)
    gt = json.loads(C677.read_text()).get("ground_truth_true_extent", {})
    gt = {int(k): v for k, v in gt.items()}
    # #677's whole-h-set center chain (nu=1 cell) = the target-h-set center M
    center_chain = by_nu677[1][0]
    center_by_n = {r["n"]: r.get("center_image_xvx") for r in center_chain}

    cert: dict[str, Any] = {
        "task": "#678",
        "stage": "Stage 12 -- higher-resolution (nu=8,10) subdivision composed covering",
        "system": {"mu": base.MU, "C": base.C_LEVEL, "h_energy": base.H_ENERGY},
        "hset_N": {
            "center_x_xdot": [base.X0, base.VX0],
            "ru": base.RU,
            "rs": base.RS,
            "note": "FIXED size, NOT shrunk (identical to #674-#677).",
        },
        "search_window": {
            "tau_max": base.TAU_MAX,
            "n_steps": base.N_STEPS,
            "step_h": base.TAU_MAX / base.N_STEPS,
        },
        "method": (
            "Identical #675/#676/#677 sub-box mean-value chain, resolutions nu=8,10 "
            "(nu<=5 reused verbatim from #677). Covering ratio = #676's sep/edgew with "
            "edges from extreme cells' outer faces, whole-set from the sub-box union. M "
            "and ground-truth extent identical to #677."
        ),
        "resolutions": {},
        "ground_truth_true_extent": gt,
    }

    all_nu = sorted(by_nu.keys())
    for nu in all_nu:
        cellsd = by_nu[nu]
        res_rows: list[dict[str, Any]] = []
        for n in range(1, N_RETURNS + 1):
            cell_imgs: list[dict[str, Any]] = []
            walled = None
            for ci in range(nu):
                rows = cellsd.get(ci, [])
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
            box_boxes = []
            for c in cell_imgs:
                bx = c["img"]["box"]
                box_boxes.append([base._iv2(iv, bx[0]), base._iv2(iv, bx[1])])
            union = vti.union_interval_boxes(iv, box_boxes)
            box_u = {
                "x": [base._f(union[0].a), base._f(union[0].b)],
                "vx": [base._f(union[1].a), base._f(union[1].b)],
            }
            left = {"x": cell_imgs[0]["img"]["left"][0], "vx": cell_imgs[0]["img"]["left"][1]}
            right = {"x": cell_imgs[-1]["img"]["right"][0], "vx": cell_imgs[-1]["img"]["right"][1]}
            cxvx = center_by_n.get(n)
            if cxvx is None:
                res_rows.append({"n": n, "result": "WALL", "reason": "global center missing"})
                continue
            cov = base._cover_from_union(iv, cxvx[0], cxvx[1], left, right, box_u)
            union_hw = max(0.5 * (union[0].b - union[0].a), 0.5 * (union[1].b - union[1].a))
            row: dict[str, Any] = {
                "n": n,
                **cov,
                "union_image_max_halfwidth": base._f(union_hw),
                "max_cell_box_ic_halfwidth": max(c["box_hw"] for c in cell_imgs),
                "max_composed_jac_width": max(c["jacw"] for c in cell_imgs),
            }
            if n in gt:
                true_hw = max(gt[n]["x_halfwidth"], gt[n]["vx_halfwidth"])
                row["union_over_true_extent"] = (
                    (base._f(union_hw) / true_hw) if true_hw > 0 else None
                )
            res_rows.append(row)
        cert["resolutions"][f"nu{nu}"] = res_rows

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
    for n in (1, 2, 3):
        parts = [f"nu{nu}:{ratio(nu, n):.4f}" for nu in all_nu if ratio(nu, n) is not None]
        walls = [f"nu{nu}:WALL" for nu in all_nu if ratio(nu, n) is None]
        if parts or walls:
            lines.append(f"N={n} ratio {{{', '.join(parts + walls)}}}")
    cert["verdict"] = " || ".join(lines) if lines else "insufficient data"
    cert["certified_returns"] = {
        f"nu{nu}": [r["n"] for r in cert["resolutions"].get(f"nu{nu}", []) if r.get("covers")]
        for nu in all_nu
    }
    CERT.write_text(json.dumps(cert, indent=2) + "\n")
    print("[#678] ===== ASSEMBLED =====")
    for ln in lines:
        print(f"[#678] {ln}")
    print(f"[#678] certified_returns={cert['certified_returns']}")
    print(f"[#678] certificate -> {CERT.relative_to(ROOT)}")


def main() -> None:
    if "--assemble" in sys.argv:
        assemble()
    else:
        run_chunk()


if __name__ == "__main__":
    main()
