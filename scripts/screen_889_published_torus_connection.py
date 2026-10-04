"""#889 -- published positive control: the SIAM 2025 Jupiter-Europa PERTBP 3:4 -> 5:6 connection.

Driver for :mod:`cyclerfinder.search.pertbp_strob_889`. Each ``--stage`` finishes in a
few minutes, prints timestamped progress and checkpoints to
``data/found/889_published_torus_connection/<stage>[_<variant>].json`` (``.npz`` for
arrays). Later stages read earlier checkpoints.

Variants (``--variant``) of the rotation numbers, both run, neither chosen in advance:
  table1   omega = 4 pi^2 / T of the CNSNS Table 1 orbits (C = 3.0024);
  printed  omega = the SIAM Sec. 8.3 printed values 1.558039 and 1.030011.

Stages
------
po        CNSNS Table 1 orbits at both published mass ratios; the fixed-period orbits
          of the ``printed`` variant.
tori      circles at e = 0, continued in e to 0.0094 at fixed omega (3:4 N = 511,
          5:6 N = 1023): residuals, reversibility, resolution.
bundles   constant-multiplier bundles, unit normalisation at theta = 0, second-order
          terms; Fig. 9 / Fig. 11 by-eye comparisons.
forward   W_u(2.39703, -77.73428) and W_s(1.83093, 202.62277) evaluated directly and
          compared with the printed point (no inversion involved).
refine    the connection refined from the printed parameters at two chart radii a factor
          of ten apart; conditioning (SVD) and decomposition of the difference from the
          printed point along the near-null direction; re-seeding test.
verify    #882 review criteria: contraction rates, junction agreement, perturbation test.

Run:  OPENBLAS_NUM_THREADS=4 uv run python scripts/screen_889_published_torus_connection.py \
          --stage tori --variant table1
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.search.pertbp_strob_889 as pm  # noqa: E402
from cyclerfinder.search.ccr4bp_strob_connection import fourier_eval, node_angles  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "889_published_torus_connection"
T_START = time.time()
MU = pm.MU_JE_CNSNS
TWO_PI = 2.0 * math.pi
E = pm.ECC_EUROPA
N34 = 511
N56 = 1023


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{stamp} +{time.time() - T_START:7.1f}s] {msg}", flush=True)


def _jsonable(o: Any) -> Any:
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.floating | np.integer):
        return o.item()
    raise TypeError(type(o))


def save(name: str, payload: dict[str, Any]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / f"{name}.json").write_text(json.dumps(payload, indent=2, default=_jsonable))
    log(f"wrote {name}.json")


def load(name: str) -> dict[str, Any]:
    out: dict[str, Any] = json.loads((OUT_DIR / f"{name}.json").read_text())
    return out


# ---------------------------------------------------------------------------


def stage_po(_: str) -> None:
    out: dict[str, Any] = {}
    for label, mu in (("mu_cnsns", pm.MU_JE_CNSNS), ("mu_aas21651", pm.MU_JE_AAS21651)):
        for name, tab in (("3:4", pm.TABLE1_3_4), ("5:6", pm.TABLE1_5_6)):
            po = pm.po_at_fixed_x(mu, tab["x"], tab["ydot"], tab["T"])
            row = {
                "mu": mu,
                "x0": po.state[0],
                "ydot0": po.state[3] - po.state[0],
                "period": po.period,
                "lam_s": po.lam_s,
                "lam_u": po.lam_u,
                "jacobi": po.jacobi,
                "omega": po.omega,
                "d_ydot": po.state[3] - po.state[0] - tab["ydot"],
                "d_period": po.period - tab["T"],
                "rel_d_lam_u": po.lam_u / tab["lam_u"] - 1.0,
                "rel_d_lam_s": po.lam_s / tab["lam_s"] - 1.0,
                "residual": po.residual,
            }
            out[f"{label} {name}"] = row
            log(
                f"{label} {name}: dT {row['d_period']:.2e} dydot {row['d_ydot']:.2e} "
                f"rel dlam_u {row['rel_d_lam_u']:.2e} C {po.jacobi:.12f} omega {po.omega:.10f}"
            )
    # fixed-period orbits at the printed rotation numbers (variant 'printed')
    for name, tab, om in (
        ("3:4", pm.TABLE1_3_4, pm.PRINTED_OMEGA_U),
        ("5:6", pm.TABLE1_5_6, pm.PRINTED_OMEGA_S),
    ):
        x_guess = tab["x"]
        vy_guess = tab["ydot"]
        po = pm.po_at_fixed_period(MU, x_guess, vy_guess, pm.period_from_omega(om))
        out[f"printed {name}"] = {
            "x0": po.state[0],
            "ydot0": po.state[3] - po.state[0],
            "period": po.period,
            "omega": po.omega,
            "lam_s": po.lam_s,
            "lam_u": po.lam_u,
            "jacobi": po.jacobi,
            "residual": po.residual,
        }
        log(f"printed-omega {name}: T {po.period:.10f} x0 {po.state[0]:.12f} C {po.jacobi:.10f}")
    save("po", out)


def _po_for(variant: str, name: str) -> pm.PeriodicOrbit:
    tab = pm.TABLE1_3_4 if name == "3:4" else pm.TABLE1_5_6
    if variant == "table1":
        return pm.po_at_fixed_x(MU, tab["x"], tab["ydot"], tab["T"])
    om = pm.PRINTED_OMEGA_U if name == "3:4" else pm.PRINTED_OMEGA_S
    return pm.po_at_fixed_period(MU, tab["x"], tab["ydot"], pm.period_from_omega(om))


def _circle_path(variant: str, name: str) -> Path:
    return OUT_DIR / f"circle_{variant}_{name.replace(':', '_')}.npz"


def load_circle(variant: str, name: str) -> pm.Circle:
    d = np.load(_circle_path(variant, name))
    return pm.Circle(
        mu=float(d["mu"]),
        e=float(d["e"]),
        omega=float(d["omega"]),
        nodes=d["nodes"],
        residual=float(d["residual"]),
        converged=True,
    )


def stage_tori(variant: str) -> None:
    if variant == "table1_hi":
        stage_tori_upsample()
        return
    out: dict[str, Any] = {}
    for name, n_nodes in (("3:4", N34), ("5:6", N56)):
        po = _po_for(variant, name)
        log(f"{name}: seed circle N = {n_nodes}, omega = {po.omega:.12f}, T = {po.period:.12f}")
        c0 = pm.seed_circle(po, n_nodes)
        on0, off0 = pm.circle_residual(c0)
        log(f"  e = 0 residual on/off-grid {on0:.2e} / {off0:.2e}; tail {c0.tail():.2e}")
        chain = pm.continue_circle_in_e(c0, E, verbose=True)
        c = chain[-1]
        if abs(c.e - E) > 1e-15:
            log(f"  FAILED: stopped at e = {c.e}")
            out[name] = {"failed_at_e": c.e}
            continue
        on, off = pm.circle_residual(c)
        rev = pm.reversibility_defect(c)
        dk0 = c.dstate(0.0)
        row = {
            "n_nodes": n_nodes,
            "omega": c.omega,
            "period_parent": po.period,
            "jacobi_parent": po.jacobi,
            "e_steps": len(chain) - 1,
            "residual_on_grid": on,
            "residual_off_grid": off,
            "reversibility_defect": rev,
            "fourier_tail": c.tail(),
            "K0": c.nodes[0],
            "px_K0_symmetry": c.nodes[0, 2],
            "max_shift_from_e0": float(np.max(np.abs(c.nodes - c0.nodes))),
            "dK0": dk0,
        }
        out[name] = row
        log(
            f"  e = {E}: on/off {on:.2e}/{off:.2e} rev {rev:.2e} tail {c.tail():.2e} "
            f"K0 {np.array2string(c.nodes[0], precision=10)}"
        )
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(
            _circle_path(variant, name),
            mu=c.mu,
            e=c.e,
            omega=c.omega,
            nodes=c.nodes,
            residual=on,
        )
    save(f"tori_{variant}", out)


def stage_tori_upsample() -> None:
    """Resolution check: the ``table1`` circles interpolated to 2N+1 nodes and
    re-corrected at e = 0.0094 (variant ``table1_hi``)."""
    out: dict[str, Any] = {}
    for name in ("3:4", "5:6"):
        c = load_circle("table1", name)
        n_new = 2 * c.n_nodes + 1
        nodes = fourier_eval(c.nodes, node_angles(n_new))
        guess = pm.Circle(mu=c.mu, e=c.e, omega=c.omega, nodes=nodes)
        log(f"{name}: upsample {c.n_nodes} -> {n_new}, correct at e = {c.e}")
        ch = pm.correct_circle(guess, c.e, tol=1e-12, verbose=True)
        on, off = pm.circle_residual(ch)
        rev = pm.reversibility_defect(ch)
        out[name] = {
            "n_nodes": n_new,
            "residual_on_grid": on,
            "residual_off_grid": off,
            "reversibility_defect": rev,
            "fourier_tail": ch.tail(),
            "max_change_vs_table1": float(np.max(np.abs(ch.nodes - nodes))),
            "history": ch.history,
        }
        log(
            f"  on/off {on:.2e}/{off:.2e} rev {rev:.2e} "
            f"change {out[name]['max_change_vs_table1']:.2e}"
        )
        np.savez(
            _circle_path("table1_hi", name),
            mu=ch.mu,
            e=ch.e,
            omega=ch.omega,
            nodes=ch.nodes,
            residual=on,
        )
    save("tori_table1_hi", out)


def _bundle_path(variant: str, name: str) -> Path:
    return OUT_DIR / f"bundles_{variant}_{name.replace(':', '_')}.npz"


def load_bundles(variant: str, name: str) -> pm.Bundles:
    d = np.load(_bundle_path(variant, name))
    return pm.Bundles(
        lam_u=float(d["lam_u"]),
        lam_s=float(d["lam_s"]),
        v_u=d["v_u"],
        v_s=d["v_s"],
        w_u=d["w_u"],
        w_s=d["w_s"],
        offgrid_u=float(d["offgrid_u"]),
        offgrid_s=float(d["offgrid_s"]),
        tail_u=float(d["tail_u"]),
        tail_s=float(d["tail_s"]),
    )


def _second_order_check(c: pm.Circle, b: pm.Bundles, branch: str, sig: float) -> float:
    """``max |F(W(th, sig)) - W(th + omega, lam sig)|`` at midpoints (third-order defect)."""
    wh = pm.Whisker(c, b, branch)
    th = node_angles(c.n_nodes)[::8] + math.pi / c.n_nodes
    pts = np.array([wh.local(t, sig)[0] for t in th])
    fp, _ = pm.strob(pts, 1, c.mu, c.e)
    tgt = np.array([wh.local(t + c.omega, wh.lam * sig)[0] for t in th])
    return float(np.max(np.abs(fp - tgt)))


def stage_bundles(variant: str) -> None:
    out: dict[str, Any] = {}
    for name in ("3:4", "5:6"):
        c = load_circle(variant, name)
        log(f"{name}: bundles (N = {c.n_nodes})")
        b = pm.hyperbolic_bundles(c)
        row: dict[str, Any] = {
            "lam_u": b.lam_u,
            "lam_s": b.lam_s,
            "lam_u_times_lam_s": b.lam_u * b.lam_s,
            "offgrid_u": b.offgrid_u,
            "offgrid_s": b.offgrid_s,
            "tail_u": b.tail_u,
            "tail_s": b.tail_s,
            "v_u0": b.v_u[0],
            "v_s0": b.v_s[0],
            "w_u0": b.w_u[0],
            "reversor_v_u_vs_v_s": float(
                np.max(
                    np.abs(
                        pm.REVERSOR * fourier_eval(b.v_u, -node_angles(31))
                        - fourier_eval(b.v_s, node_angles(31))
                    )
                )
            ),
        }
        for branch in ("unstable", "stable"):
            for sig in (1e-3, 1e-4):
                row[f"order2_defect_{branch}_{sig:g}"] = _second_order_check(c, b, branch, sig)
        out[name] = row
        log(
            f"  lam_u {b.lam_u:.10f} lam_s {b.lam_s:.10f} prod {b.lam_u * b.lam_s:.3e} "
            f"offgrid u/s {b.offgrid_u:.1e}/{b.offgrid_s:.1e} v_u(0) "
            f"{np.array2string(b.v_u[0], precision=6)}"
        )
        for k, v in row.items():
            if k.startswith("order2"):
                log(f"  {k}: {v:.2e}")
        np.savez(
            _bundle_path(variant, name),
            lam_u=b.lam_u,
            lam_s=b.lam_s,
            v_u=b.v_u,
            v_s=b.v_s,
            w_u=b.w_u,
            w_s=b.w_s,
            offgrid_u=b.offgrid_u,
            offgrid_s=b.offgrid_s,
            tail_u=b.tail_u,
            tail_s=b.tail_s,
        )
    out["paper_fig9_lam_u_by_eye"] = "about 3.0 to 3.1 near omega_u = 1.558"
    out["paper_fig11_v_u0_x_by_eye"] = "about 0.1160 to 0.1162 near omega_u = 1.558"
    save(f"bundles_{variant}", out)


def whiskers(variant: str) -> tuple[pm.Whisker, pm.Whisker]:
    cu, cs = load_circle(variant, "3:4"), load_circle(variant, "5:6")
    bu, bs = load_bundles(variant, "3:4"), load_bundles(variant, "5:6")
    return pm.Whisker(cu, bu, "unstable"), pm.Whisker(cs, bs, "stable")


def stage_forward(variant: str) -> None:
    """Forward evaluation at the printed parameters, read with the labels AS PRINTED
    (theta_u, s_u on the 3:4 unstable manifold; theta_s, s_s on the 5:6 stable one),
    both signs, our theta origin. A diagnostic of conventions only."""
    wu, ws = whiskers(variant)
    x_star = pm.PRINTED_POINT
    out: dict[str, Any] = {"printed_point": x_star}
    for su in (pm.PRINTED_S_U, -pm.PRINTED_S_U):
        for ss in (pm.PRINTED_S_S, -pm.PRINTED_S_S):
            pu, _, _, _ = wu.evaluate(pm.PRINTED_THETA_U, su, wu.n_maps(su, 1e-5))
            ps, _, _, _ = ws.evaluate(pm.PRINTED_THETA_S, ss, ws.n_maps(ss, 1e-5))
            key = f"labels_as_printed su={su:+.5f} ss={ss:+.5f}"
            out[key] = {
                "Wu": pu,
                "Ws": ps,
                "dist_Wu": float(np.linalg.norm(pu - x_star)),
                "dist_Ws": float(np.linalg.norm(ps - x_star)),
            }
            log(f"{key}: |Wu-X*| {out[key]['dist_Wu']:.3e} |Ws-X*| {out[key]['dist_Ws']:.3e}")
    save(f"forward_{variant}", out)


def _closest_convention_free(
    wh: pm.Whisker, x_star: np.ndarray, ks: range, other: pm.Whisker
) -> dict[str, Any]:
    """Closest point of one manifold to X*, seeded only from X* itself (no printed
    parameter): map X* k periods towards the circle, read the local chart, globalise
    back, minimise; keep the best over k."""
    best: dict[str, Any] | None = None
    for k in ks:
        if wh.branch == "unstable":
            z, _ = pm.seed_from_point(wh, other, x_star, k, 1)
            th0, s0 = float(z[0]), float(z[1])
        else:
            z, _ = pm.seed_from_point(other, wh, x_star, 1, k)
            th0, s0 = float(z[2]), float(z[3])
        m = wh.n_maps(s0, 1e-5)
        th, s, w, d = pm.closest_on_whisker(wh, x_star, th0, s0, m)
        log(f"  {wh.branch} k={k}: seed ({th0:.4f}, {s0:.3f}) -> ({th:.6f}, {s:.5f}) d {d:.3e}")
        if best is None or d < best["distance"]:
            best = {"k": k, "theta": th, "s": s, "point": w, "distance": d, "m": m}
    assert best is not None
    return best


def stage_refine(variant: str) -> None:
    wu, ws = whiskers(variant)
    x_star = pm.PRINTED_POINT
    out: dict[str, Any] = {"printed_point": x_star}
    # A. convention-free: is X* on both manifolds of OUR model?
    log("A. closest points to the printed point, seeded from the point only")
    cu = _closest_convention_free(wu, x_star, range(5, 11), ws)
    cs = _closest_convention_free(ws, x_star, range(4, 10), wu)
    out["closest_on_Wu_34"] = cu
    out["closest_on_Ws_56"] = cs
    log(
        f"  W_u(3:4): (theta, s) = ({cu['theta']:.6f}, {cu['s']:.5f}) distance "
        f"{cu['distance']:.3e}; W_s(5:6): ({cs['theta']:.6f}, {cs['s']:.5f}) distance "
        f"{cs['distance']:.3e}"
    )
    # B. refine the connection from those parameters at two chart radii
    z0 = np.array([cu["theta"], cu["s"], cs["theta"], cs["s"]])
    conns: dict[float, pm.Connection] = {}
    for delta in (1e-5, 1e-6):
        log(f"B. refine, chart radius delta = {delta:g}")
        cn = pm.refine_connection(wu, ws, z0, delta=delta, verbose=True)
        conns[delta] = cn
        diff = cn.point - x_star
        out[f"delta={delta:g}"] = {
            "z": cn.z,
            "point": cn.point,
            "point_rounded_5dp": np.round(cn.point, 5),
            "residual": cn.residual,
            "m_u": cn.m_u,
            "m_s": cn.m_s,
            "singular_values_normalised_columns": cn.singular_values,
            "det_normalised_columns": float(np.prod(cn.singular_values)),
            "null_vector_param": cn.null_vector / np.linalg.norm(cn.null_vector),
            "point_minus_printed": diff,
            "max_abs_point_minus_printed": float(np.max(np.abs(diff))),
            "agrees_to_printed_5dp": bool(np.all(np.round(cn.point, 5) == x_star)),
            "history": cn.history,
        }
        log(
            f"  z {np.array2string(cn.z, precision=6)} residual {cn.residual:.2e} "
            f"sv {np.array2string(cn.singular_values, precision=3)}"
        )
        log(
            f"  point {np.array2string(cn.point, precision=8)} minus printed "
            f"{np.array2string(diff, precision=2)}"
        )
    a, b = conns[1e-5], conns[1e-6]
    out["junction_agreement_point"] = float(np.max(np.abs(a.point - b.point)))
    out["junction_agreement_z"] = float(np.max(np.abs(a.z - b.z)))
    log(
        f"junction agreement delta 1e-5 vs 1e-6: point {out['junction_agreement_point']:.2e} "
        f"z {out['junction_agreement_z']:.2e}"
    )
    # C. the valley: signed mismatch along the near-degenerate direction, theta_u fixed
    log("C. valley scan around the zero")
    rows = []
    for direction in (1.0, -1.0):
        ths = b.z[0] + direction * np.arange(0.0, 0.0301, 0.003)
        for v in pm.valley_scan(wu, ws, b.z, ths, m_u=b.m_u, m_s=b.m_s):
            rows.append(
                {
                    "theta_u": v.z[0],
                    "z": v.z,
                    "mismatch": v.mismatch,
                    "norm": v.norm,
                    "point_minus_printed": v.point - x_star,
                }
            )
    rows.sort(key=lambda r: r["theta_u"])
    out["valley"] = rows
    for r in rows:
        log(
            f"  theta_u {r['theta_u']:.5f} |f| {r['norm']:.2e} f_x {r['mismatch'][0]:+.2e} "
            f"max|P-X*| {np.max(np.abs(r['point_minus_printed'])):.2e}"
        )
    # D. re-seed displaced along the near-null direction: must return to the same zero
    nv = b.null_vector / np.linalg.norm(b.null_vector)
    for k, step in enumerate((3e-3, -3e-3)):
        cn = pm.refine_connection(wu, ws, b.z + step * nv, delta=1e-6)
        out[f"reseed_{k}"] = {
            "offset": step,
            "z": cn.z,
            "residual": cn.residual,
            "max_abs_z_minus_base": float(np.max(np.abs(cn.z - b.z))),
            "max_abs_point_minus_base": float(np.max(np.abs(cn.point - b.point))),
        }
        log(
            f"D. reseed {step:+g}: residual {cn.residual:.2e} |dz| "
            f"{out[f'reseed_{k}']['max_abs_z_minus_base']:.2e} |dP| "
            f"{out[f'reseed_{k}']['max_abs_point_minus_base']:.2e}"
        )
    # E. the printed parameter pairs against ours (labels as printed and swapped)
    v_norm_shift = float(np.linalg.norm(fourier_eval(wu.bundles.v_u, TWO_PI / 1024)))
    out["parameters"] = {
        "ours_34_unstable": [b.z[0], b.z[1]],
        "ours_56_stable": [b.z[2], b.z[3]],
        "printed_first_pair": [pm.PRINTED_THETA_U, pm.PRINTED_S_U],
        "printed_second_pair": [pm.PRINTED_THETA_S, pm.PRINTED_S_S],
        "ours_34_in_hypothesised_paper_convention": [b.z[0] - TWO_PI / 1024, b.z[1] * v_norm_shift],
        "norm_v_u_at_2pi_over_1024": v_norm_shift,
    }
    log(
        f"E. ours 3:4 W_u ({b.z[0]:.5f}, {b.z[1]:.5f}); ours 5:6 W_s ({b.z[2]:.5f}, "
        f"{b.z[3]:.5f}); printed ({pm.PRINTED_THETA_U}, {pm.PRINTED_S_U}), "
        f"({pm.PRINTED_THETA_S}, {pm.PRINTED_S_S}); ours 3:4 with theta origin shifted by "
        f"2 pi/1024 and s rescaled: ({b.z[0] - TWO_PI / 1024:.5f}, {b.z[1] * v_norm_shift:.5f})"
    )
    save(f"refine_{variant}", out)


def _approach_rows(
    circle: pm.Circle, bundles: pm.Bundles, start: np.ndarray, theta0: float, n: int, sign: int
) -> list[dict[str, float]]:
    """Map ``start`` ``n`` times (``sign`` = +1 forward, -1 backward) and record the
    symplectic hyperbolic coordinates about ``K(theta0 + sign * j * omega)``."""
    rows = []
    cur = start.copy()
    for j in range(n + 1):
        th = theta0 + sign * j * circle.omega
        c_s, c_u = pm.hyperbolic_coordinates(circle, bundles, th, cur)
        d, _ = pm.distance_to_circle(circle, cur)
        rows.append({"j": j, "c_s": c_s, "c_u": c_u, "dist": d})
        cur, _ = pm.strob(cur, sign, MU, E)
    return rows


def _longest_run(ratios: list[float], target: float, tol: float = 0.02) -> tuple[int, int]:
    """Longest run of consecutive ratios within ``tol`` (relative) of ``target``."""
    best = (0, -1)
    start = None
    for k, r in enumerate([*ratios, float("nan")]):
        ok = math.isfinite(r) and abs(r / target - 1.0) <= tol
        if ok and start is None:
            start = k
        if not ok and start is not None:
            if k - start > best[0]:
                best = (k - start, start)
            start = None
    return best


def stage_verify(variant: str) -> None:
    """#882 review criteria (Sec. 4 of the #882 results note) applied to the connection."""
    wu, ws = whiskers(variant)
    ref = load(f"refine_{variant}")
    r6 = ref["delta=1e-06"]
    z = np.array(r6["z"])
    x = np.array(r6["point"])
    m_u, m_s = int(r6["m_u"]), int(r6["m_s"])
    cu, cs = wu.circle, ws.circle
    bu, bs = wu.bundles, ws.bundles
    out: dict[str, Any] = {
        "z": z,
        "junction": x,
        "m_u": m_u,
        "m_s": m_s,
        "lam_s_56": bs.lam_s,
        "inv_lam_u_34": 1.0 / bu.lam_u,
    }
    n = m_u + 6
    legs = {
        # the single trajectory through the junction (unstable-side point), both ways
        "junction_forward_to_56": (_approach_rows(cs, bs, x, z[2], m_s + 6, 1), "c_s", bs.lam_s),
        "junction_backward_to_34": (_approach_rows(cu, bu, x, z[0], n, -1), "c_u", 1 / bu.lam_u),
    }
    ps, _, _, _ = ws.evaluate(z[2], z[3], m_s)
    legs["stable_leg_as_computed_forward_to_56"] = (
        _approach_rows(cs, bs, ps, z[2], m_s + 6, 1),
        "c_s",
        bs.lam_s,
    )
    for name, (rows, key, target) in legs.items():
        ratios = [rows[k + 1][key] / rows[k][key] for k in range(len(rows) - 1)]
        other = "c_u" if key == "c_s" else "c_s"
        run, at = _longest_run(ratios, target)
        jmin = int(np.argmin([r["dist"] for r in rows]))
        out[name] = {
            "rows": rows,
            "ratios": ratios,
            "target_ratio": target,
            "longest_run_within_2pct": run,
            "run_starts_at_j": at,
            "closest_j": jmin,
            "closest_distance": rows[jmin]["dist"],
            "other_over_this_at_closest": abs(rows[jmin][other] / rows[jmin][key]),
        }
        log(
            f"{name}: {run} consecutive periods within 2% of {target:.5f} (from j={at}); "
            f"closest {rows[jmin]['dist']:.2e} at j={jmin}, |{other}/{key}| there "
            f"{out[name]['other_over_this_at_closest']:.1e}"
        )
        for k, r in enumerate(ratios):
            log(f"   j {k:2d}: ratio {r:+.5f}  dist {rows[k]['dist']:.2e}")
    # random 1e-6 perturbations of the junction must fail to reach either circle
    rng = np.random.default_rng(889)
    pert = rng.normal(size=(8, 4))
    pert = 1e-6 * pert / np.linalg.norm(pert, axis=1)[:, None]
    jf = out["junction_forward_to_56"]["closest_j"]
    jb = out["junction_backward_to_34"]["closest_j"]
    af, _ = pm.strob(x[None, :] + pert, jf, MU, E)
    ab, _ = pm.strob(x[None, :] + pert, -jb, MU, E)
    df = [pm.distance_to_circle(cs, a)[0] for a in af]
    db = [pm.distance_to_circle(cu, a)[0] for a in ab]
    out["perturbation_1e-6"] = {
        "forward_maps": jf,
        "backward_maps": jb,
        "forward_dist_true": out["junction_forward_to_56"]["closest_distance"],
        "backward_dist_true": out["junction_backward_to_34"]["closest_distance"],
        "forward_dist_perturbed": df,
        "backward_dist_perturbed": db,
    }
    log(
        f"perturbation 1e-6: forward {jf} maps true "
        f"{out['junction_forward_to_56']['closest_distance']:.2e} perturbed min {min(df):.2e}; "
        f"backward {jb} maps true {out['junction_backward_to_34']['closest_distance']:.2e} "
        f"perturbed min {min(db):.2e}"
    )
    out["junction_agreement_delta_1e-5_vs_1e-6"] = ref["junction_agreement_point"]
    save(f"verify_{variant}", out)


STAGES = {
    "po": stage_po,
    "tori": stage_tori,
    "bundles": stage_bundles,
    "forward": stage_forward,
    "refine": stage_refine,
    "verify": stage_verify,
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=sorted(STAGES))
    ap.add_argument("--variant", default="table1", choices=["table1", "printed", "table1_hi"])
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    pm.set_threads(args.threads)
    log(f"#889 stage {args.stage} variant {args.variant}")
    STAGES[args.stage](args.variant)
    log("done")


if __name__ == "__main__":
    main()
