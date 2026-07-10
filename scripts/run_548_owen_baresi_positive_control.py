#!/usr/bin/env python3
"""Task #548: Owen & Baresi positive-control gate for the qp_tori / qp_torus_heteroclinic
linking-number heteroclinic-detection pipeline.

Reframed per the #548 OUTSTANDING entry: NOT a reproduction of Owen & Baresi's exact
"latitudinal frequency" numbers (correctly abandoned as impractical by #534), but a test of
whether the linking-number scan sees ANY sign changes across a BAND of isoenergetic Earth-Moon
L1 x L2 quasi-halo torus pairs at a common energy below both collinear necks -- the regime where
Owen & Baresi (Astrodynamics 8, 2024) demonstrate 4 connections.

Two departures from #534, both forced by primary evidence found this session:
  1. ENERGY. #534's committed NRHO seeds sit at C=3.045 (verified: cr3bp.jacobi_constant of the
     corrected seeds), NOT the C=3.15 that #547 recorded as "already correct". Owen & Baresi's
     demo is at C=3.15, but exactly-3.15 isoenergetic L1+L2 quasi-halo tori are impractical here:
     the EM L1 halo family bifurcates from planar Lyapunov at C~3.146 (so C=3.15 is at/above the
     L1 quasi-halo regime) and the L2 halo family reached via NRHO continuation tops out ~C=3.087.
     Heteroclinic connections are ISOENERGETIC, so BOTH tori must share one C. We therefore build
     the positive control at the highest COMMON energy both families robustly reach (C in
     [3.05, 3.087]) -- the same physics (both necks open, unstable quasi-halo pairs), a slightly
     deeper energy than the paper's single 3.15 demo. Documented honestly; the adjudicator weighs
     the caveat.
  2. BRANCH CLASSIFICATION. Uses #547/#548's EMPIRICAL transit-branch grid builder
     (genome.qp_torus_transit.transit_torus_manifold_grid) instead of qp_torus_manifold's
     untested vec[0]*sign eigenvector heuristic -- at every torus point it propagates BOTH signed
     perturbations and keeps whichever reaches the surface of section first.

Pre-registered kill criterion (#548): zero sign changes across the whole sweep -> the
linking-number pipeline is SHELVED (postmortem + downgrade of #534/#536/#546).
"""

from __future__ import annotations

import json
import math
import pathlib
import time

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus, evaluate_torus
from cyclerfinder.genome.qp_torus_heteroclinic import (
    scan_linking_number,
)
from cyclerfinder.genome.qp_torus_transit import transit_torus_manifold_grid
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.nrho_continuation import correct_symmetric_nrho

MU = 0.012153643
SYS = cr3bp.CR3BPSystem(mu=MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0)
SURFACE_X = 1.0 - MU
LOG = pathlib.Path(
    "/tmp/claude-1000/-home-bruce-dev-cyclers/"
    "e8a086b8-fae2-4e77-b340-1425b9d3c532/scratchpad/run548_results.json"
)

# Converged halo-family samples (x0, z0, ydot0, T, C) -- seeds for the parent-C bisection.
L1_TABLE = [
    (0.836314, -0.14569, 0.25490, 2.7513, 3.04504),
    (0.83331, -0.13270, 0.24559, 2.7759, 3.06198),
    (0.83031, -0.11615, 0.23135, 2.7870, 3.08325),
    (0.82731, -0.09474, 0.20992, 2.7839, 3.10944),
    (0.82431, -0.06058, 0.17154, 2.7647, 3.14523),
]
L2_TABLE = [
    (1.023731, 0.18327, -0.10696, 1.5336, 3.04483),
    (1.01173, 0.17383, -0.07970, 1.3730, 3.05816),
    (1.00473, 0.16604, -0.06220, 1.2688, 3.06875),
    (0.99703, 0.15253, -0.04099, 1.1180, 3.08719),
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parent_halo_at_c(table: list[tuple[float, ...]], target_c: float):
    """Return a converged SymmetricNRHO halo at (approximately) target Jacobi via
    secant on x0, seeded by linear interpolation of the nearest table rows."""
    # nearest two rows by C
    rows = sorted(table, key=lambda r: abs(r[4] - target_c))[:2]
    (x0a, z0a, yda, ta, ca), (x0b, z0b, ydb, tb, cb) = rows
    if abs(cb - ca) < 1e-9:
        x0 = x0a
        z0g, ydg, tg = z0a, yda, ta
    else:
        w = (target_c - ca) / (cb - ca)
        x0 = x0a + w * (x0b - x0a)
        z0g = z0a + w * (z0b - z0a)
        ydg = yda + w * (ydb - yda)
        tg = ta + w * (tb - ta)
    # secant on x0 to hit target C
    prev = None
    for _ in range(12):
        r = correct_symmetric_nrho(SYS, float(x0), z0g, ydg, tg, with_monodromy=False)
        if not r.converged or abs(r.z0) < 1e-4:
            return None
        z0g, ydg, tg = r.z0, r.ydot0, r.T_TU
        err = r.jacobi - target_c
        if abs(err) < 2e-4:
            return r
        if prev is None:
            x0 += -math.copysign(0.001, err) * (1.0 if table is L1_TABLE else -1.0)
        else:
            px0, perr = prev
            denom = err - perr
            if abs(denom) < 1e-12:
                break
            x0 = x0 - err * (x0 - px0) / denom
        prev = (float(r.x0), err)
    return r


def build_torus(halo, amp: float, side: str) -> QPTorus | None:
    state0 = np.array([halo.x0, 0.0, halo.z0, 0.0, halo.ydot0, 0.0])
    mono = monodromy(SYS, state0, halo.T_TU)
    eigs = floquet_multipliers(mono)
    cands = [e for e in eigs if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    if len(cands) < 2:
        return None
    # Pick the complex pair nearest the unit circle (the center / quasi-halo mode)
    # and choose k so a primitive k-th root of unity best matches its phase.
    lam = min(cands, key=lambda e: abs(abs(e) - 1.0))
    phi = abs(math.atan2(lam.imag, lam.real))
    best_k, best_d = 5, math.inf
    for kk in range(3, 21):
        for j in range(1, kk):
            if math.gcd(j, kk) != 1:
                continue
            d = abs(phi - 2.0 * math.pi * j / kk)
            if d < best_d:
                best_d, best_k = d, kk
    k = best_k
    try:
        torus = correct_qp_torus(
            SYS,
            state0,
            halo.T_TU,
            (cands[0], cands[1]),
            k=int(k),
            n_trans=4,
            initial_torus_amplitude=amp,
            tol=1e-7,
            max_iter=40,
        )
    except (ValueError, RuntimeError) as e:
        log(f"    torus build failed ({side} amp={amp}): {e}")
        return None
    return torus


def torus_c(torus: QPTorus) -> float:
    return float(cr3bp.jacobi_constant(evaluate_torus(torus, 0.0, 0.0), MU))


SCAN_SPECS = [
    ("z", ("y", "ydot", "zdot")),
    ("z", ("x", "y", "ydot")),
    ("y", ("z", "ydot", "zdot")),
    ("zdot", ("y", "z", "ydot")),
]


def scan_pair(u_grid, s_grid, comp_idx: int, tag: str, results: list) -> int:
    """Run all scan specs on one (unstable=L1, stable=L2) grid pair; return #sign changes."""
    total = 0
    for scan_comp, curve_comps in SCAN_SPECS:
        ci = {"x": 0, "y": 1, "z": 2, "xdot": 3, "ydot": 4, "zdot": 5}[scan_comp]
        fu = u_grid.endpoints[:, :, ci]
        fs = s_grid.endpoints[:, :, ci]
        fu_f = fu[np.isfinite(fu)]
        fs_f = fs[np.isfinite(fs)]
        if fu_f.size < 4 or fs_f.size < 4:
            continue
        lo = max(fu_f.min(), fs_f.min())
        hi = min(fu_f.max(), fs_f.max())
        if not (lo < hi):
            continue
        dvals = np.linspace(lo, hi, 60)
        res = scan_linking_number(
            s_grid,
            u_grid,
            scanning_component=scan_comp,
            curve_components=curve_comps,
            d_values=dvals,
        )
        changes = res.sign_change_locations()
        n_nonzero = int(np.sum(res.linking_numbers != 0))
        entry = {
            "tag": tag,
            "scan_comp": scan_comp,
            "curve_comps": curve_comps,
            "overlap": [float(lo), float(hi)],
            "n_nonzero_lk": n_nonzero,
            "linking_numbers": res.linking_numbers.tolist(),
            "sign_changes": changes,
        }
        results.append(entry)
        if changes:
            log(f"    *** SIGN CHANGE {tag} scan={scan_comp} curves={curve_comps}: {changes}")
            log(f"        linking numbers: {res.linking_numbers.tolist()}")
        total += len(changes)
    return total


def main() -> None:
    log("=== #548 Owen-Baresi positive-control sweep ===")
    target_cs = [3.05, 3.06, 3.07, 3.08]
    # NOTE: torus frequency ratio is set by the parent halo (energy), NOT by the seed
    # amplitude (confirmed in-run: L1 ratio 0.3668 identical at amp 5e-4 and 1.2e-3). The
    # ENERGY sweep is what brackets the published L1 0.2739 / L2 0.02163 ratio; a single
    # representative amplitude per side therefore suffices, and keeps the run tractable.
    amps = [5e-4]
    n_long, n_lat = 20, 16
    t_max_u, t_max_s = 18.0, 24.0
    all_results: list = []
    grand_total_changes = 0
    geometry_ok_pairs = 0

    for tc in target_cs:
        log(f"--- target common C = {tc} ---")
        h1 = parent_halo_at_c(L1_TABLE, tc)
        h2 = parent_halo_at_c(L2_TABLE, tc)
        if h1 is None or h2 is None:
            log(f"  parent halo unavailable (L1={h1 is not None}, L2={h2 is not None}); skip")
            continue
        log(
            f"  L1 halo C={h1.jacobi:.5f} x0={h1.x0:.5f} z0={h1.z0:.5f}; "
            f"L2 halo C={h2.jacobi:.5f} x0={h2.x0:.5f} z0={h2.z0:.5f}"
        )

        l1_tori = []
        for amp in amps:
            t = build_torus(h1, amp, "L1")
            if t is not None and t.invariance_residual < 1e-4:
                l1_tori.append((amp, t))
                log(
                    f"  L1 torus amp={amp}: C={torus_c(t):.5f} invres={t.invariance_residual:.2e} "
                    f"omega_long={t.omega_long:.4f} omega_trans={t.omega_trans:.4f} "
                    f"ratio={abs(t.omega_trans / t.omega_long):.4f}"
                )
        l2_tori = []
        for amp in amps:
            t = build_torus(h2, amp, "L2")
            if t is not None and t.invariance_residual < 1e-4:
                l2_tori.append((amp, t))
                log(
                    f"  L2 torus amp={amp}: C={torus_c(t):.5f} invres={t.invariance_residual:.2e} "
                    f"omega_long={t.omega_long:.4f} omega_trans={t.omega_trans:.4f} "
                    f"ratio={abs(t.omega_trans / t.omega_long):.4f}"
                )

        if not l1_tori or not l2_tori:
            log("  no converged tori this energy; skip")
            continue

        # Build empirical-transit grids once per torus.
        l1_grids = []
        for amp, t in l1_tori:
            g, signs = transit_torus_manifold_grid(
                t,
                n_long=n_long,
                n_lat=n_lat,
                branch="unstable",
                surface_x=SURFACE_X,
                t_max=t_max_u,
                eps=1e-5,
            )
            ncross = int(np.sum(np.isfinite(g.endpoints[:, :, 0])))
            sign_bal = (int(np.sum(signs == 1)), int(np.sum(signs == -1)))
            log(
                f"  L1 unstable grid amp={amp}: crossings={ncross}/{n_long * n_lat} "
                f"sign(+/-)={sign_bal}"
            )
            l1_grids.append((amp, g, ncross))
        l2_grids = []
        for amp, t in l2_tori:
            g, signs = transit_torus_manifold_grid(
                t,
                n_long=n_long,
                n_lat=n_lat,
                branch="stable",
                surface_x=SURFACE_X,
                t_max=t_max_s,
                eps=1e-5,
            )
            ncross = int(np.sum(np.isfinite(g.endpoints[:, :, 0])))
            sign_bal = (int(np.sum(signs == 1)), int(np.sum(signs == -1)))
            log(
                f"  L2 stable grid amp={amp}: crossings={ncross}/{n_long * n_lat} "
                f"sign(+/-)={sign_bal}"
            )
            l2_grids.append((amp, g, ncross))

        for a1, gu, ncu in l1_grids:
            for a2, gs, ncs in l2_grids:
                if ncu < 4 or ncs < 4:
                    continue
                geometry_ok_pairs += 1
                tag = f"C{tc}_L1a{a1}_L2a{a2}"
                grand_total_changes += scan_pair(gu, gs, 0, tag, all_results)

    log(
        f"=== SWEEP DONE: geometry-usable pairs={geometry_ok_pairs}, "
        f"TOTAL sign changes across whole sweep = {grand_total_changes} ==="
    )
    LOG.write_text(
        json.dumps(
            {
                "grand_total_changes": grand_total_changes,
                "geometry_ok_pairs": geometry_ok_pairs,
                "results": all_results,
            },
            indent=1,
        )
    )
    log(f"results written to {LOG}")
    if grand_total_changes == 0:
        log("KILL CRITERION MET: zero sign changes across the swept band.")
    else:
        log("SIGN CHANGES FOUND: candidate connection(s) -- refine with deflated_newton.")


if __name__ == "__main__":
    main()
