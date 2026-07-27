"""#729 -- epoch(+torus-point) robustness scan for `#726`'s N=5 CRNBP torus
real-ephemeris consistency check (the `#705` precedent, N=5 analogue).

`#726` found that the idealized (circular-coplanar) N=5 CRNBP Jupiter-Io-
Europa-Ganymede torus (`#720`/`#723`/`#724`'s confirmed narrow-novelty object)
does NOT survive real Europa/Io/Ganymede ephemeris at ONE tested torus point
(theta1=theta2=0) and ONE tested epoch (2030-01-01, one Ganymede-synodic
forcing period): a real, VERIFIED (noise-floor-checked, not integrator
artifact) generic collapse, ``pos_gap_km=3.68e5`` -- comparable to Ganymede's
own orbital radius. A small, explicitly non-systematic sample of 2 other
torus points and 2x/4x window multiples at that SAME single epoch showed
similarly large gaps.

This EXACTLY mirrors `#704`'s own first-pass result for the Uranus Umbriel-
Titania case (median mismatch ~125,000-136,000 km across a synodic period at
its own single tested epoch) -- which `#705` then showed was NOT
representative: a dense per-epoch synodic-period scan, repeated at 10
well-separated epochs across 2000-2083, found a comparably tight near-miss
window recurring at ALL 10 epochs. This script runs the N=5 analogue of that
same scan, to settle whether `#726`'s single-epoch collapse is similarly
non-representative, or a genuine, consistent negative for this object.

Reuses `#726`'s own top-level driver
(:func:`cyclerfinder.search.crnbp_real_ephemeris_consistency.check_torus_survives_real_ephemeris`)
UNMODIFIED as a library, exactly as `#705` reused `#704`'s own
``check_connection_survives_real_ephemeris`` -- no changes to the underlying
consistency-check machinery.

Epoch axis (mirrors `#705`'s own convention exactly)
-----------------------------------------------------
10 epochs, evenly spaced across ``[2000, 2083]`` (``frac = i/9``), same
~9.22-year spacing `#705` used. That window was originally chosen for the
Uranus case as a `#312`-precedented real-ephemeris validity window; for
Jupiter, ``jup365.bsp``'s own SPICE coverage is actually far broader
(1600-2200, checked directly via ``spice.spkcov`` during this task's own
build), so nothing here is forced narrower by kernel coverage -- the SAME
2000-2083/10-epoch convention is kept anyway, for direct apples-to-apples
comparability with the `#705` precedent rather than because the physics
requires it.

At each epoch, a dense ``N_SYNODIC_PRIMARY``-point scan across one real
Ganymede-synodic forcing period (the SAME quantity `torus.period` already
uses as `#726`'s own default ``t_window_tu``) is run, then the best
(smallest ``pos_gap_km``) point is refined with the SAME 28-iteration,
3-point-sampling bisection `#705`'s own script used, to locate the true
local minimum within that epoch's own period (not just the coarse-grid best
sample).

Torus-point axis (adaptation beyond `#705`'s own scope, justified below)
-------------------------------------------------------------------------
`#705`'s object (a manifold connection) had no free "torus point" parameter
-- one connection, one departure/target pair. `#726`'s object (a torus) DOES
have one (``theta1``, ``theta2``), and `#726`'s own report flagged that a
small, non-systematic sample of other torus points showed similarly large
gaps -- an open question `#705`'s own precedent never had to answer.
Whether to scan it densely, per the dispatch's own instruction to use
judgment: measured directly during this task's own build, one
``check_torus_survives_real_ephemeris`` call costs ~0.065s (real N-body
propagation over one ~7-day window, single spacecraft-mass body) -- #705's
own single-axis choice was NOT a methodological necessity but a cost
tradeoff specific to its own (irrelevantly larger, from a different scan
axis) budget; that tradeoff does not bind here. This script therefore adds a
genuinely systematic SECOND axis: the SAME per-epoch dense-scan-plus-
bisection methodology, repeated at 4 additional, symmetrically chosen torus
points (``{0, pi} x {0, pi}`` plus the diagonal midpoint ``(pi/2, pi/2)``)
across the SAME 10 epochs -- a reduced ``N_SYNODIC_SECONDARY`` per point
(60 vs. the primary axis's 300) keeps the total added cost to a few extra
minutes rather than 4x the primary axis's own runtime, while still giving
a genuinely systematic (not ad hoc) picture of torus-point dependence.

"Narrow near-miss" bar
------------------------
`#726` (unlike `#704`) has no own committed narrow-near-miss reference value
to scale a "comparable" threshold against -- it found only the generic
collapse. This script instead borrows `#705`'s OWN established bar: `#705`'s
own dense per-epoch scan found narrow windows characterized as "sub-5000km"
(directly inherited from `#704`'s own original characterization, "3 distinct
sub-5000km-mismatch windows per synodic period"), and its own worst-case
per-epoch local minimum across all 10 tested epochs was 142.82 km -- both
comfortably under 5,000 km. ``NARROW_NEAR_MISS_KM = 5000.0`` is used here as
the same, sourced (not invented) bar. Full duty-cycle counts at the SAME
absolute threshold list `#705` used (500 km through 50,000 km) are reported
regardless, so the reader can judge independently of this one chosen bar.

No catalogue writeback; no schema design (out of scope per this task's own
dispatch, a separate, not-yet-registered future task gated on a clearly
positive verdict here).

Run:  uv run python scripts/run_729_epoch_torus_robustness_scan.py
Outputs -> data/found/729_crnbp_epoch_torus_robustness_scan/result.json
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.ccr4bp as ccr4bp  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.core.crnbp as crnbp  # noqa: E402
import cyclerfinder.search.crnbp_real_ephemeris_consistency as rec5  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
import cyclerfinder.search.variational_crnbp_torus as vc  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402
from cyclerfinder.search.ccr4bp_real_ephemeris_consistency import tu_to_seconds  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "729_crnbp_epoch_torus_robustness_scan"

N_EPOCHS = 10
START_YEAR, END_YEAR = 2000, 2083
N_SYNODIC_PRIMARY = 300
N_SYNODIC_SECONDARY = 60
DUTY_CYCLE_THRESHOLDS_KM = [500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0]

# Sourced from #705's own established "narrow near-miss" scale (see module
# docstring) -- NOT re-derived here, since #726 has no own committed
# near-miss reference the way #704 did for #705.
NARROW_NEAR_MISS_KM = 5000.0

HEADLINE_TORUS_POINT: tuple[float, float] = (0.0, 0.0)  # #726's own headline point
SECONDARY_TORUS_POINTS: list[tuple[float, float]] = [
    (np.pi, 0.0),
    (0.0, np.pi),
    (np.pi, np.pi),
    (np.pi / 2.0, np.pi / 2.0),
]


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical scaffolding to `#724`'s/`#726`'s own reproduction scripts."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


def _log(msg: str, t0: float) -> None:
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


def rebuild_724_final_torus(t0: float) -> vc.CRNBPTorusVariationalResult:
    """Reproduce `#720`/`#723`/`#724`'s own delivered N=5 torus bit-for-bit --
    identical pipeline to `#726`'s own test module's
    ``_rebuild_724_final_torus`` (``tests/search/test_crnbp_real_ephemeris_consistency.py``),
    duplicated here (not imported from a test module) exactly as `#705`'s own
    script duplicated `#701`'s reconstruction scaffolding rather than
    importing it from a test file.
    """
    system4 = ccr4bp.jupiter_europa_ganymede_default()
    target = crnbp.jupiter_europa_io_ganymede_default()
    s0, period, res = _resonant_symmetric_orbit(system4.mu, 3, 4)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system4,
        s0,
        period,
        n1=2,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    _log(f"phys (N=4) torus rebuilt: residual_rms={phys.residual_rms:.3e}", t0)
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=target.perturbers[0].a,
        omega_io=target.perturbers[0].omega,
        theta_io0=target.perturbers[0].theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    _log(f"N=5 seed (mu_io=0) built: residual_rms={seed.residual_rms:.3e}", t0)
    steps = vc.continue_crnbp_torus_mu_io(
        seed,
        target,
        n_steps=8,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    final = steps[-1]
    _log(
        f"final N=5 torus reconstructed: residual_rms={final.residual_rms:.3e}, "
        f"period={final.period:.6f} TU",
        t0,
    )
    return final


def _duty_cycle(pos_arr: np.ndarray) -> dict[str, Any]:
    """Identical counting logic to `#705`'s own ``_duty_cycle``
    (``scripts/run_705_epoch_robustness_scan.py``): contiguous runs of
    below-threshold samples, non-circular (wrap-adjacent runs counted
    separately). Duplicated (not imported) for import-independence between
    sibling scan scripts -- see this module's own docstring."""
    duty_cycle: dict[str, Any] = {}
    for thresh in DUTY_CYCLE_THRESHOLDS_KM:
        below = pos_arr < thresh
        frac = float(below.mean())
        n_windows = 0
        in_window = False
        for b in below:
            if b and not in_window:
                n_windows += 1
                in_window = True
            elif not b:
                in_window = False
        duty_cycle[f"below_{int(thresh)}km"] = {
            "fraction_of_period": frac,
            "n_distinct_windows": n_windows,
        }
    return duty_cycle


def _slim_result(phase_frac: float, r: Any) -> dict[str, Any]:
    return {
        "phase_frac": phase_frac,
        "epoch_utc": r.epoch0_utc,
        "pos_gap_km": r.pos_gap_km,
        "vel_gap_km_s": r.vel_gap_km_s,
        "propagation_success": r.propagation_success,
        "closest_approach_km": r.closest_approach_km,
    }


def _full_result(r: Any) -> dict[str, Any]:
    return {
        "epoch0_utc": r.epoch0_utc,
        "et0": r.et0,
        "t_u_tu": r.t_u_tu,
        "t_u_seconds": r.t_u_seconds,
        "r0_km": r.r0_km.tolist(),
        "v0_km_s": r.v0_km_s.tolist(),
        "r_f_km": r.r_f_km.tolist(),
        "v_f_km_s": r.v_f_km_s.tolist(),
        "r_target_km": r.r_target_km.tolist(),
        "v_target_km_s": r.v_target_km_s.tolist(),
        "pos_gap_km": r.pos_gap_km,
        "vel_gap_km_s": r.vel_gap_km_s,
        "propagation_success": r.propagation_success,
        "closest_approach_km": r.closest_approach_km,
        "notes": r.notes,
    }


def scan_one_epoch_one_torus_point(
    torus: vc.CRNBPTorusVariationalResult,
    theta1_0: float,
    theta2_0: float,
    base_dt: _dt.datetime,
    synodic_period_days: float,
    n_synodic: int,
    t0: float,
    label: str,
) -> dict[str, Any]:
    """Dense scan across one real Ganymede-synodic period starting at
    ``base_dt``, at fixed torus point ``(theta1_0, theta2_0)``, followed by a
    28-iteration bisection refine of the best (smallest ``pos_gap_km``)
    sample -- mirrors `#705`'s own ``_scan_one_epoch`` exactly, adapted to
    call `#726`'s own torus-based driver directly (no departure/target
    precompute needed -- the torus + theta point IS the departure
    specification `check_torus_survives_real_ephemeris` already takes)."""
    synodic_scan: list[dict[str, Any]] = []

    def _at_offset(offset_days: float) -> Any:
        epoch_dt = base_dt + _dt.timedelta(days=offset_days)
        epoch_utc = epoch_dt.strftime("%Y-%m-%dT%H:%M:%S")
        return rec5.check_torus_survives_real_ephemeris(epoch_utc, torus, theta1_0, theta2_0)

    for i in range(n_synodic):
        frac = i / n_synodic
        offset_days = frac * synodic_period_days
        result = _at_offset(offset_days)
        synodic_scan.append(_slim_result(frac, result))

    pos_arr = np.array([e["pos_gap_km"] for e in synodic_scan])
    duty_cycle = _duty_cycle(pos_arr)

    best_idx = int(np.argmin(pos_arr))
    best_offset_days = synodic_scan[best_idx]["phase_frac"] * synodic_period_days
    step_days = synodic_period_days / n_synodic

    lo, hi = best_offset_days - step_days, best_offset_days + step_days
    best_offset = best_offset_days
    best_result = _at_offset(best_offset_days)
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        g_lo = _at_offset(lo).pos_gap_km
        r_mid = _at_offset(mid)
        g_mid = r_mid.pos_gap_km
        g_hi = _at_offset(hi).pos_gap_km
        candidates = [(g_lo, lo), (g_mid, mid), (g_hi, hi)]
        candidates.sort(key=lambda t: t[0])
        best_offset = candidates[0][1]
        if candidates[0][1] == mid:
            best_result = r_mid
        width = hi - lo
        lo, hi = best_offset - width / 4.0, best_offset + width / 4.0
    best_result = _at_offset(best_offset)

    narrow = best_result.pos_gap_km < NARROW_NEAR_MISS_KM and best_result.propagation_success
    _log(
        f"{label} base={base_dt.strftime('%Y-%m-%d')}: "
        f"synodic min={pos_arr.min():.1f} km, median={float(np.median(pos_arr)):.1f} km, "
        f"local-min pos_gap={best_result.pos_gap_km:.3f} km "
        f"vel_gap={best_result.vel_gap_km_s * 1000:.4f} m/s "
        f"[{'NARROW near-miss' if narrow else 'not narrow (collapse-scale)'}]",
        t0,
    )

    return {
        "label": label,
        "theta1_0": theta1_0,
        "theta2_0": theta2_0,
        "base_epoch_utc": base_dt.strftime("%Y-%m-%dT%H:%M:%S"),
        "n_synodic_points": n_synodic,
        "synodic_scan_pos_gap_km_min": float(pos_arr.min()),
        "synodic_scan_pos_gap_km_max": float(pos_arr.max()),
        "synodic_scan_pos_gap_km_median": float(np.median(pos_arr)),
        "duty_cycle_dense_synodic_scan": duty_cycle,
        "local_minimum": _full_result(best_result),
        "narrow_near_miss": narrow,
        "synodic_phase_scan": synodic_scan,
    }


def main() -> dict[str, Any]:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    torus = rebuild_724_final_torus(t0)
    l_km = rec5.L_KM
    v_unit_km_s = rec5.v_unit_km_s()
    synodic_period_days = torus.period * tu_to_seconds(l_km, v_unit_km_s) / 86400.0
    _log(f"real Ganymede-synodic forcing period = {synodic_period_days:.4f} days", t0)
    _log(f"'narrow near-miss' bar (sourced from #705's own scale): {NARROW_NEAR_MISS_KM} km", t0)

    epoch_years = [
        START_YEAR + (i / (N_EPOCHS - 1)) * (END_YEAR - START_YEAR) for i in range(N_EPOCHS)
    ]
    epoch_base_dts = [
        _dt.datetime(int(year), 1, 1, tzinfo=_dt.UTC)
        + _dt.timedelta(days=(year - int(year)) * 365.25)
        for year in epoch_years
    ]

    # --- Primary axis: headline torus point, all 10 epochs, full resolution --- #
    primary_scans: list[dict[str, Any]] = []
    for i, base_dt in enumerate(epoch_base_dts):
        scan = scan_one_epoch_one_torus_point(
            torus,
            *HEADLINE_TORUS_POINT,
            base_dt,
            synodic_period_days,
            N_SYNODIC_PRIMARY,
            t0,
            label=f"primary epoch {i + 1}/{N_EPOCHS}",
        )
        scan["epoch_year"] = epoch_years[i]
        primary_scans.append(scan)

    # --- Secondary axis: 4 additional torus points, all 10 epochs, reduced resolution --- #
    secondary_scans: list[dict[str, Any]] = []
    for tp_idx, (th1, th2) in enumerate(SECONDARY_TORUS_POINTS):
        for i, base_dt in enumerate(epoch_base_dts):
            scan = scan_one_epoch_one_torus_point(
                torus,
                th1,
                th2,
                base_dt,
                synodic_period_days,
                N_SYNODIC_SECONDARY,
                t0,
                label=f"secondary torus_pt {tp_idx + 1}/{len(SECONDARY_TORUS_POINTS)} "
                f"epoch {i + 1}/{N_EPOCHS}",
            )
            scan["epoch_year"] = epoch_years[i]
            secondary_scans.append(scan)

    n_narrow_primary = sum(1 for e in primary_scans if e["narrow_near_miss"])
    n_narrow_secondary = sum(1 for e in secondary_scans if e["narrow_near_miss"])
    primary_best = [e["local_minimum"]["pos_gap_km"] for e in primary_scans]
    secondary_best = [e["local_minimum"]["pos_gap_km"] for e in secondary_scans]

    summary = {
        "narrow_near_miss_threshold_km": NARROW_NEAR_MISS_KM,
        "n_epochs_tested": N_EPOCHS,
        "primary_torus_point": HEADLINE_TORUS_POINT,
        "n_primary_synodic_points_per_epoch": N_SYNODIC_PRIMARY,
        "n_primary_epochs_narrow": n_narrow_primary,
        "fraction_primary_epochs_narrow": n_narrow_primary / N_EPOCHS,
        "primary_per_epoch_local_minimum_pos_gap_km": primary_best,
        "primary_per_epoch_year": epoch_years,
        "primary_overall_best_pos_gap_km": min(primary_best),
        "primary_overall_worst_pos_gap_km": max(primary_best),
        "secondary_torus_points": SECONDARY_TORUS_POINTS,
        "n_secondary_synodic_points_per_epoch": N_SYNODIC_SECONDARY,
        "n_secondary_epoch_torus_combos_tested": len(secondary_scans),
        "n_secondary_combos_narrow": n_narrow_secondary,
        "fraction_secondary_combos_narrow": (
            n_narrow_secondary / len(secondary_scans) if secondary_scans else float("nan")
        ),
        "secondary_overall_best_pos_gap_km": (
            min(secondary_best) if secondary_best else float("nan")
        ),
        "secondary_overall_worst_pos_gap_km": (
            max(secondary_best) if secondary_best else float("nan")
        ),
    }
    _log(f"summary: {summary}", t0)

    result_all = {
        "task": "#729",
        "system": (
            "Jupiter-Europa-Io-Ganymede N=5 CRNBP epoch(+torus-point) robustness scan "
            "of #726's generic-collapse finding"
        ),
        "source_torus": "#720/#723/#724's final continued N=5 torus (mu_io physical)",
        "torus_period_tu": torus.period,
        "torus_residual_rms": torus.residual_rms,
        "l_km_europa": l_km,
        "v_unit_europa_km_s": v_unit_km_s,
        "real_ganymede_synodic_period_days": synodic_period_days,
        "primary_scans": primary_scans,
        "secondary_scans": secondary_scans,
        "summary": summary,
        "force_model_notes": (
            "Identical to #726's own force model: Jupiter point-mass (system-GM) "
            "central term + Europa + Io + Ganymede REAL SPICE (jup365.bsp) "
            "third-body perturbations, Jupiter-centred J2000 inertial frame. "
            "See src/cyclerfinder/search/crnbp_real_ephemeris_consistency.py "
            "module docstring for the full central-GM-convention justification."
        ),
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result_all, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result_all


if __name__ == "__main__":
    main()
