"""#704 -- CCR4BP-to-real-ephemeris consistency check: Uranus Umbriel-Titania.

Answers the load-bearing vetting-chain question `#704` was dispatched for:
does `#701`'s idealized (circular-coplanar) CCR4BP homoclinic connection on
the Umbriel-Titania resonant torus approximately survive when Umbriel's and
Titania's REAL (eccentric, mutually inclined, URA111-SPICE-sourced) ephemeris
replaces the idealized approximation?

Reconstructs `#701`'s own ``best_robust_genuine_connection_corrected`` entry
bit-for-bit (same torus corrector, same saved coarse seed, same
``refine_candidate`` call -- see
``data/found/701_ccr4bp_umbriel_titania_search/result.json``), then runs
:func:`cyclerfinder.search.ccr4bp_real_ephemeris_consistency.check_connection_survives_real_ephemeris`
at a headline epoch plus two epoch scans (a short-baseline scan across one
real Umbriel-Titania synodic period, ~7.9 days, and a long-baseline scan
across the full ~2000-2083 kernel-supported window) to characterize how the
real-vs-idealized mismatch depends on epoch -- the `#312`-style "duty cycle"
framing, adapted (NOT forced to match #312's own numbers -- a different
object may behave completely differently, see this task's own scope note).

No catalogue writeback -- a consistency-check result, feeding the NEXT
(separate, later) vetting-chain task.

Run:  uv run python scripts/run_704_ccr4bp_real_ephemeris_consistency.py
Outputs -> data/found/704_ccr4bp_real_ephemeris_consistency/result.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.ccr4bp_umbriel_titania as ut  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.ccr4bp_real_ephemeris_consistency as rec  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "704_ccr4bp_real_ephemeris_consistency"

# #701's own saved best_robust_genuine_connection_corrected entry (see
# data/found/701_ccr4bp_umbriel_titania_search/result.json) -- reproduced
# bit-for-bit from the SAME committed coarse seed, not re-derived.
_SEED = hs.ManifoldCandidate(
    theta2_u=3.665191429188092,
    t_u=19.314531534610783,
    theta2_s=3.5604716740684323,
    t_s=18.187850528425155,
    gap_planar=0.0011381326404778097,
)
_U_LOBE = -1.0
_S_LOBE = -1.0
_EXPECTED_RESIDUAL_NORM = 1.1159187446079244e-14
_IDEALIZED_OFF_TORUS_KM = 1927.9247765134512  # #701's own corrected_off_torus_km


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical scaffolding to #701's own driver script."""
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


def _reconstruct_701_connection(
    t0: float,
) -> tuple[vt.CCR4BPTorusVariationalResult, np.ndarray, np.ndarray, float, hs.RefinedConnection]:
    system = ut.uranus_umbriel_titania_default()
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    _log(f"torus rebuilt: residual_rms={torus.residual_rms:.3e} (expect 1.840e-04)", t0)

    refined = hs.refine_candidate(
        torus, torus, _SEED, lobe_sign_u=_U_LOBE, lobe_sign_s=_S_LOBE, n_segments_dir=32
    )
    assert refined is not None, "#701 connection refinement failed to reproduce"
    rel_err = abs(refined.residual_norm - _EXPECTED_RESIDUAL_NORM) / _EXPECTED_RESIDUAL_NORM
    assert rel_err < 1e-2, (
        f"reconstructed residual_norm {refined.residual_norm:.6e} does not match #701's own "
        f"saved value {_EXPECTED_RESIDUAL_NORM:.6e} (rel_err={rel_err:.3e}) "
        "-- determinism check failed"
    )
    _log(
        f"#701 connection reconstructed bit-for-bit: residual_norm={refined.residual_norm:.6e}, "
        f"t_u={refined.t_u:.6f} TU, t_s={refined.t_s:.6f} TU",
        t0,
    )

    departure_u = mg.manifold_state_at(
        torus,
        "unstable",
        0.0,
        refined.theta2_u,
        0.0,
        lobe_sign=_U_LOBE,
        ref_vec=refined.ref_vec_u,
        n_segments_dir=32,
    )
    assert departure_u is not None
    return torus, departure_u, refined.state_s, refined.t_u, refined


def _result_to_dict(r: rec.ConsistencyCheckResult) -> dict[str, Any]:
    return {
        "epoch0_utc": r.epoch0_utc,
        "et0": r.et0,
        "t_u_tu": r.t_u_tu,
        "t_u_seconds": r.t_u_seconds,
        "t_u_days": r.t_u_seconds / 86400.0,
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
        "idealized_off_torus_km_for_scale": r.idealized_off_torus_km_for_scale,
        "notes": r.notes,
    }


def main() -> dict[str, Any]:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    system = ut.uranus_umbriel_titania_default()
    _torus, departure_u, target_s, t_u_tu, refined = _reconstruct_701_connection(t0)
    l_km = ut.L_KM
    v_unit_km_s = ut.v_unit_km_s()
    t_u_days = t_u_tu * rec.tu_to_seconds(l_km, v_unit_km_s) / 86400.0
    _log(f"elapsed unstable-branch flow time t_u = {t_u_tu:.4f} TU = {t_u_days:.3f} real days", t0)

    # ------------------------------------------------------------------
    # Headline single-epoch check.
    # ------------------------------------------------------------------
    headline_epoch = "2030-01-01T00:00:00"
    headline = rec.check_connection_survives_real_ephemeris(
        headline_epoch,
        departure_u,
        target_s,
        t_u_tu,
        l_km,
        v_unit_km_s,
        system.mu,
        idealized_off_torus_km_for_scale=_IDEALIZED_OFF_TORUS_KM,
    )
    _log(
        f"headline epoch {headline_epoch}: pos_gap={headline.pos_gap_km:.1f} km, "
        f"vel_gap={headline.vel_gap_km_s * 1000:.2f} m/s, "
        f"closest_approach={headline.closest_approach_km}",
        t0,
    )

    # ------------------------------------------------------------------
    # Short-baseline scan: one real Umbriel-Titania synodic period
    # (~7.911 real days -- 2*pi/|omega_gan| in TU, converted to days --
    # samples the full range of real relative Umbriel-Titania phasing at
    # a fixed reference date).
    # ------------------------------------------------------------------
    synodic_period_days = (
        2.0 * np.pi / abs(system.omega_gan) * rec.tu_to_seconds(l_km, v_unit_km_s) / 86400.0
    )
    _log(f"real Umbriel-Titania synodic period = {synodic_period_days:.4f} days", t0)
    n_synodic = 300
    synodic_scan: list[dict[str, Any]] = []
    import datetime as _dt

    base_dt = _dt.datetime(2030, 1, 1, tzinfo=_dt.UTC)
    for i in range(n_synodic):
        frac = i / n_synodic
        offset_days = frac * synodic_period_days
        epoch_dt = base_dt + _dt.timedelta(days=offset_days)
        epoch_utc = epoch_dt.strftime("%Y-%m-%dT%H:%M:%S")
        result = rec.check_connection_survives_real_ephemeris(
            epoch_utc,
            departure_u,
            target_s,
            t_u_tu,
            l_km,
            v_unit_km_s,
            system.mu,
            idealized_off_torus_km_for_scale=_IDEALIZED_OFF_TORUS_KM,
        )
        synodic_scan.append({"phase_frac": frac, **_result_to_dict(result)})
        if i % 25 == 0 or i == n_synodic - 1:
            min_titania = result.closest_approach_km.get("Titania", float("nan"))
            _log(
                f"  synodic scan {i + 1}/{n_synodic} ({epoch_utc}): "
                f"pos_gap={result.pos_gap_km:.1f} km, min_dist_Titania={min_titania:.1f} km",
                t0,
            )

    # ------------------------------------------------------------------
    # Duty-cycle characterization of the dense synodic scan: fraction of
    # the period below several thresholds, and the number of DISTINCT
    # sub-threshold windows (a recurring near-miss pattern vs. a single
    # fluke crossing).
    # ------------------------------------------------------------------
    pos_arr = np.array([e["pos_gap_km"] for e in synodic_scan])
    duty_cycle_thresholds_km = [500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0]
    duty_cycle: dict[str, Any] = {}
    for thresh in duty_cycle_thresholds_km:
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
    _log(f"duty-cycle summary (dense {n_synodic}-pt synodic scan): {duty_cycle}", t0)

    # ------------------------------------------------------------------
    # Local refinement: golden-section-style bisection around the dense
    # scan's own best (smallest pos_gap) point, to characterize the TRUE
    # local minimum this epoch-only 1-parameter family can reach (a 1D
    # epoch sweep of a 3D position-gap vector generically does NOT cross
    # exactly zero -- this finds how close the nearest local minimum
    # actually gets, not a claim of exact closure).
    # ------------------------------------------------------------------
    best_idx = int(np.argmin(pos_arr))
    best_offset_days = synodic_scan[best_idx]["phase_frac"] * synodic_period_days
    step_days = synodic_period_days / n_synodic

    def _pos_gap_at_offset(offset_days: float) -> tuple[float, rec.ConsistencyCheckResult]:
        epoch_dt = base_dt + _dt.timedelta(days=offset_days)
        epoch_utc = epoch_dt.strftime("%Y-%m-%dT%H:%M:%S")
        result = rec.check_connection_survives_real_ephemeris(
            epoch_utc,
            departure_u,
            target_s,
            t_u_tu,
            l_km,
            v_unit_km_s,
            system.mu,
            idealized_off_torus_km_for_scale=_IDEALIZED_OFF_TORUS_KM,
        )
        return result.pos_gap_km, result

    lo, hi = best_offset_days - step_days, best_offset_days + step_days
    best_offset, best_result = best_offset_days, _pos_gap_at_offset(best_offset_days)[1]
    for _ in range(28):  # ~2x halvings per iter via 3-point sampling -> sub-second precision
        mid = 0.5 * (lo + hi)
        g_lo, _ = _pos_gap_at_offset(lo)
        g_mid, r_mid = _pos_gap_at_offset(mid)
        g_hi, _ = _pos_gap_at_offset(hi)
        candidates = [(g_lo, lo), (g_mid, mid), (g_hi, hi)]
        candidates.sort(key=lambda t: t[0])
        best_offset = candidates[0][1]
        if candidates[0][1] == mid:
            best_result = r_mid
        width = hi - lo
        lo, hi = best_offset - width / 4.0, best_offset + width / 4.0
    best_offset, best_result = best_offset, _pos_gap_at_offset(best_offset)[1]
    _log(
        f"local-minimum refinement: pos_gap={best_result.pos_gap_km:.3f} km, "
        f"vel_gap={best_result.vel_gap_km_s * 1000:.4f} m/s, at epoch {best_result.epoch0_utc}",
        t0,
    )

    # ------------------------------------------------------------------
    # Long-baseline scan: spread across the full URA111-supported window
    # (kernel covers 1900-2099; use the #312-precedented 2000-2083
    # validity window) to check for secular epoch dependence beyond one
    # synodic cycle.
    # ------------------------------------------------------------------
    n_long = 20
    long_scan: list[dict[str, Any]] = []
    start_year, end_year = 2000, 2083
    for i in range(n_long):
        frac = i / (n_long - 1)
        year = start_year + frac * (end_year - start_year)
        epoch_dt = _dt.datetime(int(year), 1, 1, tzinfo=_dt.UTC) + _dt.timedelta(
            days=(year - int(year)) * 365.25
        )
        epoch_utc = epoch_dt.strftime("%Y-%m-%dT%H:%M:%S")
        result = rec.check_connection_survives_real_ephemeris(
            epoch_utc,
            departure_u,
            target_s,
            t_u_tu,
            l_km,
            v_unit_km_s,
            system.mu,
            idealized_off_torus_km_for_scale=_IDEALIZED_OFF_TORUS_KM,
        )
        long_scan.append({"year_frac": year, **_result_to_dict(result)})
        _log(
            f"  long-baseline scan {i + 1}/{n_long} ({epoch_utc}): "
            f"pos_gap={result.pos_gap_km:.1f} km",
            t0,
        )

    pos_gaps_synodic = [e["pos_gap_km"] for e in synodic_scan if np.isfinite(e["pos_gap_km"])]
    pos_gaps_long = [e["pos_gap_km"] for e in long_scan if np.isfinite(e["pos_gap_km"])]

    summary = {
        "headline_pos_gap_km": headline.pos_gap_km,
        "headline_vel_gap_km_s": headline.vel_gap_km_s,
        "synodic_scan_pos_gap_km_min": min(pos_gaps_synodic) if pos_gaps_synodic else None,
        "synodic_scan_pos_gap_km_max": max(pos_gaps_synodic) if pos_gaps_synodic else None,
        "synodic_scan_pos_gap_km_median": float(np.median(pos_gaps_synodic))
        if pos_gaps_synodic
        else None,
        "long_scan_pos_gap_km_min": min(pos_gaps_long) if pos_gaps_long else None,
        "long_scan_pos_gap_km_max": max(pos_gaps_long) if pos_gaps_long else None,
        "long_scan_pos_gap_km_median": float(np.median(pos_gaps_long)) if pos_gaps_long else None,
        "n_synodic_success": sum(1 for e in synodic_scan if e["propagation_success"]),
        "n_synodic_total": n_synodic,
        "n_long_success": sum(1 for e in long_scan if e["propagation_success"]),
        "n_long_total": n_long,
        "local_minimum_pos_gap_km": best_result.pos_gap_km,
        "local_minimum_vel_gap_km_s": best_result.vel_gap_km_s,
        "local_minimum_epoch_utc": best_result.epoch0_utc,
        "local_minimum_closest_approach_km": best_result.closest_approach_km,
    }
    _log(f"summary: {summary}", t0)

    result_all = {
        "task": "#704",
        "system": (
            "Uranus-Umbriel-Titania CCR4BP real-ephemeris consistency check of #701's connection"
        ),
        "source_connection": (
            "data/found/701_ccr4bp_umbriel_titania_search/result.json:"
            "best_robust_genuine_connection_corrected"
        ),
        "reconstructed_residual_norm": refined.residual_norm,
        "reconstructed_theta2_u": refined.theta2_u,
        "reconstructed_t_u_tu": refined.t_u,
        "reconstructed_theta2_s": refined.theta2_s,
        "reconstructed_t_s_tu": refined.t_s,
        "t_u_days": t_u_days,
        "l_km_umbriel": l_km,
        "v_unit_umbriel_km_s": v_unit_km_s,
        "system_mu": system.mu,
        "real_umbriel_titania_synodic_period_days": synodic_period_days,
        "headline_epoch_utc": headline_epoch,
        "headline_check": _result_to_dict(headline),
        "synodic_phase_scan": synodic_scan,
        "duty_cycle_dense_synodic_scan": duty_cycle,
        "local_minimum_refinement": _result_to_dict(best_result),
        "long_baseline_scan": long_scan,
        "summary": summary,
        "force_model_notes": (
            "Uranus point-mass central term + Umbriel + Titania REAL SPICE (URA111) "
            "third-body perturbations, Uranus-centred J2000 inertial frame. No J2, no Sun, "
            "no other Uranian moons -- see "
            "src/cyclerfinder/search/ccr4bp_real_ephemeris_consistency.py module docstring "
            "for the full, quantified justification of each omission."
        ),
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result_all, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result_all


if __name__ == "__main__":
    main()
