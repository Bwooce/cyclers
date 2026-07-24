"""#705 -- epoch-robustness scan for `#704`'s near-miss window.

`#704` found that `#701`'s idealized (circular-coplanar) CCR4BP Umbriel-
Titania homoclinic connection does NOT survive real Umbriel/Titania ephemeris
generically (median position mismatch ~125,000-136,000 km across a synodic
period -- comparable to Titania's own orbital radius), but WITHIN a dense
300-point scan across one real ~7.909-day Umbriel-Titania synodic period
(centred near 2030-01-01), it found a real, narrow, non-trivial near-miss
structure: 3 distinct sub-5000km-mismatch windows per synodic period (~1%
duty cycle), tightening at its own best point (2030-01-07T22:25:06 UTC) to
``pos_gap_km ~= 84.460`` km / ``vel_gap_km_s ~= 0.005905`` km/s.

`#704`'s OWN long-baseline scan (20 points spread across 2000-2083, spacing
~4.4 years) cannot confirm or rule out whether this same narrow structure
recurs at OTHER multi-year epochs, because a sample spacing of ~4.4 years is
vastly coarser than the ~7.909-day synodic period the narrow window lives
inside -- by construction it can miss the phenomenon at every one of its 20
sample points even if the phenomenon is present there too.

This script closes that gap: it reruns `#704`'s OWN dense-synodic-period-scan
methodology (unmodified -- imports
:mod:`cyclerfinder.search.ccr4bp_real_ephemeris_consistency` as a library,
exactly as that module's own docstring intends), but applies it at SEVERAL
well-separated multi-year epochs spread across the full 2000-2083
`#312`-precedented real-ephemeris validity window (see
``umbriel-titania-1-1-uranian-quasi-cycler-2026``'s own catalogue row), not
just the single epoch `#704` happened to check.

Epoch choice
------------
10 epochs, evenly spaced across [2000, 2083] (``frac = i/9``, ``i=0..9``),
giving a ~9.22-year spacing -- comparable in spirit to `#704`'s own 20-point
long-baseline scan's ~4.4-year spacing order of magnitude, but MUCH sparser
in COUNT because here each epoch itself gets `#704`'s own full 300-point
dense synodic-period sub-scan (a 10x-epoch-count blow-up of the expensive
part would have cost ~10x #704's own total runtime for little extra
information -- 10 well-separated epochs is enough to distinguish "recurs at
every tested epoch" from "occurred at one epoch only" from "occurs at some
but not all," which is the actual question this task asks).

At EACH epoch, this script:

1. Runs a dense 300-point scan across one real synodic period starting at
   that epoch (identical resolution and methodology to `#704`'s own
   ``synodic_phase_scan`` -- same thresholds, same duty-cycle window-counting
   logic).
2. Locates that scan's own best (smallest ``pos_gap_km``) point and refines
   it with the SAME golden-section-style bisection `#704`'s own script used
   (28 iterations, 3-point sampling) to find the TRUE local minimum within
   that epoch's own synodic period, not just the coarse-grid best sample.

"Comparable to #704's own 2030 near-miss" is judged against an EXPLICIT,
stated threshold: within one order of magnitude (10x) of `#704`'s own
committed best point (``pos_gap_km=84.460...`` -> threshold ``845 km``;
``vel_gap_km_s=0.005905...`` -> threshold ``0.059 km/s``), sourced directly
from `#704`'s own committed
``data/found/704_ccr4bp_real_ephemeris_consistency/result.json``.

No catalogue writeback (out of scope, per `#705`'s own dispatch). Does NOT
modify `#704`'s own module or any of `#689`-`#694`'s/`#701`'s modules --
consumes :mod:`cyclerfinder.search.ccr4bp_real_ephemeris_consistency`'s
public API exactly as `#704`'s own driver script did, and duplicates ONLY
the reconstruction scaffolding `#704`'s own script duplicated from `#701`'s
(an already-precedented pattern in this exact arc -- see `#704`'s own
``_resonant_symmetric_orbit`` docstring: "Identical scaffolding to #701's
own driver script.").

Run:  uv run python scripts/run_705_epoch_robustness_scan.py
Outputs -> data/found/705_ccr4bp_epoch_robustness_scan/result.json
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

import cyclerfinder.core.ccr4bp_umbriel_titania as ut  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.ccr4bp_real_ephemeris_consistency as rec  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "705_ccr4bp_epoch_robustness_scan"

# #701's own saved best_robust_genuine_connection_corrected entry (see
# data/found/701_ccr4bp_umbriel_titania_search/result.json), reproduced
# bit-for-bit -- identical scaffolding/constants to #704's own driver script
# (scripts/run_704_ccr4bp_real_ephemeris_consistency.py), NOT re-derived.
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

# #704's own committed best-point reference (data/found/
# 704_ccr4bp_real_ephemeris_consistency/result.json:summary), used ONLY as
# the "comparable tightness" scale below -- not recomputed here.
_REF_704_BEST_POS_GAP_KM = 84.46019822482435
_REF_704_BEST_VEL_GAP_KM_S = 0.005904706092575937
_ORDER_OF_MAGNITUDE_FACTOR = 10.0
_COMPARABLE_POS_GAP_THRESHOLD_KM = _ORDER_OF_MAGNITUDE_FACTOR * _REF_704_BEST_POS_GAP_KM
_COMPARABLE_VEL_GAP_THRESHOLD_KM_S = _ORDER_OF_MAGNITUDE_FACTOR * _REF_704_BEST_VEL_GAP_KM_S

N_EPOCHS = 10
N_SYNODIC = 300
DUTY_CYCLE_THRESHOLDS_KM = [500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0]
START_YEAR, END_YEAR = 2000, 2083


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical scaffolding to #701's/#704's own driver scripts."""
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


def _slim_result(phase_frac: float, r: rec.ConsistencyCheckResult) -> dict[str, Any]:
    """Slim per-point record (no full state vectors) -- keeps the
    10-epoch x 300-point scan's own result.json a reasonable size while
    retaining everything needed to recompute duty cycle / min / median."""
    return {
        "phase_frac": phase_frac,
        "epoch_utc": r.epoch0_utc,
        "pos_gap_km": r.pos_gap_km,
        "vel_gap_km_s": r.vel_gap_km_s,
        "propagation_success": r.propagation_success,
        "closest_approach_km": r.closest_approach_km,
    }


def _full_result(r: rec.ConsistencyCheckResult) -> dict[str, Any]:
    """Full per-point record (with state vectors) -- used only for the
    single best/local-minimum point per epoch."""
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


def _duty_cycle(pos_arr: np.ndarray) -> dict[str, Any]:
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


def _scan_one_epoch(
    base_dt: _dt.datetime,
    synodic_period_days: float,
    departure_u: np.ndarray,
    target_s: np.ndarray,
    t_u_tu: float,
    l_km: float,
    v_unit_km_s: float,
    mu: float,
    t0: float,
    epoch_idx: int,
) -> dict[str, Any]:
    synodic_scan: list[dict[str, Any]] = []
    all_results: list[rec.ConsistencyCheckResult] = []
    for i in range(N_SYNODIC):
        frac = i / N_SYNODIC
        offset_days = frac * synodic_period_days
        epoch_dt = base_dt + _dt.timedelta(days=offset_days)
        epoch_utc = epoch_dt.strftime("%Y-%m-%dT%H:%M:%S")
        result = rec.check_connection_survives_real_ephemeris(
            epoch_utc, departure_u, target_s, t_u_tu, l_km, v_unit_km_s, mu
        )
        synodic_scan.append(_slim_result(frac, result))
        all_results.append(result)
    pos_arr = np.array([e["pos_gap_km"] for e in synodic_scan])
    duty_cycle = _duty_cycle(pos_arr)

    best_idx = int(np.argmin(pos_arr))
    best_offset_days = synodic_scan[best_idx]["phase_frac"] * synodic_period_days
    step_days = synodic_period_days / N_SYNODIC

    def _pos_gap_at_offset(offset_days: float) -> tuple[float, rec.ConsistencyCheckResult]:
        epoch_dt = base_dt + _dt.timedelta(days=offset_days)
        epoch_utc = epoch_dt.strftime("%Y-%m-%dT%H:%M:%S")
        result = rec.check_connection_survives_real_ephemeris(
            epoch_utc, departure_u, target_s, t_u_tu, l_km, v_unit_km_s, mu
        )
        return result.pos_gap_km, result

    lo, hi = best_offset_days - step_days, best_offset_days + step_days
    best_offset, best_result = best_offset_days, _pos_gap_at_offset(best_offset_days)[1]
    for _ in range(28):
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

    comparable = (
        best_result.pos_gap_km < _COMPARABLE_POS_GAP_THRESHOLD_KM
        and best_result.vel_gap_km_s < _COMPARABLE_VEL_GAP_THRESHOLD_KM_S
    )
    _log(
        f"epoch {epoch_idx + 1}/{N_EPOCHS} base={base_dt.strftime('%Y-%m-%d')}: "
        f"synodic min={pos_arr.min():.1f} km, median={float(np.median(pos_arr)):.1f} km, "
        f"local-min pos_gap={best_result.pos_gap_km:.3f} km "
        f"vel_gap={best_result.vel_gap_km_s * 1000:.4f} m/s at {best_result.epoch0_utc} "
        f"[{'COMPARABLE' if comparable else 'not comparable'} to #704's 2030 window]",
        t0,
    )

    return {
        "base_epoch_utc": base_dt.strftime("%Y-%m-%dT%H:%M:%S"),
        "synodic_scan_pos_gap_km_min": float(pos_arr.min()),
        "synodic_scan_pos_gap_km_max": float(pos_arr.max()),
        "synodic_scan_pos_gap_km_median": float(np.median(pos_arr)),
        "duty_cycle_dense_synodic_scan": duty_cycle,
        "local_minimum": _full_result(best_result),
        "comparable_to_704_2030_window": comparable,
        "synodic_phase_scan": synodic_scan,
    }


def main() -> dict[str, Any]:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    system = ut.uranus_umbriel_titania_default()
    _torus, departure_u, target_s, t_u_tu, refined = _reconstruct_701_connection(t0)
    l_km = ut.L_KM
    v_unit_km_s = ut.v_unit_km_s()

    synodic_period_days = (
        2.0 * np.pi / abs(system.omega_gan) * rec.tu_to_seconds(l_km, v_unit_km_s) / 86400.0
    )
    _log(f"real Umbriel-Titania synodic period = {synodic_period_days:.4f} days", t0)
    _log(
        f"comparable-tightness thresholds (10x #704's own best point): "
        f"pos_gap < {_COMPARABLE_POS_GAP_THRESHOLD_KM:.1f} km, "
        f"vel_gap < {_COMPARABLE_VEL_GAP_THRESHOLD_KM_S * 1000:.2f} m/s",
        t0,
    )

    epoch_scans: list[dict[str, Any]] = []
    for i in range(N_EPOCHS):
        frac = i / (N_EPOCHS - 1)
        year = START_YEAR + frac * (END_YEAR - START_YEAR)
        base_dt = _dt.datetime(int(year), 1, 1, tzinfo=_dt.UTC) + _dt.timedelta(
            days=(year - int(year)) * 365.25
        )
        scan = _scan_one_epoch(
            base_dt,
            synodic_period_days,
            departure_u,
            target_s,
            t_u_tu,
            l_km,
            v_unit_km_s,
            system.mu,
            t0,
            i,
        )
        scan["epoch_year"] = year
        epoch_scans.append(scan)

    n_comparable = sum(1 for e in epoch_scans if e["comparable_to_704_2030_window"])
    best_per_epoch = [e["local_minimum"]["pos_gap_km"] for e in epoch_scans]
    summary = {
        "n_epochs_tested": N_EPOCHS,
        "n_synodic_points_per_epoch": N_SYNODIC,
        "comparable_pos_gap_threshold_km": _COMPARABLE_POS_GAP_THRESHOLD_KM,
        "comparable_vel_gap_threshold_km_s": _COMPARABLE_VEL_GAP_THRESHOLD_KM_S,
        "reference_704_best_pos_gap_km": _REF_704_BEST_POS_GAP_KM,
        "reference_704_best_vel_gap_km_s": _REF_704_BEST_VEL_GAP_KM_S,
        "n_epochs_comparable_to_704": n_comparable,
        "fraction_epochs_comparable": n_comparable / N_EPOCHS,
        "per_epoch_local_minimum_pos_gap_km": best_per_epoch,
        "per_epoch_local_minimum_vel_gap_km_s": [
            e["local_minimum"]["vel_gap_km_s"] for e in epoch_scans
        ],
        "per_epoch_year": [e["epoch_year"] for e in epoch_scans],
        "overall_best_pos_gap_km": min(best_per_epoch),
        "overall_worst_pos_gap_km": max(best_per_epoch),
    }
    _log(f"summary: {summary}", t0)

    result_all = {
        "task": "#705",
        "system": (
            "Uranus-Umbriel-Titania CCR4BP epoch-robustness scan of #704's near-miss window"
        ),
        "source_connection": (
            "data/found/701_ccr4bp_umbriel_titania_search/result.json:"
            "best_robust_genuine_connection_corrected"
        ),
        "reconstructed_residual_norm": refined.residual_norm,
        "t_u_tu": t_u_tu,
        "l_km_umbriel": l_km,
        "v_unit_umbriel_km_s": v_unit_km_s,
        "system_mu": system.mu,
        "real_umbriel_titania_synodic_period_days": synodic_period_days,
        "epoch_scans": epoch_scans,
        "summary": summary,
        "force_model_notes": (
            "Identical to #704's own force model: Uranus point-mass central term + "
            "Umbriel + Titania REAL SPICE (URA111) third-body perturbations, "
            "Uranus-centred J2000 inertial frame. No J2, no Sun, no other Uranian "
            "moons -- see "
            "src/cyclerfinder/search/ccr4bp_real_ephemeris_consistency.py module "
            "docstring for the full, quantified justification of each omission."
        ),
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result_all, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result_all


if __name__ == "__main__":
    main()
