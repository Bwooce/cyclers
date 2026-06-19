"""#393 V4 real-ephemeris (DE440) cross-check for the low-amplitude floquet cycler.

V4 question (spec §14 V4): does the V3 bounded-drift signature survive an
*independent codebase + real ephemeris* n-body realization?

Usage:
    uv run python scripts/run_393_v4_verify.py
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import rebound
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.time import Time

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.reachable_representatives import braik_ross_system

V1_VERDICT_PATH = Path("data/branch_393_v1_verdict.jsonl")
V3_VERDICT_PATH = Path("data/branch_393_v3_verdict.jsonl")
PHASE_LABEL = "393_gauntlet_v4"
V4_AGREEMENT_FLOOR_KMS = 50000.0
EXPECTED_INTEGRATOR = "REBOUND IAS15"
IAS15_EPSILON = 1e-12
EPHEMERIS = "de440"
DEFAULT_EPOCH_UTC = "2000-01-15T00:00:00"
G_KM3_KG_S2 = 6.67430e-20

_GM_KM3_S2 = {
    "earth": float(PRIMARIES["Earth"]) - float(SATELLITES["Moon"].mu_km3_s2),
    "moon": float(SATELLITES["Moon"].mu_km3_s2),
    "sun": 1.32712440018e11,
    "mars": 4.282837e4,
    "jupiter": 1.26686534e8,
}
_PERTURBER_BODIES = ("sun", "mars", "jupiter")


def _load_v1() -> tuple[str, np.ndarray, float]:
    with V1_VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "v1_verdict_cr3bp_periodic":
                candidate_id = row["candidate_id"]
                state = np.array(row["state0_corrected_nondim"], dtype=np.float64)
                period = float(row["period_corrected_TU"])
                return candidate_id, state, period
    raise AssertionError(f"V1 verdict row not found in {V1_VERDICT_PATH}")


def _load_v3_drifts() -> dict[int, list[float]]:
    drifts: dict[int, list[float]] = {}
    with V3_VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "v3_verdict_cr3bp_periodic":
                drifts[int(row["n_cycles_requested"])] = [
                    float(x) for x in row["per_cycle_v3_drift_kms"]
                ]
    return drifts


def _structural_diagnostic(
    state0_rot_nondim: np.ndarray, system: cr3bp.CR3BPSystem
) -> dict[str, float]:
    l_km = system.l_km
    r_amp_km = float(np.linalg.norm(state0_rot_nondim[:3]) * l_km)
    gm_earth_total = float(PRIMARIES["Earth"])
    gm_sun = 1.32712440018e11
    a_earth_sun_km = 1.495978707e8
    r_hill_km = a_earth_sun_km * (gm_earth_total / (3.0 * gm_sun)) ** (1.0 / 3.0)
    a_earth_g = gm_earth_total / r_amp_km**2
    a_sun_tide = 2.0 * gm_sun * r_amp_km / a_earth_sun_km**3
    return {
        "orbit_amplitude_km": r_amp_km,
        "earth_sun_hill_radius_km": r_hill_km,
        "amplitude_fraction_of_hill": r_amp_km / r_hill_km,
        "earth_gravity_accel_km_s2": a_earth_g,
        "solar_tide_accel_km_s2": a_sun_tide,
        "solar_tide_to_earth_gravity_ratio": a_sun_tide / a_earth_g,
    }


def _rotating_to_inertial(
    state_rot_nondim: np.ndarray, t_nondim: float, mu: float
) -> tuple[np.ndarray, np.ndarray]:
    psi = t_nondim
    cos_psi, sin_psi = np.cos(psi), np.sin(psi)
    rmat = np.array(
        [[cos_psi, -sin_psi, 0.0], [sin_psi, cos_psi, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    r_rot = state_rot_nondim[:3]
    v_rot = state_rot_nondim[3:6]
    omega_cross_r = np.array([-r_rot[1], r_rot[0], 0.0], dtype=np.float64)
    r_in = rmat @ r_rot
    v_in = rmat @ (v_rot + omega_cross_r)
    return r_in, v_in


def _ephem_state_rel_emb(body: str, t: Time) -> tuple[np.ndarray, np.ndarray]:
    emb_p, emb_v = get_body_barycentric_posvel("earth-moon-barycenter", t)
    p, v = get_body_barycentric_posvel(body, t)
    r = (p.xyz - emb_p.xyz).to(u.km).value
    vel = (v.xyz - emb_v.xyz).to(u.km / u.s).value
    return np.asarray(r, dtype=np.float64), np.asarray(vel, dtype=np.float64)


def _real_rotating_basis(t: Time) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    earth_r, earth_v = _ephem_state_rel_emb("earth", t)
    moon_r, moon_v = _ephem_state_rel_emb("moon", t)
    r_em = moon_r - earth_r
    v_em = moon_v - earth_v
    xhat = r_em / np.linalg.norm(r_em)
    h = np.cross(r_em, v_em)
    zhat = h / np.linalg.norm(h)
    yhat = np.cross(zhat, xhat)
    return xhat, yhat, zhat


def _run_ias15_n_cycles(
    state0_rot_nondim: np.ndarray,
    period_nondim: float,
    n_cycles: int,
    system: cr3bp.CR3BPSystem,
    epoch_utc: str,
    bodies: tuple[str, ...] = ("earth", "moon", "sun", "mars", "jupiter"),
) -> tuple[list[float], list[bool], dict[str, Any]]:
    l_km = system.l_km
    t_s = system.t_s
    v_scale = l_km / t_s
    mu_nondim = system.mu

    solar_system_ephemeris.set(EPHEMERIS)
    t0 = Time(epoch_utc, scale="tdb")

    xhat0, yhat0, zhat0 = _real_rotating_basis(t0)
    rot0 = np.column_stack([xhat0, yhat0, zhat0])

    r_sc_nd, v_sc_nd = _rotating_to_inertial(state0_rot_nondim, 0.0, mu_nondim)
    r_sc_cr3bp_km = r_sc_nd * l_km
    v_sc_cr3bp_kms = v_sc_nd * v_scale
    sc_r_km = rot0 @ r_sc_cr3bp_km
    sc_v_kms = rot0 @ v_sc_cr3bp_kms

    sim = rebound.Simulation()
    sim.G = G_KM3_KG_S2
    sim.integrator = "ias15"
    sim.integrator.epsilon = IAS15_EPSILON

    seeded: dict[str, dict[str, float]] = {}
    for body in bodies:
        r, v = _ephem_state_rel_emb(body, t0)
        m_kg = _GM_KM3_S2[body] / G_KM3_KG_S2
        sim.add(
            m=float(m_kg),
            x=float(r[0]),
            y=float(r[1]),
            z=float(r[2]),
            vx=float(v[0]),
            vy=float(v[1]),
            vz=float(v[2]),
        )
        seeded[body] = {
            "gm_km3_s2": _GM_KM3_S2[body],
            "x_km": float(r[0]),
            "y_km": float(r[1]),
            "z_km": float(r[2]),
            "vx_kms": float(v[0]),
            "vy_kms": float(v[1]),
            "vz_kms": float(v[2]),
        }
    sc_idx = len(bodies)
    sim.add(
        m=0.0,
        x=float(sc_r_km[0]),
        y=float(sc_r_km[1]),
        z=float(sc_r_km[2]),
        vx=float(sc_v_kms[0]),
        vy=float(sc_v_kms[1]),
        vz=float(sc_v_kms[2]),
    )
    sim.t = 0.0

    ic_r_rot_km = rot0.T @ sc_r_km

    per_cycle_drift_km: list[float] = []
    per_cycle_converged: list[bool] = []

    for k in range(1, n_cycles + 1):
        t_target_sec = float(k * period_nondim * t_s)
        try:
            sim.integrate(t_target_sec)
        except Exception:
            per_cycle_drift_km.append(float("inf"))
            per_cycle_converged.append(False)
            continue
        sc_p = sim.particles[sc_idx]
        sc_r_now = np.array([sc_p.x, sc_p.y, sc_p.z], dtype=np.float64)
        if not np.all(np.isfinite(sc_r_now)):
            per_cycle_drift_km.append(float("inf"))
            per_cycle_converged.append(False)
            continue
        t_now = t0 + (t_target_sec / 86400.0) * u.day
        xh, yh, zh = _real_rotating_basis(t_now)
        rot_now = np.column_stack([xh, yh, zh])
        sc_r_rot_km = rot_now.T @ sc_r_now
        drift = float(np.linalg.norm(sc_r_rot_km - ic_r_rot_km))
        per_cycle_drift_km.append(drift)
        per_cycle_converged.append(True)

    meta = {
        "ephemeris": EPHEMERIS,
        "epoch_utc": epoch_utc,
        "epoch_scale": "tdb",
        "bodies": list(bodies),
        "perturber_bodies": [b for b in bodies if b in _PERTURBER_BODIES],
        "seeded_bodies": seeded,
        "rotating_basis_xhat0": [float(v) for v in xhat0],
        "rotating_basis_zhat0": [float(v) for v in zhat0],
        "gm_km3_s2": dict(_GM_KM3_S2),
    }
    return per_cycle_drift_km, per_cycle_converged, meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_393_v4_verdict.jsonl"),
    )
    parser.add_argument("--epoch", type=str, default=DEFAULT_EPOCH_UTC)
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    candidate_id, state0, period = _load_v1()
    system = braik_ross_system()
    v3_drifts = _load_v3_drifts()

    print(f"[V4-verify] candidate_id={candidate_id}")
    print(f"[V4-verify] period_TU={period:.15f}")
    print(f"[V4-verify] ephemeris={EPHEMERIS}, epoch={args.epoch}")
    print(f"[V4-verify] integrator: {EXPECTED_INTEGRATOR}, epsilon={IAS15_EPSILON}")

    verdicts: dict[int, dict[str, Any]] = {}
    meta_record: dict[str, Any] = {}
    all_pass = True

    passed_v3_cycles = [n for n, drifts in v3_drifts.items() if max(drifts) < 50000.0]
    passed_v3_cycles.sort()

    for n in passed_v3_cycles:
        print(f"[V4-verify] running DE440 IAS15 at n_cycles={n}...")
        per_cycle_v4_km, per_cycle_conv, meta = _run_ias15_n_cycles(
            state0, period, n, system, args.epoch
        )
        meta_record = meta
        v3_drift_n = v3_drifts.get(n, [])
        agreement = [
            float(abs(v4 - v3)) for v4, v3 in zip(per_cycle_v4_km, v3_drift_n, strict=False)
        ]
        max_agreement = max(agreement) if agreement else float("inf")
        max_v4_drift = max(per_cycle_v4_km) if per_cycle_v4_km else float("inf")
        all_legs_converged = all(per_cycle_conv)
        passes_v4 = bool(all_legs_converged and max_agreement < V4_AGREEMENT_FLOOR_KMS)
        all_pass = all_pass and passes_v4
        print(
            f"  max_v4_drift={max_v4_drift:.4e} km, "
            f"max_v4_vs_v3_agreement={max_agreement:.4e} km, "
            f"all_converged={all_legs_converged}, passes_v4={passes_v4}"
        )
        verdicts[n] = {
            "kind": "v4_verdict_cr3bp_periodic",
            "candidate_id": candidate_id,
            "n_cycles_requested": int(n),
            "n_cycles_propagated": len([c for c in per_cycle_conv if c]),
            "per_cycle_v4_drift_kms": [float(x) for x in per_cycle_v4_km],
            "per_cycle_v3_drift_kms": [float(x) for x in v3_drift_n],
            "per_cycle_v4_minus_v3_kms": [float(x) for x in agreement],
            "max_v4_drift_kms": float(max_v4_drift),
            "drift_agreement_kms": float(max_agreement),
            "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
            "integrator": EXPECTED_INTEGRATOR,
            "ias15_epsilon": IAS15_EPSILON,
            "ephemeris": EPHEMERIS,
            "epoch_utc": args.epoch,
            "converged_at_each_cycle": all_legs_converged,
            "passes_v4": passes_v4,
        }

    n_ctrl = max(passed_v3_cycles) if passed_v3_cycles else 1
    print(f"[V4-verify] running Earth+Moon-only control at n_cycles={n_ctrl}...")
    ctrl_drift_km, ctrl_conv, _ = _run_ias15_n_cycles(
        state0, period, n_ctrl, system, args.epoch, bodies=("earth", "moon")
    )
    print(
        "[V4-verify] EM-only control max drift "
        f"{max(ctrl_drift_km) if ctrl_drift_km else float('inf'):.4e} km"
    )

    diag = _structural_diagnostic(state0, system)
    print(
        f"[V4-verify] structural diag: amplitude {diag['orbit_amplitude_km']:.3e} km "
        f"= {diag['amplitude_fraction_of_hill']:.3f} of Earth-Sun Hill; "
        f"solar tide / Earth gravity = {diag['solar_tide_to_earth_gravity_ratio']:.3f}"
    )

    elapsed = time.time() - t_start
    iso_end = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write(
            json.dumps(
                {
                    "kind": "header",
                    "candidate_id": candidate_id,
                    "phase": PHASE_LABEL,
                    "iso_start": iso_start,
                    "iso_end": iso_end,
                    "elapsed_seconds": elapsed,
                    "v4_agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                    "n_cycles_list": list(passed_v3_cycles),
                    "integrator_target": EXPECTED_INTEGRATOR,
                    "ias15_epsilon": IAS15_EPSILON,
                    "ephemeris": EPHEMERIS,
                    "epoch_utc": args.epoch,
                    "v4_setup_meta": meta_record,
                    "discipline": (
                        "Spec §14 V4 real-ephemeris cross-check: REBOUND IAS15 "
                        "n-body seeded from JPL DE440 (Earth + Moon + Sun + Mars "
                        "+ Jupiter as massive bodies; spacecraft a massless test "
                        "particle)."
                    ),
                }
            )
            + "\n"
        )
        for n in passed_v3_cycles:
            fh.write(json.dumps(verdicts[n]) + "\n")
        fh.write(
            json.dumps(
                {
                    "kind": "control_em_only",
                    "candidate_id": candidate_id,
                    "bodies": ["earth", "moon"],
                    "n_cycles_requested": int(n_ctrl),
                    "per_cycle_drift_kms": [float(x) for x in ctrl_drift_km],
                    "max_drift_kms": float(max(ctrl_drift_km)) if ctrl_drift_km else None,
                    "converged_at_each_cycle": bool(all(ctrl_conv)),
                    "note": (
                        "Earth+Moon-only (DE440 real eccentric/inclined lunar orbit, "
                        "NO Sun/Mars/Jupiter)."
                    ),
                }
            )
            + "\n"
        )
        fh.write(
            json.dumps(
                {
                    "kind": "structural_diagnostic",
                    "candidate_id": candidate_id,
                    **diag,
                    "note": ("Structural diagnostics to check susceptibility to solar tides."),
                }
            )
            + "\n"
        )
        fh.write(
            json.dumps(
                {
                    "kind": "footer",
                    "candidate_id": candidate_id,
                    "phase": PHASE_LABEL,
                    "iso_end": iso_end,
                    "all_pass": all_pass,
                }
            )
            + "\n"
        )
    print(f"[V4-verify] verdict written to {args.output}")
    print(f"[V4-verify] all V4 gates pass: {all_pass}")


if __name__ == "__main__":
    main()
