"""#392 V3 REBOUND IAS15 n-body cross-check for the low-amplitude floquet cycler.

V3 question: does the V2 round-off-floor bounded-drift signature survive an
independent integrator + a more realistic dynamical model?

Usage:
    uv run python scripts/run_392_v3_verify.py
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import rebound

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.search.reachable_representatives import braik_ross_system

V1_VERDICT_PATH = Path("data/branch_392_v1_verdict.jsonl")
V2_VERDICT_PATH = Path("data/branch_392_v2_verdict.jsonl")
PHASE_LABEL = "392_gauntlet_v3"
V3_AGREEMENT_FLOOR_KMS = 100.0  # spec §14 V3 V3-vs-V2 agreement floor
EXPECTED_INTEGRATOR = "REBOUND IAS15"
IAS15_EPSILON = 1e-12


def _load_v1() -> tuple[str, np.ndarray, float]:
    """Use the V1-corrected state + period (the V1-converged closure)."""
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


def _inertial_to_rotating(
    r_in_nondim: np.ndarray, v_in_nondim: np.ndarray, t_nondim: float
) -> np.ndarray:
    psi = t_nondim
    cos_psi, sin_psi = np.cos(psi), np.sin(psi)
    rmat_t = np.array(
        [[cos_psi, sin_psi, 0.0], [-sin_psi, cos_psi, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    r_rot = rmat_t @ r_in_nondim
    v_rot = rmat_t @ v_in_nondim - np.array([-r_rot[1], r_rot[0], 0.0], dtype=np.float64)
    return np.concatenate([r_rot, v_rot])


def _run_ias15_n_cycles(
    state0_rot_nondim: np.ndarray,
    period_nondim: float,
    n_cycles: int,
    system: cr3bp.CR3BPSystem,
) -> tuple[list[float], list[bool], str, dict[str, float]]:
    l_km = system.l_km
    t_s = system.t_s
    v_scale = l_km / t_s

    omega = 1.0 / t_s
    total_gm = omega * omega * l_km**3
    moon_gm = total_gm * system.mu
    earth_gm = total_gm * (1.0 - system.mu)
    g_km3_kg_s2 = 6.67430e-20
    m_moon_kg = moon_gm / g_km3_kg_s2
    m_earth_kg = earth_gm / g_km3_kg_s2

    real_earth_gm = PRIMARIES["Earth"]
    gm_relative_diff = (total_gm - real_earth_gm) / real_earth_gm

    mu_nondim = system.mu
    earth_r_km_init = np.array([-mu_nondim * l_km, 0.0, 0.0], dtype=np.float64)
    moon_r_km_init = np.array([(1.0 - mu_nondim) * l_km, 0.0, 0.0], dtype=np.float64)
    earth_v_kms = np.array([0.0, -mu_nondim * omega * l_km, 0.0], dtype=np.float64)
    moon_v_kms = np.array([0.0, (1.0 - mu_nondim) * omega * l_km, 0.0], dtype=np.float64)

    r_sc_nondim, v_sc_nondim = _rotating_to_inertial(state0_rot_nondim, 0.0, mu_nondim)
    sc_r_km = r_sc_nondim * l_km
    sc_v_kms = v_sc_nondim * v_scale

    sim = rebound.Simulation()
    sim.G = g_km3_kg_s2
    sim.integrator = "ias15"
    sim.integrator.epsilon = IAS15_EPSILON
    sim.add(
        m=m_earth_kg,
        x=float(earth_r_km_init[0]),
        y=float(earth_r_km_init[1]),
        z=float(earth_r_km_init[2]),
        vx=float(earth_v_kms[0]),
        vy=float(earth_v_kms[1]),
        vz=float(earth_v_kms[2]),
    )
    sim.add(
        m=m_moon_kg,
        x=float(moon_r_km_init[0]),
        y=float(moon_r_km_init[1]),
        z=float(moon_r_km_init[2]),
        vx=float(moon_v_kms[0]),
        vy=float(moon_v_kms[1]),
        vz=float(moon_v_kms[2]),
    )
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

    per_cycle_drift_km: list[float] = []
    per_cycle_converged: list[bool] = []
    state0_rot_pos = state0_rot_nondim[:3] * l_km

    for k in range(1, n_cycles + 1):
        t_target_sec = float(k * period_nondim * t_s)
        try:
            sim.integrate(t_target_sec)
        except Exception:
            per_cycle_drift_km.append(float("inf"))
            per_cycle_converged.append(False)
            continue
        sc_p = sim.particles[2]
        sc_r_km_now = np.array([sc_p.x, sc_p.y, sc_p.z], dtype=np.float64)
        sc_v_kms_now = np.array([sc_p.vx, sc_p.vy, sc_p.vz], dtype=np.float64)
        if not np.all(np.isfinite(sc_r_km_now)) or not np.all(np.isfinite(sc_v_kms_now)):
            per_cycle_drift_km.append(float("inf"))
            per_cycle_converged.append(False)
            continue

        sc_r_nondim_now = sc_r_km_now / l_km
        sc_v_nondim_now = sc_v_kms_now / v_scale
        t_nondim = k * period_nondim
        rot_state = _inertial_to_rotating(sc_r_nondim_now, sc_v_nondim_now, t_nondim)
        rot_r_km = rot_state[:3] * l_km
        drift = float(np.linalg.norm(rot_r_km - state0_rot_pos))
        per_cycle_drift_km.append(drift)
        per_cycle_converged.append(True)

    meta = {
        "cr3bp_implied_total_gm_km3s2": float(total_gm),
        "registry_earth_total_gm_km3s2": float(real_earth_gm),
        "gm_relative_diff": float(gm_relative_diff),
        "omega_rad_per_sec": float(omega),
        "m_earth_kg": float(m_earth_kg),
        "m_moon_kg": float(m_moon_kg),
    }
    return per_cycle_drift_km, per_cycle_converged, EXPECTED_INTEGRATOR, meta


def _load_v2_drifts() -> dict[int, list[float]]:
    drifts: dict[int, list[float]] = {}
    with V2_VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "v2_verdict_cr3bp_periodic":
                drifts[int(row["n_cycles_requested"])] = [
                    float(x) for x in row["per_cycle_drift_kms"]
                ]
    return drifts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_392_v3_verdict.jsonl"),
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    candidate_id, state0, period = _load_v1()
    system = braik_ross_system()
    v2_drifts = _load_v2_drifts()

    print(f"[V3-verify] candidate_id={candidate_id}")
    print(f"[V3-verify] period_TU={period:.15f}")

    verdicts: dict[int, dict[str, Any]] = {}
    meta_record: dict[str, float] = {}
    all_pass = True

    # We only run V3 on cycles that passed V2. For this branch, it's 3 and 5.
    passed_v2_cycles = [n for n, drifts in v2_drifts.items() if max(drifts) < 50000.0]
    passed_v2_cycles.sort()

    for n in passed_v2_cycles:
        print(f"[V3-verify] running IAS15 at n_cycles={n}...")
        per_cycle_v3_km, per_cycle_conv, integ_label, meta = _run_ias15_n_cycles(
            state0, period, n, system
        )
        meta_record = meta
        v2_drift_n = v2_drifts.get(n, [])
        agreement = [
            float(abs(v3 - v2)) for v3, v2 in zip(per_cycle_v3_km, v2_drift_n, strict=False)
        ]
        max_agreement = max(agreement) if agreement else float("inf")
        max_v3_drift = max(per_cycle_v3_km) if per_cycle_v3_km else float("inf")
        all_legs_converged = all(per_cycle_conv)
        passes_v3 = bool(all_legs_converged and max_agreement < V3_AGREEMENT_FLOOR_KMS)
        all_pass = all_pass and passes_v3
        print(
            f"  max_v3_drift={max_v3_drift:.4e} km, "
            f"max_v3_vs_v2_agreement={max_agreement:.4e} km, "
            f"all_converged={all_legs_converged}, passes_v3={passes_v3}"
        )
        verdicts[n] = {
            "kind": "v3_verdict_cr3bp_periodic",
            "candidate_id": candidate_id,
            "n_cycles_requested": int(n),
            "n_cycles_propagated": len([c for c in per_cycle_conv if c]),
            "per_cycle_v3_drift_kms": [float(x) for x in per_cycle_v3_km],
            "per_cycle_v2_drift_kms": [float(x) for x in v2_drift_n],
            "per_cycle_v3_minus_v2_kms": [float(x) for x in agreement],
            "max_v3_drift_kms": float(max_v3_drift),
            "drift_agreement_kms": float(max_agreement),
            "agreement_floor_kms": V3_AGREEMENT_FLOOR_KMS,
            "integrator": integ_label,
            "ias15_epsilon": IAS15_EPSILON,
            "converged_at_each_cycle": all_legs_converged,
            "passes_v3": passes_v3,
        }

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
                    "v3_agreement_floor_kms": V3_AGREEMENT_FLOOR_KMS,
                    "n_cycles_list": list(passed_v2_cycles),
                    "integrator_target": EXPECTED_INTEGRATOR,
                    "cr3bp_setup_meta": meta_record,
                }
            )
            + "\n"
        )
        for n in passed_v2_cycles:
            fh.write(json.dumps(verdicts[n]) + "\n")
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
    print(f"[V3-verify] verdict written to {args.output}")
    print(f"[V3-verify] all V3 gates pass: {all_pass}")


if __name__ == "__main__":
    main()
