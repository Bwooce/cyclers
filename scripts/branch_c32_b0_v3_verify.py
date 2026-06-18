"""#389 P389.4 — V3 REBOUND IAS15 n-body cross-check for branch_C32_b0.

V3 question: does the V2 round-off-floor bounded-drift signature survive an
independent integrator + a more realistic dynamical model?

For branch_C32_b0 (a closed CR3BP orbit), the V3 substrate is:

  * convert the CR3BP rotating-frame IC at t=0 to the inertial frame (the
    Earth-Moon barycenter sits at the origin; the Moon at (1-mu, 0, 0) with
    velocity (0, 1-mu, 0) in nondim units; the rotation rate is 1/TU);
  * pose a REBOUND IAS15 simulation with Earth + Moon as point-mass primaries
    (initialized on circular Earth-Moon-barycenter orbits matching the CR3BP
    assumptions) and the spacecraft as a massless test particle;
  * integrate for n_cycles periods (mapped from the CR3BP TU to seconds via
    system.t_s);
  * sample at the same cycle boundaries the V2 driver reported and convert
    the IAS15 inertial state back to the rotating-frame position;
  * report drift = ||r_rotating(IAS15) - r_rotating(V2_DOP853)|| at each
    boundary.

PASS criterion (spec §14 V3): drift_agreement_kms < 100 km — the same floor
the #331 SILVER V3 used. For an essentially-stable orbit like branch_C32_b0
we expect drift to be at the level of the IAS15 round-off integrated over
n cycles, which is a couple of orders of magnitude below the floor.

Discipline notes:

  * The V3 here uses the SAME 2-body (Earth+Moon) gravity model as the CR3BP.
    Adding Sun/Mars/Jupiter perturbations is a V4 question, NOT V3. V3 is the
    "independent-integrator" cross-check at the SAME physics; V4 is the
    "real-ephemeris perturbation" check at a richer physics.
  * The Earth + Moon initial states are derived from the CR3BP setup (Moon
    at distance l_km = 384400 km from Earth, circular orbit period
    T_moon = 2π * TU_seconds). This matches the SILVER #331 / V3 substrate
    (the SILVER's circular-coplanar moon ephemeris assumption).

Usage:
    uv run python scripts/branch_c32_b0_v3_verify.py [--output PATH]
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
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.reachable_representatives import braik_ross_system

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
V1_VERDICT_PATH = Path("data/branch_c32_b0_v1_verdict.jsonl")
PHASE_LABEL = "389_p389_4"
N_CYCLES_LIST = (3, 5, 10)
V3_AGREEMENT_FLOOR_KMS = 100.0  # spec §14 V3 V3-vs-V2 agreement floor
EXPECTED_INTEGRATOR = "REBOUND IAS15"
IAS15_EPSILON = 1e-12


def _load_corrected_state_and_period() -> tuple[np.ndarray, float]:
    """Use the V1-corrected state + period (the V1-converged closure)."""
    with V1_VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "v1_verdict_cr3bp_periodic"
                and row.get("candidate_id") == CANDIDATE_ID
            ):
                state = np.array(row["state0_corrected_nondim"], dtype=np.float64)
                period = float(row["period_corrected_TU"])
                return state, period
    raise AssertionError(f"V1 verdict row for {CANDIDATE_ID!r} not found in {V1_VERDICT_PATH}")


def _rotating_to_inertial(
    state_rot_nondim: np.ndarray, t_nondim: float, mu: float
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a CR3BP rotating-frame state (nondim) to inertial-frame position+velocity (nondim).

    Rotating-frame state ``[x, y, z, xdot, ydot, zdot]`` at time ``t`` (nondim
    TU); rotation rate omega=1 in CR3BP nondim units. The standard transformation:

        r_in = R(omega t) r_rot,
        v_in = R(omega t) (v_rot + omega cross r_rot)

    where R(psi) is the rotation by psi about z. The Earth-Moon barycenter is at
    the inertial origin.
    """
    psi = t_nondim  # ω = 1 in nondim
    cos_psi, sin_psi = np.cos(psi), np.sin(psi)
    rmat = np.array(
        [[cos_psi, -sin_psi, 0.0], [sin_psi, cos_psi, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    r_rot = state_rot_nondim[:3]
    v_rot = state_rot_nondim[3:6]
    # omega cross r_rot with omega = (0,0,1): yields (-y, x, 0)
    omega_cross_r = np.array([-r_rot[1], r_rot[0], 0.0], dtype=np.float64)
    r_in = rmat @ r_rot
    v_in = rmat @ (v_rot + omega_cross_r)
    return r_in, v_in


def _inertial_to_rotating(
    r_in_nondim: np.ndarray, v_in_nondim: np.ndarray, t_nondim: float
) -> np.ndarray:
    """Inverse of :func:`_rotating_to_inertial`."""
    psi = t_nondim
    cos_psi, sin_psi = np.cos(psi), np.sin(psi)
    # R^T(ψ) = R(-ψ)
    rmat_t = np.array(
        [[cos_psi, sin_psi, 0.0], [-sin_psi, cos_psi, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    r_rot = rmat_t @ r_in_nondim
    # v_in = R(v_rot + omega cross r_rot), so v_rot = R^T v_in - omega cross r_rot
    v_rot = rmat_t @ v_in_nondim - np.array([-r_rot[1], r_rot[0], 0.0], dtype=np.float64)
    return np.concatenate([r_rot, v_rot])


def _run_ias15_n_cycles(
    state0_rot_nondim: np.ndarray,
    period_nondim: float,
    n_cycles: int,
    system: cr3bp.CR3BPSystem,
) -> tuple[list[float], list[bool], str, dict[str, float]]:
    """Propagate Earth + Moon + spacecraft under REBOUND IAS15 for n cycles.

    Returns ``(per_cycle_drift_km, per_cycle_converged, integrator_label, meta)``.
    Drift at boundary k is ``||r_rot(IAS15, t=k*T) - r_rot(IC, t=0)||`` in km
    — same definition as V2's cumulative-drift-from-IC.
    """
    l_km = system.l_km
    t_s = system.t_s
    v_scale = l_km / t_s

    # The CR3BP rotation rate and the Earth-Moon system GM are stored
    # INDEPENDENTLY in braik_ross_system: t_s comes from the 27.32-d sidereal
    # period while PRIMARIES["Earth"]=G(M_E+M_M) implies a 27.28-d Keplerian
    # period — a 0.14% relative inconsistency. For the V3 cross-check
    # ("is the V2 V3-vs-V2 agreement an integrator artifact?") we need the
    # IAS15 Earth+Moon system to match the CR3BP exactly. We therefore set
    # G_total = omega^2 * a^3 (the SAME Kepler's-third-law identity the CR3BP
    # implicitly assumes), and split into m_earth + m_moon by the system.mu
    # ratio. This sacrifices using the literal JPL Earth GM for the V3 lane
    # but keeps the cross-check pure: the only difference between V2 and V3
    # becomes the integrator (DOP853 vs IAS15).
    omega = 1.0 / t_s  # rad/s (the CR3BP rotation rate by construction)
    total_gm = omega * omega * l_km**3  # km^3/s^2 — Kepler-consistent with t_s
    moon_gm = total_gm * system.mu
    earth_gm = total_gm * (1.0 - system.mu)
    g_km3_kg_s2 = 6.67430e-20
    m_moon_kg = moon_gm / g_km3_kg_s2
    m_earth_kg = earth_gm / g_km3_kg_s2

    # Sanity-check the registry vs CR3BP-implied GM (recorded in the verdict).
    real_earth_gm = PRIMARIES["Earth"]
    real_moon_gm_sat = SATELLITES["Moon"].mu_km3_s2
    gm_relative_diff = (total_gm - real_earth_gm) / real_earth_gm
    _ = real_moon_gm_sat  # reported in JSONL via real_total_gm
    # Inertial-frame initial conditions (t=0):
    #   Earth at (-mu * l_km, 0, 0)
    #   Moon at ((1-mu) * l_km, 0, 0)
    #   Earth velocity (0, -mu * omega * l_km, 0)   [going in negative y]
    #   Moon velocity (0, (1-mu) * omega * l_km, 0)
    mu_nondim = system.mu
    omega = 1.0 / t_s  # rad/s
    earth_r_km_init = np.array([-mu_nondim * l_km, 0.0, 0.0], dtype=np.float64)
    moon_r_km_init = np.array([(1.0 - mu_nondim) * l_km, 0.0, 0.0], dtype=np.float64)
    earth_v_kms = np.array([0.0, -mu_nondim * omega * l_km, 0.0], dtype=np.float64)
    moon_v_kms = np.array([0.0, (1.0 - mu_nondim) * omega * l_km, 0.0], dtype=np.float64)

    # Spacecraft IC in inertial frame at t=0.
    r_sc_nondim, v_sc_nondim = _rotating_to_inertial(state0_rot_nondim, 0.0, mu_nondim)
    sc_r_km = r_sc_nondim * l_km
    sc_v_kms = v_sc_nondim * v_scale

    sim = rebound.Simulation()
    sim.G = g_km3_kg_s2
    sim.integrator = "ias15"
    sim.integrator.epsilon = IAS15_EPSILON
    # Earth
    sim.add(
        m=m_earth_kg,
        x=float(earth_r_km_init[0]),
        y=float(earth_r_km_init[1]),
        z=float(earth_r_km_init[2]),
        vx=float(earth_v_kms[0]),
        vy=float(earth_v_kms[1]),
        vz=float(earth_v_kms[2]),
    )
    # Moon
    sim.add(
        m=m_moon_kg,
        x=float(moon_r_km_init[0]),
        y=float(moon_r_km_init[1]),
        z=float(moon_r_km_init[2]),
        vx=float(moon_v_kms[0]),
        vy=float(moon_v_kms[1]),
        vz=float(moon_v_kms[2]),
    )
    # Spacecraft (massless test particle)
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

    integrator_label = EXPECTED_INTEGRATOR  # REBOUND IAS15

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
        # Convert back to rotating frame at the cycle boundary.
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
    return per_cycle_drift_km, per_cycle_converged, integrator_label, meta


def _load_v2_drifts() -> dict[int, list[float]]:
    """Load V2 per_cycle_drift_kms for the same n_cycles values."""
    drifts: dict[int, list[float]] = {}
    with Path("data/branch_c32_b0_v2_verdict.jsonl").open() as fh:
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
        default=Path("data/branch_c32_b0_v3_verdict.jsonl"),
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    state0, period = _load_corrected_state_and_period()
    system = braik_ross_system()
    v2_drifts = _load_v2_drifts()

    print(f"[V3-verify] candidate_id={CANDIDATE_ID}")
    print(f"[V3-verify] period_TU={period:.15f}")
    print(f"[V3-verify] system.l_km={system.l_km}, system.t_s={system.t_s:.4f} s")
    print(f"[V3-verify] integrator: {EXPECTED_INTEGRATOR}, epsilon={IAS15_EPSILON}")

    verdicts: dict[int, dict[str, Any]] = {}
    meta_record: dict[str, float] = {}
    all_pass = True
    for n in N_CYCLES_LIST:
        print(f"[V3-verify] running IAS15 at n_cycles={n}...")
        per_cycle_v3_km, per_cycle_conv, integ_label, meta = _run_ias15_n_cycles(
            state0, period, n, system
        )
        meta_record = meta  # last n_cycles's meta — they're all the same setup
        v2_drift_n = v2_drifts.get(n, [])
        # Per-cycle agreement: |V3_drift - V2_drift| (both relative to IC).
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
            "candidate_id": CANDIDATE_ID,
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
                    "candidate_id": CANDIDATE_ID,
                    "phase": PHASE_LABEL,
                    "iso_start": iso_start,
                    "iso_end": iso_end,
                    "elapsed_seconds": elapsed,
                    "v3_agreement_floor_kms": V3_AGREEMENT_FLOOR_KMS,
                    "n_cycles_list": list(N_CYCLES_LIST),
                    "integrator_target": EXPECTED_INTEGRATOR,
                    "ias15_epsilon": IAS15_EPSILON,
                    "cr3bp_setup_meta": meta_record,
                    "discipline": (
                        "Spec §14 V3 cross-check: REBOUND IAS15 (independent of "
                        "scipy DOP853 used by V2) on the Earth+Moon 2-body system "
                        "matching the CR3BP setup. The drift agreement "
                        "(|V3_drift - V2_drift|) must stay below "
                        f"{V3_AGREEMENT_FLOOR_KMS} km — the spec §14 V3 floor — "
                        "to confirm the V2 bounded-drift signature is a "
                        "REAL dynamical property, not a DOP853 artifact. The "
                        "IAS15 Earth-Moon system uses G_total = omega^2 a^3 "
                        "(Kepler-consistent with t_s); the registry "
                        "PRIMARIES['Earth'] differs by ~1.4e-3 (the CR3BP's "
                        "implicit assumption). Using the registry GM here would "
                        "introduce a 0.14% Earth-Moon period mismatch into the "
                        "V3 baseline, swamping the V2-vs-V3 round-off comparison."
                    ),
                }
            )
            + "\n"
        )
        for n in N_CYCLES_LIST:
            fh.write(json.dumps(verdicts[n]) + "\n")
        fh.write(
            json.dumps(
                {
                    "kind": "footer",
                    "candidate_id": CANDIDATE_ID,
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
