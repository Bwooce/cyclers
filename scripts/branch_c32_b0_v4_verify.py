"""#389 P389.5 — V4 real-ephemeris (DE440) cross-check for branch_C32_b0.

V4 question (spec §14 V4): does the V3 bounded-drift signature survive an
*independent codebase + real ephemeris* n-body realization? This is the
Earth-Moon analog of the #335 V4-strict SILVER check (which used URA111 for the
Uranian moons). Here we use the JPL **DE440** planetary ephemeris (the modern
successor to the DE421 shipped with GMAT) to seed a full REBOUND IAS15 n-body
simulation:

  * Earth (399) + Moon (301) at their DE440 inertial state relative to the
    Earth-Moon barycenter at the launch epoch, with the JPL registry GM masses
    (so the Earth-Moon two-body is the REAL eccentric/inclined lunar orbit, not
    the CR3BP circular idealization V3 used);
  * Sun (10) + Mars barycenter (4) + Jupiter barycenter (5) added as massive
    third-body perturbers, also seeded from DE440 — the perturbations the CR3BP
    and the V3 2-body model both omit;
  * the spacecraft as a massless test particle, its CR3BP rotating-frame IC
    converted to the inertial frame and then *registered onto the real DE440
    Earth-Moon geometry* (rotated so the CR3BP +x axis — the Earth->Moon line —
    aligns with the real Earth->Moon direction at the epoch, and the CR3BP
    angular-momentum axis aligns with the real Earth-Moon orbit normal).

We then integrate for n_cycles CR3BP periods, sample at each cycle boundary,
convert the spacecraft state back to the (instantaneous, DE440-defined)
rotating frame, and report drift = ||r_rot(V4) - r_rot(IC)|| (the same
cumulative-drift-from-IC metric V2/V3 used) plus the V4-vs-V3 agreement.

PASS criterion (spec §14 V4): the V4-vs-V3 drift agreement stays below the
same 50 000 km same-model agreement floor the #335 SILVER V4-strict used
(``driver_floors.agreement_floor_kms = 50000.0``). For a deeply stable orbit
like branch_C32_b0 we expect the *real* lunar eccentricity + solar/planetary
tides to introduce a genuine — but bounded — perturbation, far below that floor.

Discipline notes:

  * V3 isolated the INTEGRATOR (IAS15 vs DOP853) at the SAME idealized physics.
    V4 isolates the PHYSICS (real DE440 ephemeris + Sun/Mars/Jupiter tides) at
    an independent integrator/codebase. The two together close the
    "is the bounded drift real?" question.
  * DE440 is fetched + cached by astropy; the kernel identity + the seeded
    body states are recorded in the verdict for full reproducibility.
  * The rotating-frame at each cycle boundary is defined by the INSTANTANEOUS
    DE440 Earth->Moon direction at that time (a real, breathing, inclined
    rotating frame), not the constant-rate CR3BP frame — so the drift metric
    is honest about the real geometry the spacecraft sees.

Usage:
    uv run python scripts/branch_c32_b0_v4_verify.py [--output PATH]
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

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
V1_VERDICT_PATH = Path("data/branch_c32_b0_v1_verdict.jsonl")
V3_VERDICT_PATH = Path("data/branch_c32_b0_v3_verdict.jsonl")
PHASE_LABEL = "389_p389_5"
N_CYCLES_LIST = (3, 5, 10)
# Spec §14 V4 / #335 SILVER V4-strict same-model agreement floor.
V4_AGREEMENT_FLOOR_KMS = 50000.0
EXPECTED_INTEGRATOR = "REBOUND IAS15"
IAS15_EPSILON = 1e-12
EPHEMERIS = "de440"
# Default launch epoch — matches the #335 SILVER's first scan epoch.
DEFAULT_EPOCH_UTC = "2000-01-15T00:00:00"
G_KM3_KG_S2 = 6.67430e-20

# Bodies seeded as massive REBOUND particles (NAIF/astropy names + GM km^3/s^2).
# Earth + Moon use the JPL registry GM (the REAL split, not the CR3BP mu split).
_GM_KM3_S2 = {
    "earth": float(PRIMARIES["Earth"]) - float(SATELLITES["Moon"].mu_km3_s2),  # Earth body alone
    "moon": float(SATELLITES["Moon"].mu_km3_s2),
    "sun": 1.32712440018e11,
    "mars": 4.282837e4,  # Mars-system GM (barycenter)
    "jupiter": 1.26686534e8,  # Jupiter-system GM (barycenter)
}
_PERTURBER_BODIES = ("sun", "mars", "jupiter")


def _load_corrected_state_and_period() -> tuple[np.ndarray, float]:
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
    """Quantify why a far-amplitude EM orbit is or is not stable under real solar tides.

    Computes the orbit's max Earth-distance amplitude as a fraction of the
    Earth-Sun Hill radius, and the solar tidal acceleration at that amplitude as
    a fraction of Earth's gravity. A near-unity Hill fraction / large tide ratio
    is the structural signature of an orbit the CR3BP (Sun-free) cannot model.
    """
    l_km = system.l_km
    r_amp_km = float(np.linalg.norm(state0_rot_nondim[:3]) * l_km)
    gm_earth_total = float(PRIMARIES["Earth"])  # G(M_E + M_M), km^3/s^2
    gm_sun = 1.32712440018e11  # km^3/s^2
    a_earth_sun_km = 1.495978707e8  # 1 AU
    # Earth-Sun Hill radius R_H = a (mu_E / (3 mu_Sun))^(1/3).
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
    """CR3BP rotating-frame state (nondim) -> inertial position+velocity (nondim)."""
    psi = t_nondim  # omega = 1 in nondim
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
    """DE440 inertial position (km) + velocity (km/s) of ``body`` relative to the
    Earth-Moon barycenter at time ``t``."""
    emb_p, emb_v = get_body_barycentric_posvel("earth-moon-barycenter", t)
    p, v = get_body_barycentric_posvel(body, t)
    r = (p.xyz - emb_p.xyz).to(u.km).value
    vel = (v.xyz - emb_v.xyz).to(u.km / u.s).value
    return np.asarray(r, dtype=np.float64), np.asarray(vel, dtype=np.float64)


def _real_rotating_basis(t: Time) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Instantaneous DE440 rotating-frame basis at ``t``.

    ``xhat`` = unit Earth->Moon direction; ``zhat`` = unit Earth-Moon orbital
    angular-momentum direction; ``yhat`` = zhat x xhat. Returns ``(xhat, yhat, zhat)``.
    """
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
    """Propagate the DE440-seeded n-body + spacecraft over ``bodies``.

    ``bodies`` selects which massive bodies are seeded from DE440 (the full V4
    model is Earth+Moon+Sun+Mars+Jupiter; pass ("earth", "moon") for the
    no-perturber control). Returns ``(per_cycle_drift_km, per_cycle_converged,
    meta)``. Drift at boundary k is ``||r_rot(V4, t=k*T) - r_rot(IC, t=0)||`` in
    km (same definition as V2/V3).
    """
    l_km = system.l_km
    t_s = system.t_s
    v_scale = l_km / t_s
    mu_nondim = system.mu

    solar_system_ephemeris.set(EPHEMERIS)
    t0 = Time(epoch_utc, scale="tdb")

    # --- Real DE440 rotating-frame basis at the launch epoch ---
    xhat0, yhat0, zhat0 = _real_rotating_basis(t0)
    rot0 = np.column_stack([xhat0, yhat0, zhat0])  # columns are the rotating axes

    # --- Spacecraft IC: CR3BP rotating -> CR3BP inertial -> register onto real frame ---
    # In the CR3BP, at t=0 the rotating and inertial frames coincide and the
    # +x axis is the Earth->Moon line, +z the orbit normal. We map the CR3BP
    # inertial state (nondim, scaled by l_km / v_scale) onto the real DE440
    # Earth-Moon basis: r_real = rot0 @ r_cr3bp_inertial.
    r_sc_nd, v_sc_nd = _rotating_to_inertial(state0_rot_nondim, 0.0, mu_nondim)
    r_sc_cr3bp_km = r_sc_nd * l_km
    v_sc_cr3bp_kms = v_sc_nd * v_scale
    sc_r_km = rot0 @ r_sc_cr3bp_km
    sc_v_kms = rot0 @ v_sc_cr3bp_kms

    # --- Seed REBOUND with DE440 massive bodies ---
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
    # Spacecraft (massless test particle), index == number of massive bodies.
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

    # IC position expressed in the launch-epoch rotating frame (== CR3BP r at t=0).
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
        # Instantaneous DE440 rotating frame at this boundary.
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
        default=Path("data/branch_c32_b0_v4_verdict.jsonl"),
    )
    parser.add_argument("--epoch", type=str, default=DEFAULT_EPOCH_UTC)
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    state0, period = _load_corrected_state_and_period()
    system = braik_ross_system()
    v3_drifts = _load_v3_drifts()

    print(f"[V4-verify] candidate_id={CANDIDATE_ID}")
    print(f"[V4-verify] period_TU={period:.15f}")
    print(f"[V4-verify] ephemeris={EPHEMERIS}, epoch={args.epoch}")
    print(f"[V4-verify] integrator: {EXPECTED_INTEGRATOR}, epsilon={IAS15_EPSILON}")

    verdicts: dict[int, dict[str, Any]] = {}
    meta_record: dict[str, Any] = {}
    all_pass = True
    for n in N_CYCLES_LIST:
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
            "candidate_id": CANDIDATE_ID,
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

    # --- Earth+Moon-only control (no Sun/Mars/Jupiter): isolates the perturber
    # contribution. If the EM-only drift is small but the full-model drift is
    # huge, the failure is attributable to the real solar/planetary tides, not a
    # seeding/registration artifact. ---
    n_ctrl = max(N_CYCLES_LIST)
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
                    "candidate_id": CANDIDATE_ID,
                    "phase": PHASE_LABEL,
                    "iso_start": iso_start,
                    "iso_end": iso_end,
                    "elapsed_seconds": elapsed,
                    "v4_agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                    "n_cycles_list": list(N_CYCLES_LIST),
                    "integrator_target": EXPECTED_INTEGRATOR,
                    "ias15_epsilon": IAS15_EPSILON,
                    "ephemeris": EPHEMERIS,
                    "epoch_utc": args.epoch,
                    "v4_setup_meta": meta_record,
                    "discipline": (
                        "Spec §14 V4 real-ephemeris cross-check: REBOUND IAS15 "
                        "n-body seeded from JPL DE440 (Earth + Moon + Sun + Mars "
                        "+ Jupiter as massive bodies; spacecraft a massless test "
                        "particle). Unlike V3 (which isolated the INTEGRATOR at "
                        "the idealized CR3BP 2-body), V4 isolates the PHYSICS: "
                        "real eccentric/inclined lunar orbit + solar/Mars/Jupiter "
                        "tides. The V4-vs-V3 drift agreement must stay below the "
                        f"{V4_AGREEMENT_FLOOR_KMS} km same-model floor (the #335 "
                        "SILVER V4-strict driver floor) for the bounded-drift "
                        "signature to be confirmed REAL under real ephemerides. "
                        "The drift metric uses the INSTANTANEOUS DE440 Earth->Moon "
                        "rotating frame at each cycle boundary."
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
                    "kind": "control_em_only",
                    "candidate_id": CANDIDATE_ID,
                    "bodies": ["earth", "moon"],
                    "n_cycles_requested": int(n_ctrl),
                    "per_cycle_drift_kms": [float(x) for x in ctrl_drift_km],
                    "max_drift_kms": float(max(ctrl_drift_km)) if ctrl_drift_km else None,
                    "converged_at_each_cycle": bool(all(ctrl_conv)),
                    "note": (
                        "Earth+Moon-only (DE440 real eccentric/inclined lunar orbit, "
                        "NO Sun/Mars/Jupiter). Isolates the perturber contribution: a "
                        "small EM-only drift vs a huge full-model drift attributes the "
                        "V4 result to real solar/planetary tides, not a seeding artifact."
                    ),
                }
            )
            + "\n"
        )
        fh.write(
            json.dumps(
                {
                    "kind": "structural_diagnostic",
                    "candidate_id": CANDIDATE_ID,
                    **diag,
                    "note": (
                        "branch_C32_b0 is a far-amplitude EM orbit whose max "
                        "Earth-distance reaches a large fraction of the Earth-Sun "
                        "Hill radius, where the solar tidal acceleration is a "
                        "non-negligible fraction of Earth's gravity. The CR3BP (and "
                        "V1-V3, which live in or near it) ignore the Sun entirely; "
                        "V4 with real DE440 solar gravity is the first gate that sees "
                        "this perturbation."
                    ),
                }
            )
            + "\n"
        )
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
    print(f"[V4-verify] verdict written to {args.output}")
    print(f"[V4-verify] all V4 gates pass: {all_pass}")


if __name__ == "__main__":
    main()
