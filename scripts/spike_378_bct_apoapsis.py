"""#378 Phase 0.2 — BCR4BP apoapsis-reach spike (Hiten-signature gate).

Belbruno-style honesty: measure before building. This spike decides (design
draft R1 / plan Phase 0.2) whether the project's *incoherent* BCR4BP
(`core/bcr4bp.propagate_bcr4bp`) can shape a Sun-driven high apoapsis at all
— the geometric vehicle a Hiten-class exterior WSB transfer requires
(Q_a ~ 1.5e6 km ~ 3.9 lunar distances from Earth).

It asserts NOTHING (mirrors the #412 reach-spike style): it sweeps a coarse
grid of Earth-departure states and prints a table of max Earth-relative
apoapsis reached. The GATE is read off the printed table:

  * PASS  -- at least one (t0, |V0|) reaches max apoapsis >~ 3 LD (~1.1e6 km).
  * KILL  -- max apoapsis < 2 LD across the whole sweep (R1 fires; BCT
            construction blocks on the coherent-QBCP acquisition and only
            `core/wsb.py` ships).

Units: BCR4BP is nondimensional with the Earth-Moon distance = 1, so
1 LD = 1 nondim length unit. Earth sits at (-mu, 0, 0). We launch from a
low-Earth periapsis Q0 (r_13 = r_E + 200 km in nondim units) on the
Earth-Sun line side, sweep the departure speed |V0| (km/s, converted to
nondim), the flight-path angle gamma0 (angle of V0 off the local-horizontal),
and the Sun phase t0, then forward-propagate ~120 days and record the maximum
distance from Earth.

Run:  uv run python scripts/spike_378_bct_apoapsis.py
"""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp

# Physical scales (Earth-Moon). EM mean distance and the EM-system time unit.
EM_DISTANCE_KM = 384_400.0  # nondim length unit (1 LD)
EARTH_RADIUS_KM = 6378.137
LEO_PERIAPSIS_KM = EARTH_RADIUS_KM + 200.0  # r_E + 200 km (Belbruno capture problem)
# EM synodic time unit: TU = 1 / n_moon, with n_moon the sidereal mean motion.
# sidereal month ~ 27.32166 d => angular rate 2*pi / that; TU = period / (2*pi).
SIDEREAL_MONTH_DAYS = 27.321661
TU_DAYS = SIDEREAL_MONTH_DAYS / (2.0 * math.pi)
KM_PER_TU_OVER_S = EM_DISTANCE_KM / (TU_DAYS * 86400.0)  # km/s per nondim velocity


def kms_to_nondim_v(v_kms: float) -> float:
    """Convert a speed in km/s to BCR4BP nondim velocity units."""
    return v_kms / KM_PER_TU_OVER_S


def earth_periapsis_state(mu: float, v_kms: float, gamma_deg: float) -> np.ndarray:
    """Build a 6-state at a low-Earth periapsis on the +x (anti-Sun-ish) side.

    Earth at (-mu, 0, 0). We place Q0 at r_13 = LEO_PERIAPSIS_KM beyond Earth
    on the +x axis, with velocity in the +y direction rotated by gamma_deg
    toward radial (flight-path angle off the local horizontal).
    """
    r13_nd = LEO_PERIAPSIS_KM / EM_DISTANCE_KM
    x = -mu + r13_nd
    y = 0.0
    z = 0.0
    v_nd = kms_to_nondim_v(v_kms)
    gamma = math.radians(gamma_deg)
    # Local horizontal at +x periapsis is +y; radial outward is +x.
    vx = v_nd * math.sin(gamma)
    vy = v_nd * math.cos(gamma)
    vz = 0.0
    return np.array([x, y, z, vx, vy, vz], dtype=np.float64)


def max_earth_apoapsis_ld(
    system: bcr4bp.BCR4BPSystem,
    state0: np.ndarray,
    t_days: float,
    n_samples: int = 600,
) -> float:
    """Forward-propagate and return max distance-from-Earth in LD (= nondim)."""
    mu = system.mu
    t_nd = t_days / TU_DAYS
    earth = np.array([-mu, 0.0, 0.0])
    max_r = 0.0
    # March in chunks so a blow-up (close approach) is caught and skipped.
    n_chunk = n_samples
    times = np.linspace(t_nd / n_chunk, t_nd, n_chunk)
    cur = state0.copy()
    t_prev = 0.0
    for tnext in times:
        try:
            arc = bcr4bp.propagate_bcr4bp(
                system, cur, tnext - t_prev, t0=t_prev, rtol=1e-10, atol=1e-10
            )
        except RuntimeError:
            break
        cur = arc.state_f
        t_prev = tnext
        r = float(np.linalg.norm(cur[:3] - earth))
        if not math.isfinite(r):
            break
        max_r = max(max_r, r)
    return max_r


def main() -> None:
    system = bcr4bp.andreu_default()
    mu = system.mu
    sun_period_days = system.sun_period_tu * TU_DAYS

    # Sweep grid. |V0| spans the band that, from a 200 km LEO, gives apoapsis
    # from EM-distance scale up past the Sun-shaped 3-4 LD target. Circular LEO
    # speed ~7.78 km/s; escape ~11.0 km/s. The WSB exterior transfer departs
    # just below escape so the apoapsis is huge but still bound to Earth.
    v_grid = [10.70, 10.80, 10.85, 10.90, 10.95, 10.98, 11.00, 11.02]
    gamma_grid = [0.0, 2.0, 5.0]
    # Sun phase t0 maps to theta_sun0; sweep across one synodic month.
    t0_phase_grid = [0.0, 0.25, 0.5, 0.75]  # fraction of Sun synodic period
    tof_days = 140.0

    print("#378 Phase 0.2 BCR4BP apoapsis-reach spike")
    print(
        f"  TU = {TU_DAYS:.6f} d, 1 LD = {EM_DISTANCE_KM:.0f} km, "
        f"km/s per nondim v = {KM_PER_TU_OVER_S:.5f}"
    )
    print(f"  Sun synodic period = {sun_period_days:.3f} d, TOF = {tof_days:.0f} d")
    print(
        f"  LEO periapsis = {LEO_PERIAPSIS_KM:.1f} km "
        f"(r13_nd = {LEO_PERIAPSIS_KM / EM_DISTANCE_KM:.6f})"
    )
    print()
    print(
        f"  {'|V0| km/s':>10} {'gamma deg':>10} {'t0/Tsun':>8} "
        f"{'max_apo LD':>11} {'max_apo km':>13}"
    )
    print("  " + "-" * 56)

    best = 0.0
    best_cfg = None
    for v_kms in v_grid:
        for gamma in gamma_grid:
            for frac in t0_phase_grid:
                theta0 = frac * 2.0 * math.pi
                sys_t = bcr4bp.BCR4BPSystem(
                    mu=system.mu,
                    mu_sun=system.mu_sun,
                    a_sun_nondim=system.a_sun_nondim,
                    omega_sun_nondim=system.omega_sun_nondim,
                    theta_sun0=theta0,
                )
                st = earth_periapsis_state(mu, v_kms, gamma)
                apo_ld = max_earth_apoapsis_ld(sys_t, st, tof_days)
                apo_km = apo_ld * EM_DISTANCE_KM
                print(
                    f"  {v_kms:>10.2f} {gamma:>10.1f} {frac:>8.2f} {apo_ld:>11.3f} {apo_km:>13.0f}"
                )
                if apo_ld > best:
                    best = apo_ld
                    best_cfg = (v_kms, gamma, frac)

    print("  " + "-" * 56)
    print(
        f"  BEST max apoapsis = {best:.3f} LD ({best * EM_DISTANCE_KM:.0f} km) "
        f"at (|V0|, gamma, t0/Tsun) = {best_cfg}"
    )
    print()
    if best >= 3.0:
        print(
            "  GATE: PASS (>= 3 LD) -- BCR4BP can shape a Hiten-signature apoapsis. "
            "Proceed with Phase 2/3 BCT construction."
        )
    elif best >= 2.0:
        print(
            "  GATE: MARGINAL (2-3 LD) -- apoapsis reachable but below the 3.9 LD "
            "Hiten target; proceed with caution, document the gap."
        )
    else:
        print(
            "  GATE: KILL (< 2 LD) -- R1 fires. The incoherent BCR4BP cannot hold "
            "a Hiten-signature apoapsis; BCT construction blocks on the coherent "
            "QBCP acquisition. Only core/wsb.py ships."
        )


if __name__ == "__main__":
    main()
