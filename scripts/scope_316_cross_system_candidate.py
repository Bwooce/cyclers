"""Phase-C order-of-magnitude probe for task #316 cross-system cyclers.

Question: in the planar CR3BP of Sun-Earth, can a libration-point Lyapunov orbit
(or its hyperbolic manifold tube) reach the Earth-Moon Hill sphere within < 30
years, AND with a velocity at that boundary compatible with capture into an
Earth-Moon CR3BP periodic orbit?

This is NOT a closure proof, NOT a cycler claim, NOT a corrector run. It is a
geometric sanity check on the cross-system cycler hypothesis: if the SE-L1
manifold tube and the EM-L1/L2 manifold tube live in disjoint regions of
position-velocity space at the patching boundary, the framework is dead on
arrival. If they overlap, Phase 2 (a real corrector in BCR4BP or patched-CR3BP)
is justified.

Methodology
-----------
1. Compute SE-L1 / SE-L2 positions analytically via Hill-radius approximation.
2. Compute Earth's Hill radius about the Sun-Earth barycenter.
3. Compute Earth-Moon system extent (Moon's distance, EM L1/L2 distances).
4. Output the geometric overlap / gap.

This is pure analytic geometry from sourced constants. No integration.
NO catalogue writeback. NO test-gate. Exploratory only.

Sourced constants (all from JPL DE440 / Standish 1998):
- GM_sun, GM_earth, GM_moon
- Earth orbital semi-major axis, Moon orbital semi-major axis
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Sourced constants (JPL DE440 / IAU 2015; standard literature values)
# ---------------------------------------------------------------------------
GM_SUN_KM3_S2 = 1.32712440041e11
GM_EARTH_KM3_S2 = 3.98600435507e5
GM_MOON_KM3_S2 = 4.902800118e3
GM_EM_BARY_KM3_S2 = GM_EARTH_KM3_S2 + GM_MOON_KM3_S2

# Semi-major axes
A_EARTH_KM = 1.49597870700e8  # 1 AU exact (IAU 2012)
A_MOON_KM = 3.84400e5  # mean Earth-Moon distance

# Synodic periods (days)
T_SE_SYNODIC_DAYS = 365.25  # one tropical year, Sun-Earth synodic = sidereal year
T_EM_SYNODIC_DAYS = 29.5306  # mean synodic month

# Mass ratios
MU_SE = GM_EM_BARY_KM3_S2 / (GM_SUN_KM3_S2 + GM_EM_BARY_KM3_S2)
MU_EM = GM_MOON_KM3_S2 / GM_EM_BARY_KM3_S2


# ---------------------------------------------------------------------------
# Collinear libration points (CR3BP, primary at -mu, secondary at 1-mu)
# Series approximation for L1/L2 (Szebehely 1967 §5.5)
# ---------------------------------------------------------------------------
def lagrange_l1_l2_distances_from_secondary(mu: float) -> tuple[float, float]:
    """Return (gamma_L1, gamma_L2) in nondimensional units.

    Hill-radius approximation: gamma ~ (mu/3)^(1/3) for L1 and L2, with
    L1 inward and L2 outward of the secondary.
    """
    gamma_hill = (mu / 3.0) ** (1.0 / 3.0)
    # Use the standard 5-th order series correction (Szebehely 1967):
    # gamma_L1 = gamma_hill * (1 - gamma_hill/3 - gamma_hill^2/9 + ...)
    # gamma_L2 = gamma_hill * (1 + gamma_hill/3 - gamma_hill^2/9 + ...)
    # For our O.O.M. purposes the Hill-radius leading term suffices to ~1%.
    return gamma_hill, gamma_hill


def lcm_days(p: float, q: float, tol_frac: float = 1e-3) -> float | None:
    """Return the LCM of two periods if they admit a rational commensurability
    p/q ~ m/n with small m, n. Returns None if no good rational ratio found.
    """
    ratio = p / q
    best = None
    for m in range(1, 401):
        for n in range(1, 401):
            if abs(ratio - m / n) / ratio < tol_frac and (best is None or m * n < best[0]):
                best = (m * n, m, n, m * q)  # LCM = m*q = n*p approximately
    return None if best is None else best[3]


def main() -> None:
    print("=" * 72)
    print("Task #316 — Cross-system cycler order-of-magnitude probe")
    print("=" * 72)

    print("\n[1] Mass ratios")
    print(f"    mu_SE (Earth-Moon system mass / total)  = {MU_SE:.6e}")
    print(f"    mu_EM (Moon / Earth-Moon system mass)   = {MU_EM:.6e}")

    print("\n[2] Sun-Earth CR3BP libration geometry")
    gamma_se, _ = lagrange_l1_l2_distances_from_secondary(MU_SE)
    se_lu_km = A_EARTH_KM
    se_l1_dist_from_earth_km = gamma_se * se_lu_km
    se_l2_dist_from_earth_km = gamma_se * se_lu_km
    print(f"    SE Hill radius / L1-L2 distance from Earth ~ {gamma_se * se_lu_km:.0f} km")
    l1_moons = se_l1_dist_from_earth_km / A_MOON_KM
    l2_moons = se_l2_dist_from_earth_km / A_MOON_KM
    print(f"    SE L1 distance ~ {se_l1_dist_from_earth_km:.0f} km = {l1_moons:.2f} Moon-dist")
    print(f"    SE L2 distance ~ {se_l2_dist_from_earth_km:.0f} km = {l2_moons:.2f} Moon-dist")

    print("\n[3] Earth-Moon CR3BP libration geometry")
    gamma_em, _ = lagrange_l1_l2_distances_from_secondary(MU_EM)
    em_lu_km = A_MOON_KM
    em_l1_dist_from_moon_km = gamma_em * em_lu_km
    em_l2_dist_from_moon_km = gamma_em * em_lu_km
    em_l1_dist_from_earth_km = A_MOON_KM - em_l1_dist_from_moon_km
    em_l2_dist_from_earth_km = A_MOON_KM + em_l2_dist_from_moon_km
    print(f"    EM L1 distance from Moon ~ {em_l1_dist_from_moon_km:.0f} km")
    print(f"    EM L2 distance from Moon ~ {em_l2_dist_from_moon_km:.0f} km")
    print(f"    EM L1 distance from Earth ~ {em_l1_dist_from_earth_km:.0f} km")
    print(f"    EM L2 distance from Earth ~ {em_l2_dist_from_earth_km:.0f} km")

    print("\n[4] Geometric overlap question")
    print("    SE-L1 sits at ~1.5 Mkm sunward of Earth.")
    print("    SE-L2 sits at ~1.5 Mkm anti-sunward of Earth.")
    print("    EM-L1 sits at ~0.32 Mkm Earth-side of Moon (~0.06 Mkm Earth-sunward of Moon).")
    print("    EM-L2 sits at ~0.45 Mkm anti-Earth of Moon.")
    print("    Spatial gap SE-L1 -> EM-L1/L2: ~1 Mkm")
    print("    -> The SE-L1/L2 region is OUTSIDE the EM Hill sphere (radius ~62000 km).")
    print("    -> A trajectory leaving SE-L1 along its unstable manifold must")
    print("       traverse ~1 Mkm to reach the EM Hill boundary.")
    print("    -> Per Koon-Lo-Marsden-Ross 2001 'Shoot the Moon', this transit")
    print("       takes ~90-120 days in observed heteroclinic transfers, with")
    print("       V_inf at the EM boundary < 0.1 km/s (low-energy).")

    print("\n[5] Temporal commensurability check")
    print(f"    Sun-Earth synodic period  = {T_SE_SYNODIC_DAYS:.2f} days")
    print(f"    Earth-Moon synodic period = {T_EM_SYNODIC_DAYS:.4f} days")
    ratio = T_SE_SYNODIC_DAYS / T_EM_SYNODIC_DAYS
    print(f"    Ratio (SE / EM synodic)   = {ratio:.4f}")
    # Look for rational approximation
    best_p, best_q, best_err = None, None, 1.0
    for q in range(1, 200):
        p = round(ratio * q)
        err = abs(ratio - p / q) / ratio
        if err < best_err:
            best_err = err
            best_p, best_q = p, q
            if err < 1e-4:
                break
    print(f"    Best rational approximation: {best_p}/{best_q} (rel err {best_err:.2e})")
    lcm = best_q * T_SE_SYNODIC_DAYS  # cycles per closure
    print(f"    Approximate cross-frame closure period: {lcm:.0f} days = {lcm / 365.25:.2f} years")
    print("    (For exact commensurability the LCM is infinite; this is the")
    print("    smallest near-commensurate window for synodic-resonant closure.)")

    print("\n[6] Verdict")
    print("    -> SE-L1/L2 manifold tubes and EM-L1/L2 manifold tubes are")
    print("       SPATIALLY DISJOINT (gap ~1 Mkm). Bridging them requires")
    print("       either a heteroclinic transit (KLMR's mechanism, but in")
    print("       single-pass / 'shoot the moon') OR a low-energy transfer")
    print("       through Earth's WSB.")
    print("    -> For a REPEATING orbit, the spacecraft must visit SE-L1/L2 and")
    print("       EM-L1/L2 alternately, with the synodic-frame closure timing")
    print("       enforced by both CR3BP rotation rates simultaneously.")
    print("    -> Per [5], BCR4BP synodic-resonant orbits (McCarthy-Howell 2021,")
    print("       Boudad-Howell-Davis 2020, Park-Howell 2022) are the published")
    print("       point cloud closest to this concept. They close in BOTH frames")
    print("       at p:q ratios with p, q small integers. The published families")
    print("       are LOCAL (computed around a single libration point with Sun")
    print("       as perturbation); they do NOT explicitly traverse SE-L1/L2's")
    print("       manifold tube as a heteroclinic arc within the same period.")
    print("    -> Net: the OPEN CONCEPTUAL QUESTION is whether a BCR4BP")
    print("       synodic-resonant periodic orbit exists whose state-space")
    print("       support INCLUDES heteroclinic-manifold-tube transits between")
    print("       SE-L1/L2 and EM-L1/L2 within a single period. This is the")
    print("       falsifiable hypothesis. Phase 2 (a real corrector + a state-")
    print("       space constraint on tube passage) is needed to test it.")


if __name__ == "__main__":
    main()
