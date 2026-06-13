"""Reproduce the Tito et al. (IEEE Aerospace 2013) 2018 Mars free-return.

Cross-check our real-ephemeris Lambert/flyby stack against the published
Earth-Mars-Earth ballistic free-return in Tables III/IV of:

    Tito, Anderson, Carrico, et al., "Feasibility Analysis for a Manned
    Mars Free-Return Mission in 2018," IEEE Aerospace Conference, 2013.

EXPECTED side = the Astrogator (full-force, JPL DE421) values transcribed
verbatim from Tables III/IV. ACTUAL side = our DE440 Ephemeris + Lambert
solve on the SAME published epochs. We acknowledge (not fix) the DE421 vs
DE440 modelling difference.

Method (REPRODUCE-BEFORE-WIRING — published numbers are the target):
  1. Convert the three published UTCG epochs to our t_sec (TDB s since J2000).
  2. Heliocentric Earth state at departure, Mars state at flyby, Earth state
     at return from Ephemeris("astropy") (DE440).
  3. Leg 1: Lambert(r_Earth_dep -> r_Mars_flyby, tof1). v_inf_dep = |v1 -
     v_Earth_dep|; C3 = v_inf_dep**2. v_inf_arr_Mars = |v2 - v_Mars_flyby|.
  4. Leg 2: Lambert(r_Mars_flyby -> r_Earth_ret, tof2). v_inf_dep_Mars =
     |v1 - v_Mars_flyby|; v_inf_arr_Earth = |v2 - v_Earth_ret|; return C3.
  5. Flyby ballistic-continuity check at Mars: speed match of the two Mars
     v_inf vectors and the required bend vs the achievable cone at ~100 km.
  6. Report residuals against Tables III/IV.

Run: uv run python scripts/tito_free_return_repro.py
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import bend_angle, max_bend
from cyclerfinder.core.lambert import lambert

# --------------------------------------------------------------------------- #
# Published values (Tito et al., IEEE Aerospace 2013), Tables III and IV.
# Astrogator full-force solution, JPL DE421. Epochs are UTCG.
# --------------------------------------------------------------------------- #
DEPART_EARTH_UTC = "2018-01-05 07:00:00.000"
FLYBY_MARS_UTC = "2018-08-20 08:18:19.619"
RETURN_EARTH_UTC = "2019-05-21 13:52:48.012"

# Table III flight times (days).
TOF1_DAYS_PUB = 227.05439374  # Earth -> Mars
TOF2_DAYS_PUB = 274.23227306  # Mars -> Earth
TOTAL_DAYS_PUB = 501.2866668

# Table IV (km/s and km^2/s^2). Departure / arrival per leg.
# Leg 1: Earth departure, Mars arrival.
VINF_DEP_EARTH_PUB = 6.232
C3_DEP_EARTH_PUB = 38.835
VINF_ARR_MARS_PUB = 5.417
C3_ARR_MARS_PUB = 29.344
# Leg 2: Mars departure, Earth arrival.
VINF_DEP_MARS_PUB = 5.417
VINF_ARR_EARTH_PUB = 8.837
C3_ARR_EARTH_PUB = 78.094

# Fig.3 text: Mars periapsis ~100 km altitude, peri speed ~7.27 km/s.
MARS_PERIAPSIS_ALT_KM_PUB = 100.0


def utc_to_tsec(utc_str: str) -> float:
    """UTCG calendar string -> our t_sec (TDB seconds since J2000(TDB))."""
    from astropy.time import Time

    t = Time(utc_str, scale="utc")
    return float((t.tdb.jd - 2451545.0) * 86400.0)


def main() -> None:
    eph = Ephemeris("astropy")  # DE440

    t_dep = utc_to_tsec(DEPART_EARTH_UTC)
    t_fly = utc_to_tsec(FLYBY_MARS_UTC)
    t_ret = utc_to_tsec(RETURN_EARTH_UTC)

    tof1 = t_fly - t_dep
    tof2 = t_ret - t_fly
    tof1_days = tof1 / 86400.0
    tof2_days = tof2 / 86400.0

    r_e_dep, v_e_dep = eph.state("E", t_dep)
    r_m_fly, v_m_fly = eph.state("M", t_fly)
    r_e_ret, v_e_ret = eph.state("E", t_ret)

    # --- Leg 1: Earth -> Mars (short-way prograde, single rev). --- #
    sols1 = lambert(r_e_dep, r_m_fly, tof1, prograde=True, max_revs=0)
    s1 = sols1[0]
    vinf_dep_vec = s1.v1 - v_e_dep
    vinf_arr_mars_vec_in = s1.v2 - v_m_fly
    vinf_dep = float(np.linalg.norm(vinf_dep_vec))
    vinf_arr_mars = float(np.linalg.norm(vinf_arr_mars_vec_in))
    c3_dep = vinf_dep**2

    # --- Leg 2: Mars -> Earth. --- #
    sols2 = lambert(r_m_fly, r_e_ret, tof2, prograde=True, max_revs=0)
    s2 = sols2[0]
    vinf_dep_mars_vec_out = s2.v1 - v_m_fly
    vinf_arr_earth_vec = s2.v2 - v_e_ret
    vinf_dep_mars = float(np.linalg.norm(vinf_dep_mars_vec_out))
    vinf_arr_earth = float(np.linalg.norm(vinf_arr_earth_vec))
    c3_arr_earth = vinf_arr_earth**2

    # --- Mars flyby ballistic-continuity check. --- #
    mars = PLANETS["M"]
    rp_min = mars.radius_eq_km + MARS_PERIAPSIS_ALT_KM_PUB
    speed_mismatch = abs(vinf_arr_mars - vinf_dep_mars)
    bend_required = bend_angle(vinf_arr_mars_vec_in, vinf_dep_mars_vec_out)
    v_mean = 0.5 * (vinf_arr_mars + vinf_dep_mars)
    bend_cone = max_bend(mars.mu_km3_s2, rp_min, v_mean)

    def line(label: str, ours: float, pub: float, unit: str) -> str:
        d = ours - pub
        pct = 100.0 * d / pub if pub else float("nan")
        return f"  {label:<28} ours={ours:10.4f}  pub={pub:10.4f}  d={d:+8.4f} {unit} ({pct:+.2f}%)"

    print("=" * 78)
    print("Tito et al. (IEEE Aerospace 2013) 2018 Mars free-return reproduction")
    print("OUR ephemeris: DE440 (astropy).  PUBLISHED: DE421 (Astrogator, Tables III/IV).")
    print("=" * 78)
    print("\nEpochs (published UTCG, fixed as targets):")
    print(f"  depart Earth : {DEPART_EARTH_UTC}  -> t_sec={t_dep:.3f}")
    print(f"  flyby  Mars  : {FLYBY_MARS_UTC}  -> t_sec={t_fly:.3f}")
    print(f"  return Earth : {RETURN_EARTH_UTC}  -> t_sec={t_ret:.3f}")
    print("\nFlight times (days):")
    print(line("tof1 Earth->Mars", tof1_days, TOF1_DAYS_PUB, "d"))
    print(line("tof2 Mars->Earth", tof2_days, TOF2_DAYS_PUB, "d"))
    print(line("total", tof1_days + tof2_days, TOTAL_DAYS_PUB, "d"))

    print("\nLeg 1 (Earth depart, Mars arrive):")
    print(line("vinf depart Earth (km/s)", vinf_dep, VINF_DEP_EARTH_PUB, "km/s"))
    print(line("C3 depart Earth (km2/s2)", c3_dep, C3_DEP_EARTH_PUB, "km2/s2"))
    print(line("vinf arrive Mars (km/s)", vinf_arr_mars, VINF_ARR_MARS_PUB, "km/s"))

    print("\nLeg 2 (Mars depart, Earth arrive):")
    print(line("vinf depart Mars (km/s)", vinf_dep_mars, VINF_DEP_MARS_PUB, "km/s"))
    print(line("vinf arrive Earth (km/s)", vinf_arr_earth, VINF_ARR_EARTH_PUB, "km/s"))
    print(line("C3 arrive Earth (km2/s2)", c3_arr_earth, C3_ARR_EARTH_PUB, "km2/s2"))

    print("\nMars flyby continuity (ballistic free-return constraint):")
    print(f"  vinf speed mismatch (in vs out)   : {speed_mismatch:.4f} km/s")
    print(f"  required bend angle               : {np.degrees(bend_required):.3f} deg")
    print(
        f"  achievable cone @ {MARS_PERIAPSIS_ALT_KM_PUB:.0f} km alt   : "
        f"{np.degrees(bend_cone):.3f} deg"
    )
    print(f"  ballistic-feasible (bend<=cone)?  : {bend_required <= bend_cone}")


if __name__ == "__main__":
    main()
