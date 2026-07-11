"""#571 -- independent recomputation of the Fable plan-review's analytic-empty numbers.

The 2026-07-12 Fable plan review of #571 (``data/OUTSTANDING.md``) found that
Titan-Mimas / Titan-Enceladus / Titan-Tethys / Titan-Dione are analytically
empty under this project's own two-sided #324 physical max-bend gate: each
small moon's own minimum-achievable ``V_inf`` (over ANY conic reaching
Titan's orbital radius) exceeds the ``V_inf`` at which that moon's own 5 deg
bend floor (``search/physical_sanity.py::DEFAULT_MIN_USEFUL_BEND_DEG``) is
still clearable -- so no basin can ever exist at any grid point, at any
sweep resolution.

This script independently RECOMPUTES those two numbers per moon (rather than
copying Fable's summary figures), using:

* ``core/flyby.py::max_bend`` for the ballistic bend-angle formula (the same
  function ``search/physical_sanity.py::flyby_is_useful`` calls at gate
  time), inverted by bisection to find the ``V_inf`` at which max_bend
  crosses the 5 deg floor (the "ceiling" -- ABOVE this V_inf, even a
  grazing periapsis pass cannot bend 5 deg; feasibility requires
  ``min-achievable V_inf <= ceiling``);
* ``core/satellites.py``'s sourced GM/radius/sma/safe_alt values (JPL SSD
  gm_de440 + SAT441, cited in-file) for each moon and Titan;
* a Hohmann-transfer vis-viva calculation for the "min-achievable V_inf at
  the small moon, for ANY conic reaching Titan's radius" bound -- the
  tangential minimum-energy ellipse connecting the two circular orbit radii
  minimizes ``|v_transfer - v_circular|`` at the departure radius among all
  conics connecting the two radii (any non-tangential/non-Hohmann transfer
  through the same two radii has a larger velocity-vector mismatch with the
  local circular velocity at departure), so it is the tightest (most
  favorable-to-feasibility) analytic floor -- if even this best-case V_inf
  cannot clear the moon's own gate, no Lambert-arc grid point in the actual
  #558-style sweep (which only reaches non-tangential-Hohmann, generally
  WORSE, geometries) can either.

Run as::

    uv run python scripts/verify_571_gate_analytics.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.flyby import max_bend  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402

MU_SATURN = PRIMARIES["Saturn"]
TITAN = SATELLITES["Titan"]

MOONS_TO_CHECK: tuple[str, ...] = (
    "Mimas",
    "Enceladus",
    "Tethys",
    "Dione",
    "Rhea",
    "Iapetus",
    "Hyperion",
)


def bend_ceiling_vinf(body: str, target_bend_deg: float = 5.0) -> float:
    """Solve ``max_bend(mu, rp, vinf) == target_bend_deg`` for ``vinf`` (bisection).

    ``max_bend`` is monotonically decreasing in ``vinf`` (faster flybys bend
    less), so bisection on a wide bracket is exact to float precision.
    """
    sat = SATELLITES[body]
    mu = sat.mu_km3_s2
    rp = sat.radius_eq_km + sat.safe_alt_km
    target_rad = math.radians(target_bend_deg)

    lo, hi = 1e-6, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if max_bend(mu, rp, mid) > target_rad:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def hohmann_vinf_at_r1(r1_km: float, r2_km: float, mu_primary: float) -> float:
    """Minimum-achievable ``V_inf`` at radius ``r1`` for ANY conic connecting
    circular orbits at ``r1`` and ``r2`` about ``mu_primary`` -- realized by
    the tangential (Hohmann) minimum-energy transfer ellipse.
    """
    a_t = 0.5 * (r1_km + r2_km)
    v_circ = math.sqrt(mu_primary / r1_km)
    v_transfer = math.sqrt(mu_primary * (2.0 / r1_km - 1.0 / a_t))
    return abs(v_transfer - v_circ)


def bend_at_vinf(body: str, vinf: float) -> float:
    sat = SATELLITES[body]
    return math.degrees(max_bend(sat.mu_km3_s2, sat.radius_eq_km + sat.safe_alt_km, vinf))


def main() -> int:
    print(
        f"Titan: GM={TITAN.mu_km3_s2} km^3/s^2, R={TITAN.radius_eq_km} km, "
        f"sma={TITAN.sma_km} km, safe_alt={TITAN.safe_alt_km} km"
    )
    print()
    header = (
        f"{'moon':10s} {'GM':>10s} {'sma_km':>11s} {'safe_alt':>9s}  "
        f"{'5deg-ceiling-vinf':>18s} {'min-achievable-vinf':>20s} "
        f"{'bend@min':>10s} {'feasible':>9s}"
    )
    print(header)
    for moon in MOONS_TO_CHECK:
        sat = SATELLITES[moon]
        ceiling = bend_ceiling_vinf(moon, 5.0)
        min_vinf = hohmann_vinf_at_r1(sat.sma_km, TITAN.sma_km, MU_SATURN)
        bend_at_min = bend_at_vinf(moon, min_vinf)
        feasible = min_vinf <= ceiling
        print(
            f"{moon:10s} {sat.mu_km3_s2:10.4f} {sat.sma_km:11.1f} {sat.safe_alt_km:9.1f}  "
            f"{ceiling:15.4f} km/s {min_vinf:17.4f} km/s {bend_at_min:9.4f} deg {feasible!s:>9s}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
