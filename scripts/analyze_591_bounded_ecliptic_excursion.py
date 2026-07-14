"""#591: rank already-validated BOUNDED 3D Sun-Earth structures by out-of-ecliptic
excursion (max |z|).

This is Fable's salvageable reframe of #589 (design review NO-GO'd unconstrained
z^2-maximization search, see docs/notes/2026-07-14-589-fable-design-review.md):
"maximize out-of-ecliptic excursion SUBJECT TO remaining bounded/returning" --
computed here as a derived-metric ranking over orbits ALREADY known to be
bounded, not a new search. Pool: the 5 known Gurfil-Kasdin (2002) 3D families
with nonzero z (J, K, L, M, N -- from TABLE34 in run_581_gurfil_reproduction.py)
plus #588's cluster 43 (confirmed bounded to 1000+ years by #590's long-horizon
check). ``characterize()`` in run_581 computes rmin/rmax (radial distance from
Earth) but not z alone, so this adds a small max|z| extraction over the same
5-year horizon convention the published families were themselves validated at.

No new GA search, no new fitness function, no novelty claim -- see #591 in
data/OUTSTANDING.md.

Usage:
    uv run python scripts/analyze_591_bounded_ecliptic_excursion.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core.er3bp import er3bp_eom
from cyclerfinder.core.er3bp_geocentric import (
    A_AU_KM_GURFIL_KASDIN,
    E_SUN_EARTH_GURFIL_KASDIN,
    MU_SUN_EARTH_GURFIL_KASDIN,
    geocentric_to_barycentric,
    table_interleaved_to_state,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEDUP_FILE = (
    REPO_ROOT / "data" / "found" / "583_widened_search" / "unmatched_20_for_adjudication.json"
)
OUT_FILE = REPO_ROOT / "data" / "found" / "591_bounded_ecliptic_excursion.json"

# Known Gurfil-Kasdin 3D families with nonzero z (TABLE34 in run_581_gurfil_reproduction.py).
# (ic_interleaved, theta0, type). Type matters: DEO ("Delayed Escape Orbit") is
# Gurfil-Kasdin's own label for NOT-bounded/escaping orbits -- L and N are DEO
# by definition, so they are reported but excluded from the "bounded" ranking
# below (confirmed: both actually escape within the 5yr window at t=19.8/28.7
# rad, i.e. ~3.15yr/4.57yr -- consistent with their own type label, not a bug).
KNOWN_3D_FAMILIES: dict[str, tuple[list[float], float, str]] = {
    "J": (
        [
            0.03348096835548,
            -0.00046191606162,
            0.00774766945226,
            -0.06652559750991,
            0.03675673393090,
            -0.00902011692574,
        ],
        0.0,
        "3D DRO",
    ),
    "K": (
        [
            0.00583817607709,
            0.00021213377191,
            0.00005468845617,
            -0.02914895441397,
            -0.00090918109209,
            0.00123524681998,
        ],
        0.0,
        "3D DRO",
    ),
    "L": (
        [
            0.00385929404386,
            0.00845656433011,
            0.00385933090466,
            0.00001057634437,
            0.00949726128195,
            0.00000748585020,
        ],
        0.0,
        "3D DEO",
    ),
    "M": (
        [0.00386035674766, 0, 0.00385944561146, 0, 0.00386063960861, 0],
        0.0,
        "3D ERO",
    ),
    "N": (
        [0.00766121767444, 0, 0, 0, 0.00668449197861, 0],
        0.0,
        "3D DEO",
    ),
}

FIVE_YEAR_SPAN = 10.0 * math.pi  # matches run_581's characterize() "5yr" convention


def max_abs_z_km(ic_interleaved: list[float], theta0: float) -> dict[str, float | bool]:
    """Propagate for 5 years (paper's own 'practically stable' horizon) and
    return the maximum |z| reached (km and AU), plus whether escape/collision
    terminated the run early (in which case the reported max is a lower bound,
    not the true 5yr maximum)."""
    state0 = table_interleaved_to_state(np.array(ic_interleaved, dtype=float))
    bary0 = geocentric_to_barycentric(state0, MU_SUN_EARTH_GURFIL_KASDIN)
    offset = 1.0 - MU_SUN_EARTH_GURFIL_KASDIN

    def collision(_f: float, s: np.ndarray, *_a: float) -> float:
        dx = s[0] - offset
        return float(dx * dx + s[1] ** 2 + s[2] ** 2 - (6378.0 / A_AU_KM_GURFIL_KASDIN) ** 2)

    collision.terminal = True  # type: ignore[attr-defined]

    def escape(_f: float, s: np.ndarray, *_a: float) -> float:
        dx = s[0] - offset
        return float(dx * dx + s[1] ** 2 + s[2] ** 2 - 0.25)  # 0.5 AU

    escape.terminal = True  # type: ignore[attr-defined]

    sol = solve_ivp(
        er3bp_eom,
        (theta0, theta0 + FIVE_YEAR_SPAN),
        bary0,
        args=(MU_SUN_EARTH_GURFIL_KASDIN, E_SUN_EARTH_GURFIL_KASDIN),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        max_step=0.01,
        events=(collision, escape),
    )
    z_au = sol.y[2]
    max_z_au = float(np.abs(z_au).max())
    return {
        "max_abs_z_au": max_z_au,
        "max_abs_z_km": max_z_au * A_AU_KM_GURFIL_KASDIN,
        "terminated_early": bool(sol.status == 1),
        "theta_end": float(sol.t[-1] - theta0),
    }


def main() -> None:
    results: list[dict[str, object]] = []

    for name, (ic, theta0, type_str) in KNOWN_3D_FAMILIES.items():
        r = max_abs_z_km(ic, theta0)
        results.append({"source": "known-family", "id": name, "type": type_str, **r})

    if DEDUP_FILE.exists():
        data = json.loads(DEDUP_FILE.read_text())
        c43 = next((c for c in data if c["cluster_id"] == 43), None)
        if c43 is not None:
            ic = c43["characterization"]["ic_interleaved"]
            theta0 = c43["characterization"]["theta0"]
            r = max_abs_z_km(ic, theta0)
            results.append({"source": "588-cluster", "id": "43", "type": "3D DRO", **r})

    # DEO ("Delayed Escape Orbit") is Gurfil-Kasdin's own label for a NOT-bounded
    # type -- excluded from the "bounded structures" ranking #591 is actually
    # about, reported separately instead of silently mixed in.
    bounded = [r for r in results if r["type"] != "3D DEO"]
    deo = [r for r in results if r["type"] == "3D DEO"]

    def _sort_key(r: dict[str, object]) -> float:
        return float(r["max_abs_z_km"])  # type: ignore[arg-type]

    bounded.sort(key=_sort_key, reverse=True)
    deo.sort(key=_sort_key, reverse=True)

    header = (
        f"{'id':>6}  {'type':8s}  {'source':12s}  "
        f"{'max|z| (km)':>14}  {'max|z| (AU)':>12}  {'terminated_early':>16}"
    )

    def _print_table(rows: list[dict[str, object]]) -> None:
        print(header)
        for r in rows:
            print(
                f"{r['id']:>6}  {r['type']:8s}  {r['source']:12s}  {r['max_abs_z_km']:>14,.0f}  "
                f"{r['max_abs_z_au']:>12.5f}  {r['terminated_early']!s:>16}"
            )

    print("=== BOUNDED 3D structures, ranked by out-of-ecliptic excursion ===")
    _print_table(bounded)
    print("\n=== 3D DEO (NOT bounded by type definition), reported separately ===")
    _print_table(deo)

    results = bounded + deo

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {OUT_FILE}")


if __name__ == "__main__":
    main()
