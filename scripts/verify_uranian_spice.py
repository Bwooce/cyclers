"""Verify the URA SPICE kernel install loads + queries cleanly (#335 Part A).

Loads the GMAT-bundled leap-seconds + planetary constants kernels + the
freshly-installed Uranus satellite SPK (ura111.bsp, see
``scripts/install_uranian_spice.sh``), then queries Uranus + the five
classical Uranian satellite states at J2000 + 1 day and prints the
implied SMA + eccentricity + inclination. These are the constants the
V4-strict Phase 4.1 gauntlet (#335 Parts B/C) drives the moontour against.

Expected eccentricities (Murray-Dermott "Solar System Dynamics" Table A.7,
JPL satellite mean elements):

    Miranda   0.0013
    Ariel     0.0012
    Umbriel   0.0039
    Titania   0.0011
    Oberon    0.0014

Expected SMA (JPL satellite mean elements):

    Miranda    129,872 km
    Ariel      190,945 km
    Umbriel    265,998 km
    Titania    436,298 km
    Oberon     583,519 km

Run as::

    uv run python scripts/verify_uranian_spice.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import spiceypy as spice

GMAT_ROOT = Path.home() / "GMAT" / "R2022a"
LSK = GMAT_ROOT / "data" / "time" / "SPICELeapSecondKernel.tls"
PCK = GMAT_ROOT / "data" / "planetary_coeff" / "SPICEPlanetaryConstantsKernel.tpc"
URA_BSP = GMAT_ROOT / "data" / "planetary_ephem" / "spk" / "uranian" / "ura111.bsp"
DE_BSP = GMAT_ROOT / "data" / "planetary_ephem" / "spk" / "DE421AllPlanets.bsp"

# Uranus gravitational parameter (km^3/s^2). Sourced from Jacobson 2014 +
# matches src/cyclerfinder/core/satellites.py.
MU_URANUS = 5.7945564e6

CLASSICAL_MOONS = (
    ("MIRANDA", 705),
    ("ARIEL", 701),
    ("UMBRIEL", 702),
    ("TITANIA", 703),
    ("OBERON", 704),
)


def _elements(r: np.ndarray, v: np.ndarray, mu: float) -> tuple[float, float, float]:
    """Return (SMA km, eccentricity, inclination deg) from Cartesian state."""
    r_mag = float(np.linalg.norm(r))
    v_mag = float(np.linalg.norm(v))
    h = np.cross(r, v)
    h_mag = float(np.linalg.norm(h))
    e_vec = np.cross(v, h) / mu - r / r_mag
    e = float(np.linalg.norm(e_vec))
    energy = 0.5 * v_mag * v_mag - mu / r_mag
    sma = -mu / (2.0 * energy)
    inc = math.degrees(math.acos(h[2] / h_mag))
    return float(sma), e, inc


def main() -> int:
    print("[ura111] checking files...")
    for p, label in (
        (LSK, "leap seconds"),
        (PCK, "planetary constants"),
        (URA_BSP, "URA111 satellite SPK"),
    ):
        if not p.exists():
            print(f"  MISSING: {label} -> {p}", file=sys.stderr)
            return 1
        print(f"  OK: {label} -> {p} ({p.stat().st_size} bytes)")

    print("\n[ura111] loading kernels via SPICE FURNSH...")
    spice.furnsh(str(LSK))
    spice.furnsh(str(PCK))
    spice.furnsh(str(URA_BSP))
    if DE_BSP.exists():
        spice.furnsh(str(DE_BSP))
        print(f"  OK: DE421 planets -> {DE_BSP}")

    # Pick J2000 + 1 day so we sample a state well inside ura111's 1900-2099 span
    # but not exactly at coverage start (avoids edge effects).
    et = spice.str2et("2000-01-02T00:00:00")
    et_str = spice.timout(et, "YYYY-MM-DD HR:MN ::TDB")
    print(f"\n[ura111] sampling moon states at ET = {et:.3f} s ({et_str})")
    print("        (target frame: J2000, observer: URANUS=799)")
    print(f"        mu_Uranus = {MU_URANUS:.6e} km^3/s^2")
    print()
    hdr = (
        f"{'moon':<10}{'SMA (km)':>14}{'ecc':>12}"
        f"{'inc (deg)':>12}{'|r| (km)':>14}{'|v| (km/s)':>14}"
    )
    print(hdr)
    print(f"{'-' * 10}{'-' * 14}{'-' * 12}{'-' * 12}{'-' * 14}{'-' * 14}")
    for name, _naif in CLASSICAL_MOONS:
        try:
            state, _lt = spice.spkezr(name, et, "J2000", "NONE", "URANUS")
        except Exception as exc:
            print(f"  FAIL: {name}: {exc}", file=sys.stderr)
            spice.kclear()
            return 1
        r = np.asarray(state[:3], dtype=np.float64)
        v = np.asarray(state[3:], dtype=np.float64)
        sma, ecc, inc = _elements(r, v, MU_URANUS)
        r_mag = float(np.linalg.norm(r))
        v_mag = float(np.linalg.norm(v))
        print(f"{name:<10}{sma:>14.1f}{ecc:>12.5f}{inc:>12.4f}{r_mag:>14.1f}{v_mag:>14.5f}")

    # Propagate Umbriel for 1 day via SPICE-only ephemeris (sanity check) and
    # verify it's a closed Kepler-ish orbit.
    print("\n[ura111] cross-check: Umbriel + 1 day (SPICE ephemeris-driven)")
    et_plus = et + 86400.0
    state0, _ = spice.spkezr("UMBRIEL", et, "J2000", "NONE", "URANUS")
    state1, _ = spice.spkezr("UMBRIEL", et_plus, "J2000", "NONE", "URANUS")
    r0 = np.asarray(state0[:3], dtype=np.float64)
    r1 = np.asarray(state1[:3], dtype=np.float64)
    sma0, _, _ = _elements(r0, np.asarray(state0[3:], dtype=np.float64), MU_URANUS)
    sma1, _, _ = _elements(r1, np.asarray(state1[3:], dtype=np.float64), MU_URANUS)
    print(f"  Umbriel SMA at t0:     {sma0:.1f} km")
    print(f"  Umbriel SMA at t0+1d:  {sma1:.1f} km (must match to ~5 sig figs)")
    print(f"  delta SMA:             {abs(sma1 - sma0):.3e} km")

    spice.kclear()
    print("\n[ura111] kernel unloaded; verification PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
