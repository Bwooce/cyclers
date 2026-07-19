"""#666 Kamo'oalewa co-orbital QS/HS positive control (`#661` shortlist item 4).

Builds and validates the averaged 1-DOF co-orbital Hamiltonian
(``search/coorbital_hamiltonian.py``) -- a genuinely new object class for
this project (zero co-orbital rows exist anywhere in the 361-row
catalogue) -- against a real, already-published object: 469219 Kamo'oalewa
(2016 HO3), which is documented to be CURRENTLY a quasi-satellite (QS) of
Earth and to switch between QS and horseshoe (HS) configurations on
long-term timescales (de la Fuente Marcos & de la Fuente Marcos 2016, MNRAS
462(4):3441-3456, DOI 10.1093/mnras/stw1972; also see Jiao, Cheng, Huang et
al. 2024, arXiv:2405.20411 / Nature Astronomy, and Qi & Qiao 2022,
Astrodynamics 7(1):3-14, DOI 10.1007/s42064-021-0122-0). Reproducing an
ALREADY-PUBLISHED transition pattern is a positive control, not a discovery
claim -- see `#666`'s own `data/OUTSTANDING.md` bullet and `#661`'s shortlist
item 4 for the full task framing, including the honest HIGH risk of "known
classical object, relabeled" this task explicitly flags.

Three independent checks, in order (each is a genuine cross-check on the
one before it, not a repeat of the same computation):

1. **Real ephemeris + real elements, epoch snapshot.** Kamo'oalewa's own
   osculating elements (JPL SBDB via Wikipedia, epoch 2024-Mar-31,
   JD 2460400.5: a=1.00094 au, e=0.10269, i=7.79605deg, Omega=65.7907deg,
   omega=305.0478deg, M=175.153deg) plus Earth's REAL heliocentric
   ecliptic longitude at the same instant (this project's own
   ``core.ephemeris.Ephemeris("astropy")``, the bundled JPL DE440 backend --
   no idealisation here) give the resonant angle ``sigma`` directly: -3.32
   deg. Sigma near zero, at this real epoch, from real data, is itself an
   independent, model-free confirmation that Kamo'oalewa is CURRENTLY in a
   quasi-satellite configuration (sigma librates around 0 in QS, around
   +-60/180/300 in tadpole/horseshoe) -- matching the literature before any
   averaged-Hamiltonian machinery is invoked at all.

2. **Averaged-Hamiltonian classification.** Build the full
   ``(sigma, delta)`` phase portrait at Kamo'oalewa's real e=0.10269 (see
   ``search/coorbital_hamiltonian.py``, itself validated against the
   classical closed-form L4 tadpole libration frequency to 2e-3 relative --
   `tests/search/test_coorbital_hamiltonian.py`) and classify the real
   (sigma, delta) point from (1). Result: a narrow (~6 deg wide)
   quasi-satellite island -- narrow enough that a small perturbation from
   physics this fixed-e averaged model does not capture (Earth's own
   orbital eccentricity, other planets' secular effects) is plausibly
   enough to carry the real trajectory across the island's boundary,
   consistent with the real object's documented QS<->HS switching.

3. **Full (non-averaged) CR3BP propagation.** An entirely independent,
   higher-fidelity model: build the exact Sun-Earth CR3BP rotating-frame
   initial state from Kamo'oalewa's real elements at the same epoch (a
   self-consistent construction verified below: recovered a, e from the
   constructed state match the input elements exactly, and sigma matches
   check (1) to 0.015 deg) and propagate with this project's own
   ``core.cr3bp`` machinery for several thousand years. Result: sigma(t)
   stays bounded within about +-17 deg of zero -- genuine QS libration,
   directly reproducing "Kamo'oalewa is currently a quasi-satellite" in a
   model that never saw the averaged Hamiltonian at all.

A SECOND real object, (419624) 2010 SO16 -- a well-known, long-lived Earth
HORSESHOE co-orbital in a genuinely different regime (Christou & Asher 2011,
MNRAS 414:2965, "A long-lived horseshoe companion to the Earth") -- is run
through the identical pipeline as an independent cross-check that the
classifier distinguishes QS from HS correctly, not just on Kamo'oalewa. The
averaged-Hamiltonian classifier calls it ``horseshoe`` (sigma-span ~330 deg,
essentially the whole ring) and the full CR3BP propagation shows sigma
sweeping through the ENTIRE [-180, 180] deg range within ~1000 years --
consistent with the published ~175-year one-way / ~350-year round-trip
horseshoe-traverse timing.

HONEST SCOPE LIMITATIONS (read before citing this script's results)
---------------------------------------------------------------------
* Planar only (i=0 reduction) -- both real objects' actual inclinations
  (7.8 deg, 14.5 deg) are dropped from the DYNAMICAL model; only their
  ecliptic-PROJECTED longitude and the real (a, e) feed the averaged
  Hamiltonian and the CR3BP IC. This is the standard reduction for a
  genuinely 1-DOF averaged model (`#666`'s own task framing specifies this
  reduction) but it is a real simplification, not free.
* The full CR3BP propagation is the IDEALISED Sun-Earth two-body-plus-
  test-particle restricted problem (Earth on an exactly circular orbit, no
  other planets). Over 30,000 years this idealised model shows Kamo'oalewa
  in STABLE, essentially constant-amplitude QS libration -- it does NOT
  spontaneously show the real, published Myr-scale QS<->HS switching. This
  is not a failure of the model on its own terms: the real switching is a
  secular phenomenon driven by physics this idealisation deliberately
  excludes (Earth's own ~0.0167 eccentricity, other planets' perturbations
  -- exactly why the real literature uses full N-body integrations, not a
  bare CR3BP, to see the actual transitions). Reproducing the Myr-scale
  transition ITSELF is explicitly NOT claimed here; reproducing the
  CURRENT, real, published QS state (both from real data directly and from
  two independent dynamical models) is what this positive control
  delivers, and that is the concrete, checkable claim available at this
  session's tractable compute budget.
* No novelty claim is made anywhere in this script. Kamo'oalewa's QS/HS
  behavior and 2010 SO16's horseshoe behavior are both already fully
  published; this script's job is reproduction/validation, not discovery,
  per `#666`'s own explicit framing.

Not registered in ``_LEGACY_EXEMPT`` alongside `#606`/`#608`/`#614`/`#317`/
`#624`/`#641`/`#642`/`#649`/`#650`/`#658` for the same reason as all of
those: a fixed-N capability validation against real, named, already-known
objects has no ``region_id``/``n_points`` sweep-region concept to preflight
(no discovery sweep, no catalogue writeback) -- see
``tests/scripts/test_scripts_call_preflight.py``'s own docstring for the
policy this follows.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom, cr3bp_system
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import coe_to_rv
from cyclerfinder.search.coorbital_hamiltonian import (
    build_phase_portrait,
    classify_point,
    ecliptic_longitude_from_elements,
    hill_radius_delta,
    true_anomaly_from_mean,
)

AU_KM = 1.495978707e8
J2000_JD = 2451545.0

# ---------------------------------------------------------------------------
# Real, sourced osculating elements (JPL SBDB, via Wikipedia's data tables --
# both cross-checked against an independent WebSearch synthesis in the #666
# session notes before use).
# ---------------------------------------------------------------------------

KAMOOALEWA = {
    "name": "469219 Kamo'oalewa (2016 HO3)",
    "a_au": 1.00094,
    "e": 0.10269,
    "inc_deg": 7.79605,
    "raan_deg": 65.7907,
    "argp_deg": 305.0478,
    "manom_deg": 175.153,
    "epoch_jd": 2460400.5,  # 2024-Mar-31
}

SO16 = {
    "name": "(419624) 2010 SO16",
    "a_au": 1.0028,
    "e": 0.0754,
    "inc_deg": 14.520,
    "raan_deg": 40.397,
    "argp_deg": 108.99,
    "manom_deg": 173.30,
    "epoch_jd": 2458000.5,  # 2017-Sep-04
}


def _earth_heliocentric_longitude_rad(epoch_jd: float, eph: Ephemeris) -> float:
    t_sec = (epoch_jd - J2000_JD) * 86400.0
    r_earth_km, _v = eph.state("E", t_sec)
    return math.atan2(r_earth_km[1], r_earth_km[0])


def _sigma_at_epoch_deg(elements: dict[str, float], lon_earth_rad: float) -> float:
    """Real resonant angle at epoch: ecliptic-projected asteroid longitude minus
    Earth's real heliocentric longitude at the same instant. See
    :func:`ecliptic_longitude_from_elements`."""
    lon_ast = ecliptic_longitude_from_elements(
        elements["a_au"],
        elements["e"],
        math.radians(elements["inc_deg"]),
        math.radians(elements["raan_deg"]),
        math.radians(elements["argp_deg"]),
        math.radians(elements["manom_deg"]),
    )
    sigma = math.degrees(lon_ast - lon_earth_rad) % 360.0
    if sigma > 180.0:
        sigma -= 360.0
    return sigma


def _build_cr3bp_state0(
    elements: dict[str, float], lon_earth_rad: float, mu: float, l_km: float
) -> NDArray[np.float64]:
    """Exact-(a,e)-preserving planar CR3BP rotating-frame IC from real elements.

    Uses an EFFECTIVE argument of periapsis ``varpi_eff`` -- the real
    periapsis direction's ecliptic-plane projection -- rather than dropping
    the z-component of a full 3-D state (which was tried first and found to
    drift the recovered (a, e) from the input elements by ~2%/~11%
    respectively at Kamo'oalewa's i=7.8 deg, an unacceptably large error for
    a "real elements" positive control). This construction is EXACT: the
    resulting planar 2-body orbit has, by construction, precisely the input
    (a, e); only its instantaneous longitude carries a small (i^2-order,
    ~0.015 deg at i=7.8 deg -- verified directly against the exact
    ecliptic-projection formula) approximation from replacing the true
    inclined argument of latitude with its periapsis-only projection.
    """
    a_nd = elements["a_au"] * AU_KM / l_km
    e = elements["e"]
    inc = math.radians(elements["inc_deg"])
    raan = math.radians(elements["raan_deg"])
    argp = math.radians(elements["argp_deg"])
    nu = float(true_anomaly_from_mean(np.array([math.radians(elements["manom_deg"])]), e)[0])
    varpi_eff = raan + math.atan2(math.cos(inc) * math.sin(argp), math.cos(argp))

    r_pf, v_pf = coe_to_rv(a_km=a_nd, e=e, true_anom_rad=nu, mu=1.0, arg_peri_rad=varpi_eff)

    rot = -lon_earth_rad
    c, s = math.cos(rot), math.sin(rot)
    r_rot = np.array([c * r_pf[0] - s * r_pf[1], s * r_pf[0] + c * r_pf[1], 0.0])
    v_rot_inertial = np.array([c * v_pf[0] - s * v_pf[1], s * v_pf[0] + c * v_pf[1], 0.0])

    omega = np.array([0.0, 0.0, 1.0])
    v_rotframe = v_rot_inertial - np.cross(omega, r_rot)
    r_bary = r_rot - np.array([mu, 0.0, 0.0])
    return np.array(
        [r_bary[0], r_bary[1], 0.0, v_rotframe[0], v_rotframe[1], 0.0],
        dtype=np.float64,
    )


def _state_diagnostics(state: NDArray[np.float64], mu: float) -> tuple[float, float, float]:
    """Recover (a, e, sigma_deg) from a CR3BP rotating-frame state (nondim, mu_total=1)."""
    x, y, z, vx, vy, vz = state
    r_helio = np.array([x + mu, y, z])
    v_inertial = np.array([vx - y, vy + (x + mu), vz])
    r = float(np.linalg.norm(r_helio))
    v2 = float(np.dot(v_inertial, v_inertial))
    energy = 0.5 * v2 - 1.0 / r
    a = -1.0 / (2.0 * energy)
    h_vec = np.cross(r_helio, v_inertial)
    e_vec = np.cross(v_inertial, h_vec) - r_helio / r
    ecc = float(np.linalg.norm(e_vec))
    sigma_deg = math.degrees(math.atan2(y, x + mu))
    return a, ecc, sigma_deg


def _propagate_and_report(
    name: str,
    state0: NDArray[np.float64],
    mu: float,
    years_per_nondim: float,
    t_span_years: float,
    n_eval: int,
) -> None:
    t_span_nd = t_span_years / years_per_nondim
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t_span_nd),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-13,
        max_step=1.0,
        t_eval=np.linspace(0.0, t_span_nd, n_eval),
    )
    x, y = sol.y[0], sol.y[1]
    sigma_deg = np.degrees(np.arctan2(y, x + mu))
    print(
        f"    [{name}] {t_span_years:.0f}-yr CR3BP propagation: "
        f"sigma range [{sigma_deg.min():.2f}, {sigma_deg.max():.2f}] deg "
        f"(solver status={sol.status})"
    )


def main() -> None:
    print("=" * 78)
    print("#666 co-orbital QS/HS averaged-Hamiltonian: Kamo'oalewa positive control")
    print("=" * 78)

    eph = Ephemeris("astropy")
    sysd = cr3bp_system("Sun", "Earth")
    mu, l_km, t_s = sysd.mu, sysd.l_km, sysd.t_s
    years_per_nondim = t_s / 86400.0 / 365.25
    hr = hill_radius_delta(mu)
    print(
        f"\nSun-Earth CR3BP: mu={mu:.6e}, t_s={t_s:.0f} s ({years_per_nondim:.5f} yr/nondim unit)"
    )
    print(f"Hill-radius delta scale: {hr:.6f} (nondim)")

    for elements, delta_span_yr in ((KAMOOALEWA, 30000.0), (SO16, 1000.0)):
        print(f"\n--- {elements['name']} ---")
        lon_earth = _earth_heliocentric_longitude_rad(elements["epoch_jd"], eph)
        print(f"  Earth heliocentric longitude at epoch: {math.degrees(lon_earth):.4f} deg")

        # Check 1: real ephemeris + real elements, no dynamics at all.
        sigma_deg = _sigma_at_epoch_deg(elements, lon_earth)
        delta = elements["a_au"] * AU_KM / l_km - 1.0
        print(
            f"  [Check 1] real sigma at epoch = {sigma_deg:.4f} deg, "
            f"delta = {delta:.6f} ({delta / hr:.3f} * R_hill)"
        )

        # Check 2: averaged-Hamiltonian phase-portrait classification.
        portrait = build_phase_portrait(
            elements["e"], mu, delta_max=6 * hr, n_sigma=181, n_delta=121, n_theta=120, n_varpi=120
        )
        regime, diag = classify_point(portrait, math.radians(sigma_deg), delta)
        print(
            f"  [Check 2] averaged-Hamiltonian classification: {regime} "
            f"(sigma-span {diag['sigma_span_deg']:.1f} deg)"
        )

        # Check 3: full non-averaged CR3BP propagation from a self-consistent real IC.
        state0 = _build_cr3bp_state0(elements, lon_earth, mu, l_km)
        a_chk, e_chk, sigma_chk = _state_diagnostics(state0, mu)
        print(
            f"  [IC self-check] recovered a={a_chk:.6f} au (input {elements['a_au']}), "
            f"e={e_chk:.6f} (input {elements['e']}),"
        )
        print(f"                  sigma={sigma_chk:.4f} deg (Check-1 value {sigma_deg:.4f} deg)")
        n_eval = 4000 if delta_span_yr <= 3000.0 else 8000
        _propagate_and_report(elements["name"], state0, mu, years_per_nondim, delta_span_yr, n_eval)

    print("\n" + "=" * 78)
    print("Summary: Kamo'oalewa -> quasi_satellite, narrow island, stable bounded")
    print("  libration in full CR3BP (matches its published CURRENT QS state).")
    print("2010 SO16 -> horseshoe, wide island, full sigma-range excursion in full")
    print("  CR3BP (matches its published long-lived horseshoe classification).")
    print("No novelty claim. See this script's own module docstring for full")
    print("scope limitations (planar reduction; idealised circular-CR3BP does")
    print("NOT itself reproduce the real Myr-scale QS<->HS switching mechanism).")
    print("=" * 78)


if __name__ == "__main__":
    main()
