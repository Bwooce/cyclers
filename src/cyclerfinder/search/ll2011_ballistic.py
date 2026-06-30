"""Ballistic construction for Lynam-Longuski 2011 Laplace-resonant IEG triple cyclers.

Two members from Lynam & Longuski (2011), "Laplace-resonant triple-cyclers for missions
to Jupiter," Acta Astronautica 69(3-4), pp.158-167, DOI 10.1016/j.actaastro.2011.03.011.

SINGLE-PERIOD IEG (``lynam-longuski-2011-ieg-single-period``)
-------------------------------------------------------------
Sequence E-(perijove)-I-G-(apojove)-E, one Laplace period (7.055 d sourced / 7.004 d
ideal model). This is the EIGE topology (Europa-Io-Ganymede-Europa) from Hernandez et al.
2017 (AAS 17-608) with explicit perijove/apojove nodes described. The #480 EIGE work
already characterised it fully:

* In the ideal circular-coplanar model a feasible ballistic EIGE member exists
  (``eige_ballistic.feasible_ballistic_eige()``), total flyby ΔV ~ 0.
* The L-L 2011 member (patched-conic ephemeris / MALTO) requires total ΔV ≈ 11 m/s
  powered at Europa.  The ballistic version would need a -175 km (sub-surface) Europa
  flyby — stated explicitly by L-L themselves.  This is consistent with our ideal-model
  finding that the ballistic EIGE at higher V∞ forces sub-surface Europa flybys (the
  seam's binding constraint in the 2D model; see ``eggie_ballistic`` Gate-B).

GIPEIPE (``lynam-longuski-2011-gipeipe``)
-----------------------------------------
Sequence G-I-(perijove)-E-I-(perijove)-E, orbital period 3.5 d (1:2 resonance with
Ganymede), full sequence repeats every Laplace period 7.055 d.  This is a NEW topology
not characterised by our prior #480 work (EGGIE / EIGE have 1:1 and 4:5 resonances).

Construction approach (same philosophy as ``eggie_ballistic`` / ``eige_ballistic``):
Five Lambert legs G→I₁→E₁→I₂→E₂→G closing back on Ganymede after one Laplace period
T_Gan ≈ 7.004 d (ideal model).

* Free variables: φ_Io (Io phase at t=0), φ_Eur (Europa phase at t=0), four leg ToFs
  (leg-5 ToF derived so total = T_Gan). Ganymede phase fixed at 0 (gauge).
* Core ballistic residuals: equal-in/out |V∞| at I₁, E₁, I₂, E₂ (4 constraints).
* Seam residual: |V∞_in(G₂)| = |V∞_out(G₁)| — the periodicity closure at Ganymede.
* Degrees of freedom: 5 free vars - 5 constraints = 0 DOF (well-determined).

The I→E legs pass THROUGH perijove — the spacecraft's arc goes inward past Io's orbit,
around Jupiter's closest approach, and back outward to Europa. In Lambert terms, these
arcs use ``prograde=False`` (retrograde direction selects the through-perijove branch
when Io and Europa are in a certain geometric relationship; see ``build_legs`` notes).

Sourced anchors (``feedback_golden_tests_sourced_only``): Lynam-Longuski 2011 Table 1
(as mined in docs/notes/2026-06-30-490-russell-strange-lynam-mining.md):
* GIPEIPE orbital period 3.5 d; sequence period = Laplace period 7.055 d.
* Single-period IEG: total ΔV ≈ 11 m/s (patched-conic/MALTO), Europa flyby -175 km
  (ballistic version sub-surface).
* The L-L construction uses the ideal circular-coplanar Galilean model.

References
----------
* docs/notes/2026-06-30-digest-lynam-longuski-2011-laplace-resonant-triple-cyclers.md
* docs/notes/2026-06-30-490-russell-strange-lynam-mining.md
* docs/notes/2026-06-30-480-eige-ballistic-construction-verdict.md (EIGE topology)
* src/cyclerfinder/search/eige_ballistic.py (single-period IEG = EIGE)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.lambert import lambert
from cyclerfinder.nbody.jovian import flyby_maneuver_dv
from cyclerfinder.search.eggie_ballistic import MU, SECONDS_PER_DAY, moon_state
from cyclerfinder.search.eige_ballistic import (
    ALT_MAX_KM,
    ALT_MIN_KM,
    EIGE_RESONANT_SMA_KM,
    feasible_ballistic_eige,
)
from cyclerfinder.search.resonant_conic import (
    ideal_moon_smas,
    ideal_t_syn,
    resonant_sma,
)

Vec3 = NDArray[np.float64]

# --- Sourced period constants ------------------------------------------------

#: Sourced Laplace (Ganymede) period from L-L 2011 (d).
LL2011_LAPLACE_PERIOD_DAYS: float = 7.055

#: Sourced GIPEIPE spacecraft orbital period (d).
LL2011_GIPEIPE_ORBIT_PERIOD_DAYS: float = 3.5

#: Sourced single-period IEG total ΔV (m/s), patched-conic/MALTO.
LL2011_IEG_DV_MS: float = 11.0

#: Sourced single-period IEG ballistic-mode Europa periapsis altitude (km).
#: Negative = sub-surface (infeasible).
LL2011_IEG_BALLISTIC_EUROPA_ALT_KM: float = -175.0

# --- Ideal-model derived constants ------------------------------------------

#: Ideal Ganymede orbital period (s) — the ideal Laplace period.
T_LAPLACE_S: float = ideal_t_syn()

#: Ideal Laplace period (d).
T_LAPLACE_D: float = T_LAPLACE_S / SECONDS_PER_DAY

#: GIPEIPE resonant spacecraft semi-major axis (km), n_syn=1, n_rev=2.
GIPEIPE_SMA_KM: float = resonant_sma(1, 2, T_LAPLACE_S)

#: GIPEIPE spacecraft orbital period (s).
GIPEIPE_PERIOD_S: float = 2.0 * math.pi * math.sqrt(GIPEIPE_SMA_KM**3 / MU)

#: GIPEIPE spacecraft orbital period (d).
GIPEIPE_PERIOD_D: float = GIPEIPE_PERIOD_S / SECONDS_PER_DAY

#: GIPEIPE eccentricity — apojove at ideal Ganymede orbit.
_SMAS: dict[str, float] = ideal_moon_smas()
GIPEIPE_ECC: float = _SMAS["Ganymede"] / GIPEIPE_SMA_KM - 1.0

#: GIPEIPE flyby sequence (6 nodes, 5 legs, closes back at Ganymede).
GIPEIPE_SEQUENCE: tuple[str, ...] = ("Ganymede", "Io", "Europa", "Io", "Europa", "Ganymede")

#: Single-period IEG sequence (same topology as Hernandez 2017 EIGE).
IEG_SINGLE_SEQUENCE: tuple[str, ...] = ("Europa", "Io", "Ganymede", "Europa")


def _nrm(v: Vec3) -> float:
    return float(np.linalg.norm(v))


# --- Period characterisation (read-only, no construction needed) --------------


def single_period_ieg_summary() -> dict[str, float]:
    """Return ideal-model invariants for the L-L 2011 single-period IEG.

    This cycler is the EIGE topology (Europa-Io-Ganymede-Europa, one Laplace period).
    The period and ballistic feasibility are characterised by the existing #480 EIGE
    work; this function reports the key numbers vs the L-L 2011 sourced values.

    Returns a dict with keys:
    * ``period_ideal_d``: ideal Ganymede period (ideal Laplace period, days)
    * ``period_sourced_d``: sourced L-L period (days)
    * ``period_delta_pct``: fractional error (ideal vs sourced, %)
    * ``resonant_sma_km``: 1:1 resonant spacecraft SMA (km)
    * ``eige_total_dv_ms``: total flyby ΔV for the feasible ideal-model ballistic
      EIGE (m/s) — expected ~0 (ballistic)
    * ``eige_all_feasible``: True if all altitudes in the 25-70,000 km window
    * ``sourced_powered_dv_ms``: L-L sourced powered ΔV (m/s)
    * ``sourced_ballistic_europa_alt_km``: L-L sourced ballistic Europa altitude (km)
    """
    eige = feasible_ballistic_eige()
    period_d = T_LAPLACE_D
    return {
        "period_ideal_d": period_d,
        "period_sourced_d": LL2011_LAPLACE_PERIOD_DAYS,
        "period_delta_pct": (
            100.0 * (period_d - LL2011_LAPLACE_PERIOD_DAYS) / LL2011_LAPLACE_PERIOD_DAYS
        ),
        "resonant_sma_km": EIGE_RESONANT_SMA_KM,
        "eige_total_dv_ms": eige.total_dv_ms,
        "eige_all_feasible": eige.all_feasible,
        "sourced_powered_dv_ms": LL2011_IEG_DV_MS,
        "sourced_ballistic_europa_alt_km": LL2011_IEG_BALLISTIC_EUROPA_ALT_KM,
    }


def gipeipe_period_summary() -> dict[str, float]:
    """Return ideal-model period invariants for the GIPEIPE cycler.

    Returns a dict with keys:
    * ``orbit_period_ideal_d``: ideal 1:2-resonant spacecraft orbital period (d)
    * ``orbit_period_sourced_d``: sourced L-L orbital period (d)
    * ``orbit_period_delta_pct``: fractional error (%)
    * ``sequence_period_ideal_d``: full sequence period (2 orbits = Laplace period, d)
    * ``sequence_period_sourced_d``: sourced Laplace period (d)
    * ``resonant_sma_km``: GIPEIPE resonant S/C SMA (km)
    * ``eccentricity``: conic eccentricity for apojove at Ganymede
    """
    return {
        "orbit_period_ideal_d": GIPEIPE_PERIOD_D,
        "orbit_period_sourced_d": LL2011_GIPEIPE_ORBIT_PERIOD_DAYS,
        "orbit_period_delta_pct": (
            100.0
            * (GIPEIPE_PERIOD_D - LL2011_GIPEIPE_ORBIT_PERIOD_DAYS)
            / LL2011_GIPEIPE_ORBIT_PERIOD_DAYS
        ),
        "sequence_period_ideal_d": T_LAPLACE_D,
        "sequence_period_sourced_d": LL2011_LAPLACE_PERIOD_DAYS,
        "resonant_sma_km": GIPEIPE_SMA_KM,
        "eccentricity": GIPEIPE_ECC,
    }


# --- GIPEIPE Lambert construction -------------------------------------------


def _build_gipeipe_legs(
    phi_io: float,
    phi_eur: float,
    tofs_s: tuple[float, ...] | list[float],
    ie_prograde: bool = False,
) -> tuple[list[Vec3 | None], list[Vec3 | None]] | None:
    """Build 5 Lambert legs for the GIPEIPE G→I→E→I→E→G tour.

    Ganymede phase is fixed at 0 (gauge). ``phi_io`` / ``phi_eur`` are the Io /
    Europa initial phases (rad). ``tofs_s`` contains 5 ToFs (seconds); the 5th
    must satisfy ``sum(tofs_s) ≈ T_LAPLACE_S`` for periodicity (caller enforces).

    The I→E legs (legs 1 and 3) go through perijove — the spacecraft leaves Io
    heading inward, passes perijove, and rises to Europa. Setting ``ie_prograde=False``
    selects the retrograde Lambert direction (the through-perijove branch when Io and
    Europa are on the same side of Jupiter).

    Returns ``(vinf_out, vinf_in)`` lists indexed by node 0..5 (None at unused slots),
    or None if any leg fails (non-positive ToF / no Lambert solution).
    """
    t = [0.0]
    for tof in tofs_s:
        t.append(t[-1] + float(tof))

    phases = {"Ganymede": 0.0, "Io": phi_io, "Europa": phi_eur}
    rpos: list[Vec3] = []
    rvel: list[Vec3] = []
    for k, moon in enumerate(GIPEIPE_SEQUENCE):
        r, v = moon_state(moon, phases[moon], t[k])
        rpos.append(r)
        rvel.append(v)

    vinf_out: list[Vec3 | None] = [None] * 6
    vinf_in: list[Vec3 | None] = [None] * 6

    for k in range(5):
        leg_tof = t[k + 1] - t[k]
        if leg_tof <= 0.0:
            return None
        # Legs 1 and 3 are I→E through-perijove: use ie_prograde direction.
        # Legs 0, 2, 4 are G→I, E→I, E→G: prograde direct arcs.
        is_ie_leg = k in (1, 3)
        prog = (not ie_prograde) if is_ie_leg else True
        try:
            sols = lambert(rpos[k], rpos[k + 1], leg_tof, mu=MU, prograde=prog, max_revs=0)
        except Exception:
            return None
        if not sols:
            return None
        sol = sols[0]  # single-rev solution
        vinf_out[k] = np.asarray(sol.v1, dtype=np.float64) - rvel[k]
        vinf_in[k + 1] = np.asarray(sol.v2, dtype=np.float64) - rvel[k + 1]

    return vinf_out, vinf_in


def _gipeipe_residual(x: NDArray[np.float64], ie_prograde: bool = False) -> NDArray[np.float64]:
    """Ballistic-cycler residuals for GIPEIPE (km/s).

    ``x = [phi_io, phi_eur, tof1, tof2, tof3, tof4]`` (ToFs in seconds).
    ``tof5 = T_LAPLACE_S - tof1 - tof2 - tof3 - tof4`` (periodicity constraint).

    Residuals (5 constraints):
    0: equal-in/out |V∞| at I₁ (Io first encounter, leg 0→1)
    1: equal-in/out |V∞| at E₁ (Europa first encounter, leg 1→2)
    2: equal-in/out |V∞| at I₂ (Io second encounter, leg 2→3)
    3: equal-in/out |V∞| at E₂ (Europa second encounter, leg 3→4)
    4: Ganymede seam — |V∞_in(G₂)| = |V∞_out(G₁)| (periodicity closure)
    """
    phi_io, phi_eur = x[0], x[1]
    tofs4 = list(x[2:6])
    tof5 = T_LAPLACE_S - sum(tofs4)
    if tof5 <= 0.0:
        return np.full(5, 20.0)
    tofs = [*tofs4, tof5]
    built = _build_gipeipe_legs(phi_io, phi_eur, tofs, ie_prograde=ie_prograde)
    if built is None:
        return np.full(5, 20.0)
    vo, vi = built
    return np.array(
        [
            _nrm(vi[1]) - _nrm(vo[1]),  # type: ignore[arg-type]  # I₁ ballistic
            _nrm(vi[2]) - _nrm(vo[2]),  # type: ignore[arg-type]  # E₁ ballistic
            _nrm(vi[3]) - _nrm(vo[3]),  # type: ignore[arg-type]  # I₂ ballistic
            _nrm(vi[4]) - _nrm(vo[4]),  # type: ignore[arg-type]  # E₂ ballistic
            _nrm(vi[5]) - _nrm(vo[0]),  # type: ignore[arg-type]  # Ganymede seam
        ],
        dtype=np.float64,
    )


def _gipeipe_resnorm(x: NDArray[np.float64], ie_prograde: bool = False) -> float:
    return float(np.linalg.norm(_gipeipe_residual(x, ie_prograde=ie_prograde)))


@dataclass(frozen=True)
class BallisticGipeipe:
    """A constructed GIPEIPE tour evaluated against the ballistic constraints."""

    phi_io_rad: float
    phi_eur_rad: float
    tofs_days: tuple[float, ...]
    ie_prograde: bool
    #: Departure |V∞| at G₁ and arrival |V∞| at each encounter (km/s).
    vinf_kms: dict[str, float]
    #: Per-flyby maneuver ΔV (m/s): keys G1, I1, E1, I2, E2.
    flyby_dv_ms: dict[str, float]
    #: Per-flyby periapsis altitude (km).
    flyby_alt_km: dict[str, float]
    total_dv_ms: float
    ballistic_resnorm_kms: float
    seam_defect_kms: float
    all_feasible: bool


def _evaluate_gipeipe(
    phi_io: float, phi_eur: float, tofs_s: list[float], ie_prograde: bool
) -> BallisticGipeipe | None:
    """Evaluate a GIPEIPE tour: V∞, flyby ΔV/altitudes, closure."""
    tofs5 = [*tofs_s[:4], T_LAPLACE_S - sum(tofs_s[:4])]
    if tofs5[-1] <= 0.0:
        return None
    built = _build_gipeipe_legs(phi_io, phi_eur, tofs5, ie_prograde=ie_prograde)
    if built is None:
        return None
    vo, vi = built

    flyby_nodes = {
        "G1": (vi[0], vo[0], "Ganymede"),  # initial Ganymede departure (seam)
        "I1": (vi[1], vo[1], "Io"),
        "E1": (vi[2], vo[2], "Europa"),
        "I2": (vi[3], vo[3], "Io"),
        "E2": (vi[4], vo[4], "Europa"),
        "G2": (vi[5], vo[5], "Ganymede"),  # closing Ganymede arrival (seam)
    }
    dv_ms: dict[str, float] = {}
    alt_km: dict[str, float] = {}
    feas: list[bool] = []

    for key, (vin, vout, moon) in flyby_nodes.items():
        if vin is None or vout is None:
            dv_ms[key] = 999.0
            alt_km[key] = -999.0
            feas.append(False)
            continue
        dv, alt, ok = flyby_maneuver_dv(vin, vout, moon, alt_min=ALT_MIN_KM, alt_max=ALT_MAX_KM)
        dv_ms[key] = dv
        alt_km[key] = alt
        feas.append(ok)

    vinf: dict[str, float] = {}
    for i, (vo_i, vi_i) in enumerate(zip(vo, vi, strict=False)):
        if vo_i is not None:
            vinf[f"dep_{GIPEIPE_SEQUENCE[i]}_{i}"] = _nrm(vo_i)
        if vi_i is not None:
            vinf[f"arr_{GIPEIPE_SEQUENCE[i]}_{i}"] = _nrm(vi_i)

    x = np.array([phi_io, phi_eur, *tofs5[:4]])
    return BallisticGipeipe(
        phi_io_rad=float(phi_io),
        phi_eur_rad=float(phi_eur),
        tofs_days=tuple(t / SECONDS_PER_DAY for t in tofs5),
        ie_prograde=ie_prograde,
        vinf_kms=vinf,
        flyby_dv_ms=dv_ms,
        flyby_alt_km=alt_km,
        total_dv_ms=float(sum(dv_ms.values())),
        ballistic_resnorm_kms=_gipeipe_resnorm(x, ie_prograde=ie_prograde),
        seam_defect_kms=(
            abs(_nrm(vi[5]) - _nrm(vo[0])) if vi[5] is not None and vo[0] is not None else 999.0
        ),
        all_feasible=all(feas),
    )


def refine_gipeipe(
    seed_x: NDArray[np.float64] | list[float],
    *,
    ie_prograde: bool = False,
    max_nfev: int = 4000,
) -> BallisticGipeipe | None:
    """Differential-correct GIPEIPE from ``seed_x`` (trf least-squares).

    ``seed_x = [phi_io, phi_eur, tof1, tof2, tof3, tof4]`` (ToFs in seconds).
    Returns None if the corrector cannot build any legs from the converged point.
    """
    lb = np.array([-np.inf, -np.inf, 60.0, 60.0, 60.0, 60.0])
    ub = np.array([np.inf, np.inf, T_LAPLACE_S, T_LAPLACE_S, T_LAPLACE_S, T_LAPLACE_S])
    x0 = np.clip(np.asarray(seed_x, dtype=np.float64), lb, ub)
    sol = least_squares(
        _gipeipe_residual,
        x0,
        args=(ie_prograde,),
        method="trf",
        bounds=(lb, ub),
        max_nfev=max_nfev,
        xtol=1e-14,
        ftol=1e-14,
    )
    return _evaluate_gipeipe(sol.x[0], sol.x[1], list(sol.x[2:]), ie_prograde=ie_prograde)


# --- Resonant-conic seed for GIPEIPE -----------------------------------------


def _gipeipe_conic_seed() -> NDArray[np.float64]:
    """Generate a rough seed for GIPEIPE from the 1:2 resonant-conic geometry.

    The spacecraft orbit has a=GIPEIPE_SMA_KM, e=GIPEIPE_ECC, apojove at a_Ganymede.
    True anomalies at each moon crossing are derived analytically; the timing for each
    leg estimates the eccentric-anomaly time-of-flight.
    """
    a = GIPEIPE_SMA_KM
    e = GIPEIPE_ECC
    mu = MU

    def nu_at_r(r: float) -> float:
        """True anomaly where the conic crosses radius r."""
        p = a * (1.0 - e * e)
        cos_nu = (p / r - 1.0) / e
        cos_nu = max(-1.0, min(1.0, cos_nu))
        return math.acos(cos_nu)

    def tof_nu1_to_nu2(nu1: float, nu2: float) -> float:
        """Time of flight from true anomaly nu1 to nu2 (prograde, seconds)."""

        def eccentric(nu: float) -> float:
            return 2.0 * math.atan2(
                math.sqrt(1.0 - e) * math.sin(nu / 2.0), math.sqrt(1.0 + e) * math.cos(nu / 2.0)
            )

        def mean(nu: float) -> float:
            e_anom = eccentric(nu)
            return e_anom - e * math.sin(e_anom)

        n = math.sqrt(mu / a**3)
        d_m = mean(nu2) - mean(nu1)
        if d_m < 0.0:
            d_m += 2.0 * math.pi
        return d_m / n

    # Key true anomalies on the spacecraft's reference conic
    # G at apoapsis (outbound arc at nu=pi), I inbound at -nu_Io (i.e. 2pi-nu_Io going ccw)
    nu_io = nu_at_r(_SMAS["Io"])  # ~87 deg (outbound crossing)
    nu_eur = nu_at_r(_SMAS["Europa"])  # ~126 deg (outbound crossing)

    # G at apoapsis (nu=pi). GIPEIPE starts at G going to I (inbound arc G->I_in).
    nu_i_inb = 2.0 * math.pi - nu_io  # ~273 deg inbound Io crossing (ccw from pi)

    # Leg 0: G(pi) -> I_in(2pi-nu_io)
    tof0 = tof_nu1_to_nu2(math.pi, nu_i_inb)

    # Leg 1 (I_in->P->E_out): through-perijove; tof(I_in, 2pi) + tof(0, nu_eur)
    tof1 = tof_nu1_to_nu2(nu_i_inb, 2.0 * math.pi) + tof_nu1_to_nu2(0.0, nu_eur)

    # Leg 2: E_out -> I_in (next half orbit via apoapsis)
    tof2 = tof_nu1_to_nu2(nu_eur, math.pi) + tof_nu1_to_nu2(math.pi, nu_i_inb)

    # Leg 3: same as leg 1 (I->P->E via perijove)
    tof3 = tof1

    # (Leg 4 ToF not used: the residual derives leg-5 to close the Laplace period.)

    # Moon phases: Io at 0, Europa at pi/2 as rough starting guess.
    phi_io_seed = 0.0
    phi_eur_seed = math.pi / 2.0  # 90 degrees offset as starting guess

    return np.array([phi_io_seed, phi_eur_seed, tof0, tof1, tof2, tof3])


# --- Multi-seed grid search --------------------------------------------------


def search_gipeipe(n_seeds: int = 80) -> BallisticGipeipe | None:
    """Search for a closed GIPEIPE by refining a grid of seeds.

    Tries both ie_prograde=True/False and a grid of Io/Europa initial phase angles.
    Returns the best (lowest ballistic_resnorm_kms) result, or None if none converges.

    This is the discovery function; the documented converged seed is stored in
    :data:`GIPEIPE_SEED` once found.
    """
    best: BallisticGipeipe | None = None
    best_norm = 1e9

    # Base ToF estimates from the conic-seed geometry
    conic = _gipeipe_conic_seed()
    base_tofs = list(conic[2:])

    n_phase = max(2, int(math.sqrt(n_seeds / 2)))
    for ie_prog in (False, True):
        for i_phi_io in range(n_phase):
            for i_phi_eur in range(n_phase):
                phi_io = 2.0 * math.pi * i_phi_io / n_phase
                phi_eur = 2.0 * math.pi * i_phi_eur / n_phase
                # Try multiple ToF distributions
                for scale in (0.7, 1.0, 1.3):
                    t1, t2, t3, t4 = [s * scale for s in base_tofs]
                    # Make sure total <= T_LAPLACE (leave room for leg 5)
                    total4 = t1 + t2 + t3 + t4
                    if total4 >= T_LAPLACE_S * 0.95:
                        t1 = t2 = t3 = t4 = T_LAPLACE_S * 0.95 / 4.0
                    seed = np.array([phi_io, phi_eur, t1, t2, t3, t4])
                    result = refine_gipeipe(seed, ie_prograde=ie_prog, max_nfev=800)
                    if result is not None and result.ballistic_resnorm_kms < best_norm:
                        best_norm = result.ballistic_resnorm_kms
                        best = result

    return best


# --- Documented converged seed (#493 search result) ---

#: Converged GIPEIPE seed from the #493 grid search (phi_io, phi_eur, tof1..4 in seconds).
#: Produced by ``search_gipeipe(n_seeds=80)`` — deterministic from this seed.
#: Result: resnorm ~ 6e-14 km/s (fully closed), seam ~ 3e-14 km/s (fully closed),
#: all interior flyby altitudes sub-surface (characterised negative).
GIPEIPE_SEED: NDArray[np.float64] = np.array(
    [
        -1.693676708103792,  # phi_io (rad)
        -4.174112350949415,  # phi_eur (rad)
        123216.05223477464,  # tof1 G→I (s)
        57184.344880417615,  # tof2 I→E (s)
        214326.82013619927,  # tof3 E→I (s)
        61259.32954571807,  # tof4 I→E (s)
    ]
)


__all__ = [
    "GIPEIPE_ECC",
    "GIPEIPE_PERIOD_D",
    "GIPEIPE_PERIOD_S",
    "GIPEIPE_SEED",
    "GIPEIPE_SEQUENCE",
    "GIPEIPE_SMA_KM",
    "IEG_SINGLE_SEQUENCE",
    "LL2011_GIPEIPE_ORBIT_PERIOD_DAYS",
    "LL2011_IEG_BALLISTIC_EUROPA_ALT_KM",
    "LL2011_IEG_DV_MS",
    "LL2011_LAPLACE_PERIOD_DAYS",
    "T_LAPLACE_D",
    "T_LAPLACE_S",
    "BallisticGipeipe",
    "construct_gipeipe",
    "gipeipe_period_summary",
    "refine_gipeipe",
    "search_gipeipe",
    "single_period_ieg_summary",
]


def construct_gipeipe() -> BallisticGipeipe:
    """Construct the documented GIPEIPE from the converged #493 seed — deterministic.

    Refines :data:`GIPEIPE_SEED` to a fully closed (resnorm ~ 6e-14 km/s) GIPEIPE.
    Result: all interior flyby altitudes sub-surface (characterised negative in the
    ideal 2D coplanar model; consistent with L-L's MALTO-optimised powered solution).
    """
    result = refine_gipeipe(GIPEIPE_SEED, ie_prograde=False, max_nfev=4000)
    if result is None:
        raise ValueError("GIPEIPE refine from documented seed failed")
    return result
