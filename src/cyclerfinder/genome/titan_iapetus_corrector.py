"""Titan-Iapetus 3D eccentric-Keplerian closure corrector (Saturn system).

Productized form of the throwaway probes ``scripts/probe_572_titan_iapetus_3d_closure.py``
(idealized inclined-circular 3D closure), ``scripts/run_573_titan_iapetus_population_closure.py``
(the 22-branch population census), and ``scripts/run_574_titan_iapetus_eccentric_kill_gate.py``
(the eccentric-Keplerian kill gate, which found 15 of those 22 branches survive real Titan/
Iapetus eccentricity). See ``data/OUTSTANDING.md`` task #574 (Stage B) for the authoritative
spec this module implements; this docstring only summarizes.

What this is
------------
A Lambert-leg 3-body-flyby closure engine for the Titan-Iapetus-Titan sequence about Saturn,
under an IDEALIZED model that is nonetheless a real step up from the project's usual
circular-coplanar moon ephemeris (:func:`cyclerfinder.search.discovery_campaign._moon_state`):

* Titan is kept in Saturn's equatorial (xy) reference plane (``inclination=0``), matching every
  #558-lineage script's "Titan stays equatorial, Iapetus goes 3D" convention.
* Iapetus is placed on a real ECCENTRIC orbit inclined by :data:`INCLINATION_DEG` (15.5 deg,
  documented conservative estimate -- ``core/satellites.py`` carries no inclination field) to
  Titan's plane, with an explicit ascending-node longitude ``Omega`` (measured from Titan's
  fixed t=0 position).
* BOTH moons carry real, non-negligible eccentricity (:data:`ECC_TITAN` ~= 0.0288,
  :data:`ECC_IAPETUS` ~= 0.028 -- JPL SSD Planetary Satellite Mean Orbital Parameters,
  ssd.jpl.nasa.gov/sats/elem/, sourced in the #574 Stage-A spec text). This is the load-bearing
  fidelity step over #572/#573's inclined-but-CIRCULAR model: Titan/Iapetus eccentricity is
  7-25x the Uranian moons' e<=0.004 that the #312 Uranian quasi_cycler family tolerated, so it
  was the one genuinely unresolved risk in the whole #571->#573 stack.

Free-parameter contract (C1 discipline -- #574 Fable correction)
------------------------------------------------------------------
A closure branch is parameterized by EXACTLY the 4 free parameters of
:class:`TitanIapetusClosureParams`, plus the fixed per-branch ``n_rev`` (and the
near-constant ``inclination_deg``/``e_titan``/``e_iapetus``):

* ``omega_deg`` -- Iapetus's ascending-node longitude (RAAN), the node-alignment free
  variable #572 introduced.
* ``tof_scale`` -- leg time-of-flight as a multiple of ``sqrt(P_Titan * P_Iapetus)``
  (both legs of the closed Titan->Iapetus->Titan cycle use the SAME ToF by construction).
* ``m0_titan_deg`` / ``m0_iapetus_deg`` -- each moon's MEAN ANOMALY AT EPOCH t=0. Every
  later encounter state (Titan's second state at t=2*tof; Iapetus's state at t=tof) is
  derived by Kepler-propagating (mean-motion) from these SAME two epoch values -- there is
  NO free "phase at the encounter" re-specification. This is deliberate: a free
  per-encounter phase is exactly the #480 EGGIE per-encounter self-consistency bug
  (project memory ``feedback_constructed_tour_per_encounter_self_consistency``), and #574's
  Fable plan review (C1) made this a MANDATORY discipline for this exact corrector.

Argument of periapsis is FIXED at 0 for BOTH moons (periapsis at the Omega=0 reference
direction / at the ascending node) -- an explicit, documented simplification (real apsidal
precession is not modeled at all in this idealized layer; the real-ephemeris SPICE V4-strict
gate in :mod:`cyclerfinder.data.validation.v4_saturn_strict` is what tests the real
non-precession-free geometry), NOT a free/searched parameter.

Positive control (C2 discipline)
---------------------------------
:func:`kepler_state_3d` MUST reduce exactly to the circular-coplanar
:func:`cyclerfinder.search.discovery_campaign._moon_state` at ``ecc=inc=0`` and to the
circular-INCLINED formula (the same ``R3(Omega).R1(inc)`` rotation #572's ``iapetus_state_3d``
used) at ``ecc=0``. Both reductions are pinned as regression tests in
``tests/genome/test_titan_iapetus_corrector.py`` (ported from #574's own
``_smoke_test_kepler_reduction``, which passed at a grid of M0/Omega/u test points, dr<1e-6 km,
dv<1e-9 km/s).

Physical gate
-------------
:func:`closure_passes_gate` reuses :func:`cyclerfinder.search.physical_sanity.
candidate_passes_physical_gate` VERBATIM (the project-wide #324 physical-gate machinery,
not reimplemented) -- a candidate closure must ALSO deliver a useful (>= 5 deg by default)
ballistic bend at Iapetus at a safe periapsis altitude, on top of the raw V_inf-continuity
residual gate.

Framing (mandatory, carried from #571-#574): any output of this module is quasi-cycler-class
evidence about our own idealized, internally-enumerated search space -- same standing as
#312's own Uranian family -- NOT a ballistic-cycler finding and NOT a novelty claim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day
from cyclerfinder.search.physical_sanity import (
    DEFAULT_MIN_USEFUL_BEND_DEG,
    FlybyPhysicalVerdict,
    candidate_passes_physical_gate,
)

PRIMARY = "Saturn"
ANCHOR = "Titan"
FLYBY = "Iapetus"
SEQUENCE: tuple[str, str, str] = (ANCHOR, FLYBY, ANCHOR)

INCLINATION_DEG = 15.5
"""Iapetus's inclination (deg) to Titan's orbital plane, held fixed (not searched).

Documented conservative estimate cited by #572 (``core/satellites.py`` carries no
inclination field for any Saturnian moon)."""

ECC_TITAN = 0.0288
"""Titan mean orbital eccentricity.

Source: JPL SSD Planetary Satellite Mean Orbital Parameters
(ssd.jpl.nasa.gov/sats/elem/), as sourced in the #574 Stage-A spec text
(``data/OUTSTANDING.md`` #574)."""

ECC_IAPETUS = 0.028
"""Iapetus mean orbital eccentricity. Same source as :data:`ECC_TITAN`."""

GATE_RESIDUAL_KMS = 0.05
"""Project-wide V_inf-continuity closure-gate residual floor (km/s).

Matches the #285 / #312 SILVER gate and every #558-lineage script's
``GATE_RESIDUAL_KMS`` (e.g. ``scripts/scan_558_uranus_all_pairs_offset_sweep.py``) --
a candidate closure counts as "closed" iff its worst V_inf-continuity residual across
all encounters (including the wrap back to the anchor) is <= this value."""


# --------------------------------------------------------------------------- #
# Free-parameter contract
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TitanIapetusClosureParams:
    """One Titan-Iapetus-Titan 3D eccentric-Keplerian closure branch.

    See the module docstring's "Free-parameter contract" section for the C1 discipline
    this dataclass encodes: exactly 4 free numbers (``omega_deg``, ``tof_scale``,
    ``m0_titan_deg``, ``m0_iapetus_deg``) plus the fixed per-branch ``n_rev`` and the
    near-constant eccentricity/inclination fields (overridable for e.g. the C2 e=0
    positive-control check, but NOT free search parameters in the production gate).
    """

    omega_deg: float
    """Iapetus ascending-node longitude (RAAN), degrees, measured from Titan's fixed
    t=0 reference phase."""
    tof_scale: float
    """Leg time-of-flight as a multiple of ``sqrt(P_Titan * P_Iapetus)`` (both legs of
    the closed cycle share this same ToF)."""
    n_rev: tuple[int, int]
    """Per-leg Lambert revolution count: ``(Titan->Iapetus, Iapetus->Titan)``."""
    m0_titan_deg: float
    """Titan's mean anomaly at epoch t=0, degrees."""
    m0_iapetus_deg: float
    """Iapetus's mean anomaly at epoch t=0, degrees."""
    e_titan: float = ECC_TITAN
    e_iapetus: float = ECC_IAPETUS
    inclination_deg: float = INCLINATION_DEG


@dataclass(frozen=True)
class ClosureResult:
    """Outcome of evaluating one :class:`TitanIapetusClosureParams` branch."""

    params: TitanIapetusClosureParams
    residual_kms: float
    """Worst V_inf-continuity residual across the two legs + the closed-cycle wrap
    (km/s); ``inf`` if ``lambert_infeasible``."""
    vinf_kms: tuple[float, float, float] | None
    """``(V_inf at Titan departure, V_inf at Iapetus, V_inf at Titan arrival)``, km/s;
    ``None`` if ``lambert_infeasible``."""
    lambert_infeasible: bool
    infeasible_reason: str | None = None
    """One of ``"geometry"``, ``"convergence"``, ``"infeasible_n_rev"``, or ``None``."""

    @property
    def closes(self) -> bool:
        """``residual_kms <= GATE_RESIDUAL_KMS`` (the raw closure gate, WITHOUT the
        #324 physical bend gate -- see :func:`closure_passes_gate` for the full gate)."""
        return (not self.lambert_infeasible) and self.residual_kms <= GATE_RESIDUAL_KMS


# --------------------------------------------------------------------------- #
# Eccentric-Keplerian 3D state propagation
# --------------------------------------------------------------------------- #


def _solve_kepler_e(
    mean_anomaly_rad: float, ecc: float, *, tol: float = 1e-13, max_iter: int = 60
) -> float:
    """Newton-Raphson solve of Kepler's equation ``E - e*sin(E) = M``."""
    m = mean_anomaly_rad % (2.0 * math.pi)
    e_anom = m if ecc < 0.8 else math.pi
    for _ in range(max_iter):
        f = e_anom - ecc * math.sin(e_anom) - m
        fp = 1.0 - ecc * math.cos(e_anom)
        d = f / fp
        e_anom -= d
        if abs(d) < tol:
            break
    return e_anom


def kepler_state_3d(
    m0_rad: float,
    n_rad_day: float,
    t_days: float,
    sma_km: float,
    mu: float,
    ecc: float,
    raan_rad: float,
    inc_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Eccentric-Keplerian moon position/velocity (Saturn frame, km, km/s).

    Mean anomaly is propagated by MEAN MOTION from the epoch value ``m0_rad`` (Kepler
    III mean motion depends only on ``sma_km``, not ``ecc``, so reusing the circular
    code's own ``n_rad_day`` is exact -- the C1 discipline: no free per-encounter phase,
    only epoch M0 + elapsed time). Argument of periapsis is fixed at 0 (periapsis at the
    ``Omega=0`` reference direction); ``raan_rad=inc_rad=0`` (Titan's case in this module)
    reduces this to a pure in-plane ellipse with no rotation, i.e. the perifocal frame IS
    the Saturn-frame xy plane.

    At ``ecc=0`` this reduces EXACTLY to the circular case (``E=nu=M``, ``r=sma_km``) --
    verified against :func:`cyclerfinder.search.discovery_campaign._moon_state` (at
    ``inc=raan=0``) and the standard circular-inclined ``R3(Omega).R1(inc)`` rotation (at
    ``inc!=0``) by ``tests/genome/test_titan_iapetus_corrector.py`` (the C2 positive
    control; ported from #574's own ``_smoke_test_kepler_reduction``, which passed at a
    grid of M0/Omega/u points, ``dr<1e-6`` km, ``dv<1e-9`` km/s).
    """
    m_t = m0_rad + n_rad_day * t_days
    e_anom = _solve_kepler_e(m_t, ecc)
    cos_e = math.cos(e_anom)
    nu = 2.0 * math.atan2(
        math.sqrt(1.0 + ecc) * math.sin(e_anom / 2.0),
        math.sqrt(1.0 - ecc) * math.cos(e_anom / 2.0),
    )
    r = sma_km * (1.0 - ecc * cos_e)
    p = sma_km * (1.0 - ecc * ecc)
    cos_nu, sin_nu = math.cos(nu), math.sin(nu)
    px = r * cos_nu
    py = r * sin_nu
    v_scale = math.sqrt(mu / max(p, 1e-9))
    vx_pf = -v_scale * sin_nu
    vy_pf = v_scale * (ecc + cos_nu)

    cos_o, sin_o = math.cos(raan_rad), math.sin(raan_rad)
    cosi, sini = math.cos(inc_rad), math.sin(inc_rad)

    def _rot(px_: float, py_: float) -> np.ndarray:
        x = cos_o * px_ - sin_o * cosi * py_
        y = sin_o * px_ + cos_o * cosi * py_
        z = sini * py_
        return np.array([x, y, z])

    pos = _rot(px, py)
    vel = _rot(vx_pf, vy_pf)
    return pos, vel


def titan_state(params: TitanIapetusClosureParams, t_days: float) -> tuple[np.ndarray, np.ndarray]:
    """Titan's (position km, velocity km/s) at ``t_days`` since epoch, Saturn frame.

    Titan stays in Saturn's equatorial plane (``raan=inc=0``), Kepler-propagated from
    ``params.m0_titan_deg`` by mean motion.
    """
    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    m0 = math.radians(params.m0_titan_deg)
    return kepler_state_3d(m0, n_a, t_days, sat_a.sma_km, mu, params.e_titan, 0.0, 0.0)


def iapetus_state(
    params: TitanIapetusClosureParams, t_days: float
) -> tuple[np.ndarray, np.ndarray]:
    """Iapetus's (position km, velocity km/s) at ``t_days`` since epoch, Saturn frame.

    Iapetus is inclined by ``params.inclination_deg`` with ascending node
    ``params.omega_deg``, Kepler-propagated from ``params.m0_iapetus_deg`` by mean motion.
    """
    mu = PRIMARIES[PRIMARY]
    sat_b = SATELLITES[FLYBY]
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    m0 = math.radians(params.m0_iapetus_deg)
    omega = math.radians(params.omega_deg)
    inc = math.radians(params.inclination_deg)
    return kepler_state_3d(m0, n_b, t_days, sat_b.sma_km, mu, params.e_iapetus, omega, inc)


def cycle_period_days(tof_scale: float) -> float:
    """Total Titan->Iapetus->Titan cycle period (2 legs of ``tof_scale *
    sqrt(P_Titan*P_Iapetus)`` each), days."""
    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    sat_b = SATELLITES[FLYBY]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof = tof_scale * math.sqrt(p_a * p_b)
    return 2.0 * tof


def leg_tof_days(tof_scale: float) -> float:
    """Single-leg time-of-flight, days (``tof_scale * sqrt(P_Titan * P_Iapetus)``)."""
    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    sat_b = SATELLITES[FLYBY]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    return tof_scale * math.sqrt(p_a * p_b)


# --------------------------------------------------------------------------- #
# Lambert-closure evaluation
# --------------------------------------------------------------------------- #


def _leg_best(
    r_a: np.ndarray,
    v_a: np.ndarray,
    r_b: np.ndarray,
    v_b: np.ndarray,
    tof_s: float,
    mu: float,
    n_rev: int,
) -> dict[str, float] | None:
    """Solve one Lambert leg for the exact requested ``n_rev``; return vinf_out/vinf_in.

    Returns ``None`` if that revolution count is infeasible (Lambert did not admit a
    solution at ``n_rev``). Ported verbatim from #572's ``_leg_best`` /
    #574's identical helper.
    """
    sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
    cands = [s for s in sols if s.n_revs == n_rev]
    if not cands:
        return None
    best = min(cands, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
    return {
        "vinf_out": float(np.linalg.norm(best.v1 - v_a)),
        "vinf_in": float(np.linalg.norm(best.v2 - v_b)),
    }


def evaluate_closure(params: TitanIapetusClosureParams) -> ClosureResult:
    """Evaluate one Titan-Iapetus-Titan closure branch under the eccentric 3D model.

    Titan's state at t=0 AND t=2*tof are BOTH Kepler-propagated from the SAME
    ``params.m0_titan_deg`` (C1 -- no free re-specification at the second encounter,
    mirroring #573's own doubled-TOF ``evaluate_point_tracked`` pattern). Iapetus's
    state at t=tof is Kepler-propagated from ``params.m0_iapetus_deg``.

    Residual is the worst V_inf-magnitude continuity defect across both legs
    (``|vinf_in(leg0) - vinf_out(leg1)|`` at the Iapetus midpoint,
    ``|vinf_out(leg0) - vinf_in(leg1)|`` at the closed-cycle Titan wrap).
    """
    mu = PRIMARIES[PRIMARY]
    tof = leg_tof_days(params.tof_scale)
    tof_s = tof * DAY_S

    r0, v0 = titan_state(params, 0.0)
    r1, v1 = iapetus_state(params, tof)
    r2, v2 = titan_state(params, 2.0 * tof)

    n0, n1 = params.n_rev
    try:
        leg0 = _leg_best(r0, v0, r1, v1, tof_s, mu, n0)
        leg1 = _leg_best(r1, v1, r2, v2, tof_s, mu, n1)
    except LambertGeometryError:
        return ClosureResult(params, float("inf"), None, True, "geometry")
    except LambertConvergenceError:
        return ClosureResult(params, float("inf"), None, True, "convergence")
    if leg0 is None or leg1 is None:
        return ClosureResult(params, float("inf"), None, True, "infeasible_n_rev")

    r_mid = abs(leg0["vinf_in"] - leg1["vinf_out"])
    r_periodic = abs(leg0["vinf_out"] - leg1["vinf_in"])
    residual = max(r_mid, r_periodic)

    vinf0 = leg0["vinf_out"]
    vinf1 = max(leg0["vinf_in"], leg1["vinf_out"])
    vinf2 = leg1["vinf_in"]

    return ClosureResult(params, residual, (vinf0, vinf1, vinf2), False, None)


def closure_passes_gate(
    result: ClosureResult,
    *,
    min_useful_bend_deg: float = DEFAULT_MIN_USEFUL_BEND_DEG,
) -> tuple[bool, list[FlybyPhysicalVerdict] | None]:
    """Full closure gate: raw residual gate AND the #324 physical bend gate.

    Reuses :func:`cyclerfinder.search.physical_sanity.candidate_passes_physical_gate`
    VERBATIM (not reimplemented) -- see the module docstring's "Physical gate" section.
    Returns ``(passes, per_encounter_verdicts)``; ``per_encounter_verdicts`` is ``None``
    if the raw residual gate already failed (the physical gate is not evaluated).
    """
    if not result.closes or result.vinf_kms is None:
        return False, None
    gate_pass, verdicts = candidate_passes_physical_gate(
        SEQUENCE, result.vinf_kms, min_useful_bend_deg=min_useful_bend_deg
    )
    return bool(gate_pass), verdicts


# --------------------------------------------------------------------------- #
# JSON-friendly export helper (for scripts / gauntlet runners)
# --------------------------------------------------------------------------- #


def params_to_jsonable(params: TitanIapetusClosureParams) -> dict[str, Any]:
    """``dict`` view of :class:`TitanIapetusClosureParams` for JSONL logging."""
    return {
        "omega_deg": params.omega_deg,
        "tof_scale": params.tof_scale,
        "n_rev": list(params.n_rev),
        "m0_titan_deg": params.m0_titan_deg,
        "m0_iapetus_deg": params.m0_iapetus_deg,
        "e_titan": params.e_titan,
        "e_iapetus": params.e_iapetus,
        "inclination_deg": params.inclination_deg,
    }


__all__ = [
    "ANCHOR",
    "ECC_IAPETUS",
    "ECC_TITAN",
    "FLYBY",
    "GATE_RESIDUAL_KMS",
    "INCLINATION_DEG",
    "PRIMARY",
    "SEQUENCE",
    "ClosureResult",
    "TitanIapetusClosureParams",
    "closure_passes_gate",
    "cycle_period_days",
    "evaluate_closure",
    "iapetus_state",
    "kepler_state_3d",
    "leg_tof_days",
    "params_to_jsonable",
    "titan_state",
]
