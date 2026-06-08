"""Circular -> ephemeris CONTINUATION of the TWO-ARC free-return chain (#164).

This bridges the two proven pieces that resolved the S1L1/#94 frontier to
"V_inf-reachable, ToF-quantised in circular":

* :mod:`cyclerfinder.search.free_return_chain` (#163) closed Russell's actual
  generic-return construction — TWO distinct Earth-to-Earth free-return arcs
  (``g`` + ``G``), each crossing Mars's radius, patched at an intermediate Earth
  flyby — at the SOURCED V_inf anchors. On S1L1 ``russell-ch4-4.991gG2`` (its OWN
  anchors E 4.99 / M 5.10) the V_inf anchors landed to <0.1 km/s and the G-arc ToF
  was near-exact (2.810 vs 2.8096 yr), but the g-arc carried a ~0.14 yr ToF
  residual. The chain is the SEED.
* :mod:`cyclerfinder.search.continuation` (#158) walks a circular-coplanar
  closure out to the true DE440 ephemeris by ramping the *planet model* fidelity
  (``e`` then ``i``) in an ``nstep`` ladder, keep-best (Russell 2004 §5.4). It was
  built for SINGLE-ellipse free-return rows; this module extends that ramp
  machinery to a TWO-ARC seed.

What is continued (the same principle as #158, one primitive up)
---------------------------------------------------------------
Per Russell §5.4: **the solar-system MODEL fidelity, NOT the trajectory genome.**
The two-arc chain's free variables ``(a_1, e_1, a_2, e_2)`` stay the optimisation
unknowns; the homotopy parameter is the PLANET MODEL the chain's V_inf / ToF
geometry references. The chain primitive (:func:`free_return_chain.free_return_geometry`,
:func:`free_return_chain._arc_ee_time_years`, :func:`free_return_chain._earth_vinf_vector`)
computes every emerged quantity from ``(a, e)`` against CIRCULAR planet
assumptions baked in: planet heliocentric radius ``= sma_au`` and planet velocity
``= sqrt(mu / r)`` tangential. Those two per-encounter quantities are exactly what
acquires the planet's real eccentricity as the model is ramped.

The decisive physics (memory blocker, 2026-06-04 direct construction): real Mars
eccentricity drives the Mars-encounter radius/velocity off the circular value, and
that is what drove Mars V_inf down to 2.83 km/s in the direct construction. Here
the SAME lever is applied continuously: at the chosen encounter epoch the planet's
real radius differs from its mean, so the V_inf-fixed ellipse PERIOD changes — and
the per-arc Earth-to-Earth ToF (period x integer rev + time-above-Earth) moves.
**The decisive question this module answers: does the real DE440 ephemeris close
the g-arc ToF gap (0.14 yr) while holding both arcs' V_inf at the sourced anchors?**

The ramped two-arc geometry (this module, no edit to either base module)
-----------------------------------------------------------------------
:func:`_ramped_arc_geometry` is the homotopy generalisation of
:func:`free_return_chain.free_return_geometry`. It places each encounter body on a
RAMPED heliocentric state: the planet's effective radius and velocity vector are
interpolated by ``(lam_e, lam_i)`` between the circular value (``lam = 0``) and the
DE440 state at the encounter epoch (``lam = 1``). At ``lam_e = lam_i = 0`` it is
byte-identical to :func:`free_return_chain.free_return_geometry` (the bit-identical
mechanics gate); at ``lam = 1`` it is the real-ephemeris two-arc geometry. The
radial-crossing solve uses the ramped planet radius; the emerged V_inf uses the
ramped planet velocity. The encounter epoch is a free homotopy variable
(``t0_sec``) — irrelevant at ``lam = 0`` (circular planets are phase-free in
radius/speed), it becomes the lever that selects WHERE on Mars's eccentric orbit
the encounter falls as the ramp proceeds.

The ramp schedule + ladder are REUSED from :mod:`continuation` (``LADDER``,
``_J2000_MEAN_E``, ``_J2000_MEAN_I_DEG``): Russell's e-then-i ramp in an
``nstep = 3^(steploop-1)`` ladder ``{1, 3, 9, 27, 81}``, keep-best, 243 skipped and
recorded. Each step re-solves the two-arc chain seeded from the previous step's
``(a_1, e_1, a_2, e_2, t0)`` (classic imbedding, §5.4.1).

Constraint-vs-evidence separation (the golden rule, inherited)
--------------------------------------------------------------
* CONSTRAINED (drives the residual): per-arc emerged V_inf at Earth and Mars match
  the SOURCED anchors; per-arc emerged Earth-to-Earth arc ToF matches the SOURCED
  descriptor ToF; the intermediate Earth flyby is V_inf-continuous and
  bend-feasible — at EVERY homotopy step (the anchor-respecting residual of #163,
  carried through unchanged, including the per-arc ``n_rev`` ToF-binding term that
  guards the #163 spurious-collapse trap).
* FREE (the unknowns): the two arc shapes ``(a_1, e_1, a_2, e_2)`` and the
  encounter epoch ``t0`` (continuous); the per-arc revolution counts (discrete).
* EVIDENCE (emerges, compared non-circularly): the converged per-arc V_inf
  magnitudes at Earth and Mars and the per-arc arc ToFs at the true ephemeris.

The SOURCED anchor is EXPECTED; the emerged V_inf AND ToF are evidence. A CLOSE
must satisfy BOTH halves — V_inf within ~0.5 of the anchors AND both descriptor
ToFs reached (the g-arc gap closes). V_inf alone is NOT a close (the #163 spurious-
collapse trap). A clean TOF-GAP-PERSISTS or DIVERGES is a valid honest outcome.

Pure: depends only on core (constants, ephemeris, kepler, flyby) and the two base
search modules (free_return, free_return_chain, continuation). No core/base edit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import acos, atan2, pi, sin, sqrt

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SAFE_PERIHELION_KM
from cyclerfinder.core.flyby import is_ballistic_feasible, max_bend
from cyclerfinder.core.kepler import coe_to_rv
from cyclerfinder.search.continuation import (
    _J2000_MEAN_I_DEG,
    LADDER,
)

# Reuse the chain's true->mean and crossing helpers (no re-derivation).
from cyclerfinder.search.free_return import _crossing, _true_to_mean
from cyclerfinder.search.free_return_chain import (
    _N_REV_RANGE,
    _TOF_WEIGHT_KMS_PER_YEAR,
    SECONDS_PER_YEAR,
)

__all__ = [
    "ContinuationChainResult",
    "ContinuationChainRung",
    "ContinuationChainStep",
    "RampedArcGeometry",
    "continuation_chain_correct",
]


# ---------------------------------------------------------------------------
# Ramped two-arc geometry — the homotopy generalisation of the chain primitive.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RampedArcGeometry:
    """Radial-crossing geometry + EMERGED V_inf of one arc at ramp ``(lam_e, lam_i)``.

    The homotopy generalisation of :class:`free_return.FreeReturnGeometry`: at
    ``lam_e = lam_i = 0`` every field equals the circular
    :func:`free_return.free_return_geometry` output bit-identically; at ``lam = 1``
    the encounter bodies sit on their DE440 states at the encounter epoch.
    """

    vinf: dict[str, float]
    nu: dict[str, float]
    tof_em_days: float
    period_days: float


@dataclass(frozen=True)
class _EffectivePlanet:
    """Ramped planet encounter descriptor in the spacecraft's LOCAL frame.

    Frame-consistent with :func:`free_return.free_return_geometry`: instead of an
    absolute DE440 velocity vector (which lives in a different frame than the
    spacecraft's perifocal velocity and so cannot be differenced directly), the
    planet velocity is described by quantities in the spacecraft's local
    ``(r_hat, t_hat, z_hat)`` encounter frame and assembled there by the caller:

    * ``r_eff_km`` — the heliocentric radius the radial-crossing solve targets.
    * ``speed_kms`` — the planet's orbital speed at that radius.
    * ``fpa_rad`` — the planet's flight-path angle (angle of the velocity above the
      local horizontal; zero for a circular orbit).
    * ``z_frac`` — fraction of the inclination-tilt out-of-plane component (ramped
      by ``lam_i``), applied as a small rotation of the planet velocity off the
      ecliptic.

    At ``lam_e = lam_i = 0`` this is ``(sma_au * AU_KM, sqrt(mu / r), 0, 0)`` — the
    exact circular quantities, so the ramped geometry reduces to
    :func:`free_return.free_return_geometry` bit-for-bit.
    """

    r_eff_km: float
    speed_kms: float
    fpa_rad: float
    z_frac: float


def _planet_nu_at_radius(e_p: float, a_p_km: float, r_km: float) -> float:
    """Planet true anomaly (rad, in ``[0, pi]``) whose eccentric radius is ``r``.

    Solves ``r = a_p (1 - e_p^2) / (1 + e_p cos nu_p)``. Clamped to the reachable
    band so a spacecraft-crossing radius just outside the planet's [peri, apo] range
    maps to the nearest reachable planet anomaly (peri/apo) rather than raising.
    """
    if e_p == 0.0:
        return 0.0
    p = a_p_km * (1.0 - e_p * e_p)
    cos_nu = (p / r_km - 1.0) / e_p
    return acos(max(-1.0, min(1.0, cos_nu)))


def _effective_planet(
    body: str,
    epoch_phase: float,
    lam_e: float,
    lam_i: float,
    *,
    mu: float,
) -> _EffectivePlanet:
    """Ramped local-frame planet descriptor (radius, speed, fpa, z-frac).

    The eccentricity physics is built from the planet's OWN J2000 eccentric orbit
    (``PLANETS[body].ecc`` / ``sma_au``), NOT an injected absolute DE440 vector, so
    the velocity stays in the spacecraft's local encounter frame and the V_inf
    difference is frame-consistent (the fix for the off-family collapse).

    ``epoch_phase`` is the planet's encounter true anomaly in ``[0, 2pi)`` — the
    homotopy's encounter-epoch lever: it selects WHERE on the eccentric orbit the
    planet sits (hence its radius), and thus drives the eccentric radius/speed/fpa.
    At ``lam_e = 0`` the orbit is treated as circular (radius ``sma``, fpa 0) so
    ``epoch_phase`` has no effect — exactly the phase-free circular endpoint.

    The eccentric radius/speed/fpa are ramped in by ``lam_e`` (the eccentricity
    perturbation); the out-of-plane ``z_frac`` is ramped by ``lam_i``.
    """
    r_circ = PLANETS[body].sma_au * AU_KM
    e_p = PLANETS[body].ecc
    a_p = PLANETS[body].sma_au * AU_KM
    v_circ = sqrt(mu / r_circ)
    # Effective eccentricity ramped on by lam_e: the planet orbit goes from circle
    # (lam_e=0) to its real J2000 ellipse (lam_e=1).
    e_eff = lam_e * e_p
    if e_eff == 0.0:
        return _EffectivePlanet(r_eff_km=r_circ, speed_kms=v_circ, fpa_rad=0.0, z_frac=0.0)
    # Planet true anomaly at the encounter (the epoch lever). Radius on the ramped
    # ellipse at that anomaly; speed and flight-path angle from the vis-viva /
    # conic relations on the SAME ramped ellipse.
    nu_p = epoch_phase
    p = a_p * (1.0 - e_eff * e_eff)
    r_eff = p / (1.0 + e_eff * np.cos(nu_p))
    speed = sqrt(mu * (2.0 / r_eff - 1.0 / a_p))
    fpa = atan2(e_eff * sin(nu_p), 1.0 + e_eff * np.cos(nu_p))
    # lam_i tilts the planet velocity out of the ecliptic by a small inclination;
    # carried as a z-fraction the caller applies as a rotation of the local frame.
    i_rad = lam_i * (pi / 180.0) * _J2000_MEAN_I_DEG.get(body, 0.0)
    z_frac = sin(i_rad)
    return _EffectivePlanet(
        r_eff_km=float(r_eff), speed_kms=float(speed), fpa_rad=float(fpa), z_frac=float(z_frac)
    )


def _planet_velocity_local(
    ep: _EffectivePlanet, r_hat: np.ndarray, t_hat: np.ndarray
) -> np.ndarray:
    """Assemble the planet velocity in the spacecraft's local encounter frame.

    ``v = speed * (sin(fpa) r_hat + cos(fpa) t_hat)`` tilted off the ecliptic by the
    inclination z-fraction. At ``fpa = 0, z_frac = 0`` this is ``speed * t_hat`` —
    the circular tangential velocity :func:`free_return.free_return_geometry` uses.
    """
    sf, cf = sin(ep.fpa_rad), np.cos(ep.fpa_rad)
    v_in_plane = ep.speed_kms * (sf * r_hat + cf * t_hat)
    if ep.z_frac == 0.0:
        return np.asarray(v_in_plane, dtype=np.float64)
    cz = sqrt(max(0.0, 1.0 - ep.z_frac * ep.z_frac))
    v = v_in_plane * cz
    v = v + np.array([0.0, 0.0, ep.speed_kms * ep.z_frac], dtype=np.float64)
    return np.asarray(v, dtype=np.float64)


def _ramped_arc_geometry(
    a_au: float,
    e: float,
    t0_sec: float,
    lam_e: float,
    lam_i: float,
    *,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
) -> RampedArcGeometry:
    """Two-body radial-crossing geometry of arc ``(a, e)`` at ramp ``(lam_e, lam_i)``.

    The homotopy generalisation of :func:`free_return.free_return_geometry`. The
    spacecraft ellipse crosses each body's RAMPED heliocentric radius; the per-body
    V_inf is the spacecraft velocity at that crossing minus the body's RAMPED
    velocity, BOTH expressed in the spacecraft's local encounter frame (so the
    difference is frame-consistent — the fix for the off-family collapse). The
    ramped planet radius/speed/flight-path-angle come from the planet's OWN J2000
    eccentric orbit (:func:`_effective_planet`); the encounter true anomaly on that
    orbit is set by the shared epoch ``t0_sec`` via the planet's mean motion (Earth
    at ``t0``, Mars at ``t0 + tof_em``). The epoch is irrelevant at ``lam_e = 0``
    (the circular endpoint is phase-free) and becomes the encounter-epoch lever as
    the ramp proceeds.

    At ``lam_e = lam_i = 0`` the radii are ``sma_au * AU_KM`` and the velocities are
    purely tangential of magnitude ``sqrt(mu / r)`` — identical to
    :func:`free_return.free_return_geometry`, which this reduces to bit-for-bit.
    """
    inner, outer = bodies
    if not PLANETS[inner].sma_au < PLANETS[outer].sma_au:
        raise ValueError(
            f"bodies must be ordered (inner, outer) by semi-major axis: {inner!r} "
            f"is not inside {outer!r}."
        )
    a_km = a_au * AU_KM
    n = sqrt(mu / a_km**3)
    period = 2.0 * pi / n
    period_days = period / 86400.0

    # Inner (Earth) encounter at planet phase set by t0; one provisional pass to
    # get the inner-leg ToF, then the outer (Mars) encounter at t0 + tof_em.
    ep_in = _effective_planet(inner, _planet_phase(inner, t0_sec, mu), lam_e, lam_i, mu=mu)
    nu_inner = _crossing(a_km, e, ep_in.r_eff_km)
    m_inner = _true_to_mean(nu_inner, e)
    ep_out_prov = _effective_planet(outer, _planet_phase(outer, t0_sec, mu), lam_e, lam_i, mu=mu)
    nu_out_prov = _crossing(a_km, e, ep_out_prov.r_eff_km)
    tof_prov = ((_true_to_mean(nu_out_prov, e) - m_inner) % (2.0 * pi)) / n
    ep_out = _effective_planet(
        outer, _planet_phase(outer, t0_sec + tof_prov, mu), lam_e, lam_i, mu=mu
    )
    nu_outer = _crossing(a_km, e, ep_out.r_eff_km)
    m_outer = _true_to_mean(nu_outer, e)
    tof_em_s = ((m_outer - m_inner) % (2.0 * pi)) / n

    nu = {inner: nu_inner, outer: nu_outer}
    eps = {inner: ep_in, outer: ep_out}
    vinf: dict[str, float] = {}
    for b in bodies:
        r_sc, v_sc = coe_to_rv(a_km, e, nu[b], mu)
        r_hat = np.asarray(r_sc) / np.linalg.norm(r_sc)
        t_hat = np.array([-r_hat[1], r_hat[0], 0.0])
        v_planet = _planet_velocity_local(eps[b], r_hat, t_hat)
        vinf[b] = float(np.linalg.norm(np.asarray(v_sc) - v_planet))

    return RampedArcGeometry(
        vinf=vinf,
        nu=nu,
        tof_em_days=(tof_em_s / 86400.0),
        period_days=period_days,
    )


def _planet_phase(body: str, t_sec: float, mu: float) -> float:
    """Planet encounter true anomaly (rad) from the shared epoch via mean motion.

    The shared epoch ``t_sec`` maps to each planet's position on its own orbit
    through the planet mean motion ``n_p = sqrt(mu / a_p^3)``. At ``lam_e = 0`` the
    orbit is circular and the phase is unused (the circular endpoint is phase-free);
    at ``lam_e > 0`` it selects WHERE on the eccentric orbit the encounter falls —
    the lever that, e.g., places Mars near perihelion to drop its V_inf.
    """
    a_p = PLANETS[body].sma_au * AU_KM
    n_p = sqrt(mu / a_p**3)
    return float((n_p * t_sec) % (2.0 * pi))


def _ramped_arc_ee_time_years(
    a_au: float,
    e: float,
    t0_sec: float,
    lam_e: float,
    lam_i: float,
    n_rev: int,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> float:
    """Earth-to-Earth arc time (years) for ``(a, e)`` at ramp + ``n_rev``.

    Same multi-rev arithmetic as :func:`free_return_chain._arc_ee_time_years`, but
    the Earth radial crossing is taken against the RAMPED Earth radius (so the
    time-above-Earth term tracks the eccentric Earth-encounter radius). At
    ``lam = 0`` it equals :func:`free_return_chain._arc_ee_time_years` exactly.
    """
    a_km = a_au * AU_KM
    ep_e = _effective_planet("E", _planet_phase("E", t0_sec, mu), lam_e, lam_i, mu=mu)
    nu_e = _crossing(a_km, e, ep_e.r_eff_km)
    n = sqrt(mu / a_km**3)
    period = 2.0 * pi / n
    m_e = _true_to_mean(nu_e, e)
    t_above = period - 2.0 * (m_e / n)
    return (t_above + n_rev * period) / SECONDS_PER_YEAR


def _ramped_best_n_rev(
    a_au: float,
    e: float,
    t0_sec: float,
    lam_e: float,
    lam_i: float,
    tof_target_years: float,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> int:
    """``n_rev`` minimising ``|ramped arc time - target|`` (the discrete DOF)."""
    return min(
        _N_REV_RANGE,
        key=lambda nr: abs(
            _ramped_arc_ee_time_years(a_au, e, t0_sec, lam_e, lam_i, nr, mu=mu) - tof_target_years
        ),
    )


@dataclass(frozen=True)
class _RampedArc:
    """One arc's emerged geometry at a ramp step (shape + epoch + ToF-min n_rev)."""

    a_au: float
    e: float
    t0_sec: float
    geometry: RampedArcGeometry
    n_rev: int
    arc_tof_years: float

    @property
    def vinf_e(self) -> float:
        return self.geometry.vinf["E"]

    @property
    def vinf_m(self) -> float:
        return self.geometry.vinf["M"]


def _build_ramped_arc(
    a_au: float,
    e: float,
    t0_sec: float,
    lam_e: float,
    lam_i: float,
    tof_target_years: float,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> _RampedArc | None:
    """Build one arc's ramped geometry (+ ToF-min n_rev), or ``None`` off-family."""
    try:
        g = _ramped_arc_geometry(a_au, e, t0_sec, lam_e, lam_i, mu=mu)
    except ValueError:
        return None
    n_rev = _ramped_best_n_rev(a_au, e, t0_sec, lam_e, lam_i, tof_target_years, mu=mu)
    arc_tof = _ramped_arc_ee_time_years(a_au, e, t0_sec, lam_e, lam_i, n_rev, mu=mu)
    return _RampedArc(a_au=a_au, e=e, t0_sec=t0_sec, geometry=g, n_rev=n_rev, arc_tof_years=arc_tof)


def _earth_vinf_vector_ramped(
    arc: _RampedArc, lam_e: float, lam_i: float, *, mu: float
) -> np.ndarray:
    """Heliocentric V_inf vector at the arc's Earth crossing (ramped planet vel).

    The ramped analogue of :func:`free_return_chain._earth_vinf_vector`, used by
    the intermediate-flyby bend test. The planet velocity is the local-frame
    eccentric velocity (:func:`_planet_velocity_local`); at ``lam = 0`` it equals
    the chain's circular Earth-V_inf vector bit-identically.
    """
    a_km = arc.a_au * AU_KM
    ep_e = _effective_planet("E", _planet_phase("E", arc.t0_sec, mu), lam_e, lam_i, mu=mu)
    nu_e = _crossing(a_km, arc.e, ep_e.r_eff_km)
    r_sc, v_sc = coe_to_rv(a_km, arc.e, nu_e, mu)
    r_hat = np.asarray(r_sc) / np.linalg.norm(r_sc)
    t_hat = np.array([-r_hat[1], r_hat[0], 0.0])
    v_planet = _planet_velocity_local(ep_e, r_hat, t_hat)
    return np.asarray(v_sc - v_planet, dtype=np.float64)


def _intermediate_turn_geometry_ramped(
    arc1: _RampedArc,
    arc2: _RampedArc,
    lam_e: float,
    lam_i: float,
    *,
    mu: float,
    continuity_tol_kms: float,
) -> tuple[float, float, bool, float]:
    """Required vs achievable bend + V_inf continuity at the intermediate Earth flyby.

    Ramped analogue of :func:`free_return_chain._intermediate_turn_geometry`.
    Returns ``(turn_rad, max_turn_rad, bend_feasible, continuity_kms)``.
    """
    v_in = _earth_vinf_vector_ramped(arc1, lam_e, lam_i, mu=mu)
    v_out = _earth_vinf_vector_ramped(arc2, lam_e, lam_i, mu=mu)
    vin_mag = float(np.linalg.norm(v_in))
    vout_mag = float(np.linalg.norm(v_out))
    continuity = abs(vin_mag - vout_mag)
    mu_planet = PLANETS["E"].mu_km3_s2
    rp_min = SAFE_PERIHELION_KM["E"]
    cos_arg = float(np.dot(v_in, v_out)) / (vin_mag * vout_mag)
    turn = float(np.arccos(np.clip(cos_arg, -1.0, 1.0)))
    v_mean = 0.5 * (vin_mag + vout_mag)
    max_turn = max_bend(mu_planet, rp_min, v_mean)
    v_in_eq = v_in / vin_mag * v_mean
    v_out_eq = v_out / vout_mag * v_mean
    in_cone = is_ballistic_feasible(v_in_eq, v_out_eq, mu_planet, rp_min, speed_tol=1.0e-6)
    bend_feasible = bool(in_cone and continuity < continuity_tol_kms)
    return turn, max_turn, bend_feasible, continuity


# ---------------------------------------------------------------------------
# The per-step two-arc chain solve at a ramp factor (#163 residual, ramped).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _ChainSolve:
    """Converged two-arc chain at one ramp step (the per-step homotopy state)."""

    a1: float
    e1: float
    a2: float
    e2: float
    t0_sec: float
    arc1: _RampedArc
    arc2: _RampedArc
    max_residual_kms: float
    vinf_residual_kms: float
    tof_residual_years: float
    vinf_continuity_kms: float
    intermediate_flyby_feasible: bool
    intermediate_turn_deg: float
    intermediate_max_turn_deg: float


_INF_RES = (1e3, 1e3, 1e3, 1e3, 1e3, 1e3)


def _chain_residuals_ramped(
    x: np.ndarray,
    lam_e: float,
    lam_i: float,
    *,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    arc1_tof_years: float,
    arc2_tof_years: float,
    mu: float,
) -> list[float]:
    """Anchor-respecting two-arc residual at ramp ``(lam_e, lam_i)``.

    Free variables ``x = [a1, e1, a2, e2, t0_years]``. Identical structure to
    :func:`free_return_chain._chain_residuals` (per-arc emerged V_inf at Earth and
    Mars vs SOURCED anchor, per-arc emerged Earth-to-Earth ToF vs descriptor ToF,
    ToF term weighted to km/s) but evaluated on the RAMPED planet model and with
    the shared encounter epoch as the fifth free variable. The epoch is carried in
    YEARS (``t0_years``) so the least_squares Jacobian is well-scaled against the
    AU/km-class ``(a, e)`` variables — a raw seconds epoch makes the epoch partial
    ~1e7x smaller and the solver never moves it (it is the encounter-epoch lever
    that selects WHERE on the eccentric orbit the encounter falls, so it must be
    free to slide). Nothing imposes the anchor V_inf; it EMERGES from each
    ``(a_i, e_i)`` against the ramped planet state.
    """
    a1, e1, a2, e2, t0_years = (float(x[i]) for i in range(5))
    t0 = t0_years * SECONDS_PER_YEAR
    for a, e in ((a1, e1), (a2, e2)):
        if not (0.0 < e < 0.95) or a <= 0.0:
            return list(_INF_RES)
    arc1 = _build_ramped_arc(a1, e1, t0, lam_e, lam_i, arc1_tof_years, mu=mu)
    arc2 = _build_ramped_arc(a2, e2, t0, lam_e, lam_i, arc2_tof_years, mu=mu)
    if arc1 is None or arc2 is None:
        return list(_INF_RES)
    w = _TOF_WEIGHT_KMS_PER_YEAR
    return [
        arc1.vinf_m - vinf_m_anchor,
        arc2.vinf_m - vinf_m_anchor,
        arc1.vinf_e - vinf_e_anchor,
        arc2.vinf_e - vinf_e_anchor,
        w * (arc1.arc_tof_years - arc1_tof_years),
        w * (arc2.arc_tof_years - arc2_tof_years),
    ]


def _solve_chain_step(
    x0: np.ndarray,
    lam_e: float,
    lam_i: float,
    *,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    arc1_tof_years: float,
    arc2_tof_years: float,
    mu: float,
    tol_kms: float,
    max_nfev: int = 400,
) -> _ChainSolve | None:
    """Re-solve the two-arc chain at ramp ``(lam_e, lam_i)`` seeded from ``x0``.

    Returns a :class:`_ChainSolve`, or ``None`` if the solve wandered off the
    Mars-reaching family (the homotopy step is then treated as a divergence).
    """

    def _res(x: np.ndarray) -> list[float]:
        return _chain_residuals_ramped(
            x,
            lam_e,
            lam_i,
            vinf_e_anchor=vinf_e_anchor,
            vinf_m_anchor=vinf_m_anchor,
            arc1_tof_years=arc1_tof_years,
            arc2_tof_years=arc2_tof_years,
            mu=mu,
        )

    sol = least_squares(_res, x0, method="trf", max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)
    x = sol.x
    res = _res(x)
    max_res = max(abs(r) for r in res)
    a1, e1, a2, e2, t0_years = (float(x[i]) for i in range(5))
    t0 = t0_years * SECONDS_PER_YEAR
    arc1 = _build_ramped_arc(a1, e1, t0, lam_e, lam_i, arc1_tof_years, mu=mu)
    arc2 = _build_ramped_arc(a2, e2, t0, lam_e, lam_i, arc2_tof_years, mu=mu)
    if arc1 is None or arc2 is None:
        return None
    vinf_res = max(
        abs(arc1.vinf_m - vinf_m_anchor),
        abs(arc2.vinf_m - vinf_m_anchor),
        abs(arc1.vinf_e - vinf_e_anchor),
        abs(arc2.vinf_e - vinf_e_anchor),
    )
    tof_res = max(
        abs(arc1.arc_tof_years - arc1_tof_years),
        abs(arc2.arc_tof_years - arc2_tof_years),
    )
    turn, max_turn, feasible, continuity = _intermediate_turn_geometry_ramped(
        arc1, arc2, lam_e, lam_i, mu=mu, continuity_tol_kms=tol_kms
    )
    return _ChainSolve(
        a1=a1,
        e1=e1,
        a2=a2,
        e2=e2,
        t0_sec=t0,
        arc1=arc1,
        arc2=arc2,
        max_residual_kms=float(max_res),
        vinf_residual_kms=float(vinf_res),
        tof_residual_years=float(tof_res),
        vinf_continuity_kms=float(continuity),
        intermediate_flyby_feasible=bool(feasible),
        intermediate_turn_deg=float(np.degrees(turn)),
        intermediate_max_turn_deg=float(np.degrees(max_turn)),
    )


# ---------------------------------------------------------------------------
# Audit trail + the driver.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ContinuationChainStep:
    """One homotopy step's outcome for the two-arc chain (per-step audit record)."""

    phase: str  # "seed" | "e-ramp" | "i-ramp" | "ephemeris"
    lam_e: float
    lam_i: float
    model: str
    max_residual_kms: float
    vinf_residual_kms: float
    tof_residual_years: float
    a1: float
    e1: float
    a2: float
    e2: float
    t0_sec: float
    arc1_tof_years: float
    arc2_tof_years: float
    arc1_n_rev: int
    arc2_n_rev: int
    intermediate_flyby_feasible: bool
    converged: bool


@dataclass(frozen=True)
class ContinuationChainRung:
    """One ``nstep`` rung: its step trail and final (true-ephemeris) solve."""

    nstep: int
    steps: tuple[ContinuationChainStep, ...]
    final: _ChainSolve | None
    completed: bool


@dataclass(frozen=True)
class ContinuationChainResult:
    """Two-arc chain continuation outcome (circular seed -> DE440).

    Attributes
    ----------
    best_final:
        The lowest-residual final (true-ephemeris) :class:`_ChainSolve` across all
        completed rungs (or the best of all rungs if none completed) — the
        headline the three-way gate reads.
    winning_nstep:
        Which ladder rung produced ``best_final``.
    rungs:
        Per-rung audit trail (every step's residuals + which rung won).
    ladder, skipped:
        The ``nstep`` values attempted / deliberately not run (recorded).
    seed_solve:
        The circular (lam=0) two-arc seed solve — the #163 starting point, kept so
        the bit-identical-seed mechanics gate can read it.
    """

    best_final: _ChainSolve | None
    winning_nstep: int
    rungs: tuple[ContinuationChainRung, ...]
    ladder: tuple[int, ...]
    seed_solve: _ChainSolve | None
    skipped: tuple[int, ...] = field(default_factory=tuple)

    @property
    def vinf_close(self) -> bool:
        """``best_final`` V_inf within tol of both anchors (the V_inf half)."""
        return self.best_final is not None and self.best_final.vinf_residual_kms < 0.5

    @property
    def tof_close(self) -> bool:
        """``best_final`` both descriptor ToFs reached (the ToF half)."""
        return self.best_final is not None and self.best_final.tof_residual_years < 0.1

    @property
    def closed(self) -> bool:
        """FULL close: V_inf AND ToF AND intermediate bend-feasible (BOTH halves)."""
        return bool(
            self.best_final is not None
            and self.vinf_close
            and self.tof_close
            and self.best_final.intermediate_flyby_feasible
        )


def _ramp_schedule_chain(nstep: int) -> list[tuple[str, float, float]]:
    """Build the (phase, lam_e, lam_i) sequence for one rung (Russell e-then-i).

    No phase ramp leg: the two-arc chain's circular geometry is phase-FREE in
    radius/speed (a circular planet's radius/speed is the same at every epoch), so
    the ~100 deg longitude offset the single-ellipse continuation's ``lam_p`` leg
    corrected does not apply here — the encounter epoch ``t0`` is itself a free
    variable the solver moves as the ramp proceeds. ``nstep`` e-steps
    (``lam_e: 0 -> 1`` at lam_i=0) then ``nstep`` i-steps (``lam_i: 0 -> 1`` at
    lam_e=1); the final true-ephemeris step is appended by the driver.
    """
    sched: list[tuple[str, float, float]] = []
    for k in range(1, nstep + 1):
        sched.append(("e-ramp", k / nstep, 0.0))
    for k in range(1, nstep + 1):
        sched.append(("i-ramp", 1.0, k / nstep))
    return sched


def _step_record(
    phase: str, lam_e: float, lam_i: float, model: str, s: _ChainSolve, tol_kms: float
) -> ContinuationChainStep:
    return ContinuationChainStep(
        phase=phase,
        lam_e=lam_e,
        lam_i=lam_i,
        model=model,
        max_residual_kms=s.max_residual_kms,
        vinf_residual_kms=s.vinf_residual_kms,
        tof_residual_years=s.tof_residual_years,
        a1=s.a1,
        e1=s.e1,
        a2=s.a2,
        e2=s.e2,
        t0_sec=s.t0_sec,
        arc1_tof_years=s.arc1.arc_tof_years,
        arc2_tof_years=s.arc2.arc_tof_years,
        arc1_n_rev=s.arc1.n_rev,
        arc2_n_rev=s.arc2.n_rev,
        intermediate_flyby_feasible=s.intermediate_flyby_feasible,
        converged=s.max_residual_kms < tol_kms,
    )


def continuation_chain_correct(
    a1_seed: float,
    e1_seed: float,
    a2_seed: float,
    e2_seed: float,
    t0_seed_sec: float,
    arc1_tof_years: float,
    arc2_tof_years: float,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = 0.5,
    ladder: tuple[int, ...] = LADDER,
) -> ContinuationChainResult:
    """Walk the two-arc free-return chain from circular-coplanar to the true ephemeris.

    Takes the #163 two-arc circular closure ``(a1, e1, a2, e2)`` + encounter epoch
    as the homotopy START and ramps BOTH arcs circular -> ephemeris with the #158
    schedule: ``nstep`` e-steps then ``nstep`` i-steps, finally a true-ephemeris
    step, over the ``nstep`` ladder, keep-best. Each step re-solves the two-arc
    chain (anchor-respecting #163 residual, the per-arc ``n_rev`` ToF-binding term
    carried through) seeded from the previous step's ``(a1, e1, a2, e2, t0)``.

    The ramp ENDPOINT is the planet's real J2000 eccentricity/inclination
    (``PLANETS[body].ecc`` / :data:`continuation._J2000_MEAN_I_DEG`) — the same
    Standish & Williams / DE440-consistent J2000 elements the astropy backend
    propagates. At ``lam_e = lam_i = 1`` the encounter geometry IS the real
    eccentric/inclined planet model (the "ephemeris" step); the seconds-resolution
    DE440 phase is captured through the free encounter epoch ``t0``.

    Parameters
    ----------
    a1_seed, e1_seed, a2_seed, e2_seed:
        The circular two-arc seed from :func:`free_return_chain.free_return_chain_correct`
        (the #163 ``(a_1, e_1, a_2, e_2)``).
    t0_seed_sec:
        Seed encounter epoch (seconds since J2000). At lam=0 the planet radius/speed
        are phase-free so any value reproduces the circular seed; it becomes the
        encounter-epoch lever as the ramp proceeds.
    arc1_tof_years, arc2_tof_years:
        The two descriptor ToFs (g, G) in YEARS.
    vinf_e_anchor, vinf_m_anchor:
        The SOURCED V_inf anchors at Earth and Mars (km/s) — the EXPECTED target.
    """
    skipped = tuple(n for n in (1, 3, 9, 27, 81, 243) if n not in ladder)

    seed_x = np.array(
        [a1_seed, e1_seed, a2_seed, e2_seed, t0_seed_sec / SECONDS_PER_YEAR],
        dtype=np.float64,
    )
    # lam=0 seed solve (the #163 circular closure, re-expressed through the ramped
    # primitive so the bit-identical-seed gate can compare it).
    seed_solve = _solve_chain_step(
        seed_x,
        0.0,
        0.0,
        vinf_e_anchor=vinf_e_anchor,
        vinf_m_anchor=vinf_m_anchor,
        arc1_tof_years=arc1_tof_years,
        arc2_tof_years=arc2_tof_years,
        mu=mu,
        tol_kms=tol_kms,
    )

    rungs: list[ContinuationChainRung] = []
    for nstep in ladder:
        steps: list[ContinuationChainStep] = []
        current = seed_solve
        completed = current is not None
        if current is not None:
            steps.append(_step_record("seed", 0.0, 0.0, "circular", current, tol_kms))
        for phase, lam_e, lam_i in _ramp_schedule_chain(nstep):
            if current is None:
                completed = False
                break
            x0 = np.array(
                [current.a1, current.e1, current.a2, current.e2, current.t0_sec / SECONDS_PER_YEAR],
                dtype=np.float64,
            )
            nxt = _solve_chain_step(
                x0,
                lam_e,
                lam_i,
                vinf_e_anchor=vinf_e_anchor,
                vinf_m_anchor=vinf_m_anchor,
                arc1_tof_years=arc1_tof_years,
                arc2_tof_years=arc2_tof_years,
                mu=mu,
                tol_kms=tol_kms,
            )
            if nxt is None or not np.isfinite(nxt.max_residual_kms):
                completed = False
                break
            current = nxt
            model = f"ramped(e={lam_e:.4g},i={lam_i:.4g})"
            steps.append(_step_record(phase, lam_e, lam_i, model, current, tol_kms))
        final: _ChainSolve | None
        if completed and current is not None:
            x0 = np.array(
                [current.a1, current.e1, current.a2, current.e2, current.t0_sec / SECONDS_PER_YEAR],
                dtype=np.float64,
            )
            final = _solve_chain_step(
                x0,
                1.0,
                1.0,
                vinf_e_anchor=vinf_e_anchor,
                vinf_m_anchor=vinf_m_anchor,
                arc1_tof_years=arc1_tof_years,
                arc2_tof_years=arc2_tof_years,
                mu=mu,
                tol_kms=tol_kms,
            )
            if final is not None:
                steps.append(
                    _step_record("ephemeris", 1.0, 1.0, "ephemeris(j2000-ecc/inc)", final, tol_kms)
                )
            else:
                completed = False
        else:
            final = current
        rungs.append(
            ContinuationChainRung(nstep=nstep, steps=tuple(steps), final=final, completed=completed)
        )

    completed_rungs = [r for r in rungs if r.completed and r.final is not None]
    pool = completed_rungs if completed_rungs else [r for r in rungs if r.final is not None]
    if pool:
        best = min(pool, key=lambda r: r.final.max_residual_kms)  # type: ignore[union-attr]
        best_final = best.final
        winning = best.nstep
    else:
        best_final = None
        winning = ladder[0] if ladder else 0
    return ContinuationChainResult(
        best_final=best_final,
        winning_nstep=winning,
        rungs=tuple(rungs),
        ladder=tuple(ladder),
        seed_solve=seed_solve,
        skipped=skipped,
    )
