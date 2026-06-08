"""CORRECTED-topology S1L1 / 4.991gG2 construction (#167).

The redemption of #164/#165. Task #166 (``docs/notes/2026-06-08-s1l1-source-dig.md``)
found two things at once:

1. **The topology was wrong.** In the descriptor ``g(1.4612) + G(2.8096)`` only the
   UPPERCASE ``G`` is the Mars-transit leg; the lowercase ``g`` is a *pure
   Earth-to-Earth free return that never reaches Mars* (DE440-verified closest
   approach ~1.05 AU, aphelion ~1.27 AU — well below Mars). The #164/#165 two-arc
   construction modelled BOTH arcs as Mars-crossing — the structural error behind the
   #165 ~2.6 AU drift. The correct per-cycle sequence is::

       E -> g(Earth-Earth free return, NO Mars) -> E flyby
         -> G(Earth-Mars-Earth transit, longitude rendezvous with true DE440 Mars) -> E

2. **The exact real-eph seed is in hand.** Russell 2004 Appendix C #83 prints the
   per-leg ``(epoch, time-start, v_inf-vector)`` block that reproduces the whole
   30-year, 7-cycle trajectory. Its Mars-transit leg reconstructs on DE440 to ~1.7 km
   with a 4-decimal v_inf match and a true longitude rendezvous (Mars at ecliptic
   longitude 201.0deg on 2027-06-13). This is the longitude-rendezvous constraint
   #165 omitted.

This module encodes that App-C block verbatim (:data:`APPC_4991GG2`) and provides the
reconstruction primitive Russell's own recipe specifies: at each leg start, the
spacecraft heliocentric velocity is ``v_planet(DE440) + v_inf``, then the arc is
propagated (two-body-Sun here; the n-body cross-check lives in the test). The
construction is therefore **seeded entirely from sourced numbers** — nothing here
fits anything; the achieved Mars miss / v_inf are evidence measured against Russell's
printed values.

Topology (verbatim from the App-C E/M column, §1 of the source-dig note)::

    E E M | E E M | E E M | E E M | E E M | E E M | E E M | E
    1 2 3   4 5 6   7 8 9  10 ...                              22

i.e. three encounters per cycle (E -> E -> M), one Mars encounter per cycle. The
``E -> E`` legs (1->2, 4->5, ...) are the ``g`` free returns (sub-Mars aphelion); the
``E -> M`` legs (2->3, 5->6, ...) are the ``G`` Mars-transit outbound arcs.

Provenance / honesty
--------------------
The EXPECTED side (per-leg v_inf, transit time, the 5.10/3.05 caveat) traces entirely
to Russell's printed App-C / Table 4.9 / Table 5.5; the App-C block below is a verbatim
transcription. The real-eph Mars v_inf BREATHES 3.2-8.0 km/s across the seven cycles
(avg 5.48) — it is NOT the coplanar 5.10 anchor nor the Rogers/CPOM 3.05. Any gate
must target the real-eph per-leg values, not those coplanar idealisations.

Pure: depends only on core/constants, core/kepler, core/ephemeris.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.constants import AU_KM, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate as kepler_propagate

Vec3 = np.ndarray

DAY_S = SECONDS_PER_DAY

# --- Russell 2004 Appendix C #83 "DATA NECESSARY TO REPRODUCE" block, verbatim ---
#
# Model: DE405, J2000-equatorial states rotated to ecliptic by obliquity
# 0.409092629205 rad (== our _J2000_OBLIQUITY_RAD to 1e-5 deg). The published v_inf
# vectors are ALREADY in the ecliptic frame, so reconstruction on our astropy/DE440
# ecliptic backend is direct: v_sc = v_planet + v_inf, no extra rotation.
#
# "EPOCH TIME (days after J2000)" — the single shared epoch; each leg's "time start"
# is days after this epoch.
APPC_EPOCH_DAYS: float = 9325.742435

# (leg_no, body, time_start_days, vinf_kms_ecliptic). Leg N's arc runs from its own
# time_start to leg (N+1)'s time_start; leg 22 is the terminal Earth node (no v_inf).
APPC_LEGS: tuple[tuple[int, str, float, tuple[float, float, float] | None], ...] = (
    (1, "E", -0.174074986e2, (-0.238703220e0, 0.579754982e1, -0.176632468e-2)),
    (2, "E", 0.519582315e3, (-0.227826290e1, 0.532198991e1, 0.573813715e0)),
    (3, "M", 0.699393400e3, (-0.499623059e1, 0.707421771e0, 0.144166283e1)),
    (4, "E", 0.152997105e4, (-0.575383361e1, 0.386418107e0, 0.628610087e-4)),
    (5, "E", 0.206677179e4, (-0.550681663e1, -0.160954336e1, -0.550229544e0)),
    (6, "M", 0.219229367e4, (-0.355470987e1, -0.672927795e1, 0.112658626e1)),
    (7, "E", 0.310727040e4, (0.130066966e1, -0.353985845e1, 0.127838602e-2)),
    (8, "E", 0.363605992e4, (0.372962888e1, -0.150726402e-1, -0.411042515e0)),
    (9, "M", 0.377449899e4, (0.463453129e1, -0.454933848e0, -0.467039188e-1)),
    (10, "E", 0.468119309e4, (0.398651382e1, 0.380698910e1, 0.599174579e-3)),
    (11, "E", 0.521702941e4, (0.195647577e1, 0.500455176e1, 0.131545295e1)),
    (12, "M", 0.542793047e4, (-0.106091993e1, 0.286375761e1, 0.950243209e0)),
    (13, "E", 0.621934578e4, (-0.346922206e1, 0.601287549e1, -0.294263108e-2)),
    (14, "E", 0.676095373e4, (-0.494762984e1, 0.487023018e1, 0.229687664e0)),
    (15, "M", 0.691560104e4, (-0.580355995e1, -0.130028398e1, 0.196236847e1)),
    (16, "E", 0.777160296e4, (-0.455750356e1, -0.237566423e1, -0.883416315e-3)),
    (17, "E", 0.830586489e4, (-0.289917665e1, -0.412851216e1, -0.880976268e0)),
    (18, "M", 0.841803274e4, (-0.474724840e0, -0.793519915e1, 0.124113571e1)),
    (19, "E", 0.933942338e4, (0.322773002e1, -0.335738552e1, 0.125353647e-2)),
    (20, "E", 0.987176523e4, (0.254872658e1, 0.229406465e1, 0.313437961e1)),
    (21, "M", 0.101019833e5, (0.212149789e1, 0.424849528e0, -0.238347620e1)),
    (22, "E", 0.109670202e5, None),
)

# Russell 2004 App-C "EARTH-TO-MARS TRANSIT LEG CHARACTERISTICS" table — the published
# (transit_days, vinf_Mars) for each Mars-transit (G) outbound leg, keyed by the
# ARRIVAL Mars leg number. EXPECTED side of the gate (sourced, never fit).
APPC_MARS_TRANSIT: dict[int, tuple[float, float]] = {
    3: (179.8, 5.248),
    6: (125.5, 7.693),
    9: (138.4, 4.657),
    12: (210.9, 3.198),
    15: (154.6, 6.263),
    18: (112.2, 8.046),
    21: (230.2, 3.219),
}

# Sourced AVERAGE row (Russell App-C characteristics table) and the coplanar caveat
# anchors. The real-eph cycler reproduces the App-C per-leg / average values, NOT the
# coplanar 5.10 (Russell Table 4.9) or 3.05 (Rogers/CPOM) Mars v_inf idealisations.
APPC_VINF_M_AVG: float = 5.475
APPC_VINF_E_AVG: float = 5.368
COPLANAR_VINF_M: float = 5.10  # Russell Table 4.9 idealisation — NOT a real-eph target
ROGERS_CPOM_VINF_M: float = 3.05  # Rogers/CPOM idealisation — NOT a real-eph target


@dataclass(frozen=True)
class LegArc:
    """One reconstructed App-C leg arc (seeded from a sourced v_inf node)."""

    leg_no: int
    body: str  # departure body code ("E" or "M")
    next_body: str  # arrival body code
    t0_sec: float  # departure epoch, TDB s since J2000
    t1_sec: float  # arrival epoch
    r0_km: Vec3
    v0_km_s: Vec3
    is_g_arc: bool  # True for the E->E free-return (g) leg
    is_mars_transit: bool  # True for the E->M Mars-transit (G) outbound leg


@dataclass(frozen=True)
class MarsEncounter:
    """Achieved Mars encounter of a reconstructed Mars-transit (G) leg."""

    arrival_leg_no: int
    cycle: int
    t_arr_sec: float
    miss_au: float
    vinf_kms: float
    sc_lon_deg: float
    mars_lon_deg: float
    pub_transit_days: float
    pub_vinf_kms: float


@dataclass(frozen=True)
class GArcClearance:
    """Closest approach of one g (E->E) free-return arc to real DE440 Mars."""

    leg_no: int
    cycle: int
    closest_mars_au: float
    aphelion_au: float
    earth_return_miss_au: float


@dataclass(frozen=True)
class ContinuousNode:
    """One encounter of the CONTINUOUS-from-one-seed chain (#169, the V4 measure).

    Unlike :class:`MarsEncounter` / :func:`reconstruct_mars_encounters` — which
    reconstruct each leg INDEPENDENTLY from its own App-C v_inf node (Russell's
    per-leg reproduction recipe, i.e. it RE-ANCHORS v_inf at every node) — this
    records one node of a SINGLE continuous trajectory propagated forward from the
    first Earth departure with NO v_inf re-anchoring between legs. At each node the
    spacecraft arrives ballistically with ``vinf_in_kms``; a real flyby rotates
    that v_inf freely toward the next leg's App-C direction (sourced) and the
    MAINTENANCE dv is the part a ballistic flyby cannot supply: the |v_inf|
    magnitude change ``dv_mag_kms`` plus, if the required bend exceeds the
    safe-periapsis maximum, the un-bent shortfall ``dv_bend_kms``.
    """

    arrival_leg_no: int
    body: str
    is_mars: bool
    miss_km: float
    vinf_in_kms: float
    vinf_appc_kms: float
    bend_deg: float
    max_bend_deg: float
    dv_mag_kms: float
    dv_bend_kms: float

    @property
    def dv_total_kms(self) -> float:
        """Maintenance dv at this node (magnitude change + un-bendable shortfall)."""
        return self.dv_mag_kms + self.dv_bend_kms


def continuous_chain(
    ephem: Ephemeris,
    *,
    perturbers: tuple[str, ...] = (),
    propagate: object | None = None,
    start_leg_no: int = 2,
) -> list[ContinuousNode]:
    """Propagate ONE continuous trajectory through the App-C nodes (#169 V4).

    Starts from the first Earth-departure App-C state (``start_leg_no``, default
    leg 2 — the 2026-12-15 Mars-transit departure used in the #167 results note)
    and walks forward node-to-node WITHOUT re-anchoring v_inf. Between nodes the
    state is propagated either by two-body-Sun Kepler (``propagate is None``) or by
    a supplied n-body propagator (a ``RestrictedNBody`` whose ``.propagate`` is
    called; ``perturbers`` then names the rails perturber bodies) — the test wires
    REBOUND/IAS15. At each node the FLYBY PATCH convention is applied: the position
    is snapped to the real DE440 planet (the few-thousand-km ballistic miss is
    ``<< SOI`` — the standard patched-conic boundary, recorded as evidence) and the
    velocity is rotated to the next leg's App-C v_inf direction (the free flyby
    turn). The maintenance dv (magnitude change + any un-bendable shortfall) is the
    measured horizon-TCM budget that the per-leg re-anchoring hides.

    EXPECTED = the App-C per-leg v_inf (sourced); ``miss_km`` / ``vinf_in_kms`` /
    the dv are EVIDENCE measured by the integrator, never imposed.
    """
    from cyclerfinder.core.flyby import max_bend

    nodes = [((APPC_EPOCH_DAYS + ts) * DAY_S, body, vinf) for (_no, body, ts, vinf) in APPC_LEGS]
    start_idx = next(i for i, leg in enumerate(APPC_LEGS) if leg[0] == start_leg_no)

    t0, body0, vinf0 = nodes[start_idx]
    r_p, v_p = ephem.state(body0, t0)
    r = np.asarray(r_p, dtype=np.float64)
    v = np.asarray(v_p, dtype=np.float64) + _vinf_vec(vinf0)

    out: list[ContinuousNode] = []
    for i in range(start_idx + 1, len(nodes)):
        t1, body_i, vinf_out = nodes[i]
        if propagate is None:
            rk, vk = kepler_propagate(r, v, t1 - t0)
        else:
            arc = propagate.propagate(  # type: ignore[attr-defined]
                r,
                v,
                t0,
                t1,
                bodies=perturbers,
                ephem=ephem if perturbers else None,
                accuracy=1e-11 if not perturbers else 1e-10,
            )
            rk = np.asarray(arc.r_km, dtype=np.float64)
            vk = np.asarray(arc.v_km_s, dtype=np.float64)

        r_pl_t, v_pl_t = ephem.state(body_i, t1)
        r_pl = np.asarray(r_pl_t, dtype=np.float64)
        v_pl = np.asarray(v_pl_t, dtype=np.float64)
        miss_km = float(np.linalg.norm(rk - r_pl))
        vinf_in = vk - v_pl
        vinf_in_mag = float(np.linalg.norm(vinf_in))

        dv_mag = 0.0
        dv_bend = 0.0
        bend_deg = 0.0
        max_bend_deg = 0.0
        vinf_appc_mag = vinf_in_mag
        if vinf_out is not None:
            vinf_out_vec = _vinf_vec(vinf_out)
            vinf_appc_mag = float(np.linalg.norm(vinf_out_vec))
            dv_mag = abs(vinf_appc_mag - vinf_in_mag)
            cos_b = float(np.dot(vinf_in, vinf_out_vec) / (vinf_in_mag * vinf_appc_mag))
            bend = float(np.arccos(max(min(cos_b, 1.0), -1.0)))
            bend_deg = float(np.degrees(bend))
            pdata = PLANETS[body_i]
            v_for_bend = min(vinf_in_mag, vinf_appc_mag)
            mb = max_bend(pdata.mu_km3_s2, pdata.radius_eq_km + pdata.safe_alt_km, v_for_bend)
            max_bend_deg = float(np.degrees(mb))
            if bend > mb:
                dv_bend = float(2.0 * v_for_bend * np.sin((bend - mb) / 2.0))
            # Maintained departure: flyby patch (snap to planet) + App-C v_inf turn.
            r = r_pl.copy()
            v = v_pl + vinf_out_vec
        else:
            r = rk
            v = vk

        out.append(
            ContinuousNode(
                arrival_leg_no=APPC_LEGS[i][0],
                body=body_i,
                is_mars=body_i == "M",
                miss_km=miss_km,
                vinf_in_kms=vinf_in_mag,
                vinf_appc_kms=vinf_appc_mag,
                bend_deg=bend_deg,
                max_bend_deg=max_bend_deg,
                dv_mag_kms=dv_mag,
                dv_bend_kms=dv_bend,
            )
        )
        t0 = t1
    return out


def _vinf_vec(vinf: tuple[float, float, float] | None) -> Vec3:
    assert vinf is not None
    return np.array(vinf, dtype=np.float64)


def build_seeded_arcs(ephem: Ephemeris) -> list[LegArc]:
    """Build the 21 reconstructed App-C leg arcs, seeded from the sourced v_inf nodes.

    Russell's reproduction recipe (App-C, verbatim): at each leg start take the planet
    velocity from the ephemeris and set the spacecraft velocity to
    ``v_planet + v_inf``; propagate to the next node. This places the departure node at
    the REAL DE440 planet position, so the construction is anchored to the true
    ephemeris — the longitude-rendezvous constraint #165 lacked.
    """
    arcs: list[LegArc] = []
    legs = APPC_LEGS
    for i in range(len(legs) - 1):
        leg_no, body, ts, vinf = legs[i]
        _next_no, next_body, next_ts, _ = legs[i + 1]
        t0 = (APPC_EPOCH_DAYS + ts) * DAY_S
        t1 = (APPC_EPOCH_DAYS + next_ts) * DAY_S
        r_p, v_p = ephem.state(body, t0)
        r0 = np.asarray(r_p, dtype=np.float64)
        v0 = np.asarray(v_p, dtype=np.float64) + _vinf_vec(vinf)
        is_g = body == "E" and next_body == "E"
        is_g_transit = body == "E" and next_body == "M"
        arcs.append(
            LegArc(
                leg_no=leg_no,
                body=body,
                next_body=next_body,
                t0_sec=t0,
                t1_sec=t1,
                r0_km=r0,
                v0_km_s=v0,
                is_g_arc=is_g,
                is_mars_transit=is_g_transit,
            )
        )
    return arcs


def _lon_deg(r: Vec3) -> float:
    return float(np.degrees(np.arctan2(r[1], r[0])) % 360.0)


def reconstruct_mars_encounters(
    ephem: Ephemeris,
    arcs: list[LegArc] | None = None,
) -> list[MarsEncounter]:
    """Two-body-Sun reconstruction of every Mars-transit (G) leg's encounter.

    Each Mars-transit arc is propagated (Kepler, Sun-only) from its sourced departure
    state to the published Mars-arrival epoch; the achieved miss / v_inf / longitudes
    against real DE440 Mars are measured (evidence). The n-body confirmation (REBOUND/
    IAS15 over DE440 perturbers) lives in the test — this primitive is the cheap
    structural check that the seeds intercept the true Mars at the true longitude.
    """
    if arcs is None:
        arcs = build_seeded_arcs(ephem)
    out: list[MarsEncounter] = []
    for arc in arcs:
        if not arc.is_mars_transit:
            continue
        arrival_no = arc.leg_no + 1
        cycle = (arrival_no - 3) // 3 + 1
        rk, vk = kepler_propagate(arc.r0_km, arc.v0_km_s, arc.t1_sec - arc.t0_sec)
        r_m, v_m = ephem.state("M", arc.t1_sec)
        r_m = np.asarray(r_m, dtype=np.float64)
        v_m = np.asarray(v_m, dtype=np.float64)
        miss_au = float(np.linalg.norm(rk - r_m) / AU_KM)
        vinf_kms = float(np.linalg.norm(vk - v_m))
        pub_tof, pub_vinf = APPC_MARS_TRANSIT[arrival_no]
        out.append(
            MarsEncounter(
                arrival_leg_no=arrival_no,
                cycle=cycle,
                t_arr_sec=arc.t1_sec,
                miss_au=miss_au,
                vinf_kms=vinf_kms,
                sc_lon_deg=_lon_deg(rk),
                mars_lon_deg=_lon_deg(r_m),
                pub_transit_days=pub_tof,
                pub_vinf_kms=pub_vinf,
            )
        )
    return out


def g_arc_clearances(
    ephem: Ephemeris,
    arcs: list[LegArc] | None = None,
    n_samples: int = 200,
) -> list[GArcClearance]:
    """Closest approach to real DE440 Mars over each g (E->E) free-return arc.

    The corrected topology's signature: the g arcs must stay FAR from Mars (sub-Mars
    aphelion ~1.27 AU, closest approach ~1.05 AU per #166). This is the direct refutation
    of the #164/#165 error that forced BOTH arcs to cross Mars's radius.
    """
    if arcs is None:
        arcs = build_seeded_arcs(ephem)
    out: list[GArcClearance] = []
    for arc in arcs:
        if not arc.is_g_arc:
            continue
        cycle = (arc.leg_no - 1) // 3 + 1
        dur = arc.t1_sec - arc.t0_sec
        closest = float("inf")
        aphelion = 0.0
        for k in range(n_samples + 1):
            dt = dur * k / n_samples
            rk, _ = kepler_propagate(arc.r0_km, arc.v0_km_s, dt)
            aphelion = max(aphelion, float(np.linalg.norm(rk) / AU_KM))
            r_m, _ = ephem.state("M", arc.t0_sec + dt)
            d = float(np.linalg.norm(rk - np.asarray(r_m)) / AU_KM)
            closest = min(closest, d)
        r_ret, _ = kepler_propagate(arc.r0_km, arc.v0_km_s, dur)
        r_e, _ = ephem.state("E", arc.t1_sec)
        earth_miss = float(np.linalg.norm(r_ret - np.asarray(r_e)) / AU_KM)
        out.append(
            GArcClearance(
                leg_no=arc.leg_no,
                cycle=cycle,
                closest_mars_au=closest,
                aphelion_au=aphelion,
                earth_return_miss_au=earth_miss,
            )
        )
    return out


__all__ = [
    "APPC_EPOCH_DAYS",
    "APPC_LEGS",
    "APPC_MARS_TRANSIT",
    "APPC_VINF_E_AVG",
    "APPC_VINF_M_AVG",
    "COPLANAR_VINF_M",
    "ROGERS_CPOM_VINF_M",
    "ContinuousNode",
    "GArcClearance",
    "LegArc",
    "MarsEncounter",
    "build_seeded_arcs",
    "continuous_chain",
    "g_arc_clearances",
    "reconstruct_mars_encounters",
]
