"""V-infinity vector extraction for the demanded-turn gate (#888).

Rebuilds the incoming and outgoing V-infinity VECTORS of stored closures so
that :func:`cyclerfinder.verify.turn_gate.demanded_turn_gate` can be applied
without re-deriving geometry by hand. Three representations are covered:

* :func:`symmetric_closure` -- the two-moon ``A -> B -> A`` symmetric closures
  of the ideal circular-coplanar enumerations (#563 Uranus, #575 / #655 Saturn,
  #576 Jupiter, #599 Neptune, #609 Mars). The arcs are rebuilt exactly as
  ``scripts/scan_558_uranus_all_pairs_offset_sweep.py::residual_at_point`` built
  them (moon A at ``phase0`` at t = 0, moon B at ``phase0 + rel_offset``, each leg
  one Lambert arc of time ``tof`` with ``max_revs = max(n_rev, 1)`` and, among the
  solutions with the requested revolution count, the one with the smallest
  ``|v1 - v_moon|``). The third leg (A -> B again, at ``2 tof``) is solved
  directly rather than inferred from symmetry, so the closing flyby at A is
  measured on the trajectory itself.
* :func:`mcconaghy_npr_cycler` -- the Earth-Mars ``nPr`` cyclers of McConaghy,
  Longuski and Byrnes (AIAA 2002-4420), Earth-to-Earth Lambert arcs in the
  circular-coplanar model; the demanded turn at Earth compares the arrival with
  the next departure, which is the same arc rotated by
  ``Delta Psi = 2 pi n S`` (their eq. 4).
* :func:`russell_strange_generic_cycler` -- ideal-model free-return cyclers of
  Russell and Strange (JGCD 32(1), 2009) built entirely from generic
  (``g``/``G``) free-return legs, decoded from their Tables 5-6 nomenclature
  ``g<N>;<theta>;<U|L>``: time of flight ``N`` flyby-body periods, spacecraft
  transfer angle ``theta`` degrees (revolutions ``floor(theta/360)``). The target
  body is massless in their ideal model, so turns are measured only at the
  flyby body, between consecutive legs including the wrap from the last leg
  to the first, in the flyby body's local frame (the ground tracks repeat each
  cycle, Russell and Strange p.149).

All frames are planet-centred (or Sun-centred) inertial, coplanar; units km,
km/s, s.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    LambertSolution,
    lambert,
)
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day
from cyclerfinder.search.discovery_campaign import _moon_state as _dc_moon_state
from cyclerfinder.verify.turn_gate import Encounter, to_body_local

Vec = NDArray[np.float64]
DAY_S = 86400.0

BranchRule = Literal["closest", "low", "high"]


def _circular_state(theta: float, sma_km: float, mu: float) -> tuple[Vec, Vec]:
    """Prograde circular state at longitude ``theta`` (same as discovery_campaign._moon_state)."""
    v = math.sqrt(mu / sma_km)
    return (
        np.array([sma_km * math.cos(theta), sma_km * math.sin(theta), 0.0]),
        np.array([-v * math.sin(theta), v * math.cos(theta), 0.0]),
    )


def _moon_state(moon: str, theta0: float, t_days: float, mu: float) -> tuple[Vec, Vec]:
    """Circular-coplanar moon state, signed for orbital sense (#599).

    The enumeration scripts' own helpers, so the rebuilt geometry is theirs.
    Note: the velocity is the prograde circular velocity even for a retrograde
    moon, exactly as in ``discovery_campaign._moon_state``.
    """
    sat = SATELLITES[moon]
    n = _mean_motion_rad_day(mu, sat.sma_km)
    if sat.retrograde:
        n = -n
    return _dc_moon_state(theta0, n, t_days, sat.sma_km, mu)


def _solve_leg(
    r1: Vec,
    v_body1: Vec,
    r2: Vec,
    tof_s: float,
    mu: float,
    n_rev: int,
    branch: BranchRule = "closest",
) -> LambertSolution | None:
    """One Lambert leg with ``n_rev`` revolutions; ``None`` if no such solution."""
    try:
        sols = lambert(r1, r2, tof_s, mu=mu, max_revs=max(n_rev, 1))
    except (LambertGeometryError, LambertConvergenceError):
        return None
    cands = [s for s in sols if s.n_revs == n_rev]
    if not cands:
        return None
    if branch == "closest" or n_rev == 0:
        return min(cands, key=lambda s: float(np.linalg.norm(s.v1 - v_body1)))
    picked = [s for s in cands if s.branch == branch]
    return picked[0] if picked else None


# ---------------------------------------------------------------------------
# Symmetric two-moon closures (#563 family)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SymmetricClosure:
    """Rebuilt ``A -> B -> A`` closure with its two flybys.

    ``stored_convention_vinf`` is ``[|out0|, max(|in1|, |out1|), |in2|]``, the
    per-encounter magnitudes the enumeration scripts wrote
    (``encounter_vinfs_kms``), for comparison with the stored record.
    """

    primary: str
    anchor: str
    flyby: str
    tof_days: float
    n_rev: tuple[int, int]
    rel_offset_deg: float
    residual_kms: float
    stored_convention_vinf: tuple[float, float, float]
    encounters: tuple[Encounter, Encounter]
    wrap_local_check_deg: float
    """Demanded turn at the anchor recomputed in the anchor's local frame from
    leg 0 (departure at t = 0) and leg 1 (arrival at 2 tof). Equals the direct
    value when the closure repeats in the rotating frame (tof = n T_syn / 2)."""
    vectors: dict[str, list[float]] = field(default_factory=dict)


def symmetric_closure(
    primary: str,
    anchor: str,
    flyby: str,
    *,
    tof_days: float,
    n_rev: tuple[int, int],
    rel_offset_deg: float,
    phase0_deg: float = 0.0,
    branches: tuple[BranchRule, BranchRule, BranchRule] = ("closest", "closest", "closest"),
    alt_floor_km: float | None = None,
) -> SymmetricClosure | None:
    """Rebuild one symmetric closure; ``None`` if a leg has no ``n_rev`` solution.

    ``branches`` selects the Lambert branch of legs 0, 1 and 2 (leg 2 is leg 0
    one cycle later and must use the same rule as leg 0 for a periodic chain).
    """
    mu = PRIMARIES[primary]
    th_a = math.radians(phase0_deg)
    th_b = th_a + math.radians(rel_offset_deg)
    tof_s = tof_days * DAY_S
    r0, w0 = _moon_state(anchor, th_a, 0.0, mu)
    r1, w1 = _moon_state(flyby, th_b, tof_days, mu)
    r2, w2 = _moon_state(anchor, th_a, 2.0 * tof_days, mu)
    r3, _w3 = _moon_state(flyby, th_b, 3.0 * tof_days, mu)
    n0, n1 = n_rev
    leg0 = _solve_leg(r0, w0, r1, tof_s, mu, n0, branches[0])
    leg1 = _solve_leg(r1, w1, r2, tof_s, mu, n1, branches[1])
    leg2 = _solve_leg(r2, w2, r3, tof_s, mu, n0, branches[2])
    if leg0 is None or leg1 is None or leg2 is None:
        return None
    out0 = leg0.v1 - w0
    in1 = leg0.v2 - w1
    out1 = leg1.v1 - w1
    in2 = leg1.v2 - w2
    out2 = leg2.v1 - w2
    m = np.linalg.norm
    residual = max(abs(float(m(in1)) - float(m(out1))), abs(float(m(out0)) - float(m(in2))))
    enc_b = Encounter.for_body(flyby, in1, out1, alt_floor_km=alt_floor_km, label="mid")
    enc_a = Encounter.for_body(anchor, in2, out2, alt_floor_km=alt_floor_km, label="wrap")
    loc_in = to_body_local(in2, r2, w2)
    loc_out = to_body_local(out0, r0, w0)
    cosang = float(loc_in @ loc_out) / (float(m(loc_in)) * float(m(loc_out)))
    wrap_local = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
    return SymmetricClosure(
        primary=primary,
        anchor=anchor,
        flyby=flyby,
        tof_days=tof_days,
        n_rev=(n0, n1),
        rel_offset_deg=rel_offset_deg,
        residual_kms=residual,
        stored_convention_vinf=(float(m(out0)), max(float(m(in1)), float(m(out1))), float(m(in2))),
        encounters=(enc_b, enc_a),
        wrap_local_check_deg=wrap_local,
        vectors={
            "out0": out0.tolist(),
            "in1": in1.tolist(),
            "out1": out1.tolist(),
            "in2": in2.tolist(),
            "out2": out2.tolist(),
        },
    )


# ---------------------------------------------------------------------------
# McConaghy, Longuski and Byrnes (2002) nPr Earth-Mars cyclers
# ---------------------------------------------------------------------------

#: Earth-Mars synodic period in years, S = 2 + 1/7 (McConaghy et al. 2002 p.2).
MCCONAGHY_S_YEARS: float = 15.0 / 7.0


@dataclass(frozen=True)
class NprCycler:
    n: int
    period_class: str
    r: int
    vinf_earth_kms: float
    aphelion_au: float
    period_years: float
    encounter: Encounter


def mcconaghy_npr_cycler(
    n: int, period_class: Literal["L", "S"], r: int, *, alt_floor_km: float | None = None
) -> NprCycler | None:
    """One ``nPr`` cycler of McConaghy et al. (2002) in their ideal model.

    Model (their pp.2-5): Earth on a circular 1 AU orbit (its year from
    ``MU_SUN`` so that R2 lands on Earth exactly), Lambert arc from
    R1 = (1, 0) to R2 = (cos 2 pi n S, sin 2 pi n S) in T = n S years with
    ``r`` complete revolutions; ``L`` / ``S`` is the longer / shorter period of
    the two multi-revolution solutions (``r = 0`` has one solution, ``U``).
    The next arc is the same arc rotated by ``2 pi n S``, so the Earth flyby
    must turn the arrival V-infinity onto the rotated departure V-infinity.
    """
    mu = MU_SUN_KM3_S2
    year_s = 2.0 * math.pi * math.sqrt(AU_KM**3 / mu)
    psi = 2.0 * math.pi * n * MCCONAGHY_S_YEARS
    r1, w1 = _circular_state(0.0, AU_KM, mu)
    r2, w2 = _circular_state(psi, AU_KM, mu)
    try:
        sols = lambert(r1, r2, n * MCCONAGHY_S_YEARS * year_s, mu=mu, max_revs=max(r, 1))
    except (LambertGeometryError, LambertConvergenceError):
        return None
    cands = [s for s in sols if s.n_revs == r]
    if not cands:
        return None

    def period_of(s: LambertSolution) -> float:
        energy = 0.5 * float(np.dot(s.v1, s.v1)) - mu / AU_KM
        a = -mu / (2.0 * energy)
        return 2.0 * math.pi * math.sqrt(a**3 / mu) / year_s

    cands.sort(key=period_of)
    sol = cands[-1] if period_class == "L" or len(cands) == 1 else cands[0]
    energy = 0.5 * float(np.dot(sol.v1, sol.v1)) - mu / AU_KM
    a = -mu / (2.0 * energy)
    h = np.cross(r1, sol.v1)
    e_vec = np.cross(sol.v1, h) / mu - r1 / AU_KM
    ecc = float(np.linalg.norm(e_vec))
    v_out0 = sol.v1 - w1
    v_in = sol.v2 - w2
    c, s_ = math.cos(psi), math.sin(psi)
    rot = np.array([[c, -s_, 0.0], [s_, c, 0.0], [0.0, 0.0, 1.0]])
    v_out_next = rot @ v_out0
    enc = Encounter.for_body(
        "E", v_in, v_out_next, alt_floor_km=alt_floor_km, label=f"{n}{period_class}{r}"
    )
    return NprCycler(
        n=n,
        period_class=period_class,
        r=r,
        vinf_earth_kms=float(np.linalg.norm(v_out0)),
        aphelion_au=a * (1.0 + ecc) / AU_KM,
        period_years=period_of(sol),
        encounter=enc,
    )


# ---------------------------------------------------------------------------
# Russell and Strange (2009) generic-leg ideal-model cyclers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RSBody:
    """Russell and Strange (2009) Table 2 body parameters (p.148)."""

    gm_km3_s2: float
    radius_km: float
    ideal_period_s: float


#: Table 2, "Body parameters" (Russell and Strange 2009, p.148).
RS_PRIMARY_GM: dict[str, float] = {"Jupiter": 126_686_535.0, "Saturn": 37_931_208.0}
RS_BODIES: dict[str, RSBody] = {
    "Io": RSBody(5959.92, 1827.0, 152_854.0),
    "Europa": RSBody(3202.74, 1561.0, 306_822.0),
    "Ganymede": RSBody(9887.83, 2634.0, 618_153.0),
    "Callisto": RSBody(7179.29, 2408.0, 1_441_931.0),
    "Titan": RSBody(8978.14, 2575.0, 1_377_684.0),
    "Enceladus": RSBody(6.95, 256.3, 118_387.0),
}


@dataclass(frozen=True)
class RSLeg:
    """One generic free-return leg: ``g<n_body_revs>;<theta_deg>;<label>``."""

    n_body_revs: float
    theta_deg: float
    label: str


@dataclass(frozen=True)
class RSCyclerRow:
    """A published generic-leg cycler: Tables 3/4 metrics and Tables 5/6 legs."""

    rs_id: str
    primary: str
    flyby_body: str
    vinf_flyby_kms: float
    period_days: float
    min_flyby_alt_km: float
    min_dist_primary_km: float
    max_dist_primary_km: float
    legs: tuple[RSLeg, ...]


@dataclass(frozen=True)
class RSCyclerRebuild:
    row: RSCyclerRow
    leg_vinf_kms: tuple[float, ...]
    leg_rp_km: tuple[float, ...]
    leg_ra_km: tuple[float, ...]
    leg_branch: tuple[str, ...]
    period_days: float
    encounters: tuple[Encounter, ...]


def russell_strange_generic_cycler(row: RSCyclerRow, *, alt_floor_km: float) -> RSCyclerRebuild:
    """Rebuild a generic-leg cycler in the Russell-Strange ideal model.

    Branch selection: among the Lambert solutions with ``floor(theta/360)``
    revolutions, the one whose V-infinity is closest to the PUBLISHED flyby-body
    V-infinity (Tables 3-4). The published altitude is never used to choose.
    """
    mu = RS_PRIMARY_GM[row.primary]
    body = RS_BODIES[row.flyby_body]
    sma = (mu * (body.ideal_period_s / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
    r1, w1 = _circular_state(0.0, sma, mu)
    legs_local: list[tuple[Vec, Vec]] = []
    vinfs: list[float] = []
    rps: list[float] = []
    ras: list[float] = []
    brs: list[str] = []
    for leg in row.legs:
        ang = 2.0 * math.pi * leg.n_body_revs
        n_rev = int(leg.theta_deg // 360.0)
        r2, w2 = _circular_state(ang, sma, mu)
        sols = [
            s
            for s in lambert(
                r1, r2, leg.n_body_revs * body.ideal_period_s, mu=mu, max_revs=max(n_rev, 1)
            )
            if s.n_revs == n_rev
        ]
        if not sols:
            raise ValueError(f"{row.rs_id}: no {n_rev}-rev Lambert solution for leg {leg}")
        sol = min(sols, key=lambda s: abs(float(np.linalg.norm(s.v1 - w1)) - row.vinf_flyby_kms))
        out_loc = to_body_local(sol.v1 - w1, r1, w1)
        in_loc = to_body_local(sol.v2 - w2, r2, w2)
        legs_local.append((out_loc, in_loc))
        vinfs.append(float(np.linalg.norm(out_loc)))
        energy = 0.5 * float(np.dot(sol.v1, sol.v1)) - mu / sma
        a = -mu / (2.0 * energy)
        e = float(np.linalg.norm(np.cross(sol.v1, np.cross(r1, sol.v1)) / mu - r1 / sma))
        rps.append(a * (1.0 - e))
        ras.append(a * (1.0 + e))
        brs.append(sol.branch)
    encs: list[Encounter] = []
    k = len(legs_local)
    for i in range(k):
        v_in = legs_local[i][1]
        v_out = legs_local[(i + 1) % k][0]
        encs.append(
            Encounter(
                body=row.flyby_body,
                mu_km3_s2=body.gm_km3_s2,
                radius_km=body.radius_km,
                alt_floor_km=alt_floor_km,
                vinf_in=v_in,
                vinf_out=v_out,
                label=f"leg{i + 1}->leg{(i + 1) % k + 1}",
            )
        )
    return RSCyclerRebuild(
        row=row,
        leg_vinf_kms=tuple(vinfs),
        leg_rp_km=tuple(rps),
        leg_ra_km=tuple(ras),
        leg_branch=tuple(brs),
        period_days=sum(leg.n_body_revs for leg in row.legs) * body.ideal_period_s / DAY_S,
        encounters=tuple(encs),
    )


def _rs(
    rs_id: str,
    primary: str,
    body: str,
    vinf: float,
    period: float,
    alt: float,
    dmin: float,
    dmax: float,
    legs: Sequence[tuple[float, float, str]],
) -> RSCyclerRow:
    return RSCyclerRow(
        rs_id, primary, body, vinf, period, alt, dmin, dmax, tuple(RSLeg(*leg) for leg in legs)
    )


#: Every Russell-Strange (2009) Table 3/4 cycler built ONLY from generic legs
#: (no resonant ``f``/``h`` legs, whose crank angle is a free parameter).
#: Columns from Table 3 (Jovian, p.149) and Table 4 (Titan-Enceladus, p.150):
#: flyby-body V-infinity, period (days), min flyby altitude at the flyby body (km),
#: min and max distance to the primary (km); legs from Tables 5 and 6 (p.151),
#: whose printed colons are decimal points in the text layer.
RS_GENERIC_CYCLERS: tuple[RSCyclerRow, ...] = (
    _rs("GanCal#5", "Jupiter", "Ganymede", 3.24, 37.6, 328, 821_915, 2_390_844,
        [(1.50425, 541.53130, "L"), (3.74691, 628.88825, "U")]),
    _rs("GanEur#5", "Jupiter", "Ganymede", 1.66, 35.3, 1819, 633_307, 1_080_067,
        [(4.92758, 2493.92898, "U")]),
    _rs("GanEur#43", "Jupiter", "Ganymede", 1.87, 14.1, 8861, 564_558, 1_072_330,
        [(1.97103, 1069.57159, "U")]),
    _rs("GanIo#53", "Jupiter", "Ganymede", 3.90, 21.2, 518, 280_283, 1_075_918,
        [(0.97232, 710.03419, "U"), (1.98423, 1434.32320, "U")]),
    _rs("GanIo#403", "Jupiter", "Ganymede", 4.29, 56.4, 540, 412_959, 1_336_522,
        [(3.72334, 1700.40403, "Ll"), (4.16078, 2217.88234, "L")]),
    _rs("TitEnc#37", "Saturn", "Titan", 2.50, 94.4, 1377, 219_772, 1_296_179,
        [(0.91247, 688.48877, "U"), (5.01017, 3963.66276, "L")]),
    _rs("TitEnc#314", "Saturn", "Titan", 3.59, 49.5, 1218, 196_276, 1_673_806,
        [(1.20284, 433.02184, "U"), (1.89950, 1403.81944, "Ls")]),
    _rs("TitEnc#510", "Saturn", "Titan", 5.03, 80.9, 1852, 196_399, 2_531_474,
        [(2.24732, 449.03519, "U"), (2.82923, 1378.52326, "L")]),
    _rs("TitEnc#552", "Saturn", "Titan", 5.16, 112.4, 3784, 189_939, 2_431_736,
        [(3.22220, 799.99089, "U"), (3.82857, 1738.28475, "L")]),
    _rs("TitEnc#586", "Saturn", "Titan", 5.43, 112.4, 3874, 208_262, 2_899_438,
        [(1.24728, 89.02223, "U"), (5.80348, 2089.25340, "L")]),
)  # fmt: skip

#: Russell and Strange (2009) p.144: the ideal-model Titan flyby floor.
RS_TITAN_MIN_ALT_KM: float = 1000.0


def earth_floor_km() -> float:
    """Project Earth flyby floor (Russell 2004 p.165; McConaghy 2002 p.5 also uses 200 km)."""
    return PLANETS["E"].safe_alt_km
