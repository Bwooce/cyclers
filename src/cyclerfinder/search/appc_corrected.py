"""Generic App-C real-eph reconstruction — scaling the S1L1 pipeline (#170).

:mod:`cyclerfinder.search.s1l1_corrected` hard-codes the Russell 2004 Appendix C #83
(``4.991gG2`` / S1L1) per-leg block, which is fully ballistic (every ``dv`` ~1e-11
km/s). This module generalises the same construction to ANY App-C "DATA NECESSARY TO
REPRODUCE" block, including the *powered* parents whose printed per-leg ``dv`` is a
real maintenance burn (e.g. ``8.049gGf2`` #188 = 420 m/s, ``8.165Gfh-f2`` #192 =
1678 m/s over 7 cycles, Table 5.5).

It is purely additive — S1L1 keeps its dedicated module; this one is the reusable
primitive the batch (#170) drives. The reconstruction recipe is Russell's own
(App-C, verbatim, p.202):

    at each leg start take the planet velocity from the ephemeris and set the
    spacecraft velocity to ``v_planet + v_inf``; propagate by Kepler's equation until
    ``time dv``, when the given Δv is applied; then propagate to the start of the next
    leg / next planet encounter.

For a faithful encounter measurement we reconstruct each leg INDEPENDENTLY from its
own sourced ``v_inf`` node (Russell's per-leg recipe — it re-anchors v∞ at every
node), exactly as S1L1 does. The Mars-arrival miss / v∞ are EVIDENCE measured against
Russell's printed transit-characteristics table (the EXPECTED side, sourced, never
fit). The ``dv`` mid-leg burn IS applied (these parents are not ballistic), so the
reconstructed arrival reflects the published powered solution.

Provenance / honesty
--------------------
Every number in an :class:`AppCBlock` traces verbatim to Russell 2004 Appendix C
(transcribed from the held dissertation text layer; the per-leg "DATA NECESSARY TO
REPRODUCE" block and the "EARTH TO MARS TRANSIT LEG CHARACTERISTICS" table). The
achieved miss / v∞ are integrator output. The real-eph Mars v∞ BREATHES across the
seven cycles (it is not a single coplanar anchor) — any gate targets the App-C
per-leg values.

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

# (leg_no, body, time_start_days, vinf_kms | None, time_dv_days | None, dv_kms | None)
AppCLeg = tuple[
    int,
    str,
    float,
    tuple[float, float, float] | None,
    float | None,
    tuple[float, float, float] | None,
]


@dataclass(frozen=True)
class AppCBlock:
    """One Russell App-C parent-cycler reproduction block (verbatim transcription)."""

    shorthand: str
    parent_number: int
    catalogue_id: str
    epoch_days: float
    total_dv_kms: float  # printed "Total delta v over NN yr" — sourced, the published maintenance
    legs: tuple[AppCLeg, ...]
    # ARRIVAL-leg-no -> (published transit days, published Mars v_inf km/s)
    mars_transit: dict[int, tuple[float, float]]
    mars_vinf_avg: float  # App-C AVERAGE row Mars v_inf
    earth_vinf_avg: float  # App-C AVERAGE row Earth v_inf


@dataclass(frozen=True)
class AppCArc:
    """One reconstructed App-C leg arc (seeded from a sourced v_inf node)."""

    leg_no: int
    body: str
    next_body: str
    t0_sec: float
    t1_sec: float
    r0_km: Vec3
    v0_km_s: Vec3
    t_dv_sec: float | None  # epoch of the mid-leg Δv (None if ~ballistic / no dv)
    dv_km_s: Vec3 | None
    is_mars_transit: bool  # True for an E->M outbound (Mars-arrival) leg


@dataclass(frozen=True)
class AppCMarsEncounter:
    """Achieved Mars encounter of a reconstructed Mars-transit leg."""

    arrival_leg_no: int
    cycle: int
    t_arr_sec: float
    miss_au: float
    vinf_kms: float
    sc_lon_deg: float
    mars_lon_deg: float
    pub_transit_days: float
    pub_vinf_kms: float


def _vinf_vec(vinf: tuple[float, float, float] | None) -> Vec3:
    assert vinf is not None
    return np.array(vinf, dtype=np.float64)


def _lon_deg(r: Vec3) -> float:
    return float(np.degrees(np.arctan2(r[1], r[0])) % 360.0)


def build_seeded_arcs(block: AppCBlock, ephem: Ephemeris) -> list[AppCArc]:
    """Reconstruct every leg arc of an App-C block from its sourced v_inf nodes.

    Russell's recipe: ``v_sc = v_planet(ephem) + v_inf`` at each leg start, placing
    the departure node at the REAL planet position (the longitude-rendezvous anchor).
    The mid-leg ``dv`` (if non-negligible) is recorded for the encounter propagation.
    """
    arcs: list[AppCArc] = []
    legs = block.legs
    for i in range(len(legs) - 1):
        leg_no, body, ts, vinf, t_dv, dv = legs[i]
        _next_no, next_body, next_ts, _v, _td, _dv = legs[i + 1]
        if vinf is None:
            continue
        t0 = (block.epoch_days + ts) * DAY_S
        t1 = (block.epoch_days + next_ts) * DAY_S
        r_p, v_p = ephem.state(body, t0)
        r0 = np.asarray(r_p, dtype=np.float64)
        v0 = np.asarray(v_p, dtype=np.float64) + _vinf_vec(vinf)
        t_dv_sec = (block.epoch_days + t_dv) * DAY_S if t_dv is not None else None
        dv_vec = np.array(dv, dtype=np.float64) if dv is not None else None
        arcs.append(
            AppCArc(
                leg_no=leg_no,
                body=body,
                next_body=next_body,
                t0_sec=t0,
                t1_sec=t1,
                r0_km=r0,
                v0_km_s=v0,
                t_dv_sec=t_dv_sec,
                dv_km_s=dv_vec,
                is_mars_transit=body == "E" and next_body == "M",
            )
        )
    return arcs


def _propagate_with_dv(arc: AppCArc) -> tuple[Vec3, Vec3]:
    """Two-body-Sun Kepler propagation of one arc to arrival, applying any mid-leg Δv."""
    if arc.t_dv_sec is not None and arc.dv_km_s is not None and arc.t_dv_sec < arc.t1_sec:
        r_dv, v_dv = kepler_propagate(arc.r0_km, arc.v0_km_s, arc.t_dv_sec - arc.t0_sec)
        v_dv = v_dv + arc.dv_km_s
        return kepler_propagate(r_dv, v_dv, arc.t1_sec - arc.t_dv_sec)
    return kepler_propagate(arc.r0_km, arc.v0_km_s, arc.t1_sec - arc.t0_sec)


def reconstruct_mars_encounters(block: AppCBlock, ephem: Ephemeris) -> list[AppCMarsEncounter]:
    """Two-body-Sun reconstruction of every Mars-transit leg's encounter (with Δv)."""
    arcs = build_seeded_arcs(block, ephem)
    out: list[AppCMarsEncounter] = []
    mars_arrivals = sorted(block.mars_transit)
    for arc in arcs:
        if not arc.is_mars_transit:
            continue
        arrival_no = arc.leg_no + 1
        if arrival_no not in block.mars_transit:
            continue
        cycle = mars_arrivals.index(arrival_no) + 1
        rk, vk = _propagate_with_dv(arc)
        r_m, v_m = ephem.state("M", arc.t1_sec)
        r_m = np.asarray(r_m, dtype=np.float64)
        v_m = np.asarray(v_m, dtype=np.float64)
        miss_au = float(np.linalg.norm(rk - r_m) / AU_KM)
        vinf_kms = float(np.linalg.norm(vk - v_m))
        pub_tof, pub_vinf = block.mars_transit[arrival_no]
        out.append(
            AppCMarsEncounter(
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


@dataclass(frozen=True)
class ContinuousNode:
    """One node of a SINGLE continuous trajectory (the #169 horizon-TCM measure).

    Mirrors :class:`cyclerfinder.search.s1l1_corrected.ContinuousNode` but generic
    over an :class:`AppCBlock`. Propagates ONE trajectory forward from the start node
    with NO v_inf re-anchoring between legs; at each node a free flyby rotates the
    arrived v_inf toward the next leg's App-C direction and the MAINTENANCE dv is the
    part a ballistic flyby cannot supply (|v_inf| magnitude change + un-bendable
    bend shortfall). This is the maintenance the per-leg re-anchoring hides.
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
    block: AppCBlock,
    ephem: Ephemeris,
    *,
    perturbers: tuple[str, ...] = (),
    propagate: object | None = None,
    start_leg_no: int | None = None,
) -> list[ContinuousNode]:
    """Propagate ONE continuous trajectory through an App-C block's nodes (#170 TCM).

    Generic version of :func:`s1l1_corrected.continuous_chain`. Starts from the first
    Earth-departure that begins a Mars-transit (or ``start_leg_no``), walks forward
    node-to-node WITHOUT re-anchoring v_inf. ``propagate is None`` => two-body-Sun
    Kepler; otherwise a supplied n-body propagator. At each node the flyby patch snaps
    the position to the real planet and rotates the velocity to the next App-C v_inf
    direction; the maintenance dv is measured. EXPECTED = App-C per-leg v_inf
    (sourced); miss / v_inf / dv are EVIDENCE.
    """
    from cyclerfinder.core.flyby import max_bend

    legs = block.legs
    nodes = [
        ((block.epoch_days + ts) * DAY_S, body, vinf) for (_no, body, ts, vinf, _td, _dv) in legs
    ]
    if start_leg_no is None:
        # first Earth-departure leg that goes to Mars
        start_leg_no = next(
            legs[i][0] for i in range(len(legs) - 1) if legs[i][1] == "E" and legs[i + 1][1] == "M"
        )
    start_idx = next(i for i, leg in enumerate(legs) if leg[0] == start_leg_no)

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
            r = r_pl.copy()
            v = v_pl + vinf_out_vec
        else:
            r = rk
            v = vk

        out.append(
            ContinuousNode(
                arrival_leg_no=legs[i][0],
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


# =====================================================================================
# Verbatim App-C transcriptions — the reachable #170 parents (sourced parent number).
# Russell 2004 Appendix C, "DATA NECESSARY TO REPRODUCE CYCLER" + transit table.
# Model DE405, ecliptic via obliquity 0.409092629205 rad (== our backend to 1e-5 deg);
# v_inf vectors already ecliptic. Each leg: (no, E/M, time_start_d, vinf, time_dv_d, dv).
# dv recorded only where the printed magnitude is a real burn (>~1e-6 km/s); the
# ~1e-9..1e-11 "numerical zero" dv components are recorded as None (ballistic node).
# =====================================================================================


def _dv_or_none(
    dv: tuple[float, float, float], thresh: float = 1e-6
) -> tuple[float, float, float] | None:
    return dv if float(np.linalg.norm(dv)) > thresh else None


# --- 8.049gGf2 (parent #188), total Δv 0.436091 km/s over 29.95 yr ---
# russell-ch4-8.049gGf2. Mars arrivals: legs 3, 7, 11, 15, 19, 23, 27.
_B188_EPOCH = 15567.931545
APPC_8049GGF2 = AppCBlock(
    shorthand="8.049gGf2",
    parent_number=188,
    catalogue_id="russell-ch4-8.049gGf2",
    epoch_days=_B188_EPOCH,
    total_dv_kms=0.436091,
    legs=(
        (
            1,
            "E",
            -0.215814313e2,
            (-0.686220461e1, 0.943299852e1, -0.133901350e-2),
            0.247738521e3,
            _dv_or_none((-0.314952839e-10, 0.508415129e-10, 0.113031458e-10)),
        ),
        (
            2,
            "E",
            0.539501803e3,
            (-0.977613095e1, 0.637258366e1, 0.593332504e0),
            0.553930798e3,
            _dv_or_none((0.191703091e-9, 0.314400609e-9, 0.258981706e-10)),
        ),
        (
            3,
            "M",
            0.635695107e3,
            (-0.110669417e2, -0.192220229e1, -0.510376717e0),
            0.726121716e3,
            _dv_or_none((-0.137230750e-11, -0.327766553e-10, -0.209752242e-10)),
        ),
        (
            4,
            "E",
            0.116761634e4,
            (-0.640028065e1, -0.587067019e1, -0.271176934e0),
            0.140137157e4,
            _dv_or_none((-0.169094176e-7, 0.553258335e-8, 0.346887151e-5)),
        ),
        (
            5,
            "E",
            0.153285889e4,
            (-0.716169990e1, -0.493035001e1, -0.261389264e-2),
            0.177429683e4,
            _dv_or_none((-0.325418750e-10, -0.607452743e-11, -0.733573053e-12)),
        ),
        (
            6,
            "E",
            0.208158148e4,
            (-0.442527052e1, -0.742368055e1, -0.901218506e0),
            0.209309238e4,
            _dv_or_none((-0.154419315e-9, -0.295062647e-10, 0.379244618e-12)),
        ),
        (
            7,
            "M",
            0.215832084e4,
            (-0.294368226e0, -0.116847671e2, 0.617971791e0),
            0.224696412e4,
            _dv_or_none((-0.508348994e-10, 0.183095582e-10, 0.221232917e-11)),
        ),
        (
            8,
            "E",
            0.274927606e4,
            (0.527829816e1, -0.178596893e1, 0.246550763e1),
            0.300131850e4,
            _dv_or_none((0.854999355e-5, -0.361705396e-4, -0.837469343e-3)),
        ),
        (
            9,
            "E",
            0.311455495e4,
            (0.540096030e1, -0.283063892e1, 0.148055420e-2),
            0.338901420e4,
            _dv_or_none((0.103733220e-10, -0.595548388e-11, 0.194776334e-11)),
        ),
        (
            10,
            "E",
            0.365271034e4,
            (0.447580351e1, 0.412847266e1, 0.407652954e-1),
            0.366802570e4,
            _dv_or_none((0.133871077e-9, -0.115650801e-9, -0.471398736e-11)),
        ),
        (
            11,
            "M",
            0.375481270e4,
            (0.664118892e1, 0.527066640e1, -0.858435743e0),
            0.383815060e4,
            _dv_or_none((0.620349126e-11, -0.385304265e-10, -0.632201812e-11)),
        ),
        (
            12,
            "E",
            0.431039870e4,
            (0.120341614e1, 0.935217422e1, 0.119325811e1),
            0.454416535e4,
            _dv_or_none((0.374633848e-5, 0.168821458e-5, -0.187880800e-3)),
        ),
        (
            13,
            "E",
            0.467565909e4,
            (0.258295956e1, 0.913055053e1, 0.377177379e-2),
            0.491855809e4,
            _dv_or_none((-0.377434165e-10, -0.557057643e-10, 0.855138246e-12)),
        ),
        (
            14,
            "E",
            0.522770227e4,
            (-0.199822502e1, 0.922154981e1, 0.112556115e1),
            0.524391356e4,
            _dv_or_none((0.386560836e-10, -0.571319721e-9, 0.414833677e-11)),
        ),
        (
            15,
            "M",
            0.533577755e4,
            (-0.719639155e1, 0.653951981e1, 0.261396054e1),
            0.561687807e4,
            _dv_or_none((0.352545758e-2, 0.260356785e-1, 0.342452071e0)),
        ),
        (
            16,
            "E",
            0.584686941e4,
            (-0.102058981e2, 0.425279807e1, 0.742192448e0),
            0.606967870e4,
            _dv_or_none((0.248194632e-6, -0.183087497e-5, 0.133166221e-3)),
        ),
        (
            17,
            "E",
            0.621213054e4,
            (-0.949165378e1, 0.569381648e1, -0.111286407e-2),
            0.644675101e4,
            _dv_or_none((0.455080315e-9, -0.552752235e-9, -0.188194168e-11)),
        ),
        (
            18,
            "E",
            0.677075071e4,
            (-0.109020320e2, 0.201198205e1, 0.196851251e0),
            0.678445900e4,
            _dv_or_none((-0.545979980e-9, -0.431994135e-8, -0.346580057e-9)),
        ),
        (
            19,
            "M",
            0.686213934e4,
            (-0.983689346e1, -0.613987917e1, 0.612910911e0),
            0.698918729e4,
            _dv_or_none((0.361098745e-11, 0.902085510e-11, -0.407589795e-11)),
        ),
        (
            20,
            "E",
            0.741452172e4,
            (-0.140224597e1, -0.691615025e1, -0.651351913e0),
            0.765925617e4,
            _dv_or_none((0.126995358e-6, 0.459286355e-7, -0.116222772e-4)),
        ),
        (
            21,
            "E",
            0.777979701e4,
            (-0.240009873e1, -0.666911450e1, 0.909572774e-3),
            0.800751231e4,
            _dv_or_none((0.117231435e-9, 0.101195278e-9, 0.798102240e-11)),
        ),
        (
            22,
            "E",
            0.832197628e4,
            (0.193311293e1, -0.670258718e1, -0.121590198e1),
            0.833330383e4,
            _dv_or_none((-0.554736294e-9, -0.174608769e-9, -0.306892599e-10)),
        ),
        (
            23,
            "M",
            0.839749327e4,
            (0.568618079e1, -0.913833239e1, 0.198941193e0),
            0.848728773e4,
            _dv_or_none((-0.705663665e-11, 0.249660762e-11, -0.489153222e-11)),
        ),
        (
            24,
            "E",
            0.899612303e4,
            (0.610880025e1, 0.199544212e1, -0.335866964e0),
            0.924448848e4,
            _dv_or_none((-0.582932951e-6, 0.466007410e-6, -0.124238527e-3)),
        ),
        (
            25,
            "E",
            0.936136634e4,
            (0.633042698e1, 0.111749571e1, -0.184004786e-2),
            0.962033609e4,
            _dv_or_none((0.499751404e-10, -0.354402546e-11, -0.136285712e-11)),
        ),
        (
            26,
            "E",
            0.990088665e4,
            (0.276463368e1, 0.576484300e1, 0.761426327e0),
            0.991776828e4,
            _dv_or_none((0.793466544e-10, 0.142521574e-10, -0.299324240e-11)),
        ),
        (
            27,
            "M",
            0.100134309e5,
            (0.174125530e1, 0.792187044e1, -0.227966635e1),
            0.101966079e5,
            _dv_or_none((0.504281672e-11, -0.710293136e-11, -0.293888745e-11)),
        ),
        (
            28,
            "E",
            0.105521869e5,
            (-0.385930219e1, 0.834630649e1, 0.326150615e1),
            0.107859567e5,
            _dv_or_none((-0.291378851e-5, -0.429281669e-5, 0.871638207e-4)),
        ),
        (29, "E", 0.109174522e5, None, None, None),
    ),
    mars_transit={
        3: (96.2, 11.244),
        7: (76.7, 11.705),
        11: (102.1, 8.522),
        15: (108.1, 10.069),
        19: (91.4, 11.612),
        23: (75.5, 10.765),
        27: (112.5, 8.425),
    },
    mars_vinf_avg=10.335,
    earth_vinf_avg=8.653,
)


# --- 8.165Gfh-f2 (parent #192), total Δv 1.677496 km/s over 30.36 yr ---
# russell-ch4-8.165Gfh-f2. Mars arrivals: legs 2, 7, 12, 17, 22, 27, 32.
_B192_EPOCH = 16114.706442
APPC_8165GFHF2 = AppCBlock(
    shorthand="8.165Gfh-f2",
    parent_number=192,
    catalogue_id="russell-ch4-8.165Gfh-f2",
    epoch_days=_B192_EPOCH,
    total_dv_kms=1.677496,
    legs=(
        (
            1,
            "E",
            -0.118341070e3,
            (-0.776652940e1, 0.104360865e1, 0.234580458e1),
            -0.915653486e2,
            _dv_or_none((-0.368386855e-9, -0.414311209e-9, -0.196223439e-10)),
        ),
        (
            2,
            "M",
            0.601637364e2,
            (-0.105977858e2, 0.361455507e0, -0.249654378e0),
            0.422470588e3,
            _dv_or_none((0.248575481e-3, -0.166932604e-3, 0.793744065e-2)),
        ),
        (
            3,
            "E",
            0.617558893e3,
            (-0.523317048e1, -0.428226388e1, -0.362177109e1),
            0.862272817e3,
            _dv_or_none((0.426024218e-4, -0.158347365e-4, -0.703260362e-3)),
        ),
        (
            4,
            "E",
            0.982803556e3,
            (0.486605302e0, -0.863697137e0, -0.762202740e1),
            0.100964107e4,
            _dv_or_none((-0.473680869e-9, 0.101938639e-8, -0.241711299e-9)),
        ),
        (
            5,
            "E",
            0.116172031e4,
            (-0.584834758e1, -0.213454766e1, 0.437989883e1),
            0.128225622e4,
            _dv_or_none((-0.307913680e-4, 0.857852898e-4, 0.245676729e-3)),
        ),
        (
            6,
            "E",
            0.152698063e4,
            (-0.358320783e1, -0.667744936e1, -0.849147109e0),
            0.153918699e4,
            _dv_or_none((0.445089255e-8, -0.211834313e-7, 0.527855911e-9)),
        ),
        (
            7,
            "M",
            0.160835635e4,
            (-0.702656511e0, -0.118152358e2, 0.680191080e0),
            0.169662079e4,
            _dv_or_none((0.318970558e-10, -0.144654438e-10, 0.111603321e-11)),
        ),
        (
            8,
            "E",
            0.219678597e4,
            (0.388150757e1, -0.152027835e1, -0.483733297e1),
            0.239402190e4,
            _dv_or_none((0.111616730e-2, -0.485230585e-3, 0.272941312e-3)),
        ),
        (
            9,
            "E",
            0.256203770e4,
            (0.326771043e0, 0.586216647e0, -0.633912387e1),
            0.258985618e4,
            _dv_or_none((-0.109000955e-7, -0.171172809e-7, -0.435492444e-8)),
        ),
        (
            10,
            "E",
            0.274749425e4,
            (0.488325769e1, -0.350324377e1, 0.163170248e1),
            0.286437963e4,
            _dv_or_none((-0.616969723e-5, -0.856922449e-5, 0.108924551e-3)),
        ),
        (
            11,
            "E",
            0.311276105e4,
            (0.491306178e1, 0.383178059e1, 0.843339928e-1),
            0.312761920e4,
            _dv_or_none((-0.424420684e-7, 0.385607033e-7, -0.361801309e-8)),
        ),
        (
            12,
            "M",
            0.321181536e4,
            (0.633697974e1, 0.543153340e1, -0.463629268e0),
            0.355582567e4,
            _dv_or_none((0.100165795e-2, 0.109132921e-3, 0.129168049e0)),
        ),
        (
            13,
            "E",
            0.376667070e4,
            (0.253221474e0, 0.708980874e1, -0.582734839e1),
            0.397852027e4,
            _dv_or_none((-0.101385847e-3, 0.623486707e-3, 0.140859196e-4)),
        ),
        (
            14,
            "E",
            0.413192857e4,
            (-0.138645922e1, 0.323812987e0, -0.906772276e1),
            0.415958333e4,
            _dv_or_none((0.159349702e-9, 0.985435830e-9, -0.786216462e-10)),
        ),
        (
            15,
            "E",
            0.431629360e4,
            (0.299990121e1, 0.644120252e1, 0.624389354e1),
            0.446604913e4,
            _dv_or_none((0.129088877e-2, -0.644338573e-4, -0.182838676e-3)),
        ),
        (
            16,
            "E",
            0.468155099e4,
            (-0.176404172e1, 0.924029643e1, 0.110275877e1),
            0.469824385e4,
            _dv_or_none((0.234353416e-7, -0.109730748e-6, -0.419092378e-9)),
        ),
        (
            17,
            "M",
            0.479283672e4,
            (-0.737412512e1, 0.599758478e1, 0.984229653e0),
            0.508620660e4,
            _dv_or_none((0.214869079e-2, 0.817308165e-2, 0.514513444e0)),
        ),
        (
            18,
            "E",
            0.530752072e4,
            (-0.690675347e1, 0.175047953e1, -0.608787555e1),
            0.549745662e4,
            _dv_or_none((0.713764026e-4, 0.751925203e-4, 0.847517346e-4)),
        ),
        (
            19,
            "E",
            0.567278208e4,
            (-0.683369119e0, -0.131973701e1, -0.926457865e1),
            0.575722568e4,
            _dv_or_none((-0.435264916e-6, 0.251917616e-6, 0.363727897e-5)),
        ),
        (
            20,
            "E",
            0.585244932e4,
            (-0.548136907e1, 0.460268341e1, 0.633708866e1),
            0.599124626e4,
            _dv_or_none((0.559253754e-3, 0.119821051e-2, -0.121839659e-3)),
        ),
        (
            21,
            "E",
            0.621770442e4,
            (-0.968537660e1, 0.169660829e1, 0.216873051e0),
            0.623200621e4,
            _dv_or_none((-0.193007848e-6, 0.173599081e-7, 0.601720793e-8)),
        ),
        (
            22,
            "M",
            0.631304973e4,
            (-0.100187302e2, -0.595162193e1, 0.573052933e0),
            0.641225116e4,
            _dv_or_none((-0.959656290e-10, 0.182850881e-9, -0.254171397e-11)),
        ),
        (
            23,
            "E",
            0.686416880e4,
            (-0.109267726e1, -0.514154702e1, -0.511950094e1),
            0.710889837e4,
            _dv_or_none((-0.877215913e-5, -0.630908647e-3, 0.106915796e-2)),
        ),
        (
            24,
            "E",
            0.722943682e4,
            (0.820045321e0, -0.329009544e0, -0.727422106e1),
            0.725649324e4,
            _dv_or_none((-0.919804149e-8, 0.102894872e-8, -0.227598274e-8)),
        ),
        (
            25,
            "E",
            0.740981297e4,
            (-0.294474881e1, -0.475821023e1, 0.442060633e1),
            0.752669655e4,
            _dv_or_none((-0.312221031e-4, 0.163839136e-4, 0.221006722e-3)),
        ),
        (
            26,
            "E",
            0.777507418e4,
            (0.197467747e1, -0.673431904e1, -0.121948337e1),
            0.778633563e4,
            _dv_or_none((-0.121673942e-7, -0.927388876e-9, -0.996284913e-9)),
        ),
        (
            27,
            "M",
            0.785015052e4,
            (0.566090058e1, -0.928037069e1, 0.251557026e0),
            0.793953229e4,
            _dv_or_none((-0.198405904e-12, -0.288349914e-11, 0.556838570e-11)),
        ),
        (
            28,
            "E",
            0.844602899e4,
            (0.430964624e1, 0.140038597e1, -0.489399160e1),
            0.868710077e4,
            _dv_or_none((-0.756865752e-3, -0.309597639e-3, -0.850720436e-4)),
        ),
        (
            29,
            "E",
            0.881128926e4,
            (-0.117736144e0, 0.747899649e0, -0.663404891e1),
            0.886164404e4,
            _dv_or_none((-0.391271455e-7, 0.196666049e-7, 0.588314013e-8)),
        ),
        (
            30,
            "E",
            0.899778842e4,
            (0.563010684e1, 0.295576619e0, 0.355243776e1),
            0.911101522e4,
            _dv_or_none((-0.618254890e-5, 0.767596609e-4, -0.577820269e-3)),
        ),
        (
            31,
            "E",
            0.936303618e4,
            (0.390936035e1, 0.530666105e1, 0.851636561e0),
            0.938123708e4,
            _dv_or_none((0.476513758e-7, 0.444872464e-7, -0.532394549e-9)),
        ),
        (
            32,
            "M",
            0.948437553e4,
            (-0.459285568e-1, 0.676530292e1, -0.627611643e0),
            0.959995466e4,
            _dv_or_none((-0.173357159e-10, -0.148820962e-10, 0.588707518e-12)),
        ),
        (
            33,
            "E",
            0.100622712e5,
            (-0.504453923e1, 0.134946798e1, -0.230144516e0),
            0.103106583e5,
            _dv_or_none((-0.261243213e-6, 0.821613821e-6, 0.201531072e-3)),
        ),
        (
            34,
            "E",
            0.104275464e5,
            (-0.170260265e0, -0.425378002e0, -0.522569766e1),
            0.104598411e5,
            _dv_or_none((-0.391091104e-10, -0.282232398e-10, 0.573715472e-10)),
        ),
        (
            35,
            "E",
            0.106069617e5,
            (-0.463455451e1, 0.234814522e1, 0.125165921e1),
            0.107201927e5,
            _dv_or_none((0.110210816e-6, 0.229946431e-6, 0.582576560e-5)),
        ),
        (36, "E", 0.109722230e5, None, None, None),
    ),
    mars_transit={
        2: (178.5, 10.607),
        7: (81.4, 11.856),
        12: (99.1, 8.359),
        17: (111.3, 9.556),
        22: (95.3, 11.667),
        27: (75.1, 10.874),
        32: (121.3, 6.795),
    },
    mars_vinf_avg=9.959,
    earth_vinf_avg=7.779,
)


REACHABLE_BLOCKS: dict[str, AppCBlock] = {
    APPC_8049GGF2.catalogue_id: APPC_8049GGF2,
    APPC_8165GFHF2.catalogue_id: APPC_8165GFHF2,
}


__all__ = [
    "APPC_8049GGF2",
    "APPC_8165GFHF2",
    "REACHABLE_BLOCKS",
    "AppCArc",
    "AppCBlock",
    "AppCLeg",
    "AppCMarsEncounter",
    "ContinuousNode",
    "build_seeded_arcs",
    "continuous_chain",
    "reconstruct_mars_encounters",
]
