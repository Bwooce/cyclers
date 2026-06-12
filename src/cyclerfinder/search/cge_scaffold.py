"""Liang et al. 2024 CGE triple-cycler reproduction scaffold (idealized model).

Reproduces the three IDEALIZED Callisto-Ganymede-Europa triple-cycler members
of Liang, Yang, Li, Bai & Qin, "Callisto-Ganymede-Europa Triple Cyclers",
*Journal of Guidance, Control, and Dynamics* (Engineering Note), 2024,
DOI 10.2514/1.G008387 — Members A (Tables 2-3), B (Tables 4-5) and
C (Tables 6-7) — in THEIR model, from THEIR printed inputs:

* circular, coplanar moon orbits about Jupiter at the Table 1 mean motions
  (Europa 1.7693, Ganymede 0.8782, Callisto 0.3765 rad/day); orbit radii
  DERIVED from mu_Jupiter + mean motion (the paper never prints mu_Jupiter
  nor the radii — the one reproduction gap; we use the in-repo registry
  Jupiter GM, :data:`cyclerfinder.core.satellites.PRIMARIES`);
* conic (Jupiter two-body) arcs between flybys, solved as fixed-time Lambert
  problems at the PRINTED per-leg times of flight (Tables 3/5/7);
* instantaneous zero-radius gravity assists (V-infinity magnitude preserved,
  direction turned) — the reproduction CHECKS the magnitude-continuity that
  a ballistic flyby requires rather than enforcing it.

Geometry conventions (Sec. III.A of the paper, validated numerically here):

* The initial spacecraft orbit has period T = T_Callisto, perijove radius
  r_p (member-specific) and departs perijove on the perifocal +x axis at
  t = 0; the Tables 2/4/6 moon phases are the moons' polar angles at t = 0.
  Because T = T_Callisto, the spacecraft semi-major axis EQUALS Callisto's
  orbit radius (Eq. 14), so the conic crosses Callisto's circle exactly at
  eccentric anomaly E = +/- pi/2; the structure flag 1 (flyby above the
  x-axis, Fig. 2) selects E = +pi/2.
* Per the Eq. 16 epoch law with n_cycle = 1, the first tabulated Callisto
  flyby (Tables 3/5/7 row 0, transit time 0) is t_c1 = T + t_c0: ONE FULL
  SPACECRAFT REVOLUTION after perijove departure, where
  t_c0 = (pi/2 - e) / n_Callisto is the perijove->crossing time. The
  encounter state is identical at t_c0 and t_c0 + T (both spacecraft and
  Callisto are exactly periodic with period T), but Ganymede's and Europa's
  phases at the downstream flyby epochs are not — anchoring at t_c0 instead
  of t_c0 + T mis-places them and the reproduction fails by ~1 km/s.
* All four transfer legs of every member resolve to ONE-revolution Lambert
  solutions (n_revs = 1), matching the paper's "all transfer arcs are forced
  to more than one revolution" (p. 9).

In-print defects honoured from data/errata.yaml (never imported):

* ``liang-2024-eq16-synodic-subscript``: Eq. 16's t_e1 prints S_G,E where the
  construction requires S_C,E. (Only the t_c1 term of Eq. 16 is used here.)
* ``liang-2024-eq13-denominator``: Eq. 13's gamma denominator prints
  nu_2^2 + nu_2^2 for nu_2^2 + nu_3^2. (Eq. 13 is not needed here at all:
  the reproduction is ballistic-by-table, no defect-Delta-v optimisation.)

The "Flyby Altitude" columns of Tables 3/5/7 are Delta-v-equivalent fictions
computed from the moon-centred hyperbola IGNORING Jupiter's gravity (the
paper's own caveat, p. 16). :func:`defect_altitude_km` reproduces exactly
that convention — informational only, never a golden anchor.

Golden tests: ``tests/search/test_liang_cge_reproduction.py``.
Mining note: ``docs/notes/2026-06-11-liang-2024-cge-triple-cyclers-mining.md``.
Results note: ``docs/notes/2026-06-13-liang-abc-reproduction.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64 — matches lambert.py

DAY_S: float = 86400.0

MU_JUPITER_KM3_S2: float = PRIMARIES["Jupiter"]
"""Jupiter GM used for the reproduction (km^3/s^2).

The paper never prints its mu_Jupiter (mining note section 9 item 4). We use
the in-repo registry value (JPL SSD gm_de440, satellites.py:49). Sensitivity
(2026-06-13 scan, results note): replacing it with any plausible alternative
(1.266e8..1.26712764e8) moves the worst V-infinity residual by ~1e-3 km/s —
an order of magnitude below the Table 1 print-quantization floor, so the
choice is immaterial at the validation tolerance.
"""

LIANG_MEAN_MOTIONS_RAD_DAY: dict[str, float] = {
    # SOURCED (expected-side input): Liang et al. 2024 Table 1 (p. 4).
    "Europa": 1.7693,
    "Ganymede": 0.8782,
    "Callisto": 0.3765,
}

MEAN_MOTION_HALF_ULP_RAD_DAY: float = 5.0e-5
"""Half-ULP of the 4-decimal Table 1 mean motions — the dominant input
quantization. See :func:`vinf_print_tolerance_kms` for how it propagates."""

CGCEC_SEQUENCE: tuple[str, ...] = ("Callisto", "Ganymede", "Callisto", "Europa", "Callisto")


@dataclass(frozen=True)
class CGEMemberSpec:
    """Printed inputs + printed anchors for one idealized member.

    All numeric fields are TRANSCRIBED from the paper (table cited per
    member below); none is a value our own code computed.
    """

    member: str  # "A" | "B" | "C"
    structure: str  # 0/1 flags for (Callisto, Ganymede, Europa), Fig. 2 convention
    phases_rad: dict[str, float]  # initial moon polar angles at t = 0 (perijove departure)
    perijove_scale: float  # r_p = perijove_scale * (r_Europa - 10000 km)
    tofs_days: tuple[float, float, float, float]  # printed per-leg ToF (C-G, G-C, C-E, E-C)
    vinf_printed_kms: tuple[float, float, float, float, float]  # printed per-flyby V-infinity
    altitudes_printed_km: tuple[float, float, float, float, float]  # printed fictions (p. 16)


LIANG_MEMBERS: dict[str, CGEMemberSpec] = {
    # Member A — 1-1-1 structure, high perijove. Tables 2 (phases, p. 14)
    # and 3 (per-flyby ToF / V_inf / "altitude", p. 15).
    "A": CGEMemberSpec(
        member="A",
        structure="1-1-1",
        phases_rad={"Europa": 0.063748, "Ganymede": 0.70579, "Callisto": 1.3550},
        perijove_scale=1.0,
        tofs_days=(31.8973, 18.1697, 29.9343, 19.9747),
        vinf_printed_kms=(5.6730, 6.9919, 5.6698, 4.6685, 5.8721),
        altitudes_printed_km=(1900851.0, 978172.0, 33839.0, 6241.0, 19765.0),
    ),
    # Member B — 1-1-0 structure, high perijove. Tables 4 (phases, p. 15)
    # and 5 (p. 17). First two flybys identical to Member A by construction.
    "B": CGEMemberSpec(
        member="B",
        structure="1-1-0",
        phases_rad={"Europa": -0.48843, "Ganymede": 0.70579, "Callisto": 1.3550},
        perijove_scale=1.0,
        tofs_days=(31.8973, 18.1697, 30.2850, 19.6834),
        vinf_printed_kms=(5.6730, 6.9919, 5.6698, 4.4853, 5.7914),
        altitudes_printed_km=(1900851.0, 978172.0, 60877.0, 10825.0, 36258.0),
    ),
    # Member C — 1-1-1 structure, LOW perijove (r_p halved). Tables 6
    # (phases, p. 18) and 7 (p. 18).
    "C": CGEMemberSpec(
        member="C",
        structure="1-1-1",
        phases_rad={"Europa": 0.94281, "Ganymede": 1.3883, "Callisto": 1.7936},
        perijove_scale=0.5,
        tofs_days=(32.2542, 17.8127, 31.1267, 18.9362),
        vinf_printed_kms=(7.6433, 10.4922, 7.6409, 12.0213, 7.7838),
        altitudes_printed_km=(1329021.0, 629226.0, 27438.0, 3636.0, 12516.0),
    ),
}


def derived_radius_km(moon: str, *, mu_jupiter: float = MU_JUPITER_KM3_S2) -> float:
    """Circular orbit radius (km) from the Table 1 mean motion + Jupiter GM.

    Kepler III: r = (mu / n^2)^(1/3) with n in rad/s. The paper's only
    printed radius anchor is "r_Eu - 10000, which is about 660988 km"
    (p. 14); this derivation gives 660993.5 km with the registry GM —
    5.5 km off, inside the paper's "about".
    """
    n_rad_s = LIANG_MEAN_MOTIONS_RAD_DAY[moon] / DAY_S
    return float((mu_jupiter / n_rad_s**2) ** (1.0 / 3.0))


@dataclass(frozen=True)
class LegRecord:
    """One solved transfer leg (fixed-time Lambert at the printed ToF)."""

    index: int  # 0..3 (C-G, G-C, C-E, E-C)
    n_revs: int
    branch: str
    tof_days: float  # printed input
    vinf_dep_kms: float  # |v1 - v_moon| at the departure flyby
    vinf_arr_kms: float  # |v2 - v_moon| at the arrival flyby
    selection_margin_kms: float  # cost gap to the second-best Lambert solution


@dataclass(frozen=True)
class FlybyRecord:
    """Residual row for one flyby of the first cycle."""

    index: int  # 0..4
    moon: str
    epoch_days: float  # days from perijove departure (t = 0)
    vinf_printed_kms: float
    vinf_in_kms: float  # inbound |V_inf| (initial conic for flyby 0)
    vinf_out_kms: float | None  # outbound |V_inf| (None at flyby 4: cycle end)
    residual_in_kms: float  # |vinf_in - printed|
    residual_out_kms: float | None
    continuity_kms: float | None  # |vinf_in - vinf_out| (ballistic => 0)
    defect_altitude_km: float | None  # paper-convention fiction, informational


@dataclass(frozen=True)
class MemberReproduction:
    """Full reproduction result for one member."""

    member: str
    mu_jupiter: float
    radii_km: dict[str, float]
    perijove_km: float
    eccentricity: float
    t_c0_days: float  # perijove -> Callisto-circle crossing (E = +pi/2)
    crossing_longitude_rad: float  # spacecraft polar angle at the crossing
    callisto_angle_residual_rad: float  # |theta_Callisto(t_c0) - crossing| (phase check)
    conic_vinf_kms: float  # |V_inf| of the initial conic at the crossing
    cycle_tof_days: float  # sum of the printed per-leg ToFs
    legs: tuple[LegRecord, ...]
    flybys: tuple[FlybyRecord, ...]
    max_vinf_residual_kms: float  # worst |V_inf - printed| over all in/out sides


def _moon_state(
    moon: str,
    t_days: float,
    phases_rad: dict[str, float],
    radii_km: dict[str, float],
    mu_jupiter: float,
) -> tuple[Vec3, Vec3]:
    """Circular-coplanar moon position/velocity at t (days from perijove)."""
    theta = phases_rad[moon] + LIANG_MEAN_MOTIONS_RAD_DAY[moon] * t_days
    r = radii_km[moon]
    v_circ = math.sqrt(mu_jupiter / r)
    pos = np.array([r * math.cos(theta), r * math.sin(theta), 0.0])
    vel = np.array([-v_circ * math.sin(theta), v_circ * math.cos(theta), 0.0])
    return pos, vel


def defect_altitude_km(vinf_in: Vec3, vinf_out: Vec3, moon: str) -> float:
    """Flyby "altitude" in the paper's own no-Jupiter hyperbola convention.

    Tables 3/5/7 print altitudes computed from the turn the moon-centred
    hyperbola must supply, ignoring Jupiter's gravity (paper's caveat,
    p. 16): sin(delta/2) = mu_moon / (r_p V_inf^2 + mu_moon) inverted for
    r_p, minus the moon radius. A near-zero turn maps to a huge fictitious
    altitude (e.g. 1.9e6 km at the first Callisto flyby — far outside the
    SOI). Reproduced ONLY to cross-check the turn-angle structure; never a
    physical periapsis, never ingested (mining note section 9 item 3).

    Moon GM / radius are the registry JPL SSD values (the paper prints
    neither), so trailing digits are constant-dependent.
    """
    vi = float(np.linalg.norm(vinf_in))
    vo = float(np.linalg.norm(vinf_out))
    cos_delta = float(np.dot(vinf_in, vinf_out) / (vi * vo))
    delta = math.acos(max(-1.0, min(1.0, cos_delta)))
    sin_half = math.sin(0.5 * delta)
    if sin_half <= 0.0:
        return float("inf")
    sat = SATELLITES[moon]
    r_p = sat.mu_km3_s2 * (1.0 / sin_half - 1.0) / (vi * vi)
    return r_p - sat.radius_eq_km


def vinf_print_tolerance_kms(epoch_days: float, moon: str, radii_km: dict[str, float]) -> float:
    """A-priori |V_inf - printed| tolerance at a flyby, from input precision.

    Tolerance rationale (derived BEFORE asserting; 2026-06-13 sensitivity
    scan in the results note backs every term):

    * DOMINANT — Table 1 mean-motion print quantization: each mean motion is
      printed to 4 decimals (half-ULP 5e-5 rad/day), so a moon's polar angle
      at epoch t is uncertain by up to 5e-5 * t rad. That rotates both the
      Lambert endpoint and the moon velocity vector; each contributes up to
      |v| * dtheta to the V-infinity difference, with |v| bounded by the moon
      circular speed for the moon term and empirically the same scale for
      the Lambert-endpoint term => 2 * v_moon * 5e-5 * t. A half-ULP
      perturbation of a single mean motion indeed swings the worst Member C
      residual between 1.6e-2 and 1.5e-1 km/s (results note), bracketing the
      observed 4.8e-2.
    * mu_Jupiter ambiguity (unprinted): velocities scale ~ mu^(1/3) at fixed
      mean motions/ToFs; the planet-vs-system GM spread (2e-4 fractional)
      gives <= 1e-3 km/s at 13 km/s — folded into the constant floor below,
      together with the Table 2/4/6 phase and Tables 3/5/7 ToF / V-infinity
      print quantizations (all <= a few 1e-4 km/s equivalent).

    The wrong-member / wrong-branch separation is ~1 km/s (anchoring the
    epochs at t_c0 instead of t_c0 + T, or picking any other Lambert branch,
    fails by >= 0.18 km/s; see LegRecord.selection_margin_kms), so this
    tolerance rejects every wrong reconstruction with >= 4x margin while
    admitting the print-precision noise the inputs cannot beat. NEVER widen
    it to make a member pass — a clean negative is a valid outcome.
    """
    v_moon = math.sqrt(MU_JUPITER_KM3_S2 / radii_km[moon])
    quantization = 2.0 * v_moon * MEAN_MOTION_HALF_ULP_RAD_DAY * epoch_days
    floor = 1.0e-3  # mu_J + phase + ToF + printed-V_inf quantization, km/s
    return quantization + floor


def reproduce_member(
    member: str,
    *,
    mu_jupiter: float = MU_JUPITER_KM3_S2,
    max_revs: int = 4,
) -> MemberReproduction:
    """Reproduce one idealized member from its printed inputs.

    Route (the direct REPRODUCTION route, not a re-search): place the moons
    by mean motion + Table 2/4/6 phase at the Eq. 16-anchored cumulative
    printed flyby epochs; solve the Jupiter-frame Lambert problem between
    consecutive flyby positions at the printed leg ToFs (enumerating
    revolution counts 0..max_revs and both multi-rev branches, prograde);
    select per leg the solution closest to the printed V-infinity pair
    (identification only — the residuals are then asserted against
    :func:`vinf_print_tolerance_kms`, so a wrong selection cannot pass).
    """
    spec = LIANG_MEMBERS[member]
    radii = {m: derived_radius_km(m, mu_jupiter=mu_jupiter) for m in LIANG_MEAN_MOTIONS_RAD_DAY}
    r_p = spec.perijove_scale * (radii["Europa"] - 10000.0)

    # Initial conic: period = T_Callisto => a equals Callisto's orbit radius
    # (Eq. 14 with T = 2*pi/n_C reduces to Kepler III for Callisto itself).
    a = radii["Callisto"]
    ecc = 1.0 - r_p / a
    n_c = LIANG_MEAN_MOTIONS_RAD_DAY["Callisto"]
    t_callisto = 2.0 * math.pi / n_c  # days

    # Structure flag 1 for Callisto: crossing above the x-axis, E = +pi/2.
    # Kepler: M = E - e*sin(E) = pi/2 - e; t_c0 = M / n_C.
    t_c0 = (0.5 * math.pi - ecc) / n_c
    b = a * math.sqrt(1.0 - ecc * ecc)
    crossing_pos = np.array([-a * ecc, b, 0.0])  # (a(cosE - e), b sinE) at E = pi/2
    crossing_lon = math.atan2(b, -a * ecc)
    # Perifocal velocity sqrt(mu*a)/r * (-sinE, sqrt(1-e^2) cosE) at E = pi/2,
    # r = a: v = sqrt(mu/a) * (-1, 0). Perifocal == inertial (perijove on +x).
    v_conic = np.array([-math.sqrt(mu_jupiter / a), 0.0, 0.0])

    # Phase-convention self-check: Callisto must BE at the crossing at t_c0.
    theta_c = spec.phases_rad["Callisto"] + n_c * t_c0
    angle_residual = abs(
        math.remainder(theta_c - crossing_lon, 2.0 * math.pi)
    )  # wrapped to [-pi, pi]

    # Flyby epochs: table row 0 is Eq. 16's t_c1 = T + t_c0 (n_cycle = 1) —
    # one full spacecraft revolution after perijove departure (module
    # docstring, geometry conventions). Later rows add the printed ToFs.
    epochs = [t_c0 + t_callisto]
    for tof in spec.tofs_days:
        epochs.append(epochs[-1] + tof)

    states = [
        _moon_state(m, t, spec.phases_rad, radii, mu_jupiter)
        for m, t in zip(CGCEC_SEQUENCE, epochs, strict=True)
    ]

    # Sanity: the encounter state at t_c0 + T equals the crossing state (both
    # the conic and Callisto are T-periodic), so flyby 0 sits at crossing_pos.
    crossing_gap_km = float(np.linalg.norm(states[0][0] - crossing_pos))
    if crossing_gap_km > 1.0e-3 * a:
        raise RuntimeError(
            f"Member {member}: Callisto at t_c1 is {crossing_gap_km:.3e} km "
            "from the conic crossing — phase convention violated."
        )

    vinf_in_vec: list[Vec3] = [v_conic - states[0][1]]
    vinf_out_vec: list[Vec3 | None] = []
    legs: list[LegRecord] = []
    for k in range(4):
        r_a, v_a = states[k]
        r_b, v_b = states[k + 1]
        sols = lambert(r_a, r_b, spec.tofs_days[k] * DAY_S, mu=mu_jupiter, max_revs=max_revs)

        def cost(s_v1: Vec3, s_v2: Vec3, k_leg: int = k) -> float:
            return abs(
                float(np.linalg.norm(s_v1 - states[k_leg][1])) - spec.vinf_printed_kms[k_leg]
            ) + abs(
                float(np.linalg.norm(s_v2 - states[k_leg + 1][1]))
                - spec.vinf_printed_kms[k_leg + 1]
            )

        ranked = sorted(sols, key=lambda s: cost(s.v1, s.v2))
        best = ranked[0]
        if len(ranked) > 1:
            margin = cost(ranked[1].v1, ranked[1].v2) - cost(best.v1, best.v2)
        else:
            margin = float("inf")
        vinf_out_vec.append(best.v1 - v_a)
        vinf_in_vec.append(best.v2 - v_b)
        legs.append(
            LegRecord(
                index=k,
                n_revs=best.n_revs,
                branch=best.branch,
                tof_days=spec.tofs_days[k],
                vinf_dep_kms=float(np.linalg.norm(best.v1 - v_a)),
                vinf_arr_kms=float(np.linalg.norm(best.v2 - v_b)),
                selection_margin_kms=margin,
            )
        )
    vinf_out_vec.append(None)  # flyby 4 ends the printed cycle

    flybys: list[FlybyRecord] = []
    max_residual = 0.0
    for k in range(5):
        vi_vec = vinf_in_vec[k]
        vo_vec = vinf_out_vec[k]
        vi = float(np.linalg.norm(vi_vec))
        printed = spec.vinf_printed_kms[k]
        res_in = abs(vi - printed)
        max_residual = max(max_residual, res_in)
        vo: float | None = None
        res_out: float | None = None
        continuity: float | None = None
        altitude: float | None = None
        if vo_vec is not None:
            vo = float(np.linalg.norm(vo_vec))
            res_out = abs(vo - printed)
            continuity = abs(vi - vo)
            altitude = defect_altitude_km(vi_vec, vo_vec, CGCEC_SEQUENCE[k])
            max_residual = max(max_residual, res_out)
        flybys.append(
            FlybyRecord(
                index=k,
                moon=CGCEC_SEQUENCE[k],
                epoch_days=epochs[k],
                vinf_printed_kms=printed,
                vinf_in_kms=vi,
                vinf_out_kms=vo,
                residual_in_kms=res_in,
                residual_out_kms=res_out,
                continuity_kms=continuity,
                defect_altitude_km=altitude,
            )
        )

    return MemberReproduction(
        member=member,
        mu_jupiter=mu_jupiter,
        radii_km=radii,
        perijove_km=r_p,
        eccentricity=ecc,
        t_c0_days=t_c0,
        crossing_longitude_rad=crossing_lon,
        callisto_angle_residual_rad=angle_residual,
        conic_vinf_kms=float(np.linalg.norm(vinf_in_vec[0])),
        cycle_tof_days=float(sum(spec.tofs_days)),
        legs=tuple(legs),
        flybys=tuple(flybys),
        max_vinf_residual_kms=max_residual,
    )
