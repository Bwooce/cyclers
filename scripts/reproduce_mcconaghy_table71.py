"""Reproduce McConaghy 2004 Table 7.1 (S1L1 outbound cycler vehicle 1) on DE440.

Second-source validation of the S1L1 Earth-Mars cycler (#94 follow-on). The #167
work (``src/cyclerfinder/search/s1l1_corrected.py``) reproduced Russell 2004
Appendix C #83 per-leg on DE440; this script applies the same per-leg
reproduction recipe to an INDEPENDENT published itinerary: McConaghy 2004
(Purdue PhD dissertation, UMI 3166673) Chapter 7 Table 7.1 — the 24-encounter,
33-year DE405 itinerary of S1L1 outbound cycler vehicle 1. Data transcribed in
``docs/notes/2026-06-10-mcconaghy-2004-dissertation-mining.md`` (S2.2).

Method (mirrors the #167 recipe, adapted to the published data shape)
----------------------------------------------------------------------
McConaghy prints per encounter: calendar date (day precision), |V_inf| (km/s),
closest-approach ALTITUDE (km; Mars-23's 1,770 km < Mars's radius pins the
altitude convention), and arriving-leg TOF (days). Unlike Russell App-C there
are NO v_inf VECTORS, so per-leg seeding is by Lambert between the published
encounter dates on DE440 (our ephemeris; McConaghy used DE405 — differences
are reported, not tuned away):

1. For each consecutive encounter pair, solve Lambert (multi-rev up to N=2;
   the S1L1 legs sweep >360 deg: McConaghy Table 6.2 gives theta = 657.97 deg
   and 522.29 deg). The DISCRETE (n_revs, branch) choice per leg is selected
   by best match to the published endpoint |V_inf| — a topology selection
   among ~3-5 well-separated conics, reported per leg with the runner-up gap.
   Within the chosen branch nothing is fit: the emerged V_inf at BOTH ends is
   evidence against the published values (published = golden, non-circular).
2. Flyby continuity at every interior encounter: |V_inf_in| vs |V_inf_out|
   must agree for a ballistic flyby (McConaghy's itinerary is ballistic).
3. Implied flyby periapsis: from the bend angle between the incoming/outgoing
   V_inf vectors, sin(delta/2) = 1/(1 + r_p v^2/mu) gives r_p; the implied
   ALTITUDE is checked against the published closest approach. This constraint
   is never used in selection — it is a pure cross-check.
4. Independent n-body: every leg's Lambert departure state is propagated by
   REBOUND/IAS15 (Sun-only, matching the patched-conic cruise model; plus
   real DE440 Mars as a continuous perturber for the E->M legs) and the
   arrival miss vs the real DE440 planet is measured. Mars encounters must
   land inside the 3-Mars-SOI band (~0.0116 AU) — the SAME band as #165/#167,
   never loosened.
5. Structural check: the E->E (g) free-return legs must stay sub-Mars
   (aphelion < 1.45 AU) and far from real DE440 Mars — the corrected-topology
   signature from #166/#167.

Conventions / honesty
---------------------
* Encounter epochs are taken at 12:00 TDB on the published calendar date
  (J2000 + whole days). Day-precision quantisation (+-0.5 d) is the dominant
  expected residual source (~0.01-0.1 km/s on V_inf).
* Published dates are INPUTS (constraints), so no date residual can emerge;
  the integer-day TOF cross-check (date difference vs printed leg TOF) is
  still asserted. What emerges: V_inf at both leg ends, flyby |V_inf|
  continuity, implied periapsis altitude, n-body miss distances.
* DE405 (McConaghy) vs DE440 (ours) differences are reported honestly.
* NO catalogue writeback. Results go to
  ``docs/notes/2026-06-10-mcconaghy-table71-reproduction.md``.

Run: ``uv run python scripts/reproduce_mcconaghy_table71.py`` (add
``--no-nbody`` to skip the REBOUND stage).

FINDING (2026-06-10) — see docs/notes/2026-06-10-mcconaghy-table71-reproduction.md:
  Encounters 1-19 (6 of ~7.7 cycles) REPRODUCE: |V_inf| residual <= 0.195 km/s
  (mean 0.05), flyby continuity <= 0.19, implied periapsis altitude within
  ~1,800 km of the printed closest approach, all 8 Mars encounters in the
  3-SOI band on IAS15 (Sun-only ~0 km by construction; Mars-perturbed
  8,000-53,000 km), every E->E g leg sub-Mars (aphelion <= 1.33 AU).
  Encounters 20-24 as transcribed DO NOT reproduce (residual growing
  0.15 -> 1.56 km/s; implied CA inconsistent at E22/M23). Adjudication:
  Russell 2004 App-C #83 (independently DE440-confirmed in #167) overlaps
  Table 7.1's dates 2027-2038 to <= 1.9 d; at those shared nodes Russell's
  printed v_inf matches OUR emerged values to <= 0.16 km/s and disagrees
  with the transcribed Table 7.1 tail by up to 1.7 km/s. The tail rows as
  transcribed are internally inconsistent (a ~1.1 km/s shift would need
  ~15-day date errors; dates/TOFs are mutually consistent) — flagged for
  re-verification against the dissertation PDF, NOT tuned away.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import date as _date

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate as kepler_propagate
from cyclerfinder.core.lambert import LambertError, LambertSolution, lambert

Vec3 = NDArray[np.float64]


def _lambert_candidates(r1: Vec3, r2: Vec3, tof: float) -> list[LambertSolution]:
    """All candidate conics for one leg: in-house solver, lamberthub fallback.

    The in-house universal-variable solver computes the single-rev branch first
    and raises on the longest (~900 d) legs where that branch sits hard against
    the z=(2*pi)^2 singularity — taking the (feasible) multi-rev branches down
    with it. Fall back to ``lamberthub.izzo2015`` (the same dev dependency the
    repo's Lambert crosscheck gate uses) to enumerate the branches there.
    """
    try:
        return lambert(r1, r2, tof, max_revs=2)
    except LambertError:
        from lamberthub import izzo2015  # type: ignore[import-untyped]

        sols: list[LambertSolution] = []
        for n in (0, 1, 2):
            for low_path in (True, False):
                if n == 0 and not low_path:
                    continue  # single-rev has one branch
                try:
                    v1, v2 = izzo2015(
                        MU_SUN_KM3_S2, r1, r2, tof, M=n, prograde=True, low_path=low_path
                    )
                except Exception:  # lamberthub raises on infeasible revolution counts
                    continue
                branch = "single" if n == 0 else ("low" if low_path else "high")
                sols.append(
                    LambertSolution(
                        n_revs=n,
                        branch=branch,
                        v1=np.asarray(v1, dtype=np.float64),
                        v2=np.asarray(v2, dtype=np.float64),
                    )
                )
        if not sols:
            raise
        return sols


# --- McConaghy 2004 Table 7.1, printed p.149 — VERBATIM (golden side) ---------
# (encounter_no, body, iso_date, vinf_kms, closest_approach_alt_km, leg_tof_days)
# closest approach is altitude (Mars-23: 1,770 km < Mars radius 3,396 km).
# leg_tof_days is the TOF of the leg ARRIVING at this encounter.
TABLE71: tuple[tuple[int, str, str, float, float | None, float | None], ...] = (
    (1, "E", "2005-08-13", 4.01, None, None),
    (2, "M", "2006-02-27", 3.02, 4816.0, 198.0),
    (3, "E", "2008-06-09", 6.89, 20130.0, 833.0),
    (4, "E", "2009-12-03", 6.90, 31110.0, 541.0),
    (5, "M", "2010-06-06", 4.31, 17710.0, 186.0),
    (6, "E", "2012-08-24", 6.42, 26490.0, 809.0),
    (7, "E", "2014-02-14", 6.43, 41520.0, 539.0),
    (8, "M", "2014-07-03", 7.14, 12190.0, 138.0),
    (9, "E", "2016-12-09", 4.01, 27730.0, 890.0),
    (10, "E", "2018-05-22", 4.03, 19920.0, 530.0),
    (11, "M", "2018-09-15", 6.47, 11580.0, 115.0),
    (12, "E", "2021-04-06", 4.61, 22990.0, 934.0),
    (13, "E", "2022-09-20", 4.59, 14780.0, 532.0),
    (14, "M", "2023-05-01", 2.77, 7601.0, 223.0),
    (15, "E", "2025-07-02", 7.08, 23860.0, 793.0),
    (16, "E", "2026-12-26", 7.09, 35120.0, 542.0),
    (17, "M", "2027-06-14", 5.27, 13840.0, 170.0),
    (18, "E", "2029-09-21", 5.80, 26850.0, 830.0),
    (19, "E", "2031-03-12", 5.80, 37520.0, 537.0),
    (20, "M", "2031-07-15", 7.85, 8802.0, 125.0),
    (21, "E", "2034-01-15", 4.21, 24870.0, 915.0),
    (22, "E", "2035-06-28", 4.20, 2756.0, 529.0),
    (23, "M", "2035-11-12", 5.87, 1770.0, 137.0),
    (24, "E", "2038-05-06", 7.23, None, 906.0),
)

# Same 3-Mars-SOI confirmation band as #165/#167 — kept IDENTICAL, never loosened.
MARS_SOI_AU: float = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
ENCOUNTER_MISS_TOL_AU: float = 3.0 * MARS_SOI_AU  # ~0.0116 AU


def _t_sec(iso: str) -> float:
    """Encounter epoch: 12:00 TDB on the published calendar date, s past J2000."""
    d = _date.fromisoformat(iso)
    return float((d - _date(2000, 1, 1)).days) * SECONDS_PER_DAY


@dataclass(frozen=True)
class LegResult:
    """One Lambert-built leg between consecutive published encounters."""

    dep_no: int
    arr_no: int
    dep_body: str
    arr_body: str
    t0_sec: float
    t1_sec: float
    tof_days: float
    n_revs: int
    branch: str
    r0_km: Vec3
    v0_km_s: Vec3  # heliocentric departure velocity (Lambert v1)
    vinf_dep_vec: Vec3
    vinf_arr_vec: Vec3
    vinf_dep_kms: float
    vinf_arr_kms: float
    pub_vinf_dep: float
    pub_vinf_arr: float
    runner_up_gap_kms: float  # score gap to the 2nd-best (n_revs, branch) choice

    @property
    def kind(self) -> str:
        return f"{self.dep_body}->{self.arr_body}"


def build_legs(ephem: Ephemeris) -> list[LegResult]:
    """Build all 23 legs by Lambert between published encounter dates on DE440.

    The (n_revs, branch) choice is discrete-topological, selected by summed
    endpoint |V_inf| match to the published values; everything continuous
    (the V_inf themselves) then EMERGES from the chosen conic.
    """
    legs: list[LegResult] = []
    for i in range(len(TABLE71) - 1):
        no1, b1, d1, pub_v1, _ca1, _tof1 = TABLE71[i]
        no2, b2, d2, pub_v2, _ca2, pub_tof2 = TABLE71[i + 1]
        t0 = _t_sec(d1)
        t1 = _t_sec(d2)
        tof_days = (t1 - t0) / SECONDS_PER_DAY
        if pub_tof2 is not None and abs(tof_days - pub_tof2) > 1.5:
            # +-1 d is day-rounding of the printed dates vs the printed TOF
            # (e.g. leg 3->4: dates give 542 d, printed TOF 541 d); anything
            # larger would be a transcription error.
            raise ValueError(
                f"leg {no1}->{no2}: date difference {tof_days:.0f} d != printed TOF "
                f"{pub_tof2:.0f} d — transcription error"
            )
        if pub_tof2 is not None and abs(tof_days - pub_tof2) > 0.5:
            print(
                f"  [TOF rounding] leg {no1}->{no2}: dates give {tof_days:.0f} d, "
                f"printed TOF {pub_tof2:.0f} d (+-1 d day-quantisation)"
            )
        r1_t, vp1_t = ephem.state(b1, t0)
        r2_t, vp2_t = ephem.state(b2, t1)
        r1 = np.asarray(r1_t, dtype=np.float64)
        vp1 = np.asarray(vp1_t, dtype=np.float64)
        r2 = np.asarray(r2_t, dtype=np.float64)
        vp2 = np.asarray(vp2_t, dtype=np.float64)

        sols = _lambert_candidates(r1, r2, t1 - t0)

        def _score(
            s: LambertSolution,
            _vp1: Vec3 = vp1,
            _vp2: Vec3 = vp2,
            _pv1: float = pub_v1,
            _pv2: float = pub_v2,
        ) -> float:
            vd = float(np.linalg.norm(np.asarray(s.v1) - _vp1))
            va = float(np.linalg.norm(np.asarray(s.v2) - _vp2))
            return abs(vd - _pv1) + abs(va - _pv2)

        ranked = sorted(sols, key=_score)
        best = ranked[0]
        gap = _score(ranked[1]) - _score(best) if len(ranked) > 1 else float("inf")
        vinf_dep_vec = np.asarray(best.v1, dtype=np.float64) - vp1
        vinf_arr_vec = np.asarray(best.v2, dtype=np.float64) - vp2
        legs.append(
            LegResult(
                dep_no=no1,
                arr_no=no2,
                dep_body=b1,
                arr_body=b2,
                t0_sec=t0,
                t1_sec=t1,
                tof_days=tof_days,
                n_revs=best.n_revs,
                branch=best.branch,
                r0_km=r1,
                v0_km_s=np.asarray(best.v1, dtype=np.float64),
                vinf_dep_vec=vinf_dep_vec,
                vinf_arr_vec=vinf_arr_vec,
                vinf_dep_kms=float(np.linalg.norm(vinf_dep_vec)),
                vinf_arr_kms=float(np.linalg.norm(vinf_arr_vec)),
                pub_vinf_dep=pub_v1,
                pub_vinf_arr=pub_v2,
                runner_up_gap_kms=gap,
            )
        )
    return legs


@dataclass(frozen=True)
class FlybyCheck:
    """Ballistic-flyby consistency at one interior encounter."""

    no: int
    body: str
    date_iso: str
    pub_vinf: float
    vinf_in: float
    vinf_out: float
    continuity_kms: float  # | |v_inf_in| - |v_inf_out| |  (0 for a perfect ballistic flyby)
    bend_deg: float
    implied_alt_km: float
    pub_ca_alt_km: float | None


def flyby_checks(legs: list[LegResult]) -> list[FlybyCheck]:
    """V_inf continuity + implied periapsis altitude at every interior encounter."""
    out: list[FlybyCheck] = []
    by_arr = {leg.arr_no: leg for leg in legs}
    by_dep = {leg.dep_no: leg for leg in legs}
    for no, body, d, pub_v, ca, _tof in TABLE71:
        if no not in by_arr or no not in by_dep:
            continue  # first/last encounter: only one adjoining leg
        leg_in = by_arr[no]
        leg_out = by_dep[no]
        v_in = leg_in.vinf_arr_vec
        v_out = leg_out.vinf_dep_vec
        m_in = float(np.linalg.norm(v_in))
        m_out = float(np.linalg.norm(v_out))
        cos_b = float(np.dot(v_in, v_out) / (m_in * m_out))
        bend = float(np.arccos(max(min(cos_b, 1.0), -1.0)))
        v_eff = 0.5 * (m_in + m_out)
        mu = PLANETS[body].mu_km3_s2
        sin_half = float(np.sin(bend / 2.0))
        if sin_half > 1.0e-12:
            r_p = (mu / v_eff**2) * (1.0 / sin_half - 1.0)
            implied_alt = float(r_p - PLANETS[body].radius_eq_km)
        else:
            implied_alt = float("inf")
        out.append(
            FlybyCheck(
                no=no,
                body=body,
                date_iso=d,
                pub_vinf=pub_v,
                vinf_in=m_in,
                vinf_out=m_out,
                continuity_kms=abs(m_in - m_out),
                bend_deg=float(np.degrees(bend)),
                implied_alt_km=implied_alt,
                pub_ca_alt_km=ca,
            )
        )
    return out


@dataclass(frozen=True)
class GArcCheck:
    """Sub-Mars structural check of one E->E (g) free-return leg."""

    dep_no: int
    arr_no: int
    aphelion_au: float
    closest_mars_au: float


def g_arc_checks(ephem: Ephemeris, legs: list[LegResult], n_samples: int = 200) -> list[GArcCheck]:
    """Aphelion + closest approach to real DE440 Mars over each E->E leg."""
    out: list[GArcCheck] = []
    for leg in legs:
        if not (leg.dep_body == "E" and leg.arr_body == "E"):
            continue
        dur = leg.t1_sec - leg.t0_sec
        aph = 0.0
        closest = float("inf")
        for k in range(n_samples + 1):
            dt = dur * k / n_samples
            rk, _ = kepler_propagate(leg.r0_km, leg.v0_km_s, dt)
            aph = max(aph, float(np.linalg.norm(rk) / AU_KM))
            r_m, _ = ephem.state("M", leg.t0_sec + dt)
            closest = min(closest, float(np.linalg.norm(rk - np.asarray(r_m)) / AU_KM))
        out.append(
            GArcCheck(
                dep_no=leg.dep_no, arr_no=leg.arr_no, aphelion_au=aph, closest_mars_au=closest
            )
        )
    return out


@dataclass(frozen=True)
class NBodyCheck:
    """Independent REBOUND/IAS15 propagation of one leg to its arrival epoch."""

    dep_no: int
    arr_no: int
    arr_body: str
    miss_sun_only_km: float
    vinf_sun_only_kms: float
    converged: bool
    energy_rel_drift: float
    miss_mars_pert_km: float | None  # E->M legs only (Mars as continuous perturber)


def nbody_checks(ephem: Ephemeris, legs: list[LegResult]) -> list[NBodyCheck]:
    """IAS15 Sun-only for every leg; + real-Mars-perturbed for the E->M legs.

    The arrival miss vs real DE440 planet is EVIDENCE from an integrator that
    shares none of the Lambert solver's machinery. Mars-departing legs cannot
    take Mars as a continuous perturber (the patched-conic start is at the
    planet position), matching the #167 convention exactly.
    """
    from cyclerfinder.nbody.forces import RailsEphemerisCache
    from cyclerfinder.nbody.propagator import RestrictedNBody

    prop = RestrictedNBody("rebound")
    out: list[NBodyCheck] = []
    for leg in legs:
        arc = prop.propagate(
            leg.r0_km,
            leg.v0_km_s,
            leg.t0_sec,
            leg.t1_sec,
            bodies=(),
            ephem=None,
            accuracy=1e-11,
            max_wall_sec=120.0,
        )
        r_pl_t, v_pl_t = ephem.state(leg.arr_body, leg.t1_sec)
        r_pl = np.asarray(r_pl_t, dtype=np.float64)
        v_pl = np.asarray(v_pl_t, dtype=np.float64)
        miss_s = float(np.linalg.norm(np.asarray(arc.r_km) - r_pl))
        vinf_s = float(np.linalg.norm(np.asarray(arc.v_km_s) - v_pl))

        miss_m: float | None = None
        if leg.dep_body == "E" and leg.arr_body == "M":
            cache = RailsEphemerisCache(
                ("M",),
                ephem,
                leg.t0_sec - 5 * SECONDS_PER_DAY,
                leg.t1_sec + 10 * SECONDS_PER_DAY,
            )
            arc_m = prop.propagate(
                leg.r0_km,
                leg.v0_km_s,
                leg.t0_sec,
                leg.t1_sec,
                bodies=("M",),
                ephem=ephem,
                cache=cache,
                accuracy=1e-10,
                max_wall_sec=120.0,
            )
            miss_m = float(np.linalg.norm(np.asarray(arc_m.r_km) - r_pl))
        out.append(
            NBodyCheck(
                dep_no=leg.dep_no,
                arr_no=leg.arr_no,
                arr_body=leg.arr_body,
                miss_sun_only_km=miss_s,
                vinf_sun_only_kms=vinf_s,
                converged=bool(arc.converged),
                energy_rel_drift=float(arc.energy_rel_drift),
                miss_mars_pert_km=miss_m,
            )
        )
    return out


@dataclass(frozen=True)
class RussellOverlapRow:
    """One Table 7.1 encounter that coincides with a Russell App-C #83 node.

    Russell 2004 Appendix C #83 (the trajectory #167 independently confirmed
    on DE440 to 4-decimal v_inf) spans 2025-2055 and its early nodes fall on
    the SAME calendar dates as McConaghy Table 7.1's encounters 15-24 — the two
    publications describe near-identical members of the same S1L1 vehicle.
    Russell's printed per-leg |v_inf| is therefore a THIRD, independent,
    already-DE440-confirmed reference for exactly the encounters where the
    transcribed Table 7.1 V_inf stops matching our Lambert reproduction.
    """

    enc_no: int
    body: str
    date_iso: str
    russell_date_offset_d: float  # Russell node date minus Table 7.1 date
    russell_vinf_kms: float  # |v_inf| of the App-C leg STARTING at this node
    pub_vinf: float  # McConaghy Table 7.1 printed value
    emerged_out: float | None  # our Lambert departure |v_inf| at this encounter


def russell_overlap(legs: list[LegResult]) -> list[RussellOverlapRow]:
    """Match Table 7.1 encounters to Russell App-C #83 nodes within +-8 days."""
    from cyclerfinder.search.s1l1_corrected import APPC_EPOCH_DAYS, APPC_LEGS

    out_by_dep = {leg.dep_no: leg for leg in legs}
    rows: list[RussellOverlapRow] = []
    for no, body, d, pub_v, _ca, _tof in TABLE71:
        t_enc_d = _t_sec(d) / SECONDS_PER_DAY
        for _leg_no, r_body, ts, vinf in APPC_LEGS:
            if vinf is None or r_body != body:
                continue
            dt_d = (APPC_EPOCH_DAYS + ts) - t_enc_d
            if abs(dt_d) > 8.0:
                continue
            leg_out = out_by_dep.get(no)
            rows.append(
                RussellOverlapRow(
                    enc_no=no,
                    body=body,
                    date_iso=d,
                    russell_date_offset_d=dt_d,
                    russell_vinf_kms=float(np.linalg.norm(np.asarray(vinf))),
                    pub_vinf=pub_v,
                    emerged_out=leg_out.vinf_dep_kms if leg_out is not None else None,
                )
            )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-nbody", action="store_true", help="skip the REBOUND/IAS15 stage")
    args = ap.parse_args()

    print("=== McConaghy 2004 Table 7.1 reproduction on DE440 (per-leg Lambert) ===")
    band_str = f"{ENCOUNTER_MISS_TOL_AU:.4f} AU = {ENCOUNTER_MISS_TOL_AU * AU_KM:.3e} km"
    print(f"3-Mars-SOI band: {band_str}")
    print()
    ephem = Ephemeris("astropy")
    legs = build_legs(ephem)
    print()
    print("--- Legs (emerged V_inf vs published; branch = discrete Lambert choice) ---")
    print(
        "leg    kind  TOF_d  Nrev/branch  "
        "vinf_dep  pub_dep  d_dep | vinf_arr  pub_arr  d_arr | runner-up gap"
    )
    for leg in legs:
        print(
            f"{leg.dep_no:>2}->{leg.arr_no:<3} {leg.kind}  {leg.tof_days:5.0f}  "
            f"{leg.n_revs}/{leg.branch:<6}     "
            f"{leg.vinf_dep_kms:7.3f}  {leg.pub_vinf_dep:6.2f}  "
            f"{leg.vinf_dep_kms - leg.pub_vinf_dep:+6.3f} | "
            f"{leg.vinf_arr_kms:7.3f}  {leg.pub_vinf_arr:6.2f}  "
            f"{leg.vinf_arr_kms - leg.pub_vinf_arr:+6.3f} | "
            f"{leg.runner_up_gap_kms:6.2f} km/s"
        )
    d_all = [
        abs(x)
        for leg in legs
        for x in (leg.vinf_dep_kms - leg.pub_vinf_dep, leg.vinf_arr_kms - leg.pub_vinf_arr)
    ]
    print(
        f"\n|V_inf| residual over {len(d_all)} leg-endpoints: "
        f"max {max(d_all):.3f}, mean {float(np.mean(d_all)):.3f} km/s"
    )

    print("\n--- Interior flybys (ballistic continuity + implied periapsis altitude) ---")
    print(
        "enc  body  date        pub_vinf  vinf_in  vinf_out  |in-out|  bend_deg"
        "  implied_alt_km  pub_CA_km"
    )
    fbs = flyby_checks(legs)
    for fb in fbs:
        ca = f"{fb.pub_ca_alt_km:9.0f}" if fb.pub_ca_alt_km is not None else "      n/a"
        print(
            f"{fb.no:>3}  {fb.body}     {fb.date_iso}  {fb.pub_vinf:6.2f}   "
            f"{fb.vinf_in:7.3f}  {fb.vinf_out:7.3f}   {fb.continuity_kms:6.3f}   "
            f"{fb.bend_deg:7.2f}   {fb.implied_alt_km:12.0f}  {ca}"
        )
    cont = [fb.continuity_kms for fb in fbs]
    print(f"\nflyby |V_inf| continuity: max {max(cont):.3f}, mean {float(np.mean(cont)):.3f} km/s")

    print("\n--- E->E (g) free-return legs: sub-Mars structural check ---")
    print("leg     aphelion_AU  closest_DE440_Mars_AU")
    gs = g_arc_checks(ephem, legs)
    for g in gs:
        print(f"{g.dep_no:>2}->{g.arr_no:<3}  {g.aphelion_au:9.3f}  {g.closest_mars_au:14.3f}")

    print("\n--- Russell App-C #83 overlap (third source, DE440-confirmed in #167) ---")
    print("enc  body  date        Russell_dt_d  Russell_vinf  McC_printed  our_emerged_out")
    ovs = russell_overlap(legs)
    for ov in ovs:
        eo = f"{ov.emerged_out:11.3f}" if ov.emerged_out is not None else "        n/a"
        print(
            f"{ov.enc_no:>3}  {ov.body}     {ov.date_iso}  {ov.russell_date_offset_d:+9.1f}   "
            f"{ov.russell_vinf_kms:10.3f}   {ov.pub_vinf:8.2f}   {eo}"
        )

    nbs: list[NBodyCheck] = []
    if not args.no_nbody:
        print("\n--- Independent n-body (REBOUND/IAS15 over DE440) ---")
        print("leg    arr  miss_SunOnly_km  vinf_SunOnly  miss_MarsPert_km  conv  e_drift")
        nbs = nbody_checks(ephem, legs)
        for nb in nbs:
            mm = (
                f"{nb.miss_mars_pert_km:13.0f}"
                if nb.miss_mars_pert_km is not None
                else ("          n/a")
            )
            print(
                f"{nb.dep_no:>2}->{nb.arr_no:<3} {nb.arr_body}   "
                f"{nb.miss_sun_only_km:12.1f}   {nb.vinf_sun_only_kms:9.3f}   {mm}   "
                f"{nb.converged!s:5} {nb.energy_rel_drift:8.1e}"
            )
        mars_nbs = [nb for nb in nbs if nb.arr_body == "M"]
        band_km = ENCOUNTER_MISS_TOL_AU * AU_KM
        in_band = all(
            nb.miss_sun_only_km < band_km
            and (nb.miss_mars_pert_km is None or nb.miss_mars_pert_km < band_km)
            for nb in mars_nbs
        )
        print(
            f"\nMars encounters in 3-SOI band ({band_km:.2e} km): "
            f"{'ALL ' + str(len(mars_nbs)) + ' IN-BAND' if in_band else 'OUT-OF-BAND PRESENT'}"
        )

    print("\n=== summary ===")
    print(
        f"max |V_inf| residual vs published: {max(d_all):.3f} km/s "
        f"(mean {float(np.mean(d_all)):.3f})"
    )
    print(f"max flyby continuity defect:       {max(cont):.3f} km/s")
    print(f"g-leg aphelion max:                {max(g.aphelion_au for g in gs):.3f} AU (Mars 1.52)")
    if nbs:
        print(f"n-body legs converged:             {sum(nb.converged for nb in nbs)}/{len(nbs)}")


if __name__ == "__main__":
    main()
