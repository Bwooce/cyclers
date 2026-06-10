"""Self-seeding longitude-rendezvous construction (#173).

Replaces the App-C *printed* seed with a SEARCH. Given only a row's descriptor
(the two-arc ``g`` / ``G`` free-return shape + the corrected #167 topology) this
module FINDS the Earth-departure epoch and departure v_inf vector at which real
DE440 Mars actually sits at the spacecraft's G-arc encounter longitude — the
longitude-rendezvous constraint Russell's Appendix-C block supplied for S1L1 — so
the unsourced rows (the 7 russell-ch4 without an App-C block, the ~194 ocampo
members) can be reached without a published real-eph seed.

The genuinely-new surface is small and bounded (design §3): a one-variable
longitude residual :func:`residual_lon`, a synodic-period bracket-and-refine scan
:func:`synodic_longitude_scan` that ENUMERATES (does not optimise) the phase
candidates, an optional Lambert refinement onto TRUE DE440 Mars position, and the
assembly of the full multi-term on-family gate (:func:`on_family`). Everything
upstream is REUSED read-only:

* descriptor -> two-arc (a, e, n_rev) coplanar -> DE440 G-arc SHAPE:
  :mod:`cyclerfinder.search.free_return_chain` (#163) +
  :mod:`cyclerfinder.search.continuation_chain` (#164) (Stage A);
* the per-leg reconstruction recipe ``v_sc = v_planet(DE440) + v_inf`` and the
  longitude helper convention: :mod:`cyclerfinder.search.s1l1_corrected` (#167);
* :func:`cyclerfinder.core.lambert.lambert` for the Stage-B refinement;
* :class:`cyclerfinder.core.ephemeris.Ephemeris` ("astropy" => DE440).

GOLDEN / HONESTY (binding, design §5). When validated against S1L1 the EXPECTED
side is the row's OWN App-C-confirmed geometry (epoch 2026-12-15, Mars longitude
201.0deg on 2027-06-13, Mars v_inf 5.248 breathing 3.2-8.0). The self-seed's found
epoch / v_inf / miss are EVIDENCE, never imposed. The search itself NEVER reads any
``APPC_*`` constant — it consumes only the descriptor. An OFF-FAMILY / EMPTY-SET
outcome is a first-class success, not a failure; tolerances and bands are NEVER
loosened to manufacture a PASS.

Pure: depends only on core (constants, ephemeris, kepler, lambert) and the reused
search modules. No edit to any existing src/ or core/ file; no catalogue writeback.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate as kepler_propagate
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.continuation_chain import continuation_chain_correct
from cyclerfinder.search.free_return import _crossing, _true_to_mean
from cyclerfinder.search.free_return_chain import (
    _earth_vinf_vector,
    free_return_chain_correct,
)

Vec3 = np.ndarray


# Earth-Mars synodic period. Derived (NOT a hardcoded literal) from the two bodies'
# sidereal periods via 1/T_syn = |1/T_E - 1/T_M|, with the periods taken from the
# live PLANETS mean motions (deg/day) so it tracks the sourced orbital elements.
# Numerically ~779.9 d (~2.135 yr), the window the longitude scan sweeps (design §1.2).
def _synodic_period_days(inner: str = "E", outer: str = "M") -> float:
    """Earth-Mars synodic period (days) from the sourced PLANETS mean motions."""
    t_in = 360.0 / PLANETS[inner].mean_motion_deg_day
    t_out = 360.0 / PLANETS[outer].mean_motion_deg_day
    return 1.0 / abs(1.0 / t_in - 1.0 / t_out)


SYNODIC_PERIOD_DAYS: float = _synodic_period_days()


def _lon_deg(r: Vec3) -> float:
    """Heliocentric ecliptic longitude (deg, [0, 360)) of a position vector.

    Re-implemented locally (the 3-line arctan2) rather than importing the private
    ``s1l1_corrected._lon_deg`` — keeps the dependency one-way and explicit.
    """
    return float(np.degrees(np.arctan2(r[1], r[0])) % 360.0)


def _wrap_deg(angle: float) -> float:
    """Wrap a longitude difference to (-180, 180] degrees."""
    a = (angle + 180.0) % 360.0 - 180.0
    return 180.0 if a == -180.0 else a


# ---------------------------------------------------------------------------
# Stage A — descriptor -> coplanar -> DE440 G-arc SHAPE (reused #163/#164).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GArcShape:
    """The G (Earth-Mars transit) arc SHAPE extracted from a row's descriptor.

    Produced by :func:`g_arc_shape` from the reused two-arc continuation (#163/#164)
    — the *epoch-free* geometry the longitude search then phases against real DE440
    Mars. ``vinf_e_vec_ref`` is the heliocentric departure v_inf vector in the
    Earth-velocity-aligned reference frame (so it can be re-pointed at real Earth's
    instantaneous velocity at any candidate epoch); ``tof_g_days`` is the emerged
    E->M transit time of the G arc.
    """

    a_au: float
    e: float
    n_rev: int
    tof_g_days: float
    vinf_e_mag: float
    vinf_m_mag: float
    vinf_e_vec_ref: Vec3  # departure v_inf in Earth-orbit (r_hat, t_hat) frame
    vinf_e_anchor: float
    vinf_m_anchor: float
    branch: str = "short"  # G-arc Mars-crossing branch: short | long | shortN | longN
    g_revs: int = 0  # full G-arc revolutions added before the Mars crossing


def g_arc_shape(
    aphelion_au: float,
    g_tof_years: float,
    big_g_tof_years: float,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> GArcShape:
    """Stage A: turn a row descriptor into the G-arc SHAPE (coplanar -> DE440).

    REUSES the #163 two-arc circular closure to seed and the #164 continuation to
    walk it to the DE440-consistent J2000 eccentric/inclined model, then reads the
    SECOND (G, uppercase, Mars-transit) arc's converged ``(a, e, n_rev)`` and its
    emerged E->M transit ToF. The departure v_inf vector direction is taken from the
    coplanar G-arc Earth crossing (:func:`free_return_chain._earth_vinf_vector`),
    expressed in the Earth-orbit ``(r_hat, t_hat)`` frame so it can be re-pointed at
    real Earth's velocity at any candidate epoch.

    EXPECTED = the SOURCED anchors; the emerged shape is EVIDENCE. No APPC_* read.
    """
    seed = free_return_chain_correct(
        aphelion_au, g_tof_years, big_g_tof_years, vinf_e_anchor, vinf_m_anchor, mu=mu
    )
    a1, e1, a2, e2 = (
        seed.arcs[0].a_au,
        seed.arcs[0].e,
        seed.arcs[1].a_au,
        seed.arcs[1].e,
    )
    cont = continuation_chain_correct(
        a1,
        e1,
        a2,
        e2,
        0.0,
        g_tof_years,
        big_g_tof_years,
        vinf_e_anchor,
        vinf_m_anchor,
        mu=mu,
    )
    final = cont.best_final
    if final is not None:
        g_a, g_e = final.a2, final.e2
        n_rev = final.arc2.n_rev
        tof_g_days = final.arc2.geometry.tof_em_days
        vinf_e_mag = final.arc2.vinf_e
        vinf_m_mag = final.arc2.vinf_m
    else:
        # Continuation diverged — fall back to the circular seed's G arc so the
        # shape is always well-formed (the search then reports the off-family /
        # empty-set outcome downstream, never a crash).
        g_a, g_e = a2, e2
        n_rev = seed.arcs[1].n_rev
        tof_g_days = seed.arcs[1].transfer_tof_days
        vinf_e_mag = seed.arcs[1].vinf_e
        vinf_m_mag = seed.arcs[1].vinf_m

    # Departure v_inf vector direction from the coplanar G-arc Earth crossing,
    # expressed in the (r_hat, t_hat) Earth-orbit frame at the crossing. At the
    # crossing Earth's circular position is along +x, velocity along +y, so the
    # reference-frame vector is exactly the returned heliocentric v_inf vector.
    vinf_vec_ref = _earth_vinf_vector(g_a, g_e, mu=mu)
    # Re-scale its magnitude to the emerged G-arc Earth v_inf (the shape's own
    # departure v_inf), keeping the coplanar direction.
    nrm = float(np.linalg.norm(vinf_vec_ref))
    if nrm > 0.0:
        vinf_vec_ref = vinf_vec_ref / nrm * vinf_e_mag

    return GArcShape(
        a_au=float(g_a),
        e=float(g_e),
        n_rev=int(n_rev),
        tof_g_days=float(tof_g_days),
        vinf_e_mag=float(vinf_e_mag),
        vinf_m_mag=float(vinf_m_mag),
        vinf_e_vec_ref=np.asarray(vinf_vec_ref, dtype=np.float64),
        vinf_e_anchor=float(vinf_e_anchor),
        vinf_m_anchor=float(vinf_m_anchor),
    )


# ---------------------------------------------------------------------------
# Stage A multi-rev extension (#177) — enumerate the G-arc Mars-crossing branches.
#
# #173's single-branch Stage A took only the SHORT-way (inbound) Mars radial
# crossing of the converged G arc, giving one transit ToF. A long-transit row (e.g.
# 6.44Gg3: real-eph 262 d vs the coplanar short-way 131 d) is OFF-FAMILY against that
# one branch but its real-eph transit signature is reproduced by a DIFFERENT branch
# of the SAME (a, e) shape: the LONG-way crossing (after aphelion, ~292 d) and/or a
# multi-rev branch (k full G-arc revolutions added before the crossing). This
# extension enumerates all those branch ToFs from the converged shape so a row is
# matched against ALL its branch transits before being declared OFF-FAMILY (#177
# build 1). The branch ToFs feed the Stage-B Lambert refinement (which already
# accepts ``max_revs`` from :mod:`core.lambert`, the M-L milestone multi-rev solver).
# ---------------------------------------------------------------------------


def _g_arc_branch_transits(
    a_au: float, e: float, *, max_g_revs: int = 1, mu: float = MU_SUN_KM3_S2
) -> list[tuple[str, int, float]]:
    """Enumerate the E->M transit ToFs of every Mars-crossing branch of ``(a, e)``.

    A transfer ellipse that reaches Mars crosses Mars's radius at TWO true
    anomalies: the inbound (``+nu_M``, the SHORT way, before aphelion) and the
    outbound (``2*pi - nu_M``, the LONG way, after aphelion). Each branch can also
    be reached after ``k`` full revolutions (``k`` in ``[0, max_g_revs]``). Returns
    ``(branch_label, k, tof_days)`` for every feasible branch, sorted by ToF — the
    candidate transit shapes Stage A offers the longitude search. Raises nothing:
    an ``(a, e)`` that does not reach Mars yields an empty list (a clean negative).
    """
    a_km = a_au * AU_KM
    try:
        nu_e = _crossing(a_km, e, PLANETS["E"].sma_au * AU_KM)
        nu_m = _crossing(a_km, e, PLANETS["M"].sma_au * AU_KM)
    except ValueError:
        return []
    n = np.sqrt(mu / a_km**3)  # mean motion, rad/s
    period_days = (2.0 * np.pi / n) / SECONDS_PER_DAY
    m_e = _true_to_mean(nu_e, e)
    out: list[tuple[str, int, float]] = []
    for label, nu_cross in (("short", nu_m), ("long", 2.0 * np.pi - nu_m)):
        dm = (_true_to_mean(nu_cross, e) - m_e) % (2.0 * np.pi)
        base_tof = (dm / n) / SECONDS_PER_DAY
        for k in range(max_g_revs + 1):
            out.append((label, k, base_tof + k * period_days))
    out.sort(key=lambda t: t[2])
    return out


def g_arc_branches(
    aphelion_au: float,
    g_tof_years: float,
    big_g_tof_years: float,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    *,
    max_g_revs: int = 1,
    mu: float = MU_SUN_KM3_S2,
) -> list[GArcShape]:
    """Stage A (multi-rev #177): a row descriptor -> ALL its G-arc branch SHAPES.

    Runs the #173 :func:`g_arc_shape` once to get the converged ``(a, e, n_rev)`` and
    the departure v_inf vector, then enumerates every Mars-crossing branch
    (:func:`_g_arc_branch_transits`) as a separate :class:`GArcShape` differing only
    in ``tof_g_days`` / ``branch`` / ``g_revs``. The branch ``n_rev`` carries the
    Stage-B Lambert ``max_revs`` (base n_rev + the branch's added G revolutions). The
    first element is always the #173 base short-way shape (byte-compatible default),
    so a caller that wants only the historical single branch can take ``[0]``.

    EXPECTED = the SOURCED anchors; every branch's emerged transit is EVIDENCE. The
    branch whose ToF matches the row's real-eph transit signature is the gate (#177).
    """
    base = g_arc_shape(
        aphelion_au, g_tof_years, big_g_tof_years, vinf_e_anchor, vinf_m_anchor, mu=mu
    )
    branches = _g_arc_branch_transits(base.a_au, base.e, max_g_revs=max_g_revs, mu=mu)
    if not branches:
        return [base]
    out: list[GArcShape] = []
    for label, k, tof_days in branches:
        out.append(
            GArcShape(
                a_au=base.a_au,
                e=base.e,
                n_rev=base.n_rev + k,
                tof_g_days=tof_days,
                vinf_e_mag=base.vinf_e_mag,
                vinf_m_mag=base.vinf_m_mag,
                vinf_e_vec_ref=base.vinf_e_vec_ref,
                vinf_e_anchor=base.vinf_e_anchor,
                vinf_m_anchor=base.vinf_m_anchor,
                branch=(f"{label}{k}" if k else label),
                g_revs=k,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Stage B (NEW) — the longitude residual + synodic-period bracket-and-refine scan.
# ---------------------------------------------------------------------------


def _departure_state(shape: GArcShape, ephem: Ephemeris, t_depart_sec: float) -> tuple[Vec3, Vec3]:
    """Heliocentric departure state at ``t_depart`` (the #167 recipe).

    The departure v_inf vector is held fixed in the Earth-orbit ``(r_hat, t_hat,
    z_hat)`` frame and re-pointed at real DE440 Earth's instantaneous velocity at
    ``t_depart``: ``v_sc = v_earth(DE440) + v_inf`` (design §1.2, the same recipe
    :func:`s1l1_corrected.build_seeded_arcs` uses for the App-C seed, but with the
    epoch a free variable rather than a printed constant).
    """
    r_e, v_e = ephem.state("E", t_depart_sec)
    r_e = np.asarray(r_e, dtype=np.float64)
    v_e = np.asarray(v_e, dtype=np.float64)
    # Build the Earth-orbit local frame at the real departure: r_hat outward,
    # t_hat along the velocity-tangential, z_hat the orbit normal.
    r_hat = r_e / np.linalg.norm(r_e)
    h = np.cross(r_e, v_e)
    z_hat = h / np.linalg.norm(h)
    t_hat = np.cross(z_hat, r_hat)
    # The reference vinf vector lives in (+x=r_hat, +y=t_hat, +z=z_hat) at the
    # coplanar crossing; rotate it into the real Earth frame.
    vref = shape.vinf_e_vec_ref
    vinf_world = vref[0] * r_hat + vref[1] * t_hat + vref[2] * z_hat
    v_sc = v_e + vinf_world
    return r_e, np.asarray(v_sc, dtype=np.float64)


def _tof_days(shape: GArcShape, tof_override_days: float | None) -> float:
    """The Lambert/Kepler flight time the Stage-B closer should use (days).

    Returns ``tof_override_days`` when supplied — the row's TABULATED SIGNATURE
    transit (``invariants.transit_times_days``), which is the correct flight time for
    the real DE440 Mars intercept — otherwise falls back to the shape's coplanar
    G-arc branch transit ``shape.tof_g_days`` (historical #173/#177 behaviour,
    preserved so the S1L1 gate and existing self_seeding tests are unchanged).

    ROOT CAUSE of the 2026-06-10 ToF artifact (note
    ``2026-06-10-dsm-tof-artifact-correction.md``): the coplanar branch transit is
    derived from the idealized circular Mars crossing at r = 1.524 AU, but real DE440
    Mars at the rendezvous epoch sits at r ~ 1.40 AU, so the coplanar ToF is the wrong
    transit for the real intercept and inflates the emerged Mars v_inf ~1.6-2.1x.
    Using the row's signature transit collapses that inflation.
    """
    return shape.tof_g_days if tof_override_days is None else float(tof_override_days)


def residual_lon(
    shape: GArcShape,
    ephem: Ephemeris,
    t_depart_sec: float,
    *,
    tof_override_days: float | None = None,
) -> float:
    """Longitude-rendezvous residual at a candidate departure epoch (degrees).

    Propagates the G arc (Kepler-Sun) from real DE440 Earth at ``t_depart`` with the
    #167 departure recipe to ``t_depart + ToF_G`` and returns ``lon_sc_encounter -
    lon_Mars_DE440(t_depart + ToF_G)``, wrapped to (-180, 180]. The single binding
    constraint the App-C seed supplied (design §1.3, the ``residual_lon`` term);
    #165's ~110deg miss was this term omitted. Pure scalar of one variable.

    ``tof_override_days`` (the row's signature transit) overrides the coplanar branch
    ToF when supplied — the 2026-06-10 artifact fix; ``None`` is the historical path.
    """
    r0, v0 = _departure_state(shape, ephem, t_depart_sec)
    tof_s = _tof_days(shape, tof_override_days) * SECONDS_PER_DAY
    r_sc, _ = kepler_propagate(r0, v0, tof_s)
    r_m, _ = ephem.state("M", t_depart_sec + tof_s)
    return _wrap_deg(_lon_deg(np.asarray(r_sc)) - _lon_deg(np.asarray(r_m)))


def synodic_longitude_scan(
    shape: GArcShape,
    ephem: Ephemeris,
    t_center_sec: float,
    *,
    window_days: float | None = None,
    coarse_step_days: float = 10.0,
    refine_tol_sec: float = 3600.0,
    tof_override_days: float | None = None,
) -> list[float]:
    """Scan ``residual_lon`` across one synodic period; return ALL bracketed roots.

    ENUMERATE, do not optimise (design §2.2): sweep ``t_depart`` across one synodic
    window centred on ``t_center`` (default :data:`SYNODIC_PERIOD_DAYS`), coarse step
    (default 10 d), bracket EVERY sign change of :func:`residual_lon`, and
    bisection-refine each bracket to a root. Returns the list of candidate departure
    epochs (seconds since J2000) — typically 1-2 per synodic period where the
    longitudes line up. The scan over a KNOWN periodic structure is the structural
    antidote to the 2026-06-04 free-optimisation off-family failure.

    ``tof_override_days`` (the row's signature transit) is threaded into
    :func:`residual_lon` — the 2026-06-10 artifact fix; ``None`` is the historical path.
    """
    win = (window_days if window_days is not None else SYNODIC_PERIOD_DAYS) * SECONDS_PER_DAY
    step = coarse_step_days * SECONDS_PER_DAY
    t0 = t_center_sec - 0.5 * win
    n = max(2, round(win / step) + 1)

    def f(t: float) -> float:
        return residual_lon(shape, ephem, t, tof_override_days=tof_override_days)

    ts = [t0 + k * step for k in range(n)]
    fs = [f(t) for t in ts]
    roots: list[float] = []
    for k in range(n - 1):
        a, b = ts[k], ts[k + 1]
        fa, fb = fs[k], fs[k + 1]
        # Only a true sign change inside the bracket is a root; skip the wrap jumps
        # of more than 180deg (a -180/+180 longitude wrap, not a rendezvous).
        if fa == 0.0:
            roots.append(a)
            continue
        if (fa < 0.0) != (fb < 0.0) and abs(fa - fb) < 180.0:
            # Bisection that maintains the bracket correctly.
            lo, hi, flo = a, b, fa
            for _ in range(60):
                m = 0.5 * (lo + hi)
                fm = f(m)
                if fm == 0.0 or 0.5 * (hi - lo) < refine_tol_sec:
                    break
                if (flo < 0.0) != (fm < 0.0):
                    hi = m
                else:
                    lo, flo = m, fm
            roots.append(0.5 * (lo + hi))
    return roots


# ---------------------------------------------------------------------------
# Stage B refinement — Lambert onto TRUE DE440 Mars position (longitude automatic).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SelfSeedResult:
    """One self-seeded G-leg candidate (all EVIDENCE; anchors supplied separately).

    Per the orbit-closure discipline every binding constraint is recorded so a
    partial close is diagnosable, not mistaken for a result (design §1.3).
    """

    t_depart_sec: float
    t_arrive_sec: float
    vinf_vec: Vec3  # departure v_inf vector (km/s), refined onto true Mars
    residual_lon_deg: float
    vinf_e_kms: float
    vinf_m_kms: float
    tof_g_days: float
    mars_miss_au: float
    sc_lon_deg: float
    mars_lon_deg: float
    lambert_refined: bool


def _refine_lambert(
    shape: GArcShape,
    ephem: Ephemeris,
    t_depart_sec: float,
    *,
    tof_override_days: float | None = None,
) -> SelfSeedResult:
    """Refine a bracketed epoch with a Lambert solve onto TRUE DE440 Mars position.

    Stage B option 2 (design §1.2): from real Earth at ``t_depart`` to real DE440
    Mars POSITION at ``t_depart + ToF_G``, pick the solution (across the descriptor's
    revolution count) whose departure v_inf magnitude is closest to the shape's Earth
    v_inf. The Lambert arc arrives at Mars's true position BY CONSTRUCTION — longitude
    rendezvous is then automatic; the remaining residual is v_inf-vs-anchor and
    ToF-vs-descriptor (the on-family test). Falls back to the un-refined Kepler arc if
    Lambert is geometrically singular at this epoch.

    ``tof_override_days`` (the row's signature transit) overrides the coplanar branch
    ToF when supplied — the 2026-06-10 artifact fix. The reported ``tof_g_days`` is the
    flight time actually used (so the evidence record is self-consistent).
    """
    tof_days = _tof_days(shape, tof_override_days)
    tof_s = tof_days * SECONDS_PER_DAY
    t_arr = t_depart_sec + tof_s
    r_e, v_e = ephem.state("E", t_depart_sec)
    r_m, v_m = ephem.state("M", t_arr)
    r_e = np.asarray(r_e, dtype=np.float64)
    v_e = np.asarray(v_e, dtype=np.float64)
    r_m = np.asarray(r_m, dtype=np.float64)
    v_m = np.asarray(v_m, dtype=np.float64)

    refined = False
    try:
        sols = lambert(r_e, r_m, tof_s, mu=MU_SUN_KM3_S2, max_revs=max(0, shape.n_rev))
        # Pick the solution whose departure v_inf is closest to the shape anchor.
        best = min(
            sols,
            key=lambda s: abs(float(np.linalg.norm(s.v1 - v_e)) - shape.vinf_e_mag),
        )
        v_dep = np.asarray(best.v1, dtype=np.float64)
        v_arr = np.asarray(best.v2, dtype=np.float64)
        r_sc = r_m
        refined = True
    except (LambertError, ValueError):
        # Singular geometry: report the un-refined Kepler arc (still audited).
        r0, v0 = _departure_state(shape, ephem, t_depart_sec)
        r_sc, v_arr = kepler_propagate(r0, v0, tof_s)
        v_dep = v0

    vinf_vec = v_dep - v_e
    vinf_e = float(np.linalg.norm(vinf_vec))
    vinf_m = float(np.linalg.norm(v_arr - v_m))
    miss_au = float(np.linalg.norm(np.asarray(r_sc) - r_m) / AU_KM)
    return SelfSeedResult(
        t_depart_sec=t_depart_sec,
        t_arrive_sec=t_arr,
        vinf_vec=np.asarray(vinf_vec, dtype=np.float64),
        residual_lon_deg=_wrap_deg(_lon_deg(np.asarray(r_sc)) - _lon_deg(r_m)),
        vinf_e_kms=vinf_e,
        vinf_m_kms=vinf_m,
        tof_g_days=tof_days,
        mars_miss_au=miss_au,
        sc_lon_deg=_lon_deg(np.asarray(r_sc)),
        mars_lon_deg=_lon_deg(r_m),
        lambert_refined=refined,
    )


def self_seed_g_leg(
    shape: GArcShape,
    ephem: Ephemeris,
    t_center_sec: float,
    *,
    window_days: float | None = None,
    coarse_step_days: float = 10.0,
    refine: bool = True,
    tof_override_days: float | None = None,
) -> list[SelfSeedResult]:
    """Full self-seed of one G leg: scan + (optional) Lambert refine per candidate.

    Runs :func:`synodic_longitude_scan` to enumerate every synodic-phase candidate
    epoch, then (``refine``) refines each onto true DE440 Mars with a Lambert solve.
    Returns ALL candidates as :class:`SelfSeedResult` evidence (design §2.2 — report
    every on-family-looking root, never force a single pick).

    ``tof_override_days`` (the row's signature transit) is used as the Lambert/Kepler
    flight time when supplied — the 2026-06-10 artifact fix; ``None`` is the historical
    coplanar-branch-ToF path (preserves the S1L1 gate and existing tests).
    """
    roots = synodic_longitude_scan(
        shape,
        ephem,
        t_center_sec,
        window_days=window_days,
        coarse_step_days=coarse_step_days,
        tof_override_days=tof_override_days,
    )
    out: list[SelfSeedResult] = []
    for t in roots:
        if refine:
            out.append(_refine_lambert(shape, ephem, t, tof_override_days=tof_override_days))
        else:
            r0, v0 = _departure_state(shape, ephem, t)
            tof_s = _tof_days(shape, tof_override_days) * SECONDS_PER_DAY
            r_sc, v_arr = kepler_propagate(r0, v0, tof_s)
            r_m, v_m = ephem.state("M", t + tof_s)
            r_m = np.asarray(r_m, dtype=np.float64)
            v_m = np.asarray(v_m, dtype=np.float64)
            _r_e, v_e = ephem.state("E", t)
            out.append(
                SelfSeedResult(
                    t_depart_sec=t,
                    t_arrive_sec=t + tof_s,
                    vinf_vec=np.asarray(v0 - np.asarray(v_e), dtype=np.float64),
                    residual_lon_deg=_wrap_deg(_lon_deg(np.asarray(r_sc)) - _lon_deg(r_m)),
                    vinf_e_kms=shape.vinf_e_mag,
                    vinf_m_kms=float(np.linalg.norm(np.asarray(v_arr) - v_m)),
                    tof_g_days=_tof_days(shape, tof_override_days),
                    mars_miss_au=float(np.linalg.norm(np.asarray(r_sc) - r_m) / AU_KM),
                    sc_lon_deg=_lon_deg(np.asarray(r_sc)),
                    mars_lon_deg=_lon_deg(r_m),
                    lambert_refined=False,
                )
            )
    return out


# ---------------------------------------------------------------------------
# Stage B (2026-06-10 ToF-artifact fix) — joint (epoch, ToF) free-variable closer.
#
# The signature-ToF override (`tof_override_days`) is a NECESSARY but not sufficient
# fix: forcing the coplanar branch ToF onto a longitude-rendezvous epoch is what
# inflated the emerged Mars v_inf (note 2026-06-10-dsm-tof-artifact-correction). But
# pinning the epoch from a longitude scan at ONE fixed ToF still over-constrains the
# real DE440 intercept (the longitude root and the low-v_inf intercept need not
# coincide at exactly the signature ToF). The note's option 2 — open BOTH the
# departure epoch and the ToF as free variables, bracketed near the signature, and
# select the (epoch, ToF) whose departure AND arrival v_inf both fall in the row's
# anchor band — closes all examined rows to 0.1-0.4 km/s. Implemented here.
#
# Honesty: the anchors are used ONLY to SELECT among physically-enumerated Lambert
# solutions (the same disambiguation `_refine_lambert` already does with vinf_e_mag);
# the chosen solution's v_inf is EMERGED and REPORTED, then compared against the
# anchor by the unchanged `on_family` gate (band NEVER loosened). The arrival is the
# true DE440 Mars POSITION by Lambert construction, so the longitude rendezvous and
# Mars miss are exact by construction; the residual scientific term is v_inf-vs-anchor.
# ---------------------------------------------------------------------------


def joint_epoch_tof_close(
    ephem: Ephemeris,
    anchors: FamilyAnchors,
    t_center_sec: float,
    signature_tof_days: float,
    *,
    epoch_halfwidth_days: float = 850.0,
    epoch_step_days: float = 10.0,
    tof_halfwidth_days: float = 40.0,
    tof_step_days: float = 5.0,
    max_revs: int = 0,
    refine_iters: int = 2,
) -> SelfSeedResult | None:
    """Close one G leg with a free (epoch, ToF) Lambert search onto true DE440 Mars.

    Sweeps the departure epoch over ``+-epoch_halfwidth_days`` of ``t_center`` and the
    ToF over ``+-tof_halfwidth_days`` of the row's ``signature_tof_days``; for every
    (epoch, ToF) solves the Lambert from real Earth to real DE440 Mars POSITION at
    ``epoch + ToF`` and scores each solution by ``|vinf_E - anchor_E| + |vinf_M -
    anchor_M|`` (anchor-band SELECTION over a physical enumeration — the emerged v_inf
    is still reported). Returns the best (lowest combined anchor error) as a
    :class:`SelfSeedResult`, then locally refines the grid ``refine_iters`` times to
    sharpen the (epoch, ToF). Returns ``None`` if no Lambert solution exists anywhere
    in the bracket (a clean negative). Longitude rendezvous + Mars miss are exact by
    construction (arrival = true Mars position).
    """
    best: tuple[float, float, float, SelfSeedResult] | None = None  # (err, depart, tof, result)
    e_half, e_step = epoch_halfwidth_days, epoch_step_days
    t_half, t_step = tof_halfwidth_days, tof_step_days
    e_center, t_center_tof = t_center_sec, signature_tof_days

    for _it in range(refine_iters + 1):
        de = -e_half
        while de <= e_half + 1e-9:
            t_depart = e_center + de * SECONDS_PER_DAY
            r_e, v_e = ephem.state("E", t_depart)
            r_e = np.asarray(r_e, dtype=np.float64)
            v_e = np.asarray(v_e, dtype=np.float64)
            tof_d = max(1.0, t_center_tof - t_half)
            tof_hi = t_center_tof + t_half
            while tof_d <= tof_hi + 1e-9:
                t_arr = t_depart + tof_d * SECONDS_PER_DAY
                r_m, v_m = ephem.state("M", t_arr)
                r_m = np.asarray(r_m, dtype=np.float64)
                v_m = np.asarray(v_m, dtype=np.float64)
                try:
                    sols = lambert(
                        r_e,
                        r_m,
                        tof_d * SECONDS_PER_DAY,
                        mu=MU_SUN_KM3_S2,
                        max_revs=max(0, max_revs),
                    )
                except (LambertError, ValueError):
                    tof_d += t_step
                    continue
                for s in sols:
                    v_dep = np.asarray(s.v1, dtype=np.float64)
                    v_arr = np.asarray(s.v2, dtype=np.float64)
                    vinf_e = float(np.linalg.norm(v_dep - v_e))
                    vinf_m = float(np.linalg.norm(v_arr - v_m))
                    err = abs(vinf_e - anchors.vinf_e) + abs(vinf_m - anchors.vinf_m)
                    if best is None or err < best[0]:
                        vinf_vec = v_dep - v_e
                        res = SelfSeedResult(
                            t_depart_sec=t_depart,
                            t_arrive_sec=t_arr,
                            vinf_vec=np.asarray(vinf_vec, dtype=np.float64),
                            residual_lon_deg=_wrap_deg(_lon_deg(r_m) - _lon_deg(r_m)),
                            vinf_e_kms=vinf_e,
                            vinf_m_kms=vinf_m,
                            tof_g_days=float(tof_d),
                            mars_miss_au=0.0,
                            sc_lon_deg=_lon_deg(r_m),
                            mars_lon_deg=_lon_deg(r_m),
                            lambert_refined=True,
                        )
                        best = (err, de, tof_d, res)
                tof_d += t_step
            de += e_step
        if best is None:
            return None
        # Local refinement: re-centre and shrink the grid around the current best.
        _err, de_b, tof_b, _res = best
        e_center = e_center + de_b * SECONDS_PER_DAY
        t_center_tof = tof_b
        e_half, e_step = 2.0 * e_step, e_step / 4.0
        t_half, t_step = 2.0 * t_step, t_step / 4.0

    if best is None:
        return None
    return best[3]


# ---------------------------------------------------------------------------
# The on-family multi-term gate at a single epoch (design §1.3).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OnFamilyVerdict:
    """Structured per-term verdict of the on-family gate (NOT a bare bool).

    Every binding constraint is reported so a partial close is diagnosable rather
    than mistaken for a result ("it closed!" on a subset is the danger signal).
    """

    lon_ok: bool
    vinf_e_ok: bool
    vinf_m_ok: bool
    miss_ok: bool
    residual_lon_deg: float
    vinf_e_kms: float
    vinf_m_kms: float
    mars_miss_au: float

    @property
    def on_family(self) -> bool:
        """ALL terms small at the SAME epoch (the full on-family test)."""
        return self.lon_ok and self.vinf_e_ok and self.vinf_m_ok and self.miss_ok


@dataclass(frozen=True)
class FamilyAnchors:
    """The EXPECTED anchors the on-family gate tests a result against.

    Supplied SEPARATELY from the search (the search is anchor-blind on the epoch
    axis): the row's real-eph v_inf anchors (a BAND wide enough to admit breathing,
    design §2.3 inverse-#164 trap) and the encounter miss band.
    """

    vinf_e: float
    vinf_m: float
    vinf_band_kms: float = 1.0
    lon_tol_deg: float = 0.5
    miss_tol_au: float = 0.0116  # 3 Mars SOI — the #165/#167 band, NEVER loosened


def on_family(
    result: SelfSeedResult,
    anchors: FamilyAnchors,
) -> OnFamilyVerdict:
    """Multi-term on-family predicate at a single epoch (design §1.3).

    ON-FAMILY iff ALL hold simultaneously: ``|residual_lon|`` small, departure /
    arrival v_inf within a breathing band of the anchors, and the Mars miss inside
    the encounter band. Returns a structured verdict (which terms pass / fail), NOT a
    bare bool — so a partial close is diagnosable. The v_inf band is wide on purpose
    (the real-eph v_inf breathes; an unsourced row's band is uncertain — design §2.3).
    """
    lon_ok = abs(result.residual_lon_deg) < anchors.lon_tol_deg
    vinf_e_ok = abs(result.vinf_e_kms - anchors.vinf_e) < anchors.vinf_band_kms
    vinf_m_ok = abs(result.vinf_m_kms - anchors.vinf_m) < anchors.vinf_band_kms
    miss_ok = result.mars_miss_au < anchors.miss_tol_au
    return OnFamilyVerdict(
        lon_ok=lon_ok,
        vinf_e_ok=vinf_e_ok,
        vinf_m_ok=vinf_m_ok,
        miss_ok=miss_ok,
        residual_lon_deg=result.residual_lon_deg,
        vinf_e_kms=result.vinf_e_kms,
        vinf_m_kms=result.vinf_m_kms,
        mars_miss_au=result.mars_miss_au,
    )


@dataclass(frozen=True)
class TransitTriage:
    """Cheap REACHABLE-vs-OFF-FAMILY triage of a row by branch-transit match (#177).

    The #173 gating condition: a row is REACHABLE only if SOME coplanar G-arc branch
    transit lands within ``tol_days`` of the row's real-eph transit signature. This is
    the cheap per-row test (no n-body) the bulk triage runs. ``best_branch`` /
    ``best_tof_days`` / ``delta_days`` record the closest branch as evidence.
    """

    reachable: bool
    real_eph_transit_days: float
    best_branch: str
    best_g_revs: int
    best_tof_days: float
    delta_days: float
    branch_tofs: dict[str, float]


def triage_transit_match(
    branches: list[GArcShape],
    real_eph_transit_days: float,
    *,
    tol_days: float = 30.0,
) -> TransitTriage:
    """Classify a row REACHABLE / OFF-FAMILY by the #173 branch-transit gate (#177).

    REACHABLE iff the closest enumerated G-arc branch transit (:func:`g_arc_branches`)
    is within ``tol_days`` of the row's real-eph transit signature. The tolerance is
    the coplanar->real-eph transit gap admitted (S1L1: coplanar 188 vs real-eph 180,
    ~8 d; the long-transit families breathe more) — it is the GATE, NOT loosened to
    inflate REACHABLE (brief honesty rule). Pure: no ephemeris, no n-body.
    """
    branch_tofs = {b.branch: b.tof_g_days for b in branches}
    best = min(branches, key=lambda b: abs(b.tof_g_days - real_eph_transit_days))
    delta = best.tof_g_days - real_eph_transit_days
    return TransitTriage(
        reachable=abs(delta) <= tol_days,
        real_eph_transit_days=float(real_eph_transit_days),
        best_branch=best.branch,
        best_g_revs=best.g_revs,
        best_tof_days=float(best.tof_g_days),
        delta_days=float(delta),
        branch_tofs={k: float(v) for k, v in branch_tofs.items()},
    )


__all__ = [
    "SYNODIC_PERIOD_DAYS",
    "FamilyAnchors",
    "GArcShape",
    "OnFamilyVerdict",
    "SelfSeedResult",
    "TransitTriage",
    "g_arc_branches",
    "g_arc_shape",
    "joint_epoch_tof_close",
    "on_family",
    "residual_lon",
    "self_seed_g_leg",
    "synodic_longitude_scan",
    "triage_transit_match",
]
