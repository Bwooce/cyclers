"""Two-arc free-return CHAIN — Russell's generic-return-arc construction (#163).

The single-arc free-return corrector (:mod:`cyclerfinder.search.free_return`) closes
ONE heliocentric ellipse whose radial crossings reproduce a sourced Russell cycler's
DERIVED V_inf. That primitive CLOSES the eight symmetric single-generic rows because
its objective is anchor-respecting: ``(a, e)`` is seeded from the sourced aphelion +
transit, and the per-body V_inf EMERGES and is checked against the anchor (the #137
breakthrough). But it is REFUSED on the genuinely *multi-arc* rows
(``test_free_return_v1_mechanics.py::NON_V1_ROWS``): a single ellipse cannot
represent TWO distinct generic-return arcs, so its forced single-ellipse return
breaks Mars V_inf continuity by ~tens of km/s.

This module builds the right primitive for those rows: Russell's actual
construction (2004 Ch.4; ``docs/notes/multi-arc-classification.md`` §1.2-§1.3, §5).
A multi-arc cycler such as ``russell-ch4-6.44Gg3`` (descriptor ``g(2.087,...) +
G(4.3191,...)``) is a *sequence of two distinct Earth-to-Earth free-return arcs*
patched by an intermediate Earth gravity-assist flyby. Each arc is its own
free-return ellipse ``(a_i, e_i)`` that crosses Mars's radius (the Mars encounter)
on the way out/back; the two arcs share the Earth departure/return.

What makes the two arcs DISTINCT (and prevents single-ellipse collapse)
----------------------------------------------------------------------
Each Russell leg descriptor is ``g(TOF_years, psi_deg, branch)``: the FIRST number
is the full Earth-to-Earth ARC time in YEARS, the SECOND is the transfer/wrap angle
``psi`` on the V_inf sphere (it EXCEEDS 360 deg — ``1111.33`` for the g-arc,
``1194.88`` for the G-arc of 6.44Gg3 — so the arc is MULTI-revolution). The g-arc
(2.087 yr) and the G-arc (4.3191 yr) have DIFFERENT ToFs, so they cannot be the same
ellipse. The chain therefore binds each arc to THREE conditions on its 2 continuous
DOF ``(a, e)`` plus an integer revolution count ``n_rev``:

* emerged V_inf at Mars  == sourced Mars anchor,
* emerged V_inf at Earth == sourced Earth anchor,
* emerged Earth-to-Earth arc ToF (= ``n_rev * period + time-above-Earth``) ==
  the descriptor ToF.

Three real conditions on two continuous DOF is OVERdetermined per arc — which is
exactly why a SINGLE arc is refused (#137) and why the chain can still legitimately
come up EMPTY-SET. The integer ``n_rev`` supplies the discrete freedom that brings
the (long, multi-rev) descriptor ToF into range; the chain picks, per arc, the
``n_rev`` that minimises the ToF residual at the current ``(a, e)``. Without the ToF
term the residual would collapse to a single ellipse that hits the V_inf anchors at
the WRONG ToF (a spurious CLOSE); the ToF term is what keeps the two arcs honest.

Why this is NOT a repeat of the dsm_leg floor (#162)
----------------------------------------------------
The chained-DSM evaluator (:mod:`cyclerfinder.search.dsm_leg`) floored at ~9 km/s
because its objective is ``Sum dV_DSM``-min, which does NOT pin V_inf to the anchor
(``docs/notes/2026-06-08-multiarc-basin-selection-results.md``). It optimises a
*budget*, so the emerged V_inf wanders off-anchor. The wrong primitive. This module
reuses the free-return primitive's anchor-respecting objective: the residual IS
``(emerged V_inf - sourced anchor)`` at every encounter. V_inf is never imposed.

Constraint-vs-evidence separation (the golden rule, inherited)
--------------------------------------------------------------
* CONSTRAINED (drives the residual): per-arc emerged V_inf at Earth and Mars match
  the SOURCED anchors; per-arc emerged arc ToF matches the SOURCED descriptor ToF;
  the intermediate Earth flyby is V_inf-continuous and bend-feasible.
* FREE (the unknowns): the two arc shapes ``(a_1, e_1, a_2, e_2)`` (continuous) and
  the per-arc revolution counts (discrete, chosen by ToF-min).
* EVIDENCE (emerges, compared non-circularly): the converged per-arc V_inf
  magnitudes at Earth and Mars and the per-arc arc ToFs.

The SOURCED anchor is EXPECTED; the emerged V_inf is evidence. A clean EMPTY-SET
(Russell's own primitive cannot hit the anchors in circular-coplanar) is the
strongest publishable negative, NOT a failure to soften.

Pure: depends only on core/constants, core/flyby, core/kepler, search/free_return.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SAFE_PERIHELION_KM,
)
from cyclerfinder.core.flyby import is_ballistic_feasible, max_bend
from cyclerfinder.core.kepler import coe_to_rv
from cyclerfinder.search.free_return import (
    FreeReturnGeometry,
    _crossing,
    _true_to_mean,
    free_return_geometry,
    seed_ae_from_aphelion_transit,
)

DAYS_PER_JULIAN_YEAR = 365.25
SECONDS_PER_YEAR = DAYS_PER_JULIAN_YEAR * 86400.0
# Revolution counts the per-arc ToF search ranges over. The descriptor ToFs of the
# target rows (1.46-4.32 yr) against ~1.3 yr free-return periods need n_rev up to ~4.
_N_REV_RANGE: tuple[int, ...] = (0, 1, 2, 3, 4, 5)


def _arc_ee_time_years(a_au: float, e: float, n_rev: int, *, mu: float = MU_SUN_KM3_S2) -> float:
    """Earth-to-Earth free-return arc time (years) for ``(a, e)`` and ``n_rev``.

    The free-return spacecraft departs Earth at the outbound radial crossing
    ``nu_E`` (in ``[0, pi]``), climbs past Mars to aphelion and returns to the
    inbound crossing ``2pi - nu_E``; the time *above Earth radius* is
    ``period - 2 * (M(nu_E) / n)`` where ``M`` is the mean anomaly. A MULTI-rev arc
    adds ``n_rev`` full periods (the descriptor ToFs of 2-4 yr against a ~1.3 yr
    free-return period are multi-rev — the ``psi > 360 deg`` wrap in the descriptor).
    """
    a_km = a_au * AU_KM
    nu_e = _crossing(a_km, e, PLANETS["E"].sma_au * AU_KM)
    n = sqrt(mu / a_km**3)
    period = 2.0 * pi / n
    m_e = _true_to_mean(nu_e, e)
    t_above = period - 2.0 * (m_e / n)
    return (t_above + n_rev * period) / SECONDS_PER_YEAR


def _best_n_rev(
    a_au: float, e: float, tof_target_years: float, *, mu: float = MU_SUN_KM3_S2
) -> int:
    """Integer ``n_rev`` minimising ``|arc_ee_time - tof_target|`` (discrete DOF)."""
    return min(
        _N_REV_RANGE,
        key=lambda n: abs(_arc_ee_time_years(a_au, e, n, mu=mu) - tof_target_years),
    )


@dataclass(frozen=True)
class FreeReturnArc:
    """One converged free-return arc of the chain (shape + emerged geometry)."""

    a_au: float
    e: float
    geometry: FreeReturnGeometry
    n_rev: int
    arc_tof_years: float

    @property
    def vinf_e(self) -> float:
        """EMERGED Earth (departure/return) V_inf magnitude, km/s."""
        return self.geometry.vinf["E"]

    @property
    def vinf_m(self) -> float:
        """EMERGED Mars (crossing) V_inf magnitude, km/s."""
        return self.geometry.vinf["M"]

    @property
    def transfer_tof_days(self) -> float:
        """EMERGED inner->outer (E->M) radial-crossing leg ToF, days."""
        return self.geometry.tof_em_days


@dataclass(frozen=True)
class FreeReturnChainResult:
    """Outcome of :func:`free_return_chain_correct`.

    Attributes
    ----------
    arcs:
        The two converged free-return arcs (``arc-1`` = g, ``arc-2`` = G).
    max_residual_kms:
        Max absolute element of the anchor-respecting residual vector (km/s):
        per-arc (emerged V_inf - sourced anchor) at Earth and Mars, per-arc
        (emerged arc ToF - descriptor ToF) scaled to km/s by Earth's orbital speed,
        plus the intermediate-Earth-flyby V_inf-magnitude continuity. The
        AUTHORITATIVE acceptance number (residual-magnitude only, by design — the
        residual IS the physics, exactly as in :mod:`free_return`).
    vinf_residual_kms:
        DIAGNOSTIC: the worst per-arc V_inf-anchor mismatch alone (the headline the
        three-way gate's "within 0.5 of both anchors" condition reads), separated
        from the ToF term.
    tof_residual_years:
        DIAGNOSTIC: the worst per-arc arc-ToF mismatch (years).
    vinf_continuity_kms:
        DIAGNOSTIC: ``| |V_inf_E arc1| - |V_inf_E arc2| |`` at the intermediate
        Earth flyby — the magnitude-continuity term, separated for the audit trail.
    intermediate_flyby_feasible:
        Whether the intermediate Earth flyby is BEND-feasible at the converged
        geometry (magnitudes match AND required turn within the achievable cone at
        the safe Earth periapsis), reusing
        :func:`cyclerfinder.core.flyby.is_ballistic_feasible`.
    intermediate_turn_deg, intermediate_max_turn_deg:
        DIAGNOSTIC: required and maximum-achievable intermediate Earth turn (deg).
    converged:
        ``max_residual_kms < tol_kms`` — residual-magnitude only, BY DESIGN.
    vinf_within_tol:
        ``vinf_residual_kms < vinf_tol_kms`` — the emerged-V_inf-vs-anchor half of
        the three-way gate, evaluated independently of the ToF term.
    bend_feasible_close:
        ``converged AND intermediate_flyby_feasible`` — the full CLOSE gate. The
        single boolean the probe keys on for the FIRST-multi-arc-closure verdict.
    solver_success, solver_nfev:
        DIAGNOSTIC: underlying ``least_squares`` outcome (audit trail only).
    """

    arcs: tuple[FreeReturnArc, FreeReturnArc]
    max_residual_kms: float
    vinf_residual_kms: float
    tof_residual_years: float
    vinf_continuity_kms: float
    intermediate_flyby_feasible: bool
    intermediate_turn_deg: float
    intermediate_max_turn_deg: float
    converged: bool
    vinf_within_tol: bool
    bend_feasible_close: bool
    solver_success: bool = True
    solver_nfev: int = 0


def _earth_vinf_vector(arc_a_au: float, arc_e: float, *, mu: float = MU_SUN_KM3_S2) -> np.ndarray:
    """Heliocentric ``V_inf`` vector at the outbound Earth radial crossing.

    Same radial/tangential reconstruction as :func:`free_return_geometry`: the
    spacecraft velocity from :func:`coe_to_rv` at the Earth crossing true anomaly
    minus the circular Earth velocity. The direction lets the intermediate-flyby
    bend be measured between the two arcs' Earth asymptotes; the magnitude equals
    the arc's emerged ``vinf_e`` by construction.
    """
    a_km = arc_a_au * AU_KM
    nu_e = _crossing(a_km, arc_e, PLANETS["E"].sma_au * AU_KM)
    r_sc, v_sc = coe_to_rv(a_km, arc_e, nu_e, mu)
    r_hat = r_sc / np.linalg.norm(r_sc)
    t_hat = np.array([-r_hat[1], r_hat[0], 0.0])
    v_planet = sqrt(mu / (PLANETS["E"].sma_au * AU_KM)) * t_hat
    return np.asarray(v_sc - v_planet, dtype=np.float64)


def _intermediate_turn_geometry(
    arc1: FreeReturnArc,
    arc2: FreeReturnArc,
    *,
    body: str = "E",
    mu: float = MU_SUN_KM3_S2,
    continuity_tol_kms: float = 0.5,
) -> tuple[float, float, bool, float]:
    """Required vs achievable BEND + V_inf continuity at the intermediate Earth flyby.

    Returns ``(turn_rad, max_turn_rad, bend_feasible, continuity_kms)``.

    The two are kept SEPARATE on purpose. ``continuity_kms = | |V_inf_in| -
    |V_inf_out| |`` is the V_inf-magnitude continuity — but that is the SAME physics
    as the two Earth-anchor residual terms (both arcs pull to the Earth anchor), so
    it is already counted in :func:`_chain_residuals` and is NOT re-charged here.
    ``bend_feasible`` is the genuinely additional patch-feasibility gate the task
    asks for: with the magnitudes matched to within ``continuity_tol_kms`` (the
    closure precision, NOT a hardcoded 1 m/s), is the required asymptote rotation
    inside the achievable ballistic cone at the safe Earth periapsis? Computed via
    :func:`max_bend` on the mean V_inf and compared to the turn between the two
    arcs' Earth asymptotes; :func:`is_ballistic_feasible` is reused for the in-cone
    test at the matched magnitude so the bend convention is shared with core/flyby.
    """
    v_in = _earth_vinf_vector(arc1.a_au, arc1.e, mu=mu)
    v_out = _earth_vinf_vector(arc2.a_au, arc2.e, mu=mu)
    vin_mag = float(np.linalg.norm(v_in))
    vout_mag = float(np.linalg.norm(v_out))
    continuity = abs(vin_mag - vout_mag)
    mu_planet = PLANETS[body].mu_km3_s2
    rp_min = SAFE_PERIHELION_KM[body]
    cos_arg = float(np.dot(v_in, v_out)) / (vin_mag * vout_mag)
    turn = float(np.arccos(np.clip(cos_arg, -1.0, 1.0)))
    v_mean = 0.5 * (vin_mag + vout_mag)
    max_turn = max_bend(mu_planet, rp_min, v_mean)
    # Bend feasibility at the matched magnitude (the continuity term is handled by
    # the residual; here we ask only whether the turn fits the cone). Reuse
    # is_ballistic_feasible with the magnitudes equalised to the mean so its
    # speed_tol never short-circuits the in-cone judgement, and require the
    # magnitudes to actually be continuous to the closure precision.
    v_in_eq = v_in / vin_mag * v_mean
    v_out_eq = v_out / vout_mag * v_mean
    in_cone = is_ballistic_feasible(v_in_eq, v_out_eq, mu_planet, rp_min, speed_tol=1.0e-6)
    bend_feasible = bool(in_cone and continuity < continuity_tol_kms)
    return turn, max_turn, bend_feasible, continuity


def _build_arc(
    a_au: float, e: float, tof_target_years: float, *, mu: float = MU_SUN_KM3_S2
) -> FreeReturnArc | None:
    """Build one arc's emerged geometry (+ ToF-min n_rev), or ``None`` off-family."""
    try:
        g = free_return_geometry(a_au, e, bodies=("E", "M"), mu=mu)
    except ValueError:
        return None
    n_rev = _best_n_rev(a_au, e, tof_target_years, mu=mu)
    arc_tof = _arc_ee_time_years(a_au, e, n_rev, mu=mu)
    return FreeReturnArc(a_au=a_au, e=e, geometry=g, n_rev=n_rev, arc_tof_years=arc_tof)


# km/s the residual charges per YEAR of arc-ToF mismatch (documented weight). The
# descriptor ToF is a hard STRUCTURAL constraint (it is what distinguishes the g- and
# G-arcs and forbids single-ellipse collapse), so it is weighted on a par with the
# V_inf terms: 1 km/s per year means a 0.5 km/s closure tol admits a 0.5 yr (~half-
# rev) ToF slack, which the discrete n_rev cannot finesse away. Chosen, not sourced;
# the per-term residuals are reported separately so the gate reads V_inf and ToF
# independently and the weight never manufactures a CLOSE.
_TOF_WEIGHT_KMS_PER_YEAR = 1.0


def _chain_residuals(
    x: np.ndarray,
    *,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    arc1_tof_years: float,
    arc2_tof_years: float,
    mu: float,
) -> list[float]:
    """Anchor-respecting residual for the two-arc chain (km/s-comparable).

    Free variables ``x = [a1, e1, a2, e2]``. Per arc: emerged V_inf at Earth and
    Mars must match the SOURCED anchors, and the emerged Earth-to-Earth arc ToF (at
    the ToF-min integer ``n_rev``) must match the descriptor ToF. The ToF term is
    weighted to km/s by :data:`_TOF_WEIGHT_KMS_PER_YEAR` so it is commensurate with
    the V_inf terms in ``max|res|``. Nothing imposes the anchor V_inf; it EMERGES
    from each ``(a_i, e_i)``. The intermediate-flyby continuity is evaluated
    separately (it is a function of the same two Earth-anchor terms).
    """
    a1, e1, a2, e2 = float(x[0]), float(x[1]), float(x[2]), float(x[3])
    for a, e in ((a1, e1), (a2, e2)):
        if not (0.0 < e < 0.95) or a <= 0.0:
            return [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]
    arc1 = _build_arc(a1, e1, arc1_tof_years, mu=mu)
    arc2 = _build_arc(a2, e2, arc2_tof_years, mu=mu)
    if arc1 is None or arc2 is None:
        return [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]
    w = _TOF_WEIGHT_KMS_PER_YEAR
    return [
        arc1.vinf_m - vinf_m_anchor,
        arc2.vinf_m - vinf_m_anchor,
        arc1.vinf_e - vinf_e_anchor,
        arc2.vinf_e - vinf_e_anchor,
        w * (arc1.arc_tof_years - arc1_tof_years),
        w * (arc2.arc_tof_years - arc2_tof_years),
    ]


def _seed_transit_days(arc_tof_years: float) -> float:
    """Clamp an arc descriptor ToF (years) to a reachable E->M transit seed (days).

    The descriptor ToF is the full (multi-rev) Earth-to-Earth arc time; the radial-
    crossing E->M transit that :func:`seed_ae_from_aphelion_transit` solves for is a
    much shorter fraction. We only need a transit value inside the bisection's
    reachable band to fix a STARTING ``e`` per arc; the solver then relaxes both
    ``(a, e)`` freely. Map the within-year fraction onto ``[120, 320]`` d — the band
    the symmetric Russell rows' E->M legs occupy (150 d for 4.991gG2, 262 d for
    6.44Gg3) — so the two descriptors seed at slightly different ``e``.
    """
    lo, hi = 120.0, 320.0
    frac = arc_tof_years - int(arc_tof_years)
    return lo + frac * (hi - lo)


def free_return_chain_correct(
    aphelion_au: float,
    arc1_tof_years: float,
    arc2_tof_years: float,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = 0.5,
    vinf_tol_kms: float = 0.5,
    seeds: tuple[tuple[float, float], tuple[float, float]] | None = None,
) -> FreeReturnChainResult:
    """Close a TWO-ARC free-return chain at the SOURCED anchors.

    Russell's generic-return construction: two distinct Earth-to-Earth free-return
    arcs (``arc-1`` = g, ``arc-2`` = G), each crossing Mars's radius, patched at an
    intermediate Earth gravity-assist flyby. Free variables are the two arc shapes
    ``(a_1, e_1, a_2, e_2)``; the per-arc V_inf at Earth and Mars and the per-arc
    Earth-to-Earth arc ToF (at the ToF-min ``n_rev``) EMERGE and are driven to the
    SOURCED anchors / descriptor ToFs. The intermediate Earth flyby is then checked
    for bend feasibility.

    Anchor-respecting (the residual IS the anchor + ToF match), in deliberate
    contrast to the dV-budget objective of :mod:`dsm_leg`. The sourced anchors are
    EXPECTED; the emerged V_inf is evidence.

    Parameters
    ----------
    aphelion_au:
        The cycler's SOURCED aphelion (catalogue ``orbit_elements.aphelion_au``),
        used to seed each arc. A CONSTRAINT; the V_inf is evidence.
    arc1_tof_years, arc2_tof_years:
        The two free-return arc ToFs in YEARS (the FIRST number of each Russell leg
        descriptor). The structural binding that distinguishes g from G.
    vinf_e_anchor, vinf_m_anchor:
        The SOURCED V_inf anchors at Earth and Mars (km/s) — the EXPECTED target.
    tol_kms:
        Closure tolerance on ``max_residual_kms`` (default 0.5 km/s — the campaign
        anchor-match band, applied to the COMBINED V_inf+ToF residual).
    vinf_tol_kms:
        Tolerance for the V_inf-only ``vinf_within_tol`` flag (default 0.5 km/s — the
        three-way gate's "within 0.5 of both anchors" condition).
    seeds:
        Optional explicit per-arc ``((a1, e1), (a2, e2))`` seed override (MBH / multi-
        start). When ``None`` each arc is seeded from ``aphelion_au`` + its ToF.
    """
    if seeds is None:
        a1_seed, e1_seed = seed_ae_from_aphelion_transit(
            aphelion_au, _seed_transit_days(arc1_tof_years), bodies=("E", "M"), mu=mu
        )
        a2_seed, e2_seed = seed_ae_from_aphelion_transit(
            aphelion_au, _seed_transit_days(arc2_tof_years), bodies=("E", "M"), mu=mu
        )
    else:
        (a1_seed, e1_seed), (a2_seed, e2_seed) = seeds

    def _res(x: np.ndarray) -> list[float]:
        return _chain_residuals(
            x,
            vinf_e_anchor=vinf_e_anchor,
            vinf_m_anchor=vinf_m_anchor,
            arc1_tof_years=arc1_tof_years,
            arc2_tof_years=arc2_tof_years,
            mu=mu,
        )

    x0 = np.array([a1_seed, e1_seed, a2_seed, e2_seed], dtype=np.float64)
    sol = least_squares(_res, x0, method="trf", max_nfev=400, xtol=1e-12, ftol=1e-12)
    x = sol.x
    res = _res(x)
    max_res = max(abs(r) for r in res)

    arc1 = _build_arc(float(x[0]), float(x[1]), arc1_tof_years, mu=mu)
    arc2 = _build_arc(float(x[2]), float(x[3]), arc2_tof_years, mu=mu)
    if arc1 is None or arc2 is None:
        # Solver wandered off the Mars-reaching family. Report a non-close with a
        # guaranteed-reachable audit arc (a known Mars-crossing ellipse) so the
        # result object is always well-formed even for off-family seeds.
        audit = _build_arc(1.30, 0.25, arc1_tof_years, mu=mu)
        assert audit is not None  # (1.30, 0.25) reaches Mars by construction
        return FreeReturnChainResult(
            arcs=(audit, audit),
            max_residual_kms=float("inf"),
            vinf_residual_kms=float("inf"),
            tof_residual_years=float("inf"),
            vinf_continuity_kms=float("inf"),
            intermediate_flyby_feasible=False,
            intermediate_turn_deg=float("nan"),
            intermediate_max_turn_deg=float("nan"),
            converged=False,
            vinf_within_tol=False,
            bend_feasible_close=False,
            solver_success=bool(sol.success),
            solver_nfev=int(sol.nfev),
        )

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
    turn, max_turn, feasible, continuity = _intermediate_turn_geometry(
        arc1, arc2, mu=mu, continuity_tol_kms=tol_kms
    )
    converged = max_res < tol_kms
    vinf_within = vinf_res < vinf_tol_kms
    return FreeReturnChainResult(
        arcs=(arc1, arc2),
        max_residual_kms=float(max_res),
        vinf_residual_kms=float(vinf_res),
        tof_residual_years=float(tof_res),
        vinf_continuity_kms=float(continuity),
        intermediate_flyby_feasible=bool(feasible),
        intermediate_turn_deg=float(np.degrees(turn)),
        intermediate_max_turn_deg=float(np.degrees(max_turn)),
        converged=converged,
        vinf_within_tol=vinf_within,
        bend_feasible_close=bool(converged and feasible),
        solver_success=bool(sol.success),
        solver_nfev=int(sol.nfev),
    )


def single_arc_degenerate(
    aphelion_au: float,
    arc_tof_years: float,
    vinf_e_anchor: float,
    vinf_m_anchor: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = 0.5,
) -> FreeReturnChainResult:
    """Degenerate single-arc case: both chain arcs equal -> reduces to free_return.

    Drives a chain with arc1 == arc2 (same descriptor ToF, same anchors). The
    intermediate-flyby continuity term is then exactly zero by construction and the
    two arcs converge to the SAME ``(a, e)`` — the mechanics gate that the chain
    reduces to a single free-return ellipse when the two arcs coincide.
    """
    return free_return_chain_correct(
        aphelion_au,
        arc_tof_years,
        arc_tof_years,
        vinf_e_anchor,
        vinf_m_anchor,
        mu=mu,
        tol_kms=tol_kms,
    )
