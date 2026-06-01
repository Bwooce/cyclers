"""Per-cycler ranking metrics — :class:`Score`, ``score()``, ``rank()`` (M4).

Spec references
---------------
* §5 step 4 — "rank seeds by total maintenance ΔV, max V∞, radial span,
  period error, taxi cost. Keep top-N."
* §6 — public interface sketch (the per-cycler ``score`` and the
  top-N ``rank`` reducer).
* §12(d) — hard constraints (``V∞ <= cap``, ``r_p >= r_p_min``,
  ``bend <= max``) enforced as inequalities, not soft regularisers.
* §16 — catalogue record schema; ``Score`` is the metric portion of
  that record once M7 wires it.

M4 boundaries (binding)
-----------------------
* ``score()`` reduces a *built* :class:`~cyclerfinder.model.Cycler` to
  the per-cycler metrics. It does **not** build cyclers (M3) and does
  **not** search the timing inside a cell (M5).
* ``Score.total_maintenance_dv_kms`` is computed via
  :func:`cyclerfinder.core.flyby.flyby_dv_for` (bend-and-magnitude
  decomposition) and therefore disagrees with the M3
  :meth:`Cycler.maintenance_dv` (naive velocity-discontinuity sum) on
  infeasible flybys. This is deliberate; M3's method is not removed,
  and existing M3 tests continue to consume it.
* ``taxi_cost_kms`` is a surrogate (max ``||vinf_in||`` over Earth
  encounters); see :func:`taxi_cost_kms` for the rationale and caveats.

Plan: ``docs/phases/m4-enumeration-scoring/plan.md`` §3.2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from cyclerfinder.core.constants import (
    DAYS_PER_JULIAN_YEAR,
    PLANETS,
    SAFE_PERIHELION_KM,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.flyby import flyby_dv_for, is_ballistic_feasible

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.model.cycler import Cycler


# ---------------------------------------------------------------------------
# Default composite weights (spec §5 step 4; tunable, not load-bearing)
# ---------------------------------------------------------------------------

DEFAULT_WEIGHTS: dict[str, float] = {
    "total_maintenance_dv_kms": 1.0,
    "max_vinf_kms": 0.1,
    "period_error_yr": 10.0,
    "taxi_cost_kms": 0.5,
}
"""Default weights for :func:`composite_score`.

Mixed units by design: each axis is weighted in its own unit so the
weighted sum is a unit-mixed scalar — purely ordinal, never reported as
a "cost". M5/M8 will revisit; tests pin orderings, not absolute values.
"""

_TAXI_BODY: str = "E"
"""Body whose ``vinf_in`` magnitudes are reduced into the taxi-cost
surrogate. Earth in M4 (per Aldrin / human-mission convention)."""


# ---------------------------------------------------------------------------
# Score dataclass (spec §5 step 4 + §12(d))
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Score:
    """Per-cycler ranking metrics, derived once and immutable.

    Mirrors spec §5 step 4 (rank seeds by total maintenance ΔV, max
    V∞, radial span, period error, taxi cost) and spec §12(d) (hard
    inequalities ``V∞ <= cap``, ``r_p >= r_p_min``, ``bend <= max``).

    Attributes
    ----------
    total_maintenance_dv_kms:
        Sum of per-encounter powered-flyby ΔV (km/s) per
        :func:`~cyclerfinder.core.flyby.flyby_dv_for`. Exactly zero on
        an idealised cycler with all-ballistic flybys.
    max_vinf_kms:
        Maximum ``||vinf||`` over all encounters, km/s. Hard-constraint
        driver and a spec §10 degeneracy guard.
    radial_span_au:
        ``(min perihelion AU, max aphelion AU)`` over all legs.
    period_error_yr:
        Absolute difference between ``cycler.period`` and the
        caller-supplied target period, in Julian years. Zero when
        ``score()`` is called without a target (typical for hand-built
        test cyclers).
    taxi_cost_kms:
        Surrogate hyperbolic-rendezvous cost — see :func:`taxi_cost_kms`.
    hard_constraints_pass:
        ``True`` iff every encounter satisfies spec §12(d):
        ``max(||vinf_in||, ||vinf_out||) <= vinf_cap`` AND the bend
        angle is achievable at the body's gravity and minimum safe
        periapsis (:func:`~cyclerfinder.core.flyby.is_ballistic_feasible`).
    """

    total_maintenance_dv_kms: float
    max_vinf_kms: float
    radial_span_au: tuple[float, float]
    period_error_yr: float
    taxi_cost_kms: float
    hard_constraints_pass: bool


# ---------------------------------------------------------------------------
# taxi_cost_kms — surrogate hyperbolic-rendezvous cost
# ---------------------------------------------------------------------------


def taxi_cost_kms(cycler: Cycler, taxi_body: str = _TAXI_BODY) -> float:
    """Surrogate hyperbolic-rendezvous cost for a taxi mission.

    Definition (M4):

        ``taxi_cost = max(||enc.vinf_in||)`` over encounters whose body
        is ``taxi_body`` (Earth by default). Returns ``0.0`` when the
        cycler has no encounter at ``taxi_body``.

    Rationale: a taxi launched from Earth pays roughly the hyperbolic
    excess speed of the cycler at the Earth encounter as the propulsive
    component of its injection ΔV (plus a fixed escape component that is
    constant across candidates and therefore drops out of the *ranking*).
    The surrogate is monotone in the actual cost, which is what a
    ranking metric needs.

    Caveats (documented to anchor future work):

    * Returns ``0.0`` for a cycler with no Earth encounter (e.g. a
      hypothetical Venus-Mars-only cycler). The caller is responsible
      for using a different surrogate in that case.
    * Ignores the taxi return leg.
    * Does **not** include the taxi's transfer-from-Earth-orbit ΔV,
      because that's departure-body and parking-orbit dependent. The
      proper full taxi-trajectory cost is M5/M6+ territory if it ever
      becomes a ranking target.

    The choice of ``max`` over Earth encounters (rather than ``mean``)
    reflects that a real taxi must close the *worst* encounter, not the
    average — the bottleneck dominates.
    """
    vinfs = [
        float(np.linalg.norm(enc.vinf_in)) for enc in cycler.encounters if enc.body == taxi_body
    ]
    if not vinfs:
        return 0.0
    return max(vinfs)


# ---------------------------------------------------------------------------
# score — single-cycler reduction
# ---------------------------------------------------------------------------


def score(
    cycler: Cycler,
    ephem: Ephemeris,
    vinf_cap: float,
    target_period_sec: float | None = None,
    rp_factors: dict[str, float] | None = None,
) -> Score:
    """Reduce a built cycler to a :class:`Score` record.

    Parameters
    ----------
    cycler:
        A ``Cycler`` (typically from
        :func:`cyclerfinder.search.construct.construct_cycler`).
    ephem:
        Accepted for API symmetry with M6 ephemeris-mode scoring; **not
        consumed** by M4.
    vinf_cap:
        Hard ceiling on V∞ at every encounter, km/s (spec §12(d)).
    target_period_sec:
        Synodic-anchored target period the cell was built against. If
        ``None`` (default), ``period_error_yr = 0.0``. M4 callers from
        :func:`rank` always pass an explicit target.
    rp_factors:
        Per-body multiplier on
        :data:`~cyclerfinder.core.constants.SAFE_PERIHELION_KM`
        for the bend-feasibility check. Default ``1.0`` everywhere.
        Lets callers tighten the safe periapsis without editing the
        constants table.

    Returns
    -------
    Score
        Frozen, fully-populated; immutable record suitable for the M7
        catalogue.
    """
    del ephem  # forward-compatibility hook; M4 does not consult it
    if rp_factors is None:
        rp_factors = {}

    # Maintenance ΔV via the bend-and-magnitude decomposition.
    total_dv = 0.0
    for enc in cycler.encounters:
        total_dv += flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out)

    # Direct delegations to M3's Cycler methods.
    max_vinf = cycler.max_vinf()
    radial_span = cycler.radial_span()

    # Target-period error (years).
    if target_period_sec is None:
        period_error_yr = 0.0
    else:
        seconds_per_year = DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY
        period_error_yr = abs(cycler.period - target_period_sec) / seconds_per_year

    taxi = taxi_cost_kms(cycler)

    # Hard-constraint check per spec §12(d).
    hard_ok = True
    for enc in cycler.encounters:
        vin_norm = float(np.linalg.norm(enc.vinf_in))
        vout_norm = float(np.linalg.norm(enc.vinf_out))
        if max(vin_norm, vout_norm) > vinf_cap:
            hard_ok = False
            break
        rp_factor = rp_factors.get(enc.body, 1.0)
        try:
            mu_planet = PLANETS[enc.body].mu_km3_s2
            rp_min = rp_factor * SAFE_PERIHELION_KM[enc.body]
        except KeyError:
            hard_ok = False
            break
        if not is_ballistic_feasible(enc.vinf_in, enc.vinf_out, mu_planet, rp_min):
            hard_ok = False
            break

    return Score(
        total_maintenance_dv_kms=total_dv,
        max_vinf_kms=max_vinf,
        radial_span_au=radial_span,
        period_error_yr=period_error_yr,
        taxi_cost_kms=taxi,
        hard_constraints_pass=hard_ok,
    )


# ---------------------------------------------------------------------------
# composite_score — single sortable scalar
# ---------------------------------------------------------------------------


def composite_score(s: Score, weights: dict[str, float] | None = None) -> float:
    """Weighted-sum scalar suitable for ascending sort (lower = better).

    Returns ``math.inf`` if ``s.hard_constraints_pass`` is False, so any
    sort that mixes feasible and infeasible scores puts infeasibles last
    automatically. :func:`rank` *also* filters them out; both layers
    exist as defence-in-depth.

    Parameters
    ----------
    s:
        The :class:`Score` to reduce.
    weights:
        Per-axis weights. ``None`` (default) ⇒ :data:`DEFAULT_WEIGHTS`.
        Unknown keys are ignored; missing keys contribute zero.

    Returns
    -------
    float
        Unit-mixed scalar; do **not** report as "cost". Used solely to
        order :class:`Score` records.
    """
    if not s.hard_constraints_pass:
        return math.inf
    effective = DEFAULT_WEIGHTS if weights is None else weights
    total = 0.0
    # Iterate over DEFAULT_WEIGHTS keys deterministically (dict iteration
    # order is insertion-ordered in Python 3.7+). Caller-supplied dicts
    # with missing keys contribute zero for the absent axis.
    for field in DEFAULT_WEIGHTS:
        w = effective.get(field, 0.0)
        value = float(getattr(s, field))
        total += w * value
    return total


# ---------------------------------------------------------------------------
# rank — spec §5 step 4 top-N reducer
# ---------------------------------------------------------------------------


def rank(
    cyclers: list[Cycler],
    ephem: Ephemeris,
    vinf_cap: float,
    n_keep: int = 20,
    target_period_sec: float | None = None,
    weights: dict[str, float] | None = None,
) -> list[tuple[Score, Cycler]]:
    """Score every cycler, filter by hard constraints, sort, keep top N.

    Mirrors spec §5 step 4: "rank seeds by total maintenance ΔV, max
    V∞, radial span, period error, taxi cost. Keep top-N." The
    composite ordering is documented in :func:`composite_score`.

    Parameters
    ----------
    cyclers:
        Candidate cyclers. Empty list is allowed and returns ``[]``.
    ephem:
        Forwarded to :func:`score`; unused in M4.
    vinf_cap:
        Hard V∞ ceiling, km/s.
    n_keep:
        Maximum number of items in the returned list. Output may be
        shorter when fewer cyclers are feasible.
    target_period_sec:
        Optional target period; forwarded to :func:`score`.
    weights:
        Optional composite weights; forwarded to :func:`composite_score`.

    Returns
    -------
    list[tuple[Score, Cycler]]
        ``(score, cycler)`` pairs sorted by ``composite_score`` ascending
        (best first). Ties are broken by the cycler's index in the
        original input (stable sort), so the output is reproducible
        across runs.
    """
    if not cyclers:
        return []
    scored: list[tuple[Score, Cycler]] = []
    for cyc in cyclers:
        s = score(cyc, ephem, vinf_cap, target_period_sec=target_period_sec)
        if s.hard_constraints_pass:
            scored.append((s, cyc))
    # Python's sort is stable; ties keep the input order.
    scored.sort(key=lambda pair: composite_score(pair[0], weights))
    return scored[:n_keep]


__all__ = [
    "DEFAULT_WEIGHTS",
    "Score",
    "composite_score",
    "rank",
    "score",
    "taxi_cost_kms",
]
