"""§14 V2-powered gate (the 2026-06-07 class-split amendment).

Spec references
---------------
* §14 V2-powered (the amended multi-cycle maintenance-periodicity gate): ≥3
  consecutive cycles where **(a)** every planned encounter is achieved within the
  documented encounter tolerance WITH the documented per-cycle maintenance
  applied, AND **(b)** intra-cycle drift versus the cycle's planned trajectory
  stays bounded (reset at each maneuver).
* §16.7.12 (schema v4.5 ``validation_level``).

Why a dedicated gate (the task #134 finding)
--------------------------------------------
The *original* single V2 gate measured cross-cycle rotating-frame-repeat drift
(:func:`cyclerfinder.verify.real_closure.verify_real_closure`). For a *powered*
cycler that metric is structurally unsatisfiable: the maintenance maneuver shapes
velocity, not where the planets are, and the powered Aldrin is retargeted every
cycle by design — so it never returns to the same place relative to the
incommensurately-breathing planets (drift ~4.14e8 km / 3 laps, ≈2072x tolerance;
``tests/verify/test_aldrin_v2_v3_campaign.py``). The amended V2-powered gate asks
the operationally meaningful questions instead, *per cycle*:

* **(a) encounter success WITH maintenance applied.** Each cycle is re-solved
  in-family at its lap-shifted real launch window
  (:func:`cyclerfinder.search.maintain.optimise_aldrin_maintenance_dv`), which
  places every encounter at the *real* planet state at the cumulative-ToF epoch
  (position met exactly by construction) and charges the turn-deficit maintenance
  ΔV at the closure flyby. The honest encounter-achievement measure is then the
  interior-flyby ``||V∞_in| - |V∞_out||`` continuity (a pure gravity-assist
  preserves |V∞|; the maintenance maneuver pays the heliocentric turn deficit at
  the *closure* encounter, not the interior one) — it must fall under
  :data:`ENCOUNTER_VINF_TOL_KMS`. The maintenance ΔV must be strictly positive
  (the cycle is genuinely powered, not a ballistic ΔV≈0 neighbour) and finite.
* **(b) bounded intra-cycle drift vs the planned trajectory.** Within each cycle,
  the per-leg forward **Kepler** re-propagation residual — each leg propagated
  from its departure state (``r_planet + v_planet + V∞_out``) to its planned
  arrival planet position, reset at each encounter/maneuver — must fall under
  :data:`INTRA_CYCLE_DRIFT_TOL_KM`. Reused verbatim from
  :func:`cyclerfinder.verify.agreement.crosscheck_code_paths` (path c).

The verdict ``v2_powered_passed`` is awarded only when **every** propagated cycle
clears both clauses. This is the per-cycle unit V3's horizon-TCM gate then sums.

Tolerances (declared from existing conventions, not invented)
-------------------------------------------------------------
* :data:`ENCOUNTER_VINF_TOL_KMS` = 0.5 km/s — the existing free-return genome
  V∞-continuity convention (``free_return_v1.VINF_CONTINUITY_TOL_KMS``), the
  documented bound separating a genuinely continuous flyby from a forced /
  discontinuous one. Combined with the documented safe-altitude flyby (the
  sourced 200 km Earth flyby; ``maintain._ALDRIN_EARTH_FLYBY_ALT_KM``), under
  which the turn-deficit maintenance is computed.
* :data:`INTRA_CYCLE_DRIFT_TOL_KM` = 1.0 km — the existing §14 V1 Kepler
  forward-reprop bound (``agreement.KEPLER_REPROP_TOL_KM``).

Golden discipline: the maintenance-ΔV *magnitude* is OUR computed value
(McConaghy 2002 defers it) and is only sanity-bounded / reported, never asserted
against a sourced target. The gated quantities (V∞ continuity, reprop residual)
are code-path self-consistency checks, not rediscovered published numbers.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Final

import numpy as np

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv
from cyclerfinder.verify.agreement import crosscheck_code_paths
from cyclerfinder.verify.dv_band_acceptance import accept_maintenance_dv, dv_band_threshold

# Encounter tolerance for clause (a): interior-flyby |V∞| continuity, km/s. The
# existing free-return genome-honesty convention (free_return_v1).
ENCOUNTER_VINF_TOL_KMS: Final[float] = 0.5

# Intra-cycle drift bound for clause (b): per-leg Kepler forward-reprop residual,
# km. The existing §14 V1 Kepler-reprop convention (agreement.KEPLER_REPROP_TOL_KM).
INTRA_CYCLE_DRIFT_TOL_KM: Final[float] = 1.0

# Upper sanity bound on the per-cycle maintenance ΔV (km/s). NOT a sourced target
# (McConaghy 2002 defers the magnitude); a loose bound that rejects an
# off-family degenerate (the ~55 km/s high-energy basin, finding #114). This is
# the generic (``dv_band=None``) per-cycle ceiling; when a row carries a
# ``dv_band`` the per-cycle acceptance is band-parameterised via
# :func:`cyclerfinder.verify.dv_band_acceptance.accept_maintenance_dv` (task
# #420) — the powered-DSM band's window has this same 3.5 km/s upper bound.
_MAINTENANCE_DV_SANITY_MAX_KMS: Final[float] = 3.5


def _maintenance_dv_ok(dv_kms: float, *, dv_band: str | None) -> bool:
    """Per-cycle maintenance-ΔV acceptance for clause (a) (task #420).

    When ``dv_band`` is ``None`` (or carries no impulsive window), this is the
    existing generic criterion verbatim: a strictly-positive, sanity-bounded
    maintenance ΔV (``0 < dv < _MAINTENANCE_DV_SANITY_MAX_KMS``) — the cycle is
    genuinely powered, not a ballistic ΔV≈0 neighbour, and not the off-family
    ~55 km/s degenerate.

    When ``dv_band`` carries a window, the per-cycle ΔV is judged against the
    band's per-7-cycle acceptance window scaled to a single cycle
    (:func:`cyclerfinder.verify.dv_band_acceptance.accept_maintenance_dv` with
    ``n_cycles=1``). The band window already encodes the powered-vs-ballistic
    floor (powered_dsm requires >= 300 m/s / 7 cycles; the ballistic bands admit
    near-zero), so the generic ``> 0`` rule is *subsumed* by the window and not
    re-applied — a ballistic band correctly accepts a ΔV≈0 close.
    """
    if dv_band_threshold(dv_band) is None:
        return 0.0 < dv_kms < _MAINTENANCE_DV_SANITY_MAX_KMS
    return accept_maintenance_dv(
        dv_kms * 1000.0,
        dv_band=dv_band,
        n_cycles=1,
        generic_max_mps=_MAINTENANCE_DV_SANITY_MAX_KMS * 1000.0,
    )


@dataclass(frozen=True)
class CycleResult:
    """Per-cycle V2-powered measurement.

    Attributes
    ----------
    converged:
        ``True`` when the in-family maintenance solve converged.
    a_au, e:
        Recovered Earth→Mars transfer ellipse elements (the sourced Aldrin
        anchors when in-family: a≈1.59 AU, e≈0.393).
    maintenance_dv_kms:
        Per-cycle maintenance ΔV (km/s). Strictly positive for a genuinely
        powered, in-family cycle.
    encounter_vinf_continuity_kms:
        Interior Mars-flyby ``||V∞_in| - |V∞_out||`` (km/s) — clause (a).
    intra_cycle_drift_km:
        Worst per-leg Kepler forward-reprop residual within the cycle (km) —
        clause (b).
    encounter_ok:
        Clause (a): ``encounter_vinf_continuity_kms <= ENCOUNTER_VINF_TOL_KMS``
        AND a strictly-positive, sanity-bounded maintenance ΔV.
    drift_ok:
        Clause (b): ``intra_cycle_drift_km <= INTRA_CYCLE_DRIFT_TOL_KM``.
    """

    converged: bool
    a_au: float
    e: float
    maintenance_dv_kms: float
    encounter_vinf_continuity_kms: float
    intra_cycle_drift_km: float
    encounter_ok: bool
    drift_ok: bool


@dataclass(frozen=True)
class V2PoweredResult:
    """§14 V2-powered verdict over ``n_cycles`` consecutive cycles.

    Attributes
    ----------
    n_cycles:
        Number of consecutive cycles evaluated (≥3 for the gate).
    per_cycle:
        The :class:`CycleResult` for each cycle.
    v2_powered_passed:
        ``True`` only when every cycle cleared both clause (a) and clause (b).
    detail:
        Human-readable note (the failure reason when not passed).
    """

    n_cycles: int
    per_cycle: tuple[CycleResult, ...]
    v2_powered_passed: bool
    detail: str


def verify_aldrin_v2_powered(
    ephem: Ephemeris,
    *,
    priority_date: datetime,
    period_years: float,
    n_cycles: int = 3,
    mu: float = MU_SUN_KM3_S2,
    dv_band: str | None = None,
) -> V2PoweredResult:
    """Run the §14 V2-powered gate on the Aldrin E→M→E cycler.

    For each of ``n_cycles`` consecutive cycles, re-phase the priority date one
    cycler period forward, re-solve the in-family maintenance cycle at that real
    launch window, and measure clauses (a) and (b). The verdict passes only when
    every cycle clears both clauses.

    Parameters
    ----------
    ephem:
        Real-ephemeris backend (``Ephemeris("astropy")`` for DE440).
    priority_date:
        The cycle-0 priority date (the catalogue literature epoch).
    period_years:
        The cycler repeat period (years) used to lap-shift the priority date
        between cycles (the Aldrin k=1 repeat is one E-M synodic ≈ 2.135 yr).
    n_cycles:
        Number of consecutive cycles (≥3 for the gate). Default 3.
    mu:
        Heliocentric gravitational parameter, km³/s².
    dv_band:
        Optional catalogue ``dv_band`` (task #420). When supplied, the per-cycle
        maintenance-ΔV half of clause (a) is judged against the band's
        per-7-cycle acceptance window (scaled to a single cycle) instead of the
        generic ``0 < dv < 3.5 km/s`` ceiling — so a ``powered_dsm`` row must
        close to its strictly-positive powered budget (>= 300 m/s / 7 cycles,
        not zero), while a ballistic band would require near-zero ΔV. ``None``
        (the default, most rows) keeps the existing generic criterion verbatim;
        a row is never promoted on a band it does not carry.

    Returns
    -------
    V2PoweredResult
    """
    if n_cycles < 3:
        raise ValueError(f"V2-powered requires n_cycles >= 3; got {n_cycles}")

    cycles: list[CycleResult] = []
    for k in range(n_cycles):
        pdate = priority_date + timedelta(days=k * period_years * DAYS_PER_JULIAN_YEAR)
        res = optimise_aldrin_maintenance_dv(ephem, real_window_priority_date=pdate)
        cyc = res.cycler

        # Clause (a): interior Mars-flyby |V∞| continuity.
        mars = cyc.encounters[1]
        vinf_in = float(np.linalg.norm(mars.vinf_in))
        vinf_out = float(np.linalg.norm(mars.vinf_out))
        continuity = abs(vinf_in - vinf_out)
        dv = float(res.maintenance_dv_kms)
        dv_ok = _maintenance_dv_ok(dv, dv_band=dv_band)
        encounter_ok = res.converged and continuity <= ENCOUNTER_VINF_TOL_KMS and dv_ok

        # Clause (b): intra-cycle Kepler forward-reprop residual (reset per leg).
        report = crosscheck_code_paths(cyc, ephem, mu=mu)
        intra_drift = float(report.kepler_reprop.max_residual_km)
        drift_ok = report.kepler_reprop.available and intra_drift <= INTRA_CYCLE_DRIFT_TOL_KM

        cycles.append(
            CycleResult(
                converged=bool(res.converged),
                a_au=float(res.a_au),
                e=float(res.e),
                maintenance_dv_kms=dv,
                encounter_vinf_continuity_kms=continuity,
                intra_cycle_drift_km=intra_drift,
                encounter_ok=encounter_ok,
                drift_ok=drift_ok,
            )
        )

    passed = all(c.encounter_ok and c.drift_ok for c in cycles)
    if passed:
        detail = (
            f"§14 V2-powered pass: {n_cycles} consecutive in-family cycles, each "
            f"achieving its encounters with maintenance applied (V∞ continuity "
            f"<= {ENCOUNTER_VINF_TOL_KMS} km/s) and bounded intra-cycle drift "
            f"(<= {INTRA_CYCLE_DRIFT_TOL_KM} km)"
        )
    else:
        bad = [
            (i, "encounter" if not c.encounter_ok else "", "drift" if not c.drift_ok else "")
            for i, c in enumerate(cycles)
            if not (c.encounter_ok and c.drift_ok)
        ]
        detail = f"§14 V2-powered fail in cycle(s): {bad}"

    return V2PoweredResult(
        n_cycles=n_cycles,
        per_cycle=tuple(cycles),
        v2_powered_passed=passed,
        detail=detail,
    )


__all__ = [
    "ENCOUNTER_VINF_TOL_KMS",
    "INTRA_CYCLE_DRIFT_TOL_KM",
    "CycleResult",
    "V2PoweredResult",
    "verify_aldrin_v2_powered",
]
