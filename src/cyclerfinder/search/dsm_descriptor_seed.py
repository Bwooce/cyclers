"""Descriptor -> Takao-DSM closure seed + per-row closer (spec 2026-06-10).

Bridges the two multi-arc mechanisms: parse a catalogue row's 2-arc g/G
free-return descriptor, use :func:`cyclerfinder.search.self_seeding.g_arc_branches`
to get the coplanar arc shape + the transit branch matching the tabulated transit,
and assemble the decision vector + bounds that the Takao eta-DSM corrector
(:func:`cyclerfinder.search.dsm_leg.dsm_chain_correct`) consumes. Pure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search import dsm_leg, self_seeding

# V_inf anchor-match tolerance (km/s) — aligned to the campaign gate used in
# free_return_chain.py (vinf_tol_kms=0.5) and the closer-sweep _LEVEL_EVIDENCE
# entries ("emerged V_inf within 0.5 km/s"). Note: V1_TOLERANCE_MPS in
# verify/crosscheck.py (1e-3 m/s) is the lamberthub-agreement precision bound,
# not the anchor-match band — a different quantity.
V1_TOLERANCE_KMS: float = 0.5

DAY_S = 86400.0
YEAR_DAYS = 365.25

# Russell 2004 SS2.1: the generic-return search caps ToF at 6 body periods; this is
# the global ceiling on how many revolutions the lane will enumerate.
RUSSELL_GENERIC_RETURN_BODY_PERIOD_CAP = 6


def _body_period_days(body: str) -> float:
    """Sidereal period of a planet (days), from its semi-major axis."""
    a_km = PLANETS[body].sma_au * AU_KM
    return float(2.0 * np.pi * np.sqrt(a_km**3 / MU_SUN_KM3_S2) / DAY_S)


@dataclass(frozen=True)
class DsmChainSeed:
    sequence: tuple[str, ...]
    x0: NDArray[np.float64]  # [t0, vinf_out0, alpha0, beta0, *tof, *eta]
    bounds: dsm_leg.DsmBounds
    arc_a_au: float  # coplanar descriptor arc shape
    arc_e: float
    transit_branch: str
    vinf_anchor_kms: float  # the row's sourced Russell-table V_inf cell
    per_leg_tof_days: tuple[float, ...] = ()  # SOURCED per-leg seed ToFs (audit)
    max_revs: int = 0  # global Lambert rev cap (Russell-sourced; see seed builder)


def _descriptor_params(
    row: dict[str, Any],
) -> tuple[float, float, float, float, float, tuple[str, ...]] | None:
    """Extract (aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m, sequence) from
    the row's g/G descriptor, or None if the row has no per-arc descriptor.

    Mirrors scripts/triage_self_seeding.py's extraction — same field paths:
    - aphelion_au  <- row["orbit_elements"]["aphelion_au"]
    - g/G ToFs     <- row["free_return_arcs"][*]["tof_years"] (needs >= 2)
    - vinf_e/m     <- row["vinf_kms_at_encounters"] keyed by body
    - sequence     <- row["sequence_canonical"] split on "-"
    """
    aph = (row.get("orbit_elements") or {}).get("aphelion_au")
    # Skip encounters with a null V_inf (some rows publish a body with no V_inf cell);
    # a row then missing its E/M anchor falls through to None below rather than raising.
    vinf_list = {
        e["body"]: float(e["vinf_kms"])
        for e in (row.get("vinf_kms_at_encounters") or [])
        if e.get("body") is not None and e.get("vinf_kms") is not None
    }
    fra = row.get("free_return_arcs") or []
    g_tofs = [a.get("tof_years") for a in fra if a.get("tof_years") is not None]
    seq_str: str | None = row.get("sequence_canonical")

    if (
        len(g_tofs) < 2
        or aph is None
        or "E" not in vinf_list
        or "M" not in vinf_list
        or not seq_str
    ):
        return None

    sequence = tuple(seq_str.split("-"))
    if len(sequence) < 2:
        return None

    return (
        float(aph),
        float(g_tofs[0]),
        float(g_tofs[1]),
        float(vinf_list["E"]),
        float(vinf_list["M"]),
        sequence,
    )


def seed_dsm_chain_from_descriptor(row: dict[str, Any]) -> DsmChainSeed | None:
    """Build a Takao-DSM closure seed from a catalogue row's g/G descriptor.

    Returns None for rows without a per-arc descriptor (ocampo members and any
    ch4 rows with fewer than two g/G arc ToFs).
    """
    params = _descriptor_params(row)
    if params is None:
        return None
    aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m, sequence = params
    # A descriptor whose arc geometry does not actually reach the body (the spec's
    # OFF-FAMILY-NO-CLOSE case, e.g. 5.30ggF3 / 5.75ggF3) raises "orbit does not
    # reach body" from g_arc_branches -> there is no USABLE seed, so honour the
    # documented contract and return None rather than propagating the error.
    try:
        branches = self_seeding.g_arc_branches(aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m)
    except ValueError:
        return None
    if not branches:
        return None
    arc = branches[0]  # base short-way shape; the gate may retry others
    n_legs = len(sequence) - 1
    # Seed ToFs: use the arc branch's tof_g_days as the transit-leg seed,
    # the complementary g-arc duration (big_g - transit) for the other legs.
    # SOURCED per-leg seed ToFs (spec 2026-06-20): a same-body resonant leg seeds at
    # its PUBLISHED arc ToF (free_return_arcs tof_years x 365.25, in list order); a
    # cross-body transit leg seeds at the row's sourced invariants.transit_times_days
    # (in order). Falls back to the computed arc transit (arc.tof_g_days) only if the
    # row has no transit_times_days entry left.
    arc_tofs_days = [
        float(a["tof_years"]) * YEAR_DAYS
        for a in (row.get("free_return_arcs") or [])
        if a.get("tof_years") is not None
    ]
    transit_days_list = [
        float(t) for t in ((row.get("invariants") or {}).get("transit_times_days") or [])
    ]
    big_g_tof_days = float(big_g_tof_yr * YEAR_DAYS)  # retained for the bounds cap below
    tof_seed_days = []
    per_leg_rev_cap = []
    arc_i = 0
    transit_i = 0
    for i in range(n_legs):
        body_a, body_b = sequence[i], sequence[i + 1]
        if body_a == body_b:
            # Same-body resonant return -> the next published arc ToF.
            if arc_i >= len(arc_tofs_days):
                return None  # more resonant legs than descriptor arcs -> cannot seed
            leg_tof = arc_tofs_days[arc_i]
            arc_i += 1
            rev_body = body_a
        else:
            # Cross-body transit -> the next sourced transit time (else computed arc).
            if transit_i < len(transit_days_list):
                leg_tof = transit_days_list[transit_i]
                transit_i += 1
            else:
                leg_tof = float(arc.tof_g_days)
            inner = body_a if PLANETS[body_a].sma_au <= PLANETS[body_b].sma_au else body_b
            rev_body = inner
        tof_seed_days.append(leg_tof)
        period_days = _body_period_days(rev_body)
        per_leg_rev_cap.append(int(np.floor(leg_tof / period_days)) + 1)
    max_revs = min(max(per_leg_rev_cap), RUSSELL_GENERIC_RETURN_BODY_PERIOD_CAP)
    eta_seed = tuple(0.0 for _ in range(n_legs))
    t0_seed_sec = 0.0
    bounds_raw = dsm_leg.sequence_keyed_bounds(
        sequence=sequence,
        t0_window_sec=(
            -float(big_g_tof_days) * DAY_S,
            float(big_g_tof_days) * DAY_S,
        ),
        vinf_out0_bounds_kms=(max(0.5, vinf_e - 2.0), vinf_e + 2.0),
        charge_flyby_continuity=True,
    )
    # Cap any infinite ToF upper bounds that arise from same-body legs (e.g.
    # E-E, M-M), where the synodic period is undefined. Use the total arc
    # duration as a finite ceiling so the corrector's x0 clipping stays valid
    # (inf - inf = nan in the eps computation breaks least_squares).
    tof_upper_cap = max(big_g_tof_days * 2.0, 2000.0)
    upper_capped = bounds_raw.upper.copy()
    for i in range(n_legs):
        ui = 4 + i  # tof slot i in the flat bounds vector
        if not np.isfinite(upper_capped[ui]):
            upper_capped[ui] = tof_upper_cap
    bounds = dsm_leg.DsmBounds(lower=bounds_raw.lower, upper=upper_capped)
    # Charged (#162 vector-residual) layout: the corrector runs with
    # charge_flyby_continuity=True (the only mode that rewards the bend-feasible
    # low-V_inf basin), so the seed carries the 2*(n_legs-1) intermediate-flyby
    # direction coords (seeded 0). The closer (Task 2) runs the same charged mode.
    x0 = dsm_leg.dsm_chain_decision_vector(
        t0_sec=t0_seed_sec,
        vinf_out0_kms=vinf_e,
        alpha0=0.0,
        beta0=0.0,
        tof_days_per_leg=tuple(tof_seed_days),
        eta_per_leg=eta_seed,
        alpha_int_per_leg=tuple(0.0 for _ in range(n_legs - 1)),
        beta_int_per_leg=tuple(0.0 for _ in range(n_legs - 1)),
    )
    return DsmChainSeed(
        sequence=sequence,
        x0=x0,
        bounds=bounds,
        arc_a_au=float(arc.a_au),
        arc_e=float(arc.e),
        transit_branch=str(arc.branch),
        vinf_anchor_kms=float(vinf_m),
        per_leg_tof_days=tuple(tof_seed_days),
        max_revs=int(max_revs),
    )


@dataclass(frozen=True)
class DsmClosureResult:
    """Outcome of a per-row DSM closure attempt on the real ephemeris.

    Attributes
    ----------
    converged:
        True iff the Takao-DSM corrector's residual fell below ``tol_kms``
        (``max(|residual_vector|) < tol_kms`` in charged mode).
    max_residual_kms:
        Worst per-component residual magnitude (km/s) from the corrector.
    dv_dsm_kms:
        Per-leg interior DSM impulse magnitudes (km/s). CONSTRAINED output.
    vinf_per_encounter_kms:
        EMERGED V_inf at each body in the sequence (km/s). Sourced from
        ``DsmChainResult.vinf_out_kms`` (departure side per encounter).
    vinf_anchor_kms:
        The row's sourced Russell-table V_inf anchor (km/s) used for the
        anchor-match test.
    anchor_match:
        True iff the minimum distance between any emerged V_inf and the
        anchor is within :data:`V1_TOLERANCE_KMS` (0.5 km/s campaign gate).
    hyperbolic_impossible:
        True iff any emerged V_inf exceeds the heliocentric elliptic ceiling
        (~71.9 km/s), flagging a physically impossible result.
    seed:
        The seed that was built (None when the row had no descriptor).
    """

    converged: bool
    max_residual_kms: float
    dv_dsm_kms: tuple[float, ...]
    vinf_per_encounter_kms: tuple[float, ...]
    vinf_anchor_kms: float
    anchor_match: bool
    hyperbolic_impossible: bool
    seed: DsmChainSeed | None


def close_row_dsm(
    row: dict[str, Any],
    ephem: Ephemeris,
    *,
    tol_kms: float = 0.1,
    gradient: str = "lambert",
) -> DsmClosureResult:
    """Attempt to close *row* on the real ephemeris using the Takao-DSM genome.

    Builds the charged seed (via :func:`seed_dsm_chain_from_descriptor`),
    then calls :func:`cyclerfinder.search.dsm_leg.dsm_chain_correct` in
    charged mode (``charge_flyby_continuity=True``) and packages the result.

    Never raises: an infeasible seed (``None`` from ``_descriptor_params``) or
    a non-converging corrector both return a ``converged=False`` result.

    Parameters
    ----------
    row:
        Raw catalogue row dict (from ``load_catalog``).
    ephem:
        Ephemeris instance (use ``Ephemeris("astropy")`` for real DE440).
    tol_kms:
        Convergence tolerance passed to the DSM corrector (km/s).
    gradient:
        Gradient backbone for the corrector (#244 opt-in). ``"lambert"`` (default)
        is the historical Lambert+finite-difference lane, byte-unchanged.
        ``"fbs-analytic"`` evaluates each leg with the Ellison-2018 FBS match-point
        corrector driven by the analytic Jacobian (the parity-sweep lane).

    Returns
    -------
    DsmClosureResult
    """
    seed = seed_dsm_chain_from_descriptor(row)
    if seed is None:
        return DsmClosureResult(
            converged=False,
            max_residual_kms=float("nan"),
            dv_dsm_kms=(),
            vinf_per_encounter_kms=(),
            vinf_anchor_kms=float("nan"),
            anchor_match=False,
            hyperbolic_impossible=False,
            seed=None,
        )
    res = dsm_leg.dsm_chain_correct(
        seed.x0,
        sequence=seed.sequence,
        ephem=ephem,
        bounds=seed.bounds,
        tol_kms=tol_kms,
        charge_flyby_continuity=True,
        gradient=gradient,
    )
    # Extract per-encounter emerged V_inf from vinf_out_kms (departure side).
    # vinf_out_kms is dict[int, float] keyed by body index in the sequence.
    vinf_out: tuple[float, ...] = tuple(
        float(res.vinf_out_kms[i]) for i in sorted(res.vinf_out_kms)
    )
    anchor = seed.vinf_anchor_kms
    best = min((abs(v - anchor) for v in vinf_out), default=float("inf"))
    hyper = any(v > 71.9 for v in vinf_out)
    return DsmClosureResult(
        converged=bool(res.converged),
        max_residual_kms=float(res.max_residual_kms),
        dv_dsm_kms=tuple(float(d) for d in res.dv_dsm_per_leg_kms),
        vinf_per_encounter_kms=vinf_out,
        vinf_anchor_kms=anchor,
        anchor_match=best <= V1_TOLERANCE_KMS,
        hyperbolic_impossible=hyper,
        seed=seed,
    )
