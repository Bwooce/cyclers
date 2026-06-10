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

from cyclerfinder.search import dsm_leg, self_seeding

DAY_S = 86400.0
YEAR_DAYS = 365.25


@dataclass(frozen=True)
class DsmChainSeed:
    sequence: tuple[str, ...]
    x0: NDArray[np.float64]  # [t0, vinf_out0, alpha0, beta0, *tof, *eta]
    bounds: dsm_leg.DsmBounds
    arc_a_au: float  # coplanar descriptor arc shape
    arc_e: float
    transit_branch: str
    vinf_anchor_kms: float  # the row's sourced Russell-table V_inf cell


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
    vinf_list = {e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])}
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
    branches = self_seeding.g_arc_branches(aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m)
    if not branches:
        return None
    arc = branches[0]  # base short-way shape; the gate may retry others
    n_legs = len(sequence) - 1
    # Seed ToFs: use the arc branch's tof_g_days as the transit-leg seed,
    # the complementary g-arc duration (big_g - transit) for the other legs.
    g_tof_days = float(g_tof_yr * YEAR_DAYS)
    big_g_tof_days = float(big_g_tof_yr * YEAR_DAYS)
    transit_days = float(arc.tof_g_days)
    # Distribute remaining time across non-transit legs; at minimum 30 d.
    slack_days = max(30.0, big_g_tof_days - transit_days)
    tof_seed_days: list[float] = []
    for i in range(n_legs):
        body_a, body_b = sequence[i], sequence[i + 1]
        # Transit leg: E->M or M->E
        if {body_a, body_b} == {"E", "M"}:
            tof_seed_days.append(transit_days)
        else:
            tof_seed_days.append(max(slack_days, g_tof_days))
    eta_seed = tuple(0.0 for _ in range(n_legs))
    t0_seed_sec = 0.0
    bounds = dsm_leg.sequence_keyed_bounds(
        sequence=sequence,
        t0_window_sec=(
            -float(big_g_tof_days) * DAY_S,
            float(big_g_tof_days) * DAY_S,
        ),
        vinf_out0_bounds_kms=(max(0.5, vinf_e - 2.0), vinf_e + 2.0),
        charge_flyby_continuity=False,
    )
    x0 = dsm_leg.dsm_chain_decision_vector(
        t0_sec=t0_seed_sec,
        vinf_out0_kms=vinf_e,
        alpha0=0.0,
        beta0=0.0,
        tof_days_per_leg=tuple(tof_seed_days),
        eta_per_leg=eta_seed,
    )
    return DsmChainSeed(
        sequence=sequence,
        x0=x0,
        bounds=bounds,
        arc_a_au=float(arc.a_au),
        arc_e=float(arc.e),
        transit_branch=str(arc.branch),
        vinf_anchor_kms=float(vinf_m),
    )
