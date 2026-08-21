"""Resonant Atlas pilot, Stage A: family recovery + eigenvalue survey (`#859`).

`#858`'s pre-commitment review (`docs/notes/2026-08-21-858-campaigns-789-790-791-
review.md`, its own recommended pilot shape at Sec. 3.4/6) narrowed the original
`#789` 15-20-system census into a 4-system pilot -- 3 novel targets plus one
published-anchor positive control -- and split it into two gated stages. This
module is Stage A only: for each system, sweep every coprime resonance ratio
``p:q`` with ``p, q <= 8``, recover the two-body-seeded resonant family via
this project's own proven family-continuation machinery
(:func:`cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed` +
:func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi` +
:func:`cyclerfinder.search.cr3bp_continuation.continue_family` -- the SAME
functions the Jupiter-Europa/`#753`, Saturn-Titan/`#765`, Neptune-Triton/`#776`-
`#777`, and Earth-Moon/`#780` family-confirmation gates all reuse, confirming
`#764`'s own "system-agnostic" finding once more), and report the dominant
Floquet/monodromy eigenvalue magnitude at each coarse Jacobi-constant step
along the family.

**Stage A is reporting/triage, NOT discovery.** No manifold/homoclinic/
heteroclinic machinery is touched here (mirroring every family-confirmation
module's own explicit "connection stage out of scope" boundary, e.g.
``neptune_triton_resonant_families.py``'s own docstring). The deliverable is a
per-``(system, p:q)`` census of which continuation-branch members land
"in-band" -- ``|lambda|`` inside :data:`IN_BAND_LOW`/:data:`IN_BAND_HIGH`, the
empirically-tractable window `#858` Sec. 3.2/7 derives from every converged
connection this project's own prior connection searches (`#754`/`#767`/`#781`/
`#786`) have actually produced (`|lambda| ~= 55` to `~= 2130`; every attempt
beyond ~4e3 failed honestly, e.g. `#759`'s lambda~=4445, `#781`'s 4:7 at
1.46e4). Stage B (connection attempts on in-band cells only) is explicitly
NOT this module's job -- `#859`'s own registration gates it on a human
dispatch decision after Stage A's report, not an auto-fire.

HONEST SCOPE CAVEAT (carried forward from `#776`'s own "two-body seed lineage"
finding, reused verbatim here rather than re-derived): the naive
:func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
construction frequently converges to a DIFFERENT topology than its own
``p:q`` label implies (documented for Neptune-Triton at (4,3)/(2,3)/(1,2)/
(2,1); the SAME qualitative failure Anderson & Lo 2011 report for Jupiter-
Europa and Vaquero 2013 reproduces for Saturn-Titan). This module's own
``p``/``q`` cell labels are therefore SEED-CONSTRUCTION labels, not verified
family-identity labels -- Stage A still recovers and eigenvalue-surveys a
REAL family at each cell that converges, just not necessarily the literal
``p:q`` resonance its own two-body construction targeted. This is fine for
Stage A's own triage purpose (the eigenvalue-band question does not depend on
which label a family happens to carry) but must not be over-read as "the
(3,2) row is verified to BE a genuine 3:2 resonant orbit" without the same
per-row topology check `#776`/`#777` ran by hand.

``x0_sign=-1`` ONLY (the seed's periapsis placed opposite the secondary at
t=0): ``x0_sign=+1`` places the seed exactly at the secondary's own orbital
radius, the documented DOP853 step-size-collapse hazard
(:func:`~cyclerfinder.search.jovian_resonant_families.survey_candidates`'s own
docstring, independently confirmed to time out in `#776`'s own scratchpad
check). Surveying ``+1`` seeds too is a real widening option for the actual
Stage A run, deliberately left out of this pilot harness's default scope.

Pure: math/numpy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.cr3bp_continuation`,
:mod:`cyclerfinder.search.jovian_resonant_families` (reused directly, per
`#764`'s own confirmed system-agnostic finding),
:mod:`cyclerfinder.search.campaign_runner` (`#788`/`#800`, the checkpointed
grid-campaign substrate this module's cells and worker plug into).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_families as jrf
from cyclerfinder.search.campaign_runner import CellOutcome, CellStatus

# ---------------------------------------------------------------------------
# The empirically-tractable |lambda| band (`#858` Sec. 3.2/7 -- see module
# docstring for the six data points this envelope is drawn from). NOT a
# claim that cells outside this band are worthless -- just that no connection
# search in this project's own history has succeeded past it, in either
# direction (near-unit-circle families are equally intractable -- `#776`'s
# 5:6-LI at lambda=1.000008, e-folding ~4.7e6 tu, discarded by its own source
# paper).
# ---------------------------------------------------------------------------

IN_BAND_LOW = 50.0
IN_BAND_HIGH = 2500.0


def is_in_band(abs_lambda: float) -> bool:
    """Is ``abs_lambda`` inside this project's own empirically-tractable
    connection-search band (`#858` Sec. 3.2/7, module docstring)?"""
    return IN_BAND_LOW <= abs_lambda <= IN_BAND_HIGH


# ---------------------------------------------------------------------------
# The 4-system pilot list (`#859` registration / `#858` Sec. 6), subject to
# the per-system literature/corpus pre-check `#858` Sec. 3.4 mandates --
# see ``docs/notes/2026-08-21-859-resonant-atlas-pilot-harness.md`` for the
# check's own results (this list survived it unchanged).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StageASystem:
    """One pilot system: a ``core.satellites`` registry (primary, secondary)
    pair plus its role in the pilot design."""

    system_key: str
    primary: str
    secondary: str
    role: str  # "positive_control" | "novel_target" | "band_edge_probe"


STAGE_A_SYSTEMS: tuple[StageASystem, ...] = (
    StageASystem("uranus-oberon", "Uranus", "Oberon", "positive_control"),
    StageASystem("jupiter-ganymede", "Jupiter", "Ganymede", "novel_target"),
    StageASystem("uranus-titania", "Uranus", "Titania", "novel_target"),
    StageASystem("saturn-rhea", "Saturn", "Rhea", "band_edge_probe"),
)

#: Default resonance-ratio ceiling (`#859`'s own registration: "all coprime
#: p:q with p,q<=8").
DEFAULT_MAX_PQ = 8

#: Default continuation depth per cell ("~10 steps", `#859`/`#858` Sec. 6) --
#: n_steps=9 plus the seed member itself = 10 members per cell when the walk
#: does not stop earlier (fold/topology-jump/gauntlet-reject all stop it
#: honestly sooner, per `cr3bp_continuation`'s own contract).
DEFAULT_N_C_STEPS = 9

#: Default Jacobi-constant step size per continuation step. In the same
#: range this project's own prior family-continuation calls use (1e-4 to
#: 5e-4, grep'd across `jovian_resonant_families.py` /
#: `neptune_triton_resonant_families.py`) -- NOT re-derived per system or
#: per resonance ratio for this pilot harness; the real Stage A run should
#: treat this as a tuning knob per the honest caveat in
#: `docs/notes/2026-08-21-859-resonant-atlas-pilot-harness.md`.
DEFAULT_D_JACOBI = 5e-4

#: Seed-construction sign (module docstring: ``+1`` is the documented
#: DOP853 step-collapse hazard, never surveyed by this pilot harness).
DEFAULT_X0_SIGN = -1


def coprime_pairs(max_pq: int) -> list[tuple[int, int]]:
    """Every ``(p, q)`` with ``1 <= p, q <= max_pq`` and ``gcd(p, q) == 1``.

    Ordered by ascending ``p`` then ``q`` (deterministic -- the campaign
    runner's own resumability contract requires a stable cell ordering, see
    ``campaign_runner.run_grid_campaign``'s own docstring).
    """
    if max_pq < 1:
        raise ValueError(f"max_pq must be >= 1, got {max_pq}")
    return [
        (p, q) for p in range(1, max_pq + 1) for q in range(1, max_pq + 1) if math.gcd(p, q) == 1
    ]


def build_stage_a_cells(
    systems: tuple[StageASystem, ...] = STAGE_A_SYSTEMS,
    *,
    max_pq: int = DEFAULT_MAX_PQ,
    n_c_steps: int = DEFAULT_N_C_STEPS,
    d_jacobi: float = DEFAULT_D_JACOBI,
    x0_sign: int = DEFAULT_X0_SIGN,
) -> list[dict[str, Any]]:
    """Build the full Stage A cell grid: one cell per ``(system, p, q)``.

    Cells are plain JSON-safe dicts (the ``campaign_runner``/``parallel_sweep``
    pickle-safety contract -- see ``tests/search/test_campaign_runner.py``'s
    own worker convention). Cell order is systems-outer, ``coprime_pairs``-
    inner, so resuming a partially-completed campaign never reorders already-
    checkpointed indices even if a caller later narrows ``systems``... NOTE:
    narrowing ``systems``/``max_pq`` between runs DOES change indices (the
    grid is positional, per ``run_grid_campaign``'s own contract) -- only
    re-invoke with the SAME arguments to resume a given ``results.jsonl``.
    """
    cells: list[dict[str, Any]] = []
    for sysdef in systems:
        for p, q in coprime_pairs(max_pq):
            cells.append(
                {
                    "system_key": sysdef.system_key,
                    "primary": sysdef.primary,
                    "secondary": sysdef.secondary,
                    "role": sysdef.role,
                    "p": p,
                    "q": q,
                    "x0_sign": x0_sign,
                    "n_c_steps": n_c_steps,
                    "d_jacobi": d_jacobi,
                }
            )
    return cells


def _crossing_index_near_half_period(system: cr3bp.CR3BPSystem, orbit: cp.SymmetricOrbit) -> int:
    """Recover the 1-based x-axis-crossing index nearest ``orbit.t_half``.

    Mirrors ``correct_symmetric_fixed_jacobi``'s own internal
    ``half_crossings=None`` auto-determination (module docstring of
    ``neptune_triton_resonant_families.py``: "the SAME logic ... itself uses
    when half_crossings=None ... verified this task to recover the EXACT
    integer [hand-picked] for all ten gate rows"). Needed because
    :func:`cyclerfinder.search.cr3bp_continuation.continue_family` requires an
    explicit ``half_crossings: int`` (it holds the crossing index fixed across
    the whole walk to stay on one branch), unlike the single-point corrector.
    """
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    t_hi = 1.25 * orbit.period
    times, _states = cp._xaxis_crossings(  # established cross-module convention (see docstring)
        system, state0, t_hi, with_stm=False, rtol=1e-12, atol=1e-12
    )
    if len(times) == 0:
        raise RuntimeError("no x-axis crossings found for a converged orbit (unexpected)")
    return int(np.argmin(np.abs(times - orbit.t_half))) + 1


def stage_a_worker(cell: dict[str, Any]) -> CellOutcome:
    """One Stage A cell: seed, converge, continue, eigenvalue-classify.

    ``"hit"`` iff at least one gauntlet-passing continuation member lands
    in-band (:func:`is_in_band`); ``"miss"`` if the family recovers cleanly
    but every member is out-of-band (or the seed itself never converges --
    a real, reportable negative, not an error); ``"error"`` only for an
    unexpected exception (the worker contract in
    ``campaign_runner``'s own module docstring: workers must not let
    exceptions propagate, so every anticipated failure mode above is a
    durable ``"miss"``, not an ``"error"``).
    """
    system_key = cell["system_key"]
    p, q, x0_sign = int(cell["p"]), int(cell["q"]), int(cell["x0_sign"])
    label = f"{system_key} {p}:{q}"
    try:
        system = cr3bp.cr3bp_system(cell["primary"], cell["secondary"])
        try:
            seed = jrf.two_body_resonant_seed(p, q, x0_sign=x0_sign)
        except ValueError as exc:
            # A DETERMINISTIC, well-understood construction limit, not an
            # anomaly: two_body_resonant_seed fixes the seed's periapsis
            # radius at r=1 (the secondary's own orbital radius); the
            # vis-viva speed at that radius is only real for
            # a = (q/p)**(2/3) >= r_peri/2 = 0.5, i.e. p/q not too far above
            # ~2.83 (e.g. (3,1): a=0.481, just below the cutoff -- observed
            # this task). Recorded as a "miss" (an expected, explainable
            # negative for this p:q's seed construction), never an "error"
            # (which campaign_runner reserves for genuinely unanticipated
            # failures and which blocks empty-region emission).
            return CellOutcome(
                status="miss",
                payload={
                    "system": system_key,
                    "p": p,
                    "q": q,
                    "label": label,
                    "reason": "seed_domain_invalid",
                    "detail": repr(exc),
                    "n_members": 0,
                    "n_in_band": 0,
                    "stop_reason": "seed_domain_invalid",
                },
            )
        state0 = np.array([seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0])
        seed_jacobi = cr3bp.jacobi_constant(state0, system.mu)

        orbit = cp.correct_symmetric_fixed_jacobi(
            system,
            seed.x0,
            seed_jacobi,
            seed.period_full,
            ydot0_sign=float(x0_sign),
            half_crossings=None,
        )
        if not orbit.converged:
            return CellOutcome(
                status="miss",
                payload={
                    "system": system_key,
                    "p": p,
                    "q": q,
                    "label": label,
                    "seed_jacobi": seed_jacobi,
                    "reason": "seed_did_not_converge",
                    "n_members": 0,
                    "n_in_band": 0,
                    "stop_reason": "seed_no_converge",
                },
            )

        half_crossings = _crossing_index_near_half_period(system, orbit)

        branch = cc.continue_family(
            system,
            orbit,
            direction=+1,
            d_jacobi=float(cell["d_jacobi"]),
            n_steps=int(cell["n_c_steps"]),
            min_jacobi=-1.0e9,
            max_jacobi=1.0e9,
            half_crossings=half_crossings,
            ydot0_sign=float(x0_sign),
            seed_label=label,
        )

        members = [
            {
                "jacobi": m.jacobi,
                "abs_lambda": m.abs_lambda,
                "stable": m.stable,
                "period": m.period,
            }
            for m in branch.members
        ]
        n_in_band = sum(1 for m in branch.members if is_in_band(m.abs_lambda))
        max_abs_lambda = max((m.abs_lambda for m in branch.members), default=0.0)
        min_abs_lambda = min((m.abs_lambda for m in branch.members), default=0.0)

        payload = {
            "system": system_key,
            "p": p,
            "q": q,
            "label": label,
            "seed_jacobi": seed_jacobi,
            "half_crossings": half_crossings,
            "n_members": len(branch.members),
            "n_rejected": branch.n_rejected,
            "stop_reason": str(branch.stop_reason),
            "n_in_band": n_in_band,
            "max_abs_lambda": max_abs_lambda,
            "min_abs_lambda": min_abs_lambda,
            "members": members,
        }
        status: CellStatus = "hit" if n_in_band > 0 else "miss"
        return CellOutcome(status=status, payload=payload)
    except Exception as exc:  # worker contract: never let an exception escape (module docstring)
        return CellOutcome(status="error", error=f"{label}: {exc!r}")


__all__ = [
    "DEFAULT_D_JACOBI",
    "DEFAULT_MAX_PQ",
    "DEFAULT_N_C_STEPS",
    "DEFAULT_X0_SIGN",
    "IN_BAND_HIGH",
    "IN_BAND_LOW",
    "STAGE_A_SYSTEMS",
    "StageASystem",
    "build_stage_a_cells",
    "coprime_pairs",
    "is_in_band",
    "stage_a_worker",
]
