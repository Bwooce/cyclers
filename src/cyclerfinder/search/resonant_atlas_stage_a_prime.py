"""Resonant Atlas pilot, Stage A': fold-turning family recovery (`#861`).

`#860`'s advice (`docs/notes/2026-08-21-860-resonant-seeding-topology-fix-advice.md`,
independently re-verified by the coordinating session) diagnosed `#859`'s Stage A
smoke-test failure -- the naive :func:`~cyclerfinder.search.jovian_resonant_families.
two_body_resonant_seed` + :func:`~cyclerfinder.search.cr3bp_continuation.continue_family`
pipeline converges cleanly at every tested Uranus-Oberon published ratio but lands on a
STABLE, near-unit-circle family every time, never the paper's own labeled UNSTABLE
saddle family -- as two compounding, well-understood causes: (a) the seed's phase
structurally selects the encounter-AVOIDING branch (fixed by
:func:`~cyclerfinder.search.jovian_resonant_families.two_body_conjugate_apse_seed`,
`#861`'s own sibling addition, NOT this module), and (b)
:func:`~cyclerfinder.search.cr3bp_continuation.continue_family`'s NATURAL-PARAMETER
Jacobi walk cannot cross a saddle-center fold -- exactly the failure mode
:mod:`cyclerfinder.search.cr3bp_jacobi_arclength` (`#249`) was already built and proven
to solve (its own docstring: "the 1-DOF perpendicular-x-crossing symmetric corrector
lands only on the stable branch; natural-parameter Jacobi continuation diverges at the
fold" -- and it already recovered the previously-unrecoverable unstable Earth-Moon
C11a/C21 members from a stable-branch start).

This module is the Stage A' worker: from EITHER seed phase (the naive opposition seed
OR the new conjugate-apse seed -- landing stable is fine as a starting point per `#860`
Sec. 4(c2)), converge one member via the existing
:func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`, then walk
the family curve in BOTH C directions via
:func:`~cyclerfinder.search.cr3bp_jacobi_arclength.continue_in_jacobi` (turning any
fold it meets) toward a target far beyond any real family's own extent, so the walk
terminates at the family's OWN natural boundary (``STEP_UNDERFLOW``/``NO_MEMBER``) or
the step budget, not an arbitrarily chosen local neighborhood -- this also supersedes
the `#859` note's own "``d_jacobi``/``n_c_steps`` under-samples the existence range"
finding (Sec. 4(c2) of the `#860` advice), a second fix from the same change.

Every recovered member is classified three ways, per `#861`'s own gate criteria:
  - ``abs_lambda`` in :data:`~cyclerfinder.search.resonant_atlas_stage_a.IN_BAND_LOW`/
    ``IN_BAND_HIGH`` (reused directly, `#858` Sec. 3.2/7's tractable band);
  - ``period_over_2pi`` vs the label's own integer ``q`` (the resonance-lineage check
    every prior family-confirmation module in this project's history uses --
    `jovian_resonant_families.py`/`neptune_triton_resonant_families.py` -- with the
    SAME "few-percent, reviewer-judgment" tolerance the `#755` reviewer ruling
    established, not the strict `TABLE1_PERIOD_REL_TOL=1e-2` mechanical gate alone);
  - ``winding_k1`` (about the PRIMARY -- :func:`~cyclerfinder.search.binary_star_search.
    winding_topology`) vs the label's own integer ``p`` (spacecraft revolutions per
    Anderson & Lo's own p:q convention, `jovian_resonant_families.py`'s module
    docstring), and the closest approach to the SECONDARY
    (:func:`~cyclerfinder.search.jovian_resonant_families.europa_closest_approach`,
    generic despite its name -- takes ``system``, not a hardcoded body) as the
    close-flyby instability-mechanism signature Anderson & Lo (p.177-178) and this
    project's own `#758` evidence both attribute genuine unstable resonant families to.

Stage A' is reporting/triage, mirroring Stage A's own explicit boundary (module
docstring of `resonant_atlas_stage_a.py`): no manifold/homoclinic machinery here.

Pure: math/numpy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.cr3bp_jacobi_arclength`,
:mod:`cyclerfinder.search.jovian_resonant_families`,
:mod:`cyclerfinder.search.binary_star_search`,
:mod:`cyclerfinder.search.resonant_atlas_stage_a` (``is_in_band`` /
``_crossing_index_near_half_period`` reused directly, not re-derived),
:mod:`cyclerfinder.search.campaign_runner` (``CellOutcome``/``CellStatus`` contract,
for compatibility with a future checkpointed dispatch -- not exercised via
``run_grid_campaign`` in this task's own Oberon-only gate script).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_jacobi_arclength as ja
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_families as jrf
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.campaign_runner import CellOutcome, CellStatus
from cyclerfinder.search.mu_continuation import MuMember
from cyclerfinder.search.resonant_atlas_stage_a import (
    IN_BAND_HIGH,
    IN_BAND_LOW,
    _crossing_index_near_half_period,
    is_in_band,
)

SeedKind = Literal["opposition", "conjugate_apse"]

#: Default period_over_2pi vs label-q tolerance for a "topology matches" verdict --
#: the "few-percent" reviewer-judgment tolerance the `#755` ruling established
#: (3:4-LO confirmed at a 2.1% offset; 5:6-LO's own 2.83% offset was the largest
#: precedent this project has accepted as evidence of genuine lineage). Distinct
#: from (looser than) `jovian_resonant_families.TABLE1_PERIOD_REL_TOL` (1e-2), which
#: is the STRICT mechanical dual-criterion gate those modules keep even for rows
#: with sourced-published IC tables; this module's own topology check is reported
#: at BOTH thresholds so a reader can see which line a given member crosses.
PERIOD_REVIEWER_REL_TOL = 0.03


@dataclass(frozen=True)
class FoldTurnMember:
    """One classified member along a fold-turned family curve."""

    jacobi: float
    abs_lambda: float
    stable: bool
    period: float
    x0: float
    ydot0: float
    period_over_2pi: float
    period_rel_err_from_q: float
    winding_k1_rotating: int
    winding_k2_rotating: int
    winding_p_inertial: float
    winding_p_rel_err: float
    closest_secondary_approach_nondim: float
    crossing_residual: float
    radau_djacobi: float
    in_band: bool
    period_matches_q: bool  # within PERIOD_REVIEWER_REL_TOL
    winding_matches_p: bool  # within PERIOD_REVIEWER_REL_TOL


@dataclass
class FoldTurnResult:
    """The full fold-turned census for one ``(system, p, q, seed_kind)`` cell."""

    system_key: str
    p: int
    q: int
    seed_kind: SeedKind
    seed_converged: bool
    seed_jacobi: float = float("nan")
    half_crossings: int = 0
    ydot0_sign: float = 1.0
    members: list[FoldTurnMember] = field(default_factory=list)
    stop_reason_up: str = ""
    stop_reason_down: str = ""
    c_min: float = float("nan")
    c_max: float = float("nan")
    n_in_band: int = 0
    n_unstable: int = 0
    n_period_matches_q: int = 0
    max_abs_lambda: float = 0.0
    best_member_index: int | None = None  # index into `members` of the max-abs_lambda one


def build_seed(p: int, q: int, seed_kind: SeedKind) -> tuple[float, float, float, float]:
    """Return ``(x0, ydot0, ydot0_sign, period_full)`` for one seed phase.

    ``"opposition"`` is the `#859` harness's own (only) phase --
    :func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
    with ``x0_sign=-1``. ``"conjugate_apse"`` is `#861`'s new encounter-phase seed
    (:func:`~cyclerfinder.search.jovian_resonant_families.two_body_conjugate_apse_seed`).
    """
    if seed_kind == "opposition":
        s = jrf.two_body_resonant_seed(p, q, x0_sign=-1)
        return s.x0, s.ydot0, -1.0, s.period_full
    if seed_kind == "conjugate_apse":
        c = jrf.two_body_conjugate_apse_seed(p, q)
        sign = 1.0 if c.ydot0 >= 0.0 else -1.0
        return c.x0, c.ydot0, sign, c.period_full
    raise ValueError(f"seed_kind must be 'opposition' or 'conjugate_apse'; got {seed_kind!r}")


def classify_member(system: cr3bp.CR3BPSystem, m: MuMember, p: int, q: int) -> FoldTurnMember:
    """Classify one fold-turn member against the ``p:q`` label's own topology.

    Two independent checks, both DERIVED (not asserted against themselves) from
    :func:`~cyclerfinder.search.binary_star_search.winding_topology`'s raw
    ROTATING-frame winding and the member's own period:

    - ``period_over_2pi`` vs ``q``: the secondary completes exactly ``period /
      (2*pi)`` of its own revolutions over one member period (its own orbital
      rate is 1 nondim by the rotating frame's own definition) -- this directly
      checks the label's own ``q``.
    - ``winding_p_inertial`` vs ``p``: Anderson & Lo's own p:q convention is
      "spacecraft revolutions : secondary revolutions = p:q" -- i.e. ``p`` is
      the spacecraft's own INERTIAL revolution count about the PRIMARY, not the
      rotating-frame winding ``winding_topology`` returns directly (verified
      empirically this task against the Neptune-Triton "4:5-saddle" table-
      verified saddle: rotating-frame ``k1=1``, but converting to inertial via
      ``w1_rotating + period/(2*pi)`` (the frame's own omega=1 rotation added
      back) gives 3.84, matching the label's own ``p=4`` to ~4% -- the SAME
      order of "few-percent, far-from-integrable" offset the period check
      itself shows for this row, not a coincidence). ``winding_p_inertial =
      w1_rotating + period/(2*pi)``.
    """
    period_over_2pi = m.period / (2.0 * math.pi)
    rel_err_q = abs(period_over_2pi - q) / q
    topo = winding_topology(system.mu, m.state0, m.period)
    winding_p_inertial = topo.w1 + m.period / (2.0 * math.pi)
    rel_err_p = abs(abs(winding_p_inertial) - p) / p
    closest = jrf.europa_closest_approach(system, m.x0, m.ydot0, m.period)
    return FoldTurnMember(
        jacobi=m.jacobi,
        abs_lambda=m.abs_lambda,
        stable=m.stable,
        period=m.period,
        x0=m.x0,
        ydot0=m.ydot0,
        period_over_2pi=period_over_2pi,
        period_rel_err_from_q=rel_err_q,
        winding_k1_rotating=topo.k1,
        winding_k2_rotating=topo.k2,
        winding_p_inertial=winding_p_inertial,
        winding_p_rel_err=rel_err_p,
        closest_secondary_approach_nondim=closest,
        crossing_residual=m.crossing_residual,
        radau_djacobi=m.radau_djacobi,
        in_band=is_in_band(m.abs_lambda),
        period_matches_q=rel_err_q < PERIOD_REVIEWER_REL_TOL,
        winding_matches_p=rel_err_p < PERIOD_REVIEWER_REL_TOL,
    )


def fold_turn_family(
    system: cr3bp.CR3BPSystem,
    p: int,
    q: int,
    seed_kind: SeedKind,
    *,
    system_key: str = "",
    c_span: float = 0.12,
    ds0: float = 6e-3,
    ds_max: float = 2.5e-2,
    ds_min: float = 1e-6,
    max_steps: int = 150,
    record_every: int = 2,
    corrector_tol: float = 1e-10,
    radau_closure_tol: float = 1e-6,
    radau_jacobi_tol: float = 1e-7,
) -> FoldTurnResult:
    """Converge one ``(p, q, seed_kind)`` seed and fold-turn its family both
    directions in C to the family's own natural boundaries, classifying every
    gauntlet-passing member.

    ``c_span`` is the (large, deliberately generous) one-sided target offset from
    the seed's own natural C -- e.g. Anderson & Kumar 2024's own printed Oberon
    family C-ranges (`#728` digest) are all <0.045 wide, so ``c_span=0.12`` is
    comfortably beyond any real family's own extent, meaning the walk should
    terminate at the family's own fold/topology boundary (``STEP_UNDERFLOW``/
    ``NO_MEMBER``/``TOPOLOGY_JUMP``) well before reaching the nominal target, not
    at an arbitrary local-neighborhood cutoff (the `#859` note's own documented
    gap this module's docstring names).
    """
    x0, ydot0, ydot0_sign, period_full = build_seed(p, q, seed_kind)
    state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    seed_jacobi = float(cr3bp.jacobi_constant(state0, system.mu))
    orbit = cp.correct_symmetric_fixed_jacobi(
        system, x0, seed_jacobi, period_full, ydot0_sign=ydot0_sign, half_crossings=None
    )
    if not orbit.converged:
        return FoldTurnResult(
            system_key=system_key,
            p=p,
            q=q,
            seed_kind=seed_kind,
            seed_converged=False,
            seed_jacobi=seed_jacobi,
        )
    half_crossings = _crossing_index_near_half_period(system, orbit)

    branch_up = ja.continue_in_jacobi(
        orbit,
        mu=system.mu,
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
        c_target=orbit.jacobi + c_span,
        ds0=ds0,
        ds_max=ds_max,
        ds_min=ds_min,
        max_steps=max_steps,
        record_every=record_every,
        corrector_tol=corrector_tol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
        label=f"{system_key} {p}:{q} {seed_kind} up",
    )
    branch_down = ja.continue_in_jacobi(
        orbit,
        mu=system.mu,
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
        c_target=orbit.jacobi - c_span,
        ds0=ds0,
        ds_max=ds_max,
        ds_min=ds_min,
        max_steps=max_steps,
        record_every=record_every,
        corrector_tol=corrector_tol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
        label=f"{system_key} {p}:{q} {seed_kind} down",
    )

    # Merge: branch_down in reverse (excluding its own seed-duplicate first member)
    # + branch_up in order (its first member IS the shared seed).
    raw_members = list(reversed(branch_down.members[1:])) + list(branch_up.members)
    members = [classify_member(system, m, p, q) for m in raw_members]

    result = FoldTurnResult(
        system_key=system_key,
        p=p,
        q=q,
        seed_kind=seed_kind,
        seed_converged=True,
        seed_jacobi=seed_jacobi,
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
        members=members,
        stop_reason_up=str(branch_up.stop_reason),
        stop_reason_down=str(branch_down.stop_reason),
    )
    if members:
        cs = [m.jacobi for m in members]
        result.c_min = min(cs)
        result.c_max = max(cs)
        result.n_in_band = sum(1 for m in members if m.in_band)
        result.n_unstable = sum(1 for m in members if not m.stable)
        result.n_period_matches_q = sum(1 for m in members if m.period_matches_q)
        best_idx = max(range(len(members)), key=lambda i: members[i].abs_lambda)
        result.best_member_index = best_idx
        result.max_abs_lambda = members[best_idx].abs_lambda
    return result


def fold_turn_member_to_dict(m: FoldTurnMember) -> dict[str, Any]:
    """JSON-safe payload for one member (``campaign_runner``/checkpoint contract)."""
    return {
        "jacobi": m.jacobi,
        "abs_lambda": m.abs_lambda,
        "stable": m.stable,
        "period": m.period,
        "x0": m.x0,
        "ydot0": m.ydot0,
        "period_over_2pi": m.period_over_2pi,
        "period_rel_err_from_q": m.period_rel_err_from_q,
        "winding_k1_rotating": m.winding_k1_rotating,
        "winding_k2_rotating": m.winding_k2_rotating,
        "winding_p_inertial": m.winding_p_inertial,
        "winding_p_rel_err": m.winding_p_rel_err,
        "closest_secondary_approach_nondim": m.closest_secondary_approach_nondim,
        "crossing_residual": m.crossing_residual,
        "radau_djacobi": m.radau_djacobi,
        "in_band": m.in_band,
        "period_matches_q": m.period_matches_q,
        "winding_matches_p": m.winding_matches_p,
    }


def fold_turn_worker(cell: dict[str, Any]) -> CellOutcome:
    """``campaign_runner``-compatible worker wrapping :func:`fold_turn_family`.

    ``cell`` keys: ``system_key``, ``primary``, ``secondary``, ``p``, ``q``,
    ``seed_kind``, plus optional ``mu`` (overrides the registry mu, matching every
    sourced-paper module's own "paper's own mu, not the registry's" convention)
    and any of :func:`fold_turn_family`'s own keyword tuning knobs.
    """
    system_key = cell["system_key"]
    p, q, seed_kind = int(cell["p"]), int(cell["q"]), cell["seed_kind"]
    label = f"{system_key} {p}:{q} {seed_kind}"
    try:
        mu = cell.get("mu")
        if mu is not None:
            base = cr3bp.cr3bp_system(cell["primary"], cell["secondary"])
            system = cr3bp.CR3BPSystem(
                mu=float(mu),
                primary=cell["primary"],
                secondary=cell["secondary"],
                l_km=base.l_km,
                t_s=base.t_s,
            )
        else:
            system = cr3bp.cr3bp_system(cell["primary"], cell["secondary"])
        kwargs = {
            k: cell[k]
            for k in (
                "c_span",
                "ds0",
                "ds_max",
                "ds_min",
                "max_steps",
                "record_every",
                "corrector_tol",
                "radau_closure_tol",
                "radau_jacobi_tol",
            )
            if k in cell
        }
        result = fold_turn_family(system, p, q, seed_kind, system_key=system_key, **kwargs)
        if not result.seed_converged:
            return CellOutcome(
                status="miss",
                payload={
                    "system": system_key,
                    "p": p,
                    "q": q,
                    "seed_kind": seed_kind,
                    "label": label,
                    "reason": "seed_did_not_converge",
                    "seed_jacobi": result.seed_jacobi,
                    "n_members": 0,
                    "n_in_band": 0,
                },
            )
        payload = {
            "system": system_key,
            "p": p,
            "q": q,
            "seed_kind": seed_kind,
            "label": label,
            "seed_jacobi": result.seed_jacobi,
            "half_crossings": result.half_crossings,
            "ydot0_sign": result.ydot0_sign,
            "n_members": len(result.members),
            "stop_reason_up": result.stop_reason_up,
            "stop_reason_down": result.stop_reason_down,
            "c_min": result.c_min,
            "c_max": result.c_max,
            "n_in_band": result.n_in_band,
            "n_unstable": result.n_unstable,
            "n_period_matches_q": result.n_period_matches_q,
            "max_abs_lambda": result.max_abs_lambda,
            "best_member_index": result.best_member_index,
            "members": [fold_turn_member_to_dict(m) for m in result.members],
        }
        status: CellStatus = "hit" if result.n_in_band > 0 else "miss"
        return CellOutcome(status=status, payload=payload)
    except Exception as exc:  # worker contract: never let an exception escape
        return CellOutcome(status="error", error=f"{label}: {exc!r}")


__all__ = [
    "IN_BAND_HIGH",
    "IN_BAND_LOW",
    "PERIOD_REVIEWER_REL_TOL",
    "FoldTurnMember",
    "FoldTurnResult",
    "SeedKind",
    "build_seed",
    "classify_member",
    "fold_turn_family",
    "fold_turn_member_to_dict",
    "fold_turn_worker",
]
