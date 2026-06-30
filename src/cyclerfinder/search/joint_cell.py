"""#318 Phase 2a — the unified real-ephemeris joint-cell evaluator.

The Phase-2 reframe (`docs/superpowers/plans/2026-06-30-318-phase2-realeph-joint-search-design.md`):
the four cycler-discovery axes — powered (A), multi-rev (B), 3D/broken-plane (C),
epoch-locked (D) — do NOT share a state model in the Phase-1 substrate (A/B/D are
heliocentric Lambert, C is a disconnected CR3BP corrector). They DO co-vary naturally in
ONE model: a **real-ephemeris trajectory**. This module is the joint-cell evaluator for
that model (the moon-tour case — the tractable short-cycler regime the design prioritises).

A :class:`JointCell` is one real-ephemeris cycler candidate; the four axes are its
explicit dimensions:

* **A — powered:** the per-flyby maneuver/defect ΔV charged to close the cycle.
* **B — multi-rev:** the per-leg ``(n_rev, branch)`` Lambert structure of the seed.
* **C — 3D / broken-plane:** INTRINSIC — the real ephemeris places the moons in their
  true 3D (inclined) positions; ``axis_c_max_abs_z_km`` is the out-of-plane extent the
  Phase-1 circular-coplanar substrate could not represent.
* **D — epoch-locked:** the departure epoch, which sets the real-ephemeris geometry.

:func:`evaluate_joint_cell` dispatches to the validated real-eph moon-tour lane
(``nbody.jovian.optimize_cycle`` — the #223-validated chained-cycle corrector) and returns
a unified :class:`JointCellVerdict`. This is the CHEAP patched-conic surrogate of the design
— the pre-filter that ranks cells before any expensive n-body shoot. The n-body shoot
(``jovian_shoot``) is the high-fidelity verdict, an opt-in follow-up on top survivors only
(the #480 compute lesson: short cyclers + analytic STM + surrogate-first).

Positive control (`tests/search/test_joint_cell.py`): the evaluator reproduces the
#223-validated Liang Member D CGCEC closure as a joint-cell — BEFORE any search is trusted
([[feedback_verify_gauntlet_with_positive_control]]).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cyclerfinder.nbody import jovian

SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class JointCell:
    """One real-ephemeris cycler candidate parameterised over the four #318 axes."""

    primary: str  # currently "Jupiter" (JUP365); other primaries are the generalisation
    sequence: tuple[str, ...]  # closed moon-tour (first == last)
    epoch_iso: str  # Axis D — departure epoch (TDB ISO)
    n_revs: tuple[int, ...]  # Axis B — per-leg revolution count
    branches: tuple[str, ...]  # Axis B — per-leg Lambert branch ("single"/"low"/"high")
    tof_seed_days: tuple[float, ...]  # per-leg ToF seed (optimised within bound_days)
    powered_min_alt_km: float = 50.0  # Axis A — flyby altitude floor for the maneuver model
    bound_days: float = 3.0  # epoch/ToF optimiser bound


@dataclass(frozen=True)
class JointCellVerdict:
    """Unified verdict for a joint-cell: the four axis coordinates + the closure metrics."""

    # --- the four axis coordinates (what the search sweeps over) ---
    axis_a_powered_dv_ms: float  # A: summed closure/maintenance defect ΔV
    axis_b_n_revs: tuple[int, ...]  # B: multi-rev structure
    axis_c_max_abs_z_km: float  # C: real-ephemeris out-of-plane extent (0 in coplanar model)
    axis_d_epoch_iso: str  # D: epoch
    # --- closure metrics ---
    closure_defect_ms: float
    feasible: bool
    cycle_tof_days: float
    vinf_kms: tuple[float, ...]
    min_alt_km: float
    altitudes_km: tuple[float, ...]
    notes: str = ""


def evaluate_joint_cell(
    cell: JointCell,
    ephem: jovian.JovianEphemeris,
    *,
    alt_max_km: float = 70000.0,
) -> JointCellVerdict:
    """Evaluate one joint-cell via the real-eph moon-tour patched-conic surrogate.

    Runs the #223-validated single-cycle corrector (``optimize_cycle``, cycle 1) at the
    cell's epoch with its multi-rev branch plan, then reads off the four axis coordinates
    and the closure metrics. ``axis_c_max_abs_z_km`` (the real 3D extent) is the max
    |z|-component of the encounter-moon positions over the cycle — the broken-plane
    information the Phase-1 coplanar substrate discarded.

    This is the cheap pre-filter (no n-body shoot). A cell is ``feasible`` iff the cycle
    closes with a positive ToF and every interior flyby altitude lies in
    ``[powered_min_alt_km, alt_max_km]``.
    """
    if len(cell.sequence) < 2 or cell.sequence[0] != cell.sequence[-1]:
        raise ValueError("JointCell.sequence must be a closed cycle (first == last)")
    n_legs = len(cell.sequence) - 1
    if not (len(cell.n_revs) == len(cell.branches) == len(cell.tof_seed_days) == n_legs):
        raise ValueError(f"n_revs/branches/tof_seed_days must each have {n_legs} entries")

    branch_plan = tuple(zip(cell.n_revs, cell.branches, strict=True))
    t0 = jovian.tdb_sec_from_iso(cell.epoch_iso)
    cyc, _ = jovian.optimize_cycle(
        t0,
        list(cell.tof_seed_days),
        ephem,
        cycle_index=1,
        vinf_in_prev=None,
        bound_days=cell.bound_days,
        min_alt_km=cell.powered_min_alt_km,
        sequence=cell.sequence,
        branch_plan=branch_plan,
    )

    # Axis C — the real-ephemeris out-of-plane extent (max |z| of the encounter moons).
    max_abs_z = 0.0
    for moon, t in zip(cell.sequence, cyc.epochs_sec, strict=True):
        r, _ = ephem.state(moon, float(t))
        max_abs_z = max(max_abs_z, abs(float(r[2])))

    finite_alts = [a for a in cyc.altitudes_km if np.isfinite(a)]
    min_alt = min(finite_alts) if finite_alts else float("nan")
    feasible = bool(
        cyc.cycle_tof_days > 0.0
        and finite_alts
        and all(cell.powered_min_alt_km <= a <= alt_max_km for a in finite_alts)
    )

    return JointCellVerdict(
        axis_a_powered_dv_ms=float(cyc.sum_defect_ms),
        axis_b_n_revs=tuple(int(n) for n in cell.n_revs),
        axis_c_max_abs_z_km=float(max_abs_z),
        axis_d_epoch_iso=cell.epoch_iso,
        closure_defect_ms=float(cyc.sum_defect_ms),
        feasible=feasible,
        cycle_tof_days=float(cyc.cycle_tof_days),
        vinf_kms=tuple(float(v) for v in cyc.vinf_kms),
        min_alt_km=float(min_alt),
        altitudes_km=tuple(float(a) for a in cyc.altitudes_km),
    )


def liang_member_d_cell() -> JointCell:
    """The #223-validated Liang Member D CGCEC cell — the Phase-2a positive control.

    Sequence Callisto-Ganymede-Callisto-Europa-Callisto, departing 2033-09-25, the
    4-leg (1-rev, high/low/high/low) branch plan from ``nbody.jovian.BRANCH_PLAN``. The
    evaluator must reproduce its near-ballistic closure (sub-m/s cycle-1 defect, all
    flybys above 100 km) — the proof the joint-cell evaluator is sound before any search.
    """
    revs = tuple(r for r, _ in jovian.BRANCH_PLAN)
    branches = tuple(b for _, b in jovian.BRANCH_PLAN)
    return JointCell(
        primary="Jupiter",
        sequence=jovian.CGCEC,
        epoch_iso="2033-09-25T18:04:43",
        n_revs=revs,
        branches=branches,
        tof_seed_days=(31.8973, 18.1697, 29.9343, 19.9747),
        powered_min_alt_km=100.0,
        bound_days=3.0,
    )


__all__ = [
    "JointCell",
    "JointCellVerdict",
    "evaluate_joint_cell",
    "liang_member_d_cell",
]
