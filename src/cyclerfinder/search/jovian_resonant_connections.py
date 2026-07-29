"""#754 Jupiter-Europa 3:4-LO homoclinic connection + Anderson & Lo 2011 Table 2 gate.

Re-scoped Task B (see `docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md`
Sec. 5) around the ONE family confirmed by `#755` (reviewer-ruled) out of the two
Anderson & Lo 2011 Table 1 rows Task B originally needed: `3:4-LO`. Table 2 of the
paper ("Homoclinic Trajectory State at Intersection", p.190) is the intersection of
3:4-LO's OWN stable and unstable manifolds with each other -- a self-connection, no
5:6 family required (`5:6-LO` remains a clean negative, independently retried by the
sibling task `#758`; this module does not touch it).

This module generalizes `cyclerfinder.genome.heteroclinic_cycle.correct_connection`
(unchanged in its own dynamics -- see that module's `x_sign` addition) from
libration-point Lyapunov-orbit nodes to a RESONANT periodic-orbit node
(:class:`ResonantNode`), built directly from a `#753`/`#755`-confirmed
:class:`~cyclerfinder.search.jovian_resonant_families.ResonantFamilyCandidate`'s
stored ``(x0, ydot0, period)`` -- no new stored fields, the saddle eigenvectors are
recomputed on the fly via the existing
:func:`~cyclerfinder.genome.heteroclinic_cycle._planar_floquet_pair` (one STM
propagation), exactly as `#757`'s own scoping note recommended (Sec. 2).

SOURCED CONSTANTS (verified directly against the PDF text layer + rendered pages,
this task -- see the module's own results note for the exact grep/page-read
provenance):

* :data:`TABLE2_STATE` = ``(x=-1.28427733, xdot=0.00000009, ydot=0.46372205)`` --
  Anderson & Lo 2011 Table 2 (p.190), "Homoclinic Trajectory State at Intersection",
  the paper's own single published row, at ``ANDERSON_LO_C_FLYBY``. Computed by the
  paper's own text (p.190) via INTERPOLATION between the closest sampled points on
  the two manifolds -- NOT a Newton-corrected exact intersection -- which is why the
  reproduction gate below (:data:`TABLE2_GATE_ABS_TOL`) is looser than the
  eigenvalue-reproduction gate in `jovian_resonant_families.py`.
* The paper's own one-sided Poincare section (p.170-171, Eq. 7, re-verified this
  task): ``{y=0}`` along the negative x-axis opposite Europa, one-sided ``ydot>0``.
  On the section ``(x, xdot)`` are the free coordinates and ``ydot`` is REDUNDANT --
  recovered from the Jacobi constant via Eq. 7,
  ``ydot = +sqrt(x^2 + 2(1-mu)/r1 + 2*mu/r2 - C - xdot^2)``
  (planar; the ``+`` branch matches the section's ``ydot>0`` convention and both
  Table 2/3 rows). :func:`ydot_from_section_eq7` implements this directly from
  :func:`cyclerfinder.core.cr3bp.jacobi_constant`'s own formula
  (``C = x^2+y^2 + 2(1-mu)/r1 + 2*mu/r2 - v^2``, y=0 here) -- so ``ydot`` is a
  DERIVED consistency check, not an independent third gate coordinate; the real
  gate is on ``(x, xdot)`` alone.

#759 UPDATE (2026-07-28, direct continuation of #754/#757/#758): with 5:6-LO
now ALSO reviewer-confirmed (docs/notes/2026-07-28-758-jupiter-europa-5-6-lo-table2-seeded-search.md
"Reviewer verdict"), Table 3 (p.191, "Heteroclinic Trajectory State at
Intersection") is buildable: it is ``Wu(3:4-LO) ∩ Ws(5:6-LO)``, "an
intersection of the unstable manifold of the 3:4 orbit and the stable
manifold of the 5:6 orbit **near the 3:4 orbit**" (p.191, verbatim -- the
paper explicitly selected the intersection nearest 3:4-LO's own section
point out of the several that exist, so a converged connection landing close
to 3:4-LO's own IC is the EXPECTED, sourced geometry here, not a red flag).
Re-verified directly against the paper's text layer (line ~1591 of the
``.txt`` sidecar) that Table 3, like Table 2, was "computed as before by
using interpolation between the closest points on the invariant manifolds in
the Poincare section" -- the SAME interpolation-not-Newton method, so the
same :data:`TABLE2_GATE_ABS_TOL`-class (``1e-4``) tolerance justification
applies verbatim; see :data:`TABLE3_GATE_ABS_TOL`.

:func:`find_heteroclinic` generalizes :func:`find_homoclinic` to two
DIFFERENT :class:`ResonantNode` endpoints (built via the same
``ResonantNode.from_candidate`` adapter, one from 3:4-LO's confirmed
candidate, one from 5:6-LO's --
:func:`~cyclerfinder.search.jovian_resonant_families.recover_758_table2_seeded_candidate`).
Since ``node_from != node_to`` there is no self-shadow trivial solution the
way homoclinic A=B has -- but :data:`HETERO_GHOST_GUARD_DELTA` still guards
against a genuinely degenerate corrector output (a "converged" crossing that
is really just the epsilon-scale seed offset with no real multi-crossing
manifold excursion), calibrated much tighter than :data:`GHOST_GUARD_DELTA`
precisely because landing close to 3:4-LO's own point is expected, not
degenerate, for this table (see that constant's own docstring for the full
reasoning and the specific numeric margin this task found).

#766 UPDATE (2026-07-29, direct continuation of #761): `#761` proved the
catalogued ``europa-3-4-crnbp-torus-jupiter-2026`` torus's own cited lineage
-- Kumar, Anderson, de la Llave & Gunter 2021's "arbitrarily chosen" seed
orbit at Jacobi constant :data:`~cyclerfinder.search.jovian_resonant_families.KUMAR_2021_C`
(= 3.0041) -- is a genuine real saddle on the SAME continuous family as the
confirmed 3:4-LO, via
:func:`~cyclerfinder.search.jovian_resonant_families.continue_34lo_to_kumar_c`.
This task builds a homoclinic SELF-connection AT that energy (NOT at
``ANDERSON_LO_C_FLYBY``, where :func:`find_homoclinic`'s own hit was already
found by `#754`) -- :func:`build_34lo_kumar_c_node` wraps the continuation's
own endpoint candidate into a :class:`ResonantNode` (reusing
``ResonantNode.from_candidate`` unchanged), and :func:`find_homoclinic` is
extended with an optional ``target``/``rank_by_residual`` pair of keyword
arguments so a caller scanning an energy with NO published Table-2-style
state (this task's honest situation -- Kumar 2021 reports no homoclinic
connection state for its own seed orbit) can rank surviving candidates by
Newton-residual tightness instead of distance to a nonexistent published
target. Existing callers (the Table-2/Table-3 gates above) are unaffected --
``target=None`` (the new default) still falls back to :data:`TABLE2_STATE`
exactly as before. See the results note
(``docs/notes/2026-07-29-766-torus-seed-homoclinic-connection-c30041.md``)
for the full scan log and self-consistency evidence (this task's own gate is
necessarily a self-consistency check -- Newton residual, ghost-guard margin,
independent Radau cross-check, forward/backward re-approach -- never a
reproduction claim, since no published number exists to reproduce here).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    _planar_floquet_pair,
    _seed_on_manifold,
    correct_connection,
)
from cyclerfinder.search.jovian_resonant_families import (
    ANDERSON_LO_C_FLYBY,
    ResonantFamilyCandidate,
    continue_34lo_to_kumar_c,
    jupiter_europa_system,
    recover_758_table2_seeded_candidate,
    recover_table1_candidate,
)

# ---------------------------------------------------------------------------
# Sourced constants -- Anderson & Lo 2011, Table 2 (p.190).
# ---------------------------------------------------------------------------

#: Table 2 (p.190), "Homoclinic Trajectory State at Intersection" -- verbatim.
#: (x, xdot, ydot); y == 0 is definitional to the section (not a free coordinate).
TABLE2_STATE: dict[str, float] = {
    "x": -1.28427733,
    "xdot": 0.00000009,
    "ydot": 0.46372205,
}

#: Gate tolerance on the Table-2 (x, xdot) match: absolute (both quantities are
#: O(1) or smaller in nondimensional units here, so absolute and relative coincide
#: in practice for x; xdot is near zero, where only an absolute tolerance is
#: meaningful). Looser than `jovian_resonant_families.TABLE1_GATE_REL_TOL` (1e-3
#: relative on an O(1e3) eigenvalue) because Table 2's own published digits are
#: INTERPOLATION-limited, not Newton-corrected (p.190: "computed using
#: interpolation with the closest points on the invariant manifolds"), and the
#: paper's own manifold offset epsilon is ~0.5e-5 (p.181-182) -- chasing agreement
#: tighter than the source itself can support would not be an honest gate. See
#: docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md Sec 5 item 3.
TABLE2_GATE_ABS_TOL = 1e-4

#: Homoclinic trivial-solution ghost-guard radius (nondimensional, in (x, xdot)
#: Euclidean norm). A "converged" A=B connection whose crossing lies within this
#: distance of the orbit's OWN qualifying section point(s) (see
#: :func:`own_section_points`) is the orbit trivially shadowing itself at
#: zero-length manifold offset, not a genuine homoclinic intersection. #757's own
#: scoping note documents a large margin here for 3:4-LO specifically (~0.146 in x
#: between the true homoclinic point and the orbit's own section point), so this
#: default is generous (not a delicate tolerance) while still being real (three
#: orders of magnitude tighter than that margin).
GHOST_GUARD_DELTA = 1e-3

#: The paper's own manifold-offset magnitude (p.181-182: "~0.5 x 10^-5"), used
#: verbatim as `correct_connection`'s `epsilon` for faithfulness to the paper's own
#: construction (this module's default `epsilon=1e-6` elsewhere in the codebase is
#: for the W-Z Sun-Jupiter-Oterma system, a different problem).
ANDERSON_LO_EPSILON = 0.5e-5

#: The paper's own one-sided section convention (p.170-171, Eq. 7): {y=0} along
#: the negative x-axis opposite Europa (x<0), ydot>0 half only.
SECTION_YDOT_SIGN = +1
SECTION_X_SIGN = -1

#: Table 3 (p.191), "Heteroclinic Trajectory State at Intersection" -- verbatim.
#: Wu(3:4-LO) ∩ Ws(5:6-LO), "near the 3:4 orbit" (p.191): the paper explicitly
#: selected the intersection nearest 3:4-LO's own section point out of the
#: several heteroclinic intersections that exist at this energy.
TABLE3_STATE: dict[str, float] = {
    "x": -1.43029175,
    "xdot": 0.00018678,
    "ydot": 0.67262261,
}

#: Gate tolerance on the Table-3 (x, xdot) match. Re-verified directly (this
#: task, #759) against the paper's own text (line ~1591 of the .txt sidecar):
#: Table 3, like Table 2, was "computed as before by using interpolation
#: between the closest points on the invariant manifolds" -- NOT a
#: Newton-corrected exact intersection -- so the SAME tolerance class as
#: :data:`TABLE2_GATE_ABS_TOL` applies, for the same reason (chasing tighter
#: agreement than the source's own interpolation-limited digits would not be
#: an honest gate).
TABLE3_GATE_ABS_TOL = 1e-4

#: Heteroclinic-connection ghost-guard radius (nondimensional, (x, xdot) Euclidean
#: norm) -- see the module's #759 docstring update for the full reasoning. MUCH
#: tighter than :data:`GHOST_GUARD_DELTA` (1e-3, calibrated for the Table-2
#: HOMOCLINIC self-connection, where the genuine intersection sits ~0.146 away
#: from 3:4-LO's own section point -- a huge margin with room for a coarse
#: guard). Table 3 connects TWO DIFFERENT orbits, and the paper's own text
#: (p.191) explicitly selected an intersection "near the 3:4 orbit" -- #757's
#: own scoping note already found 3:4-LO's own IC sits only ~1.16e-4 from the
#: PUBLISHED Table-3 x, so landing close to 3:4-LO's own point is the EXPECTED
#: geometry here, not a red flag, and a 1e-3-scale guard would risk wrongly
#: rejecting the genuine answer. This guard instead targets the much smaller,
#: genuinely-degenerate case: a "converged" crossing indistinguishable (to a
#: few manifold-epsilon) from EITHER orbit's own unperturbed section point --
#: i.e. the corrector found a trivial near-zero-net-propagation solution (the
#: seed's own epsilon offset, not a real multi-crossing manifold excursion).
#: Set to 4x the paper's own manifold offset (:data:`ANDERSON_LO_EPSILON` =
#: 0.5e-5 -> 2e-5), comfortably below the ~1.16e-4 expected genuine margin to
#: 3:4-LO's own point.
HETERO_GHOST_GUARD_DELTA = 4.0 * ANDERSON_LO_EPSILON


def ydot_from_section_eq7(system: cr3bp.CR3BPSystem, x: float, xdot: float, jacobi: float) -> float:
    """Recover ``ydot`` on the {y=0} section from ``(x, xdot)`` and the Jacobi
    constant via Anderson & Lo's own Eq. 7 (p.171), the ``ydot>0`` branch.

    Derived directly from :func:`cyclerfinder.core.cr3bp.jacobi_constant`'s own
    formula (``C = x^2 + y^2 + 2(1-mu)/r1 + 2*mu/r2 - (xdot^2+ydot^2)``) evaluated
    at ``y=0``, solved for ``ydot``:
    ``ydot = +sqrt(x^2 + 2(1-mu)/r1 + 2*mu/r2 - C - xdot^2)``. Raises ``ValueError``
    if the radicand is negative (the requested ``(x, xdot, C)`` triple is not
    physically reachable on this section).
    """
    mu = system.mu
    r1, r2 = cr3bp._r1_r2(x, 0.0, 0.0, mu)
    radicand = x * x + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - jacobi - xdot * xdot
    if radicand < 0.0:
        raise ValueError(
            f"ydot_from_section_eq7: negative radicand ({radicand:.3e}) at "
            f"x={x}, xdot={xdot}, C={jacobi} -- not a physically reachable section point"
        )
    return math.sqrt(radicand)


@dataclass(frozen=True)
class ResonantNode:
    """A resonant periodic orbit (e.g. 3:4-LO) serving as a heteroclinic-cycle node.

    Same six fields as
    :class:`cyclerfinder.genome.heteroclinic_cycle.LyapunovNode` (structurally
    satisfies that module's ``ConnectionNode`` Protocol) so
    :func:`~cyclerfinder.genome.heteroclinic_cycle.correct_connection` accepts it
    unchanged -- no new solver code, per `#752`'s own recommendation not to build a
    parallel system.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    unstable_eigvec: NDArray[np.float64]
    stable_eigvec: NDArray[np.float64]
    converged: bool

    @classmethod
    def from_candidate(
        cls,
        system: cr3bp.CR3BPSystem,
        cand: ResonantFamilyCandidate,
        *,
        eigenvalue_rel_tol: float = 1e-6,
    ) -> ResonantNode:
        """Build a node from a confirmed :class:`ResonantFamilyCandidate`.

        Recomputes the saddle Floquet pair via
        :func:`~cyclerfinder.genome.heteroclinic_cycle._planar_floquet_pair` from
        the candidate's own stored ``(x0, ydot0, period)`` -- no schema change to
        ``ResonantFamilyCandidate`` needed, the monodromy is fully determined by
        the stored IC + period.

        Raises ``ValueError`` if ``cand.is_real_unstable`` is ``False`` (the
        ``_planar_floquet_pair`` magnitude-only convention is only valid for a
        genuine positive-real saddle pair -- see that function's own docstring
        caveat) or if the recomputed unstable eigenvalue disagrees with the
        candidate's own stored ``max_eigenvalue`` by more than
        ``eigenvalue_rel_tol`` relative (a stale/mismatched candidate should fail
        loudly, not silently build a node from inconsistent data).
        """
        if not cand.is_real_unstable:
            raise ValueError(
                f"{cand.label}: not a real-saddle candidate (is_real_unstable=False, "
                f"barden_eigenvalue={cand.barden_eigenvalue!r}) -- ResonantNode requires "
                "a positive-real saddle pair (the _planar_floquet_pair magnitude-only "
                "convention is invalid for a complex/marginal eigenpair)"
            )
        state0 = np.array([cand.x0, 0.0, 0.0, 0.0, cand.ydot0, 0.0], dtype=np.float64)
        lam_u, v_u, _lam_s, v_s = _planar_floquet_pair(system, state0, cand.period)
        rel_err = abs(lam_u - cand.max_eigenvalue) / cand.max_eigenvalue
        if rel_err > eigenvalue_rel_tol:
            raise ValueError(
                f"{cand.label}: recomputed _planar_floquet_pair unstable eigenvalue "
                f"{lam_u} disagrees with the stored max_eigenvalue "
                f"{cand.max_eigenvalue} by {rel_err:.3e} relative "
                f"(> {eigenvalue_rel_tol:.0e}) -- candidate data may be stale"
            )
        return cls(
            label=cand.label,
            state0=state0,
            period=cand.period,
            jacobi=cand.jacobi,
            unstable_eigvec=v_u,
            stable_eigvec=v_s,
            converged=True,
        )


def own_section_points(
    system: cr3bp.CR3BPSystem,
    node: ResonantNode,
    *,
    ydot_sign: int = SECTION_YDOT_SIGN,
    x_sign: int = SECTION_X_SIGN,
    section_y: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[NDArray[np.float64]]:
    """The orbit's own qualifying ``{y=section_y}`` section crossings over one period.

    Used as the homoclinic ghost-guard reference set (a converged A=B connection
    landing within :data:`GHOST_GUARD_DELTA` of one of these points is the orbit
    trivially shadowing itself, not a genuine homoclinic intersection) and for
    reporting. The node's own IC ``state0`` is itself a qualifying section point by
    construction (a symmetric periodic orbit's IC is a perpendicular {y=0}
    crossing) -- included explicitly since ``solve_ivp`` events do not reliably
    register the t=0 boundary.
    """
    points: list[NDArray[np.float64]] = []
    x0, xdot0, ydot0 = float(node.state0[0]), float(node.state0[3]), float(node.state0[4])
    if int(np.sign(ydot0)) == ydot_sign and int(np.sign(x0)) == x_sign:
        points.append(np.array([x0, xdot0], dtype=np.float64))

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1] - section_y)

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, node.period),
        np.asarray(node.state0, dtype=np.float64),
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
        max_step=node.period / 500.0,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    t_floor = 1e-6 * node.period
    for t_ev, y_ev in zip(t_events, y_events, strict=False):
        if abs(float(t_ev)) <= t_floor or abs(float(t_ev) - node.period) <= t_floor:
            continue  # skip the t~0/t~period boundary duplicates of the IC itself
        if int(np.sign(float(y_ev[4]))) != ydot_sign:
            continue
        if int(np.sign(float(y_ev[0]))) != x_sign:
            continue
        points.append(np.array([float(y_ev[0]), float(y_ev[3])], dtype=np.float64))
    return points


def _ghost_distance(crossing_xv: NDArray[np.float64], own_pts: list[NDArray[np.float64]]) -> float:
    """Min Euclidean (x, xdot) distance from ``crossing_xv`` to the orbit's own points."""
    if not own_pts:
        return float("inf")
    return float(min(float(np.linalg.norm(crossing_xv - p)) for p in own_pts))


@dataclass(frozen=True)
class HomoclinicCandidate:
    """One surviving (ghost-guard-passed, converged) homoclinic self-connection."""

    connection: HeteroclinicConnection
    ghost_distance: float
    dist_to_table2: float


def find_homoclinic(
    system: cr3bp.CR3BPSystem,
    node: ResonantNode,
    *,
    epsilon: float = ANDERSON_LO_EPSILON,
    ydot_sign: int = SECTION_YDOT_SIGN,
    x_sign: int = SECTION_X_SIGN,
    ghost_guard_delta: float = GHOST_GUARD_DELTA,
    max_time_factor: float = 3.0,
    k_range: range = range(1, 7),
    branches: tuple[int, ...] = (+1, -1),
    tol: float = 1e-7,
    scan_n: int = 12,
    max_iter: int = 40,
    target: NDArray[np.float64] | None = None,
    rank_by_residual: bool = False,
) -> list[HomoclinicCandidate]:
    """Coarse scan over ``(branch_u, branch_s, k_u, k_s)`` for Wu(node) ∩ Ws(node).

    A homoclinic self-connection: ``correct_connection(system, node, node, ...)``
    is *structurally* already valid (the Jacobi-equality guard passes trivially for
    A=B) -- the paper gives no ``(branch, k)`` for the specific published
    intersection, so this scans the natural discrete space (branch signs, low
    crossing indices) documented in
    docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md Sec. 5 item 2.

    Every converged crossing is checked against the ghost guard
    (:func:`own_section_points` + :data:`GHOST_GUARD_DELTA`) BEFORE being kept --
    a converged "connection" landing on the orbit's own section point is rejected
    as the trivial non-connection, not reported as a hit. Returns ALL surviving
    candidates.

    ``target``/``rank_by_residual`` (#766): the original (Anderson & Lo 2011,
    Table 2) use of this function ranks surviving candidates by distance to a
    PUBLISHED state (:data:`TABLE2_STATE`, the default when ``target`` is
    ``None`` and ``rank_by_residual`` is ``False`` -- unchanged behaviour for
    every pre-existing caller). A caller scanning a DIFFERENT node/energy with
    no published homoclinic state to compare against (e.g. Kumar 2021's own
    seed orbit at ``C=3.0041`` -- the paper reports no such state) should pass
    ``rank_by_residual=True`` to rank by Newton-residual tightness instead
    (the only honest ranking criterion when there is nothing published to
    target); ``dist_to_table2`` is then reported as ``nan`` on every returned
    candidate rather than a distance to an irrelevant target. Passing an
    explicit ``target`` (a different published/derived 2-vector) ranks by
    distance to THAT target instead of :data:`TABLE2_STATE`.
    """
    own_pts = own_section_points(system, node, ydot_sign=ydot_sign, x_sign=x_sign)
    target_arr = (
        target
        if target is not None
        else np.array([TABLE2_STATE["x"], TABLE2_STATE["xdot"]], dtype=np.float64)
    )
    out: list[HomoclinicCandidate] = []
    for branch_u, branch_s, k_u, k_s in product(branches, branches, k_range, k_range):
        conn = correct_connection(
            system,
            node,
            node,
            k_u=k_u,
            k_s=k_s,
            epsilon=epsilon,
            branch_u=branch_u,
            branch_s=branch_s,
            ydot_sign_u=ydot_sign,
            ydot_sign_s=ydot_sign,
            x_sign_u=x_sign,
            x_sign_s=x_sign,
            max_time_factor=max_time_factor,
            scan_n=scan_n,
            tol=tol,
            max_iter=max_iter,
        )
        if not conn.converged:
            continue
        d_ghost = _ghost_distance(conn.crossing_xv, own_pts)
        if d_ghost < ghost_guard_delta:
            continue  # trivial self-shadow: orbit's own section point, not a genuine intersection
        d_t2 = (
            float("nan")
            if rank_by_residual
            else float(np.linalg.norm(conn.crossing_xv - target_arr))
        )
        out.append(
            HomoclinicCandidate(connection=conn, ghost_distance=d_ghost, dist_to_table2=d_t2)
        )
    if rank_by_residual:
        out.sort(key=lambda h: h.connection.residual)
    else:
        out.sort(key=lambda h: h.dist_to_table2)
    return out


@dataclass(frozen=True)
class Table2GateResult:
    """Result of gating a homoclinic candidate against Anderson & Lo's Table 2."""

    candidate: HomoclinicCandidate | None
    x: float
    xdot: float
    ydot_recovered: float
    x_err: float
    xdot_err: float
    ydot_err: float
    x_passed: bool
    xdot_passed: bool
    ydot_passed: bool
    passed: bool
    notes: str = ""


def gate_table2(
    system: cr3bp.CR3BPSystem,
    candidates: list[HomoclinicCandidate],
    *,
    abs_tol: float = TABLE2_GATE_ABS_TOL,
) -> Table2GateResult:
    """Gate the best (closest-to-Table-2) surviving candidate against the paper's
    own published state.

    The real gate is on ``(x, xdot)`` -- ``y=0`` is definitional to the section and
    ``ydot`` is Jacobi-redundant (Eq. 7); ``ydot_passed`` is reported for honesty
    but does NOT gate ``passed`` on its own (see module docstring). Honest FAIL if
    ``candidates`` is empty or the closest candidate misses tolerance -- never
    fudged, never loosened.
    """
    if not candidates:
        return Table2GateResult(
            candidate=None,
            x=float("nan"),
            xdot=float("nan"),
            ydot_recovered=float("nan"),
            x_err=float("inf"),
            xdot_err=float("inf"),
            ydot_err=float("inf"),
            x_passed=False,
            xdot_passed=False,
            ydot_passed=False,
            passed=False,
            notes="no surviving (converged, ghost-guard-passed) homoclinic candidate found",
        )
    best = candidates[0]
    x, xdot = float(best.connection.crossing_xv[0]), float(best.connection.crossing_xv[1])
    ydot_rec = ydot_from_section_eq7(system, x, xdot, best.connection.jacobi)
    x_err = abs(x - TABLE2_STATE["x"])
    xdot_err = abs(xdot - TABLE2_STATE["xdot"])
    ydot_err = abs(ydot_rec - TABLE2_STATE["ydot"])
    x_passed = x_err <= abs_tol
    xdot_passed = xdot_err <= abs_tol
    ydot_passed = ydot_err <= abs_tol
    passed = x_passed and xdot_passed
    return Table2GateResult(
        candidate=best,
        x=x,
        xdot=xdot,
        ydot_recovered=ydot_rec,
        x_err=x_err,
        xdot_err=xdot_err,
        ydot_err=ydot_err,
        x_passed=x_passed,
        xdot_passed=xdot_passed,
        ydot_passed=ydot_passed,
        passed=passed,
        notes="" if passed else "closest surviving candidate misses (x, xdot) tolerance",
    )


def build_3_4_lo_node(
    system: cr3bp.CR3BPSystem | None = None,
) -> tuple[cr3bp.CR3BPSystem, ResonantNode]:
    """Convenience: recover the confirmed `3:4-LO` candidate and build its node."""
    sys_ = system if system is not None else jupiter_europa_system()
    cand = recover_table1_candidate("3:4-LO", sys_)
    node = ResonantNode.from_candidate(sys_, cand)
    return sys_, node


def run_table2_gate(
    system: cr3bp.CR3BPSystem | None = None, **scan_kwargs: object
) -> Table2GateResult:
    """End-to-end: build the 3:4-LO node, scan for homoclinic candidates, gate Table 2."""
    sys_, node = build_3_4_lo_node(system)
    candidates = find_homoclinic(sys_, node, **scan_kwargs)  # type: ignore[arg-type]
    return gate_table2(sys_, candidates)


def build_34lo_kumar_c_node(
    system: cr3bp.CR3BPSystem | None = None,
) -> tuple[cr3bp.CR3BPSystem, ResonantNode, ResonantFamilyCandidate]:
    """#766: build a node from the confirmed 3:4-LO family's member AT Kumar et al.
    2021's own seed energy (``KUMAR_2021_C`` = 3.0041), via `#761`'s
    :func:`~cyclerfinder.search.jovian_resonant_families.continue_34lo_to_kumar_c`.

    Same pattern as :func:`build_3_4_lo_node`/:func:`build_5_6_lo_node`, but the
    candidate comes from a continuation endpoint rather than a hardcoded seed --
    ``continue_34lo_to_kumar_c`` re-runs its own ~8s continuation gauntlet each
    call (closure/period/equilibrium/Jacobi/Radau/fold checks all inherited for
    free), so this is NOT a cached/hardcoded shortcut: every call is a fresh,
    independently-gauntlet-validated re-derivation of the C=3.0041 member.

    Returns ``(system, node, endpoint_candidate)`` -- the candidate is returned
    alongside the node since it carries diagnostic fields (e.g.
    ``max_eigenvalue``, ``period_over_2pi``) not stored on ``ResonantNode`` itself,
    useful for reporting (e.g.
    :func:`~cyclerfinder.search.jovian_resonant_families.europa_closest_approach`
    needs the raw ``(x0, ydot0, period)``).
    """
    sys_ = system if system is not None else jupiter_europa_system()
    _branch, endpoint = continue_34lo_to_kumar_c(sys_)
    node = ResonantNode.from_candidate(sys_, endpoint)
    return sys_, node, endpoint


# ---------------------------------------------------------------------------
# #766: self-consistency evidence for a homoclinic candidate at an energy with
# NO published Table-2-style state to gate against (Kumar 2021's own C=3.0041
# seed orbit -- the paper reports no homoclinic connection state for its own
# orbit). Unlike gate_table2/gate_table3, this is honestly NOT a reproduction
# claim: it corroborates that a found candidate is a GENUINE homoclinic point
# (not a numerical artifact) via an independent Radau cross-check (reusing
# assemble_cycle/crosscheck_cycle, exactly as #754's own Table-2 gate does)
# plus a forward/backward re-approach check specific to this module: the
# found intersection state, propagated BACKWARD by the unstable leg's own
# elapsed transit time, must reproduce the epsilon-scale unstable-manifold
# seed it was found from (and symmetrically forward for the stable leg) --
# the numerical evidence that this state genuinely lies on both manifolds,
# not merely a close pass. #754's own results note did the same check with a
# fixed "3 orbital periods" horizon because its k=3 crossing developed
# quickly (strong instability, |lambda|~1036); at C=3.0041 (|lambda|~54.6,
# ~19x weaker) the manifold takes far longer to grow away from the orbit, so
# a fixed period count is NOT the right horizon here -- the elapsed transit
# time is read back from the connection's own manifold-seed derivation
# instead of guessed.
# ---------------------------------------------------------------------------


def _full_state_crossing(
    system: cr3bp.CR3BPSystem,
    seed: NDArray[np.float64],
    *,
    direction: str,
    k: int,
    max_time: float,
    ydot_sign: int | None,
    x_sign: int | None,
    section_y: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, NDArray[np.float64]] | None:
    """Like ``heteroclinic_cycle._section_crossing`` but returns
    ``(elapsed_time, full_6vector)`` at the k-th qualifying crossing instead of
    just the ``(x, xdot)`` section-plane point -- needed for the #766
    forward/backward re-approach check below, which must re-propagate the
    FULL state for exactly the time already elapsed reaching it (not just
    know where it crossed the section).
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    horizon = abs(float(max_time))
    t_span = (0.0, horizon) if direction == "unstable" else (0.0, -horizon)

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1] - section_y)

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        np.asarray(seed, float),
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
        max_step=horizon / 1000.0,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    t_floor = 1e-6 * horizon
    count = 0
    for t_ev, y_ev in zip(t_events, y_events, strict=False):
        if abs(float(t_ev)) <= t_floor:
            continue
        if ydot_sign is not None and int(np.sign(float(y_ev[4]))) != ydot_sign:
            continue
        if x_sign is not None and int(np.sign(float(y_ev[0]))) != x_sign:
            continue
        count += 1
        if count == k:
            return float(t_ev), np.asarray(y_ev, dtype=np.float64)
    return None


@dataclass(frozen=True)
class HomoclinicReapproachResult:
    """#766 forward/backward re-approach self-consistency evidence.

    ``backward_distance``: distance between (the found intersection state,
    propagated BACKWARD by the unstable leg's own elapsed transit time
    ``t_u``) and the ORIGINAL epsilon-scale unstable-manifold seed at
    ``tau_u`` -- should be small (genuine homoclinic behavior: retracing the
    unstable leg backward returns to where it started, near the orbit).
    ``forward_distance`` is the symmetric check on the stable leg (propagate
    forward by ``|t_s|``, compare to the stable-manifold seed at ``tau_s``).
    Both distances are expected to be small relative to the O(1) trajectory
    scale but NOT epsilon-scale themselves: re-propagating a hyperbolic
    trajectory over its own many-period growth accumulates roundoff amplified
    by the same Floquet factor that grew the original epsilon perturbation,
    so exact machine-precision reproduction is not expected (see the results
    note for the specific numeric interpretation).
    """

    t_u: float
    t_s: float
    backward_distance: float
    forward_distance: float


def homoclinic_reapproach_check(
    system: cr3bp.CR3BPSystem,
    node: ResonantNode,
    candidate: HomoclinicCandidate,
    *,
    epsilon: float = ANDERSON_LO_EPSILON,
    ydot_sign: int = SECTION_YDOT_SIGN,
    x_sign: int = SECTION_X_SIGN,
    max_time_factor: float = 8.0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> HomoclinicReapproachResult:
    """#766: forward/backward re-approach evidence for a homoclinic self-connection.

    Re-derives the FULL state at the found intersection (via
    :func:`_full_state_crossing`, using the candidate's own stored
    ``tau_u``/``tau_s``/``k_u``/``k_s``/``branch_u``/``branch_s`` -- not
    re-guessed), then:

    1. Propagates that state BACKWARD by exactly the unstable leg's own
       elapsed transit time ``t_u`` and compares to the ORIGINAL
       epsilon-scale unstable-manifold seed at ``tau_u`` (should be close: a
       genuine homoclinic trajectory, run backward, retraces the same path it
       was found on).
    2. Propagates FORWARD by ``|t_s|`` and compares to the stable-manifold
       seed at ``tau_s`` (the symmetric check).

    Raises ``ValueError`` if either leg fails to re-reach its own crossing
    within ``max_time_factor * node.period`` (should not happen for an
    already-converged candidate; a regression signal if it does).
    """
    conn = candidate.connection
    seed_u = _seed_on_manifold(
        system, node, tau=conn.tau_u, direction="unstable", branch=conn.branch_u, epsilon=epsilon
    )
    seed_s = _seed_on_manifold(
        system, node, tau=conn.tau_s, direction="stable", branch=conn.branch_s, epsilon=epsilon
    )
    max_time = max_time_factor * node.period
    hit_u = _full_state_crossing(
        system,
        seed_u,
        direction="unstable",
        k=conn.k_u,
        max_time=max_time,
        ydot_sign=ydot_sign,
        x_sign=x_sign,
    )
    hit_s = _full_state_crossing(
        system,
        seed_s,
        direction="stable",
        k=conn.k_s,
        max_time=max_time,
        ydot_sign=ydot_sign,
        x_sign=x_sign,
    )
    if hit_u is None or hit_s is None:
        raise ValueError(
            "homoclinic_reapproach_check: a leg did not re-reach its own crossing "
            f"within max_time_factor={max_time_factor} -- regression"
        )
    t_u, y_u = hit_u
    t_s, _y_s = hit_s
    hit_state = y_u  # matches y_s to the connection's own residual precision

    sol_back = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, -t_u),
        hit_state,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=abs(t_u) / 5000.0,
    )
    d_back = float(np.linalg.norm(sol_back.y[:, -1] - seed_u))

    sol_fwd = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, abs(t_s)),
        hit_state,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=abs(t_s) / 5000.0,
    )
    d_fwd = float(np.linalg.norm(sol_fwd.y[:, -1] - seed_s))

    return HomoclinicReapproachResult(
        t_u=t_u, t_s=t_s, backward_distance=d_back, forward_distance=d_fwd
    )


# ---------------------------------------------------------------------------
# #759: heteroclinic connection Wu(3:4-LO) ∩ Ws(5:6-LO) -- Anderson & Lo 2011
# Table 3 (p.191).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HeteroclinicCandidate:
    """One surviving (ghost-guard-passed, converged) heteroclinic connection
    Wu(node_from) -> Ws(node_to), i.e. two DIFFERENT resonant orbits."""

    connection: HeteroclinicConnection
    ghost_distance: float
    dist_to_table3: float


def find_heteroclinic(
    system: cr3bp.CR3BPSystem,
    node_from: ResonantNode,
    node_to: ResonantNode,
    *,
    epsilon: float = ANDERSON_LO_EPSILON,
    ydot_sign: int = SECTION_YDOT_SIGN,
    x_sign: int = SECTION_X_SIGN,
    ghost_guard_delta: float = HETERO_GHOST_GUARD_DELTA,
    max_time_factor: float = 3.0,
    k_range: range = range(1, 7),
    branches: tuple[int, ...] = (+1, -1),
    tol: float = 1e-7,
    scan_n: int = 12,
    max_iter: int = 40,
) -> list[HeteroclinicCandidate]:
    """Coarse scan over ``(branch_u, branch_s, k_u, k_s)`` for
    ``Wu(node_from) ∩ Ws(node_to)`` -- the heteroclinic analogue of
    :func:`find_homoclinic`, generalized to two DIFFERENT orbits
    (``node_from != node_to``, unlike the homoclinic self-connection).

    Reuses the exact same :func:`~cyclerfinder.genome.heteroclinic_cycle.correct_connection`
    machinery, the same one-sided section filters (``x_sign``, ``ydot_sign``),
    and the same discrete-scan discipline (the paper gives no ``(branch, k)``
    for the specific published intersection -- see :func:`find_homoclinic`'s
    own docstring). The ghost guard here is :data:`HETERO_GHOST_GUARD_DELTA`
    (NOT :data:`GHOST_GUARD_DELTA` -- see that constant's own docstring for
    why a homoclinic-scale guard would be wrong for this table): it checks the
    converged crossing against BOTH orbits' own qualifying section points
    (:func:`own_section_points`) and rejects only a crossing indistinguishable
    (to within the guard radius) from either orbit's own unperturbed point --
    a genuinely degenerate near-zero-propagation solution, not a real
    manifold intersection. Landing CLOSE to (but numerically distinct from)
    3:4-LO's own point is the paper's own stated expected geometry (p.191:
    "near the 3:4 orbit") and is explicitly NOT rejected by this guard at its
    calibrated radius.

    Returns ALL surviving candidates, ranked by distance to
    :data:`TABLE3_STATE` (closest first).
    """
    own_pts_from = own_section_points(system, node_from, ydot_sign=ydot_sign, x_sign=x_sign)
    own_pts_to = own_section_points(system, node_to, ydot_sign=ydot_sign, x_sign=x_sign)
    target = np.array([TABLE3_STATE["x"], TABLE3_STATE["xdot"]], dtype=np.float64)
    out: list[HeteroclinicCandidate] = []
    for branch_u, branch_s, k_u, k_s in product(branches, branches, k_range, k_range):
        conn = correct_connection(
            system,
            node_from,
            node_to,
            k_u=k_u,
            k_s=k_s,
            epsilon=epsilon,
            branch_u=branch_u,
            branch_s=branch_s,
            ydot_sign_u=ydot_sign,
            ydot_sign_s=ydot_sign,
            x_sign_u=x_sign,
            x_sign_s=x_sign,
            max_time_factor=max_time_factor,
            scan_n=scan_n,
            tol=tol,
            max_iter=max_iter,
        )
        if not conn.converged:
            continue
        d_from = _ghost_distance(conn.crossing_xv, own_pts_from)
        d_to = _ghost_distance(conn.crossing_xv, own_pts_to)
        d_ghost = min(d_from, d_to)
        if d_ghost < ghost_guard_delta:
            continue  # degenerate near-zero-propagation solution, not a genuine connection
        d_t3 = float(np.linalg.norm(conn.crossing_xv - target))
        out.append(
            HeteroclinicCandidate(connection=conn, ghost_distance=d_ghost, dist_to_table3=d_t3)
        )
    out.sort(key=lambda h: h.dist_to_table3)
    return out


@dataclass(frozen=True)
class Table3GateResult:
    """Result of gating a heteroclinic candidate against Anderson & Lo's Table 3."""

    candidate: HeteroclinicCandidate | None
    x: float
    xdot: float
    ydot_recovered: float
    x_err: float
    xdot_err: float
    ydot_err: float
    x_passed: bool
    xdot_passed: bool
    ydot_passed: bool
    passed: bool
    notes: str = ""


def gate_table3(
    system: cr3bp.CR3BPSystem,
    candidates: list[HeteroclinicCandidate],
    *,
    abs_tol: float = TABLE3_GATE_ABS_TOL,
) -> Table3GateResult:
    """Gate the best (closest-to-Table-3) surviving candidate against the paper's
    own published state. Same structure as :func:`gate_table2` -- the real gate
    is on ``(x, xdot)``; ``ydot`` is Jacobi-redundant (Eq. 7) and reported for
    honesty, not gated on its own. Honest FAIL if ``candidates`` is empty or the
    closest candidate misses tolerance -- never fudged, never loosened.
    """
    if not candidates:
        return Table3GateResult(
            candidate=None,
            x=float("nan"),
            xdot=float("nan"),
            ydot_recovered=float("nan"),
            x_err=float("inf"),
            xdot_err=float("inf"),
            ydot_err=float("inf"),
            x_passed=False,
            xdot_passed=False,
            ydot_passed=False,
            passed=False,
            notes="no surviving (converged, ghost-guard-passed) heteroclinic candidate found",
        )
    best = candidates[0]
    x, xdot = float(best.connection.crossing_xv[0]), float(best.connection.crossing_xv[1])
    ydot_rec = ydot_from_section_eq7(system, x, xdot, best.connection.jacobi)
    x_err = abs(x - TABLE3_STATE["x"])
    xdot_err = abs(xdot - TABLE3_STATE["xdot"])
    ydot_err = abs(ydot_rec - TABLE3_STATE["ydot"])
    x_passed = x_err <= abs_tol
    xdot_passed = xdot_err <= abs_tol
    ydot_passed = ydot_err <= abs_tol
    passed = x_passed and xdot_passed
    return Table3GateResult(
        candidate=best,
        x=x,
        xdot=xdot,
        ydot_recovered=ydot_rec,
        x_err=x_err,
        xdot_err=xdot_err,
        ydot_err=ydot_err,
        x_passed=x_passed,
        xdot_passed=xdot_passed,
        ydot_passed=ydot_passed,
        passed=passed,
        notes="" if passed else "closest surviving candidate misses (x, xdot) tolerance",
    )


def build_5_6_lo_node(
    system: cr3bp.CR3BPSystem | None = None,
) -> tuple[cr3bp.CR3BPSystem, ResonantNode]:
    """Convenience: recover the confirmed `5:6-LO` candidate (#758) and build its node."""
    sys_ = system if system is not None else jupiter_europa_system()
    cand = recover_758_table2_seeded_candidate(sys_)
    node = ResonantNode.from_candidate(sys_, cand)
    return sys_, node


def run_table3_gate(
    system: cr3bp.CR3BPSystem | None = None, **scan_kwargs: object
) -> Table3GateResult:
    """End-to-end: build both nodes, scan for heteroclinic candidates, gate Table 3."""
    sys_, node_34lo = build_3_4_lo_node(system)
    _sys2, node_56lo = build_5_6_lo_node(sys_)
    candidates = find_heteroclinic(sys_, node_34lo, node_56lo, **scan_kwargs)  # type: ignore[arg-type]
    return gate_table3(sys_, candidates)


__all__ = [
    "ANDERSON_LO_C_FLYBY",
    "ANDERSON_LO_EPSILON",
    "GHOST_GUARD_DELTA",
    "HETERO_GHOST_GUARD_DELTA",
    "SECTION_X_SIGN",
    "SECTION_YDOT_SIGN",
    "TABLE2_GATE_ABS_TOL",
    "TABLE2_STATE",
    "TABLE3_GATE_ABS_TOL",
    "TABLE3_STATE",
    "HeteroclinicCandidate",
    "HomoclinicCandidate",
    "HomoclinicReapproachResult",
    "ResonantNode",
    "Table2GateResult",
    "Table3GateResult",
    "build_3_4_lo_node",
    "build_5_6_lo_node",
    "build_34lo_kumar_c_node",
    "find_heteroclinic",
    "find_homoclinic",
    "gate_table2",
    "gate_table3",
    "homoclinic_reapproach_check",
    "jupiter_europa_system",
    "own_section_points",
    "run_table2_gate",
    "run_table3_gate",
    "ydot_from_section_eq7",
]
