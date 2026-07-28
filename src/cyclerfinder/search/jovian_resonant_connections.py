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
    correct_connection,
)
from cyclerfinder.search.jovian_resonant_families import (
    ANDERSON_LO_C_FLYBY,
    ResonantFamilyCandidate,
    jupiter_europa_system,
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
    candidates, ranked by distance to :data:`TABLE2_STATE` (closest first).
    """
    own_pts = own_section_points(system, node, ydot_sign=ydot_sign, x_sign=x_sign)
    target = np.array([TABLE2_STATE["x"], TABLE2_STATE["xdot"]], dtype=np.float64)
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
        d_t2 = float(np.linalg.norm(conn.crossing_xv - target))
        out.append(
            HomoclinicCandidate(connection=conn, ghost_distance=d_ghost, dist_to_table2=d_t2)
        )
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


__all__ = [
    "ANDERSON_LO_C_FLYBY",
    "ANDERSON_LO_EPSILON",
    "GHOST_GUARD_DELTA",
    "SECTION_X_SIGN",
    "SECTION_YDOT_SIGN",
    "TABLE2_GATE_ABS_TOL",
    "TABLE2_STATE",
    "HomoclinicCandidate",
    "ResonantNode",
    "Table2GateResult",
    "build_3_4_lo_node",
    "find_homoclinic",
    "gate_table2",
    "jupiter_europa_system",
    "own_section_points",
    "run_table2_gate",
    "ydot_from_section_eq7",
]
