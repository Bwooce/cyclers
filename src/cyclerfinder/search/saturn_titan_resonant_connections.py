"""#767 Saturn-Titan 3:4 resonant-orbit homoclinic self-connection.

Sibling of :mod:`cyclerfinder.search.jovian_resonant_connections` (#754/#759/
#766), retargeting its exact ``correct_connection``/``ResonantNode``/ghost-
guard machinery (reused directly, not reimplemented) to the Saturn-Titan
system's own `#765`-confirmed 3:4 resonant orbit (
:func:`~cyclerfinder.search.saturn_titan_resonant_families.recover_table41_candidate`,
Vaquero 2013 Table 4.1, C=3.010000). A new sibling module -- not an extension
of ``jovian_resonant_connections.py`` -- because that module is coupled to
Jupiter-Europa-specific constants and candidate builders (``TABLE2_STATE``,
``ANDERSON_LO_EPSILON``, ``jupiter_europa_system``, ``recover_table1_candidate``,
etc.); its own truly system-agnostic pieces (``ResonantNode``,
``own_section_points``, ``correct_connection`` itself, ``_full_state_crossing``,
``_ghost_distance``) are imported and reused here directly, exactly mirroring
how `#765`'s own ``saturn_titan_resonant_families.py`` reuses
``jovian_resonant_families.converge_candidate`` rather than reimplementing it.

CRITICAL HONESTY FRAMING (read before trusting any number below, same as
`#766`'s own docstring): Vaquero 2013 Sec. 4.3.1 describes a homoclinic
connection of the 3:4 resonant orbit (Fig. 4.9) as a FIGURE ONLY -- no state
table, unlike Anderson & Lo 2011's own digit-grade Table 2. There is nothing
published to reproduce here. Every number in this module's own results note
is therefore **self-consistency** evidence -- Newton residual, ghost-guard
margin, independent Radau cross-check, forward/backward re-approach -- never
a reproduction claim.

SECTION CONVENTION -- WHY THIS MODULE DOES NOT REUSE ANDERSON & LO'S OWN
``{y=0, x<0, ydot>0}`` CHOICE:

Anderson & Lo's own one-sided section (``jrc.SECTION_X_SIGN=-1``,
``jrc.SECTION_YDOT_SIGN=+1``) is specific to THEIR 3:4-LO orbit's own
geometry (whose IC sits at negative ``x``). This task's own 3:4 orbit sits at
POSITIVE ``x`` (``x0=1.0301663``, confirmed directly by reading
``recover_table41_candidate('3:4')`` this task) -- reusing ``x_sign=-1``
verbatim would exclude the orbit's own IC from :func:`own_section_points`'s
reference set entirely, silently breaking the ghost guard (an empty
reference set makes ``_ghost_distance`` return ``inf`` unconditionally, so
NOTHING would ever be rejected as a trivial self-shadow).

Direct inspection of this orbit's own {y=0} crossings over one period (this
task) found it crosses the axis 4 times per period, at exactly TWO
perpendicular (``xdot=0``) points -- ``x0=+1.0301663`` (the IC, ``ydot>0``)
and ``x=-1.3666368`` (the corrector's own half-period return target,
``_HALF_CROSSINGS['3:4']=2`` in ``saturn_titan_resonant_families.py``, also
``ydot>0``) -- plus a non-perpendicular mirror pair at ``x=+1.1064570``
(``ydot<0``, ``xdot=+-0.1033``). Restricting to a SINGLE ``x`` sign (as
Anderson & Lo's own convention does) keeps only ONE of this orbit's own two
perpendicular reference points and starves the ``k``-index of genuine
revisits (an initial attempt at ``x_sign=+1`` only found essentially no
qualifying crossings beyond ``k=1`` within many orbital periods -- see the
results note). This module's own section is instead ``{y=0, ydot>0}`` with
``x`` UNRESTRICTED -- the natural generalization that keeps BOTH of the
orbit's own perpendicular points as ghost-guard references
(:func:`own_section_points` below unions ``jrc.own_section_points`` at
``x_sign=+1`` and ``x_sign=-1``), and gives the manifold far more frequent
qualifying returns (confirmed directly: dozens of genuine crossings within a
10-period horizon, vs. essentially none under the single-sign restriction).
``correct_connection``'s own ``x_sign_u``/``x_sign_s`` kwargs already support
``None`` (unrestricted) natively -- no change to that module was needed.

MANIFOLD OFFSET (``EPSILON``): Vaquero 2013 does not publish a manifold
perturbation magnitude for this connection (there is no Table-2-style
numerical construction to match). This module reuses Anderson & Lo's own
``0.5e-5`` value (:data:`jrc.ANDERSON_LO_EPSILON`) as a reasonable,
already-validated small-offset magnitude for this class of nondimensional
CR3BP problem (both systems' orbits are O(1) nondim scale) -- not a
paper-sourced choice, just a sensible generic default carried over.

FORWARD/BACKWARD RE-APPROACH PRECISION -- WHY THIS MODULE DOES NOT REUSE
``jrc.homoclinic_reapproach_check`` VERBATIM: that function's own internal
``_full_state_crossing`` calls are hardcoded to ``rtol=atol=1e-12`` (not
threaded through from its own public kwargs). This orbit's unstable
eigenvalue (``|lambda|~2129.8``) is the STRONGEST instability of any orbit
this task chain has built a homoclinic self-connection for (`#754`'s
``C_flyby``: 1036; `#766`'s ``C=3.0041``: 54.6) -- direct empirical testing
this task found the SAME ``1e-12``-tolerance pipeline that worked for both
prior tasks left a residual full-state discrepancy between the unstable- and
stable-leg crossings of order ``1e-6`` (vs. `#766`'s own ``<1e-8``), which
this orbit's much stronger instability (Floquet growth rate ``ln(2129.8)/
period ~ 0.293`` per nondim time unit, vs. `#766`'s ``0.158``) amplifies over
the multi-period re-propagation horizon into a materially looser
``forward_distance``. Tightening the crossing-detection integrator to
``rtol=1e-13, atol=1e-14`` (this module's own :func:`homoclinic_reapproach_check`,
threading the tolerance through end-to-end where ``jrc``'s own version does
not) recovered roughly an order of magnitude -- reported honestly, not
fudged further; see the results note for the full numeric account and the
explanation for why ``forward_distance`` and ``backward_distance`` are
asymmetric (they are not equally well-conditioned: propagating the SAME
leg's own crossing back to its own seed is close to numerically-exact time-
reversal of the identical integration, while propagating one leg's estimate
forward to independently reproduce the OTHER leg's seed is a genuine cross-
leg consistency test that inherits the full corrector-residual floor,
amplified by the orbit's own strong instability).

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.genome.heteroclinic_cycle`,
:mod:`cyclerfinder.search.jovian_resonant_connections` (reused directly for
``ResonantNode``, ``own_section_points``, ``_ghost_distance``,
``_full_state_crossing``, ``HomoclinicReapproachResult`` -- none of these are
Jupiter-Europa-coupled), :mod:`cyclerfinder.search.saturn_titan_resonant_families`
(`#765`, the 3:4 candidate source).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from itertools import product

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_multiple_shooting as cms
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.saturn_titan_resonant_families as stf
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    _seed_on_manifold,
    correct_connection,
)
from cyclerfinder.search.cr3bp_periodic import _ubar_grad_x_at_axis

#: This orbit's own one-sided manifold-offset magnitude -- see module
#: docstring ("MANIFOLD OFFSET"): reused from Anderson & Lo 2011 (via
#: :data:`jrc.ANDERSON_LO_EPSILON`), not a Vaquero-sourced value.
EPSILON = jrc.ANDERSON_LO_EPSILON

#: Section convention: {y=0, ydot>0}, x UNRESTRICTED -- see module docstring
#: ("SECTION CONVENTION") for why this differs from Anderson & Lo's own
#: single-x-sign choice.
SECTION_YDOT_SIGN = +1

#: Homoclinic trivial-solution ghost-guard radius (nondimensional, (x, xdot)
#: Euclidean norm) -- IDENTICAL value to :data:`jrc.GHOST_GUARD_DELTA`, for
#: direct comparability across the whole task chain. This task's own two
#: genuine hits (see the results note) land at 274x and 306x this radius --
#: a real, non-delicate margin, not a borderline pass.
GHOST_GUARD_DELTA = jrc.GHOST_GUARD_DELTA


def own_section_points(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    ydot_sign: int = SECTION_YDOT_SIGN,
) -> list[NDArray[np.float64]]:
    """The orbit's own qualifying ``{y=0, ydot=ydot_sign}`` section points,
    BOTH x signs (union of :func:`jrc.own_section_points` at ``x_sign=+1``
    and ``x_sign=-1``) -- see module docstring ("SECTION CONVENTION").
    """
    return jrc.own_section_points(
        system, node, ydot_sign=ydot_sign, x_sign=+1
    ) + jrc.own_section_points(system, node, ydot_sign=ydot_sign, x_sign=-1)


def build_34_node(
    system: cr3bp.CR3BPSystem | None = None,
) -> tuple[cr3bp.CR3BPSystem, jrc.ResonantNode, stf.Table41GateRow]:
    """Convenience: recover `#765`'s confirmed Table 4.1 3:4 candidate and
    build its node. Returns ``(system, node, gate_row)`` -- the gate row
    (from :func:`~cyclerfinder.search.saturn_titan_resonant_families.gate_report`,
    filtered to the 3:4 row) carries the `#765`-confirmed eigenvalue/period
    match, useful for reporting.
    """
    sys_ = system if system is not None else stf.saturn_titan_system()
    cand = stf.recover_table41_candidate("3:4", sys_)
    node = jrc.ResonantNode.from_candidate(sys_, cand)
    rows = [r for r in stf.gate_report(sys_) if r.label == "3:4"]
    if not rows:
        raise ValueError("3:4 row missing from gate_report -- regression")
    return sys_, node, rows[0]


@dataclass(frozen=True)
class HomoclinicCandidate:
    """One surviving (ghost-guard-passed, converged) homoclinic self-
    connection. No ``dist_to_target`` field -- unlike `jrc`'s own
    ``HomoclinicCandidate``, there is no published state to rank against
    here (see module docstring); every candidate is ranked by Newton-
    residual tightness (the only honest criterion available).
    """

    connection: HeteroclinicConnection
    ghost_distance: float


def find_homoclinic(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    epsilon: float = EPSILON,
    ydot_sign: int = SECTION_YDOT_SIGN,
    ghost_guard_delta: float = GHOST_GUARD_DELTA,
    max_time_factor: float = 3.0,
    k_range: range = range(1, 7),
    branches: tuple[int, ...] = (+1, -1),
    tol: float = 1e-7,
    scan_n: int = 12,
    max_iter: int = 40,
    fd_step: float = 1e-6,
) -> list[HomoclinicCandidate]:
    """Coarse scan over ``(branch_u, branch_s, k_u, k_s)`` for Wu(node) ∩
    Ws(node) -- the Saturn-Titan analogue of :func:`jrc.find_homoclinic`,
    using THIS module's own ``x``-unrestricted section convention (see
    module docstring) instead of Anderson & Lo's single-``x``-sign one.

    Every converged crossing is checked against the ghost guard
    (:func:`own_section_points` + ``ghost_guard_delta``) before being kept.
    Returns ALL surviving candidates, ranked by Newton-residual tightness
    (never by distance to a published target -- none exists here).
    """
    own_pts = own_section_points(system, node, ydot_sign=ydot_sign)
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
            x_sign_u=None,
            x_sign_s=None,
            max_time_factor=max_time_factor,
            scan_n=scan_n,
            tol=tol,
            max_iter=max_iter,
            fd_step=fd_step,
        )
        if not conn.converged:
            continue
        d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
        if d_ghost < ghost_guard_delta:
            continue
        out.append(HomoclinicCandidate(connection=conn, ghost_distance=d_ghost))
    out.sort(key=lambda h: h.connection.residual)
    return out


def homoclinic_reapproach_check(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    candidate: HomoclinicCandidate,
    *,
    epsilon: float = EPSILON,
    ydot_sign: int = SECTION_YDOT_SIGN,
    max_time_factor: float = 8.0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> jrc.HomoclinicReapproachResult:
    """Forward/backward re-approach self-consistency evidence -- this
    module's own tighter-precision analogue of :func:`jrc.homoclinic_reapproach_check`
    (see module docstring, "FORWARD/BACKWARD RE-APPROACH PRECISION", for why
    this orbit's much stronger instability needed the crossing-detection
    integrator tightened to ``rtol=1e-13, atol=1e-14`` end-to-end, unlike
    `jrc`'s own hardcoded-``1e-12`` pipeline).

    Re-derives the FULL state at the found intersection (via
    ``jrc._full_state_crossing``, using the candidate's own stored
    ``tau_u``/``tau_s``/``k_u``/``k_s``/``branch_u``/``branch_s``), then:

    1. Propagates that state BACKWARD by the unstable leg's own elapsed
       transit time ``t_u`` and compares to the ORIGINAL epsilon-scale
       unstable-manifold seed at ``tau_u``.
    2. Propagates FORWARD by ``|t_s|`` and compares to the stable-manifold
       seed at ``tau_s`` (the symmetric check).

    Raises ``ValueError`` if either leg fails to re-reach its own crossing
    within ``max_time_factor * node.period``.
    """
    conn = candidate.connection
    seed_u = _seed_on_manifold(
        system,
        node,
        tau=conn.tau_u,
        direction="unstable",
        branch=conn.branch_u,
        epsilon=epsilon,
        rtol=rtol,
        atol=atol,
    )
    seed_s = _seed_on_manifold(
        system,
        node,
        tau=conn.tau_s,
        direction="stable",
        branch=conn.branch_s,
        epsilon=epsilon,
        rtol=rtol,
        atol=atol,
    )
    max_time = max_time_factor * node.period
    hit_u = jrc._full_state_crossing(
        system,
        seed_u,
        direction="unstable",
        k=conn.k_u,
        max_time=max_time,
        ydot_sign=ydot_sign,
        x_sign=None,
        rtol=rtol,
        atol=atol,
    )
    hit_s = jrc._full_state_crossing(
        system,
        seed_s,
        direction="stable",
        k=conn.k_s,
        max_time=max_time,
        ydot_sign=ydot_sign,
        x_sign=None,
        rtol=rtol,
        atol=atol,
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

    return jrc.HomoclinicReapproachResult(
        t_u=t_u, t_s=t_s, backward_distance=d_back, forward_distance=d_fwd
    )


# ---------------------------------------------------------------------------
# `#768`: the periodic 3:4<->6:5 "resonant chain" (Vaquero 2013 Fig. 4.9-4.10).
#
# STEP 1 FINDING (read the `#768` results note in full before trusting any
# number below): re-reading Vaquero 2013 Sec. 4.3.1 (p.114-117) directly shows
# the "resonant chain" is NOT a heteroclinic connection between two different
# orbits' manifolds -- it is a HOMOCLINIC self-connection of the 3:4 orbit
# ALONE (exactly the `#767` machinery above), whose Poincare-map crossing is
# deliberately SELECTED to sit near 6:5's own fixed point ("the intersection
# on the Poincare map for the homoclinic connection is selected to be near
# the fixed point corresponding to the 6:5 resonant orbit", p.115), so that
# the resulting trajectory shadows BOTH resonances without ever needing
# Ws(6:5)/Wu(6:5) themselves. This means the chain-proximity step below only
# needs 6:5's own IC location (:func:`resonant_chain_target_point`), NOT its
# eigenvalue/manifold structure -- `#765`'s own honest 6:5 eigenvalue near-
# miss (2.34e-3, well outside its own dual-criterion gate) is therefore NOT
# blocking here: this module only reads 6:5's own (x0, ydot0), which `#765`
# itself already confirms to <0.02% relative (as tight as 3:4/L1), and which
# this task independently re-confirms via a Radau cross-check (closure
# 6.8e-11, Jacobi drift 6.1e-14 over one period -- comparable in quality to
# 3:4's own 4.7e-10/8.9e-16). See the results note for the full Step-1
# reasoning and why this is judged NOT a genuinely ambiguous call for this
# narrower, IC-only use (distinct from -- and not overriding -- `#765`'s own
# unchanged eigenvalue-gate FAIL, which remains `#769`'s own separate scope).
# ---------------------------------------------------------------------------


def resonant_chain_target_point(
    system: cr3bp.CR3BPSystem | None = None,
) -> tuple[cr3bp.CR3BPSystem, NDArray[np.float64], stf.Table41GateRow]:
    """6:5's own ``{y=0, ydot>0}`` fixed point ``(x0, 0.0)`` -- the reference
    location Vaquero 2013 Fig. 4.9 (p.115) describes selecting the 3:4
    homoclinic connection's own map crossing "near" (see module docstring,
    "`#768`", and the `#768` results note's Step 1).

    Returns ``(system, target_xv, gate_row)`` where ``gate_row`` is `#765`'s
    own full 6:5 :class:`~cyclerfinder.search.saturn_titan_resonant_families.Table41GateRow`
    (eigenvalue AND period columns), included for transparency -- NOT to
    claim the eigenvalue gate passes (it does not, per `#765`), only so a
    caller can see exactly what is and is not confirmed about this reference
    point at a glance.

    Deliberately does NOT build a :class:`jrc.ResonantNode` for 6:5 (which
    would require its own Floquet saddle pair -- the very axis `#765`'s own
    gate honestly fails) -- only the IC location is used here, per this
    module's own Step-1 finding that the chain-proximity use case needs 6:5's
    location, not its manifolds.
    """
    sys_ = system if system is not None else stf.saturn_titan_system()
    cand = stf.recover_table41_candidate("6:5", sys_)
    rows = [r for r in stf.gate_report(sys_) if r.label == "6:5"]
    if not rows:
        raise ValueError("6:5 row missing from gate_report -- regression")
    target = np.array([cand.x0, 0.0], dtype=np.float64)
    return sys_, target, rows[0]


@dataclass(frozen=True)
class ChainProximityCandidate:
    """One :func:`find_homoclinic` hit, ranked by how close its own Poincare-
    map crossing sits to 6:5's own fixed point (:func:`resonant_chain_target_point`)
    -- NOT by residual (as :func:`find_homoclinic` itself ranks), because this
    task's own job is to find Vaquero's OWN selected crossing (Fig. 4.9: "the
    intersection...is selected to be near the fixed point corresponding to
    the 6:5 resonant orbit"), not merely the tightest-converged one.
    """

    candidate: HomoclinicCandidate
    dist_to_65: float


def rank_by_proximity_to_65(
    hits: list[HomoclinicCandidate], target_65: NDArray[np.float64]
) -> list[ChainProximityCandidate]:
    """Re-rank :func:`find_homoclinic`'s own converged, ghost-guard-passed
    hits by Euclidean ``(x, xdot)`` distance to 6:5's own fixed point, closest
    first.
    """
    ranked = [
        ChainProximityCandidate(h, float(np.linalg.norm(h.connection.crossing_xv - target_65)))
        for h in hits
    ]
    ranked.sort(key=lambda r: r.dist_to_65)
    return ranked


def _chain_ydot(x: float, xdot: float, jacobi: float, mu: float, sign: float = 1.0) -> float | None:
    """``ydot`` at ``(x, y=0, z=0)`` with the given ``xdot``, from the Jacobi
    constraint -- the general (``xdot`` need not be 0) form of
    :func:`~cyclerfinder.search.cr3bp_periodic.ydot0_from_jacobi`.
    """
    r1 = abs(x + mu)
    r2 = abs(x - 1.0 + mu)
    g = x * x + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2
    rad = g - xdot * xdot - jacobi
    if rad < 0.0:
        return None
    return sign * float(np.sqrt(rad))


def _chain_state0(
    x: float, xdot: float, jacobi: float, mu: float, sign: float = 1.0
) -> tuple[NDArray[np.float64], float] | tuple[None, None]:
    ydot = _chain_ydot(x, xdot, jacobi, mu, sign)
    if ydot is None:
        return None, None
    return np.array([x, 0.0, 0.0, xdot, ydot, 0.0], dtype=np.float64), ydot


def _chain_crossings(
    system: cr3bp.CR3BPSystem,
    x: float,
    xdot: float,
    jacobi: float,
    horizon: float,
    *,
    sign: float = 1.0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> tuple[list[tuple[float, NDArray[np.float64]]], float] | None:
    """All ``{y=0}`` crossings (with STM) of ``(x, xdot)``'s own state up to
    ``horizon``, excluding the ``t~0`` boundary. Returns ``(events, ydot0)``.
    """
    state0, ydot0 = _chain_state0(x, xdot, jacobi, system.mu, sign)
    if state0 is None or ydot0 is None:
        return None
    y0 = np.concatenate([state0, np.eye(6).reshape(36)])

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]
    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom,
        (0.0, horizon),
        y0,
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
        max_step=horizon / 3000.0,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    t_floor = 1e-6 * horizon
    events = [
        (float(t_ev), np.asarray(y_ev, dtype=np.float64))
        for t_ev, y_ev in zip(t_events, y_events, strict=False)
        if abs(float(t_ev)) > t_floor
    ]
    return events, ydot0


def _chain_map_step(
    system: cr3bp.CR3BPSystem,
    x: float,
    xdot: float,
    jacobi: float,
    crossing_index: int,
    horizon: float,
    *,
    sign: float = 1.0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> tuple[NDArray[np.float64], NDArray[np.float64], float, int] | None:
    """One Poincare-map-fixed-point Newton evaluation at ``(x, xdot)``: the
    ``{y=0}`` return at the fixed ``crossing_index``-th crossing, its section-
    projected 2x2 Jacobian (STM-based, y=0-constraint-projected exactly as
    :func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`'s
    own ``dy_dx0``/``dxdot_dx0`` projection, generalised to 2 free variables),
    the crossing time, and the total qualifying-crossing count within
    ``horizon`` (diagnostic -- see the `#768` results note for why this count
    exploding is itself the honest failure signature).

    Returns ``None`` if the Jacobi radicand goes negative, ``crossing_index``
    exceeds the number of crossings found, or the section-projection is
    singular (``|ydot| < 1e-13`` at the crossing).
    """
    res = _chain_crossings(system, x, xdot, jacobi, horizon, sign=sign, rtol=rtol, atol=atol)
    if res is None:
        return None
    events, ydot0 = res
    if crossing_index > len(events):
        return None
    t_cross, y_full = events[crossing_index - 1]
    state_f = y_full[:6]
    stm = y_full[6:].reshape(6, 6)
    dydot_dx = -_ubar_grad_x_at_axis(x, system.mu) / ydot0
    dydot_dxdot = -xdot / ydot0
    dstate0_dx = np.array([1.0, 0.0, 0.0, 0.0, dydot_dx, 0.0])
    dstate0_dxdot = np.array([0.0, 0.0, 0.0, 1.0, dydot_dxdot, 0.0])
    dy_dx = stm @ dstate0_dx
    dy_dxdot = stm @ dstate0_dxdot
    f_rhs = cr3bp.cr3bp_eom(t_cross, state_f, system.mu)
    vy = float(state_f[4])
    if abs(vy) < 1e-13:
        return None

    def _proj(dy_dv: NDArray[np.float64]) -> NDArray[np.float64]:
        dt = -dy_dv[1] / vy
        return np.asarray(dy_dv + f_rhs * dt, dtype=np.float64)

    p_dx = _proj(dy_dx)
    p_dxdot = _proj(dy_dxdot)
    x_ret, xdot_ret = float(state_f[0]), float(state_f[3])
    residual = np.array([x_ret - x, xdot_ret - xdot], dtype=np.float64)
    jac = np.array([[p_dx[0], p_dxdot[0]], [p_dx[3], p_dxdot[3]]], dtype=np.float64) - np.eye(2)
    return residual, jac, t_cross, len(events)


@dataclass(frozen=True)
class ChainClosureResult:
    """Result of :func:`attempt_chain_closure` -- this task's own bounded,
    honest attempt to correct a near-6:5 homoclinic excursion of the 3:4
    orbit into a genuine NEW periodic "resonant chain" orbit (Vaquero 2013
    Fig. 4.10), mirroring her own described method ("it is necessary to
    numerically correct this path via a single shooting scheme to obtain a
    periodic orbit", p.115). See the `#768` results note for the full account
    of why this does NOT converge with a plain STM-based 2-D Poincare-map
    Newton corrector, and why that is an honest, well-characterized negative
    rather than evidence against the chain's existence (Vaquero's own thesis
    demonstrates it exists; this module's own bounded effort with a
    single-shooting method is simply not equal to the task at this
    instability/horizon combination -- see ``notes``).

    `#773` extends this with a branch-drift guard: ``t_cross`` is the final
    iterate's own elapsed time to the fixed ``crossing_index``-th ``{y=0}``
    crossing. A seed closer to the desired excursion (e.g. the `#768` near-6:5
    homoclinic candidate's own crossing point, rather than ``node``'s plain
    IC) can converge to a Newton residual near machine precision WHILE having
    silently jumped onto a genuine but dynamically-UNRELATED, much
    shorter-period orbit family (`#773`'s own finding: seeding at
    ``(x=0.91407251, xdot=-0.091737)`` "converges" to residual ``8.4e-15`` at
    ``(x0=-0.249664, xdot0=~0)`` whose OWN ``t_cross`` is ``4.176`` nondim --
    nowhere near ``t_target=110.4996`` -- a real periodic orbit, just not the
    resonant chain). The line search below rejects any trial step whose
    ``t_cross`` drifts past ``max_t_cross_drift`` of ``t_target``, even if its
    residual is lower -- a *fast-converging* result is therefore NOT
    sufficient evidence of a genuine chain closure by itself; ``t_cross``
    must also be checked against ``t_target`` (this is why every reported
    convergence in the `#773` results note explicitly quotes ``t_cross``).
    """

    x0: float
    xdot0: float
    residual: float
    converged: bool
    n_iter: int
    n_events_seed: int
    n_events_after_first_step: int
    t_cross: float
    notes: str


def attempt_chain_closure(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    t_target: float,
    x0_guess: float | None = None,
    xdot0_guess: float = 0.0,
    ydot_sign: int = SECTION_YDOT_SIGN,
    max_iter: int = 3,
    max_backtrack: int = 8,
    tol: float = 1e-9,
    rtol: float = 1e-13,
    atol: float = 1e-14,
    max_t_cross_drift: float | None = None,
) -> ChainClosureResult:
    """Attempt an STM-based 2-D Newton (Poincare-map fixed-point) correction
    of a near-``node``-IC seed into a periodic "resonant chain" orbit that
    returns to itself after the crossing nearest ``t_target`` (the total
    unstable+stable transit time of a chosen `#768` near-6:5 homoclinic
    candidate -- see the results note for how ``t_target`` is derived).

    ``x0_guess``/``xdot0_guess`` default to ``node``'s own IC (``xdot=0``) --
    the natural seed, since the chain trajectory departs and re-approaches
    close to the 3:4 orbit itself. The fixed crossing index is determined
    ONCE from this seed (nearest ``t_target``), mirroring
    :func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`'s
    own ``half_crossings``-index-pinning discipline (kept fixed across
    iterations so Newton stays on a single branch rather than re-targeting a
    different crossing every step).

    Newton uses a damped backtracking line search (accept only if residual
    norm decreases); ``converged`` is ``True`` iff the residual norm drops
    below ``tol`` within ``max_iter``. Small ``max_iter``/``max_backtrack``
    defaults are deliberate -- see the results note: the very FIRST Newton
    step already leaves the seed's local linear regime (the qualifying-
    crossing count within the horizon explodes from a few dozen to hundreds),
    so additional iterations do not change the qualitative (honest-FAIL)
    verdict, only the wall-clock cost.

    `#773` branch-drift guard (``max_t_cross_drift``, default
    ``0.5 * node.period`` if ``None``): a trial step is accepted only if BOTH
    its residual is lower AND its own ``t_cross`` (the fixed
    ``crossing_index``-th crossing's own elapsed time) stays within
    ``max_t_cross_drift`` of ``t_target``. Without this guard, a seed close
    to the desired excursion (rather than ``node``'s own plain IC) can
    Newton-converge to residual near machine precision while having silently
    jumped onto an entirely different, genuine-but-unrelated, much
    shorter-period orbit family whose ``crossing_index``-th crossing happens
    to occur nearby in (x, xdot) but at a wildly different elapsed time --
    see the class docstring and the `#773` results note for the concrete
    numeric example this guard was built to catch.
    """
    horizon = t_target + 0.5 * node.period
    drift_cap = float(max_t_cross_drift) if max_t_cross_drift is not None else 0.5 * node.period
    x0 = float(x0_guess) if x0_guess is not None else float(node.state0[0])
    xdot0 = float(xdot0_guess)

    seed = _chain_crossings(system, x0, xdot0, node.jacobi, horizon, sign=float(ydot_sign))
    if seed is None:
        raise ValueError("attempt_chain_closure: seed state has no valid Jacobi solution")
    seed_events, _ = seed
    if not seed_events:
        raise ValueError("attempt_chain_closure: seed reaches no {y=0} crossing within horizon")
    crossing_index = int(np.argmin([abs(t - t_target) for t, _ in seed_events])) + 1
    n_events_seed = len(seed_events)
    t_cross_cur = seed_events[crossing_index - 1][0]

    x, xdot = x0, xdot0
    residual_norm = float("inf")
    n_events_after_first_step = n_events_seed
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        out = _chain_map_step(
            system,
            x,
            xdot,
            node.jacobi,
            crossing_index,
            horizon,
            sign=float(ydot_sign),
            rtol=rtol,
            atol=atol,
        )
        if out is None:
            return ChainClosureResult(
                x0=x,
                xdot0=xdot,
                residual=residual_norm,
                converged=False,
                n_iter=n_iter,
                n_events_seed=n_events_seed,
                n_events_after_first_step=n_events_after_first_step,
                t_cross=t_cross_cur,
                notes="map step failed (radicand<0, crossing index exceeded, or singular section)",
            )
        residual, jac, t_cross_cur, n_events = out
        residual_norm = float(np.linalg.norm(residual))
        if n_iter == 1:
            n_events_after_first_step = n_events
        if residual_norm < tol:
            return ChainClosureResult(
                x0=x,
                xdot0=xdot,
                residual=residual_norm,
                converged=True,
                n_iter=n_iter,
                n_events_seed=n_events_seed,
                n_events_after_first_step=n_events_after_first_step,
                t_cross=t_cross_cur,
                notes="converged",
            )
        try:
            step = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            return ChainClosureResult(
                x0=x,
                xdot0=xdot,
                residual=residual_norm,
                converged=False,
                n_iter=n_iter,
                n_events_seed=n_events_seed,
                n_events_after_first_step=n_events_after_first_step,
                t_cross=t_cross_cur,
                notes="singular Jacobian",
            )
        alpha = 1.0
        improved = False
        branch_drift_blocked = False
        for _bt in range(max_backtrack):
            x_t = x + alpha * float(step[0])
            xdot_t = xdot + alpha * float(step[1])
            out_t = _chain_map_step(
                system,
                x_t,
                xdot_t,
                node.jacobi,
                crossing_index,
                horizon,
                sign=float(ydot_sign),
                rtol=rtol,
                atol=atol,
            )
            if out_t is not None:
                residual_t, _jac_t, t_cross_t, _n_events_t = out_t
                if float(np.linalg.norm(residual_t)) < residual_norm:
                    if abs(t_cross_t - t_target) <= drift_cap:
                        x, xdot = x_t, xdot_t
                        improved = True
                        break
                    branch_drift_blocked = True
            alpha *= 0.5
        if not improved:
            notes = (
                "line search exhausted -- an improving step existed but was rejected: "
                f"its t_cross drifted past max_t_cross_drift={drift_cap:.6f} of "
                f"t_target={t_target} (see #773's branch-drift guard)"
                if branch_drift_blocked
                else "line search exhausted -- no improving step found"
            )
            return ChainClosureResult(
                x0=x,
                xdot0=xdot,
                residual=residual_norm,
                converged=False,
                n_iter=n_iter,
                n_events_seed=n_events_seed,
                n_events_after_first_step=n_events_after_first_step,
                t_cross=t_cross_cur,
                notes=notes,
            )
    return ChainClosureResult(
        x0=x,
        xdot0=xdot,
        residual=residual_norm,
        converged=False,
        n_iter=n_iter,
        n_events_seed=n_events_seed,
        n_events_after_first_step=n_events_after_first_step,
        t_cross=t_cross_cur,
        notes=f"max_iter={max_iter} exhausted without converging",
    )


# ---------------------------------------------------------------------------
# `#773` fix (b): multiple shooting, reusing the project's own #687 general
# CR3BP multiple-shooting utility (:mod:`cyclerfinder.search.cr3bp_multiple_shooting`)
# directly rather than writing a new corrector. See the `#773` results note
# for the full account of why this was tried (the single-shooting scheme
# above compounds ~1.2e14 of growth over the whole ~4.2-period loop into ONE
# STM; multiple shooting keeps each segment's own growth modest, ~57x for an
# 8-way split) and why it is ALSO an honest, non-forced partial/negative (real
# but decelerating progress, not a wild divergence, not a forced convergence).
# ---------------------------------------------------------------------------


def build_chain_multi_shooting_seed(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    x0_guess: float,
    xdot0_guess: float,
    t_target: float,
    n_segments: int,
    ydot_sign: int = SECTION_YDOT_SIGN,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[list[NDArray[np.float64]], list[float]]:
    """Build an ``n_segments``-node multiple-shooting seed for
    :func:`attempt_chain_closure_multiple_shooting` by propagating the single
    6-state derived from ``(x0_guess, xdot0_guess)`` (via :func:`_chain_state0`,
    same Jacobi-slaved-``ydot`` convention as the single-shooting corrector
    above) forward for ``t_target`` nondim time, split into ``n_segments``
    equal-duration arcs.

    By construction every INTERNAL segment's own continuity defect is ~0 (all
    nodes lie on the SAME unperturbed trajectory) -- only the final
    wrap-around segment (node ``n_segments - 1`` back to node ``0``) carries
    a genuine defect. This is the honest, most natural starting seed for
    :func:`~cyclerfinder.search.cr3bp_multiple_shooting.correct_multiple_shooting`;
    see the `#773` results note for the per-segment residual breakdown that
    confirms this (all internal blocks ``~1e-9``, all discrepancy in the last
    block).

    Raises ``ValueError`` if ``(x0_guess, xdot0_guess)`` has no valid Jacobi
    solution (see :func:`_chain_state0`).
    """
    state0, _ydot0 = _chain_state0(x0_guess, xdot0_guess, node.jacobi, system.mu, float(ydot_sign))
    if state0 is None:
        raise ValueError("build_chain_multi_shooting_seed: seed state has no valid Jacobi solution")
    seg_t = t_target / n_segments
    nodes = [state0]
    cur = state0
    for _ in range(n_segments - 1):
        arc = cr3bp.propagate(system, cur, seg_t, with_stm=False, rtol=rtol, atol=atol)
        cur = arc.state_f
        nodes.append(cur)
    return nodes, [seg_t] * n_segments


def attempt_chain_closure_multiple_shooting(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    x0_guess: float,
    xdot0_guess: float,
    t_target: float,
    n_segments: int = 8,
    ydot_sign: int = SECTION_YDOT_SIGN,
    tol: float = 1e-9,
    max_iter: int = 40,
    lm: float = 1e-8,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> cms.MultipleShootingOrbit:
    """Attempt to close the resonant-chain loop via multiple shooting instead
    of the fixed-crossing-index single-shooting scheme above -- `#773` fix
    (b). Builds the seed via :func:`build_chain_multi_shooting_seed` then
    hands it directly to
    :func:`~cyclerfinder.search.cr3bp_multiple_shooting.correct_multiple_shooting`
    (`#687`'s own general full-6D-state, free-segment-duration corrector,
    reused verbatim -- not reimplemented here).

    **Honest result (see the `#773` results note for the full numeric
    account)**: from the `#768` near-6:5 homoclinic candidate's own crossing
    state (``x0_guess=0.91407251, xdot0_guess=-0.091737``), both ``n_segments=8``
    and ``n_segments=16`` make real, monotonically-decreasing progress
    (closure residual ``1.62`` -> ``~0.49``-``0.55`` over dozens to hundreds of
    damped Newton iterations) but the backtracking step size locks onto a
    small fixed value (``alpha ~ 1e-3``-``4e-3``) and the residual decreases
    at a DECELERATING (not accelerating) rate -- an honest stall, clearly not
    on track to reach ``tol`` within any practical iteration budget, not a
    forced convergence. This is a genuinely different failure signature than
    the single-shooting corrector's own abrupt "line search exhausted" halt
    (a slow decelerating crawl vs. a sharp stop), but the same qualitative,
    honest conclusion: real progress, no closure.
    """
    nodes, seg_times = build_chain_multi_shooting_seed(
        system,
        node,
        x0_guess=x0_guess,
        xdot0_guess=xdot0_guess,
        t_target=t_target,
        n_segments=n_segments,
        ydot_sign=ydot_sign,
        rtol=rtol,
        atol=atol,
    )
    return cms.correct_multiple_shooting(
        system,
        nodes,
        seg_times,
        tol=tol,
        max_iter=max_iter,
        lm=lm,
        rtol=rtol,
        atol=atol,
    )


# ---------------------------------------------------------------------------
# `#775`: genuine CONTINUATION, per `#773`'s own final recommendation ("start
# from an ALREADY-converged...nearby periodic solution...and continue in
# whatever parameter increases the instability up to this regime, rather than
# a single/multiple-shooting Newton attempt from a cold seed"). Two avenues
# were tried in good faith (see the `#775` results note for the full
# quantitative account); both are honest, well-evidenced NEGATIVES, but
# materially SHARPER and more decisive than `#773`'s own "decelerating
# crawl"/"line search exhausted" characterizations:
#
# (1) An artificial homotopy in the periodicity-map's own residual TARGET
#     (:func:`continue_chain_closure_homotopy` below): anchor at `#767`'s own
#     already-converged, ghost-guard-passed near-6:5 homoclinic candidate
#     (whose periodicity-map residual is R0, nonzero -- it is a genuine
#     Wu/Ws matching solution, not a periodicity solution), and walk a
#     parameter `s in [0, 1]` that shrinks the REQUIRED residual from R0
#     (trivially satisfied at `s=0`, no correction needed) down to 0 (genuine
#     periodicity, `s=1`) via a sequence of small, warm-started Newton
#     corrections -- reusing :func:`_chain_map_step`'s own STM-based Jacobian
#     and the SAME `#773` t_cross branch-drift guard at every micro-step, not
#     a single Newton attempt with a shrunken tolerance.
#
#     FINDING: NOT ONE step -- at ANY tested step size from `ds=0.05` down to
#     `ds~7.5e-10` (nearly 10 orders of magnitude) -- was ever accepted. This
#     is a sharper diagnosis than `#773`'s own "map sensitive at the 8th
#     significant digit" finding: it is not merely that Newton eventually
#     stalls after real progress, it is that the periodicity map has NO
#     usable step-size window at all near this seed -- steps large enough to
#     matter (>~1e-5) already produce residual changes ~100x-700x worse than
#     the local linear (Jacobian) prediction (direct numeric check, `#775`
#     results note), while steps small enough to stay inside any plausible
#     linear regime (<~1e-7) are already swamped by this problem's own
#     `rtol=atol=1e-13/1e-14` integration-noise floor. There is no scale at
#     which this specific Newton-based local continuation can make progress.
#
# (2) A Jacobi-constant (`C`) sensitivity survey of node's OWN 3:4 family
#     (reusing :func:`~cyclerfinder.search.cr3bp_continuation.continue_family`
#     directly, no new code needed -- see the results note for the full
#     table): across the ENTIRE topologically-connected range this task
#     could reach (`C in [2.990, 3.014]`; the family folds/topology-jumps at
#     `C=3.014` in the `+C` direction -- coincidentally the EXACT value of
#     Vaquero's own claimed chain-family termination bound, flagged as a
#     suggestive structural coincidence, not proven), node's own `|lambda|`
#     only varies mildly (`2129.8` down to a shallow minimum `~1300` near
#     `C~2.992`, then back up) -- nowhere in the reachable range is this
#     orbit's own instability an order of magnitude gentler than at
#     `C=3.010000`. There is no "easier nearby C" to bootstrap a
#     less-unstable converged chain from, within this family's own connected
#     branch.
# ---------------------------------------------------------------------------


class ChainHomotopyStopReason(StrEnum):
    """Why :func:`continue_chain_closure_homotopy` stopped."""

    REACHED_TARGET = "reached_target"  # s=1 AND residual < tol -- genuine closure
    STEP_FLOOR = "step_floor"  # ds shrank below ds_min without any accepted step
    MAP_STEP_FAILED = "map_step_failed"  # _chain_map_step returned None (bad radicand etc.)
    SINGULAR_JACOBIAN = "singular_jacobian"
    MAX_OUTER_STEPS = "max_outer_steps"  # bounded-effort budget exhausted (not a physical edge)


@dataclass(frozen=True)
class ChainHomotopyStep:
    """One ACCEPTED step of :func:`continue_chain_closure_homotopy`'s own
    walk (the seed itself is step 0, ``s=0.0``).
    """

    s: float
    x0: float
    xdot0: float
    residual_norm: float
    t_cross: float
    n_events: int


@dataclass(frozen=True)
class ChainHomotopyResult:
    """Result of :func:`continue_chain_closure_homotopy`. ``steps`` is the
    full ACCEPTED sequence (always includes the seed at ``s=0.0``); a
    genuinely-stalled walk has ``len(steps) == 1`` (no step at ANY tested
    step size was ever accepted -- see the `#775` results note for why this
    is exactly what happened from `#767`'s own near-6:5 candidate).
    """

    seed_x0: float
    seed_xdot0: float
    t_target: float
    crossing_index: int
    steps: tuple[ChainHomotopyStep, ...]
    converged: bool
    stop_reason: ChainHomotopyStopReason
    s_reached: float


def continue_chain_closure_homotopy(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    t_target: float,
    x0_seed: float,
    xdot0_seed: float,
    ydot_sign: int = SECTION_YDOT_SIGN,
    ds_init: float = 0.05,
    ds_min: float = 1e-9,
    max_outer_steps: int = 200,
    max_inner_iter: int = 8,
    max_backtrack: int = 8,
    tol: float = 1e-9,
    rtol: float = 1e-13,
    atol: float = 1e-14,
    max_t_cross_drift: float | None = None,
) -> ChainHomotopyResult:
    """`#775`: genuine continuation in an artificial homotopy parameter ``s``
    that shrinks the periodicity-map's own REQUIRED residual from the seed's
    own (generally nonzero) residual ``R0`` down to ``0`` (true periodicity),
    rather than a single cold-start Newton attempt at the full target (which
    `#773` already tried, both single- and multiple-shooting, and found an
    honest stall for).

    At ``s=0`` the seed trivially satisfies "residual == R0" (no correction
    needed -- this is where the walk starts, warm from ``(x0_seed,
    xdot0_seed)``, e.g. `#767`'s own already-converged near-6:5 homoclinic
    candidate). Each outer step proposes ``s_try = s + ds``, sets
    ``target_residual = (1 - s_try) * R0``, and runs a SMALL, bounded inner
    Newton loop (reusing :func:`_chain_map_step`'s own STM-based Jacobian,
    damped backtracking, and the SAME `#773` ``max_t_cross_drift`` branch-
    guard) to re-converge onto that intermediate target from the PREVIOUS
    step's own converged state. If the inner loop converges, the step is
    accepted (``ds`` grows back toward ``ds_init``, capped); if not, ``ds``
    is halved and retried. The walk stops at ``s=1`` (genuine convergence,
    ``STOP_REASON.REACHED_TARGET``) or when ``ds`` shrinks below ``ds_min``
    without EVER accepting a step at that scale (``STOP_REASON.STEP_FLOOR`` --
    a genuine fold/stall, not a bug: see the `#775` results note for the
    direct numeric evidence that this specific seed has NO usable step-size
    window at all, from ``ds=0.05`` down past ``7.5e-10``).

    ``crossing_index`` is determined ONCE from the seed (nearest
    ``t_target``) and held fixed for the entire walk -- same single-branch
    discipline as :func:`attempt_chain_closure`.

    Raises ``ValueError`` if the seed itself has no valid Jacobi solution or
    reaches no ``{y=0}`` crossing within the horizon.
    """
    horizon = t_target + 0.5 * node.period
    drift_cap = float(max_t_cross_drift) if max_t_cross_drift is not None else 0.5 * node.period

    seed = _chain_crossings(
        system,
        x0_seed,
        xdot0_seed,
        node.jacobi,
        horizon,
        sign=float(ydot_sign),
        rtol=rtol,
        atol=atol,
    )
    if seed is None:
        raise ValueError("continue_chain_closure_homotopy: seed has no valid Jacobi solution")
    seed_events, _ = seed
    if not seed_events:
        raise ValueError(
            "continue_chain_closure_homotopy: seed reaches no {y=0} crossing within horizon"
        )
    crossing_index = int(np.argmin([abs(t - t_target) for t, _ in seed_events])) + 1

    def _map(
        x: float, xdot: float
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], float, int] | None:
        return _chain_map_step(
            system,
            x,
            xdot,
            node.jacobi,
            crossing_index,
            horizon,
            sign=float(ydot_sign),
            rtol=rtol,
            atol=atol,
        )

    out0 = _map(x0_seed, xdot0_seed)
    if out0 is None:
        return ChainHomotopyResult(
            seed_x0=x0_seed,
            seed_xdot0=xdot0_seed,
            t_target=t_target,
            crossing_index=crossing_index,
            steps=(),
            converged=False,
            stop_reason=ChainHomotopyStopReason.MAP_STEP_FAILED,
            s_reached=0.0,
        )
    r0, _jac0, t_cross0, n_events0 = out0
    steps: list[ChainHomotopyStep] = [
        ChainHomotopyStep(
            s=0.0,
            x0=x0_seed,
            xdot0=xdot0_seed,
            residual_norm=float(np.linalg.norm(r0)),
            t_cross=t_cross0,
            n_events=n_events0,
        )
    ]

    x, xdot = x0_seed, xdot0_seed
    s = 0.0
    ds = ds_init
    stop_reason = ChainHomotopyStopReason.MAX_OUTER_STEPS
    n_outer = 0
    while s < 1.0 - 1e-12 and n_outer < max_outer_steps:
        n_outer += 1
        s_try = min(1.0, s + ds)
        target_res = (1.0 - s_try) * r0
        x_try, xdot_try = x, xdot
        inner_converged = False
        last: tuple[NDArray[np.float64], float, int] | None = None
        for _inner in range(max_inner_iter):
            out = _map(x_try, xdot_try)
            if out is None:
                stop_reason = ChainHomotopyStopReason.MAP_STEP_FAILED
                break
            residual, jac, t_cross_t, n_events_t = out
            mod_res = residual - target_res
            mod_norm = float(np.linalg.norm(mod_res))
            last = (residual, t_cross_t, n_events_t)
            if mod_norm < tol:
                inner_converged = True
                break
            try:
                step = np.linalg.solve(jac, -mod_res)
            except np.linalg.LinAlgError:
                stop_reason = ChainHomotopyStopReason.SINGULAR_JACOBIAN
                break
            alpha = 1.0
            improved = False
            for _bt in range(max_backtrack):
                x_cand = x_try + alpha * float(step[0])
                xdot_cand = xdot_try + alpha * float(step[1])
                out_c = _map(x_cand, xdot_cand)
                if out_c is not None:
                    residual_c, _jac_c, t_cross_c, _n_c = out_c
                    mod_c = float(np.linalg.norm(residual_c - target_res))
                    if mod_c < mod_norm and abs(t_cross_c - t_target) <= drift_cap:
                        x_try, xdot_try = x_cand, xdot_cand
                        improved = True
                        break
                alpha *= 0.5
            if not improved:
                break
        if inner_converged and last is not None:
            s = s_try
            x, xdot = x_try, xdot_try
            residual, t_cross_t, n_events_t = last
            steps.append(
                ChainHomotopyStep(
                    s=s,
                    x0=x,
                    xdot0=xdot,
                    residual_norm=float(np.linalg.norm(residual)),
                    t_cross=t_cross_t,
                    n_events=n_events_t,
                )
            )
            ds = min(ds_init, ds * 1.5)
        else:
            ds *= 0.5
            if ds < ds_min:
                stop_reason = ChainHomotopyStopReason.STEP_FLOOR
                break

    converged = s >= 1.0 - 1e-12 and steps[-1].residual_norm < tol
    if converged:
        stop_reason = ChainHomotopyStopReason.REACHED_TARGET
    return ChainHomotopyResult(
        seed_x0=x0_seed,
        seed_xdot0=xdot0_seed,
        t_target=t_target,
        crossing_index=crossing_index,
        steps=tuple(steps),
        converged=converged,
        stop_reason=stop_reason,
        s_reached=s,
    )


# ---------------------------------------------------------------------------
# `#782`: reopening `#774` with a genuinely new technique (Parker, Davis &
# Born 2010, "Chaining periodic three-body orbits in the Earth-Moon system",
# Acta Astronautica 67:623-638 -- Vaquero's own ref [69]), per the digest at
# `docs/notes/2026-08-07-parker-davis-born-2010-chaining-orbits-digest.md`:
# patchpoints at NATURAL DYNAMICAL WAYPOINTS (orbit x-axis crossings + the
# crossings nearest a theoretical connection), not uniform time subdivision
# (which `#773`'s own `n_segments=8->16` test already found does not help).
# Two avenues, full account in the `#782` results note:
#
# (a) `build_chain_natural_seed`/`attempt_chain_closure_natural_multiple_shooting`
#     -- the technique the dispatch note names explicitly. Natural {y=0}
#     crossings of `#767`/`#775`'s own near-6:5 candidate's own trajectory
#     (unequal, dynamically-placed segment durations, spanning ~1.1-11.7
#     nondim time each, vs `#773`'s own uniform ~8.5-nondim slices) fed to
#     `#687`'s existing `cr3bp_multiple_shooting.correct_multiple_shooting`
#     UNCHANGED. Genuine, better-than-`#773` progress (closure residual
#     1.028 -> 0.322 over 80 warm-started iterations, vs `#773`'s own
#     uniform-seed 1.62 -> 0.49-0.55 over hundreds) confirming the natural
#     segmentation is genuinely better-conditioned (per-segment ||STM|| ~
#     20-260, vs the single-arc 1.2e14) -- BUT the progress is INVALIDATED
#     as a periodicity result by Jacobi drift (2.9715 at iteration 80, vs
#     the required 3.010000 exactly): `#687`'s own corrector does not
#     constrain C, and the 7N-vs-6N free-variable slack is spent partly
#     sliding downhill in energy rather than toward genuine periodicity at
#     Vaquero's own energy. An honest non-result -- real progress, wrong
#     question. This also surfaced a genuine, narrowly-scoped LATENT BUG in
#     `#687`'s own `cr3bp_multiple_shooting.py` (a `ValueError` from
#     `scipy`'s own event-bracketing on a short natural segment, escaping
#     the module's own documented "propagation failures become
#     non-convergence" contract) -- fixed there directly (see that module's
#     own `_is_event_bracketing_failure`), a real, useful, narrowly-scoped
#     deliverable independent of this avenue's own non-result.
#
# (b) `find_symmetric_chain_seed`/`attempt_chain_closure_symmetric` -- the
#     result that actually closes the loop. Extending `find_homoclinic`'s
#     own scan to include EQUAL crossing-index (`k_u == k_s`) combinations
#     -- never tried by `#767`/`#773`/`#775` (their own `near65_crossing_xv`
#     fixture only ever checked mismatched `(k_u, k_s)` in `{(4,5), (5,4)}`)
#     -- surfaces a genuinely NEW, materially CLOSER (`dist_to_65=0.0145`,
#     vs the previously-used `0.094`) homoclinic self-connection that is
#     also essentially EXACTLY perpendicular (`xdot ~ -2.7e-9`) -- the
#     signature of an x-axis-SYMMETRIC connecting trajectory. Vaquero 2013
#     Fig. 4.10(a) itself shows the resonant chain as a visibly
#     x-axis-symmetric star/zigzag shape (its own "Perpendicular Crossings"
#     label marks several vertices), consistent with the chain being a
#     genuine symmetric periodic orbit -- NOT requiring `#773`'s own general
#     2-free-variable, C-unconstrained Poincare-fixed-point Newton over the
#     WHOLE ~110.5-nondim loop. Seeding this project's own EXISTING
#     `cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
#     (the SAME corrector that already recovers 3:4/6:5 themselves, `#765`
#     -- reused directly, not reimplemented) at this near-perpendicular
#     candidate, with `half_crossings=4` (one HALF of the loop, ~56 nondim
#     time instead of ~110.5), converges in 3 Newton iterations to
#     `crossing_residual=3.8e-13`, `C` held EXACTLY at `3.010000` by
#     construction (no Jacobi-drift escape hatch -- `ydot0` is re-derived
#     from the Jacobi constraint every iteration). Independently confirmed:
#     a Radau cross-check (closure `2.6e-7`, Jacobi drift over the full
#     period `2.5e-14`), Barden-vs-full-monodromy eigenvalue agreement
#     (`2.4e-8` relative, `max|eig|=4.77e7`, real saddle), a symplectically-
#     consistent full 6-eigenvalue spectrum (reciprocal pairs confirmed to
#     `<1e-5` relative), ghost-guard margin `80.9x` past `GHOST_GUARD_DELTA`,
#     and a genuine-distinctness check against `node`'s own trajectory in
#     physical `(x,y)` space (max separation `0.081` nondim, `~99,000 km` --
#     the new orbit is NOT a re-discovery of `node` itself, though it does
#     spend most of its own period within a geometrically tight tube of
#     `node`'s own path, consistent with Fig. 4.10(a)'s own nested
#     "Exterior 3:4"/"Interior 6:5" loop geometry). This RECONCILES
#     Vaquero's own "it is necessary to numerically correct this path via a
#     single shooting scheme" claim (p.115): her single shooting was very
#     likely this SAME symmetric, C-pinned, one-free-variable formulation --
#     not `#773`'s own two-free-variable, C-unconstrained, full-loop one.
#     Same words, a materially easier, well-conditioned problem; this does
#     NOT contradict `#773`'s or `#775`'s own numbers, which remain valid
#     reports on a genuinely harder formulation of the same physical
#     question. See the `#782` results note for the full verification
#     bundle and honesty framing (self-consistency evidence only -- Fig.
#     4.10 has no digit-grade state table, same caveat as `#767`'s own
#     homoclinic-connection result).
# ---------------------------------------------------------------------------


def build_chain_natural_seed(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    x0_guess: float,
    xdot0_guess: float,
    t_target: float,
    ydot_sign: int = SECTION_YDOT_SIGN,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> tuple[list[NDArray[np.float64]], list[float], int]:
    """`#782`: Parker, Davis & Born 2010's own patchpoint-selection strategy
    (see module docstring, "`#782`") -- multiple-shooting nodes placed at
    this seed trajectory's own NATURAL ``{y=0}`` x-axis crossings (unequal
    segment durations, following the trajectory's own dynamics), NOT
    `#773`'s own UNIFORM time-sliced :func:`build_chain_multi_shooting_seed`.

    Reuses :func:`_chain_crossings` directly (no new propagation code) to
    find the natural crossing list, then assembles ``crossing_index`` nodes:
    node 0 is the seed's own state at ``t=0``; nodes ``1..crossing_index-1``
    are the FULL 6-states at the first ``crossing_index-1`` natural
    crossings; the final (wrap-around) segment duration is the time from the
    ``(crossing_index-1)``-th to the ``crossing_index``-th natural crossing
    -- i.e. the corrector's own periodicity constraint
    (``phi(node[-1], seg[-1]) == node[0]``) is EXACTLY the same "does the
    ``crossing_index``-th return match the seed" condition
    :func:`attempt_chain_closure` targets, just broken into
    ``crossing_index`` naturally-placed segments instead of one.

    Returns ``(nodes, seg_times, crossing_index)``.

    Raises ``ValueError`` if the seed state has no valid Jacobi solution, or
    reaches fewer than 2 qualifying ``{y=0}`` crossings before ``t_target``
    (need at least one internal natural crossing besides the final one).
    """
    horizon = t_target + 0.5 * node.period
    res = _chain_crossings(
        system,
        x0_guess,
        xdot0_guess,
        node.jacobi,
        horizon,
        sign=float(ydot_sign),
        rtol=rtol,
        atol=atol,
    )
    if res is None:
        raise ValueError("build_chain_natural_seed: seed state has no valid Jacobi solution")
    events, _ydot0 = res
    if not events:
        raise ValueError("build_chain_natural_seed: seed reaches no {y=0} crossing within horizon")
    crossing_index = int(np.argmin([abs(t - t_target) for t, _ in events])) + 1
    if crossing_index < 2:
        raise ValueError(
            "build_chain_natural_seed: crossing_index < 2 -- need at least one internal "
            "natural crossing before the target return"
        )
    state0, _ = _chain_state0(x0_guess, xdot0_guess, node.jacobi, system.mu, float(ydot_sign))
    if state0 is None:
        raise ValueError("build_chain_natural_seed: seed state has no valid Jacobi solution")
    nodes = [state0] + [events[i][1][:6] for i in range(crossing_index - 1)]
    seg_times = [events[0][0]]
    for i in range(1, crossing_index - 1):
        seg_times.append(events[i][0] - events[i - 1][0])
    seg_times.append(events[crossing_index - 1][0] - events[crossing_index - 2][0])
    return nodes, seg_times, crossing_index


def attempt_chain_closure_natural_multiple_shooting(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    x0_guess: float,
    xdot0_guess: float,
    t_target: float,
    ydot_sign: int = SECTION_YDOT_SIGN,
    tol: float = 1e-9,
    max_iter: int = 80,
    lm: float = 1e-8,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> cms.MultipleShootingOrbit:
    """`#782` (the technique the dispatch note names explicitly): builds the
    natural-crossing seed via :func:`build_chain_natural_seed` (Parker,
    Davis & Born 2010's own patchpoint strategy -- see module docstring),
    then hands it directly to
    :func:`~cyclerfinder.search.cr3bp_multiple_shooting.correct_multiple_shooting`
    (`#687`'s own general corrector, reused verbatim -- not reimplemented).

    **Honest result (see the `#782` results note for the full numeric
    account)**: from `#767`/`#775`'s own near-6:5 candidate, this makes
    genuinely BETTER progress than `#773`'s own uniform-time-sliced seed
    (closure residual ``1.028 -> 0.322`` over 80 iterations here, vs
    `#773`'s own ``1.62 -> 0.49-0.55`` over hundreds) -- per-segment
    ``||STM||`` stays ``~20-260`` (vs the single-arc ``1.2e14``), confirming
    the natural, non-uniform segmentation is genuinely better-conditioned.
    BUT this module's own `#687` corrector does not constrain the Jacobi
    constant, and the reported progress is INVALIDATED as a periodicity
    result by Jacobi drift (``2.9715`` at iteration 80, vs the required
    ``3.010000`` exactly) -- part of the residual decrease is the corrector
    sliding downhill in energy, not converging toward genuine periodicity at
    Vaquero's own C. Callers MUST check
    ``cyclerfinder.core.cr3bp.jacobi_constant(orbit.nodes[0], system.mu)``
    against ``node.jacobi`` before treating ANY residual number from this
    function as progress toward the chain -- an honest non-result, not a
    positive (see :func:`attempt_chain_closure_symmetric` for the technique
    that actually closes the loop at the correct energy).
    """
    nodes, seg_times, _crossing_index = build_chain_natural_seed(
        system,
        node,
        x0_guess=x0_guess,
        xdot0_guess=xdot0_guess,
        t_target=t_target,
        ydot_sign=ydot_sign,
        rtol=1e-13,
        atol=1e-14,
    )
    return cms.correct_multiple_shooting(
        system, nodes, seg_times, tol=tol, max_iter=max_iter, lm=lm, rtol=rtol, atol=atol
    )


def find_symmetric_chain_seed(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    target_65: NDArray[np.float64],
    *,
    branches: tuple[int, ...] = (-1,),
    k_range: range = range(4, 6),
    perp_tol: float = 1e-3,
    epsilon: float = EPSILON,
    ydot_sign: int = SECTION_YDOT_SIGN,
    ghost_guard_delta: float = GHOST_GUARD_DELTA,
    max_time_factor: float = 3.0,
    tol: float = 1e-9,
    scan_n: int = 12,
    max_iter: int = 40,
    fd_step: float = 1e-6,
) -> ChainProximityCandidate:
    """`#782`: scan for an EQUAL-crossing-index (``k_u == k_s``) homoclinic
    self-connection of ``node`` -- a genuine on-axis (perpendicular,
    ``xdot ~ 0``) intersection, the signature of an x-axis-SYMMETRIC
    connecting trajectory (see module docstring, "`#782`"). Never tried by
    `#767`/`#773`/`#775` -- their own scans/fixtures only ever checked
    MISMATCHED ``(k_u, k_s)`` pairs like ``(4,5)``/``(5,4)``.

    Returns the single closest-to-6:5, near-perpendicular
    (``|xdot| < perp_tol``) candidate found by :func:`find_homoclinic` over
    ``branches``/``k_range`` (default: ``branches=(-1,)``, ``k_range=range(4,
    6)`` -- a narrow, cheap window that already contains the confirmed
    `#782` candidate; a wider scan is a possible but NOT required follow-up,
    see the results note).

    Raises ``ValueError`` if no near-perpendicular candidate is found in
    that window.
    """
    hits = find_homoclinic(
        system,
        node,
        epsilon=epsilon,
        ydot_sign=ydot_sign,
        ghost_guard_delta=ghost_guard_delta,
        max_time_factor=max_time_factor,
        k_range=k_range,
        branches=branches,
        tol=tol,
        scan_n=scan_n,
        max_iter=max_iter,
        fd_step=fd_step,
    )
    ranked = rank_by_proximity_to_65(hits, target_65)
    perp = [r for r in ranked if abs(float(r.candidate.connection.crossing_xv[1])) < perp_tol]
    if not perp:
        raise ValueError(
            "find_symmetric_chain_seed: no near-perpendicular (|xdot|<"
            f"{perp_tol}) homoclinic candidate found in branches={branches}, "
            f"k_range={k_range}"
        )
    perp.sort(key=lambda r: r.dist_to_65)
    return perp[0]


@dataclass(frozen=True)
class SymmetricChainClosureResult:
    """Result of :func:`attempt_chain_closure_symmetric` -- `#782`'s own
    symmetric single-shooting closure of the resonant chain (see module
    docstring, "`#782`").

    ``branch_ok`` is this task's own branch-drift guard (the symmetric-orbit
    analogue of `#773`'s ``max_t_cross_drift``, mandated by the dispatch
    note): the demonstrated basin sensitivity of this map (a perturbation of
    ``1e-4`` in ``x0_seed`` can converge to a period-``5.0`` or period-``72``
    orbit instead of this one -- see the results note) means a "converged"
    result is only trustworthy if ``x0`` did not drift far from the seed.
    ``branch_ok`` is ``False`` (even when ``converged`` is ``True``) if
    ``|x0 - x0_seed| > max_x0_drift``.
    """

    x0: float
    ydot0: float
    period: float
    jacobi: float
    converged: bool
    crossing_residual: float
    n_iter: int
    half_crossings: int
    x0_seed: float
    x0_drift: float
    branch_ok: bool
    n_crossings_per_period: int
    max_eigenvalue: float
    is_real_unstable: bool
    dist_to_65_at_seed: float
    notes: str


def attempt_chain_closure_symmetric(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    x0_seed: float,
    half_crossings: int,
    target_65: NDArray[np.float64],
    period_guess: float,
    ydot0_sign: float = 1.0,
    max_x0_drift: float = 0.01,
    tol: float = 1e-12,
    max_iter: int = 40,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> SymmetricChainClosureResult:
    """`#782`: symmetric single-shooting closure of the resonant chain (see
    module docstring, "`#782`") -- reuses
    :func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
    directly (the SAME corrector that already recovers 3:4/6:5 themselves,
    `#765`), NOT reimplemented. ``ydot0`` (hence the Jacobi constant) is
    re-derived from the Jacobi constraint every Newton iteration -- ``C`` is
    held EXACTLY at ``node.jacobi`` throughout, unlike
    :func:`attempt_chain_closure_natural_multiple_shooting`'s own
    Jacobi-drift failure mode.

    Also computes the full-period monodromy eigenvalue (real-saddle check,
    same discipline as `#765`'s own Table 4.1 gate), the natural ``{y=0}``
    crossing count over the recovered period (a coarse branch-explosion
    sanity check, in the spirit of `#773`'s own ``t_cross``/
    ``n_events_after_first_step`` diagnostics), and ``x0_drift``/
    ``branch_ok`` (see the class docstring).
    """
    dist65 = float(np.linalg.norm(np.array([x0_seed, 0.0], dtype=np.float64) - target_65))
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_seed,
        node.jacobi,
        period_guess,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        tol=tol,
        max_iter=max_iter,
        rtol=rtol,
        atol=atol,
    )
    x0_drift = abs(orbit.x0 - x0_seed)
    branch_ok = x0_drift <= max_x0_drift

    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    arc = cr3bp.propagate(system, state0, orbit.period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    eigs = np.linalg.eigvals(arc.stm)
    idx = int(np.argmax(np.abs(eigs)))
    lead = eigs[idx]
    max_eigenvalue = float(np.abs(lead))
    is_real_unstable = bool(
        abs(lead.imag) < 1e-6 * max(1.0, abs(lead.real)) and max_eigenvalue > 1.0 + 1e-9
    )

    times_x, _states_x = cp._xaxis_crossings(
        system, state0, orbit.period * 1.05, with_stm=False, rtol=rtol, atol=atol
    )
    n_crossings = int(np.sum(np.asarray(times_x) <= orbit.period * 1.0001))

    if orbit.converged and branch_ok:
        notes = "converged"
    elif orbit.converged:
        notes = "converged but branch-drift-rejected (x0_drift exceeds max_x0_drift)"
    else:
        notes = "did not converge"

    return SymmetricChainClosureResult(
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=orbit.converged,
        crossing_residual=orbit.crossing_residual,
        n_iter=orbit.n_iter,
        half_crossings=half_crossings,
        x0_seed=x0_seed,
        x0_drift=x0_drift,
        branch_ok=branch_ok,
        n_crossings_per_period=n_crossings,
        max_eigenvalue=max_eigenvalue,
        is_real_unstable=is_real_unstable,
        dist_to_65_at_seed=dist65,
        notes=notes,
    )


__all__ = [
    "EPSILON",
    "GHOST_GUARD_DELTA",
    "SECTION_YDOT_SIGN",
    "ChainClosureResult",
    "ChainHomotopyResult",
    "ChainHomotopyStep",
    "ChainHomotopyStopReason",
    "ChainProximityCandidate",
    "HomoclinicCandidate",
    "SymmetricChainClosureResult",
    "attempt_chain_closure",
    "attempt_chain_closure_multiple_shooting",
    "attempt_chain_closure_natural_multiple_shooting",
    "attempt_chain_closure_symmetric",
    "build_34_node",
    "build_chain_multi_shooting_seed",
    "build_chain_natural_seed",
    "continue_chain_closure_homotopy",
    "find_homoclinic",
    "find_symmetric_chain_seed",
    "homoclinic_reapproach_check",
    "own_section_points",
    "rank_by_proximity_to_65",
    "resonant_chain_target_point",
]
