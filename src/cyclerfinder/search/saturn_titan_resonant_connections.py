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
from itertools import product

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
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
    """

    x0: float
    xdot0: float
    residual: float
    converged: bool
    n_iter: int
    n_events_seed: int
    n_events_after_first_step: int
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
    """
    horizon = t_target + 0.5 * node.period
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
                notes="map step failed (radicand<0, crossing index exceeded, or singular section)",
            )
        residual, jac, _t_cross, n_events = out
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
                notes="singular Jacobian",
            )
        alpha = 1.0
        improved = False
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
            if out_t is not None and float(np.linalg.norm(out_t[0])) < residual_norm:
                x, xdot = x_t, xdot_t
                improved = True
                break
            alpha *= 0.5
        if not improved:
            return ChainClosureResult(
                x0=x,
                xdot0=xdot,
                residual=residual_norm,
                converged=False,
                n_iter=n_iter,
                n_events_seed=n_events_seed,
                n_events_after_first_step=n_events_after_first_step,
                notes="line search exhausted -- no improving step found",
            )
    return ChainClosureResult(
        x0=x,
        xdot0=xdot,
        residual=residual_norm,
        converged=False,
        n_iter=n_iter,
        n_events_seed=n_events_seed,
        n_events_after_first_step=n_events_after_first_step,
        notes=f"max_iter={max_iter} exhausted without converging",
    )


__all__ = [
    "EPSILON",
    "GHOST_GUARD_DELTA",
    "SECTION_YDOT_SIGN",
    "ChainClosureResult",
    "ChainProximityCandidate",
    "HomoclinicCandidate",
    "attempt_chain_closure",
    "build_34_node",
    "find_homoclinic",
    "homoclinic_reapproach_check",
    "own_section_points",
    "rank_by_proximity_to_65",
    "resonant_chain_target_point",
]
