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


__all__ = [
    "EPSILON",
    "GHOST_GUARD_DELTA",
    "SECTION_YDOT_SIGN",
    "HomoclinicCandidate",
    "build_34_node",
    "find_homoclinic",
    "homoclinic_reapproach_check",
    "own_section_points",
]
