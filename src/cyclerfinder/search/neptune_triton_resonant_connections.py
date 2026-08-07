"""`#781` Neptune-Triton 4:5 resonant-orbit homoclinic self-connection.

Sibling of :mod:`cyclerfinder.search.jovian_resonant_connections` (`#754`/
`#759`/`#766`) and :mod:`cyclerfinder.search.saturn_titan_resonant_connections`
(`#767`), retargeting the SAME ``correct_connection``/``ResonantNode``/ghost-
guard machinery (reused directly, not reimplemented) to the Neptune-Triton
system's own `#776`-confirmed "4:5-saddle" resonant orbit
(:func:`~cyclerfinder.search.neptune_triton_resonant_families.recover_esm_candidate`,
Miceli & Bosanac 2026 ESM4 ``Res45+x+h``, C=2.987089791658, ``max_eigenvalue``
= -105.05423683377265) -- the `#776`/`#777` notes' own recommended first
connection-stage target, "comfortably inside the demonstrated Newton-tractable
band" of this whole task chain.

CRITICAL HONESTY FRAMING (read before trusting any number below, same
discipline as `#766`'s and `#767`'s own docstrings): unlike Anderson & Lo
2011 (a digit-grade published homoclinic state, Table 2) or even Vaquero 2013
(a figure-only published connection, Fig. 4.9), Miceli & Bosanac 2026
publishes **no homoclinic/heteroclinic connection content of any kind** for
this system -- confirmed independently THIS task by re-reading the paper's
own body text (grepped for "homoclinic"/"heteroclinic": zero hits in either
the paper or the companion Miceli 2025 dissertation) and its own Sec. 3.2
"motion primitive-informed graph" method, which connects coarse
configuration-space-adjacent ARC SEGMENTS (a discretized approximate graph
search over pre-computed manifold-arc "motion primitives"), categorically
NOT a Newton-corrected exact Poincare-section ``Wu ∩ Ws`` intersection the
way this module (and its two siblings) compute one. See the module's own
"LITERATURE NOVELTY GATE" section below for the full account, including the
``search/literature_check.py`` run -- if a genuine, ghost-guard-passed,
independently-cross-checked connection is found here, it is the first NOVEL
(not reproduction, not self-consistency-only) result of this whole
homoclinic-connection task chain: `#776`'s own note already flagged this
explicitly ("a future connection-stage task here would be genuinely
NOVEL-territory work").

SECTION CONVENTION -- WHY THIS MODULE GENERALIZES EVEN FURTHER THAN
SATURN-TITAN'S OWN ``x``-UNRESTRICTED CHOICE:

Direct inspection of the 4:5-saddle orbit's own {y=0} crossings over one
period (this task, ``rtol=atol=1e-13``) found it crosses the axis 6 times
(excluding the ``t=0``/``t=period`` boundary duplicate), at exactly TWO
perpendicular (``xdot=0``) points -- the IC itself
(``x0=-0.969056422015``, ``ydot0=-0.126767119783``, i.e. ``ydot<0``) and the
half-period return (``x=0.981105...``, ``ydot=+0.188417...``,
``half_crossings=3`` in ``neptune_triton_resonant_families.py``) -- plus TWO
non-perpendicular mirror pairs, one at ``x~=1.122044`` (``ydot<0``,
``xdot=+-0.160239``) and one at ``x~=-1.143372`` (``ydot>0``,
``xdot=+-0.153317``).

Unlike Saturn-Titan's own 3:4 orbit (whose two perpendicular points BOTH sit
at ``ydot>0`` -- only the ``x`` sign needed unioning), this orbit's own two
perpendicular points sit at OPPOSITE ``ydot`` signs (IC: ``ydot<0``;
half-period return: ``ydot>0``). Fixing ANY single ``ydot_sign`` (Saturn-
Titan's own remaining restriction) would silently drop ONE of the two -- and
if that dropped point happens to be the orbit's OWN IC, the ghost guard's
reference set would miss the single most important trivial-shadow point
entirely (confirmed directly this task: ``ydot_sign=+1`` drops the IC;
``ydot_sign=-1`` drops the half-period return -- neither single choice is
safe). This module therefore generalizes ONE STEP FURTHER than `stc`'s own
partial (``x``-only) relaxation: BOTH ``x_sign`` and ``ydot_sign`` are fully
UNRESTRICTED (``None``) in :func:`find_homoclinic`'s own search (natively
supported by
:func:`~cyclerfinder.genome.heteroclinic_cycle.correct_connection`'s
``x_sign_u``/``x_sign_s``/``ydot_sign_u``/``ydot_sign_s`` kwargs -- no change
to that module needed), and :func:`own_section_points` below unions ALL FOUR
sign combinations (``x_sign, ydot_sign in {+1,-1} x {+1,-1}``) rather than
just two, confirmed this task to recover exactly the 6 genuine crossings
(the full per-period crossing set) with no double-counting (each crossing's
own sign pair is unique to one of the four buckets).

NEGATIVE-REAL SADDLE EIGENVALUE -- PHYSICAL MEANING AND A REAL BUG FIX
(`#781`): `#776`'s own gate table records this orbit's ``max_eigenvalue`` as
``-105.05`` -- the ONLY negative entry among all ten `#776` gate rows (every
sibling orbit this task chain has ever built a node from, Jovian/Saturn-
Titan/Neptune-Triton's own Lyapunov rows alike, has a POSITIVE real
``max_eigenvalue``). This is not a transcription artifact: independently
re-derived this task via
:func:`~cyclerfinder.genome.heteroclinic_cycle._planar_floquet_pair`
directly (``lam_u=+105.05423683428756``, matching ``abs(max_eigenvalue)`` to
5e-9 relative) -- a genuinely NEGATIVE real monodromy eigenvalue: the
period-map maps the unstable eigenvector to ``-105.05 x`` itself, i.e. a
point perturbed along the unstable direction does not merely grow each
period, it also FLIPS to the opposite side of the orbit. This is still a
perfectly good real saddle (``is_real_unstable=True`` per `#776`'s own
gate) -- ``_planar_floquet_pair``'s own magnitude-only eigenvector-direction
convention is unaffected by the sign -- but it DOES matter for how a
"branch" is read: ``branch_u=+1``/``branch_u=-1`` in this module pick a sign
of the initial epsilon-offset along the SAME unstable eigenvector, and
because the underlying map itself flips sign every period, the two branches
are NOT simply "the two sides of the manifold at a fixed phase" the way they
are for a positive-eigenvalue orbit -- they are related by the map's own
period-to-period flip. This does not require any special-case code (the
existing ``correct_connection``/``_seed_on_manifold`` machinery just
integrates from a signed epsilon-offset seed, agnostic to why the sign
matters), but the "mirror pair" characterization in the results note is
therefore interpreted with this caveat, not assumed to mean the same thing
`#766`'s/`#767`'s own positive-eigenvalue mirror pairs meant.

This ALSO exposed a real (if previously unexercised) bug in
:meth:`~cyclerfinder.search.jovian_resonant_connections.ResonantNode.from_candidate`'s
own staleness guard: its ``rel_err`` was computed as
``abs(lam_u - cand.max_eigenvalue) / cand.max_eigenvalue`` -- correct only
when ``cand.max_eigenvalue`` is positive (true for every previously-built
node), but silently NEGATIVE (never triggering the guard, regardless of
actual disagreement) for a genuinely negative ``max_eigenvalue`` like this
orbit's own. `#781` fixed this in ``jovian_resonant_connections.py`` itself
(comparing ``lam_u`` against ``abs(cand.max_eigenvalue)`` -- byte-for-byte
identical arithmetic for every already-positive case, newly-correct for the
negative one) -- see that module's own updated docstring for the full
before/after account. This module's node-builder (:func:`build_45_node`)
relies on the FIXED version.

MANIFOLD OFFSET (``EPSILON``) / GHOST GUARD (``GHOST_GUARD_DELTA``): reused
VERBATIM from :mod:`jovian_resonant_connections` (``jrc.ANDERSON_LO_EPSILON``,
``jrc.GHOST_GUARD_DELTA``) -- same generic-default reasoning `stc`'s own
docstring already gives (both systems are O(1) nondim-scale CR3BP problems;
neither source publishes its own manifold-offset magnitude).

LITERATURE NOVELTY GATE (mandatory, per the `#781` dispatch note -- read
before trusting any "novel" framing in the results note):

1. Re-read Miceli & Bosanac 2026's own body text (Sec. 3.1.1 "Computing
   Periodic Orbits and Stable or Unstable Manifolds" and Sec. 3.2 "Motion
   Primitive-Informed Graph") plus the companion Miceli 2025 dissertation,
   this task: zero occurrences of "homoclinic" or "heteroclinic" in either
   document's own text layer. The paper's own manifold use is confined to
   (a) sampling arcs along a SINGLE orbit's stable/unstable manifold for
   "motion primitives" (Sec. 3.1.4) and (b) a discretized graph search that
   connects primitives whose SEGMENTS lie near each other in configuration
   space ("adding edges only between nodes representing segments of
   primitives that exist nearby in the configuration space", Sec. 3.2) --
   an approximate, coarse-grained adjacency test, not an exact Poincare-
   section manifold intersection. Neither is the object this module
   computes.
2. ``search/literature_check.py``'s own module docstring explicitly scopes
   it to CYCLER-trajectory vocabulary ("NOT FOR RAW (non-cycler) CR3BP
   PERIODIC ORBITS") -- a raw homoclinic manifold self-connection of a
   resonant orbit is further still from that scope than the raw periodic
   orbits its own docstring already excludes. Run anyway (per the dispatch
   note's own mandatory-floor instruction) with the real ``WebSearch`` tool
   wired -- see the results note for the exact queries and verdict; its
   "not-found" (if that is what it returns) is read as weak, NOT as
   dispositive, additional corroboration on top of item 1 above (the
   load-bearing evidence), never as the primary novelty grounding on its
   own.

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.genome.heteroclinic_cycle`,
:mod:`cyclerfinder.search.jovian_resonant_connections` (reused directly for
``ResonantNode``, ``_ghost_distance``, ``_full_state_crossing``,
``HomoclinicReapproachResult`` -- none of these are Jupiter-Europa-coupled),
:mod:`cyclerfinder.search.neptune_triton_resonant_families` (`#776`/`#777`,
the 4:5-saddle candidate source).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.neptune_triton_resonant_families as ntf
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    _seed_on_manifold,
    correct_connection,
)

#: This orbit's own one-sided manifold-offset magnitude -- see module
#: docstring ("MANIFOLD OFFSET"): reused from Anderson & Lo 2011 (via
#: :data:`jrc.ANDERSON_LO_EPSILON`), not a Miceli & Bosanac-sourced value.
EPSILON = jrc.ANDERSON_LO_EPSILON

#: Homoclinic trivial-solution ghost-guard radius (nondimensional, (x, xdot)
#: Euclidean norm) -- IDENTICAL value to :data:`jrc.GHOST_GUARD_DELTA`, for
#: direct comparability across the whole task chain.
GHOST_GUARD_DELTA = jrc.GHOST_GUARD_DELTA

#: The four sign-combination buckets this module's own section convention
#: unions over -- see module docstring ("SECTION CONVENTION").
_SIGN_COMBOS: tuple[tuple[int, int], ...] = tuple(product((+1, -1), (+1, -1)))


def own_section_points(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
) -> list[NDArray[np.float64]]:
    """The orbit's own qualifying ``{y=0}`` section points, ALL FOUR
    ``(x_sign, ydot_sign)`` combinations unioned -- see module docstring
    ("SECTION CONVENTION") for why this orbit needs a fuller relaxation than
    `stc.own_section_points`'s own ``x``-only union.
    """
    out: list[NDArray[np.float64]] = []
    for x_sign, ydot_sign in _SIGN_COMBOS:
        out.extend(jrc.own_section_points(system, node, ydot_sign=ydot_sign, x_sign=x_sign))
    return out


def build_45_node(
    system: cr3bp.CR3BPSystem | None = None,
) -> tuple[cr3bp.CR3BPSystem, jrc.ResonantNode, ntf.GateRow]:
    """Convenience: recover `#776`'s confirmed "4:5-saddle" candidate and
    build its node. Returns ``(system, node, gate_row)`` -- the gate row
    (from
    :func:`~cyclerfinder.search.neptune_triton_resonant_families.gate_report`,
    filtered to the "4:5-saddle" row) carries the `#776`-confirmed
    periodicity/reproduction/crosscheck triple, useful for reporting.
    """
    sys_ = system if system is not None else ntf.neptune_triton_system()
    cand = ntf.recover_esm_candidate("4:5-saddle", sys_)
    node = jrc.ResonantNode.from_candidate(sys_, cand)
    rows = [r for r in ntf.gate_report(sys_) if r.label == "4:5-saddle"]
    if not rows:
        raise ValueError("4:5-saddle row missing from gate_report -- regression")
    return sys_, node, rows[0]


@dataclass(frozen=True)
class HomoclinicCandidate:
    """One surviving (ghost-guard-passed, converged) homoclinic self-
    connection. No ``dist_to_target`` field -- like `stc`'s own
    ``HomoclinicCandidate`` (and unlike `jrc`'s Table-2-gated one), there is
    no published state to rank against here (see module docstring); every
    candidate is ranked by Newton-residual tightness (the only honest
    criterion available).
    """

    connection: HeteroclinicConnection
    ghost_distance: float


def find_homoclinic(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    epsilon: float = EPSILON,
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
    Ws(node) -- the Neptune-Triton analogue of :func:`stc.find_homoclinic`,
    using THIS module's own FULLY-unrestricted section convention (both
    ``x_sign`` and ``ydot_sign`` -- see module docstring) instead of a
    single-axis relaxation.

    Every converged crossing is checked against the ghost guard
    (:func:`own_section_points` + ``ghost_guard_delta``) before being kept.
    Returns ALL surviving candidates, ranked by Newton-residual tightness
    (never by distance to a published target -- none exists here, per the
    module docstring's literature-novelty-gate section).
    """
    own_pts = own_section_points(system, node)
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
            ydot_sign_u=None,
            ydot_sign_s=None,
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
    max_time_factor: float = 8.0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> jrc.HomoclinicReapproachResult:
    """Forward/backward re-approach self-consistency evidence -- this
    module's own analogue of :func:`stc.homoclinic_reapproach_check`
    (threading the SAME tightened ``rtol=1e-13, atol=1e-14`` integrator
    end-to-end, though this orbit's own instability, ``|lambda|~=105``, is
    far milder than Saturn-Titan's 3:4 orbit's ``|lambda|~=2129.8`` --
    closer to `#766`'s own C=3.0041 case, ``|lambda|~=54.6``, where `jrc`'s
    own un-tightened defaults were already sufficient; this module keeps the
    tighter integrator anyway, purely for direct comparability across the
    whole task chain, not because it is strictly required here).

    Re-derives the FULL state at the found intersection (via
    ``jrc._full_state_crossing``, using the candidate's own stored
    ``tau_u``/``tau_s``/``k_u``/``k_s``/``branch_u``/``branch_s``), then:

    1. Propagates that state BACKWARD by the unstable leg's own elapsed
       transit time ``t_u`` and compares to the ORIGINAL epsilon-scale
       unstable-manifold seed at ``tau_u``.
    2. Propagates FORWARD by ``|t_s|`` and compares to the stable-manifold
       seed at ``tau_s`` (the symmetric check).

    Both legs use the module's own fully-unrestricted section filters
    (``ydot_sign=None, x_sign=None``) -- see module docstring.

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
        ydot_sign=None,
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
        ydot_sign=None,
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
    "HomoclinicCandidate",
    "build_45_node",
    "find_homoclinic",
    "homoclinic_reapproach_check",
    "own_section_points",
]
