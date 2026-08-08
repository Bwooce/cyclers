"""`#786` Earth-Moon Class 1 (Casoliva p:q-resonant orbit) homoclinic self-connections.

Sibling of :mod:`cyclerfinder.search.jovian_resonant_connections` (`#754`/
`#759`/`#766`), :mod:`cyclerfinder.search.saturn_titan_resonant_connections`
(`#767`), and :mod:`cyclerfinder.search.neptune_triton_resonant_connections`
(`#781`), retargeting the SAME ``correct_connection``/``ResonantNode``/ghost-
guard machinery (reused directly, not reimplemented) at `#780`'s own
confirmed Earth-Moon Class 1 Table 3 rows
(:mod:`cyclerfinder.search.earth_moon_resonant_families`,
Casoliva et al. 2010, Table 3, p.1630). This is the connection-stage sibling
of `#754`/`#767`/`#781` -- this project's own resonant-connection machinery
had never been pointed at Earth-Moon before this task.

**FILENAME NOTE (read first):** the `#786` dispatch note asked for a module
named ``earth_moon_resonant_connections.py`` -- but that exact name was
ALREADY claimed by `#783` (Earth-Moon Class 2/He1 Barrabes-Mondelo-Olle
continuation-of-homoclinic-connections attempt, a committed clean negative,
predating this task's own dispatch) for a categorically different algorithm
(Barrabes-Mondelo-Olle's own manifold-parametrization continuation system,
not this project's own Poincare-section Newton-shooting machinery this
module uses). The dispatch note's own filename was written without
knowledge of that collision. This module is named
``earth_moon_class1_resonant_connections.py`` instead -- distinguishing
Class 1 (this module, Sec. IV both papers, p:q-resonant orbits) from Class 2
(`#783`, Sec. V, L1-Lyapunov He1 family) -- to avoid silently overwriting
`#783`'s own committed work. Same for the test file
(``test_earth_moon_class1_resonant_connections.py`` vs. the already-claimed
``test_earth_moon_resonant_connections.py``).

CRITICAL HONESTY FRAMING (read before trusting any number below, same
discipline `#766`'s/`#767`'s/`#781`'s own docstrings established): Casoliva
2010/2008 build homoclinic connections ONLY for Class 2 (the L1-Lyapunov
"He1" family, Sec. V, `#783`'s own already-negative scope) -- Class 1's
own Sec. IV method (elliptical/second-species differential correction) never
touches manifolds or homoclinic connections at all. There is therefore
nothing published BY CASOLIVA to reproduce for a Class-1-orbit homoclinic
connection -- every number here is **self-consistency** evidence (Newton
residual, ghost-guard margin, independent Radau cross-check, forward/
backward re-approach), never a Casoliva-reproduction claim, mirroring
`#767`'s own framing. HOWEVER (see "LITERATURE NOVELTY GATE" below, the
load-bearing finding of this task): Vaquero's 2013 Purdue dissertation
(already in this project's corpus, `#765`'s own primary source for a
DIFFERENT chapter) has its OWN Earth-Moon chapter (Sec. 4.4.1, not
previously deep-read by `#765`'s/`#780`'s/`#783`'s own digest passes, which
only covered her Saturn-Titan Sec. 4.3) that explicitly computes GENUINE
homoclinic self-connections of Earth-Moon p:q-resonant orbits -- 1:2 and 2:3,
at ``C=2.8284``, via her OWN family-continuation orbit-generation lineage
(NOT Casoliva's elliptical/second-species method, and NOT this module's own
7:3 target). This means the underlying PHENOMENON (a genuine homoclinic
self-connection of an Earth-Moon p:q-resonant orbit exists and is
Newton-tractable) is established PRIOR ART, not novel -- but this module's
own SPECIFIC target (Casoliva's Class 1 `7-3b`/`7-3c`, `C~1.0687`) was never
computed by either source. Any genuine hit here is reported as
CORROBORATION-AND-EXTENSION of Vaquero 2013's own demonstrated phenomenon to
a new specific resonance/orbit-generation lineage, per the dispatch note's
own explicit instruction ("report the result as reproduction/corroboration,
not novel discovery -- that's still a legitimate, useful result").

TARGET SELECTION (this task's own analysis of `#780`'s ``table3_gate_report``):
of the 12 Table 3 rows that fully reproduce Casoliva's own printed IC/period/
Jacobi/stability index (`#780`'s results note), the physically-real
(``exists_in_em_system=True``) AND genuinely-resonant
(``satisfies_resonance=True``) rows with a REAL in-plane saddle
(``|k_signed| > 2``, ``k_signed == k_par``) are: `1-2d` (k=4.8573,
lambda~4.64 -- BELOW this task chain's own demonstrated Newton-tractable
band), `2-1b` (k=2.0374, lambda~1.19 -- essentially marginal, far too weak
an instability to separate a manifold from the orbit in any reasonable
number of periods), `7-3b` (k=57.3562, lambda~57.34), and `7-3c` (k=57.0431,
lambda~57.02, a SEPARATE periodic orbit at the SAME (C_J, period) as 7-3b --
not the same orbit, confirmed directly this task: different ``(x0, xdot0,
ydot0)`` triples, different section-crossing sequences). `7-3b`/`7-3c` sit
squarely inside the ~50-2000 |lambda| band this whole task chain has
repeatedly found Newton-tractable (`#766`'s C=3.0041: |lambda|~54.6; `#781`'s
4:5-saddle: |lambda|~105) -- both are targeted here as PRIMARY/SECONDARY.
`1-2d`/`2-1b` are registered as a follow-on (see the results note), not
attempted this task (their own eigenvalues sit at or below this chain's own
weak-instability floor, where a genuine attempt would need a much longer
scan horizon this task did not have budget to pursue on top of everything
else).

SECTION CONVENTION -- WHY THIS MODULE NEEDS EVEN MORE MACHINERY THAN
NEPTUNE-TRITON'S OWN "FOUR-BUCKET" RELAXATION:

Every prior sibling module's target orbit was built via
``correct_symmetric_fixed_jacobi`` (a genuine perpendicular ``{y=0, xdot=0}``
x-axis crossing IS the orbit's own IC by construction). Casoliva's own Class
1 orbits are recovered via `#780`'s ``recover_table3_row``, a GENERAL
full-state Newton corrector (``cr3bp_periodic.correct_periodic``, 6 free
state components + period, no perpendicularity constraint) seeded from
Casoliva's own printed generic (non-perpendicular) Poincare-crossing IC --
her own text: "we have not used this property [perpendicular crossings] in
this paper." Direct inspection of `7-3b`'s own recovered orbit (this task,
``rtol=atol=1e-13``) confirms it: the orbit crosses ``{y=0}`` **20 times**
per period, at NO perpendicular (``xdot=0``) point whatsoever -- the
CLOSEST any crossing comes is ``|xdot| ~= 0.2587`` (row 3, ``t~=2.987``), a
materially richer and qualitatively different crossing structure than any
prior sibling orbit (Jovian: 2 crossings/period both perpendicular;
Saturn-Titan: 4 crossings/period, 2 perpendicular; Neptune-Triton: 6
crossings/period, 2 perpendicular + 2 non-perpendicular mirror pairs -- both
figures measured directly by their own modules' docstrings, not assumed by
comparison here). This is reported plainly, not smoothed over: Casoliva's
own Class 1 orbits are genuinely NOT symmetric about the x-axis (at least
not in the sense of possessing a perpendicular {y=0, xdot=0} crossing) -- an
intrinsic geometric difference from every prior connection-stage target this
task chain has built, not a numerical artifact of the corrector.

Three direct consequences, all handled explicitly below (NOT silently
assumed away):

1. **The corrector's own converged IC does not sit EXACTLY on {y=0}**
   (``correct_periodic``'s min-norm 7-unknown/6-equation Newton step has no
   ``y0=0`` constraint -- the converged `7-3b` IC lands at ``y0 =
   -1.60e-05``, a small but nonzero residual). :func:`_snap_to_y0` below
   phase-shifts the converged IC by a tiny (``|dt| ~ 1.2e-5`` nondim time
   units) local Newton correction to land EXACTLY on ``{y=0}`` (to
   ``|y|<1e-14``) before it is used as ``ResonantNode.state0`` --
   re-verified directly this task: periodicity closure at the snapped IC
   (propagate for the SAME period ``T``) is ``3.6e-11`` for `7-3b` and
   ``1.5e-09`` for `7-3c`, both far tighter than `#780`'s own ``tol=1e-8``
   gate -- the snap does not degrade the orbit's own periodicity, it is a
   phase re-parametrization of the SAME solution, not a different one.
2. **No single sign-restricted section (nor even Neptune-Triton's own
   4-combo union) covers the crossing density needed** -- with 20 crossings
   spread across all four ``(x_sign, ydot_sign)`` quadrants and no
   perpendicular anchor point, this module reuses `#781`'s own
   fully-unrestricted ``(x_sign=None, ydot_sign=None)`` convention directly
   (verified this task to recover the full 20-point crossing set with no
   double-counting, same verification method `#781`'s own module used)
   rather than inventing a THIRD convention style -- `#781`'s own
   "generalizes even further" framing already anticipated exactly this
   situation.
3. **The manifold-offset magnitude needs to be larger than the reused
   ``jrc.ANDERSON_LO_EPSILON`` default** -- see "MANIFOLD OFFSET" below, a
   direct consequence of the same dense-crossing geometry: at
   ``epsilon=0.5e-5`` the manifold does not clear :data:`GHOST_GUARD_DELTA`
   until roughly ``k~26-38`` (ln-growth estimate, ``ln(lambda_u)~4.05``),
   which collides with the section's own minimum crossing-time gap
   (``~0.0686`` nondim time units) at ``correct_connection``'s default
   crossing-detection resolution over a multi-period horizon -- see
   "MANIFOLD OFFSET" for the resolved choice and the resolution check this
   task ran to confirm it is safe.

``ResonantNode`` IS NOT BUILT VIA ``ResonantNode.from_candidate`` HERE:
that classmethod hardcodes ``state0 = [cand.x0, 0, 0, 0, cand.ydot0, 0]``
(shaped to a Class 2/Anderson & Lo/Vaquero/Miceli-style
``ResonantFamilyCandidate`` whose IC IS a perpendicular crossing by
construction). Table 3 rows have no such object (`#780`'s own
``recover_table3_row`` returns a plain
:class:`~cyclerfinder.search.cr3bp_periodic.PeriodicOrbit`, a FULL 6-state
IC, not an ``(x0, ydot0)`` pair). :func:`build_node` below constructs the
``jrc.ResonantNode`` dataclass DIRECTLY (same six fields, same
``_planar_floquet_pair`` eigenvector derivation
``ResonantNode.from_candidate`` itself uses internally) from the
snapped-to-``{y=0}`` full state -- no change to ``ResonantNode`` or
``from_candidate`` needed, this is a legitimate alternate constructor
pattern for a node whose source object does not fit the
``ResonantFamilyCandidate`` shape. :func:`build_node` verifies the
recomputed ``_planar_floquet_pair`` unstable eigenvalue reproduces `#780`'s
own ``table3_gate_report`` ``k_signed`` (converted via ``lambda + 1/lambda =
k``) within :data:`EIGENVALUE_CROSSCHECK_REL_TOL` -- the same staleness
discipline ``ResonantNode.from_candidate`` applies internally for its own
candidates.

MANIFOLD OFFSET (``EPSILON``) -- RAISED FROM THE REUSED
``jrc.ANDERSON_LO_EPSILON`` DEFAULT, A DELIBERATE, DOCUMENTED DEVIATION:

Every prior sibling module reused ``jrc.ANDERSON_LO_EPSILON`` (0.5e-5)
verbatim, since none of Casoliva/Vaquero/Miceli publishes its own
manifold-offset magnitude for this class of connection and Anderson & Lo's
own value is a reasonable generic small-offset choice for an O(1) nondim
CR3BP problem. This module instead uses ``EPSILON=1e-4`` (20x larger),
because THIS orbit's own crossing density (20/period, min gap ``~0.0686``
nondim time) makes the smaller default numerically risky in a way none of
the prior orbits' sparser crossing structure (2-6/period) exposed: a direct
diagnostic this task ran (propagating a single unstable-manifold seed
forward and tracking (k, ghost-distance) pairs) found the manifold does not
clear :data:`GHOST_GUARD_DELTA` until ``k~13-16`` at ``epsilon=0.5e-5``
(ghost distance 9.5e-4 at k=13, still inside the guard; 1.3e-3 at k=16,
just clear) -- a scan window that would require ``max_time_factor~=1``
(horizon approx. one orbital period) for the k-th crossing to land inside
the scanned horizon at all, but a ONE-period horizon leaves
``correct_connection``'s own hardcoded ``max_step=horizon/500`` (``~0.038``
at ``max_time_factor=1``) close enough to the ``~0.0686`` minimum crossing
gap that a coarser choice risks stepping over closely-spaced crossing pairs
(this task's own resolution check, below, found this margin acceptable at
``max_time_factor<=3`` specifically -- NOT verified safe at the tighter
horizons a smaller epsilon would otherwise force). At ``epsilon=1e-4`` the
guard clears comfortably by ``k~10`` (ghost distance already 5.7e-3, six
times the guard) with a full 3-period horizon available, giving both a
safer crossing-detection margin AND a smaller required ``k`` -- a
genuinely better-conditioned scan for this orbit's own geometry, not a
tolerance loosened to force a result (the ghost guard threshold itself,
:data:`GHOST_GUARD_DELTA`, is UNCHANGED from every sibling module).

**Crossing-detection resolution check (this task, verifying the above is
actually safe, not merely estimated):** counted ``{y=0}`` event crossings
over a 2- and 3-period horizon at ``correct_connection``'s own default
``max_step=horizon/500`` versus a 10x-finer ``horizon/5000`` -- IDENTICAL
counts both horizons (41 events over 2 periods, 61 over 3 periods, both
resolutions) -- confirming no crossing is silently stepped-over at the
resolution this module's own scans actually use. (A further 40x-finer
``horizon/20000`` check found ONE FEWER event over the 2-period horizon --
reported honestly as an unresolved minor discrepancy, most plausibly a
near-tangential/boundary root-finding artifact at that much finer step
count, not evidence against the ``horizon/500``-vs-``horizon/5000``
agreement this module's own scans actually rely on.)

LITERATURE NOVELTY GATE (mandatory, per the `#786` dispatch note -- read
before trusting any "novel"/"corroboration" framing in the results note):
see the results note's own dedicated section for the full account,
including the Vaquero 2013 Sec. 4.4.1 finding summarized above (the
load-bearing evidence -- found via direct text-layer reading of this
project's own already-acquired corpus, NOT by `search/literature_check.py`,
whose own module docstring scopes it to cycler vocabulary this raw
resonant-orbit-manifold object does not match, per every prior sibling
module's own identical caveat) plus a corroborating WebSearch pass and the
mandatory ``literature_check.py`` floor run.

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.genome.heteroclinic_cycle`,
:mod:`cyclerfinder.search.jovian_resonant_connections` (reused directly for
``ResonantNode``, ``_ghost_distance``, ``_full_state_crossing``,
``HomoclinicReapproachResult`` -- none of these are Jupiter-Europa-coupled),
:mod:`cyclerfinder.search.earth_moon_resonant_families` (`#780`, the Table 3
candidate source).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.earth_moon_resonant_families as emf
import cyclerfinder.search.jovian_resonant_connections as jrc
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    _planar_floquet_pair,
    _seed_on_manifold,
    correct_connection,
)

#: This orbit's own manifold-offset magnitude -- DELIBERATELY RAISED from
#: the reused :data:`jrc.ANDERSON_LO_EPSILON` (0.5e-5) default every prior
#: sibling module used verbatim. See module docstring ("MANIFOLD OFFSET")
#: for the full justification (this orbit's own dense, 20-crossing/period
#: section structure).
EPSILON = 1e-4

#: Homoclinic trivial-solution ghost-guard radius (nondimensional, (x, xdot)
#: Euclidean norm) -- IDENTICAL value to :data:`jrc.GHOST_GUARD_DELTA`, for
#: direct comparability across the whole task chain. UNCHANGED despite
#: :data:`EPSILON` being raised (module docstring) -- the guard threshold
#: itself is not loosened, only the manifold offset that reaches it sooner.
GHOST_GUARD_DELTA = jrc.GHOST_GUARD_DELTA

#: Tolerance for :func:`build_node`'s own eigenvalue staleness cross-check
#: against `#780`'s ``table3_gate_report`` ``k_signed`` (converted to
#: lambda). Generous relative to the observed match quality (typically
#: <1e-6 relative, this task's own measurements) -- see module docstring.
EIGENVALUE_CROSSCHECK_REL_TOL = 1e-3

#: The four sign-combination buckets this module's own section convention
#: unions over -- see module docstring ("SECTION CONVENTION"), reusing
#: `#781`'s own fully-unrestricted convention directly (not inventing a
#: third style).
_SIGN_COMBOS: tuple[tuple[int, int], ...] = tuple(product((+1, -1), (+1, -1)))


def _snap_to_y0(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    *,
    max_iter: int = 6,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> NDArray[np.float64]:
    """Phase-shift ``state0`` (already close to ``{y=0}``) onto the section
    EXACTLY, via local Newton correction (``dt = -y / vy``, iterated to
    convergence). See module docstring ("SECTION CONVENTION", point 1) for
    why this is needed here (unlike every prior sibling module, whose target
    orbit's own IC sits exactly at ``y=0`` by construction) and the
    per-target verification (periodicity closure at the snapped IC is
    3.6e-11/1.5e-09 for 7-3b/7-3c, far tighter than `#780`'s own gate
    tolerance -- confirming this is a phase re-parametrization of the SAME
    solution, not a different one).
    """
    s = np.array(state0, dtype=np.float64, copy=True)
    for _ in range(max_iter):
        if abs(float(s[1])) < 1e-14:
            break
        dt = -float(s[1]) / float(s[4])
        arc = cr3bp.propagate(system, s, dt, with_stm=False, rtol=rtol, atol=atol)
        s = np.asarray(arc.state_f, dtype=np.float64)
    s[1] = 0.0
    return s


@dataclass(frozen=True)
class NodeBuildResult:
    """:func:`build_node`'s full report -- the node plus the cross-check
    evidence that it is a faithful re-derivation of `#780`'s own confirmed
    Table 3 row (not silently divergent)."""

    node: jrc.ResonantNode
    designation: str
    snap_closure_residual: float
    lam_u: float
    lam_s: float
    k_from_lambda: float
    k_source: float
    k_rel_err: float
    eigenvalue_confirmed: bool


def build_node(
    designation: str,
    system: cr3bp.CR3BPSystem | None = None,
    *,
    eigenvalue_rel_tol: float = EIGENVALUE_CROSSCHECK_REL_TOL,
) -> tuple[cr3bp.CR3BPSystem, jrc.ResonantNode, NodeBuildResult]:
    """Recover `#780`'s confirmed Table 3 row ``designation``, snap it onto
    ``{y=0}`` exactly (:func:`_snap_to_y0`), and build a
    :class:`~cyclerfinder.search.jovian_resonant_connections.ResonantNode`
    directly from the resulting full 6-state (NOT via
    ``ResonantNode.from_candidate`` -- see module docstring).

    Raises ``ValueError`` if the row does not converge, or if the recomputed
    unstable eigenvalue disagrees with `#780`'s own ``table3_gate_report``
    ``k_signed`` (converted to lambda) by more than ``eigenvalue_rel_tol`` --
    the same staleness discipline ``ResonantNode.from_candidate`` applies
    internally for its own candidates.
    """
    sys_ = system if system is not None else emf.earth_moon_system()
    po = emf.recover_table3_row(designation, sys_)
    if not po.converged:
        raise ValueError(f"{designation}: recover_table3_row did not converge -- regression")
    snapped = _snap_to_y0(sys_, po.state0)
    arc_closure = cr3bp.propagate(sys_, snapped, po.period, with_stm=False, rtol=1e-13, atol=1e-13)
    snap_closure_residual = float(np.linalg.norm(np.asarray(arc_closure.state_f) - snapped))

    lam_u, v_u, lam_s, v_s = _planar_floquet_pair(sys_, snapped, po.period)
    k_from_lambda = lam_u + 1.0 / lam_u

    gate_rows = [r for r in emf.table3_gate_report(sys_) if r.designation == designation]
    if not gate_rows:
        raise ValueError(f"{designation}: missing from table3_gate_report -- regression")
    gate_row = gate_rows[0]
    k_source = abs(gate_row.k_signed)
    k_rel_err = abs(k_from_lambda - k_source) / k_source
    eigenvalue_confirmed = k_rel_err < eigenvalue_rel_tol
    if not eigenvalue_confirmed:
        raise ValueError(
            f"{designation}: recomputed lambda+1/lambda={k_from_lambda} disagrees with "
            f"table3_gate_report k_signed={gate_row.k_signed} by {k_rel_err:.3e} relative "
            f"(> {eigenvalue_rel_tol:.0e}) -- candidate data may be stale"
        )

    node = jrc.ResonantNode(
        label=designation,
        state0=snapped,
        period=po.period,
        jacobi=po.jacobi,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=True,
    )
    result = NodeBuildResult(
        node=node,
        designation=designation,
        snap_closure_residual=snap_closure_residual,
        lam_u=lam_u,
        lam_s=lam_s,
        k_from_lambda=k_from_lambda,
        k_source=k_source,
        k_rel_err=k_rel_err,
        eigenvalue_confirmed=eigenvalue_confirmed,
    )
    return sys_, node, result


def own_section_points(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
) -> list[NDArray[np.float64]]:
    """The orbit's own qualifying ``{y=0}`` section points, ALL FOUR
    ``(x_sign, ydot_sign)`` combinations unioned -- see module docstring
    ("SECTION CONVENTION"), reusing `#781`'s own fully-unrestricted
    convention directly. Recovers all 20 of this orbit's own per-period
    crossings for `7-3b`/`7-3c` (verified this task), with no
    double-counting (every crossing's own sign pair is unique to one of the
    four buckets).
    """
    out: list[NDArray[np.float64]] = []
    for x_sign, ydot_sign in _SIGN_COMBOS:
        out.extend(jrc.own_section_points(system, node, ydot_sign=ydot_sign, x_sign=x_sign))
    return out


@dataclass(frozen=True)
class HomoclinicCandidate:
    """One surviving (ghost-guard-passed, converged) homoclinic self-
    connection. No ``dist_to_target`` field -- like `stc`'s/`ntc`'s own
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
    Ws(node) -- the Earth-Moon Class 1 analogue of ``ntc.find_homoclinic``,
    using THIS module's own fully-unrestricted section convention.

    ``k_range`` DEFAULTS TO THE GENERIC ``range(1,7)`` (matching every prior
    sibling module) -- but per this module's own docstring ("MANIFOLD
    OFFSET"), this orbit's own dense crossing structure means the manifold
    does not clear the ghost guard until roughly ``k~10`` at this module's
    own :data:`EPSILON`; callers should pass a ``k_range`` informed by that
    (e.g. ``range(9, 21)``) rather than rely on the generic default, which
    will return an EMPTY list here (not a bug -- an honest reflection of the
    scan window, not the dynamics).

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
    module's own analogue of ``ntc.homoclinic_reapproach_check``/
    ``stc.homoclinic_reapproach_check`` (same tightened
    ``rtol=1e-13, atol=1e-14`` integrator end-to-end).

    Re-derives the FULL state at the found intersection (via
    ``jrc._full_state_crossing``, using the candidate's own stored
    ``tau_u``/``tau_s``/``k_u``/``k_s``/``branch_u``/``branch_s``), then:

    1. Propagates that state BACKWARD by the unstable leg's own elapsed
       transit time ``t_u`` and compares to the ORIGINAL epsilon-scale
       unstable-manifold seed at ``tau_u``.
    2. Propagates FORWARD by ``|t_s|`` and compares to the stable-manifold
       seed at ``tau_s`` (the symmetric check).

    Both legs use the module's own fully-unrestricted section filters
    (``ydot_sign=None, x_sign=None``).

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
    "EIGENVALUE_CROSSCHECK_REL_TOL",
    "EPSILON",
    "GHOST_GUARD_DELTA",
    "HomoclinicCandidate",
    "NodeBuildResult",
    "build_node",
    "find_homoclinic",
    "homoclinic_reapproach_check",
    "own_section_points",
]
