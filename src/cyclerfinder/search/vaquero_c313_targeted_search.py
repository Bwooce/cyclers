"""`#839` targeted Wu(3:1) <-> Ws(2:1) Earth-Moon heteroclinic search AT C=3.13.

`vaquero-31-c313-em-resonant-po-2013` sits at Jacobi constant C=3.13, INSIDE
Kumar, Rawat, Rosengren & Ross (2026)'s published ``C_J in [3.00, 3.15]`` band
for a Wu(3:1) <-> Ws(2:1) Earth-Moon CR3BP heteroclinic connection family
(reproduced digit-grade at their own SEVEN printed Table-5 points by `#827`,
:mod:`cyclerfinder.search.kumar_em_resonant_heteroclinics`). C=3.13 is NOT one
of those seven printed points (2.54, 2.70, 2.86, 3.00, 3.05, 3.10, 3.15), so
`#827` could not cover it -- there is no printed digit-grade target at this C
for anyone to reproduce. This module is a FRESH SEARCH, not a reproduction:
whether a genuine connection exists at THIS row's own C, seeded from this
project's own converged nodes, verified by `#822`'s self-consistency battery
alone (no digit comparison gate, because no print exists to compare against).

WHY THIS MATTERS: `vaquero-31-c313-em-resonant-po-2013` is classed
``resonant_po`` (not ``cycler``) specifically because `#453`/`#811` found "no
demonstrated transport utility" for it (schema v4.9). A genuine transport
connection touching this row at its own C would be evidence relevant to that
``orbit_class`` question -- but this module is the COMPUTE half only; any
``orbit_class`` change is a separate adjudication task (`#839`'s own
registration, following the `#822`/`#828` and `#827`/`#854` split precedent).

MU: the catalogue row's own ``mass_ratio`` is THIS PROJECT'S registry
Earth-Moon mu (0.01215058439469525, ``cr3bp.cr3bp_system('Earth','Moon').mu``)
-- NOT Kumar's own printed mu (1.2150584270572e-2, `#827`'s
``KUMAR_MU``, differing at the ~1.24e-10 absolute level). Since we are
searching around THIS row's own converged state (not reproducing Kumar's
printed digits), every node and manifold in this module is built at the
project's own registry mu, matching the catalogue row.

THE 3:1 NODE: the catalogue row itself. Its ``state_nd``/``period_nd`` are
`#799`/`#811`'s own DERIVED reproduction of Vaquero (2013) Sec. 4.4.7's 3:1
family top C endpoint (``VAQUERO_C_RANGE_31 = (2.54, 3.13)`` --
:data:`cyclerfinder.search.vaquero_em_cyclers.VAQUERO_C_RANGE_31` --
C=3.13 IS her family's own printed top endpoint). :func:`build_node31_c313`
re-converges from the catalogue's own IC guess with the same fixed-Jacobi
symmetric corrector (``half_crossings=3``,
:data:`~cyclerfinder.search.vaquero_em_cyclers.VAQUERO_HALF_CROSSINGS`) and
requires agreement with the catalogue row's own printed digits to
:data:`NODE31_IC_ABS_TOL` -- the same reproduce-before-trust discipline
`#827`'s ``build_kumar_node`` uses against Kumar's Table 6 print, here anchored
on this project's own already-independently-verified catalogue value instead.

THE 2:1 NODE: Vaquero's OWN 2:1 family does not reach C=3.13 at all
(``VAQUERO_C_RANGE_21 = (1.98, 2.66)``) -- so there is no catalogued or
Vaquero-sourced 2:1 member here. Kumar's own Table 6 does print 2:1 family ICs
up to C=3.15 (rows at 3.10 and 3.15), a wider continuation of the SAME
resonant family (confirmed by comparing the catalogue's own 3:1 state to
Kumar's Table 6 3:1 rows: linearly interpolating Kumar's printed
``(3, 3.10)``/``(3, 3.15)`` rows to C=3.13 lands within ~2e-3 of the catalogue
row's own ``state_nd`` -- consistent with the SAME physical family under two
independent constructions and near-identical mu, not a coincidence).
:func:`build_node21_c313` therefore builds a genuinely NEW node: converges
Kumar's own Table-6 2:1 guesses at C=3.10 and C=3.15 AT THIS PROJECT'S OWN MU
(an external published bracket, re-converged, not trusted raw), then
step-continues in 0.001 Jacobi increments (:data:`NODE21_CONTINUATION_STEP` --
a coarser 0.01 step was tried first and measured to jump BRANCHES, see that
constant's own docstring) from the C=3.10 anchor to C=3.13 (a genuine
continuation, not a single interpolated jump), checking the trajectory lands
close to a linear interpolation of the two own-mu bracket endpoints (a
continuity/family-identity check, not a correctness gate on its own).

THE SEARCH: `#822`'s generic Poincare-section Newton connection machinery
(:mod:`cyclerfinder.search.vaquero_em_cycler_connections`, reused UNCHANGED --
its ``manifold_section_crossings``/``find_connection_seeds``/
``refine_connection``/``verify_connection`` are node-generic), driven by
:func:`find_c313_connection` (structurally `#822`'s own ``find_free_transfer``,
generalized off the ``OVERLAP_GRID_ICS`` table lookup to accept arbitrary
nodes). UNSEEDED: unlike `#827`, there is no printed state at this C to select
which manifold intersection to converge, so every same-ydot-sign close
approach between the two manifolds' whole crossing sets is a candidate, in
distance order (the same search-strategy `#822`'s own C=2.60 discovery run
used). Direction Wu(3:1) -> Ws(2:1), matching Kumar's own Table-5 convention;
Kumar's own finding that "Type 1 [short] intersections cease to exist for
C >= 3.09" (Section 5.2.1, digested in `#827`'s module) means C=3.13 is
expected to need a "Type 2" (long) transfer if one exists at all -- the same
qualitative regime as Kumar's own printed 3.10/3.15 rows, which `#827`
reproduced at the SAME ``n_periods=5.0`` default used here (no special-casing
needed, the manifold whole-crossing-set construction does not distinguish
Type 1 from Type 2).

EPSILON: :data:`~cyclerfinder.search.kumar_em_resonant_heteroclinics.KUMAR_EPSILON`
(1e-4), reused verbatim and for the SAME reason `#827` derived it: the 2:1
node's saddle in this C range is extreme (Kumar's own Table 6 rows put
``|lambda|`` in the hundreds near C=3.10-3.15), which amplifies the sibling
default epsilon (0.5e-5) past the forward-reapproach gate; 1e-4 was measured
(not guessed) to bring the whole `#822` battery back under its unchanged
ceilings at this SAME family branch. Chosen for this reason before running the
search here, not loosened after seeing a result.

VERIFICATION: the full, UNCHANGED `#822` battery
(:func:`~cyclerfinder.search.vaquero_em_cycler_connections.verify_connection`)
-- full planar-4-state gap, ydot-sign hard gate, ghost guard against both
nodes' own section points, independent-Radau re-derivation, forward/backward
manifold re-approach, Jacobi drift on both legs. No digit-comparison gate (no
print exists at this C); ``passed`` alone is the acceptance criterion. An
honest negative (no seed both converges and verifies) is reported as such, not
loosened into a pass.

LITERATURE-NOVELTY: NOT this module's own scope -- `#822`'s own live-WebSearch
mandatory-floor run against this exact paper's connection concept already
returned "published" (confidence 0.95); a genuine hit at C=3.13 is a NEW POINT
in an already-published family, and any ``orbit_class``-relevant claim is
deferred to `#839`'s registered adjudication follow-up, not decided here.

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic` (fixed-Jacobi corrector),
:mod:`cyclerfinder.search.vaquero_em_cyclers` (``VAQUERO_HALF_CROSSINGS``,
``half_crossing_index``), :mod:`cyclerfinder.search.vaquero_em_cycler_connections`
(`#822` machinery, reused unchanged),
:mod:`cyclerfinder.search.kumar_em_resonant_heteroclinics` (``KUMAR_EPSILON``),
:mod:`cyclerfinder.search.jovian_resonant_connections` (``ResonantNode``),
:mod:`cyclerfinder.genome.heteroclinic_cycle` (Floquet pair).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.vaquero_em_cycler_connections as vcc
import cyclerfinder.search.vaquero_em_cyclers as vem
from cyclerfinder.genome.heteroclinic_cycle import HeteroclinicConnection, _planar_floquet_pair
from cyclerfinder.search.kumar_em_resonant_heteroclinics import KUMAR_EPSILON

#: Target Jacobi constant -- the catalogue row's own recorded value
#: (``vaquero-31-c313-em-resonant-po-2013.orbit_elements.cr3bp.jacobi_constant``,
#: 1e-14-scale grid-float offset from 3.13, not a recovered-vs-target gap).
C313 = 3.129999999999993

#: The catalogue row's own DERIVED state/period (`#799`/`#811`, already
#: independently Radau-cross-checked in the catalogue's own writeback) --
#: PRINT provenance here is ``data/catalogue.yaml``'s
#: ``vaquero-31-c313-em-resonant-po-2013`` row, transcribed verbatim.
CATALOGUE_STATE_ND_31: tuple[float, float, float, float, float, float] = (
    0.8135643069819515,
    0.0,
    0.0,
    0.0,
    -0.25304820538336525,
    0.0,
)
CATALOGUE_PERIOD_ND_31 = 6.45496522207971

#: Node-31 reproduce-before-trust ceiling: agreement of the re-converged
#: ``(x0, ydot0)`` with the catalogue row's own recorded state, absolute.
#: Looser than `#827`'s ``TABLE6_IC_ABS_TOL`` (1e-9) because the catalogue's
#: own recorded ``jacobi_constant`` (3.129999999999993) carries a ~7e-15
#: grid-float offset from the literal target 3.13 that this module targets
#: instead (see :data:`C313`'s docstring) -- 1e-8 comfortably separates a
#: genuine re-derivation (expected ~1e-9..1e-11, the corrector's own
#: convergence floor) from a mis-transcription of the row's printed digits.
NODE31_IC_ABS_TOL = 1e-8

#: Kumar et al. (2026) Table 6 (Appendix 8.2) 2:1-family rows bracketing
#: C=3.13 -- PRINT provenance, transcribed verbatim (same source `#827`'s
#: ``KUMAR_TABLE6_ICS`` vendors; repeated here rather than imported so this
#: module's own external-bracket role at THIS project's mu is explicit).
#: ``C -> (x0, ydot0)`` at their own mu (1.2150584270572e-2); used ONLY as a
#: continuation-seed guess, re-converged at THIS project's own registry mu.
KUMAR_TABLE6_21_BRACKET: dict[float, tuple[float, float]] = {
    3.10: (0.878280334961, -0.334629543419),
    3.15: (0.864401165205, -0.219058827351),
}

#: Continuation step size (Jacobi) walking the 2:1 node from the C=3.10
#: bracket anchor to C=3.13. MEASURED (not guessed): a 0.01 step lands the
#: fixed-``half_crossings``-index corrector on a DIFFERENT branch (observed
#: directly: T jumps 6.8 -> 11.8 -> 6.1 -> 3.5 nd across three 0.005-scale
#: steps -- a topology confusion, not a converged continuation) because
#: ``dx0/dC`` is O(1) here and a 0.01-wide step overshoots the corrector's
#: basin of the SAME branch. 0.001 walks smoothly and monotonically
#: (x0/ydot0/period all vary smoothly step to step, verified directly over
#: the full 3.10->3.13 walk) -- this module's own actually-used step.
NODE21_CONTINUATION_STEP = 0.001

#: Family-identity continuity check ceiling: the C=3.13 2:1 node's own
#: ``(x0, ydot0)`` must land within this absolute distance of a LINEAR
#: interpolation between the two own-mu-reconverged C=3.10/C=3.15 bracket
#: endpoints -- a smooth, non-bifurcating family produces an interpolation
#: error of O(step^2) in C (~1e-3 scale over a 0.05-wide bracket for a family
#: this mildly curved, measured directly at build time); 0.02 is generous
#: above that floor and far below the ~0.1+ scale a genuinely different
#: branch would show.
NODE21_INTERP_CHECK_TOL = 0.02

#: Manifold-offset magnitude -- reused verbatim from `#827`'s own module (see
#: this module's docstring, "EPSILON").
EPSILON = KUMAR_EPSILON


def em_system() -> cr3bp.CR3BPSystem:
    """This project's canonical Earth-Moon CR3BP system (registry mu/l*/t*) --
    the SAME mu as the catalogue row's own ``mass_ratio``."""
    return cr3bp.cr3bp_system("Earth", "Moon")


def build_node31_c313(
    system: cr3bp.CR3BPSystem,
    *,
    tol: float = 1e-12,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    ic_abs_tol: float = NODE31_IC_ABS_TOL,
) -> jrc.ResonantNode:
    """Re-derive ``vaquero-31-c313-em-resonant-po-2013``'s own orbit as a
    connection node.

    Re-converges at exactly :data:`C313` from the catalogue row's own IC guess
    (the same fixed-Jacobi symmetric corrector `#799`/`#811` used,
    ``half_crossings=3``), then REQUIRES the converged ``(x0, ydot0)`` to agree
    with the catalogue row's own recorded ``state_nd`` to
    :data:`NODE31_IC_ABS_TOL` -- reproduce-before-trust, anchored on the
    catalogue's own already-independently-verified value.
    """
    x0_guess = CATALOGUE_STATE_ND_31[0]
    ydot0_target = CATALOGUE_STATE_ND_31[4]
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_guess,
        C313,
        CATALOGUE_PERIOD_ND_31,
        ydot0_sign=-1.0,
        half_crossings=vem.VAQUERO_HALF_CROSSINGS,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"node31 corrector failed at C={C313}: residual {orbit.crossing_residual:.3e}"
        )
    dx0 = abs(orbit.x0 - x0_guess)
    dydot0 = abs(orbit.ydot0 - ydot0_target)
    if dx0 > ic_abs_tol or dydot0 > ic_abs_tol:
        raise ValueError(
            f"node31 re-converged IC disagrees with the catalogue row's own state_nd: "
            f"|dx0|={dx0:.3e}, |dydot0|={dydot0:.3e} (> {ic_abs_tol:.0e})"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    lam_u, v_u, _lam_s, v_s = _planar_floquet_pair(system, state0, orbit.period)
    if lam_u <= 1.05:
        raise ValueError(f"node31 at C={C313}: |lambda|={lam_u:.4f} is not materially unstable")
    return jrc.ResonantNode(
        label="vaquero-31-c313",
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=True,
    )


@dataclass(frozen=True)
class Node21Provenance:
    """Diagnostic record of :func:`build_node21_c313`'s continuation build."""

    bracket_310: tuple[float, float]  # own-mu re-converged (x0, ydot0) at C=3.10
    bracket_315: tuple[float, float]  # own-mu re-converged (x0, ydot0) at C=3.15
    bracket_310_lambda: float
    bracket_315_lambda: float
    trace: list[tuple[float, float, float]]  # (C, x0, ydot0) at every continuation step
    interp_x0: float  # linear interpolation of the two bracket endpoints, at C313
    interp_ydot0: float
    interp_distance: float  # |own_final - interpolated| in (x0, ydot0)


def build_node21_c313(
    system: cr3bp.CR3BPSystem,
    *,
    step: float = NODE21_CONTINUATION_STEP,
    tol: float = 1e-12,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    interp_check_tol: float = NODE21_INTERP_CHECK_TOL,
) -> tuple[jrc.ResonantNode, Node21Provenance]:
    """Build the 2:1 family's C=3.13 member at THIS project's own mu.

    Neither Vaquero's own 2:1 family (``VAQUERO_C_RANGE_21`` tops out at 2.66)
    nor this project's catalogue covers a 2:1 member anywhere near C=3.13.
    Kumar et al.'s own Table 6 prints 2:1 rows at C=3.10 and C=3.15 (their own
    mu) -- both re-converged here at THIS project's registry mu as an external
    bracket, then step-continued (:data:`NODE21_CONTINUATION_STEP`) from the
    C=3.10 anchor to C=3.13. The final member's ``(x0, ydot0)`` is checked
    against a linear interpolation of the two own-mu bracket endpoints
    (:data:`NODE21_INTERP_CHECK_TOL`) -- a continuity/family-identity sanity
    check, not a substitute for the continuation itself.
    """
    o310 = cp.correct_symmetric_fixed_jacobi(
        system,
        KUMAR_TABLE6_21_BRACKET[3.10][0],
        3.10,
        6.4,
        ydot0_sign=-1.0,
        half_crossings=None,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not o310.converged:
        raise RuntimeError(
            f"node21 bracket corrector failed at C=3.10: residual {o310.crossing_residual:.3e}"
        )
    o315 = cp.correct_symmetric_fixed_jacobi(
        system,
        KUMAR_TABLE6_21_BRACKET[3.15][0],
        3.15,
        6.4,
        ydot0_sign=-1.0,
        half_crossings=None,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not o315.converged:
        raise RuntimeError(
            f"node21 bracket corrector failed at C=3.15: residual {o315.crossing_residual:.3e}"
        )

    half_idx = vem.half_crossing_index(system, o310, rtol=rtol, atol=atol)
    half_idx_315 = vem.half_crossing_index(system, o315, rtol=rtol, atol=atol)
    if half_idx != half_idx_315:
        raise RuntimeError(
            f"node21 bracket topology mismatch: half-crossing index {half_idx} (C=3.10) "
            f"vs {half_idx_315} (C=3.15) -- not the same continuous family branch"
        )

    c = 3.10
    x0 = o310.x0
    period_guess = o310.period
    trace: list[tuple[float, float, float]] = [(c, o310.x0, o310.ydot0)]
    final = o310
    while c < C313 - 1e-9:
        c_next = min(round(c + step, 10), C313)
        orbit = cp.correct_symmetric_fixed_jacobi(
            system,
            x0,
            c_next,
            period_guess,
            ydot0_sign=-1.0,
            half_crossings=half_idx,
            tol=tol,
            rtol=rtol,
            atol=atol,
        )
        if not orbit.converged:
            raise RuntimeError(
                f"node21 continuation step to C={c_next} failed: "
                f"residual {orbit.crossing_residual:.3e}"
            )
        x0 = orbit.x0
        period_guess = orbit.period
        c = c_next
        final = orbit
        trace.append((c, orbit.x0, orbit.ydot0))

    frac = (C313 - 3.10) / (3.15 - 3.10)
    interp_x0 = o310.x0 + frac * (o315.x0 - o310.x0)
    interp_ydot0 = o310.ydot0 + frac * (o315.ydot0 - o310.ydot0)
    interp_distance = math.hypot(final.x0 - interp_x0, final.ydot0 - interp_ydot0)
    if interp_distance > interp_check_tol:
        raise ValueError(
            f"node21 at C={C313}: continuation result (x0={final.x0}, ydot0={final.ydot0}) "
            f"disagrees with the linear own-mu bracket interpolation "
            f"(x0={interp_x0}, ydot0={interp_ydot0}) by {interp_distance:.3e} "
            f"(> {interp_check_tol:.0e}) -- possible branch jump"
        )

    state0 = np.array([final.x0, 0.0, 0.0, 0.0, final.ydot0, 0.0], dtype=np.float64)
    lam_u, v_u, _lam_s, v_s = _planar_floquet_pair(system, state0, final.period)
    if lam_u <= 1.05:
        raise ValueError(f"node21 at C={C313}: |lambda|={lam_u:.4f} is not materially unstable")

    lam_310, _, _, _ = _planar_floquet_pair(
        system, np.array([o310.x0, 0.0, 0.0, 0.0, o310.ydot0, 0.0]), o310.period
    )
    lam_315, _, _, _ = _planar_floquet_pair(
        system, np.array([o315.x0, 0.0, 0.0, 0.0, o315.ydot0, 0.0]), o315.period
    )

    node = jrc.ResonantNode(
        label="kumar-21-c313-own-mu",
        state0=state0,
        period=final.period,
        jacobi=final.jacobi,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=True,
    )
    prov = Node21Provenance(
        bracket_310=(o310.x0, o310.ydot0),
        bracket_315=(o315.x0, o315.ydot0),
        bracket_310_lambda=lam_310,
        bracket_315_lambda=lam_315,
        trace=trace,
        interp_x0=interp_x0,
        interp_ydot0=interp_ydot0,
        interp_distance=interp_distance,
    )
    return node, prov


@dataclass(frozen=True)
class C313SearchResult:
    """End-to-end result of :func:`find_c313_connection`."""

    jacobi: float
    connection: HeteroclinicConnection | None
    evidence: vcc.ConnectionEvidence | None
    n_seeds: int
    n_refined: int
    n_converged: int
    notes: str = ""


def find_c313_connection(
    system: cr3bp.CR3BPSystem,
    node31: jrc.ResonantNode,
    node21: jrc.ResonantNode,
    *,
    n_tau: int = 48,
    n_periods: float = 5.0,
    max_seed_distance: float = 0.02,
    max_refine: int = 12,
    seed_diversity_radius: float = 0.05,
    epsilon: float = EPSILON,
    max_time_factor: float = 10.0,
    tol: float = 1e-9,
) -> C313SearchResult:
    """UNSEEDED Wu(3:1) -> Ws(2:1) connection search at C=3.13.

    Structurally `#822`'s own ``find_free_transfer`` (this module's docstring,
    "THE SEARCH"), generalized off the ``OVERLAP_GRID_ICS`` table to accept
    arbitrary already-built nodes. Every same-ydot-sign close approach between
    the two manifolds' whole ``{y=0}`` crossing sets is a Newton candidate, in
    distance order; the FIRST that both converges AND passes `#822`'s full
    :func:`~cyclerfinder.search.vaquero_em_cycler_connections.verify_connection`
    battery is returned. An honest negative (``connection=None``) otherwise --
    conditional on these settings, never a fabricated pass.
    """
    crossings_u: dict[int, list[vcc.ManifoldCrossing]] = {}
    crossings_s: dict[int, list[vcc.ManifoldCrossing]] = {}
    for branch in (+1, -1):
        crossings_u[branch] = vcc.manifold_section_crossings(
            system,
            node31,
            direction="unstable",
            branch=branch,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )
        crossings_s[branch] = vcc.manifold_section_crossings(
            system,
            node21,
            direction="stable",
            branch=branch,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )

    seeds: list[vcc.ConnectionSeed] = []
    for bu in (+1, -1):
        for bs in (+1, -1):
            seeds.extend(
                vcc.find_connection_seeds(
                    crossings_u[bu],
                    crossings_s[bs],
                    branch_u=bu,
                    branch_s=bs,
                    max_distance=max_seed_distance,
                )
            )
    seeds.sort(key=lambda s: s.distance)

    n_refined = 0
    n_converged = 0
    tried_xv: list[tuple[float, float]] = []
    for seed in seeds:
        if n_refined >= max_refine:
            break
        if any(
            math.hypot(seed.x - x, seed.xdot - xd) < seed_diversity_radius for x, xd in tried_xv
        ):
            continue
        tried_xv.append((seed.x, seed.xdot))
        n_refined += 1
        conn = vcc.refine_connection(
            system,
            node31,
            node21,
            seed,
            epsilon=epsilon,
            max_time_factor=max_time_factor,
            tol=tol,
        )
        if not conn.converged:
            continue
        n_converged += 1
        ev = vcc.verify_connection(
            system, node31, node21, conn, epsilon=epsilon, max_time_factor=max_time_factor
        )
        if ev.passed:
            return C313SearchResult(
                jacobi=float(node31.jacobi),
                connection=conn,
                evidence=ev,
                n_seeds=len(seeds),
                n_refined=n_refined,
                n_converged=n_converged,
            )
    return C313SearchResult(
        jacobi=float(node31.jacobi),
        connection=None,
        evidence=None,
        n_seeds=len(seeds),
        n_refined=n_refined,
        n_converged=n_converged,
        notes=(
            "no verified connection among refined seeds -- honest negative at these "
            "C=3.13 settings (conditional on n_tau/n_periods/max_refine, not proof of absence)"
        ),
    )


__all__ = [
    "C313",
    "CATALOGUE_PERIOD_ND_31",
    "CATALOGUE_STATE_ND_31",
    "EPSILON",
    "KUMAR_TABLE6_21_BRACKET",
    "NODE21_CONTINUATION_STEP",
    "NODE21_INTERP_CHECK_TOL",
    "NODE31_IC_ABS_TOL",
    "C313SearchResult",
    "Node21Provenance",
    "build_node21_c313",
    "build_node31_c313",
    "em_system",
    "find_c313_connection",
]
