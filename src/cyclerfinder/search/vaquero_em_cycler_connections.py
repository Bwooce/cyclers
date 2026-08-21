"""`#822` Vaquero 2:1 <-> 3:1 Earth-Moon cycler-family free transfer (heteroclinic).

Computes an actual FREE TRANSFER -- an unstable-to-unstable, same-Jacobi-constant
heteroclinic invariant-manifold connection -- between Vaquero (2013)'s own 2:1 and
3:1 Earth-Moon "Periodic Cycler" families (Purdue Ph.D. dissertation, Sec. 4.4.7,
pp. 169-172), in the two families' printed Jacobi-constant OVERLAP band
``C in [2.54, 2.66]``.

WHAT THE SOURCE ACTUALLY CLAIMS (and does not): Vaquero asserts the connection's
EXISTENCE in prose (pp. 171-172; including the "2.66 < C < 2.54" typesetting-slip
reading `#787`'s digest documents) -- both families are unstable across the overlap
band, which is her own stated precondition for a free (manifold-to-manifold)
transfer -- but she prints NO transfer trajectory, NO intersection state, and NO
delta-v anywhere in the section. There is therefore no digit-grade published target
to gate against (unlike Anderson & Lo 2011's Tables 2/3): success here is a
GENUINE, independently-cross-checked connection, evidenced by this module's own
self-consistency battery (Newton residual, full-4-state gap including the ydot-sign
check below, ghost-guard margin, independent-Radau re-derivation, forward/backward
re-approach), exactly the honest framing `#766`/`#767`/`#781` established for
targets without published digits. This is a REPRODUCTION of a published qualitative
claim -- nothing here is novel and nothing is claimed novel (the connection's
existence is Vaquero's own assertion; this module supplies the missing numbers).

RELATION TO `#783` (:mod:`cyclerfinder.search.earth_moon_resonant_connections`):
that module's own Earth-Moon connection attempt was a CLEAN NEGATIVE, but on a
categorically different problem: a Barrabes-Mondelo-Olle-style HOMOCLINIC
reproduction of Casoliva's own single hardest He1 energy point (|lambda|~=108.6,
legs of ~4.1 periods each, a cold start at her own 5h20m-CPU point -- see that
module's "CONDITIONING WALL"). This module is NOT that algorithm class: it uses
this project's own established Poincare-section Newton machinery
(:func:`~cyclerfinder.genome.heteroclinic_cycle.correct_connection`, the
`#754`/`#759`/`#767`/`#781` lane), which solves a 2-D phase problem
(``tau_u, tau_s``) on the {y=0} section rather than shooting the full state chain
-- and the target instabilities here (|lambda| 3.1-5.7 / 11.3-12.8) are far milder
than He1's. Method fit was treated as genuinely open (per the `#822` dispatch), and
in the event the machinery converges cleanly.

NODES: both family members at each shared 0.01-grid Jacobi constant are re-derived
at call time by the same fixed-Jacobi symmetric corrector that built them
(:func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`,
``half_crossings=3`` topology asserted via
:func:`~cyclerfinder.search.vaquero_em_cyclers.half_crossing_index`), seeded from
`#799`'s own archived 129-member record
(``data/found/799_vaquero_em_cycler_families/results.json``; the in-band subset is
vendored below as DERIVE-provenance guesses, NOT trusted -- each node build
re-converges to tol=1e-12 and re-checks the Floquet magnitude against the archived
``abs_lambda`` to 1e-6 relative, the same staleness discipline as
``ResonantNode.from_candidate``).

NEGATIVE-REAL SADDLE (2:1 family): across the whole overlap band the 2:1 family's
planar saddle eigenvalue is genuinely NEGATIVE real (Barden nu in [-2.95, -1.73],
lambda in [-5.73, -3.15] -- re-derived at node-build time; the full-period map
flips the unstable eigenvector's sign each period), the same situation `#781`'s
4:5-saddle documented: ``_planar_floquet_pair``'s magnitude-only eigenvector
convention is unaffected, but ``branch_u=+1/-1`` are related by the map's own
period-to-period flip rather than being fixed "sides" of the manifold. The 3:1
family's saddle is positive real (+11.3..+12.8) throughout. No special-case code
is needed (both branches are scanned); recorded so the branch labels are read
correctly.

SECTION CONVENTION + THE ydot-SIGN GATE (this module's own addition to the sibling
pattern): like `#781`, both ``x_sign`` and ``ydot_sign`` are fully UNRESTRICTED
(the two orbits' own perpendicular {y=0} crossings sit at mixed signs, and the
connection's location was unknown a priori), and the ghost guard unions all four
sign buckets of BOTH orbits' own section points. But an (x, xdot) match on {y=0}
at equal Jacobi determines only |ydot|, NOT its sign -- the 2-D section residual
alone can in principle "converge" onto a spurious pair of crossings with OPPOSITE
ydot. This module therefore (a) pre-filters candidate seeds to same-ydot-sign
close approaches (:func:`find_connection_seeds`) and (b) gates every verified
connection on the FULL planar 4-state gap at the matched crossing, which contains
the ydot difference explicitly (:func:`verify_connection` -- ``ydot_signs_match``
is a hard gate, never assumed).

MANIFOLD OFFSET / GHOST GUARD: ``EPSILON``/``GHOST_GUARD_DELTA`` reused verbatim
from :mod:`cyclerfinder.search.jovian_resonant_connections`
(``ANDERSON_LO_EPSILON`` = 0.5e-5, ``GHOST_GUARD_DELTA`` = 1e-3), the same
generic-default reasoning as `#767`/`#781`: Vaquero publishes no manifold-offset
magnitude of her own, and both problems are O(1) nondim-scale CR3BP systems.

SEARCH STRATEGY (promoted from `#781`'s scratchpad method into first-class module
code): a brute-force ``(branch, k_u, k_s)`` product scan is intractable here (the
manifolds need ~7-8 orbit periods to traverse the inter-family gap at this
epsilon, i.e. crossing indices k~30-50; a 50x50 k-grid of grid-scanned
``correct_connection`` calls would take days). Instead
:func:`manifold_section_crossings` builds each manifold's WHOLE crossing set (all
k) in one propagation per phase-grid point -- O(n_tau) integrations per (node,
branch) -- :func:`find_connection_seeds` pairs the two sets' same-ydot-sign close
approaches in (x, xdot), and :func:`refine_connection` feeds each seed's exact
``(tau_u, tau_s, k_u, k_s, branch_u, branch_s)`` to ``correct_connection`` with
explicit ``tau0`` values (skipping its internal grid scan entirely).

LITERATURE-NOVELTY GATE: run (mandatory floor, per the `#822` dispatch -- never
skipped reflexively even for an explicit reproduction); expected and actual
verdict is "published" anchored on Vaquero's own dissertation lineage -- which is
exactly consistent with this module's own framing (the claim is HERS; this module
reproduces it). See the results note
(``docs/notes/2026-08-11-822-vaquero-em-free-transfer.md``) for the query trail.

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.genome.heteroclinic_cycle` (``correct_connection`` and manifold
primitives, reused unchanged), :mod:`cyclerfinder.search.jovian_resonant_connections`
(``ResonantNode``/ghost machinery), :mod:`cyclerfinder.search.cr3bp_periodic`
(fixed-Jacobi corrector), :mod:`cyclerfinder.search.vaquero_em_cyclers` (`#799`,
the family-reproduction source of the vendored node guesses).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.neptune_triton_resonant_connections as ntc
import cyclerfinder.search.vaquero_em_cyclers as vem
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    _planar_floquet_pair,
    _section_crossing,
    _seed_on_manifold,
    correct_connection,
)

#: Jacobi-constant overlap band of the two families -- DERIVED as the
#: intersection of the two SOURCED printed ranges
#: (:data:`~cyclerfinder.search.vaquero_em_cyclers.VAQUERO_C_RANGE_21` and
#: :data:`~cyclerfinder.search.vaquero_em_cyclers.VAQUERO_C_RANGE_31`), and
#: itself the band Vaquero's own free-transfer prose names (pp. 171-172, via the
#: `#787` typesetting-slip reading of "2.66 < C < 2.54").
OVERLAP_C_RANGE: tuple[float, float] = (
    max(vem.VAQUERO_C_RANGE_21[0], vem.VAQUERO_C_RANGE_31[0]),
    min(vem.VAQUERO_C_RANGE_21[1], vem.VAQUERO_C_RANGE_31[1]),
)

#: Manifold-offset magnitude -- reused verbatim from the sibling connection
#: modules (module docstring, "MANIFOLD OFFSET / GHOST GUARD").
EPSILON = jrc.ANDERSON_LO_EPSILON

#: Trivial-solution ghost-guard radius -- identical to the sibling modules'.
GHOST_GUARD_DELTA = jrc.GHOST_GUARD_DELTA

#: `#799`-archived overlap-band members, ``(p, C) -> (x0, ydot0, period,
#: abs_lambda)`` -- DERIVE provenance (this project's own `#799` continuation,
#: ``data/found/799_vaquero_em_cycler_families/results.json``, NOT a Vaquero
#: print: no digit-grade IC exists anywhere in the source, see `#799`'s note).
#: Used ONLY as corrector guesses + staleness references by
#: :func:`build_vaquero_overlap_node` -- every node is re-converged at call time.
OVERLAP_GRID_ICS: dict[tuple[int, float], tuple[float, float, float, float]] = {
    # 2:1 family (p=2)
    (2, 2.54): (1.0905363960533268, -0.8231863180949408, 5.941227735609639, 3.149761131186353),
    (2, 2.55): (1.0868460759391638, -0.8212421258340358, 5.927264995240593, 3.4330710975402035),
    (2, 2.56): (1.0830550210674108, -0.8198841221290316, 5.912382475093488, 3.7256995552301233),
    (2, 2.57): (1.0791538388990847, -0.8192250242938581, 5.896480445261298, 4.027080902654688),
    (2, 2.58): (1.0751307804963242, -0.8194079602350062, 5.8794394046326985, 4.33539742190552),
    (2, 2.59): (1.0709707652405038, -0.8206186904669224, 5.8611128174423905, 4.6469640252474305),
    (2, 2.60): (1.0666538160256032, -0.8231049641748824, 5.8413155891331, 4.955136449486326),
    (2, 2.61): (1.0621524066177133, -0.8272089050995878, 5.819804706776062, 5.248338690717803),
    (2, 2.62): (1.057426638201889, -0.8334251509810034, 5.796244297747047, 5.506301194875569),
    (2, 2.63): (1.0524146115870048, -0.8425154854724315, 5.7701363344029275, 5.69224625551723),
    (2, 2.64): (1.0470105122376612, -0.855766761773147, 5.740663856668667, 5.7344291788645565),
    (2, 2.65): (1.0410035136327993, -0.8757021946207459, 5.706256620814711, 5.472604026372235),
    (2, 2.66): (1.0338302047346954, -0.9089334377051435, 5.662843584779122, 4.428157223146986),
    # 3:1 family (p=3)
    (3, 2.54): (0.9013301668020125, -0.8462249954358775, 6.269604424886022, 11.261165367653621),
    (3, 2.55): (0.9006624567792414, -0.839240465648521, 6.270721836375987, 11.37942217203676),
    (3, 2.56): (0.8999745154574539, -0.8321876135983043, 6.271849664970181, 11.49908487813),
    (3, 2.57): (0.8992661458401348, -0.825065915797768, 6.272988850083939, 11.620208160949336),
    (3, 2.58): (0.8985371479741101, -0.8178748231328931, 6.274140414722454, 11.74284439015189),
    (3, 2.59): (0.8977873189089175, -0.8106137578440186, 6.275305469377563, 11.867043382039416),
    (3, 2.60): (0.8970164526223031, -0.803282110321528, 6.276485216101558, 11.992852378826967),
    (3, 2.61): (0.8962243399082888, -0.7958792357057801, 6.277680952792588, 12.120314881178183),
    (3, 2.62): (0.8954107682223105, -0.788404450277243, 6.278894077681616, 12.249470446490303),
    (3, 2.63): (0.8945755214795653, -0.7808570276244745, 6.280126094062637, 12.380352792877849),
    (3, 2.64): (0.8937183798006902, -0.7732361945737354, 6.281378615266163, 12.512989784358279),
    (3, 2.65): (0.892839119199225, -0.7655411268634741, 6.282653369903586, 12.647402406617443),
    (3, 2.66): (0.8919375112041409, -0.7577709445440712, 6.283952207405823, 12.783603246229166),
}


def build_vaquero_overlap_node(
    system: cr3bp.CR3BPSystem,
    p: int,
    c: float,
    *,
    tol: float = 1e-12,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    eigenvalue_rel_tol: float = 1e-6,
) -> jrc.ResonantNode:
    """Re-derive one overlap-band family member and wrap it as a connection node.

    Re-converges the member at exactly Jacobi ``c`` with the same fixed-Jacobi
    corrector `#799` used (guess from :data:`OVERLAP_GRID_ICS`, never trusted
    raw), asserts the family's own ``half_crossings=3`` crossing topology, then
    recomputes the planar Floquet saddle pair and checks its magnitude against
    the `#799`-archived ``abs_lambda`` to ``eigenvalue_rel_tol`` relative (the
    same staleness discipline as ``ResonantNode.from_candidate``; the archived
    value's SIGN is not compared -- the 2:1 family's saddle is genuinely
    negative real, module docstring, and ``_planar_floquet_pair`` returns
    magnitudes).

    Raises ``KeyError`` for a ``(p, c)`` outside the vendored overlap grid and
    ``ValueError``/``RuntimeError`` on any convergence/topology/staleness
    failure -- loud, never a silently-degraded node.
    """
    key = (p, round(float(c), 4))
    if key not in OVERLAP_GRID_ICS:
        raise KeyError(
            f"({p}, {c}) is not an overlap-band grid member; valid C grid: "
            f"[{OVERLAP_C_RANGE[0]}, {OVERLAP_C_RANGE[1]}] step 0.01"
        )
    x0_g, _ydot0_g, period_g, abs_lambda_ref = OVERLAP_GRID_ICS[key]
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_g,
        float(c),
        period_g,
        ydot0_sign=-1.0,
        half_crossings=vem.VAQUERO_HALF_CROSSINGS,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"node corrector failed for p={p} C={c}: residual {orbit.crossing_residual:.3e}"
        )
    idx = vem.half_crossing_index(system, orbit, rtol=rtol, atol=atol)
    if idx != vem.VAQUERO_HALF_CROSSINGS:
        raise RuntimeError(
            f"unexpected half-crossing topology for p={p} C={c}: measured {idx}, "
            f"expected {vem.VAQUERO_HALF_CROSSINGS}"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    lam_u, v_u, _lam_s, v_s = _planar_floquet_pair(system, state0, orbit.period)
    rel_err = abs(lam_u - abs_lambda_ref) / abs_lambda_ref
    if rel_err > eigenvalue_rel_tol:
        raise ValueError(
            f"p={p} C={c}: recomputed |lambda|={lam_u} disagrees with #799's archived "
            f"{abs_lambda_ref} by {rel_err:.3e} relative (> {eigenvalue_rel_tol:.0e})"
        )
    return jrc.ResonantNode(
        label=f"vaquero-{p}1-c{round(float(c) * 100):d}",
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=True,
    )


@dataclass(frozen=True)
class ManifoldCrossing:
    """One {y=0} crossing of one manifold trajectory (all signs kept)."""

    tau: float
    k: int  # 1-based crossing index along the trajectory (all signs counted)
    t: float  # elapsed manifold-transit time to this crossing (signed: <0 stable)
    x: float
    xdot: float
    ydot: float


def manifold_section_crossings(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    *,
    direction: str,
    branch: int,
    n_tau: int,
    n_periods: float,
    epsilon: float = EPSILON,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[ManifoldCrossing]:
    """All {y=0} crossings of ``n_tau`` manifold trajectories of ``node``.

    One propagation per phase-grid point (``tau = T*(i+0.5)/n_tau``), collecting
    EVERY crossing (both ``x``/``ydot`` signs, k unbounded) out to
    ``n_periods * node.period`` -- the O(n_tau) whole-crossing-set construction
    the module docstring's "SEARCH STRATEGY" describes. ``t`` is signed (negative
    for the stable/backward direction), matching the propagation clock.
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    if branch not in (+1, -1):
        raise ValueError(f"branch must be +1 or -1; got {branch!r}")
    max_time = float(n_periods) * node.period
    t_span = (0.0, max_time) if direction == "unstable" else (0.0, -max_time)

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]

    out: list[ManifoldCrossing] = []
    t_floor = 1e-6 * max_time
    for i in range(n_tau):
        tau = node.period * (i + 0.5) / n_tau
        seed = _seed_on_manifold(
            system, node, tau=tau, direction=direction, branch=branch, epsilon=epsilon
        )
        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            t_span,
            seed,
            args=(system.mu,),  # type: ignore[call-overload]
            method="DOP853",
            rtol=rtol,
            atol=atol,
            events=_y_event,
            max_step=max_time / 1000.0,
        )
        t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
        y_events = sol.y_events[0] if sol.y_events is not None else []
        k = 0
        for t_ev, y_ev in zip(t_events, y_events, strict=False):
            if abs(float(t_ev)) <= t_floor:
                continue  # the t~0 root at the seed itself
            k += 1
            out.append(
                ManifoldCrossing(
                    tau=float(tau),
                    k=k,
                    t=float(t_ev),
                    x=float(y_ev[0]),
                    xdot=float(y_ev[3]),
                    ydot=float(y_ev[4]),
                )
            )
    return out


@dataclass(frozen=True)
class ConnectionSeed:
    """One same-ydot-sign close approach between the two manifolds' crossing
    sets -- a Newton starting point for :func:`refine_connection`."""

    distance: float  # (x, xdot) Euclidean gap at the close approach
    branch_u: int
    branch_s: int
    k_u: int
    k_s: int
    tau_u: float
    tau_s: float
    x: float
    xdot: float
    ydot: float


def find_connection_seeds(
    crossings_u: list[ManifoldCrossing],
    crossings_s: list[ManifoldCrossing],
    *,
    branch_u: int,
    branch_s: int,
    max_distance: float = 0.02,
) -> list[ConnectionSeed]:
    """Pair the unstable and stable crossing sets' close approaches in
    ``(x, xdot)``, keeping ONLY same-ydot-sign pairs (module docstring, "THE
    ydot-SIGN GATE" -- an opposite-sign pair is not a connection candidate at
    all, however close its (x, xdot) gap). For each unstable crossing the single
    nearest stable crossing is taken; pairs with gap > ``max_distance`` are
    dropped. Sorted by gap, closest first.
    """
    if not crossings_u or not crossings_s:
        return []
    arr_s = np.array([(c.x, c.xdot) for c in crossings_s])
    sign_s = np.sign([c.ydot for c in crossings_s])
    out: list[ConnectionSeed] = []
    for cu in crossings_u:
        dd = np.hypot(arr_s[:, 0] - cu.x, arr_s[:, 1] - cu.xdot)
        dd = np.where(sign_s == np.sign(cu.ydot), dd, np.inf)
        j = int(np.argmin(dd))
        if not np.isfinite(dd[j]) or float(dd[j]) > max_distance:
            continue
        cs = crossings_s[j]
        out.append(
            ConnectionSeed(
                distance=float(dd[j]),
                branch_u=branch_u,
                branch_s=branch_s,
                k_u=cu.k,
                k_s=cs.k,
                tau_u=cu.tau,
                tau_s=cs.tau,
                x=cu.x,
                xdot=cu.xdot,
                ydot=cu.ydot,
            )
        )
    out.sort(key=lambda s: s.distance)
    return out


def refine_connection(
    system: cr3bp.CR3BPSystem,
    node_from: jrc.ResonantNode,
    node_to: jrc.ResonantNode,
    seed: ConnectionSeed,
    *,
    epsilon: float = EPSILON,
    max_time_factor: float = 10.0,
    tol: float = 1e-9,
    max_iter: int = 60,
    fd_step: float = 1e-7,
) -> HeteroclinicConnection:
    """Newton-refine one :class:`ConnectionSeed` via ``correct_connection`` with
    explicit ``tau0`` values (no internal grid scan) and fully-unrestricted
    section sign filters (module docstring, "SECTION CONVENTION")."""
    return correct_connection(
        system,
        node_from,
        node_to,
        k_u=seed.k_u,
        k_s=seed.k_s,
        epsilon=epsilon,
        branch_u=seed.branch_u,
        branch_s=seed.branch_s,
        tau_u0=seed.tau_u,
        tau_s0=seed.tau_s,
        ydot_sign_u=None,
        ydot_sign_s=None,
        x_sign_u=None,
        x_sign_s=None,
        max_time_factor=max_time_factor,
        tol=tol,
        max_iter=max_iter,
        fd_step=fd_step,
    )


@dataclass(frozen=True)
class ConnectionEvidence:
    """Full verification battery for one converged connection (module
    docstring: the honest evidence bundle replacing a nonexistent published
    target). ``passed`` is the conjunction of every hard gate below."""

    full_state_gap: float  # planar 4-state |y_u - y_s| at the matched crossing
    ydot_signs_match: bool  # hard gate (module docstring, "THE ydot-SIGN GATE")
    ydot_u: float
    ydot_s: float
    ghost_distance_from: float  # min (x, xdot) distance to node_from's own points
    ghost_distance_to: float  # min (x, xdot) distance to node_to's own points
    radau_gap: float  # independent-Radau re-derivation: section gap
    radau_vs_dop853: float  # worst |radau crossing - stored crossing|
    backward_distance: float  # crossing state, propagated back t_u, vs unstable seed
    forward_distance: float  # crossing state, propagated fwd |t_s|, vs stable seed
    t_u: float  # unstable-leg manifold transit time (nd, > 0)
    t_s: float  # stable-leg manifold transit time (nd, < 0)
    crossing_state: NDArray[np.float64]  # full 6-state at the crossing (u-leg)
    jacobi_drift_u: float  # |C(crossing) - C(seed)| -- true integration drift, u-leg
    jacobi_drift_s: float
    seed_jacobi_offset_u: float  # |C(seed) - node jacobi| -- the O(eps^2) construction offset
    seed_jacobi_offset_s: float
    passed: bool
    notes: str = ""


#: Gate ceilings for :func:`verify_connection`, calibrated GENEROUSLY relative
#: to the observed evidence floor at the primary C=2.60 hit (observed:
#: full-state gap 5.2e-6, Radau gap 3.7e-6, backward 5.3e-8, forward 3.5e-2)
#: while still far below any coincidental-close-pass scale. The e-5..e-6 floor
#: itself is integrator-tolerance-limited, not corrector-limited: a ~1e-11..1e-12
#: local error amplified by the leg's own Floquet growth exp(t*ln|lambda|/T)
#: ~ 1e5 over these ~40-nd legs lands exactly in this band (same analysis as
#: `#783`'s conditioning-wall arithmetic, here benign because the legs are
#: re-derived, not shot through).
FULL_STATE_GAP_TOL = 1e-4
RADAU_GAP_TOL = 1e-4
BACKWARD_TOL = 1e-5
FORWARD_TOL = 0.5
JACOBI_DRIFT_TOL = 1e-9

#: Ceiling factor for the SEED's own Jacobi offset from the node's energy
#: surface, gated at ``max(JACOBI_DRIFT_TOL, SEED_JACOBI_OFFSET_FACTOR * eps^2)``.
#: A linear manifold seed ``x0 + eps*v`` sits O(eps^2) OFF the energy surface by
#: construction -- the saddle eigenvector is energy-orthogonal to FIRST order
#: (energy invariance of the linearized period map: ``grad(C) . M v = grad(C) . v``,
#: so ``(lambda - 1) grad(C) . v = 0``), leaving only the quadratic term. Measured
#: directly this task: ~6e-12 at eps=0.5e-5 and ~2.6e-9 at eps=1e-4 -- the clean
#: 400x = (20x)^2 quadratic scaling -- so this is a known construction property,
#: NOT integration error, and is gated separately from the true per-leg drift
#: ``|C(crossing) - C(seed)|`` (which stays at :data:`JACOBI_DRIFT_TOL`).
SEED_JACOBI_OFFSET_FACTOR = 100.0


def verify_connection(
    system: cr3bp.CR3BPSystem,
    node_from: jrc.ResonantNode,
    node_to: jrc.ResonantNode,
    conn: HeteroclinicConnection,
    *,
    epsilon: float = EPSILON,
    max_time_factor: float = 10.0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> ConnectionEvidence:
    """Run the full self-consistency battery on a converged connection.

    (a) FULL-4-STATE gap at the matched crossing (both legs re-derived at
    tightened tolerances via ``jrc._full_state_crossing`` with the connection's
    own stored ``tau``/``k``/``branch`` -- includes the ydot components the 2-D
    section residual cannot see, and the ydot-sign hard gate);
    (b) ghost-guard distances against BOTH orbits' own all-sign section points;
    (c) independent-Radau re-derivation of both legs' crossings (this project's
    standing second-integrator discipline);
    (d) forward/backward re-approach (the crossing state retraces the unstable
    leg back to its epsilon-scale seed, and runs forward to the stable leg's);
    (e) Jacobi drift along both legs, measured against each SEED's own Jacobi
    (the trajectory's true conserved value), with the seed's own O(eps^2)
    energy-surface offset gated separately
    (:data:`SEED_JACOBI_OFFSET_FACTOR` -- see that constant's docstring).

    Raises ``ValueError`` if ``conn`` is unconverged or a leg fails to re-reach
    its own crossing (a regression signal, never silently absorbed).
    """
    if not conn.converged:
        raise ValueError("verify_connection requires a converged connection")
    max_time = max_time_factor * max(node_from.period, node_to.period)

    seed_u = _seed_on_manifold(
        system,
        node_from,
        tau=conn.tau_u,
        direction="unstable",
        branch=conn.branch_u,
        epsilon=epsilon,
        rtol=rtol,
        atol=atol,
    )
    seed_s = _seed_on_manifold(
        system,
        node_to,
        tau=conn.tau_s,
        direction="stable",
        branch=conn.branch_s,
        epsilon=epsilon,
        rtol=rtol,
        atol=atol,
    )
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
        raise ValueError("verify_connection: a leg did not re-reach its own crossing")
    t_u, y_u = hit_u
    t_s, y_s = hit_s
    planar = [0, 1, 3, 4]
    full_gap = float(np.linalg.norm((y_u - y_s)[planar]))
    ydot_u, ydot_s = float(y_u[4]), float(y_s[4])
    signs_match = bool(np.sign(ydot_u) == np.sign(ydot_s))

    own_from = ntc.own_section_points(system, node_from)
    own_to = ntc.own_section_points(system, node_to)
    d_from = jrc._ghost_distance(conn.crossing_xv, own_from)
    d_to = jrc._ghost_distance(conn.crossing_xv, own_to)

    pu_r = _section_crossing(
        system,
        seed_u,
        direction="unstable",
        k=conn.k_u,
        max_time=max_time,
        ydot_sign=None,
        x_sign=None,
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    ps_r = _section_crossing(
        system,
        seed_s,
        direction="stable",
        k=conn.k_s,
        max_time=max_time,
        ydot_sign=None,
        x_sign=None,
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    if pu_r is None or ps_r is None:
        radau_gap = float("inf")
        radau_vs = float("inf")
    else:
        radau_gap = float(np.linalg.norm(pu_r - ps_r))
        radau_vs = max(
            float(np.linalg.norm(pu_r - conn.crossing_xv)),
            float(np.linalg.norm(ps_r - conn.crossing_xv)),
        )

    sol_b = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, -t_u),
        y_u,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=abs(t_u) / 5000.0,
    )
    d_back = float(np.linalg.norm(sol_b.y[:, -1] - seed_u))
    sol_f = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, abs(t_s)),
        y_u,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=abs(t_s) / 5000.0,
    )
    d_fwd = float(np.linalg.norm(sol_f.y[:, -1] - seed_s))

    c_seed_u = float(cr3bp.jacobi_constant(seed_u, system.mu))
    c_seed_s = float(cr3bp.jacobi_constant(seed_s, system.mu))
    dj_u = abs(float(cr3bp.jacobi_constant(y_u, system.mu)) - c_seed_u)
    dj_s = abs(float(cr3bp.jacobi_constant(y_s, system.mu)) - c_seed_s)
    off_u = abs(c_seed_u - node_from.jacobi)
    off_s = abs(c_seed_s - node_to.jacobi)
    offset_gate = max(JACOBI_DRIFT_TOL, SEED_JACOBI_OFFSET_FACTOR * epsilon * epsilon)

    gates = {
        "full_state_gap": full_gap <= FULL_STATE_GAP_TOL,
        "ydot_signs_match": signs_match,
        "ghost_from": d_from > GHOST_GUARD_DELTA,
        "ghost_to": d_to > GHOST_GUARD_DELTA,
        "radau_gap": radau_gap <= RADAU_GAP_TOL,
        "radau_vs_dop853": radau_vs <= RADAU_GAP_TOL,
        "backward": d_back <= BACKWARD_TOL,
        "forward": d_fwd <= FORWARD_TOL,
        "jacobi_u": dj_u <= JACOBI_DRIFT_TOL,
        "jacobi_s": dj_s <= JACOBI_DRIFT_TOL,
        "seed_jacobi_offset_u": off_u <= offset_gate,
        "seed_jacobi_offset_s": off_s <= offset_gate,
    }
    failed = [name for name, ok in gates.items() if not ok]
    return ConnectionEvidence(
        full_state_gap=full_gap,
        ydot_signs_match=signs_match,
        ydot_u=ydot_u,
        ydot_s=ydot_s,
        ghost_distance_from=d_from,
        ghost_distance_to=d_to,
        radau_gap=radau_gap,
        radau_vs_dop853=radau_vs,
        backward_distance=d_back,
        forward_distance=d_fwd,
        t_u=t_u,
        t_s=t_s,
        crossing_state=np.asarray(y_u, dtype=np.float64),
        jacobi_drift_u=dj_u,
        jacobi_drift_s=dj_s,
        seed_jacobi_offset_u=off_u,
        seed_jacobi_offset_s=off_s,
        passed=not failed,
        notes="" if not failed else "failed gates: " + ", ".join(failed),
    )


@dataclass(frozen=True)
class FreeTransferResult:
    """End-to-end result of :func:`find_free_transfer` at one shared C."""

    jacobi: float
    connection: HeteroclinicConnection | None
    evidence: ConnectionEvidence | None
    n_seeds: int
    n_refined: int
    n_converged: int
    notes: str = ""


def _find_free_transfer_between(
    system: cr3bp.CR3BPSystem,
    node_from: jrc.ResonantNode,
    node_to: jrc.ResonantNode,
    c: float,
    *,
    n_tau: int,
    n_periods: float,
    max_seed_distance: float,
    max_refine: int,
    seed_diversity_radius: float,
    epsilon: float,
    max_time_factor: float,
    tol: float,
) -> FreeTransferResult:
    """Shared search core: ``Wu(node_from) -> Ws(node_to)`` at Jacobi ``c``.

    Factored out of :func:`find_free_transfer` (`#822`) so the REVERSE
    direction (:func:`find_free_transfer_reverse`, `#840`) reuses the exact
    same crossing-set / seed-pairing / Newton-refine / verify pipeline with
    only the two nodes' roles swapped -- never a re-implementation that could
    silently drift from the forward direction's own logic.
    """
    crossings_u: dict[int, list[ManifoldCrossing]] = {}
    crossings_s: dict[int, list[ManifoldCrossing]] = {}
    for branch in (+1, -1):
        crossings_u[branch] = manifold_section_crossings(
            system,
            node_from,
            direction="unstable",
            branch=branch,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )
        crossings_s[branch] = manifold_section_crossings(
            system,
            node_to,
            direction="stable",
            branch=branch,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )

    seeds: list[ConnectionSeed] = []
    for bu in (+1, -1):
        for bs in (+1, -1):
            seeds.extend(
                find_connection_seeds(
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
        conn = refine_connection(
            system,
            node_from,
            node_to,
            seed,
            epsilon=epsilon,
            max_time_factor=max_time_factor,
            tol=tol,
        )
        if not conn.converged:
            continue
        n_converged += 1
        ev = verify_connection(
            system, node_from, node_to, conn, epsilon=epsilon, max_time_factor=max_time_factor
        )
        if ev.passed:
            return FreeTransferResult(
                jacobi=float(c),
                connection=conn,
                evidence=ev,
                n_seeds=len(seeds),
                n_refined=n_refined,
                n_converged=n_converged,
            )
    return FreeTransferResult(
        jacobi=float(c),
        connection=None,
        evidence=None,
        n_seeds=len(seeds),
        n_refined=n_refined,
        n_converged=n_converged,
        notes=(
            "no verified connection among refined seeds -- honest negative at this "
            "C/settings (conditional on n_tau/n_periods/max_refine, not proof of absence)"
        ),
    )


def find_free_transfer(
    system: cr3bp.CR3BPSystem,
    c: float,
    *,
    n_tau: int = 36,
    n_periods: float = 9.0,
    max_seed_distance: float = 0.02,
    max_refine: int = 8,
    seed_diversity_radius: float = 0.05,
    epsilon: float = EPSILON,
    max_time_factor: float = 10.0,
    tol: float = 1e-9,
) -> FreeTransferResult:
    """End-to-end free-transfer search Wu(2:1) -> Ws(3:1) at one shared ``c``.

    Builds both nodes, both manifolds' whole crossing sets (both branches),
    pairs same-ydot-sign close approaches, Newton-refines up to ``max_refine``
    seeds (closest first), and returns the FIRST candidate that both converges
    AND passes the full :func:`verify_connection` battery. An honest negative
    (``connection=None``) if none does -- never a fabricated pass.

    ``seed_diversity_radius``: a seed whose ``(x, xdot)`` lies within this
    radius of an already-refined seed's is skipped -- without this, the
    distance-ranked top of the seed list is often ``n_tau``-adjacent
    near-duplicates of ONE close approach, and a fragile geometry there burns
    the whole ``max_refine`` budget on the same failure (observed directly at
    C=2.61: 1174 seeds, 8 refined, 0 converged, all eight the same crossing
    neighborhood). Radius 0.05 is well above the seed-pair gap scale (<= 0.02)
    and well below the distinct-intersection separation scale (>~ 0.1 in this
    problem's own converged hits).
    """
    node21 = build_vaquero_overlap_node(system, 2, c)
    node31 = build_vaquero_overlap_node(system, 3, c)
    return _find_free_transfer_between(
        system,
        node21,
        node31,
        c,
        n_tau=n_tau,
        n_periods=n_periods,
        max_seed_distance=max_seed_distance,
        max_refine=max_refine,
        seed_diversity_radius=seed_diversity_radius,
        epsilon=epsilon,
        max_time_factor=max_time_factor,
        tol=tol,
    )


def find_free_transfer_reverse(
    system: cr3bp.CR3BPSystem,
    c: float,
    *,
    n_tau: int = 36,
    n_periods: float = 9.0,
    max_seed_distance: float = 0.02,
    max_refine: int = 8,
    seed_diversity_radius: float = 0.05,
    epsilon: float = EPSILON,
    max_time_factor: float = 10.0,
    tol: float = 1e-9,
) -> FreeTransferResult:
    """End-to-end free-transfer search Wu(3:1) -> Ws(2:1) at one shared ``c``
    -- the REVERSE of :func:`find_free_transfer` (`#840`).

    CR3BP time-reversal symmetry ``(x, y, xdot, ydot, t) -> (x, -y, -xdot,
    ydot, -t)`` maps ``Wu(A) ∩ Ws(B)`` connections into ``Wu(B) ∩
    Ws(A)`` ones, guaranteeing this direction's existence given a forward
    hit at the same ``c`` -- demonstrated directly here (never merely
    asserted), the same discipline `#822` Sec. 4's own C=2.60 reverse demo
    established. Identical pipeline to :func:`find_free_transfer`
    (:func:`_find_free_transfer_between`, shared verbatim) with the two
    nodes' roles swapped; an honest negative (``connection=None``) if no
    refined seed passes, never a forced pass.
    """
    node21 = build_vaquero_overlap_node(system, 2, c)
    node31 = build_vaquero_overlap_node(system, 3, c)
    return _find_free_transfer_between(
        system,
        node31,
        node21,
        c,
        n_tau=n_tau,
        n_periods=n_periods,
        max_seed_distance=max_seed_distance,
        max_refine=max_refine,
        seed_diversity_radius=seed_diversity_radius,
        epsilon=epsilon,
        max_time_factor=max_time_factor,
        tol=tol,
    )


def predict_reverse_crossing(forward_crossing_xv: tuple[float, float]) -> tuple[float, float]:
    """CR3BP time-reversal symmetry image of a forward ``{y=0}`` crossing.

    The map ``(x, y, xdot, ydot, t) -> (x, -y, -xdot, ydot, -t)`` sends a
    forward ``Wu(2:1) -> Ws(3:1)`` crossing at ``(x, xdot)`` (with ``y=0`` on
    the section) to its REVERSE ``Wu(3:1) -> Ws(2:1)`` counterpart's own
    crossing at ``(x, -xdot)`` -- both Vaquero overlap-band orbits are
    already x-axis symmetric (built by
    :func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`),
    so this prediction needs no search, only the already-computed forward
    connection. `#840` used this to RANK candidate seeds for the
    reverse-direction search at C=2.54/2.66, rather than a blind ``(x,
    xdot)``-close-approach scan. At C=2.54 the verified reverse crossing
    matched this prediction to ~6e-9, confirming it is genuinely the
    time-reversal image of the SAME forward connection. At C=2.66 the
    verified reverse crossing is a DIFFERENT intersection point entirely
    (~0.26 away from the prediction, module docstring "the intersection
    structure is rich" `#822` already documented) -- the prediction only
    biased the search toward a promising region, it did not guarantee
    landing on the specific mirror-image crossing; existence in the
    reverse direction is still a genuine, independently verified connection
    at the same C and mu either way.
    """
    x, xdot = forward_crossing_xv
    return (float(x), -float(xdot))


def find_free_transfer_reverse_near_prediction(
    system: cr3bp.CR3BPSystem,
    c: float,
    predicted_xv: tuple[float, float],
    *,
    n_tau: int = 36,
    n_periods: float = 9.0,
    epsilons: tuple[float, ...] = (1.0e-4, 2.0e-4, 3.0e-4, 5.0e-4, 1.0e-3),
    max_try_per_epsilon: int = 5,
    max_time_factor: float = 10.0,
    tol: float = 1e-9,
) -> FreeTransferResult:
    """Reverse-direction (Wu(3:1) -> Ws(2:1)) search at ``c``, ranking seeds
    by proximity to ``predicted_xv`` (:func:`predict_reverse_crossing`)
    instead of by native ``(x, xdot)`` close-approach distance, and sweeping
    ``epsilons`` in order.

    `#840`'s own diagnosed pathology at C=2.66: the blind
    :func:`find_free_transfer_reverse` scan at the sibling modules' default
    ``epsilon=1e-4`` converged 6/11 refined candidates, but EVERY ONE failed
    ``verify_connection``'s forward-reapproach gate (the SAME evidence-quality
    tradeoff `#822`'s own docstring documents for the forward direction --
    "Manifold-offset ε as an evidence-quality control"), including the
    closest-to-predicted candidate itself (``forward_distance=1.58`` vs the
    ``0.5`` ceiling). Raising ``epsilon`` shortens each leg (fewer Floquet
    amplification periods); at C=2.66 a candidate among the ``epsilon=2e-4``
    round's top-``max_try_per_epsilon`` predicted-ranked seeds then passes
    cleanly (the sibling default ``1e-4`` already sufficed at C=2.54, on the
    very first/closest-to-predicted candidate). This function promotes that
    diagnosed fix into reusable code rather than leaving it as a one-off
    retry.

    An honest negative (``connection=None``) if no candidate across the
    swept ``epsilons`` passes the full unmodified :func:`verify_connection`
    battery -- gates are never loosened, only the SEARCH is retargeted.
    """
    node21 = build_vaquero_overlap_node(system, 2, c)
    node31 = build_vaquero_overlap_node(system, 3, c)
    pred_x, pred_xdot = predicted_xv
    total_seeds_tried = 0
    total_converged = 0
    for epsilon in epsilons:
        crossings_u: dict[int, list[ManifoldCrossing]] = {}
        crossings_s: dict[int, list[ManifoldCrossing]] = {}
        for branch in (+1, -1):
            crossings_u[branch] = manifold_section_crossings(
                system,
                node31,
                direction="unstable",
                branch=branch,
                n_tau=n_tau,
                n_periods=n_periods,
                epsilon=epsilon,
            )
            crossings_s[branch] = manifold_section_crossings(
                system,
                node21,
                direction="stable",
                branch=branch,
                n_tau=n_tau,
                n_periods=n_periods,
                epsilon=epsilon,
            )
        seeds: list[ConnectionSeed] = []
        for bu in (+1, -1):
            for bs in (+1, -1):
                seeds.extend(
                    find_connection_seeds(
                        crossings_u[bu], crossings_s[bs], branch_u=bu, branch_s=bs
                    )
                )
        seeds.sort(key=lambda s: math.hypot(s.x - pred_x, s.xdot - pred_xdot))
        for seed in seeds[:max_try_per_epsilon]:
            total_seeds_tried += 1
            conn = refine_connection(
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
            total_converged += 1
            ev = verify_connection(
                system, node31, node21, conn, epsilon=epsilon, max_time_factor=max_time_factor
            )
            if ev.passed:
                return FreeTransferResult(
                    jacobi=float(c),
                    connection=conn,
                    evidence=ev,
                    n_seeds=total_seeds_tried,
                    n_refined=total_seeds_tried,
                    n_converged=total_converged,
                )
    return FreeTransferResult(
        jacobi=float(c),
        connection=None,
        evidence=None,
        n_seeds=total_seeds_tried,
        n_refined=total_seeds_tried,
        n_converged=total_converged,
        notes=(
            "no verified connection near the symmetry-predicted crossing across the "
            "swept epsilon range -- honest negative at these settings"
        ),
    )


__all__ = [
    "BACKWARD_TOL",
    "EPSILON",
    "FORWARD_TOL",
    "FULL_STATE_GAP_TOL",
    "GHOST_GUARD_DELTA",
    "JACOBI_DRIFT_TOL",
    "OVERLAP_C_RANGE",
    "OVERLAP_GRID_ICS",
    "RADAU_GAP_TOL",
    "SEED_JACOBI_OFFSET_FACTOR",
    "ConnectionEvidence",
    "ConnectionSeed",
    "FreeTransferResult",
    "ManifoldCrossing",
    "build_vaquero_overlap_node",
    "find_connection_seeds",
    "find_free_transfer",
    "find_free_transfer_reverse",
    "find_free_transfer_reverse_near_prediction",
    "manifold_section_crossings",
    "predict_reverse_crossing",
    "refine_connection",
    "verify_connection",
]
