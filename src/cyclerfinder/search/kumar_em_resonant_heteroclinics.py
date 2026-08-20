"""`#827` Kumar et al. (2026) Table-5 3:1 -> 2:1 Earth-Moon heteroclinic reproduction.

Digit-grade reproduction of the PRINTED 3:1 -> 2:1 heteroclinic manifold-
intersection states of

    B. Kumar, A. Rawat, A. J. Rosengren, S. D. Ross (2026). "Cislunar Resonant
    Transport and Heteroclinic Pathways: From 3:1 to 2:1 to L1," *Advances in
    Space Research* 77(3):3815-3845, DOI 10.1016/j.asr.2025.12.005
    (= arXiv:2509.12675v2; corpus ``kumar-2025-arxiv-2509.12675.pdf``; digest
    ``docs/notes/2026-06-20-digest-kumar-2025.md``),

using `#822`'s Poincare-section Newton machinery
(:mod:`cyclerfinder.search.vaquero_em_cycler_connections`, reused UNCHANGED --
its search/verify functions are node-generic) generalized to Kumar's own family
branches: their 3:1/2:1 unstable resonant orbits at ``C_J in {3.00, 3.05, 3.10,
3.15}``, OUTSIDE Vaquero's ``[2.54, 2.66]`` overlap band.

WHAT THE SOURCE PRINTS (unlike Vaquero 2013's prose-only existence claim):
Table 5 (Appendix 8.2) prints EXACT intersection states ``x = (x, y, xdot,
ydot)`` of Wu(3:1) and Ws(2:1) -- "accurately-computed manifold intersections"
from their Section-3.4 segment-intersection + bisection algorithm on their
perigee Poincare section -- at seven Jacobi constants (Type 1 "short" at
C in {2.54, 2.70, 2.86, 3.00, 3.05}; Type 2 "long", via a 5:2 intermediary
segment, at C in {3.10, 3.15}). Table 6 prints the exact converged ``(x,
ydot)`` (with ``y = xdot = 0``) of every parent resonant orbit. Their mass
ratio is printed in Section 2.1: ``mu = 1.2150584270572e-2`` (differs from this
project's Earth-Moon ``mu`` at the ~1e-10 absolute level, which matters at
digit grade, so :func:`kumar_system` builds the system with THEIR value).

THE SECTION MISMATCH AND HOW THE GATE IS DEFINED: their intersection states
live on the PERIGEE section (Eq. 3: ``sigma = (x + mu)*xdot + y*ydot``
crossing negative -> positive; equivalently ``rdot . r = 0`` about the Earth
with the osculating mean anomaly at 0, i.e. geocentric-distance minimum) --
NOT on this project's ``{y=0}`` Newton section. A true manifold intersection
is an epsilon-convention-independent POINT of the actual heteroclinic orbit,
so the reproduction gate is section-independent: converge the connection with
`#822`'s ``{y=0}`` machinery, then propagate the verified matched crossing
state along both legs, collect EVERY perigee-section crossing, and require
the minimum planar-4-state distance to Kumar's printed ``x`` to be at most
:data:`KUMAR_MATCH_TOL` -- with the runner-up crossing's distance far above
it (a specific-point match, not a "somewhere near the tube" coincidence).

TARGETED SEEDING (the honest use of the published digits): Kumar's printed
``x`` selects WHICH intersection to converge -- his state is propagated to its
first few ``{y=0}`` crossings, each is matched to the nearest same-ydot-sign
member of the two manifolds' whole crossing sets (`#822`'s
``manifold_section_crossings``), and the resulting ``(tau, k, branch)`` pairs
seed ``correct_connection``. The Newton residual is this project's own section
gap -- the published state is never part of the residual, so the subsequent
digit comparison is a genuine two-sided check (a wrong printed state could not
"converge to itself").

NODES: re-derived at call time from Table 6's own printed ICs by the same
fixed-Jacobi symmetric corrector the sibling modules use, at Kumar's own
``mu``; the corrector's converged IC must agree with the PRINT to
:data:`TABLE6_IC_ABS_TOL` (measured ~1e-13 -- the print is the staleness
reference here, a STRONGER anchor than `#822`'s archived-run guesses). Node
saddles here are strong: |lambda| ~ 12-17 (3:1) and ~180-333 (2:1), so
manifolds peel off within a few periods and ``n_periods = 5`` horizons
suffice (vs `#822`'s 9) -- but the 2:1's extreme saddle also amplifies the
verification legs' own error floor, which is why this module's default
manifold offset is :data:`KUMAR_EPSILON` (1e-4), not the sibling 0.5e-5
(see that constant's docstring).

LITERATURE-NOVELTY GATE: not re-run live (no catalogue writeback here, and
`#822`'s own live-WebSearch mandatory-floor run against this exact paper's
connection concept already returned "published", confidence 0.95) -- this is
an explicit REPRODUCTION of Kumar et al.'s own printed table; nothing here is
novel and nothing is claimed novel. See the results note
(``docs/notes/2026-08-21-827-kumar-table5-reproduction.md``) for the full
account, including why `#839`'s C=3.13 is NOT covered by this module's
printed-table rows.

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.vaquero_em_cycler_connections` (`#822` machinery,
reused unchanged), :mod:`cyclerfinder.search.jovian_resonant_connections`
(``ResonantNode``), :mod:`cyclerfinder.search.cr3bp_periodic` (fixed-Jacobi
corrector), :mod:`cyclerfinder.genome.heteroclinic_cycle` (Floquet pair).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.vaquero_em_cycler_connections as vcc
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    _planar_floquet_pair,
)

#: Kumar et al.'s own Earth-Moon mass ratio -- PRINT provenance (Section 2.1,
#: "we use mu = 1.2150584270572e-2"; repeated verbatim in Section 3.2). NOT this
#: project's mass-data-derived Earth-Moon mu (0.01215058439469525) -- the two
#: differ by ~1.24e-10 absolute, which is above the digit-grade comparison
#: floor, so every reproduction in this module runs at THEIR value.
KUMAR_MU = 1.2150584270572e-2

#: Table 5 (Appendix 8.2), the "3:1 to 2:1" block -- PRINT provenance,
#: transcribed verbatim from the corpus PDF (all seven printed intersection
#: rows; every row is labeled point type ``x``, an "accurately-computed
#: manifold intersection"). ``C -> (transfer_type, (x, y, xdot, ydot))``.
#: Rows 2.54-3.05 are the paper's Type 1 (short, direct) transfers; 3.10 and
#: 3.15 are Type 2 (long, via a 5:2 intermediary segment) -- Type 1
#: intersections cease to exist for C >= 3.09 (their Section 5.2.1).
KUMAR_TABLE5_31_TO_21: dict[float, tuple[int, tuple[float, float, float, float]]] = {
    2.54: (1, (-0.02495763718, -0.02013755597, 7.56003664579, -4.80802086624)),
    2.70: (1, (-0.03022725673, 0.042561712914, -5.82508502316, -2.47401118974)),
    2.86: (1, (0.06136961660, -0.00485126788, 0.322456526264, 4.886777873750)),
    3.00: (1, (-0.07187159167, -0.09858383084, 3.22059393928, -1.95100061370)),
    3.05: (1, (0.029351172720, 0.320801654785, -1.77006016626, 0.228990735501)),
    3.10: (2, (0.173726091985, 0.150836454723, -1.44182292916, 1.77676712380)),
    3.15: (2, (0.203387734653, 0.165418603215, -1.25080421783, 1.62978185751)),
}

#: The `#827` reproduction target set: the four C values of the paper's own
#: Table-1 worked summary (C = 3.00, 3.05, 3.10, 3.15) -- the same four the
#: task registration names. Table 6 prints BOTH parent-orbit ICs at exactly
#: these C values (it does NOT at 2.70/2.86, whose 3:1/2:1 columns it also
#: covers -- so all seven rows above are reproducible in principle; the four
#: below are the registered digit-grade gate).
KUMAR_REPRODUCTION_CS: tuple[float, ...] = (3.00, 3.05, 3.10, 3.15)

#: Table 6 (Appendix 8.2) -- PRINT provenance, transcribed verbatim from the
#: corpus PDF (the complete printed table: 4:1 x2, 3:1 x7, 2:1 x8 rows).
#: ``(m, C) -> (x, ydot)`` for the m:1 unstable resonant orbit at Jacobi C;
#: all rows have ``y = 0, xdot = 0`` (the table's own header note).
KUMAR_TABLE6_ICS: dict[tuple[int, float], tuple[float, float]] = {
    # 4:1 family
    (4, 2.85): (0.770658112628, -0.616166626758),
    (4, 3.15): (0.737385941470, -0.355888785284),
    # 3:1 family
    (3, 2.54): (0.901330167142, -0.846224994375),
    (3, 2.70): (0.888102851739, -0.725919236562),
    (3, 2.86): (0.868787031779, -0.584481100346),
    (3, 3.00): (0.845409258968, -0.434953051104),
    (3, 3.05): (0.834875967543, -0.37200499129),
    (3, 3.10): (0.822429022871, -0.300987128481),
    (3, 3.15): (0.806804382814, -0.218228896847),
    # 2:1 family
    (2, 2.54): (0.964310487216, -1.20233213851),
    (2, 2.628): (0.959722043439, -1.090831532610),
    (2, 2.70): (0.954308268237, -0.989690975727),
    (2, 2.86): (0.934036118105, -0.743014738588),
    (2, 3.00): (0.904862874177, -0.515850946602),
    (2, 3.05): (0.892104772476, -0.429597033392),
    (2, 3.10): (0.878280334961, -0.334629543419),
    (2, 3.15): (0.864401165205, -0.219058827351),
}

#: Node re-derivation must land back on Table 6's own printed digits to this
#: absolute tolerance in BOTH (x0, ydot0) -- measured directly at build time:
#: ~1e-13 at all eight 3:1/2:1 targets (the corrector at Kumar's own mu
#: re-converges onto their print at machine-noise level). The ceiling is set
#: 1e-9: far above the measured floor, far below the print's last digit scale
#: for a MIS-transcribed row (>= 1e-9 would only pass if the transcription
#: itself is faithful to ~9 digits).
TABLE6_IC_ABS_TOL = 1e-9

#: Manifold-offset magnitude for THIS family -- `#822`'s own documented
#: "epsilon as an evidence-quality control" calibration (its note, Sec. 2
#: finding 1), re-derived here: at the sibling default 0.5e-5 the verified
#: C=3.00 connection's honest 2.07e-6 crossing floor, amplified by the 2:1
#: node's extreme saddle (|lambda| ~ 180-333) over the ~2.4-period stable leg,
#: puts the forward re-approach at 0.87 -- OVER the unchanged 0.5 gate
#: (2.07e-6 * 180^2.4 ~ 0.5, the same arithmetic, measured directly this
#: task). At 1e-4 the legs shorten (~1.8 stable periods), the whole battery
#: floor drops to 5.6e-8..5.8e-4, and every gate passes with margin. The
#: gates themselves are `#822`'s, unchanged.
KUMAR_EPSILON = 1e-4

#: Digit-grade match ceiling for the perigee-state comparison (planar 4-state
#: Euclidean norm ours-vs-printed). Calibration: the printed states' own
#: self-consistency floor is ~1e-8 (their Jacobi vs the labeled C differs by
#: 7.5e-10..7.0e-8 across the four targets -- print-truncation scale), and this
#: module's own evidence floor is the integrator/Floquet-amplification floor of
#: the verified legs (same arithmetic as `#822` Sec. 2). Measured achieved
#: distances at the four targets: 1.19e-7 (C=3.00), 1.35e-6 (3.05), 1.75e-6
#: (3.10), 4.49e-7 (3.15); runner-up perigee crossings sit at 4.8e-2..1.6, so
#: 1e-4 is >~ 50x above the achieved floor and >~ 500x below the nearest
#: wrong-point scale -- a pass is a specific-point identification, not a
#: proximity coincidence. The ACHIEVED distance is recorded per row.
KUMAR_MATCH_TOL = 1e-4

#: How many of Kumar's-state {y=0} crossings (forward + backward) to use as
#: matching-point candidates, and how far out to look for them.
_N_TARGET_CROSSINGS_FWD = 4
_N_TARGET_CROSSINGS_BWD = 2

#: Seed-gap ceiling: a target {y=0} crossing whose nearest same-ydot-sign
#: manifold-set crossing is farther than this (in (x, xdot)) is not usable as
#: a Newton seed pair (the phase grid missed the neighborhood -- raise n_tau).
_MAX_SEED_GAP = 0.05


def kumar_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP system at Kumar et al.'s OWN printed mass ratio.

    ``l_km``/``t_s`` (used only for unit conversion in reporting, never in the
    nondimensional dynamics) are taken from this project's own Earth-Moon
    system record.
    """
    base = cr3bp.cr3bp_system("Earth", "Moon")
    return cr3bp.CR3BPSystem(
        mu=KUMAR_MU, primary=base.primary, secondary=base.secondary, l_km=base.l_km, t_s=base.t_s
    )


def kumar_table5_state6(c: float) -> NDArray[np.float64]:
    """Kumar's printed Table-5 intersection state at ``c``, lifted to a planar
    6-state ``(x, y, 0, xdot, ydot, 0)``. Raises ``KeyError`` off the table."""
    key = round(float(c), 4)
    if key not in KUMAR_TABLE5_31_TO_21:
        raise KeyError(
            f"C={c} has no printed Table-5 3:1->2:1 intersection; printed rows: "
            f"{sorted(KUMAR_TABLE5_31_TO_21)}"
        )
    _ttype, (x, y, xd, yd) = KUMAR_TABLE5_31_TO_21[key]
    return np.array([x, y, 0.0, xd, yd, 0.0], dtype=np.float64)


def build_kumar_node(
    system: cr3bp.CR3BPSystem,
    m: int,
    c: float,
    *,
    tol: float = 1e-12,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    ic_abs_tol: float = TABLE6_IC_ABS_TOL,
) -> jrc.ResonantNode:
    """Re-derive one Table-6 m:1 unstable resonant orbit and wrap it as a node.

    Re-converges at exactly Jacobi ``c`` with the fixed-Jacobi symmetric
    corrector (guess = Table 6's own printed ``x``; period guess 6.4, the
    family band's own ~2*pi scale -- the paper's Figures 7/9/12 print periods
    6.31-6.56 TU), then REQUIRES the converged ``(x0, ydot0)`` to agree with
    the PRINTED Table-6 row to ``ic_abs_tol`` (:data:`TABLE6_IC_ABS_TOL`) --
    the reproduce-before-trust gate, anchored on the source itself. The
    Floquet saddle pair is recomputed from the full-period planar monodromy;
    a non-saddle (``|lambda| <= 1.05``) is rejected loudly.

    Raises ``KeyError`` off the printed table, ``RuntimeError``/``ValueError``
    on any convergence/print-agreement/saddle failure -- never a silently
    degraded node.
    """
    key = (m, round(float(c), 4))
    if key not in KUMAR_TABLE6_ICS:
        raise KeyError(
            f"({m}:1, C={c}) is not a printed Table-6 row; rows: {sorted(KUMAR_TABLE6_ICS)}"
        )
    x0_print, ydot0_print = KUMAR_TABLE6_ICS[key]
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_print,
        float(c),
        6.4,
        ydot0_sign=-1.0,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"node corrector failed for {m}:1 C={c}: residual {orbit.crossing_residual:.3e}"
        )
    dx0 = abs(orbit.x0 - x0_print)
    dydot0 = abs(orbit.ydot0 - ydot0_print)
    if dx0 > ic_abs_tol or dydot0 > ic_abs_tol:
        raise ValueError(
            f"{m}:1 C={c}: re-converged IC disagrees with Kumar Table 6's print: "
            f"|dx0|={dx0:.3e}, |dydot0|={dydot0:.3e} (> {ic_abs_tol:.0e})"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    lam_u, v_u, _lam_s, v_s = _planar_floquet_pair(system, state0, orbit.period)
    if lam_u <= 1.05:
        raise ValueError(f"{m}:1 C={c}: |lambda|={lam_u:.4f} is not materially unstable")
    return jrc.ResonantNode(
        label=f"kumar-{m}1-c{round(float(c) * 100):d}",
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=True,
    )


def _trajectory_events(
    system: cr3bp.CR3BPSystem,
    state6: NDArray[np.float64],
    t_end: float,
    event_fn: object,
    *,
    method: str = "DOP853",
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[tuple[float, NDArray[np.float64]]]:
    """All ``event_fn`` crossings of the trajectory from ``state6`` over
    ``(0, t_end)`` (``t_end`` may be negative), as ``(t, state6)`` rows."""
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, float(t_end)),
        np.asarray(state6, dtype=np.float64),
        args=(system.mu,),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        events=event_fn,
        max_step=abs(float(t_end)) / 1000.0,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    return [
        (float(t), np.asarray(y, dtype=np.float64))
        for t, y in zip(t_events, y_events, strict=False)
    ]


def y0_section_crossings(
    system: cr3bp.CR3BPSystem,
    state6: NDArray[np.float64],
    t_end: float,
    **kwargs: object,
) -> list[tuple[float, NDArray[np.float64]]]:
    """All ``{y=0}`` crossings (both signs) over ``(0, t_end)``, excluding any
    root at the start point itself."""

    def _ev(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _ev.terminal = False  # type: ignore[attr-defined]
    _ev.direction = 0.0  # type: ignore[attr-defined]
    rows = _trajectory_events(system, state6, t_end, _ev, **kwargs)  # type: ignore[arg-type]
    t_floor = 1e-9 * abs(float(t_end))
    return [(t, y) for t, y in rows if abs(t) > t_floor]


def perigee_section_crossings(
    system: cr3bp.CR3BPSystem,
    state6: NDArray[np.float64],
    t_end: float,
    **kwargs: object,
) -> list[tuple[float, NDArray[np.float64]]]:
    """All PERIGEE-section crossings over ``(0, t_end)`` -- Kumar's own section:
    ``sigma = (x + mu)*xdot + y*ydot`` (Eq. 3, the geocentric ``rdot . r``)
    crossing negative -> positive along the trajectory's own time arrow (the
    event direction is flipped for a backward integration, where solve_ivp
    sees the sign pattern reversed)."""
    mu = system.mu

    def _ev(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return (float(y[0]) + mu) * float(y[3]) + float(y[1]) * float(y[4])

    _ev.terminal = False  # type: ignore[attr-defined]
    _ev.direction = 1.0 if float(t_end) > 0.0 else -1.0  # type: ignore[attr-defined]
    return _trajectory_events(system, state6, t_end, _ev, **kwargs)  # type: ignore[arg-type]


@dataclass(frozen=True)
class Table5Reproduction:
    """One C value's digit-grade reproduction outcome."""

    jacobi: float
    transfer_type: int  # Kumar's own Type 1 (short) / Type 2 (long) label
    connection: HeteroclinicConnection | None
    evidence: vcc.ConnectionEvidence | None
    match_distance: float  # min planar-4-state |ours - printed| over perigee crossings
    matched_state: NDArray[np.float64] | None  # our perigee state at the match (6-state)
    matched_t: float  # leg time of the matched perigee crossing (from the {y=0} crossing)
    runner_up_distance: float  # next-nearest perigee crossing's distance (separation scale)
    n_candidates: int
    n_refined: int
    n_converged: int
    matched: bool  # match_distance <= KUMAR_MATCH_TOL and battery passed
    notes: str = ""


def _match_against_print(
    system: cr3bp.CR3BPSystem,
    evidence: vcc.ConnectionEvidence,
    x_print6: NDArray[np.float64],
) -> tuple[float, NDArray[np.float64], float, float]:
    """Min planar-4-state distance from the verified connection's perigee
    crossings (both legs, re-derived from the matched ``{y=0}`` crossing
    state) to the printed intersection state. Returns ``(best_distance,
    best_state6, best_t, runner_up_distance)``."""
    planar = [0, 1, 3, 4]
    y_cross = evidence.crossing_state
    rows = perigee_section_crossings(
        system, y_cross, -abs(evidence.t_u)
    ) + perigee_section_crossings(system, y_cross, abs(evidence.t_s))
    if not rows:
        return float("inf"), y_cross, 0.0, float("inf")
    dists = [float(np.linalg.norm((y - x_print6)[planar])) for _t, y in rows]
    order = np.argsort(dists)
    i = int(order[0])
    runner_up = float(dists[int(order[1])]) if len(rows) > 1 else float("inf")
    return dists[i], rows[i][1], rows[i][0], runner_up


def reproduce_table5_intersection(
    system: cr3bp.CR3BPSystem,
    c: float,
    *,
    n_tau: int = 48,
    n_periods: float = 5.0,
    epsilon: float = KUMAR_EPSILON,
    max_refine: int = 8,
    tol: float = 1e-9,
    max_time_factor: float = 10.0,
) -> Table5Reproduction:
    """Digit-grade reproduction of Kumar's printed intersection at one ``c``.

    Pipeline (module docstring, "TARGETED SEEDING"): build both Table-6 nodes;
    build Wu(3:1)/Ws(2:1) whole ``{y=0}`` crossing sets (both branches, `#822`
    machinery); propagate the PRINTED state to its first few ``{y=0}``
    crossings; pair each with the nearest same-ydot-sign crossing of each
    manifold set; Newton-refine candidates in seed-gap order
    (``correct_connection`` -- the printed state is NOT in the residual);
    verify with `#822`'s full battery; and gate the winner's perigee-section
    states against the print (:data:`KUMAR_MATCH_TOL`). Honest negative
    (``matched=False`` with the achieved distances recorded) if nothing both
    verifies and matches -- never a fabricated pass.
    """
    key = round(float(c), 4)
    ttype, _vals = KUMAR_TABLE5_31_TO_21[key]
    x_print6 = kumar_table5_state6(key)
    node31 = build_kumar_node(system, 3, key)
    node21 = build_kumar_node(system, 2, key)

    sets_u: dict[int, list[vcc.ManifoldCrossing]] = {}
    sets_s: dict[int, list[vcc.ManifoldCrossing]] = {}
    for branch in (+1, -1):
        sets_u[branch] = vcc.manifold_section_crossings(
            system,
            node31,
            direction="unstable",
            branch=branch,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )
        sets_s[branch] = vcc.manifold_section_crossings(
            system,
            node21,
            direction="stable",
            branch=branch,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )

    targets = (
        y0_section_crossings(system, x_print6, 2.5 * node21.period)[:_N_TARGET_CROSSINGS_FWD]
        + y0_section_crossings(system, x_print6, -1.5 * node31.period)[:_N_TARGET_CROSSINGS_BWD]
    )

    candidates: list[tuple[float, vcc.ConnectionSeed]] = []
    for _t0, y0 in targets:
        tgt_x, tgt_xd = float(y0[0]), float(y0[3])
        sgn = np.sign(float(y0[4]))
        for bu in (+1, -1):
            cu_list = [cr for cr in sets_u[bu] if np.sign(cr.ydot) == sgn]
            if not cu_list:
                continue
            cu = min(cu_list, key=lambda cr: np.hypot(cr.x - tgt_x, cr.xdot - tgt_xd))
            du = float(np.hypot(cu.x - tgt_x, cu.xdot - tgt_xd))
            for bs in (+1, -1):
                cs_list = [cr for cr in sets_s[bs] if np.sign(cr.ydot) == sgn]
                if not cs_list:
                    continue
                cs = min(cs_list, key=lambda cr: np.hypot(cr.x - tgt_x, cr.xdot - tgt_xd))
                ds = float(np.hypot(cs.x - tgt_x, cs.xdot - tgt_xd))
                gap = max(du, ds)
                if gap > _MAX_SEED_GAP:
                    continue
                candidates.append(
                    (
                        gap,
                        vcc.ConnectionSeed(
                            distance=gap,
                            branch_u=bu,
                            branch_s=bs,
                            k_u=cu.k,
                            k_s=cs.k,
                            tau_u=cu.tau,
                            tau_s=cs.tau,
                            x=cu.x,
                            xdot=cu.xdot,
                            ydot=cu.ydot,
                        ),
                    )
                )
    candidates.sort(key=lambda t: t[0])

    n_refined = 0
    n_converged = 0
    best_dist = float("inf")
    seen: set[tuple[int, int, int, int]] = set()
    for _gap, seed in candidates:
        if n_refined >= max_refine:
            break
        sig = (seed.branch_u, seed.branch_s, seed.k_u, seed.k_s)
        if sig in seen:
            continue
        seen.add(sig)
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
        if not ev.passed:
            continue
        dist, state_m, t_m, runner_up = _match_against_print(system, ev, x_print6)
        best_dist = min(best_dist, dist)
        if dist <= KUMAR_MATCH_TOL:
            return Table5Reproduction(
                jacobi=key,
                transfer_type=ttype,
                connection=conn,
                evidence=ev,
                match_distance=dist,
                matched_state=state_m,
                matched_t=t_m,
                runner_up_distance=runner_up,
                n_candidates=len(candidates),
                n_refined=n_refined,
                n_converged=n_converged,
                matched=True,
            )
    return Table5Reproduction(
        jacobi=key,
        transfer_type=ttype,
        connection=None,
        evidence=None,
        match_distance=best_dist,
        matched_state=None,
        matched_t=0.0,
        runner_up_distance=float("inf"),
        n_candidates=len(candidates),
        n_refined=n_refined,
        n_converged=n_converged,
        matched=False,
        notes=(
            "no verified+matched connection among refined candidates -- honest negative "
            "at these settings (conditional on n_tau/n_periods/max_refine, not proof "
            "the printed state is wrong)"
        ),
    )


__all__ = [
    "KUMAR_EPSILON",
    "KUMAR_MATCH_TOL",
    "KUMAR_MU",
    "KUMAR_REPRODUCTION_CS",
    "KUMAR_TABLE5_31_TO_21",
    "KUMAR_TABLE6_ICS",
    "TABLE6_IC_ABS_TOL",
    "Table5Reproduction",
    "build_kumar_node",
    "kumar_system",
    "kumar_table5_state6",
    "perigee_section_crossings",
    "reproduce_table5_intersection",
    "y0_section_crossings",
]
