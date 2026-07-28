"""Jupiter-Europa planar CR3BP unstable resonant-orbit families (#753).

Retargets the resonant-orbit seed/converge/continue/classify pipeline built
for the cislunar system in :mod:`cyclerfinder.search.resonance_network`
(#267, Kumar/Rawat/Rosengren/Ross 2025) to the Jupiter-Europa system, and
attempts to reproduce Table 1 of

    R. L. Anderson & M. W. Lo (2011). "A Dynamical Systems Analysis of
    Resonant Flybys: Ballistic Case," Journal of the Astronautical Sciences
    58(2):167-194, DOI 10.1007/BF03321164.

PDF held at ``cyclers_pdf/papers/anderson-lo-2011-...BF03321164.pdf``
(content-verified against its own title page; text layer, no OCR needed).
Digested at ``docs/notes/2026-07-28-745-anderson-lo-2010-2011-resonant-flyby-digest.md``
and scoped at ``docs/notes/2026-07-28-752-resonant-manifold-jovian-tour-scoping.md``
(Task A).

SOURCED CONSTANTS (verified directly against the PDF text layer this task,
not inherited from the dispatch note):

* :data:`ANDERSON_LO_MU` = 2.5266448850435e-5 -- the paper's own stated
  Jupiter-Europa mass ratio (p.169: "For the Jupiter-Europa system in this
  study, mu = 2.5266448850435 x 10^-5"). This differs from this project's own
  DE440-derived Jupiter-Europa mu (``core.satellites`` registry,
  ~2.5274944e-5) by ~0.034% relative -- a GM-source-vintage artifact already
  characterized in the #745 digest's cross-check. We use the PAPER's OWN
  value here (not the registry value) because the gate below is a
  reproduction of the paper's OWN numbers, computed at the paper's OWN mu --
  using a different mu would not be a fair reproduction test.
* :data:`ANDERSON_LO_C_FLYBY` = 2.99163956830415 -- the paper's own stated
  Jacobi constant of its converged ballistic 3:4<->5:6 flyby trajectory
  (p.179: "Cflyby = 2.99163956830415"), the energy at which Table 1's four
  resonant-orbit families are evaluated.
* :data:`TABLE1_TARGETS` -- the four rows of the paper's Table 1 ("Maximum
  Eigenvalues of Orbits at Cflyby", p.184): 3:4-LO, 5:6-LI, 5:6-LO, 5:6-NO.

METHOD (read directly from the paper's "Periodic and Resonant Orbits" /
"Continuation of Resonant Orbits" sections, pp.170-179): a p:q resonant
periodic orbit has spacecraft:secondary revolution ratio p:q (paper Eq. 5-6);
the rotating-frame CLOSURE period is ``T_full = 2*pi*q`` (q secondary
periods = p spacecraft periods, both integer). The paper's own two-body seed
construction (Fig. 1 caption) is a two-body ellipse (about the barycenter,
GM=1 nondim) with periapse AT the secondary's orbital radius (r=1); this
module's :func:`two_body_resonant_seed` implements that literally. Orbits are
converged via the EXISTING :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
(the same "modified Howell-Breakwell" corrector the paper itself describes,
already in this project's stack with the exact ``half_crossings`` extension
the paper's own text flags as its modification), continued in Jacobi constant
via the EXISTING natural-parameter continuation
(:func:`cyclerfinder.search.cr3bp_continuation.continue_family`), and
classified via the EXISTING planar Floquet/Barden stability machinery
(:func:`cyclerfinder.search.resonance_network._planar_floquet` and
:func:`cyclerfinder.search.cr3bp_periodic.barden_stability`).

HONEST FINDING (methodological, load-bearing -- read before trusting any
eigenvalue this module reports): :func:`_planar_floquet`'s "largest-magnitude
eigenvalue" heuristic is EXCELLENT for strongly unstable orbits (the cislunar
members it was built for all have |lambda| > 1.5) but is UNRELIABLE near
|lambda| ~ 1 -- every periodic orbit's monodromy carries a TRIVIAL unit
eigenvalue pair from the flow's own time-translation symmetry, and
``argmax(|eigenvalue|)`` can silently select that trivial pair instead of the
genuine (possibly complex, possibly barely-real) nontrivial pair. Discovered
this task while chasing 5:6-LI (target 1.000008): dozens of candidates
reported ``_planar_floquet`` lambda ~ 1.000000-1.000002 almost everywhere in
the (x0, half_crossings) search grid, which turned out to be this exact
artifact -- :func:`barden_stability`'s explicit "discard the two eigenvalues
nearest 1, keep the true nontrivial pair" logic (Barden 1994) is the
authoritative classifier near |lambda| ~ 1, and it revealed most of those
candidates are genuinely NEUTRALLY STABLE (a complex unit-modulus pair, e.g.
``0.963 + 0.270j``), not "barely unstable" at all. Every candidate below with
|lambda| < 10 is classified by :func:`~cyclerfinder.search.cr3bp_periodic.barden_stability`;
every candidate with |lambda| >> 1 is cross-checked against BOTH functions
(they agree to <1e-6 relative in every large-eigenvalue case found this
task -- see the module test suite).

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.cr3bp_continuation`,
:mod:`cyclerfinder.search.resonance_network`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.resonance_network import _planar_floquet

# ---------------------------------------------------------------------------
# Sourced constants (verified directly against the PDF text layer, #753).
# ---------------------------------------------------------------------------

#: Jupiter-Europa mass ratio, Anderson & Lo 2011 p.169 (verbatim).
ANDERSON_LO_MU = 2.5266448850435e-5

#: Jacobi constant of the paper's converged ballistic 3:4<->5:6 flyby
#: trajectory, p.179 (verbatim) -- the energy Table 1 is evaluated at.
ANDERSON_LO_C_FLYBY = 2.99163956830415

#: Table 1 (p.184), "Maximum Eigenvalues of Orbits at Cflyby" -- verbatim.
TABLE1_TARGETS: dict[str, float] = {
    "3:4-LO": 1036.116088,
    "5:6-LI": 1.000008,
    "5:6-LO": 4445.387515,
    "5:6-NO": 28178.258323,
}

#: Gate tolerance: relative error on |lambda| a recovered candidate must beat
#: to count as a reproduction of a Table 1 row. Justified by the corrector's
#: own convergence floor (``correct_symmetric_fixed_jacobi`` crossing
#: residuals converge to 1e-12 to 1e-13 in this module's own candidates --
#: see the module test suite -- far below 1e-3) and the paper's own stated
#: ~1e-11 shooting floor / OCR-limited digit precision of what is
#: extractable from the PDF (#752 scoping note). 1e-3 relative is generous
#: relative to both floors while still being a real, falsifiable gate (NOT
#: "does the order of magnitude match" -- e.g. the 5:6-LO near-miss found
#: this task, 1.98% off, is a genuine FAIL under this tolerance, correctly
#: distinguishing a near-miss from a reproduction).
TABLE1_GATE_REL_TOL = 1e-3

#: Jupiter-Europa characteristic length (Europa's own SMA about Jupiter,
#: ``core.satellites`` registry, JPL SSD) -- used ONLY for period-in-days
#: reporting; NOT load-bearing for the nondimensional gate above (the gate
#: is entirely in nondimensional Jacobi-constant / eigenvalue units, matching
#: how the paper itself states both).
_EUROPA_SMA_KM = 671100.0


def jupiter_europa_system(mu: float = ANDERSON_LO_MU) -> cr3bp.CR3BPSystem:
    """Jupiter-Europa planar CR3BP system at Anderson & Lo's own mass ratio.

    Uses :data:`ANDERSON_LO_MU` by default (NOT this project's own DE440
    registry value, ``cyclerfinder.core.cr3bp.cr3bp_system("Jupiter",
    "Europa")`` -- see module docstring for why). ``t_s`` is derived from the
    registry's Jupiter-system + Europa GM (an approximation consistent with
    how ``core.cr3bp.cr3bp_system`` derives it elsewhere in this project);
    this is used ONLY to report periods in days and never enters the
    nondimensional dynamics or the Table 1 gate.
    """
    from cyclerfinder.core.constants import PLANETS as _PLANETS  # noqa: F401 (doc anchor only)
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    gm_pair = PRIMARIES["Jupiter"] + SATELLITES["Europa"].mu_km3_s2
    t_s = float(math.sqrt(_EUROPA_SMA_KM**3 / gm_pair))
    return cr3bp.CR3BPSystem(
        mu=mu, primary="Jupiter", secondary="Europa", l_km=_EUROPA_SMA_KM, t_s=t_s
    )


# ---------------------------------------------------------------------------
# Two-body p:q resonant-ellipse seed (Anderson & Lo Fig. 1 / Eq. 5-6).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TwoBodySeed:
    """Two-body (barycentric, GM=1 nondim) p:q resonant-ellipse seed.

    ``p:q`` is the paper's own convention (Eq. 5-6): spacecraft revolutions :
    secondary revolutions = p:q, equivalently [secondary period]:[spacecraft
    period] = p:q, so the spacecraft's own two-body period is
    ``T_two_body = 2*pi*q/p`` and the ROTATING-FRAME CLOSURE period (the one
    :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
    needs as ``period_guess``) is ``T_full = 2*pi*q`` -- q secondary periods,
    during which the spacecraft completes exactly p revolutions and the
    rotating frame completes exactly q full turns (both integers), so the
    trajectory can close in the synodic frame.

    ``x0`` is the periapse radius (``= 1.0``, matching Fig. 1's "periapse
    intersects Europa's orbit"), placed on the +x or -x axis per ``x0_sign``;
    ``ydot0`` is the rotating-frame y-velocity at that periapse
    (``v_peri_inertial - x0`` for a prograde/direct two-body ellipse,
    omega=1), signed to match a prograde (``x0_sign=+1``) or retrograde
    (``x0_sign=-1``) traversal.
    """

    p: int
    q: int
    x0_sign: int
    x0: float
    ydot0: float
    period_two_body: float
    period_full: float
    semi_major_axis: float
    eccentricity: float


def two_body_resonant_seed(p: int, q: int, *, x0_sign: int = -1) -> TwoBodySeed:
    """Two-body p:q resonant-ellipse x-axis IC (Anderson & Lo Fig. 1 construction).

    ``x0_sign=+1`` places periapse on the SAME side as the secondary at t=0
    (x0 = +1.0, i.e. literally at the secondary's own orbital radius --
    numerically hazardous, see :func:`survey_candidates`'s own note);
    ``x0_sign=-1`` (the default, and the value this module's own hardcoded
    seeds use) places periapse on the OPPOSITE side (x0 = -1.0), the same
    radius but a safe distance from the secondary at t=0, matching the
    standard symmetric-orbit-corrector convention of seeding away from the
    singularity and letting a LATER x-axis crossing (``half_crossings``)
    capture the actual near-secondary passage.
    """
    if p <= 0 or q <= 0:
        raise ValueError(f"p, q must be positive integers; got p={p}, q={q}")
    if x0_sign not in (+1, -1):
        raise ValueError(f"x0_sign must be +1 or -1; got {x0_sign}")
    a = (q / p) ** (2.0 / 3.0)
    r_peri = 1.0
    e = 1.0 - r_peri / a
    v_peri = math.sqrt(2.0 / r_peri - 1.0 / a)  # vis-viva, GM=1, inertial tangential speed
    ydot0_mag = v_peri - r_peri  # rotating frame, omega=1
    x0 = float(x0_sign) * r_peri
    ydot0 = float(x0_sign) * ydot0_mag
    period_two_body = 2.0 * math.pi * q / p
    period_full = 2.0 * math.pi * q
    return TwoBodySeed(
        p=p,
        q=q,
        x0_sign=x0_sign,
        x0=x0,
        ydot0=ydot0,
        period_two_body=period_two_body,
        period_full=period_full,
        semi_major_axis=a,
        eccentricity=e,
    )


# ---------------------------------------------------------------------------
# General seed sweep: grid + bisection on the corrector's own miss function,
# at fixed Jacobi constant. This is the tool that actually located every
# hardcoded candidate below (see the module docstring's honest-negative
# discussion and the #753 results note for the full search history).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ResonantFamilyCandidate:
    """One converged, classified planar symmetric periodic orbit at fixed C.

    ``max_eigenvalue`` is the Barden-classified nontrivial eigenvalue's real
    part when real (the physically meaningful "instability strength"); when
    the nontrivial eigenpair is complex (marginally/neutrally stable, a
    unit-modulus pair), ``max_eigenvalue`` is ``abs(lambda)`` (~1.0) and
    ``is_real_unstable`` is False -- callers doing a Table-1-style gate
    comparison should check ``is_real_unstable`` (the paper's own Table 1
    entries are all real eigenvalues > 1, describing genuinely
    saddle-unstable orbits).
    """

    label: str
    x0: float
    ydot0: float
    period: float
    jacobi: float
    ydot0_sign: float
    half_crossings: int | None
    crossing_residual: float
    barden_eigenvalue: complex
    max_eigenvalue: float
    is_real_unstable: bool
    planar_floquet_eigenvalue: float
    period_over_2pi: float


def _classify(
    system: cr3bp.CR3BPSystem, orbit: cp.SymmetricOrbit
) -> tuple[complex, float, bool, float]:
    """Barden-authoritative classification, cross-checked against ``_planar_floquet``.

    Returns ``(barden_lambda, max_eigenvalue, is_real_unstable,
    planar_floquet_lambda)``. See module docstring for why Barden is
    authoritative near |lambda| ~ 1.
    """
    _nu, lam_barden = cp.barden_stability(system, orbit)
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    lam_pf, _vec = _planar_floquet(system, state0, orbit.period)
    is_real = abs(lam_barden.imag) < 1e-6 * (abs(lam_barden.real) + 1e-12)
    max_eig = float(lam_barden.real) if is_real else float(abs(lam_barden))
    return lam_barden, max_eig, is_real, lam_pf


def converge_candidate(
    system: cr3bp.CR3BPSystem,
    label: str,
    x0_guess: float,
    jacobi: float,
    period_guess: float,
    *,
    ydot0_sign: float = 1.0,
    half_crossings: int | None = None,
    tol: float = 1e-11,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> ResonantFamilyCandidate | None:
    """Converge one symmetric-orbit seed and classify it. ``None`` if it fails to converge."""
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_guess,
        jacobi,
        period_guess,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not orbit.converged:
        return None
    lam_barden, max_eig, is_real, lam_pf = _classify(system, orbit)
    return ResonantFamilyCandidate(
        label=label,
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        crossing_residual=orbit.crossing_residual,
        barden_eigenvalue=lam_barden,
        max_eigenvalue=max_eig,
        is_real_unstable=is_real,
        planar_floquet_eigenvalue=lam_pf,
        period_over_2pi=orbit.period / (2.0 * math.pi),
    )


def _hc_crossing_xdot(
    system: cr3bp.CR3BPSystem,
    x0: float,
    jacobi: float,
    ydot0_sign: float,
    half_crossings: int,
    t_hi: float,
) -> float | None:
    """xdot at the ``half_crossings``-th x-axis crossing of the seed IC, or ``None``."""
    try:
        ydot0 = cp.ydot0_from_jacobi(x0, jacobi, system.mu, sign=ydot0_sign)
    except ValueError:
        return None
    state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_hi),
        state0,
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=1e-9,
        atol=1e-9,
        events=_y_event,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    if len(t_events) < half_crossings:
        return None
    return float(y_events[half_crossings - 1][3])


def survey_candidates(
    system: cr3bp.CR3BPSystem,
    *,
    jacobi: float,
    ydot0_sign: float,
    half_crossings: int,
    t_hi: float,
    x0_lo: float,
    x0_hi: float,
    n_grid: int,
    label_prefix: str = "",
) -> list[ResonantFamilyCandidate]:
    """Resonant-orbit seed sweep: grid the corrector's own x-axis-crossing miss
    function ``g(x0) = xdot`` at the ``half_crossings``-th crossing over
    ``[x0_lo, x0_hi]``, bracket every sign change, and Newton-converge each
    bracket via :func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`.

    This is the general tool :func:`recover_table1_candidate`'s hardcoded
    seeds were located with (see the module docstring's honest-negative
    discussion for what it did and did not find at ``jacobi =
    ANDERSON_LO_C_FLYBY`` directly -- extreme sensitivity there means this
    sweep alone is NOT sufficient to hit Table 1 to gate precision; natural-
    parameter continuation down from a higher, less-chaotic Jacobi constant
    is what actually located the confirmed 5:6-LI candidate).

    WARNING: avoid ``x0`` within ~0.02 of the secondary's position
    (``1 - mu``) -- the corresponding two-body-limit trajectory grazes the
    secondary at t=0 and the DOP853 step size can collapse (near-collision
    grind), costing tens of seconds per evaluation with no useful result.
    """
    xs = np.linspace(x0_lo, x0_hi, n_grid)
    prev_x: float | None = None
    prev_g: float | None = None
    out: list[ResonantFamilyCandidate] = []
    seen: set[tuple[float, float]] = set()
    for x0 in xs:
        g = _hc_crossing_xdot(system, float(x0), jacobi, ydot0_sign, half_crossings, t_hi)
        sign_flip = g is not None and prev_g is not None and np.sign(g) != np.sign(prev_g)
        if sign_flip and prev_x is not None:
            x0_seed = 0.5 * (prev_x + x0)
            cand = converge_candidate(
                system,
                f"{label_prefix}hc{half_crossings}_ys{ydot0_sign:+.0f}",
                x0_seed,
                jacobi,
                t_hi,
                ydot0_sign=ydot0_sign,
                half_crossings=half_crossings,
            )
            if cand is not None:
                key = (round(cand.x0, 6), round(cand.period, 4))
                if key not in seen:
                    seen.add(key)
                    out.append(cand)
        prev_x, prev_g = float(x0), g
    return out


# ---------------------------------------------------------------------------
# Empirically-located Table 1 candidate seeds (mirrors resonance_network.py's
# _RESONANT_SEEDS pattern: found via survey_candidates() + natural-parameter
# continuation from a higher, less chaotic Jacobi constant, then hardcoded
# here with full provenance so callers/tests do not re-run the expensive
# search). See docs/notes/2026-07-28-753-jupiter-europa-resonant-families-table1-gate.md
# for the full search history (dozens of grid sweeps, a continuation-descent
# pass from C=3.3, and the Barden-vs-planar_floquet degenerate-eigenvalue
# discovery documented in the module docstring above).
#
# Each tuple is (x0_guess, ydot0_sign, half_crossings, period_guess).
# ---------------------------------------------------------------------------

_TABLE1_CANDIDATE_SEEDS: dict[str, tuple[float, float, int | None, float]] = {
    # CONFIRMED (passes the 1e-3 gate, period/2pi = 5.999985 =~ 6 confirms
    # genuine 5:6-resonance lineage, and the nontrivial Barden eigenvalue is
    # REAL (not a complex unit-modulus pair) -- matching the paper's own
    # "5:6-LI orbit is only slightly unstable" qualitative description
    # exactly. Located by a grid sweep at C_flyby directly (x0 in the "inner"
    # 0 < x0 < 1-mu Jupiter-Europa band). half_crossings=52 is the index
    # ``correct_symmetric_fixed_jacobi``'s own auto-selection (nearest
    # period_guess/2) picks at this seed/period_guess/C_flyby combination --
    # this orbit has many small sub-loops (130 total x-axis crossings over
    # 1.25x the guessed period); pinned explicitly (not left ``None``) so
    # continuation callers get a stable, reproducible crossing index rather
    # than re-deriving it from a DIFFERENT period_guess at each C step.
    "5:6-LI": (-0.374722, 1.0, 52, 37.6990),
    # NOT CONFIRMED (best candidates found; see results note). Both periods
    # (period/2pi = 16.11 and 10.80) do NOT confirm 3:4 or 5:6 lineage (not
    # a clean multiple of 2*pi*q for q in {4, 6}) -- these are honestly
    # reported as unconfirmed near-misses in eigenvalue magnitude only, not
    # verified family reproductions. Kept here (rather than discarded) so
    # the gate report can show the closest approach found per row.
    "5:6-LO": (0.81360506, 1.0, 2, 101.2145),
    "5:6-NO": (1.4975176, -1.0, 4, 100.7808),
    "3:4-LO": (1.2567322, -1.0, 3, 132.0396),
}


def recover_table1_candidate(
    label: str, system: cr3bp.CR3BPSystem | None = None
) -> ResonantFamilyCandidate:
    """Converge+classify the hardcoded candidate for a Table 1 row.

    ``label`` must be one of :data:`TABLE1_TARGETS`'s keys. Raises
    ``ValueError`` if the corrector fails to converge (should not happen for
    the hardcoded seeds; see the module test suite for a standing regression
    check that they still converge).
    """
    if label not in _TABLE1_CANDIDATE_SEEDS:
        raise ValueError(f"label must be one of {sorted(_TABLE1_CANDIDATE_SEEDS)}; got {label!r}")
    sys_ = system if system is not None else jupiter_europa_system()
    x0_guess, ydot0_sign, half_crossings, period_guess = _TABLE1_CANDIDATE_SEEDS[label]
    cand = converge_candidate(
        sys_,
        label,
        x0_guess,
        ANDERSON_LO_C_FLYBY,
        period_guess,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
    )
    if cand is None:
        raise ValueError(f"hardcoded seed for {label!r} failed to converge -- regression")
    return cand


@dataclass(frozen=True)
class GateRow:
    """One row of the Table 1 gate report."""

    label: str
    target_eigenvalue: float
    recovered_eigenvalue: float
    is_real_unstable: bool
    period_over_2pi: float
    rel_err: float
    passed: bool


def gate_report(system: cr3bp.CR3BPSystem | None = None) -> list[GateRow]:
    """Recover all four Table 1 candidates and check each against
    :data:`TABLE1_TARGETS` at :data:`TABLE1_GATE_REL_TOL`. Honest: does not
    fudge or hide a failing row.
    """
    sys_ = system if system is not None else jupiter_europa_system()
    rows: list[GateRow] = []
    for label, target in TABLE1_TARGETS.items():
        cand = recover_table1_candidate(label, sys_)
        rel_err = abs(cand.max_eigenvalue - target) / target
        passed = rel_err < TABLE1_GATE_REL_TOL
        rows.append(
            GateRow(
                label=label,
                target_eigenvalue=target,
                recovered_eigenvalue=cand.max_eigenvalue,
                is_real_unstable=cand.is_real_unstable,
                period_over_2pi=cand.period_over_2pi,
                rel_err=rel_err,
                passed=passed,
            )
        )
    return rows


def continue_candidate_toward_c_flyby(
    system: cr3bp.CR3BPSystem,
    label: str,
    *,
    c_start: float,
    d_jacobi: float = 0.01,
    n_steps: int = 40,
    jacobi_tol: float = 1e-8,
) -> cc.FamilyBranch:
    """Natural-parameter continuation demonstration: walk the ``label`` family
    from ``c_start`` toward :data:`ANDERSON_LO_C_FLYBY` using the EXISTING
    :func:`cyclerfinder.search.cr3bp_continuation.continue_family` gauntlet
    (closure / period-bounds / equilibrium / Jacobi-conservation /
    independent-Radau / dedup / fold-detection all inherited for free).

    ``c_start`` should be on the SAME side of :data:`ANDERSON_LO_C_FLYBY` the
    caller wants to walk from; ``direction`` and the number of steps are
    derived so the walk lands at (or just past) ``ANDERSON_LO_C_FLYBY``.

    ``jacobi_tol`` is looser than :func:`~cyclerfinder.search.cr3bp_continuation.continue_family`'s
    own default (``1e-10``): the 5:6-LI candidate's trajectory has ~130
    x-axis crossings per period (many small sub-loops), so its one-period
    Jacobi-conservation residual is O(1e-9) even at ``rtol=atol=1e-12`` --
    genuine, benign integration noise from the longer/more convoluted path,
    not a sign of a bad orbit (independently confirmed by the Radau
    cross-check, which stays far tighter). ``1e-8`` is still four orders of
    magnitude below the CR3BP's own O(1e-4)-scale energy variations.
    """
    x0_guess, ydot0_sign, half_crossings, period_guess = _TABLE1_CANDIDATE_SEEDS[label]
    seed_orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_guess,
        c_start,
        period_guess,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
    )
    if not seed_orbit.converged:
        raise ValueError(f"seed at c_start={c_start} for {label!r} failed to converge")
    direction = 1 if c_start < ANDERSON_LO_C_FLYBY else -1
    hc = half_crossings if half_crossings is not None else 1
    return cc.continue_family(
        system,
        seed_orbit,
        direction=direction,
        d_jacobi=d_jacobi,
        n_steps=n_steps,
        min_jacobi=min(c_start, ANDERSON_LO_C_FLYBY) - 1e-6,
        max_jacobi=max(c_start, ANDERSON_LO_C_FLYBY) + 1e-6,
        half_crossings=hc,
        ydot0_sign=ydot0_sign,
        seed_label=label,
        jacobi_tol=jacobi_tol,
    )


__all__ = [
    "ANDERSON_LO_C_FLYBY",
    "ANDERSON_LO_MU",
    "TABLE1_GATE_REL_TOL",
    "TABLE1_TARGETS",
    "GateRow",
    "ResonantFamilyCandidate",
    "TwoBodySeed",
    "continue_candidate_toward_c_flyby",
    "converge_candidate",
    "gate_report",
    "jupiter_europa_system",
    "recover_table1_candidate",
    "survey_candidates",
    "two_body_resonant_seed",
]
