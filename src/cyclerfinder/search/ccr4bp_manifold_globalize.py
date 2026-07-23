"""CCR4BP manifold globalization (`#694`, Stage B's final sub-task per `#686`'s plan).

`#691` extracts a trustworthy stable/unstable manifold DIRECTION at a single
torus phase (the segment-anchored discrete-QR/CLV eigenvector,
:func:`cyclerfinder.search.ccr4bp_whisker.manifold_direction_segmented_clv`).
This module GLOBALIZES that direction into a genuine manifold: a discretized
"tube" of trajectories, one per phase on a grid of the torus's free angle
``theta2`` (``theta1`` is held at a fixed reference "departure phase" --
see below), each obtained by stepping a small distance ``eps`` off the torus
along the extracted direction and propagating the TRUE nonlinear CCR4BP flow
(`#689`'s :func:`cyclerfinder.core.ccr4bp.ccr4bp_eom`) forward (unstable
branch) or backward (stable branch) in time.

Why ``theta1`` is fixed at departure, not swept
------------------------------------------------
The torus is 2D (``theta1`` x ``theta2``); ``theta1`` is locked to Ganymede's
forcing clock at a FIXED rate ``omega1`` (`#690`'s own construction), so the
whole torus structure -- and therefore its whiskers -- repeats identically
every Ganymede synodic period (``torus.period``). Sweeping a 1D "generating
circle" over ``theta2`` at one representative ``theta1`` (default 0, i.e. the
corrector's own phase origin) and then letting TIME itself run is therefore
sufficient to reach every phase of the torus's own forcing cycle "for free"
as the propagation proceeds -- exactly as one sweeps only the single angle
around a periodic orbit (not also the phase along its own period) when
globalizing a classical CR3BP periodic-orbit manifold. This is a direct
extension of the existing, already-shipped precedent
:func:`cyclerfinder.genome.qp_torus_manifold.torus_manifold_grid`, which
sweeps the SAME kind of "one torus angle held fixed" 1D generating grid for a
QBCP/CR3BP quasi-periodic torus and uses the SAME default perturbation size
(``eps = 1e-6`` nondimensional length/velocity units) adopted here for direct
comparability across this project's manifold-globalization modules.

Dimension count and what "the manifold tube" actually samples
----------------------------------------------------------------
`#686`'s own strategy note frames the eventual search as living in "the 4D
stroboscopic space" and states that two 2D manifolds intersect generically in
POINTS there. That statement is exactly right for the FULL manifold (each
manifold is 2D: 1D from the torus's own ``theta2`` family, 1D from the flow
TIME since departure) intersecting the augmented 5D extended phase space
(``x, y, vx, vy, theta1 mod 2*pi``) restricted to a section: a codimension-1
section of a 5D space leaves a 4D space, and two 2D manifolds inside a 4D
space intersect generically in 0D (points) -- see the dimension count worked
through in this task's own commit-message discussion. The PRACTICAL
consequence: fixing the flow time to an INTEGER multiple of the forcing
period (a literal stroboscopic-map snapshot) freezes out one of the two free
manifold dimensions and leaves each branch's crossing set only 1D at that
instant, so two 1D curves in the full 4D map space would NOT be expected to
meet generically. The correct search therefore keeps flow time CONTINUOUS
(not stroboscopically discretized) on both branches -- this module samples
each departure phase's trajectory continuously (dense ODE output) over a
bounded time window, and the intersection search
(:mod:`cyclerfinder.search.ccr4bp_heteroclinic_search`) treats the elapsed
flow times ``t_u`` (unstable, forward) and ``t_s`` (stable, backward) as
genuinely continuous free unknowns alongside the two phases ``theta2_u``,
``theta2_s`` -- a well-posed 4-unknown / 4-residual (position + velocity)
match, not an artificially frozen stroboscopic slice.

Step-size (``eps``) choice
---------------------------
``eps = 1e-6`` (nondimensional; 1 length unit = Europa's semi-major axis,
671,100 km, so ``eps`` is roughly 0.67 km in physical position units) is the
SAME default this project's existing QBCP/CR3BP torus-manifold globalizer
(:func:`cyclerfinder.genome.qp_torus_manifold.torus_manifold_grid`) already
uses and has run against real published cases -- not a fresh, unvalidated
choice. It is additionally checked HERE against the CCR4BP's own achieved
torus invariance residual (`#690`'s own reported floor, ~1.6e-4 for the
physical-mass 3:4 torus): ``eps`` is roughly two orders of magnitude BELOW
that floor, so the perturbation is not swamped by the torus's own
representation noise at departure, while
``test_globalize_step_size_linear_regime`` in the test module directly
confirms the propagated tube position at a fixed elapsed time scales
linearly in ``eps`` (doubling ``eps`` doubles the displacement from the
unperturbed torus trajectory to within a few percent) at the phases and
horizons this module's positive control actually uses -- i.e. the chosen
``eps`` is verified to still be in the LINEAR regime at the sampled elapsed
times, not just asserted by analogy.

Scope
-----
Pure post-hoc module: does NOT modify `#689`'s ``core/ccr4bp.py``, `#690`'s
``search/variational_ccr4bp_torus.py``, or `#691`'s
``search/ccr4bp_whisker.py`` -- only calls their public functions. No
catalogue writeback (a capability-proof + discovery-attempt module, not a
vetted result).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import OdeSolution, solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
from cyclerfinder.search.ccr4bp_whisker import manifold_direction_segmented_clv
from cyclerfinder.search.variational_ccr4bp_torus import CCR4BPTorusVariationalResult

# The project-wide default manifold-departure step size (matches
# genome/qp_torus_manifold.torus_manifold_grid's own eps default -- see
# module docstring "Step-size (eps) choice").
DEFAULT_EPS = 1e-6


def _embed6(u4: NDArray[np.float64]) -> NDArray[np.float64]:
    """Embed the CCR4BP torus's planar 4-state into the full 6-state (z=vz=0).

    Same convention as ``ccr4bp_whisker._embed6`` (duplicated locally, not
    imported, to keep this module decoupled from that module's private API --
    `#691`'s module is not modified and its private helpers are not a
    supported cross-module contract)."""
    return np.array([u4[0], u4[1], 0.0, u4[2], u4[3], 0.0], dtype=np.float64)


@dataclass(frozen=True)
class ManifoldTube:
    """A discretized manifold "tube": one departure per ``theta2`` grid point,
    each propagated continuously in time and sampled at ``t_samples``.

    ``states[i, j]`` is the full 6-state at elapsed flow time ``t_samples[j]``
    since departure from torus phase ``theta2_grid[i]`` (forward in time for
    ``branch="unstable"``, backward for ``branch="stable"`` -- so
    ``t_samples`` is always a NON-NEGATIVE "elapsed time since departure"
    grid regardless of branch; the sign of actual propagation direction is
    branch-dependent and handled internally). ``base_states[i]`` /
    ``directions[i]`` are the torus state and unit CLV direction the
    departure point was built from (``base_states[i] + lobe_sign * eps *
    directions[i]``). ``valid[i]`` is ``False`` where the CLV extraction or
    the propagation failed (no data at that ``theta2`` -- ``states[i]`` is
    filled with NaN).
    """

    torus: CCR4BPTorusVariationalResult
    branch: str
    lobe_sign: float
    theta1_section: float
    eps: float
    n_segments_dir: int
    theta2_grid: NDArray[np.float64]
    t_samples: NDArray[np.float64]
    states: NDArray[np.float64]
    base_states: NDArray[np.float64]
    directions: NDArray[np.float64]
    valid: NDArray[np.bool_]


def _direction_at(
    torus: CCR4BPTorusVariationalResult,
    theta1_section: float,
    theta2: float,
    branch: str,
    ref_vec: NDArray[np.float64] | None,
    *,
    n_segments_dir: int,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """One `#691` segment-anchored CLV direction extraction, dropping diagnostics."""
    out = manifold_direction_segmented_clv(
        torus,
        theta1_section,
        theta2,
        branch,
        ref_vec,
        n_segments=n_segments_dir,
        rtol=rtol,
        atol=atol,
    )
    if out is None:
        return None
    return out[0], out[1]


def globalize_manifold_tube(
    torus: CCR4BPTorusVariationalResult,
    branch: str,
    *,
    n_theta2: int = 40,
    t_max_periods: float = 2.0,
    n_time: int = 120,
    eps: float = DEFAULT_EPS,
    lobe_sign: float = 1.0,
    theta1_section: float = 0.0,
    n_segments_dir: int = 24,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> ManifoldTube:
    """Globalize a CCR4BP torus's stable/unstable manifold into a discretized tube.

    For each of ``n_theta2`` departure phases (a grid over ``theta2`` at the
    fixed ``theta1_section``), extracts the `#691` segmented-CLV direction
    (with sign continuity threaded across the grid so the tube does not
    randomly flip lobes), perturbs by ``lobe_sign * eps`` along it, and
    propagates the TRUE nonlinear flow continuously (dense ODE output, one
    ``solve_ivp`` call per departure -- no re-anchoring, unlike `#691`'s
    direction-extraction segments) for ``t_max_periods * torus.period``,
    forward for ``branch="unstable"``, backward for ``branch="stable"``.
    Sampled at ``n_time`` points evenly spaced in elapsed flow time.

    A departure whose CLV extraction or propagation fails is marked invalid
    (:attr:`ManifoldTube.valid`) with NaN-filled states rather than raising --
    a single bad phase (e.g. a locally non-hyperbolic point) should not abort
    the whole tube.
    """
    if branch not in ("unstable", "stable"):
        raise ValueError(f"branch must be 'unstable' or 'stable', got {branch!r}")
    if n_theta2 < 1:
        raise ValueError(f"n_theta2 must be >= 1, got {n_theta2}")
    if n_time < 2:
        raise ValueError(f"n_time must be >= 2, got {n_time}")

    theta2_grid = 2.0 * np.pi * np.arange(n_theta2) / n_theta2
    t_horizon = t_max_periods * torus.period
    t_samples = np.linspace(0.0, t_horizon, n_time)
    t0 = (theta1_section % (2.0 * np.pi)) / torus.omega1
    sign_dir = 1.0 if branch == "unstable" else -1.0

    states = np.full((n_theta2, n_time, 6), np.nan, dtype=np.float64)
    base_states = np.full((n_theta2, 6), np.nan, dtype=np.float64)
    directions = np.full((n_theta2, 6), np.nan, dtype=np.float64)
    valid = np.zeros(n_theta2, dtype=np.bool_)

    ref_vec: NDArray[np.float64] | None = None
    for i, theta2 in enumerate(theta2_grid):
        extracted = _direction_at(
            torus,
            theta1_section,
            float(theta2),
            branch,
            ref_vec,
            n_segments_dir=n_segments_dir,
            rtol=rtol,
            atol=atol,
        )
        if extracted is None:
            continue
        state0, vec = extracted
        ref_vec = vec
        base_states[i] = state0
        directions[i] = vec
        perturbed = state0 + lobe_sign * eps * vec
        sol = solve_ivp(
            ccr4bp.ccr4bp_eom,
            (t0, t0 + sign_dir * t_horizon),
            perturbed,
            args=(torus.system,),
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=True,
        )
        if not sol.success or sol.sol is None:
            continue
        dense: OdeSolution = sol.sol
        for j, s in enumerate(t_samples):
            states[i, j] = dense(t0 + sign_dir * s)
        valid[i] = True

    return ManifoldTube(
        torus=torus,
        branch=branch,
        lobe_sign=lobe_sign,
        theta1_section=theta1_section,
        eps=eps,
        n_segments_dir=n_segments_dir,
        theta2_grid=theta2_grid,
        t_samples=t_samples,
        states=states,
        base_states=base_states,
        directions=directions,
        valid=valid,
    )


def manifold_state_at(
    torus: CCR4BPTorusVariationalResult,
    branch: str,
    theta1_section: float,
    theta2: float,
    t_flow: float,
    *,
    eps: float = DEFAULT_EPS,
    lobe_sign: float = 1.0,
    n_segments_dir: int = 24,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    ref_vec: NDArray[np.float64] | None = None,
) -> NDArray[np.float64] | None:
    """Full-precision single-point re-evaluation of the manifold state at
    CONTINUOUS ``(theta2, t_flow)`` -- the primitive the continuous
    differential-correction refinement in
    :mod:`cyclerfinder.search.ccr4bp_heteroclinic_search` calls repeatedly.

    Unlike :func:`globalize_manifold_tube` (a grid, built once and reused),
    this re-extracts the CLV direction and re-propagates from scratch at the
    EXACT requested ``theta2``/``t_flow`` -- the continuous primitive a
    gradient-based refinement needs. ``t_flow`` is elapsed time since
    departure (forward for ``branch="unstable"``, backward for
    ``branch="stable"``, same convention as :class:`ManifoldTube`).

    ``ref_vec`` fixes the extracted eigenvector's SIGN by dot-product
    continuity against a caller-supplied reference (e.g. a nearby
    :class:`ManifoldTube`'s own threaded ``directions`` entry, or a fixed
    anchor captured once at the start of a refinement loop). Without it, the
    raw eigen-decomposition sign is whatever ``numpy.linalg.eig`` happens to
    return at this exact phase -- deterministic for a FIXED ``theta2``, but
    with NO continuity guarantee against a neighbouring phase's convention
    (unlike :func:`globalize_manifold_tube`'s own grid, which threads
    continuity along ``theta2`` internally). A caller matching a specific
    tube's lobe convention (rather than treating this as a fresh, arbitrary
    branch pick) MUST pass ``ref_vec``, or ``lobe_sign`` could silently flip
    which of the two manifold lobes is actually stepped along. Returns
    ``None`` on any extraction or propagation failure.
    """
    if branch not in ("unstable", "stable"):
        raise ValueError(f"branch must be 'unstable' or 'stable', got {branch!r}")
    extracted = _direction_at(
        torus,
        theta1_section,
        theta2,
        branch,
        ref_vec,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
    )
    if extracted is None:
        return None
    state0, vec = extracted
    perturbed = state0 + lobe_sign * eps * vec
    t0 = (theta1_section % (2.0 * np.pi)) / torus.omega1
    sign_dir = 1.0 if branch == "unstable" else -1.0
    if t_flow == 0.0:
        return perturbed
    arc = ccr4bp.propagate_ccr4bp(
        torus.system, perturbed, sign_dir * t_flow, t0=t0, rtol=rtol, atol=atol
    )
    return arc.state_f
