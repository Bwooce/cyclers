"""BCR4BP natural-parameter continuation along ``mu_sun`` (#303 / #292 Phase 2).

Continues a periodic orbit from CR3BP (``mu_sun = 0``) into the standard
incoherent BCR4BP at the Andreu / Rosales-Jorba parameter value
(``mu_sun ~ 328900.5``). The seed at ``mu_sun = 0`` is a sourced CR3BP family
member (e.g. an Earth-Moon planar L1 Lyapunov at a sourced Jacobi level); the
target is the published BCR4BP regime. The intermediate family members are
OUR computation -- they are NOT catalogue-promotable. Per the orbit-closure
discipline (`feedback_orbit_closure_discipline`) every accepted member
satisfies an independent Radau closure cross-check inside the corrector.

Wrapper around :func:`cyclerfinder.genome.bcr4bp_genome.correct_bcr4bp_periodic`
plus a simple natural-parameter predictor-corrector loop:

  * predictor: linear extrapolation of the free unknowns from the last two
    converged members (or zero-th-order from the last one on step 1);
  * corrector: ``correct_bcr4bp_periodic`` at the new ``mu_sun`` value, with
    ``T`` left free so the family member can drift toward its own Sun-
    commensurate value;
  * stability: integrate the 6x6 STM along the full period to get the
    monodromy; report Floquet eigenvalues + a coarse stability tag.

The stepping strategy is **logarithmic** in ``mu_sun + offset`` -- the
dynamical response to ``mu_sun`` is roughly multiplicative (the Sun's
effective tidal scale is ``~ mu_sun / a_sun^3``, which grows in lockstep with
``mu_sun`` at fixed ``a_sun``). Linear stepping from 0 to ~3e5 dumps almost
all members into the strong-perturbation regime; geometric / log stepping
keeps step-size proportional to family curvature.

Discipline
----------
  * READ-ONLY on ``bcr4bp.py`` and ``bcr4bp_genome.py``; this is a wrapper.
  * No catalogue writeback. No novelty claims.
  * Independent cross-check mandatory (inherited from the corrector).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    BCR4BPPeriodicOrbit,
    correct_bcr4bp_periodic,
)

__all__ = [
    "BCR4BPFamily",
    "BCR4BPFamilyMember",
    "continue_bcr4bp_family_in_musun",
]


# ---------------------------------------------------------------------------
# Dataclasses (frozen; the family is the immutable output of a continuation).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BCR4BPFamilyMember:
    """One converged member of a BCR4BP mu_sun-continuation family.

    Attributes
    ----------
    orbit :
        The full :class:`BCR4BPPeriodicOrbit` returned by the corrector. The
        ``system`` field inside carries the actual ``mu_sun`` value used.
    mu_sun_value :
        Convenience accessor: ``orbit.system.mu_sun`` echoed for indexing
        without unpacking the system.
    monodromy :
        6x6 STM at ``t = period``. ``None`` if STM integration failed (e.g.
        very stiff member where DOP853 trips).
    floquet :
        Complex eigenvalues of ``monodromy``. ``None`` if monodromy is None.
    stability_tag :
        Coarse human-readable tag derived from Floquet eigenvalues:

        * ``"monodromy_failed"`` -- could not compute monodromy.
        * ``"marginal"`` -- max ``|lambda|`` within 1e-3 of 1 (close to
          neutral; bookkeeping for the trivial unit eigenvalues).
        * ``"stable"`` -- max ``|lambda|`` < 1.001 (the trivial unit
          eigenvalues sit right at 1 in symplectic systems).
        * ``"unstable"`` -- max ``|lambda|`` > 1.001.
        * ``"hyperbolic_pair"`` -- a real reciprocal pair with magnitudes
          >> 1 (typical of L1-substitute orbits, which are saddle x center
          in their normal form).

        These labels are EVIDENCE, not catalogue-grade stability indices --
        the Phase 5 V-pipeline would compute a sourced numerical stability
        criterion against published tables.
    """

    orbit: BCR4BPPeriodicOrbit
    mu_sun_value: float
    monodromy: NDArray[np.float64] | None
    floquet: NDArray[np.complex128] | None
    stability_tag: str


@dataclass(frozen=True)
class BCR4BPFamily:
    """A converged BCR4BP mu_sun-continuation family.

    Attributes
    ----------
    members :
        Tuple of converged members in continuation order (low ``mu_sun`` -> high).
    mu_sun_extent :
        ``(min, max)`` of ``mu_sun_value`` across converged members. ``(0, 0)``
        for an empty family.
    seed_mu_sun :
        The ``mu_sun`` value of the input seed orbit (typically 0.0).
    walk_notes :
        Free-form provenance string: stepping strategy, n_steps_attempted,
        n_steps_converged, any structured failure modes encountered.
    """

    members: tuple[BCR4BPFamilyMember, ...]
    mu_sun_extent: tuple[float, float]
    seed_mu_sun: float
    walk_notes: str


# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------


def _mu_sun_schedule(
    mu_sun_seed: float,
    mu_sun_target: float,
    n_steps: int,
    method: Literal["geometric", "linear", "log"],
) -> list[float]:
    """Build the schedule of ``mu_sun`` values for the continuation.

    Returns the list of NEW ``mu_sun`` values to step into (does NOT include
    the seed value itself). Length == ``n_steps``.

    ``method="geometric"`` is identical to ``method="log"`` for monotonically
    increasing schedules: both use geometric spacing in ``(mu_sun + offset)``
    space with a small offset so the seed at zero is well-defined. The
    distinction is bookkeeping.
    """
    if n_steps <= 0:
        raise ValueError(f"n_steps must be positive; got {n_steps}")
    if mu_sun_seed == mu_sun_target:
        raise ValueError("mu_sun_seed == mu_sun_target; nothing to continue.")
    if method == "linear":
        return [
            mu_sun_seed + (mu_sun_target - mu_sun_seed) * (i + 1) / n_steps for i in range(n_steps)
        ]
    if method in ("geometric", "log"):
        # Offset so log(0) is defined; pick small relative to target.
        offset = max(1.0, 1e-6 * abs(mu_sun_target - mu_sun_seed))
        lo = math.log(mu_sun_seed + offset)
        hi = math.log(mu_sun_target + offset)
        return [math.exp(lo + (hi - lo) * (i + 1) / n_steps) - offset for i in range(n_steps)]
    raise ValueError(f"step_method must be 'linear', 'geometric', or 'log'; got {method!r}")


def _build_unknown_vec(
    state0: NDArray[np.float64], period: float, free_vars: tuple[int, ...]
) -> NDArray[np.float64]:
    """Extract the free-variable vector from a state + period."""
    vec = np.empty(len(free_vars), dtype=np.float64)
    for col, unknown in enumerate(free_vars):
        vec[col] = period if unknown == IDX_T else float(state0[unknown])
    return vec


def _apply_unknown_vec(
    state0: NDArray[np.float64],
    period: float,
    free_vars: tuple[int, ...],
    vec: NDArray[np.float64],
) -> tuple[NDArray[np.float64], float]:
    """Inverse of :func:`_build_unknown_vec`: write ``vec`` back to state/period."""
    new_state = state0.copy()
    new_period = period
    for col, unknown in enumerate(free_vars):
        if unknown == IDX_T:
            new_period = float(vec[col])
        else:
            new_state[unknown] = float(vec[col])
    return new_state, new_period


def _compute_monodromy(
    system: bcr4bp.BCR4BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64] | None:
    """Propagate the 6x6 STM over the full period; returns the monodromy.

    Returns ``None`` on integrator failure. Independent of the corrector --
    if the corrector succeeded but the monodromy integration fails (e.g. a
    very stiff member), the family member is still admitted, but its
    ``stability_tag`` becomes ``"monodromy_failed"``.
    """
    try:
        arc = bcr4bp.propagate_bcr4bp(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    except RuntimeError:
        return None
    return arc.stm


def _classify_floquet(eigs: NDArray[np.complex128]) -> str:
    """Return a coarse stability tag from the Floquet eigenvalues.

    See :class:`BCR4BPFamilyMember` for the tag taxonomy.
    """
    if eigs is None or len(eigs) == 0:
        return "monodromy_failed"
    mags = np.abs(eigs)
    max_mag = float(np.max(mags))
    # Look for a hyperbolic pair: at least one eigenvalue with |lambda| >> 1
    # AND a reciprocal with magnitude close to 1 / max_mag.
    has_hyper_pair = False
    if max_mag > 1.1:
        for i, lam_i in enumerate(eigs):
            mag_i = float(abs(lam_i))
            if mag_i < 1.01:
                continue
            # Look for a near-reciprocal partner.
            for j, lam_j in enumerate(eigs):
                if i == j:
                    continue
                mag_j = float(abs(lam_j))
                if mag_j > 0 and abs(mag_i * mag_j - 1.0) < 0.05:
                    has_hyper_pair = True
                    break
            if has_hyper_pair:
                break
    if has_hyper_pair:
        return "hyperbolic_pair"
    if max_mag > 1.001:
        return "unstable"
    if abs(max_mag - 1.0) < 1e-3:
        return "marginal"
    return "stable"


# ---------------------------------------------------------------------------
# Top-level continuation driver.
# ---------------------------------------------------------------------------


def continue_bcr4bp_family_in_musun(
    seed_orbit: BCR4BPPeriodicOrbit,
    seed_mu_sun: float = 0.0,
    target_mu_sun: float = 328900.5423094043,
    *,
    n_steps: int = 50,
    step_method: Literal["geometric", "linear", "log"] = "geometric",
    corrector_tol: float = 1e-10,
    closure_tol: float = 1e-6,
    free_vars: tuple[int, ...] = (IDX_X, IDX_YDOT, IDX_T),
    residual_indices: tuple[int, ...] = (IDX_Y, IDX_XDOT, IDX_ZDOT),
    is_half_period_residual: bool = True,
    monodromy: bool = True,
    sun_commensurate_n: int = 1,
    state_step_cap: float = 0.05,
    period_step_cap_frac: float = 0.1,
    require_monotone_decrease: bool = True,
    max_iter: int = 80,
    on_step: Callable[[int, BCR4BPFamilyMember], None] | None = None,
) -> BCR4BPFamily:
    """Continue a BCR4BP periodic orbit in the natural parameter ``mu_sun``.

    Parameters
    ----------
    seed_orbit :
        A converged BCR4BP periodic orbit (typically the CR3BP-limit case
        at ``mu_sun = 0`` -- equivalently a CR3BP periodic orbit returned
        by the BCR4BP corrector). The seed's ``system`` field carries the
        ``mu_sun`` value; we pass that explicitly via ``seed_mu_sun`` to
        avoid floating-point ambiguity.
    seed_mu_sun :
        The seed's mu_sun value (typically 0.0 for the CR3BP-limit anchor).
    target_mu_sun :
        The target mu_sun value. Default is the Rosales-Jorba 2023 Table 3
        Andreu value 328900.5423094043.
    n_steps :
        Number of CONTINUATION STEPS (the family will have at most
        ``n_steps + 1`` members including the seed -- ``n_steps`` new
        members are appended).
    step_method :
        ``"geometric"`` / ``"log"`` (equivalent) or ``"linear"``. The Sun
        perturbation strength is multiplicative in ``mu_sun``, so geometric
        stepping is the natural choice.
    corrector_tol, closure_tol :
        Newton convergence tolerance and independent (Radau) closure tolerance
        passed down to :func:`correct_bcr4bp_periodic`.
    free_vars, residual_indices, is_half_period_residual :
        Newton mask. Default is the symmetric / perpendicular-crossing mask
        ``(x, vy, T)`` free with ``(y, vx, vz)`` residual at ``T/2`` --
        matches the Phase 1 corrector test and is the natural mask for a
        planar libration-substitute family.
    monodromy :
        If True (default), integrate the 6x6 STM over the full period at
        each converged member to get the monodromy. Drives the
        stability_tag classification.
    sun_commensurate_n :
        Commensurability bookkeeping. The seed's period is recorded as
        ``n=sun_commensurate_n`` Sun-synodic periods; subsequent members
        inherit the same n. If ``T`` is in ``free_vars`` the corrector
        lets the period drift -- ``sun_phase_drift`` then becomes nonzero
        and is reported per-member.
    state_step_cap, period_step_cap_frac :
        Per-component caps inside the corrector's Newton step. Tighter
        than the corrector default since the seed-to-member step is small
        in a properly-tuned continuation.
    require_monotone_decrease :
        Passed to the corrector. True is safer for routine steps; the
        caller may set False for known-stiff segments.
    max_iter :
        Newton iteration cap per step.
    on_step :
        Optional callback ``on_step(step_idx, member)`` invoked after each
        converged member is appended. Use for progress reporting.

    Returns
    -------
    BCR4BPFamily
        The continuation family. ``walk_notes`` records the stepping
        strategy and the converged-vs-attempted count.

    Notes
    -----
    On a divergent step the continuation BREAKS (does not silently skip
    ahead). The returned family then has fewer members than ``n_steps``
    and ``walk_notes`` documents where it broke. The seed is NOT included
    in ``members`` -- it's already in the caller's hand -- but its
    mu_sun and identity are echoed in ``seed_mu_sun`` and ``walk_notes``.
    """
    if seed_orbit.system.mu_sun != seed_mu_sun:
        raise ValueError(
            "seed_orbit.system.mu_sun != seed_mu_sun "
            f"({seed_orbit.system.mu_sun} vs {seed_mu_sun}) -- "
            "the seed orbit must be converged at the declared seed_mu_sun."
        )
    if not seed_orbit.converged:
        raise ValueError("seed_orbit.converged is False -- refusing to continue.")

    free_vars = tuple(sorted(set(int(v) for v in free_vars)))
    residual_indices = tuple(sorted(set(int(v) for v in residual_indices)))

    schedule = _mu_sun_schedule(seed_mu_sun, target_mu_sun, n_steps, step_method)
    # History of (mu_sun, unknown_vec) for the linear predictor.
    history_mu: list[float] = [seed_mu_sun]
    history_vec: list[NDArray[np.float64]] = [
        _build_unknown_vec(seed_orbit.state_initial, seed_orbit.period_nondim, free_vars)
    ]
    # Carry the FULL state (including any fixed components) so the corrector
    # has a consistent IC; the predictor only moves the free components.
    last_full_state = seed_orbit.state_initial.copy()
    last_period = seed_orbit.period_nondim

    members: list[BCR4BPFamilyMember] = []
    n_attempted = 0
    break_reason: str | None = None

    for step_idx, mu_sun_new in enumerate(schedule):
        n_attempted += 1
        # Predictor: linear in mu_sun.
        if len(history_mu) >= 2:
            mu_prev2, mu_prev1 = history_mu[-2], history_mu[-1]
            vec_prev2, vec_prev1 = history_vec[-2], history_vec[-1]
            if mu_prev1 != mu_prev2:
                slope = (vec_prev1 - vec_prev2) / (mu_prev1 - mu_prev2)
                pred_vec = vec_prev1 + slope * (mu_sun_new - mu_prev1)
            else:
                pred_vec = vec_prev1.copy()
        else:
            pred_vec = history_vec[-1].copy()

        guess_state, guess_period = _apply_unknown_vec(
            last_full_state, last_period, free_vars, pred_vec
        )

        sys_new = bcr4bp.BCR4BPSystem(
            mu=seed_orbit.system.mu,
            mu_sun=float(mu_sun_new),
            a_sun_nondim=seed_orbit.system.a_sun_nondim,
            omega_sun_nondim=seed_orbit.system.omega_sun_nondim,
            theta_sun0=seed_orbit.system.theta_sun0,
        )

        try:
            corrected = correct_bcr4bp_periodic(
                sys_new,
                guess_state,
                guess_period,
                sun_commensurate_n=sun_commensurate_n,
                free_vars=free_vars,
                residual_indices=residual_indices,
                is_half_period_residual=is_half_period_residual,
                tol=corrector_tol,
                independent_tol=closure_tol,
                state_step_cap=state_step_cap,
                period_step_cap_frac=period_step_cap_frac,
                max_iter=max_iter,
                require_monotone_decrease=require_monotone_decrease,
                notes=f"#303 mu_sun-continuation step {step_idx + 1}/{n_steps}",
            )
        except (RuntimeError, ValueError) as exc:
            break_reason = (
                f"corrector raised at step {step_idx + 1} (mu_sun={mu_sun_new:.6e}): "
                f"{type(exc).__name__}: {exc}"
            )
            break

        if not corrected.converged:
            break_reason = (
                f"corrector did NOT converge at step {step_idx + 1} "
                f"(mu_sun={mu_sun_new:.6e}); "
                f"corrector_residual={corrected.corrector_residual:.3e}, "
                f"independent_closure={corrected.independent_closure_residual:.3e}"
            )
            break

        # Monodromy + stability tag.
        mono: NDArray[np.float64] | None = None
        eigs: NDArray[np.complex128] | None = None
        tag = "monodromy_skipped"
        if monodromy:
            mono = _compute_monodromy(sys_new, corrected.state_initial, corrected.period_nondim)
            if mono is not None:
                try:
                    eigs = np.linalg.eigvals(mono).astype(np.complex128)
                    tag = _classify_floquet(eigs)
                except np.linalg.LinAlgError:
                    eigs = None
                    tag = "monodromy_failed"
            else:
                tag = "monodromy_failed"

        member = BCR4BPFamilyMember(
            orbit=corrected,
            mu_sun_value=float(mu_sun_new),
            monodromy=mono,
            floquet=eigs,
            stability_tag=tag,
        )
        members.append(member)
        if on_step is not None:
            on_step(step_idx, member)

        history_mu.append(float(mu_sun_new))
        history_vec.append(
            _build_unknown_vec(corrected.state_initial, corrected.period_nondim, free_vars)
        )
        last_full_state = corrected.state_initial.copy()
        last_period = corrected.period_nondim

    if members:
        mu_sun_extent = (
            float(min(m.mu_sun_value for m in members)),
            float(max(m.mu_sun_value for m in members)),
        )
    else:
        mu_sun_extent = (float(seed_mu_sun), float(seed_mu_sun))

    walk_notes = (
        f"step_method={step_method}; n_steps_attempted={n_attempted}; "
        f"n_steps_converged={len(members)}; seed_mu_sun={seed_mu_sun}; "
        f"target_mu_sun={target_mu_sun}"
    )
    if break_reason is not None:
        walk_notes += f"; break_reason: {break_reason}"

    return BCR4BPFamily(
        members=tuple(members),
        mu_sun_extent=mu_sun_extent,
        seed_mu_sun=float(seed_mu_sun),
        walk_notes=walk_notes,
    )
