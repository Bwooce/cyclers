"""#582 (stage 3a of #581): asymmetric/spatial-isolated 3D CR3BP GA fitness.

See ``data/OUTSTANDING.md`` **#582** for the full spec. Short version: the
existing isolated-family seed base (``search/er3bp_isolated_seeds.py``, #440)
generates only SYMMETRIC (perpendicular x-z-plane-crossing) circular seeds at
mu=0.001 for the five tabulated interior mean-motion resonances (MMRs) --
structurally it cannot generate asymmetric or genuinely-spatial seeds, because
it uses ``correct_symmetric_fixed_jacobi``. This module supplies the fitness
function half of a niching-GA layer (``search/niching_ga.py``, #581, unmodified
and reused exactly) that searches the FULL asymmetric 3D initial-condition
space directly, so a subsequent run of the existing asymmetric corrector
(``search/cr3bp_general_periodic_3d.py::correct_general_periodic_3d``, #291)
can converge candidates the symmetric seed base can never reach.

**System: mu=0.001.** A generic Sun-planet-mass system, matching #440's own
convention and the corrected Antoniadou & Libert 2019 anchor in
``genome/known_corpus_3d.py`` (mu=0.001; NOT Earth-Moon mu=0.0121506683 --
see #582's Fable correction in OUTSTANDING.md).

**Scope: 5 INTERIOR MMRs only** (``er3bp_isolated_seeds.MMR_SEMI_MAJOR_AXES``).
Exterior MMR bands (e.g. the exterior 1:2) are explicitly OUT OF SCOPE -- #440's
own docstring flags the exterior 1:2 as a known family-selection trap, and no
exterior semi-major axis is tabulated anywhere in this codebase.

Genome layout ``(x0, z0, xdot0, ydot0, zdot0, T)`` -- ``y0`` is pinned to 0
(the standard x-z-plane-crossing convention every corrector in this codebase
already uses for periapsis-referenced ICs). ``T`` is a genome variable, bounded
to roughly +-50% of the resonance's own linearized period
``T0 = 2*pi / (a1**-1.5 - 1)`` (#440's own formula).

Fitness (MAXIMIZED, bounded in (0, 1]) combines a periodicity-defect term
``||X(T) - X(0)||^2`` (full 6D state, since an asymmetric candidate need not
return through the same perpendicular crossing) with a SOFT target-Jacobi-band
penalty, in the same bounded-reciprocal form Gurfil & Kasdin (2002) Eq. 15
uses for exactly the same reason (a smooth, positive, differentiable basin
around the target rather than a hard constraint the GA could get stuck
against). The Jacobi term is SOFT by design -- #440's own note is that a
converged member sits at a "nearby, finite-amplitude-shifted" state relative
to the circular seed, and an asymmetric branch may shift the Jacobi constant
further still; this fitness must not exclude that.

Death-penalty guards (return the GA's floor fitness, 0.0):

  * ``T <= T0/2`` -- #582's own minimum-period floor, excluding degenerate
    near-equilibrium loops (mirrors Gurfil-Kasdin's own period sanity floor).
  * A primary or secondary collision at ANY point during the propagated
    arc (not just the endpoint) -- mirrors Gurfil-Kasdin's own Eq. 17-style
    collision-exclusion constraint. The exclusion radii are UNSOURCED generic
    safety margins (see :data:`DEFAULT_PRIMARY_EXCLUSION_RADIUS` /
    :data:`DEFAULT_SECONDARY_EXCLUSION_RADIUS` docstrings) -- mu=0.001 is a
    generic system with no defined physical body radii, unlike the
    Earth-anchored ``core/er3bp_geocentric.py`` collision constraint (which
    uses the real Earth radius).
  * A non-finite defect, Jacobi constant, or integrator failure.

A GA fitness peak is a BASIN INDICATOR, not a converged orbit
(per [[feedback_orbit_closure_discipline]] and #581's own counterweight) --
every survivor MUST be refined through the existing asymmetric corrector
before any claim; see ``search/isolated_3d_asymmetric_pipeline.py``.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp

# Genome index aliases: (x0, z0, xdot0, ydot0, zdot0, T). y0 is pinned to 0.
IDX_X0 = 0
IDX_Z0 = 1
IDX_XDOT0 = 2
IDX_YDOT0 = 3
IDX_ZDOT0 = 4
IDX_T = 5
GENOME_LEN = 6

#: Generic, UNSOURCED singularity-exclusion radius about the primary
#: (nondim, l_km=1 units per #440's ``CR3BPSystem(mu=0.001, ..., l_km=1.0,
#: t_s=1.0)`` convention). mu=0.001 is a generic planetary system with no
#: catalogued physical body radius -- unlike ``core/er3bp_geocentric.py``'s
#: collision constraint, which uses the real Earth radius (Eq. 17 of
#: Gurfil-Kasdin 2002). This value is a conservative engineering safety
#: margin only (order 1% of the smallest tabulated interior-MMR semi-major
#: axis, a1=0.3419 for 5:1): large enough to guard the EOM singularity,
#: small enough not to exclude any of the five tabulated resonant seeds
#: (whose crossing radius from the primary is a1, always >> this radius).
DEFAULT_PRIMARY_EXCLUSION_RADIUS = 0.01

#: Same convention as :data:`DEFAULT_PRIMARY_EXCLUSION_RADIUS`, about the
#: secondary at ``(1-mu, 0, 0)``.
DEFAULT_SECONDARY_EXCLUSION_RADIUS = 0.01


def mmr_t0(a1: float) -> float:
    """#440's own linearized resonant period formula: ``2*pi / (a1**-1.5 - 1)``."""
    if a1 <= 0.0:
        raise ValueError(f"a1 must be > 0; got {a1}")
    return float(2.0 * math.pi / (a1**-1.5 - 1.0))


def mmr_bounds(
    a1: float,
    *,
    mu: float = 0.001,
    x0_frac: float = 0.15,
    ydot0_frac: float = 0.35,
    z0_abs: float = 0.05,
    xdot0_abs: float = 0.05,
    zdot0_abs: float = 0.05,
    t_frac: float = 0.5,
) -> tuple[list[tuple[float, float]], float, float, float]:
    """GA bounds box for one interior MMR, centered on the analytic e=0 guess.

    Returns ``(bounds, x0_guess, ydot0_guess, t0)`` where ``bounds`` is in
    genome order ``(x0, z0, xdot0, ydot0, zdot0, T)``.

    ``t_frac`` implements #582's own required period bound: +-50% of
    ``T0 = mmr_t0(a1)`` (default ``t_frac=0.5``).

    The state-component half-widths (``x0_frac``, ``ydot0_frac``, the
    ``*_abs`` z/xdot/zdot widths) are NOT sourced from any paper -- there is
    no published search box for the asymmetric family this fitness function
    targets (that is the literature-open gap #582 exists to probe). They are
    a documented engineering choice, calibrated empirically against ALL FIVE
    known #440 circular members (``all_mmr_seeds()``): 4/5 land within ~1% of
    the analytic guess on every component, but the 3:2 member (#440's own
    "most eccentric" case) lands ``ydot0`` ~20% off the linearized guess --
    hence ``ydot0_frac`` defaults wide (0.35, ~1.75x headroom over the worst
    observed case) while ``x0_frac`` (worst observed -4.1%) can stay tighter.
    Widths are still narrow enough that a GA population explores THIS
    resonance's basin rather than drifting into a neighboring MMR or the
    exterior-1:2 family-selection trap #440 documents.
    """
    if a1 <= 0.0:
        raise ValueError(f"a1 must be > 0; got {a1}")
    x0_guess = -mu + a1
    ydot0_guess = math.sqrt((1.0 - mu) / a1) - a1
    t0 = mmr_t0(a1)
    bounds = [
        (x0_guess * (1.0 - x0_frac), x0_guess * (1.0 + x0_frac)),
        (-z0_abs, z0_abs),
        (-xdot0_abs, xdot0_abs),
        (ydot0_guess * (1.0 - ydot0_frac), ydot0_guess * (1.0 + ydot0_frac)),
        (-zdot0_abs, zdot0_abs),
        (t0 * (1.0 - t_frac), t0 * (1.0 + t_frac)),
    ]
    return bounds, x0_guess, ydot0_guess, t0


@dataclass(frozen=True)
class IsolatedAsymmetricFitnessConfig:
    """Constants for :func:`isolated_3d_asymmetric_fitness`."""

    mu: float
    t0: float
    """Linearized resonant period (#440 formula); backs the ``T > T0/2`` guard."""
    jacobi_target: float
    """Soft target Jacobi constant (typically the circular seed's own C)."""
    jacobi_soft_width: float = 0.05
    """Scale of the soft Jacobi penalty (nondim Jacobi units)."""
    w_period: float = 1.0
    w_jacobi: float = 1.0
    primary_radius: float = DEFAULT_PRIMARY_EXCLUSION_RADIUS
    secondary_radius: float = DEFAULT_SECONDARY_EXCLUSION_RADIUS
    rtol: float = 1e-11
    atol: float = 1e-11
    n_samples: int = 400
    """Collision-event sampling density: caps ``max_step`` to ``T/n_samples``."""


def genome_to_state0(genome: NDArray[np.float64]) -> NDArray[np.float64]:
    """``(x0, z0, xdot0, ydot0, zdot0, T)`` -> full 6D state ``(x,0,z,xdot,ydot,zdot)``."""
    g = np.asarray(genome, dtype=np.float64)
    return np.array(
        [g[IDX_X0], 0.0, g[IDX_Z0], g[IDX_XDOT0], g[IDX_YDOT0], g[IDX_ZDOT0]],
        dtype=np.float64,
    )


def _collision_event_factory(r1_excl: float, r2_excl: float) -> Callable[..., float]:
    """Terminal scipy event: fires when either exclusion radius is breached.

    3-arg signature ``(t, y, mu)`` matches ``solve_ivp``'s convention of
    forwarding the same ``args=`` tuple to every event callable (mirrors
    ``core.bcr4bp.propagate_bcr4bp``'s ``_collision_event``).
    """

    def _event(t: float, y: NDArray[np.float64], mu: float) -> float:
        r1_sq = (y[0] + mu) ** 2 + y[1] ** 2 + y[2] ** 2
        r2_sq = (y[0] - 1.0 + mu) ** 2 + y[1] ** 2 + y[2] ** 2
        if r1_sq < r1_excl * r1_excl or r2_sq < r2_excl * r2_excl:
            return 0.0
        return 1.0

    _event.terminal = True  # type: ignore[attr-defined]
    return _event


def propagate_defect(
    state0: NDArray[np.float64],
    period: float,
    config: IsolatedAsymmetricFitnessConfig,
) -> tuple[float, bool]:
    """Propagate one period; return ``(||X(T) - X(0)||^2, collided_or_failed)``.

    Collision-exclusion is enforced THROUGHOUT the propagation via a terminal
    scipy event, not merely checked at the endpoint -- an asymmetric trial
    orbit can graze a primary mid-arc and still land far from it at ``t=T``.
    """
    if period <= 0.0 or not math.isfinite(period):
        return float("inf"), True
    event = _collision_event_factory(config.primary_radius, config.secondary_radius)
    try:
        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, period),
            np.asarray(state0, dtype=np.float64),
            args=(config.mu,),
            method="DOP853",
            rtol=config.rtol,
            atol=config.atol,
            max_step=max(period / config.n_samples, 1e-6),
            events=event,
        )
    except Exception:
        return float("inf"), True
    if not sol.success or sol.status == 1:  # integration failure OR collision event fired
        return float("inf"), True
    diff = sol.y[:, -1] - np.asarray(state0, dtype=np.float64)
    if not np.all(np.isfinite(diff)):
        return float("inf"), True
    return float(np.dot(diff, diff)), False


def isolated_3d_asymmetric_fitness(
    genome: NDArray[np.float64],
    *,
    config: IsolatedAsymmetricFitnessConfig,
) -> float:
    """GA objective (MAXIMIZE): plain vector-in/scalar-out, ``niching_ga`` contract.

    See the module docstring for the fitness formula and death-penalty guards.
    Bounded in ``(0, 1]``; the death-penalty floor is exactly ``0.0``.
    """
    g = np.asarray(genome, dtype=np.float64)
    if g.shape != (GENOME_LEN,):
        raise ValueError(f"genome must have shape ({GENOME_LEN},); got {g.shape}")
    t = float(g[IDX_T])
    if not math.isfinite(t) or t <= config.t0 / 2.0:
        return 0.0

    state0 = genome_to_state0(g)
    defect, collided = propagate_defect(state0, t, config)
    if collided or not math.isfinite(defect):
        return 0.0

    try:
        jacobi = cr3bp.jacobi_constant(state0, config.mu)
    except (ValueError, ZeroDivisionError):
        return 0.0
    if not math.isfinite(jacobi):
        return 0.0

    jacobi_dev = (jacobi - config.jacobi_target) / config.jacobi_soft_width
    penalty = config.w_period * defect + config.w_jacobi * jacobi_dev * jacobi_dev
    if not math.isfinite(penalty) or penalty < 0.0:
        return 0.0
    return float(1.0 / (1.0 + penalty))
