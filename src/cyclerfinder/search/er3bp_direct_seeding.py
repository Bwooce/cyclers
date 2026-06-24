"""#436 Task 1: CR3BP-independent direct-e>0 seed grid.

Hunts for a novel ER3BP "e>0-only" periodic family (one with NO circular
3-body limit). The seeds here are deliberately INDEPENDENT of any CR3BP
family: a blind symmetric x-axis-crossing initial-condition grid placed
directly at the target eccentricity.

A seed extrapolated from a CR3BP orbit would, by construction, continue back
to that CR3BP ancestor and so could never be an e>0-only discovery. Building
the grid directly in (x0, ydot0) at the target ``system.e`` avoids that
ancestry entirely.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_continuation import continue_er3bp_family_in_e_partial
from cyclerfinder.genome.er3bp_periodic import (
    ER3BPPeriodicOrbit,
    correct_er3bp_periodic,
)


@dataclass(frozen=True)
class DirectEr3bpSeed:
    """A CR3BP-independent symmetric IC seed for the ER3BP corrector."""

    label: str
    system: ER3BPSystem  # target ER3BP system; system.e is the target eccentricity
    state0: NDArray[np.float64]  # symmetric x-axis-crossing IC, shape (6,)
    period_f: float  # true-anomaly span guess for the (half-)period
    is_half_period_residual: bool
    target_e: float
    source: str  # provenance string (must NOT reference any CR3BP family)


def direct_e_seed_grid(
    system: ER3BPSystem,
    x0_range: tuple[float, float],
    ydot0_range: tuple[float, float],
    n_x: int,
    n_ydot: int,
    period_f: float,
    is_half_period_residual: bool = True,
) -> list[DirectEr3bpSeed]:
    """Build a blind symmetric-IC grid directly at ``system.e``.

    Each seed is a symmetric x-axis crossing ``state0 = [x0, 0, 0, 0, ydot0, 0]``.
    ``x0`` spans ``x0_range`` over ``n_x`` linspace points and ``ydot0`` spans
    ``ydot0_range`` over ``n_ydot`` points; the returned list has ``n_x * n_ydot``
    seeds ordered with ``x0`` as the outer (slowest-varying) index.

    The grid is CR3BP-independent by construction: ICs are placed directly at
    the target eccentricity rather than extrapolated from a circular-problem
    family, so the ``source`` string names only the grid.
    """
    source = (
        f"#436 CR3BP-independent symmetric-IC grid "
        f"x0∈{x0_range} ydot0∈{ydot0_range} ({n_x}x{n_ydot})"
    )
    x0_values = np.linspace(x0_range[0], x0_range[1], n_x)
    ydot0_values = np.linspace(ydot0_range[0], ydot0_range[1], n_ydot)

    seeds: list[DirectEr3bpSeed] = []
    for i, x0 in enumerate(x0_values):
        for j, ydot0 in enumerate(ydot0_values):
            state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0], dtype=np.float64)
            seeds.append(
                DirectEr3bpSeed(
                    label=f"#436-direct-grid-x{i}-yd{j}",
                    system=system,
                    state0=state0,
                    period_f=period_f,
                    is_half_period_residual=is_half_period_residual,
                    target_e=system.e,
                    source=source,
                )
            )
    return seeds


def converge_direct_seed(
    seed: DirectEr3bpSeed,
    *,
    tol: float = 1e-10,
) -> ER3BPPeriodicOrbit | None:
    """Forward-converge a direct seed at its target eccentricity.

    Runs the strict-periodicity corrector on ``seed.state0`` in ``seed.system``
    (whose ``e`` is the target eccentricity). The corrector's ``period_f`` is
    the integration span: when ``is_half_period_residual`` is set this is the
    half-period, but ``seed.period_f`` stores the FULL period, so it is halved
    here (matching ``er3bp_discovery.continue_and_monitor``).

    Returns the converged orbit if both the orbit's own ``converged`` state
    (no exception raised) and ``corrector_residual < tol`` hold, else ``None``.
    """
    integration_f = seed.period_f / 2.0 if seed.is_half_period_residual else seed.period_f
    try:
        orbit = correct_er3bp_periodic(
            system=seed.system,
            state_guess=seed.state0,
            period_f=integration_f,
            is_half_period_residual=seed.is_half_period_residual,
            tol=tol,
        )
    except Exception:
        return None
    if orbit.corrector_residual < tol:
        return orbit
    return None


def classify_no_cr3bp_limit(
    converged_orbit: ER3BPPeriodicOrbit,
    system: ER3BPSystem,
    *,
    n_steps: int = 30,
    death_floor: float = 1e-3,
) -> dict[str, object]:
    """Reverse-continue a converged e>0 orbit toward e=0 to test CR3BP ancestry.

    Continues the family from ``system.e`` down to ``e_target=0.0``. The
    corrector's integration span is ``converged_orbit.period_f / 2`` for a
    half-period residual (the orbit stores the FULL period), matching the
    convention in :func:`converge_direct_seed`.

    Returns one of:
      * ``{"status": "cr3bp_continuous", "min_e": <final e>}`` if the
        continuation reaches ~0 (a CR3BP ancestor exists -> negative result);
      * ``{"status": "e_only_candidate", "death_e": <e>}`` if it dies at
        ``death_e > death_floor`` (no CR3BP limit -> candidate);
      * ``{"status": "inconclusive", "death_e": <e>}`` if it dies at
        ``death_e <= death_floor`` (numerical noise near e=0).
    """
    # The continuation corrector uses a half-period symmetric residual, so its
    # integration span is the half-period. ``converged_orbit.period_f`` is the
    # FULL period (the corrector returns ``period_f * 2`` when given a
    # half-period span), so halve it.
    integration_f = converged_orbit.period_f / 2.0

    orbits, death_e = continue_er3bp_family_in_e_partial(
        system,
        converged_orbit.state0,
        integration_f,
        0.0,
        n_steps,
        is_half_period_residual=True,
    )

    if death_e is None:
        final_e = orbits[-1].e if orbits else 0.0
        return {"status": "cr3bp_continuous", "min_e": final_e}

    if abs(death_e) > death_floor:
        return {"status": "e_only_candidate", "death_e": death_e}

    return {"status": "inconclusive", "death_e": death_e}
