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
