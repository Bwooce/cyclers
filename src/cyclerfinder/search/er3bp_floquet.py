"""#432 Floquet monitor for the ER3BP discovery campaign.

The monodromy is the full-period (f in [0, period_f]) STM from propagate_er3bp;
its eigenvalues classify stability and flag bifurcations (eigenvalue on the unit
circle). Conventions mirror the #347 Floquet framework.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp

_UNIT_CIRCLE_TOL = 1.0e-3  # |λ| within this of 1.0 counts as "on the unit circle"


def er3bp_monodromy(
    state0: NDArray[np.float64], period_f: float, system: ER3BPSystem
) -> NDArray[np.float64]:
    """Full-period monodromy (6x6 STM over f in [0, period_f]) via propagate_er3bp."""
    _f, _hist, stm = propagate_er3bp(
        np.asarray(state0, dtype=np.float64),
        (0.0, period_f),
        system,
        with_stm=True,
    )
    return np.asarray(stm, dtype=np.float64)


@dataclass(frozen=True)
class FloquetResult:
    eigenvalues: tuple[complex, ...]
    stability_tag: str  # "stable" | "unstable" | "marginal"
    on_unit_circle: bool  # a non-trivial eigenvalue sits on the unit circle


def floquet_classify(
    monodromy: NDArray[np.float64], *, unit_circle_tol: float = _UNIT_CIRCLE_TOL
) -> FloquetResult:
    """Classify stability from the monodromy eigenvalues.

    stable: all |λ| <= 1 + tol. unstable: some |λ| > 1 + tol. marginal: max |λ|
    is within tol of 1. on_unit_circle flags a non-trivial eigenvalue (not the
    trivial λ=1 pair) sitting on the unit circle — the bifurcation signal.
    """
    eig = np.linalg.eigvals(np.asarray(monodromy, dtype=np.float64))
    mags = np.abs(eig)
    max_mag = float(mags.max())
    if max_mag > 1.0 + unit_circle_tol:
        tag = "unstable"
    elif max_mag < 1.0 - unit_circle_tol:
        tag = "stable"
    else:
        tag = "marginal"
    on_uc = bool(
        np.any((np.abs(mags - 1.0) <= unit_circle_tol) & (np.abs(eig - 1.0) > unit_circle_tol))
    )
    return FloquetResult(
        eigenvalues=tuple(complex(v) for v in eig),
        stability_tag=tag,
        on_unit_circle=on_uc,
    )
