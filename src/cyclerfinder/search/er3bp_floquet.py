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
from cyclerfinder.core.floquet_classes import PLANAR_INDICES, classify_planar_monodromy

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


def _is_planar(m: NDArray[np.float64], tol: float = 1.0e-9) -> bool:
    """True if the z, z' block of a 6x6 monodromy is decoupled from the planar block."""
    if m.shape != (6, 6):
        return False
    zi = (2, 5)
    cross = np.abs(m[np.ix_(zi, PLANAR_INDICES)]).max(), np.abs(m[np.ix_(PLANAR_INDICES, zi)]).max()
    return bool(max(cross) <= tol)


@dataclass(frozen=True)
class FloquetResult:
    eigenvalues: tuple[complex, ...]
    stability_tag: str  # "stable" | "unstable" | "marginal"
    on_unit_circle: bool  # a non-trivial eigenvalue sits on the unit circle


def floquet_classify(
    monodromy: NDArray[np.float64], *, unit_circle_tol: float = _UNIT_CIRCLE_TOL
) -> FloquetResult:
    """Classify linear stability of an ER3BP periodic orbit from its 6x6 monodromy.

    Tags: "stable" (every multiplier on the unit circle, away from a critical point),
    "unstable" (a real-hyperbolic pair or a complex quartet), "marginal" (at a critical point:
    a pair at k = +-1 or the Delta = 0 Hamiltonian-Hopf collision). A symplectic monodromy of a
    linearly stable orbit has max |lambda| = 1 exactly, so the modulus alone can never say
    "stable" (it was the bug, #931), and near a Hopf transition the moduli leave the circle only as
    the square root of the parameter.

    A planar orbit (the z block decoupled) is decided by the 4x4 reciprocal-polynomial invariants
    of ``core.floquet_classes`` (sign of Delta, size of the pair indices k). A spatial orbit has no
    such closed test here: it is "unstable" if any |lambda| > 1 + unit_circle_tol, "marginal" if a
    multiplier is within unit_circle_tol of +1 or -1, else "stable" (a complex quartet whose
    modulus excess is below the tolerance would read stable; use a planar orbit or a smaller
    tolerance to resolve it). on_unit_circle flags a multiplier (other than 1) on the circle.
    """
    m = np.asarray(monodromy, dtype=np.float64)
    eig = np.linalg.eigvals(m)
    mags = np.abs(eig)
    if _is_planar(m):
        regime = classify_planar_monodromy(m).regime
        tag = {"stable": "stable", "boundary": "marginal", "critical": "marginal"}.get(
            regime, "unstable"
        )
    elif float(mags.max()) > 1.0 + unit_circle_tol:
        tag = "unstable"
    elif np.any(np.abs(eig - 1.0) <= unit_circle_tol) or np.any(
        np.abs(eig + 1.0) <= unit_circle_tol
    ):
        tag = "marginal"
    else:
        tag = "stable"
    on_uc = bool(
        np.any((np.abs(mags - 1.0) <= unit_circle_tol) & (np.abs(eig - 1.0) > unit_circle_tol))
    )
    return FloquetResult(
        eigenvalues=tuple(complex(v) for v in eig),
        stability_tag=tag,
        on_unit_circle=on_uc,
    )
