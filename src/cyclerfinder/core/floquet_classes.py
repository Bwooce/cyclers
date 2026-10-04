"""Planar 4x4 Floquet classification by the characteristic-polynomial invariants (#931).

For the planar block M (4x4) of a monodromy matrix of a Hamiltonian flow, the characteristic
polynomial is reciprocal::

    lambda^4 + alpha lambda^3 + beta lambda^2 + alpha lambda + 1 = 0
    alpha = -trace M,  beta = sum of the six principal 2x2 minors of M

and factors as (lambda^2 + b1 lambda + 1)(lambda^2 + b2 lambda + 1) with
b = (alpha +- sqrt(Delta)) / 2, Delta = alpha^2 - 4 (beta - 2) (Hadjidemetriou 1975b eqs. 39-47).
Broucke's (a1, a2) of NASA TR 32-1360 (eqs. 165-174) and AIAA J. 7:1003 (eqs. 34-45) are the same
(alpha, beta). His pair indices are k_B = lambda + 1/lambda = -b; this module reports the
normalised index ``k = -b / 2`` (cos of the rotation angle on the unit circle, critical at
|k| = 1, as in Hitzl & Henon 1977). Broucke's sign is -a1: k_B = [-a1 +- sqrt(Delta)] / 2 (TR
eq. 167, AIAA eq. 38; AIAA eq. 39 prints +a1, a sign slip that swaps regions 4/5 and 6/7).

A pair is stable iff |b| < 2 (|k| < 1). The regime is decided by the sign of Delta and the size and
sign of the b values, never by an eigenvalue modulus: at a Hamiltonian-Hopf transition the moduli
leave the unit circle only as the square root of the parameter (Olle, Pacha & Villanueva 2005), so
a modulus test is blind there, and a stable orbit has max |lambda| = 1 exactly.

Broucke's seven regions (Delta >= 0 unless stated; k_B = 2 k):

1. both pairs stable (|k_B| < 2)
2. Delta < 0, complex quartet
3. both pairs real-hyperbolic, k_B of opposite sign
4. both pairs real-hyperbolic, k_B1 > 2 and k_B2 > 2
5. both pairs real-hyperbolic, k_B1 < -2 and k_B2 < -2
6. one pair stable, the other hyperbolic with k_B > 2
7. one pair stable, the other hyperbolic with k_B < -2

``regime == "critical"`` marks Delta = 0 to a relative ``delta_tol``: the two pairs coincide. With
|k| < 1 that is the Hamiltonian-Hopf transition (Jorba & Olle 2004 map Ts at K = -1, L = 1/4:
alpha = -3, beta = 4 + L, Delta = 1 - 4L, multipliers exp(+-i arccos(3/4))); the linear data cannot
say whether it is the direct or the inverse case, and the note says so.

Orbits within ``boundary_tol`` of |k| = 1 are reported as ``regime == "boundary"`` with
``region is None``: the region is decided in the fourth decimal there (Broucke's 7P and 7A rows sit
on the line beta = -2 alpha - 2) and printed seven-digit states change region in 9 of 15 7P rows,
so classify re-corrected orbits only.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

#: indices of (x, y, x', y') inside a 6x6 state-transition matrix of a planar orbit.
PLANAR_INDICES = (0, 1, 3, 4)

DEFAULT_BOUNDARY_TOL = 1.0e-3  # on |k| = 1
DEFAULT_DELTA_TOL = 1.0e-9  # relative, on Delta = 0

HOPF_NOTE = "Hamiltonian-Hopf critical: direct or inverse undetermined from the linear data"


@dataclass(frozen=True)
class PlanarFloquetClass:
    alpha: float  # -trace; Broucke a1
    beta: float  # sum of principal 2x2 minors; Broucke a2
    delta: float  # alpha^2 - 4 (beta - 2)
    b1: complex  # (alpha + sqrt(Delta)) / 2; complex conjugate of b2 when Delta < 0
    b2: complex
    k1: complex  # -b1 / 2
    k2: complex  # -b2 / 2
    regime: str  # stable | saddle_centre | saddle_saddle | complex_quartet | boundary | critical
    region: int | None  # Broucke region 1..7; None on the boundary
    boundary: bool  # a pair within boundary_tol of |k| = 1
    reciprocity_error: float  # departure of the monodromy from a reciprocal polynomial
    note: str = ""  # for "critical": which collision, and that direct/inverse is undetermined

    @property
    def stable(self) -> bool:
        """Linearly stable (both pairs on the unit circle, away from the boundary)."""
        return self.regime == "stable"


def planar_block(monodromy: NDArray[np.float64]) -> NDArray[np.float64]:
    """The 4x4 planar block (x, y, x', y') of a 6x6 monodromy, or the matrix itself if 4x4."""
    m = np.asarray(monodromy, dtype=np.float64)
    if m.shape == (4, 4):
        return m
    if m.shape != (6, 6):
        raise ValueError(f"monodromy must be 4x4 or 6x6; got {m.shape}")
    return m[np.ix_(PLANAR_INDICES, PLANAR_INDICES)]


def planar_invariants(monodromy: NDArray[np.float64]) -> tuple[float, float, float]:
    """(alpha, beta, reciprocity_error) of a planar monodromy (4x4, or the block of a 6x6)."""
    coeffs = np.poly(planar_block(monodromy))  # [1, alpha, beta, gamma, delta]
    alpha, beta, gamma, det = (float(c) for c in coeffs[1:])
    rec = max(abs(gamma - alpha) / (1.0 + abs(alpha)), abs(det - 1.0))
    return alpha, beta, rec


def classify_invariants(
    alpha: float,
    beta: float,
    *,
    boundary_tol: float = DEFAULT_BOUNDARY_TOL,
    delta_tol: float = DEFAULT_DELTA_TOL,
    reciprocity_error: float = 0.0,
) -> PlanarFloquetClass:
    """Classify from (alpha, beta) of the reciprocal polynomial (see the module docstring)."""
    delta = alpha * alpha - 4.0 * (beta - 2.0)
    sq = math.sqrt(delta) if delta >= 0.0 else 1j * math.sqrt(-delta)
    b1 = complex((alpha + sq) / 2.0)
    b2 = complex((alpha - sq) / 2.0)
    k1, k2 = -b1 / 2.0, -b2 / 2.0
    if abs(delta) <= delta_tol * (1.0 + alpha * alpha + abs(beta)):
        # Delta = 0: the two pairs coincide (a double multiplier pair, a single Jordan block in
        # the generic case). On the unit circle this is the Hamiltonian-Hopf transition.
        note = HOPF_NOTE if abs(k1.real) < 1.0 else "double hyperbolic pair (saddle collision)"
        return PlanarFloquetClass(
            alpha, beta, delta, b1, b2, k1, k2, "critical", None, False, reciprocity_error, note
        )
    if delta < 0.0:
        return PlanarFloquetClass(
            alpha, beta, delta, b1, b2, k1, k2, "complex_quartet", 2, False, reciprocity_error
        )
    r1, r2 = k1.real, k2.real
    boundary = abs(abs(r1) - 1.0) <= boundary_tol or abs(abs(r2) - 1.0) <= boundary_tol
    if boundary:
        return PlanarFloquetClass(
            alpha, beta, delta, b1, b2, k1, k2, "boundary", None, True, reciprocity_error
        )
    stable1, stable2 = abs(r1) < 1.0, abs(r2) < 1.0
    if stable1 and stable2:
        regime, region = "stable", 1
    elif stable1 or stable2:
        hyper = r2 if stable1 else r1
        regime, region = "saddle_centre", (6 if hyper > 0.0 else 7)
    else:
        regime = "saddle_saddle"
        region = 3 if (r1 > 0.0) != (r2 > 0.0) else (4 if r1 > 0.0 else 5)
    return PlanarFloquetClass(
        alpha, beta, delta, b1, b2, k1, k2, regime, region, False, reciprocity_error
    )


def classify_planar_monodromy(
    monodromy: NDArray[np.float64],
    *,
    boundary_tol: float = DEFAULT_BOUNDARY_TOL,
    delta_tol: float = DEFAULT_DELTA_TOL,
) -> PlanarFloquetClass:
    """Classify a planar monodromy (4x4, or the planar block of a 6x6)."""
    alpha, beta, rec = planar_invariants(monodromy)
    return classify_invariants(
        alpha, beta, boundary_tol=boundary_tol, delta_tol=delta_tol, reciprocity_error=rec
    )


@dataclass(frozen=True)
class FamilyTransition:
    """A change of sign of Delta, of P(+1) (k = +1) or of P(-1) (k = -1) between two members."""

    kind: str  # "delta" | "k=+1" | "k=-1"
    index: int  # between members index and index + 1
    lo: float  # parameter of member index (the index itself when no params are given)
    hi: float
    located: float | None  # bisected parameter, when a callback was given
    region_before: int | None
    region_after: int | None


def _crossing_functions(alpha: float, beta: float) -> dict[str, float]:
    """Delta and the characteristic polynomial at +1 and -1: P(+-1) = (2 -+ ...) (2 +- b1)(2 +- b2).

    P(+1) = 2 + 2 alpha + beta = (2 + b1)(2 + b2) vanishes when a pair reaches k = +1;
    P(-1) = 2 - 2 alpha + beta = (2 - b1)(2 - b2) when a pair reaches k = -1. Both are single
    valued in (alpha, beta), so sign changes need no pairing of b1 and b2 between members.
    """
    return {
        "delta": alpha * alpha - 4.0 * (beta - 2.0),
        "k=+1": 2.0 + 2.0 * alpha + beta,
        "k=-1": 2.0 - 2.0 * alpha + beta,
    }


def find_family_transitions(
    monodromies: Sequence[NDArray[np.float64]],
    params: Sequence[float] | None = None,
    *,
    monodromy_at: Callable[[float], NDArray[np.float64]] | None = None,
    bisect_tol: float = 1.0e-9,
    max_bisect: int = 80,
) -> list[FamilyTransition]:
    """Detect sign changes of Delta, P(+1) and P(-1) along an ordered family.

    Stable windows of width of order mu are easy to step over (Hitzl & Henon 1977b), and a modulus
    test cannot see a Hamiltonian-Hopf transition (moduli move as sqrt of the parameter), so the
    scan uses the smooth invariants. A pair of crossings between two samples cancels in the sign
    test; refine the step when the family is sampled coarsely.

    With ``monodromy_at`` (parameter -> monodromy) and ``params``, each bracket is bisected to
    ``bisect_tol`` in the parameter and the located value is returned.
    """
    n = len(monodromies)
    if params is not None and len(params) != n:
        raise ValueError("params must match monodromies in length")
    if monodromy_at is not None and params is None:
        raise ValueError("monodromy_at needs params")
    xs = [float(p) for p in params] if params is not None else [float(i) for i in range(n)]
    invs = [planar_invariants(m) for m in monodromies]
    funcs = [_crossing_functions(a, b) for a, b, _ in invs]
    classes = [classify_invariants(a, b, boundary_tol=0.0, delta_tol=0.0) for a, b, _ in invs]
    out: list[FamilyTransition] = []
    for i in range(n - 1):
        for kind in ("delta", "k=+1", "k=-1"):
            g0, g1 = funcs[i][kind], funcs[i + 1][kind]
            if (g0 > 0.0) == (g1 > 0.0):  # an exact zero counts as non-positive
                continue
            located: float | None = None
            if monodromy_at is not None:
                lo, hi, glo = xs[i], xs[i + 1], g0
                for _ in range(max_bisect):
                    if abs(hi - lo) <= bisect_tol:
                        break
                    mid = 0.5 * (lo + hi)
                    a, b, _ = planar_invariants(monodromy_at(mid))
                    gm = _crossing_functions(a, b)[kind]
                    if (gm > 0.0) == (glo > 0.0):
                        lo, glo = mid, gm
                    else:
                        hi = mid
                located = 0.5 * (lo + hi)
            out.append(
                FamilyTransition(
                    kind, i, xs[i], xs[i + 1], located, classes[i].region, classes[i + 1].region
                )
            )
    return out
