"""Passing and testing a turning point (fold) of a one-parameter shooting problem (#1055 sec. 12.5).

Generalised from the #1025 / #1034 monitor-variable continuation and the #1004 two-solution check,
for any formulation exposing ``rj(z, s) -> (r, J)`` at parameter s (the moon-mass factor sigma).
No SVD tangent is used: sigma becomes an unknown and the state component that changed most over the
last secant (column-scaled) is stepped and held fixed by an extra row (the #1034 tangent rule).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Arr = NDArray[np.float64]
SigmaResJac = Callable[[Arr, float], tuple[Arr, Arr]]


def sigma_derivative(rj: SigmaResJac, z: Arr, s: float, rel: float = 1e-6) -> Arr:
    """Central difference of the residual in s at fixed z."""
    rp, _ = rj(z, s * (1.0 + rel))
    rm, _ = rj(z, s * (1.0 - rel))
    return np.asarray((rp - rm) / (2.0 * rel * s), dtype=np.float64)


@dataclass
class MonitorStep:
    z: Arr
    s: float
    ok: bool
    ratio: float
    iterations: int
    component: int


def monitor_step(
    rj: SigmaResJac,
    za: Arr,
    sa: float,
    zb: Arr,
    sb: float,
    ok: Callable[[Arr, float], bool],
    ratio: float = 1.0,
    halvings: int = 7,
    iters: int = 15,
) -> MonitorStep:
    """One monitor-variable step past (zb, sb) along the secant from (za, sa).

    The component j with the largest column-scaled change over the secant is advanced to
    target = zb[j] + ratio (zb[j] - za[j]) and held by an extra row; (z, s) solve the residual plus
    that row by damped Newton (column-scaled least squares, backtracking). ``ok(z, s)`` is the
    acceptance test of a converged point (the lane floors or the noise-floor rule). The step grows
    x1.5 after a convergence in at most 4 iterations (cap 4), and halves on failure (``halvings``).
    """
    _, jz = rj(zb, sb)
    d = np.linalg.norm(jz, axis=0)
    d[d == 0.0] = 1.0
    j = int(np.argmax(np.abs((zb - za) * d)))
    for _ in range(halvings):
        target = zb[j] + ratio * (zb[j] - za[j])
        y = np.concatenate([zb + ratio * (zb - za), [sb + ratio * (sb - sa)]])
        for it in range(iters):
            s = float(y[-1])
            try:
                r, jzz = rj(y[:-1], s)
            except RuntimeError:
                break
            if ok(y[:-1], s) and abs(y[j] - target) < 1e-9 * max(1.0, abs(target)):
                nxt = min(4.0, 1.5 * ratio) if it <= 4 else ratio
                return MonitorStep(y[:-1].copy(), s, True, nxt, it, j)
            try:
                js = sigma_derivative(rj, y[:-1], s)
            except RuntimeError:
                break
            row = np.zeros(len(y))
            row[j] = 1.0
            jac = np.vstack([np.hstack([jzz, js[:, None]]), row[None, :]])
            f = np.concatenate([r, [y[j] - target]])
            if not (np.all(np.isfinite(jac)) and np.all(np.isfinite(f))):
                break
            dc = np.linalg.norm(jac, axis=0)
            dc[dc == 0.0] = 1.0
            try:
                step = np.linalg.lstsq(jac / dc, -f, rcond=None)[0] / dc
            except np.linalg.LinAlgError:
                break
            f0 = float(np.linalg.norm(f))
            for alpha in [0.5**i for i in range(8)]:
                yt = y + alpha * step
                try:
                    rt, _ = rj(yt[:-1], float(yt[-1]))
                except RuntimeError:
                    continue
                ft = np.concatenate([rt, [yt[j] - target]])
                if np.all(np.isfinite(ft)) and float(np.linalg.norm(ft)) < f0:
                    y = yt
                    break
            else:
                break
        ratio *= 0.5
    return MonitorStep(zb.copy(), sb, False, ratio, 0, j)


def fold_supported(
    lower: list[tuple[float, float]] | None,
    upper: list[tuple[float, float]] | None,
    sigma: float,
    dv_min: float = 0.02,
    drp_min: float = 100.0,
) -> tuple[bool, list[tuple[float, float]]]:
    """The #1004 / #1025 two-solution rule. ``lower`` / ``upper`` are per-node (V_inf km/s, r_p km)
    of the two fixed-sigma solutions (None if one did not converge). Supported if both converged and
    some node differs by more than ``dv_min`` in V_inf and ``drp_min`` in r_p / sigma."""
    if lower is None or upper is None:
        return False, []
    dif = [(abs(a[0] - b[0]), abs(a[1] - b[1]) / sigma) for a, b in zip(lower, upper, strict=True)]
    return any(dv > dv_min and drp > drp_min for dv, drp in dif), dif
