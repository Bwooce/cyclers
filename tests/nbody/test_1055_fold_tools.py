"""#1055 sec. 12.5: the generic monitor step and two-solution rule on an analytic fold at
production magnitudes (km near 1e5, a 1e3 velocity-like weight). Expected values are the toy's
analytic branches, set by construction."""

from __future__ import annotations

import numpy as np

from cyclerfinder.nbody.fold_tools import fold_supported, monitor_step

SF = 0.1639  # the fold, as gc-2's (#1004)
A = 2.0e5  # km per sqrt(sigma) on the branches


def toy(z: np.ndarray, s: float) -> tuple[np.ndarray, np.ndarray]:
    """Branches x = +- A sqrt(SF - s) (km); y = 1e-3 x (km/s, weighted 1e3)."""
    x, y = z
    r = np.array([(x / A) ** 2 + s - SF, 1e3 * (y - 1e-3 * x)])
    j = np.array([[2.0 * x / A**2, 0.0], [-1.0, 1e3]])
    return r, j


def ok(z: np.ndarray, s: float) -> bool:
    r, _ = toy(z, s)
    return bool(abs(r[0]) < 1e-12 and abs(r[1]) < 1e-9)


def point(s: float, sign: float) -> np.ndarray:
    x = sign * A * np.sqrt(SF - s)
    return np.array([x, 1e-3 * x])


def test_monitor_passes_the_fold_and_turns_back() -> None:
    za, sa = point(0.160, -1.0), 0.160
    zb, sb = point(0.163, -1.0), 0.163
    ratio = 1.0
    sig = [sb]
    for _ in range(25):
        st = monitor_step(toy, za, sa, zb, sb, ok, ratio=ratio)
        assert st.ok
        za, sa, zb, sb, ratio = zb, sb, st.z, st.s, st.ratio
        sig.append(sb)
        if zb[0] > 0 and sb < 0.15:
            break
    assert max(sig) <= SF + 1e-12  # never beyond the fold
    assert zb[0] > 0  # on the other branch
    assert sb < 0.15  # sigma has turned back
    assert abs(zb[0] - A * np.sqrt(SF - sb)) < 1e-6 * A  # exactly on the upper branch


def test_two_solution_rule() -> None:
    s = 0.150
    lo = [(3.5756, 11507.0 * s), (2.9815, 11507.0 * s)]
    up = [(3.5548, 10849.0 * s), (2.9340, 10849.0 * s)]
    sup, dif = fold_supported(lo, up, s)
    assert sup  # #1004 sec. 8.3 numbers: 0.047 km/s and 658 km
    assert dif[1][0] > 0.02 and dif[1][1] > 100.0
    assert not fold_supported(lo, lo, s)[0]
    assert not fold_supported(lo, None, s)[0]
