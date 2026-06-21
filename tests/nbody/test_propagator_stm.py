from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from cyclerfinder.nbody.propagator import NBodyArc, RestrictedNBody


def _final_vec(arc: NBodyArc) -> npt.NDArray[np.float64]:
    return np.concatenate([np.asarray(arc.r_km), np.asarray(arc.v_km_s)])


def _fd_stm(
    prop: RestrictedNBody,
    r0: npt.NDArray[np.float64],
    v0: npt.NDArray[np.float64],
    t0: float,
    t1: float,
    *,
    h: float | tuple[float, ...] = (1.0, 1.0, 1.0, 1e-3, 1e-3, 1e-3),
    **kw: Any,
) -> npt.NDArray[np.float64]:
    # Per-column step. The position columns (0..2) tolerate a 1 km probe, but a
    # 1 km/s velocity probe over a multi-month arc is deep in the nonlinear
    # regime (the v->r block scales like the time-of-flight, ~1e7 s), so the
    # velocity columns (3..5) need a far smaller probe for a fair FD reference.
    steps = (h,) * 6 if isinstance(h, (int, float)) else tuple(h)
    x0 = np.concatenate([r0, v0])
    f0 = _final_vec(prop.propagate(x0[:3], x0[3:], t0, t1, **kw))
    jac = np.empty((6, 6))
    for j in range(6):
        xp = x0.copy()
        xp[j] += steps[j]
        jac[:, j] = (_final_vec(prop.propagate(xp[:3], xp[3:], t0, t1, **kw)) - f0) / steps[j]
    return jac


def test_with_stm_false_is_unchanged() -> None:
    prop = RestrictedNBody()
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    arc = prop.propagate(r0, v0, 0.0, 1.0e7)
    assert getattr(arc, "stm", None) is None


def test_sun_only_stm_matches_fd() -> None:
    prop = RestrictedNBody()
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    t0, t1 = 0.0, 120.0 * 86400.0
    arc = prop.propagate(r0, v0, t0, t1, with_stm=True)
    assert arc.stm is not None and arc.stm.shape == (6, 6)
    fd = _fd_stm(prop, r0, v0, t0, t1)
    rel = np.linalg.norm(arc.stm - fd) / np.linalg.norm(fd)
    assert rel < 1e-3, f"Sun-only STM vs FD rel={rel}"
