from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
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


def _fd_stm_colsteps(
    prop: RestrictedNBody,
    r0: npt.NDArray[np.float64],
    v0: npt.NDArray[np.float64],
    t0: float,
    t1: float,
    **kw: Any,
) -> npt.NDArray[np.float64]:
    """Per-column FD STM mirroring Task 1's fair reference (pos h=1 km, vel h=1e-3 km/s)."""
    return _fd_stm(prop, r0, v0, t0, t1, h=(1.0, 1.0, 1.0, 1e-3, 1e-3, 1e-3), **kw)


@pytest.mark.slow
def test_perturbed_stm_matches_fd() -> None:
    """Perturbed-parity gate (#388) on a WEAKLY-perturbed 1-AU heliocentric leg.

    On this geometry the s/c sits ~1 AU from the Sun and far from Earth/Mars for
    the whole 200 d, so the rails perturber shifts the terminal state by only
    ~5e3 km (on a ~1.5e8 km arc) and changes Phi itself by only ~2e-5 relative.
    The variational STM therefore matches FD to ~1.3e-4 here — but that does NOT
    prove REBOUND's variation includes the custom additional_forces gradient; it
    only proves the gradient is negligible on THIS leg. The decisive stress test
    is ``test_perturbed_stm_earthclose_omits_rails_gradient`` below, which forces
    the perturber gradient to be large and exposes that the variation omits it.
    """
    prop = RestrictedNBody()
    ephem = Ephemeris("astropy")
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    t0, t1 = 0.0, 200.0 * 86400.0
    kw: dict[str, Any] = {"bodies": ("E", "M"), "ephem": ephem}
    arc = prop.propagate(r0, v0, t0, t1, with_stm=True, **kw)
    assert arc.stm is not None
    # column-appropriate FD reference (pos cols h=1 km, vel cols h=1e-3 km/s),
    # matching the fair FD reference Task 1 established.
    fd = _fd_stm_colsteps(prop, r0, v0, t0, t1, **kw)
    rel = np.linalg.norm(arc.stm - fd) / np.linalg.norm(fd)
    # Per-column rel errors: reveals whether the velocity columns / the
    # perturbed-direction columns are the ones that diverge.
    col_rel = np.linalg.norm(arc.stm - fd, axis=0) / np.linalg.norm(fd, axis=0)
    print(f"PERTURBED STM vs FD rel = {rel:.3e}")
    print("PERTURBED STM per-column rel = " + ", ".join(f"{c:.3e}" for c in col_rel))
    assert rel < 5e-3, (
        f"perturbed STM vs FD rel={rel} — REBOUND variation likely omits "
        "the rails additional_forces gradient"
    )


@pytest.mark.slow
def test_perturbed_stm_earthclose_includes_rails_gradient() -> None:
    """Earth-close leg where the rails gradient is large — the analytic fix.

    Seeded ~2e6 km from Earth (well outside the SOI softening clamp) so the Earth
    third-body term materially shapes the STM (||Phi_pert - Phi_sun|| / ||Phi_sun||
    ~ 4e-2). REBOUND's native first-order variation differentiates only the in-sim
    Sun gravity, NOT the Python ``additional_forces`` rails callback. The callback
    therefore applies the perturber gravity-gradient tensor
    ``G_p = mu_p (3 d d^T/|d|^5 - I/|d|^3)`` (``d = r_p - r_sc``) to the variational
    particles, so the co-integrated variational equations carry the full dynamics
    and the STM matches a full-dynamics FD reference inside the 5e-3 gate — the
    leg that previously diverged by ~7.5e-2 with the gradient omitted.
    """
    prop = RestrictedNBody()
    ephem = Ephemeris("astropy")
    r_earth, _ = ephem.state("E", 0.0)
    r0 = r_earth + np.array([2.0e6, 0.0, 0.0])
    v_hat = np.array([-r0[1], r0[0], 0.0])
    v_hat /= np.linalg.norm(v_hat)
    v0 = 29.0 * v_hat
    t0, t1 = 0.0, 200.0 * 86400.0
    kw: dict[str, Any] = {"bodies": ("E", "M"), "ephem": ephem}
    arc = prop.propagate(r0, v0, t0, t1, with_stm=True, **kw)
    assert arc.stm is not None
    fd = _fd_stm_colsteps(prop, r0, v0, t0, t1, **kw)
    rel = np.linalg.norm(arc.stm - fd) / np.linalg.norm(fd)
    col_rel = np.linalg.norm(arc.stm - fd, axis=0) / np.linalg.norm(fd, axis=0)
    print(f"EARTH-CLOSE STM vs FD rel = {rel:.3e}")
    print("EARTH-CLOSE STM per-column rel = " + ", ".join(f"{c:.3e}" for c in col_rel))
    assert rel < 5e-3, (
        f"Earth-close perturbed STM vs FD rel={rel} — REBOUND variation omits "
        "the rails additional_forces gradient (STM is Sun-gravity-only)"
    )
