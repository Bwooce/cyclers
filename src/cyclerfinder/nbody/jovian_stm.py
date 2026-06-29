"""Self-contained state+STM co-integrator for the Jupiter-central model (#480).

The Jovian propagator :class:`cyclerfinder.nbody.jovian.JovianRestrictedNBody` is
REBOUND/IAS15 with the moons applied through a Python ``additional_forces``
callback. REBOUND's order-1 variational particles do NOT differentiate that
Python callback (``memory/reference_rebound_variation_custom_force_gotcha``), so a
REBOUND-variational STM would be Jupiter-gravity-only — wrong exactly where the
EGGIE corrector needs it (near a Galilean flyby). The finite-difference Jacobian of
the same propagator is noise-limited and is what walls the Stage-2 close
(``docs/notes/2026-06-29-480-eggie-stage2-nbody-verdict.md``).

This module sidesteps both by co-integrating the spacecraft state AND the 6x6 state
transition matrix with an ANALYTIC gravity gradient — no REBOUND variational
particles. The dynamics are Jupiter-central with the Galilean moons as point masses
on the injected ephemeris rails (the same model the corrector residual uses):

* a(r,t) = -mu_J r/|r|^3 + sum_m mu_m [ (r_m - r)/|r_m - r|^3 - r_m/|r_m|^3 ]
  (the bracketed second term is the Jupiter-frame indirect term — the moon's pull
  on Jupiter, which the central frame must subtract; it does not depend on r so it
  contributes nothing to the gradient).
* Gravity gradient G(r,t) = da/dr = -mu_J/|r|^3 (I - 3 rhat rhat^T)
  + sum_m -mu_m/|d_m|^3 (I - 3 dhat_m dhat_m^T),  d_m = r_m - r.
  (d_m = r_m - r so d(d_m)/dr = -I; the two sign flips cancel, leaving the moon
  gradient with the SAME outer-product form as the central term.)
* Moon-surface softening (matching :class:`JovianRestrictedNBody`): below the moon
  surface the distance is clamped to ``radius_eq_km`` so a below-surface dive stays
  finite (a divergent-seed signal, never a NaN). In that softened regime the moon
  acceleration is mu_m (d/r_surf^3 - r_m/|r_m|^3), whose gradient is -mu_m/r_surf^3 I
  (a constant cube, no outer-product term). A real >= ~25 km flyby is above the
  surface and integrates exactly.

The STM ODE is dPhi/dt = F(t) Phi with F = [[0, I],[G(r,t), 0]], Phi(t0) = I_6,
co-integrated with the state via scipy ``solve_ivp`` (DOP853). The moon positions
r_m(t) are read from the injected ephemeris at the integrator's own time t (not a
coarse rails cache) for gradient fidelity.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.nbody.jovian import GALILEAN, MU_JUPITER_KM3_S2

Vec3 = NDArray[np.float64]
Mat6 = NDArray[np.float64]


class _EphemLike(Protocol):
    """Duck-typed moon ephemeris (``JovianEphemeris`` / ``IdealJovianEphemeris``)."""

    def state(self, moon: str, t_sec: float) -> tuple[Vec3, Vec3]: ...


def _accel_and_gradient(
    r: Vec3,
    t_sec: float,
    *,
    ephem: _EphemLike,
    moons: Sequence[str],
    mus: dict[str, float],
    surf: dict[str, float],
) -> tuple[Vec3, Mat6]:
    """Return ``(a, G)``: Jupiter-central acceleration and its 3x3 gravity gradient.

    ``a`` includes the central, moon-direct, and Jupiter-frame indirect terms;
    ``G = da/dr`` includes the central + moon-direct gradients (the indirect term is
    r-independent). Moon-surface softening is applied consistently to both.
    """
    eye3 = np.eye(3)
    rn = float(np.linalg.norm(r))
    rhat = r / rn
    a = -MU_JUPITER_KM3_S2 * r / rn**3
    grad = -MU_JUPITER_KM3_S2 / rn**3 * (eye3 - 3.0 * np.outer(rhat, rhat))
    for m in moons:
        r_m = np.asarray(ephem.state(m, float(t_sec))[0], dtype=np.float64)
        d = r_m - r
        dn = float(np.linalg.norm(d))
        rm3 = float(np.linalg.norm(r_m)) ** 3
        r_surf = surf[m]
        if dn >= r_surf:
            d_eff3 = dn**3
            dhat = d / dn
            grad = grad - mus[m] / dn**3 * (eye3 - 3.0 * np.outer(dhat, dhat))
        else:
            # Below the moon surface: distance clamped to r_surf (constant cube),
            # so a = mu_m (d/r_surf^3 - ...) and da/dr = -mu_m/r_surf^3 I.
            d_eff3 = r_surf**3
            grad = grad - mus[m] / r_surf**3 * eye3
        a = a + mus[m] * (d / d_eff3 - r_m / rm3)
    return a, grad


def propagate_with_stm(
    r0: Vec3,
    v0: Vec3,
    t0_sec: float,
    t1_sec: float,
    *,
    ephem: _EphemLike,
    moons: Sequence[str] = GALILEAN,
    rtol: float = 1e-11,
    atol: float = 1e-9,
) -> tuple[Vec3, Vec3, Mat6]:
    """Co-integrate the state and the 6x6 STM over one Jupiter-central leg.

    Parameters mirror :meth:`JovianRestrictedNBody.propagate` (Jupiter-central,
    moons as point masses at ``ephem.state(moon, t)[0]``). Returns ``(rf, vf, Phi)``
    where ``Phi = d(rf, vf)/d(r0, v0)`` is the 6x6 state transition matrix.

    The ephemeris is evaluated at the integrator's own time ``t`` (not a coarse
    rails cache) so the gravity gradient matches the propagated trajectory exactly.
    Moon-surface softening matches :class:`JovianRestrictedNBody`: a real flyby
    above the moon surface integrates exactly (the STM is valid there).
    """
    moons = tuple(moons)
    mus = {m: SATELLITES[m].mu_km3_s2 for m in moons}
    surf = {m: SATELLITES[m].radius_eq_km for m in moons}
    r0 = np.asarray(r0, dtype=np.float64)
    v0 = np.asarray(v0, dtype=np.float64)

    def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        r = y[0:3]
        v = y[3:6]
        phi = y[6:].reshape(6, 6)
        a, grad = _accel_and_gradient(r, t, ephem=ephem, moons=moons, mus=mus, surf=surf)
        f = np.zeros((6, 6), dtype=np.float64)
        f[0:3, 3:6] = np.eye(3)
        f[3:6, 0:3] = grad
        dphi = f @ phi
        dy = np.empty(42, dtype=np.float64)
        dy[0:3] = v
        dy[3:6] = a
        dy[6:] = dphi.ravel()
        return dy

    y0 = np.concatenate([r0, v0, np.eye(6, dtype=np.float64).ravel()])
    sol = solve_ivp(
        rhs,
        (float(t0_sec), float(t1_sec)),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not sol.success:  # pragma: no cover - defensive; DOP853 rarely fails here
        raise RuntimeError(f"propagate_with_stm integration failed: {sol.message}")
    yf = sol.y[:, -1]
    rf = np.asarray(yf[0:3], dtype=np.float64)
    vf = np.asarray(yf[3:6], dtype=np.float64)
    phi = np.asarray(yf[6:].reshape(6, 6), dtype=np.float64)
    return rf, vf, phi


__all__ = ["propagate_with_stm"]
