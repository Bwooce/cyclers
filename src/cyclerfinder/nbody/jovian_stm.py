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
from typing import TYPE_CHECKING, Protocol

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.nbody.jovian import _W_VEL, GALILEAN, MU_JUPITER_KM3_S2

if TYPE_CHECKING:
    from cyclerfinder.nbody.jovian_ideal import SubarcSeed
    from cyclerfinder.nbody.shooter import ShootingSeed

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


def jovian_stm_jacobian(
    seed: ShootingSeed,
    x: NDArray[np.float64],
    *,
    ephem: _EphemLike,
    moons: Sequence[str] = GALILEAN,
    rtol: float = 1e-11,
    atol: float = 1e-9,
) -> NDArray[np.float64]:
    """Analytic block-bidiagonal Jacobian of :func:`jovian_defect_residual`.

    Jupiter-central analogue of :func:`cyclerfinder.nbody.shooter._stm_jacobian`,
    assembled from per-leg state-transition matrices (one analytic co-integrated
    :func:`propagate_with_stm` per leg) instead of the ``6*n_nodes+1`` REBOUND
    finite-difference re-propagations. The free variables are the per-node full
    Cartesian states (6 each, the :func:`cyclerfinder.nbody.shooter._states_to_x`
    packing); the residual is, in order, the ``n-1`` leg continuity defects (6
    each), the ``n-2`` interior flyby hinges (1 each), and the 6-component
    periodicity wrap — exactly the layout of
    :func:`cyclerfinder.nbody.jovian.jovian_defect_residual`.

    The crucial difference from the heliocentric Jacobian: the Jovian residual
    weights every velocity row by ``_W_VEL`` (position km vs velocity km/s on
    comparable least-squares scales). So both the per-leg ``Phi`` block AND the
    ``-I`` / ``+I`` coupling blocks carry a velocity-row scaling
    ``diag(1,1,1,W,W,W)``:

    - Leg defect ``c_i = propagate(node_i, leg_i) - node_{i+1}`` gives
      ``dc_i/dnode_i = R_W @ Phi_i`` and ``dc_i/dnode_{i+1} = -R_W`` (the ``-I_6``
      with velocity rows scaled), where ``R_W = diag(1,1,1,W,W,W)``.
    - The flyby hinges read ``seed.vinf_in/out`` (carried constants), so their rows
      are identically zero — matching what FD sees.
    - The wrap residual ``(node_{n-1} - r_wrap_pl) - (node_0 - r_home)`` (velocity
      rows ``*W``) gives ``-R_W`` on node 0 and ``+R_W`` on node ``n-1`` (the moon
      states are epoch-fixed constants).

    A leg whose analytic propagation fails leaves its ``Phi`` block zero (the
    ``-R_W`` coupling is kept) — the honest local linearisation of the residual's
    divergence sentinel (a constant in ``x`` there), matching the FD oracle.
    """
    from cyclerfinder.nbody.shooter import _STATE_DIM, _x_to_states

    moons = tuple(moons)
    n = len(seed.sequence)
    states = _x_to_states(x, n)

    row_w = np.array([1.0, 1.0, 1.0, _W_VEL, _W_VEL, _W_VEL], dtype=np.float64)
    rw = np.diag(row_w)  # velocity-row weighting (matches jovian_defect_residual)

    n_leg = (n - 1) * _STATE_DIM
    n_hinge = max(0, n - 2)
    n_rows = n_leg + n_hinge + _STATE_DIM
    n_cols = n * _STATE_DIM
    jac = np.zeros((n_rows, n_cols), dtype=np.float64)

    for i in range(n - 1):
        s_i = states[i]
        r0 = np.asarray(s_i[:3], dtype=np.float64)
        v0 = np.asarray(s_i[3:], dtype=np.float64)
        rows = slice(i * _STATE_DIM, (i + 1) * _STATE_DIM)
        try:
            _, _, phi = propagate_with_stm(
                r0,
                v0,
                seed.epochs[i],
                seed.epochs[i + 1],
                ephem=ephem,
                moons=moons,
                rtol=rtol,
                atol=atol,
            )
        except Exception:
            phi = None
        if phi is not None and np.all(np.isfinite(phi)):
            jac[rows, i * _STATE_DIM : (i + 1) * _STATE_DIM] = rw @ phi
        jac[rows, (i + 1) * _STATE_DIM : (i + 2) * _STATE_DIM] = -rw

    # Hinge rows (n_leg : n_leg + n_hinge) are constant in x -> left zero.

    # Periodicity wrap rows.
    wrap = slice(n_leg + n_hinge, n_leg + n_hinge + _STATE_DIM)
    jac[wrap, 0:_STATE_DIM] = -rw
    jac[wrap, (n - 1) * _STATE_DIM : n * _STATE_DIM] = rw

    return jac


def subarc_stm_jacobian(
    sub: SubarcSeed,
    x: NDArray[np.float64],
    *,
    ephem: _EphemLike,
    moons: Sequence[str] = GALILEAN,
    rtol: float = 1e-11,
    atol: float = 1e-9,
) -> NDArray[np.float64]:
    """Analytic block-bidiagonal Jacobian of :func:`subarc_defect_residual` (#480 Stage 4).

    Sub-arc generalisation of :func:`jovian_stm_jacobian`: the free variables are the
    per-node full states over ALL ``m`` nodes (encounter + interior continuity nodes,
    :func:`cyclerfinder.nbody.shooter._states_to_x` packing); the residual is, in
    order, the ``m-1`` sub-arc continuity defects (6 each), the ``n_enc-2`` interior
    flyby hinges (constant in ``x`` -> zero rows), and the 6-component periodicity wrap
    between the first and last encounter nodes — exactly the layout of
    :func:`cyclerfinder.nbody.jovian_ideal.subarc_defect_residual`.

    Each sub-arc ``j`` (node ``j`` -> node ``j+1``) contributes
    ``dc_j/dnode_j = R_W @ Phi_j`` and ``dc_j/dnode_{j+1} = -R_W`` with the same
    velocity-row weighting ``R_W = diag(1,1,1,W,W,W)`` as the one-node-per-leg case
    (``Phi_j`` from one co-integrated :func:`propagate_with_stm` over the sub-arc). A
    sub-arc whose analytic propagation fails leaves its ``Phi`` block zero (keeping the
    ``-R_W`` coupling) — the honest local linearisation of the divergence sentinel.
    With ``sub.n_subarcs == 1`` this is identical to :func:`jovian_stm_jacobian`.
    """
    from cyclerfinder.nbody.shooter import _STATE_DIM, _x_to_states

    moons = tuple(moons)
    m = len(sub.node_states)
    n_enc = len(sub.encounter_idx)
    states = _x_to_states(x, m)

    row_w = np.array([1.0, 1.0, 1.0, _W_VEL, _W_VEL, _W_VEL], dtype=np.float64)
    rw = np.diag(row_w)  # velocity-row weighting (matches subarc_defect_residual)

    n_leg = (m - 1) * _STATE_DIM
    n_hinge = max(0, n_enc - 2)
    n_rows = n_leg + n_hinge + _STATE_DIM
    n_cols = m * _STATE_DIM
    jac = np.zeros((n_rows, n_cols), dtype=np.float64)

    for j in range(m - 1):
        s_j = states[j]
        r0 = np.asarray(s_j[:3], dtype=np.float64)
        v0 = np.asarray(s_j[3:], dtype=np.float64)
        rows = slice(j * _STATE_DIM, (j + 1) * _STATE_DIM)
        try:
            _, _, phi = propagate_with_stm(
                r0,
                v0,
                sub.epochs[j],
                sub.epochs[j + 1],
                ephem=ephem,
                moons=moons,
                rtol=rtol,
                atol=atol,
            )
        except Exception:
            phi = None
        if phi is not None and np.all(np.isfinite(phi)):
            jac[rows, j * _STATE_DIM : (j + 1) * _STATE_DIM] = rw @ phi
        jac[rows, (j + 1) * _STATE_DIM : (j + 2) * _STATE_DIM] = -rw

    # Hinge rows (n_leg : n_leg + n_hinge) are constant in x -> left zero.

    # Periodicity wrap rows (between the first and last encounter nodes).
    i0 = sub.encounter_idx[0]
    i_last = sub.encounter_idx[-1]
    wrap = slice(n_leg + n_hinge, n_leg + n_hinge + _STATE_DIM)
    jac[wrap, i0 * _STATE_DIM : (i0 + 1) * _STATE_DIM] = -rw
    jac[wrap, i_last * _STATE_DIM : (i_last + 1) * _STATE_DIM] = rw

    return jac


__all__ = ["jovian_stm_jacobian", "propagate_with_stm", "subarc_stm_jacobian"]
