"""Seedless 2D pseudospectral quasi-periodic-torus corrector for the QBCP (#617).

Motivation and the wall this crosses
-------------------------------------
This module combines the two precedents `#611` and `#612` built this session
into the one corrector `#544`'s deepest diagnosis named as the real fix for the
`#538` cislunar-cycler blocker: a torus corrector that is BOTH (a) genuinely 2D
quasi-periodic (two angles, a transverse rotation number -- like `#612`'s CR3BP
corrector ``search.variational_qp_torus``) AND (b) adapted to the QBCP's
non-autonomous, time-periodic, fixed-period, first-order-6-canonical-state
setting (like `#611`'s periodic-orbit corrector
``search.variational_periodic_orbit_qbcp``).

The concrete blocker is that ``genome.qbcp_torus.correct_qbcp_torus`` (the
classical GMOS / Gomez-Mondelo-Olikara-Scheeres invariant-circle corrector)
resolves the LONGITUDINAL torus angle by stroboscopic FLOW INTEGRATION:
``qbcp_torus_residual`` propagates every invariant-circle sample forward through
the true nonlinear ``qbcp_eom`` for one full Sun-synodic period ``T_s`` and
matches the rotated circle. This is fine for the mildly-unstable Sun-Earth L2
torus (converges cleanly, invariance residual ~1.5e-5 to 3e-5) but structurally
hopeless for the violently-unstable Earth-Moon L1/L2 region: `#544` and `#611`
both measured a frozen-time linearization rate of ~2-3 over ``T_s ~ 6.79``,
implying a one-period amplification ``exp(rate * T_s) ~ 1e6-1e8`` -- so every
single-period propagation (and every finite-difference column of its Jacobian)
blows any model-instance / roundoff offset up to O(1). `#544`'s own conclusion,
verbatim: "The real fix is a multiple-shooting GMOS corrector ... NOT more
Fourier modes, a different target orbit, or coefficient hunting." `#611`/`#612`
showed a DIFFERENT fix than multiple-shooting also works, and is cleaner: an
integration-free pseudospectral (harmonic-balance / collocation) formulation in
which the unstable flow's amplification never enters the residual or its
Jacobian at all.

The single most important design decision: omega1 is FIXED, not free
-------------------------------------------------------------------
A quasi-periodic invariant 2-torus of the QBCP is a map ``u : T^2 -> R^6`` from
two angles ``(theta1, theta2)`` to the canonical (PM) 6-state, on which the flow
is linear: ``theta1' = omega1``, ``theta2' = omega2``. The QBCP is
NON-AUTONOMOUS and time-periodic with period ``T_s`` (its ``qbcp_eom`` depends
on absolute time ``t`` ONLY through the Sun phase ``theta_sun = theta_sun0 +
omega_s * t``, an explicit ``T_s``-periodic forcing). This forces a specific,
asymmetric treatment of the two angles that differs from BOTH precedents:

* ``theta1`` is the model's own time-phase, LOCKED to the forcing:
  ``theta1 = omega_s * t`` (with the standard ``theta_sun0 = 0`` epoch, so
  ``theta_sun = theta1`` exactly, or ``= period_multiple * theta1`` for a
  subharmonic). Its frequency ``omega1 = omega_s = 2*pi/T_s`` is a FIXED input,
  NOT a solved-for unknown -- this is EXACTLY `#611`'s reasoning for the
  periodic-orbit case ("shifting a candidate trajectory in time does NOT produce
  another valid solution unless the epoch is shifted with it, so there is no
  continuous time-shift symmetry to gauge-fix; the period becomes a fixed input,
  not an unknown"), lifted to the first torus angle. Consequently the
  LONGITUDINAL phase gauge that `#612` needed (pinning ``sin(theta1)`` of ``x``
  to remove the CR3BP's autonomous time-shift degeneracy) has NO analogue here
  and is DROPPED -- theta1 has no phase freedom to fix, its origin is pinned by
  the Sun epoch.
* ``theta2`` is the genuine internal quasi-periodic angle (the libration /
  center-manifold oscillation the torus wraps). Its frequency ``omega2`` is
  FREE, exactly like `#612`'s transverse frequency. The rotation number
  ``rho = omega2 / omega1`` (equivalently the GMOS ``rho_strob = omega2 * T_s``,
  the advance of theta2 per stroboscopic period) selects the family member.

This is precisely what the existing GMOS ``correct_qbcp_torus`` already assumes
(``omega_long = 2*pi/T_s`` fixed; ``rho`` the free rotation number) -- so the
formulation here is consistent with the validated corrector, it only replaces
the shooting residual with an integration-free one. It differs from `#612`'s
CR3BP torus (where BOTH ``omega1`` and ``omega2`` were free, because the CR3BP
is autonomous -- neither angle is locked to an external clock) in dropping
``omega1`` from the unknowns entirely.

The invariance PDE (no integration anywhere in the search)
----------------------------------------------------------
A trajectory ``x(t) = u(omega1 t + phi1, omega2 t + phi2)`` satisfies
``x' = f(t, x)`` (the QBCP canonical RHS :func:`cyclerfinder.core.qbcp.qbcp_eom`)
if and only if ``u`` solves the quasi-periodic invariance PDE (chain rule)

    omega1 * du/dtheta1  +  omega2 * du/dtheta2  =  F(theta1, u(theta1, theta2))

for all ``(theta1, theta2) in T^2``, where ``F(theta1, u)`` is ``qbcp_eom``
evaluated with the Sun phase set to ``theta_sun = period_multiple * theta1``
(i.e. at time ``t = theta1 / omega1``). The explicit non-autonomous
time-dependence of the RHS thus becomes an explicit dependence of ``F`` on the
FIRST torus angle -- clean and exact. Crucially the PDE contains NO time
integration: ``F(theta1, u)`` is evaluated pointwise and the angle-derivatives
are taken SPECTRALLY, so the ~1e6-1e8 per-period amplification never enters the
residual or its Jacobian. The residual's Jacobian is the local operator
``omega1 d/dtheta1 + omega2 d/dtheta2 - D_u F(theta1, u)``, well-conditioned
regardless of how unstable the flow is -- the whole reason `#611`/`#612` beat
shooting, now lifted to the QBCP 2-torus.

Each of the six canonical state components (x, y, z, px, py, pz) is a real
tensor-product Fourier series over both angles (mirroring `#612`, one extra
angle vs `#611`):

    u_c(theta1, theta2) = sum_{a,b} C[c,a,b] * phi_a(theta1) * phi_b(theta2)

with the real 1D Fourier basis ``phi_0 = 1, phi_k = cos(k .), phi_{N+k} =
sin(k .)`` (reality automatic). Angle-derivatives are analytic term-by-term. The
RHS ``F`` is the vectorized ``qbcp_eom`` (see :func:`_qbcp_rhs_grid`); its
6x6 state-Jacobian is exactly ``qbcp_stm_eom``'s ``jac_a`` (see
:func:`_qbcp_jacobian_grid`), evaluated pointwise with the theta1-dependent
alphas.

Free unknowns and gauge fixing (three rows, one fewer than `#612`)
------------------------------------------------------------------
Free unknowns: the coefficients ``C`` (shape ``(6, 2N1+1, 2N2+1)``) and ONLY
``omega2`` (``omega1`` is fixed). Because ``theta1`` is locked, only ONE
continuous phase-shift symmetry survives (``theta2 -> theta2 + d2``), plus the
family-collapse degeneracy -- so THREE scalar gauge/anchor rows, not `#612`'s
four (the longitudinal phase gauge is dropped):

* **Transverse phase gauge** (``theta2``): pin the ``sin(theta2)``,
  ``theta1``-constant coefficient of a chosen coordinate ``c*`` (the one with
  the largest transverse-mode content in the guess) to zero -- `#612`'s
  transverse gauge, kept verbatim.
* **Transverse amplitude anchor**: fix the L2 norm of all ``theta2``-first-
  harmonic (``k2 = 1``) coefficients to a caller value. Without it the optimizer
  collapses the transverse structure to zero -- i.e. slides back to the parent
  PERIODIC orbit (`#611`'s object, the torus's zero-amplitude center). This is
  the "how big a torus" family-selection choice, not information from a target.
* **Rotation-number pin**: ``omega2 = rho_target * omega1``. `#612` found this
  load-bearing (without it the least-squares escapes to a degenerate
  ``omega2 -> 0`` resonant branch that reaches machine precision and always
  beats the truncation-limited genuine torus). Since ``omega1`` is fixed here,
  this simply pins ``omega2`` to the seed's rotation number; ``omega2`` stays a
  formal unknown so the corrector can make tiny adjustments and the Jacobian
  structure matches `#612`.

Bootstrap (integration allowed ONLY in the initial guess)
---------------------------------------------------------
Like `#611`/`#612` the SEARCH never integrates, but a rough guess is needed. For
the SE-L2 positive control the honest bootstrap is `#612`'s "one known family
member" pattern: sample a converged GMOS :class:`~genome.qbcp_torus.QBCPTorus`
on the 2D grid via ``evaluate_qbcp_torus`` (this propagates -- fine, it is only
the guess), convert PV->PM, least-squares project onto the 2D Fourier basis, and
hand off. See :func:`project_gmos_qbcp_torus_to_2d` and
:func:`discover_qbcp_torus_from_gmos`.

Independent closure check (NOT circular with the residual)
----------------------------------------------------------
The minimized quantity is the algebraic PDE residual on the collocation grid.
The independent check propagates a torus point ``u(theta1, theta2)`` (a PM state
at time ``t0 = theta1/omega1``) through the TRUE nonlinear ``qbcp_eom`` for a
SHORT time ``dt`` and compares against ``u(theta1 + omega1 dt, theta2 + omega2
dt)`` -- genuinely independent (it uses the integrator, not the spectral
series). ``dt`` is a small fraction of the period so the violent-instability
amplification stays bounded (see `#612`'s identical reasoning).

Discipline / scope
------------------
The returned torus is OUR computation; no novelty claim, no catalogue writeback
-- a capability build (a corrector crossing a documented convergence wall), not
a discovery result. This module deliberately does NOT re-run `#538`'s full
cross-system SE<->EM boundary-value chain (portfolio-parked 2026-07-10); its
scope is build + validate only.

References
----------
* Jorba, A. (2001). Nonlinearity 14, 943-976 (invariant curves of maps).
* Jorba, A., & Villanueva, J. (1997). On the persistence of lower-dimensional
  invariant tori under quasi-periodic perturbations. J. Nonlinear Sci.
* Olikara, Z., & Scheeres, D. (2010). Numerical method for computing quasi-
  periodic orbits (the GMOS invariant-circle reduction replaced here).
* Gimeno, J., & Jorba, A. (2018) / Rosales, J., & Jorba, A. (2023) -- the QBCP
  model instance this codebase implements in ``core.qbcp``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.genome.qbcp_torus import QBCPTorus, evaluate_qbcp_torus
from cyclerfinder.search.outcome_log import log_outcome

_N_STATE = 6  # canonical PM state: x, y, z, px, py, pz


@dataclass(frozen=True)
class QBCPTorusVariationalResult:
    """A QBCP quasi-periodic 2-torus found by seedless 2D pseudospectral collocation.

    ``coeffs`` is the ``(6, 2*n1+1, 2*n2+1)`` real tensor-product Fourier
    coefficient array (PM canonical state; see module docstring for the basis);
    ``omega1`` is FIXED at ``period_multiple * omega_s = 2*pi/period`` (the Sun's
    synodic frequency -- NOT solved for), ``omega2`` is the free transverse
    frequency (rad/TU), ``rotation_number = omega2 / omega1``, and
    ``rho_strob = omega2 * period`` is the advance of theta2 per stroboscopic
    period (the GMOS ``rho``). ``residual_rms`` is the RMS of the quasi-periodic
    invariance-PDE residual at the 2D collocation points (the quantity
    minimized). ``closure_residual`` is the INDEPENDENT short-time flow check
    (propagate ``u(theta1,theta2)`` for ``closure_dt`` through the true nonlinear
    ``qbcp_eom``, compare to ``u(theta1+omega1*dt, theta2+omega2*dt)`` -- NOT
    circular with ``residual_rms``). ``transverse_amplitude`` is the anchored L2
    norm of the ``k2=1`` coefficients.
    """

    system: qbcp.QBCPSystem
    coeffs: NDArray[np.float64]
    omega1: float
    omega2: float
    rotation_number: float
    rho_strob: float
    period: float
    n1: int
    n2: int
    m1: int
    m2: int
    period_multiple: int
    transverse_amplitude: float
    residual_rms: float
    closure_residual: float
    converged: bool
    n_iter: int
    notes: str = ""
    extras: dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Real tensor-product Fourier basis (identical convention to #612).
# ---------------------------------------------------------------------------


def _basis_matrices(
    n_modes: int, thetas: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ``(P, Pd)`` of shape ``(len(thetas), 2*n_modes+1)``: the real 1D
    Fourier basis ``[1, cos(1.), .., cos(N.), sin(1.), .., sin(N.)]`` at
    ``thetas`` and its theta-derivative.
    """
    m = thetas.size
    size = 2 * n_modes + 1
    mat_p = np.zeros((m, size))
    mat_pd = np.zeros((m, size))
    mat_p[:, 0] = 1.0
    for k in range(1, n_modes + 1):
        ck, sk = np.cos(k * thetas), np.sin(k * thetas)
        mat_p[:, k] = ck
        mat_pd[:, k] = -k * sk
        mat_p[:, n_modes + k] = sk
        mat_pd[:, n_modes + k] = k * ck
    return mat_p, mat_pd


def _n_free(n1: int, n2: int) -> int:
    """Free-variable count: all coefficients plus omega2 (omega1 is fixed)."""
    return _N_STATE * (2 * n1 + 1) * (2 * n2 + 1) + 1


def _pack(coeffs: NDArray[np.float64], omega2: float) -> NDArray[np.float64]:
    return np.concatenate([coeffs.reshape(-1), [omega2]])


def _unpack(z: NDArray[np.float64], n1: int, n2: int) -> tuple[NDArray[np.float64], float]:
    n_coef = _N_STATE * (2 * n1 + 1) * (2 * n2 + 1)
    coeffs = z[:n_coef].reshape(_N_STATE, 2 * n1 + 1, 2 * n2 + 1)
    return coeffs, float(z[-1])


def _k2_first_harmonic_cols(n2: int) -> tuple[int, int]:
    """Column indices in the theta2 basis for the k2=1 cos and sin terms."""
    return 1, n2 + 1


def evaluate_torus_state(
    result: QBCPTorusVariationalResult,
    theta1: float | NDArray[np.float64],
    theta2: float | NDArray[np.float64],
) -> NDArray[np.float64]:
    """Evaluate the 2D Fourier torus PM state ``u(theta1, theta2)`` (6-vector).

    Scalars in -> ``(6,)`` out. Uses the analytic Fourier series directly, no
    propagation.
    """
    t1 = np.atleast_1d(np.asarray(theta1, dtype=np.float64))
    t2 = np.atleast_1d(np.asarray(theta2, dtype=np.float64))
    p1, _ = _basis_matrices(result.n1, t1)
    p2, _ = _basis_matrices(result.n2, t2)
    out = np.empty((t1.size, _N_STATE))
    for c in range(_N_STATE):
        full = p1 @ result.coeffs[c] @ p2.T
        out[:, c] = np.diagonal(full) if t1.size == t2.size else full[:, 0]
    if np.isscalar(theta1) and np.isscalar(theta2):
        return out[0]
    return out


# ---------------------------------------------------------------------------
# Vectorized QBCP RHS and its state-Jacobian on the (theta1, theta2) grid.
# ---------------------------------------------------------------------------


def _alphas_on_theta1(
    theta1: NDArray[np.float64], omega1: float, system: qbcp.QBCPSystem
) -> NDArray[np.float64]:
    """Evaluate the 8 QBCP alpha functions at each theta1 grid value.

    ``theta1`` locks to the Sun phase via ``t = theta1 / omega1`` so that
    ``evaluate_alphas`` sees ``theta_sun = theta_sun0 + omega_s * t =
    theta_sun0 + (omega_s/omega1) * theta1`` (``= theta1`` for
    ``period_multiple=1`` and the default ``theta_sun0=0``). Returns shape
    ``(8, m1)`` in order (a1, a2, a3, a4, a5, a6, xs=alpha7, ys=alpha8).
    """
    m1 = theta1.size
    out = np.empty((8, m1), dtype=np.float64)
    for i in range(m1):
        t = float(theta1[i]) / omega1
        al = qbcp.evaluate_alphas(t, system)
        out[:, i] = al[1:9]
    return out


def _qbcp_rhs_grid(
    u: NDArray[np.float64], alphas: NDArray[np.float64], mu: float, mu_sun: float
) -> NDArray[np.float64]:
    """Vectorized ``qbcp_eom`` on the grid: ``u`` shape ``(6, m1, m2)`` ->
    ``(6, m1, m2)``. ``alphas`` shape ``(8, m1)`` (per theta1 row, broadcast
    across theta2). Matches :func:`cyclerfinder.core.qbcp.qbcp_eom` exactly,
    including the Sun-only 1/alpha6 convention.
    """
    x, y, z, px, py, pz = (u[i] for i in range(6))
    a1, a2, a3, a4, a5, a6, xs, ys = (alphas[i][:, None] for i in range(8))
    dx = a1 * px + a2 * x + a3 * y
    dy = a1 * py + a2 * y - a3 * x
    dz = a1 * pz + a2 * z
    rpe2 = (x + mu) ** 2 + y * y + z * z
    rpm2 = (x - 1.0 + mu) ** 2 + y * y + z * z
    rps2 = (x - xs) ** 2 + (y - ys) ** 2 + z * z
    rpe3 = rpe2 * np.sqrt(rpe2)
    rpm3 = rpm2 * np.sqrt(rpm2)
    rps3 = rps2 * np.sqrt(rps2)
    msun_a6 = mu_sun / a6
    pot_x = (1.0 - mu) * (x + mu) / rpe3 + mu * (x - 1.0 + mu) / rpm3 + msun_a6 * (x - xs) / rps3
    pot_y = (1.0 - mu) * y / rpe3 + mu * y / rpm3 + msun_a6 * (y - ys) / rps3
    pot_z = (1.0 - mu) * z / rpe3 + mu * z / rpm3 + msun_a6 * z / rps3
    dpx = -a2 * px + a3 * py - a4 - pot_x
    dpy = -a2 * py - a3 * px - a5 - pot_y
    dpz = -a2 * pz - pot_z
    return np.stack([dx, dy, dz, dpx, dpy, dpz], axis=0)


def _qbcp_jacobian_grid(
    u: NDArray[np.float64], alphas: NDArray[np.float64], mu: float, mu_sun: float
) -> NDArray[np.float64]:
    """Vectorized 6x6 state-Jacobian ``dF/du`` at every grid point.

    ``u`` shape ``(6, m1, m2)``, ``alphas`` shape ``(8, m1)``; returns
    ``(6, 6, m1, m2)``. Exactly ``qbcp_stm_eom``'s ``jac_a`` matrix (kinematic
    a1/a2/a3 blocks + potential-Hessian block from
    ``qbcp_potential_second_derivatives``), evaluated pointwise with the
    theta1-dependent alphas -- the analytic Jacobian of the PDE's ``F(u)`` term.
    """
    x, y, z = u[0], u[1], u[2]
    a1, a2, a3, _a4, _a5, a6, xs, ys = (alphas[i][:, None] for i in range(8))
    rpe2 = (x + mu) ** 2 + y * y + z * z
    rpm2 = (x - 1.0 + mu) ** 2 + y * y + z * z
    rps2 = (x - xs) ** 2 + (y - ys) ** 2 + z * z
    rpe3 = rpe2 * np.sqrt(rpe2)
    rpm3 = rpm2 * np.sqrt(rpm2)
    rps3 = rps2 * np.sqrt(rps2)
    rpe5 = rpe3 * rpe2
    rpm5 = rpm3 * rpm2
    rps5 = rps3 * rps2
    om1 = 1.0 - mu
    msun_a6 = mu_sun / a6
    uxx = (
        -om1 * (1.0 / rpe3 - 3.0 * (x + mu) ** 2 / rpe5)
        - mu * (1.0 / rpm3 - 3.0 * (x - 1.0 + mu) ** 2 / rpm5)
        - msun_a6 * (1.0 / rps3 - 3.0 * (x - xs) ** 2 / rps5)
    )
    uyy = (
        -om1 * (1.0 / rpe3 - 3.0 * y * y / rpe5)
        - mu * (1.0 / rpm3 - 3.0 * y * y / rpm5)
        - msun_a6 * (1.0 / rps3 - 3.0 * (y - ys) ** 2 / rps5)
    )
    uzz = (
        -om1 * (1.0 / rpe3 - 3.0 * z * z / rpe5)
        - mu * (1.0 / rpm3 - 3.0 * z * z / rpm5)
        - msun_a6 * (1.0 / rps3 - 3.0 * z * z / rps5)
    )
    uxy = (
        3.0 * om1 * (x + mu) * y / rpe5
        + 3.0 * mu * (x - 1.0 + mu) * y / rpm5
        + 3.0 * msun_a6 * (x - xs) * (y - ys) / rps5
    )
    uxz = (
        3.0 * om1 * (x + mu) * z / rpe5
        + 3.0 * mu * (x - 1.0 + mu) * z / rpm5
        + 3.0 * msun_a6 * (x - xs) * z / rps5
    )
    uyz = 3.0 * om1 * y * z / rpe5 + 3.0 * mu * y * z / rpm5 + 3.0 * msun_a6 * (y - ys) * z / rps5

    m1, m2 = x.shape
    jf = np.zeros((6, 6, m1, m2))
    ones = np.ones((m1, m2))
    a1b, a2b, a3b = a1 * ones, a2 * ones, a3 * ones
    jf[0, 0], jf[0, 1], jf[0, 3] = a2b, a3b, a1b
    jf[1, 0], jf[1, 1], jf[1, 4] = -a3b, a2b, a1b
    jf[2, 2], jf[2, 5] = a2b, a1b
    jf[3, 0], jf[3, 1], jf[3, 2], jf[3, 3], jf[3, 4] = uxx, uxy, uxz, -a2b, a3b
    jf[4, 0], jf[4, 1], jf[4, 2], jf[4, 3], jf[4, 4] = uxy, uyy, uyz, -a3b, -a2b
    jf[5, 0], jf[5, 1], jf[5, 2], jf[5, 5] = uxz, uyz, uzz, -a2b
    return jf


def _transverse_amplitude(coeffs: NDArray[np.float64], n2: int) -> float:
    cos1, sin1 = _k2_first_harmonic_cols(n2)
    return float(np.sqrt(np.sum(coeffs[:, :, cos1] ** 2 + coeffs[:, :, sin1] ** 2)))


# ---------------------------------------------------------------------------
# Invariance-PDE residual and its analytic Jacobian (no integration).
# ---------------------------------------------------------------------------


def _residual(
    z: NDArray[np.float64],
    omega1: float,
    mu: float,
    mu_sun: float,
    n1: int,
    n2: int,
    p1: NDArray[np.float64],
    p1d: NDArray[np.float64],
    p2: NDArray[np.float64],
    p2d: NDArray[np.float64],
    alphas: NDArray[np.float64],
    phase2_coord: int,
    amplitude_anchor: float,
    gauge_weight: float,
    rho_target: float,
    rho_weight: float,
) -> NDArray[np.float64]:
    """Real residual: QBCP quasi-periodic invariance PDE on the 2D grid + 3 gauges.

    Gauge rows (see module docstring): transverse phase, transverse amplitude
    anchor, rotation-number pin ``omega2 = rho_target * omega1``. NO longitudinal
    phase gauge (theta1 is locked to the Sun epoch, no phase freedom).
    """
    coeffs, omega2 = _unpack(z, n1, n2)
    if not np.isfinite(omega2):
        return np.full(p1.shape[0] * p2.shape[0] * _N_STATE + 3, 1e6)

    m1, m2 = p1.shape[0], p2.shape[0]
    u = np.empty((_N_STATE, m1, m2))
    du1 = np.empty((_N_STATE, m1, m2))
    du2 = np.empty((_N_STATE, m1, m2))
    for c in range(_N_STATE):
        cc = coeffs[c]
        u[c] = p1 @ cc @ p2.T
        du1[c] = p1d @ cc @ p2.T
        du2[c] = p1 @ cc @ p2d.T

    rhs = _qbcp_rhs_grid(u, alphas, mu, mu_sun)
    lhs = omega1 * du1 + omega2 * du2
    pde_res = (lhs - rhs).reshape(-1)

    _, sin1_col_t2 = _k2_first_harmonic_cols(n2)
    phase2 = coeffs[phase2_coord, 0, sin1_col_t2]  # const(theta1) x sin(theta2)
    amp = _transverse_amplitude(coeffs, n2)
    gauge = np.array(
        [
            gauge_weight * phase2,
            gauge_weight * (amp - amplitude_anchor),
            rho_weight * (omega2 - rho_target * omega1),
        ]
    )
    return np.concatenate([pde_res, gauge])


def _jacobian(
    z: NDArray[np.float64],
    omega1: float,
    mu: float,
    mu_sun: float,
    n1: int,
    n2: int,
    p1: NDArray[np.float64],
    p1d: NDArray[np.float64],
    p2: NDArray[np.float64],
    p2d: NDArray[np.float64],
    alphas: NDArray[np.float64],
    phase2_coord: int,
    amplitude_anchor: float,
    gauge_weight: float,
    rho_target: float,
    rho_weight: float,
) -> NDArray[np.float64]:
    """Analytic Jacobian of :func:`_residual` (exact; no finite differencing).

    Row block: PDE residual (6*m1*m2) then 3 gauge rows. Columns: vec(coeffs)
    (6*A1*A2) then [omega2]. Per module docstring:
    ``dR_c/dC_d = delta_cd (w1 P1d P2 + w2 P1 P2d) - Jf_cd (P1 P2)``,
    ``dR_c/domega2 = du2_c``.
    """
    coeffs, omega2 = _unpack(z, n1, n2)
    a1n, a2n = 2 * n1 + 1, 2 * n2 + 1
    m1, m2 = p1.shape[0], p2.shape[0]
    n_coef = _N_STATE * a1n * a2n
    n_free = n_coef + 1
    n_pde = _N_STATE * m1 * m2
    jac = np.zeros((n_pde + 3, n_free))

    b_val = np.einsum("ia,jb->ijab", p1, p2)
    b_d1 = np.einsum("ia,jb->ijab", p1d, p2)
    b_d2 = np.einsum("ia,jb->ijab", p1, p2d)
    lhs_basis = (omega1 * b_d1 + omega2 * b_d2).reshape(m1 * m2, a1n * a2n)
    val_basis = b_val.reshape(m1 * m2, a1n * a2n)

    u = np.empty((_N_STATE, m1, m2))
    du2 = np.empty((_N_STATE, m1, m2))
    for c in range(_N_STATE):
        u[c] = p1 @ coeffs[c] @ p2.T
        du2[c] = p1 @ coeffs[c] @ p2d.T
    jf = _qbcp_jacobian_grid(u, alphas, mu, mu_sun)  # (6,6,m1,m2)

    for c in range(_N_STATE):
        rows = slice(c * m1 * m2, (c + 1) * m1 * m2)
        cols_c = slice(c * a1n * a2n, (c + 1) * a1n * a2n)
        jac[rows, cols_c] = lhs_basis
        for d in range(_N_STATE):
            cols_d = slice(d * a1n * a2n, (d + 1) * a1n * a2n)
            jf_cd = jf[c, d].reshape(m1 * m2)
            jac[rows, cols_d] -= jf_cd[:, None] * val_basis
        jac[rows, n_coef] = du2[c].reshape(m1 * m2)

    _, sin1_col_t2 = _k2_first_harmonic_cols(n2)
    cos1_col_t2, _ = _k2_first_harmonic_cols(n2)

    def _coef_index(c: int, a: int, b: int) -> int:
        return c * a1n * a2n + a * a2n + b

    jac[n_pde + 0, _coef_index(phase2_coord, 0, sin1_col_t2)] = gauge_weight
    amp = _transverse_amplitude(coeffs, n2)
    if amp > 0:
        for c in range(_N_STATE):
            for a in range(a1n):
                jac[n_pde + 1, _coef_index(c, a, cos1_col_t2)] = (
                    gauge_weight * coeffs[c, a, cos1_col_t2] / amp
                )
                jac[n_pde + 1, _coef_index(c, a, sin1_col_t2)] = (
                    gauge_weight * coeffs[c, a, sin1_col_t2] / amp
                )
    jac[n_pde + 2, n_coef] = rho_weight
    return jac


# ---------------------------------------------------------------------------
# Bootstrap from a converged GMOS QBCP torus.
# ---------------------------------------------------------------------------


def project_gmos_qbcp_torus_to_2d(
    torus: QBCPTorus, n1: int, n2: int, m1: int, m2: int
) -> NDArray[np.float64]:
    """Least-squares project a converged GMOS :class:`QBCPTorus` onto the 2D real
    Fourier basis -> initial ``(6, 2*n1+1, 2*n2+1)`` PM coefficient array.

    Samples the GMOS torus on the ``(theta1, theta2)`` grid via
    ``evaluate_qbcp_torus`` (which propagates -- allowed, this is only the
    guess; ``theta_long`` is our ``theta1``, ``theta_trans`` our ``theta2``),
    converts each PV sample to PM at ``t = theta1 / omega_long``, then solves the
    two 1D least-squares fits ``U = P1 C P2^T`` for ``C``.
    """
    omega1 = torus.omega_long
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, _ = _basis_matrices(n1, t1)
    p2, _ = _basis_matrices(n2, t2)
    p1_pinv = np.linalg.pinv(p1)
    p2_pinv = np.linalg.pinv(p2)
    u_grid = np.empty((_N_STATE, m1, m2))
    for i, a in enumerate(t1):
        for j, b in enumerate(t2):
            pv = evaluate_qbcp_torus(torus, float(a), float(b))
            pm = qbcp.state_pv_to_pm(pv, float(a) / omega1, torus.system)
            u_grid[:, i, j] = pm
    coeffs = np.empty((_N_STATE, 2 * n1 + 1, 2 * n2 + 1))
    for c in range(_N_STATE):
        coeffs[c] = p1_pinv @ u_grid[c] @ p2_pinv.T
    return coeffs


# ---------------------------------------------------------------------------
# The corrector.
# ---------------------------------------------------------------------------


def correct_qbcp_torus_pseudospectral(
    system: qbcp.QBCPSystem,
    coeffs0: NDArray[np.float64],
    omega2_0: float,
    *,
    n1: int,
    n2: int,
    amplitude_anchor: float,
    period_multiple: int = 1,
    rho_target: float | None = None,
    rho_weight: float = 1.0,
    m1: int | None = None,
    m2: int | None = None,
    tol: float = 1e-6,
    closure_tol: float = 1e-4,
    max_nfev: int = 200,
    gauge_weight: float = 1.0,
    tr_solver: str = "exact",
    closure_dt_frac: float = 0.02,
    n_closure_samples: int = 12,
    notes: str = "",
) -> QBCPTorusVariationalResult:
    """Refine a 2D-Fourier PM torus guess by driving the QBCP invariance-PDE
    residual to zero -- no forward integration in the search (see module
    docstring).

    ``coeffs0`` shape ``(6, 2*n1+1, 2*n2+1)``; ``omega2_0`` the initial
    transverse frequency. ``omega1`` is fixed internally at ``period_multiple *
    system.omega_sun_nondim``. ``amplitude_anchor`` fixes the transverse (k2=1)
    coefficient L2 norm. Collocation grid defaults to ``m1 = 2*n1+3``, ``m2 =
    2*n2+3`` (mild oversampling). Uses ``scipy.optimize.least_squares`` with
    ``method="trf"`` and an exact analytic Jacobian (``method="trf"`` respects
    the ``nfev`` budget far better than ``"lm"`` on ill-conditioned problems --
    the `#611` lesson).
    """
    if coeffs0.shape != (_N_STATE, 2 * n1 + 1, 2 * n2 + 1):
        raise ValueError(
            f"coeffs0 shape {coeffs0.shape} != expected {(_N_STATE, 2 * n1 + 1, 2 * n2 + 1)}"
        )
    if period_multiple < 1:
        raise ValueError(f"period_multiple must be >= 1, got {period_multiple}")
    mu = system.mu
    mu_sun = system.mu_sun
    omega1 = period_multiple * system.omega_sun_nondim
    period = 2.0 * np.pi / omega1
    if rho_target is None:
        rho_target = omega2_0 / omega1
    m1 = m1 if m1 is not None else 2 * n1 + 3
    m2 = m2 if m2 is not None else 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    alphas = _alphas_on_theta1(t1, omega1, system)

    cos1, sin1 = _k2_first_harmonic_cols(n2)
    trans_by_coord = np.sqrt(np.sum(coeffs0[:, :, cos1] ** 2 + coeffs0[:, :, sin1] ** 2, axis=1))
    phase2_coord = int(np.argmax(trans_by_coord))

    z0 = _pack(coeffs0, omega2_0)
    solver_args = (
        omega1,
        mu,
        mu_sun,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        alphas,
        phase2_coord,
        amplitude_anchor,
        gauge_weight,
        rho_target,
        rho_weight,
    )
    sol = least_squares(
        _residual,
        z0,
        jac=_jacobian,
        args=solver_args,
        method="trf",
        tr_solver=tr_solver,
        x_scale="jac",
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=max_nfev,
    )
    coeffs, omega2 = _unpack(sol.x, n1, n2)
    n_res = m1 * m2 * _N_STATE
    residual_rms = float(np.sqrt(np.sum(sol.fun[:n_res] ** 2) / n_res))
    amp = _transverse_amplitude(coeffs, n2)
    rot = omega2 / omega1

    result = QBCPTorusVariationalResult(
        system=system,
        coeffs=coeffs,
        omega1=omega1,
        omega2=omega2,
        rotation_number=rot,
        rho_strob=omega2 * period,
        period=period,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        period_multiple=period_multiple,
        transverse_amplitude=amp,
        residual_rms=residual_rms,
        closure_residual=float("nan"),
        converged=False,
        n_iter=int(sol.nfev),
        notes=notes,
    )
    closure = _independent_closure(result, closure_dt_frac, n_closure_samples)
    converged = (residual_rms < tol) and (closure < closure_tol)
    result = QBCPTorusVariationalResult(
        system=system,
        coeffs=coeffs,
        omega1=omega1,
        omega2=omega2,
        rotation_number=rot,
        rho_strob=omega2 * period,
        period=period,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        period_multiple=period_multiple,
        transverse_amplitude=amp,
        residual_rms=residual_rms,
        closure_residual=closure,
        converged=converged,
        n_iter=int(sol.nfev),
        notes=notes,
    )
    log_outcome(
        solver="variational_qbcp_torus.correct_qbcp_torus_pseudospectral",
        inputs={
            "mu": float(mu),
            "mu_sun": float(mu_sun),
            "n1": int(n1),
            "n2": int(n2),
            "period_multiple": int(period_multiple),
            "amplitude_anchor": float(amplitude_anchor),
        },
        outcome={
            "converged": bool(converged),
            "residual_rms": residual_rms,
            "closure_residual": closure,
            "omega2": float(omega2),
            "rotation_number": float(rot),
            "rho_strob": float(omega2 * period),
            "transverse_amplitude": float(amp),
        },
        meta={"model": "qbcp"},
    )
    return result


def _independent_closure(
    result: QBCPTorusVariationalResult, dt_frac: float, n_samples: int
) -> float:
    """Short-time flow-consistency check (independent of the algebraic residual).

    For several base angles, propagate the PM state ``u(theta1, theta2)`` (at
    absolute time ``t0 = theta1/omega1``) for ``dt = dt_frac * period`` through
    the true nonlinear ``qbcp_eom`` and compare to ``u(theta1 + omega1*dt,
    theta2 + omega2*dt)``. Short ``dt`` keeps the unstable-flow amplification
    bounded (see module docstring). Returns the max over samples.
    """
    dt = dt_frac * result.period
    rng = np.random.default_rng(0xC0FFEE)
    max_err = 0.0
    for _ in range(n_samples):
        th1 = float(rng.uniform(0, 2 * np.pi))
        th2 = float(rng.uniform(0, 2 * np.pi))
        u0 = evaluate_torus_state(result, th1, th2)
        t0 = th1 / result.omega1
        sol = solve_ivp(
            lambda t, y: qbcp.qbcp_eom(t, y, result.system),
            (t0, t0 + dt),
            u0,
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
        )
        if not sol.success:
            return float("inf")
        u_target = evaluate_torus_state(result, th1 + result.omega1 * dt, th2 + result.omega2 * dt)
        max_err = max(max_err, float(np.linalg.norm(sol.y[:, -1] - u_target)))
    return max_err


def discover_qbcp_torus_from_gmos(
    system: qbcp.QBCPSystem,
    seed_gmos_torus: QBCPTorus,
    *,
    n1: int = 8,
    n2: int = 3,
    m1: int | None = None,
    m2: int | None = None,
    rho_target: float | None = None,
    rho_weight: float = 1.0,
    gauge_weight: float = 1.0,
    tr_solver: str = "exact",
    tol: float = 1e-6,
    closure_tol: float = 1e-4,
    max_nfev: int = 200,
    closure_dt_frac: float = 0.02,
    notes: str = "",
) -> QBCPTorusVariationalResult:
    """Build a 2D pseudospectral QBCP torus from a converged GMOS seed.

    Projects ``seed_gmos_torus`` onto the 2D PM Fourier basis (bootstrap; uses
    propagation for the guess only), sets ``omega2 = omega_trans``, anchors the
    transverse amplitude at the seed's own value, and refines via
    :func:`correct_qbcp_torus_pseudospectral` (no integration in the solve). This
    is the entry point for the SE-L2 positive control.
    """
    m1 = m1 if m1 is not None else 2 * n1 + 3
    m2 = m2 if m2 is not None else 2 * n2 + 3
    coeffs0 = project_gmos_qbcp_torus_to_2d(seed_gmos_torus, n1, n2, m1, m2)
    amp_anchor = _transverse_amplitude(coeffs0, n2)
    return correct_qbcp_torus_pseudospectral(
        system,
        coeffs0,
        seed_gmos_torus.omega_trans,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        amplitude_anchor=amp_anchor,
        rho_target=rho_target,
        rho_weight=rho_weight,
        gauge_weight=gauge_weight,
        tr_solver=tr_solver,
        tol=tol,
        closure_tol=closure_tol,
        max_nfev=max_nfev,
        closure_dt_frac=closure_dt_frac,
        notes=notes or "bootstrapped_from_gmos",
    )


def continue_qbcp_torus_amplitude(
    result: QBCPTorusVariationalResult,
    target_amplitude: float,
    *,
    n_steps: int = 8,
    rho_weight: float = 1.0,
    gauge_weight: float = 1.0,
    tr_solver: str = "exact",
    tol: float = 1e-6,
    closure_tol: float = 1e-4,
    max_nfev: int = 400,
    closure_dt_frac: float = 0.02,
) -> list[QBCPTorusVariationalResult]:
    """Natural-parameter continuation in transverse amplitude, using ONLY this
    module's own integration-free solves (`#612`'s pattern).

    Steps the amplitude anchor linearly from ``result.transverse_amplitude`` to
    ``target_amplitude`` in ``n_steps`` steps, warm-starting each solve from the
    previous torus and holding the rotation number pinned at the running
    value. Returns the list of per-step results (the last is the
    target-amplitude torus). Stops early if a step fails to converge.
    """
    steps: list[QBCPTorusVariationalResult] = []
    current = result
    amps = np.linspace(current.transverse_amplitude, target_amplitude, n_steps + 1)[1:]
    for amp in amps:
        nxt = correct_qbcp_torus_pseudospectral(
            current.system,
            current.coeffs,
            current.omega2,
            n1=current.n1,
            n2=current.n2,
            m1=current.m1,
            m2=current.m2,
            period_multiple=current.period_multiple,
            amplitude_anchor=float(amp),
            rho_target=current.rotation_number,
            rho_weight=rho_weight,
            gauge_weight=gauge_weight,
            tr_solver=tr_solver,
            tol=tol,
            closure_tol=closure_tol,
            max_nfev=max_nfev,
            closure_dt_frac=closure_dt_frac,
            notes=f"continuation_amp={amp:.5f}",
        )
        steps.append(nxt)
        current = nxt
        if not nxt.converged:
            break
    return steps
