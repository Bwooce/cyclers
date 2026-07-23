"""Pseudospectral quasi-periodic-torus corrector for the CCR4BP (#690).

The "EOM swap" of `#617`'s QBCP torus corrector
-----------------------------------------------
This module is `#617`'s seedless 2D pseudospectral quasi-periodic-torus
corrector (``search.variational_qbcp_torus``) with ONE thing changed: the
underlying equations of motion + state Jacobian are `#689`'s
:mod:`cyclerfinder.core.ccr4bp` (planar Jupiter-Europa-Ganymede CCR4BP) instead
of the QBCP's. `#686`'s strategy note (section 3, Stage B item 2) flagged this as
the legitimate "EOM swap" case precisely because the two problems have the
IDENTICAL shape:

* both are NON-AUTONOMOUS, time-periodic, forced by a single external clock;
* one torus angle ``theta1`` is LOCKED to that clock -- for the QBCP the Sun's
  synodic phase, here Ganymede's synodic phase ``theta_gan(t)`` -- with a FIXED
  frequency ``omega1`` (not a solved-for unknown);
* the other torus angle ``theta2`` is the genuine internal (orbit / libration)
  angle, with a FREE frequency ``omega2``; the rotation number
  ``rho = omega2 / omega1`` selects the family member.

So the whole pseudospectral machinery -- the real tensor-product Fourier
discretization, the quasi-periodic invariance-PDE residual on a 2D collocation
grid, the exact analytic Jacobian, the three gauge rows (transverse phase,
transverse amplitude anchor, rotation-number pin), the ``least_squares`` /
``trf`` Newton solve, and the independent short-time-flow closure check --
transfers VERBATIM from `#617`. What genuinely changed (documented so the swap is
auditable):

1. **State dimension 6 -> 4.** The QBCP corrector carries the 6-canonical
   (PM) state ``(x, y, z, px, py, pz)``. The CCR4BP is PLANAR and worked in the
   ordinary rotating-frame velocity state ``(x, y, vx, vy)`` (``_N_STATE = 4``).
   No position/momentum (PV<->PM) conversion is needed anywhere -- the CCR4BP RHS
   is already in the state the corrector stores, which is strictly SIMPLER than
   the QBCP's canonical bookkeeping.
2. **RHS + state-Jacobian backend.** :func:`_ccr4bp_rhs_grid` is the vectorized
   :func:`cyclerfinder.core.ccr4bp.ccr4bp_eom` (planar), replacing the QBCP's
   ``_qbcp_rhs_grid``; :func:`_ccr4bp_jacobian_grid` is the vectorized 4x4
   Jacobian block of :func:`cyclerfinder.core.ccr4bp.ccr4bp_stm_eom`, replacing
   ``_qbcp_jacobian_grid``. Both are checked pointwise to machine precision
   against the scalar ``core.ccr4bp`` functions in the test module -- the same
   discipline `#617` used against ``qbcp_eom`` / ``qbcp_stm_eom``.
3. **Forcing clock = Ganymede's synodic angle.** The QBCP evaluated 8 ``alpha``
   forcing functions at each ``theta1``; here :func:`_ganymede_on_theta1`
   evaluates Ganymede's synodic-frame position ``(gx, gy)`` at each ``theta1``,
   the only time-dependent input the CCR4BP RHS has. ``omega1 = period_multiple *
   |omega_gan|`` (POSITIVE; Ganymede regresses in the Europa-synodic frame so
   ``omega_gan < 0``, and ``theta1 = omega1 * t`` runs forward with the Ganymede
   phase mapped as ``theta_gan = theta_gan0 + omega_gan * (theta1/omega1) =
   theta_gan0 - theta1/period_multiple``). This mirrors the QBCP's
   ``omega1 = omega_s`` exactly in structure, only the sign convention of the
   forcing rate differs -- handled once, in :func:`_ganymede_on_theta1`.

Everything else -- the ``theta1``-locked / ``theta2``-free asymmetry, dropping
`#612`'s longitudinal phase gauge (``theta1`` has no phase freedom, its origin is
pinned by the Ganymede epoch), the three surviving gauge rows, and the
integration-free residual whose Jacobian never sees the unstable flow's
per-period amplification -- is `#617`'s reasoning unchanged. See that module's
docstring for the full derivation; only the deltas above are repeated here.

The physical target (positive control)
--------------------------------------
Kumar-Anderson-de la Llave-Gunter 2021 (arXiv:2109.14815) computed CCR4BP
quasi-periodic analogues of the Jupiter-Europa 3:4 and Jupiter-Ganymede 3:2 / 7:5
resonant periodic orbits at PHYSICAL Galilean masses (Europa:Ganymede period
ratio 2.014, `#689`'s / this module's default -- NOT an exact-2:1 idealization).
Their object is exactly this corrector's object: an unstable resonant periodic
orbit of the base CRTBP that persists as a 2D invariant torus under the fourth
body's periodic forcing, with rotation number ``omega = 2*pi*Omega1/|Omega_p-1|``
in their notation (our ``rho_strob = omega2 * period``). The positive control
(:func:`discover_ccr4bp_torus_from_resonant_orbit` + mass continuation) targets
the Jupiter-Europa 3:4 member of that family. The published paper does NOT list
initial conditions / energy / rotation number to pixel-reproduce, so the control
is honest-partial in the `#664` sense: the ``mu_gan -> 0`` regression is a
machine-precision check (the corrector must recover the base resonant orbit as a
``theta1``-flat torus), and the physical-mass torus is reported at its achieved
invariance residual + rotation number, matching the object CLASS not a specific
published torus.

Discipline / scope (`#690`)
---------------------------
Build + validate ONLY the torus corrector and its positive control. Explicitly
NOT the whisker/manifold globalization or the mesh-intersection heteroclinic
search (later `#686` Stage-B sub-tasks). The returned torus is OUR computation;
no novelty claim, no catalogue writeback.

References
----------
* Kumar, B., Anderson, R. L., de la Llave, R., Gunter, B. (2021). AAS 21-651 /
  arXiv:2109.14815 (the CCR4BP resonant tori this reproduces).
* Olikara, Z., & Scheeres, D. (2010); Jorba, A. (2001) -- see
  ``search.variational_qbcp_torus`` for the full method lineage this reuses.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.ccr4bp as ccr4bp
from cyclerfinder.search.outcome_log import log_outcome

# Reuse the state-dimension-agnostic Fourier-basis machinery from the QBCP
# corrector VERBATIM (this is the "don't rebuild the pseudospectral machinery"
# point of the EOM swap): the real tensor-product basis and the k2=1 harmonic
# column indices are identical -- they depend only on the mode counts, not on the
# physics or the state dimension.
from cyclerfinder.search.variational_qbcp_torus import (
    _basis_matrices,
    _k2_first_harmonic_cols,
)

_N_STATE = 4  # planar rotating-frame state: x, y, vx, vy


@dataclass(frozen=True)
class CCR4BPTorusVariationalResult:
    """A CCR4BP quasi-periodic 2-torus found by 2D pseudospectral collocation.

    ``coeffs`` is the ``(4, 2*n1+1, 2*n2+1)`` real tensor-product Fourier
    coefficient array (planar rotating-frame state ``(x, y, vx, vy)``; basis as in
    :mod:`cyclerfinder.search.variational_qbcp_torus`). ``omega1`` is FIXED at
    ``period_multiple * |omega_gan|`` (Ganymede's synodic frequency -- NOT solved
    for); ``omega2`` is the free internal frequency (rad/TU);
    ``rotation_number = omega2 / omega1`` and ``rho_strob = omega2 * period`` is
    the advance of ``theta2`` per Ganymede synodic period (the stroboscopic-map
    rotation number). ``residual_rms`` is the RMS of the quasi-periodic
    invariance-PDE residual at the 2D collocation points (the quantity
    minimized). ``closure_residual`` is the INDEPENDENT short-time flow check
    (propagate ``u(theta1,theta2)`` for ``closure_dt`` through the true nonlinear
    ``ccr4bp_eom``, compare to ``u(theta1+omega1*dt, theta2+omega2*dt)`` -- NOT
    circular with ``residual_rms``). ``transverse_amplitude`` is the anchored L2
    norm of the ``k2=1`` coefficients.
    """

    system: ccr4bp.CCR4BPSystem
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
# Pack / unpack and torus evaluation (4-state analogues of the QBCP corrector).
# ---------------------------------------------------------------------------


def _pack(coeffs: NDArray[np.float64], omega2: float) -> NDArray[np.float64]:
    return np.concatenate([coeffs.reshape(-1), [omega2]])


def _unpack(z: NDArray[np.float64], n1: int, n2: int) -> tuple[NDArray[np.float64], float]:
    n_coef = _N_STATE * (2 * n1 + 1) * (2 * n2 + 1)
    coeffs = z[:n_coef].reshape(_N_STATE, 2 * n1 + 1, 2 * n2 + 1)
    return coeffs, float(z[-1])


def evaluate_torus_state(
    result: CCR4BPTorusVariationalResult,
    theta1: float | NDArray[np.float64],
    theta2: float | NDArray[np.float64],
) -> NDArray[np.float64]:
    """Evaluate the 2D Fourier torus state ``u(theta1, theta2)`` (4-vector).

    Scalars in -> ``(4,)`` out (the analytic Fourier series, no propagation).
    Identical algorithm to the QBCP corrector's ``evaluate_torus_state``.
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
        return np.asarray(out[0], dtype=np.float64)
    return out


# ---------------------------------------------------------------------------
# Vectorized CCR4BP RHS and its state-Jacobian on the (theta1, theta2) grid.
# ---------------------------------------------------------------------------


def _ganymede_on_theta1(
    theta1: NDArray[np.float64], omega1: float, system: ccr4bp.CCR4BPSystem
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Ganymede synodic-frame position ``(gx, gy)`` at each ``theta1`` grid value.

    ``theta1`` locks to the Ganymede synodic phase via ``t = theta1 / omega1``
    (``omega1 = period_multiple * |omega_gan| > 0``), so
    :func:`cyclerfinder.core.ccr4bp._ganymede_position` sees
    ``theta_gan = theta_gan0 + omega_gan * t = theta_gan0 - theta1/period_multiple``
    (``omega_gan < 0``). Returns two ``(m1,)`` arrays. This is the CCR4BP analogue
    of the QBCP corrector's ``_alphas_on_theta1``: the only time-dependent input
    the RHS has.
    """
    m1 = theta1.size
    gx = np.empty(m1, dtype=np.float64)
    gy = np.empty(m1, dtype=np.float64)
    for i in range(m1):
        t = float(theta1[i]) / omega1
        px, py, _ = ccr4bp._ganymede_position(t, system)
        gx[i] = px
        gy[i] = py
    return gx, gy


def _ccr4bp_rhs_grid(
    u: NDArray[np.float64],
    gx: NDArray[np.float64],
    gy: NDArray[np.float64],
    mu: float,
    mu_gan: float,
    a_gan: float,
) -> NDArray[np.float64]:
    """Vectorized planar ``ccr4bp_eom`` on the grid: ``u`` shape ``(4, m1, m2)`` ->
    ``(4, m1, m2)``. ``gx, gy`` shape ``(m1,)`` (Ganymede position per theta1 row,
    broadcast across theta2). Matches
    :func:`cyclerfinder.core.ccr4bp.ccr4bp_eom` (planar) exactly.
    """
    x, y, vx, vy = u[0], u[1], u[2], u[3]
    gxb, gyb = gx[:, None], gy[:, None]
    r1c = ((x + mu) ** 2 + y * y) ** 1.5
    r2c = ((x - 1.0 + mu) ** 2 + y * y) ** 1.5
    om1 = 1.0 - mu
    base_ax = x + 2.0 * vy - om1 * (x + mu) / r1c - mu * (x - 1.0 + mu) / r2c
    base_ay = y - 2.0 * vx - om1 * y / r1c - mu * y / r2c
    dxg = x - gxb
    dyg = y - gyb
    dg3 = (dxg * dxg + dyg * dyg) ** 1.5
    a_gan3 = a_gan**3
    gax = -mu_gan * dxg / dg3 - mu_gan * gxb / a_gan3
    gay = -mu_gan * dyg / dg3 - mu_gan * gyb / a_gan3
    return np.stack([vx, vy, base_ax + gax, base_ay + gay], axis=0)


def _ccr4bp_jacobian_grid(
    u: NDArray[np.float64],
    gx: NDArray[np.float64],
    gy: NDArray[np.float64],
    mu: float,
    mu_gan: float,
    a_gan: float,
) -> NDArray[np.float64]:
    """Vectorized 4x4 state-Jacobian ``dF/du`` at every grid point.

    ``u`` shape ``(4, m1, m2)``; returns ``(4, 4, m1, m2)``. Exactly the planar
    4x4 block of :func:`cyclerfinder.core.ccr4bp.ccr4bp_stm_eom`'s A matrix
    (CR3BP pseudo-potential second derivatives + Ganymede direct-term Hessian +
    Coriolis), the analytic Jacobian of the PDE's ``F(u)`` term. The indirect
    Ganymede term is state-independent, so it drops here (as it does in
    ``core.ccr4bp``).
    """
    x, y = u[0], u[1]
    gxb, gyb = gx[:, None], gy[:, None]
    r1c = ((x + mu) ** 2 + y * y) ** 1.5
    r2c = ((x - 1.0 + mu) ** 2 + y * y) ** 1.5
    r1f = ((x + mu) ** 2 + y * y) ** 2.5
    r2f = ((x - 1.0 + mu) ** 2 + y * y) ** 2.5
    om1 = 1.0 - mu
    uxx = (
        1.0
        - om1 / r1c
        - mu / r2c
        + 3.0 * om1 * (x + mu) ** 2 / r1f
        + 3.0 * mu * (x - 1.0 + mu) ** 2 / r2f
    )
    uyy = 1.0 - om1 / r1c - mu / r2c + 3.0 * om1 * y * y / r1f + 3.0 * mu * y * y / r2f
    uxy = 3.0 * om1 * (x + mu) * y / r1f + 3.0 * mu * (x - 1.0 + mu) * y / r2f

    dxg = x - gxb
    dyg = y - gyb
    dg2 = dxg * dxg + dyg * dyg
    dg3 = dg2**1.5
    dg5 = dg2**2.5
    gxx = -mu_gan * (1.0 / dg3 - 3.0 * dxg * dxg / dg5)
    gyy = -mu_gan * (1.0 / dg3 - 3.0 * dyg * dyg / dg5)
    gxy = 3.0 * mu_gan * dxg * dyg / dg5

    m1, m2 = x.shape
    jf = np.zeros((_N_STATE, _N_STATE, m1, m2))
    ones = np.ones((m1, m2))
    jf[0, 2] = ones
    jf[1, 3] = ones
    jf[2, 0], jf[2, 1], jf[2, 3] = uxx + gxx, uxy + gxy, 2.0 * ones
    jf[3, 0], jf[3, 1], jf[3, 2] = uxy + gxy, uyy + gyy, -2.0 * ones
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
    mu_gan: float,
    a_gan: float,
    n1: int,
    n2: int,
    p1: NDArray[np.float64],
    p1d: NDArray[np.float64],
    p2: NDArray[np.float64],
    p2d: NDArray[np.float64],
    gx: NDArray[np.float64],
    gy: NDArray[np.float64],
    phase2_coord: int,
    amplitude_anchor: float,
    gauge_weight: float,
    rho_target: float,
    rho_weight: float,
) -> NDArray[np.float64]:
    """Real residual: CCR4BP quasi-periodic invariance PDE on the 2D grid + 3 gauges.

    ``omega1*du/dtheta1 + omega2*du/dtheta2 - F(theta1, u) = 0`` at every
    collocation point, then the transverse phase gauge, the transverse amplitude
    anchor, and the rotation-number pin ``omega2 = rho_target*omega1``. No
    longitudinal phase gauge (``theta1`` is locked to the Ganymede epoch).
    Structurally identical to the QBCP corrector's ``_residual``.
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

    rhs = _ccr4bp_rhs_grid(u, gx, gy, mu, mu_gan, a_gan)
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
    mu_gan: float,
    a_gan: float,
    n1: int,
    n2: int,
    p1: NDArray[np.float64],
    p1d: NDArray[np.float64],
    p2: NDArray[np.float64],
    p2d: NDArray[np.float64],
    gx: NDArray[np.float64],
    gy: NDArray[np.float64],
    phase2_coord: int,
    amplitude_anchor: float,
    gauge_weight: float,
    rho_target: float,
    rho_weight: float,
) -> NDArray[np.float64]:
    """Analytic Jacobian of :func:`_residual` (exact; no finite differencing).

    Same block structure as the QBCP corrector's ``_jacobian``:
    ``dR_c/dC_d = delta_cd (w1 P1d P2 + w2 P1 P2d) - Jf_cd (P1 P2)``,
    ``dR_c/domega2 = du2_c``, plus the three gauge rows.
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
    jf = _ccr4bp_jacobian_grid(u, gx, gy, mu, mu_gan, a_gan)  # (4,4,m1,m2)

    for c in range(_N_STATE):
        rows = slice(c * m1 * m2, (c + 1) * m1 * m2)
        cols_c = slice(c * a1n * a2n, (c + 1) * a1n * a2n)
        jac[rows, cols_c] = lhs_basis
        for d in range(_N_STATE):
            cols_d = slice(d * a1n * a2n, (d + 1) * a1n * a2n)
            jf_cd = jf[c, d].reshape(m1 * m2)
            jac[rows, cols_d] -= jf_cd[:, None] * val_basis
        jac[rows, n_coef] = du2[c].reshape(m1 * m2)

    cos1_col_t2, sin1_col_t2 = _k2_first_harmonic_cols(n2)

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
# Bootstrap from a base-CRTBP resonant periodic orbit (theta1-flat torus).
# ---------------------------------------------------------------------------


def project_base_orbit_to_2d(
    orbit_state0: NDArray[np.float64],
    orbit_period: float,
    mu: float,
    n1: int,
    n2: int,
    m1: int,
    m2: int,
) -> NDArray[np.float64]:
    """Project a base-CR3BP (planar) periodic orbit onto a ``theta1``-FLAT 2D
    Fourier coefficient array ``(4, 2*n1+1, 2*n2+1)``.

    The orbit is sampled over one period as ``theta2`` in ``[0, 2*pi)`` (mapping
    ``t = theta2 * orbit_period / (2*pi)``, propagated through the bare CR3BP,
    ``mu_gan = 0``), then least-squares fit onto the ``theta2`` basis; all
    ``theta1`` structure is zero (no Ganymede forcing yet -- the corrector builds
    it). This is the CCR4BP analogue of the QBCP corrector's
    ``project_gmos_qbcp_torus_to_2d`` bootstrap (integration is allowed in the
    GUESS only), and the direct analogue of `#618`'s "seed the bare CR3BP orbit,
    let the corrector build the forcing structure" pattern.
    """
    t2 = 2.0 * np.pi * np.arange(m2) / m2
    p2, _ = _basis_matrices(n2, t2)
    p2_pinv = np.linalg.pinv(p2)
    base = ccr4bp.CCR4BPSystem(mu=mu, mu_gan=0.0, a_gan=2.0, omega_gan=-0.5)
    # Sample the orbit: state (x,y,z,vx,vy,vz); keep planar components (x,y,vx,vy).
    u_row = np.empty((_N_STATE, m2))
    for j, b in enumerate(t2):
        tj = float(b) * orbit_period / (2.0 * np.pi)
        if tj == 0.0:
            s = np.asarray(orbit_state0, dtype=np.float64)
        else:
            arc = ccr4bp.propagate_ccr4bp(base, np.asarray(orbit_state0, dtype=np.float64), tj)
            s = arc.state_f
        u_row[0, j], u_row[1, j] = s[0], s[1]
        u_row[2, j], u_row[3, j] = s[3], s[4]
    row0 = np.empty((_N_STATE, 2 * n2 + 1))
    for c in range(_N_STATE):
        row0[c] = p2_pinv @ u_row[c]
    coeffs = np.zeros((_N_STATE, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[:, 0, :] = row0  # theta1-constant row only (flat in theta1)
    return coeffs


# ---------------------------------------------------------------------------
# The corrector.
# ---------------------------------------------------------------------------


def correct_ccr4bp_torus_pseudospectral(
    system: ccr4bp.CCR4BPSystem,
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
    tr_solver: Literal["exact", "lsmr"] = "exact",
    closure_dt_frac: float = 0.02,
    n_closure_samples: int = 12,
    notes: str = "",
) -> CCR4BPTorusVariationalResult:
    """Refine a 2D-Fourier planar torus guess by driving the CCR4BP invariance-PDE
    residual to zero -- no forward integration in the search.

    ``coeffs0`` shape ``(4, 2*n1+1, 2*n2+1)``; ``omega2_0`` the initial internal
    frequency. ``omega1`` is fixed internally at
    ``period_multiple * |system.omega_gan|``. ``amplitude_anchor`` fixes the
    ``k2=1`` coefficient L2 norm. Everything below is the QBCP corrector's
    ``correct_qbcp_torus_pseudospectral`` with the CCR4BP RHS/Jacobian/forcing
    substituted.
    """
    if coeffs0.shape != (_N_STATE, 2 * n1 + 1, 2 * n2 + 1):
        raise ValueError(
            f"coeffs0 shape {coeffs0.shape} != expected {(_N_STATE, 2 * n1 + 1, 2 * n2 + 1)}"
        )
    if period_multiple < 1:
        raise ValueError(f"period_multiple must be >= 1, got {period_multiple}")
    mu = system.mu
    mu_gan = system.mu_gan
    a_gan = system.a_gan
    omega1 = period_multiple * abs(system.omega_gan)
    period = 2.0 * np.pi / omega1
    if rho_target is None:
        rho_target = omega2_0 / omega1
    m1 = m1 if m1 is not None else 2 * n1 + 3
    m2 = m2 if m2 is not None else 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    gx, gy = _ganymede_on_theta1(t1, omega1, system)

    cos1, sin1 = _k2_first_harmonic_cols(n2)
    trans_by_coord = np.sqrt(np.sum(coeffs0[:, :, cos1] ** 2 + coeffs0[:, :, sin1] ** 2, axis=1))
    phase2_coord = int(np.argmax(trans_by_coord))

    z0 = _pack(coeffs0, omega2_0)
    solver_args = (
        omega1,
        mu,
        mu_gan,
        a_gan,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        gx,
        gy,
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

    result = CCR4BPTorusVariationalResult(
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
    result = replace(result, closure_residual=closure, converged=converged)
    log_outcome(
        solver="variational_ccr4bp_torus.correct_ccr4bp_torus_pseudospectral",
        inputs={
            "mu": float(mu),
            "mu_gan": float(mu_gan),
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
        meta={"model": "ccr4bp"},
    )
    return result


def _independent_closure(
    result: CCR4BPTorusVariationalResult, dt_frac: float, n_samples: int
) -> float:
    """Short-time flow-consistency check (independent of the algebraic residual).

    For several base angles, propagate the state ``u(theta1, theta2)`` (at
    absolute time ``t0 = theta1/omega1``) for ``dt = dt_frac * period`` through
    the true nonlinear ``ccr4bp_eom`` and compare to ``u(theta1 + omega1*dt,
    theta2 + omega2*dt)``. The planar 4-state is embedded into the 6-state
    ``(x, y, 0, vx, vy, 0)`` for propagation. Returns the max over samples.
    """
    dt = dt_frac * result.period
    rng = np.random.default_rng(0xC0FFEE)
    max_err = 0.0
    for _ in range(n_samples):
        th1 = float(rng.uniform(0, 2 * np.pi))
        th2 = float(rng.uniform(0, 2 * np.pi))
        u0 = evaluate_torus_state(result, th1, th2)
        s6 = np.array([u0[0], u0[1], 0.0, u0[2], u0[3], 0.0], dtype=np.float64)
        t0 = th1 / result.omega1
        sol = solve_ivp(
            ccr4bp.ccr4bp_eom,
            (t0, t0 + dt),
            s6,
            args=(result.system,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
        )
        if not sol.success:
            return float("inf")
        u_target = evaluate_torus_state(result, th1 + result.omega1 * dt, th2 + result.omega2 * dt)
        end = sol.y[:, -1]
        planar_end = np.array([end[0], end[1], end[3], end[4]], dtype=np.float64)
        max_err = max(max_err, float(np.linalg.norm(planar_end - u_target)))
    return max_err


def discover_ccr4bp_torus_from_resonant_orbit(
    system: ccr4bp.CCR4BPSystem,
    orbit_state0: NDArray[np.float64],
    orbit_period: float,
    *,
    n1: int = 1,
    n2: int = 20,
    m1: int | None = None,
    m2: int | None = None,
    rho_target: float | None = None,
    rho_weight: float = 100.0,
    gauge_weight: float = 30.0,
    tr_solver: Literal["exact", "lsmr"] = "exact",
    tol: float = 1e-6,
    closure_tol: float = 1e-4,
    max_nfev: int = 600,
    closure_dt_frac: float = 0.02,
    notes: str = "",
) -> CCR4BPTorusVariationalResult:
    """Build a 2D pseudospectral CCR4BP torus from a base-CRTBP resonant orbit.

    Projects the (planar) base periodic orbit onto a ``theta1``-flat 2D Fourier
    guess (bootstrap; propagation for the guess only), sets ``omega2 =
    2*pi/orbit_period``, anchors the transverse amplitude at the seed's own value,
    and refines via :func:`correct_ccr4bp_torus_pseudospectral`. This is the entry
    point for the Jupiter-Europa 3:4 positive control (call with
    ``system.mu_gan = 0`` for the machine-precision regression floor, then
    continue the mass up with :func:`continue_ccr4bp_torus_mass`).

    Defaults are tuned for the CCR4BP's WEAK forcing (Ganymede ``mu_gan ~ 8e-5``):
    a SINGLE forcing harmonic (``n1 = 1`` -- the physical response to a
    single-frequency perturber is first-harmonic-dominated) and STRONG gauge /
    rotation-number weights (30 / 100), because at such weak forcing the
    ``theta1``-shear direction is near-degenerate and, with many ``theta1`` modes
    or weak gauges, the trust-region can wander into spurious large ``theta1``
    structure or off-branch rotation numbers (`#612`'s escape mode). ``n2`` is
    large (20) because the reachable symmetric resonant orbit is eccentric and
    needs many internal-angle modes to resolve its periapsis passage.
    """
    m1 = m1 if m1 is not None else 2 * n1 + 3
    m2 = m2 if m2 is not None else 2 * n2 + 3
    coeffs0 = project_base_orbit_to_2d(orbit_state0, orbit_period, system.mu, n1, n2, m1, m2)
    amp_anchor = _transverse_amplitude(coeffs0, n2)
    omega2_0 = 2.0 * np.pi / orbit_period
    return correct_ccr4bp_torus_pseudospectral(
        system,
        coeffs0,
        omega2_0,
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
        notes=notes or "bootstrapped_from_resonant_orbit",
    )


def continue_ccr4bp_torus_mass(
    result: CCR4BPTorusVariationalResult,
    target_system: ccr4bp.CCR4BPSystem,
    *,
    n_steps: int = 8,
    rho_weight: float = 1.0,
    gauge_weight: float = 1.0,
    tr_solver: Literal["exact", "lsmr"] = "exact",
    tol: float = 1e-6,
    closure_tol: float = 1e-4,
    max_nfev: int = 400,
    closure_dt_frac: float = 0.02,
    divergence_factor: float = 1.0e3,
) -> list[CCR4BPTorusVariationalResult]:
    """Homotopy in Ganymede's mass from ``result.system.mu_gan`` to
    ``target_system.mu_gan``, using ONLY integration-free solves.

    Steps ``mu_gan`` linearly (holding ``omega_gan``, ``a_gan``, ``mu`` at the
    target's -- turning UP the forcing amplitude at the physical forcing rate;
    the final step IS the physical system), warm-starting each solve from the
    previous torus. Returns the list of per-step results (the last is the
    target-mass torus). The rotation number is held pinned at the running value
    (the family is selected by the seed orbit, not re-chosen mid-continuation).

    Stops early only on genuine DIVERGENCE -- a non-finite residual or one that
    grows past ``divergence_factor`` times the seed's -- NOT on the strict
    ``converged`` flag (whose ``tol=1e-6`` invariance gate is below the
    truncation floor of the eccentric resonant orbit, so it is expected to read
    ``False`` at a perfectly good ~1e-5 torus; the returned per-step residuals
    are the honest measure).
    """
    steps: list[CCR4BPTorusVariationalResult] = []
    current = result
    floor = max(result.residual_rms, 1.0e-6)
    mu_gans = np.linspace(current.system.mu_gan, target_system.mu_gan, n_steps + 1)[1:]
    for mg in mu_gans:
        sys_step = replace(target_system, mu_gan=float(mg))
        nxt = correct_ccr4bp_torus_pseudospectral(
            sys_step,
            current.coeffs,
            current.omega2,
            n1=current.n1,
            n2=current.n2,
            m1=current.m1,
            m2=current.m2,
            period_multiple=current.period_multiple,
            amplitude_anchor=current.transverse_amplitude,
            rho_target=current.rotation_number,
            rho_weight=rho_weight,
            gauge_weight=gauge_weight,
            tr_solver=tr_solver,
            tol=tol,
            closure_tol=closure_tol,
            max_nfev=max_nfev,
            closure_dt_frac=closure_dt_frac,
            notes=f"mass_continuation_mu_gan={mg:.6e}",
        )
        steps.append(nxt)
        current = nxt
        if not np.isfinite(nxt.residual_rms) or nxt.residual_rms > divergence_factor * floor:
            break
    return steps
