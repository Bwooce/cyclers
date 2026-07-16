"""Seedless 2D pseudospectral quasi-periodic-torus corrector for the CR3BP (#612).

Motivation and the wall this crosses
-------------------------------------
`#606`/`#611` built seedless *periodic-orbit* correctors (a truncated Fourier
series in ONE angle, driven to satisfy the equations of motion at collocation
points, with NO forward integration in the search) and used them to cross two
documented "shooting-fragility" walls -- `#556`'s CR3BP L1 quasi-halo periodic
orbit and `#538`/`#544`'s QBCP EM-L1/L2 periodic orbit. `#611`'s own docstring
flagged the natural next step and declined it as out of scope: the genuinely
QUASI-PERIODIC 2-TORUS (two angles, a free rotation number), which is
`#556`'s actual parked target. This module builds exactly that, for the CR3BP.

The concrete blocker `#556` names is that ``genome.qp_tori.correct_qp_torus``
(the GMOS / Olikara-Scheeres invariant-circle corrector) cannot converge the
EM L1 quasi-halo torus above amplitude ~0.01. `#612`'s diagnosis (live,
2026-07-16) confirmed *why*, and it is the same fragility `#606`/`#611`
attacked: the GMOS residual resolves the LONGITUDINAL angle by **stroboscopic
flow integration** -- ``_gmos_residual`` propagates each invariant-circle
sample point forward by ``t_strob = 2*pi/omega_long`` and matches the result
to the rotated circle. The parent L1 halo at C=3.15 has a **monodromy spectral
radius of ~1540** (measured), so that one-period propagation -- and every
finite-difference column of its Jacobian -- is amplified ~1540x. As the torus
amplitude grows and its sample points spread away from the (violently
unstable) center, the propagation and its Jacobian become catastrophically
ill-conditioned; the measured invariance residual degrades monotonically
(5.5e-7 at amp 5e-4 -> 2.5e-4 at amp 1e-2) and a single amp>=0.015 build does
not finish in 90-250 s. This is a *shooting* fragility wearing an
invariant-circle disguise, not a conceptual limit of the torus itself.

The 2D pseudospectral formulation (no integration anywhere in the search)
-------------------------------------------------------------------------
A quasi-periodic invariant 2-torus is a map ``u : T^2 -> R^6`` from two angles
``(theta1, theta2) = (theta_long, theta_trans)`` to the CR3BP state, on which
the flow is linear: ``theta1' = omega1``, ``theta2' = omega2``. A trajectory
``x(t) = u(omega1 t + phi1, omega2 t + phi2)`` then satisfies ``x' = f(x)``
(the CR3BP rotating-frame RHS, :func:`cyclerfinder.core.cr3bp.cr3bp_eom`) if
and only if ``u`` solves the **quasi-periodic invariance PDE**

    omega1 * du/dtheta1  +  omega2 * du/dtheta2  =  f(u(theta1, theta2))     (*)

for all ``(theta1, theta2) in T^2`` (chain rule). This is the standard
"big-system" torus equation (Jorba 2001; Schilder-Osinga-Vogt 2005;
Olikara-Scheeres 2010 give the invariant-circle reduction ``correct_qp_torus``
implements instead). Crucially, unlike the invariant-circle reduction, the
full PDE (*) contains **no time integration**: ``f(u)`` is evaluated pointwise
and the angle-derivatives are taken SPECTRALLY, so the ~1540x monodromy
amplification never enters the residual or its Jacobian. The residual's
Jacobian is the local operator ``omega1 d/dtheta1 + omega2 d/dtheta2 - Df(u)``,
well-conditioned regardless of how unstable the flow is -- exactly the reason
`#606`/`#611` beat shooting on their walls, now lifted to the 2-torus.

Represent each of the six state components as a real tensor-product Fourier
series (mirroring `#606`'s ``_eval_series``, one extra angle):

    u_c(theta1, theta2) = sum_{a,b} C[c,a,b] * phi_a(theta1) * phi_b(theta2)

where ``phi_0 = 1``, ``phi_{k} = cos(k .)`` for ``k=1..N``, ``phi_{N+k} =
sin(k .)`` for ``k=1..N`` is the real 1D Fourier basis of size ``2N+1`` (so
reality is automatic -- no complex conjugate bookkeeping). Angle-derivatives
are analytic term-by-term. The residual (*) is evaluated on a uniform 2D
collocation grid; ``f`` is applied to all six grid-sampled components at once.

Free unknowns and gauge fixing
------------------------------
Free unknowns: the coefficients ``C`` (shape ``(6, 2N1+1, 2N2+1)``) and BOTH
frequencies ``(omega1, omega2)`` (their ratio ``omega2/omega1`` is the
rotation number). Two continuous phase-shift symmetries
``(theta1, theta2) -> (theta1 + d1, theta2 + d2)`` and one family-selection /
collapse degeneracy must be removed -- THREE scalar conditions, appended as
residual rows (the ``correct_qp_torus`` convention, not `#606`'s
variable-elimination convention; kept soft so the 2D packing stays a plain
dense array):

* **Longitudinal phase gauge** (``theta1``): pin the ``sin(theta1)``,
  ``theta2``-constant coefficient of ``x`` to zero -- the 2-angle analogue of
  `#606`'s ``sin_x1 = 0`` (pins the time origin to an x-extremum of the mean
  orbit). This is a NULL direction of the PDE residual, so pinning it costs the
  fit nothing.
* **Transverse phase gauge** (``theta2``): pin the ``sin(theta2)``,
  ``theta1``-constant coefficient of a chosen coordinate ``c*`` (the one with
  the largest transverse-mode content in the initial guess -- the analogue of
  ``correct_qp_torus``'s ``phase_pin_idx``) to zero.
* **Transverse amplitude anchor**: fix the L2 norm of all ``theta2``-first-
  harmonic (``k2 = 1``) coefficients to a caller-specified value. Without it
  the optimizer can collapse the transverse structure to zero -- i.e. slide
  back to the parent PERIODIC orbit (a zero-``theta2``-amplitude torus, with
  ``omega2`` then undetermined). This is ``correct_qp_torus``'s
  ``amplitude_pin`` / `#606`'s amplitude anchor: a legitimate "how big a torus"
  family-selection choice, not information leaked from a target.

Bootstrap and continuation (integration allowed ONLY in the initial guess)
--------------------------------------------------------------------------
Like `#606`, the SEARCH never integrates, but a rough initial guess is needed.
The cleanest, most honest bootstrap is `#606`'s "one known family member"
pattern: take a SMALL-amplitude torus from the existing GMOS corrector (which
converges cleanly below amp~0.01, where the stroboscopic amplification is
survivable), sample it on the 2D grid via ``evaluate_torus`` (this uses
propagation -- fine, it is only the guess), least-squares project onto the 2D
Fourier basis, and hand off to the pseudospectral solve. To cross the wall,
step the amplitude anchor upward and warm-start each solve from the previous
converged 2D torus -- a natural-parameter continuation done ENTIRELY with this
module's own (integration-free) solves, reaching amplitudes the GMOS corrector
cannot. See :func:`discover_qp_torus` and :func:`continue_qp_torus_amplitude`.

Independent closure check (NOT circular with the residual)
----------------------------------------------------------
The minimized quantity is the algebraic PDE residual on the collocation grid.
The independent check propagates a torus point ``u(theta1, theta2)`` through
the TRUE nonlinear flow for a SHORT time ``dt`` and compares against
``u(theta1 + omega1 dt, theta2 + omega2 dt)`` -- genuinely independent (it uses
the integrator, not the spectral series). ``dt`` is deliberately a small
fraction of ``t_strob``: over the full stroboscopic period the ~1540x
amplification would swamp the check, so a short hop (amplification ``exp(lambda
dt)`` with ``lambda ~ ln(1540)/T``) keeps it meaningful at large amplitude
while still exercising the real dynamics off the collocation grid.

Discipline / scope
------------------
* The returned torus is OUR computation; no novelty claim, no catalogue
  writeback. This is a capability build (a corrector that crosses a documented
  convergence wall), not a discovery result.
* Crossing the *convergence* wall (converging at amp>=0.02) is a DIFFERENT
  thing from reaching Owen & Baresi's L1 latitudinal frequency 0.2739: `#555`
  established that at C=3.15 the L1 quasi-halo rotation number is ENERGY-pinned
  near 0.074, flat in amplitude -- a physical family fact, independent of which
  corrector is used, re-confirmed here (see the test file). This module removes
  the corrector limitation `#556` named; it does not (and cannot) move the
  rotation number that the energy fixes.

References
----------
* Jorba, A. (2001). "Numerical computation of the normal behaviour of
  invariant curves of n-dimensional maps." Nonlinearity 14, 943-976.
* Schilder, F., Osinga, H. M., & Vogt, W. (2005). "Continuation of quasi-
  periodic invariant tori." SIAM J. Appl. Dyn. Syst. 4(3), 459-488.
* Olikara, Z., & Scheeres, D. (2010). "Numerical Method for Computing Quasi-
  Periodic Orbits and Their Stability in the Restricted Three-Body Problem."
* Moore, C. (1993), PRL 70, 3675 / Chenciner & Montgomery (2000), Ann. Math.
  152, 881 -- the harmonic-balance/action-collocation precedent `#606` cites.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, evaluate_torus
from cyclerfinder.search.outcome_log import log_outcome

_N_STATE = 6


@dataclass(frozen=True)
class QPTorusVariationalResult:
    """A CR3BP quasi-periodic 2-torus found by seedless 2D pseudospectral collocation.

    ``coeffs`` is the ``(6, 2*n1+1, 2*n2+1)`` real tensor-product Fourier
    coefficient array (see module docstring for the basis); ``omega1``/
    ``omega2`` are the longitudinal/transverse frequencies (rad/TU) and
    ``rotation_number = omega2 / omega1``. ``residual_rms`` is the RMS of the
    quasi-periodic invariance-PDE residual at the 2D collocation points (the
    quantity minimized). ``closure_residual`` is the INDEPENDENT short-time
    flow check (propagate ``u(theta1,theta2)`` for ``closure_dt`` through the
    true nonlinear CR3BP EOM, compare to ``u(theta1+omega1*dt, theta2+omega2*dt)``
    -- NOT circular with ``residual_rms``). ``transverse_amplitude`` is the
    anchored L2 norm of the ``k2=1`` coefficients.
    """

    system: cr3bp.CR3BPSystem
    coeffs: NDArray[np.float64]
    omega1: float
    omega2: float
    rotation_number: float
    n1: int
    n2: int
    m1: int
    m2: int
    transverse_amplitude: float
    residual_rms: float
    closure_residual: float
    converged: bool
    n_iter: int
    jacobi: float
    notes: str = ""


# ---------------------------------------------------------------------------
# Real tensor-product Fourier basis.
# ---------------------------------------------------------------------------


def _basis_matrices(
    n_modes: int, thetas: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ``(P, Pd)`` of shape ``(len(thetas), 2*n_modes+1)``: the real 1D
    Fourier basis ``[1, cos(1.), .., cos(N.), sin(1.), .., sin(N.)]`` evaluated
    at ``thetas``, and its theta-derivative.
    """
    m = thetas.size
    size = 2 * n_modes + 1
    mat_p = np.zeros((m, size))
    mat_pd = np.zeros((m, size))
    mat_p[:, 0] = 1.0  # constant mode -> derivative 0
    for k in range(1, n_modes + 1):
        ck, sk = np.cos(k * thetas), np.sin(k * thetas)
        mat_p[:, k] = ck
        mat_pd[:, k] = -k * sk
        mat_p[:, n_modes + k] = sk
        mat_pd[:, n_modes + k] = k * ck
    return mat_p, mat_pd


def _n_free(n1: int, n2: int) -> int:
    return _N_STATE * (2 * n1 + 1) * (2 * n2 + 1) + 2


def _pack(coeffs: NDArray[np.float64], omega1: float, omega2: float) -> NDArray[np.float64]:
    return np.concatenate([coeffs.reshape(-1), [omega1, omega2]])


def _unpack(z: NDArray[np.float64], n1: int, n2: int) -> tuple[NDArray[np.float64], float, float]:
    n_coef = _N_STATE * (2 * n1 + 1) * (2 * n2 + 1)
    coeffs = z[:n_coef].reshape(_N_STATE, 2 * n1 + 1, 2 * n2 + 1)
    return coeffs, float(z[-2]), float(z[-1])


def evaluate_torus_state(
    result: QPTorusVariationalResult,
    theta1: float | NDArray[np.float64],
    theta2: float | NDArray[np.float64],
) -> NDArray[np.float64]:
    """Evaluate the 2D Fourier torus state ``u(theta1, theta2)`` (6-vector).

    Scalars in -> ``(6,)`` out. Uses the analytic Fourier series directly, no
    propagation.
    """
    t1 = np.atleast_1d(np.asarray(theta1, dtype=np.float64))
    t2 = np.atleast_1d(np.asarray(theta2, dtype=np.float64))
    p1, _ = _basis_matrices(result.n1, t1)
    p2, _ = _basis_matrices(result.n2, t2)
    # u_c = p1 @ C_c @ p2.T -> take diagonal pairing (t1[i], t2[i])
    out = np.empty((t1.size, _N_STATE))
    for c in range(_N_STATE):
        full = p1 @ result.coeffs[c] @ p2.T  # (len t1, len t2)
        out[:, c] = np.diagonal(full) if t1.size == t2.size else full[:, 0]
    if np.isscalar(theta1) and np.isscalar(theta2):
        return out[0]
    return out


# ---------------------------------------------------------------------------
# Invariance-PDE residual (no integration).
# ---------------------------------------------------------------------------


def _k2_first_harmonic_cols(n2: int) -> tuple[int, int]:
    """Column indices in the theta2 basis for the k2=1 cos and sin terms."""
    return 1, n2 + 1  # cos(1 theta2), sin(1 theta2)


def _residual(
    z: NDArray[np.float64],
    mu: float,
    n1: int,
    n2: int,
    p1: NDArray[np.float64],
    p1d: NDArray[np.float64],
    p2: NDArray[np.float64],
    p2d: NDArray[np.float64],
    phase1_coord: int,
    phase2_coord: int,
    amplitude_anchor: float,
    gauge_weight: float,
    rho_target: float,
    rho_weight: float,
) -> NDArray[np.float64]:
    """Real residual: quasi-periodic invariance PDE on the 2D grid + 4 gauges.

    Gauge rows: longitudinal phase, transverse phase, transverse amplitude
    anchor, and a rotation-number pin ``omega2 = rho_target * omega1``. The
    rotation-number pin (Jorba-style rho fixing) is essential: without it the
    least-squares escapes to a DEGENERATE spurious branch (``omega2 -> 0``, a
    resonant "tube of periodic orbits" that reaches machine-precision residual
    and so always beats the truncation-limited genuine torus). Pinning rho is
    also physically correct for the L1 quasi-halo target, whose rotation number
    `#555` showed is energy-pinned (flat in amplitude) at C=3.15. Both omegas
    remain free unknowns in ``z``; only their ratio is constrained.
    """
    coeffs, omega1, omega2 = _unpack(z, n1, n2)
    if not (np.isfinite(omega1) and np.isfinite(omega2)) or omega1 <= 0.0:
        return np.full(p1.shape[0] * p2.shape[0] * _N_STATE + 4, 1e6)

    # u and its angle-derivatives on the grid: shape (m1, m2) per coordinate.
    m1, m2 = p1.shape[0], p2.shape[0]
    u = np.empty((_N_STATE, m1, m2))
    du1 = np.empty((_N_STATE, m1, m2))
    du2 = np.empty((_N_STATE, m1, m2))
    for c in range(_N_STATE):
        cc = coeffs[c]
        u[c] = p1 @ cc @ p2.T
        du1[c] = p1d @ cc @ p2.T
        du2[c] = p1 @ cc @ p2d.T

    x, y, zc, vx, vy, vz = (u[i] for i in range(6))
    r1 = np.sqrt((x + mu) ** 2 + y * y + zc * zc)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y * y + zc * zc)
    r1c, r2c = r1**3, r2**3
    ax = x + 2.0 * vy - (1.0 - mu) * (x + mu) / r1c - mu * (x - 1.0 + mu) / r2c
    ay = y - 2.0 * vx - (1.0 - mu) * y / r1c - mu * y / r2c
    az = -(1.0 - mu) * zc / r1c - mu * zc / r2c
    rhs = np.stack([vx, vy, vz, ax, ay, az], axis=0)  # (6, m1, m2)

    lhs = omega1 * du1 + omega2 * du2
    pde_res = (lhs - rhs).reshape(-1)

    # Gauge / anchor rows.
    sin1_col_t1 = n1 + 1  # sin(1 theta1)
    _, sin1_col_t2 = _k2_first_harmonic_cols(n2)
    phase1 = coeffs[phase1_coord, sin1_col_t1, 0]  # sin(theta1) x const(theta2)
    phase2 = coeffs[phase2_coord, 0, sin1_col_t2]  # const(theta1) x sin(theta2)
    cos1_col_t2, _ = _k2_first_harmonic_cols(n2)
    trans = coeffs[:, :, cos1_col_t2] ** 2 + coeffs[:, :, sin1_col_t2] ** 2
    amp = float(np.sqrt(np.sum(trans)))
    gauge = np.array(
        [
            gauge_weight * phase1,
            gauge_weight * phase2,
            gauge_weight * (amp - amplitude_anchor),
            rho_weight * (omega2 - rho_target * omega1),
        ]
    )
    return np.concatenate([pde_res, gauge])


def _transverse_amplitude(coeffs: NDArray[np.float64], n2: int) -> float:
    cos1, sin1 = _k2_first_harmonic_cols(n2)
    return float(np.sqrt(np.sum(coeffs[:, :, cos1] ** 2 + coeffs[:, :, sin1] ** 2)))


def _cr3bp_jacobian_grid(u: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Vectorized 6x6 CR3BP RHS Jacobian ``df/dstate`` at every grid point.

    ``u`` has shape ``(6, m1, m2)``; returns ``(6, 6, m1, m2)``. The matrix is
    exactly ``cr3bp.cr3bp_stm_eom``'s ``A`` (pseudo-potential Hessian + Coriolis
    + velocity-selection block), evaluated pointwise -- the analytic Jacobian of
    the invariance-PDE's ``f(u)`` term.
    """
    x, y, z = u[0], u[1], u[2]
    m1, m2 = x.shape
    r1 = np.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    r1c, r2c, r1f, r2f = r1**3, r2**3, r1**5, r2**5
    om1 = 1.0 - mu
    uxx = (
        1 - om1 / r1c - mu / r2c + 3 * om1 * (x + mu) ** 2 / r1f + 3 * mu * (x - 1 + mu) ** 2 / r2f
    )
    uyy = 1 - om1 / r1c - mu / r2c + 3 * om1 * y * y / r1f + 3 * mu * y * y / r2f
    uzz = -om1 / r1c - mu / r2c + 3 * om1 * z * z / r1f + 3 * mu * z * z / r2f
    uxy = 3 * om1 * (x + mu) * y / r1f + 3 * mu * (x - 1 + mu) * y / r2f
    uxz = 3 * om1 * (x + mu) * z / r1f + 3 * mu * (x - 1 + mu) * z / r2f
    uyz = 3 * om1 * y * z / r1f + 3 * mu * y * z / r2f
    jf = np.zeros((6, 6, m1, m2))
    jf[0, 3] = jf[1, 4] = jf[2, 5] = 1.0
    jf[3, 0], jf[3, 1], jf[3, 2] = uxx, uxy, uxz
    jf[4, 0], jf[4, 1], jf[4, 2] = uxy, uyy, uyz
    jf[5, 0], jf[5, 1], jf[5, 2] = uxz, uyz, uzz
    jf[3, 4], jf[4, 3] = 2.0, -2.0
    return jf


def _jacobian(
    z: NDArray[np.float64],
    mu: float,
    n1: int,
    n2: int,
    p1: NDArray[np.float64],
    p1d: NDArray[np.float64],
    p2: NDArray[np.float64],
    p2d: NDArray[np.float64],
    phase1_coord: int,
    phase2_coord: int,
    amplitude_anchor: float,
    gauge_weight: float,
    rho_target: float,
    rho_weight: float,
) -> NDArray[np.float64]:
    """Analytic Jacobian of :func:`_residual` (exact; no finite differencing).

    Row block: PDE residual (6*m1*m2) then 4 gauge rows. Columns: vec(coeffs)
    (6*A1*A2) then [omega1, omega2]. See module docstring for the derivation
    ``dR_c/dC_d = delta_cd (w1 P1d P2 + w2 P1 P2d) - Jf_cd (P1 P2)``,
    ``dR_c/dw1 = du1_c``, ``dR_c/dw2 = du2_c``.
    """
    coeffs, omega1, omega2 = _unpack(z, n1, n2)
    a1, a2 = 2 * n1 + 1, 2 * n2 + 1
    m1, m2 = p1.shape[0], p2.shape[0]
    n_coef = _N_STATE * a1 * a2
    n_free = n_coef + 2
    n_pde = _N_STATE * m1 * m2
    jac = np.zeros((n_pde + 4, n_free))

    # Basis tensors on the grid: B[i,j,a,b] = P1[i,a] P2[j,b], etc.
    b_val = np.einsum("ia,jb->ijab", p1, p2)
    b_d1 = np.einsum("ia,jb->ijab", p1d, p2)
    b_d2 = np.einsum("ia,jb->ijab", p1, p2d)
    lhs_basis = (omega1 * b_d1 + omega2 * b_d2).reshape(m1 * m2, a1 * a2)
    val_basis = b_val.reshape(m1 * m2, a1 * a2)

    # u and its angle derivatives (for the omega columns).
    u = np.empty((_N_STATE, m1, m2))
    du1 = np.empty((_N_STATE, m1, m2))
    du2 = np.empty((_N_STATE, m1, m2))
    for c in range(_N_STATE):
        u[c] = p1 @ coeffs[c] @ p2.T
        du1[c] = p1d @ coeffs[c] @ p2.T
        du2[c] = p1 @ coeffs[c] @ p2d.T
    jf = _cr3bp_jacobian_grid(u, mu)  # (6,6,m1,m2)

    for c in range(_N_STATE):
        rows = slice(c * m1 * m2, (c + 1) * m1 * m2)
        # Diagonal LHS block (same coord).
        cols_c = slice(c * a1 * a2, (c + 1) * a1 * a2)
        jac[rows, cols_c] = lhs_basis
        # -Jf coupling across all coords d.
        for d in range(_N_STATE):
            cols_d = slice(d * a1 * a2, (d + 1) * a1 * a2)
            jf_cd = jf[c, d].reshape(m1 * m2)  # per grid point
            jac[rows, cols_d] -= jf_cd[:, None] * val_basis
        # omega columns.
        jac[rows, n_coef] = du1[c].reshape(m1 * m2)
        jac[rows, n_coef + 1] = du2[c].reshape(m1 * m2)

    # Gauge rows.
    sin1_col_t1 = n1 + 1
    _, sin1_col_t2 = _k2_first_harmonic_cols(n2)
    cos1_col_t2, _ = _k2_first_harmonic_cols(n2)

    def _coef_index(c: int, a: int, b: int) -> int:
        return c * a1 * a2 + a * a2 + b

    jac[n_pde + 0, _coef_index(phase1_coord, sin1_col_t1, 0)] = gauge_weight
    jac[n_pde + 1, _coef_index(phase2_coord, 0, sin1_col_t2)] = gauge_weight
    amp = _transverse_amplitude(coeffs, n2)
    if amp > 0:
        for c in range(_N_STATE):
            for a in range(a1):
                jac[n_pde + 2, _coef_index(c, a, cos1_col_t2)] = (
                    gauge_weight * coeffs[c, a, cos1_col_t2] / amp
                )
                jac[n_pde + 2, _coef_index(c, a, sin1_col_t2)] = (
                    gauge_weight * coeffs[c, a, sin1_col_t2] / amp
                )
    # Rotation-number pin row: rho_weight * (omega2 - rho_target * omega1).
    jac[n_pde + 3, n_coef] = -rho_weight * rho_target
    jac[n_pde + 3, n_coef + 1] = rho_weight
    return jac


# ---------------------------------------------------------------------------
# Initial-guess construction from a small-amplitude GMOS torus.
# ---------------------------------------------------------------------------


def project_gmos_torus_to_2d(
    torus: QPTorus, n1: int, n2: int, m1: int, m2: int
) -> NDArray[np.float64]:
    """Least-squares project a converged GMOS :class:`QPTorus` onto the 2D real
    Fourier basis -> initial ``(6, 2*n1+1, 2*n2+1)`` coefficient array.

    Samples the GMOS torus on the ``(theta1, theta2)`` grid via
    ``evaluate_torus`` (which propagates -- allowed, this is only the guess),
    then solves the two 1D least-squares fits ``U = P1 C P2^T`` for ``C``.
    """
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, _ = _basis_matrices(n1, t1)
    p2, _ = _basis_matrices(n2, t2)
    p1_pinv = np.linalg.pinv(p1)
    p2_pinv = np.linalg.pinv(p2)
    coeffs = np.empty((_N_STATE, 2 * n1 + 1, 2 * n2 + 1))
    u_grid = np.empty((_N_STATE, m1, m2))
    for i, a in enumerate(t1):
        for j, b in enumerate(t2):
            u_grid[:, i, j] = evaluate_torus(torus, float(a), float(b))
    for c in range(_N_STATE):
        # C = P1^+ U P2^+^T  (solves U = P1 C P2^T in least squares)
        coeffs[c] = p1_pinv @ u_grid[c] @ p2_pinv.T
    return coeffs


# ---------------------------------------------------------------------------
# The corrector.
# ---------------------------------------------------------------------------


def correct_qp_torus_pseudospectral(
    system: cr3bp.CR3BPSystem,
    coeffs0: NDArray[np.float64],
    omega1_0: float,
    omega2_0: float,
    *,
    n1: int,
    n2: int,
    amplitude_anchor: float,
    rho_target: float | None = None,
    rho_weight: float = 1.0,
    m1: int | None = None,
    m2: int | None = None,
    tol: float = 1e-6,
    closure_tol: float = 1e-5,
    max_nfev: int = 200,
    gauge_weight: float = 1.0,
    closure_dt_frac: float = 0.02,
    n_closure_samples: int = 12,
    notes: str = "",
) -> QPTorusVariationalResult:
    """Refine a 2D-Fourier torus guess by driving the invariance-PDE residual to
    zero -- no forward integration in the search (see module docstring).

    ``coeffs0`` shape ``(6, 2*n1+1, 2*n2+1)``; ``omega1_0``/``omega2_0`` the
    initial frequencies. ``amplitude_anchor`` fixes the transverse (k2=1)
    coefficient L2 norm. Collocation grid defaults to ``m1 = 2*n1+3``,
    ``m2 = 2*n2+3`` (mild oversampling). Uses ``scipy.optimize.least_squares``
    with ``method="trf"`` and a bounded ``max_nfev`` (finite-difference
    Jacobian; ``method="trf"`` respects the ``nfev`` budget far better than
    ``"lm"`` on ill-conditioned problems -- the `#611` lesson).
    """
    if coeffs0.shape != (_N_STATE, 2 * n1 + 1, 2 * n2 + 1):
        raise ValueError(
            f"coeffs0 shape {coeffs0.shape} != expected {(_N_STATE, 2 * n1 + 1, 2 * n2 + 1)}"
        )
    mu = system.mu
    if rho_target is None:
        rho_target = omega2_0 / omega1_0
    m1 = m1 if m1 is not None else 2 * n1 + 3
    m2 = m2 if m2 is not None else 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)

    # Phase-pin coordinates: theta1 gauge on x (coord 0); theta2 gauge on the
    # coordinate with the largest transverse (k2=1) content in the guess.
    cos1, sin1 = _k2_first_harmonic_cols(n2)
    trans_by_coord = np.sqrt(np.sum(coeffs0[:, :, cos1] ** 2 + coeffs0[:, :, sin1] ** 2, axis=1))
    phase2_coord = int(np.argmax(trans_by_coord))
    phase1_coord = 0

    z0 = _pack(coeffs0, omega1_0, omega2_0)
    solver_args = (
        mu,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        phase1_coord,
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
        x_scale="jac",
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=max_nfev,
    )
    coeffs, omega1, omega2 = _unpack(sol.x, n1, n2)
    n_res = m1 * m2 * _N_STATE
    residual_rms = float(np.sqrt(np.sum(sol.fun[:n_res] ** 2) / n_res))
    amp = _transverse_amplitude(coeffs, n2)
    rot = omega2 / omega1 if omega1 != 0 else float("nan")

    result = QPTorusVariationalResult(
        system=system,
        coeffs=coeffs,
        omega1=omega1,
        omega2=omega2,
        rotation_number=rot,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        transverse_amplitude=amp,
        residual_rms=residual_rms,
        closure_residual=float("nan"),
        converged=False,
        n_iter=int(sol.nfev),
        jacobi=float("nan"),
        notes=notes,
    )
    closure = _independent_closure(result, closure_dt_frac, n_closure_samples)
    jacobi = float(cr3bp.jacobi_constant(evaluate_torus_state(result, 0.0, 0.0), mu))
    converged = (residual_rms < tol) and (closure < closure_tol)
    result = QPTorusVariationalResult(
        system=system,
        coeffs=coeffs,
        omega1=omega1,
        omega2=omega2,
        rotation_number=rot,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        transverse_amplitude=amp,
        residual_rms=residual_rms,
        closure_residual=closure,
        converged=converged,
        n_iter=int(sol.nfev),
        jacobi=jacobi,
        notes=notes,
    )
    log_outcome(
        solver="variational_qp_torus.correct_qp_torus_pseudospectral",
        inputs={
            "mu": float(mu),
            "n1": int(n1),
            "n2": int(n2),
            "amplitude_anchor": float(amplitude_anchor),
        },
        outcome={
            "converged": bool(converged),
            "residual_rms": residual_rms,
            "closure_residual": closure,
            "omega1": float(omega1),
            "omega2": float(omega2),
            "rotation_number": float(rot),
            "transverse_amplitude": float(amp),
        },
        meta={"primary": system.primary, "secondary": system.secondary},
    )
    return result


def _independent_closure(result: QPTorusVariationalResult, dt_frac: float, n_samples: int) -> float:
    """Short-time flow-consistency check (independent of the algebraic residual).

    For several base angles, propagate ``u(theta1, theta2)`` for
    ``dt = dt_frac * (2*pi/omega1)`` through the true CR3BP EOM and compare to
    ``u(theta1 + omega1*dt, theta2 + omega2*dt)``. Short ``dt`` keeps the
    unstable-flow amplification bounded (see module docstring). Returns the max
    over samples.
    """
    t_strob = 2 * np.pi / result.omega1
    dt = dt_frac * t_strob
    rng = np.random.default_rng(0xC0FFEE)
    max_err = 0.0
    for _ in range(n_samples):
        th1 = float(rng.uniform(0, 2 * np.pi))
        th2 = float(rng.uniform(0, 2 * np.pi))
        u0 = evaluate_torus_state(result, th1, th2)
        try:
            arc = cr3bp.propagate(result.system, u0, dt, with_stm=False)
        except RuntimeError:
            return float("inf")
        u_target = evaluate_torus_state(result, th1 + result.omega1 * dt, th2 + result.omega2 * dt)
        max_err = max(max_err, float(np.linalg.norm(arc.state_f - u_target)))
    return max_err


def discover_qp_torus(
    system: cr3bp.CR3BPSystem,
    seed_gmos_torus: QPTorus,
    *,
    n1: int = 12,
    n2: int = 4,
    m1: int | None = None,
    m2: int | None = None,
    rho_target: float | None = None,
    rho_weight: float = 1.0,
    tol: float = 1e-6,
    closure_tol: float = 1e-5,
    max_nfev: int = 200,
    closure_dt_frac: float = 0.02,
    notes: str = "",
) -> QPTorusVariationalResult:
    """Build a 2D pseudospectral torus from a small-amplitude GMOS seed.

    Projects ``seed_gmos_torus`` onto the 2D Fourier basis (bootstrap; uses
    propagation for the guess only), sets ``(omega1, omega2) = (omega_long,
    omega_trans)``, anchors the transverse amplitude at the seed's own value,
    and refines via :func:`correct_qp_torus_pseudospectral` (no integration in
    the solve). This is the entry point for the positive control; walk to
    higher amplitude with :func:`continue_qp_torus_amplitude`.
    """
    m1 = m1 if m1 is not None else 2 * n1 + 3
    m2 = m2 if m2 is not None else 2 * n2 + 3
    coeffs0 = project_gmos_torus_to_2d(seed_gmos_torus, n1, n2, m1, m2)
    amp_anchor = _transverse_amplitude(coeffs0, n2)
    return correct_qp_torus_pseudospectral(
        system,
        coeffs0,
        seed_gmos_torus.omega_long,
        seed_gmos_torus.omega_trans,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        amplitude_anchor=amp_anchor,
        rho_target=rho_target,
        rho_weight=rho_weight,
        tol=tol,
        closure_tol=closure_tol,
        max_nfev=max_nfev,
        closure_dt_frac=closure_dt_frac,
        notes=notes or "bootstrapped_from_gmos",
    )


def continue_qp_torus_amplitude(
    result: QPTorusVariationalResult,
    target_amplitude: float,
    *,
    n_steps: int = 8,
    tol: float = 1e-6,
    closure_tol: float = 1e-5,
    max_nfev: int = 200,
    closure_dt_frac: float = 0.02,
) -> list[QPTorusVariationalResult]:
    """Natural-parameter continuation in transverse amplitude, using ONLY this
    module's own integration-free solves.

    Steps the amplitude anchor linearly from ``result.transverse_amplitude`` to
    ``target_amplitude`` in ``n_steps`` steps, warm-starting each solve from the
    previous converged torus. Returns the list of per-step results (the last is
    the target-amplitude torus). Stops early if a step fails to converge.
    """
    steps: list[QPTorusVariationalResult] = []
    current = result
    amps = np.linspace(current.transverse_amplitude, target_amplitude, n_steps + 1)[1:]
    for amp in amps:
        nxt = correct_qp_torus_pseudospectral(
            current.system,
            current.coeffs,
            current.omega1,
            current.omega2,
            n1=current.n1,
            n2=current.n2,
            m1=current.m1,
            m2=current.m2,
            amplitude_anchor=float(amp),
            rho_target=current.rotation_number,
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
