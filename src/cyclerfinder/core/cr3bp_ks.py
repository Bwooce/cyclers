"""Moon-centred KS-regularised propagator and transition matrix for close passes (#928).

Formulation: Stiefel & Scheifele 1971, total-energy form (9,53) (digest
``docs/notes/2026-10-04-digest-stiefel-scheifele-1971-linear-regular-celestial-mechanics.md``,
sections 3.3, 7.2, 11, 12, 16.1), with the KS map of :mod:`cyclerfinder.core.ks`. The
regularisation centre is one body with gravitational parameter ``K^2``; every other
acceleration is a perturbation. With ``x`` the position relative to the centre, ``r = |x|``,
``u`` the KS vector (``x = L(u) u``), ``w = du/ds`` and ``dt = r ds``:

    u'' + (h/2) u = -(1/4) d/du (|u|^2 V) + (|u|^2/2) L^T P,
    h'  = -|u|^2 dV/dt - 2 (u', L^T P),
    t'  = |u|^2,

where ``V(t, x)`` is the perturbing potential (force ``-grad V``), ``P`` the remaining force and
``h = K^2/r - |xdot|^2/2 - V`` the NEGATIVE total energy (the book's sign: h > 0 bound about the
centre; a hyperbolic pass has h < 0 and the equations hold unchanged). The state is the
10-vector ``y = (u1..u4, w1..w4, h, t)``. With ``-(1/4) d/du(|u|^2 V) = -(1/2) u V -
(r/2) L3^T grad V`` and the Coriolis force of a frame rotating at ``omega`` about +z,
``P = 2 omega (vy, -vx, 0)``, written as ``(r/2) L3^T P = 2 omega L3^T J L3 w`` (no division by
r), every term is regular at ``u = 0`` provided ``V`` and ``grad V`` are finite at the centre.

Model (:class:`MoonCentredCR3BP`): the project's CR3BP of :mod:`cyclerfinder.core.cr3bp`
(rotating frame, unit angular velocity, mu = m2/(m1 + m2), primary of mass 1 - mu at
(-mu, 0, 0), secondary of mass mu at (1 - mu, 0, 0)) centred on the SECONDARY: ``K^2 = mu``,
``x = (X - (1 - mu), Y, Z)``, and ``V`` = centrifugal plus the primary's direct attraction with
the book's interior-rule constant so that ``V(0) = 0`` (S&S section 22, comment 2 p.117):

    V = -(1/2)[(x1 + 1 - mu)^2 + x2^2] - (1 - mu)/r1 + (1 - mu)^2/2 + (1 - mu),
    r1 = |x + e1| (the primary sits at x = -e1).

Hence ``h = C/2 - (1 - mu)^2/2 - (1 - mu)`` with ``C`` the project's ``jacobi_constant`` (which
omits the constant mu(1 - mu)); ``h`` is constant (Coriolis does no work) and the model declares
``h' = 0`` exactly. ``V`` and ``grad V`` are evaluated in a form free of cancellation near the
centre (both are O(r^2) and O(r) there). ``mu = 1`` gives the mu = 0 rotating Kepler problem of
Llibre 1982 about the origin. Passes near the PRIMARY are not regularised here (one centre only,
S&S p.125); the force term is a model object so the elliptic or bicircular problems can be added
as further :class:`KSModel` subclasses (with ``h`` then integrated, ``autonomous = False``).

Monitors: the bilinear relation ``l(u, w) = 0`` (KS (7)) and the regularised Hamiltonian
``K = 2|w|^2 - K^2 + r (h + V) = 0`` (the energy relation (9,73) times r; Peters 1968's
identically-zero K), both evaluated at every accepted step.

Transition matrix (S&S digest section 12; the book prints none): integrate the 10 x 10
variational equations of the system above from the identity, lift with the Jacobian ``J0`` of
the initial KS state with respect to the physical state, project with the Jacobian of
``z = (x, xdot)`` with respect to ``y`` at the end, and subtract ``zdot_f (dt_f/dz0)`` (the fixed
physical-time correction; without it the matrix is wrong at O(1)). The result is the 6 x 6
Jacobian of the physical flow at fixed physical time in the output frame. It is invariant under
the choice of fibre angle and fibre gauge of the lift.

Hand-over (not automated here): Aarseth 1971's criterion (perturbing to central force ratio
gamma = 0.01) places the plain-to-KS switch at about 0.04 lunar distances from the Moon (digest
``docs/notes/2026-10-04-digest-aarseth-1971-direct-integration-n-body.md`` section 7, item 5).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.ks import KS_BASIS, ks_lift, ks_matrix, ks_velocity_lift

if TYPE_CHECKING:
    from scipy.integrate._ivp.ivp import OdeResult

FloatArray = NDArray[np.float64]

_E3 = np.ascontiguousarray(KS_BASIS[:, :3, :])  # _E3[k] = d L3 / d u_k (3 x 4)


class KSModel:
    """Perturbations of a KS-regularised two-body problem (base class: pure Kepler, V = 0).

    Attributes
    ----------
    k2 : gravitational parameter of the regularisation centre (book's K^2).
    omega : rotation rate of the frame about +z; the Coriolis force 2 omega (vy, -vx, 0) is
        applied in regular form. The centrifugal term, if any, belongs in the potential.
    centre : position of the regularisation centre in the output frame (3,).
    autonomous : False when ``V`` depends explicitly on time (then ``potential_dt`` enters
        ``h'``). With ``autonomous`` True and no ``force``, ``h' = 0`` exactly and the
        transition matrix is available; a ``force`` always adds its work ``-2 (w, L^T P)``.

    A subclass overrides ``potential_and_grad`` (and ``potential_hess`` for the transition
    matrix), and optionally ``force`` and ``potential_dt``.
    """

    k2: float = 1.0
    omega: float = 0.0
    centre: FloatArray = np.zeros(3)
    autonomous: bool = True

    def __init__(self, k2: float = 1.0) -> None:
        if not (math.isfinite(k2) and k2 > 0.0):
            raise ValueError(f"KSModel: k2 must be positive, got {k2}")
        self.k2 = float(k2)
        self.centre = np.zeros(3)

    def potential_and_grad(
        self, t: float, x1: float, x2: float, x3: float
    ) -> tuple[float, float, float, float]:
        """``(V, dV/dx1, dV/dx2, dV/dx3)`` at the centre-relative point (the hot path).

        Subclasses override this one method (and ``potential_hess`` for the transition
        matrix); ``potential`` and ``potential_grad`` are derived from it.
        """
        return 0.0, 0.0, 0.0, 0.0

    def potential(self, t: float, x: FloatArray) -> float:
        """Perturbing potential V(t, x) (force -grad V), x relative to the centre."""
        return self.potential_and_grad(t, float(x[0]), float(x[1]), float(x[2]))[0]

    def potential_grad(self, t: float, x: FloatArray) -> FloatArray:
        _, g1, g2, g3 = self.potential_and_grad(t, float(x[0]), float(x[1]), float(x[2]))
        return np.array([g1, g2, g3])

    def potential_hess(self, t: float, x: FloatArray) -> FloatArray:
        return np.zeros((3, 3))

    def potential_dt(self, t: float, x: FloatArray) -> float:
        return 0.0

    def force(self, t: float, x: FloatArray) -> FloatArray | None:
        """Extra non-potential force P(t, x) per unit mass, or None (Coriolis excluded)."""
        return None


class MoonCentredCR3BP(KSModel):
    """The project's CR3BP centred on the secondary (mass mu at (1 - mu, 0, 0)); 0 < mu <= 1."""

    def __init__(self, mu: float) -> None:
        if not (0.0 < mu <= 1.0):
            raise ValueError(f"MoonCentredCR3BP: mu must be in (0, 1], got {mu}")
        super().__init__(k2=mu)
        self.mu = float(mu)
        self.a = 1.0 - self.mu  # primary mass, and the barycentre's distance from the centre
        self.omega = 1.0
        self.centre = np.array([1.0 - self.mu, 0.0, 0.0])
        self.autonomous = True

    def potential_and_grad(
        self, t: float, x1: float, x2: float, x3: float
    ) -> tuple[float, float, float, float]:
        # V = -(x1^2 + x2^2)/2 - a [x1 + 1/r1 - 1] and
        # grad V = -(x1 + a, x2, 0) + a (x + e1)/r1^3, with q = r1^2 - 1 = 2 x1 + |x|^2 exact and
        # the O(1) cancellations removed: x1 + 1/r1 - 1 = x1 q (r1 + 2)/(r1 (1 + r1)^2)
        # - |x|^2/(r1 (1 + r1)) and 1/r1^3 - 1 = -q (r1^2 + r1 + 1)/((1 + r1) r1^3).
        a = self.a
        xx = x1 * x1 + x2 * x2 + x3 * x3
        q = 2.0 * x1 + xx
        r1 = math.sqrt(1.0 + q)
        op = 1.0 + r1
        bracket = x1 * q * (r1 + 2.0) / (r1 * op * op) - xx / (r1 * op)
        pot = -0.5 * (x1 * x1 + x2 * x2) - a * bracket
        r13 = r1 * r1 * r1
        c = a / r13
        g1 = -x1 + c * x1 - a * q * (r1 * r1 + r1 + 1.0) / (op * r13)
        g2 = -x2 + c * x2
        g3 = c * x3
        return pot, g1, g2, g3

    def potential_hess(self, t: float, x: FloatArray) -> FloatArray:
        d = x + np.array([1.0, 0.0, 0.0])
        r1 = float(np.linalg.norm(d))
        hess = (self.a / r1**3) * (np.eye(3) - 3.0 * np.outer(d, d) / (r1 * r1))
        hess[0, 0] -= 1.0
        hess[1, 1] -= 1.0
        return hess

    def energy_constant(self) -> float:
        """c in ``h = C/2 - c`` (C the project's ``jacobi_constant``)."""
        return 0.5 * self.a * self.a + self.a


# ------------------------------------------------------------------------------------------
# KS state <-> physical state


def ks_energy(model: KSModel, t: float, x: FloatArray, v: FloatArray) -> float:
    """Book energy ``h = K^2/r - |v|^2/2 - V(t, x)`` (x relative to the centre)."""
    r = float(np.linalg.norm(x))
    return float(model.k2 / r - 0.5 * float(v @ v) - model.potential(t, x))


def ks_state_from_physical(
    model: KSModel, state6: ArrayLike, t: float = 0.0, fibre_angle: float = 0.0
) -> FloatArray:
    """The 10-vector ``(u, w, h, t)`` of an output-frame physical state (lift (9,69)-(9,72))."""
    z = np.asarray(state6, dtype=np.float64)
    x = z[:3] - model.centre
    v = z[3:6].copy()
    u = ks_lift(x, phi=fibre_angle)
    y = np.empty(10)
    y[:4] = u
    y[4:8] = ks_velocity_lift(u, v)
    y[8] = ks_energy(model, t, x, v)
    y[9] = t
    return y


def physical_from_ks_state(model: KSModel, y: ArrayLike) -> FloatArray:
    """Output-frame ``(x, xdot)`` of a KS 10-vector; undefined (raises) at ``u = 0``."""
    ya = np.asarray(y, dtype=np.float64)
    u, w = ya[:4], ya[4:8]
    r = float(u @ u)
    if r == 0.0:
        raise ZeroDivisionError("physical_from_ks_state: velocity undefined at collision (u = 0)")
    l3 = ks_matrix(u)[:3]
    out = np.empty(6)
    out[:3] = l3 @ u + model.centre
    out[3:] = 2.0 * (l3 @ w) / r
    return out


def ks_monitors(model: KSModel, y: ArrayLike) -> tuple[float, float]:
    """``(K, l)``: regularised Hamiltonian ``2|w|^2 - K^2 + r (h + V)`` and bilinear ``l(u, w)``.

    Both vanish identically along a motion; regular at ``u = 0``.
    """
    ya = np.asarray(y, dtype=np.float64)
    u, w, h, t = ya[:4], ya[4:8], float(ya[8]), float(ya[9])
    r = float(u @ u)
    x = ks_matrix(u)[:3] @ u
    k_res = 2.0 * float(w @ w) - model.k2 + r * (h + model.potential(t, x))
    bil = u[3] * w[0] - u[2] * w[1] + u[1] * w[2] - u[0] * w[3]
    return float(k_res), float(bil)


def physical_acceleration(model: KSModel, t: float, x: FloatArray, v: FloatArray) -> FloatArray:
    """Physical acceleration of the model (x relative to the centre, singular at x = 0)."""
    r = float(np.linalg.norm(x))
    a = -model.k2 * x / r**3 - model.potential_grad(t, x)
    if model.omega != 0.0:
        a = a + 2.0 * model.omega * np.array([v[1], -v[0], 0.0])
    f = model.force(t, x)
    if f is not None:
        a = a + f
    return np.asarray(a, dtype=np.float64)


# ------------------------------------------------------------------------------------------
# Right-hand sides


def ks_rhs(s: float, y: FloatArray, model: KSModel) -> FloatArray:
    """d y / d s for ``y = (u, w, h, t)`` (S&S (9,53) with Coriolis in regular form)."""
    u1, u2, u3, u4, w1, w2, w3, w4, h, t = y.tolist()
    r = u1 * u1 + u2 * u2 + u3 * u3 + u4 * u4
    x1 = u1 * u1 - u2 * u2 - u3 * u3 + u4 * u4
    x2 = 2.0 * (u1 * u2 - u3 * u4)
    x3 = 2.0 * (u1 * u3 + u2 * u4)
    pot, g1, g2, g3 = model.potential_and_grad(t, x1, x2, x3)
    # force vector f (3) whose image (r/2) L3^T f enters the u-equation: -grad V, Coriolis
    # 2 omega J xdot = 2 omega J (2 L3 w / r) (so (r/2) L3^T of it is 2 omega L3^T J L3 w), and
    # the model's extra force.
    hr = 0.5 * r
    f1, f2, f3 = -hr * g1, -hr * g2, -hr * g3
    if model.omega != 0.0:
        lw1 = u1 * w1 - u2 * w2 - u3 * w3 + u4 * w4
        lw2 = u2 * w1 + u1 * w2 - u4 * w3 - u3 * w4
        two_om = 2.0 * model.omega
        f1 += two_om * lw2
        f2 -= two_om * lw1
    hdot = 0.0
    extra = model.force(t, np.array([x1, x2, x3]))
    if extra is not None:
        p1, p2, p3 = (float(c) for c in extra)
        # h' = -2 (w, L3^T P)
        hdot -= 2.0 * (
            p1 * (u1 * w1 - u2 * w2 - u3 * w3 + u4 * w4)
            + p2 * (u2 * w1 + u1 * w2 - u4 * w3 - u3 * w4)
            + p3 * (u3 * w1 + u4 * w2 + u1 * w3 + u2 * w4)
        )
        f1 += hr * p1
        f2 += hr * p2
        f3 += hr * p3
    if not model.autonomous:
        hdot -= r * model.potential_dt(t, np.array([x1, x2, x3]))
    k = -0.5 * (h + pot)
    out = np.empty(10)
    out[0], out[1], out[2], out[3] = w1, w2, w3, w4
    # L3^T f = f1 (u1, -u2, -u3, u4) + f2 (u2, u1, -u4, -u3) + f3 (u3, u4, u1, u2)
    out[4] = k * u1 + f1 * u1 + f2 * u2 + f3 * u3
    out[5] = k * u2 - f1 * u2 + f2 * u1 + f3 * u4
    out[6] = k * u3 - f1 * u3 - f2 * u4 + f3 * u1
    out[7] = k * u4 + f1 * u4 - f2 * u3 + f3 * u2
    out[8] = hdot
    out[9] = r
    return out


_JC = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])  # P_cor = 2 omega JC v


def ks_jacobian(y: FloatArray, model: KSModel) -> FloatArray:
    """Analytic 10 x 10 Jacobian of :func:`ks_rhs` (autonomous models only)."""
    u = y[:4]
    w = y[4:8]
    h = float(y[8])
    t = float(y[9])
    r = float(u @ u)
    l3 = ks_matrix(u)[:3]
    x = l3 @ u
    if not model.autonomous or model.force(t, x) is not None:
        raise NotImplementedError(
            "ks_jacobian: only autonomous potential models (plus Coriolis) are supported"
        )
    pot = model.potential(t, x)
    grad = model.potential_grad(t, x)
    hess = model.potential_hess(t, x)
    ltg = l3.T @ grad
    # N_g[:, k] = E3_k^T grad  (d/du of L3(u)^T grad at fixed grad)
    n_g = np.einsum("kij,i->jk", _E3, grad)
    j_uu = (
        -(0.5 * h + 0.5 * pot) * np.eye(4)
        - np.outer(u, ltg)
        - np.outer(ltg, u)
        - 0.5 * r * n_g
        - r * (l3.T @ hess @ l3)
    )
    jac = np.zeros((10, 10))
    jac[:4, 4:8] = np.eye(4)
    if model.omega != 0.0:
        jl = _JC @ l3
        lw = l3 @ w
        jlw = _JC @ lw
        # d/du_k [L3^T JC L3 w] = E3_k^T JC L3 w + L3^T JC E3_k w
        term1 = np.einsum("kij,i->jk", _E3, jlw)
        term2 = l3.T @ _JC @ np.einsum("kij,j->ik", _E3, w)
        j_uu = j_uu + 2.0 * model.omega * (term1 + term2)
        jac[4:8, 4:8] = 2.0 * model.omega * (l3.T @ jl)
    jac[4:8, :4] = j_uu
    jac[4:8, 8] = -0.5 * u
    jac[9, :4] = 2.0 * u
    return jac


def ks_stm_rhs(s: float, y: FloatArray, model: KSModel) -> FloatArray:
    """State plus 10 x 10 variational equations, ``y = (state10, Phi.ravel())``."""
    out = np.empty(110)
    out[:10] = ks_rhs(s, y[:10], model)
    out[10:] = (ks_jacobian(y[:10], model) @ y[10:].reshape(10, 10)).ravel()
    return out


# ------------------------------------------------------------------------------------------
# Lift and projection Jacobians


def lift_jacobian(model: KSModel, y0: FloatArray, gauge: ArrayLike | None = None) -> FloatArray:
    """``J0 = d y0 / d z0`` (10 x 6) of the lift at ``y0``; ``gauge`` (6,) adds a fibre component.

    Position: ``du0 = L3^T dx / (2r) + tangent * (gauge . dz)``; velocity from
    ``w0 = (1/2) L3(u0)^T v``; energy from ``h = K^2/r - v^2/2 - V``; ``t0`` fixed.
    """
    u = y0[:4]
    t = float(y0[9])
    r = float(u @ u)
    l3 = ks_matrix(u)[:3]
    x = l3 @ u
    v = 2.0 * (l3 @ y0[4:8]) / r
    du = np.zeros((4, 6))
    du[:, :3] = l3.T / (2.0 * r)
    if gauge is not None:
        tangent = np.array([-u[3], u[2], -u[1], u[0]])
        du = du + np.outer(tangent, np.asarray(gauge, dtype=np.float64))
    m_v = np.einsum("kij,i->jk", _E3, v)  # d (L3(u)^T v) / du at fixed v
    dw = 0.5 * (m_v @ du)
    dw[:, 3:] += 0.5 * l3.T
    dh = np.zeros(6)
    dh[:3] = -model.k2 * x / r**3 - model.potential_grad(t, x)
    dh[3:] = -v
    jac = np.zeros((10, 6))
    jac[:4] = du
    jac[4:8] = dw
    jac[8] = dh
    return jac


def projection_jacobian(y: FloatArray) -> FloatArray:
    """``d z / d y`` (6 x 10) of ``z = (L3 u, 2 L3 w / r)`` (the centre offset is constant)."""
    u = y[:4]
    w = y[4:8]
    r = float(u @ u)
    l3 = ks_matrix(u)[:3]
    g_w = np.einsum("kij,j->ik", _E3, w)  # d (L3(u) w) / du at fixed w
    lw = l3 @ w
    jac = np.zeros((6, 10))
    jac[:3, :4] = 2.0 * l3
    jac[3:, :4] = (2.0 / r) * g_w - (4.0 / (r * r)) * np.outer(lw, u)
    jac[3:, 4:8] = (2.0 / r) * l3
    return jac


# ------------------------------------------------------------------------------------------
# Propagation


@dataclass(frozen=True)
class KSArc:
    """Result of :func:`propagate_ks` (physical state in the model's output frame)."""

    state: FloatArray  # (6,) at t_final
    t: float
    s: float
    y: FloatArray  # (10,) KS state at the end
    stm: FloatArray | None  # (6, 6) d state(t_final) / d state0, fixed physical time
    energy_drift: float  # integrated h minus algebraic h at the end
    hamiltonian_max: float  # max |K| over accepted steps
    bilinear_max: float  # max |l(u, w)| over accepted steps
    r_min: float  # closest approach to the centre on the arc
    nfev: int


def integrate_ks(
    model: KSModel,
    y0: ArrayLike,
    s_span: tuple[float, float],
    *,
    rtol: float = 1e-12,
    atol: float = 1e-14,
    with_stm: bool = False,
    events: object = None,
    dense_output: bool = False,
    max_step: float = np.inf,
) -> OdeResult:
    """Integrate the KS system (and optionally its variational equations) in ``s`` (DOP853).

    Raw entry point for starts at ``u = 0`` (ejection) and section events; returns scipy's
    ``OdeResult``. With ``with_stm`` the state is ``(y10, Phi.ravel())`` starting from the
    identity.
    """
    y0a = np.asarray(y0, dtype=np.float64)
    if with_stm and y0a.size == 10:
        y0a = np.concatenate([y0a, np.eye(10).ravel()])
    rhs = ks_stm_rhs if with_stm else ks_rhs
    sol: OdeResult = solve_ivp(  # type: ignore[call-overload]
        rhs,
        s_span,
        y0a,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        args=(model,),
        events=events,
        dense_output=dense_output,
        max_step=max_step,
    )
    return sol


def _min_distance_event(s: float, y: FloatArray, model: KSModel) -> float:
    return float(y[:4] @ y[4:8])  # r'/2


_min_distance_event.direction = 0.0  # type: ignore[attr-defined]


def propagate_ks(
    model: KSModel,
    state0: ArrayLike,
    t_final: float,
    *,
    t0: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-14,
    with_stm: bool = False,
    fibre_angle: float = 0.0,
    gauge: ArrayLike | None = None,
    s_max: float = 1e12,
) -> KSArc:
    """Propagate an output-frame state from ``t0`` to the physical time ``t_final``.

    The arc is integrated in ``s`` to a terminal event ``t(s) = t_final``; the end point is then
    recomputed by integrating from the last accepted step (not interpolated) and polished by
    Newton steps ``ds = (t_final - t)/r``. ``fibre_angle`` and ``gauge`` choose the lift; the
    physical result does not depend on them. ``with_stm`` adds the 6 x 6 fixed-time transition
    matrix (autonomous models only).
    """
    if t_final == t0:
        raise ValueError("propagate_ks: t_final equals t0")
    sign = 1.0 if t_final > t0 else -1.0
    y0 = ks_state_from_physical(model, state0, t=t0, fibre_angle=fibre_angle)
    y_start = np.concatenate([y0, np.eye(10).ravel()]) if with_stm else y0

    def _time_event(s: float, y: FloatArray, model: KSModel) -> float:
        return float(y[9] - t_final)

    _time_event.terminal = True  # type: ignore[attr-defined]
    sol = integrate_ks(
        model,
        y_start,
        (0.0, sign * s_max),
        rtol=rtol,
        atol=atol,
        with_stm=with_stm,
        events=[_time_event, _min_distance_event],
    )
    t_events = sol.t_events
    y_events = sol.y_events
    if sol.status != 1 or t_events is None or y_events is None or t_events[0].size != 1:
        raise RuntimeError(f"propagate_ks: did not reach t_final ({sol.message})")
    s_grid = sol.t
    y_grid = np.asarray(sol.y, dtype=np.float64)
    nfev = int(sol.nfev)
    s_event = float(t_events[0][0])
    # recompute the end from the last accepted step before the event
    k_prev = len(s_grid) - 2
    s_cur = float(s_grid[k_prev])
    y_cur: FloatArray = y_grid[:, k_prev].copy()
    for _ in range(6):
        if s_event == s_cur:
            break
        seg = integrate_ks(model, y_cur, (s_cur, s_event), rtol=rtol, atol=atol, with_stm=with_stm)
        nfev += int(seg.nfev)
        s_cur, y_cur = s_event, np.asarray(seg.y[:, -1], dtype=np.float64)
        dt = t_final - float(y_cur[9])
        if abs(dt) <= 4e-16 * max(1.0, abs(t_final)):
            break
        s_event = s_cur + dt / float(y_cur[:4] @ y_cur[:4])
    y_end = y_cur[:10]
    # monitors over the accepted steps
    k_max, b_max = 0.0, 0.0
    for j in range(y_grid.shape[1] - 1):
        k_res, bil = ks_monitors(model, y_grid[:10, j])
        k_max, b_max = max(k_max, abs(k_res)), max(b_max, abs(bil))
    k_res, bil = ks_monitors(model, y_end)
    k_max, b_max = max(k_max, abs(k_res)), max(b_max, abs(bil))
    r_vals = [float(y_start[:4] @ y_start[:4]), float(y_end[:4] @ y_end[:4])]
    for ye in y_events[1]:
        r_vals.append(float(ye[:4] @ ye[:4]))
    state = physical_from_ks_state(model, y_end)
    x_rel = state[:3] - model.centre
    h_alg = ks_energy(model, float(y_end[9]), x_rel, state[3:])
    stm = None
    if with_stm:
        phi_y = y_cur[10:].reshape(10, 10)
        chain = phi_y @ lift_jacobian(model, y0, gauge)
        phi_s = projection_jacobian(y_end) @ chain
        zdot = np.concatenate(
            [state[3:], physical_acceleration(model, float(y_end[9]), x_rel, state[3:])]
        )
        stm = phi_s - np.outer(zdot, chain[9])
    return KSArc(
        state=state,
        t=float(y_end[9]),
        s=s_cur,
        y=y_end.copy(),
        stm=stm,
        energy_drift=float(y_end[8]) - h_alg,
        hamiltonian_max=k_max,
        bilinear_max=b_max,
        r_min=min(r_vals),
        nfev=nfev,
    )
