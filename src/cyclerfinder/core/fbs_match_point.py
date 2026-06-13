r"""Ellison 2018 forward-backward-shooting (FBS) match point — massless DSM leg.

Path-B re-transcription of the DSM lane (#226; mining note
``docs/notes/2026-06-10-ellison-2018-analytic-gradients-mining.md`` §3, §7). A
single two-body leg carries ONE interior impulse Δv (the MGAnDSMs n=1 case,
Ellison Sec. II.D). The leg is propagated FORWARD from its left boundary and
BACKWARD from its right boundary to a common match point placed AT the impulse;
the position/velocity mismatch there (the defect, Ellison Eq. 2 with the mass row
dropped — our cycler legs are massless) must vanish for the leg to be dynamically
consistent. Unlike the Takao/Lambert DSM leg (``search/dsm_leg.py``) the
post-impulse velocity is ``v_fwd + Δv`` with Δv an explicit decision variable —
there is NO Lambert solve anywhere (the structural payoff of two-sided shooting,
mining note §1).

Massless simplification (mining note §7, Path B): Ellison's maneuver transition
matrix (Eq. 29) collapses to the identity, and the impulse enters the forward
match state only through the Eq. 42 velocity-slot connection ``∂X_k^+/∂Δv =
[0; I]``. The only non-trivial matrices in the Eq. 31-32 match-point chain are
therefore the analytic two-body STMs (:func:`cyclerfinder.core.kepler_stm.shepperd_stm`).

Phase 1 scope: the single-leg defect and its analytic Jacobian w.r.t. Δv and the
boundary velocities. Phase-TOF partials (Pitkin Eqs. 43-44), the multi-leg chain
(Eqs. 31-32), and the corrector wiring are later phases (see the plan
``docs/superpowers/plans/2026-06-13-ellison-fbs-path-b-plan.md``).

Source: D. H. Ellison, B. A. Conway, J. A. Englander, M. T. Ozimek, "Analytic
Gradient Computation for Bounded-Impulse Trajectory Models Using Two-Sided
Shooting," *Journal of Guidance, Control, and Dynamics*, Vol. 41, No. 7, 2018,
pp. 1449-1462, doi:10.2514/1.G003077 (Eqs. 2, 16, 17, 29, 31-32, 42).

CONSISTENCY DISCIPLINE: Ellison prints no unit-level numeric gradient (mining
note §6); the Jacobian is validated against central differences of the defect,
never a sourced golden. The A1-A6 NLP-scaling caveat noted in the
``nbody/flyby_gradients.py`` docstring carries over.

Units: km, km/s, s, km^3/s^2. State ordered ``[r; v]`` (position rows/columns
first), matching :mod:`cyclerfinder.core.kepler_stm`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.kepler_stm import shepperd_stm

Vec3 = NDArray[np.float64]  # (3,), dtype float64
Vec6 = NDArray[np.float64]  # (6,), dtype float64
Mat6 = NDArray[np.float64]  # (6, 6), dtype float64


class FbsMatchPointError(Exception):
    """Errors raised by the FBS match-point model."""


@dataclass(frozen=True)
class FbsLeg:
    r"""A single massless two-body leg with one interior impulse (MGAnDSMs n=1).

    Attributes
    ----------
    r0, v0:
        Left-boundary inertial state, ``(3,)`` float64, km / km/s.
    rf, vf:
        Right-boundary inertial state, ``(3,)`` float64, km / km/s.
    tof_s:
        Total leg time of flight, s. Must be positive.
    alpha:
        Burn-index fraction in ``(0, 1)`` (Ellison p. 6 ``alpha_1``; our Takao
        ``eta`` is exactly this for n=1, mining note §2). The impulse is applied
        at ``t_burn = alpha * tof_s`` after the leg start; the match point sits
        there.
    mu:
        Central-body gravitational parameter, km^3/s^2 (heliocentric default).
    """

    r0: Vec3
    v0: Vec3
    rf: Vec3
    vf: Vec3
    tof_s: float
    alpha: float
    mu: float = MU_SUN_KM3_S2

    def __post_init__(self) -> None:
        if self.tof_s <= 0.0:
            raise FbsMatchPointError(f"tof_s must be positive, got {self.tof_s}")
        if not (0.0 < self.alpha < 1.0):
            raise FbsMatchPointError(f"alpha must lie in (0, 1), got {self.alpha}")

    @property
    def t_burn_s(self) -> float:
        """Forward time from the leg start to the impulse, s (``alpha * tof_s``)."""
        return self.alpha * self.tof_s


def _two_body_accel(r: Vec3, mu: float) -> Vec3:
    """Two-body gravitational acceleration ``-mu r / |r|^3`` (km/s^2)."""
    r_arr = np.asarray(r, dtype=np.float64)
    r_n = float(np.linalg.norm(r_arr))
    return -mu * r_arr / (r_n * r_n * r_n)


def _forward_match_state(leg: FbsLeg, dv: Vec3) -> tuple[Vec6, Vec3, Vec3, Mat6]:
    """Forward state at the match point (just AFTER the impulse) + STM pieces.

    Returns ``(X_fwd_post (6,), r_fwd, v_fwd_pre, phi_fwd)`` where ``phi_fwd`` is
    the 6x6 STM of the forward coast ``propagate(r0, v0, t_burn)``. The impulse
    adds ``Δv`` to the velocity slot (massless MTM = identity + Eq. 42), which is
    a state shift only — it does not alter ``phi_fwd``.
    """
    r_fwd, v_fwd_pre, phi_fwd = shepperd_stm(leg.r0, leg.v0, leg.t_burn_s, leg.mu)
    x_post = np.empty(6, dtype=np.float64)
    x_post[0:3] = r_fwd
    x_post[3:6] = v_fwd_pre + np.asarray(dv, dtype=np.float64)
    return x_post, r_fwd, v_fwd_pre, phi_fwd


def _backward_match_state(leg: FbsLeg) -> tuple[Vec6, Mat6]:
    """Backward state at the match point + the backward-coast STM.

    Coasts ``(rf, vf)`` by ``-(tof_s - t_burn)`` back to the impulse point.
    Returns ``(X_bwd (6,), phi_bwd)`` with ``phi_bwd = ∂X_bwd/∂(rf, vf)``.
    """
    dt_back = -(leg.tof_s - leg.t_burn_s)
    r_bwd, v_bwd, phi_bwd = shepperd_stm(leg.rf, leg.vf, dt_back, leg.mu)
    x_bwd = np.empty(6, dtype=np.float64)
    x_bwd[0:3] = r_bwd
    x_bwd[3:6] = v_bwd
    return x_bwd, phi_bwd


def match_point_defect(leg: FbsLeg, dv: Vec3) -> Vec6:
    r"""Match-point defect 6-vector ``c_mp = X^B - X^F`` (Ellison Eq. 2, massless).

    ``[r^B - r^F; v^B - v^F]`` — position mismatch (km, 3) and velocity mismatch
    (km/s, 3) between the backward- and forward-propagated states at the impulse
    point. Zero ⇒ the leg is dynamically consistent. No mass row (massless leg),
    no Lambert solve (Δv is the explicit decision variable).
    """
    arr = np.asarray(dv, dtype=np.float64)
    if arr.shape != (3,):
        raise FbsMatchPointError(f"dv must have shape (3,), got {arr.shape}")
    x_fwd, _, _, _ = _forward_match_state(leg, arr)
    x_bwd, _ = _backward_match_state(leg)
    return x_bwd - x_fwd


def match_point_defect_jacobian(
    leg: FbsLeg, dv: Vec3, *, include_phase: bool = False
) -> NDArray[np.float64]:
    r"""Analytic Jacobian of the match-point defect (Ellison Eqs. 16, 31-32, 43-44).

    Columns are ordered ``[∂c/∂Δv (3) | ∂c/∂v0 (3) | ∂c/∂vf (3)]`` for the
    decision sub-vector ``x = [Δv; v0; vf]`` (holds ``r0, rf, tof_s, alpha``
    fixed; ``v0`` / ``vf`` are the boundary-velocity / v∞ slots — the v∞ → state
    map is an additive constant, so ``∂c/∂v0 = ∂c/∂v∞_out`` etc.).

    For the massless single-DSM leg the Eq. 31-32 chain is one STM per side and
    the maneuver connection is the Eq. 42 velocity slot ``[0; I]``:

    * ``∂c/∂Δv  = -[0; I]``           (only ``X^F`` sees Δv, via Eq. 42)
    * ``∂c/∂v0  = -Φ_fwd[:, 3:6]``    (only ``X^F`` sees v0, via the fwd-coast STM)
    * ``∂c/∂vf  = +Φ_bwd[:, 3:6]``    (only ``X^B`` sees vf, via the bwd-coast STM)

    Phase-TOF columns (``include_phase=True``; Pitkin Eqs. 43-44, Ellison Eq. 58)
    append ``[∂c/∂tof_s (1) | ∂c/∂alpha (1)]`` giving a 6x11 Jacobian for the
    decision vector ``x = [Δv; v0; vf; tof_s; alpha]``. The defect depends on
    these only through the propagation intervals ``t_burn = alpha*tof_s`` (forward)
    and ``dt_back = -(1-alpha)*tof_s`` (backward); the time-derivative of a coasted
    state w.r.t. its own propagation interval is the state's flow
    ``[v; a]`` with ``a = -mu r/|r|^3`` (Pitkin's Lagrange-coefficient time
    derivatives reduce to this), so with ``∂t_burn/∂tof = alpha``,
    ``∂dt_back/∂tof = -(1-alpha)``, ``∂t_burn/∂alpha = ∂dt_back/∂alpha = tof``:

    * ``∂c/∂tof   = -(1-alpha)·[v^B; a^B] - alpha·[v^F; a^F]``
    * ``∂c/∂alpha = tof·([v^B; a^B] - [v^F; a^F])``

    (the impulse adds a constant to the forward velocity slot, so it drops out of
    the time derivatives). Validated FD-vs-analytic (consistency; Ellison
    publishes no numeric gradient).
    """
    arr = np.asarray(dv, dtype=np.float64)
    if arr.shape != (3,):
        raise FbsMatchPointError(f"dv must have shape (3,), got {arr.shape}")
    x_fwd, r_fwd, v_fwd_pre, phi_fwd = _forward_match_state(leg, arr)
    x_bwd, phi_bwd = _backward_match_state(leg)

    n_cols = 11 if include_phase else 9
    jac = np.zeros((6, n_cols), dtype=np.float64)
    # ∂c/∂Δv = -[0; I]  (Eq. 42 velocity-slot connection; sign from c = X^B - X^F)
    jac[3:6, 0:3] = -np.eye(3, dtype=np.float64)
    # ∂c/∂v0 = -Φ_fwd[:, 3:6]  (forward-coast STM velocity columns)
    jac[:, 3:6] = -phi_fwd[:, 3:6]
    # ∂c/∂vf = +Φ_bwd[:, 3:6]  (backward-coast STM velocity columns)
    jac[:, 6:9] = phi_bwd[:, 3:6]

    if include_phase:
        # Flow (state time-derivative) at each match state. The forward velocity
        # slot carries +dv, but d/dt of a constant shift is zero, so the flow uses
        # the COAST velocity v_fwd_pre and the coast acceleration at r_fwd.
        flow_fwd = np.empty(6, dtype=np.float64)
        flow_fwd[0:3] = v_fwd_pre
        flow_fwd[3:6] = _two_body_accel(r_fwd, leg.mu)
        flow_bwd = np.empty(6, dtype=np.float64)
        flow_bwd[0:3] = x_bwd[3:6]
        flow_bwd[3:6] = _two_body_accel(x_bwd[0:3], leg.mu)

        one_minus_a = 1.0 - leg.alpha
        # ∂c/∂tof = ∂X^B/∂tof - ∂X^F/∂tof
        jac[:, 9] = -one_minus_a * flow_bwd - leg.alpha * flow_fwd
        # ∂c/∂alpha = ∂X^B/∂alpha - ∂X^F/∂alpha = tof*(flow_bwd - flow_fwd)
        jac[:, 10] = leg.tof_s * (flow_bwd - flow_fwd)

    _ = x_fwd  # forward match state retained for clarity / future use
    return jac
