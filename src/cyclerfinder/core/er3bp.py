"""Elliptic Restricted 3-Body Problem (ER3BP) core dynamics.

Pulsating-rotating frame with true anomaly `f` as the independent variable.
Equations of motion derived from Szebehely (1967) Ch. 10.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import CR3BPSystem, StmMode, _r1_r2


@dataclass(frozen=True)
class ER3BPSystem:
    """ER3BP system definition (mass ratio and eccentricity)."""

    mu: float
    e: float
    primary_name: str
    secondary_name: str

    @classmethod
    def from_cr3bp(cls, cr3bp: CR3BPSystem, e: float) -> ER3BPSystem:
        return cls(
            mu=cr3bp.mu,
            e=e,
            primary_name=cr3bp.primary,
            secondary_name=cr3bp.secondary,
        )


def er3bp_eom(f: float, state6: NDArray[np.float64], mu: float, e: float) -> NDArray[np.float64]:
    """ER3BP equations of motion.

    Independent variable `f` is true anomaly.
    State is [x, y, z, x', y', z'] where prime denotes d/df.

    The pulsating potential is:
    omega = (1 / (1 + e*cos(f))) * [ 0.5*(x^2 + y^2 - e*z^2*cos(f)) + (1-mu)/r1 + mu/r2 ]
    """
    x, y, z, xprime, yprime, zprime = (float(v) for v in state6)

    r1, r2 = _r1_r2(x, y, z, mu)
    r1_3 = r1**3
    r2_3 = r2**3

    # Gravitational terms
    grav_x = (1.0 - mu) * (x + mu) / r1_3 + mu * (x - 1.0 + mu) / r2_3
    grav_y = (1.0 - mu) * y / r1_3 + mu * y / r2_3
    grav_z = (1.0 - mu) * z / r1_3 + mu * z / r2_3

    # Scale factor
    scale = 1.0 / (1.0 + e * math.cos(f))

    # Accelerations
    xdoubleprime = 2.0 * yprime + scale * (x - grav_x)
    ydoubleprime = -2.0 * xprime + scale * (y - grav_y)
    zdoubleprime = -scale * (e * math.cos(f) * z + grav_z)

    return np.array(
        [
            xprime,
            yprime,
            zprime,
            xdoubleprime,
            ydoubleprime,
            zdoubleprime,
        ]
    )


def er3bp_stm_eom(
    f: float, state42: NDArray[np.float64], mu: float, e: float
) -> NDArray[np.float64]:
    """ER3BP variational equations of motion (state + STM)."""
    x, y, z, xprime, yprime, zprime = (float(v) for v in state42[:6])

    r1, r2 = _r1_r2(x, y, z, mu)
    r1_3 = r1**3
    r2_3 = r2**3
    r1_5 = r1**5
    r2_5 = r2**5

    # EOM evaluation
    grav_x = (1.0 - mu) * (x + mu) / r1_3 + mu * (x - 1.0 + mu) / r2_3
    grav_y = (1.0 - mu) * y / r1_3 + mu * y / r2_3
    grav_z = (1.0 - mu) * z / r1_3 + mu * z / r2_3

    scale = 1.0 / (1.0 + e * math.cos(f))

    xdoubleprime = 2.0 * yprime + scale * (x - grav_x)
    ydoubleprime = -2.0 * xprime + scale * (y - grav_y)
    zdoubleprime = -scale * (e * math.cos(f) * z + grav_z)

    # State derivative
    dstate = np.array(
        [
            xprime,
            yprime,
            zprime,
            xdoubleprime,
            ydoubleprime,
            zdoubleprime,
        ]
    )

    # Variational terms (Jacobian of EOMs w.r.t state)
    # Let w_xx = d(omega_x)/dx etc., where omega is the pulsating pseudo-potential
    # omega_x = scale * (x - grav_x)

    u_xx = scale * (
        1.0
        - (1.0 - mu) / r1_3
        - mu / r2_3
        + 3.0 * (1.0 - mu) * (x + mu) ** 2 / r1_5
        + 3.0 * mu * (x - 1.0 + mu) ** 2 / r2_5
    )
    u_yy = scale * (
        1.0
        - (1.0 - mu) / r1_3
        - mu / r2_3
        + 3.0 * (1.0 - mu) * y**2 / r1_5
        + 3.0 * mu * y**2 / r2_5
    )
    u_zz = scale * (
        -e * math.cos(f)
        - (1.0 - mu) / r1_3
        - mu / r2_3
        + 3.0 * (1.0 - mu) * z**2 / r1_5
        + 3.0 * mu * z**2 / r2_5
    )

    u_xy = scale * (3.0 * (1.0 - mu) * (x + mu) * y / r1_5 + 3.0 * mu * (x - 1.0 + mu) * y / r2_5)
    u_xz = scale * (3.0 * (1.0 - mu) * (x + mu) * z / r1_5 + 3.0 * mu * (x - 1.0 + mu) * z / r2_5)
    u_yz = scale * (3.0 * (1.0 - mu) * y * z / r1_5 + 3.0 * mu * y * z / r2_5)

    a_mat = np.zeros((6, 6))
    a_mat[0, 3] = 1.0
    a_mat[1, 4] = 1.0
    a_mat[2, 5] = 1.0

    a_mat[3, 0] = u_xx
    a_mat[3, 1] = u_xy
    a_mat[3, 2] = u_xz
    a_mat[3, 4] = 2.0

    a_mat[4, 0] = u_xy
    a_mat[4, 1] = u_yy
    a_mat[4, 2] = u_yz
    a_mat[4, 3] = -2.0

    a_mat[5, 0] = u_xz
    a_mat[5, 1] = u_yz
    a_mat[5, 2] = u_zz

    phi = state42[6:].reshape((6, 6))
    phidot = a_mat @ phi

    return np.concatenate((dstate, phidot.flatten()))


def propagate_er3bp(
    state6: NDArray[np.float64],
    f_span: tuple[float, float],
    sys: ER3BPSystem,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    with_stm: bool = False,
    stm_mode: StmMode = "variable",
    method: Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"] = "DOP853",
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Propagate the ER3BP state.

    Returns:
        (f_eval, state_history, STM_final) if with_stm=True (state_history is 6xN, STM_final is 6x6)
        (f_eval, state_history, empty_array) if with_stm=False
    """
    if not with_stm:
        sol = solve_ivp(
            fun=er3bp_eom,
            t_span=f_span,
            y0=state6,
            args=(sys.mu, sys.e),
            method=method,
            rtol=rtol,
            atol=atol,
        )
        return sol.t, sol.y, np.zeros((6, 6))

    y0_stm = np.concatenate((state6, np.eye(6).flatten()))

    if stm_mode == "variable":
        sol = solve_ivp(
            fun=er3bp_stm_eom,
            t_span=f_span,
            y0=y0_stm,
            args=(sys.mu, sys.e),
            method=method,
            rtol=rtol,
            atol=atol,
        )
        return sol.t, sol.y[:6, :], sol.y[6:, -1].reshape((6, 6))

    elif stm_mode == "fixed_path":
        # 1. State-only propagation to fix the grid
        sol_state = solve_ivp(
            fun=er3bp_eom,
            t_span=f_span,
            y0=state6,
            args=(sys.mu, sys.e),
            method=method,
            rtol=rtol,
            atol=atol,
        )
        f_eval = sol_state.t

        # 2. STM propagation along the fixed grid
        sol_stm = solve_ivp(
            fun=er3bp_stm_eom,
            t_span=(f_span[0], f_span[1]),
            t_eval=f_eval,
            y0=y0_stm,
            args=(sys.mu, sys.e),
            method=method,
            rtol=rtol,
            atol=atol,
        )
        return sol_stm.t, sol_stm.y[:6, :], sol_stm.y[6:, -1].reshape((6, 6))

    else:
        raise ValueError(f"Unknown stm_mode: {stm_mode}")
