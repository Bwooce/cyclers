"""Circular-restricted three-body problem (CR3BP) dynamics core (spec 2026-06-10).

Nondimensional rotating frame; mu = m2/(m1+m2); primary at (-mu,0,0), secondary at
(1-mu,0,0). See the plan's "CR3BP equations" block. Pure: math/numpy/scipy +
core.satellites only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES


def _r1_r2(x: float, y: float, z: float, mu: float) -> tuple[float, float]:
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return r1, r2


def jacobi_constant(state6: NDArray[np.float64], mu: float) -> float:
    """C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2 (conserved)."""
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    r1, r2 = _r1_r2(x, y, z, mu)
    return float(
        (x * x + y * y) + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - (vx * vx + vy * vy + vz * vz)
    )


def cr3bp_eom(t: float, state6: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Rotating-frame equations of motion (state [x,y,z,vx,vy,vz])."""
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    r1, r2 = _r1_r2(x, y, z, mu)
    r1c, r2c = r1**3, r2**3
    ax = x + 2.0 * vy - (1.0 - mu) * (x + mu) / r1c - mu * (x - 1.0 + mu) / r2c
    ay = y - 2.0 * vx - (1.0 - mu) * y / r1c - mu * y / r2c
    az = -(1.0 - mu) * z / r1c - mu * z / r2c
    return np.array([vx, vy, vz, ax, ay, az], dtype=np.float64)


@dataclass(frozen=True)
class CR3BPSystem:
    mu: float
    primary: str
    secondary: str
    l_km: float  # characteristic length (secondary SMA about primary)
    t_s: float  # characteristic time = sqrt(l^3 / (G(m1+m2)))


@dataclass(frozen=True)
class CR3BPArc:
    state_f: NDArray[np.float64]
    stm: NDArray[np.float64] | None
    t: float


def cr3bp_stm_eom(t: float, y42: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """State (6) + flattened 6x6 STM (36) variational EOM."""
    s = y42[:6]
    phi = y42[6:].reshape(6, 6)
    x, y, z = float(s[0]), float(s[1]), float(s[2])
    r1, r2 = _r1_r2(x, y, z, mu)
    r1c, r2c = r1**3, r2**3
    r1f, r2f = r1**5, r2**5
    om1 = 1.0 - mu
    # Pseudo-potential second derivatives (Uxx etc.) for the A matrix.
    uxx = (
        1 - om1 / r1c - mu / r2c + 3 * om1 * (x + mu) ** 2 / r1f + 3 * mu * (x - 1 + mu) ** 2 / r2f
    )
    uyy = 1 - om1 / r1c - mu / r2c + 3 * om1 * y * y / r1f + 3 * mu * y * y / r2f
    uzz = -om1 / r1c - mu / r2c + 3 * om1 * z * z / r1f + 3 * mu * z * z / r2f
    uxy = 3 * om1 * (x + mu) * y / r1f + 3 * mu * (x - 1 + mu) * y / r2f
    uxz = 3 * om1 * (x + mu) * z / r1f + 3 * mu * (x - 1 + mu) * z / r2f
    uyz = 3 * om1 * y * z / r1f + 3 * mu * y * z / r2f
    mat_a = np.zeros((6, 6))
    mat_a[0:3, 3:6] = np.eye(3)
    mat_a[3, 0], mat_a[3, 1], mat_a[3, 2] = uxx, uxy, uxz
    mat_a[4, 0], mat_a[4, 1], mat_a[4, 2] = uxy, uyy, uyz
    mat_a[5, 0], mat_a[5, 1], mat_a[5, 2] = uxz, uyz, uzz
    mat_a[3, 4], mat_a[4, 3] = 2.0, -2.0  # Coriolis
    dphi = mat_a @ phi
    return np.concatenate([cr3bp_eom(t, s, mu), dphi.reshape(36)])


def propagate(
    system: CR3BPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    with_stm: bool = False,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> CR3BPArc:
    """Propagate a state (and optionally the STM) for nondimensional time ``t``."""
    if with_stm:
        y0 = np.concatenate([np.asarray(state6, float), np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp_stm_eom,
            (0.0, t),
            y0,
            args=(system.mu,),
            rtol=rtol,
            atol=atol,
            method="DOP853",
            dense_output=False,
        )
        yf = sol.y[:, -1]
        return CR3BPArc(state_f=yf[:6], stm=yf[6:].reshape(6, 6), t=t)
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t),
        np.asarray(state6, float),
        args=(system.mu,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
    )
    return CR3BPArc(state_f=sol.y[:, -1], stm=None, t=t)


def cr3bp_system(primary: str, secondary: str) -> CR3BPSystem:
    """Build a CR3BPSystem from the registry: mu = GM2/(GM1+GM2), scales from SMA."""
    gm1 = PRIMARIES[primary]
    gm2 = SATELLITES[secondary].mu_km3_s2
    mu = gm2 / (gm1 + gm2)
    l_km = SATELLITES[secondary].sma_km
    t_s = float(math.sqrt(l_km**3 / (gm1 + gm2)))
    return CR3BPSystem(mu=mu, primary=primary, secondary=secondary, l_km=l_km, t_s=t_s)
