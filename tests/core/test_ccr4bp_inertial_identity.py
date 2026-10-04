"""#893: the concentric circular restricted four-body model against an independent derivation.

The model in ``core/ccr4bp.py`` is written in the frame rotating with the planet and the base
moon. The same physics written in a NON-rotating frame has no sign convention to get wrong: the
planet and the base moon move counter-clockwise on circles about their barycentre at rate 1, the
perturbing moon moves counter-clockwise on a circle about that barycentre at its own inertial
rate ``1 + omega_gan`` (slower than 1 for an outer moon, so ``omega_gan < 0`` and the perturber
regresses in the rotating frame), and it acts on the particle directly and through the
barycentre's acceleration. The two must agree to integration accuracy.

This is the check that `core/bcr4bp.py` lacked until 2026-10-04 (#891), when its perturber was
found to move the wrong way round.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp

FloatArray = NDArray[np.float64]


def _rot_to_inertial(state: FloatArray, t: float) -> FloatArray:
    c, s = math.cos(t), math.sin(t)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    r = rot @ state[:3]
    v = rot @ state[3:] + np.cross([0.0, 0.0, 1.0], r)
    return np.concatenate([r, v])


def _inertial_to_rot(state: FloatArray, t: float) -> FloatArray:
    c, s = math.cos(t), math.sin(t)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    r = rot.T @ state[:3]
    v = rot.T @ (state[3:] - np.cross([0.0, 0.0, 1.0], state[:3]))
    return np.concatenate([r, v])


def _inertial_rhs(t: float, state: FloatArray, system: ccr4bp.CCR4BPSystem) -> FloatArray:
    mu, mu_p, a_p = system.mu, system.mu_gan, system.a_gan
    r = state[:3]
    c, s = math.cos(t), math.sin(t)
    r_planet = np.array([-mu * c, -mu * s, 0.0])
    r_base = np.array([(1.0 - mu) * c, (1.0 - mu) * s, 0.0])
    th = system.theta_gan0 + (1.0 + system.omega_gan) * t
    r_pert = a_p * np.array([math.cos(th), math.sin(th), 0.0])
    acc = (
        -(1.0 - mu) * (r - r_planet) / np.linalg.norm(r - r_planet) ** 3
        - mu * (r - r_base) / np.linalg.norm(r - r_base) ** 3
        - mu_p * (r - r_pert) / np.linalg.norm(r - r_pert) ** 3
        - mu_p * r_pert / a_p**3
    )
    return np.concatenate([state[3:], acc])


def _difference(system: ccr4bp.CCR4BPSystem, x0: FloatArray, t_end: float) -> float:
    rotating = solve_ivp(
        ccr4bp.ccr4bp_eom, (0.0, t_end), x0, args=(system,), method="DOP853", rtol=1e-13, atol=1e-13
    ).y[:, -1]
    inertial = solve_ivp(
        _inertial_rhs,
        (0.0, t_end),
        _rot_to_inertial(x0, 0.0),
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    ).y[:, -1]
    return float(np.linalg.norm(rotating - _inertial_to_rot(inertial, t_end)))


_X0 = np.array([0.6, 0.1, 0.02, 0.05, 0.75, 0.01])
_T_END = 8.0


@pytest.mark.parametrize("theta0", [0.0, 1.3])
def test_outer_perturber_matches_the_inertial_derivation(theta0: float) -> None:
    system = dataclasses.replace(ccr4bp.jupiter_europa_ganymede_default(), theta_gan0=theta0)
    assert system.omega_gan < 0.0  # an outer moon regresses in the base moon's frame
    assert _difference(system, _X0, _T_END) < 1e-9


def test_inner_perturber_matches_the_inertial_derivation() -> None:
    """A perturber inside the base moon's orbit advances in the rotating frame."""
    base = ccr4bp.jupiter_europa_ganymede_default()
    a_inner = 0.62
    system = dataclasses.replace(
        base,
        a_gan=a_inner,
        omega_gan=ccr4bp.two_body_synodic_rate(base.mu, base.mu_gan, a_inner),
        theta_gan0=0.4,
    )
    assert system.omega_gan > 0.0
    x0 = np.array([1.25, 0.1, 0.02, 0.05, -0.35, 0.01])
    assert _difference(system, x0, _T_END) < 1e-9


def test_the_comparison_detects_a_reversed_perturber() -> None:
    """Control: with the perturber's rate reversed in sign the two frames disagree, and with the
    perturber removed they agree, so the comparison itself is sound and discriminating."""
    system = ccr4bp.jupiter_europa_ganymede_default()
    reversed_rate = dataclasses.replace(system, omega_gan=-system.omega_gan)
    rotating = solve_ivp(
        ccr4bp.ccr4bp_eom,
        (0.0, _T_END),
        _X0,
        args=(reversed_rate,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    ).y[:, -1]
    inertial = solve_ivp(
        _inertial_rhs,
        (0.0, _T_END),
        _rot_to_inertial(_X0, 0.0),
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    ).y[:, -1]
    assert float(np.linalg.norm(rotating - _inertial_to_rot(inertial, _T_END))) > 1e-6
    no_perturber = dataclasses.replace(system, mu_gan=0.0)
    assert _difference(no_perturber, _X0, _T_END) < 1e-9
