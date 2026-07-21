"""Task #665: cannonball solar-radiation-pressure (SRP) augmented real-binary
(k1,k2)-cycler dynamics.

Motivation (see `#665`'s own `data/OUTSTANDING.md` bullet and `#661`'s
strategy pass that proposed it): `#549`/`#657`/`#659`/`#660` found GRAVITY-ONLY
negatives for six real binary-asteroid systems (Patroclus-Menoetius, Didymos-
Dimorphos, Orcus-Vanth, Eris-Dysnomia, Sila-Nunam, Lempo-Hiisi). Per
`[[project_negative_results_registry]]`'s own "empty is always conditional on
the method" rule, every one of those negatives is GRAVITY-ONLY-conditional: at
small-body GMs, solar radiation pressure is a real, non-negligible
perturbation for spacecraft-scale area-to-mass ratios. This module adds a
cannonball SRP force term to the existing `real_binary_kk_sweep.py` CR3BP
genome/corrector so those six systems can be re-swept under a genuinely new
(not just re-tuned) model axis.

Physics: the simplest, most defensible model
----------------------------------------------
SRP on a real binary asteroid is caused by the SUN, a THIRD body external to
the binary's own mutual orbit. The Sun-to-binary-barycenter direction, in an
inertial frame, changes on the binary's HELIOCENTRIC period (years) -- far
slower than the mutual orbital period this module's CR3BP rotating frame spins
at (hours to days for every `REAL_BINARY_SYSTEMS` entry). Expressed IN that
fast-rotating frame, though, an inertially-near-fixed Sun direction appears to
sweep a full circle once per mutual-orbit period -- i.e. genuinely FAST, not
slowly-varying, on the frame's own timescale. A fully faithful treatment would
add an explicitly time-periodic forcing term (breaking the autonomous
structure this whole corrector's fixed-Jacobi machinery depends on -- see
`cr3bp_periodic.correct_symmetric_fixed_jacobi`).

This module instead implements the "simplest model" explicitly named in the
positive-control literature survey (Voyatzis, Gkolias, Gaitanas & Tsiganis
2025, CMDA 138:2, abstract: "the small body as a point mass with constant SRP
magnitude and direction, is integrable" -- see
`docs/notes/2026-07-21-srp-binary-asteroid-corpus.md`): SRP as a CONSTANT
vector, fixed magnitude and direction, in the ROTATING frame -- a snapshot
approximation valid for one instant of the binary's heliocentric orbit (one
particular Sun-binary geometry angle `phi0`). This keeps the augmented system
AUTONOMOUS: a uniform, velocity-independent force is still derivable from a
time-independent potential in the rotating frame, so a Jacobi-LIKE conserved
quantity still exists (derived below), and the existing
`correct_symmetric_fixed_jacobi`-style fixed-Jacobi differential correction
still applies almost unchanged -- only the effective potential gains a linear
term.

IMPORTANT constraint on ``phi0``: only 0 or pi are valid
----------------------------------------------------------
The half-period "symmetric orbit" shooting trick this whole corrector
architecture is built on (Ross & Roberts-Tsoukkas 2025, AAS 25-621, Eqs.
9-12; `cr3bp_periodic.correct_symmetric_fixed_jacobi`'s own docstring) relies
on the standard CR3BP's MIRROR symmetry: if ``(x(t), y(t), vx(t), vy(t))`` is
a solution, so is ``(x(-t), -y(-t), vx(-t), -vy(-t))`` -- which is what lets a
single half-period crossing condition (``y=0, xdot=0``) certify FULL-period
closure without ever propagating the second half. A constant SRP force
``beta_nd*(cos(phi0), sin(phi0))`` breaks that symmetry UNLESS
``sin(phi0)=0`` (the y-component of the added force must vanish, since an
added CONSTANT does not flip sign under ``y -> -y`` the way the gravity terms
do). For any other ``phi0`` the corrector still numerically "converges" (the
half-period crossing residual goes to zero) but the result is NOT a genuine
periodic orbit -- confirmed directly: a full-period independent-Radau
crosscheck (:func:`crosscheck_periodic_srp`) closes to ~1e-8 at
``phi0 in {0, pi}`` (matching the beta=0 gravity-only baseline) but only to
~1e-3--1e-4 (four to five orders of magnitude worse) at
``phi0 in {0.3, pi/2}`` in a direct numerical check -- a real, load-bearing
finding, not a hypothetical caveat. **Every caller of
:func:`correct_symmetric_fixed_jacobi_srp`/`beta_step_to_target`/
`c_sweep_find_nu_zero_srp` (and `real_binary_kk_sweep.sweep_family_srp`/
`sweep_family_grid_srp`) MUST restrict ``phi0`` to ``{0.0, math.pi}``** --
SRP exactly along the primary-secondary axis (pointing from primary toward
secondary, or the reverse) -- for the returned orbit to be trustworthy; the
``crosscheck_ok``/``crosscheck_dj`` fields on the returned `SweepResult`
remain the authoritative, always-checked guard against silently trusting a
non-periodic "convergence" at any other angle. A fully general (non-mirror-
symmetric) SRP direction would need a full-period (not half-period) multiple-
shooting corrector, out of this task's scope.

Derivation of the augmented Jacobi integral
--------------------------------------------
Standard CR3BP pseudo-potential (rotating frame): ``Ubar(x,y,z) =
-1/2(x^2+y^2) - (1-mu)/r1 - mu/r2``, with ``a = -grad(Ubar) + 2*(ydot,-xdot,0)``
(Coriolis) and Jacobi ``C = -2*Ubar - v^2`` conserved (`cr3bp.jacobi_constant`).
A constant SRP force ``f_srp = beta_nd*(cos(phi0), sin(phi0), 0)`` (ND units)
is exactly ``-grad(U_srp)`` for ``U_srp(x,y) = -beta_nd*(x*cos(phi0) +
y*sin(phi0))`` -- a time-independent addition to the rotating-frame potential,
so the SAME derivation that gives the ordinary Jacobi integral (Coriolis does
no work; every other force is conservative in the rotating frame) gives:

    Ubar_tot = Ubar + U_srp
    C_srp    = -2*Ubar_tot - v^2 = C_std + 2*beta_nd*(x*cos(phi0) + y*sin(phi0))

conserved along trajectories of :func:`cr3bp_srp_eom`. This is checked
directly (a cheap, rigorous, paper-independent internal-consistency test) in
`tests/search/test_665_real_binary_srp.py` and again in the mandatory
positive-control script `scripts/run_665_srp_positive_control.py`.

Because the added force is CONSTANT (independent of position and velocity),
the Jacobian (A-matrix) of the equations of motion is IDENTICAL to the
gravity-only case (`cr3bp.cr3bp_stm_eom`'s own A-matrix derivation is
unaffected -- ``d(a_srp)/d(state) = 0`` identically). :func:`cr3bp_srp_stm_eom`
therefore reuses `cr3bp.cr3bp_stm_eom`'s variational block verbatim, only
splicing in the SRP-augmented state derivative.

Cannonball SRP magnitude
-------------------------
Standard cannonball model (e.g. Montenbruck & Gill 2000, *Satellite Orbits*,
Springer, Sec. 3.4 -- widely cited baseline in astrodynamics texts):
``a_srp = P_1AU * C_R * (A/m) / d_AU^2`` where ``P_1AU`` is the solar
radiation pressure at 1 AU for a fully-absorbing surface (``P_1AU =
Phi_1AU/c``, with the ~1367 W/m^2 solar constant and ``c`` the speed of
light -- ``P_1AU = 4.56e-6 N/m^2``, the commonly-quoted textbook value),
``C_R`` is the dimensionless reflectivity coefficient (1.0 fully absorbing,
2.0 fully specularly reflecting; bare rock/regolith is intermediate, ``C_R ~
1.0-1.4`` is a standard assumed range for an uncoated rocky/metallic body),
and ``A/m`` is the area-to-mass ratio (``beta`` in mission-design usage,
m^2/kg).

Bare-rock beta range (user-decided 2026-07-21, see `#665`'s own
`data/OUTSTANDING.md` bullet): ``1e-4`` to ``1e-3`` m^2/kg -- the range
physically defensible for an uncoated rocky/metallic cycler spacecraft body
itself, NOT a sail/balloon-augmented design. Sweeping outside this range
without separately flagging it as a more speculative extension is explicitly
out of this task's scope.
"""

from __future__ import annotations

import math
from typing import Final

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.binary_star_search import Topology, collinear_lpoints

__all__ = [
    "SOLAR_RADIATION_PRESSURE_1AU_N_M2",
    "barden_stability_srp",
    "beta_step_to_target",
    "c_sweep_find_nu_zero_srp",
    "cannonball_srp_accel_m_s2",
    "correct_symmetric_fixed_jacobi_srp",
    "cr3bp_srp_eom",
    "cr3bp_srp_stm_eom",
    "crosscheck_periodic_srp",
    "jacobi_constant_srp",
    "min_body_clearance_km_srp",
    "srp_beta_nd",
    "winding_topology_srp",
    "ydot0_from_jacobi_srp",
]

#: Solar radiation pressure at 1 AU for a fully-absorbing (C_R=1) surface,
#: N/m^2. Standard textbook cannonball-SRP baseline value (Montenbruck & Gill
#: 2000, *Satellite Orbits*, Springer, Sec. 3.4; equivalently
#: Phi_1AU/c = 1367 W/m^2 / 2.998e8 m/s = 4.56e-6 N/m^2, widely reproduced
#: across the astrodynamics literature -- not tied to any one paywalled
#: source).
SOLAR_RADIATION_PRESSURE_1AU_N_M2: Final[float] = 4.56e-6


def cannonball_srp_accel_m_s2(
    area_to_mass_m2_per_kg: float, c_r: float, sun_distance_au: float
) -> float:
    """Cannonball SRP acceleration magnitude, m/s^2, on a body at
    ``sun_distance_au`` AU from the Sun.

    ``a_srp = P_1AU * C_R * (A/m) / d_AU^2`` -- see module docstring.
    """
    return (
        SOLAR_RADIATION_PRESSURE_1AU_N_M2
        * c_r
        * area_to_mass_m2_per_kg
        / (sun_distance_au * sun_distance_au)
    )


def srp_beta_nd(
    area_to_mass_m2_per_kg: float,
    c_r: float,
    sun_distance_au: float,
    system: cr3bp.CR3BPSystem,
) -> float:
    """Convert a physical area-to-mass ratio into the nondimensional CR3BP
    perturbing-acceleration parameter ``beta_nd`` used by :func:`cr3bp_srp_eom`
    (ND acceleration unit = ``system.l_km*1000 / system.t_s**2`` m/s^2, by the
    same length/time-unit convention `cr3bp.CR3BPSystem` already uses)."""
    a_srp = cannonball_srp_accel_m_s2(area_to_mass_m2_per_kg, c_r, sun_distance_au)
    l_m = system.l_km * 1000.0
    return a_srp * system.t_s * system.t_s / l_m


# ---------------------------------------------------------------------------
# SRP-augmented dynamics core (constant-direction, autonomous-frame model)
# ---------------------------------------------------------------------------


def cr3bp_srp_eom(
    t: float, state6: NDArray[np.float64], mu: float, beta_nd: float, phi0: float
) -> NDArray[np.float64]:
    """CR3BP EOM plus a constant SRP acceleration ``beta_nd*(cos(phi0),
    sin(phi0), 0)`` (ND, rotating frame). Reduces exactly to
    :func:`cr3bp.cr3bp_eom` at ``beta_nd=0``."""
    base = cr3bp.cr3bp_eom(t, state6, mu)
    out = base.copy()
    out[3] += beta_nd * math.cos(phi0)
    out[4] += beta_nd * math.sin(phi0)
    return out


def cr3bp_srp_stm_eom(
    t: float, y42: NDArray[np.float64], mu: float, beta_nd: float, phi0: float
) -> NDArray[np.float64]:
    """State (6) + flattened STM (36) variational EOM under
    :func:`cr3bp_srp_eom`.

    The added SRP force is CONSTANT (independent of state), so its
    contribution to the Jacobian (A-matrix) of the RHS is exactly zero --
    the variational block is IDENTICAL to `cr3bp.cr3bp_stm_eom`'s own; only
    the state-derivative part changes (see module docstring)."""
    base_full = cr3bp.cr3bp_stm_eom(t, y42, mu)
    s = y42[:6]
    ds_srp = cr3bp_srp_eom(t, s, mu, beta_nd, phi0)
    return np.concatenate([ds_srp, base_full[6:]])


def jacobi_constant_srp(
    state6: NDArray[np.float64], mu: float, beta_nd: float, phi0: float
) -> float:
    """Augmented (SRP-inclusive) Jacobi-like conserved quantity -- see module
    docstring for the derivation. Reduces exactly to `cr3bp.jacobi_constant`
    at ``beta_nd=0``."""
    c0 = cr3bp.jacobi_constant(state6, mu)
    x, y = float(state6[0]), float(state6[1])
    return c0 + 2.0 * beta_nd * (x * math.cos(phi0) + y * math.sin(phi0))


def _ubar_x_at_axis_srp(x0: float, mu: float, beta_nd: float, phi0: float) -> float:
    """SRP-augmented analog of `cr3bp_periodic._ubar_x_at_axis`: ``-2*Ubar_tot``
    evaluated on the x-axis (y=0)."""
    return cp._ubar_x_at_axis(x0, mu) + 2.0 * beta_nd * x0 * math.cos(phi0)


def ydot0_from_jacobi_srp(
    x0: float, jacobi: float, mu: float, beta_nd: float, phi0: float, *, sign: float = 1.0
) -> float:
    """SRP-augmented analog of `cr3bp_periodic.ydot0_from_jacobi`: solve
    ``C_srp(x0,0,0,0,ydot0,0) = jacobi`` for ``ydot0``."""
    rad = _ubar_x_at_axis_srp(x0, mu, beta_nd, phi0) - jacobi
    if rad < 0.0:
        raise ValueError(
            f"ydot0_from_jacobi_srp: negative radicand {rad:.3e} at x0={x0:.6f}, "
            f"C_srp={jacobi:.6f}, beta_nd={beta_nd:.3e}, phi0={phi0:.4f}"
        )
    return float(sign) * math.sqrt(rad)


def _ubar_grad_x_at_axis_srp(x0: float, mu: float, beta_nd: float, phi0: float) -> float:
    """SRP-augmented analog of `cr3bp_periodic._ubar_grad_x_at_axis`:
    ``dUbar_tot/dx`` at ``(x0,0,0)``."""
    return cp._ubar_grad_x_at_axis(x0, mu) - beta_nd * math.cos(phi0)


def _xaxis_crossings_srp(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t_hi: float,
    beta_nd: float,
    phi0: float,
    *,
    with_stm: bool,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], list[NDArray[np.float64]]]:
    """SRP-augmented analog of `cr3bp_periodic._xaxis_crossings`."""

    def _y_event(t: float, y: NDArray[np.float64], *_args: float) -> float:
        return float(y[1])

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]

    rhs = cr3bp_srp_stm_eom if with_stm else cr3bp_srp_eom
    if with_stm:
        y0 = np.concatenate([np.asarray(state0, float), np.eye(6).reshape(36)])
    else:
        y0 = np.asarray(state0, float)

    sol = solve_ivp(
        rhs,
        (0.0, t_hi),
        y0,
        args=(system.mu, beta_nd, phi0),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events: list[NDArray[np.float64]] = list(sol.y_events[0]) if sol.y_events is not None else []
    t_lo = 1e-6 * t_hi
    pairs = [(t, y) for t, y in zip(t_events, y_events, strict=True) if t > t_lo]
    if not pairs:
        return np.array([]), []
    times = np.array([t for t, _ in pairs])
    states = [y for _, y in pairs]
    return times, states


def correct_symmetric_fixed_jacobi_srp(
    system: cr3bp.CR3BPSystem,
    x0_guess: float,
    jacobi: float,
    period_guess: float,
    beta_nd: float,
    phi0: float,
    *,
    ydot0_sign: float = 1.0,
    half_crossings: int | None = None,
    tol: float = 1e-8,
    max_iter: int = 30,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    x0_bounds: tuple[float, float] = (-2.0, 2.0),
) -> cp.SymmetricOrbit:
    """SRP-augmented analog of
    `cr3bp_periodic.correct_symmetric_fixed_jacobi` -- identical fixed-Jacobi
    single-shooting Newton algorithm (Ross & Roberts-Tsoukkas 2025, AAS 25-621,
    Eqs. 9-12), generalized with the SRP-augmented Jacobi/EOM/STM helpers
    above. At ``beta_nd=0`` this reduces to calling the gravity-only
    corrector (checked directly in the test suite).

    Raises
    ------
    ValueError
        If ``phi0`` is not (to within ``1e-9``) an integer multiple of pi --
        see the module docstring's "IMPORTANT constraint on phi0" section:
        the half-period symmetric-orbit trick this corrector relies on is
        only valid when the SRP force has zero y-component in the rotating
        frame.
    """
    if abs(math.sin(phi0)) > 1e-9:
        raise ValueError(
            f"correct_symmetric_fixed_jacobi_srp: phi0={phi0!r} has sin(phi0)="
            f"{math.sin(phi0):.3e} != 0 -- the half-period symmetric-orbit corrector is only "
            "valid for phi0 in {0, pi} (SRP along the primary-secondary axis); see the module "
            "docstring for why any other angle silently breaks full-period closure"
        )
    x0 = float(x0_guess)
    t_half_guess = 0.5 * float(period_guess)
    t_hi = 1.25 * float(period_guess)
    n_target = half_crossings
    crossing_res = float("inf")
    n_iter = 0
    t_half = t_half_guess
    yf: NDArray[np.float64] = np.zeros(6)
    lo, hi = x0_bounds
    for n_iter in range(1, max_iter + 1):  # noqa: B007
        ydot0 = ydot0_from_jacobi_srp(x0, jacobi, system.mu, beta_nd, phi0, sign=ydot0_sign)
        state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
        times, states = _xaxis_crossings_srp(
            system, state0, t_hi, beta_nd, phi0, with_stm=True, rtol=rtol, atol=atol
        )
        if len(times) == 0:
            break
        if n_target is None:
            n_target = int(np.argmin(np.abs(times - t_half_guess))) + 1
        idx = int(np.argmin(np.abs(times - t_half))) if n_target > len(times) else n_target - 1
        ystm = states[idx]
        yf = ystm[:6]
        stm = ystm[6:].reshape(6, 6)
        t_half = float(times[idx])
        xdot1 = float(yf[3])
        crossing_res = abs(xdot1)
        if crossing_res < tol:
            break
        ydot1 = float(yf[4])
        f1 = cr3bp_srp_eom(t_half, yf, system.mu, beta_nd, phi0)
        xddot1 = float(f1[3])
        dydot0_dx0 = -_ubar_grad_x_at_axis_srp(x0, system.mu, beta_nd, phi0) / ydot0
        dy_dx0 = float(stm[1, 0]) + float(stm[1, 4]) * dydot0_dx0
        dxdot_dx0 = float(stm[3, 0]) + float(stm[3, 4]) * dydot0_dx0
        if abs(ydot1) < 1e-14:
            break
        dxdot_total = dxdot_dx0 - xddot1 * dy_dx0 / ydot1
        if abs(dxdot_total) < 1e-14:
            break
        dx0 = -xdot1 / dxdot_total
        max_step = 0.2
        if abs(dx0) > max_step:
            dx0 = math.copysign(max_step, dx0)
        x0 = x0 + dx0
        x0 = min(max(x0, lo), hi)
    ydot0_final = ydot0_from_jacobi_srp(x0, jacobi, system.mu, beta_nd, phi0, sign=ydot0_sign)
    period = 2.0 * t_half
    converged = crossing_res < tol
    return cp.SymmetricOrbit(
        x0=x0,
        ydot0=ydot0_final,
        jacobi=jacobi_constant_srp(
            np.array([x0, 0.0, 0.0, 0.0, ydot0_final, 0.0]), system.mu, beta_nd, phi0
        ),
        t_half=t_half,
        period=period,
        converged=converged,
        crossing_residual=crossing_res,
        n_iter=n_iter,
    )


def barden_stability_srp(
    system: cr3bp.CR3BPSystem,
    orbit: cp.SymmetricOrbit,
    beta_nd: float,
    phi0: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, complex]:
    """SRP-augmented analog of `cr3bp_periodic.barden_stability` -- identical
    half-period-monodromy Barden stability parameter, propagated under
    :func:`cr3bp_srp_stm_eom`."""
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    y0 = np.concatenate([state0, np.eye(6).reshape(36)])
    sol = solve_ivp(
        cr3bp_srp_stm_eom,
        (0.0, orbit.t_half),
        y0,
        args=(system.mu, beta_nd, phi0),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"barden_stability_srp: propagation failed: {sol.message}")
    phi_half = sol.y[6:, -1].reshape(6, 6)
    idx = [0, 1, 3, 4]
    phi4 = phi_half[np.ix_(idx, idx)]
    g4 = np.diag([1.0, -1.0, -1.0, 1.0])
    monodromy = g4 @ np.linalg.inv(phi4) @ g4 @ phi4
    eigs = np.linalg.eigvals(monodromy)
    order = np.argsort(np.abs(eigs - 1.0))
    nontrivial = eigs[order[2:]]
    lam = complex(nontrivial[np.argmax(np.abs(nontrivial))])
    nu = 0.5 * (lam + 1.0 / lam)
    return float(nu.real), lam


def winding_topology_srp(
    mu: float,
    state0: NDArray[np.float64],
    period: float,
    beta_nd: float,
    phi0: float,
    *,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    n_samples: int = 4000,
) -> Topology:
    """SRP-augmented analog of `binary_star_search.winding_topology`."""
    sol = solve_ivp(
        cr3bp_srp_eom,
        (0.0, period),
        np.asarray(state0, float),
        args=(mu, beta_nd, phi0),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=period / n_samples,
    )
    x, y = sol.y[0], sol.y[1]

    def wind(px: float, py: float) -> float:
        th = np.unwrap(np.arctan2(y - py, x - px))
        return float((th[-1] - th[0]) / (2.0 * np.pi))

    w1 = wind(-mu, 0.0)
    w2 = wind(1.0 - mu, 0.0)
    l1, _l2, _l3 = collinear_lpoints(mu)
    return Topology(
        k1=round(abs(w1)),
        k2=round(abs(w2)),
        w1=w1,
        w2=w2,
        prograde=(w1 > 0.0 and w2 > 0.0),
        x_min=float(x.min()),
        x_max=float(x.max()),
        reaches_secondary=bool(x.max() > l1),
    )


def crosscheck_periodic_srp(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    beta_nd: float,
    phi0: float,
    *,
    method: str = "Radau",
    rtol: float = 1e-11,
    atol: float = 1e-11,
    closure_tol: float = 1e-6,
    jacobi_tol: float = 1e-8,
) -> tuple[bool, float]:
    """SRP-augmented analog of `cr3bp_periodic.crosscheck_periodic`: an
    independent (different-integrator) full-period closure + Jacobi_srp
    preservation check."""
    s0 = np.asarray(state0, float)
    c0 = jacobi_constant_srp(s0, system.mu, beta_nd, phi0)
    sol = solve_ivp(
        cr3bp_srp_eom,
        (0.0, period),
        s0,
        args=(system.mu, beta_nd, phi0),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    sf = sol.y[:, -1]
    closure = float(np.linalg.norm(sf - s0))
    dj = abs(jacobi_constant_srp(sf, system.mu, beta_nd, phi0) - c0)
    ok = closure < closure_tol and dj < jacobi_tol
    return ok, dj


def min_body_clearance_km_srp(
    target_system: cr3bp.CR3BPSystem,
    x0: float,
    ydot0: float,
    period: float,
    beta_nd: float,
    phi0: float,
    *,
    n_samples: int = 4000,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> tuple[float, float]:
    """SRP-augmented analog of
    `real_binary_kk_sweep.min_body_clearance_km` -- #660's physical
    body-clearance gate still applies unchanged under SRP (SRP does not
    relax the physical-collision problem)."""
    state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    sol = solve_ivp(
        cr3bp_srp_eom,
        (0.0, period),
        state0,
        args=(target_system.mu, beta_nd, phi0),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=period / n_samples,
    )
    x, y = sol.y[0], sol.y[1]
    d_primary_nd = float(np.hypot(x - (-target_system.mu), y).min())
    d_secondary_nd = float(np.hypot(x - (1.0 - target_system.mu), y).min())
    return d_primary_nd * target_system.l_km, d_secondary_nd * target_system.l_km


# ---------------------------------------------------------------------------
# Beta-continuation (mirrors mu_step_to_system's structure) and C_srp-sweep
# (mirrors pluto_charon_kk_sweep.c_sweep_find_nu_zero)
# ---------------------------------------------------------------------------


def beta_step_to_target(
    system: cr3bp.CR3BPSystem,
    x0_seed: float,
    jacobi_seed: float,
    period_seed: float,
    beta_target: float,
    phi0: float,
    *,
    hc: int | None,
    sign: float = -1.0,
    n_steps: int = 40,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> cp.SymmetricOrbit | None:
    """Step ``beta_nd`` from 0 (the gravity-only seed) to ``beta_target``,
    holding the (SRP-augmented) Jacobi constant fixed at each step -- the
    exact continuation pattern `real_binary_kk_sweep.mu_step_to_system` uses
    for mu, applied here to beta instead. The gravity-only seed's ordinary
    Jacobi IS `jacobi_constant_srp(..., beta_nd=0, ...)` (the two formulas
    coincide at beta=0), so this starts from the seed's own gravity-only
    corrector result unchanged.

    Returns ``None`` if any step fails to converge (a genuine continuation
    branch fold, the same failure mode `mu_step_to_system` reports)."""
    betas = np.linspace(0.0, beta_target, n_steps + 1)[1:]
    x0_cur, t_cur = x0_seed, period_seed
    jacobi = jacobi_seed  # == jacobi_constant_srp(..., beta_nd=0, phi0) by construction
    for beta_next in betas:
        try:
            o = correct_symmetric_fixed_jacobi_srp(
                system,
                x0_cur,
                jacobi,
                t_cur,
                float(beta_next),
                phi0,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=tol,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            return None
        if not o.converged:
            return None
        x0_cur, t_cur = o.x0, o.period
    return correct_symmetric_fixed_jacobi_srp(
        system,
        x0_cur,
        jacobi,
        t_cur,
        beta_target,
        phi0,
        ydot0_sign=sign,
        half_crossings=hc,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )


def c_sweep_find_nu_zero_srp(
    system: cr3bp.CR3BPSystem,
    x0_start: float,
    jacobi_start: float,
    period_start: float,
    beta_nd: float,
    phi0: float,
    *,
    hc: int | None,
    sign: float = -1.0,
    c_lo: float,
    c_hi: float,
    n_coarse: int = 60,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    nu_tol: float = 1e-10,
) -> cp.SymmetricOrbit | None:
    """SRP-augmented analog of `pluto_charon_kk_sweep.c_sweep_find_nu_zero`:
    C_srp-sweep the branch at fixed ``beta_nd``/``phi0``, locate the first
    Barden ``nu``-sign-change bracket, and brentq to the ``nu=0`` midpoint.
    Returns ``None`` if no stable window exists in ``[c_lo, c_hi]``."""
    from scipy.optimize import brentq

    x0_cur, t_cur = x0_start, period_start
    if jacobi_start != c_lo:
        for c_walk in np.linspace(jacobi_start, c_lo, 20)[1:]:
            try:
                o = correct_symmetric_fixed_jacobi_srp(
                    system,
                    x0_cur,
                    c_walk,
                    t_cur,
                    beta_nd,
                    phi0,
                    ydot0_sign=sign,
                    half_crossings=hc,
                    tol=tol,
                    rtol=rtol,
                    atol=atol,
                )
            except ValueError:
                break
            if o.converged:
                x0_cur, t_cur = o.x0, o.period

    c_grid = np.linspace(c_lo, c_hi, n_coarse)
    nu_prev: float | None = None
    orbit_prev: cp.SymmetricOrbit | None = None
    bracket: tuple[float, float, float, float] | None = None
    x0_sweep, t_sweep = x0_cur, t_cur

    for i, c_val in enumerate(c_grid):
        try:
            o = correct_symmetric_fixed_jacobi_srp(
                system,
                x0_sweep,
                c_val,
                t_sweep,
                beta_nd,
                phi0,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=tol,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            nu_prev = None
            orbit_prev = None
            continue
        if not o.converged:
            nu_prev = None
            orbit_prev = None
            continue
        nu, _ = barden_stability_srp(system, o, beta_nd, phi0, rtol=rtol, atol=atol)
        if nu_prev is not None and nu_prev * nu < 0.0:
            bracket = (c_grid[i - 1], nu_prev, c_val, nu)
            break
        nu_prev = nu
        orbit_prev = o
        x0_sweep, t_sweep = o.x0, o.period

    if bracket is None:
        return None

    c_a, _nu_a, c_b, _nu_b = bracket
    assert orbit_prev is not None
    x0_brent, t_brent = orbit_prev.x0, orbit_prev.period

    def _nu_at(c_val_inner: float) -> float:
        try:
            o = correct_symmetric_fixed_jacobi_srp(
                system,
                x0_brent,
                c_val_inner,
                t_brent,
                beta_nd,
                phi0,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=1e-11,
                rtol=1e-13,
                atol=1e-13,
            )
        except ValueError:
            return float("nan")
        if not o.converged:
            return float("nan")
        nu_inner, _ = barden_stability_srp(system, o, beta_nd, phi0, rtol=1e-13, atol=1e-13)
        return float(nu_inner)

    try:
        c_mid = brentq(_nu_at, c_a, c_b, xtol=nu_tol, rtol=nu_tol, maxiter=60)
    except ValueError:
        return None

    try:
        o_mid = correct_symmetric_fixed_jacobi_srp(
            system,
            x0_brent,
            c_mid,
            t_brent,
            beta_nd,
            phi0,
            ydot0_sign=sign,
            half_crossings=hc,
            tol=1e-11,
            rtol=1e-13,
            atol=1e-13,
        )
    except ValueError:
        return None
    if not o_mid.converged:
        return None
    return o_mid
