"""Bicircular Restricted 4-Body Problem (BCR4BP) Sun-Earth-Moon synodic dynamics.

#292 Phase 1 (Track-A capability build). Extends the Earth-Moon CR3BP one
fidelity rung by adding a Sun perturbation. The Sun is treated as a 4th body in
circular motion around the Sun-Earth-Moon barycenter; the Earth-Moon barycenter
therefore moves on a small circle around that same barycenter. The model is the
**standard (incoherent) BCR4BP** of Simo / Jorba / Gomez, parameterised by the
four constants in the Andreu (1998) digest (also tabulated in Rosales-Jorba 2023
Table 3 and Gimeno-Jorba 2018 Table 3):

  - mu      Earth-Moon mass ratio              ~ 0.012150581600000
  - mu_S    Sun mass (Earth+Moon = 1 units)    ~ 328900.5423094043
  - a_S     Sun distance (EM-distance units)   ~ 388.8111430233511
  - omega_S Sun angular freq in synodic frame  ~ 0.925195985520347 rad/TU

NOTE on Andreu's coherent QBCP vs the standard BCR4BP
-----------------------------------------------------
Andreu's (1998) **Quasi-BiCircular Problem (QBCP)** is the *coherent* refinement
of the BCR4BP and uses eight 2pi-periodic Fourier functions alpha_i(theta_S)
that encode the slight deviation of the Sun-Earth-Moon primaries from a true
bicircular configuration. Those alpha tables are NOT in the in-repo digest -- per
the digest itself, harvesting them from Gimeno-Jorba 2018 Table 4 is a
documented *future* step. Phase 1 therefore implements the standard incoherent
BCR4BP (the four-constant model), which:

  * is one fidelity rung above CR3BP, as required;
  * reduces EXACTLY to CR3BP when ``mu_S -> 0`` (structural sanity test);
  * uses *only* digest-sourced constants (no fabricated alpha values);
  * has a small but nonzero residual against the QBCP POL1/POL2 substitute ICs
    (the O(eps^2) deviation noted in the digest section "The Sun's perturbation
    is O(eps^2)"). This is the EXPECTED model gap and is documented in the
    Phase 1 genome tests: the BCR4BP corrector starting from the POL1 IC
    re-converges to a NEARBY BCR4BP periodic orbit, not the QBCP POL1 itself.

A future QBCP module (separate task) would carry the alpha_i Fourier tables
once digested.

Frame conventions (match cr3bp.py)
----------------------------------
  * Earth-Moon synodic rotating frame, nondim units (l = EM distance, t = 1/n).
  * Primaries: Earth at (-mu, 0, 0), Moon at (1 - mu, 0, 0).
  * Synodic frame angular rate = 1 (nondim).
  * Sun position (synodic): (a_S cos(theta_S), a_S sin(theta_S), 0)
    with ``theta_S = theta_S0 + omega_S * t``. omega_S is the Sun's frequency
    *relative to the synodic frame* (so the synodic-frame Sun period is
    ``2*pi/omega_S ~ 6.79 TU ~ 29.5 d``, the lunar synodic month).
  * Standard incoherent direct+indirect Sun acceleration:

        a_Sun = -mu_S * (r - r_Sun) / |r - r_Sun|^3
                - mu_S * r_Sun / a_S^3

    The second term is the "indirect" piece accounting for the Earth-Moon
    barycenter's own acceleration toward the Sun (so that a test particle at
    r=0 feels zero net Sun force on average -- which is what the synodic
    frame demands). Equivalently this is the standard BCR4BP Hamiltonian with
    the Sun's distance treated as a *constant* a_S (not a function of time)
    inside the indirect-term normalisation.

References
----------
  * Andreu (1998) PhD thesis (filed as cyclers_pdf, not in-session readable).
  * Gimeno, Jorba et al. (2018), Frontiers AMS 4:32. Table 3 (parameters).
  * Rosales, Jorba et al. (2023), CeMDA 135:15. Table 3 (parameters), Table 4
    (POL1/POL2 dynamical substitute ICs for the QBCP -- used as the BCR4BP
    Phase 1 corrector SEEDS, not as closure goldens).
  * Simo, Gomez, Jorba, Masdemont (1995), in "The dynamical behaviour of our
    solar system" -- standard BCR4BP formulation.

Discipline
----------
  * Pure module: math/numpy/scipy only.
  * Does NOT modify cr3bp.py. Reduces to it exactly at mu_S = 0.
  * Sourced goldens only: the four constants trace to Rosales-Jorba 2023
    Table 3 (cross-checked vs Gimeno-Jorba 2018 to printed precision -- see
    the digest's parameter table).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp

# Module-level sourced constants from the Andreu digest, traceable to
# Rosales-Jorba (2023) Table 3 / Gimeno-Jorba (2018) Table 3. These values
# are PUBLISHED; do NOT replace them with values computed by our own code.
_ANDREU_MU_EM: float = 0.012150581600000
"""Earth-Moon mass ratio (Rosales-Jorba 2023 Table 3)."""

_ANDREU_MU_S: float = 328900.5423094043
"""Sun mass in Earth+Moon = 1 units (Rosales-Jorba 2023 Table 3)."""

_ANDREU_A_S: float = 388.8111430233511
"""Sun semi-major axis in EM-distance units (Rosales-Jorba 2023 Table 3)."""

_ANDREU_OMEGA_S: float = 0.925195985520347
"""Sun angular frequency in the EM synodic frame (Rosales-Jorba 2023 Table 3)."""


@dataclass(frozen=True)
class BCR4BPSystem:
    """Bicircular Restricted 4-Body Problem parameters (Sun-Earth-Moon synodic).

    Attributes
    ----------
    mu :
        Earth-Moon mass ratio (Moon mass / (Earth + Moon) mass).
    mu_sun :
        Sun mass in units where Earth + Moon = 1.
    a_sun_nondim :
        Sun semi-major axis in EM-distance units.
    omega_sun_nondim :
        Sun angular frequency in the synodic frame, rad/TU.
    theta_sun0 :
        Sun phase angle at t=0, rad. Default 0 means the Sun starts on the
        +x axis at t=0 in the synodic frame. Phase choice does not affect
        which orbits *exist* but does set the time anchor for the periodicity
        commensurability constraint.

    Notes
    -----
    Constructed via :func:`andreu_default` for the digest-sourced parameter
    set, or directly for sensitivity / continuation studies. Set ``mu_sun = 0``
    to recover CR3BP exactly (structural test). Attribute names are lower-
    case for ruff N815 compliance; in published formulas the Sun-related
    quantities are conventionally written with a subscript S (mu_S, a_S,
    omega_S) -- the digest uses the latter, here read mu_sun = mu_S etc.
    """

    mu: float
    mu_sun: float
    a_sun_nondim: float
    omega_sun_nondim: float
    theta_sun0: float = 0.0

    @property
    def sun_period_tu(self) -> float:
        """Synodic-frame Sun period = 2*pi / omega_S (TU)."""
        return 2.0 * math.pi / self.omega_sun_nondim


def andreu_default() -> BCR4BPSystem:
    """BCR4BP with the four digest-sourced Andreu / Rosales-Jorba constants.

    Returns the canonical Sun-Earth-Moon BCR4BP parameter set:

      mu       = 0.012150581600000
      mu_sun   = 328900.5423094043
      a_sun    = 388.8111430233511 (EM-distance units)
      omega_sun= 0.925195985520347 (rad / EM time unit)
      theta_sun0 = 0

    These are the Rosales-Jorba 2023 Table 3 values, cross-checked against
    Gimeno-Jorba 2018 Table 3 (Andreu's published 1998 constants).
    """
    return BCR4BPSystem(
        mu=_ANDREU_MU_EM,
        mu_sun=_ANDREU_MU_S,
        a_sun_nondim=_ANDREU_A_S,
        omega_sun_nondim=_ANDREU_OMEGA_S,
        theta_sun0=0.0,
    )


def _sun_position(t: float, system: BCR4BPSystem) -> tuple[float, float, float]:
    """Sun position in synodic frame at nondim time ``t``."""
    theta = system.theta_sun0 + system.omega_sun_nondim * t
    a_sun = system.a_sun_nondim
    return a_sun * math.cos(theta), a_sun * math.sin(theta), 0.0


def _sun_acceleration(
    x: float, y: float, z: float, t: float, system: BCR4BPSystem
) -> tuple[float, float, float]:
    """Sun direct + indirect acceleration on the test particle (synodic frame).

    Returns ``(ax_sun, ay_sun, az_sun)`` to be added to the CR3BP RHS.

    Direct term: ``-mu_sun * (r - r_Sun) / |r - r_Sun|^3``.
    Indirect term: ``-mu_sun * r_Sun / a_sun^3`` (compensates for the EM
    barycenter's own Sun-ward acceleration, since the synodic frame is
    EM-barycenter-centred).
    """
    if system.mu_sun == 0.0:
        return 0.0, 0.0, 0.0
    sx, sy, _ = _sun_position(t, system)  # sz = 0 (Sun on planar circle)
    dx = x - sx
    dy = y - sy
    dz = z  # - sz, but sz = 0
    d2 = dx * dx + dy * dy + dz * dz
    d3 = d2 * math.sqrt(d2)
    a_sun3 = system.a_sun_nondim**3
    ax = -system.mu_sun * dx / d3 - system.mu_sun * sx / a_sun3
    ay = -system.mu_sun * dy / d3 - system.mu_sun * sy / a_sun3
    # Sun is planar (z=0); indirect z-term vanishes. Direct contributes -mu_sun*z/d3.
    az = -system.mu_sun * dz / d3
    return ax, ay, az


def bcr4bp_eom(t: float, state6: NDArray[np.float64], system: BCR4BPSystem) -> NDArray[np.float64]:
    """BCR4BP equations of motion in the Sun-Earth-Moon synodic frame.

    ``state6 = (x, y, z, vx, vy, vz)``. Returns d(state)/dt.

    Reduces to :func:`cyclerfinder.core.cr3bp.cr3bp_eom` exactly when
    ``system.mu_S == 0`` (the indirect term also vanishes by inspection).
    """
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    # CR3BP portion (reused so the limit test is structural, not coincidental).
    base = cr3bp.cr3bp_eom(t, np.asarray(state6, dtype=np.float64), system.mu)
    ax_sun, ay_sun, az_sun = _sun_acceleration(x, y, z, t, system)
    return np.array(
        [vx, vy, vz, float(base[3]) + ax_sun, float(base[4]) + ay_sun, float(base[5]) + az_sun],
        dtype=np.float64,
    )


def _cr3bp_uxx_block(
    x: float, y: float, z: float, mu: float
) -> tuple[float, float, float, float, float, float]:
    """Reproduce the CR3BP pseudo-potential second-derivative entries.

    Returns ``(uxx, uyy, uzz, uxy, uxz, uyz)``. Pulled out from cr3bp_stm_eom
    so the BCR4BP STM can add Sun contributions without duplicating the CR3BP
    math. Matches cr3bp.py lines 110-117 exactly.
    """
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    r1c, r2c = r1**3, r2**3
    r1f, r2f = r1**5, r2**5
    om1 = 1.0 - mu
    uxx = (
        1 - om1 / r1c - mu / r2c + 3 * om1 * (x + mu) ** 2 / r1f + 3 * mu * (x - 1 + mu) ** 2 / r2f
    )
    uyy = 1 - om1 / r1c - mu / r2c + 3 * om1 * y * y / r1f + 3 * mu * y * y / r2f
    uzz = -om1 / r1c - mu / r2c + 3 * om1 * z * z / r1f + 3 * mu * z * z / r2f
    uxy = 3 * om1 * (x + mu) * y / r1f + 3 * mu * (x - 1 + mu) * y / r2f
    uxz = 3 * om1 * (x + mu) * z / r1f + 3 * mu * (x - 1 + mu) * z / r2f
    uyz = 3 * om1 * y * z / r1f + 3 * mu * y * z / r2f
    return uxx, uyy, uzz, uxy, uxz, uyz


def _sun_second_deriv_block(
    x: float, y: float, z: float, t: float, system: BCR4BPSystem
) -> tuple[float, float, float, float, float, float]:
    """Second derivatives of the Sun direct acceleration w.r.t. (x, y, z).

    The indirect term ``-mu_sun * r_Sun / a_sun^3`` is independent of (x, y, z)
    and contributes nothing to the variational Jacobian.

    The DIRECT term is ``a_dir = -mu_sun * (r - r_S) / d^3`` with
    ``d = |r - r_S|``. Its Jacobian w.r.t. (x, y, z) is the standard inverse-
    square-law variational matrix:

        d a_dir / d r = -mu_sun * (I/d^3 - 3 (r - r_S)(r - r_S)^T / d^5)

    The returned entries (sxx, syy, szz, sxy, sxz, syz) are ADDED to the CR3BP
    (uxx, ..., uyz) block when assembling the BCR4BP STM A matrix.
    """
    if system.mu_sun == 0.0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    sx, sy, _ = _sun_position(t, system)
    dx = x - sx
    dy = y - sy
    dz = z
    d2 = dx * dx + dy * dy + dz * dz
    d3 = d2 * math.sqrt(d2)
    d5 = d3 * d2
    inv_d3 = 1.0 / d3
    inv_d5 = 1.0 / d5
    # d a_x / d x = -mu_sun * (1/d^3 - 3 dx^2 / d^5)
    sxx = -system.mu_sun * (inv_d3 - 3.0 * dx * dx * inv_d5)
    syy = -system.mu_sun * (inv_d3 - 3.0 * dy * dy * inv_d5)
    szz = -system.mu_sun * (inv_d3 - 3.0 * dz * dz * inv_d5)
    sxy = -system.mu_sun * (-3.0 * dx * dy * inv_d5)
    sxz = -system.mu_sun * (-3.0 * dx * dz * inv_d5)
    syz = -system.mu_sun * (-3.0 * dy * dz * inv_d5)
    return sxx, syy, szz, sxy, sxz, syz


def bcr4bp_stm_eom(
    t: float, state_and_stm: NDArray[np.float64], system: BCR4BPSystem
) -> NDArray[np.float64]:
    """BCR4BP variational EOM: state (6) + flattened 6x6 STM (36).

    A matrix layout matches CR3BP (cr3bp_stm_eom, line 118 onward): identity
    on the upper-right 3x3 (positions -> velocities), Coriolis on the
    lower-middle, and pseudo-potential second derivatives on the lower-left.
    The Sun direct-term contribution adds to the (uxx, uyy, uzz, uxy, uxz,
    uyz) entries.

    Reduces to :func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom` exactly when
    ``mu_S = 0``.
    """
    s = state_and_stm[:6]
    phi = state_and_stm[6:].reshape(6, 6)
    x, y, z = float(s[0]), float(s[1]), float(s[2])
    uxx, uyy, uzz, uxy, uxz, uyz = _cr3bp_uxx_block(x, y, z, system.mu)
    sxx, syy, szz, sxy, sxz, syz = _sun_second_deriv_block(x, y, z, t, system)
    uxx += sxx
    uyy += syy
    uzz += szz
    uxy += sxy
    uxz += sxz
    uyz += syz
    mat_a = np.zeros((6, 6))
    mat_a[0:3, 3:6] = np.eye(3)
    mat_a[3, 0], mat_a[3, 1], mat_a[3, 2] = uxx, uxy, uxz
    mat_a[4, 0], mat_a[4, 1], mat_a[4, 2] = uxy, uyy, uyz
    mat_a[5, 0], mat_a[5, 1], mat_a[5, 2] = uxz, uyz, uzz
    mat_a[3, 4], mat_a[4, 3] = 2.0, -2.0  # Coriolis (frame-fixed)
    dphi = mat_a @ phi
    return np.concatenate([bcr4bp_eom(t, s, system), dphi.reshape(36)])


@dataclass(frozen=True)
class BCR4BPArc:
    """Result of a BCR4BP propagation. Mirrors :class:`cr3bp.CR3BPArc`."""

    state_f: NDArray[np.float64]
    stm: NDArray[np.float64] | None
    t: float


def propagate_bcr4bp(
    system: BCR4BPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    with_stm: bool = False,
    t0: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> BCR4BPArc:
    """Propagate a BCR4BP state (and optionally the STM) from ``t0`` to ``t0+t``.

    Mirrors :func:`cyclerfinder.core.cr3bp.propagate` (DOP853 integrator). The
    Sun phase advances along the integration so the absolute time origin
    ``t0`` matters; supply it when stitching arcs (e.g. inside a
    multi-shooter). Independent-integrator cross-checks (Radau) are invoked
    directly via ``scipy.integrate.solve_ivp`` by the genome-level corrector,
    not exposed as a parameter here.

    Raises
    ------
    RuntimeError
        If the integrator fails (e.g. close-approach to a primary driving the
        step size below floating-point resolution).
    """
    state0 = np.asarray(state6, dtype=np.float64)
    if with_stm:
        y0 = np.concatenate([state0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            bcr4bp_stm_eom,
            (t0, t0 + t),
            y0,
            args=(system,),
            rtol=rtol,
            atol=atol,
            method="DOP853",
        )
        if not sol.success:
            raise RuntimeError(f"BCR4BP STM propagation failed at t={sol.t[-1]}: {sol.message}")
        yf = sol.y[:, -1]
        return BCR4BPArc(state_f=yf[:6], stm=yf[6:].reshape(6, 6), t=float(sol.t[-1]) - t0)
    sol = solve_ivp(
        bcr4bp_eom,
        (t0, t0 + t),
        state0,
        args=(system,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
    )
    if not sol.success:
        raise RuntimeError(f"BCR4BP propagation failed at t={sol.t[-1]}: {sol.message}")
    return BCR4BPArc(state_f=sol.y[:, -1], stm=None, t=float(sol.t[-1]) - t0)


def sun_commensurate_period(omega_sun: float, n: int) -> float:
    """Smallest period commensurate with ``n`` Sun revolutions: T = 2*pi*n/omega_sun.

    Strict periodicity in the BCR4BP requires the orbit period to satisfy
    ``omega_sun * T = 2*pi * n`` for some positive integer ``n`` (the Sun must
    return to the same synodic-frame phase after one orbit period). Smaller
    ``n`` => longer orbit; n=1 gives one Sun-synodic-period revolution
    (T ~ 6.79 TU ~ 29.5 d in EM units).
    """
    if n <= 0:
        raise ValueError(f"sun_commensurate_period: n must be positive integer; got {n}")
    return 2.0 * math.pi * n / omega_sun
