"""Concentric Circular Restricted 4-Body Problem (CCR4BP): Jupiter-Europa-Ganymede.

#689 Stage B sub-task 1 of #686's CCR4BP discovery plan (see
``docs/notes/2026-07-22-686-nbody-discovery-strategy-pass.md`` section 3 and the
in-corpus digests of Kumar-Anderson-de la Llave-Gunter 2021, arXiv:2109.14815,
AAS 21-651, and Kumar-Anderson-de la Llave 2023, arXiv:2309.06073, AAS 23-397).

This is the STRUCTURAL analogue of :mod:`cyclerfinder.core.bcr4bp` (the
Sun-Earth-Moon BCR4BP), one fidelity rung above the Jupiter-Europa CR3BP: a
massless spacecraft under Jupiter + Europa (the base CR3BP) PLUS a periodic
gravitational perturbation from Ganymede, a fourth body on a concentric circular
coplanar orbit. The two moons exert NO mutual force (that is the CCR4BP
idealisation that makes the model tractable and TIME-PERIODIC in a rotating
frame -- exactly as BCR4BP's incoherent Sun forcing does).

Model (per the arXiv:2109.14815 / arXiv:2309.06073 digests)
-----------------------------------------------------------
  * Jupiter (m1) + Europa (m2) + Ganymede (m3) + massless spacecraft.
  * Europa and Ganymede on CONCENTRIC CIRCLES of radii r12, r13 about Jupiter,
    at mean motions n2, n3, with no mutual perturbation.
  * Mass ratios mu = m2/(m1+m2) (Europa), and the Ganymede-frame ratio
    mu3 = m3/(m1+m3) used for the Jupiter-Ganymede reduction.
  * KEY STRUCTURE (the tractability hinge, and the structural test this module
    implements): the CCR4BP is simultaneously a TIME-PERIODIC PERTURBATION of
    BOTH the Jupiter-Europa CRTBP AND the Jupiter-Ganymede CRTBP:
        - at mu_Ganymede -> 0 it reduces EXACTLY to the Jupiter-Europa PCRTBP
          (trivial, pointwise -- like bcr4bp's mu_S -> 0 -> CR3BP);
        - at mu_Europa -> 0 it reduces EXACTLY to the Jupiter-Ganymede PCRTBP,
          but ONLY after a frame/units reinterpretation (see below) -- the base
          frame is Europa-synchronised, so the Ganymede limit is not a naive
          pointwise mu -> 0.

Frame convention (base = Jupiter-Europa synodic, mirrors bcr4bp exactly)
------------------------------------------------------------------------
  * Jupiter-Europa synodic rotating frame, nondimensional units:
    length = r12 (Europa SMA about Jupiter), time = 1/n2, frame rate = 1.
    Jupiter at (-mu, 0, 0), Europa at (1-mu, 0, 0). Origin = Jupiter-Europa
    barycentre. This is IDENTICAL to :func:`cyclerfinder.core.cr3bp.cr3bp_eom`
    at the Jupiter-Europa mass ratio.
  * Ganymede is the periodic perturber, playing the exact structural role
    BCR4BP's Sun plays. In the rotating frame Ganymede sits on a circle of
    radius ``a_gan = r13/r12`` and its synodic angle advances LINEARLY:

        theta_gan(t) = theta_gan0 + omega_gan * t

    where ``omega_gan = n3/n2 - 1`` is Ganymede's synodic angular rate (its
    inertial rate n3 minus the frame rate n2, nondimensionalised by n2). With
    the outer moon slower than the inner, ``omega_gan < 0`` (Ganymede regresses
    in the Europa-synodic frame). The forcing is time-periodic with period
    ``2*pi/|omega_gan|`` -- Ganymede's synodic period -- REGARDLESS of exact
    commensurability, exactly as the Sun's synodic month makes bcr4bp periodic.
    (The Europa:Ganymede mean-motion ratio is physically ~2.014, near but not
    exactly the 2:1 Laplace resonance; see the docstring of
    :func:`jupiter_europa_ganymede_default`. The 2:1 idealisation is a one-value
    change of ``omega_gan`` -- the structure is identical either way.)
  * Standard incoherent direct + indirect Ganymede acceleration (same shape as
    bcr4bp's Sun term, with Ganymede's nondim mass ``mu_gan = m3/(m1+m2)``):

        a_gan = -mu_gan * (r - r_gan) / |r - r_gan|^3
                - mu_gan * r_gan / a_gan^3

    The second (indirect) term compensates for the Jupiter-Europa barycentre
    origin's own acceleration toward Ganymede, so a particle co-located with the
    barycentre feels zero net synodic-frame Ganymede force. This indirect term
    is precisely what makes the ``mu_Europa -> 0`` limit reduce to the
    Jupiter-Ganymede *barycentric* CR3BP (not a Jupiter-centred variant); see
    ``tests/core/test_ccr4bp.py`` for the explicit frame transformation.

Discipline (mirrors bcr4bp)
---------------------------
  * Pure module: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`.
  * Does NOT modify cr3bp.py. Reduces to it exactly at mu_gan = 0 (and, under
    the documented transformation, at mu = 0).
  * Sourced constants only: the default parameter set is derived by
    :func:`jupiter_europa_ganymede_default` from the SAME JPL SSD
    ``core.satellites`` registry that #688's ``genome/composed_moon_map.py``
    uses -- no independently re-sourced constant set (cross-checked exactly in
    the tests).
  * This task builds ONLY the EOM + STM + structural validation. The
    pseudospectral torus corrector, whisker/manifold machinery, and
    mesh-intersection heteroclinic search are LATER #686 Stage-B sub-tasks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES


@dataclass(frozen=True)
class CCR4BPSystem:
    """Concentric Circular Restricted 4-Body Problem parameters (Jupiter-Europa-Ganymede).

    All quantities are in the Jupiter-Europa synodic nondimensional units
    (length = Europa SMA r12, time = 1/n2, base frame rate = 1).

    Attributes
    ----------
    mu :
        Europa mass ratio m2/(m1+m2) -- the base Jupiter-Europa CR3BP mu.
    mu_gan :
        Ganymede's nondimensional gravitational mass, m3/(m1+m2) (i.e.
        GM_Ganymede in units where G(m1+m2) = 1). This is the coefficient of the
        direct + indirect Ganymede acceleration -- NOT the Ganymede mass RATIO
        m3/(m1+m3) (that is :attr:`mu3_reduction`). Set ``mu_gan = 0`` to recover
        the Jupiter-Europa CR3BP exactly (structural test).
    a_gan :
        Ganymede orbital radius in Europa-SMA units, r13/r12 (> 1).
    omega_gan :
        Ganymede synodic angular rate in the Europa-synodic frame, n3/n2 - 1
        (nondim; < 0 for the slower outer moon).
    theta_gan0 :
        Ganymede synodic phase at t = 0, rad. Default 0 puts Ganymede on the
        +x axis at t = 0. Sets the time anchor for the periodicity
        commensurability constraint (does not affect which orbits exist).

    Notes
    -----
    Constructed via :func:`jupiter_europa_ganymede_default` for the
    registry-sourced Galilean parameter set, or directly for sensitivity /
    idealised-commensurability (exact 2:1) studies. Attribute names are
    lower-case for ruff N815 compliance.
    """

    mu: float
    mu_gan: float
    a_gan: float
    omega_gan: float
    theta_gan0: float = 0.0

    @property
    def ganymede_synodic_period(self) -> float:
        """Synodic-frame Ganymede period = 2*pi / |omega_gan| (TU)."""
        return 2.0 * math.pi / abs(self.omega_gan)

    @property
    def mu3_reduction(self) -> float:
        """Jupiter-Ganymede mass ratio of the ``mu -> 0`` reduced system.

        In the Europa-mass-vanishing limit Jupiter's nondimensional mass is
        ``1 - mu -> 1`` and Ganymede's is ``mu_gan``, so the reduced
        Jupiter-Ganymede CR3BP has mass ratio ``mu_gan / (1 + mu_gan)``. This is
        the value the ``mu -> 0`` structural reduction test targets; it equals
        the physical Ganymede mass ratio m3/(m1+m3) to the model's precision.
        """
        return self.mu_gan / (1.0 + self.mu_gan)


def two_body_synodic_rate(mu: float, mu_gan: float, a_gan: float) -> float:
    """Ganymede's synodic angular rate ``omega_gan = n3/n2 - 1`` in the frame.

    Uses the physically-correct two-body mean motions (Kepler's third law with
    the actual nondimensional masses), NOT a geometric approximation:

        n2 = sqrt(G(m1+m2)/r12^3) = 1 (the frame rate, by unit definition),
        n3 = sqrt(G(m1+m3)/r13^3),
        => omega_gan = n3/n2 - 1 = sqrt((1 - mu + mu_gan) / a_gan^3) - 1

    where ``1 - mu + mu_gan = (m1+m3)/(m1+m2)`` in the nondimensional mass unit
    ``G(m1+m2) = 1``. This exact two-body form is what makes the ``mu -> 0``
    Jupiter-Ganymede reduction EXACT: the reduced system's time-rescaled frame
    rate ``k = 1 + omega_gan = n3/n2`` matches the Jupiter-Ganymede CR3BP's own
    two-body-consistent unit frame rate. (A geometric ``sqrt(1/a_gan^3) - 1``,
    ignoring the moon masses, would leave an O(m/M) ~ 3e-5 residual in that
    reduction -- adequate for #688's exterior-map SCREEN but not for an exact
    structural reduction test.)
    """
    return math.sqrt((1.0 - mu + mu_gan) / a_gan**3) - 1.0


def jupiter_europa_ganymede_default() -> CCR4BPSystem:
    """CCR4BP parameters for Jupiter-Europa-Ganymede from the JPL SSD registry.

    Reuses the SAME sourced constants as #688's ``genome/composed_moon_map.py``
    (the ``core.satellites`` ``SATELLITES`` / ``PRIMARIES`` gm_de440 registry),
    read by the identical formulae so the two modules cannot drift
    (cross-checked in ``tests/core/test_ccr4bp.py``):

      * mu     = GM_Europa   / (GM_Jupiter_sys + GM_Europa)   (== composed's Europa mu)
      * mu_gan = GM_Ganymede / (GM_Jupiter_sys + GM_Europa)
      * a_gan  = SMA_Ganymede / SMA_Europa                    (~1.595)
      * omega_gan = two-body synodic rate (see :func:`two_body_synodic_rate`)  (~ -0.5035)

    The physical Europa:Ganymede mean-motion ratio is ``n2/n3 ~ 2.014`` -- near
    but NOT exactly the 2:1 Laplace resonance. The model is time-periodic in the
    Europa-synodic frame at Ganymede's synodic period regardless; the exact-2:1
    idealisation (``omega_gan = -0.5``) is a later-task option, not imposed here,
    so this default stays faithful to the sourced physical constants.

    ``omega_gan`` uses the two-body mean motion (:func:`two_body_synodic_rate`),
    which differs from #688's screen-grade geometric ``sqrt(GM_J/r^3)`` value by
    ~1.3e-5 (the moons' own masses in Kepler's law); this is the physically
    correct rate and is what makes the ``mu -> 0`` reduction exact.

    ``GM_Jupiter_sys`` is the JPL *system* GM (Jupiter + all Galilean moons);
    the residual folding of the other moons' mass into the Jupiter term is the
    same <= 2.4e-4 approximation ``cr3bp_system`` / ``composed_moon_map`` already
    accept, and is documented there.
    """
    gm_j = PRIMARIES["Jupiter"]
    gm_e = SATELLITES["Europa"].mu_km3_s2
    gm_g = SATELLITES["Ganymede"].mu_km3_s2
    sma_e = SATELLITES["Europa"].sma_km
    sma_g = SATELLITES["Ganymede"].sma_km
    denom = gm_j + gm_e  # mass unit G(m1+m2), matches composed_moon_map's Europa-mu denominator
    mu = gm_e / denom
    mu_gan = gm_g / denom
    a_gan = sma_g / sma_e
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)


def _ganymede_position(t: float, system: CCR4BPSystem) -> tuple[float, float, float]:
    """Ganymede position in the Jupiter-Europa synodic frame at nondim time ``t``."""
    theta = system.theta_gan0 + system.omega_gan * t
    a = system.a_gan
    return a * math.cos(theta), a * math.sin(theta), 0.0


def _ganymede_acceleration(
    x: float, y: float, z: float, t: float, system: CCR4BPSystem
) -> tuple[float, float, float]:
    """Ganymede direct + indirect acceleration on the spacecraft (synodic frame).

    Direct term: ``-mu_gan * (r - r_gan) / |r - r_gan|^3``.
    Indirect term: ``-mu_gan * r_gan / a_gan^3`` (compensates for the
    Jupiter-Europa barycentre origin's own acceleration toward Ganymede).

    Structurally identical to bcr4bp's ``_sun_acceleration``; returns
    ``(ax, ay, az)`` to be ADDED to the CR3BP RHS.
    """
    if system.mu_gan == 0.0:
        return 0.0, 0.0, 0.0
    gx, gy, _ = _ganymede_position(t, system)  # gz = 0 (Ganymede on planar circle)
    dx = x - gx
    dy = y - gy
    dz = z  # - gz, gz = 0
    d2 = dx * dx + dy * dy + dz * dz
    d3 = d2 * math.sqrt(d2)
    a_gan3 = system.a_gan**3
    ax = -system.mu_gan * dx / d3 - system.mu_gan * gx / a_gan3
    ay = -system.mu_gan * dy / d3 - system.mu_gan * gy / a_gan3
    # Ganymede is planar (z=0); indirect z-term vanishes. Direct contributes -mu_gan*z/d3.
    az = -system.mu_gan * dz / d3
    return ax, ay, az


def ccr4bp_eom(t: float, state6: NDArray[np.float64], system: CCR4BPSystem) -> NDArray[np.float64]:
    """CCR4BP equations of motion in the Jupiter-Europa synodic frame.

    ``state6 = (x, y, z, vx, vy, vz)``. Returns d(state)/dt.

    Reduces to :func:`cyclerfinder.core.cr3bp.cr3bp_eom` at the Europa mass ratio
    exactly when ``system.mu_gan == 0`` (the indirect term also vanishes by
    inspection). Under the frame/units transformation documented in the module
    docstring it reduces to the Jupiter-Ganymede CR3BP when ``system.mu == 0``.
    """
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    # Jupiter-Europa CR3BP portion (reused so the mu_gan -> 0 limit is structural).
    base = cr3bp.cr3bp_eom(t, np.asarray(state6, dtype=np.float64), system.mu)
    ax_g, ay_g, az_g = _ganymede_acceleration(x, y, z, t, system)
    return np.array(
        [vx, vy, vz, float(base[3]) + ax_g, float(base[4]) + ay_g, float(base[5]) + az_g],
        dtype=np.float64,
    )


def _cr3bp_uxx_block(
    x: float, y: float, z: float, mu: float
) -> tuple[float, float, float, float, float, float]:
    """CR3BP pseudo-potential second-derivative entries ``(uxx, uyy, uzz, uxy, uxz, uyz)``.

    Replicated from :func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom` (and matching
    bcr4bp's own copy) so the CCR4BP STM can add Ganymede contributions without
    duplicating the CR3BP math.
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


def _ganymede_second_deriv_block(
    x: float, y: float, z: float, t: float, system: CCR4BPSystem
) -> tuple[float, float, float, float, float, float]:
    """Second derivatives of the Ganymede DIRECT acceleration w.r.t. (x, y, z).

    The indirect term ``-mu_gan * r_gan / a_gan^3`` is independent of (x, y, z)
    and contributes nothing to the variational Jacobian. The direct term's
    Jacobian is the standard inverse-square-law variational matrix

        d a_dir / d r = -mu_gan * (I/d^3 - 3 (r - r_gan)(r - r_gan)^T / d^5)

    evaluated at Ganymede's TIME-VARYING position r_gan(t) (so this block is not
    constant along a trajectory). Entries are ADDED to the CR3BP
    ``(uxx, ..., uyz)`` block when assembling the CCR4BP STM A matrix.
    """
    if system.mu_gan == 0.0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    gx, gy, _ = _ganymede_position(t, system)
    dx = x - gx
    dy = y - gy
    dz = z
    d2 = dx * dx + dy * dy + dz * dz
    d3 = d2 * math.sqrt(d2)
    d5 = d3 * d2
    inv_d3 = 1.0 / d3
    inv_d5 = 1.0 / d5
    sxx = -system.mu_gan * (inv_d3 - 3.0 * dx * dx * inv_d5)
    syy = -system.mu_gan * (inv_d3 - 3.0 * dy * dy * inv_d5)
    szz = -system.mu_gan * (inv_d3 - 3.0 * dz * dz * inv_d5)
    sxy = -system.mu_gan * (-3.0 * dx * dy * inv_d5)
    sxz = -system.mu_gan * (-3.0 * dx * dz * inv_d5)
    syz = -system.mu_gan * (-3.0 * dy * dz * inv_d5)
    return sxx, syy, szz, sxy, sxz, syz


def ccr4bp_stm_eom(
    t: float, state_and_stm: NDArray[np.float64], system: CCR4BPSystem
) -> NDArray[np.float64]:
    """CCR4BP variational EOM: state (6) + flattened 6x6 STM (36).

    A-matrix layout matches CR3BP / bcr4bp: identity on the upper-right 3x3
    (positions -> velocities), Coriolis on the lower-middle, and pseudo-potential
    second derivatives on the lower-left. The Ganymede direct-term contribution
    (evaluated at Ganymede's time-varying position) adds to the
    ``(uxx, uyy, uzz, uxy, uxz, uyz)`` entries.

    Reduces to :func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom` at the Europa mass
    ratio exactly when ``mu_gan = 0``.
    """
    s = state_and_stm[:6]
    phi = state_and_stm[6:].reshape(6, 6)
    x, y, z = float(s[0]), float(s[1]), float(s[2])
    uxx, uyy, uzz, uxy, uxz, uyz = _cr3bp_uxx_block(x, y, z, system.mu)
    gxx, gyy, gzz, gxy, gxz, gyz = _ganymede_second_deriv_block(x, y, z, t, system)
    uxx += gxx
    uyy += gyy
    uzz += gzz
    uxy += gxy
    uxz += gxz
    uyz += gyz
    mat_a = np.zeros((6, 6))
    mat_a[0:3, 3:6] = np.eye(3)
    mat_a[3, 0], mat_a[3, 1], mat_a[3, 2] = uxx, uxy, uxz
    mat_a[4, 0], mat_a[4, 1], mat_a[4, 2] = uxy, uyy, uyz
    mat_a[5, 0], mat_a[5, 1], mat_a[5, 2] = uxz, uyz, uzz
    mat_a[3, 4], mat_a[4, 3] = 2.0, -2.0  # Coriolis (frame-fixed)
    dphi = mat_a @ phi
    return np.concatenate([ccr4bp_eom(t, s, system), dphi.reshape(36)])


@dataclass(frozen=True)
class CCR4BPArc:
    """Result of a CCR4BP propagation. Mirrors :class:`cr3bp.CR3BPArc`."""

    state_f: NDArray[np.float64]
    stm: NDArray[np.float64] | None
    t: float


def propagate_ccr4bp(
    system: CCR4BPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    with_stm: bool = False,
    t0: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> CCR4BPArc:
    """Propagate a CCR4BP state (and optionally the STM) from ``t0`` to ``t0+t``.

    Mirrors :func:`cyclerfinder.core.bcr4bp.propagate_bcr4bp` (DOP853). Ganymede's
    synodic phase advances along the integration, so the absolute time origin
    ``t0`` matters; supply it when stitching arcs (e.g. inside a multi-shooter).

    Raises
    ------
    RuntimeError
        If the integrator fails (e.g. close approach to a primary driving the
        step size below floating-point resolution).
    """
    state0 = np.asarray(state6, dtype=np.float64)
    if with_stm:
        y0 = np.concatenate([state0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            ccr4bp_stm_eom,
            (t0, t0 + t),
            y0,
            args=(system,),
            rtol=rtol,
            atol=atol,
            method="DOP853",
        )
        if not sol.success:
            raise RuntimeError(f"CCR4BP STM propagation failed at t={sol.t[-1]}: {sol.message}")
        yf = sol.y[:, -1]
        return CCR4BPArc(state_f=yf[:6], stm=yf[6:].reshape(6, 6), t=float(sol.t[-1]) - t0)
    sol = solve_ivp(
        ccr4bp_eom,
        (t0, t0 + t),
        state0,
        args=(system,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
    )
    if not sol.success:
        raise RuntimeError(f"CCR4BP propagation failed at t={sol.t[-1]}: {sol.message}")
    return CCR4BPArc(state_f=sol.y[:, -1], stm=None, t=float(sol.t[-1]) - t0)


def ganymede_commensurate_period(omega_gan: float, n: int) -> float:
    """Smallest period commensurate with ``n`` Ganymede synodic revolutions.

    ``T = 2*pi*n / |omega_gan|``. Strict periodicity in the CCR4BP requires the
    orbit period to satisfy ``omega_gan * T = 2*pi * n`` for some positive
    integer ``n`` (Ganymede must return to the same synodic-frame phase after one
    orbit period). Mirrors bcr4bp's :func:`sun_commensurate_period`.
    """
    if n <= 0:
        raise ValueError(f"ganymede_commensurate_period: n must be positive integer; got {n}")
    return 2.0 * math.pi * n / abs(omega_gan)
