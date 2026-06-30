"""Ross-Scheeres (2007) Keplerian map for gravity-assist moon-tour design (#500).

Implements the 2D symplectic twist map on (omega, K) state space — periapsis
angle in the rotating frame and Keplerian energy — for a spacecraft on an
EXTERIOR orbit (periapsis outside the moon's orbital radius) in the Planar
Circular Restricted 3-Body Problem (PCR3BP).  The map advances the state by
one periapsis passage of the spacecraft without integrating the CR3BP ODEs,
making it an extremely fast surrogate for resonant gravity-assist dynamics.

The controlled version (Grover & Ross 2009) adds a scalar control input u_n
that shifts the periapsis angle before the kick is applied, modelling an
impulsive velocity change at periapsis.  A greedy coarse controller steers
the semimajor axis along a decreasing (or increasing) path through resonance
space with bounded total DeltaV.

References
----------
Ross, S. D. and Scheeres, D. J., "Multiple Gravity Assists, Capture, and
Escape in the Restricted Three-Body Problem," SIAM Journal on Applied Dynamical
Systems (SIADS), Vol. 6, No. 3, pp. 576-596, 2007.
DOI: 10.1137/06065195X  [hereafter RS07]

Grover, P. and Ross, S. D., "Designing Trajectories in a Planet-Moon
Environment Using the Controlled Keplerian Map," Journal of Guidance, Control,
and Dynamics, Vol. 32, No. 2, March-April 2009, pp. 436-443.
DOI: 10.2514/1.38320  [hereafter GR09]

Normalisation (CR3BP canonical units for a given planet-moon pair)
------------------------------------------------------------------
  - length unit  : moon orbital semi-major axis (= 1)
  - time unit    : moon orbital period / (2 pi)  (=> moon mean motion = 1)
  - GM_primary ~ 1  (exact when mu << 1)

In these units:
  K = -1/(2a)        Keplerian energy  (a in moon-radii)
  omega              argument of periapsis in the rotating frame [rad]
  T_sc = 2pi a^{3/2} spacecraft orbital period

Sign convention for omega (RS07)
---------------------------------
  omega = 0   : periapsis points TOWARD the moon
  omega > 0   : periapsis slightly AHEAD of the moon in its orbit
                => moon decelerates spacecraft => K decreases => a decreases
  omega < 0   : periapsis slightly BEHIND the moon
                => moon accelerates spacecraft => K increases => a increases
  f(omega) is ODD: f(-omega) = -f(omega); f(0) = 0.

Map equations (RS07 eq. 4.2)
------------------------------
  K_{n+1}     = K_n + mu * f(omega_n + u_n)       [kick, with optional control]
  omega_{n+1} = (omega_n + u_n) - 2pi*(-2*K_{n+1})^{-3/2}  (mod 2pi)

where f(omega) is the kick function (precomputed by Picard-first-order
quadrature along the unperturbed Keplerian orbit).

Validity
--------
The map is valid for the EXTERIOR periapsis regime (periapsis > moon orbit).
Validity requires: a(1-e) >= 1 (spacecraft periapsis outside the moon).
The approximation degrades near the Hill sphere (r_peri very close to 1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

# ---------------------------------------------------------------------------
# Physical-conversion helpers
# ---------------------------------------------------------------------------

_TWO_PI: float = 2.0 * math.pi


def eccentricity_from_tisserand(K: float, C_J: float) -> float:
    """Eccentricity e for a Keplerian orbit with energy K at Jacobi constant C_J.

    Uses the Tisserand (Jacobi-constant) approximation (valid for mu << 1):

        C_J ~= 1/a + 2*sqrt(a*(1-e^2))

    Returns NaN when no real eccentricity is consistent with K and C_J.

    Source: RS07 eq. 2.5 / Tisserand-parameter approximation.
    """
    if K >= 0.0:
        return math.nan  # unbound orbit
    a = -0.5 / K
    rhs = (C_J - 1.0 / a) / 2.0  # = sqrt(a*(1-e^2))
    if rhs < 0.0:
        return math.nan
    e_sq = 1.0 - rhs * rhs / a
    if e_sq < 0.0:
        return 0.0  # round-off; treat as circular
    return math.sqrt(e_sq)


def semimajor_from_K(K: float) -> float:
    """Semimajor axis from Keplerian energy K = -1/(2a)."""
    return -0.5 / K


def periapsis_radius(K: float, C_J: float) -> float:
    """Periapsis radius r_peri = a*(1-e) in moon-radius units."""
    a = semimajor_from_K(K)
    e = eccentricity_from_tisserand(K, C_J)
    if math.isnan(e):
        return math.nan
    return a * (1.0 - e)


# ---------------------------------------------------------------------------
# Kick-function: first-order Picard quadrature
# ---------------------------------------------------------------------------


def _eccentric_anomaly(theta: float, e: float) -> float:
    """Eccentric anomaly E from true anomaly theta.

    Uses atan2 form for stability at theta near +-pi.
    """
    return 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(theta / 2.0),
        math.sqrt(1.0 + e) * math.cos(theta / 2.0),
    )


def _theta_to_time(theta: float, e: float, n_mm: float) -> float:
    """Time from periapsis (theta=0) via Kepler's equation.

    Parameters
    ----------
    theta : true anomaly [rad], range (-pi, pi)
    e     : eccentricity
    n_mm  : spacecraft mean motion = a^{-3/2} = (-2K)^{3/2} [rad/time]
    """
    E = _eccentric_anomaly(theta, e)
    M = E - e * math.sin(E)  # mean anomaly
    return M / n_mm


def _kick_integrand(
    theta: float,
    omega: float,
    a: float,
    e: float,
    h: float,
    n_mm: float,
) -> float:
    """Integrand of the Picard energy-kick integral at true anomaly theta.

    Returns   (F_moon . v_sc) * (r^2/h)

    which when integrated over theta in [-pi, pi] gives f(omega) = DeltaK/mu.

    Physical picture
    ----------------
    The spacecraft follows its unperturbed Keplerian orbit.  At true anomaly
    theta, we evaluate the work done on the spacecraft by the moon's gravity
    (F_moon . v_sc), weighted by the Jacobian dt/dtheta = r^2/h.  Integrating
    over one full orbit (-pi to pi) gives the net energy change per passage.

    Coordinate system
    -----------------
    Inertial frame with the primary at the origin.  At t=0 (periapsis), the
    moon is at (1, 0) on the +x axis and the spacecraft periapsis is in the
    direction (cos omega, sin omega).  The moon moves counterclockwise on a
    unit circle: r_moon(t) = (cos t, sin t).
    """
    # Spacecraft radius and inertial position
    r = a * (1.0 - e * e) / (1.0 + e * math.cos(theta))
    phi = omega + theta  # spacecraft inertial angle
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    x_sc = r * cos_phi
    y_sc = r * sin_phi

    # Time from periapsis (Kepler's equation)
    t = _theta_to_time(theta, e, n_mm)

    # Moon inertial position (circular orbit at unit radius, passes (1,0) at t=0)
    x_m = math.cos(t)
    y_m = math.sin(t)

    # Spacecraft-to-moon vector and cubed distance
    dx = x_m - x_sc
    dy = y_m - y_sc
    r2_sq = dx * dx + dy * dy
    if r2_sq < 1.0e-24:
        return 0.0  # degenerate: at moon (should not happen for exterior map)
    r2 = math.sqrt(r2_sq)
    r2_3 = r2_sq * r2

    # Gravitational force per unit mu (dimensionless in normalised units)
    f_x = dx / r2_3
    f_y = dy / r2_3

    # Spacecraft inertial velocity for Keplerian orbit (GM = 1)
    #   v_radial    = (e sinθ) / h      (positive outward)
    #   v_transverse = h / r             (positive counterclockwise)
    v_r = e * math.sin(theta) / h
    v_t = h / r
    v_x = v_r * cos_phi - v_t * sin_phi
    v_y = v_r * sin_phi + v_t * cos_phi

    # Power per unit mu: F_moon . v_sc
    power = f_x * v_x + f_y * v_y

    # Jacobian: dt/dtheta = r^2/h  (Kepler's second law)
    return power * (r * r / h)


def compute_kick(omega: float, K: float, C_J: float, n_quad: int = 400) -> float:
    """Kick function f(omega) = DeltaK/mu for one exterior periapsis passage.

    Numerically integrates the Picard first-order energy-kick integral:

        f(omega) = ∫_{-pi}^{pi} (F_moon . v_sc) * (r^2/h) d(theta)

    along the unperturbed Keplerian orbit with Keplerian energy K and
    eccentricity derived from K and C_J via the Tisserand relation.

    Parameters
    ----------
    omega : periapsis angle in rotating frame [rad]
    K     : Keplerian energy (< 0) in CR3BP normalised units
    C_J   : Jacobi constant (fixes eccentricity given K)
    n_quad: scipy.integrate.quad 'limit' parameter (subdivisions)

    Returns
    -------
    f such that DeltaK_per_passage = mu * f(omega).
    """
    if K >= 0.0:
        return 0.0
    a = semimajor_from_K(K)
    e = eccentricity_from_tisserand(K, C_J)
    if math.isnan(e) or e >= 1.0:
        return 0.0
    h = math.sqrt(a * (1.0 - e * e))
    n_mm = (-2.0 * K) ** 1.5  # mean motion = a^{-3/2}

    def integrand(theta: float) -> float:
        return _kick_integrand(theta, omega, a, e, h, n_mm)

    # Split at theta = 0 (periapsis, near-peak when omega ~ 0) for robustness
    result, _err = quad(
        integrand,
        -math.pi,
        math.pi,
        limit=n_quad,
        points=[0.0],
        epsabs=1.0e-7,
        epsrel=1.0e-7,
    )
    return result


# ---------------------------------------------------------------------------
# KeplerianMap class
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MapState:
    """State of the spacecraft in the Keplerian map phase space."""

    omega: float  # periapsis angle in rotating frame [rad], in (-pi, pi]
    K: float  # Keplerian energy = -1/(2a), < 0
    C_J: float  # Jacobi constant (approximately conserved along a trajectory)

    @property
    def a(self) -> float:
        """Semimajor axis in moon-radius units."""
        return semimajor_from_K(self.K)

    @property
    def e(self) -> float:
        """Eccentricity (from Tisserand relation with C_J)."""
        return eccentricity_from_tisserand(self.K, self.C_J)

    @property
    def r_peri(self) -> float:
        """Periapsis radius in moon-radius units."""
        a = self.a
        e = self.e
        if math.isnan(e):
            return math.nan
        return a * (1.0 - e)

    @property
    def r_apo(self) -> float:
        """Apoapsis radius in moon-radius units."""
        a = self.a
        e = self.e
        if math.isnan(e):
            return math.nan
        return a * (1.0 + e)


class KeplerianMap:
    """Keplerian periapsis map (RS07 eq. 4.2) with optional coarse control (GR09).

    The map precomputes the kick function f(omega) on a grid of omega values
    (using first-order Picard quadrature) and uses a cubic spline for fast
    interpolation during orbit propagation.

    Parameters
    ----------
    mu : mass ratio of the perturbing moon (m_moon/(m_primary + m_moon))
    C_J : Jacobi constant (conserved along a trajectory; fixes e given a)
    K_ref : reference Keplerian energy for the kick table (sets e for quadrature)
    n_grid : number of omega grid points in (-pi, pi] for the kick table
    n_quad : scipy.integrate.quad limit per grid point (controls accuracy)
    """

    def __init__(
        self,
        mu: float,
        C_J: float,
        K_ref: float,
        n_grid: int = 101,
        n_quad: int = 300,
    ) -> None:
        self.mu = mu
        self.C_J = C_J
        self.K_ref = K_ref

        # Precompute f(omega) on a uniform grid in (-pi, pi]
        # Use n_grid+1 points from -pi to pi (inclusive) for periodicity
        omega_arr = np.linspace(-math.pi, math.pi, n_grid)
        kick_arr = np.array([compute_kick(w, K_ref, C_J, n_quad) for w in omega_arr])

        self._omega_grid: np.ndarray = omega_arr
        self._kick_grid: np.ndarray = kick_arr

        # Periodic cubic spline (f is 2pi-periodic)
        # Enforce f(-pi) = f(pi) exactly (they should be equal by periodicity)
        kick_arr[-1] = kick_arr[0]
        self._spline = CubicSpline(omega_arr, kick_arr, bc_type="periodic")

    # ------------------------------------------------------------------
    # Kick function
    # ------------------------------------------------------------------

    def kick(self, omega: float) -> float:
        """Interpolated kick function f(omega). Valid for omega in (-pi, pi]."""
        # Wrap to (-pi, pi]
        omega = ((omega + math.pi) % _TWO_PI) - math.pi
        return float(self._spline(omega))

    # ------------------------------------------------------------------
    # Map step
    # ------------------------------------------------------------------

    def step(
        self,
        omega: float,
        K: float,
        u: float = 0.0,
    ) -> tuple[float, float]:
        """One periapsis passage (RS07 eq. 4.2, with optional GR09 control).

        Map equations:
            K_{n+1}     = K_n + mu * f(omega_n + u_n)
            omega_{n+1} = (omega_n + u_n) - 2pi*(-2*K_{n+1})^{-3/2}  (mod 2pi)

        Control convention (GR09): u_n shifts the periapsis angle before the
        kick, equivalent to applying an impulsive DeltaV before closest approach.
        u = 0 is the uncontrolled map.

        Parameters
        ----------
        omega : current periapsis angle [rad]
        K     : current Keplerian energy
        u     : control input (angle shift, radians); positive u = shift ahead

        Returns
        -------
        (omega_new, K_new) after one periapsis passage.
        """
        omega_eff = omega + u  # apply control before kick

        # Energy kick
        f = self.kick(omega_eff)
        K_new = K + self.mu * f

        if K_new >= 0.0:
            # Orbit became unbound
            K_new = K  # freeze (caller should check)

        # Periapsis angle advance: rotating frame has moved by T_sc during one orbit
        T_sc = _TWO_PI * (-2.0 * K_new) ** (-1.5)
        omega_new = omega_eff - T_sc

        # Wrap to (-pi, pi]
        omega_new = ((omega_new + math.pi) % _TWO_PI) - math.pi

        return omega_new, K_new

    # ------------------------------------------------------------------
    # Trajectory propagation
    # ------------------------------------------------------------------

    def propagate(
        self,
        omega0: float,
        K0: float,
        n_steps: int,
        u_sequence: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Propagate the map for n_steps periapsis passages.

        Parameters
        ----------
        omega0, K0  : initial state
        n_steps     : number of passages
        u_sequence  : optional array of control inputs (length n_steps);
                      None => uncontrolled (u = 0 at every step)

        Returns
        -------
        (omegas, Ks) : arrays of shape (n_steps+1,) including the initial state
        """
        omegas = np.empty(n_steps + 1)
        Ks = np.empty(n_steps + 1)
        omegas[0] = omega0
        Ks[0] = K0

        for i in range(n_steps):
            u = 0.0 if u_sequence is None else float(u_sequence[i])
            omegas[i + 1], Ks[i + 1] = self.step(omegas[i], Ks[i], u)

        return omegas, Ks

    # ------------------------------------------------------------------
    # Greedy coarse controller (GR09 §III coarse-control algorithm)
    # ------------------------------------------------------------------

    def coarse_control(
        self,
        omega0: float,
        K0: float,
        K_target: float,
        u_max: float,
        n_max_steps: int = 200,
        lookahead: int = 5,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Greedy energy-descent coarse controller (GR09 §III).

        At each step, applies control u in {0, +u_max, -u_max} that best
        steers K toward K_target over the next `lookahead` steps (similar to
        the GR09 coarse-control look-ahead).

        Parameters
        ----------
        omega0, K0  : initial state
        K_target    : target Keplerian energy (< K0 for energy descent)
        u_max       : maximum control shift per step [rad]
        n_max_steps : maximum number of periapsis passages
        lookahead   : look-ahead horizon for control selection

        Returns
        -------
        (omegas, Ks, u_applied) : trajectory arrays (length n_max_steps+1 or shorter)
        """
        omegas = [omega0]
        Ks = [K0]
        u_applied = []

        omega = omega0
        K = K0

        for _ in range(n_max_steps):
            if abs(K - K_target) < 1e-6:
                break  # close enough

            # Try u = 0, +u_max, -u_max; pick the one that moves K most toward target
            best_u = 0.0
            best_score = -math.inf
            for u_cand in (0.0, u_max, -u_max):
                # Look ahead 'lookahead' steps with this control
                w, k = omega, K
                for _j in range(lookahead):
                    w, k = self.step(w, k, u_cand if _j == 0 else 0.0)
                    if k >= 0.0:
                        break
                # Score: net progress toward K_target
                score = K - k if K_target < K else k - K
                if score > best_score:
                    best_score = score
                    best_u = u_cand

            omega, K = self.step(omega, K, best_u)
            omegas.append(omega)
            Ks.append(K)
            u_applied.append(best_u)

            if K >= 0.0:
                break  # unbound

        return np.array(omegas), np.array(Ks), np.array(u_applied)


# ---------------------------------------------------------------------------
# Physical-unit helpers
# ---------------------------------------------------------------------------


def periapsis_velocity_norm(K: float, C_J: float) -> float:
    """Periapsis velocity in CR3BP normalised units (units of moon orbital speed).

    v_peri = h / r_peri where h = sqrt(a*(1-e^2)).
    """
    a = semimajor_from_K(K)
    e = eccentricity_from_tisserand(K, C_J)
    if math.isnan(e):
        return math.nan
    h = math.sqrt(a * (1.0 - e * e))
    r_p = a * (1.0 - e)
    return h / r_p


def control_dv_m_s(u_rad: float, K: float, C_J: float, v_moon_km_s: float) -> float:
    """Convert a control angle shift u_rad [rad] to physical DeltaV [m/s].

    The control u is approximated as a tangential impulse at periapsis:
        DeltaV [m/s] ~ |u_rad| * v_peri_norm * v_moon_km_s * 1000

    Parameters
    ----------
    u_rad       : control angle shift [rad]
    K, C_J      : current map state (to compute periapsis velocity)
    v_moon_km_s : moon orbital speed [km/s] (physical units conversion)
    """
    v_peri = periapsis_velocity_norm(K, C_J)
    return abs(u_rad) * v_peri * v_moon_km_s * 1.0e3  # km/s -> m/s
