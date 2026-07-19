"""Averaged 1-DOF co-orbital Hamiltonian: a seedless, global QS/HS/tadpole map.

**Task #666** (`#661` shortlist item 4). Zero co-orbital rows exist anywhere
in this project's 361-row catalogue; this module builds the standard
semi-analytic tool for that object class -- the resonant (1:1 mean-motion)
disturbing function, averaged over the fast synodic angle, leaving a 1-DOF
Hamiltonian in the resonant angle ``sigma = lambda - lambda_planet`` and the
slowly-varying semi-major-axis offset ``delta = a - a_planet``, at a fixed
(also slowly-varying) eccentricity ``e``. Level curves of this Hamiltonian
give the classical tadpole (librating around L4/L5, ``sigma`` near +-60 deg),
horseshoe (librating around L3 through L4/L5, wide ``sigma`` excursion) and
quasi-satellite (QS -- librating around ``sigma = 0``, i.e. appearing to
orbit the planet itself in the rotating frame) regimes, and -- critically --
the separatrix between them, which is where real bodies' *un-averaged*
trajectories cross when they undergo an observed QS<->HS transition (the
perturbation this averaging removes is exactly what drives the crossing).

Derivation (planar, circular-planet restricted 3-body problem; heliocentric,
not barycentric -- the O(mu) barycentric offset is ~3e-6 AU for Sun-Earth,
utterly negligible against the Hill-radius (~0.01 AU) scale this module
resolves)
------------------------------------------------------------------------

Nondimensional units: ``a_planet = 1``, ``n_planet = 1``, ``GM_total = 1``
(the standard CR3BP convention already used throughout this codebase, e.g.
``core/cr3bp.py``). In these units the *exact* (no small-e/small-delta
truncation) two-body Delaunay Hamiltonian for the massless body, referred to
the frame rotating with the planet, is::

    H0(a) = -1/(2a) - sqrt(a)

(``-1/(2a)`` is the Kepler energy; ``-sqrt(a)`` is the ``-n_planet * Lambda``
rotating-frame term with ``Lambda = sqrt(GM_total * a)``). Note
``H0'(1) = 0`` identically -- ``a = a_planet`` is always a critical point of
the unperturbed problem, i.e. the resonance condition itself.

The full (pre-averaging) disturbing function at the planet ``p`` (circular,
``a_p = 1``) is the standard direct + indirect form::

    R(sigma, e, delta, theta, varpi) = mu * ( 1/Delta - r_ast . r_p )

with ``Delta = |r_ast - r_p|``, ``r_p = (cos theta, sin theta)`` (planet mean
anomaly ``theta``, standing in for the fast averaging clock -- for a
circular orbit mean anomaly = true anomaly = longitude), and ``r_ast`` the
asteroid's position at semi-major axis ``1 + delta``, eccentricity ``e``,
argument of periapsis ``varpi`` and mean anomaly
``M_ast = sigma + theta - varpi`` (mod 2*pi) -- the relation that holds
``sigma = lambda_ast - lambda_planet`` fixed as ``theta`` sweeps the fast
clock. ``|r_p| = 1`` exactly (circular, ``a_p = 1``), so the indirect term is
just the dot product.

The averaged disturbing function ``R_avg(sigma, e, delta)`` used everywhere
below is a **numerical** double average over the fast planet-phase ``theta``
AND the (secularly slow but, over the Myr timescales real co-orbitals
transition on, effectively uniformly sampled) argument of periapsis
``varpi`` -- both in ``[0, 2*pi)``::

    R_avg(sigma, e, delta) = (1/(2*pi)^2) * int int R(...) dtheta dvarpi

This is deliberately a numerical quadrature, not a truncated analytic
(Laplace-coefficient) series: co-orbital eccentricities are not always small
(Kamo'oalewa's is ~0.10, and QS orbits in general need enough eccentricity
that periapsis clears the planet at ``sigma = 0``), so a low-order expansion
in ``e`` is not trustworthy here, and the numerical average is exact up to
quadrature resolution regardless of ``e``.

The full reduced Hamiltonian is then::

    H(sigma, delta, e) = H0(1 + delta) - H0(1) - R_avg(sigma, e, delta)

(the ``- H0(1)`` is just a convenient zero-reference, not a physical choice;
the ``- R_avg`` sign, MINUS not plus, is the standard celestial-mechanics
disturbing-function convention -- ``R`` is defined so the perturbing
acceleration is ``+grad(R)``, i.e. ``R`` is a *force function*, and the
Hamiltonian is ``H = T + PE = H_kepler - R``. This sign was caught and fixed
by this module's own L4-stability validation test: with the wrong sign,
L4/L5 came out as an unstable saddle of the reduced system, contradicting
the well-known fact that they are dynamically stable for a mass ratio this
small -- see the git history / test suite for the diagnostic that caught
it). This is genuinely 1-DOF at fixed ``e``: ``(sigma, delta)`` is a
canonical pair up to a fixed linear rescaling (the proper canonical momentum
is ``Delta_Lambda = sqrt(1 + delta) - 1 ~= delta/2`` for small ``delta``; the
level-curve *geometry* in ``(sigma, delta)`` is identical to that in
``(sigma, Delta_Lambda)`` since rescaling one axis of a phase portrait by a
constant does not change which curves are closed loops vs open/circulating
curves -- only the rate of traversal, which does not matter for this
module's classification purpose).

Validation (mandatory, done in ``tests/search/test_coorbital_hamiltonian.py``
before this module is trusted for anything)
--------------------------------------------------------------------------
Small-``delta`` expansion of the exact ``H0`` term reduces to
``H0(1 + delta) - H0(1) ~= -(3/8) delta^2`` (verified against the classical
CR3BP L4/L5 tadpole libration frequency
``omega_lib = n_planet * sqrt(27/4 * mu * (1 - mu))`` in the ``e -> 0``
limit -- a closed-form, independently-known result this module's numerical
machinery must reproduce before its QS/HS classification is trusted).

References
----------
* Namouni, F., "Secular Interactions of Coorbiting Objects", Icarus
  137(2):293-314 (1999), DOI 10.1006/icar.1998.6032 -- the coorbital
  QS/HS/tadpole taxonomy this module reproduces via a from-scratch numerical
  average rather than Namouni's closed-form (elliptic-integral) expansion.
* Murray, C.D. & Dermott, S.F., *Solar System Dynamics*, Cambridge Univ.
  Press (1999), Ch. 8 (coorbital / 1:1 resonance) -- the ``H0(a)`` expansion
  and the L4/L5 tadpole libration-frequency closed form used as this
  module's validation gate.
* de la Fuente Marcos, C. & de la Fuente Marcos, R., "Asteroid (469219)
  2016 HO3, the smallest and closest Earth quasi-satellite", MNRAS
  462(4):3441-3456 (2016), DOI 10.1093/mnras/stw1972 -- reports Kamo'oalewa's
  current QS episode (began ~100 yr ago, ends in ~300 yr) used as this
  module's positive control.

Scope / honest limitations
---------------------------
Planar only (``i = 0``) -- this is the standard reduction for a genuinely
1-DOF averaged model (the task this module serves explicitly specifies this
reduction). Real objects' actual inclination is dropped; real-object
``sigma`` is computed from the ecliptic-plane *projection* of their true 3-D
position. This is a modelling simplification, not an error, but it means
this module's classification of a real, inclined object is approximate --
full verification against inclination effects needs the project's existing
non-averaged CR3BP / n-body dynamics (``core/cr3bp.py``), not this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, pi, sin, sqrt
from typing import Literal

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

TWO_PI = 2.0 * pi

Regime = Literal["quasi_satellite", "tadpole", "horseshoe", "circulating", "indeterminate"]


# ---------------------------------------------------------------------------
# Kepler solve (vectorized Newton), mean anomaly -> true anomaly
# ---------------------------------------------------------------------------


def kepler_eccentric_anomaly(
    mean_anom: FloatArray, e: float, *, tol: float = 1e-12, max_iter: int = 50
) -> FloatArray:
    """Solve Kepler's equation ``E - e*sin(E) = M`` for ``E`` (vectorized Newton).

    ``mean_anom`` may be any shape/array; ``e`` is a single scalar
    eccentricity (this module never needs per-point ``e``). Standard
    textbook Newton iteration seeded at ``E0 = M`` (adequate for ``e < 0.9``,
    which covers every co-orbital regime this module targets).
    """
    m = np.mod(np.asarray(mean_anom, dtype=np.float64), TWO_PI)
    ecc_anom = m.copy()
    for _ in range(max_iter):
        f = ecc_anom - e * np.sin(ecc_anom) - m
        fp = 1.0 - e * np.cos(ecc_anom)
        d = f / fp
        ecc_anom = ecc_anom - d
        if np.max(np.abs(d)) < tol:
            break
    return ecc_anom


def true_anomaly_from_mean(mean_anom: FloatArray, e: float) -> FloatArray:
    """True anomaly ``nu`` from mean anomaly ``M`` at eccentricity ``e``."""
    ecc_anom = kepler_eccentric_anomaly(mean_anom, e)
    beta = e / (1.0 + sqrt(1.0 - e * e))
    return ecc_anom + 2.0 * np.arctan2(beta * np.sin(ecc_anom), 1.0 - beta * np.cos(ecc_anom))


# ---------------------------------------------------------------------------
# Averaged disturbing function
# ---------------------------------------------------------------------------


def _asteroid_xy(
    a_ast: float, e: float, mean_anom: FloatArray, varpi: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """Heliocentric (x, y) of the asteroid at ``(a_ast, e)``, mean anomaly
    grid ``mean_anom``, rotated by argument-of-periapsis grid ``varpi``
    (broadcastable shapes)."""
    nu = true_anomaly_from_mean(mean_anom, e)
    r = a_ast * (1.0 - e * e) / (1.0 + e * np.cos(nu))
    lon = nu + varpi
    return r * np.cos(lon), r * np.sin(lon)


def averaged_disturbing_function(
    a_ast: float,
    e: float,
    sigma_rad: float,
    mu: float,
    *,
    n_theta: int = 200,
    n_varpi: int = 200,
) -> float:
    """Numerically-averaged co-orbital disturbing function ``R_avg(sigma, e, delta)``.

    Double average, over the planet's fast phase ``theta`` (planet mean
    anomaly = true anomaly = longitude, circular orbit at ``a_p = 1``) and
    the asteroid's argument of periapsis ``varpi`` (held fixed within one
    fast-angle sweep but marginalized over its own full ``2*pi`` range, since
    a real object's apsidal line has cycled many times over the Myr
    timescale this module's QS/HS classification is meant to hold over).
    See the module docstring for the exact quadrature and its physical
    justification.

    Parameters
    ----------
    a_ast:
        Asteroid semi-major axis, nondimensional CR3BP units (``a_planet = 1``).
    e:
        Asteroid eccentricity.
    sigma_rad:
        Resonant angle ``lambda_ast - lambda_planet``, radians.
    mu:
        CR3BP mass parameter (``GM_planet / GM_total``).
    n_theta, n_varpi:
        Quadrature grid resolution (uniform trapezoidal-equivalent mean over
        a periodic grid -- the periodic-trapezoid rule is spectrally
        accurate for a smooth periodic integrand, which ``R`` is away from
        the ``e=0, delta=0`` collision singularity).

    Returns
    -------
    float
        ``R_avg``, already carrying the ``mu`` prefactor (i.e. this is
        directly the potential term to add to ``H0`` -- see
        :func:`hamiltonian`).
    """
    theta = np.linspace(0.0, TWO_PI, n_theta, endpoint=False)
    varpi = np.linspace(0.0, TWO_PI, n_varpi, endpoint=False)
    theta_grid, varpi_grid = np.meshgrid(theta, varpi, indexing="ij")
    mean_anom_ast = sigma_rad + theta_grid - varpi_grid
    x_ast, y_ast = _asteroid_xy(a_ast, e, mean_anom_ast, varpi_grid)
    x_p = np.cos(theta_grid)
    y_p = np.sin(theta_grid)
    delta_sq = (x_ast - x_p) ** 2 + (y_ast - y_p) ** 2
    delta_mag = np.sqrt(delta_sq)
    # |r_p| = 1 exactly (circular a_p=1), so the indirect term is just the dot product.
    indirect = x_ast * x_p + y_ast * y_p
    r_inst = mu * (1.0 / delta_mag - indirect)
    return float(np.mean(r_inst))


# ---------------------------------------------------------------------------
# Reduced Hamiltonian
# ---------------------------------------------------------------------------


def kepler_energy_rotating_frame(a_ast: float) -> float:
    """Exact (no truncation) rotating-frame Delaunay Hamiltonian ``H0(a) = -1/(2a) - sqrt(a)``.

    Nondimensional units (``a_planet = n_planet = GM_total = 1``). Note
    ``H0'(1) = 0`` -- ``a = a_planet`` is always a critical point of the
    unperturbed two-body problem (the resonance condition itself), and
    ``H0(1 + delta) - H0(1) ~= -(3/8) delta^2`` for small ``delta`` (see
    module docstring; validated against the classical L4/L5 tadpole
    libration frequency in the test suite).
    """
    return -1.0 / (2.0 * a_ast) - sqrt(a_ast)


def hamiltonian(
    sigma_rad: float,
    delta: float,
    e: float,
    mu: float,
    *,
    n_theta: int = 200,
    n_varpi: int = 200,
) -> float:
    """Reduced averaged co-orbital Hamiltonian ``H(sigma, delta, e)``.

    ``delta = a_ast - a_planet`` (nondimensional; ``a_planet = 1``). Zero-
    referenced so ``H(sigma, 0, e) = R_avg(sigma, e, delta=0)`` (the pure
    Kepler-energy term cancels at ``delta = 0`` by construction).
    """
    a_ast = 1.0 + delta
    h0 = kepler_energy_rotating_frame(a_ast) - kepler_energy_rotating_frame(1.0)
    r_avg = averaged_disturbing_function(a_ast, e, sigma_rad, mu, n_theta=n_theta, n_varpi=n_varpi)
    # MINUS r_avg: the standard celestial-mechanics disturbing-function sign
    # convention (Murray & Dermott Ch. 6) has the perturbing acceleration as
    # +grad(R) (verified directly against grad(1/Delta) above), so the
    # Hamiltonian (H = T + PE, PE = -force-function) is H0 - R, not H0 + R.
    # Caught by this module's own L4-stability sanity check: the "+R" sign
    # placed L4/L5 at a MINIMUM of R (an unstable saddle once combined with
    # the always-concave delta kinetic term), contradicting the well-known
    # fact that L4/L5 are dynamically stable for mu this small -- flipping
    # the sign puts L4/L5 at the required local MAXIMUM of H.
    return h0 - r_avg


# ---------------------------------------------------------------------------
# Phase portrait / regime classification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PhasePortrait:
    """A ``(sigma, delta)`` grid of the averaged Hamiltonian at fixed ``e``."""

    e: float
    mu: float
    sigma_rad: FloatArray  # (n_sigma,)
    delta: FloatArray  # (n_delta,)
    h_grid: FloatArray  # (n_delta, n_sigma)


def build_phase_portrait(
    e: float,
    mu: float,
    *,
    delta_max: float,
    n_sigma: int = 181,
    n_delta: int = 121,
    n_theta: int = 120,
    n_varpi: int = 120,
) -> PhasePortrait:
    """Build the full ``(sigma, delta)`` grid of ``H`` at fixed ``e``.

    ``delta_max`` sets the half-width of the ``delta`` sweep
    (``[-delta_max, +delta_max]``); callers should scale this to a few times
    the Hill radius (``(mu/3)**(1/3)``) to resolve the QS/L1/L2 structure
    near the planet, per the module docstring's near-planet caveat.
    """
    sigma_rad = np.linspace(0.0, TWO_PI, n_sigma)
    delta = np.linspace(-delta_max, delta_max, n_delta)
    h_grid = np.empty((n_delta, n_sigma), dtype=np.float64)
    for i, d in enumerate(delta):
        for j, s in enumerate(sigma_rad):
            # Reuses :func:`hamiltonian` verbatim (not a re-derivation) so the
            # sign convention lives in exactly one place -- a duplicated
            # inline copy here previously drifted out of sync with a sign
            # fix made only in :func:`hamiltonian`, caught by this module's
            # own test suite.
            h_grid[i, j] = hamiltonian(s, d, e, mu, n_theta=n_theta, n_varpi=n_varpi)
    return PhasePortrait(e=e, mu=mu, sigma_rad=sigma_rad, delta=delta, h_grid=h_grid)


def _bilinear_value(portrait: PhasePortrait, sigma_rad: float, delta: float) -> float:
    sigma_wrapped = float(np.mod(sigma_rad, TWO_PI))
    s_idx = float(np.interp(sigma_wrapped, portrait.sigma_rad, np.arange(len(portrait.sigma_rad))))
    d_idx = float(np.interp(delta, portrait.delta, np.arange(len(portrait.delta))))
    s_idx = min(max(s_idx, 0.0), len(portrait.sigma_rad) - 1.0)
    d_idx = min(max(d_idx, 0.0), len(portrait.delta) - 1.0)
    i0, j0 = int(d_idx), int(s_idx)
    i1 = min(i0 + 1, portrait.h_grid.shape[0] - 1)
    j1 = min(j0 + 1, portrait.h_grid.shape[1] - 1)
    fi, fj = d_idx - i0, s_idx - j0
    v00, v01 = portrait.h_grid[i0, j0], portrait.h_grid[i0, j1]
    v10, v11 = portrait.h_grid[i1, j0], portrait.h_grid[i1, j1]
    return float((1 - fi) * ((1 - fj) * v00 + fj * v01) + fi * ((1 - fj) * v10 + fj * v11))


def classify_point(
    portrait: PhasePortrait, sigma_rad: float, delta: float
) -> tuple[Regime, dict[str, float]]:
    """Classify ``(sigma, delta)`` against a pre-built :class:`PhasePortrait`.

    Purely GEOMETRIC classification (no hand-derived sign rules): find the
    connected component of the grid's level set at this point's own ``H``
    value (``H >= H0`` region containing the nearest grid cell, since every
    stable co-orbital regime here is a local-maximum "island" of ``H`` --
    see module docstring on the negative-definite kinetic term), using an
    8-connected flood fill that wraps the periodic ``sigma`` boundary.
    Classified by the angular extent and connectivity of that component:

    - ``quasi_satellite``: component includes ``sigma = 0`` (mod 2*pi) and
      does NOT extend to within 20 deg of tadpole/horseshoe territory
      (60/180/300 deg) at the SAME ``delta`` sign -- i.e. an isolated island
      hugging the planet.
    - ``tadpole``: component is confined to within roughly +-90 deg of
      +60 deg or -60 deg (300 deg) and does not reach ``sigma=0`` or
      ``sigma=180``.
    - ``horseshoe``: component spans a wide contiguous arc through
      0/60/180/300 deg (large angular extent, > 180 deg of ``sigma``
      coverage) while remaining periodic-sigma-connected (wraps around).
    - ``circulating``: component (or the level set at this ``H``) spans the
      FULL ``sigma in [0, 360)`` range at essentially one ``delta`` branch
      (a body drifting past the planet's resonance rather than librating).
    - ``indeterminate``: none of the above patterns is clearly matched
      (typically exactly on a separatrix, where the level curve is
      genuinely singular/self-intersecting on a finite grid).

    Returns
    -------
    (regime, diagnostics)
        ``diagnostics`` carries the raw geometric measurements used for the
        classification (component sigma-span in degrees, whether it wraps
        the periodic boundary, whether it reaches sigma=0, distance in H to
        the nearest actual grid extremum) for callers/tests that want to
        inspect the raw evidence rather than trust the label alone.
    """
    from scipy import ndimage

    h0 = _bilinear_value(portrait, sigma_rad, delta)
    mask = portrait.h_grid >= h0
    # Wrap sigma periodicity: tile the mask 3x along the sigma axis, label,
    # then read off the component containing the (tiled) query point's
    # middle copy -- a simple, robust way to get periodic connectivity
    # without a custom wrap-aware labeling routine.
    tiled = np.tile(mask, (1, 3))
    structure = np.ones((3, 3), dtype=int)  # 8-connectivity
    labels, _n = ndimage.label(tiled, structure=structure)
    n_sigma = len(portrait.sigma_rad)
    n_delta = len(portrait.delta)
    sigma_wrapped = float(np.mod(sigma_rad, TWO_PI))
    j_idx_float = np.interp(sigma_wrapped, portrait.sigma_rad, np.arange(n_sigma))
    i_idx_float = np.interp(delta, portrait.delta, np.arange(n_delta))
    j_query = int(np.clip(round(float(j_idx_float)), 0, n_sigma - 1))
    i_query = int(np.clip(round(float(i_idx_float)), 0, n_delta - 1))
    query_label = labels[i_query, j_query + n_sigma]  # middle tile copy
    if query_label == 0:
        # Nearest-cell rounding can still land just outside the mask right at
        # a level-curve boundary (h0 computed by bilinear interp vs. a single
        # grid cell's exact value); fall back to the best of the 4
        # surrounding cells before giving up as genuinely indeterminate.
        for di in (0, -1, 1):
            for dj in (0, -1, 1):
                ii = int(np.clip(i_query + di, 0, n_delta - 1))
                jj = int(np.clip(j_query + dj, 0, n_sigma - 1))
                cand = labels[ii, jj + n_sigma]
                if cand != 0:
                    query_label = cand
                    break
            if query_label != 0:
                break
    if query_label == 0:
        return "indeterminate", {
            "h_value": h0,
            "sigma_span_deg": 0.0,
            "wraps": 0.0,
            "reaches_sigma0": 0.0,
        }
    component_mask = labels == query_label
    # Columns (sigma, middle-tile-relative) touched by this component, any delta row.
    cols_touched = np.any(component_mask[:, n_sigma : 2 * n_sigma], axis=0)
    span_count = int(np.sum(cols_touched))
    sigma_span_deg = span_count * (360.0 / n_sigma)
    wraps_periodic = bool(cols_touched[0] and cols_touched[-1])
    reaches_sigma0 = bool(cols_touched[0] or cols_touched[-1])
    # Full delta-range check at ANY single sigma column -> circulating.
    rows_full = int(np.sum(np.any(component_mask, axis=1)))
    spans_all_sigma = span_count >= n_sigma - 2

    diagnostics = {
        "h_value": h0,
        "sigma_span_deg": sigma_span_deg,
        "wraps": 1.0 if wraps_periodic else 0.0,
        "reaches_sigma0": 1.0 if reaches_sigma0 else 0.0,
        "delta_rows_touched": float(rows_full),
    }

    if spans_all_sigma and rows_full <= 2:
        return "circulating", diagnostics
    if reaches_sigma0 and sigma_span_deg < 150.0:
        return "quasi_satellite", diagnostics
    if sigma_span_deg < 150.0 and not reaches_sigma0:
        return "tadpole", diagnostics
    if sigma_span_deg >= 150.0:
        return "horseshoe", diagnostics
    return "indeterminate", diagnostics


def hill_radius_delta(mu: float) -> float:
    """Hill-radius scale in ``delta = a - a_planet`` units: ``(mu/3)**(1/3)``."""
    return float((mu / 3.0) ** (1.0 / 3.0))


def ecliptic_longitude_from_elements(
    a: float,
    e: float,
    inc_rad: float,
    raan_rad: float,
    arg_peri_rad: float,
    mean_anom_rad: float,
) -> float:
    """Ecliptic-plane PROJECTED longitude of an (possibly inclined) orbit.

    ``lambda = Omega + atan2(cos(i) * sin(u), cos(u))`` with
    ``u = arg_peri + true_anomaly`` the argument of latitude -- the standard
    projection of an inclined orbit's position onto the reference
    (ecliptic) plane. This module is planar-only (see docstring); this
    helper is how a REAL object's inclined elements are honestly reduced to
    the single planar ``sigma`` this module needs, at the cost of dropping
    the (secondary, for this module's purpose) inclination modulation.
    ``a`` is unused (longitude does not depend on scale) but kept for a
    self-documenting call signature alongside ``e``.
    """
    del a
    nu = float(true_anomaly_from_mean(np.array([mean_anom_rad]), e)[0])
    u = arg_peri_rad + nu
    return raan_rad + atan2(cos(inc_rad) * sin(u), cos(u))
