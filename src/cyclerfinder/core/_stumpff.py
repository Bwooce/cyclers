"""Stumpff functions C(z) and S(z) for universal-variable orbit mechanics.

Used by both :mod:`cyclerfinder.core.kepler` and :mod:`cyclerfinder.core.lambert`.
Kept in a private module so any future numerical improvement (better series
cutoff, vectorisation, etc.) lands in one place.

References
----------
Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th ed.,
Microcosm Press, 2013, §2.2 (universal variables and Stumpff functions),
eqs. 2-93 through 2-98.

Definitions
-----------
For the universal variable formulation::

    C(z) = (1 - cos(sqrt(z))) / z                  for z > 0
         = (cosh(sqrt(-z)) - 1) / (-z)             for z < 0
         = 1/2                                     for z = 0

    S(z) = (sqrt(z) - sin(sqrt(z))) / sqrt(z)^3    for z > 0
         = (sinh(sqrt(-z)) - sqrt(-z)) / sqrt(-z)^3 for z < 0
         = 1/6                                     for z = 0

Both functions are analytic everywhere; the closed forms above lose catastrophic
precision near ``z = 0`` so a Maclaurin series is used inside ``|z| < 1e-3``.
"""

from __future__ import annotations

from math import cos, cosh, sin, sinh, sqrt

# Series cutoff per the module docstring. Below this magnitude the closed-form
# expressions lose precision to subtraction of near-equal terms.
_SERIES_CUTOFF: float = 1.0e-3


def stumpff_c(z: float) -> float:
    """Return the Stumpff C(z) function.

    Parameters
    ----------
    z:
        Universal-variable argument, dimensionless (= chi**2 * alpha).

    Returns
    -------
    float
        C(z), always positive and finite for any real z.
    """
    if abs(z) < _SERIES_CUTOFF:
        # Maclaurin series: C(z) = 1/2 - z/24 + z**2/720 - z**3/40320 + ...
        # Truncated at z**3 the error is below 1e-18 for |z| < 1e-3.
        return 0.5 - z / 24.0 + z * z / 720.0 - (z * z * z) / 40320.0
    if z > 0.0:
        sqrt_z = sqrt(z)
        return (1.0 - cos(sqrt_z)) / z
    # z < 0: hyperbolic regime
    sqrt_mz = sqrt(-z)
    return (cosh(sqrt_mz) - 1.0) / (-z)


def stumpff_s(z: float) -> float:
    """Return the Stumpff S(z) function.

    Parameters
    ----------
    z:
        Universal-variable argument, dimensionless (= chi**2 * alpha).

    Returns
    -------
    float
        S(z), always positive and finite for any real z.
    """
    if abs(z) < _SERIES_CUTOFF:
        # Maclaurin series: S(z) = 1/6 - z/120 + z**2/5040 - z**3/362880 + ...
        return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - (z * z * z) / 362880.0
    if z > 0.0:
        sqrt_z = sqrt(z)
        return (sqrt_z - sin(sqrt_z)) / (sqrt_z * sqrt_z * sqrt_z)
    # z < 0: hyperbolic regime
    sqrt_mz = sqrt(-z)
    return (sinh(sqrt_mz) - sqrt_mz) / (sqrt_mz * sqrt_mz * sqrt_mz)
