# src/cyclerfinder/search/leveraging_leg.py
"""Phase-full VILM single-leg evaluator (spec 2026-06-09, Component 1).

The phase-FULL counterpart of the phase-free :mod:`cyclerfinder.search.vilm`
(Campagnola & Russell, "The Endgame Problem"). A VILM leg departs a moon M with
V∞_in on an orbit resonant with M, applies a deep-space impulse at the apse, and
returns to M with a changed V∞_out — the apse impulse IS the leveraging maneuver.

Canonical units about the primary (see plan "Canonical units"). Coplanar, circular
moon. Pure: math/scipy + core.satellites + search.vilm only.
"""

from __future__ import annotations

import math

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES


def v_m_kms(moon: str) -> float:
    """Moon circular velocity about its primary, km/s (the canonical V_M)."""
    sat = SATELLITES[moon]
    return math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)


def resonant_sma(n_moon_revs: int, m_sc_revs: int) -> float:
    """Canonical SC semimajor axis for an n:m resonance (a = (n/m)**(2/3), a_M=1)."""
    return float((n_moon_revs / m_sc_revs) ** (2.0 / 3.0))


def tisserand_vinf(*, a: float, e: float) -> float:
    """Adimensional V∞ at the moon for an SC orbit (a, e). nan if no real V∞."""
    val = 3.0 - 1.0 / a - 2.0 * math.sqrt(a * (1.0 - e * e))
    if abs(val) < 1e-15:
        return 0.0
    return math.sqrt(val) if val > 0.0 else float("nan")


def eccentricity_from_vinf(*, a: float, vinf: float) -> float:
    """Eccentricity of the orbit with semimajor axis ``a`` and V∞ ``vinf``.

    Inverts :func:`tisserand_vinf`. nan if no real moon-crossing orbit.
    """
    h = (3.0 - 1.0 / a - vinf * vinf) / 2.0
    if h < 0.0:
        return float("nan")
    ratio = (h * h) / a
    if ratio > 1.0:
        return float("nan")
    e = math.sqrt(1.0 - ratio)
    if not (a * (1.0 - e) <= 1.0 + 1e-12 and a * (1.0 + e) + 1e-12 >= 1.0):
        return float("nan")
    return e
