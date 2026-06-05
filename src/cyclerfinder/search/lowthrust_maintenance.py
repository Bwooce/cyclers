r"""Low-thrust cycler-maintenance evaluator (Phase 4 of the v2 low-thrust scope).

**Machinery only.** This module delivers the v2-goal machinery for evaluating a
cycler's per-synodic maintenance manoeuvre under a *thrust-bounded low-thrust*
model and reporting a ``propellant_mass_fraction`` — the comparability metric
flagged in ``docs/v2-future-references.md`` §1 ("``model/score.py`` adds
``propellant_mass_fraction`` so low-thrust candidates are comparable to
ballistic ones"). It does **not** add catalogue rows or schema fields: no
published source supplies a powered low-thrust cycler row we hold with
extractable numbers, so fabricating one is forbidden by the project's
golden-test discipline. See the plan's "Execution deviation (Phase 4)" note.

The evaluator is the bridge a future *sourced* powered row (or the Forge) would
consume: given a maintenance ``Delta V`` — whether from the impulsive
turn-deficit surrogate in :mod:`cyclerfinder.search.maintain` or from a
distributed Sims-Flanagan low-thrust solve
(:mod:`cyclerfinder.search.lowthrust`) — it attaches the Tsiolkovsky propellant
cost so low-thrust candidates can be ranked alongside ballistic ones.

Golden discipline: the wrapped ``Delta V`` magnitudes are *our own computed
values* (the literature publishes no powered cycler maintenance ΔV; McConaghy
2002 defers Aldrin's). Nothing here is a golden EXPECTED. The propellant
fraction is a source-free physics identity (the rocket equation), and the only
cross-check against the existing stack is an *internal-consistency* one (the
powered model re-expresses the impulsive maintenance ΔV unchanged).

Plan: ``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md`` (Phase 4).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp

from cyclerfinder.core.constants import STANDARD_GRAVITY_KM_S2


def propellant_mass_fraction(maintenance_dv_kms: float, isp_s: float) -> float:
    r"""Propellant mass fraction for a maintenance ``Delta V`` at ``Isp``.

    The Tsiolkovsky identity ``m_p / m_0 = 1 - exp(-Delta V / (g0 * Isp))``: the
    fraction of the pre-manoeuvre wet mass spent as propellant to deliver
    ``maintenance_dv_kms``. A source-free physics invariant — this is the
    comparability metric that lets a powered low-thrust maintenance cost be
    ranked against a ballistic one.

    Parameters
    ----------
    maintenance_dv_kms:
        Per-synodic maintenance ``Delta V``, km/s. Must be non-negative.
    isp_s:
        Specific impulse, seconds. Must be positive.

    Returns
    -------
    float
        Propellant mass fraction in ``[0, 1)``.
    """
    if maintenance_dv_kms < 0.0:
        raise ValueError(f"maintenance_dv_kms must be non-negative, got {maintenance_dv_kms}")
    if isp_s <= 0.0:
        raise ValueError(f"isp_s must be positive, got {isp_s}")
    return 1.0 - exp(-maintenance_dv_kms / (STANDARD_GRAVITY_KM_S2 * isp_s))


@dataclass(frozen=True)
class PoweredMaintenanceResult:
    r"""Powered (low-thrust) maintenance cost for one synodic period.

    Attributes
    ----------
    maintenance_dv_kms:
        The per-synodic maintenance ``Delta V`` modelled, km/s. Carried through
        unchanged from the supplied cost (impulsive surrogate or distributed
        low-thrust solve). **Our own computed value**, never source-attested.
    isp_s:
        Specific impulse assumed for the propellant accounting, seconds.
    propellant_mass_fraction:
        Tsiolkovsky propellant fraction of the pre-manoeuvre wet mass (see
        :func:`propellant_mass_fraction`).
    propellant_mass_kg:
        Propellant mass spent per synodic period, kg, for the supplied
        post-manoeuvre (dry) mass: ``m_0 - m_dry`` where
        ``m_0 = m_dry / (1 - propellant_mass_fraction)``.
    dry_mass_kg:
        Post-manoeuvre spacecraft mass the accounting is referenced to, kg.
    """

    maintenance_dv_kms: float
    isp_s: float
    propellant_mass_fraction: float
    propellant_mass_kg: float
    dry_mass_kg: float


def powered_maintenance_from_dv(
    maintenance_dv_kms: float,
    isp_s: float,
    dry_mass_kg: float,
) -> PoweredMaintenanceResult:
    r"""Wrap a maintenance ``Delta V`` into a powered-maintenance result.

    Attaches the Tsiolkovsky propellant accounting to a maintenance ``Delta V``
    so a low-thrust maintenance cost is comparable to a ballistic one. The
    ``Delta V`` is carried through unchanged — this routine does not re-solve the
    trajectory; it expresses an already-computed maintenance cost (from the
    impulsive surrogate in :mod:`cyclerfinder.search.maintain` or a distributed
    low-thrust solve in :mod:`cyclerfinder.search.lowthrust`) in propellant
    terms.

    Parameters
    ----------
    maintenance_dv_kms:
        Per-synodic maintenance ``Delta V``, km/s. Non-negative.
    isp_s:
        Specific impulse, seconds. Positive.
    dry_mass_kg:
        Post-manoeuvre (dry) spacecraft mass the propellant mass is referenced
        to, kg. Positive.

    Returns
    -------
    PoweredMaintenanceResult
    """
    if dry_mass_kg <= 0.0:
        raise ValueError(f"dry_mass_kg must be positive, got {dry_mass_kg}")
    frac = propellant_mass_fraction(maintenance_dv_kms, isp_s)
    # m_dry = m_0 * (1 - frac)  =>  m_0 = m_dry / (1 - frac); propellant = m_0 - m_dry.
    wet_mass = dry_mass_kg / (1.0 - frac)
    propellant_mass = wet_mass - dry_mass_kg
    return PoweredMaintenanceResult(
        maintenance_dv_kms=maintenance_dv_kms,
        isp_s=isp_s,
        propellant_mass_fraction=frac,
        propellant_mass_kg=propellant_mass,
        dry_mass_kg=dry_mass_kg,
    )


__all__ = [
    "PoweredMaintenanceResult",
    "powered_maintenance_from_dv",
    "propellant_mass_fraction",
]
