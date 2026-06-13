"""Recover the Braik-Ross 2026 common-energy representative orbits at C_J=3.1294.

The paper (Table 2) lists thirteen planar Earth-Moon periodic orbits, one per
family, all re-selected at a single common Jacobi constant ``C_J = 3.1294``, but
publishes ONLY period-in-days and the Floquet rate sigma -- no state vectors. To
score the accessibility network we must first *recover the member* at that
energy, then confirm it against the sourced period before trusting it.

Two recovery routes (mining note Q3 / proposed task 2):

* Lyapunov (LL1, LL2), distant prograde (DPO), and the resonant families are
  standard JPL 3-Body Periodic-Orbit families. We pull the member nearest
  ``C_J = 3.1294`` from the JPL oracle
  (:func:`cyclerfinder.verify.jpl_periodic_orbits.query`) as a seed, then
  re-correct it under OUR mass ratio with the fixed-Jacobi symmetric corrector so
  the recovered orbit sits exactly on ``C_J = 3.1294`` (the JPL mu differs from
  ours by ~1e-7; see ``verify.jpl_periodic_orbits`` CONVENTION RECONCILIATION).
* The four cyclers (C11a, C11b, C21, C32) are NOT in the JPL DB. We recover them
  with our own fixed-Jacobi symmetric corrector from the Ross & Roberts-Tsoukkas
  2025 (AAS 25-621) family seed regions (x0 region + half-crossing index), driven
  to the sourced ``C_J = 3.1294`` period.

SOURCED-CONFIRMATION DISCIPLINE: every recovered member is only *trusted* once its
period matches the Braik-Ross Table-2 sourced period (in days) within tolerance.
A member whose period does not match is reported as unconfirmed and must NOT enter
the scored network -- selecting members by our own criteria instead of a
published anchor would make the validation gate circular.

This module performs a network call (JPL) for the non-cycler families; the
cyclers are recovered offline. Pure-corrector helpers are network-free.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.verify.jpl_periodic_orbits import query

#: Braik-Ross common Jacobi constant (Table 2 header / Sec. 2).
C_J_BRAIK_ROSS = 3.1294

#: Ross & Roberts-Tsoukkas 2025 Earth-Moon mass ratio (AAS 25-621, p. 3). Used
#: for the cycler recovery so the corrector matches the family-defining paper.
ROSS_MU = 1.2150584270572e-2

#: 1 nondimensional time unit in days (T_EM = 27.321661 d, TU = T_EM / 2pi).
TU_DAYS = 27.321661 / (2.0 * math.pi)


def braik_ross_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP system in the Braik-Ross / Ross-RT nondimensional scales.

    Uses the AAS 25-621 mass ratio (identical to Braik-Ross Table 1 to all
    printed digits) and the standard a_M = 384400 km length / T_EM/2pi time
    scales, so periods convert to days via :data:`TU_DAYS`.
    """
    return cr3bp.CR3BPSystem(
        mu=ROSS_MU,
        primary="Earth",
        secondary="Moon",
        l_km=384400.0,
        t_s=TU_DAYS * 86400.0,
    )


@dataclass(frozen=True)
class Representative:
    """One recovered common-energy representative orbit.

    ``state0`` is the planar IC ``(x0, 0, 0, 0, ydot0, 0)``; ``period`` is
    nondimensional (multiply by :data:`TU_DAYS` for days). ``sourced_period_days``
    is the Braik-Ross Table-2 value; ``confirmed`` is True iff the recovered
    period matches it within ``tol_days``.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    sourced_period_days: float
    period_days: float
    confirmed: bool
    converged: bool


# Braik-Ross Table 2 sourced periods (days) at C_J = 3.1294 for the families we
# can source-confirm. (Resonant R31/R52 and the U-branches carry no period in the
# data available to us, so they are intentionally excluded from the gate set --
# including them would assert un-cross-checkable members; see module docstring.)
SOURCED_PERIODS_DAYS: dict[str, float] = {
    "LL1": 12.811,
    "LL2": 15.117,
    "DPO": 11.184,
    "R21-S": 26.500,
    "C11a": 42.140,
    "C11b": 55.995,
    "C21": 84.533,
    "C32": 78.613,
}

# Ross & Roberts-Tsoukkas 2025 cycler family seed regions: (x0 seed, ydot0 sign,
# half-crossing index). Source for the x0 region / half-crossing index: AAS
# 25-621 Table 3 seeds (carried in tests/search/test_cr3bp_ross_families.py) and
# this module's recovery scan at C_J=3.1294.
#
# RECOVERY STATUS at C_J=3.1294 with the 1-DOF perpendicular-x-crossing symmetric
# corrector (empirically determined; see the results note):
#   * C11b: recovers EXACTLY to the sourced 55.995 d (confirmed).
#   * C32:  recovers to ~79.50 d, ~1.1% above the sourced 78.613 d -- the
#           perpendicular-x-crossing branch intersection of the (3,2) family at
#           this OFF-STABLE common energy. Marginal; flagged, not silently
#           accepted (confirm tol below admits it only at a relaxed tolerance).
#   * C11a (42.140 d) and C21 (84.533 d): NOT recoverable with this corrector at
#           this energy (no perpendicular-crossing member within 1.5 d of the
#           sourced period across the scanned seed/half-crossing grid). They are
#           a known recovery gap (a multi-segment / 2-D shooting corrector would
#           be needed); excluded from the scored set rather than faked.
_CYCLER_SEEDS: dict[str, tuple[float, float, int]] = {
    "C11b": (-0.7682140805, -1.0, 3),
    "C32": (-0.4000000000, -1.0, 3),
}


def _jpl_seed_near_cj(
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    cj: float = C_J_BRAIK_ROSS,
    stable: bool | None = None,
) -> NDArray[np.float64]:
    """JPL member nearest ``cj`` (optionally filtered to the stable branch).

    ``stable=True`` keeps members with stability index near 1 (|stab| < ~1+eps,
    bounded/stable in JPL's reduced index convention); ``stable=False`` keeps the
    strongly unstable members; ``None`` takes the global nearest-in-C member.
    """
    _constants, orbits = query("earth-moon", family, libr=libr, branch=branch)
    js = np.array([o.jacobi for o in orbits])
    stab = np.array([abs(o.stability) for o in orbits])
    mask = np.ones(len(orbits), dtype=bool)
    if stable is True:
        mask = stab <= 1.0 + 1e-6
    elif stable is False:
        mask = stab > 1.0 + 1e-6
    if not mask.any():
        mask = np.ones(len(orbits), dtype=bool)
    cand = np.where(mask)[0]
    i = int(cand[np.argmin(np.abs(js[cand] - cj))])
    return np.asarray(orbits[i].state0, dtype=np.float64)


def recover_from_seed(
    system: cr3bp.CR3BPSystem,
    label: str,
    x0_seed: float,
    period_guess_days: float,
    *,
    ydot0_sign: float,
    half_crossings: int,
    tol_days: float = 0.5,
    corrector_tol: float = 1e-10,
) -> Representative:
    """Correct a symmetric member to ``C_J=3.1294`` and confirm its period.

    Wraps :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
    (Jacobi held at :data:`C_J_BRAIK_ROSS`) and compares the recovered period to
    the Braik-Ross sourced period. ``period_guess_days`` seeds the corrector's
    period; the sourced period is the confirmation target.
    """
    sourced = SOURCED_PERIODS_DAYS[label]
    period_guess = period_guess_days / TU_DAYS
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_seed,
        C_J_BRAIK_ROSS,
        period_guess,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        tol=corrector_tol,
    )
    period_days = orbit.period * TU_DAYS
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    confirmed = orbit.converged and abs(period_days - sourced) <= tol_days
    return Representative(
        label=label,
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        sourced_period_days=sourced,
        period_days=period_days,
        confirmed=confirmed,
        converged=orbit.converged,
    )


def recover_jpl_family(
    system: cr3bp.CR3BPSystem,
    label: str,
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    half_crossings: int = 1,
    stable: bool | None = None,
    tol_days: float = 0.5,
) -> Representative:
    """Recover a JPL-DB family member at ``C_J=3.1294`` (seed from JPL, re-correct).

    Pulls the member nearest the common energy from the JPL oracle, then
    re-corrects it under our mass ratio with the fixed-Jacobi symmetric corrector
    so it sits exactly on ``C_J = 3.1294``; confirms against the sourced period.
    """
    seed = _jpl_seed_near_cj(family, libr=libr, branch=branch, stable=stable)
    x0_seed = float(seed[0])
    ydot0_sign = math.copysign(1.0, float(seed[4]))
    return recover_from_seed(
        system,
        label,
        x0_seed,
        SOURCED_PERIODS_DAYS[label],
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        tol_days=tol_days,
    )
