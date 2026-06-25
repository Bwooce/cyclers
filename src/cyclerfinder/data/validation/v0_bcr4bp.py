"""V0 internal-consistency BCR4BP gauntlet (#305 Part B).

Spec reinterpretation for BCR4BP periodic orbits
------------------------------------------------
Spec §14 V0 is "internal consistency: hard constraints met; closure residual
<= tol (idealized)". For a BCR4BP periodic orbit the hard constraints are:

  * the corrector converged (its masked residual <= the corrector floor);
  * the IC and the trajectory over one period are finite;
  * the minimum periapsis radius to Earth and to Moon over one period is at or
    above the sourced per-body flyby floor (the orbit does not pass inside a
    primary);
  * ``sun_commensurate_n`` is a positive integer (a well-formed commensurability
    label);
  * the periodicity is *conditional on Sun-commensurability* — V0 SURFACES the
    ``sun_phase_drift`` rather than hiding it. An orbit whose drift exceeds the
    labelled convention is tagged ``quasi_periodic=True`` (it does not silently
    pass as a strict periodic orbit).

Periapsis floors (sourced)
--------------------------
Earth: ``radius_eq_km + safe_alt_km`` = 6378.137 + 200 (Russell-Ocampo 2003
design floor, ``data/flyby_altitude_references.yaml`` / ``constants.py``).
Moon:  1737.4 + 100 (engineering default, ``satellites.py`` / refs table).
These are read live from :data:`cyclerfinder.core.constants.PLANETS` and
:data:`cyclerfinder.core.satellites.SATELLITES`, never hard-coded here, per the
digest-not-adoption discipline.

The phase-drift threshold is a LABELLED CONVENTION
--------------------------------------------------
``V0_BCR4BP_PHASE_DRIFT_CONVENTION = 1e-6 rad`` is a JUDGMENT-CALL convention
for "Sun-commensurate" (``|omega_sun*T - 2*pi*n| < threshold``), NOT a sourced
physical constant. It is labelled as such (§5 of the design draft). V0 does not
FAIL a non-commensurate orbit — it TAGS it ``quasi_periodic`` so V2+ can refuse
a strict-period lap on it.

Discipline
----------
* READ-ONLY on the BCR4BP genome (wrap, never re-solve).
* Periapsis floors read live from the sourced constants tables.
* NO catalogue writeback.

References
----------
* spec §14 V0; ``data/flyby_altitude_references.yaml`` (#428 flyby floors).
* Rosales-Jorba 2023 (BCR4BP commensurability).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.validation.v1_bcr4bp import SEM_L_KM
from cyclerfinder.genome.bcr4bp_genome import BCR4BPPeriodicOrbit

V0_BCR4BP_CORRECTOR_FLOOR: Final[float] = 1e-10
"""V0 corrector-residual floor (nondim). Matches the corrector's ``tol`` default."""

V0_BCR4BP_PHASE_DRIFT_CONVENTION: Final[float] = 1e-6
"""LABELLED CONVENTION (not sourced): ``sun_phase_drift`` below this (rad) is
treated as Sun-commensurate. Above it the orbit is tagged ``quasi_periodic``.
A judgment call per the design draft §5, NOT a physical constant."""

V0_BCR4BP_PERIAPSIS_SAMPLES: Final[int] = 2000
"""Number of dense samples over one period for the min-periapsis scan."""


def _earth_flyby_floor_km() -> float:
    """Earth periapsis floor (km): equatorial radius + sourced safe-altitude."""
    e = PLANETS["E"]
    return float(e.radius_eq_km) + float(e.safe_alt_km)


def _moon_flyby_floor_km() -> float:
    """Moon periapsis floor (km): mean radius + sourced safe-altitude."""
    m = SATELLITES["Moon"]
    return float(m.radius_eq_km) + float(m.safe_alt_km)


@dataclass(frozen=True)
class V0VerdictBCR4BP:
    """Frozen V0 verdict for a BCR4BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id :
        Identifier carried for the audit trail.
    converged_corrector :
        Whether ``corrector_residual <= corrector_floor``.
    state_finite :
        Whether the IC and the sampled trajectory over one period are finite.
    min_periapsis_earth_km / min_periapsis_moon_km :
        Minimum distance to Earth / Moon over one period (km). ``inf`` if the
        period propagation failed.
    earth_floor_km / moon_floor_km :
        The sourced flyby floors this verdict was held against.
    periapsis_ok :
        Whether both min-periapsis values are at or above their floors.
    sun_commensurate_n :
        Echoed commensurability number.
    sun_commensurate_n_is_positive_int :
        Whether ``sun_commensurate_n`` is a well-formed positive integer.
    sun_phase_drift :
        Echoed ``|omega_sun*T - 2*pi*n|`` (rad).
    phase_drift_convention :
        The labelled convention threshold this verdict used.
    quasi_periodic :
        ``True`` iff ``sun_phase_drift > phase_drift_convention`` — the orbit is
        NOT Sun-commensurate (quasi-periodic). Surfaced, not hidden.
    passes_v0_bcr4bp :
        ``converged_corrector AND state_finite AND periapsis_ok AND
        sun_commensurate_n_is_positive_int``. Headline boolean.

        NOTE: ``quasi_periodic`` does NOT veto V0 (a quasi-periodic orbit can
        still be internally consistent over one nominal period) — but a passing
        V0 with ``quasi_periodic=True`` is a flag that V2+ must refuse a strict
        multi-lap period (see :mod:`v2_bcr4bp`).
    notes :
        Free-form audit string.
    """

    candidate_id: str
    converged_corrector: bool
    state_finite: bool
    min_periapsis_earth_km: float
    min_periapsis_moon_km: float
    earth_floor_km: float
    moon_floor_km: float
    periapsis_ok: bool
    sun_commensurate_n: int
    sun_commensurate_n_is_positive_int: bool
    sun_phase_drift: float
    phase_drift_convention: float
    quasi_periodic: bool
    passes_v0_bcr4bp: bool
    notes: str = ""


def _min_periapsis_over_period(
    orbit: BCR4BPPeriodicOrbit,
    *,
    n_samples: int,
    rtol: float,
    atol: float,
) -> tuple[float, float, bool]:
    """Scan min distance to Earth and Moon over one full period (km).

    Returns ``(min_earth_km, min_moon_km, all_finite)``. The Earth sits at
    ``(-mu, 0, 0)`` and the Moon at ``(1-mu, 0, 0)`` in the synodic frame; the
    scan propagates the FULL period (not the half) so a symmetric-residual
    orbit still gets its whole geometry checked.
    """
    mu = float(orbit.system.mu)
    earth = np.array([-mu, 0.0, 0.0])
    moon = np.array([1.0 - mu, 0.0, 0.0])
    period = float(orbit.period_nondim)
    times = np.linspace(0.0, period, max(int(n_samples), 2))

    min_e = float("inf")
    min_m = float("inf")
    cur = np.asarray(orbit.state_initial, dtype=np.float64).copy()
    if not np.all(np.isfinite(cur)):
        return float("inf"), float("inf"), False
    t_prev = 0.0
    for t in times:
        if t > 0.0:
            try:
                arc = bcr4bp.propagate_bcr4bp(
                    orbit.system, cur, t - t_prev, with_stm=False, t0=t_prev, rtol=rtol, atol=atol
                )
            except RuntimeError:
                return min_e, min_m, False
            cur = arc.state_f
            t_prev = t
        if not np.all(np.isfinite(cur)):
            return min_e, min_m, False
        pos = cur[:3]
        de = float(np.linalg.norm(pos - earth))
        dm = float(np.linalg.norm(pos - moon))
        min_e = min(min_e, de)
        min_m = min(min_m, dm)
    return min_e * SEM_L_KM, min_m * SEM_L_KM, True


def run_v0_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    corrector_floor: float = V0_BCR4BP_CORRECTOR_FLOOR,
    phase_drift_convention: float = V0_BCR4BP_PHASE_DRIFT_CONVENTION,
    earth_floor_km: float | None = None,
    moon_floor_km: float | None = None,
    n_periapsis_samples: int = V0_BCR4BP_PERIAPSIS_SAMPLES,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    notes: str = "",
) -> V0VerdictBCR4BP:
    """Run the V0 internal-consistency BCR4BP gauntlet on a periodic orbit.

    Pipeline:
      1. Re-verify the corrector residual against ``corrector_floor``.
      2. Check the IC + sampled trajectory over one period are finite.
      3. Scan min periapsis to Earth and Moon over one period; gate against the
         sourced flyby floors.
      4. Check ``sun_commensurate_n`` is a positive integer.
      5. Surface the ``sun_phase_drift``; tag ``quasi_periodic`` if it exceeds
         the labelled convention.

    Parameters
    ----------
    candidate_id :
        Identifier carried into the verdict.
    orbit :
        A :class:`BCR4BPPeriodicOrbit` from the genome. V0 re-verifies it.
    corrector_floor :
        Bar for the corrector residual (default :data:`V0_BCR4BP_CORRECTOR_FLOOR`).
    phase_drift_convention :
        Labelled commensurability threshold (default
        :data:`V0_BCR4BP_PHASE_DRIFT_CONVENTION`).
    earth_floor_km / moon_floor_km :
        Override the sourced flyby floors (default: read live from the
        constants tables).
    n_periapsis_samples :
        Dense samples over one period for the periapsis scan.
    rtol, atol :
        Integrator tolerances for the periapsis scan.
    notes :
        Free-form audit note.

    Returns
    -------
    V0VerdictBCR4BP
        ``passes_v0_bcr4bp`` is the headline boolean.

    Notes
    -----
    A V0 PASS does NOT admit to the catalogue. ``quasi_periodic=True`` on a
    passing V0 is a propagated flag, not a veto — V2+ uses it to refuse a
    strict multi-lap period.
    """
    if not isinstance(orbit, BCR4BPPeriodicOrbit):
        raise TypeError(f"orbit must be a BCR4BPPeriodicOrbit; got {type(orbit).__name__}")
    if corrector_floor <= 0.0:
        raise ValueError(f"corrector_floor must be > 0; got {corrector_floor}")
    if phase_drift_convention <= 0.0:
        raise ValueError(f"phase_drift_convention must be > 0; got {phase_drift_convention}")

    e_floor = _earth_flyby_floor_km() if earth_floor_km is None else float(earth_floor_km)
    m_floor = _moon_flyby_floor_km() if moon_floor_km is None else float(moon_floor_km)

    converged_corrector = bool(orbit.corrector_residual <= corrector_floor)

    ic_finite = bool(np.all(np.isfinite(np.asarray(orbit.state_initial, dtype=np.float64))))
    min_e_km, min_m_km, traj_finite = _min_periapsis_over_period(
        orbit, n_samples=n_periapsis_samples, rtol=rtol, atol=atol
    )
    state_finite = bool(ic_finite and traj_finite)

    periapsis_ok = bool(state_finite and min_e_km >= e_floor and min_m_km >= m_floor)

    n = orbit.sun_commensurate_n
    n_is_pos_int = bool(isinstance(n, int) and n >= 1)

    drift = float(orbit.sun_phase_drift)
    quasi_periodic = bool(drift > phase_drift_convention)

    passes = bool(converged_corrector and state_finite and periapsis_ok and n_is_pos_int)

    return V0VerdictBCR4BP(
        candidate_id=candidate_id,
        converged_corrector=converged_corrector,
        state_finite=state_finite,
        min_periapsis_earth_km=float(min_e_km),
        min_periapsis_moon_km=float(min_m_km),
        earth_floor_km=float(e_floor),
        moon_floor_km=float(m_floor),
        periapsis_ok=periapsis_ok,
        sun_commensurate_n=int(n),
        sun_commensurate_n_is_positive_int=n_is_pos_int,
        sun_phase_drift=drift,
        phase_drift_convention=float(phase_drift_convention),
        quasi_periodic=quasi_periodic,
        passes_v0_bcr4bp=passes,
        notes=notes,
    )


__all__ = [
    "V0_BCR4BP_CORRECTOR_FLOOR",
    "V0_BCR4BP_PERIAPSIS_SAMPLES",
    "V0_BCR4BP_PHASE_DRIFT_CONVENTION",
    "V0VerdictBCR4BP",
    "run_v0_bcr4bp",
]
