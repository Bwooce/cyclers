"""Verification subpackage (M6a + M6b).

Hosts the spec §14 V0-V5 validation gauntlet:

* :mod:`cyclerfinder.verify.propagate` (M6a) — spec §14 V2 multi-lap
  periodicity gate. Provides
  :func:`~cyclerfinder.verify.propagate.verify_long_term_stability`,
  the spec §12 entry point that converts an idealized cycler into a
  checkable bounded-drift statement on a real ephemeris over
  ``n_laps`` laps.
* :mod:`cyclerfinder.verify.real_closure` (M6b) — spec §14 V2-real
  multi-cycle real-ephemeris closure gate. Provides
  :func:`~cyclerfinder.verify.real_closure.verify_real_closure`, the
  M6b entry point that promotes a catalogue entry from V1 idealised
  reproduction to V2-real ephemeris instance.

* :mod:`cyclerfinder.verify.crosscheck` (M7) — spec §14 V1 single-leg
  Lambert cross-check against ``lamberthub``. Provides
  :func:`~cyclerfinder.verify.crosscheck.crosscheck_cycler`.

Future occupants (per the M6a/M6b plans' non-goals):

* ``gmat_bridge.py`` (V4 GMAT bridge, stretch per spec §7).

Importing this subpackage is cheap — no astropy import is triggered;
callers pass an :class:`~cyclerfinder.core.ephemeris.Ephemeris`
instance into the verifier so the astropy backend is constructed at
the caller's discretion.
"""

from __future__ import annotations

from cyclerfinder.verify.crosscheck import (
    V1_TOLERANCE_MPS,
    LambertCrosscheckResult,
    crosscheck_cycler,
    crosscheck_leg,
)
from cyclerfinder.verify.fidelity import (
    FidelityRungUnavailableError,
    FidelitySolution,
    PersistenceClass,
    PersistenceReport,
    fidelity_persistence,
    solve_at_fidelity,
)
from cyclerfinder.verify.propagate import (
    DRIFT_TOLERANCE_KM,
    StabilityReport,
    lap_to_lap_drift,
    multi_lap_propagation,
    propagate_lap,
    verify_long_term_stability,
)
from cyclerfinder.verify.real_closure import (
    EXPECTED_SKIPS,
    N_CYCLES_DEFAULT,
    REAL_DRIFT_TOLERANCE_KM,
    RealClosureConstructionError,
    RealClosureResult,
    construct_real_ephemeris_cycler,
    verify_real_closure,
)

__all__ = [
    "DRIFT_TOLERANCE_KM",
    "EXPECTED_SKIPS",
    "N_CYCLES_DEFAULT",
    "REAL_DRIFT_TOLERANCE_KM",
    "V1_TOLERANCE_MPS",
    "FidelityRungUnavailableError",
    "FidelitySolution",
    "LambertCrosscheckResult",
    "PersistenceClass",
    "PersistenceReport",
    "RealClosureConstructionError",
    "RealClosureResult",
    "StabilityReport",
    "construct_real_ephemeris_cycler",
    "crosscheck_cycler",
    "crosscheck_leg",
    "fidelity_persistence",
    "lap_to_lap_drift",
    "multi_lap_propagation",
    "propagate_lap",
    "solve_at_fidelity",
    "verify_long_term_stability",
    "verify_real_closure",
]
