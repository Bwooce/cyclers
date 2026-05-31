"""Verification subpackage (M6a+).

Hosts the spec §14 V0-V5 validation gauntlet:

* :mod:`cyclerfinder.verify.propagate` (M6a) — spec §14 V2 multi-lap
  periodicity gate. Provides :func:`~cyclerfinder.verify.propagate.\
verify_long_term_stability`, the spec §12 entry point that converts
  an idealized cycler into a checkable bounded-drift statement on a
  real ephemeris over ``n_laps`` laps.

Future M6b/M7 occupants (per the M6a plan §1.5 non-goals):

* ``crosscheck.py`` (V1 Lambert cross-check, M7 stretch).
* ``gmat_bridge.py`` (V4 GMAT bridge, stretch per spec §7).

Importing this subpackage is cheap — no astropy import is triggered;
callers pass an :class:`~cyclerfinder.core.ephemeris.Ephemeris`
instance into :func:`~cyclerfinder.verify.propagate.\
verify_long_term_stability` so the astropy backend is constructed at
the caller's discretion.
"""

from __future__ import annotations
