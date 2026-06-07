"""Restricted n-body harness — VALIDATION INFRASTRUCTURE ONLY.

This package is a cross-check / validation layer, not part of the production
trajectory pipeline. Nothing here is consumed by the construct / score / verify
of catalogue rows (seeds-not-tracks is intact — mirror
``verify/spice_kernels.py:5-8``). It exists to:

1. propagate held SILVER candidates and the Jones / Aldrin trajectories through
   an *independent integrator* (REBOUND / IAS15) over the *same* DE440 ephemeris
   astropy already caches, and
2. measure the n-body correction ΔV needed to restore periodicity from a
   patched-conic seed.

The whole design rests on reusing ``core/ephemeris.py``'s time/frame convention
(TDB seconds-since-J2000 axis, heliocentric J2000-ecliptic frame) and #129's
shared-DE440 kernel machinery (``verify/spice_kernels.py``) rather than
reinventing either — see ``docs/superpowers/specs/2026-06-06-nbody-harness-design.md``
§0 (the "~64.184 s" TDB↔TT/UTC trap).

REBOUND is an optional dependency (the ``validation`` extra, design Q2). Modules
that need it import lazily and tests skip cleanly when it is absent.
"""

from __future__ import annotations

__all__: list[str] = []
