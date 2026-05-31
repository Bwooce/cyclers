"""Cycler model layer — frozen dataclasses representing a result trajectory.

Re-exports the three structures M3 defines so callers write
``from cyclerfinder.model import Cycler, Leg, Encounter`` rather than
reaching into the implementation module.
"""

from __future__ import annotations

from cyclerfinder.model.cycler import Cycler, Encounter, Leg

__all__ = ["Cycler", "Encounter", "Leg"]
