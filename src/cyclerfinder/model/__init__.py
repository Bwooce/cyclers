"""Cycler model layer — frozen dataclasses representing a result trajectory.

Re-exports the M3 structures (``Cycler``, ``Leg``, ``Encounter``) and the
M4 :class:`Score` so callers write ``from cyclerfinder.model import
Cycler, Score`` rather than reaching into the implementation modules.
"""

from __future__ import annotations

from cyclerfinder.model.cycler import Cycler, Encounter, Leg
from cyclerfinder.model.score import Score

__all__ = ["Cycler", "Encounter", "Leg", "Score"]
