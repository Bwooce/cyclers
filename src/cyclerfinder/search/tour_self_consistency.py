"""Per-encounter self-consistency invariant for constructed multi-encounter tours.

An "encounter" only exists if the spacecraft is actually AT the body. In an idealized
(non-ephemeris) tour CONSTRUCTION — where body positions are built from phases / conic
crossings rather than queried from a global ephemeris — it is easy to satisfy the FIRST
occurrence of each body and let REPEATED occurrences silently drift off (e.g. #480 EGGIE:
the 2nd Ganymede landed 391,000 km from the conic crossing, ≫ its ~31,700 km SOI, giving
garbage V∞ and masquerading as a physics "wall"). ``cge_scaffold`` already guards its
single Callisto crossing this way; this lifts that into a shared, reusable check.

Use it on any constructed tour BEFORE trusting its V∞ / feeding it to a corrector:
every encounter — INCLUDING repeats — must place the spacecraft within (a fraction of)
the body's sphere of influence of the body.

See ``docs/notes/2026-06-29-480-eggie-ballistic-construction-verdict.md`` and the memory
``feedback_constructed_tour_per_encounter_self_consistency``.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

Vec3 = NDArray[np.float64]


class TourSelfConsistencyError(ValueError):
    """A constructed tour places the spacecraft outside a body's SOI at an encounter."""


def soi_km(body: str) -> float:
    """Laplace sphere-of-influence radius (km) of a moon about its primary.

    ``r_SOI = a · (mu_moon / (3 · mu_primary))**(1/3)`` — the standard patched-conic SOI.
    """
    sat = SATELLITES[body]
    mu_primary = PRIMARIES[sat.primary]
    return float(sat.sma_km * (sat.mu_km3_s2 / (3.0 * mu_primary)) ** (1.0 / 3.0))


def encounter_gaps_km(
    spacecraft_positions: Sequence[Vec3], body_positions: Sequence[Vec3]
) -> list[float]:
    """Per-encounter ``|spacecraft - body|`` (km). Both sequences in encounter order."""
    if len(spacecraft_positions) != len(body_positions):
        raise ValueError(
            f"length mismatch: {len(spacecraft_positions)} spacecraft vs "
            f"{len(body_positions)} body positions"
        )
    return [
        float(np.linalg.norm(np.asarray(sc, dtype=np.float64) - np.asarray(b, dtype=np.float64)))
        for sc, b in zip(spacecraft_positions, body_positions, strict=True)
    ]


def assert_encounters_self_consistent(
    spacecraft_positions: Sequence[Vec3],
    body_positions: Sequence[Vec3],
    bodies: Sequence[str],
    *,
    fraction: float = 1.0,
    context: str = "",
) -> None:
    """Raise :class:`TourSelfConsistencyError` if any encounter body is too far from the node.

    For each encounter ``i`` the spacecraft must be within ``fraction · SOI(bodies[i])`` of
    ``body_positions[i]``. This catches the repeated-encounter-drift bug class: a body whose
    phase was pinned only at its first encounter and silently drifts at later ones.

    ``bodies`` may repeat (e.g. ``("Europa","Ganymede","Ganymede","Io","Europa")``) — every
    occurrence is checked independently, which is the whole point.
    """
    gaps = encounter_gaps_km(spacecraft_positions, body_positions)
    if not (len(gaps) == len(bodies)):
        raise ValueError(f"length mismatch: {len(gaps)} gaps vs {len(bodies)} bodies")
    failures: list[str] = []
    for i, (gap, body) in enumerate(zip(gaps, bodies, strict=True)):
        limit = fraction * soi_km(body)
        if gap > limit:
            failures.append(
                f"  encounter {i} ({body}): gap {gap:.0f} km > {fraction:g}·SOI {limit:.0f} km"
            )
    if failures:
        raise TourSelfConsistencyError(
            (f"{context}: " if context else "")
            + "constructed tour is not self-consistent — spacecraft is outside the body SOI "
            + "at "
            + str(len(failures))
            + " encounter(s) (repeated-encounter drift?):\n"
            + "\n".join(failures)
        )


__all__ = [
    "TourSelfConsistencyError",
    "assert_encounters_self_consistent",
    "encounter_gaps_km",
    "soi_km",
]
