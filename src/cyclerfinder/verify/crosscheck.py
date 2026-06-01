"""V1 Lambert cross-check (M7 sub-deliverable B) — spec §14 V1.

Spec reference
--------------
* §14 V1 — "every leg re-solved with **lamberthub izzo + gooding**,
  agreement < 1e-3 m/s; full trajectory re-propagated with the
  **Kepler** propagator (not the Lambert that built it), planet
  positions met < tol."

This module ships the lamberthub-agreement half of V1. The
Kepler-re-propagation half is the V2 path: M6a's
:func:`cyclerfinder.verify.propagate.propagate_lap` re-propagates the
whole trajectory with :func:`cyclerfinder.core.kepler.propagate`, not
the Lambert solver that built it, so V1+V2 together cover the spec §14
V1 wording. M7's V0/V1 batch runner (:mod:`cyclerfinder.data.discover`)
consumes :func:`crosscheck_cycler`.

The agreement is between the in-house single-rev :func:`lambert` and
``lamberthub``'s ``izzo2015`` / ``gooding1990`` solvers on the same
``(r1, r2, tof)`` endpoints — an independent-implementation check that
the in-house solver did not converge to a subtly wrong velocity. The
M1 gate already validates this on synthetic legs; M7 wires it to the
catalogue's real cycler legs.

Plan: ``docs/phases/m7-catalogue-novelty-matching/plan.md`` §3.4.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.model.cycler import Cycler, Leg

V1_TOLERANCE_MPS: Final[float] = 1.0e-3
"""Spec §14 V1 agreement bound: the worst per-leg velocity disagreement
between the in-house Lambert and either lamberthub solver, in metres
per second. Spec-fixed; NOT test-tunable."""


@dataclass(frozen=True)
class LambertCrosscheckResult:
    """Per-leg V1 cross-check record (spec §14 V1).

    Attributes
    ----------
    leg_index:
        Position of the leg within the cycler (0-based).
    mine_v1_kms:
        Departure velocity (km/s) from the in-house :func:`lambert`.
    lamberthub_izzo_v1_kms:
        Departure velocity (km/s) from ``lamberthub.izzo2015``.
    lamberthub_gooding_v1_kms:
        Departure velocity (km/s) from ``lamberthub.gooding1990``.
    max_diff_mps:
        Worst of ``||mine - izzo||`` / ``||mine - gooding||`` over the
        departure velocity, in metres per second.
    passed:
        ``max_diff_mps < V1_TOLERANCE_MPS``. Named ``passed`` rather
        than ``pass`` because ``pass`` is a Python keyword.
    """

    leg_index: int
    mine_v1_kms: tuple[float, float, float]
    lamberthub_izzo_v1_kms: tuple[float, float, float]
    lamberthub_gooding_v1_kms: tuple[float, float, float]
    max_diff_mps: float
    passed: bool


def _leg_endpoints(
    leg: Leg,
    cycler: Cycler,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Heliocentric endpoint positions for ``leg`` from the cycler's encounters.

    Matches the leg's departure / arrival epochs against the cycler's
    encounter records and returns the planet positions ``(r1, r2)``.
    Raises :class:`ValueError` if no encounter matches an endpoint
    epoch (a malformed cycler).
    """
    r1: NDArray[np.float64] | None = None
    r2: NDArray[np.float64] | None = None
    for enc in cycler.encounters:
        if r1 is None and enc.body == leg.from_body and enc.t == leg.t_depart:
            r1 = np.asarray(enc.r, dtype=np.float64)
        if r2 is None and enc.body == leg.to_body and enc.t == leg.t_arrive:
            r2 = np.asarray(enc.r, dtype=np.float64)
    if r1 is None or r2 is None:
        raise ValueError(
            f"leg {leg.from_body}->{leg.to_body} "
            f"({leg.t_depart}..{leg.t_arrive}) has no matching encounter "
            f"endpoints in the cycler."
        )
    return r1, r2


def crosscheck_leg(
    leg: Leg,
    cycler: Cycler,
    ephem: Ephemeris,
    *,
    leg_index: int = 0,
    mu: float = MU_SUN_KM3_S2,
) -> LambertCrosscheckResult:
    """Re-solve ``leg`` with the in-house and lamberthub Lambert solvers.

    Compares the in-house :func:`lambert` departure velocity against
    ``lamberthub.izzo2015`` and ``lamberthub.gooding1990`` on the same
    ``(r1, r2, tof)`` endpoints. The ``ephem`` argument is accepted for
    interface symmetry with the V2 path (the endpoints come from the
    cycler's encounter records, which already carry the real planet
    positions); it is currently unused.

    Single-rev only: the in-house solver is single-rev (M1), so this
    cross-check is meaningful for ``leg.n_revs == 0`` legs. Multi-rev
    legs raise :class:`NotImplementedError` so callers route them to
    the multi-rev blocker rather than silently comparing against a
    single-rev in-house result.
    """
    del ephem  # endpoints come from the cycler; kept for interface symmetry
    if leg.n_revs > 0:
        raise NotImplementedError(
            f"crosscheck_leg is single-rev only; leg {leg_index} has "
            f"n_revs={leg.n_revs}. Multi-rev Lambert is out of M1/M7 scope."
        )

    from lamberthub import gooding1990, izzo2015  # type: ignore[import-untyped]

    r1, r2 = _leg_endpoints(leg, cycler)
    tof_sec = leg.t_arrive - leg.t_depart

    mine = lambert(r1, r2, tof_sec, mu=mu, prograde=True)[0]
    v1_izzo, _v2_izzo = izzo2015(mu, r1, r2, tof_sec, M=0, prograde=True)
    v1_gooding, _v2_gooding = gooding1990(mu, r1, r2, tof_sec, M=0, prograde=True)

    mine_v1 = np.asarray(mine.v1, dtype=np.float64)
    izzo_v1 = np.asarray(v1_izzo, dtype=np.float64)
    gooding_v1 = np.asarray(v1_gooding, dtype=np.float64)

    diff_izzo = float(np.linalg.norm(mine_v1 - izzo_v1))
    diff_gooding = float(np.linalg.norm(mine_v1 - gooding_v1))
    max_diff_mps = 1000.0 * max(diff_izzo, diff_gooding)

    return LambertCrosscheckResult(
        leg_index=leg_index,
        mine_v1_kms=(float(mine_v1[0]), float(mine_v1[1]), float(mine_v1[2])),
        lamberthub_izzo_v1_kms=(float(izzo_v1[0]), float(izzo_v1[1]), float(izzo_v1[2])),
        lamberthub_gooding_v1_kms=(
            float(gooding_v1[0]),
            float(gooding_v1[1]),
            float(gooding_v1[2]),
        ),
        max_diff_mps=max_diff_mps,
        passed=max_diff_mps < V1_TOLERANCE_MPS,
    )


def crosscheck_cycler(
    cycler: Cycler,
    ephem: Ephemeris,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> tuple[LambertCrosscheckResult, ...]:
    """Run :func:`crosscheck_leg` across every single-rev leg of ``cycler``.

    Multi-rev legs are skipped (no in-house single-rev result to
    compare against). The aggregated V1 pass is
    ``all(r.passed for r in result)`` AND every single-rev leg is
    represented; a cycler with only multi-rev legs returns an empty
    tuple, which the caller must treat as "V1 not applicable" rather
    than "V1 passed".
    """
    results: list[LambertCrosscheckResult] = []
    for j, leg in enumerate(cycler.legs):
        if leg.n_revs > 0:
            continue
        results.append(crosscheck_leg(leg, cycler, ephem, leg_index=j, mu=mu))
    return tuple(results)


__all__ = [
    "V1_TOLERANCE_MPS",
    "LambertCrosscheckResult",
    "crosscheck_cycler",
    "crosscheck_leg",
]
