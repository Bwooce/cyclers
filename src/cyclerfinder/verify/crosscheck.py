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

Endpoint independence (task #197, the #180/63-s false-consensus class)
----------------------------------------------------------------------
The velocity-side independence above is real, but the Lambert endpoints
themselves must NOT be trusted from the artifact under test: a cycler
whose embedded encounter position/epoch is wrong upstream would be
re-solved on the same wrong ``(r1, r2)`` by every solver, and all three
would happily agree (false consensus). So by default
:func:`crosscheck_leg` re-derives ``r1``/``r2`` independently from the
passed ephemeris at the leg's own epochs, checks them against the
encounter-embedded positions to within
:data:`POSITION_CONSISTENCY_TOL_KM`, and runs the Lambert re-solve on
the *re-queried* endpoints. A mismatch fails the leg
(``passed = False``) regardless of solver agreement.

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

POSITION_CONSISTENCY_TOL_KM: Final[float] = 1.0
"""Endpoint-consistency bound (km) between a cycler's encounter-embedded
planet position and the independent ephemeris re-query at the same
epoch. Every in-tree constructor reads ``Encounter.r`` straight from
``ephem.state(body, t)``, so a clean cycler checked against the same
ephemeris matches to float round-off; 1 km corresponds to ~0.03 s of
planetary motion — far above round-off, far below any real
encounter-position or epoch bug. Module constant; NOT test-tunable."""


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
    endpoint_mismatch_km:
        Worst position difference (km) between the encounter-embedded
        endpoints and the independent ephemeris re-query at the leg's
        epochs. ``None`` when the check was skipped
        (``independent_endpoints=False``).
    passed:
        ``max_diff_mps < V1_TOLERANCE_MPS`` AND, when the endpoint
        check ran, ``endpoint_mismatch_km <=
        POSITION_CONSISTENCY_TOL_KM``. Named ``passed`` rather than
        ``pass`` because ``pass`` is a Python keyword.
    """

    leg_index: int
    mine_v1_kms: tuple[float, float, float]
    lamberthub_izzo_v1_kms: tuple[float, float, float]
    lamberthub_gooding_v1_kms: tuple[float, float, float]
    max_diff_mps: float
    passed: bool
    endpoint_mismatch_km: float | None = None


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
    independent_endpoints: bool = True,
) -> LambertCrosscheckResult:
    """Re-solve ``leg`` with the in-house and lamberthub Lambert solvers.

    Compares the in-house :func:`lambert` departure velocity against
    ``lamberthub.izzo2015`` and ``lamberthub.gooding1990`` on the same
    ``(r1, r2, tof)`` endpoints.

    Endpoint independence (default ON): with
    ``independent_endpoints=True`` the endpoints are re-derived from
    ``ephem.state(body, t)`` at the leg's own epochs, checked against
    the cycler's encounter-embedded positions to within
    :data:`POSITION_CONSISTENCY_TOL_KM` (recorded as
    ``endpoint_mismatch_km``), and the re-solve runs on the
    *re-queried* endpoints — so a poisoned upstream encounter position
    fails the leg instead of feeding the same wrong ``(r1, r2)`` to
    every solver (false consensus, the #180/63-s bug class). All
    validation paths must keep this ON. ``independent_endpoints=False``
    is the documented escape hatch for the one legitimate exception — a
    caller whose cycler epochs are not on the passed ephemeris's time
    base (no such caller exists in-tree); it restores the
    endpoints-from-artifact behaviour and leaves
    ``endpoint_mismatch_km`` as ``None``.

    Single- and multi-rev: the in-house solver matches ``leg.n_revs`` and
    ``leg.branch`` ("low" maps to ``lamberthub``'s ``low_path=True``,
    "high" to ``low_path=False``); the reference solvers are called with
    ``M=leg.n_revs`` and the corresponding ``low_path`` so the comparison
    is on the same branch.
    """
    from lamberthub import gooding1990, izzo2015  # type: ignore[import-untyped]

    r1_enc, r2_enc = _leg_endpoints(leg, cycler)
    endpoint_mismatch_km: float | None = None
    endpoints_consistent = True
    if independent_endpoints:
        r1 = np.asarray(ephem.state(leg.from_body, leg.t_depart)[0], dtype=np.float64)
        r2 = np.asarray(ephem.state(leg.to_body, leg.t_arrive)[0], dtype=np.float64)
        endpoint_mismatch_km = max(
            float(np.linalg.norm(r1 - r1_enc)),
            float(np.linalg.norm(r2 - r2_enc)),
        )
        endpoints_consistent = endpoint_mismatch_km <= POSITION_CONSISTENCY_TOL_KM
    else:
        r1, r2 = r1_enc, r2_enc
    tof_sec = leg.t_arrive - leg.t_depart

    sols = lambert(r1, r2, tof_sec, mu=mu, prograde=True, max_revs=leg.n_revs)
    mine = next(s for s in sols if s.n_revs == leg.n_revs and s.branch == leg.branch)
    low_path = leg.branch != "high"
    v1_izzo, _v2_izzo = izzo2015(
        mu, r1, r2, tof_sec, M=leg.n_revs, prograde=True, low_path=low_path
    )
    v1_gooding, _v2_gooding = gooding1990(
        mu, r1, r2, tof_sec, M=leg.n_revs, prograde=True, low_path=low_path
    )

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
        passed=(max_diff_mps < V1_TOLERANCE_MPS) and endpoints_consistent,
        endpoint_mismatch_km=endpoint_mismatch_km,
    )


def crosscheck_cycler(
    cycler: Cycler,
    ephem: Ephemeris,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> tuple[LambertCrosscheckResult, ...]:
    """Run :func:`crosscheck_leg` across every leg of ``cycler``.

    Single- and multi-rev legs are both crosschecked (each on its own
    ``n_revs``/``branch``). The aggregated V1 pass is
    ``all(r.passed for r in result)`` with every leg represented.
    """
    results: list[LambertCrosscheckResult] = []
    for j, leg in enumerate(cycler.legs):
        results.append(crosscheck_leg(leg, cycler, ephem, leg_index=j, mu=mu))
    return tuple(results)


__all__ = [
    "V1_TOLERANCE_MPS",
    "LambertCrosscheckResult",
    "crosscheck_cycler",
    "crosscheck_leg",
]
