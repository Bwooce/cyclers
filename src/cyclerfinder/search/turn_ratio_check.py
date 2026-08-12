"""#833 — Russell's published TURN RATIO as a reusable sourced cross-check.

Russell-Ocampo 2003 p.13 (digest ``docs/notes/2026-06-17-digest-russell-ocampo-2003.md``):

    "Turn Ratio TR = max physically allowable turn angle / max required turn
     angle (delta_MAX)"; TR > 1 => all flybys physically attainable; the max
     allowable is based on a **200 km altitude** Earth flyby.

Our :func:`cyclerfinder.search.correct._max_bend_deg` uses
``PLANETS[body].safe_alt_km = 200.0`` km for Earth and Mars — itself sourced from
Russell 2004 p.165 ``r_p,min`` — so

    ``min`` over a closure's intermediate flybys of ``max_bend / required_bend``

is a **like-for-like reproduction of Russell's own TR**, and a catalogue row's
published ``invariants.turn_ratio`` is a legitimate GOLDEN expected side for it
(``[[feedback_golden_tests_sourced_only]]``: the expected value is published, the
measured value emerges from the converged trajectory).

Why this is worth a named instrument
------------------------------------
`#826` used this check as throwaway adjudication code and it **caught an
off-family closure that both the closure residual and the V_inf anchors passed**
(``russell-ch4-5.30ggF3``: residual 0.000, emerged V_inf within 0.03 km/s of the
sourced anchors, yet an Earth node demanding 141.2 deg where 86.6 deg is
available). Flyby TURN GEOMETRY is an evidence axis independent of both closure
continuity and V_inf magnitude — the classic
``[[feedback_orbit_closure_discipline]]`` "it closed!" blind spot. It reproduced
the published TR to **0.001-0.024 on 7 of 8** `#820` rows.

Standing use: wherever a closure is produced for a catalogue row that publishes
``invariants.turn_ratio``, measure it and compare. A measured TR far BELOW the
published one means the closure is turning somewhere the sourced cycler does not
(off-family). A measured TR < 1 is **not automatically a defect**: for a row whose
PUBLISHED TR is < 1 (Russell Table 4.13's near-ballistic rows, e.g.
``russell-ch4-6.44Gg3`` at 0.95) reproducing TR < 1 is *evidence the closure is on
the sourced family* — the published cycler itself needs a small powered nudge.
That distinction is exactly why the raw ``bend_feasible`` boolean is not a
sufficient gate on its own (`#829`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import (
    BallisticClosureResult,
    _bend_deg,
    _max_bend_deg,
    _vinf_nodes,
)

DAY_S = 86400.0

# Agreement tolerance for "the measured TR reproduces the published TR".
# Russell prints TR to 2 decimals; #826's observed spread over the 7 reproducing
# rows was 0.001-0.024, so 0.05 has teeth (it catches a wrong node, a wrong flyby
# body or a mis-posed genome) without pinning an exact computed value.
TURN_RATIO_TOL = 0.05

# A node the closure turns through by less than this is UNCONSTRAINED: it imposes
# no turn at all, so it can never be the binding flyby Russell's delta_MAX refers
# to. (Under the #820 designated-arc posing the Mars node is exactly this: legs 0
# and 1 are the same conic split at the Mars encounter, so required_bend ~ 0.)
UNCONSTRAINED_BEND_DEG = 1e-9


@dataclass(frozen=True)
class FlybyTurn:
    """One intermediate flyby's required vs maximum ballistic turn."""

    index: int
    """Position of the flyby in the encounter ``sequence``."""
    body: str
    required_bend_deg: float
    """Angle between the incoming and outgoing V_inf vectors."""
    max_bend_deg: float
    """V_inf-limited maximum turn at the body's ``safe_alt_km`` periapsis."""
    ratio: float
    """``max_bend / required_bend`` — Russell's TR for this node.
    ``inf`` when the node is unconstrained (required turn ~ 0)."""
    unconstrained: bool
    """True when the required turn is below :data:`UNCONSTRAINED_BEND_DEG`."""

    @property
    def feasible(self) -> bool:
        """The turn fits inside the ballistic maximum (``ratio >= 1``)."""
        return self.required_bend_deg <= self.max_bend_deg


@dataclass(frozen=True)
class TurnRatioReport:
    """Measured turn ratio of a closure, per node and overall."""

    flybys: tuple[FlybyTurn, ...]
    turn_ratio: float
    """``min`` over the flybys of ``max_bend / required_bend`` — the measured
    reproduction of Russell's published TR. ``inf`` if no node imposes a turn."""
    binding_index: int | None
    """``sequence`` index of the flyby setting :attr:`turn_ratio` (``None`` if
    every node is unconstrained)."""

    @property
    def all_feasible(self) -> bool:
        """Every intermediate flyby's turn is ballistically attainable.

        Mirrors :attr:`BallisticClosureResult.bend_feasible` — asserting the two
        agree is a cheap cross-check that this instrument sees the same nodes the
        corrector's own post-hoc feasibility check saw.
        """
        return all(f.feasible for f in self.flybys)

    def agrees_with_published(self, published: float, *, tol: float = TURN_RATIO_TOL) -> bool:
        """Does the measured TR reproduce a row's published ``turn_ratio``?"""
        return abs(self.turn_ratio - float(published)) <= tol

    def summary(self) -> str:
        """One-line per-node digest, e.g. ``M1:inf E2:0.613 E3:1.255``."""
        return " ".join(
            f"{f.body}{f.index}:" + ("inf" if f.unconstrained else f"{f.ratio:.3f}")
            for f in self.flybys
        )


def measure_turn_ratio(
    nodes: dict[str, np.ndarray],
    sequence: tuple[str, ...],
    *,
    rp_factors: dict[str, float] | None = None,
) -> TurnRatioReport:
    """Measure Russell's turn ratio from a closure's V_inf node vectors.

    ``nodes`` is the ``b{i}_in`` / ``b{i}_out`` mapping produced by
    :func:`cyclerfinder.search.correct._vinf_nodes`; only the INTERMEDIATE
    encounters (``1 .. len(sequence) - 2``) are flybys — the chain's ends carry a
    single leg each and are joined by the periodicity residual, not by a turn.
    """
    flybys: list[FlybyTurn] = []
    for i in range(1, len(sequence) - 1):
        v_in, v_out = nodes[f"b{i}_in"], nodes[f"b{i}_out"]
        required = _bend_deg(v_in, v_out)
        max_turn = _max_bend_deg(float(np.linalg.norm(v_in)), sequence[i], rp_factors)
        unconstrained = required < UNCONSTRAINED_BEND_DEG
        flybys.append(
            FlybyTurn(
                index=i,
                body=sequence[i],
                required_bend_deg=required,
                max_bend_deg=max_turn,
                ratio=float("inf") if unconstrained else max_turn / required,
                unconstrained=unconstrained,
            )
        )
    binding: int | None = None
    ratio = float("inf")
    for f in flybys:
        if f.ratio < ratio:
            ratio, binding = f.ratio, f.index
    return TurnRatioReport(flybys=tuple(flybys), turn_ratio=ratio, binding_index=binding)


def closure_turn_ratio(
    result: BallisticClosureResult,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    slack_leg: int,
    period_sec: float,
    ephem: Ephemeris,
    rp_factors: dict[str, float] | None = None,
) -> TurnRatioReport:
    """:func:`measure_turn_ratio` for a :class:`BallisticClosureResult`.

    The result object carries only per-encounter V_inf MAGNITUDES, so the node
    vectors are re-extracted from its converged ``(t0, tof_days)`` with the same
    genome the corrector ran. Raises whatever ``_vinf_nodes`` raises if the
    converged point is not re-solvable (a caller wanting a soft failure should
    check ``result.converged`` first).
    """
    free_tof = tuple(float(v) for i, v in enumerate(result.tof_days) if i != slack_leg)
    nodes = _vinf_nodes(
        sequence=sequence,
        per_leg_revs=per_leg_revs,
        per_leg_branch=per_leg_branch,
        t0_sec=float(result.t0_sec),
        free_tof_days=free_tof,
        slack_leg=slack_leg,
        period_days=period_sec / DAY_S,
        ephem=ephem,
    )
    return measure_turn_ratio(nodes, sequence, rp_factors=rp_factors)
