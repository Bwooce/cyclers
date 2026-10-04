"""Demanded-turn wiring for the moon-tour validation lane (#888).

The V2 / V3 / V4 / V4-strict moon-tour tiers re-solve one Lambert leg per
encounter pair and, before #888, compared only V-infinity MAGNITUDES at the
encounters. Every leg restarted from its own Lambert departure velocity, so no
incoming V-infinity vector was ever compared with the outgoing one, and six
catalogued Uranian rows passed every tier while demanding 1.8 to 28 times the
bend the moons can supply (``docs/notes/2026-10-04-888-demanded-turn-gate.md``).

This module is the glue between those tiers and the gate in
:mod:`cyclerfinder.verify.turn_gate` (which it does not modify):

* each tier's ``_cycle_*`` function, when handed a :class:`CycleTurnRecord`,
  stores every leg's departure and arrival V-infinity VECTORS (spacecraft
  Lambert velocity minus the moon's velocity at that encounter, in the tier's
  own planet-centred frame) and the departure V-infinity of the NEXT cycle's
  first leg, i.e. the leg that leaves the anchor after the closing flyby;
* :func:`chain_encounters` turns that into one :class:`Encounter` per
  intermediate flyby plus the anchor WRAP flyby (arrival of the last leg
  against the next cycle's first departure, at the same epoch and the same
  moon state, so no frame rotation is needed);
* :func:`gate_cycle_record` applies :func:`demanded_turn_gate` with the
  registry altitude floors (or a uniform override) and stores the report;
* :func:`summarize_turn_gate` folds the per-cycle reports into the verdict
  fields every tier now carries (``turn_feasible``, ``worst_turn_ratio``,
  ``turn_failure_reason``), and every tier's ``passes_*`` requires
  ``turn_feasible``.

The verdict used is the gate's ``turn_feasible`` (demanded turn no larger than
the available bend at every flyby), NOT its stricter ``ballistic`` flag: the
V-infinity magnitude continuity is already held to each tier's own closure
residual, and ``ballistic`` at its 1e-3 km/s default would reject published
ballistic cyclers whose printed parameters carry a 0.011 km/s spread (#888 note,
section 7).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.verify.turn_gate import (
    Encounter,
    EncounterTurn,
    TurnGateReport,
    demanded_turn_gate,
)

Vec = NDArray[np.float64]

#: ``turn_failure_reason`` of a verdict built without the gate being run.
TURN_GATE_NOT_EVALUATED = "turn gate not evaluated"


@dataclass(frozen=True)
class LegVinf:
    """V-infinity vectors of one leg (km/s, the tier's planet-centred frame).

    ``depart``: spacecraft departure velocity minus the departure moon's
    velocity. ``arrive``: spacecraft arrival velocity minus the arrival moon's
    velocity.
    """

    depart: Vec
    arrive: Vec


@dataclass
class CycleTurnRecord:
    """Collector a tier's ``_cycle_*`` function fills when one is passed in.

    ``legs`` holds one :class:`LegVinf` per leg of the cycle, in order;
    ``wrap_depart`` is the departure V-infinity of the next cycle's first leg
    (``None`` until known); ``report`` is set by :func:`gate_cycle_record`;
    ``error`` explains why the gate could not be applied (no flyable wrap leg,
    a zero V-infinity, ...), in which case the cycle counts as NOT turn-feasible.
    """

    legs: list[LegVinf] = field(default_factory=list)
    wrap_depart: Vec | None = None
    report: TurnGateReport | None = None
    error: str = ""

    @property
    def turn_feasible(self) -> bool:
        return self.report is not None and self.report.turn_feasible

    @property
    def encounters(self) -> tuple[EncounterTurn, ...]:
        return () if self.report is None else self.report.encounters


def leg_vinf(v1_sc: Vec, v_moon_depart: Vec, v2_sc: Vec, v_moon_arrive: Vec) -> LegVinf:
    """:class:`LegVinf` from a leg's spacecraft velocities and the two moon velocities."""
    return LegVinf(
        depart=np.asarray(v1_sc, dtype=np.float64) - np.asarray(v_moon_depart, dtype=np.float64),
        arrive=np.asarray(v2_sc, dtype=np.float64) - np.asarray(v_moon_arrive, dtype=np.float64),
    )


def chain_encounters(
    sequence: Sequence[str],
    legs: Sequence[LegVinf],
    wrap_depart: Sequence[float] | Vec,
    *,
    alt_floor_km: float | None = None,
) -> list[Encounter]:
    """Every flyby of one closed cycle ``sequence[0] -> ... -> sequence[-1]``.

    Encounters ``1 .. n_legs - 1`` compare leg ``k - 1``'s arrival with leg
    ``k``'s departure at ``sequence[k]`` (label ``"e<k>"``); the closing
    encounter at ``sequence[-1]`` compares the last leg's arrival with
    ``wrap_depart``, the departure of the next cycle's first leg (label
    ``"wrap"``). ``alt_floor_km=None`` uses each body's registry floor.
    """
    n_legs = len(sequence) - 1
    if len(legs) != n_legs:
        raise ValueError(f"{len(legs)} legs for a {len(sequence)}-encounter sequence")
    if sequence[0] != sequence[-1]:
        raise ValueError(f"sequence must be closed (first == last); got {tuple(sequence)!r}")
    out = [
        Encounter.for_body(
            sequence[k],
            legs[k - 1].arrive,
            legs[k].depart,
            alt_floor_km=alt_floor_km,
            label=f"e{k}",
        )
        for k in range(1, n_legs)
    ]
    out.append(
        Encounter.for_body(
            sequence[-1], legs[-1].arrive, wrap_depart, alt_floor_km=alt_floor_km, label="wrap"
        )
    )
    return out


def chain_turn_report(
    sequence: Sequence[str],
    legs: Sequence[LegVinf],
    wrap_depart: Sequence[float] | Vec,
    *,
    alt_floor_km: float | None = None,
) -> TurnGateReport:
    """:func:`demanded_turn_gate` over :func:`chain_encounters`."""
    return demanded_turn_gate(
        chain_encounters(sequence, legs, wrap_depart, alt_floor_km=alt_floor_km)
    )


def gate_cycle_record(
    record: CycleTurnRecord,
    sequence: Sequence[str],
    *,
    alt_floor_km: float | None = None,
) -> None:
    """Apply the gate to a filled record (sets ``record.report`` or ``record.error``)."""
    if record.error:
        return
    if len(record.legs) != len(sequence) - 1:
        record.error = f"only {len(record.legs)} of {len(sequence) - 1} legs were solved"
        return
    if record.wrap_depart is None:
        record.error = "no flyable leg leaves the anchor after the closing flyby (wrap leg)"
        return
    try:
        record.report = chain_turn_report(
            sequence, record.legs, record.wrap_depart, alt_floor_km=alt_floor_km
        )
    except ValueError as exc:  # e.g. a zero V-infinity: no flyby is defined
        record.error = f"turn gate not applicable: {exc}"


@dataclass(frozen=True)
class TurnGateSummary:
    """Verdict-level fold of the per-cycle turn gates."""

    turn_feasible: bool
    worst_turn_ratio: float
    turn_failure_reason: str


def _describe(e: EncounterTurn, cycle_index: int) -> str:
    return (
        f"demanded turn exceeds the available bend at {e.body} ({e.label}) in cycle "
        f"{cycle_index}: {e.demanded_turn_deg:.1f} deg demanded, {e.available_bend_deg:.1f} "
        f"deg available at the {e.alt_floor_km:.0f} km floor (ratio {e.ratio:.2f}, "
        f"required periapsis altitude {e.required_alt_km:.0f} km)"
    )


def summarize_turn_gate(records: Sequence[CycleTurnRecord]) -> TurnGateSummary:
    """Fold the gated records of the COMPLETED cycles into verdict fields.

    ``turn_feasible`` requires at least one record and every record feasible.
    The reason names the encounter with the largest demanded/available ratio
    (or the first record that could not be gated).
    """
    if not records:
        return TurnGateSummary(False, math.inf, "no completed cycle to gate")
    for k, rec in enumerate(records):
        if rec.report is None:
            return TurnGateSummary(False, math.inf, f"cycle {k}: {rec.error or 'not gated'}")
    worst: tuple[float, int, EncounterTurn] | None = None
    for k, rec in enumerate(records):
        assert rec.report is not None
        b = rec.report.binding
        if b is not None and (worst is None or b.ratio > worst[0]):
            worst = (b.ratio, k, b)
    feasible = all(rec.turn_feasible for rec in records)
    if worst is None:
        return TurnGateSummary(feasible, 0.0, "")
    ratio, k, enc = worst
    return TurnGateSummary(feasible, float(ratio), "" if feasible else _describe(enc, k))


__all__ = [
    "TURN_GATE_NOT_EVALUATED",
    "CycleTurnRecord",
    "LegVinf",
    "TurnGateSummary",
    "chain_encounters",
    "chain_turn_report",
    "gate_cycle_record",
    "leg_vinf",
    "summarize_turn_gate",
]
