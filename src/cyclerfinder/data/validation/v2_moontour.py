"""V2 long-span bounded-drift moontour gauntlet (#306 Phase 2 / #330).

Spec reference
--------------
* §14 V2-ballistic — "≥3 continuous laps; **bounded** drift in the dynamic
  rotating frame (tolerant of geometric breathing), evaluated **in the row's
  defining model**".

For a moontour (a Lambert-leg patched-conic tour through moons of a primary,
e.g. the #327 Umbriel-Oberon-Umbriel SILVER), the row's defining model is the
circular-coplanar Keplerian moon ephemeris + planet-frame Lambert legs. There
is NO (state0, period) periodic IC — the candidate is described by
``(sequence, V_inf-tuple-at-each-encounter, ToF-per-leg, rel_offset_deg)``.

Phase 1's V2 (``v2_3d.py``) is for periodic CR3BP orbits: it propagates a
single 6D IC for ``n_cycles * period`` and gates the position drift at each
cycle boundary. That ``run_v2_3d`` does NOT apply to moontours — there is no
single 6D IC to propagate.

Phase 2's V2 (this module) re-solves the Lambert legs over ``n_cycles``
consecutive cycles with moon longitudes advanced through their natural
Keplerian motion across the cycle. The "drift" is the inter-cycle rendezvous
defect: how much the encounter geometry / V_inf-continuity of cycle k departs
from cycle 0.

Discipline distinction vs Phase 1
---------------------------------
* Phase 1: ``state0 + period`` -> propagate the 6D IC for ``n_cycles``.
* Phase 2 (this): ``sequence + V_inf-tuple + ToFs + rel_offset_deg`` ->
  re-solve Lambert legs ``x n_cycles``, with moon ephemerides advanced each
  cycle.

PASS criterion: ``n_cycles >= n_cycles_min`` complete (every leg's Lambert
converged in every cycle) AND every cycle's V_inf-continuity residual stays
within ``closure_floor_kms`` AND the max inter-cycle rendezvous drift stays
below ``drift_floor_kms``.

The drift floor (50,000 km) mirrors the Phase-1 same-model floor
:data:`cyclerfinder.data.validation.v2_3d.V2_DRIFT_FLOOR_KMS`. The per-cycle
closure floor (0.05 km/s) matches the #285 / #312 SILVER closure-gate
threshold.

Composition map
---------------
This module composes directly on
:func:`cyclerfinder.core.lambert.lambert` (the same Lambert solver the #327
SILVER closed under) and the circular-coplanar moon-state primitive
:func:`cyclerfinder.search.discovery_campaign._moon_state`. The leg geometry
mirrors :meth:`RepeatedMoonTarget._close_one_phasing` verbatim; the cycle
loop simply repeats with the moon-longitude initial conditions advanced by
the cumulative ToF of the previous cycles.

Discipline
----------
* NO catalogue writeback. A V2 pass alone DOES NOT admit a row to
  ``catalogue.yaml``. V3 + V4 + #329 lit-check still gate.
* The 3-cycle minimum (``n_cycles_min=3``) is spec §14, NOT test-tunable.
* The V2 verdict is whatever the math says — if a SILVER candidate's drift
  blows past the floor at cycle 2, the verdict is FAIL and the candidate is
  honestly recorded as unstable on the 3-cycle horizon (the
  ``feedback_orbit_closure_discipline`` "clean negative is success" rule).
* The Lambert geometry is the same machinery #327 closed under: the inputs
  here are OUR computation, and what's being asserted is CONVERGENCE +
  DRIFT-BOUND, never a specific number.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day, _moon_state

V2_MOONTOUR_N_CYCLES_MIN: Final[int] = 3
"""Spec §14 V2-ballistic minimum: ``>= 3`` continuous laps. Spec-fixed."""

V2_MOONTOUR_DRIFT_FLOOR_KMS: Final[float] = 50_000.0
"""Default same-model inter-cycle rendezvous drift floor in km.

Mirrors :data:`cyclerfinder.data.validation.v2_3d.V2_DRIFT_FLOOR_KMS` so the
moontour bar is the same as the 3D-periodic-orbit bar. For a moontour the
"drift" is the final-encounter position offset between cycle k and cycle 0
(the patched-conic equivalent of the 6D CR3BP cumulative position drift)."""

V2_MOONTOUR_CLOSURE_FLOOR_KMS: Final[float] = 0.05
"""Per-cycle V_inf-continuity residual floor (matches #285 / #312 SILVER gate)."""


@dataclass(frozen=True)
class MoontourCycleVerdict:
    """Per-cycle verdict for the moontour V2 Lambert-relegs gauntlet.

    Attributes
    ----------
    cycle_index:
        Zero-indexed cycle number (0, 1, 2, ...).
    converged_legs:
        Number of Lambert legs that closed inside this cycle.
    n_legs:
        Total Lambert legs in one cycle (``len(sequence) - 1``).
    rendezvous_drift_kms:
        Position offset of THIS cycle's final encounter vs the SAME encounter
        in cycle 0 (km). Captures how far the moons have drifted out of the
        cycle-0 phasing by cycle k. 0.0 for ``cycle_index == 0`` by
        construction.
    rendezvous_drift_seconds:
        Timing offset of this cycle vs cycle 0 — the cumulative ToF up to
        this cycle's final encounter minus the cumulative ToF at cycle 0.
        Always 0.0 here (per-cycle ToFs are held fixed at the SILVER's
        stored ToFs); kept on the dataclass for forward-compat with a
        time-retargeted variant.
    closure_residual_kms:
        This cycle's V_inf-continuity residual (max over interior flybys +
        closed-cycle anchor wrap, km/s). Same definition as
        :meth:`RepeatedMoonTarget._close_one_phasing`'s ``worst``.
    notes:
        Free-form audit string.
    """

    cycle_index: int
    converged_legs: int
    n_legs: int
    rendezvous_drift_kms: float
    rendezvous_drift_seconds: float
    closure_residual_kms: float
    notes: str = ""


@dataclass(frozen=True)
class V2MoontourVerdict:
    """Frozen V2 verdict for a moontour Lambert-relegs gauntlet candidate.

    Attributes
    ----------
    candidate_id:
        Identifier of the candidate (carried for the audit trail).
    sequence:
        Body sequence of one cycle (e.g. ``("Umbriel", "Oberon", "Umbriel")``).
    n_cycles_requested:
        Cycles the caller asked for.
    n_cycles_completed:
        Cycles where every Lambert leg converged. Equals
        ``n_cycles_requested`` unless a leg failed mid-flight.
    per_cycle:
        One :class:`MoontourCycleVerdict` per attempted cycle.
    max_drift_kms:
        Max of ``per_cycle[k].rendezvous_drift_kms`` for ``k >= 1``.
    max_drift_seconds:
        Max of ``per_cycle[k].rendezvous_drift_seconds`` for ``k >= 1``.
    max_closure_residual_kms:
        Max of ``per_cycle[k].closure_residual_kms`` across all completed cycles.
    drift_floor_kms:
        The bar the rendezvous drift was held against.
    closure_floor_kms:
        The bar each cycle's V_inf-continuity residual was held against.
    n_cycles_min:
        Spec §14 minimum cycles (3). Stored for audit.
    passes_v2:
        ``n_cycles_completed >= n_cycles_min AND max_drift_kms <=
        drift_floor_kms AND max_closure_residual_kms <= closure_floor_kms``.
        The headline boolean.
    notes:
        Free-form audit string.
    """

    candidate_id: str
    sequence: tuple[str, ...]
    n_cycles_requested: int
    n_cycles_completed: int
    per_cycle: tuple[MoontourCycleVerdict, ...]
    max_drift_kms: float
    max_drift_seconds: float
    max_closure_residual_kms: float
    drift_floor_kms: float
    closure_floor_kms: float
    n_cycles_min: int
    passes_v2: bool
    notes: str = ""


def _resolve_primary(system: cr3bp.CR3BPSystem | None, sequence: tuple[str, ...]) -> str:
    """Resolve the primary body name from the system or from the moon sequence.

    The moontour CR3BP system carries the primary as ``system.primary``. If
    the caller passes ``None``, fall back to the moon registry: the sequence's
    first moon's ``SATELLITES[moon].primary`` is authoritative.
    """
    if system is not None and system.primary:
        # ``CR3BPSystem.primary`` is stored as a lowercase string in some
        # callers; the registry uses titlecase. Normalise to the registry's
        # canonical form by matching case-insensitively.
        target = system.primary.strip().lower()
        for name in PRIMARIES:
            if name.lower() == target:
                return name
        raise ValueError(f"unknown primary {system.primary!r} from system")
    if not sequence:
        raise ValueError("empty sequence; cannot resolve primary")
    head = sequence[0]
    if head not in SATELLITES:
        raise ValueError(f"unknown moon {head!r} in sequence")
    return SATELLITES[head].primary


def _moon_constants(primary: str, sequence: tuple[str, ...]) -> dict[str, tuple[float, float]]:
    """``{moon: (sma_km, mean_motion_rad_day)}`` for every moon in the sequence."""
    mu = PRIMARIES[primary]
    out: dict[str, tuple[float, float]] = {}
    for moon in sequence:
        if moon not in SATELLITES:
            raise ValueError(f"unknown moon {moon!r}")
        sat = SATELLITES[moon]
        if sat.primary != primary:
            raise ValueError(
                f"moon {moon!r} orbits {sat.primary!r}, not the resolved primary {primary!r}"
            )
        out[moon] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))
    return out


def _cycle_residual(
    *,
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    theta_base: dict[str, float],
    t_cycle_offset_days: float,
    consts: dict[str, tuple[float, float]],
    mu: float,
    n_revs: tuple[int, ...] | None,
) -> tuple[bool, float, list[tuple[np.ndarray, np.ndarray]]]:
    """Re-solve all Lambert legs of one cycle; return (converged, residual, states).

    Mirrors :meth:`RepeatedMoonTarget._close_one_phasing` verbatim for the
    geometry; only the moon-longitude phase is shifted by
    ``n * t_cycle_offset_days`` to advance the ephemeris through cycle k.

    Returns:
        ``converged``: every leg's Lambert had a solution at the requested n_rev
        (or the lowest-energy solution at default n_rev=0 if ``n_revs`` is None).
        ``residual``: worst per-flyby V_inf-magnitude continuity defect
        (km/s, including the closed-cycle anchor wrap).
        ``states``: planet-frame (pos km, vel km/s) at every encounter.
    """
    n_legs = len(sequence) - 1
    if n_revs is None:
        n_revs_used: tuple[int, ...] = tuple(0 for _ in range(n_legs))
    else:
        if len(n_revs) != n_legs:
            raise ValueError(f"n_revs length {len(n_revs)} != n_legs {n_legs}")
        n_revs_used = tuple(n_revs)

    # Cumulative epochs at each encounter, in days, measured from cycle k's start.
    epochs_days = [0.0]
    for tof in leg_tofs_days:
        epochs_days.append(epochs_days[-1] + tof)

    # Each moon's longitude at this encounter: theta_base + n * (t_cycle_offset + t).
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for moon, t in zip(sequence, epochs_days, strict=True):
        sma, n_rad_day = consts[moon]
        states.append(_moon_state(theta_base[moon], n_rad_day, t_cycle_offset_days + t, sma, mu))

    vinf_in: list[float | None] = [None] * len(sequence)
    vinf_out: list[float | None] = [None] * len(sequence)
    for k in range(n_legs):
        r_a, v_a = states[k]
        r_b, _ = states[k + 1]
        nrev = max(0, n_revs_used[k])
        sols = _lambert(r_a, r_b, leg_tofs_days[k] * DAY_S, mu=mu, max_revs=nrev)
        wanted = [s for s in sols if s.n_revs == n_revs_used[k]]
        if not wanted:
            return (False, math.inf, states)
        v_a_captured = v_a
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a_captured)))
        _, v_b_moon = states[k + 1]
        vinf_out[k] = float(np.linalg.norm(best.v1 - v_a))
        vinf_in[k + 1] = float(np.linalg.norm(best.v2 - v_b_moon))

    worst = 0.0
    for k in range(len(sequence)):
        vi = vinf_in[k]
        vo = vinf_out[k]
        if vi is not None and vo is not None:
            worst = max(worst, abs(vi - vo))
    # Closed-cycle anchor wrap (same definition as #259 fix C).
    wrap_out = vinf_out[0]
    wrap_in = vinf_in[-1]
    if wrap_out is not None and wrap_in is not None:
        worst = max(worst, abs(wrap_out - wrap_in))
    return (True, worst, states)


def run_v2_moontour(
    candidate_id: str,
    sequence: tuple[str, ...],
    vinf_tuple_kms: tuple[float, ...],
    leg_tofs_days: tuple[float, ...],
    rel_offset_deg: float,
    system: cr3bp.CR3BPSystem | None,
    *,
    n_cycles: int = V2_MOONTOUR_N_CYCLES_MIN,
    drift_floor_kms: float = V2_MOONTOUR_DRIFT_FLOOR_KMS,
    closure_floor_kms: float = V2_MOONTOUR_CLOSURE_FLOOR_KMS,
    n_revs: tuple[int, ...] | None = None,
    phase0_deg: float = 0.0,
    notes: str = "",
) -> V2MoontourVerdict:
    """Run V2 for a moontour: re-solve Lambert legs over ``n_cycles``.

    Pipeline:
      1. Build the moon ephemerides (circular-coplanar Keplerian) for every
         moon in ``sequence``, anchored at ``phase0_deg`` for the FIRST moon
         and with relative offset ``rel_offset_deg`` between the SILVER's
         two distinct moons (mirrors the #327 SILVER's deterministic offset
         convention).
      2. For each cycle ``k = 0, 1, ..., n_cycles - 1``:
           a. Re-solve all Lambert legs with the moon longitudes advanced by
              ``k * cycle_period_days`` where ``cycle_period_days = sum(leg_tofs_days)``.
           b. Record the per-cycle V_inf-continuity residual.
           c. Record the rendezvous drift = ``||r_final_k - r_final_0||``
              (km), where ``r_final_k`` is the planet-frame position of the
              cycle's final encounter (closes the loop back at the anchor).
      3. Verdict: PASS iff every cycle converged AND
         ``max_drift <= drift_floor`` AND ``max_closure <= closure_floor``.

    Parameters
    ----------
    candidate_id:
        Identifier carried into the verdict.
    sequence:
        Body sequence of one cycle (``(M_0, M_1, ..., M_0)``); first and last
        must match (closed cycle).
    vinf_tuple_kms:
        Stored V_inf tuple at each encounter (carried for the audit trail
        and the ``notes`` field; the V2 driver does NOT enforce it as a
        constraint — V2 re-solves the legs from scratch).
    leg_tofs_days:
        Per-leg time-of-flight in days. Length = ``len(sequence) - 1``.
    rel_offset_deg:
        Relative-offset between the two distinct moons in the SILVER's basin
        (the #327 gate-passing offset is 180°). Used to set
        ``theta0[moon_b] = phase0 + rel_offset_deg``. For sequences with more
        than 2 distinct moons the offset is applied to the second distinct
        moon in registry-sorted order.
    system:
        CR3BP system. Used ONLY to resolve the primary name (the moontour
        runs in the planet frame, not the CR3BP rotating frame). May be
        ``None``; in that case the primary is resolved from the first moon's
        registry entry.
    n_cycles:
        Cycles to attempt. Must be ``>= V2_MOONTOUR_N_CYCLES_MIN``. Default 3.
    drift_floor_kms:
        Bar against which ``max_drift_kms`` is held. Default 50,000 km.
    closure_floor_kms:
        Bar against which each cycle's V_inf-continuity residual is held.
        Default 0.05 km/s (matches the #285 / #312 SILVER gate).
    n_revs:
        Per-leg revolution count. If ``None``, every leg uses ``n_revs=0``
        (the SILVER's stored ``n_rev=(1,1)`` should be passed explicitly).
    phase0_deg:
        Initial longitude of the first moon at cycle 0, degrees. Default 0°.
        For reproducing a SILVER row precisely, callers pass the row's stored
        ``phase0_deg`` (the #327 basin-floor record stores 29.999...°).
    notes:
        Free-form audit note.

    Returns
    -------
    V2MoontourVerdict
        ``passes_v2`` is the headline.

    Notes
    -----
    The cycle ToFs are held FIXED at ``leg_tofs_days`` for every cycle (this
    is the SILVER's stored geometry — re-solving the time grid each cycle
    would be a different gate: V2-powered-retarget, not V2-ballistic). The
    moons advance through their natural ephemerides; the rendezvous drift
    measures whether the moons remain in cycle-0 phasing after ``n_cycles``.

    A V2 PASS does NOT admit to the catalogue. V3 + V4 + #329 follow.
    """
    if not sequence:
        raise ValueError("empty sequence")
    if sequence[0] != sequence[-1]:
        raise ValueError(f"moontour sequence must be CLOSED (first == last); got {sequence!r}")
    n_legs = len(sequence) - 1
    if len(leg_tofs_days) != n_legs:
        raise ValueError(f"leg_tofs_days length {len(leg_tofs_days)} != n_legs {n_legs}")
    if len(vinf_tuple_kms) != len(sequence):
        raise ValueError(
            f"vinf_tuple_kms length {len(vinf_tuple_kms)} != len(sequence) {len(sequence)}"
        )
    if any(tof <= 0.0 for tof in leg_tofs_days):
        raise ValueError(f"leg_tofs_days must be positive; got {leg_tofs_days!r}")
    if n_cycles < V2_MOONTOUR_N_CYCLES_MIN:
        raise ValueError(
            f"V2-moontour requires n_cycles >= {V2_MOONTOUR_N_CYCLES_MIN} "
            f"(spec §14); got {n_cycles}"
        )
    if drift_floor_kms <= 0.0:
        raise ValueError(f"drift_floor_kms must be > 0; got {drift_floor_kms}")
    if closure_floor_kms <= 0.0:
        raise ValueError(f"closure_floor_kms must be > 0; got {closure_floor_kms}")

    primary = _resolve_primary(system, sequence)
    consts = _moon_constants(primary, sequence)
    mu = PRIMARIES[primary]

    # Anchor moon longitudes for cycle 0. The convention mirrors the #327
    # SILVER: phase0_deg sets the FIRST distinct moon's longitude, and the
    # SECOND distinct moon (in registry-sorted order) is offset by
    # ``rel_offset_deg``. Any further distinct moons get evenly-spaced
    # additional offsets so a 3-moon sequence still phases deterministically.
    phase0_rad = math.radians(phase0_deg)
    rel_off_rad = math.radians(rel_offset_deg)
    distinct_moons = tuple(sorted({m for m in sequence}))
    if len(distinct_moons) == 1:
        raise ValueError(f"moontour requires >= 2 distinct moons; got {distinct_moons!r}")
    theta_base: dict[str, float] = {}
    for j, moon in enumerate(distinct_moons):
        if j == 0:
            theta_base[moon] = phase0_rad
        elif j == 1:
            theta_base[moon] = phase0_rad + rel_off_rad
        else:
            # Spread any extra distinct moons evenly around (preserves
            # determinism while not collapsing to colocation).
            theta_base[moon] = (
                phase0_rad + rel_off_rad + 2.0 * math.pi * (j - 1) / len(distinct_moons)
            )

    cycle_period_days = float(sum(leg_tofs_days))

    per_cycle: list[MoontourCycleVerdict] = []
    cycle_zero_final_pos_km: np.ndarray | None = None
    n_completed = 0
    max_drift_kms = 0.0
    max_drift_seconds = 0.0
    max_closure = 0.0

    for k in range(n_cycles):
        t_offset_days = k * cycle_period_days
        converged, residual, states = _cycle_residual(
            sequence=sequence,
            leg_tofs_days=leg_tofs_days,
            theta_base=theta_base,
            t_cycle_offset_days=t_offset_days,
            consts=consts,
            mu=mu,
            n_revs=n_revs,
        )
        if not converged:
            per_cycle.append(
                MoontourCycleVerdict(
                    cycle_index=k,
                    converged_legs=0,
                    n_legs=n_legs,
                    rendezvous_drift_kms=float("inf"),
                    rendezvous_drift_seconds=0.0,
                    closure_residual_kms=float("inf"),
                    notes="Lambert failed at least one leg in this cycle",
                )
            )
            break

        final_pos_km = states[-1][0]
        if k == 0:
            cycle_zero_final_pos_km = final_pos_km.copy()
            drift_kms = 0.0
        else:
            assert cycle_zero_final_pos_km is not None
            drift_kms = float(np.linalg.norm(final_pos_km - cycle_zero_final_pos_km))
            max_drift_kms = max(max_drift_kms, drift_kms)
        max_closure = max(max_closure, residual)

        per_cycle.append(
            MoontourCycleVerdict(
                cycle_index=k,
                converged_legs=n_legs,
                n_legs=n_legs,
                rendezvous_drift_kms=drift_kms,
                rendezvous_drift_seconds=0.0,
                closure_residual_kms=residual,
            )
        )
        n_completed += 1

    passes_v2 = bool(
        n_completed >= V2_MOONTOUR_N_CYCLES_MIN
        and max_drift_kms <= drift_floor_kms
        and max_closure <= closure_floor_kms
    )

    return V2MoontourVerdict(
        candidate_id=candidate_id,
        sequence=tuple(sequence),
        n_cycles_requested=int(n_cycles),
        n_cycles_completed=int(n_completed),
        per_cycle=tuple(per_cycle),
        max_drift_kms=float(max_drift_kms),
        max_drift_seconds=float(max_drift_seconds),
        max_closure_residual_kms=float(max_closure),
        drift_floor_kms=float(drift_floor_kms),
        closure_floor_kms=float(closure_floor_kms),
        n_cycles_min=V2_MOONTOUR_N_CYCLES_MIN,
        passes_v2=passes_v2,
        notes=notes,
    )


__all__ = [
    "V2_MOONTOUR_CLOSURE_FLOOR_KMS",
    "V2_MOONTOUR_DRIFT_FLOOR_KMS",
    "V2_MOONTOUR_N_CYCLES_MIN",
    "MoontourCycleVerdict",
    "V2MoontourVerdict",
    "run_v2_moontour",
]
