"""Phase-matching: idealised cycler → real launch-window dates.

This module implements a slice of spec §12.1 "the idealized→ephemeris
bridge". Given an idealised cycler (orbit defined in the circular-coplanar
model OR a catalogue YAML entry), it finds the next N calendar dates on the
real JPL ephemeris where the planetary geometry matches what the cycler
requires for a ballistic injection.

What this module does:

- :func:`phase_signature` extracts a cycler's required leg ToFs and per-leg
  V∞ targets — the "fingerprint" of geometry the cycler needs.
- :func:`phase_signature_from_catalogue_entry` does the same starting from
  the catalogue YAML dict shape (so the site repo can use it without
  reconstructing a full :class:`~cyclerfinder.model.cycler.Cycler`).
- :func:`find_real_windows` grid-scans an astropy ephemeris across a date
  range, solving Lambert at each candidate departure date and computing the
  V∞ mismatch vs the cycler's signature. Returns the N lowest-mismatch
  dates with their actual V∞ at departure.

What this module does NOT do (those remain in the full M6):

- Multi-lap propagation through the rotating frame (closure-residual
  budget — M6a).
- Ephemeris-mode TCM minimisation over a 3-5 lap horizon (M6b).
- C₃ and time-of-flight cost columns on the launch-window page.

In practical terms: this module produces *honest real dates* at which a
ballistic injection geometry exists. It does not produce a TCM-budgeted
mission design — that requires the M6b optimisation slice that wasn't
brought forward in this 2026-06-01 slice.

References
----------
Spec §12.1 (phase-match bridge), §12(a) (idealised vs ephemeris modes).
Plan: docs/phases/m4-enumeration-scoring/plan.md describes the M4 boundary;
this module's launch-windows slice was implemented ahead of the full M6
plan.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from math import inf
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError, lambert
from cyclerfinder.data.catalog import _segments_as_legs

# Reference epoch for converting between datetime and seconds-since-J2000.
_J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def _dt_to_t_sec(dt: datetime) -> float:
    """Convert a UTC datetime to seconds since the J2000 reference epoch."""
    if dt.tzinfo is None:
        raise ValueError(f"datetime must be timezone-aware; got naive {dt!r}")
    return (dt - _J2000_EPOCH).total_seconds()


def _t_sec_to_dt(t_sec: float) -> datetime:
    """Convert seconds-since-J2000 back to a UTC-aware datetime."""
    return _J2000_EPOCH + timedelta(seconds=t_sec)


@dataclass(frozen=True)
class PhaseSignature:
    """A cycler's geometric fingerprint: what real-ephemeris geometry it needs.

    Attributes
    ----------
    bodies:
        Ordered tuple of body codes the cycler visits, e.g. ``("E", "M")`` for
        Aldrin or ``("E", "V", "M")`` for VEM. The first element is the
        "departure" body (a launch window is for departure from this body).
    leg_durations_s:
        Time of flight for each leg in seconds. ``len() == len(bodies) - 1``
        for an open chain; full-cycle closure isn't required here.
    vinf_target_kms:
        Per-encounter |V∞| targets, km/s. Same length as ``bodies``. Used to
        score each candidate date's Lambert solution.
    primary:
        Gravitational primary; defaults to ``"Sun"``. Future Earth-Moon /
        Jovian launch-window work will read this; the current
        :func:`find_real_windows` only handles heliocentric (``"Sun"``).
    """

    bodies: tuple[str, ...]
    leg_durations_s: tuple[float, ...]
    vinf_target_kms: tuple[float, ...]
    primary: str = "Sun"

    def __post_init__(self) -> None:
        if len(self.leg_durations_s) != len(self.bodies) - 1:
            raise ValueError(
                f"leg_durations_s length {len(self.leg_durations_s)} must equal "
                f"len(bodies) - 1 = {len(self.bodies) - 1}"
            )
        if len(self.vinf_target_kms) != len(self.bodies):
            raise ValueError(
                f"vinf_target_kms length {len(self.vinf_target_kms)} must equal "
                f"len(bodies) = {len(self.bodies)}"
            )


@dataclass(frozen=True)
class LaunchWindow:
    """A real-ephemeris date at which the cycler's required geometry exists.

    Attributes
    ----------
    departure_date:
        UTC-aware ``datetime`` of departure from the cycler's first body.
    mismatch_kms:
        Sum of per-encounter |V∞_actual - V∞_target| at this date. Lower is
        better; the search returns the ``n`` lowest-mismatch dates.
    vinf_actual_kms:
        Tuple of per-encounter actual |V∞| at this date, in the same order as
        :attr:`PhaseSignature.bodies`. Lets the caller present "you'd actually
        leave Earth at 6.4 km/s" rather than just the target.
    """

    departure_date: datetime
    mismatch_kms: float
    vinf_actual_kms: tuple[float, ...]


def phase_signature_from_catalogue_entry(entry: dict[str, Any]) -> PhaseSignature:
    """Build a :class:`PhaseSignature` directly from a catalogue YAML dict.

    Site and other consumers that haven't built a full :class:`Cycler`
    instance can call this. Skips entries that don't carry enough data —
    raises :class:`ValueError` if required fields are null/missing.
    """
    bodies_raw = entry.get("bodies")
    if not bodies_raw:
        raise ValueError(f"entry {entry.get('id')!r} has empty bodies")
    bodies = tuple(bodies_raw)

    # Schema v3 (spec §16.6.2): read trajectory.segments when migrated,
    # else the legacy flat legs[]. Segments reuse the from/to/tof_days keys.
    legs_raw = _segments_as_legs(entry)
    if len(legs_raw) < len(bodies) - 1:
        raise ValueError(
            f"entry {entry.get('id')!r} has {len(legs_raw)} legs; need {len(bodies) - 1}"
        )
    leg_durations_s: list[float] = []
    for leg in legs_raw[: len(bodies) - 1]:
        tof = leg.get("tof_days")
        if tof is None:
            raise ValueError(f"entry {entry.get('id')!r} leg {leg!r} has null tof_days")
        leg_durations_s.append(float(tof) * SECONDS_PER_DAY)

    vinf_raw = entry.get("vinf_kms_at_encounters") or []
    if len(vinf_raw) < len(bodies):
        raise ValueError(
            f"entry {entry.get('id')!r} has {len(vinf_raw)} V∞ entries; need {len(bodies)}"
        )
    vinf_target_kms: list[float] = []
    for v in vinf_raw[: len(bodies)]:
        val = v.get("vinf_kms")
        if val is None:
            raise ValueError(f"entry {entry.get('id')!r} V∞ entry {v!r} has null vinf_kms")
        vinf_target_kms.append(float(val))

    return PhaseSignature(
        bodies=bodies,
        leg_durations_s=tuple(leg_durations_s),
        vinf_target_kms=tuple(vinf_target_kms),
        primary=entry.get("primary", "Sun"),
    )


def _vinf_at_lambert(
    r1: NDArray[np.float64],
    v1_planet: NDArray[np.float64],
    r2: NDArray[np.float64],
    v2_planet: NDArray[np.float64],
    tof_s: float,
) -> tuple[float, float] | None:
    """Single-rev Lambert from r1→r2 in tof_s; return (|v∞_depart|, |v∞_arrive|) or None.

    Returns None if the Lambert solver fails for any reason (no solution,
    Newton non-convergence, near-singular geometry). The grid-scan caller
    treats None as "skip this date" rather than propagating the exception —
    real-ephemeris geometries occasionally produce Lambert problems where
    Newton stalls; one bad date shouldn't kill the whole window search.
    """
    try:
        solutions = lambert(r1, r2, tof_s, mu=MU_SUN_KM3_S2, max_revs=0)
    except (LambertConvergenceError, LambertGeometryError):
        return None
    if not solutions:
        return None
    sol = solutions[0]  # single-rev only
    vinf_depart = float(np.linalg.norm(sol.v1 - v1_planet))
    vinf_arrive = float(np.linalg.norm(sol.v2 - v2_planet))
    return vinf_depart, vinf_arrive


def _mismatch_at_date(
    signature: PhaseSignature,
    ephem: Ephemeris,
    departure_dt: datetime,
) -> tuple[float, tuple[float, ...]] | None:
    """Compute mismatch + per-encounter actual V∞ at a candidate departure date.

    Returns (mismatch_sum, (vinf_actual_per_encounter,)). Returns None if
    Lambert fails on any leg (geometry singular).
    """
    t_sec = _dt_to_t_sec(departure_dt)
    encounter_states: list[tuple[NDArray[np.float64], NDArray[np.float64]]] = []
    t_running = t_sec
    encounter_states.append(ephem.state(signature.bodies[0], t_running))
    for i, leg_tof in enumerate(signature.leg_durations_s):
        t_running += leg_tof
        encounter_states.append(ephem.state(signature.bodies[i + 1], t_running))

    # Per-encounter actual V∞: for encounter i (1..N-1) we have both an inbound
    # (arrival from leg i-1) and outbound (departure on leg i if exists); we
    # compare each to the target. For boundaries (encounter 0, encounter N-1)
    # we only have one direction.
    vinf_actual_components: list[float] = []
    mismatch = 0.0
    for i in range(len(signature.bodies)):
        if i == 0:
            # Outbound only from encounter 0
            r1, v1_planet = encounter_states[0]
            r2, v2_planet = encounter_states[1]
            result = _vinf_at_lambert(r1, v1_planet, r2, v2_planet, signature.leg_durations_s[0])
            if result is None:
                return None
            vinf_d, _ = result
            vinf_actual_components.append(vinf_d)
            mismatch += abs(vinf_d - signature.vinf_target_kms[0])
        elif i == len(signature.bodies) - 1:
            # Inbound only at terminal encounter
            r1, v1_planet = encounter_states[i - 1]
            r2, v2_planet = encounter_states[i]
            result = _vinf_at_lambert(
                r1, v1_planet, r2, v2_planet, signature.leg_durations_s[i - 1]
            )
            if result is None:
                return None
            _, vinf_a = result
            vinf_actual_components.append(vinf_a)
            mismatch += abs(vinf_a - signature.vinf_target_kms[i])
        else:
            # Both inbound (leg i-1) and outbound (leg i); report average mismatch
            r_prev, v_prev_p = encounter_states[i - 1]
            r_here, v_here_p = encounter_states[i]
            r_next, v_next_p = encounter_states[i + 1]
            res_in = _vinf_at_lambert(
                r_prev, v_prev_p, r_here, v_here_p, signature.leg_durations_s[i - 1]
            )
            res_out = _vinf_at_lambert(
                r_here, v_here_p, r_next, v_next_p, signature.leg_durations_s[i]
            )
            if res_in is None or res_out is None:
                return None
            vinf_in = res_in[1]
            vinf_out = res_out[0]
            # Cycler convention: V∞ magnitude is preserved across a ballistic
            # flyby; report the average as the "actual" V∞ at this encounter.
            vinf_actual_components.append(0.5 * (vinf_in + vinf_out))
            mismatch += abs(vinf_in - signature.vinf_target_kms[i])
            mismatch += abs(vinf_out - signature.vinf_target_kms[i])

    return mismatch, tuple(vinf_actual_components)


def find_real_windows(
    signature: PhaseSignature,
    ephem: Ephemeris,
    date_range: tuple[datetime, datetime],
    n: int = 5,
    step_days: float = 5.0,
    mismatch_cap_kms: float = 5.0,
) -> list[LaunchWindow]:
    """Grid-scan ``ephem`` across ``date_range`` for the N best-matching dates.

    Algorithm: simple grid scan at ``step_days`` resolution; track local
    minima (mismatch is below both neighbours). Sort by mismatch, drop any
    above ``mismatch_cap_kms``, return the N lowest.

    The ``"astropy"`` :class:`Ephemeris` is required for real dates; passing a
    ``"circular"`` ephemeris gives idealised geometry that's epoch-agnostic
    and not what callers asking for "launch windows" want. (Function does
    not enforce this — caller's choice.)

    Parameters
    ----------
    signature:
        The cycler's :class:`PhaseSignature`.
    ephem:
        Ephemeris backend. Use ``Ephemeris("astropy")`` for real launch dates.
    date_range:
        ``(start, end)`` UTC-aware datetimes bounding the search.
    n:
        Maximum number of windows to return (default 5).
    step_days:
        Grid resolution in days (default 5.0). Tighter = more accurate but
        slower; 5d is a good default for Earth-Mars launches (synodic
        ≈ 779d, so 5d captures < 1° of synodic phase per step).
    mismatch_cap_kms:
        Reject candidates whose mismatch exceeds this (km/s). Default 5.0;
        a real geometric match for a ballistic cycler is typically < 1 km/s.

    Returns
    -------
    Up to ``n`` :class:`LaunchWindow` instances, sorted by ascending mismatch.
    Empty list if no candidate beats ``mismatch_cap_kms``.
    """
    if signature.primary != "Sun":
        raise NotImplementedError(
            f"find_real_windows only handles heliocentric (primary='Sun'); "
            f"got primary={signature.primary!r}. Lunar/Jovian launch-window "
            "support is future work."
        )
    if ephem.model != "astropy":
        # Soft warning via docstring; programmatic enforcement would block
        # legitimate test use of the circular backend. Leave as caller's call.
        pass

    start, end = date_range
    total_days = (end - start).total_seconds() / SECONDS_PER_DAY
    n_steps = max(2, int(total_days / step_days) + 1)

    # Grid scan: compute mismatch at each step. Store (mismatch, dt, actual_vinf).
    samples: list[tuple[float, datetime, tuple[float, ...]]] = []
    for i in range(n_steps):
        dt = start + timedelta(days=i * step_days)
        if dt > end:
            break
        result = _mismatch_at_date(signature, ephem, dt)
        if result is None:
            samples.append((inf, dt, ()))
            continue
        mismatch, vinf_actual = result
        samples.append((mismatch, dt, vinf_actual))

    # Find local minima — date i is a local min if mismatch[i] < mismatch[i-1]
    # and mismatch[i] <= mismatch[i+1] (or at boundaries).
    local_mins: list[LaunchWindow] = []
    for i, (mismatch, dt, vinf_actual) in enumerate(samples):
        if mismatch > mismatch_cap_kms:
            continue
        left_ok = (i == 0) or (mismatch < samples[i - 1][0])
        right_ok = (i == len(samples) - 1) or (mismatch <= samples[i + 1][0])
        if left_ok and right_ok:
            local_mins.append(LaunchWindow(dt, mismatch, vinf_actual))

    local_mins.sort(key=lambda w: w.mismatch_kms)
    return local_mins[:n]


def leg_duration_seeds(
    bodies: tuple[str, ...],
    primary_leg_durations_s: tuple[float, ...],
    vinf_target_kms: tuple[float, ...],
    period_s: float,
    *,
    perturb_fracs: Sequence[float] = (0.0, 0.10, -0.10, 0.20, -0.20),
    min_leg_days: float = 30.0,
    max_leg_days_frac: float = 0.95,
) -> list[PhaseSignature]:
    """Generate asymmetric leg-duration perturbation seeds for window search.

    Given a primary leg-duration vector, produce one
    :class:`PhaseSignature` per entry in ``perturb_fracs`` by rescaling the
    first leg by ``(1 + f)`` and redistributing the change across the
    remaining legs so the total period is conserved (see below). This lets
    the robust epoch resolver fan a single literature signature out into a
    family of asymmetric basins rather than betting everything on the
    equispaced (symmetric) seed.

    Period conservation
    -------------------
    For a 2-leg family ``[d1, d2]`` and fraction ``f`` the new legs are
    ``[d1*(1+f), d2 - d1*f]`` — exact conservation of ``d1 + d2``. For an
    N-leg family the perturbation ``f*d0`` applied to leg 0 is spread as
    ``-f*d0/(N-1)`` across each of the remaining legs. Conservation holds
    exactly unless a clip triggers (see below).

    Clipping
    --------
    Each resulting leg is clipped to
    ``[min_leg_days * SECONDS_PER_DAY, max_leg_days_frac * period_s]`` to
    avoid Lambert-degenerate durations. Clipping a leg breaks exact period
    conservation for that seed; this is intentional and documented (the
    downstream window search tolerates a slightly off-period seed because it
    only uses the seed's per-leg ToFs to *probe* the ephemeris, not to assert
    closure).

    Parameters
    ----------
    bodies:
        Ordered body codes; passed straight through to each signature.
    primary_leg_durations_s:
        The literature/primary per-leg ToFs in seconds. ``len() ==
        len(bodies) - 1``.
    vinf_target_kms:
        Per-encounter |V∞| targets; passed straight through. ``len() ==
        len(bodies)``.
    period_s:
        Total period in seconds used as the conservation target and the
        ``max_leg_days_frac`` reference. Typically
        ``sum(primary_leg_durations_s)``.
    perturb_fracs:
        Fractional perturbations applied to leg 0. ``0.0`` reproduces the
        primary seed. Default ``(0.0, ±0.10, ±0.20)``.
    min_leg_days, max_leg_days_frac:
        Clip bounds (lower in days, upper as a fraction of ``period_s``).

    Returns
    -------
    A de-duplicated list of :class:`PhaseSignature` (duplicate leg vectors,
    e.g. produced when clipping collapses two fractions to the same bound,
    are dropped). Always contains at least the primary (``f=0``) seed.
    """
    n_legs = len(primary_leg_durations_s)
    if n_legs < 1:
        raise ValueError("primary_leg_durations_s must have >= 1 entry")
    min_leg_s = min_leg_days * SECONDS_PER_DAY
    # Upper clip protects against a leg consuming nearly the whole multi-leg
    # period. For a single-leg (open-chain) signature the leg *is* the period,
    # so floor the cap at the longest primary leg — otherwise the primary
    # (f=0) seed would itself be clipped below its literature value.
    max_leg_s = max(max_leg_days_frac * period_s, max(primary_leg_durations_s))
    d0 = primary_leg_durations_s[0]

    seen: set[tuple[int, ...]] = set()
    seeds: list[PhaseSignature] = []
    for f in perturb_fracs:
        legs = list(primary_leg_durations_s)
        delta = f * d0
        legs[0] = d0 + delta
        if n_legs > 1:
            share = delta / (n_legs - 1)
            for j in range(1, n_legs):
                legs[j] = primary_leg_durations_s[j] - share
        # Clip each leg to the Lambert-safe band.
        legs = [min(max(leg, min_leg_s), max_leg_s) for leg in legs]
        # De-duplicate on a coarse second-resolution key so float jitter and
        # clip collapses don't yield identical signatures twice.
        key = tuple(round(leg) for leg in legs)
        if key in seen:
            continue
        seen.add(key)
        seeds.append(
            PhaseSignature(
                bodies=bodies,
                leg_durations_s=tuple(legs),
                vinf_target_kms=vinf_target_kms,
            )
        )
    return seeds


def find_candidate_windows(
    signatures: Sequence[PhaseSignature],
    ephem: Ephemeris,
    date_range: tuple[datetime, datetime],
    n: int = 10,
    step_days: float = 5.0,
    mismatch_cap_kms: float = 20.0,
    dedup_window_days: float = 30.0,
) -> list[LaunchWindow]:
    """Fan multiple signatures through :func:`find_real_windows`, merge, rank.

    Each signature in ``signatures`` is scanned independently over
    ``date_range`` via :func:`find_real_windows`; the resulting windows are
    pooled, de-duplicated (windows whose departure dates are within
    ``dedup_window_days`` of each other collapse to the lower-mismatch one),
    sorted by ascending ``mismatch_kms``, and capped at ``n`` total.

    This is the multi-seed, ranked-by-mismatch replacement for the single
    signature / calendar-proximity selection that the old
    ``_resolve_real_t_start`` performed. The candidate pool is ranked
    purely by V∞ mismatch; calendar proximity is left to the caller (which
    uses ``date_range`` only to centre/bound the search).

    Parameters
    ----------
    signatures:
        Seeds to fan out (e.g. the output of :func:`leg_duration_seeds`).
    ephem:
        Ephemeris backend. Use ``Ephemeris("astropy")`` for real dates.
    date_range:
        ``(start, end)`` UTC-aware datetimes bounding every seed's scan.
    n:
        Maximum windows to return after merging.
    step_days:
        Grid resolution passed to each :func:`find_real_windows` scan.
    mismatch_cap_kms:
        Reject windows above this mismatch (km/s). Default 20.0 — looser than
        :func:`find_real_windows`'s own default so asymmetric seeds can land
        their basin; the downstream optimiser re-filters on the V∞ cap.
    dedup_window_days:
        Two windows within this many days are duplicates (lower mismatch
        kept).

    Returns
    -------
    Up to ``n`` :class:`LaunchWindow` sorted by ascending ``mismatch_kms``.
    """
    pool: list[LaunchWindow] = []
    for sig in signatures:
        pool.extend(
            find_real_windows(
                sig,
                ephem,
                date_range,
                n=n,
                step_days=step_days,
                mismatch_cap_kms=mismatch_cap_kms,
            )
        )

    # Sort by mismatch first so the O(N^2) dedup sweep always retains the
    # lowest-mismatch representative of each calendar cluster.
    pool.sort(key=lambda w: w.mismatch_kms)
    kept: list[LaunchWindow] = []
    for w in pool:
        is_dup = any(
            abs((w.departure_date - k.departure_date).total_seconds())
            < dedup_window_days * SECONDS_PER_DAY
            for k in kept
        )
        if not is_dup:
            kept.append(w)

    return kept[:n]


__all__ = [
    "LaunchWindow",
    "PhaseSignature",
    "find_candidate_windows",
    "find_real_windows",
    "leg_duration_seeds",
    "phase_signature_from_catalogue_entry",
]
