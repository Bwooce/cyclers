"""Physical-sanity flyby gate — reject patched-conic flybys with unusable bend.

Motivation (task #324, 2026-06-16). The full automated guard chain of #312
(closure + cross-check + NN + lit-corpus + ML) admitted a SILVER survivor —
the Umbriel-Oberon-Umbriel (1,1) Uranian moon-tour candidate (commit ``bebaeaf``
/ ``d0d5898``) — where the V_inf at Umbriel was substantially above what the
moon's gravity can usefully bend at the minimum safe periapsis. The candidate
was *V_inf-magnitude-continuous* across the encounter at numerical tolerance,
so every magnitude / closure guard passed; but the encounter delivered
essentially **zero** asymptote rotation, so the "flyby" was not in physical
fact a gravity assist. The agent caught it manually; the automated pipeline
didn't. This module is the missing gate.

Physics
-------
For a hyperbolic flyby at the minimum-safe periapsis radius
``r_p = r_body + min_safe_altitude_km``, the maximum ballistic deflection
angle is the patched-conic bound (Bate-Mueller-White §6.4; also documented in
:mod:`cyclerfinder.core.flyby`):

.. math::

    \\delta_\\max = 2 \\arcsin\\!\\left(
        \\frac{\\mu_\\text{body}}{\\mu_\\text{body} + r_p V_\\infty^2}
    \\right).

A useful gravity-assist must rotate the velocity asymptote by *some* finite
amount. The floor adopted here is **5°** — a judgment-call threshold (NOT a
sourced physical constant). Rationale:

  * Galileo's Earth flybys (V_inf ~6.2 km/s) achieve ~75° max bend; Cassini's
    Venus flybys (V_inf ~7 km/s) ~61°; Aldrin Mars flybys (V_inf ~5.5 km/s)
    ~32°. The flybys engineered into real cycler / tour missions sit far above
    5° — see e.g. Strange & Longuski JSR 2002 tour designs, Sims-Flanagan
    JPL planners.
  * Below ~5° the flyby contributes deflection at the level of typical
    targeting / TCM noise; calling it a "gravity assist" is essentially
    sleight-of-hand. The Umbriel case at the prompt's worst-case 2.27 km/s
    is 2.7° (below floor → reject); the actual #312 SILVER row at 0.92 km/s
    is 14.7° (above floor → admit, and the candidate proceeds to other guards).
  * A 5° floor is **stricter** than zero (catches the pathological case) but
    **looser** than any operationally interesting flyby — so admitting a flyby
    via this gate does NOT certify it as useful, only as not-pathological.

The threshold is parametric (``min_useful_bend_deg``) so callers can tighten
it. No catalogue writeback or novelty claim depends on it; gate-passing
candidates remain subject to lit-check + ML + gauntlet, per task #324 Phase 1
discipline.

Implementation
--------------
Thin wrapper around :func:`cyclerfinder.core.flyby.max_bend` (the patched-conic
formula above). No new physics; no modification of ``core/flyby.py``. The
gate is exposed at two levels:

  * :func:`flyby_is_useful` — per-encounter scalar verdict + structured record
    (``FlybyPhysicalVerdict``).
  * :func:`candidate_passes_physical_gate` — sequence-level: run the per-encounter
    check at every flyby body in a tour and reject if *any* fails.

Body lookup tries :data:`cyclerfinder.core.constants.PLANETS` first (so
heliocentric flybys at V/E/M/J/S/U/N/Me work) and then
:data:`cyclerfinder.core.satellites.SATELLITES` (so planetocentric moon-tour
flybys at Io/Europa/Titan/Umbriel/Oberon/… also work). Unknown body raises
``KeyError`` (the gate must NEVER silently pass an unknown body).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import degrees
from typing import Final

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.satellites import SATELLITES

# Default useful-bend floor (deg). Judgment threshold (see module docstring);
# stricter than zero, looser than any operationally interesting flyby.
DEFAULT_MIN_USEFUL_BEND_DEG: Final[float] = 5.0


@dataclass(frozen=True)
class FlybyPhysicalVerdict:
    """Per-encounter physical-sanity verdict.

    Attributes
    ----------
    body:
        Body code or full name as supplied to :func:`flyby_is_useful` (planet
        code like ``"E"`` or moon name like ``"Umbriel"``).
    vinf_kms:
        Hyperbolic excess speed at the flyby, km/s.
    min_safe_altitude_km:
        Periapsis-altitude floor used (effective; either caller-supplied or
        the body-default from the registry).
    periapsis_radius_km:
        Periapsis radius corresponding to the altitude floor:
        ``radius_eq_km + min_safe_altitude_km``.
    max_bend_deg:
        Maximum ballistic deflection at this V_inf and periapsis, deg.
    is_useful:
        ``True`` iff ``max_bend_deg >= min_useful_bend_deg``.
    notes:
        Free-form short note (e.g. ``"V_inf below 5 deg floor"``).
    """

    body: str
    vinf_kms: float
    min_safe_altitude_km: float
    periapsis_radius_km: float
    max_bend_deg: float
    is_useful: bool
    notes: str = ""


def _resolve_body(
    body: str,
) -> tuple[float, float, float]:
    """Return ``(mu_km3_s2, radius_eq_km, default_safe_alt_km)`` for ``body``.

    Tries :data:`PLANETS` first then :data:`SATELLITES`. Raises ``KeyError`` if
    unknown — the gate must never silently pass an unrecognised body.
    """
    if body in PLANETS:
        p = PLANETS[body]
        return p.mu_km3_s2, p.radius_eq_km, p.safe_alt_km
    if body in SATELLITES:
        s = SATELLITES[body]
        return s.mu_km3_s2, s.radius_eq_km, s.safe_alt_km
    raise KeyError(
        f"Unknown body {body!r}; not in PLANETS or SATELLITES registries. "
        "The physical-sanity gate refuses to silently admit an unknown body."
    )


def flyby_is_useful(
    body: str,
    vinf_kms: float,
    *,
    min_safe_altitude_km: float | None = None,
    min_useful_bend_deg: float = DEFAULT_MIN_USEFUL_BEND_DEG,
) -> FlybyPhysicalVerdict:
    """Check whether a flyby at ``body`` with ``V_inf=vinf_kms`` can usefully bend.

    Computes the patched-conic max-bend at the supplied (or registry-default)
    safe-altitude floor via :func:`cyclerfinder.core.flyby.max_bend`, and
    compares against ``min_useful_bend_deg``. A flyby with ``max_bend`` below
    the floor is *V_inf-magnitude-continuous-but-geometrically-vacuous* —
    formally a "flyby" in the patched-conic accounting, but with negligible
    asymptote rotation.

    Parameters
    ----------
    body:
        Planet code (V/E/M/...) or full moon name (Umbriel, Europa, ...).
        Lookup is :data:`PLANETS` then :data:`SATELLITES`; unknown → ``KeyError``.
    vinf_kms:
        Hyperbolic excess speed at the encounter, km/s. Must be non-negative.
    min_safe_altitude_km:
        Periapsis altitude floor, km. If ``None`` (default), uses the registry
        default (``safe_alt_km`` on :class:`PlanetData` / :class:`SatelliteData`).
    min_useful_bend_deg:
        Floor on ``max_bend`` for the flyby to count as useful, deg. Default
        :data:`DEFAULT_MIN_USEFUL_BEND_DEG` (5.0).

    Returns
    -------
    FlybyPhysicalVerdict
    """
    if vinf_kms < 0.0:
        raise ValueError(f"vinf_kms must be non-negative, got {vinf_kms}")
    if min_useful_bend_deg < 0.0:
        raise ValueError(f"min_useful_bend_deg must be non-negative, got {min_useful_bend_deg}")

    mu, radius_km, default_alt = _resolve_body(body)
    alt = float(default_alt if min_safe_altitude_km is None else min_safe_altitude_km)
    if alt < 0.0:
        raise ValueError(f"min_safe_altitude_km must be non-negative, got {alt}")
    rp = radius_km + alt

    bend_rad = max_bend(mu, rp, vinf_kms)
    bend_deg = degrees(bend_rad)
    useful = bend_deg >= min_useful_bend_deg

    if useful:
        notes = ""
    else:
        notes = (
            f"max_bend {bend_deg:.4f} deg below {min_useful_bend_deg:.2f} deg floor "
            f"(V_inf {vinf_kms:.4f} km/s at {body} too high for usable bend at "
            f"r_p={rp:.1f} km)"
        )

    return FlybyPhysicalVerdict(
        body=body,
        vinf_kms=float(vinf_kms),
        min_safe_altitude_km=alt,
        periapsis_radius_km=float(rp),
        max_bend_deg=float(bend_deg),
        is_useful=bool(useful),
        notes=notes,
    )


def candidate_passes_physical_gate(
    sequence: tuple[str, ...],
    vinf_kms_per_encounter: tuple[float, ...],
    *,
    min_useful_bend_deg: float = DEFAULT_MIN_USEFUL_BEND_DEG,
    per_body_min_safe_altitude_km: dict[str, float] | None = None,
) -> tuple[bool, list[FlybyPhysicalVerdict]]:
    """Run :func:`flyby_is_useful` at every encounter and reject on any fail.

    A multi-leg patched-conic tour passes the gate iff *every* encounter
    delivers at least ``min_useful_bend_deg`` of ballistic bend at the safe
    periapsis. The gate is conservative: even one unphysical flyby fails the
    candidate. The full per-encounter verdict list is returned so the caller
    can log / triage either way.

    Parameters
    ----------
    sequence:
        Tuple of body codes / moon names, one per encounter. Length must equal
        ``len(vinf_kms_per_encounter)``.
    vinf_kms_per_encounter:
        Tuple of V_inf magnitudes (km/s) at each encounter in ``sequence``.
    min_useful_bend_deg:
        Threshold passed through to :func:`flyby_is_useful` (deg).
    per_body_min_safe_altitude_km:
        Optional per-body override of the safe-altitude floor (km). Bodies
        not present fall back to the registry default. Useful for sweeping
        the gate (e.g. "what if we accept 50 km Umbriel periapsis?").

    Returns
    -------
    (passed, verdicts):
        ``passed`` is ``True`` iff all verdicts are useful. ``verdicts`` is a
        list of length ``len(sequence)`` with one :class:`FlybyPhysicalVerdict`
        per encounter (in order).
    """
    if len(sequence) != len(vinf_kms_per_encounter):
        raise ValueError(
            f"sequence (len {len(sequence)}) and vinf_kms_per_encounter "
            f"(len {len(vinf_kms_per_encounter)}) must have the same length"
        )
    overrides = per_body_min_safe_altitude_km or {}

    verdicts: list[FlybyPhysicalVerdict] = []
    for body, vinf in zip(sequence, vinf_kms_per_encounter, strict=True):
        alt = overrides.get(body)
        verdicts.append(
            flyby_is_useful(
                body,
                vinf,
                min_safe_altitude_km=alt,
                min_useful_bend_deg=min_useful_bend_deg,
            )
        )

    passed = all(v.is_useful for v in verdicts)
    return passed, verdicts
