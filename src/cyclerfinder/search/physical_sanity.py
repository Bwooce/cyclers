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

Task #818 adds the complementary PASSIVE-node gate
(:func:`passive_node_is_self_consistent`): where #324 rejects a node that is
*claimed to work but physically cannot*, #818 rejects a node that is *claimed
to do nothing but physically must do something* — see the inline section
comment above :data:`DEFAULT_MAX_PARASITIC_TURN_FRACTION`.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import degrees, inf
from typing import Final

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# Default useful-bend floor (deg). Judgment threshold (see module docstring);
# stricter than zero, looser than any operationally interesting flyby.
DEFAULT_MIN_USEFUL_BEND_DEG: Final[float] = 5.0

# --- Passive-node self-consistency gate (task #818, 2026-08-10) -------------
#
# A node whose REQUIRED turn is ~0 ("passive node": the closure treats the
# body's gravity as doing nothing) is only self-consistent if the body's OWN
# unavoidable patched-conic deflection is negligible against the trajectory's
# working turn budget. The #817 adjudication (docs/notes/
# 2026-08-10-817-passive-oberon-node-adjudication.md) did this check by hand
# for #816's Ariel-Oberon object: Oberon's parasitic deflection at its own
# Laplace SOI boundary is 1.397 deg = 33% of the 4.233 deg Ariel working turn,
# an unmodelled contradiction (the node geometry demands an encounter, the
# node dynamics demand no encounter). This gate automates that check,
# generalizing #480's geometric SOI-containment rule
# (search/tour_self_consistency.py) to the DYNAMICAL side.

# A node counts as "passive" below this required turn (deg). Same value as
# scan_816_unequal_tof_discrete_roots.py's TURN_TRIVIAL_DEG classifier
# threshold, adopted here so gate and classifier agree on what "passive" is.
PASSIVE_NODE_TURN_MAX_DEG: Final[float] = 0.05

# Max tolerable parasitic deflection as a fraction of the working turn budget.
# Judgment threshold (NOT a sourced physical constant), grounded in this
# project's own precedents rather than invented fresh:
#   * #817's adjudication used "under 2% of the working turn" as its own
#     negligibility line (its 160,000-km standoff row), and treated 6.5%
#     as still within the disqualifying discussion.
#   * Same spirit as the #324 5-deg floor above: ~2% of a working budget is
#     the level of typical targeting / TCM noise.
#   * Calibration against the two anchor cases: a Russell-Strange-2009-style
#     genuinely-negligible passive target (Enceladus at ~4 km/s vs a
#     tens-of-degrees Titan working turn) sits well below 1%; #817's Oberon
#     case is 33%. The two populations are separated by ~2 orders of
#     magnitude, so the verdict is insensitive to the exact value over
#     roughly 0.5%-10%.
DEFAULT_MAX_PARASITIC_TURN_FRACTION: Final[float] = 0.02


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


# ---------------------------------------------------------------------------
# Passive-node self-consistency gate (task #818)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PassiveNodeVerdict:
    """Self-consistency verdict for a near-zero-required-turn ("passive") node.

    Attributes
    ----------
    body:
        Body at the passive node (planet code or full moon name).
    vinf_kms:
        Hyperbolic excess speed at the passive node, km/s.
    working_turn_budget_deg:
        Sum of the required turns at the trajectory's WORKING (non-passive)
        nodes, deg — the denominator the parasitic deflection is ranked
        against.
    soi_radius_km:
        Laplace sphere-of-influence radius ``a * (mu / mu_primary)**(2/5)``,
        km. The patched-conic encounter boundary; the verdict is taken here.
    hill_radius_km:
        Hill radius ``a * (mu / (3 * mu_primary))**(1/3)``, km (same formula
        as :func:`cyclerfinder.search.tour_self_consistency.soi_km`).
        Informational: the more lenient boundary; reported so a caller can see
        both bounds, as #817's own table did.
    parasitic_deflection_deg:
        ``max_bend`` at the Laplace-SOI radius, deg — the MINIMUM deflection
        the body imparts on ANY pass that counts as an encounter (any real
        periapsis inside the SOI bends MORE). Verdict-bearing.
    parasitic_deflection_hill_deg:
        Same lower bound taken at the Hill radius instead, deg. Informational.
    parasitic_fraction:
        ``parasitic_deflection_deg / working_turn_budget_deg`` (``inf`` when
        the budget is non-positive).
    is_self_consistent:
        ``True`` iff ``parasitic_fraction <= max_parasitic_fraction``.
    notes:
        Free-form short note on rejection.
    """

    body: str
    vinf_kms: float
    working_turn_budget_deg: float
    soi_radius_km: float
    hill_radius_km: float
    parasitic_deflection_deg: float
    parasitic_deflection_hill_deg: float
    parasitic_fraction: float
    is_self_consistent: bool
    notes: str = ""


def _resolve_orbit(body: str) -> tuple[float, float, float]:
    """Return ``(mu_km3_s2, sma_km, mu_primary_km3_s2)`` for ``body``.

    Planets orbit the Sun (``sma_au * AU_KM``, :data:`MU_SUN_KM3_S2`);
    satellites orbit their registry primary (:data:`PRIMARIES`). Raises
    ``KeyError`` if unknown — same never-silently-pass rule as
    :func:`_resolve_body`.
    """
    if body in PLANETS:
        p = PLANETS[body]
        return p.mu_km3_s2, p.sma_au * AU_KM, MU_SUN_KM3_S2
    if body in SATELLITES:
        s = SATELLITES[body]
        return s.mu_km3_s2, s.sma_km, PRIMARIES[s.primary]
    raise KeyError(
        f"Unknown body {body!r}; not in PLANETS or SATELLITES registries. "
        "The passive-node gate refuses to silently admit an unknown body."
    )


def laplace_soi_km(body: str) -> float:
    """Laplace sphere-of-influence radius (km): ``a * (mu / mu_primary)**(2/5)``.

    The standard patched-conic encounter boundary (Bate-Mueller-White §7.4).
    NOTE this is NOT the same as :func:`cyclerfinder.search.
    tour_self_consistency.soi_km`, which computes the (larger) HILL radius
    ``a * (mu / (3 * mu_primary))**(1/3)`` — for Oberon: 9,678 km Laplace vs
    13,288 km Hill. #817's verdict-bearing row used the Laplace SOI.
    """
    mu, sma_km, mu_primary = _resolve_orbit(body)
    return float(sma_km * (mu / mu_primary) ** 0.4)


def hill_radius_km(body: str) -> float:
    """Hill radius (km): ``a * (mu / (3 * mu_primary))**(1/3)``.

    Same formula as :func:`cyclerfinder.search.tour_self_consistency.soi_km`
    (guarded by a parity test), extended to planets about the Sun.
    """
    mu, sma_km, mu_primary = _resolve_orbit(body)
    return float(sma_km * (mu / (3.0 * mu_primary)) ** (1.0 / 3.0))


def passive_node_is_self_consistent(
    body: str,
    vinf_kms: float,
    working_turn_budget_deg: float,
    *,
    max_parasitic_fraction: float = DEFAULT_MAX_PARASITIC_TURN_FRACTION,
) -> PassiveNodeVerdict:
    """Check whether a ~zero-required-turn node at ``body`` is self-consistent.

    A closure that routes the trajectory THROUGH a body's position while
    requiring ~0 deg of turn there asserts a contradiction: in patched-conic
    terms an encounter is a pass inside the body's SOI, and ANY pass inside
    the SOI (periapsis ``r_p <= r_SOI``) deflects the asymptote by at least

    .. math::

        \\delta_\\text{parasitic} = 2 \\arcsin\\!\\left(
            \\frac{1}{1 + r_\\text{SOI} V_\\infty^2 / \\mu}\\right),

    the same law as :func:`cyclerfinder.core.flyby.max_bend` (reused verbatim
    — #817 verified this reproduces the stored ``max_bend_deg_per_encounter``
    values bit-for-bit) evaluated at the SOI boundary instead of the safe
    periapsis. That unmodelled deflection is *parasitic*: the closure was
    solved without it. The node is self-consistent only if the parasitic
    lower bound is negligible against the trajectory's working turn budget:

        ``parasitic_deflection_deg <= max_parasitic_fraction *
        working_turn_budget_deg``

    Threshold rationale: see :data:`DEFAULT_MAX_PARASITIC_TURN_FRACTION`.
    Whether a node counts as "passive" in the first place is the caller's
    classification (:data:`PASSIVE_NODE_TURN_MAX_DEG` is the project
    convention, matching #816's ``TURN_TRIVIAL_DEG``).

    Parameters
    ----------
    body:
        Planet code (V/E/M/...) or full moon name at the passive node.
        Unknown body raises ``KeyError``.
    vinf_kms:
        Hyperbolic excess speed at the passive node, km/s. Non-negative.
    working_turn_budget_deg:
        Total required turn at the trajectory's working node(s), deg. A
        non-positive budget can never absorb a parasitic deflection, so the
        verdict is automatically inconsistent (fraction ``inf``).
    max_parasitic_fraction:
        Tolerable ``parasitic / budget`` ratio. Must be positive. Default
        :data:`DEFAULT_MAX_PARASITIC_TURN_FRACTION` (0.02).

    Returns
    -------
    PassiveNodeVerdict
    """
    if vinf_kms < 0.0:
        raise ValueError(f"vinf_kms must be non-negative, got {vinf_kms}")
    if max_parasitic_fraction <= 0.0:
        raise ValueError(f"max_parasitic_fraction must be positive, got {max_parasitic_fraction}")

    mu, _sma_km, _mu_primary = _resolve_orbit(body)
    r_soi = laplace_soi_km(body)
    r_hill = hill_radius_km(body)
    parasitic_deg = degrees(max_bend(mu, r_soi, vinf_kms))
    parasitic_hill_deg = degrees(max_bend(mu, r_hill, vinf_kms))

    fraction = parasitic_deg / working_turn_budget_deg if working_turn_budget_deg > 0.0 else inf
    consistent = fraction <= max_parasitic_fraction

    if consistent:
        notes = ""
    elif working_turn_budget_deg <= 0.0:
        notes = (
            f"non-positive working turn budget {working_turn_budget_deg:.4f} deg cannot "
            f"absorb the {parasitic_deg:.4f} deg parasitic deflection at {body}"
        )
    else:
        notes = (
            f"parasitic deflection {parasitic_deg:.4f} deg at {body}'s SOI "
            f"(r={r_soi:.0f} km, V_inf {vinf_kms:.4f} km/s) is "
            f"{100.0 * fraction:.1f}% of the {working_turn_budget_deg:.4f} deg working "
            f"turn budget, above the {100.0 * max_parasitic_fraction:.1f}% floor — the "
            f"closure models this node as doing nothing but the encounter it demands "
            f"cannot do nothing"
        )

    return PassiveNodeVerdict(
        body=body,
        vinf_kms=float(vinf_kms),
        working_turn_budget_deg=float(working_turn_budget_deg),
        soi_radius_km=float(r_soi),
        hill_radius_km=float(r_hill),
        parasitic_deflection_deg=float(parasitic_deg),
        parasitic_deflection_hill_deg=float(parasitic_hill_deg),
        parasitic_fraction=float(fraction),
        is_self_consistent=bool(consistent),
        notes=notes,
    )
