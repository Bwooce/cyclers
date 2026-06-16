"""Tisserand-Poincaré graph enumerator for epoch-windowed MGA chains.

Phase 2 of task #289 (Track-A epoch-locked trajectory substrate). Phase 1
(:mod:`cyclerfinder.genome.epoch_aware_genome`) delivered the
``EpochLockedTrajectory`` value type and the per-leg single-rev prograde
Lambert + ballistic-continuity + Kepler-cross-check closure driver
:func:`close_epoch_locked`. The Tito 2018 free-return closes through it at
0.063 km/s worst residual and 0.0105 km/s flyby continuity.

This module *proposes* candidate chains; Phase 1 *closes* them. The
enumerator's responsibility is pure geometry — it never solves a Lambert
itself, and it never writes the catalogue.

Methodology
-----------
The Tisserand parameter ``T_p`` at a flyby body conserves the spacecraft
hyperbolic excess ``V_inf`` across that flyby (Tisserand 1896). A
constant-:math:`V_\\infty` flyby at body ``p`` lives on a constant-:math:`T_p`
curve in ``(a, e)`` space; the Strange-Russell 2007 "V_inf globe" lifts this
to a (V_inf, declination) sphere where geometric continuity between two
flybys is exactly the Tisserand contour intersection (Strange & Longuski JSR
39(1):9-16, 2002; Campagnola & Russell JGCD 33(2):476-486, 2010). The
Petropoulos-Longuski pump-tour combinatorics map that geometric pre-screen
onto a directed-graph BFS over (planet, V_inf-bin) nodes.

Concretely:

  * **Nodes** are ``(planet_code, vinf_kms_bin)`` pairs. Each node names a
    constant-:math:`V_\\infty` family at one body — the family of spacecraft
    orbits that can encounter ``planet_code`` at ``vinf_kms_bin`` ± a small
    binning tolerance.
  * **Edges** are V_inf-preserving geometric continuations between two
    nodes. For two *different* planets the edge is admissible iff the
    inclined Tisserand predicate
    :func:`cyclerfinder.search.tisserand.linkable_3d` returns True at this
    common V_inf — i.e. there exists a heliocentric ``(a, e, i_sc)`` orbit
    reachable from both bodies at this :math:`V_\\infty`. For the *same*
    planet the edge is admissible if a resonance band exists whose
    spacecraft period puts the TOF inside ``tof_box_days_per_leg``.
  * **TOF per leg** is estimated geometrically (Hohmann half-period for
    hetero hops; resonance ratio for same-planet hops) and constrained to
    the caller's TOF box. The estimate is a *proposal*; the real
    flight time at the actual launch epoch falls out of the Phase-1 Lambert
    closure.
  * **Chain score** is the sum of (number of legs) + a small chain-length
    penalty + a tiny per-V_inf-spread term — it is a heuristic for ordering
    candidates so the most-likely-to-close come first.

The enumerator is deliberately pessimistic: it over-prunes (anything the
3-D Tisserand predicate rejects at the requested V_inf grid is dropped) and
under-claims (every surviving candidate must still pass the Phase-1
closure to be even SILVER-grade). It is *search* infrastructure, not a
catalogue.

Galileo VEEGA integration test
------------------------------
The 1989-October Galileo trajectory (V, E, E, J) is the canonical
``mga_tour`` archetype (Diehl-Belbruno-Roberts AAS 1986;
KNOWN_CORPUS at corpus rev ``568d8a4``). The integration test sweeps
``planet_set=("V","E","J")``, ``launch_window=(1989-10-01, 1989-11-15)``,
``max_legs=4`` and verifies that a (V, E, E, J) candidate surfaces with a
launch epoch within ±4 weeks of the published 1989-10-18T16:53:40 UTC
liftoff. If the closure is not within tolerance the test reports the
calibration gap honestly — it does not silently pass.

Scope (Phase 2)
---------------
* Ballistic only (no DSMs); Phase 3+ wraps :mod:`dsm_leg` around the
  surviving candidates.
* Single-rev prograde Lambert closure (inherited from Phase 1).
* No catalogue writeback — surviving candidates go to a discovery JSONL.
* No "novel" framing — :mod:`literature_check` (#272) is mandatory before
  any such claim.

References
----------
* Tisserand, F., *Traité de Mécanique Céleste*, vol. 4, 1896.
* Strange, N. J. & Longuski, J. M., "Graphical Method for Gravity-Assist
  Trajectory Design", *J. Spacecraft and Rockets*, 39(1):9-16, 2002.
* Campagnola, S. & Russell, R. P., "Endgame Problem Part 2: Multibody
  Technique and the Tisserand-Poincaré Graph", *JGCD*, 33(2):476-486, 2010.
* Petropoulos, A. E. & Longuski, J. M., "Shape-Based Algorithm for the
  Automated Design of Low-Thrust, Gravity Assist Trajectories",
  AAS/AIAA Astrodynamics Specialist Conf., 2000.
* Diehl, R. E., Belbruno, E. A., Roberts, R. M., "Galileo Mission Plan
  for the (Venus, Earth, Earth, Jupiter) Trajectory", AAS 1986 (sourced).
* See ``docs/notes/2026-06-16-298-289-phase2-tisserand-enumerator.md`` for
  the phase plan and the Galileo VEEGA reproduction probe.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SAFE_PERIHELION_KM,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.genome.epoch_aware_genome import (
    EpochLockedClosure,
    EpochLockedTrajectory,
    close_epoch_locked,
)
from cyclerfinder.search.tisserand import linkable_3d, vinf_to_tisserand

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris


# --------------------------------------------------------------------------- #
# Public dataclass
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MGAChainCandidate:
    """A geometrically-feasible MGA chain proposal — pre-Lambert.

    Produced by :func:`find_mga_chains`. The fields name a sequence of body
    encounters, a per-encounter ``V_inf`` (km/s), a per-leg TOF (days), and
    a nominal launch epoch (ISO-8601 UTC); the chain score is a heuristic
    ordering hint (lower is better) and the Tisserand parameter is the
    V_inf-invariant of the chain at the first body (sanity hook for
    downstream cross-checks).

    Candidates are *proposals* — they still need to pass
    :func:`validate_chain_candidate` (which wraps them in an
    :class:`EpochLockedTrajectory` and calls Phase 1's
    :func:`close_epoch_locked`) to be even SILVER-grade.
    """

    sequence: tuple[str, ...]
    """One-letter body codes in encounter order (e.g. ``("V","E","E","J")``)."""

    vinf_tuple_kms: tuple[float, ...]
    """``|V_inf|`` (km/s) at each encounter, one per body. Length matches
    ``sequence``. Same-body hops keep V_inf bit-identical (Tisserand
    conservation); hetero hops carry the per-body V_inf bin (which may
    differ slightly across the chain — small ΔV_inf is tolerated by the
    binning width). Used as the EXPECTED side downstream."""

    leg_tofs_days: tuple[float, ...]
    """Geometric TOF estimate per leg (days). Length = ``len(sequence) - 1``.
    Inside the caller's ``tof_box_days_per_leg``. The real flight time at
    the actual launch epoch is fixed by Phase 1's Lambert solve, but this
    proposal is the *starting point* for that solve and a ranking signal."""

    launch_epoch_utc: str
    """ISO-8601 UTC of the first body encounter."""

    tisserand_parameter: float
    """``T_p`` at the FIRST body of the chain (chain-invariant under
    Tisserand conservation). Sanity hook — every downstream encounter's
    implied ``T_p`` should equal this within the V_inf binning width."""

    chain_score: float
    """Heuristic chain quality (lower is better). Tracks (# legs) +
    chain-V_inf-spread + safe-altitude-margin proxies. Used to order the
    enumerator output; downstream callers may cap by ``chain_score_threshold``
    in :func:`find_mga_chains`."""


# --------------------------------------------------------------------------- #
# Internal helpers
# --------------------------------------------------------------------------- #


def _planet_period_days(body: str) -> float:
    """Heliocentric orbital period of ``body`` in days, from Kepler's third law.

    Read from :data:`PLANETS` and the canonical ``MU_SUN_KM3_S2``. Matches the
    ``mean_motion_deg_day`` derivation in :mod:`cyclerfinder.core.constants`.
    """
    a_km = PLANETS[body].sma_au * AU_KM
    period_s = 2.0 * math.pi * math.sqrt(a_km**3 / MU_SUN_KM3_S2)
    return period_s / SECONDS_PER_DAY


def _hohmann_tof_days(body_a: str, body_b: str) -> float:
    """Hohmann half-period (days) between the orbits of ``body_a`` and ``body_b``.

    The transfer ellipse has ``a_t = (a_a + a_b) / 2`` and TOF = half its
    period. This is the *geometric* TOF estimate the enumerator gives the
    downstream Lambert as a seed; the actual flight time at the real launch
    epoch may differ by the geometric-vs-real-ephemeris phasing.

    For same-body hops (``body_a == body_b``) the Hohmann TOF degenerates to
    the body's own orbital period; callers should branch on body equality and
    use :func:`_resonant_tof_days` instead.
    """
    a_a_km = PLANETS[body_a].sma_au * AU_KM
    a_b_km = PLANETS[body_b].sma_au * AU_KM
    a_t_km = 0.5 * (a_a_km + a_b_km)
    period_s = 2.0 * math.pi * math.sqrt(a_t_km**3 / MU_SUN_KM3_S2)
    return 0.5 * period_s / SECONDS_PER_DAY


def _resonant_tof_days_candidates(
    body: str,
    vinf_kms: float,
    tof_box_days: tuple[float, float],
) -> list[tuple[float, int, int]]:
    """Resonance bands for a same-body (``body -> body``) hop at this V_inf.

    A V_inf-resonant return places the spacecraft on a heliocentric ellipse
    whose period satisfies ``M * T_sc = N * T_p`` for integers ``(M, N)``;
    the spacecraft makes ``M`` revolutions while the body makes ``N``. The
    spacecraft semi-major axis is set by the Tisserand identity at this V_inf
    (one of the two branches of the contour); the TOF of the leg is
    ``M * T_sc`` days. We enumerate ``(M, N)`` with ``M, N <= 5`` and return
    those whose TOF lies inside ``tof_box_days``.

    Returns a list of ``(tof_days, M, N)`` triples sorted by ``tof_days``
    ascending. Empty list if the body cannot resonate inside the TOF box at
    this V_inf (e.g. inner-planet hop with TOF box too narrow).

    The Tisserand identity in coplanar form is
    ``T_p = a_p/a + 2 sqrt(a/a_p * (1 - e^2))``; with ``T_p`` known from V_inf,
    each ``e`` gives two ``a`` solutions. The resonance fixes ``a`` directly
    via ``a / a_p = (M / N)**(2/3)`` (Kepler's third law). We accept that ``a``
    only if it produces a real eccentricity ``e in [0, 1)``; otherwise that
    resonance is geometrically infeasible at this V_inf.
    """
    t_p = vinf_to_tisserand(body, vinf_kms)
    period_p_days = _planet_period_days(body)
    out: list[tuple[float, int, int]] = []
    for m_rev in range(1, 6):
        for n_p in range(1, 6):
            a_ratio = (float(n_p) / float(m_rev)) ** (2.0 / 3.0)
            # Tisserand back-solve for e^2:
            # T_p = 1/a_ratio + 2 sqrt(a_ratio * (1-e^2))
            # -> sqrt(a_ratio * (1-e^2)) = (T_p - 1/a_ratio) / 2
            rhs = 0.5 * (t_p - 1.0 / a_ratio)
            if rhs <= 0.0:
                continue
            one_minus_e2 = (rhs * rhs) / a_ratio
            if not (0.0 < one_minus_e2 <= 1.0):
                continue
            # Spacecraft period in days = m_rev * planet revs / n_p? No:
            # n_p * T_p = m_rev * T_sc, so T_sc = (n_p / m_rev) * T_p.
            # Total TOF spent over m_rev sc revolutions:
            t_sc_days = (float(n_p) / float(m_rev)) * period_p_days
            tof_days = m_rev * t_sc_days
            if tof_box_days[0] <= tof_days <= tof_box_days[1]:
                out.append((tof_days, m_rev, n_p))
    out.sort(key=lambda triple: triple[0])
    return out


def delta_vinf_max_kms(
    body: str,
    vinf_in_kms: float,
    *,
    periapsis_altitude_km: float | None = None,
) -> float:
    """Strange-Longuski 2002 §12 gravity-assist pump envelope at body ``body``.

    The maximum heliocentric ``|V_inf|`` change a single flyby at ``body``
    can impart at the **next** body, given an incoming hyperbolic excess
    ``vinf_in_kms``. This is the pump-tour increment (Strange-Longuski 2002
    §12; pump-tour combinatorics in Petropoulos-Longuski 2000).

    Sourced golden. For Earth at V_inf = 9 km/s the published figure (this
    repository's KNOWN_CORPUS anchor at corpus rev ``568d8a4``,
    Strange-Longuski 2002 §12) is ~2.7 km/s per flyby. See
    ``docs/notes/2026-06-16-298-289-phase2-tisserand-enumerator.md`` §Phase
    3 hand-off.

    Derivation. At a flyby, the body-frame :math:`|v_\\infty|` is
    conserved and the :math:`v_\\infty` vector rotates by at most
    :math:`\\delta_\\text{max}` (cone half-angle at safe periapsis,
    :func:`max_bend`). The change in heliocentric speed (scalar) that drives
    the change in V_inf at the **next** encounter follows the
    Strange-Longuski 2002 §12 leveraging identity:

    .. math::

        \\Delta V_\\infty^{\\text{next}}_\\text{max}
        \\;\\approx\\; V_\\infty^{\\text{in}}\\;
                       \\bigl(1 - \\cos\\delta_\\text{max}\\bigr).

    This is the heliocentric speed-magnitude change picked up by the
    spacecraft when the v_inf asymptote is rotated through the maximum cone
    half-angle (the perpendicular-to-V_p case). For Earth at V_inf = 9 km/s
    and safe periapsis (~300 km altitude) this evaluates to ~3.24 km/s; the
    sourced 2.7 km/s figure is recovered with a realised-vs-geometric factor
    of ~0.83 (the rotation usually doesn't reach the geometric maximum
    because it must point at the next body, not just anywhere on the V_inf
    globe). The function returns the GEOMETRIC bound; multi-shell callers
    that want the realised-shift bound pass ``pump_envelope_factor=0.83``
    to :func:`find_mga_chains`.

    Discrepancy honesty. The textbook closed-form
    :math:`2 V_p \\sin(\\delta/2)` (where :math:`V_p` is the body's
    heliocentric circular speed) gives ~25 km/s for Earth at V_inf = 9 —
    that is the maximum INSTANTANEOUS heliocentric velocity-vector change
    at the flyby, NOT the maximum change in **V_inf at the next body**.
    The two quantities are distinct: the instantaneous heliocentric ΔV is
    bounded by the body's circular speed; the per-encounter V_inf change at
    the NEXT body is bounded by the energy / Tisserand-shift envelope —
    which is the smaller :math:`V_\\infty(1 - \\cos\\delta)` expression
    used here. Source: Strange-Longuski 2002 §12.

    Parameters
    ----------
    body:
        One-letter body code from :data:`PLANETS`.
    vinf_in_kms:
        Incoming hyperbolic excess speed, km/s. Must be non-negative.
    periapsis_altitude_km:
        Optional override on the flyby altitude. ``None`` uses
        :data:`SAFE_PERIHELION_KM` (the engineering safe periapsis radius).

    Returns
    -------
    float
        Maximum per-flyby V_inf change at ``body`` for an incoming
        :math:`|v_\\infty|` = ``vinf_in_kms``, km/s. Non-negative.
    """
    if vinf_in_kms < 0.0:
        raise ValueError(f"vinf_in_kms must be non-negative, got {vinf_in_kms}")
    if body not in PLANETS:
        raise ValueError(f"unknown body code {body!r}; valid: {tuple(PLANETS)!r}")
    if vinf_in_kms == 0.0:
        return 0.0
    planet = PLANETS[body]
    mu = planet.mu_km3_s2
    if periapsis_altitude_km is None:
        rp = SAFE_PERIHELION_KM[body]
    else:
        rp = planet.radius_eq_km + float(periapsis_altitude_km)
    # Deflection cone half-angle at safe periapsis for this V_inf.
    delta_max = max_bend(mu, rp, vinf_in_kms)
    # Strange-Longuski 2002 §12 leveraging identity: per-flyby change in
    # V_inf at the NEXT body, evaluated at the maximum-effect geometry.
    return vinf_in_kms * (1.0 - math.cos(delta_max))


def _snap_to_grid(v: float, grid: tuple[float, ...]) -> float:
    """Snap a continuous V_inf value to its nearest grid point."""
    return min(grid, key=lambda g: abs(g - v))


def _edge_admissible(
    body_a: str,
    body_b: str,
    vinf_kms: float,
    *,
    a_range_au: tuple[float, float],
) -> bool:
    """Inclined-Tisserand admissibility of a (body_a -> body_b) hop at V_inf.

    Same-body edges are always Tisserand-trivial (the spacecraft simply
    revisits the same body on a resonant return). Hetero-body edges defer to
    :func:`linkable_3d` — the 3-D extension of the M2 ``linkable`` predicate
    that allows non-zero spacecraft inclination.

    Returns False on any predicate exception so the enumerator never crashes
    on a degenerate (V_inf, body, body) corner.
    """
    if body_a == body_b:
        return True
    try:
        return bool(linkable_3d(body_a, body_b, vinf_kms, a_range_au=a_range_au))
    except Exception:
        return False


def _adjacency(
    planet_set: tuple[str, ...],
    vinf_grid_kms: tuple[float, ...],
    *,
    a_range_au: tuple[float, float],
) -> dict[tuple[str, float], list[tuple[str, float]]]:
    """Build the Tisserand-Poincaré (body, V_inf) adjacency dict.

    For every node ``(body_a, vinf_a)`` in the cartesian product
    ``planet_set x vinf_grid_kms``, compute the set of successor nodes
    ``(body_b, vinf_a)`` reachable in one Tisserand-admissible flyby. We
    enforce V_inf conservation across the flyby so the successor's V_inf
    bin is identical to the predecessor's.

    Returns a dict mapping each node to its successor list. Empty successor
    lists are kept (for diagnostic enumeration).
    """
    adj: dict[tuple[str, float], list[tuple[str, float]]] = {}
    for body_a in planet_set:
        for vinf in vinf_grid_kms:
            node = (body_a, float(vinf))
            successors: list[tuple[str, float]] = []
            for body_b in planet_set:
                if _edge_admissible(body_a, body_b, vinf, a_range_au=a_range_au):
                    successors.append((body_b, float(vinf)))
            adj[node] = successors
    return adj


def _multi_shell_adjacency(
    planet_set: tuple[str, ...],
    vinf_grid_kms: tuple[float, ...],
    *,
    a_range_au: tuple[float, float],
    pump_envelope_factor: float = 1.0,
) -> dict[tuple[str, float], list[tuple[str, float]]]:
    """Multi-shell Tisserand-Poincaré adjacency: allow V_inf shifts within the
    Strange-Longuski 2002 §12 pump envelope at every hetero-body flyby.

    Each node ``(body_a, vinf_a)`` may transition to ``(body_b, vinf_b)``
    where ``|vinf_b - vinf_a| <= delta_vinf_max_kms(body_a, vinf_a) *
    pump_envelope_factor``. Same-body hops stay V_inf-conservative (the
    Tisserand identity is bit-exact at the same flyby body).

    ``pump_envelope_factor`` defaults to 1.0 (the full geometric envelope);
    pass <1.0 to be more conservative (the realised shift is usually a
    fraction of the geometric maximum because the V_inf rotation must also
    point at the next body, not just anywhere on the V_inf globe).
    """
    if pump_envelope_factor <= 0.0:
        raise ValueError(
            f"pump_envelope_factor must be positive, got {pump_envelope_factor}",
        )
    adj: dict[tuple[str, float], list[tuple[str, float]]] = {}
    for body_a in planet_set:
        for vinf_a in vinf_grid_kms:
            v_a = float(vinf_a)
            node = (body_a, v_a)
            successors: list[tuple[str, float]] = []
            # The pump envelope at this body and V_inf.
            envelope_kms = delta_vinf_max_kms(body_a, v_a) * pump_envelope_factor
            for body_b in planet_set:
                if body_a == body_b:
                    # Same-body resonant return: V_inf is Tisserand-conserved.
                    if _edge_admissible(body_a, body_b, v_a, a_range_au=a_range_au):
                        successors.append((body_b, v_a))
                else:
                    # Hetero body: try every neighbouring V_inf bin inside the
                    # envelope. The next-body V_inf is constrained both by the
                    # pump envelope (geometric reachability) AND by the 3-D
                    # Tisserand predicate at the *new* V_inf.
                    for vinf_b in vinf_grid_kms:
                        v_b = float(vinf_b)
                        if abs(v_b - v_a) > envelope_kms:
                            continue
                        # The Tisserand admissibility is at the *outgoing*
                        # V_inf (the V_inf the spacecraft carries on the leg
                        # to body_b). For a hetero hop the outgoing V_inf at
                        # body_a equals the body_a frame's V_inf at the moment
                        # of flyby; we use v_a (V_inf is body-frame conserved
                        # AT the flyby), but the leg's heliocentric geometry
                        # only matters for the predicate at v_b.
                        if _edge_admissible(body_a, body_b, v_b, a_range_au=a_range_au):
                            successors.append((body_b, v_b))
            adj[node] = successors
    return adj


def _enumerate_multi_shell_sequences(
    adj: dict[tuple[str, float], list[tuple[str, float]]],
    planet_set: tuple[str, ...],
    vinf_grid_kms: tuple[float, ...],
    *,
    max_legs: int,
    start_body_filter: Iterable[str] | None = None,
) -> Iterator[tuple[tuple[str, ...], tuple[float, ...]]]:
    """BFS over the multi-shell adjacency; V_inf may shift across hetero flybys.

    Generalises :func:`_enumerate_sequences` (single-shell) by letting V_inf
    take per-node values rather than carrying one shell value across the
    whole chain. The Galileo pattern (V_inf ≈ 4 → 9 → 9 → 5.6 km/s) is now
    representable.
    """
    start_filter: set[str] | None = None if start_body_filter is None else set(start_body_filter)
    for vinf in vinf_grid_kms:
        v = float(vinf)
        for start_body in planet_set:
            if start_filter is not None and start_body not in start_filter:
                continue
            frontier: list[list[tuple[str, float]]] = [[(start_body, v)]]
            while frontier:
                next_frontier: list[list[tuple[str, float]]] = []
                for path in frontier:
                    last = path[-1]
                    if len(path) >= 2:
                        seq = tuple(p[0] for p in path)
                        vinf_tuple = tuple(p[1] for p in path)
                        yield seq, vinf_tuple
                    if len(path) - 1 >= max_legs:
                        continue
                    for succ in adj.get(last, []):
                        next_frontier.append([*path, succ])
                frontier = next_frontier


def _enumerate_sequences(
    adj: dict[tuple[str, float], list[tuple[str, float]]],
    planet_set: tuple[str, ...],
    vinf_grid_kms: tuple[float, ...],
    *,
    max_legs: int,
    start_body_filter: Iterable[str] | None = None,
) -> Iterator[tuple[tuple[str, ...], tuple[float, ...]]]:
    """BFS over the Tisserand adjacency, yielding (sequence, vinf_tuple) pairs.

    The traversal builds sequences of length 2..(max_legs + 1) bodies (i.e.
    1..max_legs legs). V_inf is preserved across edges (each flyby conserves
    it), so the V_inf tuple for a sequence is just the per-encounter V_inf
    bin — for ballistic chains under the Tisserand model, every encounter on
    the chain shares the same V_inf bin (the chain is one "row" of the
    Strange-Russell V_inf globe).

    For Phase 2 we cap chain length at ``max_legs`` and emit every prefix of
    length >= 2. The downstream Lambert closure deduplicates: two sequences
    that share a prefix close to the same trajectory.
    """
    start_filter: set[str] | None = None if start_body_filter is None else set(start_body_filter)
    for vinf in vinf_grid_kms:
        v = float(vinf)
        for start_body in planet_set:
            if start_filter is not None and start_body not in start_filter:
                continue
            # BFS frontier: each entry is a path of (body, vinf) tuples.
            frontier: list[list[tuple[str, float]]] = [[(start_body, v)]]
            while frontier:
                next_frontier: list[list[tuple[str, float]]] = []
                for path in frontier:
                    last = path[-1]
                    if len(path) >= 2:
                        seq = tuple(p[0] for p in path)
                        vinf_tuple = tuple(p[1] for p in path)
                        yield seq, vinf_tuple
                    if len(path) - 1 >= max_legs:
                        continue
                    for succ in adj.get(last, []):
                        next_frontier.append([*path, succ])
                frontier = next_frontier


def _leg_tof_proposal(
    body_a: str,
    body_b: str,
    vinf_kms: float,
    tof_box_days: tuple[float, float],
) -> float | None:
    """Pick a single proposal TOF (days) for the (a -> b) leg at this V_inf.

    * Same-body (``body_a == body_b``): the lowest-TOF resonance whose
      period sits inside ``tof_box_days``. Returns None if none exists at
      this V_inf.
    * Hetero (``body_a != body_b``): the Hohmann half-period if inside
      ``tof_box_days``, otherwise the box edge nearest the Hohmann value.
      Hetero legs always have *some* proposal TOF (the Lambert downstream
      adapts), so we never return None for hetero in Phase 2.
    """
    if body_a == body_b:
        cands = _resonant_tof_days_candidates(body_a, vinf_kms, tof_box_days)
        if not cands:
            return None
        return cands[0][0]
    hoh = _hohmann_tof_days(body_a, body_b)
    if tof_box_days[0] <= hoh <= tof_box_days[1]:
        return hoh
    # Clamp into the box. The Lambert solver will re-fit; the proposal is a
    # ranking-and-validity-check signal, not the final answer.
    if hoh < tof_box_days[0]:
        return tof_box_days[0]
    return tof_box_days[1]


def _chain_score(
    sequence: tuple[str, ...],
    vinf_tuple: tuple[float, ...],
    leg_tofs_days: tuple[float, ...],
) -> float:
    """Heuristic score: lower is better.

    Combines:
      * number of legs (we prefer shorter chains, all else equal),
      * V_inf spread (we prefer chains where every encounter sits on the same
        V_inf level set; spread > 0 indicates the binning approximated a
        slightly-mismatched chain),
      * total TOF (per-encounter normalised; shorter chains rank higher).

    Pure ranking — does NOT gate admission. Use ``chain_score_threshold`` on
    :func:`find_mga_chains` to cap.
    """
    n_legs = len(sequence) - 1
    vinf_spread = max(vinf_tuple) - min(vinf_tuple) if vinf_tuple else 0.0
    total_tof = sum(leg_tofs_days)
    return float(n_legs + vinf_spread + 1.0e-4 * total_tof)


def _add_days_utc(utc_iso: str, days: float) -> str:
    """Add ``days`` to an ISO-UTC string. Astropy-backed, matches Phase 1."""
    from astropy.time import Time, TimeDelta

    t = Time(utc_iso, scale="utc") + TimeDelta(days * 86400.0, format="sec")
    return str(t.utc.isot) + "Z"


# --------------------------------------------------------------------------- #
# Public surface
# --------------------------------------------------------------------------- #


def find_mga_chains(
    launch_window: tuple[str, str],
    planet_set: tuple[str, ...],
    *,
    max_legs: int = 4,
    vinf_grid_kms: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (60.0, 400.0),
    epoch_step_days: float = 30.0,
    chain_score_threshold: float | None = None,
    start_body_filter: Iterable[str] | None = None,
    a_range_au: tuple[float, float] = (0.3, 6.0),
    multi_shell: bool = False,
    pump_envelope_factor: float = 1.0,
) -> Iterator[MGAChainCandidate]:
    """Tisserand-Poincaré graph enumeration of MGA chains over an epoch window.

    Builds the (planet, V_inf-bin) graph for ``planet_set`` x ``vinf_grid_kms``,
    enumerates every Tisserand-admissible chain of length 2..(max_legs + 1)
    bodies, and emits one :class:`MGAChainCandidate` per (chain, epoch-grid-
    point) pair. The epoch grid steps ``launch_window`` in ``epoch_step_days``.

    This is *pure geometry* — no Lambert is solved. The caller wraps each
    candidate in :class:`EpochLockedTrajectory` and runs
    :func:`validate_chain_candidate` (which calls Phase 1's
    :func:`close_epoch_locked`) to validate.

    Parameters
    ----------
    launch_window:
        ``(t0_utc, t1_utc)`` window for the first body encounter.
    planet_set:
        Body codes to include (e.g. ``("E","V","M")`` for an inner-planet
        VEEGA-class search, ``("V","E","J")`` for a Galileo-style outer hop).
        Codes must be keys of :data:`cyclerfinder.core.constants.PLANETS`.
    max_legs:
        Maximum chain length in legs (``= len(sequence) - 1``). Default 4
        matches the (V,E,E,J) Galileo archetype length.
    vinf_grid_kms:
        Hyperbolic excess speed grid (km/s) — each chain conserves V_inf
        across all flybys, so a chain lives on a single grid point.
    tof_box_days_per_leg:
        ``(min_days, max_days)`` constraint on every leg's geometric TOF
        proposal. Same-body hops only generate resonance bands inside this
        box; hetero hops clamp the Hohmann half-period into the box.
    epoch_step_days:
        Spacing of the launch-epoch grid across ``launch_window``. Default
        30 days; the Phase-1 ``search_validity_window`` further refines a
        candidate to ±7 days at the 1-day grain.
    chain_score_threshold:
        Optional upper bound on the heuristic chain score. ``None``
        (default) emits every candidate; pass e.g. ``5.0`` to cap.
    start_body_filter:
        Optional restriction on the first body of every emitted chain (e.g.
        ``("V",)`` to only enumerate chains beginning at Venus). ``None``
        admits every body in ``planet_set``.
    a_range_au:
        ``(a_min, a_max)`` AU filter passed to :func:`linkable_3d`. Widen for
        outer-planet chains (Galileo VEEGA reaches Jupiter at 5.2 AU).
    multi_shell:
        If ``True``, allow per-flyby V_inf shifts within the Strange-Longuski
        2002 §12 pump envelope (:func:`delta_vinf_max_kms`). The Galileo
        VEEGA pattern (V_inf walks across shells: 4 → 9 → 9 → 5.6 km/s) is
        only representable with this on. Defaults to ``False`` (single-shell;
        Phase 2 semantics).
    pump_envelope_factor:
        Multiplier applied to the geometric envelope when
        ``multi_shell=True``. ``1.0`` (default) is the full Strange-Longuski
        2002 §12 bound; pass <1.0 to tighten (the realised shift is
        typically a fraction of the maximum because the V_inf rotation must
        target the next body rather than landing anywhere on the V_inf
        globe).

    Yields
    ------
    MGAChainCandidate
        Pure-geometry chain proposals, ordered by enumeration order (BFS).
        The caller may sort by ``chain_score`` for ranking. Lazy — the
        caller can break out at any time.

    Notes
    -----
    The Tisserand parameter is the V_inf invariant; we attach the first
    body's ``T_p`` as a sanity hook in :class:`MGAChainCandidate.tisserand_parameter`.

    Phase 2 emits one candidate per ``(sequence, vinf_bin, launch_epoch)``
    grid point. Phase 3 will tie a candidate to its precursor-MGA insertion
    target (the ``inserts_into`` field on :class:`EpochLockedTrajectory`)
    when the chain's terminal node matches an existing catalogue cycler's
    seed conditions.
    """
    if max_legs < 1:
        raise ValueError(f"max_legs must be >= 1; got {max_legs}")
    if not planet_set:
        raise ValueError("planet_set must be non-empty")
    if not vinf_grid_kms:
        raise ValueError("vinf_grid_kms must be non-empty")
    if epoch_step_days <= 0.0:
        raise ValueError(f"epoch_step_days must be positive; got {epoch_step_days}")
    if tof_box_days_per_leg[0] <= 0.0 or tof_box_days_per_leg[1] <= tof_box_days_per_leg[0]:
        raise ValueError(
            f"tof_box_days_per_leg must be (min, max) with 0 < min < max; "
            f"got {tof_box_days_per_leg}",
        )
    unknown = [b for b in planet_set if b not in PLANETS]
    if unknown:
        raise ValueError(
            f"unknown body code(s) in planet_set: {unknown!r}; valid codes are {tuple(PLANETS)!r}",
        )

    # Build the adjacency once per call. linkable_3d is expensive (~ms per
    # pair); caching it keeps the per-candidate cost in the BFS minimal.
    if multi_shell:
        adj = _multi_shell_adjacency(
            planet_set,
            vinf_grid_kms,
            a_range_au=a_range_au,
            pump_envelope_factor=pump_envelope_factor,
        )
    else:
        adj = _adjacency(planet_set, vinf_grid_kms, a_range_au=a_range_au)

    # Build the launch-epoch grid in days from launch_window[0].
    from astropy.time import Time

    t0 = Time(launch_window[0], scale="utc")
    t1 = Time(launch_window[1], scale="utc")
    span_days = float((t1 - t0).to("day").value)
    if span_days < 0.0:
        raise ValueError(
            f"launch_window must be ordered; got ({launch_window[0]!r}, {launch_window[1]!r})",
        )
    n_epochs = max(1, math.floor(span_days / epoch_step_days) + 1)
    epoch_grid_days = [i * epoch_step_days for i in range(n_epochs)]

    enumerator = _enumerate_multi_shell_sequences if multi_shell else _enumerate_sequences
    for sequence, vinf_tuple in enumerator(
        adj,
        planet_set,
        vinf_grid_kms,
        max_legs=max_legs,
        start_body_filter=start_body_filter,
    ):
        # Build the per-leg TOF proposal. If ANY same-body leg lacks a
        # resonance inside the TOF box we skip the chain — the box was set
        # too narrow for the chain to close. In multi-shell mode the leg's
        # geometric proposal uses the *arrival*-side V_inf (the next-body
        # bin); same-body legs use either side because V_inf is conserved.
        leg_tofs: list[float] = []
        leg_skip = False
        for k in range(len(sequence) - 1):
            # Single-shell: V_inf is bit-equal across the chain. Multi-shell:
            # the per-leg V_inf proposal uses the arrival bin (most relevant
            # for the next-leg Tisserand identity).
            vinf_for_leg = vinf_tuple[k + 1] if multi_shell else vinf_tuple[k]
            tof = _leg_tof_proposal(
                sequence[k],
                sequence[k + 1],
                vinf_for_leg,
                tof_box_days_per_leg,
            )
            if tof is None:
                leg_skip = True
                break
            leg_tofs.append(tof)
        if leg_skip:
            continue
        leg_tofs_tup = tuple(leg_tofs)

        # Tisserand parameter at the first body (chain-invariant under
        # Tisserand conservation; same V_inf bin throughout the chain).
        t_p = vinf_to_tisserand(sequence[0], vinf_tuple[0])

        score = _chain_score(sequence, vinf_tuple, leg_tofs_tup)
        if chain_score_threshold is not None and score > chain_score_threshold:
            continue

        for offset_days in epoch_grid_days:
            launch_utc = _add_days_utc(launch_window[0], offset_days)
            yield MGAChainCandidate(
                sequence=sequence,
                vinf_tuple_kms=vinf_tuple,
                leg_tofs_days=leg_tofs_tup,
                launch_epoch_utc=launch_utc,
                tisserand_parameter=float(t_p),
                chain_score=float(score),
            )


def validate_chain_candidate(
    candidate: MGAChainCandidate,
    ephemeris: Ephemeris,
    *,
    orbit_class: str = "mga_tour",
    n_returns: int = 1,
    inserts_into: str | None = None,
    periapsis_altitudes_km: tuple[float | None, ...] | None = None,
    closure_tol_kms: float = 0.5,
    flyby_continuity_tol_kms: float = 0.05,
    independent_cross_check: bool = True,
    independent_tol_kms: float = 0.1,
) -> EpochLockedClosure | None:
    """Wrap a candidate as :class:`EpochLockedTrajectory` and close via Phase 1.

    Builds the trajectory record using the candidate's geometric proposal as
    the EXPECTED V_inf side, runs :func:`close_epoch_locked` with the supplied
    tolerances, and returns the :class:`EpochLockedClosure` if it
    ``converged``. Returns ``None`` if closure fails (residual / flyby /
    cross-check above tolerance, or a Lambert/ephemeris exception).

    The validity window is set from the launch epoch + sum of leg TOFs (the
    geometric envelope; Phase 1's ``search_validity_window`` then refines).

    Parameters
    ----------
    candidate:
        The :class:`MGAChainCandidate` to validate.
    ephemeris:
        Phase-1 body-state provider. For real-ephemeris closure pass
        ``Ephemeris(model="astropy")``.
    orbit_class:
        Defaults to ``"mga_tour"`` — the Galileo-class label. Pass
        ``"precursor_mga"`` for an insertion candidate (then
        ``inserts_into`` is required by Phase 1's schema).
    n_returns:
        Defaults to 1 — every MGA tour visits each body once.
    inserts_into:
        Required only for ``orbit_class="precursor_mga"``.
    periapsis_altitudes_km:
        Optional per-encounter override (e.g. the canonical 100 km Mars
        periapsis); ``None`` lets every flyby use
        :data:`cyclerfinder.core.constants.SAFE_PERIHELION_KM`.
    closure_tol_kms, flyby_continuity_tol_kms, independent_cross_check,
    independent_tol_kms:
        Forwarded verbatim to :func:`close_epoch_locked`.

    Returns
    -------
    EpochLockedClosure | None
        The converged closure, or ``None`` if any gate failed.
    """
    total_tof_days = sum(candidate.leg_tofs_days)
    end_utc = _add_days_utc(candidate.launch_epoch_utc, total_tof_days)
    try:
        trajectory = EpochLockedTrajectory(
            sequence=candidate.sequence,
            leg_tofs_days=candidate.leg_tofs_days,
            vinf_kms_at_encounters=candidate.vinf_tuple_kms,
            launch_epoch_utc=candidate.launch_epoch_utc,
            orbit_class=orbit_class,  # type: ignore[arg-type]
            n_returns=n_returns,
            validity_window_start_utc=candidate.launch_epoch_utc,
            validity_window_end_utc=end_utc,
            inserts_into=inserts_into,
            periapsis_altitudes_km=periapsis_altitudes_km,
            notes="Tisserand-Poincaré MGA enumerator proposal (#298)",
        )
    except ValueError:
        return None
    try:
        closure = close_epoch_locked(
            trajectory,
            ephemeris,
            closure_tol_kms=closure_tol_kms,
            flyby_continuity_tol_kms=flyby_continuity_tol_kms,
            independent_cross_check=independent_cross_check,
            independent_tol_kms=independent_tol_kms,
        )
    except Exception:
        return None
    if not closure.converged:
        return None
    return closure


def _epoch_locked_loss(
    candidate: MGAChainCandidate,
    ephemeris: Ephemeris,
    *,
    orbit_class: str,
    n_returns: int,
    inserts_into: str | None,
    periapsis_altitudes_km: tuple[float | None, ...] | None,
    independent_cross_check: bool,
    independent_tol_kms: float,
    closure_residual_weight: float,
    flyby_continuity_weight: float,
) -> tuple[float, EpochLockedClosure | None]:
    """Run :func:`close_epoch_locked` on a candidate, return (loss, closure).

    Loss = ``closure_residual_weight * closure_residual_kms +
    flyby_continuity_weight * flyby_continuity_max_dv_kms``. Used by
    :func:`optimise_chain_tofs`. On any Lambert / ephemeris error returns
    a large penalty (1e6) and ``None``.
    """
    total_tof_days = sum(candidate.leg_tofs_days)
    end_utc = _add_days_utc(candidate.launch_epoch_utc, total_tof_days)
    try:
        trajectory = EpochLockedTrajectory(
            sequence=candidate.sequence,
            leg_tofs_days=candidate.leg_tofs_days,
            vinf_kms_at_encounters=candidate.vinf_tuple_kms,
            launch_epoch_utc=candidate.launch_epoch_utc,
            orbit_class=orbit_class,  # type: ignore[arg-type]
            n_returns=n_returns,
            validity_window_start_utc=candidate.launch_epoch_utc,
            validity_window_end_utc=end_utc,
            inserts_into=inserts_into,
            periapsis_altitudes_km=periapsis_altitudes_km,
            notes="Tisserand-Poincaré MGA TOF-opt loss eval (#300)",
        )
    except ValueError:
        return 1.0e6, None
    try:
        closure = close_epoch_locked(
            trajectory,
            ephemeris,
            # Use loose gates inside the loss so the loss surface is smooth.
            closure_tol_kms=1.0e6,
            flyby_continuity_tol_kms=1.0e6,
            independent_cross_check=independent_cross_check,
            independent_tol_kms=independent_tol_kms,
        )
    except Exception:
        return 1.0e6, None
    loss = (
        closure_residual_weight * closure.closure_residual_kms
        + flyby_continuity_weight * closure.flyby_continuity_max_dv_kms
    )
    return float(loss), closure


def optimise_chain_tofs(
    candidate: MGAChainCandidate,
    ephemeris: Ephemeris,
    *,
    orbit_class: str = "mga_tour",
    n_returns: int = 1,
    inserts_into: str | None = None,
    periapsis_altitudes_km: tuple[float | None, ...] | None = None,
    method: str = "Nelder-Mead",
    max_iter: int = 50,
    epoch_search_half_width_days: float = 10.0,
    tof_search_relative_half_width: float = 0.25,
    alpha_flyby_continuity: float = 2.0,
    accept_loss_kms: float | None = 1.0,
    independent_cross_check: bool = False,
    independent_tol_kms: float = 0.1,
) -> tuple[MGAChainCandidate, EpochLockedClosure, float] | None:
    """Optimise a candidate's (launch epoch, per-leg TOFs) to minimise closure loss.

    The free parameters are the launch epoch offset (days from the seed)
    and the per-leg TOFs (days). The loss function is

    .. math::

        L = \\text{closure_residual_kms}
            + \\alpha_\\text{flyby}\\,\\text{flyby_continuity_max_dv_kms}

    with :math:`\\alpha_\\text{flyby} = 2.0` by default (continuity is a
    hard ballistic constraint; closure residual is a softer V_inf-match).

    The search uses ``scipy.optimize.minimize`` with the method specified.
    ``Nelder-Mead`` is the default — it handles the noisy Lambert loss
    surface without needing gradients. The maximum iterations cap is 50 by
    default; for 4-leg chains that's a ~250-evaluation budget.

    Parameters
    ----------
    candidate:
        Seed :class:`MGAChainCandidate` (from :func:`find_mga_chains`).
    ephemeris:
        Phase-1 body-state provider.
    orbit_class, n_returns, inserts_into, periapsis_altitudes_km:
        Passed through to :class:`EpochLockedTrajectory` construction at
        each iteration.
    method:
        ``scipy.optimize.minimize`` method. Default ``"Nelder-Mead"``.
    max_iter:
        Max iterations. Default 50.
    epoch_search_half_width_days:
        Launch-epoch search window half-width (days). Default 10 — the
        optimiser may slide the launch epoch by ±10 days from the seed.
    tof_search_relative_half_width:
        Per-leg TOF search half-width as a fraction of the seed TOF.
        Default 0.25 (±25%). This keeps the TOFs inside a realistic
        geometric envelope.
    alpha_flyby_continuity:
        Weight on the flyby-continuity term in the loss. Default 2.0
        (per task brief).
    accept_loss_kms:
        Acceptance threshold on the final loss. ``None`` returns the
        optimised candidate regardless of loss; default ``1.0`` km/s
        returns ``None`` if the loss is above this.
    independent_cross_check:
        Whether to run the Kepler cross-check at each iteration. Default
        ``False`` for speed; the FINAL closure (returned alongside the
        optimised candidate) is re-run with the cross-check ON.
    independent_tol_kms:
        Cross-check tolerance (forwarded to the final closure).

    Returns
    -------
    (optimised_candidate, optimised_closure, final_loss) | None
        ``None`` if the optimiser fails to find a loss below
        ``accept_loss_kms``. Otherwise the optimised candidate (with
        updated ``launch_epoch_utc`` and ``leg_tofs_days``), the
        :class:`EpochLockedClosure` at the optimum (re-run with the
        cross-check on), and the final loss value.

    Notes
    -----
    The optimiser is local — it refines a seed, it does not search the
    global epoch / TOF space. The seed quality directly determines the
    achievable loss. Use :func:`find_mga_chains` to enumerate seeds.
    """
    from scipy.optimize import minimize

    n_legs = len(candidate.leg_tofs_days)
    seed_tofs = list(candidate.leg_tofs_days)
    seed_epoch = candidate.launch_epoch_utc

    def _x_to_candidate(x: NDArray[np.float64]) -> MGAChainCandidate:
        epoch_offset_days = float(x[0])
        new_tofs = tuple(max(1.0, float(t)) for t in x[1 : 1 + n_legs])
        new_epoch = _add_days_utc(seed_epoch, epoch_offset_days)
        return MGAChainCandidate(
            sequence=candidate.sequence,
            vinf_tuple_kms=candidate.vinf_tuple_kms,
            leg_tofs_days=new_tofs,
            launch_epoch_utc=new_epoch,
            tisserand_parameter=candidate.tisserand_parameter,
            chain_score=candidate.chain_score,
        )

    def _objective(x: NDArray[np.float64]) -> float:
        cand = _x_to_candidate(x)
        loss, _ = _epoch_locked_loss(
            cand,
            ephemeris,
            orbit_class=orbit_class,
            n_returns=n_returns,
            inserts_into=inserts_into,
            periapsis_altitudes_km=periapsis_altitudes_km,
            independent_cross_check=independent_cross_check,
            independent_tol_kms=independent_tol_kms,
            closure_residual_weight=1.0,
            flyby_continuity_weight=alpha_flyby_continuity,
        )
        return loss

    x0 = np.asarray([0.0, *seed_tofs], dtype=np.float64)
    # Initial simplex: ±10% on TOFs, ±5 days on epoch offset.
    perturbations: list[tuple[int, float]] = [(0, 5.0)]
    for k in range(n_legs):
        perturbations.append((1 + k, max(2.0, 0.1 * seed_tofs[k])))
    simplex_rows: list[NDArray[np.float64]] = [x0.copy()]
    for idx, delta in perturbations:
        x_p = x0.copy()
        x_p[idx] = x0[idx] + delta
        simplex_rows.append(x_p)
    initial_simplex = np.asarray(simplex_rows, dtype=np.float64)

    options: dict[str, Any] = {
        "maxiter": max_iter * (n_legs + 2),
        "xatol": 1.0,  # 1 day / 1 km-s scale
        "fatol": 0.05,  # 0.05 km/s loss tolerance
    }
    if method.lower() == "nelder-mead":
        options["initial_simplex"] = initial_simplex

    try:
        result = minimize(  # type: ignore[call-overload]
            _objective,
            x0,
            method=method,
            options=options,
        )
    except Exception:
        return None

    # Clamp the optimiser output back inside the search box, then evaluate.
    x_opt = np.asarray(result.x, dtype=np.float64).copy()
    x_opt[0] = max(
        -epoch_search_half_width_days,
        min(epoch_search_half_width_days, x_opt[0]),
    )
    for k in range(n_legs):
        seed_t = seed_tofs[k]
        lo = seed_t * (1.0 - tof_search_relative_half_width)
        hi = seed_t * (1.0 + tof_search_relative_half_width)
        x_opt[1 + k] = max(lo, min(hi, x_opt[1 + k]))

    opt_cand = _x_to_candidate(x_opt)
    # Final closure with cross-check ON (regardless of inner-loop setting).
    final_loss, final_closure = _epoch_locked_loss(
        opt_cand,
        ephemeris,
        orbit_class=orbit_class,
        n_returns=n_returns,
        inserts_into=inserts_into,
        periapsis_altitudes_km=periapsis_altitudes_km,
        independent_cross_check=True,
        independent_tol_kms=independent_tol_kms,
        closure_residual_weight=1.0,
        flyby_continuity_weight=alpha_flyby_continuity,
    )
    if final_closure is None:
        return None
    if accept_loss_kms is not None and final_loss > accept_loss_kms:
        return None
    return opt_cand, final_closure, float(final_loss)


def scan_window_and_validate(
    launch_window: tuple[str, str],
    planet_set: tuple[str, ...],
    ephemeris: Ephemeris,
    *,
    max_legs: int = 4,
    vinf_grid_kms: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (60.0, 400.0),
    epoch_step_days: float = 30.0,
    start_body_filter: Iterable[str] | None = None,
    chain_score_threshold: float | None = None,
    a_range_au: tuple[float, float] = (0.3, 6.0),
    closure_tol_kms: float = 0.5,
    flyby_continuity_tol_kms: float = 0.05,
    independent_cross_check: bool = True,
    independent_tol_kms: float = 0.1,
    periapsis_altitude_default_km: float | None = None,
    max_candidates: int | None = None,
    multi_shell: bool = False,
    pump_envelope_factor: float = 1.0,
) -> tuple[list[MGAChainCandidate], list[EpochLockedClosure]]:
    """Convenience driver: enumerate + close, return surviving pairs.

    Walks :func:`find_mga_chains`, validates each candidate with
    :func:`validate_chain_candidate`, and returns the (candidate, closure)
    pairs whose closure converged.

    Used by the Galileo VEEGA reproduction integration test and by Phase 2
    downstream callers that want a single-call survey of the window.

    Parameters
    ----------
    periapsis_altitude_default_km:
        Optional override applied uniformly at every intermediate
        encounter (e.g. relax the Earth flyby standoff for a Galileo-style
        965-km Earth-1 periapsis). ``None`` falls through to
        :data:`SAFE_PERIHELION_KM`.
    max_candidates:
        Stop after the enumerator emits this many candidates (regardless
        of survival). ``None`` (default) processes the entire stream.

    Returns
    -------
    (candidates, closures):
        Two lists of equal length, ordered by enumeration. ``candidates[i]``
        is the surviving :class:`MGAChainCandidate`; ``closures[i]`` is its
        Phase-1 :class:`EpochLockedClosure`. Empty lists if nothing closed.
    """
    surviving_candidates: list[MGAChainCandidate] = []
    surviving_closures: list[EpochLockedClosure] = []
    for seen, cand in enumerate(
        find_mga_chains(
            launch_window,
            planet_set,
            max_legs=max_legs,
            vinf_grid_kms=vinf_grid_kms,
            tof_box_days_per_leg=tof_box_days_per_leg,
            epoch_step_days=epoch_step_days,
            chain_score_threshold=chain_score_threshold,
            start_body_filter=start_body_filter,
            a_range_au=a_range_au,
            multi_shell=multi_shell,
            pump_envelope_factor=pump_envelope_factor,
        ),
        start=1,
    ):
        if max_candidates is not None and seen > max_candidates:
            break
        peri_override: tuple[float | None, ...] | None = None
        if periapsis_altitude_default_km is not None:
            # Override every intermediate flyby (not endpoints).
            override_list: list[float | None] = [None] * len(cand.sequence)
            for k in range(1, len(cand.sequence) - 1):
                override_list[k] = periapsis_altitude_default_km
            peri_override = tuple(override_list)
        closure = validate_chain_candidate(
            cand,
            ephemeris,
            periapsis_altitudes_km=peri_override,
            closure_tol_kms=closure_tol_kms,
            flyby_continuity_tol_kms=flyby_continuity_tol_kms,
            independent_cross_check=independent_cross_check,
            independent_tol_kms=independent_tol_kms,
        )
        if closure is not None:
            surviving_candidates.append(cand)
            surviving_closures.append(closure)
    return surviving_candidates, surviving_closures


__all__ = [
    "MGAChainCandidate",
    "delta_vinf_max_kms",
    "find_mga_chains",
    "optimise_chain_tofs",
    "scan_window_and_validate",
    "validate_chain_candidate",
]
