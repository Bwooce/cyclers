"""Planet-pair synodic periods and multi-body integer commensurabilities.

The cycler period bank: each candidate cycler period is an integer multiple
``k`` of a synodic period (the time between successive same-geometry
alignments of two planets), and a multi-body cycler must satisfy several
such commensurabilities simultaneously. The natural Venus-Earth-Mars (VEM)
beat is ``3 * T_syn(E,M) ~ 4 * T_syn(E,V) ~ 6.4 yr`` (spec §3, §9).

All periods are returned in days; year conversions use
:data:`cyclerfinder.core.constants.DAYS_PER_JULIAN_YEAR` (365.25, the Julian
year, matching ``astropy.units.yr``).

Reference: spec §3 (synodic background), §5 step 1 (resonance pipeline),
§9 (gate anchors), §13 (search prioritisation).

Plan: ``docs/phases/m2-flyby-maps/plan.md`` §3.3.
"""

from __future__ import annotations

from itertools import product

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, PLANETS

# ---------------------------------------------------------------------------
# Pairwise synodic periods
# ---------------------------------------------------------------------------


def _orbital_period_days(body: str) -> float:
    """Sidereal orbital period (days) of ``body`` from its mean motion."""
    n = PLANETS[body].mean_motion_deg_day  # deg/day
    return 360.0 / n


def synodic_period_days(body_a: str, body_b: str) -> float:
    """Synodic period of ``body_a`` and ``body_b`` in days.

    ``1 / |1/T_a - 1/T_b|``, where ``T_x`` is the sidereal orbital period
    derived from :data:`PLANETS[x].mean_motion_deg_day`.

    Parameters
    ----------
    body_a, body_b:
        One-letter planet codes from :data:`PLANETS`.

    Returns
    -------
    float
        Synodic period in days.

    Raises
    ------
    ValueError
        If ``body_a == body_b`` (the synodic period is undefined for
        identical orbits; ``1/T - 1/T = 0`` would yield infinity).
    KeyError
        If either code is not in :data:`PLANETS`.
    """
    if body_a == body_b:
        raise ValueError(f"synodic_period_days requires distinct bodies; got both = {body_a!r}")
    t_a = _orbital_period_days(body_a)
    t_b = _orbital_period_days(body_b)
    return 1.0 / abs(1.0 / t_a - 1.0 / t_b)


def synodic_period_years(body_a: str, body_b: str) -> float:
    """Synodic period in Julian years (``synodic_period_days / 365.25``).

    Convenience for direct comparison with the spec §9 anchors
    (E-M = 2.135 yr, E-V = 1.599 yr).
    """
    return synodic_period_days(body_a, body_b) / DAYS_PER_JULIAN_YEAR


def k_synodic_periods_days(body_a: str, body_b: str, k_max: int) -> list[float]:
    """Candidate cycler periods (days) at the body pair: ``[k * T_syn for k in 1..k_max]``.

    This is the period bank M3 uses for the Aldrin reproduction (``k=1`` on
    E-M ~ 2.135 yr) and M4 uses to enumerate periods per cell.

    Parameters
    ----------
    body_a, body_b:
        Distinct planet codes.
    k_max:
        Maximum integer multiple (inclusive). Must be ``>= 1``.

    Returns
    -------
    list[float]
        ``[T_syn, 2*T_syn, ..., k_max * T_syn]`` in days.

    Raises
    ------
    ValueError
        If ``k_max < 1`` or if ``body_a == body_b``.
    """
    if k_max < 1:
        raise ValueError(f"k_max must be >= 1; got {k_max}")
    t_syn = synodic_period_days(body_a, body_b)
    return [k * t_syn for k in range(1, k_max + 1)]


# ---------------------------------------------------------------------------
# Multi-body beats
# ---------------------------------------------------------------------------


def _non_ref_bodies(bodies: list[str], ref: str) -> list[str]:
    """Return ``bodies`` minus the (first) occurrence of ``ref``, preserving order."""
    out: list[str] = []
    seen_ref = False
    for b in bodies:
        if b == ref and not seen_ref:
            seen_ref = True
            continue
        out.append(b)
    return out


def _choose_reference(bodies: list[str]) -> str:
    """Reference-body rule per plan §3.3.

    For a 2-body set the rule is degenerate (any of the two works; pick
    the first). For a 3-body set the middle element is the reference
    (so ``[V, E, M]`` uses E, giving co-resonances ``T_syn(V,E)`` and
    ``T_syn(E,M)``). For N >= 4 the same "middle element" rule applies
    pragmatically; M2 only validates the 3-body case and the N >= 4
    extension is M8 territory.
    """
    n = len(bodies)
    if n < 2:
        raise ValueError(f"multi-body resonance requires >= 2 bodies; got {n}")
    if n == 2:
        return bodies[0]
    return bodies[n // 2]


def multi_body_beat_days(
    bodies: list[str],
    k_max: int = 6,
    tol_frac: float = 0.02,
) -> list[tuple[int, ...]]:
    """Find integer commensurabilities of pairwise synodic periods.

    For a body set ``[b_0, ..., b_{N-1}]`` with reference body chosen per
    plan §3.3 (middle element for N=3; first for N=2), search integer
    tuples ``(k_1, ..., k_{N-1})`` with each ``k_i in 1..k_max`` such that
    every ``k_i * T_syn(b_i, ref)`` agrees to within ``tol_frac``
    (fractional) of every other.

    The tuple ordering follows ``[b for b in bodies if b != ref]`` (input
    order, reference removed). For ``["V", "E", "M"]`` the natural beat
    ``3 * T_syn(E,M) ~ 4 * T_syn(E,V)`` returns as ``(4, 3)`` because the
    non-reference bodies in input order are ``[V, M]``. The same tuple
    fed to :func:`beat_period_days` reproduces ~6.406 yr.

    Returns
    -------
    list[tuple[int, ...]]
        Tuples within the tolerance, ranked by ascending fractional
        mismatch. May be empty; usually a single dominant entry but
        near-misses are kept for transparency.

    Raises
    ------
    ValueError
        On fewer than 2 bodies or ``k_max < 1``.
    """
    if k_max < 1:
        raise ValueError(f"k_max must be >= 1; got {k_max}")
    ref = _choose_reference(bodies)
    non_ref = _non_ref_bodies(bodies, ref)

    # Pre-compute each non-reference body's synodic period to the reference.
    syn = [synodic_period_days(b, ref) for b in non_ref]

    matches: list[tuple[float, tuple[int, ...]]] = []
    for tup in product(range(1, k_max + 1), repeat=len(non_ref)):
        scaled = [k * s for k, s in zip(tup, syn, strict=True)]
        mean = sum(scaled) / len(scaled)
        max_dev = max(abs(x - mean) for x in scaled)
        frac = max_dev / mean
        if frac <= tol_frac:
            matches.append((frac, tup))
    matches.sort(key=lambda kv: kv[0])
    return [tup for _frac, tup in matches]


def beat_period_days(bodies: list[str], k_tuple: tuple[int, ...]) -> float:
    """Beat period (days) implied by a commensurability tuple.

    Given a ``k_tuple`` returned by :func:`multi_body_beat_days` for the
    same ``bodies`` argument, returns the mean of
    ``k_i * T_syn(b_i, ref)`` across the non-reference bodies. This is the
    period the tuple corresponds to (well-defined when the tuple is within
    the fractional tolerance; "approximately equal" times reduce to a mean).

    Raises
    ------
    ValueError
        If ``len(k_tuple) != len(bodies) - 1``.
    """
    ref = _choose_reference(bodies)
    non_ref = _non_ref_bodies(bodies, ref)
    if len(k_tuple) != len(non_ref):
        raise ValueError(
            f"k_tuple length {len(k_tuple)} != non-reference body count {len(non_ref)}"
        )
    scaled = [k * synodic_period_days(b, ref) for k, b in zip(k_tuple, non_ref, strict=True)]
    return sum(scaled) / len(scaled)
