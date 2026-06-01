"""Real-ephemeris minimum-ΔV periodic optimiser for the Aldrin E-M cycler.

The classic Aldrin Earth-Mars cycler (McConaghy/Longuski/Byrnes 2002, AIAA
2002-4420, Table 4 row "1L1") is a *powered* cycler: the **Earth** (geocentric)
flyby cannot deliver the required turn ballistically (≈84° required vs ≈72°
achievable at a 200 km Earth flyby), so a non-zero maintenance ΔV is unavoidable
even in the idealised model. McConaghy's dissertation ("Global Search and
Optimization for Free-Return Earth-Mars Cyclers") tabulates the *Required
Geocentric Turning Angle at each Flyby*, bases the maximum allowed turn on a
"200 km altitude Earth flyby", and chooses the turn-ratio cutoff TRMIN = 0.85
specifically to include the Aldrin cycler (72 / 84.7 ≈ 0.85). No ΔV magnitude is
published anywhere in the literature — McConaghy 2002 explicitly defers it — so
the ΔV this optimiser computes is **our value only** and must never be used as a
golden-test assertion target.

What *is* source-attested are the orbital anchors (Rogers et al. 2012 Table 1;
Russell 2004 Table 3.4; McConaghy 2002 Table 4 row "1L1"):

* semi-major axis a = 1.60 AU
* eccentricity e = 0.393
* V∞ at Earth = 6.5 km/s ; V∞ at Mars = 9.7-9.75 km/s
* Earth→Mars leg time-of-flight = 146 days

This module finds the periodic E→M→E slice that minimises the per-synodic
maintenance ΔV over real planet positions, and reports the computed anchors so
the caller can validate them against the published values.

Formulation (locked)
---------------------
Free variables ``x = [t0_sec, tof1_days, tof2_days]`` — the Earth-departure
epoch and the two leg times-of-flight. For a trial ``x`` we build encounter
epochs ``[t0, t0 + tof1, t0 + tof1 + tof2]`` with sequence ``["E", "M", "E"]``,
Lambert-solve each leg (single-rev, single branch) over the supplied
ephemeris, and charge

    f(x) = flyby_dv_for("M", vinf_in_M, vinf_out_M)
         + flyby_dv_for("E", vinf_in_E_arrival, vinf_out_E_departure)

The Earth term closes the loop: the arrival V∞ at the final Earth encounter
must be bent to match the departure V∞ at the first Earth encounter
(periodicity). This sum is the per-synodic maintenance ΔV.

The minimiser mirrors :func:`cyclerfinder.search.optimize.optimise_cell_idealized`
— scipy ``differential_evolution`` for global coverage followed by an SLSQP
polish — but operates on the focused ``[t0, tof1, tof2]`` parameterisation
directly rather than through the ``Cell`` abstraction (a catalogue cycler is an
impedance mismatch for ``Cell``; the ephemeris-mode ``Cell`` stub stays locked).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import acos, degrees, sqrt
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import differential_evolution, minimize

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SAFE_PERIHELION_KM,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import dv_from_turn_deficit, flyby_dv_for, max_bend
from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError
from cyclerfinder.model import Cycler
from cyclerfinder.model.cycler import orbit_elements_au
from cyclerfinder.search.construct import build_aldrin_seed, construct_cycler
from cyclerfinder.search.resonance import synodic_period_days

# ---------------------------------------------------------------------------
# Module-level constants / tolerances
# ---------------------------------------------------------------------------

_PATHOLOGICAL_OBJECTIVE_KMS: float = 1.0e6
"""Penalty returned when the trial cycler cannot be built (Lambert pathology
or non-monotonic epochs). Large enough to be unambiguously worse than any
feasible objective, finite so DE / SLSQP do not choke on ``inf`` / ``nan``."""

_TOF1_BOUNDS_DAYS: tuple[float, float] = (100.0, 250.0)
"""Earth→Mars leg ToF bounds (days). Brackets the published 146 d."""

_TOF2_BOUNDS_DAYS: tuple[float, float] = (400.0, 900.0)
"""Mars→Earth return-leg ToF bounds (days). Brackets ~634 d (the E-M synodic
period of ~780 d minus the ~146 d outbound leg)."""

_T0_WINDOW_SYNODIC_FRAC: float = 0.15
"""Half-width of the ``t0`` search window, in E-M synodic periods, around the
``t0_guess``.

This is deliberately *narrow*. The Earth-departure epoch ``t0`` selects the
synodic launch *phase* — i.e. *which* periodic E→M→E cycler family you are on.
The pure minimum-ΔV objective is degenerate over a wide ``t0`` window: several
distinct families are ballistically feasible (ΔV ≈ 0) under the powered-flyby
surrogate, so a wide global search slides off the published Aldrin anchors onto
a neighbouring lower-V∞ ballistic family (empirically a ≈ 1.587 AU, V∞_E ≈ 5.1
km/s at ToF1 ≈ 127 d). The published Aldrin is the family at the phase the seed
inverts from the 132° transfer geometry; ±0.15 synodic periods keeps the search
inside that launch phase while still leaving the two leg ToFs fully free. See
the module docstring and the ``test_maintain`` detuned-start test, which
exercises the wider degeneracy explicitly."""

_T0_TIE_BREAK_KMS: float = 1.0e-3
"""When two candidate optima have maintenance ΔV within this band, prefer the
one closest to ``t0_guess``. The min-ΔV objective is flat (≈ 0) across the
Aldrin basin, so without a tie-break the polish could wander along the plateau;
this pins the result to the Aldrin launch phase deterministically."""

_DE_MAXITER: int = 60
_DE_POPSIZE: int = 12
_DE_TOL: float = 1.0e-5
_SLSQP_MAXITER: int = 200
_SLSQP_FTOL: float = 1.0e-8

_CONVERGENCE_OBJECTIVE_KMS: float = 1.0e3
"""A finite optimum below this is considered a real (non-pathological)
solution; at or above it every trial hit the Lambert penalty and the run is
flagged ``converged=False``."""

_ALDRIN_EARTH_FLYBY_ALT_KM: float = 200.0
"""Sourced Earth-flyby periapsis altitude (km) for the classic Aldrin cycler.

McConaghy's dissertation bases the maximum allowable geocentric turn on a
"200 km altitude Earth flyby"; the published ≈72° achievable turn (and the
TRMIN = 0.85 cutoff chosen to include Aldrin) follow from it. We adopt the same
sourced altitude rather than the conservative 300 km default so the recovered
turn deficit reproduces the literature geometry. The value is also recorded
per-body in the catalogue entry's ``flyby_mechanics`` so it is configurable per
orbit."""


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MaintenanceOptimResult:
    """Outcome of :func:`optimise_aldrin_maintenance_dv`.

    Attributes
    ----------
    cycler:
        The optimised real-ephemeris ``["E", "M", "E"]`` cycler.
    t0_sec:
        Optimised Earth-departure epoch, seconds since J2000.
    leg_tofs_days:
        Optimised leg times-of-flight (Earth→Mars, Mars→Earth), days.
    maintenance_dv_kms:
        **Computed** per-synodic maintenance ΔV (km/s) — the turn-deficit ΔV
        the return flyby (Earth for Aldrin) must supply to keep the recovered
        orbit repeating (see :attr:`turn_deficit`). This is our own number; no
        published counterpart exists (McConaghy 2002 defers it), so it must be
        reported / sanity-bounded only, never asserted against a sourced
        target. The *turn angles*, by contrast, are source-traceable.
    per_encounter_dv_kms:
        ``(body, flyby_dv)`` breakdown summing to ``maintenance_dv_kms``.
    converged:
        ``True`` when a real (non-penalty) optimum was found.
    a_au:
        Computed semi-major axis of the Earth→Mars transfer ellipse, AU.
    e:
        Computed eccentricity of the Earth→Mars transfer ellipse.
    vinf_kms_at_encounters:
        Computed ``(body, |V∞|)`` per encounter, km/s.
    turn_deficit:
        The :class:`FlybyTurnDeficit` at the return body (Earth for Aldrin),
        computed from the recovered ``(a, e)``. Its ``turn_required_deg`` /
        ``turn_max_deg`` are the per-orbit turn numbers (≈84° / ≈72° for
        Aldrin) — the source-traceable validation target. ``None`` if the
        orbit does not cross the return body's orbit (or the run did not
        converge).
    """

    cycler: Cycler
    t0_sec: float
    leg_tofs_days: tuple[float, ...]
    maintenance_dv_kms: float
    per_encounter_dv_kms: tuple[tuple[str, float], ...]
    converged: bool
    a_au: float
    e: float
    vinf_kms_at_encounters: tuple[tuple[str, float], ...]
    turn_deficit: FlybyTurnDeficit | None


# ---------------------------------------------------------------------------
# Sourced-geometry turn deficit (the faithful "powered" evidence)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FlybyTurnDeficit:
    """Turn a return-body flyby must supply each synodic period, vs the most it
    can supply ballistically.

    Computed analytically from the *idealised* coplanar cycler orbit
    (semi-major axis ``a``, eccentricity ``e``) and the return body's circular
    heliocentric orbit. The cycler crosses the body's orbit radius twice per
    revolution — once outbound (radial velocity > 0) and once inbound (< 0) —
    with identical ``|V∞|``. For the orbit to repeat, the flyby must rotate the
    inbound-arrival ``V∞`` onto the outbound-departure ``V∞``; the angle between
    them is :attr:`turn_required_deg`. The body can ballistically deliver at
    most :attr:`turn_max_deg` (see :func:`cyclerfinder.core.flyby.max_bend`).

    For the classic Aldrin E-M cycler the binding body is **Earth**: this
    reproduces McConaghy/Longuski/Byrnes 2002 (AIAA 2002-4420, Table 4 "1L1")
    — 84° required against 72° achievable — directly from the sourced
    ``a = 1.60 AU, e = 0.393`` anchors. The required turn depends only on the
    sourced elements (so it is itself source-traceable); the achievable turn
    additionally depends on our assumed safe periapsis. The resulting
    :attr:`dv_kms` is **our computed value** (McConaghy defers the ΔV) and must
    never be a golden-test target.
    """

    body: str
    vinf_kms: float
    turn_required_deg: float
    turn_max_deg: float
    deficit_deg: float
    dv_kms: float
    ballistically_feasible: bool


def idealized_flyby_turn_deficit(
    a_au: float,
    e: float,
    body: str,
    *,
    mu_sun: float = MU_SUN_KM3_S2,
    flyby_alt_km: float | None = None,
) -> FlybyTurnDeficit | None:
    """Turn deficit at ``body`` for a coplanar cycler orbit ``(a_au, e)``.

    Returns ``None`` when the orbit does not cross the body's circular orbit
    (the body is never encountered), otherwise a :class:`FlybyTurnDeficit`.

    The geometry is evaluated in the orbital plane using a (radial,
    transverse) decomposition of the heliocentric velocity at the body's
    orbital radius; the body's velocity is purely transverse at circular
    speed. See :class:`FlybyTurnDeficit` for the physical interpretation and
    sourcing.

    Parameters
    ----------
    flyby_alt_km:
        Flyby periapsis altitude above the body's equatorial radius (km),
        which sets the maximum ballistic turn (smaller altitude → tighter
        flyby → larger achievable turn). Configurable per orbit and per body;
        a catalogue entry records the source's value (the classic Aldrin
        cycler uses a 200 km Earth flyby per McConaghy's dissertation). ``None``
        falls back to the conservative per-body default
        (:attr:`~cyclerfinder.core.constants.PlanetData.safe_alt_km`, 300 km).
    """
    r_body = PLANETS[body].sma_au * AU_KM
    a_km = a_au * AU_KM
    r_peri = a_km * (1.0 - e)
    r_apo = a_km * (1.0 + e)
    if not (r_peri <= r_body <= r_apo):
        return None

    v = sqrt(mu_sun * (2.0 / r_body - 1.0 / a_km))
    h = sqrt(mu_sun * a_km * (1.0 - e * e))
    v_t = h / r_body
    v_r = sqrt(max(0.0, v * v - v_t * v_t))
    v_circ = sqrt(mu_sun / r_body)

    # In/out V∞ in the (radial, transverse) frame; both share |V∞|.
    dv_t = v_t - v_circ
    vinf = sqrt(v_r * v_r + dv_t * dv_t)
    if vinf == 0.0:
        return None
    cos_turn = (-(v_r * v_r) + dv_t * dv_t) / (vinf * vinf)
    turn_req = acos(max(-1.0, min(1.0, cos_turn)))
    if flyby_alt_km is None:
        r_peri_flyby = SAFE_PERIHELION_KM[body]
    else:
        r_peri_flyby = PLANETS[body].radius_eq_km + flyby_alt_km
    turn_max = max_bend(PLANETS[body].mu_km3_s2, r_peri_flyby, vinf)
    dv = dv_from_turn_deficit(vinf, turn_req, turn_max)

    return FlybyTurnDeficit(
        body=body,
        vinf_kms=vinf,
        turn_required_deg=degrees(turn_req),
        turn_max_deg=degrees(turn_max),
        deficit_deg=max(0.0, degrees(turn_req) - degrees(turn_max)),
        dv_kms=dv,
        ballistically_feasible=turn_req <= turn_max,
    )


# ---------------------------------------------------------------------------
# Forward map + objective
# ---------------------------------------------------------------------------


def _build_chain(
    x: NDArray[np.float64], sequence: Sequence[str], ephem: Ephemeris
) -> Cycler | None:
    """Build the ``sequence`` cycler from ``x = [t0, tof_1, …, tof_{n-1}]``.

    ``x[0]`` is the departure epoch (seconds since J2000); the remaining
    ``len(sequence) - 1`` elements are per-leg times-of-flight (days). Encounter
    epochs are the running cumulative sum. Returns ``None`` on any Lambert /
    construction pathology so the objective can convert it into a finite penalty
    rather than propagating ``inf`` / ``nan`` into scipy.
    """
    t0 = float(x[0])
    tofs_sec = [float(v) * SECONDS_PER_DAY for v in x[1:]]
    if any(tof <= 0.0 for tof in tofs_sec):
        return None
    times = [t0]
    for tof in tofs_sec:
        times.append(times[-1] + tof)
    try:
        return construct_cycler(
            sequence=list(sequence),
            encounter_times_sec=times,
            ephem=ephem,
        )
    except (ValueError, LambertConvergenceError, LambertGeometryError):
        return None


def _maintenance_dv_chain(cycler: Cycler) -> float:
    """Per-synodic maintenance ΔV (km/s) for an arbitrary closed cycler.

    Sums the turn-deficit ΔV at each *intermediate* encounter (bending its
    in-V∞ onto its out-V∞) plus the closure term that bends the final arrival
    V∞ onto the first departure V∞ (periodicity). For the Aldrin ``E-M-E``
    sequence this is exactly ``dv_mars + dv_earth``: the single intermediate
    Mars flyby plus the Earth loop-closure.
    """
    total = 0.0
    for enc in cycler.encounters[1:-1]:
        total += flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out)
    first = cycler.encounters[0]
    last = cycler.encounters[-1]
    total += flyby_dv_for(first.body, last.vinf_in, first.vinf_out)
    return total


def _objective(x: NDArray[np.float64], sequence: Sequence[str], ephem: Ephemeris) -> float:
    """Per-synodic maintenance ΔV (km/s); finite penalty on construction
    failure."""
    cycler = _build_chain(x, sequence, ephem)
    if cycler is None:
        return _PATHOLOGICAL_OBJECTIVE_KMS
    return _maintenance_dv_chain(cycler)


def _first_distinct_pair(sequence: Sequence[str]) -> tuple[str, str]:
    """First ``(body_a, body_b)`` pair of distinct bodies in ``sequence``.

    Used to size the ``t0`` search window by their synodic period when the
    caller does not name an explicit ``synodic_pair``. For ``E-M-E`` this is
    ``("E", "M")``.
    """
    first = sequence[0]
    for body in sequence[1:]:
        if body != first:
            return (first, body)
    raise ValueError(f"sequence {list(sequence)!r} has no two distinct bodies")


# ---------------------------------------------------------------------------
# Public optimiser
# ---------------------------------------------------------------------------


def _default_t0_guess(em_tof_days: float) -> float:
    """Earth-departure epoch (s) that produces the Aldrin 132° transfer arc.

    Mirrors :func:`cyclerfinder.search.construct.build_aldrin_seed`'s inversion
    of the phase equation by building the seed and reading its first-encounter
    epoch — this keeps the single source of truth for the seed phase in
    ``construct`` rather than duplicating the algebra here.
    """
    seed = build_aldrin_seed(Ephemeris("circular"), em_tof_days=em_tof_days)
    return float(seed.encounters[0].t)


def optimise_maintenance_dv(
    sequence: Sequence[str],
    ephem: Ephemeris,
    *,
    t0_guess_sec: float,
    tof_days_guesses: Sequence[float],
    tof_bounds_days: Sequence[tuple[float, float]],
    synodic_pair: tuple[str, str] | None = None,
    closure_body: str | None = None,
    closure_flyby_alt_km: float | None = None,
    t0_window_synodic_frac: float = _T0_WINDOW_SYNODIC_FRAC,
    tof_jitter_half_days: Sequence[float] | None = None,
    n_starts: int = 5,
    seed: int = 0,
    seed_cycler_factory: Callable[[], Cycler] | None = None,
) -> MaintenanceOptimResult:
    """Find the minimum-ΔV periodic cycler for an arbitrary closed ``sequence``.

    Body-agnostic generalisation of :func:`optimise_aldrin_maintenance_dv`. The
    free variables are ``x = [t0, tof_1, …, tof_{n-1}]`` (departure epoch in
    seconds plus ``len(sequence) - 1`` leg times-of-flight in days). A global
    ``differential_evolution`` pass over a narrow ``t0`` window (sized by the
    ``synodic_pair`` synodic period) followed by SLSQP polishes from the seed
    and jittered restarts; near-ties in the flat ΔV plateau break toward
    ``t0_guess_sec``.

    Parameters
    ----------
    sequence:
        Body-code encounter sequence; first and last codes must match (closed
        cycler), e.g. ``["E", "M", "E"]``.
    ephem:
        Planet-state provider (``"circular"`` for the test surface, ``"astropy"``
        for real DE440 states).
    t0_guess_sec:
        Departure-epoch guess (seconds since J2000) — selects the launch phase.
    tof_days_guesses:
        ``len(sequence) - 1`` initial leg ToFs (days).
    tof_bounds_days:
        ``len(sequence) - 1`` ``(lo, hi)`` per-leg ToF bounds (days).
    synodic_pair:
        Body pair whose synodic period sizes the ``t0`` window. ``None`` derives
        the first distinct pair from ``sequence``.
    closure_body:
        Body whose idealised turn-deficit ΔV is reported as the maintenance
        cost (the source-traceable "powered" evidence). ``None`` reports the
        Lambert closure ΔV from :func:`_maintenance_dv_chain` instead.
    closure_flyby_alt_km:
        Flyby periapsis altitude (km) for the turn-deficit computation. ``None``
        uses the per-body conservative default.
    t0_window_synodic_frac:
        Half-width of the ``t0`` search window in synodic periods.
    tof_jitter_half_days:
        ``len(sequence) - 1`` per-leg jitter half-widths (days) for the
        multi-start restarts. ``None`` defaults to zero jitter on every leg.
    n_starts:
        Number of SLSQP polish starts (seed plus ``n_starts - 1`` jittered),
        in addition to the global DE pass.
    seed:
        RNG seed for ``differential_evolution`` and the multi-start jitter.
    seed_cycler_factory:
        Builds a fallback cycler if every trial is pathological even at the
        seed. ``None`` raises ``RuntimeError`` in that case.

    Returns
    -------
    MaintenanceOptimResult
        With the optimised cycler, computed anchors, and the (computed,
        not source-attested) maintenance ΔV.
    """
    n_legs = len(sequence) - 1
    if len(tof_days_guesses) != n_legs:
        raise ValueError(f"tof_days_guesses needs {n_legs} entries, got {len(tof_days_guesses)}")
    if len(tof_bounds_days) != n_legs:
        raise ValueError(f"tof_bounds_days needs {n_legs} entries, got {len(tof_bounds_days)}")
    if tof_jitter_half_days is None:
        tof_jitter_half_days = tuple(0.0 for _ in range(n_legs))
    elif len(tof_jitter_half_days) != n_legs:
        raise ValueError(
            f"tof_jitter_half_days needs {n_legs} entries, got {len(tof_jitter_half_days)}"
        )
    if synodic_pair is None:
        synodic_pair = _first_distinct_pair(sequence)

    syn_sec = synodic_period_days(*synodic_pair) * SECONDS_PER_DAY
    t0_lo = t0_guess_sec - t0_window_synodic_frac * syn_sec
    t0_hi = t0_guess_sec + t0_window_synodic_frac * syn_sec
    bounds: list[tuple[float, float]] = [(t0_lo, t0_hi)]
    bounds.extend((float(lo), float(hi)) for lo, hi in tof_bounds_days)

    # --- Global pass: differential_evolution (mirrors optimise_cell_idealized).
    de_any = cast(Any, differential_evolution)
    candidate_starts: list[NDArray[np.float64]] = []
    try:
        de_result: Any = de_any(
            _objective,
            bounds,
            args=(sequence, ephem),
            seed=seed,
            polish=False,
            maxiter=_DE_MAXITER,
            popsize=_DE_POPSIZE,
            tol=_DE_TOL,
        )
        candidate_starts.append(np.asarray(de_result.x, dtype=np.float64))
    except (ValueError, RuntimeError):
        pass

    # --- Multi-start seeds: the phase-inverted guess plus jittered variants.
    seed_x = np.array([t0_guess_sec, *tof_days_guesses], dtype=np.float64)
    candidate_starts.append(seed_x)
    rng = np.random.default_rng(seed)
    for _ in range(max(0, n_starts - 1)):
        # Draw order is load-bearing for bit-reproducibility: t0 jitter first,
        # then each leg jitter in sequence order.
        jitter_values = [rng.uniform(-0.25 * syn_sec, 0.25 * syn_sec)]
        for half in tof_jitter_half_days:
            jitter_values.append(rng.uniform(-half, half))
        candidate_starts.append(seed_x + np.array(jitter_values, dtype=np.float64))

    # --- SLSQP polish from every start; keep the best feasible optimum,
    # breaking near-ties in the (flat, ≈0) ΔV plateau by proximity to the
    # launch phase ``t0_guess`` so the result is deterministic and does not
    # drift along the plateau onto a neighbouring cycler family.
    best_x: NDArray[np.float64] | None = None
    best_f = np.inf
    lower = np.array([b[0] for b in bounds], dtype=np.float64)
    upper = np.array([b[1] for b in bounds], dtype=np.float64)
    for start in candidate_starts:
        x0 = np.clip(start, lower, upper)
        try:
            res: Any = minimize(
                _objective,
                x0,
                args=(sequence, ephem),
                method="SLSQP",
                bounds=bounds,
                options={"maxiter": _SLSQP_MAXITER, "ftol": _SLSQP_FTOL},
            )
        except (ValueError, RuntimeError):
            continue
        x_final = np.asarray(res.x, dtype=np.float64)
        f_final = _objective(x_final, sequence, ephem)
        if best_x is None:
            best_f, best_x = f_final, x_final
            continue
        improves = f_final < best_f - _T0_TIE_BREAK_KMS
        ties = abs(f_final - best_f) <= _T0_TIE_BREAK_KMS
        closer = abs(x_final[0] - t0_guess_sec) < abs(best_x[0] - t0_guess_sec)
        if improves or (ties and closer):
            best_f, best_x = f_final, x_final

    # Fallback to the raw seed if every polish failed to build a cycler.
    if best_x is None or best_f >= _PATHOLOGICAL_OBJECTIVE_KMS:
        best_x = np.clip(seed_x, lower, upper)
        best_f = _objective(best_x, sequence, ephem)

    cycler = _build_chain(best_x, sequence, ephem)
    converged = cycler is not None and best_f < _CONVERGENCE_OBJECTIVE_KMS
    leg_tofs_days = tuple(float(v) for v in best_x[1:])
    if cycler is None:
        # Pathological even at the seed — surface the failure honestly.
        if seed_cycler_factory is None:
            raise RuntimeError(
                f"sequence {list(sequence)!r} is pathological at the seed and no "
                "seed_cycler_factory was provided to build a fallback cycler"
            )
        cycler = seed_cycler_factory()
        return MaintenanceOptimResult(
            cycler=cycler,
            t0_sec=float(best_x[0]),
            leg_tofs_days=leg_tofs_days,
            maintenance_dv_kms=float(best_f),
            per_encounter_dv_kms=(),
            converged=False,
            a_au=float("nan"),
            e=float("nan"),
            vinf_kms_at_encounters=(),
            turn_deficit=None,
        )

    a_au, e = orbit_elements_au(cycler.encounters[0].r, cycler.legs[0].v_depart, MU_SUN_KM3_S2)
    vinf_per_enc = tuple(
        (enc.body, float(np.linalg.norm(enc.vinf_in))) for enc in cycler.encounters
    )

    # The faithful maintenance ΔV is the turn deficit the *closure* flyby
    # (Earth for Aldrin) must make up to keep the recovered orbit repeating — a
    # geometric property of the sourced (a, e). The free-ToF Lambert closure
    # term (``_maintenance_dv_chain``) instead slides onto a cheaper neighbouring
    # ballistic geometry (a ~32° Earth turn, ΔV≈0) and is therefore NOT used as
    # the reported cost; it survives only as the optimiser objective that pins
    # the orbit anchors.
    turn_deficit: FlybyTurnDeficit | None = None
    if closure_body is not None:
        turn_deficit = idealized_flyby_turn_deficit(
            float(a_au), float(e), closure_body, flyby_alt_km=closure_flyby_alt_km
        )
    if turn_deficit is not None:
        maintenance_dv = turn_deficit.dv_kms
        per_enc_dv: tuple[tuple[str, float], ...] = ((turn_deficit.body, turn_deficit.dv_kms),)
    else:
        maintenance_dv = _maintenance_dv_chain(cycler)
        per_enc_dv = ()

    return MaintenanceOptimResult(
        cycler=cycler,
        t0_sec=float(best_x[0]),
        leg_tofs_days=leg_tofs_days,
        maintenance_dv_kms=float(maintenance_dv),
        per_encounter_dv_kms=per_enc_dv,
        converged=converged,
        a_au=float(a_au),
        e=float(e),
        vinf_kms_at_encounters=vinf_per_enc,
        turn_deficit=turn_deficit,
    )


def optimise_aldrin_maintenance_dv(
    ephem: Ephemeris,
    *,
    t0_guess_sec: float | None = None,
    em_tof_days_guess: float = 146.0,
    me_tof_days_guess: float = 634.0,
    n_starts: int = 5,
    seed: int = 0,
) -> MaintenanceOptimResult:
    """Find the minimum-ΔV periodic Aldrin E→M→E cycler over ``ephem``.

    Thin Aldrin-specific wrapper over :func:`optimise_maintenance_dv`: pins the
    ``E-M-E`` sequence, the ``("E", "M")`` synodic window, the Earth closure
    flyby at the sourced 200 km altitude, and the ``(20, 60)`` day leg jitters.
    The parameterisation and default guesses are unchanged from the original
    Aldrin-only optimiser, so results are bit-for-bit identical.

    Parameters
    ----------
    ephem:
        Planet-state provider. ``Ephemeris("circular")`` is fast and used by
        the test surface; ``Ephemeris("astropy")`` gives real DE440 states.
    t0_guess_sec:
        Earth-departure epoch guess (s). ``None`` derives it from the Aldrin
        seed phase (see :func:`_default_t0_guess`).
    em_tof_days_guess, me_tof_days_guess:
        Initial leg ToFs (days) for the Earth→Mars and Mars→Earth legs.
    n_starts:
        Number of SLSQP polish starts (the seed plus ``n_starts - 1`` jittered
        starts), in addition to the global DE pass.
    seed:
        RNG seed for ``differential_evolution`` and the multi-start jitter, so
        results are reproducible.

    Returns
    -------
    MaintenanceOptimResult
        With the optimised cycler, computed anchors, and the (computed,
        not source-attested) maintenance ΔV.
    """
    if t0_guess_sec is None:
        t0_guess_sec = _default_t0_guess(em_tof_days_guess)

    return optimise_maintenance_dv(
        ["E", "M", "E"],
        ephem,
        t0_guess_sec=t0_guess_sec,
        tof_days_guesses=(em_tof_days_guess, me_tof_days_guess),
        tof_bounds_days=(_TOF1_BOUNDS_DAYS, _TOF2_BOUNDS_DAYS),
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        tof_jitter_half_days=(20.0, 60.0),
        n_starts=n_starts,
        seed=seed,
        seed_cycler_factory=lambda: build_aldrin_seed(Ephemeris("circular")),
    )


__all__ = [
    "FlybyTurnDeficit",
    "MaintenanceOptimResult",
    "idealized_flyby_turn_deficit",
    "optimise_aldrin_maintenance_dv",
    "optimise_maintenance_dv",
]
