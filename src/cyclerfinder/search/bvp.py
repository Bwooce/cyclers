"""Phase B — flyby-constrained periodic boundary-value solver for *powered*
cyclers (task #71).

Motivation
----------
The classic Aldrin E-M cycler is **powered**: its Earth flyby needs an ~84°
heliocentric turn but can ballistically deliver only ~72° at a 200 km flyby
(McConaghy 2002 Table 4 "1L1"; see
:func:`tests/verify/test_real_closure.py::test_aldrin_powered_turn_deficit_gate`).
Because of that turn deficit no purely ballistic Lambert chain closes the
cycler on real ephemeris — empirically it drifts ~10^8 km over two cycles
(see :func:`test_aldrin_ballistic_closure_fails_because_powered`). Closing a
powered cycler requires solving a periodic boundary-value problem (BVP):

    find the launch epoch and per-leg times-of-flight such that

      * the spacecraft heliocentric state repeats after one cycler period
        with both planets realigned at the synodic repeat period,
      * ``|V_inf|`` is preserved across each gravity-assist flyby, and
      * the turn the geometry demands beyond the achievable ballistic cone is
        paid as a single periapsis maintenance maneuver,

    minimising the summed maintenance ``Delta V`` over the cycle.

This is the ephemeris-mode construction the M6b plan (§3.1 Path A) deferred and
that the closure xfail reasons point to. It is the *positive* counterpart to
the ballistic-closure negative test.

Implementation
--------------
The solver resolves a phase-matched launch epoch from the catalogue signature
(or uses a caller-supplied epoch) and then delegates the periodic slice to
:func:`cyclerfinder.search.maintain.optimise_aldrin_maintenance_dv`, which
finds the minimum-ΔV E→M→E family at that launch phase on the real ephemeris
and charges the maintenance ΔV as the sourced return-flyby turn deficit
(Earth for Aldrin; ≈84° required vs ≈72° achievable). The returned cycler is
the genuine *powered* periodic cycler with a strictly positive maintenance ΔV.

Closure caveat
--------------
Producing the cycler is **not** the same as the rotating-frame drift metric
falling under :data:`~cyclerfinder.verify.real_closure.REAL_DRIFT_TOLERANCE_KM`
(200,000 km). For the k=1 Aldrin cycler that bound is physically unreachable on
DE440: the drift propagator pins each leg-start to the *real* planet position
at the lap-shifted epoch, and Mars's heliocentric radius breathes ≈0.117 AU
(≈1.75e7 km) per 2.135 yr cycle because the cycler period is not commensurate
with Mars's 1.881 yr orbit. The empirical drift floor is ≈4.14e8 km (≈2072x
the tolerance, #134) regardless of the maneuver — the maneuver shapes velocity, not
where Mars is. This is exactly why the real Aldrin cycler needs a per-cycle
retargeting maneuver; the 200,000 km rotating-frame-repeat criterion is a
circular-ephemeris idealisation that eccentric Mars cannot satisfy.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv
from cyclerfinder.search.phase_match import phase_signature_from_catalogue_entry


@dataclass(frozen=True)
class PoweredCyclerSolution:
    """Result of :func:`solve_powered_periodic_cycler`.

    Attributes
    ----------
    cycler:
        The maneuvered, flyby-constrained periodic :class:`Cycler` built on the
        real ephemeris — ready to feed to
        :func:`cyclerfinder.verify.real_closure.verify_real_closure`.
    t_start_sec:
        Inertial-frame launch epoch (seconds since J2000) the solution was
        solved at.
    total_maintenance_dv_kms:
        Summed maintenance ``Delta V`` (km/s) over one cycle — the quantity the
        BVP minimises. Strictly positive for a genuinely powered cycler.
    per_encounter_dv_kms:
        ``(body_code, dv_kms)`` pairs giving the maneuver paid at each
        encounter. Sums to ``total_maintenance_dv_kms``.
    """

    cycler: Cycler
    t_start_sec: float
    total_maintenance_dv_kms: float
    per_encounter_dv_kms: tuple[tuple[str, float], ...]


def solve_powered_periodic_cycler(
    catalogue_entry: dict[str, object],
    ephem: Ephemeris,
    *,
    t_start_sec: float | None = None,
    signature_priority_date: datetime | None = None,
) -> PoweredCyclerSolution:
    """Solve the flyby-constrained periodic BVP for a powered cycler.

    Parameters
    ----------
    catalogue_entry:
        A catalogue YAML entry (e.g. ``aldrin-classic-em-k1-outbound``) whose
        ``trajectory_regime`` is ``"powered"``.
    ephem:
        Ephemeris backend supplying real planet states (``Ephemeris("astropy")``
        for DE440).
    t_start_sec:
        Optional fixed launch epoch (seconds since J2000). When ``None`` the
        solver resolves a phase-matched epoch near ``signature_priority_date``.
    signature_priority_date:
        Priority date used to resolve the launch epoch when ``t_start_sec`` is
        ``None``.

    Returns
    -------
    PoweredCyclerSolution

    Raises
    ------
    ValueError
        If no launch epoch is supplied and none can be resolved (neither
        ``t_start_sec`` nor a ``signature_priority_date`` that yields a
        phase-matched window).
    """
    if t_start_sec is None:
        if signature_priority_date is None:
            raise ValueError(
                "solve_powered_periodic_cycler needs a launch epoch: pass "
                "t_start_sec, or signature_priority_date so the epoch can be "
                "resolved from the catalogue signature."
            )
        # Local import: verify.real_closure depends on search machinery, so a
        # module-level import here would risk a construction-time cycle.
        from cyclerfinder.verify.real_closure import _resolve_real_t_start

        signature = phase_signature_from_catalogue_entry(catalogue_entry)
        resolved = _resolve_real_t_start(signature, ephem, signature_priority_date)
        if resolved is None:
            raise ValueError(
                "could not resolve a phase-matched launch epoch within the "
                f"search window of {signature_priority_date!r}; supply "
                "t_start_sec explicitly."
            )
        t_start_sec = resolved

    em_guess, me_guess = _leg_tof_guesses(catalogue_entry)
    result = optimise_aldrin_maintenance_dv(
        ephem,
        t0_guess_sec=t_start_sec,
        em_tof_days_guess=em_guess,
        me_tof_days_guess=me_guess,
    )

    return PoweredCyclerSolution(
        cycler=result.cycler,
        t_start_sec=result.t0_sec,
        total_maintenance_dv_kms=result.maintenance_dv_kms,
        per_encounter_dv_kms=result.per_encounter_dv_kms,
    )


def _leg_tof_guesses(catalogue_entry: dict[str, object]) -> tuple[float, float]:
    """Initial ``(E→M, M→E)`` leg ToF guesses (days) for the optimiser.

    Reads the catalogue legs' ``tof_days`` when present, falling back to the
    classic Aldrin 146 d / 634 d when a leg omits it. The optimiser leaves both
    leg ToFs fully free, so these only seed the multi-start polish.
    """
    em_default, me_default = 146.0, 634.0
    legs = catalogue_entry.get("legs")
    if not isinstance(legs, list) or len(legs) < 2:
        return em_default, me_default

    def _tof(leg: object, default: float) -> float:
        if isinstance(leg, dict):
            value = leg.get("tof_days")
            if isinstance(value, (int, float)) and value > 0:
                return float(value)
        return default

    return _tof(legs[0], em_default), _tof(legs[1], me_default)
