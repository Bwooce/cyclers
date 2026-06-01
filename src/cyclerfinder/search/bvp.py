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

Status
------
**Not yet implemented.** :func:`solve_powered_periodic_cycler` raises
:class:`NotImplementedError`. The positive closure gate
(:func:`test_aldrin_powered_cycler_closes_on_de440`) is ``xfail(strict=True)``
against this stub and flips to XPASS — failing CI and prompting removal of the
marker — the moment the solver lands and closes the Aldrin cycler.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler


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
    NotImplementedError
        Always, until the BVP solver is implemented (task #71).
    """
    del catalogue_entry, ephem, t_start_sec, signature_priority_date
    raise NotImplementedError(
        "solve_powered_periodic_cycler (Phase B / task #71) is not yet "
        "implemented. It must solve the flyby-constrained periodic BVP "
        "(|V_inf| continuity at each flyby + planets realigning at the synodic "
        "period + maintenance maneuver paying the turn deficit, minimising "
        "summed Delta V). Until then, a powered cycler cannot be closed on real "
        "ephemeris; see test_aldrin_ballistic_closure_fails_because_powered for "
        "the ballistic negative result."
    )
