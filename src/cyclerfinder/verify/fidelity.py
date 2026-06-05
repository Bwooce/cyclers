"""Axis B — the fidelity-ladder gate (Forge phase 1).

A candidate cycler can be *solved* at several model fidelities. A trustworthy
solution must, across that ladder, either keep a tracked quantity stable
(``PERSISTS``) or move it in a *source-documented* direction/band
(``SHIFTS_DOCUMENTED``). A quantity that moves with no documented reason is the
red flag (``SHIFTS_UNDOCUMENTED``) the gauntlet is built to catch — it is the
cross-fidelity confusion class that produced the S1L1 5.65-vs-4.99 episode.

Two pieces, mirroring the spec §16.7 Axis-B sketch:

* :func:`solve_at_fidelity` — a thin wrapper dispatching a :class:`Cell` to the
  existing machinery per fidelity rung. It does **not** reimplement any physics:

  - ``"circular-coplanar"`` →
    :func:`~cyclerfinder.search.resonant_construct.construct_resonant_cycler`
    from the cell's *sourced* ``(a, e)`` (the resonance-anchored coplanar
    construction);
  - ``"real-de440"`` →
    :func:`~cyclerfinder.search.optimize.optimise_cell_ephemeris` on a DE440
    astropy ephemeris.
  - ``"analytic-ephemeris"`` is a documented **extension point**: no in-house
    analytic-ephemeris backend exists (the Jones AAS 17-577 multisets are an
    external tabulation, and VEM real-ephemeris convergence is M-ED's job), so
    this rung raises :class:`FidelityRungUnavailableError` rather than silently
    standing in for one of the other two.

* :func:`fidelity_persistence` — the classifier. Given a tracked quantity's
  value at each available rung plus a *documented expectation*, it returns a
  frozen :class:`PersistenceReport` labelling the shift.

Golden discipline
------------------
This module computes solutions; it never supplies the EXPECTED side of a golden
check. The documented expectation passed to :func:`fidelity_persistence`
(a tolerance, a direction, or a sourced band) is the caller's responsibility to
source from the catalogue — see ``tests/verify/test_fidelity_gate.py``, where
both the coplanar idealization ToF (146 d) and the STOUR band ([161, 172] d)
are read from the Aldrin row's ``trajectory.segments`` (schema v4.2).
"""

from __future__ import annotations

import enum
from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.provenance import Fidelity, is_fidelity
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.resonant_construct import construct_resonant_cycler
from cyclerfinder.search.sequence import Cell

# ---------------------------------------------------------------------------
# solve_at_fidelity — the ladder wrapper
# ---------------------------------------------------------------------------


class FidelityRungUnavailableError(NotImplementedError):
    """Raised when a :data:`~cyclerfinder.data.provenance.Fidelity` rung has no
    in-house backend wired today.

    Currently only ``"analytic-ephemeris"`` is unavailable — there is no
    in-house analytic/mean-element ephemeris optimiser (the documented middle
    rung of the ladder). Raising rather than approximating keeps the ladder
    honest: a cross-fidelity comparison must never silently substitute a
    neighbouring tier (the S1L1 mismatch class).
    """


@dataclass(frozen=True)
class FidelitySolution:
    """One cell solved at one fidelity rung.

    All fields are COMPUTED outputs of the wrapped machinery — never a golden
    EXPECTED value. ``outbound_tof_days`` and ``vinf_kms`` are the two tracked
    quantities Axis B ladders.

    Attributes
    ----------
    fidelity:
        The rung this solution was produced at.
    cell_id:
        ``cell.id`` echoed for traceability.
    outbound_tof_days:
        Time of flight of the *first* leg (the outbound transit), days.
        For ``circular-coplanar`` this is the resonance construction's first
        inter-crossing ToF; for ``real-de440`` it is the optimised first leg's
        ``t_arrive - t_depart``.
    vinf_kms:
        Per-body V∞ at each encounter, km/s, keyed by body code. For
        ``circular-coplanar`` the construction yields one value per body; for
        ``real-de440`` the per-encounter departure V∞ magnitudes are reduced to
        a per-body value (first encounter of each body).
    converged:
        Whether the underlying solve reports a feasible converged result.
        ``circular-coplanar`` is a closed-form construction and is always
        ``True``; ``real-de440`` follows the optimiser's
        ``converged ∧ constraints_satisfied`` predicate.
    """

    fidelity: Fidelity
    cell_id: str
    outbound_tof_days: float
    vinf_kms: dict[str, float]
    converged: bool


def solve_at_fidelity(
    cell: Cell,
    fidelity: Fidelity,
    *,
    a_au: float | None = None,
    e: float | None = None,
    ephem: Ephemeris | None = None,
    vinf_cap: float = 12.0,
    vinf_targets_kms: dict[str, float] | None = None,
    priority_date_iso: str | None = None,
    n_starts: int = 5,
    seed: int = 0,
) -> FidelitySolution:
    """Solve *cell* at one fidelity rung by dispatching to existing machinery.

    Parameters
    ----------
    cell:
        The discrete structural cell to solve.
    fidelity:
        One of the :data:`~cyclerfinder.data.provenance.Fidelity` tiers.
    a_au, e:
        Sourced spacecraft heliocentric ``(a, e)`` — **required** for the
        ``circular-coplanar`` rung (the resonance construction's inputs). The
        ``real-de440`` rung ignores them (it discovers the geometry).
    ephem:
        Planet-state provider for the ``real-de440`` rung. ``None`` ⇒ an
        ``Ephemeris("astropy")`` is constructed. Ignored by ``circular-coplanar``.
    vinf_cap:
        Hard V∞ ceiling for the ``real-de440`` optimiser, km/s.
    vinf_targets_kms, priority_date_iso:
        Phase-match anchors for the ``real-de440`` epoch resolution — both
        required there (blind discovery is out of scope, per
        :func:`optimise_cell_ephemeris`).
    n_starts, seed:
        Multi-start count and RNG seed for the ``real-de440`` optimiser
        (reproducible).

    Returns
    -------
    FidelitySolution
        Frozen record of the two tracked quantities at this rung.

    Raises
    ------
    FidelityRungUnavailableError
        If ``fidelity == "analytic-ephemeris"`` (no in-house backend).
    ValueError
        If ``fidelity`` is not a known tier, or required inputs are missing for
        the chosen rung.
    """
    if not is_fidelity(fidelity):
        raise ValueError(f"unknown fidelity tier {fidelity!r}")

    if fidelity == "circular-coplanar":
        return _solve_coplanar(cell, a_au=a_au, e=e)
    if fidelity == "analytic-ephemeris":
        raise FidelityRungUnavailableError(
            "no in-house analytic-ephemeris backend: the middle ladder rung is "
            "an explicit extension point (Jones AAS 17-577 multisets are an "
            "external tabulation; VEM real-ephemeris convergence is M-ED's job)",
        )
    # real-de440
    return _solve_real_de440(
        cell,
        ephem=ephem,
        vinf_cap=vinf_cap,
        vinf_targets_kms=vinf_targets_kms,
        priority_date_iso=priority_date_iso,
        n_starts=n_starts,
        seed=seed,
    )


def _solve_coplanar(cell: Cell, *, a_au: float | None, e: float | None) -> FidelitySolution:
    """Circular-coplanar rung via the resonance-anchored construction."""
    if a_au is None or e is None:
        raise ValueError(
            "circular-coplanar rung requires sourced (a_au, e) — they are the "
            "resonance construction's inputs",
        )
    # Encounter bodies along the trajectory, de-duplicating the closing body so
    # the construction sees each distinct body once (E-M-E -> ("E", "M")).
    bodies: list[str] = []
    for b in cell.sequence:
        if b not in bodies:
            bodies.append(b)
    rc = construct_resonant_cycler(a_au, e, bodies=tuple(bodies))
    # Outbound leg ToF: first inter-crossing in body order.
    first_from, first_to = cell.sequence[0], cell.sequence[1]
    outbound_tof_days = rc.leg_tofs_days[f"{first_from}->{first_to}"]
    return FidelitySolution(
        fidelity="circular-coplanar",
        cell_id=cell.id,
        outbound_tof_days=outbound_tof_days,
        vinf_kms=dict(rc.vinf_kms),
        converged=True,
    )


def _solve_real_de440(
    cell: Cell,
    *,
    ephem: Ephemeris | None,
    vinf_cap: float,
    vinf_targets_kms: dict[str, float] | None,
    priority_date_iso: str | None,
    n_starts: int,
    seed: int,
) -> FidelitySolution:
    """Real-DE440 rung via the general ephemeris cell optimiser."""
    if ephem is None:
        ephem = Ephemeris(model="astropy")
    result = optimise_cell_ephemeris(
        cell,
        ephem,
        vinf_cap=vinf_cap,
        priority_date_iso=priority_date_iso,
        vinf_targets_kms=vinf_targets_kms,
        n_starts=n_starts,
        seed=seed,
    )
    cyc = result.best_cycler
    first_leg = cyc.legs[0]
    outbound_tof_days = (first_leg.t_arrive - first_leg.t_depart) / SECONDS_PER_DAY
    # Per-body V∞: first encounter of each body, departure-side magnitude.
    vinf_kms: dict[str, float] = {}
    for enc in cyc.encounters:
        if enc.body not in vinf_kms:
            vinf_kms[enc.body] = float(np.linalg.norm(enc.vinf_out))
    return FidelitySolution(
        fidelity="real-de440",
        cell_id=cell.id,
        outbound_tof_days=outbound_tof_days,
        vinf_kms=vinf_kms,
        converged=bool(result.converged and result.constraints_satisfied),
    )


# ---------------------------------------------------------------------------
# fidelity_persistence — the classifier
# ---------------------------------------------------------------------------


class PersistenceClass(enum.Enum):
    """How a tracked quantity behaves across the fidelity ladder.

    * :attr:`PERSISTS` — the value stays within ``abs_tol`` across all rungs.
      A fidelity-robust quantity; no documented shift is needed.
    * :attr:`SHIFTS_DOCUMENTED` — the value moves *beyond* tolerance but in the
      direction (or into the band) the source documents. The expected,
      explained fidelity sensitivity (e.g. the Aldrin outbound ToF rising from
      the coplanar idealization toward the STOUR ephemeris band).
    * :attr:`SHIFTS_UNDOCUMENTED` — the value moves beyond tolerance with no
      documented reason (wrong direction, or outside any documented band). The
      red flag: a cross-fidelity discrepancy that has not been explained, the
      class the gauntlet exists to surface.
    """

    PERSISTS = "persists"
    SHIFTS_DOCUMENTED = "shifts-documented"
    SHIFTS_UNDOCUMENTED = "shifts-undocumented"


@dataclass(frozen=True)
class PersistenceReport:
    """Result of :func:`fidelity_persistence` for ONE tracked quantity.

    Attributes
    ----------
    classification:
        The :class:`PersistenceClass` label.
    quantity:
        Human-readable name of the tracked quantity (e.g. ``"outbound_tof_days"``).
    low_value, high_value:
        The quantity's value at the lower-fidelity and higher-fidelity rung
        respectively (the two rungs being compared).
    delta:
        ``high_value - low_value`` (signed). Positive ⇒ the higher-fidelity
        value is larger.
    within_tol:
        Whether ``abs(delta) <= abs_tol`` (i.e. the quantity PERSISTS).
    """

    classification: PersistenceClass
    quantity: str
    low_value: float
    high_value: float
    delta: float
    within_tol: bool


def fidelity_persistence(
    quantity: str,
    low_value: float,
    high_value: float,
    *,
    abs_tol: float,
    expected_direction: int | None = None,
    documented_band: tuple[float, float] | None = None,
) -> PersistenceReport:
    """Classify a tracked quantity's behaviour from low → high fidelity.

    Pure function (no catalogue read, no physics): it classifies a *shift*
    against a caller-supplied *documented expectation*. The expectation
    (``expected_direction`` and/or ``documented_band``) is the caller's; the
    caller must source it from the catalogue (golden discipline) — this function
    only compares.

    Parameters
    ----------
    quantity:
        Name of the tracked quantity (carried into the report).
    low_value:
        The quantity at the lower-fidelity rung.
    high_value:
        The quantity at the higher-fidelity rung.
    abs_tol:
        Absolute tolerance below which the quantity is deemed to PERSIST.
    expected_direction:
        ``+1`` ⇒ the higher-fidelity value is documented to *increase*;
        ``-1`` ⇒ documented to *decrease*; ``None`` ⇒ no direction is
        documented (then only ``documented_band`` can justify a shift).
    documented_band:
        ``(lo, hi)`` sourced band the higher-fidelity value is documented to
        move *toward or into*. A shift is documented if the high value either
        lands inside the band, or moves *closer* to the band than the low value
        did (toward-the-band counts, per the plan's "in/toward" wording).
        ``None`` ⇒ no band documented.

    Returns
    -------
    PersistenceReport
        Frozen classification.

    Notes
    -----
    Classification order:

    1. If ``abs(delta) <= abs_tol`` ⇒ :attr:`PersistenceClass.PERSISTS`.
    2. Else if the shift agrees with ``expected_direction`` *or* moves
       toward/into ``documented_band`` ⇒
       :attr:`PersistenceClass.SHIFTS_DOCUMENTED`.
    3. Else ⇒ :attr:`PersistenceClass.SHIFTS_UNDOCUMENTED`.

    With *neither* ``expected_direction`` nor ``documented_band`` supplied, any
    beyond-tolerance shift is necessarily ``SHIFTS_UNDOCUMENTED`` — there is no
    documentation to justify it (the conservative default).
    """
    delta = high_value - low_value
    within_tol = abs(delta) <= abs_tol
    if within_tol:
        return PersistenceReport(
            classification=PersistenceClass.PERSISTS,
            quantity=quantity,
            low_value=low_value,
            high_value=high_value,
            delta=delta,
            within_tol=True,
        )

    documented = False
    if expected_direction is not None:
        sign = 1 if delta > 0 else -1
        if sign == (1 if expected_direction > 0 else -1):
            documented = True
    if not documented and documented_band is not None:
        documented = _moves_toward_band(low_value, high_value, documented_band)

    classification = (
        PersistenceClass.SHIFTS_DOCUMENTED if documented else PersistenceClass.SHIFTS_UNDOCUMENTED
    )
    return PersistenceReport(
        classification=classification,
        quantity=quantity,
        low_value=low_value,
        high_value=high_value,
        delta=delta,
        within_tol=False,
    )


def _moves_toward_band(low_value: float, high_value: float, band: tuple[float, float]) -> bool:
    """True iff ``high_value`` is inside *band* or closer to it than ``low_value``."""
    lo, hi = (band[0], band[1]) if band[0] <= band[1] else (band[1], band[0])

    def dist_to_band(x: float) -> float:
        if x < lo:
            return lo - x
        if x > hi:
            return x - hi
        return 0.0

    d_high = dist_to_band(high_value)
    if d_high == 0.0:
        return True
    return d_high < dist_to_band(low_value)


__all__ = [
    "FidelityRungUnavailableError",
    "FidelitySolution",
    "PersistenceClass",
    "PersistenceReport",
    "fidelity_persistence",
    "solve_at_fidelity",
]
