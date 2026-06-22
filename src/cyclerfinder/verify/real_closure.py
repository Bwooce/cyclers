"""Real-ephemeris closure verification (M6b).

Spec references
---------------
* §12.1 (idealised → ephemeris bridge — :func:`find_real_windows` is
  the shipping entry point M6b composes).
* §12.2 (three-representation framework — M6b promotes a catalogue
  entry from V1 idealised reproduction to V2-real ephemeris instance,
  it does not flatten one representation into another).
* §14 V2-real (multi-cycle real-ephemeris periodicity gate — the
  binding gate M6b stands up).
* §16.1 (catalogue schema v2 fields ``validation.gates.V2.*`` that
  M7's writer fills from :class:`RealClosureResult`).

Purpose
-------
Given a catalogue entry (V1 reproduction passed M6a's
:func:`verify_long_term_stability` on circular ephemeris), this module
verifies the corresponding **real-ephemeris instance** closes on JPL
DE440 over ``n_cycles`` consecutive periods within
:data:`REAL_DRIFT_TOLERANCE_KM` of geometric drift in the dynamic
rotating frame.

Architectural template
----------------------
Pascarella et al. 2024 (Solar System Pony Express) separates
patched-conic, medium-fidelity, and high-fidelity GMAT stages. M6b is
the **medium-fidelity output** stage: real-ephemeris ballistic Lambert
closure without TCMs. TCM ΔV budgeting (the operational maintenance
cost) is V3 of the §14 gauntlet and lands in M7.

Composition map
---------------
M6b deliberately reuses M6a / M6-slice infrastructure rather than
duplicating it:

* :func:`verify_long_term_stability` (M6a) — the propagator +
  drift-measurement core. :func:`verify_real_closure` delegates to
  this; the composition is asserted by
  ``test_real_closure_uses_m6a_machinery``.
* :func:`find_real_windows` (M6 slice) — picks the real launch epoch
  that matches the catalogue's signature.
* :func:`lambert` — solves each leg's velocity vectors over real
  planet positions, single- or multi-rev per the leg's ``n_revs`` and
  ``branch``.

Plan: ``docs/phases/m6b-real-ephemeris-closure/plan.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from typing import Any, Final

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import (
    DAYS_PER_JULIAN_YEAR,
    MU_SUN_KM3_S2,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)
from cyclerfinder.data.catalog import _segments_as_legs
from cyclerfinder.model.cycler import Cycler, Encounter, Leg
from cyclerfinder.search.phase_match import (
    LaunchWindow,
    PhaseSignature,
    find_candidate_windows,
    leg_duration_seeds,
    phase_signature_from_catalogue_entry,
)
from cyclerfinder.verify.propagate import (
    StabilityReport,
    verify_long_term_stability,
)

REAL_DRIFT_TOLERANCE_KM: Final[float] = 200_000.0
"""Maximum permissible cycle-to-cycle drift on real DE440 ephemeris
over ``N >= 2`` cycles. Derivation in plan §4.3:

* 4x M6a's idealised 50,000 km tolerance to absorb real eccentricity
  breathing on E-M trajectories (Earth e ≈ 0.017, Mars e ≈ 0.093).
* Calibrated against the Aldrin classic k=1 cycler: Pascarella 2024
  and Russell 2004 report cycle-to-cycle drift on real ephemeris in
  the 100,000-200,000 km range before TCMs.
* Tight enough to reject the spec §10 degenerate-closure basin
  (V∞ > ~11 km/s), which produces AU/cycle drift.

M7 may tighten this if its TCM-budget computation reports post-TCM
drift consistently below 20,000 km."""

N_CYCLES_DEFAULT: Final[int] = 2
"""Default cycle count for the M6b binding gate. Trade-off (plan §4.3):

* ``n_cycles=2`` covers ≈ 4.3 yr for k=1 cyclers, ≈ 8.5 yr for k=2.
* ``n_cycles=3`` would double CI runtime per gate test.
* ``n_cycles=5`` (spec §12(a) horizon) is M7's batch concern.
"""

_J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
_DAY_S = SECONDS_PER_DAY
_YEAR_S = DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY


# ---------------------------------------------------------------------------
# Exception types
# ---------------------------------------------------------------------------


class RealClosureConstructionError(Exception):
    """Raised when single-rev Lambert fails on a real-ephemeris leg.

    Carries the leg index and the underlying Lambert exception so
    diagnostic output can show the geometry that broke.
    """

    def __init__(
        self,
        catalogue_id: str | None,
        leg_index: int,
        cause: Exception,
    ) -> None:
        self.catalogue_id = catalogue_id
        self.leg_index = leg_index
        self.cause = cause
        super().__init__(
            f"catalogue entry {catalogue_id!r} leg {leg_index} Lambert "
            f"construction failed: {type(cause).__name__}: {cause}"
        )


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RealClosureResult:
    """Result of :func:`verify_real_closure`.

    All fields are immutable. The ``horizon_tcm_mps`` and
    ``per_cycle_tcm_mps`` fields are zeros at M6b; M7 populates them
    when the catalogue ingest writer consumes this result.

    Spec references: §12.1, §12.2, §14 V2-real, §16.1
    ``validation.gates.V2.*``.

    Attributes
    ----------
    cycler_id:
        Catalogue entry id (e.g. ``"aldrin-classic-em-k1-outbound"``)
        if the result was produced from a catalogue cycler, else
        ``None``.
    n_cycles_propagated:
        Number of cycles the propagation actually completed. Equals
        the ``n_cycles`` argument unless an early-termination tripped
        (multi-rev Lambert blocker, propagator divergence, no real
        launch window).
    max_drift_km:
        Maximum consecutive-cycle-pair drift across the
        ``range(n_cycles - 1)`` pairs in the dynamic rotating frame.
        Sourced from :attr:`StabilityReport.max_drift_km`.
    per_cycle_drift_km:
        Cumulative drift at each cycle boundary, sourced from
        :attr:`StabilityReport.per_lap_drift_km`.
    per_encounter_vinf_mismatch_kms:
        At each interior encounter, ``||V_inf_in| - |V_inf_out||`` in
        km/s — diagnostic measure of how far the real-ephemeris
        construction departs from the published ballistic ideal.
        Empty tuple if construction failed or the cycler has no
        interior encounter (e.g. Aldrin).
    closes:
        ``max_drift_km < REAL_DRIFT_TOLERANCE_KM``. The headline
        boolean for spec §14 V2-real (the M6b binding gate).
    v3_status:
        One of ``"v3-real-closure-pass"``, ``"v3-real-closure-fail"``,
        ``"v3-no-real-window"``, ``"v3-construction-error"``. M7's
        catalogue writer records this in
        ``validation.gates.V2.status``.
    horizon_tcm_mps:
        **M6b: 0.0**. M7 populates with the summed TCM ΔV over the
        horizon (m/s). Locked here so M7 does not reshape the
        dataclass.
    per_cycle_tcm_mps:
        **M6b: zero-tuple of length** ``n_cycles_propagated``. M7
        populates per-cycle TCM ΔV (m/s).
    frame_used:
        ``"dynamic"`` always at M6b — the M6a-locked dynamic frame.
        Field exists for forward compatibility with M7 batch runs.
    t_start_sec:
        The inertial-frame launch epoch (seconds since J2000) the
        result was computed against. Either passed in by the caller
        or derived from :func:`_resolve_real_t_start`.
    """

    cycler_id: str | None
    n_cycles_propagated: int
    max_drift_km: float
    per_cycle_drift_km: tuple[float, ...]
    per_encounter_vinf_mismatch_kms: tuple[float, ...]
    closes: bool
    v3_status: str
    horizon_tcm_mps: float
    per_cycle_tcm_mps: tuple[float, ...]
    frame_used: str
    t_start_sec: float | None


# ---------------------------------------------------------------------------
# EXPECTED_SKIPS registry
# ---------------------------------------------------------------------------


EXPECTED_SKIPS: Final[dict[str, str]] = {
    "s1l1-2syn-em-cpom": (
        "wrong topology, not missing data: the M->E direct return leg is "
        "undefined by construction. In the McConaghy/Longuski/Byrnes "
        "nomenclature the S and L labels denote consecutive Earth-to-Earth "
        "resonant intervals (S1 short, L1 long), not Earth<->Mars transits; "
        "Mars is a flyby target of opportunity on the outbound (~154-d) arc. "
        "A single S1L1 vehicle provides only one-way fast transit, so there "
        "is no tabulated direct Mars->Earth crewed return ToF to source — the "
        "return requires a mirrored companion cycler (the conjugate L1S1). "
        "Closing s1l1 as an E->M->E round trip therefore models a leg that "
        "does not exist; reframing the trajectory as outbound E->M plus the "
        "S1/L1 Earth-to-Earth resonant intervals is the correct fix and is "
        "tracked in the entry's data_gaps. See the catalogue entry."
    ),
    "mcconaghy-2006-em-k2": (
        "incomplete leg data: catalogue entry has only 2 legs (E-M, M-E) "
        "totalling 306 d but advertises a 2-synodic period of 4.27 yr; the "
        "intermediate Earth/Mars encounters and their legs are not "
        "tabulated in the published abstract source. Full-cycle "
        "reconstruction is M7's catalogue-completion concern."
    ),
    "russell-ocampo-2.1.1+2-case2": (
        "incomplete leg data: catalogue entry has 1 leg (E-M) of 207 d "
        "but advertises a 2-synodic 4.27 yr period; the return / "
        "intermediate-encounter legs are not in Russell 2004 Table 3.4 "
        "(only AR/TR/aphelion-ratio summarised). Full-cycle "
        "reconstruction is M7's concern."
    ),
    "russell-ocampo-2.5.1+0": (
        "incomplete leg data: catalogue entry has 1 leg (E-M) of 94 d "
        "but advertises a 2-synodic 4.27 yr period; same Russell-table "
        "summarisation as 2.1.1+2-case2. Full-cycle reconstruction is "
        "M7's concern."
    ),
    "jones-2017-vem-triple-family": ("VEM 3-body real-closure is M8 scope (VEM campaign + viz)."),
    "vem-emeeve-3syn": ("VEM 3-body real-closure is M8 scope."),
}
"""Registry of catalogue ids M6b deliberately does not verify.

Why this exists: M6b's regression set runs over a hand-picked subset
of literature-anchored cyclers. The Aldrin outbound/inbound entries
carry complete 2-leg cycle data and are the binding gate. The other
regression candidates have only the outbound leg published, or (S1L1
CPOM) have no direct M->E return leg to verify at all — the S/L labels
denote Earth-to-Earth resonant intervals, not Earth<->Mars transits, so
the round-trip topology must be re-modelled rather than back-filled.
Full-cycle reconstruction is catalogue-completion work that belongs to
M7 alongside the catalogue writer. The non-Sun primary cases are filtered at load time by
:func:`tests.data._catalogue_loader_m6b.load_m6b_entries`; the VEM
entries here document the 3-body deferral. Every contributor reading
this list knows what is outstanding."""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _dt_to_t_sec(dt: datetime) -> float:
    """Convert a UTC-aware datetime to seconds since J2000."""
    if dt.tzinfo is None:
        raise ValueError(f"datetime must be timezone-aware; got naive {dt!r}")
    return (dt - _J2000_EPOCH).total_seconds()


def _parse_priority_date(value: Any) -> datetime | None:
    """Parse a catalogue ``priority_date`` field to a UTC-aware datetime.

    Accepts ISO-8601 ``"YYYY-MM-DD"`` strings, ``datetime`` instances,
    and ``date`` instances. Returns ``None`` for ``None`` or unparseable
    values; callers default to "no priority date" in that case.
    """
    if value is None:
        return None
    if isinstance(value, datetime):
        if value.tzinfo is None:
            return value.replace(tzinfo=UTC)
        return value
    if isinstance(value, str):
        try:
            parsed = datetime.fromisoformat(value)
        except ValueError:
            return None
        if parsed.tzinfo is None:
            return parsed.replace(tzinfo=UTC)
        return parsed
    # ``datetime.date`` (non-datetime) — convert via isoformat round-trip.
    iso = getattr(value, "isoformat", None)
    if callable(iso):
        try:
            parsed = datetime.fromisoformat(iso())
        except ValueError:
            return None
        if parsed.tzinfo is None:
            return parsed.replace(tzinfo=UTC)
        return parsed
    return None


def _resolve_real_t_start(
    signature: PhaseSignature,
    ephem: Ephemeris,
    priority_date: datetime,
    *,
    window_years: float = 10.0,
    n_candidates: int = 5,
    mismatch_cap_kms: float = 20.0,
) -> float | None:
    """Pick a real-ephemeris launch epoch matching the signature.

    Searches ``priority_date ± window_years`` for the best
    :func:`find_real_windows` match; returns the lowest-mismatch
    window's J2000-relative seconds-since-epoch or ``None`` if no
    window beats ``mismatch_cap_kms``.

    Parameters
    ----------
    signature:
        The cycler's :class:`PhaseSignature` (geometric fingerprint).
    ephem:
        Ephemeris backend (must be ``"astropy"`` for real dates; the
        underlying :func:`find_real_windows` does not enforce but the
        result is meaningless on ``"circular"``).
    priority_date:
        Centre of the search range — the catalogue's literature epoch.
    window_years:
        Half-width of the search range, years. Default 10 (broad enough
        to catch the next favourable synodic alignment even when the
        priority date sits between alignments).
    n_candidates:
        Maximum windows to retrieve from :func:`find_real_windows`.
    mismatch_cap_kms:
        Reject windows whose summed |V∞| mismatch exceeds this.
        Default 20 km/s is generous; the M6a xfail test uses the same
        cap.
    """
    delta = timedelta(days=window_years * DAYS_PER_JULIAN_YEAR)
    start = priority_date - delta
    end = priority_date + delta

    # STAGE 3: fan the literature signature into a family of asymmetric
    # leg-duration seeds, then rank the merged candidate pool by V_inf
    # mismatch rather than calendar proximity. The old proximity tie-break
    # biased the resolver toward a near-date window even when a far-date
    # window matched the signature far better — the degenerate-basin bug.
    period_s = sum(signature.leg_durations_s)
    seeds = leg_duration_seeds(
        bodies=signature.bodies,
        primary_leg_durations_s=signature.leg_durations_s,
        vinf_target_kms=signature.vinf_target_kms,
        period_s=period_s,
    )
    windows = find_candidate_windows(
        seeds,
        ephem,
        (start, end),
        n=n_candidates * len(seeds),
        mismatch_cap_kms=mismatch_cap_kms,
    )
    if not windows:
        return None
    # Lowest mismatch wins; find_candidate_windows already sorts ascending.
    return _dt_to_t_sec(windows[0].departure_date)


def _full_body_chain(catalogue_entry: dict[str, Any]) -> tuple[str, ...]:
    """Derive the per-leg body chain from a catalogue entry's legs.

    Returns a tuple of length ``len(legs) + 1`` where consecutive
    elements equal each leg's ``from``/``to`` bodies. Raises
    :class:`ValueError` if the legs are inconsistent (e.g. leg ``j+1``
    does not start at leg ``j``'s ``to``).

    The catalogue's top-level ``bodies`` field is the deduplicated set
    (e.g. ``["E", "M"]`` for Aldrin's 2-leg ``E→M→E`` chain), so it
    cannot be indexed directly by leg.
    """
    legs = _segments_as_legs(catalogue_entry)
    if not legs:
        raise ValueError(
            f"catalogue entry {catalogue_entry.get('id')!r} has no legs; cannot derive body chain."
        )
    chain: list[str] = [str(legs[0]["from"])]
    for j, leg in enumerate(legs):
        leg_from = str(leg.get("from", ""))
        leg_to = str(leg.get("to", ""))
        if not leg_from or not leg_to:
            raise ValueError(
                f"catalogue entry {catalogue_entry.get('id')!r} leg {j} missing from/to: {leg!r}"
            )
        if chain[-1] != leg_from:
            raise ValueError(
                f"catalogue entry {catalogue_entry.get('id')!r} leg {j} "
                f"starts at {leg_from!r} but previous leg arrived at "
                f"{chain[-1]!r}; legs are not chained."
            )
        chain.append(leg_to)
    return tuple(chain)


def _check_vinf_continuity(
    cycler: Cycler,
    ephem: Ephemeris,
) -> tuple[float, ...]:
    """Per-interior-encounter ``||V∞_in| - |V∞_out||`` magnitude residual.

    Pure ballistic flybys preserve |V∞|; M6b reports the mismatch as
    a diagnostic, not a gate (drift is the gate). Returns an empty
    tuple when the cycler has no interior encounter (e.g. Aldrin's
    2-encounter ``E-M`` chain).

    Module-internal but exercised by ``test_check_vinf_continuity_*``.
    """
    encounters = cycler.encounters
    if len(encounters) <= 2:
        return ()
    out: list[float] = []
    for enc in encounters[1:-1]:
        v_in = float(np.linalg.norm(enc.vinf_in))
        v_out = float(np.linalg.norm(enc.vinf_out))
        out.append(abs(v_in - v_out))
    return tuple(out)


def construct_real_ephemeris_cycler(
    catalogue_entry: dict[str, Any],
    ephem: Ephemeris,
    launch_window: LaunchWindow | float,
) -> Cycler:
    """Build a real-ephemeris :class:`Cycler` from a catalogue entry.

    Algorithm (plan §3.1.1): derive the per-leg body chain from
    ``legs[].from/to``; place the first encounter at ``launch_window``;
    cumulate leg ToFs to set successor encounter epochs; read planet
    positions from ``ephem`` at each epoch; Lambert-solve each leg to
    populate ``Leg`` and ``Encounter`` velocity vectors.

    Parameters
    ----------
    catalogue_entry:
        Raw catalogue YAML dict with at least ``id``, ``bodies``,
        ``legs[].{from,to,tof_days,n_revs}`` and ``period.years``.
    ephem:
        Real-ephemeris backend (typically ``Ephemeris("astropy")``).
    launch_window:
        Either a :class:`LaunchWindow` (its ``departure_date`` is the
        launch epoch) or a raw J2000-relative ``t_start`` in seconds.

    Returns
    -------
    Cycler
        Frozen cycler whose ``period`` is set from
        ``catalogue_entry["period"]["years"] * year_seconds``,
        whose ``encounters`` are placed at the real planet states at
        each cumulative-ToF epoch, and whose ``legs`` are Lambert
        solutions over those real positions. A leg's ``n_revs`` selects
        the revolution count; its optional ``branch`` ("low"/"high")
        picks the multi-rev branch, defaulting to "low" when omitted.

    Raises
    ------
    RealClosureConstructionError
        If Lambert fails on any leg (geometry singular, non-convergent,
        or the requested ``n_revs``/``branch`` solution does not exist --
        e.g. a multi-rev leg whose time of flight is below ``t_min(n)``).
    ValueError
        If the entry's legs are missing or inconsistent.
    """
    cat_id = catalogue_entry.get("id")
    chain = _full_body_chain(catalogue_entry)
    legs_meta = _segments_as_legs(catalogue_entry)

    if isinstance(launch_window, LaunchWindow):
        t_start_sec = _dt_to_t_sec(launch_window.departure_date)
    else:
        t_start_sec = float(launch_window)

    # Encounter epochs from cumulative ToFs.
    encounter_times: list[float] = [t_start_sec]
    for j, leg in enumerate(legs_meta):
        tof_days = leg.get("tof_days")
        if tof_days is None:
            raise ValueError(f"catalogue entry {cat_id!r} leg {j} has null tof_days")
        encounter_times.append(encounter_times[-1] + float(tof_days) * _DAY_S)

    # Per-encounter planet states from real ephemeris.
    planet_states: list[tuple[NDArray[np.float64], NDArray[np.float64]]] = []
    for i, body in enumerate(chain):
        r, v = ephem.state(body, encounter_times[i])
        planet_states.append(
            (
                np.asarray(r, dtype=np.float64),
                np.asarray(v, dtype=np.float64),
            )
        )

    # Per-leg Lambert solutions. Each leg carries its revolution count
    # (``n_revs``) and an optional ``branch`` ("low"/"high"); the catalogue
    # omits ``branch``, so multi-rev legs default to "low".
    leg_vels: list[tuple[NDArray[np.float64], NDArray[np.float64], int, str]] = []
    for j, leg in enumerate(legs_meta):
        n_revs = int(leg.get("n_revs", 0) or 0)
        requested_branch = leg.get("branch") or ("single" if n_revs == 0 else "low")
        r1, _v1_planet = planet_states[j]
        r2, _v2_planet = planet_states[j + 1]
        tof = encounter_times[j + 1] - encounter_times[j]
        try:
            sols = lambert(r1, r2, tof, mu=MU_SUN_KM3_S2, max_revs=n_revs)
        except (LambertConvergenceError, LambertGeometryError) as exc:
            raise RealClosureConstructionError(
                catalogue_id=cat_id if cat_id is None else str(cat_id),
                leg_index=j,
                cause=exc,
            ) from exc
        chosen = next(
            (s for s in sols if s.n_revs == n_revs and s.branch == requested_branch),
            None,
        )
        if chosen is None:
            raise RealClosureConstructionError(
                catalogue_id=cat_id if cat_id is None else str(cat_id),
                leg_index=j,
                cause=ValueError(
                    f"no Lambert solution n_revs={n_revs} branch={requested_branch!r}; "
                    f"available={[(s.n_revs, s.branch) for s in sols]}"
                ),
            )
        leg_vels.append(
            (
                np.asarray(chosen.v1, dtype=np.float64),
                np.asarray(chosen.v2, dtype=np.float64),
                chosen.n_revs,
                chosen.branch,
            )
        )

    # Build Encounters.
    encounters: list[Encounter] = []
    n_enc = len(chain)
    for i in range(n_enc):
        r_p, v_p = planet_states[i]
        if i == 0:
            v_dep = leg_vels[0][0]
            vinf_out = v_dep - v_p
            vinf_in = vinf_out
        elif i == n_enc - 1:
            v_arr = leg_vels[-1][1]
            vinf_in = v_arr - v_p
            vinf_out = vinf_in
        else:
            v_arr = leg_vels[i - 1][1]
            v_dep = leg_vels[i][0]
            vinf_in = v_arr - v_p
            vinf_out = v_dep - v_p
        encounters.append(
            Encounter(
                body=chain[i],
                t=encounter_times[i],
                r=r_p,
                v_planet=v_p,
                vinf_in=vinf_in,
                vinf_out=vinf_out,
            )
        )

    # Build Legs.
    legs: list[Leg] = []
    for j, _leg in enumerate(legs_meta):
        v_dep, v_arr, leg_n_revs, leg_branch = leg_vels[j]
        legs.append(
            Leg(
                from_body=chain[j],
                to_body=chain[j + 1],
                t_depart=encounter_times[j],
                t_arrive=encounter_times[j + 1],
                v_depart=v_dep,
                v_arrive=v_arr,
                n_revs=leg_n_revs,
                branch=leg_branch,
            )
        )

    # Period: prefer the catalogue's published value (schema-v2 canonical);
    # fall back to the leg-sum if absent (shouldn't happen for in-scope
    # entries, but the loader doesn't enforce non-null period).
    period_meta = catalogue_entry.get("period") or {}
    period_years = period_meta.get("years")
    if period_years is not None:
        period_s = float(period_years) * _YEAR_S
    else:
        period_s = encounter_times[-1] - encounter_times[0]

    return Cycler(
        bodies=list(chain),
        period=period_s,
        encounters=encounters,
        legs=legs,
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def _empty_result(
    *,
    cycler_id: str | None,
    n_cycles: int,
    v3_status: str,
    t_start_sec: float | None,
) -> RealClosureResult:
    """Build a RealClosureResult representing an early-termination outcome."""
    return RealClosureResult(
        cycler_id=cycler_id,
        n_cycles_propagated=0,
        max_drift_km=math.inf,
        per_cycle_drift_km=(),
        per_encounter_vinf_mismatch_kms=(),
        closes=False,
        v3_status=v3_status,
        horizon_tcm_mps=0.0,
        per_cycle_tcm_mps=(0.0,) * n_cycles,
        frame_used="dynamic",
        t_start_sec=t_start_sec,
    )


def _compute_horizon_tcm(
    constructed: Cycler,
    n_cycles: int,
    ephem: Ephemeris,
    *,
    tcm_perturbers: tuple[str, ...] | None = None,
) -> tuple[float, tuple[float, ...]]:
    """M7: real-ephemeris maintenance ΔV over ``n_cycles`` for a constructed cycler.

    Builds the continuous N-cycle encounter sequence from ``constructed`` (repeating
    the one-cycle encounters by the encounter-span offset so the inter-cycle HOME
    flyby is included — dropping it would undercount), then runs
    :func:`cyclerfinder.nbody.maintenance_shoot.continuous_maintenance_chain`
    (position-targeted, per-leg endpoint-excluded perturbers). Returns
    ``(horizon_tcm_mps, per_cycle_tcm_mps)``. On any non-converging leg the chain is
    honestly unmeasurable: returns ``(inf, (inf,)*n_cycles)`` — the row stays V0,
    never forced to a fake ballistic 0.0.

    Lazy import of the nbody/rebound stack keeps the module's cheap astropy-only
    import contract for the M6b path.
    """
    from cyclerfinder.nbody.maintenance_shoot import continuous_maintenance_chain
    from cyclerfinder.nbody.propagator import RestrictedNBody

    encounters = constructed.encounters
    m = len(encounters)
    if m < 2 or n_cycles < 1:
        return 0.0, (0.0,) * max(n_cycles, 0)

    # Repeat by the actual encounter span (not the published period) so the seam node
    # epoch is exact: encounters[0].t + span == encounters[-1].t, making the home
    # flyby between cycles continuous.
    span = encounters[-1].t - encounters[0].t
    node_epochs: list[float] = []
    node_bodies: list[str] = []
    for c in range(n_cycles):
        for j in range(m - 1):
            node_epochs.append(encounters[j].t + c * span)
            node_bodies.append(encounters[j].body)
    node_epochs.append(encounters[-1].t + (n_cycles - 1) * span)
    node_bodies.append(encounters[-1].body)

    # Per-leg seeds = the constructed cycler's rev-correct departure velocities
    # (construct_real_ephemeris_cycler solved Lambert at the catalogue n_revs/branch).
    # The concatenated N-cycle leg list is exactly the one-cycle legs repeated, so the
    # seeds repeat per cycle. Essential for multi-rev resonant cyclers (a single-rev
    # re-guess lands a high-energy transfer and inflates the TCM by orders of magnitude).
    leg_v_guess: list[NDArray[np.float64]] | None = None
    if len(constructed.legs) == m - 1:
        one_cycle = [np.asarray(leg.v_depart, dtype=np.float64) for leg in constructed.legs]
        leg_v_guess = one_cycle * n_cycles

    if tcm_perturbers is None:
        # System significant perturbers = the cycler's own flyby bodies. Endpoint
        # exclusion (chain default) drops each from legs it bounds; a body that is a
        # flyby target but NOT this leg's endpoint still perturbs that leg's cruise.
        tcm_perturbers = tuple(sorted(set(node_bodies)))

    prop = RestrictedNBody("rebound")
    chain = continuous_maintenance_chain(
        node_epochs,
        node_bodies,
        ephem,
        prop,
        cruise_perturbers=tcm_perturbers,
        leg_v_guess=leg_v_guess,
        n_cycles=n_cycles,
    )
    if chain.diverged:
        return float("inf"), (float("inf"),) * n_cycles
    # Per-cycle reporting: uniform average over the horizon (the position-targeted
    # model has no cross-cycle accumulation; real per-cycle variation across the
    # N-cycle walk is a Phase-2 refinement).
    per_cycle = chain.horizon_tcm_mps / n_cycles
    return chain.horizon_tcm_mps, (per_cycle,) * n_cycles


def verify_real_closure(
    cycler: Cycler | dict[str, Any],
    n_cycles: int,
    ephem: Ephemeris,
    *,
    t_start: float | None = None,
    frame_bodies: tuple[str, ...] | None = None,
    cycler_id: str | None = None,
    signature_priority_date: datetime | None = None,
    compute_tcm: bool = False,
    tcm_perturbers: tuple[str, ...] | None = None,
) -> RealClosureResult:
    """Spec §14 V2-real gate machinery; M6b's binding entry point.

    Pipeline (plan §3.1.2):

    1. Resolve ``t_start`` if ``None`` via :func:`_resolve_real_t_start`.
    2. If ``cycler`` is a dict (catalogue entry), build a real-ephemeris
       :class:`Cycler` via :func:`construct_real_ephemeris_cycler`;
       catch :class:`RealClosureConstructionError` →
       ``v3-construction-error``.
    3. Delegate propagation + drift measurement to
       :func:`verify_long_term_stability` (M6a).
       **This is the binding composition** asserted by
       ``test_real_closure_uses_m6a_machinery``.
    4. Compare ``stability.max_drift_km`` to
       :data:`REAL_DRIFT_TOLERANCE_KM`. Note: M6a's
       :attr:`StabilityReport.stable` uses its 50,000 km idealised
       tolerance; M6b's ``closes`` uses the 200,000 km real tolerance.
       The two disagree by design for cyclers in the 50k-200k km band.
    5. Populate V3 placeholders (zeros) and return.

    Parameters
    ----------
    cycler:
        A constructed :class:`Cycler` (used directly) or a catalogue
        entry dict (built into a real-ephemeris :class:`Cycler` via
        :func:`construct_real_ephemeris_cycler`).
    n_cycles:
        Number of cycles to propagate. Must be >= 2.
    ephem:
        Ephemeris backend; typically ``Ephemeris("astropy")``.
    t_start:
        J2000-relative launch epoch in seconds. If ``None``,
        :func:`_resolve_real_t_start` derives one from
        ``signature_priority_date``.
    frame_bodies:
        Optional override for the dynamic-frame anchor bodies. Forwards
        to :func:`verify_long_term_stability`.
    cycler_id:
        Catalogue id pass-through for the returned
        :class:`RealClosureResult`.
    signature_priority_date:
        Required when ``t_start`` is ``None`` and ``cycler`` is a dict
        or when the caller does not pre-pick a launch epoch.

    Returns
    -------
    RealClosureResult
        Frozen V2-real gate record. The ``closes`` boolean is the spec
        §14 V2-real headline; ``v3_status`` carries the precise reason
        for any non-pass outcome.

    Raises
    ------
    ValueError
        If ``n_cycles < 2`` (consecutive-pair drift is undefined for a
        single cycle) or if ``t_start`` and
        ``signature_priority_date`` are both ``None``.
    """
    if n_cycles < 2:
        raise ValueError(f"verify_real_closure requires n_cycles >= 2; got {n_cycles}")

    # ----- Step 1: resolve t_start -----
    resolved_t_start = t_start
    if resolved_t_start is None:
        if isinstance(cycler, dict):
            signature = phase_signature_from_catalogue_entry(cycler)
        else:
            # A constructed cycler with concrete leg ToFs; derive a
            # signature from its first encounter's body / leg ToFs /
            # |V∞| magnitudes.
            signature = _signature_from_cycler(cycler)
        priority = signature_priority_date
        if priority is None and isinstance(cycler, dict):
            priority = _parse_priority_date(cycler.get("priority_date"))
        if priority is None:
            raise ValueError(
                "verify_real_closure: t_start or signature_priority_date "
                "is required to resolve a real launch epoch."
            )
        resolved_t_start = _resolve_real_t_start(signature, ephem, priority)
        if resolved_t_start is None:
            return _empty_result(
                cycler_id=cycler_id,
                n_cycles=n_cycles,
                v3_status="v3-no-real-window",
                t_start_sec=None,
            )

    # ----- Step 2: construct cycler from catalogue dict if needed -----
    if isinstance(cycler, dict):
        try:
            constructed = construct_real_ephemeris_cycler(cycler, ephem, resolved_t_start)
        except RealClosureConstructionError:
            return _empty_result(
                cycler_id=cycler_id,
                n_cycles=n_cycles,
                v3_status="v3-construction-error",
                t_start_sec=resolved_t_start,
            )
    else:
        constructed = cycler

    # ----- Step 3: delegate to M6a -----
    stability: StabilityReport = verify_long_term_stability(
        constructed,
        n_laps=n_cycles,
        ephem=ephem,
        t_start=resolved_t_start,
        frame_bodies=frame_bodies,
        cycler_id=cycler_id,
        use_uniform_frame=False,
    )

    # ----- Step 4: M6b tolerance -----
    closes = stability.max_drift_km < REAL_DRIFT_TOLERANCE_KM
    v3_status = "v3-real-closure-pass" if closes else "v3-real-closure-fail"

    # ----- Step 5: diagnostics + result -----
    vinf_mismatch = _check_vinf_continuity(constructed, ephem)

    # ----- Step 5b (M7): real-ephemeris maintenance ΔV (horizon TCM) -----
    # Opt-in (compute_tcm): the Mars-perturbed/position-targeted shoot is minutes/row.
    # Off by default => the historical 0.0 placeholders, preserving M6b behaviour and
    # the cheap astropy-only import (the nbody/rebound import is lazy, inside the branch).
    horizon_tcm_mps = 0.0
    per_cycle_tcm_mps: tuple[float, ...] = (0.0,) * stability.n_laps_propagated
    if compute_tcm:
        horizon_tcm_mps, per_cycle_tcm_mps = _compute_horizon_tcm(
            constructed,
            stability.n_laps_propagated,
            ephem,
            tcm_perturbers=tcm_perturbers,
        )

    return RealClosureResult(
        cycler_id=cycler_id,
        n_cycles_propagated=stability.n_laps_propagated,
        max_drift_km=stability.max_drift_km,
        per_cycle_drift_km=stability.per_lap_drift_km,
        per_encounter_vinf_mismatch_kms=vinf_mismatch,
        closes=closes,
        v3_status=v3_status,
        horizon_tcm_mps=horizon_tcm_mps,
        per_cycle_tcm_mps=per_cycle_tcm_mps,
        frame_used="dynamic",
        t_start_sec=resolved_t_start,
    )


def _signature_from_cycler(cycler: Cycler) -> PhaseSignature:
    """Derive a :class:`PhaseSignature` from a constructed :class:`Cycler`.

    Maps the cycler's body chain to ``bodies``, per-leg
    ``(t_arrive - t_depart)`` to ``leg_durations_s``, and per-encounter
    ``|vinf_in|`` (or ``|vinf_out|`` at the first encounter) to
    ``vinf_target_kms``.

    Used by :func:`verify_real_closure` when the caller passes a
    :class:`Cycler` directly but does not supply ``t_start`` — the
    signature drives :func:`_resolve_real_t_start`.
    """
    bodies = tuple(cycler.bodies)
    leg_durations_s = tuple(leg.t_arrive - leg.t_depart for leg in cycler.legs)
    vinf_target_kms = tuple(
        float(np.linalg.norm(enc.vinf_in if i > 0 else enc.vinf_out))
        for i, enc in enumerate(cycler.encounters)
    )
    return PhaseSignature(
        bodies=bodies,
        leg_durations_s=leg_durations_s,
        vinf_target_kms=vinf_target_kms,
    )


__all__ = [
    "EXPECTED_SKIPS",
    "N_CYCLES_DEFAULT",
    "REAL_DRIFT_TOLERANCE_KM",
    "RealClosureConstructionError",
    "RealClosureResult",
    "construct_real_ephemeris_cycler",
    "verify_real_closure",
]
