"""Epoch-aware (non-repeating, window-bounded) trajectory data model.

Phase 1 of task #289 (Track-A capability for the four-class catalogue scope
admitted at schema v4.7, 2026-06-15). The non-cycler classes — ``quasi_cycler``,
``precursor_mga``, ``mga_tour`` — are *epoch-locked*: closure exists only inside
a (typically narrow) launch window, the spacecraft visits a finite sequence of
encounters, and there is no period to repeat. This module supplies the data
model + closure-evaluation substrate those classes need; Phase 2+ adds the
Tisserand-graph window enumerator, MGA chain optimisation, precursor-MGA
insertion matching, and V0-V5 gauntlet integration.

The substrate is deliberately thin: every dynamical primitive
(:func:`cyclerfinder.core.lambert.lambert`,
:class:`cyclerfinder.core.ephemeris.Ephemeris`,
:func:`cyclerfinder.core.flyby.bend_angle`,
:func:`cyclerfinder.core.flyby.max_bend`,
:func:`cyclerfinder.core.kepler.propagate`) already exists. The Phase-1
contribution is

  1. an ``EpochLockedTrajectory`` value type that encodes a candidate
     epoch-locked trajectory and its closure semantics,
  2. an ``EpochLockedClosure`` result type that records the residual,
     flyby continuity, and independent-cross-check disposition, and
  3. a :func:`close_epoch_locked` driver that runs per-leg Lambert against
     the real ephemeris and an independent Kepler re-propagation
     cross-check.

The Tito 2018 free-return reproduction (``scripts/tito_free_return_repro.py``)
is the reference closure pattern; this module generalises that script's
mechanics to an arbitrary ``E[-B-...-]E`` body sequence.

Phase 1 scope:

* Ballistic only (no DSMs). Matches Tito and matches the simplest
  ``mga_tour`` / ``precursor_mga`` cases. DSMs are a Phase 3+ extension.
* Single-rev prograde Lambert per leg. Multi-rev is a Phase-2+ option.
* No catalogue writeback. Phase 1 is substrate + tests only.

See ``docs/notes/2026-06-16-297-289-phase1-epoch-aware-genome.md`` for the
phase plan and ``docs/notes/2026-06-16-catalogue-scope-taxonomy.md`` for the
class semantics.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Literal

import numpy as np

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import bend_angle, max_bend
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import LambertSolution, lambert

OrbitClass = Literal["quasi_cycler", "precursor_mga", "mga_tour"]
"""Epoch-locked orbit classes admitted by this module.

The strict ``cycler`` class (period + invariants + repeating
``sequence_canonical``) is intentionally excluded: cyclers have no
``launch_epoch`` and no ``validity_window``, so they do not fit the
``EpochLockedTrajectory`` shape. See the schema v4.7 taxonomy note.
"""

_VALID_ORBIT_CLASSES: frozenset[str] = frozenset(("quasi_cycler", "precursor_mga", "mga_tour"))


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def utc_to_tsec(utc_iso: str) -> float:
    """Convert an ISO-8601 UTC string to ``t_sec`` (TDB s since J2000(TDB)).

    The astropy ephemeris backend's ``t_sec`` axis is the TDB-seconds offset
    from JD 2451545.0 TDB (see :mod:`cyclerfinder.core.ephemeris` module
    docstring). This helper matches the conversion used in
    ``scripts/tito_free_return_repro.py``.
    """
    from astropy.time import Time

    t = Time(utc_iso, scale="utc")
    return float((t.tdb.jd - 2451545.0) * 86400.0)


@dataclass(frozen=True)
class DSMSpec:
    """An optional deep-space-manoeuvre on one leg of an epoch-locked trajectory.

    A DSM splits one Lambert leg at ``fraction_along_leg`` (in [0, 1]),
    propagates the first segment ballistically up to that fraction, applies
    an instantaneous ``delta_v_kms`` vector in the heliocentric frame, then
    runs a fresh Lambert from the post-impulse position to the next body
    over the remaining fraction. Closure is the V_inf match at the FAR
    endpoint (after both Lamberts).

    The DSM is a *user-supplied* placement; Phase 4 will add automated
    placement. The data model and closure driver only know how to evaluate
    a hand-placed DSM, not where to place it.

    Class invariants enforced in :meth:`__post_init__`:

    * ``leg_index`` is non-negative.
    * ``fraction_along_leg`` is strictly inside ``(0, 1)`` (an endpoint DSM
      is degenerate — use a different ``leg_index`` or absorb into the
      flyby).
    * ``delta_v_kms`` is a 3-tuple of finite floats.
    """

    leg_index: int
    """0-based index of the leg this DSM belongs to.
    Must satisfy ``0 <= leg_index < len(leg_tofs_days)``."""

    fraction_along_leg: float
    """Fraction of the leg's TOF at which the DSM fires. In ``(0, 1)``.
    A value of 0.5 fires the DSM at the midpoint of the leg."""

    delta_v_kms: tuple[float, float, float]
    """Inertial heliocentric ``Delta V`` vector (km/s), applied to the
    spacecraft at the split point. The magnitude is the DSM cost."""

    def __post_init__(self) -> None:
        if self.leg_index < 0:
            raise ValueError(f"leg_index must be >= 0, got {self.leg_index}")
        if not (0.0 < self.fraction_along_leg < 1.0):
            raise ValueError(
                f"fraction_along_leg must be in (0, 1), got {self.fraction_along_leg}",
            )
        if len(self.delta_v_kms) != 3:
            raise ValueError(
                f"delta_v_kms must be a length-3 tuple, got {self.delta_v_kms!r}",
            )
        for k, c in enumerate(self.delta_v_kms):
            if not np.isfinite(c):
                raise ValueError(
                    f"delta_v_kms[{k}] = {c} is not finite",
                )


def _add_days_to_iso_utc(utc_iso: str, days: float) -> str:
    """Return an ISO-8601 UTC string offset from ``utc_iso`` by ``days``."""
    from astropy.time import Time, TimeDelta

    t = Time(utc_iso, scale="utc") + TimeDelta(days * 86400.0, format="sec")
    return str(t.utc.isot) + "Z"


# --------------------------------------------------------------------------- #
# Dataclasses
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class EpochLockedTrajectory:
    """A non-repeating (or finite-return) trajectory tied to a specific epoch window.

    Distinct from a CR3BP / cyclers-only :class:`Cycler` (period + invariants +
    repeating ``sequence_canonical``): has no period; has ``launch_epoch_utc``;
    closure is window-bounded; ``n_returns`` is finite by construction.

    Class invariants enforced in :meth:`__post_init__`:

    * ``orbit_class`` is one of ``quasi_cycler`` / ``precursor_mga`` /
      ``mga_tour`` (strict cyclers do not fit this shape).
    * ``n_returns`` is a positive integer.
    * ``precursor_mga`` requires a non-empty ``inserts_into`` catalogue id.
    * ``sequence`` length ``L``: ``leg_tofs_days`` has length ``L-1`` (one
      per leg) and ``vinf_kms_at_encounters`` has length ``L`` (one per
      body encounter).
    """

    sequence: tuple[str, ...]
    """Body codes in encounter order (e.g. ``("E", "M", "E")`` for Tito).
    Codes are one-letter keys from :data:`cyclerfinder.core.constants.PLANETS`.
    Sequence length must be >= 2."""

    leg_tofs_days: tuple[float, ...]
    """Time of flight per leg, days. Length = ``len(sequence) - 1``.
    Strictly positive."""

    vinf_kms_at_encounters: tuple[float, ...]
    """Published / sourced ``|V_inf|`` (km/s) at each encounter, in
    encounter order. Length = ``len(sequence)``. Used as the EXPECTED side
    of the closure residual."""

    launch_epoch_utc: str
    """ISO-8601 UTC of the trajectory start (first body encounter)."""

    orbit_class: OrbitClass
    """One of ``quasi_cycler`` / ``precursor_mga`` / ``mga_tour``. The strict
    ``cycler`` class is rejected — those have no epoch."""

    n_returns: int
    """Finite return count. ``mga_tour`` / ``precursor_mga`` typically
    ``n_returns=1``; ``quasi_cycler`` is 3-15."""

    validity_window_start_utc: str
    """Earliest closure UTC. For ``mga_tour`` this is typically equal to
    ``launch_epoch_utc``."""

    validity_window_end_utc: str
    """Latest closure UTC (terminal encounter epoch for ``mga_tour``;
    final return for ``quasi_cycler`` / ``precursor_mga``)."""

    inserts_into: str | None = None
    """Catalogue id of the cycler row this precursor inserts into. Required
    for ``precursor_mga``; ``None`` for ``mga_tour`` and ``quasi_cycler``."""

    periapsis_altitudes_km: tuple[float | None, ...] | None = None
    """Optional per-encounter periapsis altitude override (km above body
    equatorial radius). ``None`` (the default) at an entry means
    :data:`cyclerfinder.core.constants.SAFE_PERIHELION_KM` for that body is
    used. A non-``None`` entry is the *published* flyby periapsis (e.g.
    Tito's 100 km Mars periapsis vs the engineering default 300 km). If the
    field itself is ``None`` (the default), every encounter uses the safe
    default. If supplied, length must equal ``len(sequence)``."""

    dsm_specs: tuple[DSMSpec, ...] | None = None
    """Optional per-leg deep-space-manoeuvres. ``None`` (the default) is
    the ballistic case (Phase 1 / Tito reproduction). If supplied, every
    entry's ``leg_index`` must satisfy ``0 <= leg_index < len(leg_tofs_days)``.
    At most one DSM per leg is supported by the Phase-3 closure driver
    (the Phase-4 multi-DSM-per-leg extension is deferred). DSMs are
    *hand-placed* in Phase 3; automated placement is a Phase-4+ extension.

    See :class:`DSMSpec` for the per-DSM contract. The closure driver
    (:func:`close_epoch_locked`) reads this field; if it is non-``None``,
    each affected leg is split at the DSM fraction, propagated ballistically
    to the split, gets the ΔV applied, then runs a fresh Lambert from the
    post-impulse position to the next body."""

    closure_residual_kms: float | None = None
    """Worst-encounter ``|V_inf|`` residual vs ``vinf_kms_at_encounters`` after
    closure. Set by :func:`close_epoch_locked`; ``None`` before closure."""

    flyby_continuity_max_dv_kms: float | None = None
    """Max ballistic-continuity ``Delta V`` (km/s) at any intermediate flyby
    after closure. Set by :func:`close_epoch_locked`; ``None`` before closure."""

    ephemeris: str = "DE440"
    """Ephemeris kernel name. ``"DE440"`` for the astropy backend; other
    names are informational only."""

    notes: str = ""

    def __post_init__(self) -> None:
        # 1. Sequence shape.
        if len(self.sequence) < 2:
            raise ValueError(
                f"sequence must have >= 2 body codes (legs requires both endpoints); "
                f"got {self.sequence!r}"
            )
        n_legs = len(self.sequence) - 1
        if len(self.leg_tofs_days) != n_legs:
            raise ValueError(
                f"leg_tofs_days length {len(self.leg_tofs_days)} != len(sequence) - 1 = {n_legs}",
            )
        if len(self.vinf_kms_at_encounters) != len(self.sequence):
            raise ValueError(
                f"vinf_kms_at_encounters length {len(self.vinf_kms_at_encounters)} "
                f"!= len(sequence) = {len(self.sequence)}",
            )

        # 2. Leg TOF positivity.
        for i, tof in enumerate(self.leg_tofs_days):
            if not (tof > 0.0):
                raise ValueError(f"leg_tofs_days[{i}] = {tof} is not strictly positive")

        # 3. Body code validity.
        unknown = [b for b in self.sequence if b not in PLANETS]
        if unknown:
            raise ValueError(
                f"unknown body code(s) in sequence: {unknown!r}; "
                f"valid codes are {tuple(PLANETS)!r}",
            )

        # 4. Class invariants.
        if self.orbit_class not in _VALID_ORBIT_CLASSES:
            raise ValueError(
                f"orbit_class={self.orbit_class!r} not admitted by EpochLockedTrajectory; "
                f"strict 'cycler' does not fit this shape (it has no epoch). "
                f"Valid: {sorted(_VALID_ORBIT_CLASSES)!r}",
            )

        if not isinstance(self.n_returns, int) or self.n_returns < 1:
            raise ValueError(
                f"n_returns must be a positive integer (>= 1); got {self.n_returns!r}",
            )

        if self.orbit_class == "precursor_mga" and not self.inserts_into:
            raise ValueError(
                "precursor_mga requires a non-empty 'inserts_into' catalogue id "
                "(per schema v4.7 invariant; see "
                "docs/notes/2026-06-16-catalogue-scope-taxonomy.md)",
            )

        # 5. Periapsis-altitude override shape (if supplied).
        if self.periapsis_altitudes_km is not None and len(self.periapsis_altitudes_km) != len(
            self.sequence
        ):
            raise ValueError(
                f"periapsis_altitudes_km length {len(self.periapsis_altitudes_km)} "
                f"!= len(sequence) = {len(self.sequence)}",
            )

        # 6. DSM specs validation (Phase 3 extension).
        if self.dsm_specs is not None:
            seen_legs: set[int] = set()
            for k, spec in enumerate(self.dsm_specs):
                if not isinstance(spec, DSMSpec):
                    raise ValueError(
                        f"dsm_specs[{k}] is not a DSMSpec instance: {spec!r}",
                    )
                if not (0 <= spec.leg_index < n_legs):
                    raise ValueError(
                        f"dsm_specs[{k}].leg_index = {spec.leg_index} is outside "
                        f"[0, {n_legs}) (n_legs = len(sequence) - 1)",
                    )
                if spec.leg_index in seen_legs:
                    raise ValueError(
                        f"dsm_specs[{k}].leg_index = {spec.leg_index} is duplicated; "
                        "Phase 3 supports at most one DSM per leg",
                    )
                seen_legs.add(spec.leg_index)


@dataclass(frozen=True)
class EpochLockedClosure:
    """Result of closing an :class:`EpochLockedTrajectory` at a specific epoch.

    Carries the per-leg Lambert solutions, the worst-encounter ``V_inf``
    residual vs the trajectory's published expected values, the worst-flyby
    ballistic-continuity ``Delta V``, and an independent-integrator
    cross-check residual (mandatory per the orbit-closure discipline:
    ``feedback_orbit_closure_discipline``).

    A closure is **converged** iff every gate passes:

    * ``closure_residual_kms <= closure_tol_kms``
    * ``flyby_continuity_max_dv_kms <= flyby_continuity_tol_kms``
    * if ``independent_check_residual_kms is not None``,
      ``independent_check_residual_kms <= independent_tol_kms``
    """

    trajectory: EpochLockedTrajectory
    closure_residual_kms: float
    flyby_continuity_max_dv_kms: float
    per_leg_lambert_solutions: tuple[LambertSolution, ...]
    per_encounter_vinf_kms: tuple[float, ...]
    """Our DE440 ``|V_inf|`` (km/s) at each encounter, ordered matching
    ``trajectory.sequence``. For an intermediate body the speed is the mean of
    the inbound and outbound legs (the ballistic-continuity ideal); the
    individual leg speeds are recoverable from ``per_leg_lambert_solutions``.
    """
    independent_check_residual_kms: float | None
    converged: bool
    dsm_delta_v_kms_per_leg: tuple[float, ...] = ()
    """Per-DSM ``Delta V`` magnitudes (km/s), one entry per DSM in
    ``trajectory.dsm_specs`` (empty for a ballistic closure). Ordered the
    same as ``trajectory.dsm_specs``. Phase 3+; for Phase-1 / ballistic
    closures this is ``()`` (empty tuple)."""
    notes: str = ""


# --------------------------------------------------------------------------- #
# Closure driver
# --------------------------------------------------------------------------- #


def close_epoch_locked(
    trajectory_spec: EpochLockedTrajectory,
    ephemeris: Ephemeris,
    *,
    closure_tol_kms: float = 0.5,
    flyby_continuity_tol_kms: float = 0.05,
    independent_cross_check: bool = True,
    independent_tol_kms: float = 0.1,
) -> EpochLockedClosure:
    """Close a candidate epoch-locked trajectory at its ``launch_epoch_utc``.

    Per-leg method (matches ``scripts/tito_free_return_repro.py``):

      1. Convert each encounter UTC to ``t_sec`` (TDB seconds since J2000(TDB)).
         The encounter epochs are derived by cumulatively adding
         ``leg_tofs_days`` to ``launch_epoch_utc``.
      2. Real-ephemeris heliocentric body state at each encounter via
         ``ephemeris.state(body, t_sec)``.
      3. Per leg: single-rev prograde Lambert
         ``lambert(r_a, r_b, tof, prograde=True, max_revs=0)``. The selected
         solution is the first / single-rev branch.
      4. ``|V_inf|`` at each encounter:

         * first body: ``|v1_leg0 - v_body|``
         * last body: ``|v2_legN - v_body|``
         * intermediate body: mean of inbound ``|v2 - v_body|`` and outbound
           ``|v1 - v_body|`` (the ballistic-continuity ideal)

      5. Flyby continuity at every intermediate body: compute the
         speed mismatch ``abs(|v_inf_in| - |v_inf_out|)`` AND the bend deficit
         ``max(0, bend_angle - max_bend)`` paid at the published / safe
         periapsis. Both are added to a single per-flyby ``dv`` and the worst
         across all intermediates is reported.
      6. If ``independent_cross_check`` is on, each leg's Lambert
         ``(v1, tof)`` is re-propagated with the universal-variable Kepler
         propagator (algorithmically independent from the Lambert UV
         boundary-value solver); the worst ``|r2_propagated - r2_ephemeris|``
         residual divided by the leg's chord length is reported as the
         dimensionless cross-check signal, multiplied by the encounter
         ``|V_inf|`` scale to express the disagreement as an effective km/s.

    Parameters
    ----------
    trajectory_spec:
        The :class:`EpochLockedTrajectory` to close. Read-only.
    ephemeris:
        Body-state provider. For real-ephemeris closure (the typical case)
        pass ``Ephemeris(model="astropy")``. The analytic-circular backend is
        accepted for unit testing.
    closure_tol_kms:
        Max allowed ``|V_inf|`` residual vs ``vinf_kms_at_encounters`` at any
        encounter, km/s. Default ``0.5`` km/s is the Tito-reproduction-grade
        ``mga_tour`` gate (Tito reproduction's worst residual is ~0.063 km/s
        on DE440 vs DE421).
    flyby_continuity_tol_kms:
        Max allowed ballistic-continuity ``Delta V`` (km/s) at any
        intermediate flyby. Default ``0.05`` is Tito-grade (his Mars flyby
        leaves 0.011 km/s on DE440).
    independent_cross_check:
        Run the universal-variable Kepler re-propagation cross-check. Default
        ``True``; turn off only for fast unit tests that exercise the
        residual machinery alone.
    independent_tol_kms:
        Independent-cross-check disagreement tolerance, km/s. Default
        ``0.1`` km/s.

    Returns
    -------
    EpochLockedClosure
        The closure result. ``converged`` is True iff all three gates pass.

    Raises
    ------
    ValueError
        On non-positive tolerances.
    """
    if closure_tol_kms <= 0.0:
        raise ValueError(f"closure_tol_kms must be positive, got {closure_tol_kms}")
    if flyby_continuity_tol_kms <= 0.0:
        raise ValueError(
            f"flyby_continuity_tol_kms must be positive, got {flyby_continuity_tol_kms}",
        )
    if independent_tol_kms <= 0.0:
        raise ValueError(f"independent_tol_kms must be positive, got {independent_tol_kms}")

    seq = trajectory_spec.sequence
    n_legs = len(seq) - 1

    # Encounter epochs in t_sec (TDB).
    t_encounter_sec: list[float] = [utc_to_tsec(trajectory_spec.launch_epoch_utc)]
    for tof_d in trajectory_spec.leg_tofs_days:
        t_encounter_sec.append(t_encounter_sec[-1] + tof_d * 86400.0)

    # Real-ephemeris body states at each encounter.
    body_states: list[tuple[np.ndarray, np.ndarray]] = [
        ephemeris.state(seq[i], t_encounter_sec[i]) for i in range(len(seq))
    ]

    # Build a quick lookup for DSMs by leg index (Phase 3 DSM extension).
    dsm_by_leg: dict[int, DSMSpec] = {}
    if trajectory_spec.dsm_specs is not None:
        for spec in trajectory_spec.dsm_specs:
            dsm_by_leg[spec.leg_index] = spec

    # Per-leg single-rev prograde Lambert. For a leg with a DSM, we model
    # the leg as a TWO-arc trajectory that ACTUALLY EXECUTES the user's
    # prescribed ΔV:
    #
    #   * Arc 1: depart body_k with the ballistic-Lambert v1, propagate
    #     ballistically for ``fraction_along_leg * tof`` seconds; reach
    #     r_dsm, v_dsm_pre.
    #   * Apply user ΔV: v_dsm_post = v_dsm_pre + ΔV (heliocentric inertial).
    #   * Arc 2: propagate (r_dsm, v_dsm_post) for the REMAINING TOF
    #     forward to the leg's arrival time. The propagated arrival is r2_arc
    #     (may differ from r2_body when ΔV is non-zero) and v2_arc.
    #
    # The leg's "v2" used by downstream V_inf computation is v2_arc — i.e.
    # the velocity the spacecraft actually has when it reaches the arrival
    # epoch. The arrival V_inf is then v2_arc - v_body_arrival; if the DSM
    # was *targeted* (the user picked it to redirect to body_{k+1}), the
    # closure residual will pick up the magnitude of the rendezvous miss
    # via the V_inf-mismatch term. The DSM cost reported on the closure
    # record is the magnitude of the user-prescribed ΔV (the "actuator
    # cost"); the closure residual measures whether the prescribed DSM
    # actually closes the trajectory.
    per_leg_sols: list[LambertSolution] = []
    # DSM diagnostic data — empty when there are no DSMs (backward-compatible).
    per_dsm_delta_v_kms: dict[int, float] = {}
    for k in range(n_legs):
        r1, _ = body_states[k]
        r2, _ = body_states[k + 1]
        tof_s = t_encounter_sec[k + 1] - t_encounter_sec[k]
        if k not in dsm_by_leg:
            sols = lambert(r1, r2, tof_s, prograde=True, max_revs=0)
            per_leg_sols.append(sols[0])
            continue
        dsm = dsm_by_leg[k]
        t_dsm_s = tof_s * dsm.fraction_along_leg
        # Ballistic Lambert (r1 -> r2 over tof). Used to seed the leg's v1
        # (the spacecraft's departure velocity at body_k); the DSM perturbs
        # the trajectory partway through.
        sols_ref = lambert(r1, r2, tof_s, prograde=True, max_revs=0)
        v1_ref = sols_ref[0].v1
        # Arc 1: propagate ballistically to the DSM time.
        r_dsm, v_dsm_pre = propagate(r1, v1_ref, t_dsm_s)
        # Apply the user ΔV.
        dv_vec = np.asarray(dsm.delta_v_kms, dtype=np.float64)
        v_dsm_post = np.asarray(v_dsm_pre, dtype=np.float64) + dv_vec
        # Arc 2: propagate the post-impulse state forward to the leg's
        # arrival epoch. r2_arc is where the spacecraft ACTUALLY arrives; it
        # equals r2 only when the DSM happens to be the right correction.
        tof_remain_s = tof_s - t_dsm_s
        _r2_arc, v2_arc = propagate(
            np.asarray(r_dsm, dtype=np.float64),
            v_dsm_post,
            tof_remain_s,
        )
        # The leg solution reports the actual v1 / v2 the spacecraft has
        # at the leg endpoints. v1 = v1_ref (ballistic Lambert departure
        # from body_k); v2 = v2_arc (post-DSM forward propagation).
        synthetic = LambertSolution(
            n_revs=sols_ref[0].n_revs,
            branch=sols_ref[0].branch,
            v1=v1_ref,
            v2=np.asarray(v2_arc, dtype=np.float64),
        )
        per_leg_sols.append(synthetic)
        # DSM cost = magnitude of the prescribed ΔV (the spacecraft pays
        # this regardless of whether the DSM is well-targeted).
        per_dsm_delta_v_kms[k] = float(np.linalg.norm(dv_vec))

    # Per-encounter |V_inf|: at endpoints it's the leg V_inf; at intermediate
    # bodies it's the MEAN of inbound and outbound (the ballistic-continuity
    # ideal; the raw per-leg values land in per_leg_sols for diagnostics).
    per_encounter_vinf: list[float] = []
    # Endpoint 0 (departure).
    _, v_body_0 = body_states[0]
    vinf_dep_vec = per_leg_sols[0].v1 - v_body_0
    per_encounter_vinf.append(float(np.linalg.norm(vinf_dep_vec)))
    # Intermediate bodies.
    for k in range(1, len(seq) - 1):
        _, v_body_k = body_states[k]
        vinf_in_vec = per_leg_sols[k - 1].v2 - v_body_k
        vinf_out_vec = per_leg_sols[k].v1 - v_body_k
        v_in = float(np.linalg.norm(vinf_in_vec))
        v_out = float(np.linalg.norm(vinf_out_vec))
        per_encounter_vinf.append(0.5 * (v_in + v_out))
    # Endpoint N (arrival).
    _, v_body_n = body_states[-1]
    vinf_arr_vec = per_leg_sols[-1].v2 - v_body_n
    per_encounter_vinf.append(float(np.linalg.norm(vinf_arr_vec)))

    # Closure residual: worst absolute |V_inf| disagreement vs expected.
    expected = trajectory_spec.vinf_kms_at_encounters
    closure_residual = max(
        abs(ours - exp) for ours, exp in zip(per_encounter_vinf, expected, strict=True)
    )

    # Flyby continuity at each intermediate body: speed mismatch + bend deficit
    # at the safe periapsis (or the per-encounter override, if supplied). The
    # "max dv" is the worst per-flyby total.
    flyby_continuity_max_dv = 0.0
    for k in range(1, len(seq) - 1):
        body_code = seq[k]
        _, v_body_k = body_states[k]
        vinf_in_vec = per_leg_sols[k - 1].v2 - v_body_k
        vinf_out_vec = per_leg_sols[k].v1 - v_body_k
        v_in = float(np.linalg.norm(vinf_in_vec))
        v_out = float(np.linalg.norm(vinf_out_vec))
        speed_mismatch = abs(v_in - v_out)
        delta_req = bend_angle(vinf_in_vec, vinf_out_vec)
        mu = PLANETS[body_code].mu_km3_s2
        # Per-encounter periapsis-altitude override (e.g. Tito's published
        # 100 km Mars periapsis vs the engineering safe 300 km default).
        rp = SAFE_PERIHELION_KM[body_code]
        if trajectory_spec.periapsis_altitudes_km is not None:
            alt_override = trajectory_spec.periapsis_altitudes_km[k]
            if alt_override is not None:
                rp = PLANETS[body_code].radius_eq_km + float(alt_override)
        v_mean = 0.5 * (v_in + v_out)
        delta_cone = max_bend(mu, rp, v_mean)
        # Bend deficit charged as 2 * v_mean * sin(deficit / 2) -- the
        # Strange-Longuski single-impulse asymptote-rotation surrogate.
        # The flyby module's flyby_dv() composes the same two parts but
        # short-circuits to 0.0 inside the cone; we keep the parts visible
        # here so the test can introspect each gate cleanly.
        deficit = max(0.0, delta_req - delta_cone)
        bend_dv = 2.0 * v_mean * float(np.sin(0.5 * deficit))
        per_flyby_dv = speed_mismatch + bend_dv
        if per_flyby_dv > flyby_continuity_max_dv:
            flyby_continuity_max_dv = per_flyby_dv

    # Independent cross-check: re-propagate each leg's (r1, v1, tof) with the
    # universal-variable Kepler propagator and compare arrival to the
    # ephemeris position r2. The boundary-value Lambert solver and the
    # initial-value Kepler propagator are algorithmically independent
    # (different residual functions, different convergence loops); their
    # agreement on the SAME conic is the cross-check.
    independent_residual: float | None = None
    if independent_cross_check:
        worst_pos_km = 0.0
        worst_chord_km = 1.0
        for k in range(n_legs):
            r1, _ = body_states[k]
            r2_eph, _ = body_states[k + 1]
            tof_s = t_encounter_sec[k + 1] - t_encounter_sec[k]
            v1 = per_leg_sols[k].v1
            if k in dsm_by_leg:
                # Match the closure: propagate to the DSM fraction, apply
                # the ΔV, then propagate the remainder.
                dsm = dsm_by_leg[k]
                t_dsm_s = tof_s * dsm.fraction_along_leg
                r_dsm, v_dsm_pre = propagate(r1, v1, t_dsm_s)
                v_dsm_post = np.asarray(v_dsm_pre, dtype=np.float64) + np.asarray(
                    dsm.delta_v_kms,
                    dtype=np.float64,
                )
                r2_prop, _v2_prop = propagate(
                    np.asarray(r_dsm, dtype=np.float64),
                    v_dsm_post,
                    tof_s - t_dsm_s,
                )
            else:
                r2_prop, _v2_prop = propagate(r1, v1, tof_s)
            pos_err = float(np.linalg.norm(np.asarray(r2_prop) - np.asarray(r2_eph)))
            chord = float(np.linalg.norm(np.asarray(r2_eph) - np.asarray(r1)))
            if pos_err > worst_pos_km:
                worst_pos_km = pos_err
                worst_chord_km = max(chord, 1.0)
        # Translate dimensionless miss / chord into an effective km/s by
        # multiplying by the encounter V_inf scale (the orbital-speed scale
        # that maps a positional mismatch to an effective velocity defect).
        # The result is a pessimistic single-number disagreement: a true
        # solver-to-solver disagreement of x km on a chord of L km roughly
        # corresponds to (x/L) * |V_inf| km/s of equivalent velocity error.
        vinf_scale = max(per_encounter_vinf) if per_encounter_vinf else 1.0
        independent_residual = (worst_pos_km / worst_chord_km) * vinf_scale

    closure_gate = closure_residual <= closure_tol_kms
    flyby_gate = flyby_continuity_max_dv <= flyby_continuity_tol_kms
    indep_gate = (
        True if independent_residual is None else independent_residual <= independent_tol_kms
    )
    converged = closure_gate and flyby_gate and indep_gate

    # Per-DSM ΔV diagnostic, ordered to match trajectory_spec.dsm_specs
    # (so caller can zip the two).
    dsm_dv_tuple: tuple[float, ...] = ()
    if trajectory_spec.dsm_specs is not None:
        dsm_dv_tuple = tuple(
            per_dsm_delta_v_kms.get(spec.leg_index, 0.0) for spec in trajectory_spec.dsm_specs
        )

    return EpochLockedClosure(
        trajectory=trajectory_spec,
        closure_residual_kms=float(closure_residual),
        flyby_continuity_max_dv_kms=float(flyby_continuity_max_dv),
        per_leg_lambert_solutions=tuple(per_leg_sols),
        per_encounter_vinf_kms=tuple(per_encounter_vinf),
        independent_check_residual_kms=(
            None if independent_residual is None else float(independent_residual)
        ),
        converged=bool(converged),
        dsm_delta_v_kms_per_leg=dsm_dv_tuple,
    )


# --------------------------------------------------------------------------- #
# Window search
# --------------------------------------------------------------------------- #


def search_validity_window(
    trajectory_spec: EpochLockedTrajectory,
    ephemeris: Ephemeris,
    *,
    epoch_grid_step_days: float = 1.0,
    epoch_grid_padding_days: float = 7.0,
    **closure_kwargs: Any,
) -> list[EpochLockedClosure]:
    """Walk ``launch_epoch_utc`` across a grid and return converged closures.

    Walks the launch epoch over
    ``[nominal - padding, nominal + padding]`` at step ``epoch_grid_step_days``
    days, runs :func:`close_epoch_locked` at each grid point with the
    flight times held fixed, and returns the closures whose
    ``converged is True``.

    For ``mga_tour`` the function typically yields a single grid point (the
    nominal launch epoch) because the closure gate is tight and the launch
    window is geometrically narrow. For ``quasi_cycler`` it typically
    yields a series of closures spaced by the closure cycle (the
    cycle-of-opportunity pattern), and the caller is expected to widen the
    padding accordingly.

    Parameters
    ----------
    trajectory_spec:
        Candidate trajectory. The ``launch_epoch_utc`` is treated as the
        nominal centre of the search window.
    ephemeris:
        Body-state provider.
    epoch_grid_step_days:
        Grid step in days. Default 1.0 (24 h). Must be > 0.
    epoch_grid_padding_days:
        Half-width of the search window in days, on either side of the
        nominal launch epoch. Default 7.0 days (±1 week).
    **closure_kwargs:
        Forwarded to :func:`close_epoch_locked`.

    Returns
    -------
    list[EpochLockedClosure]
        Converged closures, in ascending launch-epoch order.
    """
    if epoch_grid_step_days <= 0.0:
        raise ValueError(
            f"epoch_grid_step_days must be positive, got {epoch_grid_step_days}",
        )
    if epoch_grid_padding_days < 0.0:
        raise ValueError(
            f"epoch_grid_padding_days must be non-negative, got {epoch_grid_padding_days}",
        )

    # Build the grid in days off the nominal launch epoch.
    n_each_side = int(np.floor(epoch_grid_padding_days / epoch_grid_step_days))
    offsets_days = [i * epoch_grid_step_days for i in range(-n_each_side, n_each_side + 1)]

    found: list[EpochLockedClosure] = []
    for d in offsets_days:
        launch_iso = _add_days_to_iso_utc(trajectory_spec.launch_epoch_utc, d)
        candidate = replace(trajectory_spec, launch_epoch_utc=launch_iso)
        try:
            cl = close_epoch_locked(candidate, ephemeris, **closure_kwargs)
        except Exception:
            # A grid point may go through a singular Lambert geometry; we
            # treat that as "not closed at this epoch" and move on.
            continue
        if cl.converged:
            found.append(cl)
    return found


__all__ = [
    "DSMSpec",
    "EpochLockedClosure",
    "EpochLockedTrajectory",
    "OrbitClass",
    "close_epoch_locked",
    "search_validity_window",
    "utc_to_tsec",
]
