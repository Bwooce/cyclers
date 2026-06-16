r"""#318 Phase 1 - multi-axis joint-search substrate.

A small *composition driver* that wires four independently-developed cycler-
discovery axes into a single search cell so a candidate can be parameterised
across all of them at once. The four axes are:

  A. **Powered maintenance** (#309).
     :func:`cyclerfinder.search.low_thrust_cycler_search.search_low_thrust_cyclers`
     drives the Lambert-closed periodic cycler at minimum per-synodic
     maintenance dV; the maintenance value is then wrapped in the Tsiolkovsky
     propellant accounting and tested as a Sims-Flanagan distributed-thrust
     leg.
  B. **Multi-revolution Lambert** (existing in
     :func:`cyclerfinder.core.lambert.lambert`).
     The discovery sweeps so far have left ``max_revs=0`` on every leg; the
     joint driver opens ``per_leg_revs`` to a *grid* per cell, threading
     ``per_leg_revs`` into the same ``optimise_maintenance_dv`` call the
     powered driver uses.
  C. **3D / broken-plane** (#291).
     :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
     is the full-3D CR3BP corrector. Phase 1 RECORDS the requested z0 amplitude
     on every candidate (the joint-axis index); the *core* cycler solve still
     runs in the 2D Lambert engine (composition, not rewrite). The recorded z0
     amplitude is the joint-axis Phase 2 hook - Phase 2 will pipe the 3D
     corrector in as the primary residual when the scientific case requires it.
  D. **Epoch-locked validity windows** (#289).
     :func:`cyclerfinder.genome.epoch_aware_genome.close_epoch_locked` is the
     real-ephemeris closure driver. When the joint driver is given a
     non-``None`` ``launch_epoch_utc``, it calls ``close_epoch_locked`` on the
     baseline cycler's encounter tour and records the worst V_inf residual +
     worst-flyby ballistic-continuity dV + independent-cross-check residual.

Phase 1 scope
-------------
Substrate, not a full discovery sweep:

  * :class:`JointAxisCandidate` - value type carrying the four-axis coords
    PLUS the verdict fields the existing axis-drivers already report.
  * :func:`joint_axis_search` - composes the four axes over the user's
    Cartesian product of (powered budgets x n_revs grids x z0 amplitudes x
    launch epochs) for ONE primary_sequence.
  * Phase 1 = the Cartesian product. Phase 2+ adds smarter sampling (Sobol,
    surrogate-driven, etc.). Phase 1 is the *capability*; Phase 2 is the
    efficiency.

The four-axis joint sweep has not been done in the published cycler
literature (per #328 / #302 / #309 surveys: each axis has been published in
isolation). If the joint manifold has genuinely-novel pockets, the Phase 1
substrate + a small EM probe is what surfaces the cells that warrant a
Phase 2 follow-up.

What this module does NOT do
----------------------------
  * No catalogue writeback.
  * No novelty claims (the literature-check verdict reuses #309 / #256 wiring;
    a NOT-found cleared the necessary-not-sufficient novelty gate but the
    V0-V5 gauntlet still governs).
  * No modification of any of the four axis modules. Composition only.
  * No "true cycler" verdict: a Phase 1 candidate is a CELL of the joint
    sweep, not a closed-three-period cycler. The V0-V5 floor applies.

Discipline (per ``feedback_orbit_closure_discipline``)
------------------------------------------------------
  * Independent cross-check is mandatory at the powered-axis layer (a seed
    bump must reproduce the maintenance dV) AND at the epoch-locked layer
    (Kepler re-propagation per leg).
  * The sourced golden anchor is the catalogue row
    ``aldrin-classic-em-k1-outbound``: at the (powered=0, n_revs=0, z0=0,
    epoch=None) corner the joint driver re-finds Aldrin.

References
----------
  * #309 substrate: :mod:`cyclerfinder.search.low_thrust_cycler_search`.
  * #291 substrate: :mod:`cyclerfinder.search.cr3bp_general_periodic_3d`.
  * #289 substrate: :mod:`cyclerfinder.genome.epoch_aware_genome`.
  * #318 phase plan: ``docs/notes/2026-06-318-multi-axis-joint-search-phase1.md``.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from typing import Any

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.literature_check import SearchFn
from cyclerfinder.search.low_thrust_cycler_search import (
    LowThrustCyclerCandidate,
    datetime_to_sec_since_j2000,
    search_low_thrust_cyclers,
)

# ---------------------------------------------------------------------------
# Axis defaults
# ---------------------------------------------------------------------------

_DEFAULT_POWERED_BUDGETS_KMS: tuple[float, ...] = (0.0,)
"""Default Axis-A grid: pure ballistic."""

_DEFAULT_N_REVS_GRID_PER_LEG: tuple[tuple[int, ...], ...] = ((0,),)
"""Default Axis-B grid: zero-rev (single-rev) Lambert on every leg."""

_DEFAULT_Z0_AMPLITUDES_NONDIM: tuple[float, ...] = (0.0,)
"""Default Axis-C grid: planar (z0 = 0)."""

_DEFAULT_LAUNCH_EPOCH_GRID: tuple[str | None, ...] = (None,)
"""Default Axis-D grid: epoch-blind."""


# ---------------------------------------------------------------------------
# Joint candidate data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JointAxisCandidate:
    """A cycler candidate parameterised across all four axes simultaneously.

    The structural fingerprint (sequence, V_inf tuple, leg ToFs) carries the
    trajectory identity; the four axis-coord fields tag the cell of the joint
    sweep the candidate came from; the verdict fields record per-axis closure
    outcomes.

    Phase 1 emits these to JSONL only - no catalogue writeback.
    """

    sequence: tuple[str, ...]
    vinf_tuple_kms: tuple[float, ...]
    leg_tofs_days: tuple[float, ...]

    # ----- Axis A - Powered maintenance ------------------------------------
    powered_maintenance_dv_kms_per_synodic: float
    """The maintenance dV the powered driver computed at this cell, km/s."""
    powered_budget_kms_requested: float
    """The cell's Axis-A coordinate (km/s). 0.0 = pure ballistic; >0 = the
    powered cell. The driver does NOT reject candidates that exceed the
    budget; the relation is recorded on the candidate so a Phase 2 promotion
    can apply a hard gate."""
    powered_leg_indices: tuple[int, ...]
    """Indices of legs that the powered burn is distributed across (Phase 1:
    all legs; Phase 2 will localise)."""
    sims_flanagan_feasible: bool
    """The Sims-Flanagan per-segment capability witness from the #309
    driver - does the requested dV fit in the N-segment train?"""

    # ----- Axis B - Multi-revolution Lambert -------------------------------
    n_revs_per_leg: tuple[int, ...]
    """The cell's Axis-B coordinate (per-leg n_revs). 0 = direct (single-rev);
    1+ = multi-rev."""

    # ----- Axis C - 3D / broken-plane --------------------------------------
    z0_amplitude_nondim: float
    """The cell's Axis-C coordinate (non-dim z0). 0 = planar; >0 = a 3D
    request. Phase 1: the value is RECORDED on the candidate, not yet
    driven into the 2D Lambert engine. A Phase 2 follow-up iterates over
    candidates with ``z0_amplitude_nondim > 0`` and pipes them through
    the #291 3D corrector for a same-model 3D closure."""

    # ----- Axis D - Epoch-locked validity window ---------------------------
    launch_epoch_utc: str | None
    """The cell's Axis-D coordinate. ``None`` = epoch-blind."""
    validity_window_start_utc: str | None
    """Earliest validity UTC if the epoch-locked closure ran; ``None``
    otherwise."""
    validity_window_end_utc: str | None
    """Latest validity UTC if the epoch-locked closure ran; ``None``
    otherwise."""
    epoch_locked_closure_residual_kms: float | None = None
    """Worst per-encounter V_inf residual from :func:`close_epoch_locked`
    when Axis D is active; ``None`` for an epoch-blind cell."""
    epoch_locked_flyby_continuity_max_dv_kms: float | None = None
    """Worst-flyby ballistic-continuity dV from :func:`close_epoch_locked`
    when Axis D is active; ``None`` otherwise."""
    epoch_locked_independent_check_residual_kms: float | None = None
    """Independent Kepler-propagation cross-check residual from
    :func:`close_epoch_locked`; ``None`` otherwise."""

    # ----- Shared verdict fields ------------------------------------------
    closure_residual_kms: float = math.nan
    """Lambert closure residual at the powered driver's optimum, km/s.
    This is the per-synodic maintenance dV - see the #309 docstring."""
    flyby_continuity_max_dv_kms: float = math.nan
    """Worst-flyby ballistic-continuity dV (km/s). For a Phase 1
    epoch-blind cell this falls back to the powered driver's
    ``maintenance_dv_kms``. For an epoch-locked cell this equals
    :attr:`epoch_locked_flyby_continuity_max_dv_kms`."""
    independent_cross_check_residual_kms: float = math.nan
    """Seed-bump independent cross-check from the #309 driver, km/s."""
    powered_axis_record: LowThrustCyclerCandidate | None = None
    """Full #309 candidate row (Tsiolkovsky wrap, lit-check, etc.)
    underlying this cell."""
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict[str, Any]:
        """Serialise to a JSONL-friendly dict (no NaN; floats round-trip)."""

        def _maybe(v: float | None) -> float | None:
            if v is None:
                return None
            if math.isnan(v):
                return None
            return float(v)

        return {
            "sequence": list(self.sequence),
            "vinf_tuple_kms": list(self.vinf_tuple_kms),
            "leg_tofs_days": list(self.leg_tofs_days),
            "powered_maintenance_dv_kms_per_synodic": float(
                self.powered_maintenance_dv_kms_per_synodic
            ),
            "powered_budget_kms_requested": float(self.powered_budget_kms_requested),
            "powered_leg_indices": list(self.powered_leg_indices),
            "sims_flanagan_feasible": bool(self.sims_flanagan_feasible),
            "n_revs_per_leg": list(self.n_revs_per_leg),
            "z0_amplitude_nondim": float(self.z0_amplitude_nondim),
            "launch_epoch_utc": self.launch_epoch_utc,
            "validity_window_start_utc": self.validity_window_start_utc,
            "validity_window_end_utc": self.validity_window_end_utc,
            "epoch_locked_closure_residual_kms": _maybe(self.epoch_locked_closure_residual_kms),
            "epoch_locked_flyby_continuity_max_dv_kms": _maybe(
                self.epoch_locked_flyby_continuity_max_dv_kms
            ),
            "epoch_locked_independent_check_residual_kms": _maybe(
                self.epoch_locked_independent_check_residual_kms
            ),
            "closure_residual_kms": _maybe(self.closure_residual_kms),
            "flyby_continuity_max_dv_kms": _maybe(self.flyby_continuity_max_dv_kms),
            "independent_cross_check_residual_kms": _maybe(
                self.independent_cross_check_residual_kms
            ),
            "powered_axis_record": (
                self.powered_axis_record.as_dict() if self.powered_axis_record is not None else None
            ),
            "notes": list(self.notes),
        }


# ---------------------------------------------------------------------------
# Joint search driver
# ---------------------------------------------------------------------------


def _is_planar_cell(z0_amplitude_nondim: float) -> bool:
    """True when the Axis-C coord is below a planar-equivalence epsilon."""
    return abs(float(z0_amplitude_nondim)) < 1e-9


def _add_days_to_iso_utc(utc_iso: str, days: float) -> str:
    """Add ``days`` to an ISO-8601 UTC string; thin wrapper for testability."""
    from datetime import datetime, timedelta

    if utc_iso.endswith("Z"):
        base = datetime.fromisoformat(utc_iso[:-1])
    else:
        base = datetime.fromisoformat(utc_iso)
    return (base + timedelta(days=float(days))).strftime("%Y-%m-%dT%H:%M:%S")


def _run_epoch_locked_witness(
    powered: LowThrustCyclerCandidate,
    *,
    launch_epoch_utc: str,
    ephem: Any,
) -> tuple[float | None, float | None, float | None, str]:
    """Run the #289 epoch-locked closure as a same-tour witness.

    Returns ``(closure_residual_kms, flyby_continuity_max_dv_kms,
    independent_check_residual_kms, note)``. On any failure (Lambert
    singularity, astropy unavailable, etc.) all three fields are ``None`` and
    the note explains.
    """
    try:
        from cyclerfinder.genome.epoch_aware_genome import (
            EpochLockedTrajectory,
            close_epoch_locked,
        )
    except ImportError as exc:  # pragma: no cover - environment-dependent
        return None, None, None, f"epoch_locked_unavailable:{exc}"

    leg_tofs_days = tuple(float(t) for t in powered.leg_tofs_days)
    if not leg_tofs_days:
        return None, None, None, "epoch_locked_skipped_no_legs"

    total_days = float(sum(leg_tofs_days))
    try:
        end_utc = _add_days_to_iso_utc(launch_epoch_utc, total_days)
    except Exception as exc:  # pragma: no cover
        return None, None, None, f"epoch_window_compute_failed:{exc}"

    seq = tuple(powered.sequence)
    n_enc = len(seq)
    vinf_kms = tuple(powered.vinf_per_encounter_kms)
    if len(vinf_kms) != n_enc:
        # Powered driver reports a vinf per body in the closed sequence; if
        # there's a length mismatch, pad/clip defensively.
        _pad_val = vinf_kms[-1] if vinf_kms else 0.0
        _padding = [_pad_val] * (n_enc - len(vinf_kms))
        vinf_kms = tuple(list(vinf_kms)[:n_enc] + _padding)

    try:
        spec = EpochLockedTrajectory(
            sequence=seq,
            leg_tofs_days=leg_tofs_days,
            vinf_kms_at_encounters=vinf_kms,
            launch_epoch_utc=launch_epoch_utc,
            orbit_class="quasi_cycler",
            n_returns=1,
            validity_window_start_utc=launch_epoch_utc,
            validity_window_end_utc=end_utc,
        )
    except ValueError as exc:
        return None, None, None, f"epoch_locked_spec_invalid:{exc}"

    try:
        closure = close_epoch_locked(
            spec,
            ephem,
            closure_tol_kms=10.0,  # Phase 1 witness - tolerate large
            flyby_continuity_tol_kms=10.0,  # disagreement; record, don't gate.
            independent_cross_check=True,
            independent_tol_kms=10.0,
        )
    except Exception as exc:
        return None, None, None, f"epoch_locked_closure_failed:{type(exc).__name__}:{exc}"

    return (
        float(closure.closure_residual_kms),
        float(closure.flyby_continuity_max_dv_kms),
        (
            None
            if closure.independent_check_residual_kms is None
            else float(closure.independent_check_residual_kms)
        ),
        "epoch_locked_witness_ran",
    )


def joint_axis_search(
    primary_sequence: tuple[str, ...],
    *,
    k_synodic: int = 1,
    ephem: Ephemeris | None = None,
    # Per-axis configuration
    powered_budgets_kms: Sequence[float] = _DEFAULT_POWERED_BUDGETS_KMS,
    n_revs_grid_per_leg: Sequence[Sequence[int]] | None = None,
    z0_amplitudes_nondim: Sequence[float] = _DEFAULT_Z0_AMPLITUDES_NONDIM,
    launch_epoch_grid: Sequence[str | None] = _DEFAULT_LAUNCH_EPOCH_GRID,
    # Shared - forwarded to the #309 powered driver
    leg_tof_guesses_days: Sequence[float] | None = None,
    leg_tof_bounds_days: Sequence[tuple[float, float]] | None = None,
    tof_jitter_half_days: Sequence[float] | None = None,
    synodic_pair: tuple[str, str] | None = None,
    closure_body: str | None = None,
    closure_flyby_alt_km: float | None = None,
    closure_floor_kms: float = 0.5,
    flyby_continuity_floor_kms: float = 0.05,
    # Powered-axis defaults forwarded
    isp_s: float = 3000.0,
    tmax_kn: float = 2.5e-4,
    m0_kg: float = 1.0e4,
    dry_mass_kg: float = 5.0e3,
    n_segments: int = 20,
    search: SearchFn | None = None,
    n_starts: int = 4,
    seed: int = 0,
    epoch_locked_ephem: Any | None = None,
    t0_guess_sec_for_epoch_blind: float | None = None,
) -> Iterator[JointAxisCandidate]:
    """Joint axis search over (powered x multi-rev x 3D x epoch) for one tour.

    Phase 1 = full 4D Cartesian product over the per-axis grids. Phase 2+ will
    add smarter sampling (Sobol, MCMC, surrogate-driven). Phase 1 IS the
    capability; Phase 2 is the efficiency.

    For each (powered_dv, n_revs, z0, launch_epoch) cell:

      1. Build the per-leg ``per_leg_revs`` tuple from the Axis-B grid.
      2. Run the #309 powered driver
         (:func:`search_low_thrust_cyclers`) with the cell's per-leg revs +
         the requested propulsion + Tsiolkovsky wrap. The driver computes
         the maintenance dV optimum + Sims-Flanagan feasibility witness +
         seed-bump independent cross-check.
      3. If ``z0 > 0`` (Axis-C non-trivial): RECORD the amplitude on the
         candidate; Phase 1 does NOT yet drive the 3D corrector inline.
      4. If ``launch_epoch_utc is not None`` (Axis-D non-trivial): run the
         #289 :func:`close_epoch_locked` driver as a witness.

    A cell *survives* and emits a :class:`JointAxisCandidate` iff the powered
    driver closes (returns non-``None``). The closure_floor_kms and
    flyby_continuity_floor_kms thresholds are recorded on the candidate's
    notes; Phase 1 does NOT pre-filter on them.

    Parameters
    ----------
    primary_sequence:
        Encounter tour, must close (first == last). E.g. ``("E", "M", "E")``.
    k_synodic:
        Cycler period in synodic multiples of the ``synodic_pair`` pair.
    ephem:
        Planet-state provider. Defaults to ``Ephemeris("circular")``.
    powered_budgets_kms:
        Axis-A grid. ``(0.0,)`` (the default) means pure ballistic; values
        > 0 record the cell as a powered request.
    n_revs_grid_per_leg:
        Axis-B grid. Outer tuple length = n_legs; each inner tuple is the
        revs grid for that leg. ``None`` defaults to ``((0,),) * n_legs``.
    z0_amplitudes_nondim:
        Axis-C grid. ``(0.0,)`` (the default) is planar; positive values
        record a 3D request (Phase 1: record-only).
    launch_epoch_grid:
        Axis-D grid. ``(None,)`` (the default) is epoch-blind; UTC strings
        activate the #289 epoch-locked witness for that cell.
    leg_tof_guesses_days, leg_tof_bounds_days, tof_jitter_half_days,
    synodic_pair, closure_body, closure_flyby_alt_km:
        Forwarded to the #309 powered driver.
    closure_floor_kms, flyby_continuity_floor_kms:
        Phase 1 ADVISORY floors (recorded on the candidate's notes; not used
        to gate). Phase 2 will tighten into hard gates.
    isp_s, tmax_kn, m0_kg, dry_mass_kg, n_segments:
        Powered-axis configuration; forwarded to the #309 driver.
    search:
        Optional ``SearchFn`` for the #309 literature-check.
    n_starts, seed:
        Forwarded to the #309 optimiser.
    epoch_locked_ephem:
        Ephemeris for the #289 witness. If ``None``, reuses ``ephem``. Pass
        ``Ephemeris("astropy")`` for a real-DE440 epoch-locked closure.
    t0_guess_sec_for_epoch_blind:
        For epoch-blind cells (``launch_epoch_grid`` entry is ``None``), an
        explicit departure-epoch seed (s since J2000) forwarded to the #309
        optimiser. ``None`` lets the powered driver pick its own default
        (2030-01-01). For Aldrin reproduction the value should be the
        phase-inversion seed (``maintain._default_t0_guess(em_tof_days)``).

    Yields
    ------
    JointAxisCandidate
        One candidate per *surviving* cell (powered driver converged).
        Non-convergent cells are dropped silently - non-closure is a
        non-event.
    """
    if primary_sequence[0] != primary_sequence[-1]:
        raise ValueError(
            f"primary_sequence must close (first body == last body); "
            f"got {list(primary_sequence)!r}",
        )
    if k_synodic < 1:
        raise ValueError(f"k_synodic must be >= 1; got {k_synodic}")

    eph = ephem if ephem is not None else Ephemeris("circular")
    eph_epoch = epoch_locked_ephem if epoch_locked_ephem is not None else eph

    n_legs = len(primary_sequence) - 1
    if n_revs_grid_per_leg is None:
        n_revs_grid_per_leg = tuple((0,) for _ in range(n_legs))
    n_revs_grid_per_leg = tuple(tuple(int(r) for r in g) for g in n_revs_grid_per_leg)
    if len(n_revs_grid_per_leg) != n_legs:
        raise ValueError(
            f"n_revs_grid_per_leg length {len(n_revs_grid_per_leg)} != n_legs {n_legs}",
        )

    powered_budgets = tuple(float(p) for p in powered_budgets_kms)
    z0_amps = tuple(float(z) for z in z0_amplitudes_nondim)
    launch_epochs = tuple(launch_epoch_grid)

    def _per_leg_revs_combinations(
        grids: Sequence[Sequence[int]],
    ) -> Iterator[tuple[int, ...]]:
        if not grids:
            yield ()
            return
        head = grids[0]
        for rest in _per_leg_revs_combinations(grids[1:]):
            for r in head:
                yield (int(r), *rest)

    seen_fingerprints: set[tuple[Any, ...]] = set()

    for revs_tuple in _per_leg_revs_combinations(n_revs_grid_per_leg):
        for powered_kms in powered_budgets:
            for z0_amp in z0_amps:
                for epoch_utc in launch_epochs:
                    t0_guess_sec: float | None
                    if epoch_utc is not None:
                        try:
                            from datetime import datetime as _dt

                            if epoch_utc.endswith("Z"):
                                _base = _dt.fromisoformat(epoch_utc[:-1])
                            else:
                                _base = _dt.fromisoformat(epoch_utc)
                            t0_guess_sec = datetime_to_sec_since_j2000(_base)
                        except Exception:
                            t0_guess_sec = None
                    else:
                        t0_guess_sec = t0_guess_sec_for_epoch_blind

                    try:
                        powered_record = search_low_thrust_cyclers(
                            sequence=primary_sequence,
                            k_synodic=k_synodic,
                            ephem=eph,
                            t0_guess_sec=t0_guess_sec,
                            leg_tof_guesses_days=leg_tof_guesses_days,
                            leg_tof_bounds_days=leg_tof_bounds_days,
                            per_leg_revs=revs_tuple,
                            tof_jitter_half_days=tof_jitter_half_days,
                            synodic_pair=synodic_pair,
                            closure_body=closure_body,
                            closure_flyby_alt_km=closure_flyby_alt_km,
                            isp_s=isp_s,
                            tmax_kn=tmax_kn,
                            m0_kg=m0_kg,
                            dry_mass_kg=dry_mass_kg,
                            n_segments=n_segments,
                            search=search,
                            n_starts=n_starts,
                            seed=seed,
                        )
                    except (RuntimeError, ValueError):
                        continue

                    if powered_record is None:
                        continue

                    cell_notes: list[str] = []
                    if not _is_planar_cell(z0_amp):
                        cell_notes.append(
                            f"axis_C_3d_request_recorded:z0={z0_amp:.3e}_"
                            f"phase2_pipe_to_cr3bp_3d_corrector",
                        )

                    epoch_closure_kms: float | None = None
                    epoch_flyby_kms: float | None = None
                    epoch_indep_kms: float | None = None
                    validity_start: str | None = None
                    validity_end: str | None = None
                    if epoch_utc is not None:
                        (
                            epoch_closure_kms,
                            epoch_flyby_kms,
                            epoch_indep_kms,
                            note,
                        ) = _run_epoch_locked_witness(
                            powered_record,
                            launch_epoch_utc=epoch_utc,
                            ephem=eph_epoch,
                        )
                        cell_notes.append(note)
                        validity_start = epoch_utc
                        if powered_record.leg_tofs_days:
                            try:
                                validity_end = _add_days_to_iso_utc(
                                    epoch_utc,
                                    float(sum(powered_record.leg_tofs_days)),
                                )
                            except Exception:
                                validity_end = None

                    if powered_record.closure_residual_kms > closure_floor_kms:
                        cell_notes.append(
                            f"closure_residual_above_advisory_floor:"
                            f"{powered_record.closure_residual_kms:.3f}>{closure_floor_kms:.3f}",
                        )
                    if powered_kms > 0.0 and powered_record.maintenance_dv_kms > powered_kms:
                        cell_notes.append(
                            f"maintenance_dv_exceeds_cell_budget:"
                            f"{powered_record.maintenance_dv_kms:.3f}>{powered_kms:.3f}",
                        )

                    cand = JointAxisCandidate(
                        sequence=tuple(primary_sequence),
                        vinf_tuple_kms=tuple(powered_record.vinf_per_encounter_kms),
                        leg_tofs_days=tuple(powered_record.leg_tofs_days),
                        powered_maintenance_dv_kms_per_synodic=float(
                            powered_record.maintenance_dv_kms
                        ),
                        powered_budget_kms_requested=float(powered_kms),
                        powered_leg_indices=tuple(range(n_legs)),
                        sims_flanagan_feasible=bool(powered_record.sims_flanagan_feasible),
                        n_revs_per_leg=tuple(int(r) for r in revs_tuple),
                        z0_amplitude_nondim=float(z0_amp),
                        launch_epoch_utc=epoch_utc,
                        validity_window_start_utc=validity_start,
                        validity_window_end_utc=validity_end,
                        epoch_locked_closure_residual_kms=epoch_closure_kms,
                        epoch_locked_flyby_continuity_max_dv_kms=epoch_flyby_kms,
                        epoch_locked_independent_check_residual_kms=epoch_indep_kms,
                        closure_residual_kms=float(powered_record.closure_residual_kms),
                        flyby_continuity_max_dv_kms=(
                            epoch_flyby_kms
                            if epoch_flyby_kms is not None
                            else float(powered_record.maintenance_dv_kms)
                        ),
                        independent_cross_check_residual_kms=float(
                            powered_record.independent_cross_check_residual_kms
                        ),
                        powered_axis_record=powered_record,
                        notes=tuple(cell_notes),
                    )

                    fp = (
                        tuple(revs_tuple),
                        round(float(powered_kms), 6),
                        round(float(z0_amp), 6),
                        epoch_utc,
                        tuple(round(float(v), 2) for v in cand.vinf_tuple_kms),
                        tuple(round(float(t), 0) for t in cand.leg_tofs_days),
                    )
                    if fp in seen_fingerprints:
                        continue
                    seen_fingerprints.add(fp)
                    yield cand


__all__ = [
    "JointAxisCandidate",
    "joint_axis_search",
]
