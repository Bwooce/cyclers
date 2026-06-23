"""Precursor MGA matcher — find Earth-launched chains that insert into a cycler.

Phase 4 of task #289 (precursor-MGA matcher).  This module *wraps* the Phase
1-3 substrate (epoch-aware genome, Tisserand-Poincaré BFS, multi-shell, TOF
optimisation, optional DSM) and answers a single concrete question:

    Given an existing cycler row in ``data/catalogue.yaml``, what one-shot
    Earth-launched MGA chains end with the same ``|V_inf|`` at the same
    first cycler-encounter body — i.e. could insert a spacecraft FROM Earth
    ON to that cycler?

The question matters because published cycler papers (Aldrin, Russell-Ocampo,
McConaghy, Braik-Ross, Roberts-Tsoukkas) typically document the steady-state
cycler invariants but only rarely the precursor insertion trajectory.  The
insertion trajectory is therefore a low-corpus-coverage region of the
literature where genuine literature-fresh candidates are plausible — and the
gauntlet for promoting one to the catalogue requires (a) closure under the
Phase-1 closure driver, (b) a literature-novelty check that does not find
the chain in the corpus, and (c) a separate V0-V5 round (not done here).

DISCIPLINE
----------
* **No catalogue writeback** — even literature-fresh candidates only land in
  ``data/precursor_302_*.jsonl`` (the discovery probe output).  Adapting the
  V0-V5 gauntlet to epoch-locked classes is a separate future task.
* **No novelty claims** — a literature-fresh hit per ``check_literature`` is
  necessary-not-sufficient (see ``feedback_literature_novelty_check_baseline``
  and the discipline preamble on :mod:`literature_check`).  The matcher
  reports "literature-fresh" only; the caller decides what to do with it.
* **READ-ONLY on the Phase 1-3 modules** — this module wraps them, never
  edits them.  If the matcher genuinely needs an extension to
  :mod:`epoch_aware_genome` or :mod:`tisserand_mga_window` that is a sign
  that the Phase-5 task (#303) is needed.

Methodology
-----------
1.  Read the target cycler's first encounter body and its published
    ``|V_inf|`` at that body from the catalogue row's
    ``vinf_kms_at_encounters``.  This is the "seed V_inf" that any
    insertion chain must hit.
2.  Enumerate Earth-launched MGA chains via
    :func:`cyclerfinder.search.tisserand_mga_window.find_mga_chains` whose
    terminal body is the cycler's first encounter body and whose terminal
    V_inf bin is close to the seed V_inf (within
    ``vinf_terminal_tol_kms``).
3.  Validate each chain through :func:`close_epoch_locked` (Phase 1) and
    optionally optimise the per-leg TOFs and launch epoch via
    :func:`optimise_chain_tofs` (Phase 3).  Only converged closures pass.
4.  Run :func:`check_literature` on each surviving candidate's structural
    fingerprint.  Tag the result on the :class:`PrecursorMatch` record.

The matcher returns ranked :class:`PrecursorMatch` records sorted by a
combined (V_inf-match-residual, flyby-continuity, closure-residual)
score.  Lower is better.

See ``docs/notes/2026-06-16-302-289-phase4-precursor-matcher.md`` for the
Aldrin / S1L1 scan verdicts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from cyclerfinder.genome.epoch_aware_genome import (
    EpochLockedClosure,
    EpochLockedTrajectory,
    close_epoch_locked,
)
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    LiteratureCheckResult,
    SearchFn,
    check_literature,
)
from cyclerfinder.search.tisserand_mga_window import (
    MGAChainCandidate,
    _add_days_utc,
    find_mga_chains,
    optimise_chain_tofs,
)

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.data.catalog import Catalog, CatalogueEntry


# --------------------------------------------------------------------------- #
# Public dataclass
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PrecursorMatch:
    """One precursor candidate that inserts into a catalogue cycler row.

    A surviving candidate has:

    * a closed-form epoch-locked trajectory (the ``candidate``),
    * the catalogue id of the cycler row this precursor inserts into,
    * the cycler's published seed V_inf at the FIRST encounter body (the
      EXPECTED side of the V_inf match residual),
    * the residual ``|terminal_vinf - seed_vinf|`` for this candidate,
    * the launch-epoch alignment score (a continuous number representing
      proximity of the candidate's TERMINAL epoch to the cycler's
      windowed-cadence — 0.0 for "no scoring info available", in which
      case the alignment is a heuristic only),
    * the Phase-1 closure record (residual + per-encounter V_inf +
      independent cross-check),
    * the literature-novelty verdict from :func:`check_literature`.
    """

    candidate: EpochLockedTrajectory
    cycler_id: str
    cycler_seed_vinf_kms: float
    vinf_match_residual_kms: float
    epoch_alignment_score: float
    closure: EpochLockedClosure
    literature_check: LiteratureCheckResult

    def is_literature_fresh(self) -> bool:
        """``True`` iff the literature-check verdict is a clean not-found.

        Necessary-not-sufficient (see the
        :mod:`cyclerfinder.search.literature_check` discipline preamble).
        A fresh hit is a CANDIDATE for the V0-V5 gauntlet, NOT a novelty
        claim.
        """
        return self.literature_check.status == "not-found"

    def quality_score(self) -> float:
        """Combined ranking score (lower is better).

        Combines

          * the V_inf match residual at the cycler's seed body,
          * the closure residual (worst per-encounter V_inf residual),
          * 2x the flyby-continuity max DV (matches the
            :func:`optimise_chain_tofs` loss weighting; ballistic
            continuity is a harder constraint than V_inf match).

        Used to sort the matcher output.  An optimiser's ``accept_loss_kms``
        is NOT folded into this — quality is a ranking signal, not a gate.
        """
        return float(
            self.vinf_match_residual_kms
            + self.closure.closure_residual_kms
            + 2.0 * self.closure.flyby_continuity_max_dv_kms
        )


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


_BODY_CODE_FOR_CATALOGUE = {
    # The catalogue uses one-letter codes for the inner planets, matching the
    # Phase-1 / Phase-2 substrate.  This helper is a no-op for the inner-planet
    # catalogue rows (Aldrin, S1L1, Russell-Ocampo, McConaghy) but is kept
    # explicit so a future moon-tour precursor extension has an obvious place
    # to wire body-code translation.
    "E": "E",
    "M": "M",
    "V": "V",
    "Me": "Me",
    "J": "J",
    "S": "S",
}


def _first_encounter_body_and_vinf(entry: CatalogueEntry) -> tuple[str, float]:
    """Extract the cycler's first encounter body + sourced V_inf from a row.

    The "first encounter" of a steady-state cycler in the catalogue's
    ``vinf_kms_at_encounters`` ordering is canonically the body where the
    spacecraft arrives after a precursor insertion — i.e. the natural
    insertion target.  Aldrin's first encounter is Earth (V_inf 6.5 km/s);
    S1L1's first encounter is Earth (V_inf ~4.9-5.0 km/s in the coplanar
    single-ellipse framing; the catalogue currently stores 5.65 km/s for
    s1l1-2syn-em-cpom, but #366 (2026-06-17) identified this as a
    cycler-confusion artifact — 5.65 traces to McConaghy 2004 Table 4
    row 2L3 / Byrnes 2002 Case 1, a DIFFERENT cycler.  V0→V1 promotion
    (#365) will correct the row; this matcher reads whatever value the
    row carries at runtime, so the precursor target moves with the row).
    Both are insertion-from-Earth-launch candidates.

    Raises
    ------
    ValueError
        If the row has no usable ``vinf_kms_at_encounters`` block.
    """
    if not entry.vinf_kms_at_encounters:
        raise ValueError(
            f"cycler {entry.id!r} has no vinf_kms_at_encounters block; "
            "cannot anchor a precursor V_inf match"
        )
    body_raw, vinf_raw = entry.vinf_kms_at_encounters[0]
    if vinf_raw is None:
        raise ValueError(
            f"cycler {entry.id!r}'s first encounter ({body_raw!r}) has no sourced "
            "V_inf; cannot anchor a precursor V_inf match"
        )
    body = _BODY_CODE_FOR_CATALOGUE.get(str(body_raw))
    if body is None:
        raise ValueError(
            f"cycler {entry.id!r}'s first encounter body {body_raw!r} is not in "
            f"the matcher's body-code translation table {sorted(_BODY_CODE_FOR_CATALOGUE)!r}"
        )
    return body, float(vinf_raw)


def _terminal_vinf_kms(closure: EpochLockedClosure) -> float:
    """The closure's |V_inf| at the FINAL encounter (the cycler-insertion body).

    For an MGA-tour candidate the last per-encounter V_inf is the V_inf at
    arrival at the chain's terminal body — which the matcher is comparing
    against the cycler's seed V_inf.
    """
    return float(closure.per_encounter_vinf_kms[-1])


def _epoch_alignment_score(
    launch_epoch_utc: str,
    target_phase_window_utc: tuple[str, str] | None,
) -> float:
    """Heuristic alignment score for the candidate's launch epoch.

    If ``target_phase_window_utc`` is provided (a (t0, t1) window inside which
    the cycler's first encounter would naturally occur), the score is

      * 1.0 if the launch epoch sits inside the window
      * a linear fall-off to 0 over one extra synodic period (taken as 779.8
        days for E-M, the conservative default)

    If no target window is provided (the typical case, since the catalogue
    rarely publishes the cycler's seed_epoch), the score is a constant 0.0
    — the matcher reports this as "no alignment signal available" and the
    caller must not gate on it.  The matcher's RANKING is dominated by the
    V_inf match residual + closure quality; alignment is informational
    only.
    """
    if target_phase_window_utc is None:
        return 0.0
    from astropy.time import Time

    t_launch = Time(launch_epoch_utc, scale="utc")
    t0 = Time(target_phase_window_utc[0], scale="utc")
    t1 = Time(target_phase_window_utc[1], scale="utc")
    if t0 <= t_launch <= t1:
        return 1.0
    if t_launch < t0:
        gap_days = float((t0 - t_launch).to("day").value)
    else:
        gap_days = float((t_launch - t1).to("day").value)
    # Linear fall-off over one Earth-Mars synodic period (779.8 d) — a
    # conservative default that matches the published Aldrin / S1L1 cadence.
    fall_off_days = 779.8
    return max(0.0, 1.0 - gap_days / fall_off_days)


def _phase_window_for_entry(
    entry: CatalogueEntry,
    override: tuple[str, str] | None,
) -> tuple[str, str] | None:
    """Resolve the cycler-cadence target phase window (#307 Task 3).

    An explicit ``override`` always wins. Otherwise, auto-derive the window from
    the cycler row's published ``validity_window`` (the epoch span over which the
    cycler is insertion-ready) — the natural target window for the precursor's
    terminal-Earth arrival. Returns ``None`` when neither is available (the row
    publishes no validity window — the common case; the alignment score then
    stays the informational, ungated 0.0 per :func:`_epoch_alignment_score`).

    A *period* alone cannot anchor a phase window (it gives the cadence interval,
    not the phase), so only an epoch-anchored ``validity_window`` is used.
    """
    if override is not None:
        return override
    raw = getattr(entry, "raw", None)
    if not isinstance(raw, dict):
        return None
    vw = raw.get("validity_window")
    if isinstance(vw, dict):
        start, end = vw.get("start"), vw.get("end")
        if start and end:
            return (str(start), str(end))
    return None


def _signature_from_candidate(
    candidate: EpochLockedTrajectory,
) -> CandidateSignature:
    """Build the literature-check signature for a precursor candidate."""
    return CandidateSignature(
        primary="Sun",  # heliocentric MGA precursor — Phase 4 scope
        sequence=candidate.sequence,
        period_k=None,  # epoch-locked: not a periodic cycler
        period_years=None,
        vinf_per_encounter_kms=candidate.vinf_kms_at_encounters,
        n_rev=(),
    )


def _candidate_to_trajectory(
    candidate: MGAChainCandidate,
    inserts_into: str,
    *,
    periapsis_altitudes_km: tuple[float | None, ...] | None,
) -> EpochLockedTrajectory | None:
    """Wrap an :class:`MGAChainCandidate` as a ``precursor_mga`` trajectory.

    Returns ``None`` if the resulting :class:`EpochLockedTrajectory`'s
    invariants fail (e.g. an unknown body code surfaces).
    """
    total_tof_days = sum(candidate.leg_tofs_days)
    end_utc = _add_days_utc(candidate.launch_epoch_utc, total_tof_days)
    try:
        return EpochLockedTrajectory(
            sequence=candidate.sequence,
            leg_tofs_days=candidate.leg_tofs_days,
            vinf_kms_at_encounters=candidate.vinf_tuple_kms,
            launch_epoch_utc=candidate.launch_epoch_utc,
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc=candidate.launch_epoch_utc,
            validity_window_end_utc=end_utc,
            inserts_into=inserts_into,
            periapsis_altitudes_km=periapsis_altitudes_km,
            notes="Precursor MGA matcher candidate (#302)",
        )
    except ValueError:
        return None


# --------------------------------------------------------------------------- #
# Public surface
# --------------------------------------------------------------------------- #


def find_cycler_precursors(
    cycler_id: str,
    catalogue: Catalog,
    ephemeris: Ephemeris,
    *,
    launch_window: tuple[str, str],
    max_legs: int = 4,
    intermediate_bodies: tuple[str, ...] = ("V", "E"),
    vinf_terminal_tol_kms: float = 0.5,
    vinf_terminal_post_closure_tol_kms: float | None = None,
    vinf_grid_kms: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (60.0, 600.0),
    epoch_step_days: float = 30.0,
    tof_optimise: bool = True,
    closure_tol_kms: float = 1.0,
    flyby_continuity_tol_kms: float = 0.10,
    independent_cross_check: bool = True,
    independent_tol_kms: float = 0.2,
    target_phase_window_utc: tuple[str, str] | None = None,
    multi_shell: bool = True,
    pump_envelope_factor: float = 1.0,
    a_range_au: tuple[float, float] = (0.3, 5.0),
    max_candidates_to_validate: int | None = 200,
    literature_check_search: SearchFn | None = None,
    progress_hook: object | None = None,
    max_revs: int = 0,
) -> list[PrecursorMatch]:
    """Find precursor MGA chains that insert a spacecraft into a cycler row.

    Reads the cycler's first encounter body and seed V_inf from the
    catalogue row (``catalogue.by_id[cycler_id].vinf_kms_at_encounters[0]``);
    enumerates Earth-launched chains via
    :func:`find_mga_chains` whose terminal body matches the cycler's first
    encounter body and whose terminal V_inf bin is within
    ``vinf_terminal_tol_kms`` of the seed V_inf; closes each via
    :func:`close_epoch_locked` (optionally via :func:`optimise_chain_tofs`
    when ``tof_optimise`` is on); runs :func:`check_literature` on the
    survivors.

    Parameters
    ----------
    cycler_id:
        Catalogue id of the target cycler (e.g.
        ``"aldrin-classic-em-k1-outbound"``).
    catalogue:
        Parsed :class:`Catalog` from
        :func:`cyclerfinder.data.catalog.load_catalog`.  READ-ONLY here.
    ephemeris:
        Phase-1 body-state provider.  Real-ephemeris (``Ephemeris("astropy")``)
        is the typical case.
    launch_window:
        ``(t0_utc, t1_utc)`` window for the precursor's first (Earth) departure.
        Typically one Earth-Mars synodic period (~26 months) for full coverage.
    max_legs:
        Maximum chain length in legs.  4 admits VEEGA-style triple-flyby
        precursors; lower keeps the search box small.
    intermediate_bodies:
        Body codes admissible at intermediate flybys.  The first body is
        always ``E`` (Earth launch); the last body is the cycler's first
        encounter; intermediate bodies are drawn from this tuple.  Default
        ``("V", "E")`` covers Venus and Earth resonant returns — the
        standard precursor archetype.
    vinf_terminal_tol_kms:
        Tolerance on the terminal-body V_inf match vs the cycler's seed V_inf.
        Default 0.5 km/s.
    vinf_grid_kms:
        V_inf grid used by :func:`find_mga_chains`.  Default covers 3-8 km/s
        which contains both Aldrin (~6.5) and S1L1 (~4.9-5.0 coplanar
        single-ellipse; #366 retraction note in _first_encounter_body_and_vinf
        docstring) seed V_inf.
    tof_box_days_per_leg:
        Per-leg TOF constraint window.  Default ``(60, 600)`` covers single-
        synodic Earth-Mars Hohmanns (~258 days) and short Venus-Earth hops
        (~90 days).
    epoch_step_days:
        Launch-epoch grid step.  Default 30 days; the optimiser refines.
    tof_optimise:
        If ``True``, run :func:`optimise_chain_tofs` on each Phase-1
        ballistic candidate.  Halves the closure residual at ~3-10x cost
        per candidate.
    closure_tol_kms, flyby_continuity_tol_kms, independent_cross_check,
    independent_tol_kms:
        Forwarded to :func:`close_epoch_locked`.  Defaults relax the
        Phase-1 gates slightly (1.0 vs 0.5 km/s closure, 0.10 vs 0.05 km/s
        continuity) because precursor matching is a search probe — we
        want candidates to survive into the doc with their residuals
        recorded, not to be silently dropped at a tight gate.
    target_phase_window_utc:
        Optional (t0, t1) UTC window inside which the cycler's first
        encounter would naturally occur.  Used by
        :func:`_epoch_alignment_score`.  ``None`` (default) means no
        alignment signal is available — the matcher's ranking is
        dominated by V_inf match + closure quality.
    multi_shell:
        Forwarded to :func:`find_mga_chains`.  Default ``True`` admits
        per-flyby V_inf shifts within the Strange-Longuski pump envelope
        — needed for any precursor where the Earth launch V_inf differs
        from the cycler insertion V_inf (the typical case).
    pump_envelope_factor:
        Forwarded to :func:`find_mga_chains`.
    a_range_au:
        Forwarded to :func:`find_mga_chains` and the Tisserand 3-D
        predicate.  Default ``(0.3, 5.0)`` covers inner-planet + Jupiter.
    max_candidates_to_validate:
        Cap on the number of geometric proposals to validate (per
        ``feedback_incremental_progress_reports`` — we don't want a 2h
        black-box scan).  ``None`` validates the entire stream.
    literature_check_search:
        Injected :data:`SearchFn` for :func:`check_literature`.  ``None``
        (default) skips the literature check and tags every survivor as
        ``inconclusive`` — useful for fast unit tests; the production
        scan scripts pass an offline-corpus fake_search wrapped from
        :data:`KNOWN_CORPUS`.
    progress_hook:
        Optional ``callable(stage_message: str) -> None`` invoked at
        coarse milestones (enumeration start, validation start, each
        100 candidates validated, completion).  ``None`` is silent.
        Per ``feedback_incremental_progress_reports`` — never run a black
        box.

    Returns
    -------
    list[PrecursorMatch]
        Ranked by :meth:`PrecursorMatch.quality_score` ascending (best
        first).  Empty list if no chain survives the closure + V_inf
        match gates.
    """
    if cycler_id not in catalogue.by_id:
        raise ValueError(
            f"cycler_id {cycler_id!r} not in catalogue (by_id has {len(catalogue.by_id)} entries)"
        )
    if vinf_terminal_tol_kms <= 0.0:
        raise ValueError(f"vinf_terminal_tol_kms must be positive, got {vinf_terminal_tol_kms}")
    if "E" not in intermediate_bodies and "E" != "E":
        # Earth is always implicitly the launch body, regardless of
        # intermediate_bodies; this branch is kept for parameter clarity.
        pass

    entry = catalogue.by_id[cycler_id]
    target_body, seed_vinf = _first_encounter_body_and_vinf(entry)
    # #307 Task 3: auto-derive the cycler-cadence target phase window from the
    # row's validity_window when the caller does not override it (informational
    # alignment score only — never gated; see _epoch_alignment_score).
    target_phase_window_utc = _phase_window_for_entry(entry, target_phase_window_utc)

    # The matcher's planet_set is the launch body (Earth) + intermediate
    # bodies + target body.  Earth is always included (we launch from
    # Earth); target is always included (we arrive at it).  Intermediate
    # bodies are drawn from the caller's tuple.
    planet_set_set = {"E", target_body} | set(intermediate_bodies)
    planet_set = tuple(sorted(planet_set_set))

    def _emit(msg: str) -> None:
        if callable(progress_hook):
            progress_hook(msg)

    _emit(
        f"precursor matcher: target={cycler_id} "
        f"first_encounter={target_body} seed_vinf={seed_vinf:.3f} km/s "
        f"planet_set={planet_set}"
    )

    # ------------------------------------------------------------- step 1
    # Enumerate geometric chain proposals.  We only retain chains that
    # START at Earth (precursor) and END at the cycler's first encounter
    # body within the V_inf tolerance.
    candidates_filtered: list[MGAChainCandidate] = []
    for cand in find_mga_chains(
        launch_window=launch_window,
        planet_set=planet_set,
        max_legs=max_legs,
        vinf_grid_kms=vinf_grid_kms,
        tof_box_days_per_leg=tof_box_days_per_leg,
        epoch_step_days=epoch_step_days,
        start_body_filter=("E",),
        a_range_au=a_range_au,
        multi_shell=multi_shell,
        pump_envelope_factor=pump_envelope_factor,
    ):
        if cand.sequence[-1] != target_body:
            continue
        if cand.sequence[0] != "E":
            continue
        if len(cand.sequence) < 2:
            continue
        terminal_vinf_bin = cand.vinf_tuple_kms[-1]
        if abs(terminal_vinf_bin - seed_vinf) > vinf_terminal_tol_kms:
            continue
        candidates_filtered.append(cand)
        if (
            max_candidates_to_validate is not None
            and len(candidates_filtered) >= max_candidates_to_validate
        ):
            break

    _emit(
        f"precursor matcher: enumerated {len(candidates_filtered)} geometric "
        f"proposals matching target={target_body} within "
        f"|V_inf| - {seed_vinf:.2f} <= {vinf_terminal_tol_kms} km/s"
    )

    # ------------------------------------------------------------- step 2
    # Validate each proposal via Phase 1 closure (optionally TOF-optimised).
    matches: list[PrecursorMatch] = []
    for idx, cand in enumerate(candidates_filtered):
        if idx and idx % 50 == 0:
            _emit(
                f"precursor matcher: validated {idx} / {len(candidates_filtered)} "
                f"survivors={len(matches)}"
            )

        if tof_optimise:
            opt_result = optimise_chain_tofs(
                cand,
                ephemeris,
                orbit_class="precursor_mga",
                n_returns=1,
                inserts_into=cycler_id,
                periapsis_altitudes_km=None,
                max_iter=40,
                epoch_search_half_width_days=15.0,
                tof_search_relative_half_width=0.30,
                alpha_flyby_continuity=2.0,
                accept_loss_kms=None,
                independent_cross_check=False,
                max_revs=max_revs,
            )
            if opt_result is None:
                continue
            opt_cand, closure, _opt_loss = opt_result
            cand_for_traj = opt_cand
        else:
            traj = _candidate_to_trajectory(
                cand,
                inserts_into=cycler_id,
                periapsis_altitudes_km=None,
            )
            if traj is None:
                continue
            try:
                closure = close_epoch_locked(
                    traj,
                    ephemeris,
                    closure_tol_kms=1.0e6,
                    flyby_continuity_tol_kms=1.0e6,
                    independent_cross_check=False,
                    independent_tol_kms=1.0e6,
                    max_revs=max_revs,
                )
            except Exception:
                continue
            cand_for_traj = cand

        # Phase-1 / 3 closure produces per_encounter_vinf as a derived value;
        # the V_inf-match residual is THAT (our DE440 number) vs the cycler's
        # sourced seed V_inf — not the proposal's vinf_tuple bin (which is the
        # GEOMETRIC seed).
        terminal_vinf_kms = _terminal_vinf_kms(closure)
        vinf_match_residual_kms = abs(terminal_vinf_kms - seed_vinf)

        # Gate on the loose Phase-4 thresholds (the search is a discovery
        # probe; reporting wide is the point).
        if closure.closure_residual_kms > closure_tol_kms:
            continue
        if closure.flyby_continuity_max_dv_kms > flyby_continuity_tol_kms:
            continue
        post_closure_vinf_tol = (
            vinf_terminal_post_closure_tol_kms
            if vinf_terminal_post_closure_tol_kms is not None
            else vinf_terminal_tol_kms
        )
        if vinf_match_residual_kms > post_closure_vinf_tol:
            # Re-check post-closure (proposal bin may have shifted under
            # TOF optimisation).  The proposal-time tolerance
            # ``vinf_terminal_tol_kms`` filters the BFS output; the
            # post-closure tolerance gates the survivor record.  Setting
            # ``vinf_terminal_post_closure_tol_kms`` to a larger value
            # admits candidates whose proposal-V_inf was in-tolerance but
            # whose actual closure-V_inf drifted — useful for a discovery
            # probe that wants to see the residual distribution.
            continue
        if (
            independent_cross_check
            and closure.independent_check_residual_kms is not None
            and closure.independent_check_residual_kms > independent_tol_kms
        ):
            continue

        # Wrap the (possibly optimised) candidate as a ``precursor_mga``
        # trajectory for the PrecursorMatch record.
        trajectory_record = _candidate_to_trajectory(
            cand_for_traj,
            inserts_into=cycler_id,
            periapsis_altitudes_km=None,
        )
        if trajectory_record is None:
            continue

        # Literature check.
        sig = _signature_from_candidate(trajectory_record)
        if literature_check_search is None:
            lit = LiteratureCheckResult(
                status="inconclusive",
                citation=None,
                doi=None,
                confidence=0.0,
                query_trail=[],
                notes="No literature_check_search provided; verdict deferred",
            )
        else:
            lit = check_literature(sig, search=literature_check_search)

        alignment_score = _epoch_alignment_score(
            trajectory_record.launch_epoch_utc,
            target_phase_window_utc,
        )

        matches.append(
            PrecursorMatch(
                candidate=trajectory_record,
                cycler_id=cycler_id,
                cycler_seed_vinf_kms=seed_vinf,
                vinf_match_residual_kms=vinf_match_residual_kms,
                epoch_alignment_score=alignment_score,
                closure=closure,
                literature_check=lit,
            )
        )

    matches.sort(key=lambda m: m.quality_score())
    _emit(
        f"precursor matcher: complete — {len(matches)} survivors "
        f"(literature-fresh: {sum(1 for m in matches if m.is_literature_fresh())})"
    )
    return matches


# --------------------------------------------------------------------------- #
# JSONL serialisation
# --------------------------------------------------------------------------- #


def precursor_match_to_jsonl_record(match: PrecursorMatch) -> dict[str, Any]:
    """Render a :class:`PrecursorMatch` as a JSONL-emittable dict.

    Captures every field downstream consumers need to (a) reproduce the
    closure, (b) inspect the V_inf-match quality, and (c) cross-check the
    literature verdict.  Floats are not rounded — the consumer rounds for
    display.
    """
    traj = match.candidate
    closure = match.closure
    record = {
        "cycler_id": match.cycler_id,
        "cycler_seed_vinf_kms": match.cycler_seed_vinf_kms,
        "vinf_match_residual_kms": match.vinf_match_residual_kms,
        "epoch_alignment_score": match.epoch_alignment_score,
        "quality_score": match.quality_score(),
        "candidate": {
            "sequence": list(traj.sequence),
            "leg_tofs_days": list(traj.leg_tofs_days),
            "vinf_kms_at_encounters": list(traj.vinf_kms_at_encounters),
            "launch_epoch_utc": traj.launch_epoch_utc,
            "orbit_class": traj.orbit_class,
            "inserts_into": traj.inserts_into,
            "validity_window_start_utc": traj.validity_window_start_utc,
            "validity_window_end_utc": traj.validity_window_end_utc,
        },
        "closure": {
            "closure_residual_kms": closure.closure_residual_kms,
            "flyby_continuity_max_dv_kms": closure.flyby_continuity_max_dv_kms,
            "per_encounter_vinf_kms": list(closure.per_encounter_vinf_kms),
            "independent_check_residual_kms": closure.independent_check_residual_kms,
            "converged": closure.converged,
        },
        "literature_check": {
            "status": match.literature_check.status,
            "citation": match.literature_check.citation,
            "doi": match.literature_check.doi,
            "confidence": match.literature_check.confidence,
            "matched_url": match.literature_check.matched_url,
            "notes": match.literature_check.notes,
            "n_queries": len(match.literature_check.query_trail),
        },
        "is_literature_fresh": match.is_literature_fresh(),
    }
    return record


__all__ = [
    "PrecursorMatch",
    "find_cycler_precursors",
    "precursor_match_to_jsonl_record",
]
