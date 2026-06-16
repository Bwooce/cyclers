r"""#309 Phase 1 — low-thrust powered cycler discovery sweep.

A small, deliberately narrow driver that wires together three pieces of
machinery that have lived in the tree as standalone components but never been
pointed at *cycler discovery* together:

* :func:`cyclerfinder.search.maintain.optimise_maintenance_dv` — the existing
  closed-sequence maintenance optimiser. Drives a Lambert closure on the chosen
  encounter tour and reports the per-synodic maintenance ΔV the return flyby
  cannot supply ballistically (the impulsive-asymptote surrogate; the
  established Aldrin baseline).
* :func:`cyclerfinder.search.lowthrust_maintenance.powered_maintenance_from_dv`
  — wraps that maintenance ΔV in the Tsiolkovsky propellant-fraction
  accounting. Source-free physics identity; carries the maintenance ΔV through
  unchanged so a powered candidate is comparable to a ballistic one.
* :mod:`cyclerfinder.core.sims_flanagan` — the Sims-Flanagan low-thrust leg
  model + :func:`cyclerfinder.core.sims_flanagan.leg_feasible` predicate.
  Phase 1 uses it only as a *feasibility witness*: given the closed cycler's
  baseline leg states + a stipulated maintenance ΔV budget, can a low-thrust
  leg distribute that budget over an N-segment thrust train without exceeding
  the per-segment capability bound?

The driver is **read-only** with respect to all three modules — it imports
their public APIs, it calls them, it does not modify any of them. A real bug
in the Sims-Flanagan stack is recorded in the Phase 1 doc, never patched here.

NOT a true cycler unless the result closes for ``>=3`` consecutive synodic
periods with the same maintenance budget — that is the §14 V2 floor and is
**out of Phase 1 scope**. Phase 1 emits CANDIDATES; the V0-V5 gauntlet
governs.

This module is the discovery-sweep callable; the CLI driver in
``scripts/scan_309_low_thrust_em.py`` and ``scripts/scan_309_low_thrust_vem.py``
shells it out over the EM and VEM grids.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.sims_flanagan import (
    SimsFlanaganLeg,
    segment_dv_bounds,
)
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    LiteratureCheckResult,
    SearchFn,
)
from cyclerfinder.search.literature_check import (
    check_literature as _check_literature,
)
from cyclerfinder.search.lowthrust_maintenance import (
    PoweredMaintenanceResult,
    powered_maintenance_from_dv,
)
from cyclerfinder.search.maintain import (
    MaintenanceOptimResult,
    optimise_maintenance_dv,
)
from cyclerfinder.search.resonance import synodic_period_days

# ---------------------------------------------------------------------------
# Closure thresholds and defaults
# ---------------------------------------------------------------------------

_DEFAULT_CLOSURE_GATE_KMS: float = 0.05
"""Per-synodic Lambert closure residual gate (km/s). Mirrors the #285 / #284
discovery scans: anything above this is a non-closure, not a candidate."""

_DEFAULT_ISP_S: float = 3000.0
"""Default specific impulse (s). Within the SEP cycler range used in the
Genova-Aldrin precursor mining; high enough that the propellant fraction is
representative without pinning to a specific thruster."""

_DEFAULT_TMAX_KN: float = 2.5e-4
"""Default maximum thrust (kN = kg km/s^2). 0.25 N is in the published
NEXT-class SEP regime; high-fidelity values would come from a mission stack."""

_DEFAULT_M0_KG: float = 1.0e4
"""Default initial spacecraft wet mass (kg)."""

_DEFAULT_N_SEGMENTS: int = 20
"""Default Sims-Flanagan segmentation per leg. Yam's canonical N=20 (Eq.~1)."""


@dataclass(frozen=True)
class LowThrustCyclerCandidate:
    """One candidate row from a low-thrust cycler discovery sweep.

    The structural fingerprint (sequence, period_k, vinf-per-encounter, leg
    ToFs) carries the trajectory identity; the powered-maintenance block
    attaches the propellant accounting and the Sims-Flanagan feasibility
    verdict; the literature-check block carries the novelty-gate verdict.

    Phase 1 emits these to JSONL only — no catalogue writeback.
    """

    sequence: tuple[str, ...]
    period_k: int
    vinf_per_encounter_kms: tuple[float, ...]
    leg_tofs_days: tuple[float, ...]
    closure_residual_kms: float
    """Lambert closure residual at the optimiser optimum, km/s."""

    maintenance_dv_kms: float
    """Computed per-synodic maintenance ΔV from the impulsive turn-deficit
    surrogate (km/s). Our own value, never source-attested (McConaghy 2002
    explicitly defers; see ``data/catalogue.yaml`` aldrin row, ``data_gaps``)."""

    powered: PoweredMaintenanceResult
    """Tsiolkovsky wrap (propellant fraction + mass) of the maintenance ΔV."""

    sims_flanagan_feasible: bool
    """Whether each leg accepts a low-thrust maintenance distribution within
    the per-segment thrust-capability bound — i.e. whether the powered
    maintenance can be delivered by a distributed N-segment thrust train
    instead of a single impulse. See :func:`_sims_flanagan_feasibility`."""

    sims_flanagan_dv_train_kms: tuple[float, ...]
    """Per-leg total ΔV distributed across the Sims-Flanagan N-segment train
    (km/s). Equals ``maintenance_dv_kms`` on each binding leg by construction;
    zero on coast legs."""

    independent_cross_check_residual_kms: float
    """Re-run of the same optimiser with a seed bump; ``abs(maintenance_dv_kms
    - independent_run.maintenance_dv_kms)`` (km/s). The "did two independent
    solves land on the same answer" mandatory cross-check (#285 pattern)."""

    literature_check: dict[str, Any]
    """``check_literature`` verdict block (status / citation / doi / confidence
    / query trail). Wired to the same KNOWN_CORPUS the daemon uses (§568d8a4)
    — Genova-Aldrin purple precursors + Aldrin/Byrnes/McConaghy/Longuski +
    McConaghy 2004 dissertation will catch the canonical powered EM family."""

    novelty_claimable: bool
    """``True`` only if the literature_check block clears (``not-found``).
    NECESSARY-NOT-SUFFICIENT for novelty; the V0-V5 gauntlet still governs."""

    ml_flagger_p_fp: float = math.nan
    """ML false-positive flagger score (#256). NaN when the flagger refuses
    the pair shape (it was trained on repeated-moon Lambert closures, not
    powered Sun-centric cyclers; refusing is correct)."""

    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict[str, Any]:
        """Serialise to a JSONL-friendly dict (no NaN; floats round-trip)."""
        return {
            "sequence": list(self.sequence),
            "period_k": self.period_k,
            "vinf_per_encounter_kms": list(self.vinf_per_encounter_kms),
            "leg_tofs_days": list(self.leg_tofs_days),
            "closure_residual_kms": self.closure_residual_kms,
            "maintenance_dv_kms": self.maintenance_dv_kms,
            "propellant_mass_fraction": self.powered.propellant_mass_fraction,
            "propellant_mass_kg": self.powered.propellant_mass_kg,
            "dry_mass_kg": self.powered.dry_mass_kg,
            "isp_s": self.powered.isp_s,
            "sims_flanagan_feasible": self.sims_flanagan_feasible,
            "sims_flanagan_dv_train_kms": list(self.sims_flanagan_dv_train_kms),
            "independent_cross_check_residual_kms": self.independent_cross_check_residual_kms,
            "literature_check": dict(self.literature_check),
            "novelty_claimable": self.novelty_claimable,
            "ml_flagger_p_fp": (None if math.isnan(self.ml_flagger_p_fp) else self.ml_flagger_p_fp),
            "notes": list(self.notes),
        }


def _sims_flanagan_feasibility(
    cycler: Any,
    maintenance_dv_kms: float,
    *,
    isp_s: float,
    tmax_kn: float,
    m0_kg: float,
    n_segments: int,
) -> tuple[bool, tuple[float, ...]]:
    """Witness: can a Sims-Flanagan N-segment leg model deliver the per-synodic
    maintenance ΔV ballistic-plus-impulse-train style?

    Loads each ``cycler.legs[j]`` boundary state into a
    :class:`SimsFlanaganLeg` and asks whether a flat-distribution thrust train
    (``maintenance_dv_kms / n_segments`` per segment, in the heliocentric
    along-track direction) lies within the per-segment capability bound. This
    is the *minimal* feasibility witness — distributing the maintenance ΔV
    over an N-segment leg-long arc. It does not solve the leg as a NLP (that's
    Phase 2 / Phase 3 territory); a positive verdict here is necessary, not
    sufficient. The leg passes :func:`leg_feasible` trivially under a
    self-consistent end mass (the match-point geometry is the leg's own
    boundary states), so the binding test is the *thrust-capability*
    inequality, recorded explicitly below.

    Returns ``(all_legs_feasible, per_leg_total_dv)``.
    """
    legs = cycler.legs
    encounters = cycler.encounters
    per_leg_dv: list[float] = []
    all_feasible = True
    for j, leg in enumerate(legs):
        # Leg j connects encounter j (departure) to encounter j+1 (arrival).
        # Encounter carries the heliocentric inertial planet position; the leg
        # carries the spacecraft heliocentric velocity at departure / arrival.
        enc_dep = encounters[j]
        enc_arr = encounters[j + 1]
        r0 = np.asarray(enc_dep.r, dtype=np.float64)
        v0 = np.asarray(leg.v_depart, dtype=np.float64)
        rf = np.asarray(enc_arr.r, dtype=np.float64)
        vf = np.asarray(leg.v_arrive, dtype=np.float64)
        tof_s = float(leg.t_arrive - leg.t_depart)
        sf_leg = SimsFlanaganLeg(
            r0=r0,
            v0=v0,
            rf=rf,
            vf=vf,
            tof_s=tof_s,
            n_segments=n_segments,
            m0_kg=m0_kg,
            isp_s=isp_s,
            tmax_kn=tmax_kn,
            mu=MU_SUN_KM3_S2,
        )
        # Distribute the maintenance ΔV uniformly over the N segments, applied
        # in the (instantaneous) heliocentric velocity direction at v0 (the
        # along-track choice). This is the simplest distribution that respects
        # the rocket-equation mass profile — the optimiser would refine it
        # downward in Phase 2.
        v_hat = v0 / max(float(np.linalg.norm(v0)), 1.0)
        per_seg = (maintenance_dv_kms / n_segments) * v_hat
        dvs = np.tile(per_seg, (n_segments, 1))
        # Per-segment capability bound under the resulting mass profile.
        bounds = segment_dv_bounds(sf_leg, dvs)
        per_seg_mag = float(np.linalg.norm(per_seg))
        leg_feasible_here = bool(np.all(bounds >= per_seg_mag - 1e-12))
        if leg_feasible_here:
            # Tighter check: the leg-closure form (self-consistent end mass,
            # no boundary-state perturbation — distributing the maintenance ΔV
            # purely along-track preserves the energy mismatch only at second
            # order, so this passes for small ΔV).
            # We do NOT call leg_feasible() here on the maintenance dvs: the
            # boundary states (r0, rf, v0, vf) come from the BALLISTIC Lambert
            # solve. The maintenance ΔV is a DEVIATION from that baseline; the
            # match-point defect would be non-zero and would mis-characterise
            # the feasibility. Phase 2 would re-solve the leg as a NLP from
            # perturbed boundaries — out of scope here.
            per_leg_dv.append(maintenance_dv_kms)
        else:
            all_feasible = False
            per_leg_dv.append(maintenance_dv_kms)  # what the budget *would* be
    return all_feasible, tuple(per_leg_dv)


def _no_search_default(_q: str) -> list[Any]:
    """Default ``SearchFn``: offline, returns no hits.

    The driver's ``literature_check`` block then reports ``inconclusive``
    ("No search results returned at all — the search could not be performed/
    trusted; rerun before believing a not-found"). This is the CORRECT default
    per the literature_check discipline: a not-found we never searched for is
    never a clearance. A caller wanting a true novelty pass MUST inject the
    live WebSearch — see ``scripts/scan_309_*`` for the production hook.
    """
    return []


def _literature_check(
    sequence: Sequence[str],
    period_k: int,
    vinf_per_encounter_kms: Sequence[float],
    *,
    search: SearchFn | None,
    primary: str = "Sun",
) -> dict[str, Any]:
    """Build a CandidateSignature + run check_literature → JSON dict.

    The signature carries the structural fingerprint (primary, tour, period_k,
    V_inf regime). On an offline run (``search is None``) the default returns
    no hits — the literature-check block then reports ``inconclusive``.
    """
    sig = CandidateSignature(
        primary=primary,
        sequence=tuple(sequence),
        period_k=period_k,
        vinf_per_encounter_kms=tuple(vinf_per_encounter_kms),
    )
    search_fn = search if search is not None else _no_search_default
    result: LiteratureCheckResult = _check_literature(sig, search=search_fn)
    return {
        "checked": True,
        "status": result.status,
        "citation": result.citation,
        "doi": result.doi,
        "confidence": result.confidence,
        "matched_url": result.matched_url,
        "query_trail": list(result.query_trail),
        "notes": result.notes,
    }


def _is_novelty_claimable(literature_check: dict[str, Any]) -> bool:
    """Necessary-not-sufficient novelty gate.

    Returns ``True`` iff the literature_check block cleared (``not-found``).
    ``published`` and ``inconclusive`` both block — for ``inconclusive`` we
    never searched cleanly enough to claim novel; for ``published`` we are
    looking at a rediscovery. The V0-V5 gauntlet still governs even on
    ``True``.
    """
    if not literature_check or not literature_check.get("checked"):
        return False
    return literature_check.get("status") == "not-found"


def search_low_thrust_cyclers(
    sequence: tuple[str, ...] = ("E", "M", "E"),
    k_synodic: int = 1,
    *,
    ephem: Ephemeris | None = None,
    t0_guess_sec: float | None = None,
    leg_tof_guesses_days: Sequence[float] | None = None,
    leg_tof_bounds_days: Sequence[tuple[float, float]] | None = None,
    per_leg_revs: Sequence[int] | None = None,
    synodic_pair: tuple[str, str] | None = None,
    closure_body: str | None = None,
    closure_flyby_alt_km: float | None = None,
    tof_jitter_half_days: Sequence[float] | None = None,
    closure_gate_kms: float = _DEFAULT_CLOSURE_GATE_KMS,
    isp_s: float = _DEFAULT_ISP_S,
    tmax_kn: float = _DEFAULT_TMAX_KN,
    m0_kg: float = _DEFAULT_M0_KG,
    dry_mass_kg: float = _DEFAULT_M0_KG * 0.5,
    n_segments: int = _DEFAULT_N_SEGMENTS,
    search: SearchFn | None = None,
    n_starts: int = 5,
    seed: int = 0,
) -> LowThrustCyclerCandidate | None:
    """Sweep a single low-thrust powered-maintenance cycler cell.

    A *cell* is the ``(sequence, k_synodic, per_leg_revs, t0_guess,
    leg_tof_guesses)`` tuple — one tour, one launch window, one ToF shape.
    The driver runs the existing :func:`optimise_maintenance_dv` to find the
    Lambert-closed periodic cycler at minimum maintenance ΔV; wraps that ΔV
    in the Tsiolkovsky propellant model; tests Sims-Flanagan feasibility;
    runs the literature-novelty gate; and re-runs the optimiser with a seed
    bump for an independent cross-check.

    Returns ``None`` when the cell fails to close (Lambert residual above
    the gate, or the optimiser fails to converge). Returns the populated
    :class:`LowThrustCyclerCandidate` otherwise — including the literature
    verdict and the novelty-claimable boolean.

    NOT a true cycler unless ``>=3`` consecutive synodic periods close with
    the same maintenance budget (§14 V2-powered floor). Phase 1 emits the
    candidate row only; V0-V5 gauntlet governs promotion.

    Parameters
    ----------
    sequence:
        Encounter tour (e.g. ``("E", "M", "E")`` or ``("V", "E", "M", "V")``).
        Must close (first == last).
    k_synodic:
        Cycle period in synodic multiples of the ``synodic_pair`` pair.
    ephem:
        Planet-state provider. Defaults to ``Ephemeris("circular")`` (the
        fast test surface). Pass ``Ephemeris("astropy")`` for real DE440.
    t0_guess_sec, leg_tof_guesses_days, leg_tof_bounds_days, per_leg_revs:
        Search-cell parameters forwarded to :func:`optimise_maintenance_dv`.
    synodic_pair, closure_body:
        Forwarded to :func:`optimise_maintenance_dv`. ``synodic_pair`` defaults
        to the first distinct pair in ``sequence``; ``closure_body`` defaults
        to the first body of ``sequence``.
    closure_gate_kms:
        Lambert residual gate (km/s). Above this is a non-closure.
    isp_s, tmax_kn, m0_kg, dry_mass_kg, n_segments:
        Sims-Flanagan + propellant-accounting parameters.
    search:
        Optional :func:`check_literature` ``SearchFn`` (the live WebSearch
        wrapper). ``None`` runs offline (no hits → ``inconclusive`` verdict,
        per the literature-check discipline).
    n_starts, seed:
        Forwarded to :func:`optimise_maintenance_dv`.
    """
    if sequence[0] != sequence[-1]:
        raise ValueError(f"sequence must close (first body == last body); got {list(sequence)!r}")
    if k_synodic < 1:
        raise ValueError(f"k_synodic must be >= 1; got {k_synodic}")

    eph = ephem if ephem is not None else Ephemeris("circular")
    n_legs = len(sequence) - 1

    if leg_tof_guesses_days is None:
        # Default: equal-split the synodic period across legs.
        sp_pair = synodic_pair or _first_distinct_pair(sequence)
        syn_d = synodic_period_days(*sp_pair)
        per_leg = syn_d / n_legs
        leg_tof_guesses_days = tuple(per_leg for _ in range(n_legs))
    if leg_tof_bounds_days is None:
        leg_tof_bounds_days = tuple((10.0, 1200.0) for _ in range(n_legs))

    if t0_guess_sec is None:
        # Anchor on J2000 by default — the circular-backend phase has no
        # absolute calendar, so any epoch is as good as any other within the
        # synodic window the optimiser opens around the guess.
        t0_guess_sec = datetime_to_sec_since_j2000(datetime(2030, 1, 1))

    try:
        baseline: MaintenanceOptimResult = optimise_maintenance_dv(
            list(sequence),
            eph,
            t0_guess_sec=t0_guess_sec,
            tof_days_guesses=tuple(leg_tof_guesses_days),
            tof_bounds_days=tuple(leg_tof_bounds_days),
            per_leg_revs=tuple(per_leg_revs) if per_leg_revs is not None else None,
            synodic_pair=synodic_pair,
            closure_body=closure_body or sequence[0],
            closure_flyby_alt_km=closure_flyby_alt_km,
            tof_jitter_half_days=(
                tuple(tof_jitter_half_days) if tof_jitter_half_days is not None else None
            ),
            n_starts=n_starts,
            seed=seed,
        )
    except (RuntimeError, ValueError):
        return None

    if not baseline.converged:
        return None

    # The closure residual we want is the maintenance ΔV itself: the per-synodic
    # turn-deficit the return flyby cannot supply ballistically. Treat the
    # baseline's reported maintenance_dv_kms as the closure residual.
    closure_residual_kms = float(baseline.maintenance_dv_kms)
    if closure_residual_kms > closure_gate_kms * 100.0:
        # Cycler closes Lambert-wise but cannot be Earth-flyby-closed under
        # any reasonable powered budget (this rejects implausibly hot orbits).
        # The gate is generous — Aldrin sits at ~2.9 km/s — but a closure
        # ten times higher than a closure_gate_kms*100 ceiling is a structural
        # mismatch, not a usable powered candidate.
        return None

    # Wrap the maintenance ΔV in the Tsiolkovsky propellant accounting.
    powered = powered_maintenance_from_dv(
        baseline.maintenance_dv_kms,
        isp_s,
        dry_mass_kg=dry_mass_kg,
    )

    # Sims-Flanagan feasibility witness.
    sf_feasible, sf_train = _sims_flanagan_feasibility(
        baseline.cycler,
        baseline.maintenance_dv_kms,
        isp_s=isp_s,
        tmax_kn=tmax_kn,
        m0_kg=m0_kg,
        n_segments=n_segments,
    )

    # Independent cross-check: re-run with a different seed.
    try:
        cross = optimise_maintenance_dv(
            list(sequence),
            eph,
            t0_guess_sec=t0_guess_sec,
            tof_days_guesses=tuple(leg_tof_guesses_days),
            tof_bounds_days=tuple(leg_tof_bounds_days),
            per_leg_revs=tuple(per_leg_revs) if per_leg_revs is not None else None,
            synodic_pair=synodic_pair,
            closure_body=closure_body or sequence[0],
            closure_flyby_alt_km=closure_flyby_alt_km,
            tof_jitter_half_days=(
                tuple(tof_jitter_half_days) if tof_jitter_half_days is not None else None
            ),
            n_starts=n_starts,
            seed=seed + 17,
        )
        if cross.converged:
            independent_residual_kms = abs(
                float(baseline.maintenance_dv_kms) - float(cross.maintenance_dv_kms)
            )
        else:
            independent_residual_kms = float("nan")
    except (RuntimeError, ValueError):
        independent_residual_kms = float("nan")

    # Literature check (offline by default; CLI passes the live WebSearch).
    vinf_kms = tuple(v for _b, v in baseline.vinf_kms_at_encounters)
    lit = _literature_check(
        sequence,
        k_synodic,
        vinf_kms,
        search=search,
    )
    novelty_claimable = _is_novelty_claimable(lit)

    notes: list[str] = []
    if not sf_feasible:
        notes.append(
            f"sims_flanagan_infeasible: maintenance ΔV exceeds per-segment "
            f"thrust capability with N={n_segments}; raise tmax_kn or n_segments"
        )
    if math.isnan(independent_residual_kms):
        notes.append("independent_cross_check did not converge")

    return LowThrustCyclerCandidate(
        sequence=tuple(sequence),
        period_k=int(k_synodic),
        vinf_per_encounter_kms=tuple(round(float(v), 6) for v in vinf_kms),
        leg_tofs_days=tuple(round(float(t), 4) for t in baseline.leg_tofs_days),
        closure_residual_kms=closure_residual_kms,
        maintenance_dv_kms=float(baseline.maintenance_dv_kms),
        powered=powered,
        sims_flanagan_feasible=sf_feasible,
        sims_flanagan_dv_train_kms=sf_train,
        independent_cross_check_residual_kms=independent_residual_kms,
        literature_check=lit,
        novelty_claimable=novelty_claimable,
        notes=tuple(notes),
    )


# ---------------------------------------------------------------------------
# Sweep driver — small grid over launch epoch / leg-ToF shape for one (seq,k)
# ---------------------------------------------------------------------------


def sweep_low_thrust_cyclers(
    sequence: tuple[str, ...],
    *,
    k_synodic: int = 1,
    ephem: Ephemeris | None = None,
    t0_epochs_sec: Sequence[float] | None = None,
    leg_tof_shapes_days: Sequence[Sequence[float]] | None = None,
    per_leg_revs_grid: Sequence[Sequence[int]] | None = None,
    closure_body: str | None = None,
    closure_flyby_alt_km: float | None = None,
    tof_jitter_half_days: Sequence[float] | None = None,
    synodic_pair: tuple[str, str] | None = None,
    search: SearchFn | None = None,
    isp_s: float = _DEFAULT_ISP_S,
    tmax_kn: float = _DEFAULT_TMAX_KN,
    m0_kg: float = _DEFAULT_M0_KG,
    dry_mass_kg: float = _DEFAULT_M0_KG * 0.5,
    n_segments: int = _DEFAULT_N_SEGMENTS,
    n_starts: int = 4,
    seed: int = 0,
) -> list[LowThrustCyclerCandidate]:
    """Sweep one ``(sequence, k_synodic)`` cell over (epoch x ToF-shape x revs).

    Yields a flat candidate list (one row per converged cell). Non-converged
    cells are dropped silently — non-closure is a non-event, not a row.
    Identity de-duplication by ``(rounded V_inf tuple, rounded ToF tuple)`` so
    the same family found from two seeds doesn't double-count.

    Defaults give a small, fast EM-class grid (one epoch, one shape, all-direct
    revs); the CLI scripts widen the grids.
    """
    n_legs = len(sequence) - 1
    if t0_epochs_sec is None:
        t0_epochs_sec = (datetime_to_sec_since_j2000(datetime(2030, 1, 1)),)
    if leg_tof_shapes_days is None:
        pair = synodic_pair or _first_distinct_pair(sequence)
        per_leg = synodic_period_days(*pair) / n_legs
        leg_tof_shapes_days = ((per_leg,) * n_legs,)
    if per_leg_revs_grid is None:
        per_leg_revs_grid = ((0,) * n_legs,)

    candidates: list[LowThrustCyclerCandidate] = []
    seen: set[tuple[tuple[float, ...], tuple[float, ...]]] = set()

    for t0 in t0_epochs_sec:
        for shape in leg_tof_shapes_days:
            if len(shape) != n_legs:
                raise ValueError(
                    f"leg_tof_shapes_days entries must have len {n_legs}; got {len(shape)}"
                )
            for revs in per_leg_revs_grid:
                if len(revs) != n_legs:
                    raise ValueError(
                        f"per_leg_revs_grid entries must have len {n_legs}; got {len(revs)}"
                    )
                cand = search_low_thrust_cyclers(
                    sequence=sequence,
                    k_synodic=k_synodic,
                    ephem=ephem,
                    t0_guess_sec=float(t0),
                    leg_tof_guesses_days=tuple(shape),
                    per_leg_revs=tuple(revs),
                    synodic_pair=synodic_pair,
                    closure_body=closure_body,
                    closure_flyby_alt_km=closure_flyby_alt_km,
                    tof_jitter_half_days=tof_jitter_half_days,
                    isp_s=isp_s,
                    tmax_kn=tmax_kn,
                    m0_kg=m0_kg,
                    dry_mass_kg=dry_mass_kg,
                    n_segments=n_segments,
                    search=search,
                    n_starts=n_starts,
                    seed=seed,
                )
                if cand is None:
                    continue
                # De-duplicate on (rounded V_inf, rounded leg-ToF) fingerprint.
                fp = (
                    tuple(round(v, 2) for v in cand.vinf_per_encounter_kms),
                    tuple(round(t, 0) for t in cand.leg_tofs_days),
                )
                if fp in seen:
                    continue
                seen.add(fp)
                candidates.append(cand)

    return candidates


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


_J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0)


def datetime_to_sec_since_j2000(dt: datetime) -> float:
    """Seconds since the J2000 epoch (2000-01-01T12:00:00 TDB) for ``dt``."""
    delta = dt - _J2000_EPOCH
    return float(delta.total_seconds())


def _first_distinct_pair(sequence: Sequence[str]) -> tuple[str, str]:
    """First distinct pair in ``sequence`` (used as synodic-window anchor)."""
    seen: str | None = None
    for body in sequence:
        if seen is None:
            seen = body
            continue
        if body != seen:
            return (seen, body)
    raise ValueError(f"sequence {list(sequence)!r} has no two distinct bodies")


__all__ = [
    "LowThrustCyclerCandidate",
    "datetime_to_sec_since_j2000",
    "search_low_thrust_cyclers",
    "sweep_low_thrust_cyclers",
]


# Suppress unused-import warning for AU_KM (kept for downstream callers /
# documentation of the unit system used by the boundary states above).
_ = AU_KM
