"""Phase 4 — the novelty discovery loop (the Forge).

``discover_novel`` is the construction-first novelty frontier. Where the M7
:func:`cyclerfinder.data.discover.discover` walks the idealised optimiser over a
cell stream and matches signatures, this module drives the **real-ephemeris
ballistic corrector** (:func:`cyclerfinder.search.correct.ballistic_correct` via
the parallel scan engine :func:`cyclerfinder.search.scan.scan_parallel`) over
multi-arc topologies, bridges each closed chain into a full
:class:`~cyclerfinder.model.cycler.Cycler`, computes its canonical signature,
matches it against the catalogue (supersession-aware, R1 delta 3), runs the
Axis-A code-path agreement cross-check, and folds everything through
:func:`cyclerfinder.verify.gauntlet.run_gauntlet`.

Why multi-arc E-M, not VEM single-ellipse
-----------------------------------------
The #110 dense VEM scan produced **zero** bend-feasible closures (floors
17.9 / 18.5 km/s vs the sourced Jones 2.42-7.0): VEM ballistic novelty at
coplanar-circular / single-ellipse-per-leg fidelity is empirically nil (see
``data/OUTSTANDING.md`` M-ED open-research entry). Where the corrector
*demonstrably* closes bend-feasible, sub-cap chains is the **E-M multi-arc
space** (the S1L1-prototype family, ``scripts/correct_s1l1_twoarc.py``, and the
Sanchez-regime chain, ``tests/search/test_correct_sanchez_regime.py``). That is
the novelty frontier this loop drives by default.

Routing (the human gate)
------------------------
Each survivor's verdict is the gauntlet tier:

* **GOLD** is *impossible* for a novel candidate — GOLD requires an independent
  source, and a novel candidate has none by definition. (A *rediscovery* — a
  candidate whose signature matches a sourced catalogue row — can in principle
  reach GOLD via Axis C, but only when the caller supplies that row's
  provenance; the discovery loop itself never fabricates a source.)
* **SILVER** — machine-confirmed (>= 2 code paths agree) but unsourced. This is
  the novelty holding tier: it routes to the human-review queue and is **never**
  auto-promoted (golden discipline).
* **BRONZE** — a weak, non-refuted signal (Axis A had fewer than two paths).
* **REJECTED** — falsified / adversarially refuted / a disqualifying axis
  failure.

Golden discipline
-----------------
This module asserts no computed value as a golden EXPECTED. Signatures, agreement
predicates, and verdicts are all *classifications* derived from the candidate's
own geometry; a novel candidate is held at SILVER pending a human, by design.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.catalog import (
    CanonicalSignature,
    Catalog,
    MatchResult,
    canonical_signature,
    load_catalog,
    match,
)
from cyclerfinder.model.cycler import Cycler
from cyclerfinder.search.construct import construct_cycler
from cyclerfinder.search.correct import BallisticClosureResult
from cyclerfinder.search.moon_prune import prune_topology_legs
from cyclerfinder.search.scan import (
    ScanResult,
    build_epoch_branch_grid,
    scan_parallel,
)
from cyclerfinder.verify.agreement import AgreementReport, crosscheck_code_paths
from cyclerfinder.verify.gauntlet import (
    ValidationVerdict,
    run_gauntlet,
)

DAY_S = 86400.0

MatchOutcome = Literal["known", "probable-match-NEEDS-HUMAN", "novel"]


# ---------------------------------------------------------------------------
# The multi-arc bridge: BallisticClosureResult -> Cycler
# ---------------------------------------------------------------------------


def cycler_from_closure(
    closure: BallisticClosureResult,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    ephem: Ephemeris,
    *,
    mu_sun: float = MU_SUN_KM3_S2,
) -> Cycler:
    """Reconstruct a full :class:`Cycler` from a converged ballistic closure.

    The corrector returns only ``(t0_sec, tof_days, per-encounter V_inf)``; the
    Axis-A code-path machinery and the canonical signature need a full
    :class:`Cycler` (every leg's departure/arrival velocity + every encounter's
    V_inf in/out vectors). This deterministically re-solves the same Lambert legs
    via :func:`cyclerfinder.search.construct.construct_cycler` at the converged
    encounter epochs and the same ``(n_revs, branch)`` selection — so the
    reconstructed geometry is bit-for-bit the one the corrector closed.

    Parameters
    ----------
    closure:
        A converged :class:`~cyclerfinder.search.correct.BallisticClosureResult`.
    sequence, per_leg_revs, per_leg_branch:
        The topology the corrector ran (same objects passed to
        :func:`~cyclerfinder.search.correct.ballistic_correct`).
    ephem:
        The DE440 (``"astropy"``) ephemeris.
    mu_sun:
        Heliocentric gravitational parameter.

    Returns
    -------
    Cycler
        ``sense="n/a"`` (a discovered candidate carries no catalogue sense until
        a human assigns one).
    """
    epochs = [closure.t0_sec]
    for tof in closure.tof_days:
        epochs.append(epochs[-1] + float(tof) * DAY_S)
    return construct_cycler(
        list(sequence),
        epochs,
        ephem,
        mu_sun=mu_sun,
        max_revs_per_leg=list(per_leg_revs),
        branch_per_leg=list(per_leg_branch),
    )


# ---------------------------------------------------------------------------
# Verdict classification (pure logic; supersession-aware)
# ---------------------------------------------------------------------------


def _superseded_chain(entry_raw: dict[str, Any] | None) -> tuple[str, ...]:
    """Extract a catalogue row's ``superseded_by`` chain (schema v4.3).

    Accepts the row dict's ``superseded_by`` as a string or list; returns a
    normalised tuple. ``None`` / missing -> empty tuple.
    """
    if not entry_raw:
        return ()
    raw = entry_raw.get("superseded_by")
    if raw is None:
        return ()
    if isinstance(raw, str):
        return (raw,)
    if isinstance(raw, (list, tuple)):
        return tuple(str(x) for x in raw)
    return ()


def classify_candidate_verdict(
    *,
    candidate_id: str,
    agreed: bool,
    n_paths_available: int,
    n_paths_passed: int,
    match_outcome: MatchOutcome,
    known_id: str | None,
    superseded_by: tuple[str, ...],
    falsified: bool,
    notes: str | None = None,
) -> ValidationVerdict:
    """Route a candidate to a tiered verdict via the pure gauntlet combiner.

    A novel candidate has no independent source, so Axis C is left empty
    (``provenance_tier=None``); the only way it reaches GOLD is if the caller
    supplied a sourced match — which the discovery loop never fabricates. The
    result is therefore SILVER for a machine-confirmed novel candidate, BRONZE
    for a thin one, REJECTED for a falsified one (golden discipline: no
    auto-promotion past SILVER).

    This wraps :func:`cyclerfinder.verify.gauntlet.run_gauntlet` by synthesising
    a minimal :class:`~cyclerfinder.verify.agreement.AgreementReport`-shaped
    summary from the already-computed agreement booleans; it does not recompute
    physics.
    """
    from cyclerfinder.verify.agreement import (
        ConstructionOptimiserPathResult,
        KeplerRepropPathResult,
        LamberthubPathResult,
    )

    empty_a = LamberthubPathResult(available=False, per_leg=(), max_diff_mps=0.0, passed=False)
    empty_b = ConstructionOptimiserPathResult(
        available=False,
        resonant_vinf_kms={},
        cycler_vinf_kms={},
        construction_max_diff_kms=float("inf"),
        optimiser_available=False,
        optimiser_vinf_kms={},
        optimiser_max_diff_kms=None,
        max_diff_kms=float("inf"),
        passed=False,
    )
    empty_c = KeplerRepropPathResult(
        available=False, per_leg_residual_km=(), max_residual_km=float("inf"), passed=False
    )
    agreement = AgreementReport(
        lamberthub=empty_a,
        construction_optimiser=empty_b,
        kepler_reprop=empty_c,
        n_paths_available=n_paths_available,
        n_paths_passed=n_paths_passed,
        agreed=agreed,
    )
    return run_gauntlet(
        candidate_id,
        agreement=agreement,
        persistence_reports=(),
        provenance_tier=None,
        corroboration=None,
        falsified=falsified,
        known_id=known_id if match_outcome != "novel" else None,
        superseded_by=superseded_by,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Per-candidate evaluation: closure -> signature -> agreement -> verdict
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NoveltyFinding:
    """One evaluated discovery candidate (frozen, audit-complete).

    Attributes
    ----------
    candidate_id:
        Deterministic identifier (topology + epoch seed).
    signature:
        The candidate's :class:`CanonicalSignature`, or ``None`` if the closure
        geometry could not be reduced to one.
    match_outcome:
        The catalogue match outcome (``known`` / ``probable`` / ``novel``).
    known_id:
        The matched catalogue row id (if any).
    superseded_by:
        The matched row's supersession chain (R1 delta 3), empty if none.
    verdict:
        The gauntlet :class:`ValidationVerdict`.
    agreement:
        The Axis-A :class:`AgreementReport`.
    vinf_per_encounter_kms:
        Per-encounter V_inf magnitudes (km/s) of the closed chain.
    tof_days:
        Per-leg times of flight (days).
    bend_feasible:
        Whether every intermediate flyby is bend-feasible.
    max_vinf_kms:
        Peak per-encounter V_inf (km/s).
    """

    candidate_id: str
    signature: CanonicalSignature | None
    match_outcome: MatchOutcome
    known_id: str | None
    superseded_by: tuple[str, ...]
    verdict: ValidationVerdict
    agreement: AgreementReport
    vinf_per_encounter_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    bend_feasible: bool
    max_vinf_kms: float
    # Topology + corrector seed, carried so an adversarial panel can re-run the
    # exact closure without lossy id-parsing (Phase 5 orchestrator).
    sequence: tuple[str, ...] = ()
    per_leg_revs: tuple[int, ...] = ()
    per_leg_branch: tuple[str, ...] = ()
    period_k: int = 0
    t0_sec: float = 0.0
    slack_leg: int = 0


def _is_multi_arc(cycler: Cycler, *, a_bin_au: float = 0.01, e_bin: float = 0.01) -> bool:
    """Whether the cycler's legs trace more than one distinct heliocentric ellipse.

    A single-ellipse cycler (Aldrin family) has every leg on one ``(a, e)``; a
    multi-arc cycler (S1L1 / Sanchez-regime E-M-E-E) has legs on two or more
    distinct ellipses. Binned to the catalogue signature widths so numerical
    noise does not split one ellipse into two.
    """
    from cyclerfinder.model.cycler import orbit_elements_au

    seen: set[tuple[float, float]] = set()
    for leg, enc in zip(cycler.legs, cycler.encounters[:-1], strict=True):
        a_au, e = orbit_elements_au(enc.r, leg.v_depart, MU_SUN_KM3_S2)
        if not (0.0 <= e < 1.0) or a_au <= 0.0:
            # A hyperbolic / degenerate leg is its own (non-elliptic) arc.
            return True
        seen.add((round(a_au / a_bin_au), round(e / e_bin)))
    return len(seen) >= 2


def _agreement_predicate(agreement: AgreementReport, *, multi_arc: bool) -> tuple[bool, int, int]:
    """Re-derive (agreed, n_available, n_passed) for the candidate's regime.

    For a **multi-arc real-ephemeris** candidate the resonance-construction path
    (b) is structurally inapplicable: it builds ONE circular-coplanar resonant
    ellipse from the first leg's ``(a, e)`` and asserts that single ellipse
    reproduces every encounter's V_inf. A multi-arc chain is, by definition, not
    one ellipse — and comparing a coplanar-idealised construction against a real
    DE440 multi-arc geometry is exactly the cross-fidelity confusion the Forge
    refuses (the S1L1 5.65-vs-4.99 episode). So for multi-arc candidates path (b)
    is treated as **unavailable** (not a failing veto), and the ">= 2 paths
    agree" predicate rests on the two genuinely-independent real-ephemeris
    witnesses: in-house Lambert vs lamberthub (a) and forward Kepler
    re-propagation (c).

    For single-ellipse candidates the unmodified report is returned verbatim.
    """
    if not multi_arc:
        return (agreement.agreed, agreement.n_paths_available, agreement.n_paths_passed)
    # Drop path (b) entirely: count only paths (a) and (c).
    paths = (agreement.lamberthub, agreement.kepler_reprop)
    n_available = sum(1 for p in paths if p.available)
    n_passed = sum(1 for p in paths if p.available and p.passed)
    all_available_passed = all(p.passed for p in paths if p.available)
    agreed = n_available >= 2 and all_available_passed
    return (agreed, n_available, n_passed)


def _candidate_id(
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_k: int,
    t0_sec: float,
) -> str:
    """A deterministic, sortable id for a discovery candidate."""
    seq = "-".join(sequence)
    revs = "".join(str(r) for r in per_leg_revs)
    br = "".join(b[0] for b in per_leg_branch)
    return f"novel|{seq}|k{period_k}|r{revs}|b{br}|t{round(t0_sec)}"


def evaluate_closure(
    closure: BallisticClosureResult,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_k: int,
    ephem: Ephemeris,
    catalog: Catalog | None = None,
    model_assumption: str = "analytic-ephemeris",
    falsified: bool = False,
    slack_leg: int = 0,
    mu_central: float = MU_SUN_KM3_S2,
) -> NoveltyFinding:
    """Bridge a closed chain to a Cycler, sign it, match it, gauntlet it.

    The construction-first per-candidate pipeline:

    1. :func:`cycler_from_closure` -> full :class:`Cycler`.
    2. :func:`cyclerfinder.data.catalog.canonical_signature` -> signature.
    3. :func:`cyclerfinder.data.catalog.match` -> match outcome
       (supersession surfaced from the matched row's raw YAML).
    4. :func:`cyclerfinder.verify.agreement.crosscheck_code_paths` -> Axis A.
    5. :func:`classify_candidate_verdict` -> the tiered verdict.

    ``falsified`` is the Axis-D hook: a caller's adversarial probe may set it
    True to force REJECTED (the falsification gate). It is ALSO set
    automatically when the closure is physically inadmissible — not
    bend-feasible, or busting the V_inf cap — because such a chain is not a
    realisable ballistic cycler (a flyby that cannot deliver the required turn
    is an impossible bend). A bend-infeasible closure therefore routes to
    REJECTED, never SILVER (golden discipline: only realisable candidates reach
    the human queue).
    """
    catalog = catalog if catalog is not None else load_catalog()
    cycler = cycler_from_closure(
        closure, sequence, per_leg_revs, per_leg_branch, ephem, mu_sun=mu_central
    )
    signature = canonical_signature(
        cycler, model_assumption=model_assumption, period_k=period_k, mu_central=mu_central
    )
    match_result: MatchResult = match(signature, catalog)
    known_id = match_result.entry.id if match_result.entry is not None else None
    superseded_by = (
        _superseded_chain(match_result.entry.raw) if match_result.entry is not None else ()
    )

    # Physical-admissibility falsification: a non-converged, bend-infeasible, or
    # over-cap closure is not a realisable cycler -> Axis D refutes it.
    physically_inadmissible = not (
        closure.converged and closure.bend_feasible and closure.vinf_cap_ok
    )
    falsified = falsified or physically_inadmissible

    agreement = crosscheck_code_paths(cycler, ephem, mu=mu_central)
    multi_arc = _is_multi_arc(cycler)
    agreed, n_available, n_passed = _agreement_predicate(agreement, multi_arc=multi_arc)
    cid = _candidate_id(sequence, per_leg_revs, per_leg_branch, period_k, closure.t0_sec)
    verdict = classify_candidate_verdict(
        candidate_id=cid,
        agreed=agreed,
        n_paths_available=n_available,
        n_paths_passed=n_passed,
        match_outcome=match_result.outcome,
        known_id=known_id,
        superseded_by=superseded_by,
        falsified=falsified,
        notes=f"discover_novel/{model_assumption}",
    )
    vinf = closure.vinf_per_encounter_kms
    return NoveltyFinding(
        candidate_id=cid,
        signature=signature,
        match_outcome=match_result.outcome,
        known_id=known_id,
        superseded_by=superseded_by,
        verdict=verdict,
        agreement=agreement,
        vinf_per_encounter_kms=vinf,
        tof_days=closure.tof_days,
        bend_feasible=closure.bend_feasible,
        max_vinf_kms=float(max(vinf)) if vinf else float("inf"),
        sequence=sequence,
        per_leg_revs=per_leg_revs,
        per_leg_branch=per_leg_branch,
        period_k=period_k,
        t0_sec=closure.t0_sec,
        slack_leg=slack_leg,
    )


# ---------------------------------------------------------------------------
# The discovery loop over the E-M multi-arc space
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TopologySpec:
    """One multi-arc topology to scan (sequence + per-leg rev/branch + period).

    Attributes
    ----------
    sequence:
        Encounter body sequence, e.g. ``("E","M","E","E")``.
    per_leg_revs, per_leg_branch:
        Per-leg Lambert rev count + branch.
    period_k:
        Period in synodic multiples (carried into the signature).
    period_sec:
        Target heliocentric period (seconds) the corrector pins.
    tof_seed_days:
        The *free* (slack-eliminated) ToF seed for the corrector.
    slack_leg:
        Index of the pinned (period-absorbing) leg.
    """

    sequence: tuple[str, ...]
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    period_k: int
    period_sec: float
    tof_seed_days: tuple[float, ...]
    slack_leg: int


# The S1L1-prototype / Sanchez-regime E-M-E-E multi-arc topology family. This is
# the empirically-demonstrated bend-feasible frontier (NOT VEM single-ellipse,
# which the #110 dense scan proved is ballistically nil). Period and seed trace
# to the sourced Russell 4.991gG arcs (1.4612 + 2.8096 yr); the corrector's V_inf
# output is OUR computation (never a golden).
_S1L1_PERIOD_SEC = (1.4612 + 2.8096) * DAYS_PER_JULIAN_YEAR * DAY_S


def em_multiarc_topologies() -> tuple[TopologySpec, ...]:
    """The default E-M multi-arc topology set for the novelty loop.

    The E-M-E-E two-arc family the prototype + Sanchez-regime gate close, swept
    over the two converging E-E Lambert branches (low/high at n_revs in {1,2}).
    The slack (period-absorbing) leg is the E-E phasing leg (index 2).
    """
    specs: list[TopologySpec] = []
    for ee_revs, ee_branch in ((1, "low"), (2, "low"), (1, "high"), (2, "high")):
        specs.append(
            TopologySpec(
                sequence=("E", "M", "E", "E"),
                per_leg_revs=(0, 0, ee_revs),
                per_leg_branch=("single", "single", ee_branch),
                period_k=2,
                period_sec=_S1L1_PERIOD_SEC,
                tof_seed_days=(154.0, 379.0),  # E-M, M-E free legs; E-E is slack
                slack_leg=2,
            )
        )
    return tuple(specs)


# ---------------------------------------------------------------------------
# Moon-system (Jovian) topology set + centred sweep (Forge Phase 6)
# ---------------------------------------------------------------------------

# Galilean synodic ToF seeds (days) — the Laplace-resonance spacings the moon
# orbital periods set: Io 1.769 d, Europa 3.551 d, Ganymede 7.155 d. These are
# SEEDS ONLY (NON-GOLDEN); the corrector refines them. period_sec seeds from the
# resonance multiple (the sum of the three moon periods).
_GALILEAN_PERIODS_DAYS: dict[str, float] = {
    "Io": 1.769,
    "Europa": 3.551,
    "Ganymede": 7.155,
    "Callisto": 16.689,
}


def jovian_galilean_topologies() -> tuple[TopologySpec, ...]:
    """The Galilean Jovian topology set the first campaign sweeps (plan Task 2.0).

    The empirically-closing Io-Europa-Ganymede(-Io) chain (#76 Phase 3 closes it;
    the open question is bend feasibility under VILM gating). At the Laplace-
    resonance ToF seeds these short Galilean legs admit only the direct (0-rev,
    ``"single"``) Lambert branch — a fictitious multi-rev/low-high branch would be
    rejected at reconstruction (``construct_cycler`` raises on a missing branch),
    so the physically-honest enrichment is over the *period multiple* (``period_k``
    -> the resonance-multiple period seed) and the launch-epoch grid (the
    ``n_epochs`` sweep in :func:`discover_novel_moon`), NOT over branches the
    geometry does not support. The slack (period-absorbing) leg is the Ganymede->Io
    return leg (index 2). ``period_sec`` seeds from the Laplace-resonance multiple.
    NON-GOLDEN: the corrector refines the period/ToF seeds.
    """
    seq = ("Io", "Europa", "Ganymede", "Io")
    one_synodic_days = (
        _GALILEAN_PERIODS_DAYS["Io"]
        + _GALILEAN_PERIODS_DAYS["Europa"]
        + _GALILEAN_PERIODS_DAYS["Ganymede"]
    )
    # Free (slack-eliminated) ToF seeds: the Io->Europa and Europa->Ganymede legs;
    # the Ganymede->Io leg is slack (period-absorbing).
    tof_seed_days = (
        _GALILEAN_PERIODS_DAYS["Europa"],
        _GALILEAN_PERIODS_DAYS["Ganymede"],
    )
    specs: list[TopologySpec] = []
    for period_k in (1, 2):
        specs.append(
            TopologySpec(
                sequence=seq,
                per_leg_revs=(0, 0, 0),
                per_leg_branch=("single", "single", "single"),
                period_k=period_k,
                period_sec=period_k * one_synodic_days * DAY_S,
                tof_seed_days=tof_seed_days,
                slack_leg=2,
            )
        )
    return tuple(specs)


def _moon_period_days(moon: str) -> float:
    """Orbital period (days) of a moon about its primary, Kepler III.

    DERIVED at call time from ``SATELLITES[moon].sma_km`` + the primary's mu
    (``PRIMARIES``), never a hand-copied magic number — same Kepler-III mechanism
    :func:`cyclerfinder.core.satellites.mean_motion_deg_day_about` uses. This is
    the sourced sibling of the literal ``_GALILEAN_PERIODS_DAYS`` table (the latter
    pre-dates the registry and is retained only for the original Galilean campaign).
    NON-GOLDEN here: a ToF/period *seed* the corrector refines, not a published anchor.
    """
    sat = SATELLITES[moon]
    period_s = 2.0 * np.pi * np.sqrt(sat.sma_km**3 / PRIMARIES[sat.primary])
    return float(period_s / DAY_S)


def _tour_topologies(
    seq: tuple[str, ...],
    *,
    period_ks: tuple[int, ...],
) -> tuple[TopologySpec, ...]:
    """Build period-multiple-swept VILM topologies for a closed moon tour.

    Shared construction discipline (mirrors :func:`jovian_galilean_topologies`):
    short resonant moon legs admit only the direct 0-rev ``"single"`` Lambert
    branch (a fictitious multi-rev/low-high branch is rejected at reconstruction),
    so enrichment is over the *period multiple* (``period_k``) and the launch-epoch
    grid, NOT branches the geometry cannot support. The return leg (last index) is
    the slack (period-absorbing) leg. ``period_sec`` seeds from the sum-of-periods
    synodic multiple. The free (slack-eliminated) ToF seeds are each intermediate
    leg's destination-moon period. NON-GOLDEN: the corrector refines all seeds.
    """
    if len(seq) < 3 or seq[0] != seq[-1]:
        raise ValueError(f"tour sequence must be a closed chain of >=3 bodies: {seq!r}")
    n_legs = len(seq) - 1
    one_synodic_days = sum(_moon_period_days(m) for m in seq[:-1])
    # Free ToF seeds: every leg except the slack (return) leg, seeded by the
    # destination moon's period.
    tof_seed_days = tuple(_moon_period_days(seq[i + 1]) for i in range(n_legs - 1))
    specs: list[TopologySpec] = []
    for period_k in period_ks:
        specs.append(
            TopologySpec(
                sequence=seq,
                per_leg_revs=tuple(0 for _ in range(n_legs)),
                per_leg_branch=tuple("single" for _ in range(n_legs)),
                period_k=period_k,
                period_sec=period_k * one_synodic_days * DAY_S,
                tof_seed_days=tof_seed_days,
                slack_leg=n_legs - 1,
            )
        )
    return tuple(specs)


def saturnian_titan_tour_topologies() -> tuple[TopologySpec, ...]:
    """Saturnian midsize+Titan resonant moon-tour topology set (Forge Phase 6 #178).

    Closed resonant sub-chains over the Saturnian midsize moons + Titan. Periods
    are DERIVED from :data:`cyclerfinder.core.satellites.SATELLITES` (sma + Saturn
    mu, Kepler III) via :func:`_moon_period_days` — no hardcoded magic numbers. The
    construction discipline matches :func:`jovian_galilean_topologies`: 0-rev
    ``"single"`` branches only, return leg = slack, ``period_k`` in (1, 2). Tours:

    * Dione-Rhea-Titan-Dione: the outer resonant chain anchored on Titan.
    * Enceladus-Dione-Rhea-Enceladus: an inner midsize resonant chain.
    * Enceladus-Tethys-Dione-Rhea-Titan style is too long/non-resonant for the
      no-leveraging genome, so the swept set is the two 3-leg sub-chains.

    NON-GOLDEN: the corrector refines the period/ToF seeds.
    """
    tours = (
        ("Dione", "Rhea", "Titan", "Dione"),
        ("Enceladus", "Dione", "Rhea", "Enceladus"),
    )
    specs: list[TopologySpec] = []
    for seq in tours:
        specs.extend(_tour_topologies(seq, period_ks=(1, 2)))
    return tuple(specs)


def jovian_galilean_permutation_topologies() -> tuple[TopologySpec, ...]:
    """Galilean orderings beyond the swept I-E-G-I chain (Forge Phase 6 #178).

    Alternate Galilean orderings (including Callisto) the first campaign did NOT
    sweep: Io-Ganymede-Europa-Io, Europa-Ganymede-Callisto-Europa, and
    Io-Europa-Callisto-Io. Periods are DERIVED from :data:`SATELLITES` (sma +
    Jupiter mu, Kepler III) via :func:`_moon_period_days`. Same discipline as
    :func:`jovian_galilean_topologies` (0-rev ``"single"`` branches only, return
    leg = slack), swept over ``period_k`` in (1, 2, 3).

    NON-GOLDEN: the corrector refines the period/ToF seeds.
    """
    tours = (
        ("Io", "Ganymede", "Europa", "Io"),
        ("Europa", "Ganymede", "Callisto", "Europa"),
        ("Io", "Europa", "Callisto", "Io"),
    )
    specs: list[TopologySpec] = []
    for seq in tours:
        specs.extend(_tour_topologies(seq, period_ks=(1, 2, 3)))
    return tuple(specs)


def discover_novel_moon(
    *,
    base_t0_sec: float,
    topologies: Sequence[TopologySpec] | None = None,
    center: str = "Jupiter",
    budget_kms: float = 50.0,
    n_epochs: int = 16,
    span_days: float = 8.0,
    vinf_cap: float = 14.0,
    vinf_seed_kms: float = 4.0,
    max_workers: int | None = None,
    catalog: Catalog | None = None,
    distinct_only: bool = True,
) -> Iterator[NoveltyFinding]:
    """The VILM-pruned centred Jovian novelty sweep (plan Task 2.1).

    A sibling of :func:`discover_novel` for the moon space. For each topology it
    applies the Phase-1 incremental prune (:func:`prune_topology_legs`) BEFORE
    building its scan grid — pruned topologies are SKIPPED (their per-leg reasons
    are available to the caller via :func:`prune_topology_legs`, which the
    orchestrator records for the empty-region report). Survivors run the centred
    scan (``Ephemeris(model="circular", center=center)``,
    ``mu_central=PRIMARIES[center]``, threaded via the ``center`` ScanPoint field)
    and each closure flows through :func:`evaluate_closure` verbatim — the bridge
    -> signature -> match -> agreement -> gauntlet pipeline is centre-agnostic.

    The closures read ``novel`` against the null-numeric Jovian catalogue bucket
    (design §3a); a bend-INFEASIBLE closure is auto-falsified to REJECTED inside
    :func:`evaluate_closure` (never SILVER).
    """
    ephem = Ephemeris(model="circular", center=center)
    catalog = catalog if catalog is not None else load_catalog()
    topologies = tuple(topologies) if topologies is not None else jovian_galilean_topologies()

    seen_hashes: set[str] = set()
    for spec in topologies:
        survives, _reasons = prune_topology_legs(
            spec,
            vinf_seed_kms=vinf_seed_kms,
            budget_kms=budget_kms,
            primary=center,
        )
        if not survives:
            continue
        t0_seeds = _epoch_seeds_sec(base_t0_sec=base_t0_sec, n_epochs=n_epochs, span_days=span_days)
        grid = build_epoch_branch_grid(
            sequence=spec.sequence,
            period_sec=spec.period_sec,
            vinf_cap=vinf_cap,
            t0_seeds_sec=t0_seeds,
            branch_topologies=[(spec.per_leg_revs, spec.per_leg_branch)],
            tof_seed_days=spec.tof_seed_days,
            slack_leg=spec.slack_leg,
            center=center,
        )
        results: list[ScanResult] = scan_parallel(
            grid,
            max_workers=max_workers,
            closed_only=True,
        )
        for sr in results:
            finding = evaluate_closure(
                sr.result,
                sequence=spec.sequence,
                per_leg_revs=spec.per_leg_revs,
                per_leg_branch=spec.per_leg_branch,
                period_k=spec.period_k,
                ephem=ephem,
                catalog=catalog,
                model_assumption="circular-coplanar",
                slack_leg=spec.slack_leg,
                mu_central=PRIMARIES[center],
            )
            if distinct_only and finding.signature is not None:
                if finding.signature.hash in seen_hashes:
                    continue
                seen_hashes.add(finding.signature.hash)
            yield finding


def _epoch_seeds_sec(
    *,
    base_t0_sec: float,
    n_epochs: int,
    span_days: float,
) -> list[float]:
    """``n_epochs`` launch-epoch seeds spread evenly across ``span_days``.

    Centred on ``base_t0_sec`` so the prototype's converging window is sampled.
    """
    if n_epochs <= 1:
        return [base_t0_sec]
    offsets = np.linspace(-span_days / 2.0, span_days / 2.0, n_epochs)
    return [base_t0_sec + float(o) * DAY_S for o in offsets]


def discover_novel(
    *,
    ephem: Ephemeris | None = None,
    topologies: Sequence[TopologySpec] | None = None,
    base_t0_sec: float,
    n_epochs: int = 16,
    span_days: float = 1280.0,
    vinf_cap: float = 14.0,
    max_workers: int | None = None,
    catalog: Catalog | None = None,
    distinct_only: bool = True,
) -> Iterator[NoveltyFinding]:
    """Drive the construction-first novelty loop over the E-M multi-arc space.

    For each topology, build an epoch x branch scan grid, run it in parallel
    through :func:`cyclerfinder.search.scan.scan_parallel`, and for every
    *closed* result evaluate the full pipeline (bridge -> signature -> match ->
    agreement -> gauntlet), yielding a :class:`NoveltyFinding` per survivor.

    Loop-until-dry semantics (bounded here): the deepening frontier is the
    topology set x epoch grid; this run is bounded by ``n_epochs`` x
    ``len(topologies)``. The unbounded version is the Phase 5 workflow's job.

    Parameters
    ----------
    ephem:
        DE440 (``"astropy"``) ephemeris; constructed if ``None``.
    topologies:
        Topology set; defaults to :func:`em_multiarc_topologies`.
    base_t0_sec:
        Centre launch epoch (seconds from J2000) for the scan window.
    n_epochs:
        Launch-epoch grid points per topology.
    span_days:
        Total launch-epoch span (centred on ``base_t0_sec``).
    vinf_cap:
        V_inf cap (km/s) handed to the corrector.
    max_workers:
        Worker process count for the parallel scan (``None`` -> cpu_count).
    catalog:
        Catalogue to match against; loaded once if ``None``.
    distinct_only:
        Collapse closures that bin to the same canonical signature hash to a
        single yielded finding (the distinct-family count).

    Yields
    ------
    NoveltyFinding
        One per closed candidate (deduped by signature hash when
        ``distinct_only``), in deterministic scan order.
    """
    ephem = ephem if ephem is not None else Ephemeris("astropy")
    catalog = catalog if catalog is not None else load_catalog()
    topologies = tuple(topologies) if topologies is not None else em_multiarc_topologies()

    seen_hashes: set[str] = set()
    for spec in topologies:
        t0_seeds = _epoch_seeds_sec(base_t0_sec=base_t0_sec, n_epochs=n_epochs, span_days=span_days)
        grid = build_epoch_branch_grid(
            sequence=spec.sequence,
            period_sec=spec.period_sec,
            vinf_cap=vinf_cap,
            t0_seeds_sec=t0_seeds,
            branch_topologies=[(spec.per_leg_revs, spec.per_leg_branch)],
            tof_seed_days=spec.tof_seed_days,
            slack_leg=spec.slack_leg,
        )
        results: list[ScanResult] = scan_parallel(
            grid,
            ephem_model="astropy",
            max_workers=max_workers,
            closed_only=True,
        )
        for sr in results:
            finding = evaluate_closure(
                sr.result,
                sequence=spec.sequence,
                per_leg_revs=spec.per_leg_revs,
                per_leg_branch=spec.per_leg_branch,
                period_k=spec.period_k,
                ephem=ephem,
                catalog=catalog,
                slack_leg=spec.slack_leg,
            )
            if distinct_only and finding.signature is not None:
                if finding.signature.hash in seen_hashes:
                    continue
                seen_hashes.add(finding.signature.hash)
            yield finding


__all__ = [
    "MatchOutcome",
    "NoveltyFinding",
    "TopologySpec",
    "classify_candidate_verdict",
    "cycler_from_closure",
    "discover_novel",
    "discover_novel_moon",
    "em_multiarc_topologies",
    "evaluate_closure",
    "jovian_galilean_permutation_topologies",
    "jovian_galilean_topologies",
    "saturnian_titan_tour_topologies",
]
