"""Alternating-double-cycler construction operator (#526).

Generalizes the "switched-double-cycler" construction of Liang, Yang, Li,
Bai & Qin, "Callisto-Ganymede-Europa Triple Cyclers," *Journal of Guidance,
Control, and Dynamics*, Vol. 48, No. 1, 2025 (Technical Note, published
online 12 Sep 2024), DOI 10.2514/1.G008387 -- a portion presented as Paper
IAC-23,C1,9,9,x76777 at the 74th International Astronautical Congress,
Baku, Azerbaijan, 2-6 October 2023 -- into a reusable genome operator.

Grounding (#526 explicit requirement -- see the mis-citation this project
was burned by before, ``feedback_ground_citations_against_content``)
--------------------------------------------------------------------------
The citation was independently re-verified 2026-07-03 against the AIAA/ARC
article page (title, full author list, DOI and venue match exactly:
Guoliang Liang, Hongwei Yang, Shuang Li, Xiaoli Bai and Limin Qin; JGCD
Vol. 48 No. 1; published online 12 Sep 2024), Semantic Scholar (independent
author-list cross-check) and Unpaywall (``doi:10.2514/1.g008387`` ->
``is_oa: false``, ``has_repository_copy: false`` -- confirms this is a
closed-access AIAA note with no legitimate open copy anywhere, so the
free AIAA "Read Now" first-page preview and this project's own prior
corpus mining are the only sourced content, not a substitute full read).
This rules out the same-author/same-year citation-collision failure mode
(cf. the Hernandez-Jones-Jesick 2017 VEM-vs-Galilean-triple-cycler
collision) -- there is exactly one paper at this DOI, and it says what the
task assumed.

The paper itself was already filed and deeply digested under an earlier
task (#216, mined 2026-06-11, character-by-character rescanned 2026-06-12):

* ``papers/liang-2024-callisto-ganymede-europa-triple-cyclers-JGCD.pdf``
  (private corpus; text-layer; see ``docs/notes/CORPUS_INDEX.md``).
* ``docs/notes/2026-06-11-liang-2024-cge-triple-cyclers-mining.md`` -- the
  near-resonance derivation (Eqs. 1-2), the switched-double-cycler
  strategy (Sec. III.B) and the published Tables 1-7.
* ``docs/notes/2026-06-13-liang-abc-reproduction.md`` -- same-model
  numeric reproduction results.
* :mod:`cyclerfinder.search.cge_scaffold` -- exact reproduction of the
  three idealized members (A/B/C) from their printed inputs.
* :mod:`cyclerfinder.search.moon_cycler_genome` -- the repeated-moon
  multi-revolution SEARCH genome (``MoonSystem``, ``MoonCyclerGenome``,
  ``EncounterGene``, ``LegGene``) the CGE members are encoded into, plus
  the reproduce-before-search gate against Liang's published members.

None of that machinery is duplicated here. What #526 asks for -- and what
was still missing -- is a REUSABLE operator that generalizes the paper's
construction *idea* beyond the one hard-coded Callisto-Ganymede-Europa
case, so it composes over other moon (or planet) triples.

The idea generalized (mining note Sec. 1 and Sec. 3.2)
--------------------------------------------------------
A triple-body cycler that must revisit a "hub" body H alternating with two
"partner" bodies P1 and P2 is hard to close directly when H, P1, P2 are
only *near*-resonant: the three-body relative phase configuration never
resets exactly (Callisto-Ganymede-Europa sit on a near-4:7, not exact,
synodic ratio -- Laplace resonance only binds Io-Europa-Ganymede). Liang
et al.'s fix (their "key trick", mining note Sec. 3.2) is to never close
the *triple* relationship at all: decompose the tour into two *double*
cyclers -- (H, P1) and (H, P2) -- each of which only needs an easy
two-body phase reset, and ALTERNATE between them once per near-resonance
quasi-period, with H visited at every switch. This module operationalizes
that decomposition as two independent, reusable pieces:

1. :func:`analyze_near_resonance` -- the generalized Eq. 1-2 near-
   commensurability analysis. Liang et al. compute it on the *chain* of
   consecutive synodic periods S(Callisto,Ganymede) and S(Ganymede,Europa)
   (Eq. 1: their ratio is ~7/4; Eq. 2: the per-quasi-period mismatch is
   4*S_C,G - 7*S_G,E = 0.7365 d). This function performs that same
   chain analysis for ANY three-body chain (A, B, C) and any small-integer
   denominator bound, given only the bodies' sidereal mean motions -- it
   is not hard-coded to Jupiter's Galilean moons.
2. :func:`build_alternating_double_cycler_seed` -- the switching operator
   itself. Given two ALREADY-CONSTRUCTED double-cycler encounter
   half-sequences that both start and end at the shared hub (e.g.
   ``(Callisto, Ganymede, Callisto)`` and ``(Callisto, Europa,
   Callisto)``), stitch them end-to-end (dropping the duplicated junction
   encounter) into one repeating
   :class:`~cyclerfinder.search.moon_cycler_genome.MoonCyclerGenome`
   sequence. Feeding it Liang's own two CGE halves reproduces their
   published CGCEC sequence exactly (the positive control below).

What this module deliberately does NOT do
-------------------------------------------
* It does NOT re-implement Liang et al.'s specific Eq. 16 epoch-multiplier
  law (``t_c1 = (6n-5)*T + t_c0``, ...) or their Sec. III.C sequential
  local Lambert / differential-evolution optimizer. Those are specific
  consequences of their particular choice (spacecraft period = T_Callisto,
  the 7:4 resonance, the CGCEC structure) and are already faithfully
  reproduced, unmodified, in :mod:`cyclerfinder.search.cge_scaffold`.
  Re-deriving a "generalized Eq. 16" from their specific multipliers
  without the full III.C derivation in hand would risk silently
  fabricating a formula the paper never published -- this module
  generalizes only the parts that generalize cleanly from the stated
  PRINCIPLE (near-resonance analysis + alternation), not their algebra.
* It does NOT construct a physical double cycler from scratch -- it
  composes two caller-supplied encounter half-sequences. Building those
  halves (e.g. via a Lambert/DE search, or via
  :mod:`cyclerfinder.search.moon_cycler_genome`'s corrector) is a separate,
  future concern.
* No novelty claim. No catalogue writeback. This is scaffolding that
  operationalizes a published construction *strategy* for reuse (e.g. the
  Saturnian Titan-Enceladus/Titan-Dione follow-on the mining note itself
  flags, Sec. 6 item 3, citing Russell & Strange 2009's Titan-Enceladus
  double cyclers) -- it does not claim any new cycler exists.

Positive controls (sourced-only, never circular)
---------------------------------------------------
* :func:`analyze_near_resonance` is checked against Liang et al.'s own
  printed Eq. 1-2 numbers using their exact Table 1 mean motions
  (:data:`cyclerfinder.search.cge_scaffold.LIANG_MEAN_MOTIONS_RAD_DAY`):
  ratio ~ 7/4 and mismatch = 0.7365 d, both printed in the paper.
* As an independent (non-golden) cross-check, the same analysis is run on
  this project's own registry mean motions
  (:mod:`cyclerfinder.core.satellites`, sourced from JPL SSD, NOT copied
  from the paper) and recovers the same 7:4 near-resonance with a mismatch
  within ~2% of Liang's printed value -- two independent data sources
  agreeing, without feeding the paper's numbers into the registry path.
* :func:`build_alternating_double_cycler_seed` is checked by re-deriving
  the published CGCEC sequence
  (:data:`cyclerfinder.search.moon_cycler_genome.CGE_SEQUENCE`) from its
  two double-cycler halves and confirming an exact match.
* Reusability (not a novelty claim) is demonstrated by running
  :func:`analyze_near_resonance` on a Saturnian chain
  (Enceladus-Dione-Rhea) via the registry -- the function executes
  body-agnostically and returns a well-formed (finite, positive-period)
  result; no resonance is asserted or claimed for that chain.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from cyclerfinder.search.moon_cycler_genome import (
    EncounterGene,
    LegGene,
    MoonCyclerGenome,
    MoonSystem,
)


def synodic_period_days(mean_motion_a_rad_day: float, mean_motion_b_rad_day: float) -> float:
    """Two-body synodic period (days) from sidereal mean motions (rad/day).

    ``S = 2*pi / |n_a - n_b|`` (Liang et al. 2024 Sec. II.A, the unlabeled
    equation immediately preceding Eq. 1).
    """
    dn = abs(mean_motion_a_rad_day - mean_motion_b_rad_day)
    if dn == 0.0:
        raise ValueError("equal mean motions: synodic period is undefined (co-orbital)")
    return 2.0 * math.pi / dn


def mean_motion_rad_day(system: MoonSystem, body: str) -> float:
    """Sidereal mean motion (rad/day) of ``body`` from the moon-system registry.

    Thin wrapper over :meth:`MoonSystem.period_days` (``n = 2*pi / T``) so
    :func:`analyze_near_resonance` can be driven directly from the in-repo
    registry (JPL SSD-sourced, :mod:`cyclerfinder.core.satellites`) instead
    of a hand-copied table -- the independent (non-golden) cross-check path
    described in the module docstring.
    """
    return 2.0 * math.pi / system.period_days(body)


@dataclass(frozen=True)
class NearResonanceAnalysis:
    """Generalized Eq. 1-2 near-commensurability analysis for a 3-body chain.

    For a chain ``(body_a, body_b, body_c)`` of consecutive bodies (e.g.
    Callisto-Ganymede-Europa, ordered by distance from the primary),
    ``synodic_ab_days`` / ``synodic_bc_days`` are the pairwise synodic
    periods of (a, b) and (b, c). ``ratio = synodic_ab / synodic_bc``;
    ``(p, q)`` is the lowest-order continued-fraction convergent
    ``ratio ~ p/q`` whose scale-free mismatch clears a tolerance (see
    :func:`analyze_near_resonance` / :func:`_lowest_order_convergent`), so
    that ``q * synodic_ab_days ~ p * synodic_bc_days``. ``mismatch_days`` is
    that residual (generalized Eq. 2) -- the closer to zero, the more
    exactly the a-b-c relative-phase configuration repeats every
    ``quasi_period_days = q * synodic_ab_days``.

    A small ``mismatch_days`` relative to ``quasi_period_days`` is a
    NECESSARY condition for the alternating-double-cycler construction to
    be worth attempting on this chain (mining note Sec. 1); it is not
    sufficient (nothing here checks that ballistic double cyclers actually
    exist for either pair).
    """

    body_a: str
    body_b: str
    body_c: str
    synodic_ab_days: float
    synodic_bc_days: float
    ratio: float
    p: int
    q: int
    mismatch_days: float
    quasi_period_days: float

    @property
    def mismatch_fraction(self) -> float:
        """``mismatch_days / quasi_period_days`` -- a scale-free tightness measure."""
        if self.quasi_period_days == 0.0:
            return float("inf")
        return abs(self.mismatch_days) / self.quasi_period_days


def _lowest_order_convergent(
    ratio: float,
    s_ab: float,
    s_bc: float,
    *,
    max_denominator: int,
    mismatch_tolerance: float,
) -> tuple[int, int]:
    """Smallest-denominator continued-fraction convergent of ``ratio`` whose
    ``mismatch_fraction`` clears ``mismatch_tolerance`` (falls back to the
    single best convergent within ``max_denominator`` if none does).

    A plain "closest fraction with denominator <= max_denominator" search
    (``Fraction.limit_denominator(max_denominator)`` alone) is NOT what
    Liang et al. report: for their own numbers, the globally-closest
    fraction with denominator <= 20 is 16/9 (mismatch fraction 0.09%), yet
    the paper picks the lower-order 7/4 (mismatch fraction ~1.5%) -- the
    standard resonance-identification convention of preferring the
    *lowest-order* commensurability that is already tight, rather than
    chasing an arbitrarily tighter but physically-arbitrary high-order
    match (Eq. 3, p. 3, ties T_S to small-integer multiples of ALL THREE
    orbital periods, not just the two-period ratio). Scanning convergents
    by increasing denominator and stopping at the first one under
    ``mismatch_tolerance`` reproduces the paper's own 7/4 exactly (see the
    ``liang_cge_near_resonance`` golden test) without hand-coding "4" as a
    special case.
    """
    last: tuple[int, int] | None = None
    for denominator in range(1, max_denominator + 1):
        frac = Fraction(ratio).limit_denominator(denominator)
        candidate = (frac.numerator, frac.denominator)
        if candidate == last:
            continue
        last = candidate
        p, q = candidate
        quasi_period = q * s_ab
        if quasi_period <= 0.0:
            continue
        mismatch_fraction = abs(q * s_ab - p * s_bc) / quasi_period
        if mismatch_fraction < mismatch_tolerance:
            return p, q
    assert last is not None  # max_denominator >= 1 guarantees at least one iteration
    return last


def analyze_near_resonance(
    body_a: str,
    body_b: str,
    body_c: str,
    mean_motions_rad_day: dict[str, float],
    *,
    max_denominator: int = 20,
    mismatch_tolerance: float = 0.05,
) -> NearResonanceAnalysis:
    """Generalized Eq. 1-2: near-commensurate synodic-period chain analysis.

    ``mean_motions_rad_day`` must carry entries for ``body_a``, ``body_b``
    and ``body_c`` -- caller-supplied so this works equally on the paper's
    own printed Table 1 values (the sourced-golden path,
    :func:`liang_cge_near_resonance`) or on an independent source such as
    the in-repo registry (:func:`mean_motion_rad_day`).

    ``(p, q)`` is the lowest-order continued-fraction convergent of
    ``ratio = synodic_ab / synodic_bc`` whose scale-free mismatch
    (:attr:`NearResonanceAnalysis.mismatch_fraction`) is below
    ``mismatch_tolerance`` (see :func:`_lowest_order_convergent` for why
    this -- not a naive closest-fraction search -- is what reproduces
    Liang et al.'s own 7/4 choice); if no convergent up to
    ``max_denominator`` clears the tolerance, the best convergent found is
    returned (and ``mismatch_fraction`` will reveal that it is loose).

    Raises :class:`ValueError` if ``max_denominator < 1`` (a degenerate
    search that could never find a nontrivial ``(p, q)``).
    """
    if max_denominator < 1:
        raise ValueError(f"max_denominator must be >= 1, got {max_denominator}")
    s_ab = synodic_period_days(mean_motions_rad_day[body_a], mean_motions_rad_day[body_b])
    s_bc = synodic_period_days(mean_motions_rad_day[body_b], mean_motions_rad_day[body_c])
    ratio = s_ab / s_bc
    p, q = _lowest_order_convergent(
        ratio, s_ab, s_bc, max_denominator=max_denominator, mismatch_tolerance=mismatch_tolerance
    )
    mismatch = q * s_ab - p * s_bc
    quasi_period = q * s_ab
    return NearResonanceAnalysis(
        body_a=body_a,
        body_b=body_b,
        body_c=body_c,
        synodic_ab_days=s_ab,
        synodic_bc_days=s_bc,
        ratio=ratio,
        p=p,
        q=q,
        mismatch_days=mismatch,
        quasi_period_days=quasi_period,
    )


def build_alternating_double_cycler_seed(
    system: MoonSystem,
    half_cycle_1: tuple[str, ...],
    half_cycle_2: tuple[str, ...],
    *,
    n_cycles: int = 1,
    legs: tuple[LegGene, ...] | None = None,
) -> MoonCyclerGenome:
    """Stitch two hub-sharing double-cycler halves into an alternating seed.

    This is the "key trick" of Liang et al. 2024 Sec. III.B operationalized
    as a composable operator: ``half_cycle_1`` and ``half_cycle_2`` must
    each be a moon (or planet) sequence that starts AND ends at the same
    shared hub body (e.g. ``("Callisto", "Ganymede", "Callisto")`` and
    ``("Callisto", "Europa", "Callisto")``). One repeating cycle is their
    concatenation with the duplicated junction hub encounter dropped --
    ``half_cycle_1 + half_cycle_2[1:]`` -- and ``n_cycles`` repeats that
    pattern, again dropping the duplicated hub encounter at each cycle
    boundary.

    Feeding Liang's own two CGE halves reproduces
    :data:`cyclerfinder.search.moon_cycler_genome.CGE_SEQUENCE`
    (``Callisto-Ganymede-Callisto-Europa-Callisto``) exactly -- the
    positive control for this function (see
    ``tests/genome/test_alternating_double_cycler.py``).

    ``legs`` defaults to ``LegGene(p=1, q=1, n_rev=1)`` for every leg (the
    Liang CGE convention: every leg is a one-revolution Lambert arc in the
    1:1 hub-period frame, per
    :func:`cyclerfinder.search.moon_cycler_genome.liang_member_genome`) --
    override it when the caller's double-cycler halves use a different
    resonance/rev-count bookkeeping.

    The result is a structural SEED only (all ``b_plane_angle_rad`` are
    0.0 and ``epoch_days`` is 0.0) -- it encodes the alternating TOPOLOGY,
    not a closed, dynamically-consistent trajectory; closing it is a
    separate corrector concern (e.g.
    :func:`cyclerfinder.search.moon_cycler_genome.liang_periodicity_residual`
    for the Liang CGE case).
    """
    if n_cycles < 1:
        raise ValueError(f"n_cycles must be >= 1, got {n_cycles}")
    for name, half in (("half_cycle_1", half_cycle_1), ("half_cycle_2", half_cycle_2)):
        if len(half) < 2:
            raise ValueError(f"{name} must have at least 2 bodies (hub...hub), got {half}")
        if half[0] != half[-1]:
            raise ValueError(f"{name} must start and end at the shared hub: got {half}")
    if half_cycle_1[0] != half_cycle_2[0]:
        raise ValueError(
            "half_cycle_1 and half_cycle_2 must share the same hub: "
            f"{half_cycle_1[0]!r} != {half_cycle_2[0]!r}"
        )

    one_cycle = half_cycle_1 + half_cycle_2[1:]
    sequence: list[str] = list(one_cycle)
    for _ in range(n_cycles - 1):
        sequence.extend(one_cycle[1:])

    encounters = tuple(EncounterGene(moon=m, b_plane_angle_rad=0.0) for m in sequence)
    n_legs = len(sequence) - 1
    resolved_legs = (
        legs if legs is not None else tuple(LegGene(p=1, q=1, n_rev=1) for _ in range(n_legs))
    )
    if len(resolved_legs) != n_legs:
        raise ValueError(
            f"expected {n_legs} legs for {len(sequence)} encounters, got {len(resolved_legs)}"
        )

    return MoonCyclerGenome(
        system=system,
        encounters=encounters,
        legs=resolved_legs,
        epoch_days=0.0,
        perijove_scale=1.0,
    )


# ---------------------------------------------------------------------------
# Liang CGE positive-control anchors (sourced-only; never our own output)
# ---------------------------------------------------------------------------


def liang_cge_near_resonance() -> NearResonanceAnalysis:
    """The generalized analysis run on Liang et al.'s own printed Table 1.

    Sourced golden: expect ``p == 7``, ``q == 4`` and
    ``mismatch_days`` == the paper's printed 0.7365 d (Eq. 2), both from
    :data:`cyclerfinder.search.cge_scaffold.LIANG_MEAN_MOTIONS_RAD_DAY`
    (their Table 1, p. 4) -- not a value this repo computed independently.
    """
    from cyclerfinder.search.cge_scaffold import LIANG_MEAN_MOTIONS_RAD_DAY

    return analyze_near_resonance(
        "Callisto", "Ganymede", "Europa", LIANG_MEAN_MOTIONS_RAD_DAY, max_denominator=20
    )


def liang_cge_alternating_seed() -> MoonCyclerGenome:
    """Reconstruct Liang's published CGCEC sequence from its two halves.

    Structural positive control for :func:`build_alternating_double_cycler_seed`:
    the stitched sequence must equal
    :data:`cyclerfinder.search.moon_cycler_genome.CGE_SEQUENCE` exactly.
    """
    from cyclerfinder.search.moon_cycler_genome import jupiter_system

    return build_alternating_double_cycler_seed(
        jupiter_system(),
        half_cycle_1=("Callisto", "Ganymede", "Callisto"),
        half_cycle_2=("Callisto", "Europa", "Callisto"),
    )
