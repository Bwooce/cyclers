"""#583 (stage 3b of #581): widened-domain Sun-Earth ER3BP bounded-drift search.

Widens #581 stage 2's niching-GA reproduction beyond Gurfil-Kasdin (2002)'s
own 12 published optimization-set boxes, reusing BOTH
``search/niching_ga.py::run_deterministic_crowding`` and
``core/er3bp_geocentric.py::gurfil_kasdin_fitness`` UNMODIFIED (per #583's
own scope) -- the only new work is the widened bounds/genome + the long-span
bounded-vs-divergent drift classifier that gates any survivor before it is
trusted (``data/validation/er3bp_drift_classifier.py``).

Explicit, written-down decisions (per #583's own mandate -- do not leave to
the reader to infer):

1. **Free vs fixed state dims.** Genome is 7-dimensional: the 6 interleaved
   geocentric state slots ``[x, x', y, y', z, z']`` (same convention as
   stage 2's ``TABLE2``/``table_interleaved_to_state``) PLUS ``theta0`` as a
   free 7th gene. Each of the 7 :data:`PARTITIONS` below fixes a SUBSET of
   the 6 state slots at 0 (mirroring which published families it targets)
   and frees the rest -- same "fix most components, vary a few" structure
   Gurfil-Kasdin's own Table 2 uses, generalized across their per-set
   groupings rather than redone per set.
2. **theta0 IS free** (recommended by the task spec, done here): every
   partition searches ``theta0 in [0, 2*pi)`` as a genome variable, instead
   of stage 2's per-set FIXED theta0 in {0, pi}. Leaving it fixed would not
   genuinely explore beyond the 12 published sets (theta0 is the true-
   anomaly launch phase; stage 2 only ever tested Earth-perihelion and
   Earth-aphelion launches).
3. **n_rev stays 1**, matching stage 2 and ``gurfil_kasdin_fitness``'s own
   default -- not widened here (no motivating reason given in the spec to
   widen it, and doing so would conflate two independent axes of widening
   in one run).
4. **escape_radius stays 0.5 AU** (the ``gurfil_kasdin_fitness`` default,
   ~50x Earth's Hill radius, ~1.5e6 km) -- already implicitly admits the
   full heliocentric co-orbital regime (quasi-satellites, horseshoes,
   Sun-Earth DROs) that the #583 corpus-anchor prerequisite exists to
   police. This is a stated decision, not a surprise: a bounded-drift
   survivor near that boundary is exactly the territory the 3 new
   ``literature_check.py`` anchors (Gurfil-Kasdin 2002, Sun-Earth
   co-orbital, Henon family-f) were filed to catch.
5. **Widened bounds** (the actual numbers, justified below) cover the space
   BETWEEN and AROUND the union of the 12 published sets' boxes -- NOT a
   redo of stage 2's own per-set boxes, and NOT the paper's own
   sqrt(2)/sqrt(3) norm-budgeting convention (that convention existed only
   to keep a MULTI-COMPONENT vector's overall magnitude inside
   [``LO_R``, ``HI_R``] when 2-3 components shared one budget; the widened
   search instead gives every free position/velocity component its OWN
   independent bound, which is itself part of the widening -- a published
   family can no longer "hide" behind the norm-sharing constraint).

   * ``LO_WIDE = LO_R`` (1e6 km, the paper's own floor). CORRECTED
     2026-07-13 after a Fable review of a failed positive control: an
     earlier choice of ``0.002`` AU (~299,200 km, closer-in than ``LO_R``)
     re-admitted a physically trivial, strictly fitness-dominant deep-Hill-
     sphere quasi-circular basin that Eq. 15 cannot distinguish from a
     genuine family (it rewards only annulus thinness, with no periodicity
     content) -- deterministic crowding collapsed the entire population onto
     it. See the constant's own inline comment and docs/notes/2026-07-13-
     583-corpus-anchors-and-drift-classifier.md for the full diagnosis.
   * ``HI_WIDE = 0.15`` AU (~2.244e7 km): ~2.24x FARTHER than the paper's
     own ``HI_R`` (1e7 km) -- stays at 30% of ``escape_radius`` (0.5 AU),
     leaving room for a genuinely bounded family to breathe without
     immediately tripping the escape death-penalty.
   * ``V_WIDE = 0.5`` (AU/rad, ~15 km/s): numerically the SAME magnitude as
     the paper's own ``V``, but every partition below grants it to EACH free
     velocity component INDEPENDENTLY and over the FULL SIGNED range
     ``[-V_WIDE, V_WIDE]`` -- stage 2's sets each searched only one half
     (``[-V, 0]`` or ``[0, V]``) per set and shared the budget via
     sqrt(2)/sqrt(3) when >1 component was free. This is the actual
     "between and around" widening for the velocity axis.

Positive control (the FIRST dispatch's deliverable, superseded below by a
partition redesign -- kept here for history): run partition ``P1`` (the
richest and cheapest -- 2 free state dims + theta0 -- covering 6 of the 14
published families A-F) to completion and confirm it still recovers ITS
neighboring known families under stage 2's own pre-registered non-bit-exact
match criterion (``IC_PROXIMITY_TOL``/``FEATURE_FACTOR_TOL``, imported
unmodified from ``run_581_gurfil_reproduction.py``). See
``docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md`` for the
result (1/6, root-caused to the two problems the redesign below fixes).

## 2026-07-13 partition redesign (this dispatch, replaces P1-P7)

The original 7 wide ``P1``-``P7`` partitions pooled up to 6 published
families into ONE box whose position range ran all the way to a single
global ceiling (``HI_WIDE = 0.15`` AU) regardless of where those families
actually live. Two independent Fable reviews (see the note) diagnosed why
``P1`` (A-F) only recovered 1/6: (a) every one of A-F's own published ICs
sits well inside the paper's own ``HI_R`` (1e7 km) -- the extra 2.24x of
position range past ``HI_R`` is territory NONE of the 6 targets occupy, so
it only dilutes the population's search effort -- and (b) pooling 6 (``P1``)
or 3 (old ``P5``: J, K, L, in the RICHEST 6-free-dim signature) published
families into one ``population=200`` niching-GA run gives each family far
less effective niche capacity than stage 2's own per-set runs (each covering
only 1-3 families).

The fix partitions RADIALLY into 3 bands and keeps each family-targeted box
narrow -- BUT the smoke-scale validation below (3 rounds, not 1) found the
niche-capacity/box-width story to be necessary but NOT sufficient: 2 more
hard bugs were caught and fixed along the way, and a THIRD, deeper limit
(inherent to Eq. 15 + deterministic crowding, not a bug) was confirmed
empirically and is reported honestly rather than papered over. See the full
3-round account in
``docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md``
(redesign addendum). Summary:

1. **``deep_hill`` band**, position range ``[LO_DEEP, LO_R]`` (0.002 AU ..
   1e6 km): the originally-mistaken floor's own territory (see ``LO_DEEP``'s
   comment for the full diagnosis of why this basin is trivial). Its own
   partition, ``DEEP_HILL``. Judged ONLY by :func:`classify_bounded_drift`
   plus the corpus anchors filed in ``literature_check.py`` (Henon family-f
   distant-retrograde-orbit physics, Sun-Earth co-orbital dynamics) --
   NEVER by the Table 3/4 family-match criterion, which would just always,
   uninformatively, MISS (nothing published lives in this band).
2. **``paper`` band**: position range per partition CORRECTLY matches stage
   2's own per-set convention (verified against every ``TABLE34`` IC's own
   per-component values, not just its Euclidean radial magnitude -- see the
   redesign note's table): ``[LO_R, HI_R]`` undivided for the 1-free-
   position-component families (A-F), ``[LO_R/sqrt(2), HI_R/sqrt(2)]`` for
   G/H/I/N, ``[LO_R/sqrt(3), HI_R/sqrt(3)]`` (further bracketed per-family
   for J/K/L, which would otherwise collide) for J/K/L/M -- a bug CAUGHT
   during this redesign, not inherited: an earlier draft used bare
   ``[LO_R, HI_R]`` everywhere and, checked directly against every TABLE34
   IC's per-component values, was found to EXCLUDE G/H/K/L/M's own published
   ICs entirely (their true stage-2 per-component floor was ``LO_R`` divided
   by the paper's own sqrt(2)/sqrt(3) norm-sharing divisor, not bare
   ``LO_R``). ``A``-``F`` are FULLY single-family partitions (not pooled at
   all, not even the A/B/C-together grouping stage 2's own set 1 used) --
   see point 4 below for why.
3. **``beyond_hi_r`` band**, position range ``[HI_R, HI_OUTER]`` (``HI_OUTER``
   is the old ``HI_WIDE`` ceiling, renamed for clarity now that it only
   bounds this one band): genuinely uncharted territory past every
   published family's footprint. Its own partition, ``BEYOND_HI_R``, judged
   the same way as ``deep_hill`` (drift classifier + corpus anchors, not
   Table 3/4 match) -- "no new family found here" is an expected, legitimate,
   honest outcome (this territory is already covered by the Henon family-f /
   Sun-Earth co-orbital anchors), not a pipeline failure.
4. **Empirical finding from the smoke test (reported honestly, not spun):**
   narrowing the box (fixing both bugs above, PLUS a velocity-sign split
   mirroring stage 2's DRO/DPO branches) did NOT restore multi-family
   recovery when >1 family shared a population -- 2 corrected multi-family
   attempts (a 3-family ABC/DEF pair with an accidental byte-identical-bounds
   bug, then a properly disjoint, bug-fixed 4-family ABCF box) BOTH still
   collapsed the population onto family C alone (1/N every time, N=3,4,6
   across old P1 and both redesign attempts). A follow-up diagnostic (varied
   random seed, varied box width, theta0 fixed vs free) traced this to
   something deeper than box width: Eq. 15's landscape has NO per-family
   structure once a box contains more than one family's neighborhood -- 2 of
   3 fresh seeds converged almost exactly onto family C's own IC even when
   SOLELY targeting family A, and even a box tightly bracketed around ONLY
   family A's own position range converged to the box's own edge (smaller
   x), never to A's specific published point, REGARDLESS of whether theta0
   was free or fixed at the paper's own value. This directly confirms
   Fable's own second-round diagnosis: Eq. 15 is a boundedness measure that
   saturates near 1.0 across the whole bounded continuum, so nothing in the
   objective steers a population toward a SPECIFIC published family once
   more than one is reachable. Single-family granularity for A-F (fully
   split, not stage 2's own 3-family A/B/C grouping) is adopted as the
   safest available mitigation (removes the family-vs-family competition
   pathway entirely), but is NOT proven, by this dispatch's evidence, to
   guarantee recovery of each family's SPECIFIC published point on any given
   seed -- the coordinator's eventual full sweep should run multiple seeds
   per partition and treat "recovers >=1 known family cleanly" (this
   dispatch's own already-accepted #583 gate) as the realistic per-partition
   bar, not "recovers every family in its group."

## 2026-07-14 #586 follow-up: multi-seed + cluster-everything harvesting

The redesign above (single-family partitions) proved NECESSARY but not
SUFFICIENT: even a partition targeting exactly one family did not reliably
converge to that family's own published point on every seed (module
docstring point 4, above). A Fable design review (recorded in full in
``data/OUTSTANDING.md``'s ``#586`` entry) settled the follow-up: (a) run
multiple INDEPENDENT seeds per partition, and reframe per-partition
family-match recovery as a POSITIVE-CONTROL / machinery-health signal, not
the novelty deliverable -- (b) a target-proximity fitness augmentation was
explicitly REJECTED (it would duplicate stage 2 while destroying this
track's one distinctive property, unaimed exploration); (c) the
seed-conditional limitation is documented, not engineered away.

The actual novelty deliverable is **cluster-everything harvesting**
(mirroring ``run_582_asymmetric_3d_niching_search.py``'s own
``--mode analyze``/:func:`cluster_representatives` pattern): for EACH
seed's finished population, cluster into distinct high-fitness basins (not
just the single best member, not just a check against the partition's own
named target family), then run every cluster representative through the
drift classifier (bounded vs divergent) PLUS the family-match check against
ALL 14 published families (not just the partition's own target -- a cluster
found while searching for family A could plausibly land near family C's own
footprint instead, and that is expected, not a bug, given the redesign's own
finding that Eq. 15 has no per-family structure). Any BOUNDED cluster that
matches NONE of the 14 families is populated into a
``literature_check.CandidateSignature`` (:func:`_build_er3bp_candidate_signature`)
and checked for matcher engagement
(:func:`_er3bp_literature_anchors_engaged`) -- this gets a candidate to
"bounded, unmatched, literature-matcher-engaged", which is necessary but NOT
sufficient for any "novel" claim per
[[feedback_literature_novelty_check_baseline]]; the live
``search/literature_check.py::check_literature()`` web search is explicitly
NOT run here (see :func:`harvest_seed`'s own ``literature_check_status``
field).

**Seeding convention** (:func:`_seed_for`): ``seed_idx=0`` reproduces the
ORIGINAL single-seed base seed (``583000 + partition_index``) byte-for-byte
(the historical positive-control run, e.g. the committed
``C_final.npz``/``C_checkpoint.npz``/``C_analysis_summary.json`` from the
2026-07-13 dispatch); ``seed_idx>=1`` offsets by ``SEED_STRIDE`` (a prime
comfortably larger than the 16-partition index range) so no ``(partition,
seed_idx)`` pair can ever collide. The multi-seed run/harvest paths ALWAYS
use seed-suffixed filenames (``{name}_seed{i}_*``), including seed 0 --
recomputed fresh under the new naming rather than aliased to the legacy
unsuffixed files, keeping the two artifact generations unambiguous. The
legacy unsuffixed ``run_partition``/``analyze_partition``/``--analyze`` path
is left fully intact for backward compatibility with those already-committed
single-seed artifacts.

See ``docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md``
(#586 addendum) for the small-scale validation result.

Usage (chunked, same convention as #581's stage-2 script):
    # legacy single-seed run/analyze (original #583 positive control):
    uv run python scripts/run_583_widened_bounded_drift_search.py \\
        --partition C [--max-gens 100] --workers 8
    uv run python scripts/run_583_widened_bounded_drift_search.py --analyze --partition C

    # #586 multi-seed run (default --n-seeds 3; runs seed_idx 0..N-1 in one call):
    uv run python scripts/run_583_widened_bounded_drift_search.py \\
        --partition C --n-seeds 3 [--max-gens 100] --workers 8
    # or one seed at a time (chunking across turns/coordinator dispatches):
    uv run python scripts/run_583_widened_bounded_drift_search.py \\
        --partition C --seed-index 1 [--max-gens 100] --workers 8

    # #586 cluster-everything harvest + aggregate across all --n-seeds seeds:
    uv run python scripts/run_583_widened_bounded_drift_search.py \\
        --harvest --partition C --n-seeds 3

Checkpoints/results: data/found/583_widened_search/
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from cyclerfinder.core.er3bp_geocentric import (
    SUN_EARTH_ER3BP,
    gurfil_kasdin_fitness,
    table_interleaved_to_state,
)
from cyclerfinder.data.preflight import MethodCapability, preflight_search
from cyclerfinder.data.validation.er3bp_drift_classifier import (
    N_REVS_DEFAULT,
    classify_bounded_drift,
    spot_check_theta0_robustness,
)
from cyclerfinder.search.literature_check import CandidateSignature
from cyclerfinder.search.niching_ga import (
    DeterministicCrowdingConfig,
    run_deterministic_crowding,
)

# Reuse stage 2's own match criterion + dynamics characterizer UNMODIFIED --
# the widened search must be judged by the SAME non-bit-exact bar, not a new
# one invented for this task. Same repo-root sys.path convention as other
# scripts/-to-scripts/ reuse in this project (see pyproject.toml's mypy
# overrides comment, e.g. run_392_v4_annual_sweep.py).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.run_581_gurfil_reproduction import (
    FEATURE_FACTOR_TOL,
    HI_R,
    IC_PROXIMITY_TOL,
    LO_R,
    S2,
    S3,
    TABLE34,
    characterize,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "found" / "583_widened_search"

_REGION_ID = "sun-earth-er3bp-widened-bounded-drift-583-positive-control-2026-07-13"
_METHOD = MethodCapability(
    genome=(
        "Gurfil & Kasdin (2002) geocentric ER3BP objective (Eq. 15, "
        "core/er3bp_geocentric.py) evaluated by deterministic-crowding niching GA "
        "(search/niching_ga.py) over a genuinely WIDENED 7-dim genome (6 "
        "interleaved state slots + free theta0), radially partitioned into a "
        "deep-Hill band, a paper-footprint band (single-family partitions "
        "with correct per-family position bounds), and a beyond-HI_R band"
    ),
    corrector=(
        "no corrector -- this dispatch validates the redesigned partitioning "
        "at smoke scale (3 rounds, see the module docstring's redesign "
        "section) against the old pooled-6-family P1 (1/6); found that box "
        "narrowing alone does not restore multi-family recovery under this "
        "objective (an honestly-reported empirical limit, not a bug), and "
        "adopted single-family granularity as the safest mitigation. The "
        "full novel-territory sweep across all 16 partitions is explicitly "
        "deferred to a future dispatch"
    ),
    capability_tags=frozenset(
        {
            "er3bp",
            "geocentric",
            "niching-ga",
            "widened-domain",
            "positive-control",
            "sun-earth",
        }
    ),
    git_sha="working-tree",
)

# ---------------------------------------------------------------------------
# Radial band + widened bound constants (see module docstring's 2026-07-13
# redesign section for the derivation).
# ---------------------------------------------------------------------------
LO_DEEP = 0.002  # AU (~299,200 km): the ORIGINALLY-mistaken floor from the
# first #583 attempt -- do NOT use this as a general-purpose lower bound. A
# 0.002 AU floor re-admits a physically trivial, strictly-dominant
# quasi-circular deep-Hill-sphere basin: Eq. 15 rewards only 1yr annulus
# thinness with no periodicity/family content, and deep inside the Hill
# sphere essentially ANY near-circular orbit scores ~1 (deficit ~5.6e-9 for
# the trivial basin vs ~8.7e-7 for the best genuine target family, a 150x
# fitness gap) -- deterministic crowding cannot protect a lower-fitness niche
# against a larger, smoothly-reachable, strictly-higher-fitness basin under
# child>=parent replacement, so a population searching THIS RANGE ALONE (not
# mixed with the paper footprint) collapses onto the trivial solution. See
# docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md for the
# full diagnosis. Kept here ONLY as the floor of the dedicated ``DEEP_HILL``
# partition below, which is judged by drift classification + corpus anchors,
# never by Eq. 15 rank (uninformative in this band) or Table 3/4 match
# (nothing published lives here).
# LO_R, HI_R (imported above): the paper's own bounds (1e6 km .. 1e7 km),
# CORRECT ONLY as a per-component floor/ceiling for the 1-free-position-
# component families (A-F, stage-2 sets 1-4, which TABLE2 never divided).
# For every OTHER family, stage 2's OWN Table 2 divided BOTH lo and hi by
# sqrt(2) or sqrt(3) (a shared "vector norm budget" across however many
# genes -- position AND velocity -- that set's box left simultaneously
# free): G/H/I/N used sqrt(2); J/K/L/M used sqrt(3). Checked numerically
# against every TABLE34 published IC (see the redesign note's table) --
# blindly reusing bare (LO_R, HI_R) for these components would EXCLUDE
# several families' own published state entirely (G's x=y=0.004727 AU sits
# 29% below LO_R; H's x=0.005013 sits 25% below; K's y/z components sit
# far below or slightly negative; L's x/y sit 42% below; M's x/y/z all sit
# 42% below) -- a hard correctness bug, not just a niche-capacity dilution,
# discovered during this redesign (the ORIGINAL #583 build's P2-P7 were
# NEVER actually run to completion, only P1, so this never surfaced before).
HI_OUTER = 0.15  # AU (~2.244e7 km): ~2.24x farther than the paper's HI_R;
# 30% of escape_radius (0.5 AU) -- renamed from the original HI_WIDE now that
# it bounds only the dedicated BEYOND_HI_R partition (genuinely uncharted
# territory past every published family), not a position ceiling shared by
# every partition regardless of which families it targets.
V_WIDE = 0.5  # AU/rad (~15 km/s): same magnitude as paper's V, full signed range, per-component.
# Velocity stays UNDIVIDED (full V_WIDE) even for the sqrt(2)/sqrt(3) groups:
# unlike position, no published family's own velocity component is anywhere
# near +-0.5 AU/rad, so an undivided full-width velocity bound is a genuine
# widening with no correctness risk (verified against every TABLE34 IC).
THETA0_LO, THETA0_HI = 0.0, 2.0 * math.pi

_POS_DEEP = (LO_DEEP, LO_R)  # deep_hill band: known Henon-f/co-orbital territory
_POS_S1 = (LO_R, HI_R)  # 1-free-position-component families: A-F (no division in stage 2)
_POS_S2 = (LO_R / S2, HI_R / S2)  # G, H, I, N's own stage-2 divisor
_POS_S3 = (LO_R / S3, HI_R / S3)  # J, K, L, M's own stage-2 divisor (full shared envelope)
_POS_OUTER = (HI_R, HI_OUTER)  # beyond_hi_r band: genuinely uncharted territory
_VEL = (-V_WIDE, V_WIDE)
_FIX = (0.0, 0.0)
_VEL_NEG = (-V_WIDE, 0.0)  # DRO-branch sign restriction (A/B/C/F): matches stage 2's own sets 1,4
_VEL_POS = (0.0, V_WIDE)  # DPO-branch sign restriction (D/E): matches stage 2's own sets 2,3

# J, K, L all share the SAME (0,1,2,3,4,5) all-6-state-free signature and the
# SAME sqrt(3)-divided envelope (_POS_S3) -- if all 3 partitions also shared
# the identical position range, they would be byte-identical search boxes
# (3 seeded replicates of the SAME landscape, not narrower boxes at all; this
# is exactly the byte-identical-bounds bug this dispatch's own smoke test
# caught, see the module docstring's redesign section). Bracket each
# family's OWN published per-component values (verified against TABLE34,
# see the redesign note) with a modest margin, clipped inside _POS_S3's
# envelope, so the 3 partitions are genuinely distinct (if overlapping --
# that's fine and expected, mirroring how several stage-2 families
# genuinely share territory, e.g. C/D/E/F/G/M/N all sit near LO_R and are
# discriminated by signature/velocity-sign/theta0 instead of position alone).
_POS_J = (0.006, HI_R / S3)  # J's own x=0.0335, y=0.0077, z=0.0368 (min 0.0077)
_POS_K = (-0.002, 0.015)  # K's own x=0.0058, y=0.0000547, z=-0.0009 (irregular per TABLE2 set9)
_POS_L = (0.003, 0.018)  # L's own x=y=0.0039, z=0.0095


def _bounds7(state6_bounds: list[tuple[float, float]]) -> list[tuple[float, float]]:
    return [*state6_bounds, (THETA0_LO, THETA0_HI)]


# Per-family free-dim signatures, generalized from stage 2's own TABLE2
# structure (see the redesign note for the rpos/vmag table that justifies
# each grouping). ``pos`` is the partition-specific position bound; ``vel``
# defaults to the full-signed _VEL (genuine widening) but the A-F group
# overrides it to a sign-restricted half-range matching stage 2's own DRO
# (vy<0: A/B/C/F) vs DPO (vy>0: D/E) branches -- an earlier draft that gave
# every A-F partition the SAME full-signed velocity made 2 of them
# BYTE-IDENTICAL (same position range, same free dims), i.e. not narrower at
# all, just 2 seeded replicates of the SAME landscape (the bug this
# dispatch's first smoke-test iteration caught).
def _sig_x_vy(
    pos: tuple[float, float], vel: tuple[float, float] = _VEL
) -> list[tuple[float, float]]:
    """(0,3): free x, vy -- A/B/C/F and D/E's own stage-2 signature."""
    return _bounds7([pos, _FIX, _FIX, vel, _FIX, _FIX])


def _sig_x_y(pos: tuple[float, float]) -> list[tuple[float, float]]:
    """(0,2): free x, y -- G's own stage-2 signature."""
    return _bounds7([pos, _FIX, pos, _FIX, _FIX, _FIX])


def _sig_x_vx(pos: tuple[float, float]) -> list[tuple[float, float]]:
    """(0,1): free x, vx -- H's own stage-2 signature."""
    return _bounds7([pos, _VEL, _FIX, _FIX, _FIX, _FIX])


def _sig_x_vx_y_vy(pos: tuple[float, float]) -> list[tuple[float, float]]:
    """(0,1,2,3): free x, vx, y, vy -- I's own stage-2 signature."""
    return _bounds7([pos, _VEL, pos, _VEL, _FIX, _FIX])


def _sig_full(pos: tuple[float, float]) -> list[tuple[float, float]]:
    """(0,1,2,3,4,5): all 6 state slots free -- J/K/L's own stage-2 signature."""
    return _bounds7([pos, _VEL, pos, _VEL, pos, _VEL])


def _sig_x_y_z(pos: tuple[float, float]) -> list[tuple[float, float]]:
    """(0,2,4): free x, y, z -- M's own stage-2 signature."""
    return _bounds7([pos, _FIX, pos, _FIX, pos, _FIX])


def _sig_x_z(pos: tuple[float, float]) -> list[tuple[float, float]]:
    """(0,4): free x, z -- N's own stage-2 signature."""
    return _bounds7([pos, _FIX, _FIX, _FIX, pos, _FIX])


@dataclass(frozen=True)
class PartitionSpec:
    """One search box: its genome bounds, target families, and band."""

    bounds7: list[tuple[float, float]]  # [x,vx,y,vy,z,vz,theta0]
    families: tuple[str, ...]  # TABLE34 keys this partition is judged against;
    # EMPTY for deep_hill/beyond_hi_r (judged by drift classifier only, never
    # by Table 3/4 match -- nothing published lives in those bands).
    band: str  # "paper" | "deep_hill" | "beyond_hi_r"
    state_sig: tuple[int, ...]  # documentary: which of the 6 state slots are free


# 2026-07-13 redesign: replaces the original 7 wide P1-P7 partitions (see the
# module docstring's redesign section + docs/notes/2026-07-13-583-corpus-
# anchors-and-drift-classifier.md for the full 3-round derivation: the
# position-floor correctness fix, the byte-identical-bounds bug, and the
# smoke-scale finding that multi-family pooling doesn't recover under this
# objective regardless of box width). 14 paper-footprint single-family
# partitions (A-N, each its own population -- see below for why A-F went all
# the way to single-family, not stage 2's own 3-family A/B/C grouping) + 1
# deep-Hill partition + 1 beyond-HI_R partition (both judged by drift
# classification + corpus anchors, never Table 3/4 match).
#
# A-F: FULLY split to single-family granularity. Empirically forced by this
# dispatch's OWN smoke tests: 2 corrected multi-family attempts (a 3-family
# ABC/DEF pair that turned out to have an accidental byte-identical-bounds
# bug, then a properly disjoint, velocity-sign-restricted 4-family ABCF box)
# BOTH still collapsed the ENTIRE population onto family C alone (std hit
# EXACTLY 0.0 by generation ~120 in the corrected ABCF run) -- box-narrowing
# does not, by itself, restore multi-family recovery under deterministic
# crowding on this objective. Splitting fully to single-family removes the
# family-vs-family competition pathway (the most likely fixable part of the
# problem), though a follow-up seed/theta0/box-width diagnostic (see the
# module docstring's redesign section, point 4) found even single-family,
# tightly-bracketed boxes don't reliably converge to a family's SPECIFIC
# published point either -- Eq. 15's landscape appears to have no per-family
# structure at all, a limitation of the (unmodified, per #583's own scope)
# objective itself, not of this partitioning.
PARTITIONS: dict[str, PartitionSpec] = {
    "A": PartitionSpec(_sig_x_vy(_POS_S1, _VEL_NEG), ("A",), "paper", (0, 3)),
    "B": PartitionSpec(_sig_x_vy(_POS_S1, _VEL_NEG), ("B",), "paper", (0, 3)),
    "C": PartitionSpec(_sig_x_vy(_POS_S1, _VEL_NEG), ("C",), "paper", (0, 3)),
    "D": PartitionSpec(_sig_x_vy(_POS_S1, _VEL_POS), ("D",), "paper", (0, 3)),
    "E": PartitionSpec(_sig_x_vy(_POS_S1, _VEL_POS), ("E",), "paper", (0, 3)),
    "F": PartitionSpec(_sig_x_vy(_POS_S1, _VEL_NEG), ("F",), "paper", (0, 3)),
    "G": PartitionSpec(_sig_x_y(_POS_S2), ("G",), "paper", (0, 2)),
    "H": PartitionSpec(_sig_x_vx(_POS_S2), ("H",), "paper", (0, 1)),
    "I": PartitionSpec(_sig_x_vx_y_vy(_POS_S2), ("I",), "paper", (0, 1, 2, 3)),
    "J": PartitionSpec(_sig_full(_POS_J), ("J",), "paper", (0, 1, 2, 3, 4, 5)),
    "K": PartitionSpec(_sig_full(_POS_K), ("K",), "paper", (0, 1, 2, 3, 4, 5)),
    "L": PartitionSpec(_sig_full(_POS_L), ("L",), "paper", (0, 1, 2, 3, 4, 5)),
    "M": PartitionSpec(_sig_x_y_z(_POS_S3), ("M",), "paper", (0, 2, 4)),
    "N": PartitionSpec(_sig_x_z(_POS_S2), ("N",), "paper", (0, 4)),
    "DEEP_HILL": PartitionSpec(_sig_x_vy(_POS_DEEP), (), "deep_hill", (0, 3)),
    "BEYOND_HI_R": PartitionSpec(_sig_x_vy(_POS_OUTER), (), "beyond_hi_r", (0, 3)),
}

# Known-territory notes for the 2 non-"paper" bands (printed + written into
# the analysis summary in lieu of a Table 3/4 match -- see literature_check.py
# KNOWN_CORPUS for the actual anchors these reference).
_KNOWN_TERRITORY_NOTE = {
    "deep_hill": (
        "deep_hill band [LO_DEEP, LO_R] is NOT judged against Table 3/4 "
        "(no published Gurfil-Kasdin family lives this close in) -- this is "
        "the near-circular deep-Hill-sphere continuum, physically the same "
        "territory as Henon family-f distant retrograde orbits / Sun-Earth "
        "co-orbital dynamics (literature_check.py KNOWN_CORPUS anchors). A "
        "clean 'nothing novel here' result is expected and legitimate, not a "
        "failure -- report bounded-fraction + clustering only."
    ),
    "beyond_hi_r": (
        "beyond_hi_r band [HI_R, HI_OUTER] is NOT judged against Table 3/4 "
        "(no published Gurfil-Kasdin family reaches this far out) -- this is "
        "genuinely uncharted position range past every published family's "
        "footprint, but still physically continuous with the same "
        "co-orbital/DRO continuum the corpus anchors already cover. A clean "
        "'nothing novel here' result is expected and legitimate, not a "
        "failure -- report bounded-fraction + clustering only."
    ),
}


def fitness_widened(vec7: np.ndarray) -> float:
    """GA fitness: decode the 7-gene [state6_interleaved, theta0] genome."""
    vec7 = np.asarray(vec7, dtype=float)
    state = table_interleaved_to_state(vec7[:6])
    theta0 = float(vec7[6])
    return gurfil_kasdin_fitness(
        state, theta0, SUN_EARTH_ER3BP, n_rev=1.0, rtol=1e-9, atol=1e-9, n_samples=2000
    )


def run_partition(name: str, max_gens: int | None, workers: int) -> None:
    """Legacy SINGLE-seed run (original #583 positive-control naming, unsuffixed).

    Kept unmodified for backward compatibility with the already-committed
    ``{name}_checkpoint.npz``/``{name}_final.npz`` artifacts (e.g. ``C``,
    ``P1``, ``ABCF``). New work should use :func:`run_partition_multiseed`
    (#586) instead -- see the module docstring's 2026-07-14 section.
    """
    spec = PARTITIONS[name]
    bounds, families = spec.bounds7, spec.families
    config = DeterministicCrowdingConfig(seed=583000 + sorted(PARTITIONS).index(name))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ckpt = OUT_DIR / f"{name}_checkpoint.npz"
    runlog = OUT_DIR / f"{name}_runlog.jsonl"
    t_start = time.monotonic()

    def progress(stats: dict[str, float]) -> None:
        rec = {
            "ts": datetime.now(UTC).isoformat(timespec="seconds"),
            "partition": name,
            **stats,
            "elapsed_s": round(time.monotonic() - t_start, 1),
        }
        with runlog.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        gen = int(stats["generation"])
        if gen % 20 == 0 or gen == config.generations:
            print(
                f"[{name}] gen {gen}/{config.generations} "
                f"mean={stats['fitness_mean']:.6f} max={stats['fitness_max']:.6f} "
                f"std={stats['fitness_std']:.6f} elapsed={rec['elapsed_s']}s",
                flush=True,
            )

    result = run_deterministic_crowding(
        fitness_widened,
        bounds,
        config,
        workers=workers,
        checkpoint_path=ckpt,
        max_generations_this_call=max_gens,
        progress_fn=progress,
    )
    fam_desc = families if families else f"none -- band={spec.band} (known territory)"
    print(
        f"[{name}] done through generation {result.generations_run} (target families: {fam_desc})",
        flush=True,
    )
    if result.generations_run >= config.generations:
        np.savez(
            OUT_DIR / f"{name}_final.npz",
            phenotypes=result.phenotypes,
            fitness=result.fitness,
        )
        print(f"[{name}] final population saved", flush=True)


# ---------------------------------------------------------------------------
# #586: multi-seed run support.
# ---------------------------------------------------------------------------

SEED_STRIDE = 1009
"""Per-seed-index offset (#586). A prime comfortably larger than the
16-partition index range (0-15): ``seed_idx=0`` reproduces the legacy
``583000 + partition_index`` base seed byte-for-byte, and every
``(partition_index, seed_idx)`` pair across all 16 partitions and any
realistic number of seeds maps to a distinct integer seed (no collision is
possible below ``seed_idx = SEED_STRIDE / 16 ~= 63``, far past any planned
``--n-seeds``)."""

DEFAULT_N_SEEDS = 3
"""Per #586's Fable design review: 3 independent seeds is the recommended
default -- cheap enough to run per partition (~2-4 min/seed per #583's own
timing note) while giving the cluster-everything harvest 3 independent
samples of basin space rather than 1."""


def _seed_for(name: str, seed_idx: int) -> int:
    """The GA seed for one ``(partition, seed_idx)`` pair. See :data:`SEED_STRIDE`."""
    return 583000 + sorted(PARTITIONS).index(name) + seed_idx * SEED_STRIDE


def run_partition_seed(name: str, seed_idx: int, max_gens: int | None, workers: int) -> None:
    """One independent-seed niching-GA run for ``name`` (#586).

    Always uses seed-suffixed filenames (``{name}_seed{seed_idx}_*``), even
    for ``seed_idx=0`` -- see the module docstring's 2026-07-14 section for
    why this is kept distinct from the legacy unsuffixed :func:`run_partition`
    artifacts rather than aliased to them.
    """
    spec = PARTITIONS[name]
    bounds, families = spec.bounds7, spec.families
    config = DeterministicCrowdingConfig(seed=_seed_for(name, seed_idx))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"{name}_seed{seed_idx}"
    ckpt = OUT_DIR / f"{tag}_checkpoint.npz"
    runlog = OUT_DIR / f"{tag}_runlog.jsonl"
    t_start = time.monotonic()

    def progress(stats: dict[str, float]) -> None:
        rec = {
            "ts": datetime.now(UTC).isoformat(timespec="seconds"),
            "partition": name,
            "seed_index": seed_idx,
            "seed": config.seed,
            **stats,
            "elapsed_s": round(time.monotonic() - t_start, 1),
        }
        with runlog.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        gen = int(stats["generation"])
        if gen % 20 == 0 or gen == config.generations:
            print(
                f"[{tag}] gen {gen}/{config.generations} "
                f"mean={stats['fitness_mean']:.6f} max={stats['fitness_max']:.6f} "
                f"std={stats['fitness_std']:.6f} elapsed={rec['elapsed_s']}s",
                flush=True,
            )

    result = run_deterministic_crowding(
        fitness_widened,
        bounds,
        config,
        workers=workers,
        checkpoint_path=ckpt,
        max_generations_this_call=max_gens,
        progress_fn=progress,
    )
    fam_desc = families if families else f"none -- band={spec.band} (known territory)"
    print(
        f"[{tag}] done through generation {result.generations_run} "
        f"(seed={config.seed}, target families: {fam_desc})",
        flush=True,
    )
    if result.generations_run >= config.generations:
        np.savez(
            OUT_DIR / f"{tag}_final.npz",
            phenotypes=result.phenotypes,
            fitness=result.fitness,
        )
        print(f"[{tag}] final population saved", flush=True)


def run_partition_multiseed(
    name: str,
    *,
    n_seeds: int,
    seed_index: int | None,
    max_gens: int | None,
    workers: int,
) -> None:
    """Run one or all of ``name``'s ``--n-seeds`` independent seeds (#586).

    ``seed_index`` (if given) restricts this call to that ONE seed (chunking
    across turns/dispatches, same convention as ``--max-gens``); otherwise
    every ``seed_idx`` in ``range(n_seeds)`` is run sequentially in this call
    (each partition/seed pair is fast -- ~2-4 min -- so this stays well within
    a single foreground turn for the default ``n_seeds=3``).
    """
    seed_indices = [seed_index] if seed_index is not None else list(range(n_seeds))
    for si in seed_indices:
        run_partition_seed(name, si, max_gens, workers)


# ---------------------------------------------------------------------------
# Analysis: stage 2's own match criterion, generalized to the 7-gene genome.
# ---------------------------------------------------------------------------


def match_family_in_widened_population(
    fam: str,
    partition: str,
    phen: np.ndarray,
    fitness: np.ndarray,
    *,
    precomputed_characterization: dict[str, object] | None = None,
) -> dict[str, object]:
    """Family-directed match against the 7-gene widened population.

    Same three-part criterion as stage 2's own
    ``run_581_gurfil_reproduction.py::match_family_in_population``
    (IC proximity + type match + rmin/rmax feature ratio), generalized to
    include theta0 as a compared gene (stage 2 fixed theta0 per set so it
    never entered the distance; here it is free, so a genuine match must
    also land near the published theta0).

    ``partition`` selects whose bounds normalize the IC-proximity distance --
    NOT necessarily the partition the candidate was actually found in. #586's
    cluster-everything harvest calls this with ``partition == fam`` (every
    ``fam in TABLE34`` has a same-named single-family partition, see
    :data:`PARTITIONS`) for EVERY one of the 14 families against a single
    harvested cluster representative, checking whether that candidate --
    found while searching for a possibly DIFFERENT family, or even in the
    ``deep_hill``/``beyond_hi_r`` bands -- happens to sit near a different
    family's own published point.

    ``precomputed_characterization`` (#586): :func:`characterize` is a 5-year
    high-accuracy propagation, independent of ``fam``/``partition`` -- when a
    caller already has it for this exact candidate (e.g. checking one cluster
    representative against all 14 families in a loop), passing it here avoids
    14x redundant propagation of the SAME state. ``None`` (the default)
    preserves the original behaviour exactly (compute it fresh from the
    population's own nearest-IC member).
    """
    bounds = PARTITIONS[partition].bounds7
    icv, ic_theta0, ftype, rmin_km, rmax_km = TABLE34[fam]
    target7 = np.array([*icv, ic_theta0])
    free = [k for k, (lo, hi) in enumerate(bounds) if hi > lo]
    span = np.array([bounds[k][1] - bounds[k][0] for k in free])
    lo_v = np.array([bounds[k][0] for k in free])
    target_norm = (target7[free] - lo_v) / span
    norm = (phen[:, free] - lo_v) / span
    diff = norm - target_norm
    if 6 in free:
        # theta0 (index 6) is periodic over its own bounds span -- a linear
        # normalized difference spuriously penalizes a genuine near-target
        # candidate that wrapped past 2*pi back toward 0 (found by Fable
        # review 2026-07-13). Wrap the normalized difference into [-0.5, 0.5]
        # (circular distance in unit-period coordinates) before combining.
        theta_col = free.index(6)
        diff[:, theta_col] -= np.round(diff[:, theta_col])
    dist = np.sqrt(np.mean(diff**2, axis=1))
    i = int(np.argmin(dist))
    cand_theta0 = float(phen[i, 6])
    cand = (
        precomputed_characterization
        if precomputed_characterization is not None
        else characterize(phen[i, :6], cand_theta0)
    )
    got_rmin = float(cand["rmin_km_1yr"])
    got_rmax = float(cand["rmax_km_1yr"])
    ok_ic = float(dist[i]) < IC_PROXIMITY_TOL
    ok_type = cand["type"] == ftype
    ft = FEATURE_FACTOR_TOL
    ok_feat = 1.0 / ft <= got_rmin / rmin_km <= ft and 1.0 / ft <= got_rmax / rmax_km <= ft
    return {
        "family": fam,
        "partition": partition,
        "matched": bool(ok_ic and ok_type and ok_feat),
        "ic_rms_distance": float(dist[i]),
        "ic_ok": bool(ok_ic),
        "type_expected": ftype,
        "type_got": cand["type"],
        "type_ok": bool(ok_type),
        "rmin_ratio": got_rmin / rmin_km,
        "rmax_ratio": got_rmax / rmax_km,
        "features_ok": bool(ok_feat),
        "member_fitness": float(fitness[i]),
        "candidate": cand,
    }


def analyze_partition(name: str) -> None:
    spec = PARTITIONS[name]
    families = spec.families
    fpath = OUT_DIR / f"{name}_final.npz"
    if not fpath.exists():
        print(f"[{name}] final population missing, run to completion first")
        return
    data = np.load(fpath)
    phen, fitness = data["phenotypes"], data["fitness"]

    summary: dict[str, object] = {"partition": name, "band": spec.band, "families": {}}
    if spec.band != "paper":
        # deep_hill / beyond_hi_r: no published family lives in this band, so
        # a Table 3/4 match would always, uninformatively, MISS. Judge only
        # by drift classification + corpus anchors (per Fable's own framing).
        note = _KNOWN_TERRITORY_NOTE[spec.band]
        summary["known_territory_note"] = note
        summary["reproduction_rate"] = "n/a (known territory, not judged against Table 3/4)"
        print(f"[{name}] {note}", flush=True)
    else:
        n_matched = 0
        for fam in families:
            rec = match_family_in_widened_population(fam, name, phen, fitness)
            summary["families"][fam] = rec  # type: ignore[index]
            n_matched += int(rec["matched"])
            print(
                f"[{name}] family {fam}: {'MATCH' if rec['matched'] else 'MISS '} "
                f"ic_dist={rec['ic_rms_distance']:.4f}({'ok' if rec['ic_ok'] else 'FAR'}) "
                f"type={rec['type_got']} vs {rec['type_expected']} "
                f"rmin_ratio={rec['rmin_ratio']:.3f} rmax_ratio={rec['rmax_ratio']:.3f}",
                flush=True,
            )
        summary["reproduction_rate"] = f"{n_matched}/{len(families)}"

    # Bounded-vs-divergent drift classification (#583's own new gate) on
    # every population member above a modest fitness floor -- the positive-
    # control question here is only "does the widened machinery still see
    # the known families", so this is diagnostic breadth, not a novelty claim.
    bounded_ids: list[int] = []
    checked = 0
    for idx, (row, fit) in enumerate(zip(phen, fitness, strict=True)):
        if fit < 0.9:  # skip low-fitness (unbounded/collided) genome noise
            continue
        checked += 1
        state = table_interleaved_to_state(row[:6])
        theta0 = float(row[6])
        v = classify_bounded_drift(state, theta0, n_revs=N_REVS_DEFAULT)
        if v.bounded:
            bounded_ids.append(idx)
    summary["drift_classified"] = {"checked": checked, "bounded": len(bounded_ids)}
    print(
        f"[{name}] drift classification: {len(bounded_ids)}/{checked} high-fitness "
        f"members bounded at {N_REVS_DEFAULT}yr",
        flush=True,
    )

    # Cheap theta0-robustness spot-check (per #583's own mandate) on up to 5
    # bounded survivors: re-test at 2 other phases and flag any flip. This is
    # diagnostic-only (does not itself gate anything) -- capped at 5 to keep
    # this analysis pass bounded regardless of how many members pass drift
    # classification.
    spot_checks: dict[str, object] = {}
    for idx in bounded_ids[:5]:
        row = phen[idx]
        state = table_interleaved_to_state(row[:6])
        theta0 = float(row[6])
        out = spot_check_theta0_robustness(state, theta0, n_revs=N_REVS_DEFAULT)
        flipped = not all(out.values())
        spot_checks[str(idx)] = {"per_phase_bounded": {str(k): v for k, v in out.items()}}
        print(
            f"[{name}] theta0 spot-check member {idx}: "
            f"{'FLIPPED' if flipped else 'stable'} across phases {list(out.values())}",
            flush=True,
        )
    summary["theta0_spot_checks"] = spot_checks

    out = OUT_DIR / f"{name}_analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[{name}] reproduction: {summary['reproduction_rate']}; summary -> {out}")


# ---------------------------------------------------------------------------
# #586: cluster-everything harvesting (mirrors run_582's own --mode analyze /
# cluster_representatives pattern -- see the module docstring's 2026-07-14
# section). The novelty deliverable for this track: cluster EVERY seed's
# final population into distinct high-fitness basins, drift-classify +
# family-match (against ALL 14 families) each representative, and populate a
# literature_check.CandidateSignature for any bounded-but-unmatched cluster.
# ---------------------------------------------------------------------------

DEFAULT_HARVEST_FITNESS_FLOOR = 0.9
"""Same choice, same reason, as #582's own ``DEFAULT_ANALYZE_FITNESS_FLOOR``
and this module's own legacy ``analyze_partition`` (0.9 skips low-fitness
unbounded/collided genome noise while staying well below any corrector's own
tighter closure tolerance -- there is no corrector step in this pipeline, the
drift classifier plays that role)."""

DEFAULT_HARVEST_DISTANCE_THRESHOLD = 0.1
"""Same value as #582's own ``DEFAULT_ANALYZE_DISTANCE_THRESHOLD`` -- coarse
enough to collapse near-duplicate members of one niche, tight enough to keep
genuinely distinct basins separate. Not paper-sourced (no published
clustering recipe for this landscape either); documented and adjustable via
``--distance-threshold``, mirroring #582's own disclosure."""

DEFAULT_HARVEST_MAX_CLUSTERS = 25
"""Same cap as #582's own ``DEFAULT_ANALYZE_MAX_CLUSTERS``."""


def cluster_representatives(
    phen: np.ndarray,
    fitness: np.ndarray,
    bounds: list[tuple[float, float]],
    *,
    fitness_floor: float,
    distance_threshold: float,
    max_clusters: int,
) -> list[int]:
    """Greedy deterministic-crowding-style niche clustering (#586, #582 pattern).

    Walks the population in DESCENDING fitness order (stopping once fitness
    drops below ``fitness_floor``) and accepts a member as a NEW cluster
    representative only if it is more than ``distance_threshold`` away, in
    the bounds-normalized ``[0,1]^7`` Euclidean metric, from EVERY
    representative already accepted. Ported from
    ``run_582_asymmetric_3d_niching_search.py::cluster_representatives``
    (same algorithm, same defaults, generalized to this module's 7-gene
    ``bounds7`` instead of #582's 6-gene bounds).

    Unlike #582's genome, #583's ``bounds7`` routinely has FIXED (``lo==hi``)
    columns (most partitions fix 2-4 of the 6 state slots at 0) -- the
    ``span[span == 0.0] = 1.0`` guard below is actually EXERCISED here (in
    #582 it was defensive-only, that genome has no fixed genes): a fixed
    column's normalized value is identical (``phen - lo == 0``) for every
    population member, so it contributes zero offset to every pairwise
    distance and cannot perturb clustering.
    """
    lo = np.array([b[0] for b in bounds], dtype=np.float64)
    span = np.array([b[1] - b[0] for b in bounds], dtype=np.float64)
    span[span == 0.0] = 1.0
    norm = (phen - lo) / span

    order = np.argsort(-fitness)
    reps: list[int] = []
    rep_vecs: list[np.ndarray] = []
    for idx in order:
        if float(fitness[idx]) < fitness_floor:
            break  # order is descending, so nothing further qualifies either
        v = norm[idx]
        if all(float(np.linalg.norm(v - rv)) > distance_threshold for rv in rep_vecs):
            reps.append(int(idx))
            rep_vecs.append(v)
            if len(reps) >= max_clusters:
                break
    return reps


def _build_er3bp_candidate_signature() -> CandidateSignature:
    """Populate the literature-matcher's structural fingerprint for a bounded,
    Table-3/4-unmatched Sun-Earth ER3BP cluster (#586).

    ``primary="Sun"``, ``sequence=("E",)`` reaches all 3 of the #583-filed
    corpus anchors (Gurfil-Kasdin 2002, Sun-Earth co-orbital/horseshoe, Henon
    1969 family-f -- all three declare ``body_set=frozenset({"E"})`` under
    ``primary="Sun"``, verified directly against
    ``literature_check.py``'s own anchor definitions). ``topology_label`` is
    left at its empty default DELIBERATELY: the one anchor that DOES declare
    a non-empty label (``"binary-coorbital"``, the co-orbital anchor) would
    be silently EXCLUDED by ``_candidate_anchors``'s "both non-empty and
    disjoint" filter if this signature declared any other label -- leaving it
    empty guarantees all 3 target anchors stay reachable.

    Note this also reaches every OTHER Sun-primary corpus anchor whose own
    ``body_set`` merely contains ``"E"`` as one of several bodies (e.g. the
    Earth-Mars interplanetary cycler anchors) -- ``_candidate_anchors``'s
    match is a body-SUBSET test, so a single-body candidate is structurally a
    subset of any multi-body anchor's set. This is a wide, not narrow, net;
    documented rather than worked around, since narrowing it (e.g. via
    ``topology_label``) risks excluding the anchors that actually matter (see
    above) and the live ``check_literature()`` search -- not this structural
    engagement check -- is what actually disambiguates relevance.
    """
    return CandidateSignature(primary="Sun", sequence=("E",), period_years=1.0, n_rev=(1,))


def _er3bp_literature_anchors_engaged(sig: CandidateSignature) -> list[str]:
    """Names of ``literature_check`` corpus anchors this signature would search.

    Mirrors ``isolated_3d_asymmetric_pipeline.py::literature_anchors_engaged``
    (#582's own pattern) exactly: proves the matcher actually engages
    (non-empty return), NOT the live search itself.
    """
    from cyclerfinder.search.literature_check import _candidate_anchors

    return [a.name for a in _candidate_anchors(sig)]


_HARVEST_LITERATURE_CHECK_STATUS = (
    "NOT RUN. Per [[feedback_literature_novelty_check_baseline]], "
    "search/literature_check.py::check_literature() requires a live injected "
    "SearchFn and MUST be run against every cluster with "
    "ready_for_literature_check=true before any 'novel' claim -- this harvest "
    "pass only confirms the structural matcher engages (non-empty anchor pool "
    "via _er3bp_literature_anchors_engaged()) and that the candidate matches "
    "NONE of the 14 published Table 3/4 families; it does not itself search "
    "or adjudicate novelty. Not-found is necessary, not sufficient."
)


def harvest_seed(
    name: str,
    seed_idx: int,
    *,
    fitness_floor: float = DEFAULT_HARVEST_FITNESS_FLOOR,
    distance_threshold: float = DEFAULT_HARVEST_DISTANCE_THRESHOLD,
    max_clusters: int = DEFAULT_HARVEST_MAX_CLUSTERS,
) -> dict[str, object] | None:
    """Cluster-everything harvest of ONE seed's finished population (#586).

    For each cluster representative: classify bounded-vs-divergent
    (:func:`~cyclerfinder.data.validation.er3bp_drift_classifier.classify_bounded_drift`),
    then -- ONLY for bounded members (a divergent candidate cannot be any
    published family or any genuine novel bounded orbit) -- check against
    ALL 14 published families (not just this partition's own named target;
    see :func:`match_family_in_widened_population`'s own docstring). A
    bounded cluster matching NONE of the 14 is populated into a
    :class:`CandidateSignature` and checked for literature-matcher
    engagement. Writes ``{name}_seed{seed_idx}_harvest_summary.json``.
    """
    spec = PARTITIONS[name]
    tag = f"{name}_seed{seed_idx}"
    fpath = OUT_DIR / f"{tag}_final.npz"
    if not fpath.exists():
        print(f"[{tag}] final population missing ({fpath}); run to completion first")
        return None

    data = np.load(fpath)
    phen, fitness = data["phenotypes"], data["fitness"]
    rep_idx = cluster_representatives(
        phen,
        fitness,
        spec.bounds7,
        fitness_floor=fitness_floor,
        distance_threshold=distance_threshold,
        max_clusters=max_clusters,
    )
    n_above_floor = int(np.sum(fitness >= fitness_floor))
    print(
        f"[{tag}] population={phen.shape[0]} n_above_floor({fitness_floor})={n_above_floor} "
        f"clusters_selected={len(rep_idx)}",
        flush=True,
    )

    clusters: list[dict[str, object]] = []
    own_family_matched = False
    n_bounded = 0
    n_unmatched_bounded = 0
    for rank, idx in enumerate(rep_idx):
        genome = phen[idx]
        fit = float(fitness[idx])
        state6 = genome[:6]
        theta0 = float(genome[6])
        state_geo = table_interleaved_to_state(state6)
        verdict = classify_bounded_drift(state_geo, theta0, n_revs=N_REVS_DEFAULT)
        entry: dict[str, object] = {
            "cluster_rank": rank,
            "population_index": int(idx),
            "ga_genome": [float(v) for v in genome],
            "ga_fitness": fit,
            "drift": {
                "bounded": bool(verdict.bounded),
                "growth_ratio": float(verdict.growth_ratio),
                "trend_fraction": float(verdict.trend_fraction),
                "terminated_early": bool(verdict.terminated_early),
                "termination_reason": verdict.termination_reason,
                "n_windows_complete": verdict.n_windows_complete,
            },
        }
        if not verdict.bounded:
            entry["family_matches"] = []
            entry["ready_for_literature_check"] = False
            clusters.append(entry)
            print(
                f"[{tag}] cluster {rank}: DIVERGENT fitness={fit:.6f} "
                f"({verdict.termination_reason or verdict.notes})",
                flush=True,
            )
            continue

        n_bounded += 1
        cand = characterize(state6, theta0)  # computed ONCE, reused for all 14 family checks
        hits: list[str] = []
        for fam in sorted(TABLE34):
            rec = match_family_in_widened_population(
                fam,
                fam,  # normalize against FAM's OWN home partition, not this cluster's partition
                genome.reshape(1, -1),
                np.array([fit]),
                precomputed_characterization=cand,
            )
            if rec["matched"]:
                hits.append(fam)
                if fam in spec.families:
                    own_family_matched = True
        entry["family_matches"] = hits
        entry["candidate_type"] = cand["type"]

        if hits:
            entry["ready_for_literature_check"] = False
            print(
                f"[{tag}] cluster {rank}: BOUNDED fitness={fit:.6f} type={cand['type']} "
                f"matches known family(ies) {hits}",
                flush=True,
            )
        else:
            n_unmatched_bounded += 1
            sig = _build_er3bp_candidate_signature()
            anchors = _er3bp_literature_anchors_engaged(sig)
            entry["candidate_signature"] = {
                "primary": sig.primary,
                "sequence": list(sig.sequence),
            }
            entry["literature_anchors_engaged"] = anchors
            entry["ready_for_literature_check"] = bool(anchors)
            print(
                f"[{tag}] cluster {rank}: BOUNDED fitness={fit:.6f} type={cand['type']} "
                f"UNMATCHED against all 14 families -- literature-check candidate "
                f"(anchors={anchors})",
                flush=True,
            )
        clusters.append(entry)

    summary: dict[str, object] = {
        "partition": name,
        "seed_index": seed_idx,
        "seed": _seed_for(name, seed_idx),
        "band": spec.band,
        "own_target_families": list(spec.families),
        "own_family_recovered_this_seed": own_family_matched,
        "population_size": int(phen.shape[0]),
        "fitness_floor": fitness_floor,
        "distance_threshold": distance_threshold,
        "n_above_floor": n_above_floor,
        "n_clusters": len(clusters),
        "n_bounded": n_bounded,
        "n_unmatched_bounded": n_unmatched_bounded,
        "clusters": clusters,
        "literature_check_status": _HARVEST_LITERATURE_CHECK_STATUS,
    }
    out = OUT_DIR / f"{tag}_harvest_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    print(
        f"[{tag}] harvest: {len(clusters)} cluster(s), {n_bounded} bounded, "
        f"{n_unmatched_bounded} unmatched-bounded; own target family recovered="
        f"{own_family_matched}; summary -> {out}",
        flush=True,
    )
    return summary


def harvest_partition(
    name: str,
    *,
    n_seeds: int,
    seed_index: int | None,
    fitness_floor: float = DEFAULT_HARVEST_FITNESS_FLOOR,
    distance_threshold: float = DEFAULT_HARVEST_DISTANCE_THRESHOLD,
    max_clusters: int = DEFAULT_HARVEST_MAX_CLUSTERS,
) -> None:
    """Harvest one or all of ``name``'s seeds, then aggregate (#586).

    If ``seed_index`` is given, harvests only that one seed and returns (no
    aggregate -- a partial-seed-set call has nothing meaningful to
    aggregate). Otherwise harvests every ``seed_idx`` in ``range(n_seeds)``
    and writes ``{name}_aggregate_harvest_summary.json``: per-seed recovery
    of THIS partition's own target family (the positive-control/machinery-
    health signal, per #586's Fable design review -- NOT the novelty
    deliverable), plus every unmatched-bounded candidate found across all
    seeds (the actual novelty deliverable, still gated on the live
    ``check_literature()`` call this dispatch does not run).
    """
    seed_indices = [seed_index] if seed_index is not None else list(range(n_seeds))
    per_seed = [
        s
        for si in seed_indices
        if (
            s := harvest_seed(
                name,
                si,
                fitness_floor=fitness_floor,
                distance_threshold=distance_threshold,
                max_clusters=max_clusters,
            )
        )
        is not None
    ]
    if seed_index is not None:
        return

    spec = PARTITIONS[name]
    n_recovered = sum(1 for s in per_seed if s["own_family_recovered_this_seed"])
    unmatched: list[dict[str, object]] = []
    for s in per_seed:
        for c in s["clusters"]:  # type: ignore[union-attr]
            if c.get("ready_for_literature_check"):
                unmatched.append({"seed_index": s["seed_index"], **c})

    if spec.band == "paper":
        recovery_note = (
            "Per #586's Fable design review: per-partition/per-seed recovery of "
            "this partition's OWN target family is a POSITIVE-CONTROL/machinery-"
            "health signal, not the novelty deliverable -- gurfil_kasdin_fitness "
            "(Eq. 15) saturates near 1.0 across the whole bounded continuum once "
            "the trivial deep-Hill basin is excluded, so recovery per partition is "
            "SEED-CONDITIONAL under Eq. 15 + deterministic crowding, not a "
            "guarantee on any given seed. A result below n_seeds_run recovering "
            "is NOT, by itself, a pipeline failure -- see the module docstring's "
            "2026-07-14 section."
        )
    else:
        recovery_note = _KNOWN_TERRITORY_NOTE[spec.band]

    aggregate = {
        "partition": name,
        "band": spec.band,
        "own_target_families": list(spec.families),
        "n_seeds_run": len(per_seed),
        "n_seeds_recovered_own_family": n_recovered,
        "recovery_note": recovery_note,
        "unmatched_bounded_candidates": unmatched,
        "n_unmatched_bounded_candidates": len(unmatched),
        "literature_novelty_status": (
            "NOT ADJUDICATED. Per [[feedback_literature_novelty_check_baseline]], "
            "these unmatched-bounded candidates are literature-matcher-ENGAGED "
            "(non-empty anchor pool) but the live check_literature() web search "
            "has NOT been run against any of them -- necessary, not sufficient, "
            "for any novelty claim. Live search + adjudication is a further step "
            "for the coordinator/a future adjudication dispatch, not done here."
        ),
        "per_seed_summaries": per_seed,
    }
    out = OUT_DIR / f"{name}_aggregate_harvest_summary.json"
    out.write_text(json.dumps(aggregate, indent=2, default=str))
    print(
        f"[{name}] AGGREGATE across {len(per_seed)} seed(s): "
        f"{n_recovered}/{len(per_seed)} recovered own target family; "
        f"{len(unmatched)} unmatched-bounded literature-check candidate(s); "
        f"summary -> {out}",
        flush=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--partition", choices=sorted(PARTITIONS), help="run/analyze one partition")
    ap.add_argument("--max-gens", type=int, default=None, help="generation cap this call")
    ap.add_argument("--workers", type=int, default=8, help="process-pool workers (M3: 8 cores)")
    ap.add_argument(
        "--analyze",
        action="store_true",
        help="legacy: analyze a finished SINGLE-seed {name}_final.npz (original "
        "#583 unsuffixed positive-control naming). Superseded by --harvest (#586) "
        "for new work",
    )
    ap.add_argument(
        "--legacy-run",
        action="store_true",
        help="legacy: run the SINGLE-seed unsuffixed positive-control path "
        "(run_partition, original #583 naming) instead of the #586 multi-seed "
        "path -- pairs with --analyze",
    )
    ap.add_argument(
        "--harvest",
        action="store_true",
        help="#586: cluster-everything harvest of finished {name}_seed{i}_final.npz "
        "populations across --n-seeds seeds, instead of running the GA",
    )
    ap.add_argument(
        "--n-seeds",
        type=int,
        default=DEFAULT_N_SEEDS,
        help="#586: number of independent seeds to run/harvest per partition "
        f"(default {DEFAULT_N_SEEDS} per the Fable design review)",
    )
    ap.add_argument(
        "--seed-index",
        type=int,
        default=None,
        help="#586: run/harvest only this one 0-based seed index instead of all "
        "--n-seeds (chunking across turns/dispatches, same convention as --max-gens)",
    )
    ap.add_argument(
        "--fitness-floor",
        type=float,
        default=DEFAULT_HARVEST_FITNESS_FLOOR,
        help="#586 --harvest only: minimum GA fitness to be considered for clustering",
    )
    ap.add_argument(
        "--distance-threshold",
        type=float,
        default=DEFAULT_HARVEST_DISTANCE_THRESHOLD,
        help="#586 --harvest only: bounds-normalized distance a candidate must "
        "exceed from every existing representative to start a new cluster",
    )
    ap.add_argument(
        "--max-clusters",
        type=int,
        default=DEFAULT_HARVEST_MAX_CLUSTERS,
        help="#586 --harvest only: cap on cluster representatives analyzed per seed",
    )
    args = ap.parse_args()
    preflight_search(
        task_no=583,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(PARTITIONS) * args.n_seeds,
        override_reason=(
            "this dispatch validates the 2026-07-13 partition redesign at "
            "smoke scale (3 rounds: ABC/DEF, corrected ABCF, single-family "
            "A) against the old pooled-6-family P1 (1/6), and #586's "
            "follow-up multi-seed + cluster-everything harvest at small "
            "scale (1 partition, 2-3 seeds); the full novel-territory sweep "
            "across all 16 redesigned partitions x --n-seeds is explicitly "
            "deferred to a future coordinator-run dispatch, not an "
            "unbudgeted discovery sweep run here"
        ),
    )
    if args.partition is None:
        ap.error("specify --partition NAME (with --analyze, --harvest, or --max-gens)")
    if args.analyze:
        analyze_partition(args.partition)
    elif args.harvest:
        harvest_partition(
            args.partition,
            n_seeds=args.n_seeds,
            seed_index=args.seed_index,
            fitness_floor=args.fitness_floor,
            distance_threshold=args.distance_threshold,
            max_clusters=args.max_clusters,
        )
    elif args.legacy_run:
        run_partition(args.partition, args.max_gens, args.workers)
    else:
        run_partition_multiseed(
            args.partition,
            n_seeds=args.n_seeds,
            seed_index=args.seed_index,
            max_gens=args.max_gens,
            workers=args.workers,
        )


if __name__ == "__main__":
    main()
