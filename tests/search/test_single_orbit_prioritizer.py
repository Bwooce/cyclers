"""Tests for the single-orbit prioritizer adapter (#310 Phase 1).

The adapter closes the architectural gap surfaced by #284: the #282 five-tier
prioritizer is shaped for orbit PAIRS (representative-pair mode) or
patched-conic Lambert LEGS (tier-0-only mode). A single-orbit discovery
candidate fits NEITHER.

The four smoke / contract tests below exercise the two adaptation strategies
documented in the adapter's module docstring:

1. **Smoke** -- a sourced IC (#287 spike `planar_baseline_z0eq0`, the Braik-
   Ross C11a-cycler at C=3.1294) goes in; the parallel-pipeline tiers
   (4 + 5) must populate when their scorers are provided.
2. **Surrogate-pair** -- the same IC, with the catalogue's pair-shape rows
   in scope; the nearest neighbor is found and reported (the planar Braik-
   Ross C11a row matches the IC verbatim, since the IC IS the state_nd row).
3. **Skip-tiers** -- explicit ``skip_tiers=(0, 1, 2, 3)`` yields a partial
   score with only tiers 4 + 5 populated.
4. **No-neighbor fallback** -- a far-mu candidate (a different CR3BP system)
   yields no surrogate within the bounded distance; tiers 0-3 stay None,
   tiers 4 + 5 still run.

Sourced golden discipline: the smoke IC is loaded from
``data/spike_287.jsonl`` (the #287 spike record, `case=planar_baseline_z0eq0`).
That state itself was reproduced by the #249 corrector against Braik-Ross
2026 Table 2 (T_PO=42.140 d to 0.0011%); the IC is sourced through the
documented reproduction chain. Tests assert STRUCTURE (which tiers
populated, which surrogate was matched) NOT magnitudes (the source papers
don't publish single-orbit single-pair scores for comparison).

Wall-time budget: each test runs in <30 s on the smoke FTLE grid; total
suite <90 s per the workflow gate.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.five_tier_prioritizer as ftp
import cyclerfinder.search.ftle_scorer as fs
import cyclerfinder.search.lobe_overlap_scorer as los
import cyclerfinder.search.reachable_representatives as rr
import cyclerfinder.search.single_orbit_prioritizer as sop

# ---------------------------------------------------------------------------
# Fixtures.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def em_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP system in Braik-Ross / Ross-RT nondimensional scales."""
    return rr.braik_ross_system()


@pytest.fixture(scope="module")
def spike_287_planar_ic() -> dict[str, Any]:
    """Load #287 spike's `planar_baseline_z0eq0` IC (Braik-Ross C11a).

    The state matches the catalogue row ``braik-ross-c11a-cycler-2026``
    ``orbit_elements.cr3bp.state_nd`` verbatim (#249 reproduction). Used
    as the smoke-test IC -- a SOURCED single-orbit candidate.
    """
    spike_path = Path(__file__).resolve().parents[2] / "data" / "spike_287.jsonl"
    with spike_path.open() as f:
        for line in f:
            row = json.loads(line)
            if row.get("case") == "planar_baseline_z0eq0":
                assert isinstance(row, dict)
                return row
    raise RuntimeError(
        "Could not find planar_baseline_z0eq0 row in data/spike_287.jsonl -- "
        "the #287 spike data has changed under test"
    )


@pytest.fixture(scope="module")
def smoke_ftle_field(em_system: cr3bp.CR3BPSystem) -> fs.FTLEField:
    """Coarse FTLE field at C=3.05 over the Braik-Ross interior box.

    Cribbed from ``tests/search/test_ftle_scorer.py::smoke_field``: 8x8
    grid + 2 TU horizon = ~5s. The Braik-Ross C11a IC at x~-0.81 lies
    inside this box (x_bounds=(0.6, 1.4) does NOT contain it; we use a
    wider box to ensure the candidate's position is on-grid).
    """
    return fs.compute_ftle_field(
        em_system,
        c_j=3.05,
        x_bounds=(-1.5, 1.5),
        y_bounds=(-1.0, 1.0),
        grid_shape=(8, 8),
        integration_time_tu=2.0,
    )


@pytest.fixture(scope="module")
def prioritizer_with_scorers(
    em_system: cr3bp.CR3BPSystem,
    smoke_ftle_field: fs.FTLEField,
) -> ftp.FiveTierPrioritizer:
    """FiveTierPrioritizer with Tier 4 (FTLE) and Tier 5 (lobe) attached.

    Tier 0 (NN) is the default vendored model. Tiers 1-3 are NOT attached
    (they require an expensive atlas / continuation; not needed for the
    smoke contract tests). The parallel-pipeline tiers (4 + 5) are.
    """
    return ftp.FiveTierPrioritizer(
        ftle=fs.FTLEScorer(system=em_system, c_j=3.05, ftle_field=smoke_ftle_field),
        lobe=los.LobeOverlapScorer(system=em_system),
    )


# ---------------------------------------------------------------------------
# Test 1: smoke -- the #287 spike IC gets a SingleOrbitScore back.
# ---------------------------------------------------------------------------


def test_smoke_287_spike_returns_score(
    spike_287_planar_ic: dict[str, Any],
    em_system: cr3bp.CR3BPSystem,
    prioritizer_with_scorers: ftp.FiveTierPrioritizer,
) -> None:
    """The #287 planar baseline IC produces a populated SingleOrbitScore.

    Asserts STRUCTURE (the dataclass returns, tier 4 + tier 5 populate)
    rather than magnitudes -- the source papers don't publish a single-
    orbit single-pair score for comparison.
    """
    state0 = np.array(
        [
            spike_287_planar_ic["x0"],
            0.0,
            spike_287_planar_ic.get("z0", 0.0),
            0.0,
            spike_287_planar_ic["ydot0"],
            0.0,
        ],
        dtype=np.float64,
    )
    period_nd = float(spike_287_planar_ic["T_TU"])
    jacobi = float(spike_287_planar_ic["jacobi"])

    score = sop.score_single_orbit(
        candidate_state0=state0,
        candidate_period_nondim=period_nd,
        candidate_system=em_system,
        candidate_id="spike-287-c11a",
        candidate_jacobi=jacobi,
        prioritizer=prioritizer_with_scorers,
        # No max distance -- we want to test that SOME neighbor is found;
        # the surrogate test below pins which one.
        surrogate_neighbor_max_distance=None,
    )
    assert isinstance(score, sop.SingleOrbitScore)
    assert score.candidate_id == "spike-287-c11a"
    # Tier 4 always runs (parallel pipeline, FTLE scorer attached).
    assert score.tier_4_score is not None, (
        "tier 4 (FTLE corridor strength) must populate when an FTLE scorer "
        "is attached -- parallel pipeline doesn't need a surrogate"
    )
    assert 0.0 <= score.tier_4_score <= 1.0, (
        f"tier 4 corridor strength out of [0,1]: {score.tier_4_score}"
    )
    # Tier 5 attempts to run; the lobe pipeline can return zero-flux if the
    # synthetic ResonantMember's manifold doesn't produce overlapping lobes
    # in self-graph mode, but the scorer must NOT crash and the field must
    # be a finite float (or skipped with a notes-recorded reason).
    if score.tier_5_score is not None:
        assert math.isfinite(score.tier_5_score)
    else:
        # If tier 5 was skipped, the reason must be on the notes audit.
        assert "tier5" in score.notes or "lobe" in score.notes.lower(), (
            f"tier 5 returned None without an audit entry; notes={score.notes!r}"
        )
    # Notes must record which strategies ran.
    assert "strategy-2 parallel-pipeline" in score.notes


# ---------------------------------------------------------------------------
# Test 2: surrogate-neighbor lookup -- the C11a row matches.
# ---------------------------------------------------------------------------


def test_surrogate_neighbor_matches_c11a(
    spike_287_planar_ic: dict[str, Any],
    em_system: cr3bp.CR3BPSystem,
) -> None:
    """The C11a IC's nearest catalogue neighbor IS its catalogue-row twin.

    The planar baseline IC of the #287 spike is the verbatim state for the
    ``braik-ross-c11a-cycler-2026`` catalogue row (the row's
    ``state_nd`` was BACKFILLED from the #249 reproduction that the spike
    re-runs). The tuple-distance metric over (mu, C, T) must therefore
    pick that row as the nearest neighbor at zero distance to within
    numerical noise.
    """
    period_nd = float(spike_287_planar_ic["T_TU"])
    jacobi = float(spike_287_planar_ic["jacobi"])
    found = sop.find_surrogate_neighbor(
        candidate_mu=em_system.mu,
        candidate_jacobi=jacobi,
        candidate_period_nondim=period_nd,
    )
    assert found is not None, (
        "surrogate-neighbor search returned None -- catalogue.yaml has no "
        "rows with (mass_ratio, jacobi_constant, period_nd, state_nd) "
        "populated, or the lookup is broken"
    )
    surrogate, distance = found
    # The planar Braik-Ross (1,1) family row matches the spike IC verbatim.
    # We allow any of the Braik-Ross / Ross-RT rows at mu=1.215e-2, C~3.13
    # as a satisfactory match -- the C11a row is the EXACT match.
    assert surrogate.id.startswith("braik-ross-") or surrogate.id.startswith("ross-rt-"), (
        f"nearest neighbor was {surrogate.id!r}, expected a Braik-Ross or "
        "Ross-RT Earth-Moon cycler row"
    )
    # The EXPECTED nearest is C11a (same C=3.1294 and same T=9.69107744 TU).
    assert surrogate.id == "braik-ross-c11a-cycler-2026", (
        f"nearest neighbor was {surrogate.id!r}, expected braik-ross-c11a-"
        f"cycler-2026 (the spike's catalogue-row twin); distance={distance:.3e}"
    )
    # The distance must be small (mu / C / T all match to <0.001%).
    assert distance < 1e-4, (
        f"surrogate distance to C11a is {distance:.3e}; "
        "the spike IC is supposed to be the verbatim C11a state, distance ~ 0"
    )


# ---------------------------------------------------------------------------
# Test 3: skip_tiers=(0, 1, 2, 3) -- partial score with only 4 + 5.
# ---------------------------------------------------------------------------


def test_skip_tiers_yields_partial_score(
    spike_287_planar_ic: dict[str, Any],
    em_system: cr3bp.CR3BPSystem,
    prioritizer_with_scorers: ftp.FiveTierPrioritizer,
) -> None:
    """``skip_tiers=(0, 1, 2, 3)`` collapses to the parallel-pipeline only."""
    state0 = np.array(
        [
            spike_287_planar_ic["x0"],
            0.0,
            spike_287_planar_ic.get("z0", 0.0),
            0.0,
            spike_287_planar_ic["ydot0"],
            0.0,
        ],
        dtype=np.float64,
    )
    period_nd = float(spike_287_planar_ic["T_TU"])
    jacobi = float(spike_287_planar_ic["jacobi"])

    score = sop.score_single_orbit(
        candidate_state0=state0,
        candidate_period_nondim=period_nd,
        candidate_system=em_system,
        candidate_id="spike-287-c11a-skip0123",
        candidate_jacobi=jacobi,
        prioritizer=prioritizer_with_scorers,
        skip_tiers=(0, 1, 2, 3),
    )
    assert score.tier_0_score is None
    assert score.tier_1_score is None
    assert score.tier_2_score is None
    assert score.tier_3_score is None
    # Tier 4 must still populate (FTLE scorer is attached, and tier 4 is not
    # skipped).
    assert score.tier_4_score is not None
    # When all tiers 0-3 are skipped, the surrogate lookup is BYPASSED.
    assert score.surrogate_pair_neighbor_id is None


# ---------------------------------------------------------------------------
# Test 4: no-neighbor fallback -- a tiny-mu system has no near-by catalogue row.
# ---------------------------------------------------------------------------


def test_no_neighbor_fallback_to_parallel_pipeline(
    em_system: cr3bp.CR3BPSystem,
    prioritizer_with_scorers: ftp.FiveTierPrioritizer,
) -> None:
    """A candidate with no admissible neighbor falls back to the parallel pipeline.

    We use a far-from-EM tiny-mu Sun-Earth-like system so the catalogue's
    EM-only CR3BP rows are all far in the (mu, C, T) tuple metric. With a
    tight ``max_distance``, no surrogate qualifies; tiers 0-3 stay None,
    tier 4 still populates from the parallel pipeline. The point is to
    exercise the FALLBACK PATH, not to propose a physically-meaningful run
    (the FTLE field is computed in the EM system; the candidate state is
    just a smoke 6-vector at a position in the field's range).
    """
    # Tiny-mu system (Sun-Earth-ish: mu=3e-6 vs EM mu=1.215e-2).
    far_system = cr3bp.CR3BPSystem(
        mu=3.0e-6,
        primary="Sun",
        secondary="Earth",
        l_km=1.496e8,
        t_s=5.022e6,
    )
    # State + period at a generic point; only the (mu, C, T) tuple is used for
    # neighbor lookup. Tier 4 uses (state[0], state[1]) sampled in the EM
    # FTLE field; we put the candidate at (0.7, 0.1) (interior, on-grid).
    state0 = np.array([0.7, 0.1, 0.0, 0.0, 0.5, 0.0], dtype=np.float64)
    period_nd = 1.0  # far from any EM cycler period
    jacobi = cr3bp.jacobi_constant(state0, far_system.mu)

    score = sop.score_single_orbit(
        candidate_state0=state0,
        candidate_period_nondim=period_nd,
        candidate_system=far_system,
        candidate_id="fallback-test",
        candidate_jacobi=jacobi,
        prioritizer=prioritizer_with_scorers,
        surrogate_neighbor_max_distance=0.1,
    )
    # Fallback: no surrogate id (catalogue is all EM mu=1.215e-2; relative
    # mu distance to mu=3e-6 is ~4000, well above 0.1).
    assert score.surrogate_pair_neighbor_id is None, (
        f"expected no surrogate for tiny-mu candidate; got {score.surrogate_pair_neighbor_id!r}"
    )
    # Tiers 0-3 cannot populate without a surrogate.
    assert score.tier_0_score is None
    assert score.tier_1_score is None
    assert score.tier_2_score is None
    assert score.tier_3_score is None
    # Tier 4 still runs (uses the EM FTLE field; valid on-grid sample).
    assert score.tier_4_score is not None
    # Notes record the fallback reason.
    assert "no catalogue row within max_distance" in score.notes
