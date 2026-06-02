"""Catalogue-driven rediscovery tests.

For every constructible entry in the seed catalogue (per
``tests/_catalogue_loader.py``), build the corresponding M4 cell from
its published ``(bodies, sequence_canonical, period_k)`` and run the M5
optimiser. Assert the rediscovered V∞ signature matches the entry's
published targets within tolerance.

Tolerance rationale
-------------------
The M5 gate (``tests/search/test_optimize.py``
``test_2syn_em_rediscovers_5_65_kms_earth``) uses ±0.2 km/s on the
canonical anchor (5.65 km/s Earth, 3.05 km/s Mars on the 2-syn E-M
cycler). Catalogue entries come from various circular-coplanar
references; absolute values shift by ~0.1-0.3 km/s between sources
(e.g. Russell 2004 tabulates 4.99 / 5.10 km/s for the McConaghy "Notable
Two-Synodic"; McConaghy's own 2006 abstract gives 4.7 / 5.0 km/s for the
SAME cycler). We use ±0.3 km/s to absorb these inter-source roundings
and the optimiser's finite stochasticity without losing rediscovery
fidelity.

What the test asserts per entry
-------------------------------
* ``result.constraints_satisfied`` — the optimiser found a
  feasible interior solution (V∞ ≤ cap, r_p ≥ floor at every
  encounter).
* For every body listed in the entry, the rediscovered
  ``max(||vinf_in||, ||vinf_out||)`` at any encounter of that body
  matches the entry's published ``vinf_kms`` within ``VINF_TOL_KMS``.

The ``vinf_cap`` passed to the optimiser per entry is the entry's
maximum published V∞ plus ``VINF_CAP_HEADROOM_KMS`` so the optimiser
has room to refine without tripping the cap at the published value
itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import (
    OptimisationResult,
    _target_period_sec,
    interior_epochs_from_leg_tofs,
    optimise_cell_idealized,
)
from cyclerfinder.search.sequence import Cell
from tests._catalogue_loader import (
    CatalogueEntry,
    ExclusionReason,
    classify_catalogue,
    load_constructible_entries,
)

# ---------------------------------------------------------------------------
# Tolerances and harness constants
# ---------------------------------------------------------------------------

VINF_TOL_KMS: float = 0.3
"""±0.3 km/s per encounter. See module docstring rationale."""

VINF_CAP_HEADROOM_KMS: float = 2.5
"""Headroom above the entry's max published V∞ for the optimiser's
``vinf_cap``. Without headroom the cap pins the constraint at the target
value and the optimiser has no room to refine. The 2.5 km/s value
matches the spacing used by the M5 gate tests (the Aldrin gate uses
``vinf_cap=12.0`` against a 9.7 km/s anchor — 2.3 km/s of headroom)."""

_OPTIMISER_N_STARTS: int = 5
"""Matches the M5 spec §13.4 default; same as the M5 gate test."""

_OPTIMISER_SEED: int = 0
"""Reproducibility seed shared across the suite."""

_OPTIMISER_USE_DE: bool = True
"""Enable the ``differential_evolution`` global pass per-entry.

Most catalogue entries (e.g. the Aldrin family at V∞ ≈ 9.7 km/s at
Mars) live in a basin the multi-start grid alone cannot reach from the
equispaced free-return seed; DE is required to drive SLSQP into the
correct family. Matches the M5 gate test's
``test_aldrin_regression_anchor`` setup. Per spec §13.4 / plan §4.5
the DE budget is bounded (~400 evaluations / cell) but the per-cell
wall-clock is still 5-10 minutes on CI's 2-core runners — the
parametrised test is therefore tagged as a slow-suite gate in the
project's pytest config.
"""


# ---------------------------------------------------------------------------
# Expected-skip registry — every entry here MUST have a documented reason
# ---------------------------------------------------------------------------

EXPECTED_SKIPS: dict[str, str] = {
    "aldrin-classic-em-k1-outbound": (
        "MODEL LIMITATION, not a seeding fix (re-diagnosed 2026-06-01, #52). "
        "Aldrin's 146-day E->M transit is only cheap with Mars's real "
        "eccentricity; in the circular-coplanar idealisation a 146-day leg "
        "is near-hyperbolic. Verified empirically: building the cycler at the "
        "published interior epoch (Mars at t=146 d) gives V_inf~21.5 km/s, and "
        "a 1-D sweep over the entire interior-epoch range (0.02T..0.98T) never "
        "produces a closing sub-cap solution — the lowest max-V_inf is 6.09 "
        "km/s at the T/2 midpoint with closure residual 21.8 km/s. So the "
        "catalogue-seeded warm-start (#52, now implemented and unit-tested) "
        "cannot close this gate: there is no feasible Aldrin basin in the "
        "circular-coplanar model regardless of seed. Rediscovering Aldrin's "
        "9.74 km/s signature requires the real-ephemeris (M6b) optimiser. The "
        "cell IS strict-ballistic per Russell 2004 Table 3.4 (1.0.1.-1); the "
        "v1 idealised model just cannot host its asymmetric real-ephemeris "
        "geometry."
    ),
    "aldrin-classic-em-k1-inbound": (
        "Symmetric counterpart of aldrin-classic-em-k1-outbound; same "
        "circular-coplanar model limitation. See that entry's reason."
    ),
}
"""Entries whose published V∞ signature is known not to match the
circular-coplanar Lambert reconstruction the v1 optimiser computes
(e.g. published values are ephemeris-optimised or come from a
non-Keplerian model). Document the reason per entry; reference
``data/OUTSTANDING.md`` if a pattern emerges.

Empty by default — the loader's constructibility filter is intended to
exclude exactly the entries that can't round-trip, so this dict should
stay short.
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_cell_from_entry(entry: CatalogueEntry) -> Cell:
    """Construct the v1-optimiser cell template for ``entry``.

    The loader has already filtered to entries whose
    ``sequence_canonical`` is a 2-body alternation (e.g. ``"E-M"``).
    The matching cell is the minimal closed 3-encounter cycle
    ``(body_a, body_b, body_a)`` with ``period_k`` from the entry and
    direct (single-branch, 0-rev) legs — the only ballistic regime the
    v1 Lambert solver supports.
    """
    body_a, body_b = entry.bodies
    return Cell(
        bodies=(body_a, body_b),
        sequence=(body_a, body_b, body_a),
        period_k=entry.period_k,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )


def _warm_starts_for_entry(entry: CatalogueEntry, cell: Cell) -> list[tuple[float, ...]] | None:
    """Build a catalogue-seeded warm start from the entry's leg ToFs (#52).

    Returns ``None`` unless the entry **fully tabulates** its legs — i.e.
    it carries exactly ``len(cell.sequence) - 1`` leg ToFs, one per
    transfer. Most seed-catalogue entries are *under-tabulated* (only the
    forward transit is recorded, not the return loop the closed sequence
    implies); for those the cumulative mapping cannot produce a
    correctly-sized interior-epoch vector, so we fall back to the
    grid + DE search rather than feed the optimiser a malformed seed.

    The warm start is opportunistic: it only ever *adds* a polished
    candidate to the optimiser's ranking, never removes the grid/DE
    coverage, so it cannot regress an entry that already rediscovers.
    """
    n_legs = len(cell.sequence) - 1
    if len(entry.leg_tofs_days) != n_legs:
        return None
    warm = interior_epochs_from_leg_tofs(entry.leg_tofs_days, _target_period_sec(cell))
    return [warm]


def _vinf_magnitudes_by_body(result: OptimisationResult) -> dict[str, float]:
    """Map ``body_code -> max(||vinf_in||, ||vinf_out||)`` across encounters.

    Same helper shape as ``tests/search/test_optimize.py``. ``max``
    rather than ``mean`` handles cells where the same body appears at
    multiple encounters — the largest V∞ is the binding magnitude for
    rediscovery against published values.
    """
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:
        m = max(
            float(np.linalg.norm(enc.vinf_in)),
            float(np.linalg.norm(enc.vinf_out)),
        )
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def _entry_target_by_body(entry: CatalogueEntry) -> dict[str, float]:
    """Map ``body_code -> published vinf_kms`` from the entry.

    Loader guarantees ``len(bodies) == len(vinf_targets_kms)``. If a
    body appears twice in ``bodies`` (not possible under the loader's
    2-body alternation filter, but defensive), the later occurrence
    wins — and the test will fail loudly so the loader can be tightened.
    """
    return dict(zip(entry.bodies, entry.vinf_targets_kms, strict=True))


# ---------------------------------------------------------------------------
# Sanity tests on the loader (run cheaply; protect against silent
# filter regressions as the catalogue grows)
# ---------------------------------------------------------------------------


def test_loader_returns_at_least_one_entry() -> None:
    """Sanity: the filter shouldn't return an empty list as the catalogue
    grows. Aldrin classic outbound is always present (catalogue entry 1)."""
    entries = load_constructible_entries()
    assert len(entries) > 0
    assert any(e.id == "aldrin-classic-em-k1-outbound" for e in entries)


def test_loader_filter_excludes_family_seeds() -> None:
    """Family-seed entries (null orbital data, citation-only) must be skipped."""
    entries = load_constructible_entries()
    family_seed_ids = {
        "jones-2017-vem-triple-family",
        "wittal-2022-em-cycler-family",
        "russell-strange-2009-jovian-multimoon-family",
        "russell-strange-2009-saturnian-multimoon-family",
        "hernandez-2017-jovian-ieg-triple-family",
        "mcconaghy-2005-em-snlm-broadclass-family",
    }
    found = {e.id for e in entries} & family_seed_ids
    assert not found, f"family-seed entries leaked through filter: {found}"


def test_loader_filter_excludes_non_heliocentric() -> None:
    """Lunar / Jovian entries (``primary != "Sun"``) must be skipped —
    the v1 optimiser is heliocentric."""
    entries = load_constructible_entries()
    # All current entries have ``primary`` absent or ``"Sun"``; this
    # test guards the filter against a future regression that lets a
    # non-heliocentric entry slip through.
    import yaml  # type: ignore[import-untyped]

    from tests._catalogue_loader import CATALOGUE_PATH

    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    non_helio_ids = {row["id"] for row in raw if row.get("primary") not in (None, "Sun")}
    leaked = {e.id for e in entries} & non_helio_ids
    assert not leaked, f"non-heliocentric entries leaked: {leaked}"


# ---------------------------------------------------------------------------
# Coverage audit — the V0→V1 promotion gauntlet census (task #55)
#
# The gauntlet rediscovers only ``CONSTRUCTIBLE`` entries. The vast
# majority of the catalogue is *not* constructible by the v1 optimiser,
# and that is fine — but it must never be excluded *silently*. These
# tests assert that every catalogue row is classified into exactly one
# documented :class:`ExclusionReason`, and freeze the distribution as a
# ratchet so any change to the catalogue's testability profile surfaces
# as a reviewed diff rather than a quiet drop in coverage.
# ---------------------------------------------------------------------------

EXPECTED_COVERAGE: dict[ExclusionReason, int] = {
    ExclusionReason.MULTI_ENCOUNTER_SEQUENCE: 202,
    ExclusionReason.NON_HELIOCENTRIC: 6,
    ExclusionReason.MISSING_VINF: 6,
    ExclusionReason.CONSTRUCTIBLE: 2,
    ExclusionReason.NOT_TWO_BODY: 2,
    ExclusionReason.MISSING_PERIOD: 1,
}
"""Frozen census of how the 219-row catalogue distributes across
exclusion reasons (as of 2026-06-02). This is a *ratchet*: when the
catalogue changes, this dict must be updated in the same commit, which
forces a conscious review of whether the change moved entries into or
out of the v1 gauntlet's reach.

Reasons absent from this dict are expected to have a count of zero (e.g.
``NON_BALLISTIC`` and ``MISSING_LEG_TOFS`` — no such rows today).
"""


def _census_counts() -> dict[ExclusionReason, int]:
    counts: dict[ExclusionReason, int] = {}
    for _id, reason in classify_catalogue():
        counts[reason] = counts.get(reason, 0) + 1
    return counts


def test_census_classifies_every_entry_exactly_once() -> None:
    """No catalogue row may be dropped silently: the census must cover the
    full catalogue, with one reason per id and no duplicate ids."""
    import yaml

    from tests._catalogue_loader import CATALOGUE_PATH

    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    all_ids = [row["id"] for row in raw]

    census = classify_catalogue()
    census_ids = [cid for cid, _ in census]

    assert len(census_ids) == len(all_ids), (
        f"census size {len(census_ids)} != catalogue size {len(all_ids)}"
    )
    assert len(set(census_ids)) == len(census_ids), "duplicate ids in census"
    assert set(census_ids) == set(all_ids), (
        f"census/catalogue id mismatch: "
        f"missing={set(all_ids) - set(census_ids)}, "
        f"extra={set(census_ids) - set(all_ids)}"
    )


def test_census_constructible_set_agrees_with_loader() -> None:
    """The census and the gauntlet's own loader must name the *same*
    constructible entries — one classifier, no divergence."""
    census_constructible = {
        cid for cid, reason in classify_catalogue() if reason is ExclusionReason.CONSTRUCTIBLE
    }
    loader_constructible = {e.id for e in load_constructible_entries()}
    assert census_constructible == loader_constructible, (
        f"census vs loader disagree on constructible set: "
        f"census_only={census_constructible - loader_constructible}, "
        f"loader_only={loader_constructible - census_constructible}"
    )


def test_census_breakdown_matches_frozen_ratchet() -> None:
    """The exclusion-reason distribution must match ``EXPECTED_COVERAGE``.

    A mismatch means the catalogue changed in a way that altered which
    entries the v1 gauntlet can reach. Update ``EXPECTED_COVERAGE`` in the
    same commit, having confirmed the shift is intended.
    """
    counts = _census_counts()
    # Normalise: reasons absent from EXPECTED_COVERAGE are expected to be 0.
    expected = {reason: EXPECTED_COVERAGE.get(reason, 0) for reason in ExclusionReason}
    actual = {reason: counts.get(reason, 0) for reason in ExclusionReason}
    assert actual == expected, (
        "catalogue testability profile changed.\n"
        f"  expected: { {r.value: n for r, n in expected.items() if n} }\n"
        f"  actual:   { {r.value: n for r, n in actual.items() if n} }\n"
        "If intended, update EXPECTED_COVERAGE in this commit."
    )


def test_expected_skips_are_constructible() -> None:
    """Every ``EXPECTED_SKIPS`` key must still be a CONSTRUCTIBLE entry.

    A skip documented for an entry that is no longer constructible (renamed,
    re-classified, or removed) is dead weight that hides a coverage gap —
    fail loudly so it is cleaned up or re-justified.
    """
    constructible = {
        cid for cid, reason in classify_catalogue() if reason is ExclusionReason.CONSTRUCTIBLE
    }
    stale = set(EXPECTED_SKIPS) - constructible
    assert not stale, (
        f"EXPECTED_SKIPS names entries that are not constructible (stale skips): {stale}"
    )


# ---------------------------------------------------------------------------
# Catalogue-driven rediscovery — parametrised over every constructible
# entry returned by the loader
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize(
    "entry",
    load_constructible_entries(),
    ids=lambda e: e.id,
)
def test_catalogue_entry_rediscovers(entry: CatalogueEntry) -> None:
    """Optimiser should rediscover the published V∞ signature.

    Per entry:

    1. Build the M4 cell ``(body_a, body_b, body_a)`` with the entry's
       ``period_k`` and direct (single-branch, 0-rev) legs.
    2. Run ``optimise_cell_idealized`` against
       ``Ephemeris("circular")`` with ``vinf_cap`` set to the entry's
       max published V∞ + ``VINF_CAP_HEADROOM_KMS``.
    3. Assert ``constraints_satisfied`` and that every body's
       rediscovered V∞ matches its target within ``VINF_TOL_KMS``.
    """
    if entry.id in EXPECTED_SKIPS:
        pytest.skip(EXPECTED_SKIPS[entry.id])

    cell = _build_cell_from_entry(entry)
    target_by_body = _entry_target_by_body(entry)
    vinf_cap = max(target_by_body.values()) + VINF_CAP_HEADROOM_KMS

    warm_starts = _warm_starts_for_entry(entry, cell)

    eph = Ephemeris(model="circular")
    result = optimise_cell_idealized(
        cell,
        eph,
        vinf_cap=vinf_cap,
        n_starts=_OPTIMISER_N_STARTS,
        seed=_OPTIMISER_SEED,
        use_de=_OPTIMISER_USE_DE,
        warm_starts=warm_starts,
    )

    assert result.constraints_satisfied, (
        f"{entry.id}: hard constraints violated; "
        f"max_vinf={result.best_score.max_vinf_kms}, "
        f"residual={result.closure_residual_kms}, vinf_cap={vinf_cap}"
    )
    assert math.isfinite(result.closure_residual_kms), (
        f"{entry.id}: non-finite closure residual ({result.closure_residual_kms})"
    )

    rediscovered = _vinf_magnitudes_by_body(result)
    for body, target in target_by_body.items():
        assert body in rediscovered, (
            f"{entry.id}: body {body!r} expected in rediscovered "
            f"encounters; got {sorted(rediscovered)}"
        )
        got = rediscovered[body]
        assert abs(got - target) < VINF_TOL_KMS, (
            f"{entry.id}: V_inf at {body} = {got:.3f} km/s vs target "
            f"{target:.3f} km/s (|diff| = {abs(got - target):.3f} > "
            f"tol {VINF_TOL_KMS} km/s, vinf_cap={vinf_cap:.2f})"
        )
