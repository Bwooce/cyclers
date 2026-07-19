"""Tests for `#658`'s epoch-locking pilot: the small epoch-carrying subset
of the catalogue #650's transfer-network machinery (already-built,
already-positive-controlled, unmodified here) can be pointed at.

Verifies the real-catalogue candidate-selection finding directly (this is
the actual result the task depends on, not a fixture-only unit test):
`#654`'s named candidates were "Jones 2017 VEM-triple" and "Aldrin/Byrnes"
as already-epoch-carrying rows. Checked against the live catalogue:

* Jones-lineage: 0 rows qualify (every Jones-authored row is
  `epoch_locked=false`, same as every other `cycler`-class row -- #654's
  framing was wrong for this half).
* Aldrin-lineage: 9 rows qualify, but they are the Rogers/Hughes/Longuski/
  Aldrin 2012/2015 one-shot `precursor_mga` "establishment" trajectories
  (`n_returns=1`), not the steady-state `aldrin-classic-em-k1-*` cycler rows
  (which are themselves `epoch_locked=false`).

These tests use the programmatic selection rule
(`scripts/run_658_epoch_carrying_transfer_check.epoch_carrying_candidates`)
so a future catalogue edit that changes this picture fails loudly here
rather than silently invalidating a stale hardcoded id list.
"""

from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path
from typing import Any

import pytest

from cyclerfinder.data.catalog import CATALOGUE_PATH
from cyclerfinder.data.transfer_network import (
    compute_edge,
    epoch_window_intersection,
    usable_bodies,
)

_SCRIPT_PATH = (
    Path(__file__).resolve().parent.parent.parent
    / "scripts"
    / "run_658_epoch_carrying_transfer_check.py"
)


def _load_script_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "run_658_epoch_carrying_transfer_check", _SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_MODULE = _load_script_module()


def _load_catalogue_rows() -> list[dict[str, Any]]:
    import yaml  # type: ignore[import-untyped]

    with open(CATALOGUE_PATH) as fh:
        data = yaml.safe_load(fh)
    assert isinstance(data, list)
    return data


# ---------------------------------------------------------------------------
# Candidate-selection unit tests (fixture rows, not the live catalogue).
# ---------------------------------------------------------------------------


def test_has_real_epoch_requires_epoch_locked_true_and_real_window() -> None:
    assert _MODULE._has_real_epoch(
        {"epoch_locked": True, "validity_window": {"start": "2020-01-01", "end": "2020-06-01"}}
    )
    assert not _MODULE._has_real_epoch({"epoch_locked": False, "validity_window": None})
    assert not _MODULE._has_real_epoch({"epoch_locked": True, "validity_window": None})
    assert not _MODULE._has_real_epoch(
        {"epoch_locked": True, "validity_window": {"start": None, "end": None}}
    )


def test_epoch_carrying_candidates_filters_by_author_and_real_epoch() -> None:
    rows = [
        {
            "id": "a",
            "first_published": {"authors": ["Aldrin, B."]},
            "epoch_locked": True,
            "validity_window": {"start": "2020-01-01", "end": "2020-06-01"},
        },
        {
            "id": "b",
            "first_published": {"authors": ["Aldrin, B."]},
            "epoch_locked": False,
            "validity_window": None,
        },
        {
            "id": "c",
            "first_published": {"authors": ["Jones, D. R."]},
            "epoch_locked": True,
            "validity_window": {"start": "2020-01-01", "end": "2020-06-01"},
        },
    ]
    got = _MODULE.epoch_carrying_candidates(rows, "Aldrin")
    assert [r["id"] for r in got] == ["a"]
    got_jones = _MODULE.epoch_carrying_candidates(rows, "Jones")
    assert [r["id"] for r in got_jones] == ["c"]


# ---------------------------------------------------------------------------
# Real-catalogue finding -- the actual result this task depends on.
# ---------------------------------------------------------------------------


def test_jones_lineage_has_zero_epoch_carrying_rows_in_real_catalogue() -> None:
    """#654's framing named Jones 2017 VEM-triple as epoch-carrying; verified false."""
    rows = _load_catalogue_rows()
    got = _MODULE.epoch_carrying_candidates(rows, "Jones")
    assert got == [], (
        "expected 0 Jones-lineage epoch-carrying rows (all Jones-authored rows are "
        "epoch_locked=false, cycler-class); if this now finds rows, #654's original "
        "framing may have become true after a catalogue edit -- re-verify by hand "
        "before trusting a nonzero count here"
    )


def test_aldrin_lineage_epoch_carrying_rows_are_exactly_the_expected_nine() -> None:
    """The 9 real epoch-carrying Aldrin-lineage rows are the Rogers 2012/2015
    one-shot precursor_mga establishment trajectories, not the steady-state
    aldrin-classic-em-k1-* cycler rows (which are epoch_locked=false)."""
    rows = _load_catalogue_rows()
    got = _MODULE.epoch_carrying_candidates(rows, "Aldrin")
    got_ids = {r["id"] for r in got}
    expected_ids = {
        "aldrin-4-3-2-establishment",
        "aldrin-3-2-1-establishment",
        "visit-1-4-3-2-establishment",
        "visit-2-4-3-2-establishment",
        "case-1-4-3-2-establishment",
        "case-2-4-3-2-establishment",
        "case-3-4-3-2-establishment",
        "s1l1-4-3-2-establishment",
        "u0l1-3-2-1-establishment",
    }
    assert got_ids == expected_ids
    for r in got:
        assert r["orbit_class"] == "precursor_mga"
        assert r["n_returns"] == 1
        assert r.get("inserts_into") is not None
    # The steady-state Aldrin rows are NOT in this set (they're epoch_locked=false).
    assert "aldrin-classic-em-k1-outbound" not in got_ids
    assert "aldrin-classic-em-k1-inbound" not in got_ids


def test_all_nine_aldrin_lineage_rows_share_an_earth_encounter() -> None:
    rows = {r["id"]: r for r in _load_catalogue_rows()}
    candidates = _MODULE.epoch_carrying_candidates(list(rows.values()), "Aldrin")
    for row in candidates:
        assert "Earth" in usable_bodies(row), row["id"]


def test_pairwise_raw_epoch_window_overlap_matches_expected_five_pairs() -> None:
    """Deterministic (ungated) real-calendar-window overlap check via #650's own
    `epoch_window_intersection`, unmodified -- exactly 5 of the 36 pairs overlap."""
    rows = {r["id"]: r for r in _load_catalogue_rows()}
    candidates = _MODULE.epoch_carrying_candidates(list(rows.values()), "Aldrin")
    ids = sorted(r["id"] for r in candidates)
    assert len(ids) == 9

    overlapping = []
    for id_a, id_b in itertools.combinations(ids, 2):
        if epoch_window_intersection(rows[id_a], rows[id_b]) is not None:
            overlapping.append((id_a, id_b))

    expected = {
        ("aldrin-4-3-2-establishment", "case-3-4-3-2-establishment"),
        ("aldrin-4-3-2-establishment", "s1l1-4-3-2-establishment"),
        ("aldrin-3-2-1-establishment", "visit-2-4-3-2-establishment"),
        ("case-1-4-3-2-establishment", "u0l1-3-2-1-establishment"),
        ("case-3-4-3-2-establishment", "s1l1-4-3-2-establishment"),
    }
    assert set(overlapping) == expected


def test_no_genuinely_cheap_independent_transfer_opportunity_among_candidates() -> None:
    """Headline result: of 36 pairs, exactly 2 are B0_ballistic_compatible
    (delta_vinf<=0.1 km/s at Earth). (1) case-3 <-> s1l1: overlapping windows,
    but a documented near-twin coincidence (same Rogers-2015-table launch
    date, near-identical orbital elements), not an independent-mission
    scheduling discovery. (2) case-2 <-> visit-1: ballistically cheap (both
    near-Hohmann Earth v_inf) but their real calendar windows are 17 YEARS
    apart (2029 vs 2046) -- no overlap at all, reproducing #650's own general
    finding that DeltaV-cheap does not imply a real-date-realizable transfer.
    The 4 real-distinct-mission window-overlap pairs are all B2_moderate
    (dv_hop 0.7-1.7 km/s) -- too expensive to count as cheap under #650's own
    B0/B1 bands. No edge in this candidate set is `cheap_edge=True`."""
    rows = {r["id"]: r for r in _load_catalogue_rows()}
    candidates = _MODULE.epoch_carrying_candidates(list(rows.values()), "Aldrin")
    ids = sorted(r["id"] for r in candidates)

    n_cheap = 0
    b0_pairs = []
    for id_a, id_b in itertools.combinations(ids, 2):
        row_a, row_b = rows[id_a], rows[id_b]
        shared = usable_bodies(row_a) & usable_bodies(row_b)
        for body in shared:
            edge = compute_edge(row_a, row_b, body)
            if edge.cheap_edge:
                n_cheap += 1
            if edge.band == "B0_ballistic_compatible":
                b0_pairs.append((id_a, id_b))

    assert n_cheap == 0
    assert set(b0_pairs) == {
        ("case-3-4-3-2-establishment", "s1l1-4-3-2-establishment"),
        ("case-2-4-3-2-establishment", "visit-1-4-3-2-establishment"),
    }
    # The case-2/visit-1 B0 pair has NO real calendar-window overlap (2029 vs 2046).
    assert (
        epoch_window_intersection(
            rows["case-2-4-3-2-establishment"], rows["visit-1-4-3-2-establishment"]
        )
        is None
    )


def test_case3_s1l1_b0_pair_is_the_known_documented_near_twin_not_a_new_finding() -> None:
    rows = {r["id"]: r for r in _load_catalogue_rows()}
    case3 = rows["case-3-4-3-2-establishment"]
    s1l1 = rows["s1l1-4-3-2-establishment"]
    edge = compute_edge(case3, s1l1, "Earth")
    assert edge.band == "B0_ballistic_compatible"
    assert edge.dv_hop_kms == pytest.approx(0.00981, abs=1e-4)
    # Same source launch date -- documented in the catalogue's own notes as a
    # near-twin coincidence (Rogers 2015 Table 4), not independently derived here.
    assert case3["launch_epoch"] == s1l1["launch_epoch"] == "2022-12-20T00:00:00Z"
    assert "near-twin" not in (case3.get("notes") or "")  # the caveat lives on s1l1's row
    assert "Same LD" in (s1l1.get("notes") or "")
