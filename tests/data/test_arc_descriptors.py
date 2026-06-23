"""Schema v4.1 arc-type descriptor tests (spec §16.7.7, 2026-06-03).

Tests for the `free_return_arcs` field -- Russell's Earth-to-Earth free-return
arc decomposition, distinct from encounter-segment decomposition (§16.6.2).

Source authority
----------------
All golden assertions trace to Russell 2004 dissertation Tables 4.9-4.13
(p.127-134) plus the descriptor notation explained at p.126:
- "the first number in the parenthesis is the time of flight in years"
- "the number following the colon in the full-rev strings represent the number
  of revolutions by the spacecraft" (M:N resonant orbit)
- uppercase letter = designated transit leg (transit times + Mars v-inf computed
  from this leg); lowercase = secondary leg
- letters: g/G = generic, h/H = half-rev, f/F = full-rev
"""

from __future__ import annotations

import re
from typing import Any

import pytest

from cyclerfinder.data.catalog import CatalogueEntry, _entry_from_yaml, load_catalog

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_GENERIC_HALF_RE = re.compile(
    r"""
    (?P<letter>[gGhH])         # arc-type letter (generic or half-rev)
    \(
        (?P<tof>[0-9.]+)       # TOF in years (Russell p.127: "first number is TOF in years")
        ,
        (?P<psi>[^,)]+)        # ψ angle in degrees (referencing angle on v∞ sphere)
        ,
        (?P<branch>[^)]+)      # branch code: Ll, Ls, U, L, etc.
    \)
    """,
    re.VERBOSE,
)

_FULL_REV_RE = re.compile(
    r"""
    (?P<letter>[fF])           # full-rev arc
    \(
        (?P<resonance>\d+:\d+) # M:N resonant orbit (e.g. 3:2, 2:1, 1:1)
        ,
        (?P<p2>[^,)]+)         # longitude degrees
        ,
        (?P<p3>[^)]+)          # branch or second longitude
    \)
    """,
    re.VERBOSE,
)


def _parse_descriptor(token: str) -> dict[str, Any]:
    """Parse a single Russell descriptor token into a dict.

    Returns a dict with keys: arc_type, tof_years, resonance (or None),
    raw_descriptor.

    Per Russell p.127:
    - g/h arcs: format is letter(tof_years, psi_deg, branch_code)
    - f/F arcs: format is letter(M:N, longitude_deg, branch_deg)
      The TOF is NOT explicitly listed for full-rev arcs — it is determined
      by the M:N resonance condition and is therefore None here.
    """
    raw = token.strip()
    letter = raw[0].lower()

    if letter in ("g", "h"):
        m = _GENERIC_HALF_RE.match(raw)
        if m is None:
            raise ValueError(f"Cannot parse g/h descriptor: {raw!r}")
        arc_type = "generic" if letter == "g" else "half-rev"
        return {
            "arc_type": arc_type,
            "resonance": None,
            "tof_years": float(m.group("tof")),
            "raw_descriptor": raw,
        }
    elif letter == "f":
        m = _FULL_REV_RE.match(raw)
        if m is None:
            raise ValueError(f"Cannot parse f/F descriptor: {raw!r}")
        return {
            "arc_type": "full-rev",
            "resonance": m.group("resonance"),
            "tof_years": None,  # not given explicitly; determined by M:N resonance
            "raw_descriptor": raw,
        }
    else:
        raise ValueError(f"Unknown arc type letter: {raw[0]!r}")


# ---------------------------------------------------------------------------
# Unit tests for descriptor parser
# ---------------------------------------------------------------------------


def test_parse_generic_descriptor_g() -> None:
    """g(1.4612,526.02,Ll) → generic, no resonance, tof=1.4612 yr."""
    result = _parse_descriptor("g(1.4612,526.02,Ll)")
    assert result["arc_type"] == "generic"
    assert result["resonance"] is None
    assert result["tof_years"] == pytest.approx(1.4612)
    assert result["raw_descriptor"] == "g(1.4612,526.02,Ll)"


def test_parse_generic_descriptor_uppercase_g() -> None:
    """G(2.8096,651.46,U) -> generic (uppercase = designated transit leg)."""
    result = _parse_descriptor("G(2.8096,651.46,U)")
    assert result["arc_type"] == "generic"
    assert result["resonance"] is None
    assert result["tof_years"] == pytest.approx(2.8096)
    assert result["raw_descriptor"] == "G(2.8096,651.46,U)"


def test_parse_full_rev_descriptor_f() -> None:
    """f(3:2,82.487,118.851) → full-rev, resonance=3:2, tof_years=None.

    For full-rev arcs, the first argument is the M:N resonance ratio, NOT a
    TOF in years. Per Russell p.127: "the number following the colon in the
    full-rev strings represent the number of revolutions by the spacecraft,
    thus the full-rev return is an M:N resonant orbit." The TOF is determined
    by the resonance condition and Earth's orbital period.
    """
    result = _parse_descriptor("f(3:2,82.487,118.851)")
    assert result["arc_type"] == "full-rev"
    assert result["resonance"] == "3:2"
    assert result["tof_years"] is None  # not given explicitly in full-rev descriptors


def test_parse_full_rev_descriptor_uppercase_f() -> None:
    """F(3:2,82.487,180.000) -> full-rev, resonance=3:2, uppercase=designated transit."""
    result = _parse_descriptor("F(3:2,82.487,180.000)")
    assert result["arc_type"] == "full-rev"
    assert result["resonance"] == "3:2"
    assert result["raw_descriptor"] == "F(3:2,82.487,180.000)"


def test_parse_full_rev_1_1_resonance() -> None:
    """f(1:1,74.468,-180.000) → full-rev, resonance=1:1."""
    result = _parse_descriptor("f(1:1,74.468,-180.000)")
    assert result["arc_type"] == "full-rev"
    assert result["resonance"] == "1:1"


def test_parse_full_rev_2_1_resonance() -> None:
    """F(2:1,85.196,0.000) → full-rev, resonance=2:1."""
    result = _parse_descriptor("F(2:1,85.196,0.000)")
    assert result["arc_type"] == "full-rev"
    assert result["resonance"] == "2:1"


# ---------------------------------------------------------------------------
# Loader round-trip: CatalogueEntry reads free_return_arcs
# ---------------------------------------------------------------------------


def _make_row_with_arcs(arcs: list[dict[str, Any]]) -> dict[str, Any]:
    """Minimal multi-arc YAML row with a free_return_arcs block."""
    return {
        "id": "test-arc",
        "bodies": ["E", "M"],
        "sequence_canonical": "E-E-M-M",
        "cycler_class": "multi-arc",
        "orbit_elements": {"a_au": None, "e": None},
        "free_return_arcs": arcs,
    }


def test_loader_reads_free_return_arcs_field() -> None:
    """CatalogueEntry.free_return_arcs is populated when present in YAML."""
    arcs = [
        {"arc_type": "generic", "resonance": None, "raw_descriptor": "g(1.4612,526.02,Ll)"},
        {"arc_type": "generic", "resonance": None, "raw_descriptor": "G(2.8096,651.46,U)"},
    ]
    row = _make_row_with_arcs(arcs)
    entry = _entry_from_yaml(row)
    assert entry.free_return_arcs is not None
    assert len(entry.free_return_arcs) == 2
    assert entry.free_return_arcs[0]["arc_type"] == "generic"
    assert entry.free_return_arcs[1]["arc_type"] == "generic"
    assert entry.free_return_arcs[1]["raw_descriptor"] == "G(2.8096,651.46,U)"


def test_loader_free_return_arcs_none_when_absent() -> None:
    """CatalogueEntry.free_return_arcs is None when field not in YAML."""
    row: dict[str, Any] = {
        "id": "test-no-arcs",
        "bodies": ["E", "M"],
        "sequence_canonical": "E-M",
        "cycler_class": "multi-arc",
        "orbit_elements": {"a_au": None, "e": None},
    }
    entry = _entry_from_yaml(row)
    assert entry.free_return_arcs is None


def test_loader_free_return_arcs_with_full_rev_resonance() -> None:
    """full-rev arc with resonance is round-tripped correctly."""
    arcs: list[dict[str, Any]] = [
        {"arc_type": "generic", "resonance": None, "raw_descriptor": "g(1.4646,527.25,Ll)"},
        {"arc_type": "generic", "resonance": None, "raw_descriptor": "g(1.9416,338.97,U)"},
        {"arc_type": "full-rev", "resonance": "3:2", "raw_descriptor": "F(3:2,82.487,180.000)"},
    ]
    row = _make_row_with_arcs(arcs)
    entry = _entry_from_yaml(row)
    assert entry.free_return_arcs is not None
    assert len(entry.free_return_arcs) == 3
    assert entry.free_return_arcs[2]["arc_type"] == "full-rev"
    assert entry.free_return_arcs[2]["resonance"] == "3:2"


# ---------------------------------------------------------------------------
# Golden tests for backfilled catalogue entries
# ---------------------------------------------------------------------------
# For every multi-arc row whose notes contain an explicit Russell leg descriptor,
# assert the parsed arc_type sequence + resonance match the raw descriptor string.
# This is the "descriptor IS the source" contract.
#
# Source: Russell 2004 dissertation Tables 4.9-4.13 (pp.127-134), descriptor
# notation explained p.126.


def _arc_seq(entry: CatalogueEntry) -> list[str]:
    """Extract arc_type sequence from free_return_arcs."""
    if entry.free_return_arcs is None:
        return []
    return [a["arc_type"] for a in entry.free_return_arcs]


def _resonances(entry: CatalogueEntry) -> list[str | None]:
    """Extract resonance list from free_return_arcs."""
    if entry.free_return_arcs is None:
        return []
    return [a.get("resonance") for a in entry.free_return_arcs]


def _raw_descriptors(entry: CatalogueEntry) -> list[str | None]:
    """Extract raw_descriptor list from free_return_arcs."""
    if entry.free_return_arcs is None:
        return []
    return [a.get("raw_descriptor") for a in entry.free_return_arcs]


@pytest.fixture(scope="module")
def catalog() -> Any:
    return load_catalog()


def test_4991gg2_arc_sequence(catalog: Any) -> None:
    """russell-ch4-4.991gG2: g(1.4612,526.02,Ll) G(2.8096,651.46,U) → [generic, generic].

    Source: Russell 2004 dissertation Table 4.9 p.127 row 1 (the S1L1 cycler,
    footnote 'a': documented in Ref. 15 = McConaghy 2002 AIAA 2002-4420).
    Descriptor confirmed verbatim on p.127.
    """
    entry = catalog.by_id["russell-ch4-4.991gG2"]
    assert _arc_seq(entry) == ["generic", "generic"]
    assert _resonances(entry) == [None, None]
    # Verify raw descriptors are preserved verbatim
    raw = _raw_descriptors(entry)
    assert raw[0] == "g(1.4612,526.02,Ll)"
    assert raw[1] == "G(2.8096,651.46,U)"


def test_5_30ggf3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-5.30ggF3: g g F(3:2,...) → [generic, generic, full-rev] resonance=3:2.

    Source: Russell 2004 Table 4.12 p.132 (cycler 5.30ggF, bolded entry).
    Descriptor: g(1.4646,527.25,Ll) g(1.9416,338.97,U) F(3:2,82.487,180.000).
    """
    entry = catalog.by_id["russell-ch4-5.30ggF3"]
    assert _arc_seq(entry) == ["generic", "generic", "full-rev"]
    assert _resonances(entry) == [None, None, "3:2"]
    raw = _raw_descriptors(entry)
    assert raw[2] == "F(3:2,82.487,180.000)"


def test_8049ggf2_arc_sequence(catalog: Any) -> None:
    """russell-ch4-8.049gGf2: g G f(1:1,...) → [generic, generic, full-rev] resonance=1:1.

    Source: Russell 2004 Table 4.9 p.127 row 2 (new cycler, not previously documented).
    Descriptor: g(1.4951,538.24,Ll) G(1.7757,279.24,U) f(1:1,74.468,-180.000).
    """
    entry = catalog.by_id["russell-ch4-8.049gGf2"]
    assert _arc_seq(entry) == ["generic", "generic", "full-rev"]
    assert _resonances(entry) == [None, None, "1:1"]
    raw = _raw_descriptors(entry)
    assert raw[0] == "g(1.4951,538.24,Ll)"
    assert raw[2] == "f(1:1,74.468,-180.000)"


def test_5_75ggf3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-5.75ggF3: g g F(2:1,...) → [generic, generic, full-rev] resonance=2:1.

    Source: Russell 2004 Table 4.12 p.132 (bolded cycler, short transfer time).
    Descriptor: g(1.8987,323.54,U) g(2.5074,902.67,L) F(2:1,85.196,0.000).
    """
    entry = catalog.by_id["russell-ch4-5.75ggF3"]
    assert _arc_seq(entry) == ["generic", "generic", "full-rev"]
    assert _resonances(entry) == [None, None, "2:1"]


def test_3_64ggg3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-3.64gGg3: g G f(1:1,...) → [generic, generic, full-rev] resonance=1:1.

    Source: Russell 2004 Table 4.10 p.128. The descriptor in the notes reads
    "g(2.4845,894.42,Ll) G(2.9217,691.79,U) f(1:1,82.995,-180.000)".
    The entry shorthand name says 'gGg3' but the explicit descriptor has a
    1:1 full-rev orbit as the third leg. The explicit table descriptor takes
    precedence over the shorthand label.
    """
    entry = catalog.by_id["russell-ch4-3.64gGg3"]
    assert _arc_seq(entry) == ["generic", "generic", "full-rev"]
    assert _resonances(entry) == [None, None, "1:1"]


def test_3_78gg3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-3.78Gg3: G g → [generic, generic].

    Source: Russell 2004 Table 4.10 p.128 (bolded entry, favorable turning angles).
    Descriptor: G(2.9043,685.56,U) g(3.5018,1260.65,L).
    """
    entry = catalog.by_id["russell-ch4-3.78Gg3"]
    assert _arc_seq(entry) == ["generic", "generic"]
    assert _resonances(entry) == [None, None]


def test_9_94gg3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-9.94Gg3: G g → [generic, generic].

    Source: Russell 2004 Table 4.11 (extremely short transit 82d).
    Descriptor: G(1.7025,252.9,U) g(4.7037,2053.31,Ls).
    """
    entry = catalog.by_id["russell-ch4-9.94Gg3"]
    assert _arc_seq(entry) == ["generic", "generic"]
    assert _resonances(entry) == [None, None]


def test_3_66gff3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-3.66gfF3: g f(1:1,...) F(3:2,...) → [generic, full-rev, full-rev].

    Source: Russell 2004 Table 4.12 p.132 (bolded entry 3.66gfF).
    Descriptor: g(2.4062,866.21,Ls) f(1:1,82.955,87.388) F(3:2,87.27,0.000).
    """
    entry = catalog.by_id["russell-ch4-3.66gfF3"]
    assert _arc_seq(entry) == ["generic", "full-rev", "full-rev"]
    assert _resonances(entry) == [None, "1:1", "3:2"]


def test_6_44gg3_arc_sequence(catalog: Any) -> None:
    """russell-ch4-6.44Gg3: g G → [generic, generic].

    Source: Russell 2004 Table 4.13 p.134 (near-ballistic, very long ToF).
    Descriptor: g(2.087,1111.33,L) G(4.3191,1194.88,L).
    """
    entry = catalog.by_id["russell-ch4-6.44Gg3"]
    assert _arc_seq(entry) == ["generic", "generic"]
    assert _resonances(entry) == [None, None]


def test_9_353gg2_arc_sequence(catalog: Any) -> None:
    """russell-ch4-9.353Gg2: G g → [generic, generic].

    Source: Russell 2004 Table 4.9 p.127 row 4 (high-energy 2-synodic).
    Descriptor: G(1.7238,260.58,U) g(2.5469,916.9,L).
    """
    entry = catalog.by_id["russell-ch4-9.353Gg2"]
    assert _arc_seq(entry) == ["generic", "generic"]
    assert _resonances(entry) == [None, None]


# ---------------------------------------------------------------------------
# Entries without explicit descriptors must have null free_return_arcs
# ---------------------------------------------------------------------------


def test_entries_without_descriptors_have_null_arcs(catalog: Any) -> None:
    """Multi-arc entries with no known Russell descriptor have null free_return_arcs.

    This is the 'gap, not error' contract from spec §16.7.7.
    """
    # These entries have no explicit descriptor in the current notes.
    # (russell-ocampo-4.3.1-5 removed #388 2026-06-23: its McConaghy-2005 Table 2
    # per-arc descriptor was ingested — it is now descriptor-bearing.)
    no_descriptor_ids = [
        "russell-ocampo-2.1.1+2-case2",
        "russell-ocampo-2.3.1+1-case3",
        "russell-ch4-5.66Gfh3",
    ]
    for eid in no_descriptor_ids:
        entry = catalog.by_id[eid]
        assert entry.free_return_arcs is None, (
            f"{eid}: expected null free_return_arcs (no source descriptor), "
            f"got {entry.free_return_arcs!r}"
        )


# ---------------------------------------------------------------------------
# mcconaghy-2006-em-k2 also carries descriptor from Russell Table 4.9
# ---------------------------------------------------------------------------


def test_mcconaghy_2006_arc_sequence(catalog: Any) -> None:
    """mcconaghy-2006-em-k2: same descriptor as 4.991gG2 → [generic, generic].

    Source: Russell 2004 Table 4.9 row 1 notes in the mcconaghy entry.
    The descriptor "g(1.4612,526.02,Ll) G(2.8096,651.46,U)" is quoted verbatim.
    """
    entry = catalog.by_id["mcconaghy-2006-em-k2"]
    assert _arc_seq(entry) == ["generic", "generic"]
    assert _resonances(entry) == [None, None]


# ---------------------------------------------------------------------------
# Raw-descriptor fidelity: the stored string must reproduce the original token
# ---------------------------------------------------------------------------


def test_raw_descriptor_faithfully_preserved(catalog: Any) -> None:
    """Every stored raw_descriptor must parse back to the same arc_type + resonance.

    This validates the 'raw_descriptor IS the source' contract: if the
    stored string is correct, the arc_type can be verified independently
    by parsing the descriptor letter.
    """
    # Entries with backfilled descriptors
    backfilled_ids = [
        "russell-ch4-4.991gG2",
        "russell-ch4-5.30ggF3",
        "russell-ch4-5.75ggF3",
        "russell-ch4-8.049gGf2",
        "russell-ch4-3.78Gg3",
    ]
    for eid in backfilled_ids:
        entry = catalog.by_id[eid]
        assert entry.free_return_arcs is not None, f"{eid}: expected free_return_arcs"
        for arc in entry.free_return_arcs:
            raw = arc.get("raw_descriptor")
            assert raw is not None, f"{eid}: arc missing raw_descriptor"
            # The stored arc_type must be consistent with the descriptor letter
            letter = raw[0].lower()
            expected_type = {"g": "generic", "h": "half-rev", "f": "full-rev"}[letter]
            assert arc["arc_type"] == expected_type, (
                f"{eid}: arc_type={arc['arc_type']!r} inconsistent "
                f"with descriptor letter {raw[0]!r}"
            )
