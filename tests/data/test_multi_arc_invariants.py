"""Phase 3 golden test: Russell multi-arc invariants (AR/TR/transit) in typed invariants{}.

Two-layer ratchet:
1. Note-consistency sweep: for EVERY multi-arc row, if the row's sourced notes
   state an explicit aphelion_ratio or turn_ratio, then ``invariants.aphelion_ratio``
   /``invariants.turn_ratio`` must be present and match it (pytest.approx, atol=0.005).
   This ensures the typed field never drifts from the prose source.

2. Golden spot-checks: a handful of rows whose AR/TR values are explicitly cited in
   the Russell 2004 dissertation tables, listed here as hard-coded expected values
   with the source table in an inline comment. These anchor the typed values to the
   primary literature, not to our own extraction.

Source: Russell 2004 PhD dissertation (UT Austin, hdl.handle.net/2152/1253)
  - Table 3.4: 2-synodic and 3-synodic ballistic cyclers (AR/TR columns explicit)
  - Table 4.9: 2-synodic generic-return cyclers (TR explicit; AR not tabulated)
  - Tables 4.10-4.13: 3-synodic and multi-synodic generic-return cyclers
"""

from __future__ import annotations

import re
from typing import Any

import pytest

# ---------------------------------------------------------------------------
# Helpers to extract what the notes say about AR/TR
# (mirrors the extraction logic in scripts/backfill_invariants.py)
# ---------------------------------------------------------------------------


def _collect_text(obj: Any, depth: int = 0) -> list[str]:
    parts: list[str] = []
    if depth > 10:
        return parts
    if isinstance(obj, str):
        parts.append(obj)
    elif isinstance(obj, dict):
        for v in obj.values():
            parts.extend(_collect_text(v, depth + 1))
    elif isinstance(obj, list):
        for v in obj:
            parts.extend(_collect_text(v, depth + 1))
    return parts


def _notes_ar(row: dict[str, Any]) -> float | None:
    """Return the AR value explicitly stated in the row's sourced text, or None."""
    text = "\n".join(_collect_text(row))
    matches: list[float] = []
    for pat in [
        r"\bAR\s*=\s*([0-9]+\.[0-9]+)",
        r"Aphelion Ratio\s*(?:AR\s*=?\s*|is\s+|=\s*)?([0-9]+\.[0-9]+)",
    ]:
        for m in re.finditer(pat, text):
            matches.append(float(m.group(1)))
    unique: list[float] = []
    for v in matches:
        if not any(abs(v - u) < 0.005 for u in unique):
            unique.append(v)
    if len(unique) == 1:
        return unique[0]
    return None  # absent or ambiguous → skip consistency check


def _notes_tr(row: dict[str, Any]) -> float | None:
    """Return the TR value explicitly stated in the row's sourced text, or None."""
    text = "\n".join(_collect_text(row))
    matches = [float(m.group(1)) for m in re.finditer(r"\bTR\s*=\s*([0-9]+\.[0-9]+)", text)]
    unique: list[float] = []
    for v in matches:
        if not any(abs(v - u) < 0.005 for u in unique):
            unique.append(v)
    if len(unique) == 1:
        return unique[0]
    return None  # absent or ambiguous → skip consistency check


# ---------------------------------------------------------------------------
# 1. Note-consistency sweep
# ---------------------------------------------------------------------------


# Planet-centric (moon-tour Tier-1, task #76) family-seed rows that are
# multi-arc patched-conic moon tours but carry NULL numerics: they record the
# citation + qualitative geometry only (data/README.md "family-seed null-numeric
# records"; Russell-Strange / Hernandez Jovian seeds). They have no sourced
# trajectory.segments and therefore no transit_times / AR / TR to populate —
# fabricating any would violate the golden discipline. They are exempt from the
# completeness invariants below (which assume the full-numeric heliocentric
# Russell Earth-Mars rows). The note-consistency sweeps still apply to them: if
# such a row ever DID state an AR/TR in its notes, the typed-field check fires.
_FAMILY_SEED_NULL_NUMERIC_MULTI_ARC: frozenset[str] = frozenset(
    {
        "hernandez-2017-jovian-ieg-triple-family",
        "russell-strange-2009-jovian-multimoon-family",
        # #216 (2026-06-12): Liang et al. 2024 CGE triple-cycler Member D, the
        # SPICE-ephemeris (JUP365.bsp) member — its per-flyby/per-leg ToF are
        # published as figure traces only (no numeric table; mining note §4.4), so
        # transit_times_days is a genuine data_gap (null). Exempt from the
        # transit-times completeness invariant. (Members A-C DO carry sourced
        # per-leg ToF from Tables 3/5/7 and are NOT exempt.) AR/TR are Earth-Mars
        # free-return concepts inapplicable to a Jovian moon tour -> null for all
        # four Liang members; they state no AR/TR in their notes, so the
        # note-consistency sweeps do not fire.
        "liang-2024-cgcec-ephemeris-2033",
        # Task #408: Russell 2006 cyclers added from Table 5 without trajectory segments.
        "russell-2006-54-3.768ghminus3",
        "russell-2006-55-3.768ghplus3",
        "russell-2006-111-5.219gghminus3",
        "russell-2006-112-5.219gghplus3",
        "russell-2006-117-5.225ggg3",
        "russell-2006-177-5.751ggf3",
        "russell-2006-178-5.751ggf3",
        # #491 (2026-06-30): Russell-Strange 2009 (10 Jovian + 20 Titan-Enceladus)
        # + Lynam-Longuski 2011 (2 IEG) moon cyclers added from the papers' tables
        # (citation + V_inf + qualitative geometry; no sourced per-leg trajectory
        # segments) -> family-seed null-numeric, exempt like the russell-2006 rows.
        "russell-strange-2009-eurgan-131",
        "russell-strange-2009-eurgan-159",
        "russell-strange-2009-gancal-1",
        "russell-strange-2009-gancal-5",
        "russell-strange-2009-ganeur-5",
        "russell-strange-2009-ganeur-43",
        "russell-strange-2009-ganeur-316",
        "russell-strange-2009-ganio-53",
        "russell-strange-2009-ganio-185",
        "russell-strange-2009-ganio-403",
        "russell-strange-2009-titenc-37",
        "russell-strange-2009-titenc-145",
        "russell-strange-2009-titenc-183",
        "russell-strange-2009-titenc-207",
        "russell-strange-2009-titenc-217",
        "russell-strange-2009-titenc-227",
        "russell-strange-2009-titenc-231",
        "russell-strange-2009-titenc-235",
        "russell-strange-2009-titenc-314",
        "russell-strange-2009-titenc-370",
        "russell-strange-2009-titenc-492",
        "russell-strange-2009-titenc-510",
        "russell-strange-2009-titenc-539",
        "russell-strange-2009-titenc-552",
        "russell-strange-2009-titenc-572",
        "russell-strange-2009-titenc-586",
        "russell-strange-2009-titenc-594",
        "russell-strange-2009-titenc-602",
        "russell-strange-2009-titenc-624",
        "russell-strange-2009-titenc-631",
        "lynam-longuski-2011-ieg-single-period",
        "lynam-longuski-2011-gipeipe",
    }
)


def _load_multi_arc_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [r for r in rows if r.get("cycler_class") == "multi-arc"]


def _load_full_numeric_multi_arc_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Multi-arc rows that carry full numerics (excludes family-seed null rows
    AND non-cycler orbit_class rows per schema v4.7 / task #294).

    The completeness invariants below — invariants{} present, AR/TR populated,
    transit_times_days populated — are inherently cycler concepts (aphelion
    ratio and turn ratio are *cycle-level* descriptors; transit times are
    per-flyby ToFs of a *repeating* sequence). The v4.7 scope expansion admits
    non-cycler orbit classes (quasi_cycler / precursor_mga / mga_tour) that
    legitimately carry cycler_class=multi-arc as a *structural* tag (multiple
    distinct heliocentric legs) without being cyclers — the cycler invariants
    do not apply to them. Filter them out here so the test stays focused on
    actual cyclers.
    """
    return [
        r
        for r in _load_multi_arc_rows(rows)
        if r["id"] not in _FAMILY_SEED_NULL_NUMERIC_MULTI_ARC
        and r.get("orbit_class", "cycler") == "cycler"
    ]


def test_note_ar_matches_invariants_field(catalogue_rows: list[dict[str, Any]]) -> None:
    """For every multi-arc row that states an AR in its notes, invariants.aphelion_ratio
    must be present and equal to the note-stated value (±0.005).

    This is the complete ratchet: no row's typed field may drift from its prose source.
    """
    rows = _load_multi_arc_rows(catalogue_rows)
    failures: list[str] = []
    for row in rows:
        rid = row["id"]
        note_ar = _notes_ar(row)
        if note_ar is None:
            continue  # no unambiguous AR in notes — nothing to check
        inv = row.get("invariants") or {}
        typed_ar = inv.get("aphelion_ratio")
        if typed_ar is None:
            failures.append(
                f"{rid}: notes state AR={note_ar} but invariants.aphelion_ratio is absent/null"
            )
        elif abs(float(typed_ar) - note_ar) > 0.005:
            failures.append(f"{rid}: notes AR={note_ar} != invariants.aphelion_ratio={typed_ar}")
    assert failures == [], f"Note→field AR consistency violations ({len(failures)}):\n" + "\n".join(
        failures
    )


def test_note_tr_matches_invariants_field(catalogue_rows: list[dict[str, Any]]) -> None:
    """For every multi-arc row that states a TR in its notes, invariants.turn_ratio
    must be present and equal to the note-stated value (±0.005).

    Rows where multiple distinct TR values appear in the notes (ambiguous) are skipped;
    those are handled in test_ambiguous_tr_rows_are_gapped.
    """
    rows = _load_multi_arc_rows(catalogue_rows)
    failures: list[str] = []
    for row in rows:
        rid = row["id"]
        note_tr = _notes_tr(row)
        if note_tr is None:
            continue  # absent or ambiguous — nothing to check here
        inv = row.get("invariants") or {}
        typed_tr = inv.get("turn_ratio")
        if typed_tr is None:
            failures.append(
                f"{rid}: notes state TR={note_tr} but invariants.turn_ratio is absent/null"
            )
        elif abs(float(typed_tr) - note_tr) > 0.005:
            failures.append(f"{rid}: notes TR={note_tr} != invariants.turn_ratio={typed_tr}")
    assert failures == [], f"Note→field TR consistency violations ({len(failures)}):\n" + "\n".join(
        failures
    )


def test_all_multi_arc_rows_have_invariants_block(catalogue_rows: list[dict[str, Any]]) -> None:
    """Every full-numeric multi-arc row must carry an invariants{} block.

    Planet-centric family-seed null-numeric rows (moon-tour Tier-1) are exempt —
    they record citation + qualitative geometry only (see
    _FAMILY_SEED_NULL_NUMERIC_MULTI_ARC).
    """
    rows = _load_full_numeric_multi_arc_rows(catalogue_rows)
    missing = [r["id"] for r in rows if "invariants" not in r]
    assert missing == [], (
        f"Multi-arc rows missing invariants{{}} block ({len(missing)}):\n" + "\n".join(missing[:20])
    )


def test_transit_times_present_for_all_multi_arc_rows(
    catalogue_rows: list[dict[str, Any]],
) -> None:
    """Every multi-arc row must have invariants.transit_times_days populated.

    Transit times (E→M outbound, M→E inbound) are always available from the
    trajectory.segments for every full-numeric multi-arc row in the catalogue.
    Planet-centric family-seed null-numeric rows (moon-tour Tier-1) are exempt —
    they have no sourced trajectory.segments (see
    _FAMILY_SEED_NULL_NUMERIC_MULTI_ARC).
    """
    rows = _load_full_numeric_multi_arc_rows(catalogue_rows)
    failures: list[str] = []
    for row in rows:
        rid = row["id"]
        inv = row.get("invariants") or {}
        tt = inv.get("transit_times_days")
        if tt is None or (isinstance(tt, list) and len(tt) == 0):
            failures.append(f"{rid}: invariants.transit_times_days is absent/null/empty")
    assert failures == [], f"Rows missing transit_times_days ({len(failures)}):\n" + "\n".join(
        failures[:20]
    )


# ---------------------------------------------------------------------------
# 2. Golden spot-checks anchored to Russell 2004 tables
# ---------------------------------------------------------------------------

# Format: (id, expected_ar, expected_tr)
# None means not tabulated in this Russell table → expect invariants.* = null
# Source for each value is in the inline comment.
GOLDEN_CHECKS: list[tuple[str, float | None, float | None]] = [
    # Russell 2004 Table 3.4: ballistic cyclers from the 2-synodic sweep
    # AR and TR are explicit columns in Table 3.4.
    ("russell-ocampo-2.3.1+1-case3", 1.08, 0.92),  # Table 3.4 row 2.3.1.+1
    ("russell-ocampo-3.1.1+3", 1.07, 1.19),  # Table 3.4 row 3.1.1.+3
    ("russell-ocampo-2.1.1+2-case2", 0.95, 1.11),  # Table 3.4 row 2.1.1.+2
    ("russell-ocampo-3.3.1+2", 1.19, 1.06),  # Table 3.4 row 3.3.1.+2
    ("russell-ocampo-3.5.1+1", 1.43, 1.15),  # Table 3.4 row 3.5.1.+1
    ("russell-ocampo-4.3.1-5", 0.99, 1.55),  # Table 3.4 row 4.3.1.-5 (AR=0.992)
    # Russell 2004 Table 4.9: generic-return 2-synodic cyclers
    # TR is explicit; AR is NOT tabulated in Table 4.9 (→ null).
    ("russell-ch4-4.991gG2", None, 2.65),  # Table 4.9 row 1 (= McConaghy 'Notable' S1L1)
    ("mcconaghy-2006-em-k2", None, 2.65),  # same cycler, McConaghy-convention entry
    # Russell 2004 Table 4.13: 3-synodic generic-return
    # TR=0.95 per Table 4.13 row 6.44 Gg (near-ballistic).
    ("russell-ch4-6.44Gg3", None, 0.95),  # Table 4.13 row 6.44 Gg (NEAR-BALLISTIC)
    # Russell 2004 Table 4.12: full-rev variants
    # TR=1.27 per dissertation line 5595-5597.
    ("russell-ch4-5.30ggF3", None, 1.27),  # Table 4.12 + dissertation line 5595-5597
]


@pytest.mark.parametrize("row_id,expected_ar,expected_tr", GOLDEN_CHECKS)
def test_golden_invariants(
    row_id: str,
    expected_ar: float | None,
    expected_tr: float | None,
    catalogue_rows: list[dict[str, Any]],
) -> None:
    """Spot-check: invariants{} fields match values read directly from the Russell tables.

    ``None`` means the value is not tabulated in that Russell table; the typed field
    must be null (or absent → treated as null).
    """
    rows = catalogue_rows
    row = next((r for r in rows if r["id"] == row_id), None)
    assert row is not None, f"Row {row_id!r} not found in catalogue"
    inv = row.get("invariants") or {}

    # Check aphelion_ratio
    actual_ar = inv.get("aphelion_ratio")
    if expected_ar is None:
        assert actual_ar is None, (
            f"{row_id}: expected invariants.aphelion_ratio=null (not tabulated), got {actual_ar!r}"
        )
    else:
        assert actual_ar is not None, (
            f"{row_id}: expected invariants.aphelion_ratio={expected_ar}, got null"
        )
        assert actual_ar == pytest.approx(expected_ar, abs=0.005), (
            f"{row_id}: invariants.aphelion_ratio={actual_ar} != expected {expected_ar}"
        )

    # Check turn_ratio
    actual_tr = inv.get("turn_ratio")
    if expected_tr is None:
        assert actual_tr is None, (
            f"{row_id}: expected invariants.turn_ratio=null (not tabulated), got {actual_tr!r}"
        )
    else:
        assert actual_tr is not None, (
            f"{row_id}: expected invariants.turn_ratio={expected_tr}, got null"
        )
        assert actual_tr == pytest.approx(expected_tr, abs=0.005), (
            f"{row_id}: invariants.turn_ratio={actual_tr} != expected {expected_tr}"
        )


def test_golden_transit_times(catalogue_rows: list[dict[str, Any]]) -> None:
    """Spot-check transit times for key rows against Russell table values.

    Table 4.9 row 1: t_out=150, t_in=150 (symmetric).
    McConaghy 2006 abstract: "153-day transfer times between Earth and Mars" (both dirs).
    """
    rows = catalogue_rows
    row_map = {r["id"]: r for r in rows}

    # russell-ch4-4.991gG2: t_out=150, t_in=150 per Russell 2004 Table 4.9 row 1
    inv = row_map["russell-ch4-4.991gG2"].get("invariants") or {}
    tt = inv.get("transit_times_days")
    assert tt is not None, "russell-ch4-4.991gG2: transit_times_days absent"
    assert tt[0] == pytest.approx(150, abs=1), f"t_out expected 150, got {tt[0]}"
    assert len(tt) >= 2 and tt[1] == pytest.approx(150, abs=1), (
        f"t_in expected 150, got {tt[1] if len(tt) >= 2 else 'missing'}"
    )

    # mcconaghy-2006-em-k2: 153 d both ways per McConaghy 2006 abstract
    inv2 = row_map["mcconaghy-2006-em-k2"].get("invariants") or {}
    tt2 = inv2.get("transit_times_days")
    assert tt2 is not None, "mcconaghy-2006-em-k2: transit_times_days absent"
    assert tt2[0] == pytest.approx(153, abs=1), f"t_out expected 153, got {tt2[0]}"
    assert len(tt2) >= 2 and tt2[1] == pytest.approx(153, abs=1), (
        f"t_in expected 153, got {tt2[1] if len(tt2) >= 2 else 'missing'}"
    )

    # russell-ocampo-3.1.1+3: t_out=174 per Russell 2004 Table 3.4 row 3.1.1.+3
    inv3 = row_map["russell-ocampo-3.1.1+3"].get("invariants") or {}
    tt3 = inv3.get("transit_times_days")
    assert tt3 is not None, "russell-ocampo-3.1.1+3: transit_times_days absent"
    assert tt3[0] == pytest.approx(174, abs=1), f"t_out expected 174, got {tt3[0]}"


# ---------------------------------------------------------------------------
# #831 (2026-08-12): the deliberate two-row carriage of ONE physical object.
#
# ``mcconaghy-2006-em-k2`` and ``russell-ch4-4.991gG2`` are the same physical
# cycler (Russell 2004 Table 4.9 row 1 = parent cycler 4.991gG2 #83 = the
# McConaghy "S1L1"/"Notable" cycler) carried as two rows, one per SOURCE
# REALIZATION. ``#831`` adjudicated this and chose KEEP-BOTH; see
# docs/notes/2026-08-12-831-mcconaghy-russell-4991gg2-duplicate-disposition.md.
#
# That disposition is only defensible while the split stays exactly where it is:
# the SHARED physics (one object) identical, the SOURCE-FLAVOURED anchors
# (two realizations) distinct. This test pins both halves so a future edit
# cannot silently converge the rows into a real duplicate, nor diverge their
# shared physics into two different objects without the ratchet noticing.
# ---------------------------------------------------------------------------

_S1L1_PAIR = ("mcconaghy-2006-em-k2", "russell-ch4-4.991gG2")


def test_831_s1l1_pair_shares_identical_physics(catalogue_rows: list[dict[str, Any]]) -> None:
    """The `#831` pair must agree on every field that carries the object's identity.

    These are the fields that make them ONE object. If any diverges, either a row
    was edited in error or they really are two objects — both cases need a human.
    """
    row_map = {r["id"]: r for r in catalogue_rows}
    a, b = (row_map[rid] for rid in _S1L1_PAIR)

    assert a["orbit_source"] == b["orbit_source"] == "russell-2004-t49_413"
    assert a["free_return_arcs"] == b["free_return_arcs"], (
        "#831 pair: free_return_arcs diverged — this descriptor pair "
        "g(1.4612,526.02,Ll)+G(2.8096,651.46,U) IS the object's structural fingerprint"
    )
    assert (a["invariants"]["turn_ratio"], b["invariants"]["turn_ratio"]) == (2.65, 2.65)
    assert a["orbit_elements"]["aphelion_au"] == b["orbit_elements"]["aphelion_au"] == 1.64
    assert a["sequence_canonical"] == b["sequence_canonical"] == "E-E-M-M"
    assert a["period"]["k"] == b["period"]["k"] == 2
    assert a["cycler_class"] == b["cycler_class"] == "multi-arc"


def test_831_s1l1_pair_keeps_distinct_source_anchors(catalogue_rows: list[dict[str, Any]]) -> None:
    """The two rows must NOT converge onto one anchor set — that is their whole reason to exist.

    Spec §16.2 bins V∞ to 0.05 km/s, so 4.7 vs 4.99 (~6 bins) and 5.0 vs 5.10
    (2 bins) give the two rows DIFFERENT canonical signatures. §16.3's "exact"
    branch (log a re-derivation, do not create a new entry) therefore never
    fires on this pair; it is at most a `probable-match-NEEDS-HUMAN`, and `#831`
    is that human confirmation. Convergence here would retroactively turn a
    lawful two-realization record into a genuine duplicate.
    """
    row_map = {r["id"]: r for r in catalogue_rows}
    mcc, rus = (row_map[rid] for rid in _S1L1_PAIR)

    assert mcc["vinf_source"] == "mcconaghy-2006"
    assert rus["vinf_source"] == "russell-2004-t49_413"

    mcc_vinf = [e["vinf_kms"] for e in mcc["vinf_kms_at_encounters"]]
    rus_vinf = [e["vinf_kms"] for e in rus["vinf_kms_at_encounters"]]
    assert mcc_vinf == [4.7, 5.0], f"McConaghy-2006 anchor changed: {mcc_vinf}"
    assert rus_vinf == [4.99, 5.10], f"Russell Table 4.9 anchor changed: {rus_vinf}"

    # Binned to spec §16.2's 0.05 km/s, the multisets must still differ.
    assert [round(v / 0.05) for v in mcc_vinf] != [round(v / 0.05) for v in rus_vinf], (
        "#831 pair: V∞ anchors now collide under the spec §16.2 0.05 km/s binning — "
        "the rows would share a canonical signature and become a true duplicate"
    )

    assert mcc["invariants"]["transit_times_days"] == [153, 153]
    assert rus["invariants"]["transit_times_days"] == [150, 150]


def test_831_s1l1_pair_validation_levels_are_not_laundered(
    catalogue_rows: list[dict[str, Any]],
) -> None:
    """`#826` held the McConaghy row at V0 because ITS OWN cited anchor is unreproduced.

    The probe on `mcconaghy-2006-em-k2` is byte-for-byte the `4.991gG2` probe
    (``build_genome`` reads `free_return_arcs`, which are Russell's, and never
    touches the 4.7/5.0 fields), so it supplies zero independent evidence for
    this row: emerged V∞ 5.008 vs the row's own cited 4.7 km/s. The V0/V3 split
    is therefore INFORMATION, not a duplication defect — and the single largest
    hazard of a merge would be laundering V0 evidence into a V3 row by absorption.
    """
    row_map = {r["id"]: r for r in catalogue_rows}
    assert row_map["mcconaghy-2006-em-k2"]["validation_level"] == "V0", (
        "#826/#831: the McConaghy-flavoured realization is V0 because its own "
        "4.7 km/s anchor is unreproduced — a promotion needs evidence against "
        "THAT anchor, not the sibling row's Russell closure"
    )
    assert row_map["russell-ch4-4.991gG2"]["validation_level"] == "V3"
