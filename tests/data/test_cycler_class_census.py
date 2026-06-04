"""Task 2.1: Census ratchet test for cycler_class tags on all 235 catalogue rows.

Verifies:
1. Every row has a ``cycler_class`` key (no missing tags).
2. The class distribution is exactly {single-ellipse: 28, multi-arc: 201, non-keplerian: 6}.
3. The set of ids tagged multi-arc exactly equals the 201-id MULTI_ARC_ALLOWLIST from
   docs/notes/multi-arc-classification.md §9 (frozen ratchet).
4. The 6 non-keplerian ids match the §3 list.

Source of truth: docs/notes/multi-arc-classification.md §3 and §9.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"

# ---------------------------------------------------------------------------
# Frozen ratchet: the 201-id MULTI_ARC_ALLOWLIST
# Source: docs/notes/multi-arc-classification.md §9
# (184 russell-ocampo-* + 14 russell-ch4-* + mcconaghy-2006-em-k2
#  + sanchez-net-2022-eem-cycler1 + sanchez-net-2022-em-cycler2)
# ---------------------------------------------------------------------------
MULTI_ARC_ALLOWLIST: frozenset[str] = frozenset(
    [
        # russell-ocampo family (184 rows) - Russell 2004 Ch3, Tables 3.4-3.11
        "russell-ocampo-2.1.1+2-case2",
        "russell-ocampo-2.3.1+1-case3",
        "russell-ocampo-2.5.1+0",
        "russell-ocampo-3.1.1+3",
        "russell-ocampo-3.1.2+1",
        "russell-ocampo-3.3.1+2",
        "russell-ocampo-3.5.1+1",
        "russell-ocampo-3.5.2+0",
        "russell-ocampo-3.7.1+1",
        "russell-ocampo-3.9.1+0",
        "russell-ocampo-4.0.3+1",
        "russell-ocampo-4.1.1-5",
        "russell-ocampo-4.10.1+2",
        "russell-ocampo-4.11.1-2",
        "russell-ocampo-4.12.1+1",
        "russell-ocampo-4.13.1-1",
        "russell-ocampo-4.14.1+0",
        "russell-ocampo-4.14.1-1",
        "russell-ocampo-4.3.1-4",
        "russell-ocampo-4.3.1-5",
        "russell-ocampo-4.5.1-3",
        "russell-ocampo-4.5.1-4",
        "russell-ocampo-4.5.2-2",
        "russell-ocampo-4.5.3-1",
        "russell-ocampo-4.7.1-3",
        "russell-ocampo-4.9.1-2",
        "russell-ocampo-4.9.2-1",
        "russell-ocampo-5.1.1-7",
        "russell-ocampo-5.1.2-3",
        "russell-ocampo-5.1.5-1",
        "russell-ocampo-5.2.1+7",
        "russell-ocampo-5.2.2+2",
        "russell-ocampo-5.2.5+0",
        "russell-ocampo-5.3.1-6",
        "russell-ocampo-5.3.1-7",
        "russell-ocampo-5.3.3-2",
        "russell-ocampo-5.4.1+5",
        "russell-ocampo-5.4.1+6",
        "russell-ocampo-5.4.3+1",
        "russell-ocampo-5.5.1-4",
        "russell-ocampo-5.5.1-5",
        "russell-ocampo-5.5.1-6",
        "russell-ocampo-5.5.2-2",
        "russell-ocampo-5.5.2-3",
        "russell-ocampo-5.5.4-1",
        "russell-ocampo-5.6.1+3",
        "russell-ocampo-5.6.1+4",
        "russell-ocampo-5.6.1+5",
        "russell-ocampo-5.6.2+1",
        "russell-ocampo-5.6.2+2",
        "russell-ocampo-5.6.4+0",
        "russell-ocampo-5.7.1-3",
        "russell-ocampo-5.7.1-4",
        "russell-ocampo-5.7.1-5",
        "russell-ocampo-5.8.1+2",
        "russell-ocampo-5.8.1+3",
        "russell-ocampo-5.8.1+4",
        "russell-ocampo-5.9.1-2",
        "russell-ocampo-5.9.1-3",
        "russell-ocampo-5.9.1-4",
        "russell-ocampo-5.9.2-1",
        "russell-ocampo-5.9.2-2",
        "russell-ocampo-5.9.3-1",
        "russell-ocampo-5.10.1+2",
        "russell-ocampo-5.10.1+3",
        "russell-ocampo-5.10.2+0",
        "russell-ocampo-5.10.2+1",
        "russell-ocampo-5.10.3+0",
        "russell-ocampo-5.11.1-2",
        "russell-ocampo-5.11.1-3",
        "russell-ocampo-5.11.2+1",
        "russell-ocampo-5.12.1+1",
        "russell-ocampo-5.12.1+2",
        "russell-ocampo-5.13.1-2",
        "russell-ocampo-5.13.1-3",
        "russell-ocampo-5.13.2-1",
        "russell-ocampo-5.14.1+1",
        "russell-ocampo-5.14.1+2",
        "russell-ocampo-5.14.2+0",
        "russell-ocampo-5.15.1-1",
        "russell-ocampo-5.15.1-2",
        "russell-ocampo-5.16.1+0",
        "russell-ocampo-5.16.1+1",
        "russell-ocampo-5.17.1-1",
        "russell-ocampo-5.18.1+0",
        "russell-ocampo-6.0.1+6d",
        "russell-ocampo-6.0.1+7c",
        "russell-ocampo-6.0.1+8b",
        "russell-ocampo-6.0.1+9a",
        "russell-ocampo-6.1.2-4",
        "russell-ocampo-6.1.3-3",
        "russell-ocampo-6.1.4-2",
        "russell-ocampo-6.1.6-1",
        "russell-ocampo-6.2.1+6",
        "russell-ocampo-6.2.1+7",
        "russell-ocampo-6.2.1+8",
        "russell-ocampo-6.2.2+2",
        "russell-ocampo-6.2.2+3",
        "russell-ocampo-6.2.3+1",
        "russell-ocampo-6.2.3+2",
        "russell-ocampo-6.2.4+1",
        "russell-ocampo-6.2.6+0",
        "russell-ocampo-6.3.1-9",
        "russell-ocampo-6.3.4+1",
        "russell-ocampo-6.4.1+4",
        "russell-ocampo-6.4.1+5",
        "russell-ocampo-6.4.1+6",
        "russell-ocampo-6.4.1+7",
        "russell-ocampo-6.5.1-6",
        "russell-ocampo-6.5.1-7",
        "russell-ocampo-6.5.1-8",
        "russell-ocampo-6.5.5-1",
        "russell-ocampo-6.6.1+3",
        "russell-ocampo-6.6.1+4",
        "russell-ocampo-6.6.1+5",
        "russell-ocampo-6.6.1+6",
        "russell-ocampo-6.6.2+1",
        "russell-ocampo-6.6.2+2",
        "russell-ocampo-6.6.5+0",
        "russell-ocampo-6.7.1-6",
        "russell-ocampo-6.7.1-7",
        "russell-ocampo-6.7.2+3",
        "russell-ocampo-6.7.3-2",
        "russell-ocampo-6.7.5+0",
        "russell-ocampo-6.8.1+2",
        "russell-ocampo-6.8.1+3",
        "russell-ocampo-6.8.1+4",
        "russell-ocampo-6.8.1+5",
        "russell-ocampo-6.8.1+6",
        "russell-ocampo-6.8.3+1",
        "russell-ocampo-6.9.1-4",
        "russell-ocampo-6.9.1-5",
        "russell-ocampo-6.9.1-6",
        "russell-ocampo-6.9.2-2",
        "russell-ocampo-6.9.2-3",
        "russell-ocampo-6.9.4-1",
        "russell-ocampo-6.10.1+2",
        "russell-ocampo-6.10.1+3",
        "russell-ocampo-6.10.1+4",
        "russell-ocampo-6.10.1+5",
        "russell-ocampo-6.10.2+1",
        "russell-ocampo-6.10.2+2",
        "russell-ocampo-6.10.4+0",
        "russell-ocampo-6.11.1-4",
        "russell-ocampo-6.11.1-5",
        "russell-ocampo-6.11.2+2",
        "russell-ocampo-6.12.1+2",
        "russell-ocampo-6.12.1+3",
        "russell-ocampo-6.12.1+4",
        "russell-ocampo-6.13.1-3",
        "russell-ocampo-6.13.1-4",
        "russell-ocampo-6.13.1-5",
        "russell-ocampo-6.13.1+5",
        "russell-ocampo-6.13.2-2",
        "russell-ocampo-6.13.3-1",
        "russell-ocampo-6.14.1+1",
        "russell-ocampo-6.14.1+2",
        "russell-ocampo-6.14.1+3",
        "russell-ocampo-6.14.2+0",
        "russell-ocampo-6.14.2+1",
        "russell-ocampo-6.14.3+0",
        "russell-ocampo-6.15.1-2",
        "russell-ocampo-6.15.1-3",
        "russell-ocampo-6.15.1-4",
        "russell-ocampo-6.15.1+4",
        "russell-ocampo-6.15.2+1",
        "russell-ocampo-6.16.1+1",
        "russell-ocampo-6.16.1+2",
        "russell-ocampo-6.17.1-2",
        "russell-ocampo-6.17.1-3",
        "russell-ocampo-6.17.1+3",
        "russell-ocampo-6.17.2-1",
        "russell-ocampo-6.18.1+1",
        "russell-ocampo-6.18.1+2",
        "russell-ocampo-6.18.2+0",
        "russell-ocampo-6.19.1-2",
        "russell-ocampo-6.19.1+2",
        "russell-ocampo-6.19.2+0",
        "russell-ocampo-6.20.1-4",
        "russell-ocampo-6.20.1+0",
        "russell-ocampo-6.20.1+1",
        "russell-ocampo-6.21.1-1",
        "russell-ocampo-6.21.1+1",
        "russell-ocampo-6.22.1+0",
        # russell-ch4 family (14 rows) - Russell 2004 Ch4, Tables 4.9-4.13
        "russell-ch4-4.991gG2",
        "russell-ch4-8.049gGf2",
        "russell-ch4-8.165Gfh-f2",
        "russell-ch4-9.353Gg2",
        "russell-ch4-3.64gGg3",
        "russell-ch4-3.77Gh3",
        "russell-ch4-3.78Gg3",
        "russell-ch4-5.30gGf3",
        "russell-ch4-5.66Gfh3",
        "russell-ch4-9.94Gg3",
        "russell-ch4-3.66gfF3",
        "russell-ch4-5.30ggF3",
        "russell-ch4-5.75ggF3",
        "russell-ch4-6.44Gg3",
        # resolved from §7 — same physical cycler as russell-ch4-4.991gG2 (Russell Table 4.9)
        "mcconaghy-2006-em-k2",
        # Sanchez Net 2022 EEM near-ballistic real-date patched-conic cycler (Fig. 2a)
        "sanchez-net-2022-eem-cycler1",
        # Sanchez Net 2022 EM near-ballistic real-date patched-conic cycler (Fig. 2b)
        "sanchez-net-2022-em-cycler2",
    ]
)

assert len(MULTI_ARC_ALLOWLIST) == 201, (
    f"Allowlist must have 201 entries, got {len(MULTI_ARC_ALLOWLIST)}"
)

# ---------------------------------------------------------------------------
# Frozen ratchet: the 6 non-keplerian ids
# Source: docs/notes/multi-arc-classification.md §3
# ---------------------------------------------------------------------------
NON_KEPLERIAN_IDS: frozenset[str] = frozenset(
    [
        "arenstorf-em-figure8-1963",
        "genova-aldrin-2015-em-3petal-cycler",
        "wittal-2022-em-cycler-family",
        "hernandez-2017-jovian-ieg-triple-family",
        "russell-strange-2009-jovian-multimoon-family",
        "russell-strange-2009-saturnian-multimoon-family",
    ]
)

assert len(NON_KEPLERIAN_IDS) == 6


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _load_rows() -> list[dict]:  # type: ignore[type-arg]
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_all_rows_have_cycler_class() -> None:
    """Every row in catalogue.yaml must have a cycler_class key."""
    rows = _load_rows()
    missing = [r["id"] for r in rows if "cycler_class" not in r]
    assert missing == [], f"Rows missing cycler_class ({len(missing)} rows): {missing[:10]}"


def test_census_distribution() -> None:
    """Exact class distribution: single-ellipse=28, multi-arc=200, non-keplerian=6."""
    rows = _load_rows()
    counts = Counter(r.get("cycler_class", "single-ellipse") for r in rows)
    expected = {"single-ellipse": 28, "multi-arc": 201, "non-keplerian": 6}
    assert dict(counts) == expected, (
        f"Census mismatch.\n  Expected: {expected}\n  Got:      {dict(counts)}"
    )


def test_multi_arc_ids_match_allowlist() -> None:
    """The exact set of multi-arc ids matches the 201-id MULTI_ARC_ALLOWLIST ratchet."""
    rows = _load_rows()
    actual = frozenset(r["id"] for r in rows if r.get("cycler_class") == "multi-arc")
    extra = actual - MULTI_ARC_ALLOWLIST
    missing_from_actual = MULTI_ARC_ALLOWLIST - actual
    assert actual == MULTI_ARC_ALLOWLIST, (
        f"multi-arc id mismatch.\n"
        f"  In catalogue but NOT in allowlist ({len(extra)}): {sorted(extra)}\n"
        f"  In allowlist but NOT in catalogue ({len(missing_from_actual)}):"
        f" {sorted(missing_from_actual)}"
    )


def test_non_keplerian_ids_match_ratchet() -> None:
    """The exact set of non-keplerian ids matches the 6-id NON_KEPLERIAN_IDS ratchet."""
    rows = _load_rows()
    actual = frozenset(r["id"] for r in rows if r.get("cycler_class") == "non-keplerian")
    extra = actual - NON_KEPLERIAN_IDS
    missing_from_actual = NON_KEPLERIAN_IDS - actual
    assert actual == NON_KEPLERIAN_IDS, (
        f"non-keplerian id mismatch.\n"
        f"  In catalogue but NOT in ratchet ({len(extra)}): {sorted(extra)}\n"
        f"  In ratchet but NOT in catalogue ({len(missing_from_actual)}):"
        f" {sorted(missing_from_actual)}"
    )
