"""Task 4 remainder: live-row validation-tier census ratchet.

Runs :func:`cyclerfinder.data.provenance.classify_validation` over every row of
``data/catalogue.yaml`` using the per-field provenance tags back-filled by Task 3
(``orbit_source`` / ``vinf_source`` / ``orbit_fidelity`` / ``vinf_fidelity``),
and freezes the resulting tier distribution — the validation-strength sibling of
``test_cycler_class_census.py``'s ``cycler_class`` census.

Any catalogue change that shifts a row's provenance (and hence its tier) must
update the frozen counts AND the ``CROSS_VALIDATED_IDS`` set in the SAME commit,
making the strongly-validated set visible and monotone: a row cannot silently
fall out of ``cross_validated``, and no row can silently claim it.

The classification here is computed from the catalogue's own tags only — it is
not a golden physics value, so freezing it is a coverage ratchet, not a sourced
EXPECTED check (cf. the golden-tests discipline: the EXPECTED side of a *golden*
test must trace to a published source; a tier census freezes our own
provenance bookkeeping, which is the thing under test).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.provenance import Tier, classify_validation

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"

# ---------------------------------------------------------------------------
# Frozen ratchet: the live validation-tier distribution (2026-06-05, v4.4).
# Update in the same commit as any provenance-tag change.
#
# 2026-06-07 (#142 catalogue ingest, batch 1): +15 Rall 1970 (rall-1970-*)
# multi-arc rows, each tagged orbit_source == vinf_source == rall-1970-te34 at
# one fidelity (App E circular-coplanar / App F analytic-ephemeris) =>
# CONSISTENCY_CHECKED. consistency_checked 218 -> 233.
#
# 2026-06-07 (#142 catalogue ingest, batch 2): +16 Russell 2004 Table 3.4
# circular-coplanar cyclers (russell-ocampo-* not previously catalogued), each
# tagged orbit_source == vinf_source == russell-2004-t34 at one fidelity
# (circular-coplanar) => CONSISTENCY_CHECKED. consistency_checked 233 -> 249.
#
# 2026-06-12 (#216, USER-approved writeback): +5 Ross & Roberts-Tsoukkas 2025
# (AAS 25-621) stable Earth-Moon CR3BP cycler rows (ross-rt-em-cycler-*). These
# carry no orbit_source/vinf_source provenance tags (the SOURCE_REGISTRY keys
# are heliocentric Earth-Mars / patched-conic citations; the Ross CR3BP source
# is not a registry key, and their identity is a Jacobi/period tuple, not the
# orbit/vinf-multiset pair the tier classifier reads), so they classify
# UNVALIDATED by the provenance-tag tier system — exactly as the existing CR3BP
# rows (arenstorf, wittal, genova-aldrin) do. Their gauntlet strength lives in
# the orthogonal validation_level axis (V1, _LEVEL_EVIDENCE), NOT this
# orbit/vinf provenance-tag census. unvalidated 14 -> 19.
#
# 2026-06-12 (#216, USER-approved writeback): +4 Liang et al. 2024 (JGCD doi
# 10.2514/1.G008387) Callisto-Ganymede-Europa triple-cycler members
# (liang-2024-cgcec-*). These are Jovian moon-tour multi-arc rows whose identity
# is a per-flyby ToF/V_inf table, not the heliocentric orbit/vinf-multiset pair
# the SOURCE_REGISTRY tier classifier reads (the Liang source is not a registry
# key), so they classify UNVALIDATED by the provenance-tag tier system — exactly
# as the existing Jovian seeds (hernandez, russell-strange) do. Members A-C are
# idealized-model V0; member D is ephemeris figures-only V0 (validation_level
# axis). unvalidated 19 -> 23.
#
# 2026-06-15 (#249, USER-approved writeback): +3 Braik & Ross 2026 (arXiv
# 2605.31543) common-energy Earth-Moon CR3BP cycler reproductions
# (braik-ross-c11a/c11b/c32-cycler-2026). Like the Ross-RT 2025 rows, these
# carry no orbit_source/vinf_source provenance tags (their identity is a
# Jacobi/period tuple, not the orbit/vinf-multiset pair the SOURCE_REGISTRY
# tier classifier reads), so they classify UNVALIDATED on this orthogonal
# provenance-tag axis. Their gauntlet strength lives in the validation_level
# axis (V1, _LEVEL_EVIDENCE). unvalidated 23 -> 26.
# ---------------------------------------------------------------------------
# 249 -> 256 (2026-06-17, #367): seven Rogers 2015 Table 4 precursor_mga
# rows (VISIT-1/2, Case 1/2/3, S1L1, U0L1) admitted with
# orbit_source=vinf_source=rogers-2012-t1 at fidelity=analytic-ephemeris.
# Single-source same-fidelity pair classifies as CONSISTENCY_CHECKED
# under the existing tier rules (matches the Aldrin '4:3(2)-' and
# '3:2(1)-' Rogers establishment rows already in this bucket).
# 256 -> 257 (2026-06-17, #367 wave 2): one Rogers 2015 Table 3
# circular-coplanar precursor_mga sub-variant (VISIT-1 5:4(3)-) admitted
# with orbit_source=vinf_source=rogers-2012-t1 at fidelity=circular-
# coplanar. Same CONSISTENCY_CHECKED classification.
EXPECTED_TIER_CENSUS: dict[str, int] = {
    "cross_validated": 5,
    "consistency_checked": 257,
    # 26 -> 27 (2026-06-15, #294): Tito 2018 Mars free-return admitted as
    # mga_tour with orbit_source=derived / vinf_source=derived (same fidelity)
    # which classifies as 'unvalidated' under the existing provenance rules.
    # The row's V0 evidence (the DE440 reproduction <1.5% to Tito's published
    # DE421) lives in validation_level, not the source-pair tier — separate
    # axis. See docs/notes/2026-06-16-catalogue-scope-taxonomy.md.
    # 27 -> 28 (2026-06-16, #336): Heaton-Longuski 2003 Uranian satellite tour
    # U00-01 admitted as mga_tour (second mga_tour row, after Tito 2018) with
    # orbit_source=derived / vinf_source=derived (same fidelity) — also
    # classifies as 'unvalidated' on the provenance-tag tier axis. Its V0
    # evidence (Tables 3 + 5 verbatim from the JSR paper, DOI 10.2514/2.3981)
    # lives in validation_level, not the source-pair tier — separate axis.
    # 28 -> 29 (2026-06-17, #339): umbriel-oberon-1-1-uranian-quasi-cycler-2026
    # admitted as catalogue's first computed quasi_cycler row with
    # orbit_source=derived / vinf_source=derived (same fidelity) — classifies
    # as 'unvalidated' on the provenance-tag tier axis. Its evidence
    # (validation_level=V4, gated on #335 V4-strict + #338 annual epoch
    # sweep EFFECTIVELY_CYCLIC) lives on a SEPARATE axis from this
    # source-pair census. The two axes are orthogonal by design (see
    # earlier #294/#336 comments above).
    # 29 -> 30 (2026-06-17, #356): damario-1992-galileo-veega admitted as
    # the third computed mga_tour row with orbit_source=derived /
    # vinf_source=derived (same fidelity) — classifies as 'unvalidated'
    # under the same orthogonal-axis convention. Its V0 evidence
    # (D'Amario 1992 SSR 60(1-4):23-78 Table I + Figs 3/5/8 verbatim) lives
    # in validation_level, not the source-pair tier.
    "unvalidated": 30,
}

# The exact set of CROSS_VALIDATED rows: each pairs two DIFFERENT independent
# citations at the SAME fidelity (orbit_source != vinf_source). Frozen so the
# strong set cannot silently change.
CROSS_VALIDATED_IDS: frozenset[str] = frozenset(
    [
        "aldrin-classic-em-k1-outbound",  # orbit Rogers 2012 T1 / vinf Russell 2004 T3.4
        "aldrin-classic-em-k1-inbound",  # orbit Rogers 2012 T1 / vinf Russell 2004 T3.4
        "mcconaghy-2006-em-k2",  # orbit Russell 2004 T4.9 / vinf McConaghy 2006
        "mcconaghy-2005-em-case1",  # orbit Rogers 2012 T1 / vinf McConaghy 2002
        "s1l1-2syn-em-cpom",  # orbit Rogers 2012 T1 / vinf McConaghy 2002
    ]
)


def _load_rows() -> list[dict]:  # type: ignore[type-arg]
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def _row_tier(row: dict) -> Tier:  # type: ignore[type-arg]
    """Classify a single row from its back-filled provenance tags."""
    orbit_fid = row.get("orbit_fidelity")
    vinf_fid = row.get("vinf_fidelity")
    same_fid = orbit_fid is not None and vinf_fid is not None and orbit_fid == vinf_fid
    return classify_validation(
        row.get("orbit_source"),
        row.get("vinf_source"),
        same_fidelity=same_fid,
    )


def test_tier_census_distribution() -> None:
    """Live tier distribution exactly matches the frozen census."""
    rows = _load_rows()
    counts = Counter(_row_tier(r).value for r in rows)
    assert dict(counts) == EXPECTED_TIER_CENSUS, (
        f"Validation-tier census mismatch.\n"
        f"  Expected: {EXPECTED_TIER_CENSUS}\n"
        f"  Got:      {dict(counts)}\n"
        "If this is an intended provenance change, update EXPECTED_TIER_CENSUS "
        "(and CROSS_VALIDATED_IDS) in the same commit."
    )


def test_cross_validated_ids_match_ratchet() -> None:
    """The exact set of cross_validated rows matches the frozen id set."""
    rows = _load_rows()
    actual = frozenset(r["id"] for r in rows if _row_tier(r) is Tier.CROSS_VALIDATED)
    extra = actual - CROSS_VALIDATED_IDS
    missing = CROSS_VALIDATED_IDS - actual
    assert actual == CROSS_VALIDATED_IDS, (
        f"cross_validated id set mismatch.\n"
        f"  In catalogue but NOT in ratchet ({len(extra)}): {sorted(extra)}\n"
        f"  In ratchet but NOT in catalogue ({len(missing)}): {sorted(missing)}"
    )


def test_cross_validated_rows_have_distinct_sources() -> None:
    """Every cross_validated row genuinely pairs two different sources at one
    fidelity (guards the classifier's independence rule against a tag typo)."""
    rows = _load_rows()
    byid = {r["id"]: r for r in rows}
    for rid in CROSS_VALIDATED_IDS:
        row = byid[rid]
        os_, vs = row.get("orbit_source"), row.get("vinf_source")
        assert os_ is not None and vs is not None, f"{rid}: missing a source tag"
        assert os_ != vs, f"{rid}: cross_validated but shares source {os_!r}"
        assert row.get("orbit_fidelity") == row.get("vinf_fidelity"), (
            f"{rid}: cross_validated requires matching fidelity"
        )
