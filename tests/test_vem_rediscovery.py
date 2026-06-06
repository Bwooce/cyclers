"""M8 VEM rediscovery gate (spec §8 M8 anchor, §9 line 160).

GOLDEN DISCIPLINE: the EXPECTED side of every assertion below is a *sourced*
catalogue value (``period.years`` and ``sequence_canonical``) — never a number
our own optimiser computes. The converged-geometry case is xfail pending M-ED
(real-ephemeris discovery) — see plan §4 Task 4.2.

R1 (2026-06-05, Jones AAS 17-577 ingest). The primary period round-trip anchors
on the two *member* rows (``jones-2017-vem-emevve-outbound`` and
``jones-2017-vem-meevem-inbound``), whose sourced period is **12.8 yr** and
whose ``period.pair`` is the **beat token** ``"VEM-syn"`` (NOT a body pair).
The loader maps a non-2-body-pair token to ``period_basis=None``, so the
resolver falls back to the natural VEM beat ``multi_body_beat_days[0]`` (6.406
yr); with the sourced ``k=2`` this gives ``6.406*2 = 12.81 yr``, agreeing with
the sourced 12.8 within rounding (the gate asserts that agreement).

``vem-emeeve-3syn`` (period.pair ``"E-M"``, k=3, 6.41 yr) is **UNREALIZED**
(Jones p.8: no feasible 1-synodic / 6.4-yr family member exists; realized
cyclers repeat at 12.8 yr). Its 6.41-yr value survives only as the beat
archetype and is used here purely as an *anchor-pair resolver-arithmetic*
check, not as a realized-cycler claim.
"""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.search.optimize import _target_period_sec
from cyclerfinder.search.sequence import Cell
from tests._catalogue_loader import CATALOGUE_PATH

# Single-letter body codes used to decide whether a period.pair token is a
# 2-body anchor pair (e.g. "E-M") or a beat token (e.g. "VEM-syn").
_BODY_CODES = frozenset({"V", "E", "M", "S", "J"})


def _row(entry_id: str) -> dict[str, Any]:
    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    for row in raw:
        if row["id"] == entry_id:
            return cast("dict[str, Any]", row)
    raise AssertionError(f"catalogue row {entry_id!r} not found")


def _basis_from_pair(pair: str | None) -> tuple[str, str] | None:
    """Map a ``period.pair`` token to a Cell anchor pair, or None for a beat
    token (plan §5 R1 delta 3). Only ``A-B`` over known single-letter body
    codes is a real pair; ``"VEM-syn"`` and friends fall back to None so the
    resolver uses the natural beat (and we never call
    ``synodic_period_days("VEM","syn")``)."""
    if not pair:
        return None
    parts = pair.split("-")
    if len(parts) == 2 and all(p in _BODY_CODES for p in parts):
        return (parts[0], parts[1])
    return None


# ---------------------------------------------------------------------------
# Primary anchor: the two member rows (sourced 12.8 yr, beat-token pair)
# ---------------------------------------------------------------------------

_MEMBER_ROW_IDS = (
    "jones-2017-vem-emevve-outbound",
    "jones-2017-vem-meevem-inbound",
)


@pytest.mark.parametrize("entry_id", _MEMBER_ROW_IDS)
def test_member_row_sourced_period_round_trips(entry_id: str) -> None:
    """A Jones member row's SOURCED 12.8-yr period round-trips through the
    resolver via the natural VEM beat.

    EXPECTED side = the catalogue's published ``period.years`` (sourced from
    Jones AAS 17-577). The resolver output is the side under test. ``period.pair
    = "VEM-syn"`` is a beat token, so ``period_basis=None`` and the resolver
    uses ``beat_period_days`` (6.406 yr) * sourced k=2 = 12.81 yr.
    """
    row = _row(entry_id)
    sourced_years = float(row["period"]["years"])  # 12.8, Jones AAS 17-577
    assert sourced_years == pytest.approx(12.8, abs=0.005)  # fixture sanity

    seq = tuple(row["sequence_canonical"].split("-"))
    basis = _basis_from_pair(row["period"]["pair"])  # "VEM-syn" -> None
    assert basis is None  # beat token, not an anchor pair (R1 delta 3)
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=int(row["period"]["k"]),  # 2, traceable to the YAML
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=basis,
    )
    resolved_years = _target_period_sec(cell) / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    # Sourced 12.8 yr vs beat 6.406*2 = 12.81 yr: agree within inter-source
    # rounding (Jones rounds the beat to 6.4; our constants give 6.406).
    assert resolved_years == pytest.approx(sourced_years, abs=0.05)


@pytest.mark.parametrize("entry_id", _MEMBER_ROW_IDS)
def test_member_row_sequence_is_sourced(entry_id: str) -> None:
    """The member rows carry a sourced ``sequence_canonical`` (11-encounter
    Jones itineraries) — assert it parses to a VEM-only body set with the
    sourced length, the golden-legal sequence anchor."""
    row = _row(entry_id)
    seq = tuple(row["sequence_canonical"].split("-"))
    assert len(seq) >= 6  # multi-encounter itinerary
    assert set(seq) <= {"V", "E", "M"}


# ---------------------------------------------------------------------------
# Loader-driven generalisation: every CONSTRUCTIBLE_MULTIBODY row's sourced
# period round-trips the resolver (Task 5.4)
# ---------------------------------------------------------------------------


def _multibody_entries() -> list[str]:
    from tests._catalogue_loader import ExclusionReason, classify_catalogue

    return [
        cid
        for cid, reason in classify_catalogue()
        if reason is ExclusionReason.CONSTRUCTIBLE_MULTIBODY
    ]


def test_all_multibody_rows_period_round_trips_via_loader() -> None:
    """Every CONSTRUCTIBLE_MULTIBODY row's SOURCED period round-trips through
    the resolver, using the loader-carried ``period_basis`` (no YAML re-read of
    the pair). This generalises the per-row anchors across all four VEM rows.

    EXPECTED side = each row's catalogue ``period.years``. The basis comes from
    the loader's ``CatalogueEntry.period_basis`` (``("E","M")`` for EMEEVE,
    ``None`` for the beat-token member/family rows), proving the loader feeds
    the resolver correctly.
    """
    from tests._catalogue_loader import classify_row

    ids = _multibody_entries()
    assert set(ids) == {
        "jones-2017-vem-triple-family",
        "vem-emeeve-3syn",
        "jones-2017-vem-emevve-outbound",
        "jones-2017-vem-meevem-inbound",
    }
    for entry_id in ids:
        row = _row(entry_id)
        _reason, entry = classify_row(row)
        assert entry is not None
        seq = tuple(row["sequence_canonical"].split("-"))
        n_legs = max(len(seq) - 1, 1)
        cell = Cell(
            bodies=("V", "E", "M"),
            sequence=seq if len(seq) >= 2 else ("V", "E", "M"),
            period_k=entry.period_k,
            per_leg_revs=(0,) * n_legs,
            per_leg_branch=("single",) * n_legs,
            period_basis=entry.period_basis,
        )
        resolved_years = _target_period_sec(cell) / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
        assert resolved_years == pytest.approx(entry.period_years, abs=0.05), (
            f"{entry_id}: resolved {resolved_years:.4f} vs sourced {entry.period_years}"
        )


# ---------------------------------------------------------------------------
# Resolver-arithmetic check: the UNREALIZED EMEEVE archetype (anchor pair)
# ---------------------------------------------------------------------------


def test_emeeve_unrealized_period_matches_anchor_resolver() -> None:
    """The EMEEVE row's SOURCED period (6.41 yr, catalogue.yaml) matches what
    _target_period_sec resolves from the row's anchor pair + sourced k.

    NOTE: ``vem-emeeve-3syn`` is UNREALIZED (Jones p.8: no feasible 1-synodic
    6.4-yr family member exists; realized cyclers repeat at 12.8 yr). This is a
    resolver *arithmetic* check on the beat archetype, NOT a realized-cycler
    claim. EXPECTED side = the catalogue's published period.years.
    """
    row = _row("vem-emeeve-3syn")
    sourced_years = float(row["period"]["years"])  # 6.41
    assert sourced_years == pytest.approx(6.41, abs=0.005)  # fixture sanity

    seq = tuple(row["sequence_canonical"].split("-"))  # E,M,E,E,V,E
    basis = _basis_from_pair(row["period"]["pair"])  # "E-M" -> ("E","M")
    assert basis == ("E", "M")
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=int(row["period"]["k"]),  # 3, traceable to the YAML
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=basis,
    )
    resolved_years = _target_period_sec(cell) / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    # Sourced 6.41 yr vs 3*T_syn(E,M) ~ 6.405 yr: agree within rounding.
    assert resolved_years == pytest.approx(sourced_years, abs=0.02)


def test_vem_syn_beat_token_never_parsed_as_body_pair() -> None:
    """Guard (R1 delta 3): the beat token "VEM-syn" must NOT be split into a
    body pair. Splitting it would crash synodic_period_days("VEM","syn"); the
    loader/gate map it to period_basis=None instead."""
    assert _basis_from_pair("VEM-syn") is None
    assert _basis_from_pair("E-M") == ("E", "M")
    assert _basis_from_pair(None) is None


# ---------------------------------------------------------------------------
# M-ED handoff: full ballistic convergence (documented xfail)
# ---------------------------------------------------------------------------


# M-ED HEADLINE GATE: Jones VEM ballistic rediscovery (plan Phase 5 Task 5.1).
#
# GOLDEN DISCIPLINE: EXPECTED = the catalogue's SOURCED vinf_kms_at_encounters
# (AAS 17-577 Tables 2/3). The corrector output is the side under test. No
# self-computed value is ever the EXPECTED side.
#
# RISK (spec §7): the S1L1 corrector family floors Mars V_inf ~6.4 km/s; the Jones
# Mars targets are 2.50/2.79 (EMEVVE) and ~2.42-3.85 (MEEVEM). Convergence is NOT
# assumed. Until it converges within VEM_VINF_TOL_KMS this stays xfail with the
# finding recorded; the STOP/report branch (Task 5.4) governs that outcome.

VEM_VINF_TOL_KMS = 0.5  # tied to sourced Jones rounding + model difference (spec §5)


def _sourced_vinf_multiset(entry_id: str) -> list[float]:
    return sorted(float(e["vinf_kms"]) for e in _row(entry_id)["vinf_kms_at_encounters"])


def _sourced_cycle_tofs(entry_id: str) -> list[float]:
    """First-cycle per-leg ToFs (days) from the SOURCED trajectory segments.

    One cycle has ``len(sequence_canonical) - 1`` legs; the row's segments list
    two cycles (the 11-encounter Table 2/3 span), so the first ``n_legs`` segment
    ToFs are the single-cycle seed. These are SOURCED (Table 2/3 Flight Time
    column), used only as the anchor-rung *seed* — never an EXPECTED.
    """
    row = _row(entry_id)
    n_legs = len(row["sequence_canonical"].split("-")) - 1
    segs = row["trajectory"]["segments"]
    return [float(segs[i]["tof_days"]) for i in range(n_legs)]


@pytest.mark.xfail(
    reason="M-ED HEADLINE GATE: ballistic VEM rediscovery to the sourced Jones "
    "multiset within 0.5 km/s. STOP/report outcome (plan Task 5.4): a DENSE "
    "parallel epoch x branch scan (task #110) STILL FAILS, so the xfail is NOT "
    "flipped and the tolerance is NOT loosened. Hunt 2026-06-06 "
    "(scripts/hunt_vem_ballistic.py, 256 epochs over the full 12.8-yr repeat "
    "period x 11 rev/branch topologies = 2816 points/row, parallel 16-core): "
    "EMEVVE outbound = 831 closed / 474 distinct families, BEST max-V_inf 17.86 "
    "(per-encounter [13.88,13.91,15.39,15.39,17.43,17.86]); MEEVEM inbound = "
    "1239 closed / 570 distinct families, BEST max-V_inf 18.49 (per-encounter "
    "[11.40,16.34,16.34,17.83,18.47,18.49]). ZERO bend-feasible solutions in "
    "either survey. The closed families floor ~11-18 km/s, far above the sourced "
    "Jones 2.42-7.0 km/s -- the S1L1 Mars-V_inf ~6.4 floor generalised: the "
    "single-ellipse-per-leg corrector closes a DIFFERENT, higher-V_inf, powered "
    "family than the Jones members. The sharpened hypothesis: the Jones VEM "
    "family needs 3D inclination (M-3D), real-eccentricity intermediate "
    "flybys, or a different (e.g. multi-arc per leg) topology seeding -- not "
    "more scan density. Per-cycle corrector yields 6 encounter V_inf; the "
    "sourced multiset is the 11-encounter two-cycle Table 2/3 span, so the "
    "strict-multiset compare also surfaces a length gap. Flip ONLY when a member "
    "row genuinely converges within VEM_VINF_TOL_KMS; see plan Phase 5 Task 5.1.",
    strict=False,
)
@pytest.mark.slow
@pytest.mark.parametrize("entry_id", _MEMBER_ROW_IDS)
def test_jones_vem_ballistic_rediscovers_sourced_multiset(entry_id: str) -> None:
    """HEADLINE GATE: the M-ED ballistic corrector, seeded via the sourced-anchor
    rung (the row's sourced transit/segment ToFs + per-encounter V_inf targets),
    converges to a closed ballistic chain whose per-encounter V_inf magnitudes
    match the row's SOURCED multiset within VEM_VINF_TOL_KMS.

    EXPECTED side = the published Jones multiset only (golden discipline).
    """
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import optimise_cell_ephemeris
    from cyclerfinder.search.seed_ladder import resolve_seed

    row = _row(entry_id)
    seq = tuple(row["sequence_canonical"].split("-"))
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=int(row["period"]["k"]),
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=_basis_from_pair(row["period"]["pair"]),  # "VEM-syn" -> None
    )

    # Sourced anchor: per-body V_inf targets (first occurrence of each body in
    # encounter order) + the first-cycle sourced segment ToFs.
    vinf_targets: dict[str, float] = {}
    for enc in row["vinf_kms_at_encounters"]:
        vinf_targets.setdefault(enc["body"], float(enc["vinf_kms"]))
    seed_plan = resolve_seed(
        cell,
        anchor_tofs=_sourced_cycle_tofs(entry_id),
        anchor_vinf=vinf_targets,
    )
    assert seed_plan.source == "anchor"

    # Priority epoch = the row's sourced t-zero (Table 2/3 first encounter date).
    # Drive the scan rung (task #110): a parallel epoch grid across one full
    # 12.8-yr repeat period is the density lever family selection needs (spec
    # 3.4; the prototype's main() scan). The single start lands the degenerate
    # high-V_inf basin; the dense scan samples the whole period for a member.
    priority = row["trajectory"]["epoch_tzero"]
    result = optimise_cell_ephemeris(
        cell,
        Ephemeris("astropy"),
        vinf_cap=8.0,
        priority_date_iso=priority,
        vinf_targets_kms=vinf_targets,
        tof_seed_days=list(seed_plan.tof_seed_days),
        mode="ballistic",
        scan_epochs=64,
    )
    assert result.converged and result.constraints_satisfied

    expected = _sourced_vinf_multiset(entry_id)
    got = sorted(
        max(float(np.linalg.norm(e.vinf_in)), float(np.linalg.norm(e.vinf_out)))
        for e in result.best_cycler.encounters
    )
    for g, x in zip(got, expected, strict=True):
        assert g == pytest.approx(x, abs=VEM_VINF_TOL_KMS)
