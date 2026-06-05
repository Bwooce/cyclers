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


@pytest.mark.xfail(
    reason="VEM strict ballistic closure is open research (spec §17 line 199); "
    "no sourced V_inf anchor exists to assert against and the circular-coplanar "
    "optimiser is not expected to converge a ballistic VEM cycler. Flipped by "
    "M-ED (real-ephemeris discovery). See roadmap M-N test-gate row.",
    strict=False,
)
@pytest.mark.slow
def test_emeeve_idealized_optimiser_converges_feasible() -> None:
    """ASPIRATIONAL: the idealized optimiser finds a feasible (constraints-
    satisfied) interior solution for the EMEEVE VEM cell. Expected to xfail in
    the circular-coplanar model; documents the M-ED handoff target.

    This asserts ONLY result.constraints_satisfied (a feasibility predicate,
    not a sourced number) — it never asserts a computed V_inf as golden.
    """
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import optimise_cell_idealized

    # EMEEVE with loop leg requires multi-rev (M-L). Until M-L lands this also
    # cannot construct; the xfail covers both the M-L and the convergence gaps.
    # Built as the loader would: anchor pair (E,M), sourced k=3, no k rewrite.
    seq = ("E", "M", "E", "E", "V", "E")
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=3,
        per_leg_revs=(0, 0, 1, 0, 0),  # the E-E loop leg is multi-rev
        per_leg_branch=("single", "single", "low", "single", "single"),
        period_basis=("E", "M"),
    )
    result = optimise_cell_idealized(
        cell,
        Ephemeris(model="circular"),
        vinf_cap=7.0,
        seed=0,
    )
    assert result.constraints_satisfied
