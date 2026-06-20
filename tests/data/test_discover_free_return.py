"""SnLm free-return discover-reachability gate (task #106 Phase-2 re-scope).

The ~204 ``MULTI_ENCOUNTER_SEQUENCE`` catalogue rows are NOT reachable through the
2-encounter idealised-optimiser sweep in :func:`cyclerfinder.data.discover.discover`.
The 12 that carry a ``free_return_arcs[]`` descriptor (the McConaghy/Russell-12
SnLm chains, now classified ``DESCRIPTOR_CLOSABLE``) ARE reachable through the #137
free-return (radial-crossing) genome, exposed as the sibling discover path
:func:`cyclerfinder.data.discover.discover_free_return`.

This module is the discover-reachability gate for that path. It reproduces the
#137 campaign result (``docs/notes/2026-06-07-russell12-freereturn-results.md``:
8 CLOSE-AND-MATCH, 3 CLOSE-MATCH-SYMMETRIC-ONLY, 1 CLOSE-OFF-ANCHOR, 0 NO-CLOSE;
4 of the closing rows promote to §14 V1) through the package API rather than the
campaign script.

GOLDEN DISCIPLINE (project memory ``feedback_golden_tests_sourced_only``): the
EXPECTED side of every V_inf match is the row's SOURCED catalogue anchor
(``vinf_kms_at_encounters``, Russell 2004 / McConaghy 2006). The emerged V_inf is
EVIDENCE, never imposed; nothing our code computes is used as an EXPECTED.
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.discover import FreeReturnDiscovery, discover_free_return

# The #137 campaign per-row outcome (recovered run, dense phase scan — the current
# default). The DESCRIPTOR_CLOSABLE bucket's expectation. Same as the
# campaign_russell12.py --genome free-return --model circular summary.
_EXPECTED_OUTCOME: dict[str, str] = {
    "mcconaghy-2006-em-k2": "CLOSE-AND-MATCH",
    "russell-ch4-4.991gG2": "CLOSE-AND-MATCH",
    "russell-ch4-8.049gGf2": "CLOSE-AND-MATCH",
    "russell-ch4-9.353Gg2": "CLOSE-AND-MATCH",
    "russell-ch4-3.64gGg3": "CLOSE-AND-MATCH",
    "russell-ch4-3.78Gg3": "CLOSE-AND-MATCH",
    "russell-ch4-5.30gGf3": "CLOSE-AND-MATCH",
    "russell-ch4-9.94Gg3": "CLOSE-AND-MATCH",
    "russell-ch4-3.66gfF3": "CLOSE-MATCH-SYMMETRIC-ONLY",
    "russell-ch4-5.30ggF3": "CLOSE-MATCH-SYMMETRIC-ONLY",
    "russell-ch4-5.75ggF3": "CLOSE-MATCH-SYMMETRIC-ONLY",
    "russell-ch4-6.44Gg3": "CLOSE-OFF-ANCHOR",
}

# The free-return rows that clear §14 V1 mechanics (#137 Part 1 + Part 3): a closed,
# V_inf-continuous reconstructed E->M->E arc. 8.049gGf2 joined after the #200/#205
# Lambert-accuracy fixes made its arc close V_inf-continuously (continuity 103 m/s,
# in-family with the 0.9-190.5 m/s accepted band; row is catalogue-V3) — see the
# V1_ROWS note in tests/search/test_free_return_v1_mechanics.py.
_EXPECTED_V1: frozenset[str] = frozenset(
    {
        "russell-ch4-5.30gGf3",
        "russell-ch4-9.94Gg3",
        "russell-ch4-5.75ggF3",
        "russell-ch4-9.353Gg2",
        "russell-ch4-8.049gGf2",
    }
)


@pytest.fixture(scope="module")
def _discoveries() -> dict[str, FreeReturnDiscovery]:
    """Run the free-return discover path once for the whole module (~7 s)."""
    return {d.row_id: d for d in discover_free_return(ephem=Ephemeris("circular"))}


@pytest.mark.slow
def test_discover_free_return_reaches_every_descriptor_row(
    _discoveries: dict[str, FreeReturnDiscovery],
) -> None:
    """Every DESCRIPTOR_CLOSABLE row is reached (yielded) by the discover path.

    The bucket is no longer unreachable: discover_free_return surfaces all 12
    free_return_arcs[]-bearing rows. N-agnostic — keyed on descriptor presence.
    """
    assert set(_discoveries) == set(_EXPECTED_OUTCOME), (
        f"discover path did not reach exactly the descriptor rows: "
        f"missing={set(_EXPECTED_OUTCOME) - set(_discoveries)}, "
        f"extra={set(_discoveries) - set(_EXPECTED_OUTCOME)}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("rid", sorted(_EXPECTED_OUTCOME))
def test_discover_free_return_outcome_matches_campaign(
    rid: str, _discoveries: dict[str, FreeReturnDiscovery]
) -> None:
    """Per-row outcome reproduces the #137 campaign (closes-and-matches the SOURCED
    V_inf anchor for the known-good rows; honest off-anchor / symmetric-only else).
    """
    d = _discoveries[rid]
    assert d.outcome == _EXPECTED_OUTCOME[rid], (
        f"{rid}: outcome {d.outcome!r} != expected {_EXPECTED_OUTCOME[rid]!r}; "
        f"derived V_inf={d.derived_vinf_kms} vs sourced={d.sourced_vinf_kms}"
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "rid",
    [
        r
        for r, o in _EXPECTED_OUTCOME.items()
        if o in ("CLOSE-AND-MATCH", "CLOSE-MATCH-SYMMETRIC-ONLY")
    ],
)
def test_matched_rows_emerged_vinf_within_tolerance_of_sourced(
    rid: str, _discoveries: dict[str, FreeReturnDiscovery]
) -> None:
    """For every closing-and-matching row the EMERGED per-body V_inf lands within
    the 0.5 km/s campaign tolerance of the INDEPENDENTLY SOURCED anchor.

    This is the golden assertion: EXPECTED = the row's sourced ``vinf_kms``; the
    emerged value is never imposed.
    """
    d = _discoveries[rid]
    assert d.closed
    assert d.sourced_vinf_kms, f"{rid}: no sourced V_inf anchor to compare against"
    for body, src in d.sourced_vinf_kms.items():
        ach = d.derived_vinf_kms.get(body)
        assert ach is not None, f"{rid}: no emerged V_inf for body {body}"
        assert abs(ach - src) <= 0.5, (
            f"{rid}: emerged V_inf[{body}]={ach:.3f} off sourced {src:.3f} (> 0.5 km/s)"
        )


@pytest.mark.slow
def test_discover_free_return_v1_promotions_match_campaign(
    _discoveries: dict[str, FreeReturnDiscovery],
) -> None:
    """Exactly the four #137 free-return rows clear §14 V1 mechanics on a closed,
    V_inf-continuous reconstructed arc — like-for-like circular, NOT a V3 result."""
    promoted = {rid for rid, d in _discoveries.items() if d.v1_passed}
    assert promoted == _EXPECTED_V1, (
        f"V1 promotions diverged from the #137 campaign: "
        f"missing={_EXPECTED_V1 - promoted}, extra={promoted - _EXPECTED_V1}"
    )
