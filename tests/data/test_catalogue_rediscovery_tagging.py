"""M7 BINDING GATE — spec §8 M7 anchor (rediscovery tagging).

The 2-synodic Earth-Mars cycler is *well-known* (Byrnes/McConaghy/
Longuski 2002; catalogued as ``s1l1-2syn-em-cpom``). Spec §13.4 / §16.3
require that re-deriving it end-to-end be tagged ``known`` and inherit
the published attribution — never presented as a discovery.

This gate runs the M5 finder end-to-end, reduces the best result to a
canonical signature, and matches it against the seed catalogue. **The
catalogue entry is the fixture, not the tolerance** — no per-test V∞
widening, no signature-distance loosening.

Current status: **xfail (strict=False)** under the pre-existing M5
optimiser regression (task #54): ``find_cyclers(("E","M"), k_synodic=2,
vinf_cap=7.0)`` returns 0 results (the optimiser does not converge to
the 5.65/3.05 km/s Russell cycler), so there is no best result to
signature-match. Same root cause as
``tests/search/test_optimize.py::test_2syn_em_rediscovers_5_65_kms_earth``.
Flip ``strict=True`` once #54 lands.
"""

from __future__ import annotations

import pytest

from cyclerfinder.data.catalog import canonical_signature, load_catalog, match
from cyclerfinder.search.optimize import find_cyclers


@pytest.mark.xfail(
    strict=False,
    reason=(
        "Pre-existing M5 optimiser regression (task #54): "
        "find_cyclers(('E','M'), k_synodic=2, vinf_cap=7.0, seed=0) "
        "returns 0 results (optimiser does not converge to the 5.65 km/s "
        "Russell cycler), so there is no best result to signature-match. "
        "Same root cause as test_2syn_em_rediscovers_5_65_kms_earth. "
        "Flip strict=True once #54 lands. "
        "NOTE (2026-06-04): the 5.65 km/s anchor is unverified-provenance "
        "(catalogue data_gap vinf_kms_at_encounters, s1l1-2syn-em-cpom): "
        "traces only to spec.md §9; unconfirmed in Patel 2019 / McConaghy "
        "2006 / Sanchez Net 2022 — see docs/notes/s1l1-target-topology-mining.md."
    ),
)
def test_rediscovered_2syn_em_tagged_known() -> None:
    """M7 BINDING GATE — re-derived 2-synodic E-M cycler matches the catalogue."""
    results = find_cyclers(("E", "M"), k_synodic=2, vinf_cap=7.0, seed=0)
    assert results, "M5 finder returned no feasible cyclers (task #54)"
    best = results[0]
    signature = canonical_signature(best.best_cycler, model_assumption="circular-coplanar")
    match_result = match(signature, load_catalog())
    assert match_result.outcome == "known"
    assert match_result.entry is not None
    assert match_result.entry.id == "s1l1-2syn-em-cpom"
    assert match_result.entry.priority_date == "2002-08-05"
