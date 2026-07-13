"""#585 scope item 2 -- regression coverage for the analyze-mode drift-detection check.

``scripts/run_582_asymmetric_3d_niching_search.py``'s ``--mode analyze`` now
flags a converged cluster representative as ``drifted_to_neighboring_mmr`` if
its OWN converged period implies a semi-major axis nearer a DIFFERENT
tabulated interior MMR than the one this run targeted -- widening
``mmr_bounds``'s symmetry-breaking bounds (via #585's new ``s`` parameter)
risks exactly this failure mode (#440's own documented neighboring-MMR
family-selection trap). These are property/regression tests on the pure
helper logic (:func:`nearest_mmr_by_implied_a1`, :func:`_s_tag`), not a
sourced-discovery golden.
"""

from __future__ import annotations

import pytest

import scripts.run_582_asymmetric_3d_niching_search as run
from cyclerfinder.search.er3bp_isolated_seeds import MMR_SEMI_MAJOR_AXES
from cyclerfinder.search.isolated_3d_asymmetric_fitness import mmr_t0


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_nearest_mmr_by_implied_a1_recovers_own_mmr_at_exact_t0(p: int, q: int, a1: float) -> None:
    """A period exactly at the MMR's own T0 must resolve back to that same MMR."""
    label, a1_implied, dist = run.nearest_mmr_by_implied_a1(mmr_t0(a1))
    assert label == f"{p}:{q}"
    assert a1_implied == pytest.approx(a1, rel=1e-9)
    assert dist == pytest.approx(0.0, abs=1e-9)


def test_nearest_mmr_by_implied_a1_flags_drift_to_neighbor() -> None:
    """A period roughly midway between two tabulated MMRs' T0 resolves to
    whichever tabulated a1 is numerically closest -- the exact mechanism
    ``--mode analyze`` uses to flag ``drifted_to_neighboring_mmr``."""
    t0_32 = mmr_t0(0.763143)  # 3:2
    t0_52 = mmr_t0(0.5428)  # 5:2
    # a period much closer to 5:2's own T0 than to 3:2's should resolve to 5:2.
    period_near_52 = t0_52 + 0.05 * (t0_32 - t0_52)
    label, _a1_implied, _dist = run.nearest_mmr_by_implied_a1(period_near_52)
    assert label == "5:2"


def test_s_tag_is_filesystem_safe_and_distinguishes_rungs() -> None:
    tag_1 = run._s_tag(0.15)
    tag_2 = run._s_tag(0.30)
    assert tag_1 != tag_2
    for tag in (tag_1, tag_2):
        assert "." not in tag
        assert "/" not in tag
        assert tag.startswith("s")
