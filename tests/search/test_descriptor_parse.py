"""M-ED Phase 2: free_return_arcs[] -> per-leg tuples (plan Phase 2)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
from cyclerfinder.search.descriptor import parse_free_return_arcs


def test_s1l1_two_generic_arcs() -> None:
    arcs = [
        {
            "arc_type": "generic",
            "tof_years": 1.4612,
            "resonance": None,
            "raw_descriptor": "g(1.4612,526.02,Ll)",
        },
        {
            "arc_type": "generic",
            "tof_years": 2.8096,
            "resonance": None,
            "raw_descriptor": "G(2.8096,...)",
        },
    ]
    revs, branches, seeds = parse_free_return_arcs(arcs)
    assert revs == (0, 0)
    assert branches == ("single", "single")
    assert seeds[0] == pytest.approx(1.4612 * DAYS_PER_JULIAN_YEAR)
    assert seeds[1] == pytest.approx(2.8096 * DAYS_PER_JULIAN_YEAR)
