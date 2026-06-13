"""JPL periodic_orbits client — parser tests (no network) + mu reconciliation.

The parser is exercised against a saved 4-row Earth-Moon L1-halo fixture
captured from the live API (tests/verify/fixtures/). The live ``query`` is NOT
tested here (no network in CI); the fixture IS a real response, so the parser
contract is pinned against ground truth without a network dependency.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.verify.jpl_periodic_orbits import (
    OUR_EM_MU,
    parse_payload,
    reconcile_mu,
)

_FIXTURE = Path(__file__).parent / "fixtures" / "jpl_earth_moon_l1_halo_sample.json"


@pytest.fixture
def payload() -> dict[str, Any]:
    data: dict[str, Any] = json.loads(_FIXTURE.read_text())
    return data


def test_parse_system_constants(payload: dict[str, Any]) -> None:
    constants, _ = parse_payload(payload)
    assert constants.name == "Earth-Moon"
    # JPL's reported Earth-Moon mass ratio (the value to reconcile against ours).
    assert constants.mu == pytest.approx(1.215058560962404e-2, rel=0, abs=1e-15)
    assert constants.lunit_km > 0
    assert constants.tunit_s > 0
    # L1/L2/L3 lie on the x-axis (y=z=0); L1 is between the primaries (0<x<1).
    assert "L1" in constants.libration_points
    l1 = constants.libration_points["L1"]
    assert 0.0 < l1[0] < 1.0
    assert l1[1] == 0.0 and l1[2] == 0.0


def test_parse_orbits_shape_and_fields(payload: dict[str, Any]) -> None:
    _, orbits = parse_payload(payload)
    assert len(orbits) == 4
    for o in orbits:
        assert o.state0.shape == (6,)
        assert o.state0.dtype == np.float64
        assert np.isfinite(o.state0).all()
        assert o.period > 0
        assert np.isfinite(o.jacobi)
        assert np.isfinite(o.stability)


def test_parse_is_field_name_keyed_not_positional(payload: dict[str, Any]) -> None:
    """A column re-order in the payload must not silently mis-map values.

    Shuffle ``fields`` (and each data row in lockstep) and confirm the parser
    recovers the SAME orbits — proving it keys by name, not position.
    """
    _, orbits_ref = parse_payload(payload)

    fields = list(payload["fields"])
    perm = [fields.index(c) for c in reversed(fields)]
    shuffled = {
        **payload,
        "fields": [fields[i] for i in perm],
        "data": [[row[i] for i in perm] for row in payload["data"]],
    }
    _, orbits_shuf = parse_payload(shuffled)

    assert len(orbits_shuf) == len(orbits_ref)
    for a, b in zip(orbits_ref, orbits_shuf, strict=True):
        np.testing.assert_array_equal(a.state0, b.state0)
        assert a.jacobi == b.jacobi
        assert a.period == b.period
        assert a.stability == b.stability


def test_parse_raises_on_missing_field(payload: dict[str, Any]) -> None:
    broken = {**payload, "fields": [f for f in payload["fields"] if f != "jacobi"]}
    with pytest.raises(ValueError, match="missing expected fields"):
        parse_payload(broken)


def test_reconcile_mu_reports_the_gap() -> None:
    # The convention gap that gates any JPL-vs-ours cross-check: ~1e-7 relative.
    rec = reconcile_mu(1.215058560962404e-2)
    assert rec["our_mu"] == OUR_EM_MU
    assert rec["abs_diff"] == pytest.approx(1.215e-9, abs=5e-10)
    assert 1e-8 < rec["rel_diff"] < 1e-6
