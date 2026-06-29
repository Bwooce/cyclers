"""Tests for the per-encounter self-consistency guard (#480 repeated-encounter bug class).

The decisive test is the NEGATIVE one: the guard must FLAG the real ``resonant_conic``
EGGIE construction, whose 2nd Ganymede / 2nd Europa drift hundreds of thousands of km off
the conic crossing (the bug that masqueraded as a physics wall). See
``docs/notes/2026-06-29-480-eggie-ballistic-construction-verdict.md``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2 as MU
from cyclerfinder.search.resonant_conic import (
    EGGIE_SEQUENCE,
    _eggie_conic_nodes,
    _moon_phase_ics,
    ideal_moon_smas,
)
from cyclerfinder.search.tour_self_consistency import (
    TourSelfConsistencyError,
    assert_encounters_self_consistent,
    encounter_gaps_km,
    soi_km,
)


def test_soi_km_ganymede_sane() -> None:
    """Ganymede's SOI is ~3.2e4 km (sma·(mu/3mu_J)^(1/3))."""
    assert soi_km("Ganymede") == pytest.approx(31717.0, rel=0.02)


def test_gaps_and_pass_when_within_soi() -> None:
    """A tour whose bodies sit at (or just inside SOI of) the nodes passes."""
    bodies = ("Europa", "Ganymede", "Io")
    nodes = [np.array([1e6, 0.0, 0.0]), np.array([0.0, 1.07e6, 0.0]), np.array([4.2e5, 0.0, 0.0])]
    # Bodies offset by a small fraction of each SOI -> consistent.
    bod = [n + np.array([0.3 * soi_km(b), 0.0, 0.0]) for n, b in zip(nodes, bodies, strict=True)]
    gaps = encounter_gaps_km(nodes, bod)
    assert all(g < soi_km(b) for g, b in zip(gaps, bodies, strict=True))
    assert_encounters_self_consistent(nodes, bod, bodies)  # must NOT raise


def test_catches_resonant_conic_repeated_encounter_drift() -> None:
    """The guard FLAGS the real resonant_conic EGGIE construction (the #480 bug).

    The conic crossing positions (where the spacecraft/conic is) vs the moon positions at
    the encounter epochs (where the body actually is): the FIRST encounter of each body
    matches by construction, but the 2nd Ganymede / 2nd Europa drift ~1e5-4e5 km — far
    outside their SOI. A correct guard MUST raise.
    """
    e = 0.62
    _a, _nus, times, conic_pos, _cvel = _eggie_conic_nodes(e, 0.0, MU)
    phases, anchors = _moon_phase_ics(times, conic_pos)
    smas = ideal_moon_smas()
    n_m = {m: math.sqrt(MU / smas[m] ** 3) for m in smas}

    moon_pos = []
    for k, moon in enumerate(EGGIE_SEQUENCE):
        theta = phases[moon] + n_m[moon] * (times[k] - anchors[moon])
        moon_pos.append(smas[moon] * np.array([math.cos(theta), math.sin(theta), 0.0]))

    gaps = encounter_gaps_km(conic_pos, moon_pos)
    # First-of-each-body is consistent (~0); a repeated body drifts >> SOI.
    assert gaps[0] < soi_km("Europa")  # Europa-depart (1st) is fine
    assert max(gaps) > 1.0e5  # a repeat is hundreds of thousands of km off

    with pytest.raises(TourSelfConsistencyError, match="repeated-encounter drift"):
        assert_encounters_self_consistent(conic_pos, moon_pos, EGGIE_SEQUENCE)
