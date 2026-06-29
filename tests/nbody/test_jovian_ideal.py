"""Fast tests for the ideal circular-coplanar Galilean model (#480 Stage 2).

These are the fast unit gates on the λ=0 homotopy endpoint: the ideal ephemeris is a
true circle, its moons sit exactly on the Stage-1 conic guess's node positions at the
guess epochs, and the derived :class:`ShootingSeed` carries V∞ magnitudes on the
sourced Table-4 targets. The expensive n-body close (``ideal_eggie_shoot``) is NOT in
the default suite — it lives in the Stage-2 runner / a slow-marked verdict test.

V∞ golden values are SOURCED to Hernandez-Jones-Jesick 2017 (AAS 17-608) Table 4 via
docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
(``feedback_golden_tests_sourced_only`` — never a value our own code computed).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2
from cyclerfinder.nbody.jovian_ideal import (
    EGGIE_MOONS,
    IdealJovianEphemeris,
    build_eggie_shooting_seed,
    ideal_ephemeris_from_guess,
)
from cyclerfinder.search.resonant_conic import (
    EGGIE_VINF_TARGET_KMS,
    eggie_initial_guess,
    eggie_refined_guess,
    ideal_moon_smas,
)

# SOURCED — Table 4, Hernandez 2017 (AAS 17-608); the make-or-break energy targets.
TABLE4_VINF_KMS = {"Europa": 9.12, "Ganymede": 7.07, "Io": 8.38}


def test_ideal_ephemeris_is_a_circle() -> None:
    """Each moon's state is a prograde circle: |r|=a, |v|=sqrt(mu/a), r.v=0."""
    smas = ideal_moon_smas()
    phases = {m: 0.3 * i for i, m in enumerate(EGGIE_MOONS)}
    anchors = {m: 0.0 for m in EGGIE_MOONS}
    eph = IdealJovianEphemeris(smas, phases, anchors)
    for moon in EGGIE_MOONS:
        for t in (0.0, 1.0e5, 5.0e5):
            r, v = eph.state(moon, t)
            assert r[2] == 0.0 and v[2] == 0.0
            assert float(np.linalg.norm(r)) == pytest.approx(smas[moon], rel=1e-12)
            assert float(np.linalg.norm(v)) == pytest.approx(
                math.sqrt(MU_JUPITER_KM3_S2 / smas[moon]), rel=1e-12
            )
            assert abs(float(np.dot(r, v))) / smas[moon] < 1e-6  # circular: r ⟂ v


def test_ideal_moons_sit_on_guess_nodes() -> None:
    """The corrector must integrate against moons that ARE at the seed node positions.

    At every encounter epoch the ideal moon position/velocity must coincide with the
    Stage-1 guess's node position / moon velocity — otherwise the seed defect would be
    a placement artefact, not a dynamics signal.
    """
    guess = eggie_initial_guess(e=0.620)
    eph = ideal_ephemeris_from_guess(guess)
    for k, moon in enumerate(guess.sequence):
        epoch = guess.epochs_days[k] * SECONDS_PER_DAY
        r, v = eph.state(moon, epoch)
        np.testing.assert_allclose(r, guess.node_positions[k], rtol=0, atol=1e-6)
        np.testing.assert_allclose(v, guess.moon_velocities[k], rtol=0, atol=1e-9)


def test_seed_builder_matches_guess_nodes() -> None:
    """The ShootingSeed node states echo the guess (position | spacecraft velocity)."""
    guess = eggie_initial_guess(e=0.620)
    seed = build_eggie_shooting_seed(guess)
    assert seed.sequence == guess.sequence
    assert len(seed.node_states) == len(guess.sequence)
    for k in range(len(guess.sequence)):
        np.testing.assert_allclose(seed.node_states[k][:3], guess.node_positions[k], atol=1e-9)
        np.testing.assert_allclose(seed.node_states[k][3:], guess.node_sc_velocities[k], atol=1e-9)
        assert seed.epochs[k] == guess.epochs_days[k] * SECONDS_PER_DAY


def test_refined_seed_vinf_on_table4_targets() -> None:
    """The REFINED seed's per-node V∞ all land on Table 4 (±0.5 km/s).

    The raw conic guess leaves Ganymede V∞ ~4 km/s off; the in-basin refine
    (V∞ pinned to the band) is the Stage-2 family-correct seed. Computed against the
    SAME ideal ephemeris the corrector uses; the EXPECTED side is the sourced Table-4
    multiset (±0.5 km/s, the digest's reproduction band).
    """
    guess = eggie_refined_guess()
    eph = ideal_ephemeris_from_guess(guess)
    seed = build_eggie_shooting_seed(guess)
    # Node 0 departs Europa (outgoing leg velocity); nodes 1..4 are arrivals.
    checks = [(0, "Europa")] + [(k, guess.sequence[k]) for k in range(1, len(guess.sequence))]
    for k, moon in checks:
        _, v_m = eph.state(moon, seed.epochs[k])
        vinf = float(np.linalg.norm(seed.node_states[k][3:] - v_m))
        assert abs(vinf - TABLE4_VINF_KMS[moon]) <= 0.5, f"{moon} node {k}: V∞={vinf:.3f}"
    # Sanity: the targets we assert against are exactly the module's sourced dict.
    assert EGGIE_VINF_TARGET_KMS == TABLE4_VINF_KMS


def test_refined_ideal_moons_sit_on_seed_nodes() -> None:
    """The ideal moons coincide with the refined seed's node positions at each epoch."""
    guess = eggie_refined_guess()
    eph = ideal_ephemeris_from_guess(guess)
    seed = build_eggie_shooting_seed(guess)
    for k, moon in enumerate(guess.sequence):
        r, _ = eph.state(moon, seed.epochs[k])
        np.testing.assert_allclose(r, seed.node_states[k][:3], rtol=0, atol=1e-6)
