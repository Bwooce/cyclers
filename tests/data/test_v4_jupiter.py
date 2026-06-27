"""Tests for the IEG EGGIE Lambert-real seed adapter (#480 Task 3).

Sourced golden discipline
-------------------------
The test asserts STRUCTURE (sequence, array lengths, finite values, ToF
matching the sourced digest) and a SANITY V-inf range for the Lambert-derived
seed.  It does NOT assert that the V-inf equals the digest values (9.12 / 7.07 /
7.07 / 8.38 km/s) -- the seed is a Lambert-geometry GUESS in real ephemeris at
departure ET=0; matching the paper's actual V-inf is Task 4/6's job.

Source:
  docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
"""

from __future__ import annotations

import math


def test_ieg_eggie_seed_structure() -> None:
    """EGGIE seed has the correct sequence, node count, ToFs, and period."""
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    seed = ieg_eggie_seed()
    assert seed.sequence == ("Europa", "Ganymede", "Ganymede", "Io", "Europa")
    assert len(seed.node_states) == 5
    assert len(seed.tofs) == 4
    # ToFs must match the sourced Table 4 values to within 0.01 d.
    sourced_tofs = [1.59, 8.60, 7.34, 10.69]
    assert seed.tofs == sourced_tofs or all(
        abs(a - b) < 0.01 for a, b in zip(seed.tofs, sourced_tofs, strict=True)
    )
    assert abs(seed.period_days - 28.22) < 0.05
    # Every node state is a length-6 array of finite floats.
    for st in seed.node_states:
        assert len(st) == 6 and all(math.isfinite(x) for x in st)


def test_ieg_seed_vinf_in_reasonable_range() -> None:
    """Each vinf_in vector (after node 0, which has no inbound leg) is finite and
    within a physically plausible 1-30 km/s range for Galilean moon flybys.

    This is a SANITY check only.  The Lambert seed is a geometry guess in real
    ephemeris at departure ET=0; the actual sourced golden is 7-9 km/s (Table 4)
    and is asserted only after the corrector (Task 4/6) has converged.
    """
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    seed = ieg_eggie_seed()
    import numpy as np

    # vinf_in for nodes 1..4 are the inbound V-inf vectors.
    for i in range(1, len(seed.vinf_in)):
        v = seed.vinf_in[i]
        assert len(v) == 3
        mag = float(np.linalg.norm(v))
        assert math.isfinite(mag), f"vinf_in[{i}] magnitude is not finite: {mag}"
        assert 1.0 <= mag <= 30.0, (
            f"vinf_in[{i}] magnitude {mag:.2f} km/s is outside the 1-30 km/s sanity window; "
            "seed may be degenerate (co-located positions -> Lambert failure)"
        )


def test_ieg_seed_epochs_consistent_with_tofs() -> None:
    """Epoch differences must equal the sourced ToFs (in seconds)."""
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    seed = ieg_eggie_seed(departure_et=0.0)
    spd = 86400.0
    for i, tof in enumerate(seed.tofs):
        dt_sec = seed.epochs[i + 1] - seed.epochs[i]
        assert abs(dt_sec - tof * spd) < 1.0, (
            f"epoch gap {dt_sec:.1f} s != tof[{i}] * 86400 = {tof * spd:.1f} s"
        )


def test_ieg_seed_slack_leg_is_longest() -> None:
    """slack_leg must be the index of the longest ToF (the period pin convention)."""
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    seed = ieg_eggie_seed()
    expected_slack = seed.tofs.index(max(seed.tofs))
    assert seed.slack_leg == expected_slack, (
        f"slack_leg={seed.slack_leg} but longest ToF is leg {expected_slack}"
        f" ({max(seed.tofs):.2f} d)"
    )


def test_ieg_seed_vinf_out_length_matches() -> None:
    """vinf_out has one vector per encounter, all length-3 finite."""
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    seed = ieg_eggie_seed()
    import numpy as np

    assert len(seed.vinf_out) == 5
    for i, v in enumerate(seed.vinf_out):
        assert len(v) == 3
        mag = float(np.linalg.norm(v))
        assert math.isfinite(mag), f"vinf_out[{i}] not finite"


def test_ieg_seed_departure_et_offset() -> None:
    """departure_et argument shifts all epochs by that offset in seconds."""
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    s0 = ieg_eggie_seed(departure_et=0.0)
    offset = 7.05 * 86400.0  # one synodic period
    s1 = ieg_eggie_seed(departure_et=offset)
    for e0, e1 in zip(s0.epochs, s1.epochs, strict=True):
        assert abs((e1 - e0) - offset) < 1.0, (
            f"epoch not shifted by offset: delta={(e1 - e0):.1f}, expected {offset:.1f}"
        )
