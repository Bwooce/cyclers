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


# ---------------------------------------------------------------------------
# #480 Task 4 — corrector against real Galilean ephemeris at the paper epoch.
#
# CHARACTERISED RESULT (math decided — see scripts/_ieg_epoch_scan_480.py and the
# Task-4 report). Two independent honest negatives were found:
#
#   1. CORRECTOR GAP. ``nbody.shooter.shoot`` integrates a HELIOCENTRIC central mass
#      (MU_SUN, ~1047x Jupiter's GM) with planet perturbers from ``constants.PLANETS``.
#      The IEG seed is Jupiter-centred with Galilean-moon nodes (not in PLANETS), so
#      ``shoot()`` cannot run on it: the REBOUND callback raises KeyError (swallowed),
#      the wrong central GM makes IAS15 step-collapse, and the call does not terminate
#      in a usable budget. The Jupiter-correct propagator (``nbody.jovian``) is NOT
#      wired into ``shoot()``. The prescribed ``shoot(...)`` corrector test is therefore
#      NOT runnable here; we assert the characterised result with the Jovian-correct
#      propagator instead (the only physically valid path available).
#
#   2. BASIN MISS. Even with the correct propagator, the single-revolution Lambert
#      EGGIE seed at the paper epoch sits FAR from the basin: per-epoch seed defect
#      ~2.3e6 (best in a +/-2-synodic window ~1.6e6 at ET 6.55335e8, +4.9 d), seed
#      V_inf 3-27 km/s vs the paper's 9.12/7.07/7.07/8.38. The construction does not
#      reach the EGGIE basin anywhere in the launch window.
#
# These tests assert STRUCTURE + the characterised negative only; they do NOT assert
# reproduction pass/fail (Task 6's golden).
# ---------------------------------------------------------------------------

# Best departure ET from the +/-2-synodic epoch scan (scripts/_ieg_epoch_scan_480.py;
# runlog /tmp/ieg_epoch_scan_480.jsonl). +4.9 d past the 02-Oct-2020 paper epoch.
BEST_ET: float = 655335360.0
PAPER_ET: float = 654912000.0  # 2020-OCT-02 12:00 TDB.


def test_ieg_paper_departure_et_sanity() -> None:
    """Paper-epoch ET resolves to ~6.55e8 s (NOT the historical 10x-wrong ~6.6e7)."""
    from cyclerfinder.search.ieg_seed import PAPER_DEPARTURE_ET_APPROX, paper_departure_et

    et = paper_departure_et()
    # Within one day of the hardcoded 02-Oct-2020 12:00 TDB value used by the scan.
    assert abs(et - PAPER_ET) < 86400.0, f"paper ET {et} not near {PAPER_ET}"
    # The order-of-magnitude guard: ~6.5e8, never ~6.6e7.
    assert 6.0e8 < et < 7.0e8, f"paper ET {et} 10x off (expected ~{PAPER_DEPARTURE_ET_APPROX:.2e})"


def test_ieg_multirev_seed_vinf_near_paper() -> None:
    """Multi-rev EGGIE seed V∞ improvement over single-rev (Task 4b characterised result).

    Sourced targets (Hernandez 2017, AAS 17-608, Table 4):
      Node 0 Europa departure: 9.12 km/s
      Node 1 Ganymede arrival: 7.07 km/s
      Node 2 Ganymede arrival: 7.07 km/s   <- IMPROVED by multi-rev (leg 1 n=1,high)
      Node 3 Io arrival:       8.38 km/s
      Node 4 Europa arrival:   9.12 km/s   <- IMPROVED by multi-rev (leg 3 n=1,high)

    CHARACTERISED FINDING at departure_et=0:
    - Node 2 (Ganymede, from leg 1 n=1,high): 7.04 km/s vs 7.07 target — WITHIN 2.5 km/s.
    - Node 4 (Europa,   from leg 3 n=1,high): 8.73 km/s vs 9.12 target — WITHIN 2.5 km/s.
    - Nodes 0,1: leg 0 (1.59 d) has no feasible multi-rev arc; departure_et=0 geometry
      gives ~5.4/2.4 km/s vs 9.12/7.07.  Gap is geometric, not a solver gap.
    - Node 3 (Io):  leg 2 Io arrival ~20 km/s vs 8.38 target — best multi-rev n=2,high
      still cannot match; the Keplerian Lambert arc cannot reproduce the paper's
      multi-body patched-conic encounter at this geometry/epoch.

    This test asserts the two IMPROVED nodes are within 2.5 km/s and the departure
    V∞ is also explicitly characterized (the other nodes are sanity-range checked only).
    """
    import numpy as np

    from cyclerfinder.search.ieg_seed import EGGIE_SEQUENCE, ieg_eggie_seed

    seed = ieg_eggie_seed()

    # --- Nodes that ARE improved to within 2.5 km/s by multi-rev ---
    # Node 2: Ganymede arrival from leg 1 (n=1,high): sourced 7.07 km/s.
    v_node2 = float(np.linalg.norm(seed.vinf_in[2]))
    assert abs(v_node2 - 7.07) <= 2.5, (
        f"vinf_in[2] (Ganymede, leg1 n=1,high) = {v_node2:.2f} km/s is more than 2.5 km/s "
        f"from sourced 7.07 km/s — multi-rev selection may have regressed"
    )

    # Node 4: Europa arrival from leg 3 (n=1,high): sourced 9.12 km/s.
    v_node4 = float(np.linalg.norm(seed.vinf_in[4]))
    assert abs(v_node4 - 9.12) <= 2.5, (
        f"vinf_in[4] (Europa, leg3 n=1,high) = {v_node4:.2f} km/s is more than 2.5 km/s "
        f"from sourced 9.12 km/s — multi-rev selection may have regressed"
    )

    # --- CHARACTERIZED GAP NODES: single-rev only possible / Keplerian geometry gap ---
    # These are asserted to be in a PHYSICALLY PLAUSIBLE range only (not near the paper target).
    # Node 0 (Europa departure, from leg 0 n=0 single): not within 2.5 of 9.12 at ET=0.
    v_dep = float(np.linalg.norm(seed.vinf_out[0]))
    assert 1.0 <= v_dep <= 30.0, (
        f"vinf_out[0] (Europa departure) {v_dep:.2f} km/s outside 1-30 km/s sanity range"
    )

    # Node 1 (Ganymede arrival from leg 0 single): not within 2.5 of 7.07 at ET=0.
    v_node1 = float(np.linalg.norm(seed.vinf_in[1]))
    assert 1.0 <= v_node1 <= 30.0, (
        f"vinf_in[1] (Ganymede, leg0 arrival) {v_node1:.2f} km/s outside 1-30 km/s sanity range"
    )

    # Node 3 (Io arrival from leg 2 n=2,high): characterized ~20 km/s vs 8.38 target.
    v_node3 = float(np.linalg.norm(seed.vinf_in[3]))
    assert 1.0 <= v_node3 <= 30.0, (
        f"vinf_in[3] (Io arrival) {v_node3:.2f} km/s outside 1-30 km/s sanity range"
    )

    # Sanity: all vectors are finite.
    for i in range(5):
        assert all(np.isfinite(x) for x in seed.vinf_in[i]), f"vinf_in[{i}] not finite"
        assert all(np.isfinite(x) for x in seed.vinf_out[i]), f"vinf_out[{i}] not finite"

    # Report the full V∞ picture as contextual info (not asserted as golden).
    _ = EGGIE_SEQUENCE  # suppress unused import warning


def test_ieg_multirev_seed_defect_improves() -> None:
    """Multi-rev seed defect after jovian_shoot is substantially below the single-rev 2.47e6.

    Fast structural assertion: just build the multi-rev seed and compute the seed_defect_norm
    WITHOUT running the corrector iterations.  Checks that the multi-rev V∞ give a seed
    substantially closer to the EGGIE basin than the single-rev ~2.47e6.
    """
    import math

    import numpy as np

    from cyclerfinder.nbody.jovian import (
        JovianEphemeris,
        JovianRailsCache,
        JovianRestrictedNBody,
    )
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed
    from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

    seed = ieg_eggie_seed(departure_et=0.0)
    moons = ("Io", "Europa", "Ganymede")
    jeph = JovianEphemeris(ensure_jup365_kernel())
    cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
    prop = JovianRestrictedNBody()

    res: list[float] = []
    n = len(seed.sequence)
    for i in range(n - 1):
        s_i = seed.node_states[i]
        arc = prop.propagate(
            np.asarray(s_i[:3], dtype=np.float64),
            np.asarray(s_i[3:], dtype=np.float64),
            t0_sec=seed.epochs[i],
            t1_sec=seed.epochs[i + 1],
            moons=moons,
            cache=cache,
            max_wall_sec=30.0,
        )
        assert arc.converged, f"leg {i} propagation did not converge"
        assert np.all(np.isfinite(arc.r_km)) and np.all(np.isfinite(arc.v_km_s))
        s_next = seed.node_states[i + 1]
        res.extend(float(x) for x in (arc.r_km - np.asarray(s_next[:3])))
        res.extend(float(x) for x in (arc.v_km_s - np.asarray(s_next[3:])))

    seed_defect_norm = float(np.linalg.norm(res))
    assert math.isfinite(seed_defect_norm)
    # CHARACTERISED IMPROVEMENT: multi-rev seed reduces defect from ~2.47e6 to ~1.58e6
    # (improvement ~1.57x).  The dominant residual is Leg 2 Io-arrival (~20 km/s vs
    # 8.38 km/s paper target) which no multi-rev Keplerian arc can resolve at ET=0.
    # The threshold <2.0e6 confirms the multi-rev legs are materially closer to the basin
    # than single-rev without asserting a basin that geometry cannot reach.
    assert seed_defect_norm < 2.0e6, (
        f"Multi-rev seed defect {seed_defect_norm:.3e} is NOT an improvement over "
        f"the single-rev ~2.47e6 baseline; multi-rev topology may have regressed."
    )
    # Also guard against accidental over-improvement (a sign something is wrong with the test).
    assert seed_defect_norm > 1.0e4, (
        f"Multi-rev seed defect {seed_defect_norm:.3e} is unexpectedly tiny; "
        "if the seed now closes, revisit the Task-4b characterised finding."
    )


def test_ieg_corrects_against_real_ephemeris_at_paper_epoch() -> None:
    """Characterised #480 Task-4 result against the Jovian-correct propagator.

    The prescribed ``nbody.shooter.shoot`` corrector cannot run on the Jupiter-centred
    IEG seed (heliocentric central mass + PLANETS-only perturber registry, see the
    block comment above). We therefore exercise the only valid propagator
    (``JovianRestrictedNBody``, central Jupiter, moons on rails) and assert STRUCTURE +
    finiteness of the best-epoch seed-continuity residual — math decides, no
    reproduction pass/fail asserted.
    """
    import math

    import numpy as np

    from cyclerfinder.nbody.jovian import (
        JovianEphemeris,
        JovianRailsCache,
        JovianRestrictedNBody,
    )
    from cyclerfinder.search.ieg_seed import EGGIE_SEQUENCE, ieg_eggie_seed
    from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

    seed = ieg_eggie_seed(departure_et=BEST_ET)
    assert seed.sequence == ("Europa", "Ganymede", "Ganymede", "Io", "Europa")
    assert seed.sequence == EGGIE_SEQUENCE
    assert len(seed.node_states) == 5

    moons = ("Io", "Europa", "Ganymede")
    jeph = JovianEphemeris(ensure_jup365_kernel())
    cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
    prop = JovianRestrictedNBody()

    res: list[float] = []
    n = len(seed.sequence)
    for i in range(n - 1):
        s_i = seed.node_states[i]
        arc = prop.propagate(
            np.asarray(s_i[:3], dtype=np.float64),
            np.asarray(s_i[3:], dtype=np.float64),
            t0_sec=seed.epochs[i],
            t1_sec=seed.epochs[i + 1],
            moons=moons,
            cache=cache,
            max_wall_sec=30.0,
        )
        assert arc.converged, f"leg {i} propagation did not converge"
        assert np.all(np.isfinite(arc.r_km)) and np.all(np.isfinite(arc.v_km_s))
        s_next = seed.node_states[i + 1]
        res.extend(float(x) for x in (arc.r_km - np.asarray(s_next[:3])))
        res.extend(float(x) for x in (arc.v_km_s - np.asarray(s_next[3:])))

    seed_defect_norm = float(np.linalg.norm(res))
    assert math.isfinite(seed_defect_norm)
    # CHARACTERISED NEGATIVE: the Lambert single-rev seed is far from the EGGIE basin
    # even at the best epoch (residual ~1.6e6 km/(km/s)). This is the honest result,
    # NOT a reproduction pass. If a future seed/corrector lands the basin this loosens.
    assert seed_defect_norm > 1.0e5, (
        f"seed_defect_norm={seed_defect_norm:.3e} unexpectedly small — if the seed now "
        "reaches the basin, revisit the Task-4 characterised negative."
    )
