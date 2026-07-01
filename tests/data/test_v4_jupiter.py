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

import pytest

from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

try:
    _KERNEL: str | None = ensure_jup365_kernel()
except Exception:  # jup365.bsp is local-only (~50 MB, absent in CI); the IEG seed needs it -> skip
    _KERNEL = None

pytestmark = pytest.mark.skipif(_KERNEL is None, reason="JUP365 kernel not furnished (local-only)")


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
# #480 M1 SCIENCE VERDICT — Jovian corrector against real Galilean ephemeris.
#
# Both infra fixes are now in: the Jupiter-central multiple-shooting corrector
# ``nbody.jovian.jovian_shoot`` (Task 4a) and the MULTI-REV ``ieg_eggie_seed``
# (Task 4b).  The corrector RUNS on the Jupiter-centred Galilean seed (the old
# heliocentric ``shoot()`` could not — KeyError 'Io' + wrong central GM).
#
# CHARACTERISED RESULT (math decided — scripts/_ieg_rerun_scan_480.py epoch scan +
# the corrector run; report in the #480 M1 thread).  A fine +/-1-synodic epoch scan
# at the BEST_ET below (-4.54 d from the 02-Oct-2020 paper epoch) drops the multi-rev
# seed defect to ~3.89e5 km/(km/s) — a sharp local minimum (~1.4e6 just 0.2 d away).
# From that warm start ``jovian_shoot`` cuts the defect ~5x (3.89e5 -> ~7.8e4) but does
# NOT converge: it lands in an OFF-PAPER basin (corrected V∞ ~ 11.5/8.3/8.5/4.9/8.0 km/s
# vs the paper's 9.12/7.07/7.07/8.38; flyby bends NOT all feasible) with a ~5.9 km/s
# correction ΔV — orders of magnitude above the paper's 0.70 m/s ballistic close.
# The deep solve is FD-Jacobian-bound (~20 s/eval; a single Gauss-Newton step is ~10
# min), the same compute wall recorded for the heliocentric multi-rev cyclers in MEMORY.
#
# These tests assert STRUCTURE + the corrector's defect-reduction (defect_norm <=
# seed_defect_norm); they do NOT assert reproduction pass/fail (Task 6's golden).
# ---------------------------------------------------------------------------

# Best departure ET from the +/-1-synodic fine epoch scan (scripts/_ieg_rerun_scan_480.py
# at 0.1 d, refined at 0.02 d; runlog /tmp/ieg_rerun_scan_480.jsonl).  -4.54 d from the
# 02-Oct-2020 paper epoch; multi-rev seed_defect_norm ~3.89e5 there (a sharp minimum).
BEST_ET: float = 654519744.0
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


def test_ieg_seed_defect_at_best_epoch_structure() -> None:
    """Fast (non-slow) structural check of the best-epoch multi-rev seed continuity.

    Propagates the four EGGIE legs of the multi-rev seed at ``BEST_ET`` through the
    Jovian-correct propagator and asserts STRUCTURE + a finite, characterised
    seed-continuity residual (~3.89e5 km/(km/s) — the sharp epoch-scan minimum).  This
    is the warm-start the slow corrector test below consumes; running the corrector
    itself is FD-Jacobian-bound (~10 min) and lives in the @slow test.
    """
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
    # CHARACTERISED: at the best epoch the multi-rev seed is materially closer to the
    # basin (~3.89e5 km/(km/s)) than at the paper epoch (~1.44e6) but still far from a
    # ballistic close. NOT a reproduction pass. If a future seed lands the basin this
    # loosens — math decides.
    assert 1.0e4 < seed_defect_norm < 1.0e6, (
        f"seed_defect_norm={seed_defect_norm:.3e} outside the characterised "
        "[1e4, 1e6] window at BEST_ET — re-run the epoch scan."
    )


@pytest.mark.slow
def test_ieg_jovian_shoot_at_best_epoch() -> None:
    """#480 M1 SCIENCE VERDICT: run ``jovian_shoot`` from the best-epoch multi-rev seed.

    Exercises the prescribed corrector (``nbody.jovian.jovian_shoot``, the
    Jupiter-central multiple-shooting analogue of ``nbody.shooter.shoot``) on the
    multi-rev EGGIE seed at ``BEST_ET``.  Asserts STRUCTURE (sequence, 5 nodes, finite
    defect) and that the corrector does NOT worsen the seed (defect_norm <=
    seed_defect_norm).  Does NOT assert the V∞-vs-paper reproduction tolerance — that
    is Task 6's golden.

    CHARACTERISED (math decided): from the ~3.89e5 seed the corrector cuts the defect
    ~5x (to ~7.8e4) at max_nfev=30 but does NOT converge; it lands in an off-paper basin
    (corrected V∞ ~ 11.5/8.3/8.5/4.9/8.0 km/s vs the paper 9.12/7.07/7.07/8.38; ~5.9 km/s
    correction ΔV, four orders above the paper's 0.70 m/s ballistic close).  Slow: a
    single FD-Jacobian Gauss-Newton step is ~10 min.
    """
    import numpy as np

    from cyclerfinder.nbody.jovian import jovian_shoot
    from cyclerfinder.search.ieg_seed import EGGIE_SEQUENCE, ieg_eggie_seed

    # Budget kept well under the 600 s per-test timeout: ~20 s/residual-eval, so
    # max_nfev=18 (~one FD Jacobian + a step) lands around ~350 s. The deeper
    # max_nfev=30/120 runs (recorded in the #480 M1 report) confirm the verdict is
    # stable — the corrector cuts the defect ~5x but stays in the off-paper basin.
    seed = ieg_eggie_seed(departure_et=BEST_ET)
    result = jovian_shoot(
        seed,
        moons=("Io", "Europa", "Ganymede"),
        max_nfev=18,
        max_wall_sec=45.0,
    )

    # STRUCTURE.
    assert result.sequence == EGGIE_SEQUENCE
    assert len(result.sequence) == 5
    assert len(result.vinf_per_encounter_kms) == 5
    assert math.isfinite(result.defect_norm)
    assert math.isfinite(result.seed_defect_norm)
    assert all(np.isfinite(v) for v in result.vinf_per_encounter_kms)
    assert math.isfinite(result.correction_dv_kms)

    # The corrector must not WORSEN the seed continuity.
    assert result.defect_norm <= result.seed_defect_norm, (
        f"corrector worsened the defect: {result.defect_norm:.3e} > "
        f"seed {result.seed_defect_norm:.3e}"
    )

    # CHARACTERISED NEGATIVE: it does not converge to the ballistic EGGIE basin at this
    # budget. Documented, not a reproduction pass (Task 6's golden owns V∞-vs-paper).
    assert not result.converged, (
        "jovian_shoot CONVERGED at BEST_ET — the characterised #480 M1 negative no "
        "longer holds; promote this to a reproduction assertion and revisit Task 6."
    )
