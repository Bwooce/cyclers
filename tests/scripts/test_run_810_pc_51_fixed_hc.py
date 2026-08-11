"""#810: evidence pins for the Pluto-Charon (5,1) fixed-hc sweep discovery.

The #810 run (scripts/run_810_pc_51_fixed_hc_sweep.py) held the half-crossing
index FIXED through the C-sweep (#504's (3,2) positive-control convention,
closing #656's hc=None branch-loss method gap) and found the prograde (5,1)
family's ONE stable window at Pluto-Charon mu, with the nu=0 midpoint at
C=3.167935964707279 passing every gate (winding topology, Barden |nu|<1,
independent Radau crosscheck, #660 body clearance).

These tests are the INDEPENDENT-REPRODUCTION evidence for that candidate:
each rebuilds its claim from scratch through the library primitives (fresh
corrector call, fresh winding classification, fresh Radau crosscheck, fresh
clearance propagation), pinned against the run's recorded values. They are
deliberately NOT marked slow — an evidence test CI never runs is an
unverified claim.

The expected-side numbers are this project's own #810 run record (a
discovery-evidence regression pin, like test_656's grid-seed pin), NOT
sourced golden values — the (5,1) family at PC mu appears in no publication
(Ross & Roberts-Tsoukkas 2026, arXiv:2606.29189, tabulates only
(1,1),(3,1),(3,2),(3,3) and never addresses mu~0.1087; grounded against the
source 2026-08-11).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import cyclerfinder.search.cr3bp_periodic as cp
import scripts.run_810_pc_51_fixed_hc_sweep as run
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.pluto_charon_kk_sweep import (
    CHARON_RADIUS_KM,
    PLUTO_RADIUS_KM,
    make_pluto_charon_system,
)
from cyclerfinder.search.real_binary_kk_sweep import min_body_clearance_km

# The #810 run's nu=0 stable midpoint (data/found/810_pc51_fixed_hc_sweep/
# sweep_up.jsonl, nu_zero event; see the module docstring for provenance).
MID_C = 3.167935964707279
MID_X0 = -0.7058054139668293
MID_T = 24.305715846921732
MID_HC = 5


def test_656_seed_is_a_genuine_prograde_51_orbit_with_hc5() -> None:
    """#656's recorded seed reconverges and is a prograde (5,1); its own
    perpendicular-crossing index (the one #810 holds fixed) is 5, not the
    #656 bullet's prose 'hc=4' (that was the lost (4,0) branch's index)."""
    sys_pc = make_pluto_charon_system()
    o = cp.correct_symmetric_fixed_jacobi(
        sys_pc,
        run.SEED_X0,
        run.SEED_C,
        run.SEED_T,
        ydot0_sign=-1.0,
        half_crossings=None,
        tol=1e-10,
    )
    assert o.converged
    state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
    topo = winding_topology(sys_pc.mu, state0, o.period)
    assert (topo.k1, topo.k2) == (5, 1)
    assert topo.prograde and topo.reaches_secondary
    hc, _n, _times = run.measure_half_crossing_index(sys_pc, o)
    assert hc == MID_HC


def test_stable_51_member_reconverges_standalone_with_all_gates() -> None:
    """The discovered stable member, rebuilt from scratch: fresh fixed-hc
    corrector, fresh winding topology, Barden |nu|<1, independent Radau
    crosscheck, and the #660 body-clearance gate."""
    sys_pc = make_pluto_charon_system()
    o = cp.correct_symmetric_fixed_jacobi(
        sys_pc,
        MID_X0,
        MID_C,
        MID_T,
        ydot0_sign=-1.0,
        half_crossings=MID_HC,
        tol=1e-11,
    )
    assert o.converged
    assert abs(o.x0 - MID_X0) < 1e-8
    assert abs(o.period - MID_T) < 1e-7
    assert abs(o.jacobi - MID_C) < 1e-10

    state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
    topo = winding_topology(sys_pc.mu, state0, o.period)
    assert (topo.k1, topo.k2) == (5, 1)
    assert topo.prograde and topo.reaches_secondary

    nu, _lam = cp.barden_stability(sys_pc, o, rtol=1e-13, atol=1e-13)
    assert abs(nu) < 1e-6  # run recorded 6.0e-10; the member is the nu=0 midpoint

    po = cp.PeriodicOrbit(
        state0=state0,
        period=o.period,
        jacobi=o.jacobi,
        converged=o.converged,
        closure_residual=o.crossing_residual,
    )
    ok_cc, dj = cp.crosscheck_periodic(sys_pc, po, closure_tol=1e-6, jacobi_tol=1e-8)
    assert ok_cc
    assert dj < 1e-10

    d_p, d_s = min_body_clearance_km(sys_pc, o.x0, o.ydot0, o.period)
    # Run recorded 5128.5 / 745.5 km; Charon clearance is only ~140 km above
    # the surface -- narrow but PASSING (the #659/#660 Antiope lesson is what
    # this assertion guards).
    assert d_p > PLUTO_RADIUS_KM
    assert d_s > CHARON_RADIUS_KM
    assert abs(d_p - 5128.5) < 5.0
    assert abs(d_s - 745.5) < 5.0


def test_checkpoint_record_has_exactly_one_gate_passing_candidate() -> None:
    """The committed run record contains exactly one stable topology-verified
    candidate (the nu=0 midpoint) and it passed the clearance gate."""
    cands = run._stable_candidates()
    assert len(cands) == 1
    c = cands[0]
    assert c["source"] == "sweep_up_nu_zero_mid"
    assert abs(c["c"] - MID_C) < 1e-12
    assert c["crosscheck_ok"] is True

    gates_path = Path(run.OUT_DIR) / "gates.jsonl"
    gates = [json.loads(line) for line in gates_path.read_text().splitlines() if line.strip()]
    assert len(gates) == 1
    assert gates[0]["clearance_ok"] is True
