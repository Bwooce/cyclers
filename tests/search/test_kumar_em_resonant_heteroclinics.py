"""Tests for `#827`: Kumar et al. (2026) Table-5 3:1 <-> 2:1 Earth-Moon
heteroclinic digit-grade reproduction
(docs/notes/2026-08-21-827-kumar-table5-reproduction.md).

Reproduce-before-trust discipline: the connection-evidence tests below rebuild
a result from scratch through the library primitives -- a fresh `#822` node
correction (:func:`build_kumar_node`, which itself re-derives the Table-6 IC
against the print and rejects any disagreement), a fresh Newton connection
refinement seeded ONLY by the recorded phase indices (``tau_u``, ``tau_s``,
``k_u``, ``k_s``, ``branch_u``, ``branch_s`` -- never the converged crossing
state itself, per :func:`~cyclerfinder.search.vaquero_em_cycler_connections.
refine_connection`'s own signature, which ignores a seed's ``x``/``xdot``/
``ydot``), a fresh full `#822` verification battery, and a fresh perigee-
section match against Kumar's own printed Table-5 state. Nothing here
replays the committed JSON directly -- the achieved distances are pinned
against ``data/found/827_kumar_table5_reproduction/results.json``'s own
recorded values as regression targets, following the `#810`
(``tests/scripts/test_run_810_pc_51_fixed_hc.py``) precedent. None are
``@pytest.mark.slow`` -- an evidence test CI never runs is an unverified
claim (`feedback_delegation_fresh_agent_not_fork`).
"""

from __future__ import annotations

import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.kumar_em_resonant_heteroclinics as keh
import cyclerfinder.search.vaquero_em_cycler_connections as vcc


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return keh.kumar_system()


# ---------------------------------------------------------------------------
# (1) Sourced/derived constants and table structure.
# ---------------------------------------------------------------------------


def test_kumar_mu_matches_the_paper_print() -> None:
    """Section 2.1's own printed mass ratio, NOT this project's Earth-Moon
    mu -- the two differ at the ~1e-10 absolute level (module docstring)."""
    assert keh.KUMAR_MU == 1.2150584270572e-2
    base = cr3bp.cr3bp_system("Earth", "Moon")
    assert abs(keh.KUMAR_MU - base.mu) > 1e-11


def test_kumar_system_uses_kumar_mu() -> None:
    system = keh.kumar_system()
    assert system.mu == keh.KUMAR_MU


def test_table5_has_all_seven_printed_jacobi_rows_with_correct_transfer_type() -> None:
    """Type 1 (short, direct) below the paper's own C=3.09 cutoff (Section
    5.2.1), Type 2 (long, via a 5:2 intermediary) at/above it."""
    assert sorted(keh.KUMAR_TABLE5_31_TO_21) == [2.54, 2.70, 2.86, 3.00, 3.05, 3.10, 3.15]
    for c, (ttype, _state) in keh.KUMAR_TABLE5_31_TO_21.items():
        assert ttype == (1 if c < 3.09 else 2)


def test_kumar_table5_state6_lifts_the_planar_print_correctly() -> None:
    state = keh.kumar_table5_state6(2.86)
    _ttype, (x, y, xd, yd) = keh.KUMAR_TABLE5_31_TO_21[2.86]
    assert state[0] == x
    assert state[1] == y
    assert state[2] == 0.0
    assert state[3] == xd
    assert state[4] == yd
    assert state[5] == 0.0


def test_kumar_table5_state6_rejects_a_c_not_in_the_printed_table() -> None:
    """`#839`'s C=3.13 (`vaquero-31-c313-em-resonant-po-2013`) is deliberately
    NOT one of the seven printed rows -- Kumar's Table 5 never tabulates it,
    so `#839` needs its own targeted run, out of `#827`'s scope."""
    with pytest.raises(KeyError):
        keh.kumar_table5_state6(3.13)


def test_reproduction_cs_is_a_subset_of_the_printed_table() -> None:
    assert set(keh.KUMAR_REPRODUCTION_CS) <= set(keh.KUMAR_TABLE5_31_TO_21)


# ---------------------------------------------------------------------------
# (2) Node re-derivation (never a raw trust of the Table-6 print).
# ---------------------------------------------------------------------------


def test_build_kumar_node_rederives_table6_print_at_2_86(system: cr3bp.CR3BPSystem) -> None:
    node31 = keh.build_kumar_node(system, 3, 2.86)
    node21 = keh.build_kumar_node(system, 2, 2.86)
    assert node31.converged
    assert node21.converged
    x0_31, ydot0_31 = keh.KUMAR_TABLE6_ICS[(3, 2.86)]
    x0_21, ydot0_21 = keh.KUMAR_TABLE6_ICS[(2, 2.86)]
    assert abs(node31.state0[0] - x0_31) < keh.TABLE6_IC_ABS_TOL
    assert abs(node31.state0[4] - ydot0_31) < keh.TABLE6_IC_ABS_TOL
    assert abs(node21.state0[0] - x0_21) < keh.TABLE6_IC_ABS_TOL
    assert abs(node21.state0[4] - ydot0_21) < keh.TABLE6_IC_ABS_TOL


def test_build_kumar_node_rejects_c_off_the_table6_print(system: cr3bp.CR3BPSystem) -> None:
    with pytest.raises(KeyError):
        keh.build_kumar_node(system, 3, 3.13)


# ---------------------------------------------------------------------------
# (3) Independent from-scratch connection reconstruction, pinned against the
#     committed #827 run record (data/found/827_kumar_table5_reproduction/
#     results.json) as a regression target, not as a value to replay.
# ---------------------------------------------------------------------------


def test_connection_reconverges_and_matches_print_at_2_86_type1(
    system: cr3bp.CR3BPSystem,
) -> None:
    """C=2.86 (Type 1, short transfer). Seeded ONLY by the recorded phase
    indices; the Newton residual, verification battery, and match distance
    are all freshly computed here, then pinned against the committed run."""
    c = 2.86
    node31 = keh.build_kumar_node(system, 3, c)
    node21 = keh.build_kumar_node(system, 2, c)
    seed = vcc.ConnectionSeed(
        distance=0.0,
        branch_u=-1,
        branch_s=-1,
        k_u=20,
        k_s=14,
        tau_u=1.2504446611561342,
        tau_s=5.555534605411551,
        x=0.0,
        xdot=0.0,
        ydot=0.0,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=keh.KUMAR_EPSILON)
    assert conn.converged
    assert conn.residual < 1e-8
    assert abs(conn.crossing_xv[0] - (-0.11271809747808621)) < 1e-6
    assert abs(conn.crossing_xv[1] - 1.9035816166167323) < 1e-6

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=keh.KUMAR_EPSILON)
    assert ev.passed
    assert ev.ydot_signs_match
    assert ev.full_state_gap < 1e-5

    dist, _state_m, _t_m, runner_up = keh._match_against_print(
        system, ev, keh.kumar_table5_state6(c)
    )
    assert dist < keh.KUMAR_MATCH_TOL
    assert abs(dist - 1.282737565150418e-06) < 1e-8  # regression pin (results.json)
    assert runner_up > 1e-3  # a specific-point match, not a proximity coincidence


def test_connection_reconverges_and_matches_print_at_3_15_type2(
    system: cr3bp.CR3BPSystem,
) -> None:
    """C=3.15 (Type 2, long transfer via a 5:2 intermediary segment) -- the
    other printed transfer type, at the far edge of the published band."""
    c = 3.15
    node31 = keh.build_kumar_node(system, 3, c)
    node21 = keh.build_kumar_node(system, 2, c)
    seed = vcc.ConnectionSeed(
        distance=0.0,
        branch_u=1,
        branch_s=1,
        k_u=22,
        k_s=8,
        tau_u=0.06047261430256611,
        tau_s=2.714308658676767,
        x=0.0,
        xdot=0.0,
        ydot=0.0,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=keh.KUMAR_EPSILON)
    assert conn.converged
    assert conn.residual < 1e-8
    assert abs(conn.crossing_xv[0] - (-0.43008884712315554)) < 1e-6
    assert abs(conn.crossing_xv[1] - 0.7404473027937676) < 1e-6

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=keh.KUMAR_EPSILON)
    assert ev.passed
    assert ev.ydot_signs_match

    dist, _state_m, _t_m, runner_up = keh._match_against_print(
        system, ev, keh.kumar_table5_state6(c)
    )
    assert dist < keh.KUMAR_MATCH_TOL
    assert abs(dist - 4.4851711537154614e-07) < 1e-8  # regression pin (results.json)
    assert runner_up > 1e-3
