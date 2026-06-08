"""Tests for the GMAT V4 script generator (#171).

String-templating only — NO GMAT is installed or invoked here (the GMAT run is the
manual, out-of-CI step). The sourced inputs are each row's published v_inf nodes
(Aldrin documented Mars return geometry; S1L1 App-C ``APPC_LEGS``) and the Jones
continuity tolerance; the reference ΔVs are OUR values under external check,
asserted only as present-in-text, never as an EXPECTED-from-source equality on a
computed quantity.

The generator emits one **runnable** GMAT script per Mars flyby (Mars-relative
B-plane targeting, the ``Ex_MarsBPlane`` pattern; Beeson BC at flyby periapse).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.search.s1l1_corrected import APPC_LEGS
from scripts.gmat_v4_generate import (
    JONES_VEL_TOL_KMS,
    CyclerRow,
    aldrin_row,
    generate_aldrin_script,
    generate_s1l1_scripts,
    generate_script,
    s1l1_row,
)

# --- Phase 1: Aldrin powered-periodic generator ------------------------------


def test_aldrin_script_force_model_and_epoch() -> None:
    text = generate_script(aldrin_row(epoch_iso="01 Jan 2030 00:00:00.000"))
    assert "Create ForceModel" in text
    assert "Mars" in text
    assert "Sun" in text
    assert "01 Jan 2030 00:00:00.000" in text


def test_aldrin_script_has_flyby_target_block() -> None:
    text = generate_script(aldrin_row())
    # One Target block, BC at Mars periapse; the Jones B-plane goal is recorded as
    # the approach geometry (comment), the Achieve targets the outgoing v_inf^+.
    assert text.count("Target DC") == 1
    assert "Jones B-plane approach geometry" in text
    assert "B.R =" in text and "B.T =" in text
    assert "{Sat.Mars.Periapsis}" in text
    assert "Achieve DC(Sat.VX" in text
    assert "Achieve DC(Sat.VY" in text
    # The maintenance TCM is varied (the impulse components).
    assert "TCM_MarsReturn.Element1" in text
    assert "Maneuver TCM_MarsReturn(Sat)" in text
    # Jones SOURCED velocity tolerance written into the Achieve goal.
    assert f"Tolerance = {JONES_VEL_TOL_KMS:g}" in text


def test_aldrin_script_provides_initial_guess() -> None:
    text = generate_script(aldrin_row())
    # The seed is a real Mars-relative Cartesian state (not an empty slot).
    assert "GMAT Sat.X =" in text
    assert "GMAT Sat.VX =" in text
    assert "GMAT Sat.CoordinateSystem = MarsInertial" in text


def test_aldrin_reference_dv_recorded_as_our_value() -> None:
    text = generate_script(aldrin_row())
    assert "2.9138" in text
    assert "OUR value" in text  # honesty: reference under external check


def test_aldrin_script_is_runnable_shape() -> None:
    """Sanity: GMAT requires all objects before BeginMissionSequence."""
    text = generate_script(aldrin_row())
    begin = text.index("BeginMissionSequence")
    assert "Create Spacecraft Sat" in text[:begin]
    assert "Create DifferentialCorrector DC" in text[:begin]
    assert "Create ImpulsiveBurn TCM_MarsReturn" in text[:begin]
    # The Target/Maneuver/Achieve live AFTER BeginMissionSequence.
    assert "Target DC" in text[begin:]


# --- Phase 2: S1L1 flyby-station-keep chain ----------------------------------


def test_s1l1_row_seeds_appc_leg2() -> None:
    row = s1l1_row()
    # The seed v_inf must be the App-C leg-2 value read from APPC_LEGS (not a literal).
    leg2 = next(leg for leg in APPC_LEGS if leg[0] == 2)
    expected = np.array(leg2[3], dtype=np.float64)
    np.testing.assert_allclose(row.seed_vinf_kms, expected, rtol=1e-12)


def test_s1l1_has_per_mars_flyby_targets() -> None:
    row = s1l1_row()
    n_mars_nodes = sum(1 for leg in APPC_LEGS if leg[1] == "M")
    assert n_mars_nodes == 7  # App-C has 7 Mars encounters
    assert len(row.mars_flybys) == n_mars_nodes
    # Each flyby renders to its own runnable script with one Target block targeting
    # that node's outgoing v_inf^+, BC at Mars periapse.
    for i in range(n_mars_nodes):
        text = generate_script(row, node_index=i)
        assert text.count("Target DC") == 1
        assert "Jones B-plane approach geometry" in text
        assert "Achieve DC(Sat.VX" in text
        assert "{Sat.Mars.Periapsis}" in text


def test_s1l1_mars_nodes_use_consecutive_appc_pairs() -> None:
    """Each Mars node's (v_inf^-, v_inf^+) is the consecutive App-C node pair."""
    row = s1l1_row()
    by_no = {leg[0]: leg for leg in APPC_LEGS}
    mars_nos = [leg[0] for leg in APPC_LEGS if leg[1] == "M"]
    assert len(row.mars_flybys) == len(mars_nos)
    for node, leg_no in zip(row.mars_flybys, mars_nos, strict=True):
        prev_vinf = np.array(by_no[leg_no - 1][3], dtype=np.float64)
        own_vinf = np.array(by_no[leg_no][3], dtype=np.float64)
        np.testing.assert_allclose(node.vinf_minus, prev_vinf, rtol=1e-12)
        np.testing.assert_allclose(node.vinf_plus, own_vinf, rtol=1e-12)


# --- Phase 2.3: generalisation contract --------------------------------------


def test_generator_generalises_from_row_descriptor() -> None:
    """Aldrin and S1L1 flow through the same entry point with different mode."""
    aldrin = aldrin_row()
    s1l1 = s1l1_row()
    assert isinstance(aldrin, CyclerRow)
    assert isinstance(s1l1, CyclerRow)
    assert aldrin.mode == "powered-periodic"
    assert s1l1.mode == "flyby-station-keep"
    # Same templating function renders both.
    assert "BeginMissionSequence" in generate_script(aldrin)
    assert "BeginMissionSequence" in generate_script(s1l1)


def test_write_scripts_to_disk(tmp_path) -> None:  # type: ignore[no-untyped-def]
    a = generate_aldrin_script(tmp_path / "aldrin.script")
    assert a.read_text().count("Target DC") == 1
    s_paths = generate_s1l1_scripts(tmp_path / "s1l1")
    assert len(s_paths) == 7  # one runnable script per Mars flyby
    for p in s_paths:
        assert p.read_text().count("Target DC") == 1
