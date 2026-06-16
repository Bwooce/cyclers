"""Tests for the #325 topology audit harness.

The audit harness re-verifies discovery-JSONL topology claims through
independent, sourced checkers.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.tulip import TULIP_Z_AMPLITUDE_FLOOR_NONDIM
from cyclerfinder.search.topology_audit import (
    INDEPENDENT_CLOSURE_FLOOR_NONDIM,
    REAL_AXIS_TOL,
    UNIT_CIRCLE_TOL,
    audit_topology,
    check_floquet_neimark_sacker,
    check_periodic_orbit_closure,
    check_tulip_topology,
    write_findings,
)


def test_real_eigenvalue_at_plus_one_is_not_neimark_sacker() -> None:
    """The trivial Liouville pair (+1, +1) MUST NOT classify as Neimark-Sacker."""
    is_ns, extras = check_floquet_neimark_sacker(1.0, 0.0)
    assert not is_ns
    assert extras["imag_part"] == 0.0
    assert extras["unit_circle_distance"] == 0.0


def test_real_eigenvalue_at_minus_one_is_not_neimark_sacker() -> None:
    """A real eigenvalue at -1 is period-doubling, not Neimark-Sacker."""
    is_ns, extras = check_floquet_neimark_sacker(-1.0, 0.0)
    assert not is_ns
    assert extras["imag_part"] == 0.0


def test_complex_pair_on_unit_circle_is_neimark_sacker() -> None:
    """A genuine complex pair on the unit circle IS Neimark-Sacker."""
    angle = np.pi / 4
    is_ns, extras = check_floquet_neimark_sacker(np.cos(angle), np.sin(angle))
    assert is_ns
    assert extras["imag_part"] > REAL_AXIS_TOL
    assert extras["unit_circle_distance"] < UNIT_CIRCLE_TOL


def test_complex_eig_far_off_unit_circle_is_not_neimark_sacker() -> None:
    """A complex pair OFF the unit circle (|lambda| != 1) is not NS."""
    is_ns, extras = check_floquet_neimark_sacker(0.5, 0.5)
    assert not is_ns
    assert extras["unit_circle_distance"] > UNIT_CIRCLE_TOL


def test_borderline_real_within_real_axis_tol_is_not_neimark_sacker() -> None:
    """Numerical-noise real eigenvalues (small imag from numerics) reject NS."""
    is_ns, extras = check_floquet_neimark_sacker(-1.0, 1e-5)
    assert not is_ns
    assert extras["imag_part"] < REAL_AXIS_TOL


def test_planar_orbit_below_floor_is_not_3d() -> None:
    """A planar orbit (z0 = 0) MUST fail the 3D-topology gate."""
    system = cr3bp.cr3bp_system("Earth", "Moon")
    state0 = np.array([0.836, 0.0, 0.0, 0.0, 0.05, 0.0], dtype=np.float64)
    period = 0.5
    is_3d, extras = check_tulip_topology(state0, period, system)
    assert not is_3d
    assert extras["z0"] == 0.0
    assert extras["max_abs_z"] < TULIP_Z_AMPLITUDE_FLOOR_NONDIM


def test_genuine_3d_ic_above_floor_is_3d() -> None:
    """An IC with z0 = 2 * floor MUST pass the 3D-amplitude gate."""
    system = cr3bp.cr3bp_system("Earth", "Moon")
    z0 = 2.0 * TULIP_Z_AMPLITUDE_FLOOR_NONDIM
    state0 = np.array([0.836, 0.0, z0, 0.0, 0.05, 0.0], dtype=np.float64)
    period = 0.5
    is_3d, extras = check_tulip_topology(state0, period, system)
    assert is_3d
    assert extras["max_abs_z"] >= TULIP_Z_AMPLITUDE_FLOOR_NONDIM


def test_non_periodic_ic_does_not_close() -> None:
    """A blatantly non-periodic IC must fail the Radau cross-check."""
    system = cr3bp.cr3bp_system("Earth", "Moon")
    state0 = np.array([0.836, 0.0, 0.0, 0.0, 0.05, 0.0], dtype=np.float64)
    period = 0.123
    closes, extras = check_periodic_orbit_closure(state0, period, system)
    assert not closes
    assert extras["closure_residual"] > INDEPENDENT_CLOSURE_FLOOR_NONDIM


def test_audit_topology_flags_planar_collapse_row(tmp_path: Path) -> None:
    """End-to-end: a synthetic JSONL row with z0 ~ 0 gets flagged as planar."""
    header = {
        "type": "header",
        "issue": 9999,
        "system": {
            "primary": "Earth",
            "secondary": "Moon",
            "mu": 0.012150584270572,
            "l_km": 384400.0,
            "t_s": 375699.8,
        },
    }
    row = {
        "step_index": 0,
        "state_nd": [0.836, 0.0, 1e-12, 0.0, 0.05, 0.0],
        "T_TU": 0.5,
    }
    out = tmp_path / "synthetic.jsonl"
    with out.open("w") as fh:
        fh.write(json.dumps(header) + "\n")
        fh.write(json.dumps(row) + "\n")

    findings = audit_topology(str(out))
    tulip_findings = [f for f in findings if f.genome_gate == "tulip_or_3d_periodic"]
    assert len(tulip_findings) == 1
    f = tulip_findings[0]
    assert f.genome_verdict is True
    assert f.independent_verdict is False
    assert f.discrepancy is True
    assert f.failure_mode == "planar_collapse_under_322_floor"


def test_audit_topology_skips_degenerate_planar_flagged_rows(tmp_path: Path) -> None:
    """Rows already tagged ``degenerate_planar=True`` are skipped."""
    header = {
        "type": "header",
        "issue": 9999,
        "system": {
            "primary": "Earth",
            "secondary": "Moon",
            "mu": 0.012150584270572,
            "l_km": 384400.0,
            "t_s": 375699.8,
        },
    }
    row = {
        "step_index": 138,
        "state_nd": [0.836, 0.0, 1e-12, 0.0, 0.05, 0.0],
        "T_TU": 0.5,
        "degenerate_planar": True,
    }
    out = tmp_path / "honest_planar.jsonl"
    with out.open("w") as fh:
        fh.write(json.dumps(header) + "\n")
        fh.write(json.dumps(row) + "\n")
    findings = audit_topology(str(out))
    tulip_findings = [f for f in findings if f.genome_gate == "tulip_or_3d_periodic"]
    assert tulip_findings == []


def test_audit_topology_floquet_bracket_neimark_sacker_classification(tmp_path: Path) -> None:
    """A bracket_inventory with a real eig-pair near -1 must be flagged as NOT-NS."""
    header = {
        "type": "header",
        "issue": 9999,
        "bracket_inventory": [
            {
                "k": 2,
                "classification": "neimark_sacker",
                "eig_a_re": -1.0,
                "eig_a_im": 0.0,
                "eig_b_re": -1.0,
                "eig_b_im": 0.0,
            }
        ],
    }
    out = tmp_path / "synthetic_bracket.jsonl"
    with out.open("w") as fh:
        fh.write(json.dumps(header) + "\n")

    findings = audit_topology(str(out))
    ns_findings = [f for f in findings if f.genome_gate == "neimark_sacker"]
    assert len(ns_findings) == 2
    for f in ns_findings:
        assert f.discrepancy is True
        assert f.failure_mode == "real_eig_near_-1_period_doubling_k2"


def test_audit_topology_handles_missing_file() -> None:
    """audit_topology raises FileNotFoundError on a missing path."""
    with pytest.raises(FileNotFoundError):
        audit_topology("/nonexistent/path/does_not_exist.jsonl")


def test_write_findings_round_trip(tmp_path: Path) -> None:
    """write_findings produces a JSONL that reloads to the same payload."""
    from cyclerfinder.search.topology_audit import TopologyAuditFinding

    finding = TopologyAuditFinding(
        source_jsonl="/tmp/foo.jsonl",
        row_id="0",
        genome_gate="tulip_or_3d_periodic",
        genome_verdict=True,
        independent_verdict=False,
        discrepancy=True,
        failure_mode="planar_collapse_under_322_floor",
        extras={"max_abs_z": 1e-12, "z_floor": 5e-3},
    )
    out = tmp_path / "findings.jsonl"
    write_findings([finding], str(out))
    with out.open() as fh:
        line = fh.readline()
        payload = json.loads(line)
    assert payload["genome_gate"] == "tulip_or_3d_periodic"
    assert payload["discrepancy"] is True
    assert payload["failure_mode"] == "planar_collapse_under_322_floor"
    assert payload["extras"]["max_abs_z"] == 1e-12
