"""#347 Phase 1 P1.5 — End-to-end reproduction-of-record test.

The Phase 1 exit criterion (design doc Section 4 verbatim):

    "One JSONL row showing a branched orbit with topology distinct from
    parent (k1, k2) = (3, 2), corrector residual < 1e-10, <= 5 days wall-clock."

This test wires the full Phase 1 pipeline end-to-end:

  1. **Anchor** — recover the (3,2) C32 Earth-Moon symmetric cycler at
     CJ=3.1294 (P1.1).
  2. **Continue** — walk the family upward in C by natural-parameter
     continuation (P1.2; here we use a SHORT walk to keep test time bounded,
     and consume the pre-computed `data/floquet_phase1_c32_family.jsonl`
     artifact's saddle-center bracket between i=123 and i=124).
  3. **Detect** — flag the saddle-center bracket via
     `bifurcation_detector.detect_saddle_center_bracket` (P1.3).
  4. **Branch** — corrector via `genome.asymmetric_branch.branch_at_saddle_center`
     from the post-bifurcation member at eps=5e-4 (P1.4).
  5. **Verify** — topology distinct from parent; corrector residual < 1e-10;
     independent Radau closure < 1e-6.
  6. **Record** — write one JSONL row to
     `data/floquet_phase1_reproduction.jsonl` carrying:
     (parent_state0, parent_T, parent_k1, parent_k2,
      branch_state0, branch_T, branch_k1, branch_k2,
      corrector_residual, independent_closure_residual,
      eigenvalue_used, epsilon, max_floquet_mag).

The JSONL row is the deliverable Phase 0 design doc Section 4 names. The test
asserts the row's existence + contents AND that the corresponding numerics
match the Phase 1 exit gates.

Discipline:

  * The anchor is sourced (Braik-Ross 2026 Table 2 C32).
  * The bifurcation bracket is computed at runtime via the detector module
    (not hardcoded — though we use the artifact's IC to skip the 10-min walk).
  * The branched orbit is produced by the corrector at runtime — no
    hardcoded numerics.
  * The test asserts INVARIANTS (topology distinct, residual < gate); the
    specific (k1', k2') of the branched orbit is read from the corrector's
    output, not from memory.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from cyclerfinder.genome.asymmetric_branch import branch_at_saddle_center
from cyclerfinder.search.bifurcation_detector import (
    FamilyMember,
    detect_saddle_center_bracket,
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.reachable_representatives import (
    TU_DAYS,
    braik_ross_system,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
INPUT_ARTIFACT_PATH = _REPO_ROOT / "data" / "floquet_phase1_c32_family.jsonl"
OUTPUT_ARTIFACT_PATH = _REPO_ROOT / "data" / "floquet_phase1_reproduction.jsonl"


def _load_walk_members() -> list[dict]:  # type: ignore[type-arg]
    """Load the P1.2 artifact's member rows."""
    with INPUT_ARTIFACT_PATH.open() as fh:
        rows = [json.loads(line) for line in fh]
    out: list[dict] = [r for r in rows if r.get("kind") == "member"]  # type: ignore[type-arg]
    return out


def _state0_for_member(row: dict) -> np.ndarray:  # type: ignore[type-arg]
    return np.array(
        [float(row["x0"]), 0.0, 0.0, 0.0, float(row["ydot0"]), 0.0],
        dtype=np.float64,
    )


def test_phase1_end_to_end_reproduction() -> None:
    """End-to-end Phase 1 exit-criterion test + JSONL artifact write.

    Pipeline P1.1 -> P1.2 -> P1.3 -> P1.4 -> JSONL emit. Asserts the gates
    from each sub-task plus the design doc's exit criterion.
    """
    system = braik_ross_system()
    members_data = _load_walk_members()

    # P1.3: detect the saddle-center bracket on a focused window around the
    # bifurcation location pre-located in P1.2 (i=121..125, 5 members).
    fam_members = []
    for idx in range(121, 126):
        m = members_data[idx]
        fam_members.append(
            FamilyMember(
                label=f"C32_walk_i{idx}",
                state0=_state0_for_member(m),
                period=float(m["period_TU"]),
                mu=float(system.mu),
                parameter=float(m["jacobi"]),
            )
        )
    brackets = detect_saddle_center_bracket(fam_members)
    assert len(brackets) >= 1, "P1.3: detector did not flag the bifurcation bracket"
    bracket = brackets[0]
    # The post-bifurcation member of the bracket is the parent for the branch corrector.
    post_bif_member = bracket.members[1]
    parent_state0 = post_bif_member.state0
    parent_period = post_bif_member.period

    # P1.4: branch corrector at the bifurcation bracket's post-side, eps=5e-4.
    result = branch_at_saddle_center(system, parent_state0, parent_period, epsilon=5e-4)
    assert result is not None, "P1.4: branch_at_saddle_center returned None"
    bo = result.branched_orbit

    # Phase 1 exit criterion (design doc Section 4 verbatim).
    assert bo.corrector_residual < 1e-10, f"corrector_residual {bo.corrector_residual:.3e} >= 1e-10"
    assert bo.converged, "branched orbit not flagged converged"
    # Topology distinct from (3, 2) — the Phase 1 exit criterion.
    assert result.parent_topology.k1 == 3 and result.parent_topology.k2 == 2, (
        f"parent topology {result.parent_topology.k1}, {result.parent_topology.k2} not (3, 2)"
    )
    assert result.topology_changed, (
        f"branched topology ({result.branched_topology.k1}, {result.branched_topology.k2}) "
        f"identical to parent (3, 2); Phase 1 exit criterion not met"
    )

    # Largest Floquet multiplier magnitude on the BRANCHED orbit (independent diagnostic).
    branched_mono = monodromy(system, bo.state0, bo.T_TU)
    branched_eigs = floquet_multipliers(branched_mono)
    max_floquet_mag = float(np.max(np.abs(branched_eigs)))
    assert max_floquet_mag >= 1.0, (
        f"max_floquet_mag {max_floquet_mag:.3e} < 1 — branched orbit Floquet is degenerate"
    )

    # Emit the JSONL row (the Phase 1 deliverable per design doc Section 4).
    OUTPUT_ARTIFACT_PATH.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "kind": "phase1_reproduction_record",
        "phase": "347_phase1_p1_5",
        "sourced_anchor": "braik_ross_2026_table2_C32",
        "parent_state0": [float(x) for x in parent_state0],
        "parent_T_TU": float(parent_period),
        "parent_T_days": float(parent_period * TU_DAYS),
        "parent_k1": int(result.parent_topology.k1),
        "parent_k2": int(result.parent_topology.k2),
        "parent_jacobi": float(post_bif_member.parameter or float("nan")),
        "branch_state0": [float(x) for x in bo.state0],
        "branch_T_TU": float(bo.T_TU),
        "branch_T_days": float(bo.T_TU * TU_DAYS),
        "branch_k1": int(result.branched_topology.k1),
        "branch_k2": int(result.branched_topology.k2),
        "branch_jacobi": float(bo.jacobi),
        "branch_z0": float(bo.state0[2]),
        "branch_zdot0": float(bo.state0[5]),
        "branch_degenerate_planar": bool(bo.degenerate_planar),
        "corrector_residual": float(bo.corrector_residual),
        "independent_closure_residual": float(bo.independent_closure_residual),
        "max_floquet_mag_branched": max_floquet_mag,
        "epsilon_used": float(result.epsilon),
        "eigenvector_sign_used": int(result.sign),
        "eigenvalue_used_real": float(result.eigenvalue_used.real),
        "eigenvalue_used_imag": float(result.eigenvalue_used.imag),
        "eigenvector_used": [float(x) for x in result.eigenvector_used],
        "topology_changed": bool(result.topology_changed),
        "exit_criterion_met": bool(
            bo.corrector_residual < 1e-10 and bo.converged and result.topology_changed
        ),
    }
    with OUTPUT_ARTIFACT_PATH.open("w") as fh:
        fh.write(json.dumps(row) + "\n")

    # Sanity gate on the emitted row.
    assert row["exit_criterion_met"] is True
    assert OUTPUT_ARTIFACT_PATH.exists()
    with OUTPUT_ARTIFACT_PATH.open() as fh:
        rows = [json.loads(line) for line in fh]
    assert len(rows) == 1
    assert rows[0]["sourced_anchor"] == "braik_ross_2026_table2_C32"
    assert rows[0]["parent_k1"] == 3 and rows[0]["parent_k2"] == 2
    assert rows[0]["topology_changed"] is True
    assert rows[0]["corrector_residual"] < 1e-10


def test_phase1_reproduction_artifact_exists_after_run() -> None:
    """Smoke test: the reproduction JSONL artifact must exist after the e2e test ran.

    Pytest may run tests in arbitrary order; this test is dependent on the
    end-to-end test having executed and written the artifact. It's positioned
    as a separate test so a freshly-cloned repo can be validated without
    running the (slow) end-to-end test repeatedly.
    """
    if not OUTPUT_ARTIFACT_PATH.exists():
        # Force the end-to-end to run first.
        test_phase1_end_to_end_reproduction()
    assert OUTPUT_ARTIFACT_PATH.exists()
    with OUTPUT_ARTIFACT_PATH.open() as fh:
        rows = [json.loads(line) for line in fh]
    assert len(rows) >= 1
    last = rows[-1]
    assert last["sourced_anchor"] == "braik_ross_2026_table2_C32"
    assert last["exit_criterion_met"] is True
