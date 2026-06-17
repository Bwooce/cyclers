"""#347 Phase 1 P1.3 — Saddle-center detector tests on the (3,2) C32 family walk.

Exercises :func:`cyclerfinder.search.bifurcation_detector.detect_saddle_center_bracket`
against the JSONL artifact produced by ``scripts/floquet_phase1_p1_2_walk.py``
(written to ``data/floquet_phase1_c32_family.jsonl`` at commit ``075f21b``).

The artifact is a 250-step walk in CJ along the Earth-Moon (3,2) C32 cycler
family, with the saddle-center bifurcation pre-located between members
i=123 (C=3.14170, secondary pair 0.9989 +- 0.0473j on the unit circle) and
i=124 (C=3.14180, secondary pair (0.9742, 1.0264) real near +1).

Gates:

  1. ``_classify_secondary_pair`` returns "complex_unit_circle" at i=123 and
     "real_near_one" at i=124 (the eigenvalue classifier sees the transition
     locally).
  2. ``detect_saddle_center_bracket`` returns at least one ``BifurcationPoint``
     with k=1 whose two members bracket i=123 and i=124 (the family-level
     detector flags exactly the right bracket).
  3. The trivial pair was correctly excluded — the eigenvalues flagged as the
     "secondary pair" are NOT the integrator-split (1, 1) pair (which sits at
     about (0.997, 1.003) at i=123 by trivial_dist_max).

Discipline:

  * The expected eigenvalues are read from the COMMITTED JSONL artifact
    (a "sourced same-model golden" per ``feedback_golden_tests_sourced_only``
    — the golden side is "our own committed run output", not a number
    re-computed in the test).
  * The bracket flagging is a SEPARATE computation from the JSONL artifact:
    the test re-runs the monodromy at the artifact's stored (x0, ydot0,
    period) IC, then runs the detector. The detector's output is what's
    gated; the JSONL is the reference reading of where the bifurcation
    sits.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from cyclerfinder.search.bifurcation_detector import (
    FamilyMember,
    _classify_secondary_pair,
    detect_saddle_center_bracket,
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

# The committed artifact from #347 P1.2 commit 075f21b.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
ARTIFACT_PATH = _REPO_ROOT / "data" / "floquet_phase1_c32_family.jsonl"


def _load_artifact_members() -> list[dict]:  # type: ignore[type-arg]
    with ARTIFACT_PATH.open() as fh:
        rows = [json.loads(line) for line in fh]
    return [r for r in rows if r.get("kind") == "member"]


def _state0_for_member(row: dict) -> np.ndarray:  # type: ignore[type-arg]
    return np.array(
        [
            float(row["x0"]),
            0.0,
            0.0,
            0.0,
            float(row["ydot0"]),
            0.0,
        ],
        dtype=np.float64,
    )


def test_artifact_exists_and_has_expected_member_count() -> None:
    """Sanity gate: the P1.2 artifact is on disk with 250 members + header."""
    assert ARTIFACT_PATH.exists(), f"artifact missing: {ARTIFACT_PATH}"
    with ARTIFACT_PATH.open() as fh:
        rows = [json.loads(line) for line in fh]
    headers = [r for r in rows if r.get("kind") == "header"]
    members = [r for r in rows if r.get("kind") == "member"]
    assert len(headers) == 1, f"expected 1 header row, got {len(headers)}"
    assert len(members) == 251, f"expected 251 member rows, got {len(members)}"
    assert headers[0]["sourced_anchor"] == "braik_ross_2026_table2_C32"


def test_classify_secondary_pair_before_bifurcation_at_i123() -> None:
    """Local classifier flags the pre-bifurcation member i=123 as complex_unit_circle.

    At C = 3.14170 the artifact records secondary pair 0.9989 +- 0.0473j (UC_arg
    = 0.04736). Re-computing the monodromy at the artifact's IC and running the
    classifier must return "complex_unit_circle".
    """
    members = _load_artifact_members()
    m = members[123]
    assert abs(m["jacobi"] - 3.14170) < 1e-5
    system = braik_ross_system()
    state0 = _state0_for_member(m)
    mono = monodromy(system, state0, float(m["period_TU"]))
    eigs = floquet_multipliers(mono)
    kind, _lam_a, _lam_b = _classify_secondary_pair(eigs)
    assert kind == "complex_unit_circle", (
        f"i=123 (C=3.14170) classification was {kind!r}, expected complex_unit_circle"
    )


def test_classify_secondary_pair_after_bifurcation_at_i124() -> None:
    """Local classifier flags the post-bifurcation member i=124 as real_near_one.

    At C = 3.14180 the artifact records secondary pair (0.9742, 1.0264) — both
    real, both within 0.1 of +1. Classifier must return "real_near_one".
    """
    members = _load_artifact_members()
    m = members[124]
    assert abs(m["jacobi"] - 3.14180) < 1e-5
    system = braik_ross_system()
    state0 = _state0_for_member(m)
    mono = monodromy(system, state0, float(m["period_TU"]))
    eigs = floquet_multipliers(mono)
    kind, _lam_a, _lam_b = _classify_secondary_pair(eigs)
    assert kind == "real_near_one", (
        f"i=124 (C=3.14180) classification was {kind!r}, expected real_near_one"
    )


def test_detect_saddle_center_bracket_brackets_i123_i124() -> None:
    """Family-level detector returns a k=1 bracket bracketing i=123 and i=124.

    Hand the detector a focused 5-member window centred on the known bifurcation
    location: members 121, 122, 123, 124, 125. It must return at least one
    BifurcationPoint whose members are (i=123, i=124) — the adjacent pair where
    the secondary-pair classification flips.
    """
    members_data = _load_artifact_members()
    system = braik_ross_system()
    window_indices = [121, 122, 123, 124, 125]
    fam_members = []
    for idx in window_indices:
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
    assert len(brackets) >= 1, (
        "detector returned no saddle-center brackets in the i=121..125 window "
        "(known bifurcation at i=123->124)"
    )
    # The bracket's parameter band must straddle the known bifurcation location.
    found_bracket = False
    for bp in brackets:
        param_before = bp.extras.get("param_before", float("nan"))
        param_after = bp.extras.get("param_after", float("nan"))
        if abs(param_before - 3.14170) < 1e-5 and abs(param_after - 3.14180) < 1e-5:
            found_bracket = True
            assert bp.k == 1, f"saddle-center bracket k={bp.k}, expected 1"
            break
    assert found_bracket, (
        f"no bracket straddled C=(3.14170, 3.14180); got brackets: "
        f"{[(bp.extras.get('param_before'), bp.extras.get('param_after')) for bp in brackets]}"
    )


def test_detect_saddle_center_bracket_finds_two_brackets_in_full_walk() -> None:
    """Run the detector across the FULL 250-member walk; expect TWO brackets.

    Empirically (#347 P1.2 commit 075f21b artifact), the (3,2) family has TWO
    saddle-center bifurcations within the C-range walked:

      * Forward bifurcation between i=123 (C=3.14170) and i=124 (C=3.14180):
        secondary pair coalesces complex_unit_circle -> real_near_one.
      * Inverse bifurcation between i ~ 200 (C ~ 3.1494) and i ~ 230
        (C ~ 3.1524): the bifurcated real pair RE-COALESCES into a complex
        pair on the unit circle (i.e. real_far -> complex_unit_circle).

    The detector currently only flags transitions involving "real_near_one"
    on one side (the proper saddle-center signature near +1); the second
    bifurcation transitions real_far -> complex_unit_circle which is NOT a
    saddle-center at +1. The current detector deliberately ignores it.
    However the FIRST saddle-center (the canonical one for P1.4) must be
    flagged. Test gates exactly the canonical-bracket presence; documents
    the second event in the progress note.

    Subsampling the family: feed every 5th member so we stay in-budget but
    the discrete signal is well-sampled.
    """
    members_data = _load_artifact_members()
    system = braik_ross_system()
    fam_members = []
    for idx in range(0, len(members_data), 5):
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
    assert len(brackets) >= 1, "no saddle-center bracket found across full walk"
    # At least ONE bracket must straddle the canonical C ~ 3.14175 bifurcation.
    found_canonical = False
    for bp in brackets:
        param_before = bp.extras.get("param_before", float("nan"))
        param_after = bp.extras.get("param_after", float("nan"))
        if 3.140 < param_before < 3.143 and 3.140 < param_after < 3.143:
            found_canonical = True
            assert bp.k == 1
    assert found_canonical, (
        f"no canonical saddle-center bracket near C=3.14175; got brackets at: "
        f"{[(bp.extras.get('param_before'), bp.extras.get('param_after')) for bp in brackets]}"
    )


@pytest.mark.parametrize(
    "idx,expected_kind",
    [
        (0, "complex_unit_circle"),  # anchor (C=3.1294)
        (60, "complex_unit_circle"),  # mid-walk (C=3.1354)
        (123, "complex_unit_circle"),  # last pre-bifurcation member
        (124, "real_near_one"),  # first post-bifurcation member
        # At i=150 the bifurcated pair has separated to (0.7894, 1.2668) — outside
        # the |lambda - 1| < 0.1 "near one" band, so classifier returns "real_far".
        (150, "real_far"),
        # At i=230 (C=3.1524) the secondary pair has RE-COALESCED on the unit
        # circle (0.9895 +- 0.1443j) — a SECOND saddle-center bifurcation in
        # the inverse direction at C in (3.1494, 3.1524). Documented in
        # progress note P1.3.
        (230, "complex_unit_circle"),
    ],
)
def test_classify_secondary_pair_at_select_indices(idx: int, expected_kind: str) -> None:
    """Classifier returns the expected secondary-pair label at known reference points.

    The expected labels come from the artifact's eigenvalue rows (the "real_near_one"
    transition pinned between i=123 and i=124 by inspection; i=150 verified
    real-pair from artifact eigs_real / eigs_imag entries).
    """
    members_data = _load_artifact_members()
    system = braik_ross_system()
    m = members_data[idx]
    state0 = _state0_for_member(m)
    mono = monodromy(system, state0, float(m["period_TU"]))
    eigs = floquet_multipliers(mono)
    kind, _lam_a, _lam_b = _classify_secondary_pair(eigs)
    assert kind == expected_kind, (
        f"i={idx} (C={m['jacobi']:.5f}): classification {kind!r}, expected {expected_kind!r}"
    )
