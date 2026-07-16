"""#607 -- small-body multi-moon symmetric-closure sweep tests.

Two kinds of check, matching this project's own precedent
(``test_enumerate_563_symmetric_closures.py``):

1. Registry sanity: the new `core/satellites.py` PRIMARIES/SATELLITES
   entries (Sylvia/Elektra/Eugenia/Kleopatra + their moons) resolve, give a
   positive ``n_max`` (a non-degenerate search space) for every swept pair,
   and Lempo-Paha-Hiisi remains deliberately ABSENT (structural exclusion
   guard -- see `enumerate_607_smallbody_multimoon_symmetric_closures.py`'s
   own docstring for why).
2. Regression pin for the fast (2-moon) systems: reproduces the actual
   discovery-run counts (`n_evaluated`, `n_all_gates_passed`) committed in
   `data/enumerate_607_smallbody_multimoon_symmetric_closures.jsonl`. The
   Elektra 3-moon sweep (96,768 evaluations, ~570s) is intentionally NOT
   re-run here -- only smoke-tested on one sequence with a small
   ``tof_scale_max`` override to keep this file fast; its full-scale result
   is committed data, not re-derived on every CI run.
"""

from __future__ import annotations

import json
from pathlib import Path

import scripts.enumerate_563_symmetric_closures as enum563
import scripts.enumerate_600_3moon_symmetric_closures as enum600
import scripts.enumerate_607_smallbody_multimoon_symmetric_closures as enum607
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

_ROOT = Path(__file__).resolve().parent.parent.parent
_COMMITTED_PATH = _ROOT / "data" / "enumerate_607_smallbody_multimoon_symmetric_closures.jsonl"


def test_new_primaries_registered() -> None:
    for primary in ("Sylvia", "Elektra", "Eugenia", "Kleopatra"):
        assert primary in PRIMARIES
        assert PRIMARIES[primary] > 0.0


def test_new_satellites_registered_with_correct_primary() -> None:
    expected = {
        "Romulus": "Sylvia",
        "Remus": "Sylvia",
        "ElektraBeta": "Elektra",
        "ElektraGamma": "Elektra",
        "ElektraDelta": "Elektra",
        "PetitPrince": "Eugenia",
        "EugeniaS2": "Eugenia",
        "AlexHelios": "Kleopatra",
        "CleoSelene": "Kleopatra",
    }
    for moon, primary in expected.items():
        assert moon in SATELLITES
        assert SATELLITES[moon].primary == primary
        assert SATELLITES[moon].mu_km3_s2 > 0.0
        assert SATELLITES[moon].sma_km > 0.0


def test_lempo_paha_hiisi_deliberately_absent() -> None:
    """Structural exclusion guard (Lempo:Hiisi mass ratio ~1.27:1, a genuine
    near-equal-mass binary -- see the #607 script's module docstring). If
    this ever needs to change, it must be a deliberate, reviewed decision
    (e.g. a genuine CR3BP-based treatment), not an accidental re-add via the
    same point-mass-primary construction."""
    for name in ("Lempo", "Paha", "Hiisi"):
        assert name not in PRIMARIES
        assert name not in SATELLITES


def test_pair_n_max_positive_for_every_swept_pair() -> None:
    pairs = [
        ("Romulus", "Remus", "Sylvia"),
        ("PetitPrince", "EugeniaS2", "Eugenia"),
        ("AlexHelios", "CleoSelene", "Kleopatra"),
        ("ElektraBeta", "ElektraGamma", "Elektra"),
        ("ElektraBeta", "ElektraDelta", "Elektra"),
        ("ElektraGamma", "ElektraDelta", "Elektra"),
    ]
    for a, b, primary in pairs:
        _t_syn, _p_a, _p_b, n_max = enum563.pair_n_max(a, b, primary=primary)
        assert n_max >= 1, f"{a}-{b} ({primary}) has a degenerate (empty) search space"


def test_fast_2moon_systems_reproduce_committed_clean_negative() -> None:
    """Re-runs the 3 fast (2-moon) systems fully and diffs against the
    committed discovery-run counts -- a regression pin, not a fresh claim."""
    committed: dict[str, dict[str, int]] = {}
    with _COMMITTED_PATH.open() as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("_meta"):
                for sysinfo in d["systems"]:
                    committed[sysinfo["primary"]] = sysinfo

    for primary, moons in (
        ("Sylvia", ("Romulus", "Remus")),
        ("Eugenia", ("PetitPrince", "EugeniaS2")),
        ("Kleopatra", ("AlexHelios", "CleoSelene")),
    ):
        results = enum607.run_2moon_system(primary, moons)
        n_evaluated = sum(r["n_evaluated"] for r in results)
        n_passed = sum(r["n_all_gates_passed"] for r in results)
        assert n_evaluated == committed[primary]["n_evaluated"], primary
        assert n_passed == committed[primary]["n_all_gates_passed"] == 0, primary


def test_elektra_3moon_construction_smoke() -> None:
    """Cheap smoke test only (small tof_scale_max, ONE sequence) -- the full
    6-sequence/96768-evaluation/~570s Elektra sweep is committed data, not
    re-derived here. Just confirms the #600 3-moon-chain construction runs
    end-to-end against the new Elektra registry entries without error and
    returns the expected structure."""
    res = enum600.enumerate_sequence(
        "ElektraBeta", "ElektraGamma", "ElektraDelta", primary="Elektra", tof_scale_max=2.0
    )
    assert res["anchor"] == "ElektraBeta"
    assert res["sequence"] == ["ElektraBeta", "ElektraGamma", "ElektraDelta", "ElektraBeta"]
    assert res["n_evaluated"] > 0
    assert isinstance(res["passes"], list)


def test_committed_output_shows_clean_negative_overall() -> None:
    with _COMMITTED_PATH.open() as fh:
        meta = json.loads(fh.readline())
    assert meta["total_all_gates_passed"] == 0
    assert meta["total_evaluated"] > 0
    excluded_systems = {e["system"] for e in meta["excluded"]}
    assert "Lempo-Paha-Hiisi" in excluded_systems
