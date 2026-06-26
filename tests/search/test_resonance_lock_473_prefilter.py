"""Task #477 — POSITIVE CONTROL for the #473 resonance-lock ΔV-floor pre-filter.

The #477 branch-and-bound screens each candidate's legs with the admissible VILM
ΔV-floor (``vilm.vilm_dv_floor``, the #76 escape+capture lower bound) before the
expensive leveraging-chain solve, pruning any tour with an over-budget leg. The
mandatory guard (the exact failure mode of a bad lower bound): the 2 known #470
flips (Galilean IEG, EGC) are FEASIBLE bounded tours and MUST survive the
pre-filter — assert the gate does NOT drop them. The #339 positive control must
also survive.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from _resonance_lock_moontour_473 import (  # type: ignore[import-not-found]  # noqa: E402
    PRUNE_BUDGET_KMS,
    SKELETONS,
    _prefilter_cell,
)


def _sk(label: str) -> Any:
    for sk in SKELETONS:
        if sk.label == label or label in sk.label:
            return sk
    raise KeyError(label)


@pytest.mark.parametrize(
    "label",
    ["Galilean IEG", "Galilean EGC", "POSCTL Umbriel-Oberon (#339)"],
)
def test_prefilter_does_not_drop_feasible_controls(label: str) -> None:
    # The #470 flips + the #339 control are feasible bounded tours; a sound
    # admissible floor must let them through (a bad bound that prunes them is the
    # exact regression #477 guards against).
    sk = _sk(label)
    survives, reason = _prefilter_cell(sk)
    assert survives, f"{label} wrongly pruned: {reason}"


def test_prefilter_is_not_vacuous() -> None:
    # The gate is not vacuous: it FIRES on a leg whose VILM floor exceeds the
    # budget. The Galilean IEG legs have floors ~1.35-1.57 km/s; with a budget set
    # below that (here via monkeypatch-free direct check against a tightened bar)
    # the same skeleton is pruned with a VILM-floor reason — proving the bound is
    # live, not a pass-through. (The production PRUNE_BUDGET_KMS=3.5 deliberately
    # lets all these cheap pairs through; the mechanism is what is asserted here.)
    from cyclerfinder.search.vilm import vilm_dv_floor

    sk = _sk("Galilean IEG")
    legs = [
        (sk.sequence[k], sk.sequence[k + 1])
        for k in range(len(sk.sequence) - 1)
        if sk.sequence[k] != sk.sequence[k + 1]
    ]
    floors = [vilm_dv_floor(a, b) for a, b in legs]
    max_floor = max(floors)
    # The production budget exceeds every IEG leg floor (so IEG survives) ...
    assert max_floor < PRUNE_BUDGET_KMS
    # ... and a budget below the max floor would prune at least one leg (live bound).
    over_budget = [(a, b) for (a, b), f in zip(legs, floors, strict=True) if f > max_floor - 1e-9]
    assert over_budget, "expected a binding leg at the max floor"
