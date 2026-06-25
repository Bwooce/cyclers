"""Cislunar BCT discovery driver tests (#378 Phase 3).

`search/cislunar_bct_search.py` sweeps the BCT family parameters, constructs +
classifies each as a `transfer` (precursor_mga capability output) or a
`quasi_cycler_candidate` (the return leg re-acquires W -- the only catalogue-row
candidate, design §6 / Table 1.1), runs the BCT novelty self-test, and emits
records. Per the design's critical honesty, the expected default outcome is
capability-only (transfers), with the repeating-chain quasi-cycler unproven.

Tasks:
  * 3.1 -- classify transfer vs return-leg-re-acquires-W chain.
  * 3.2 -- novel candidate runs lit-check; non-novel suppressed.
"""

from __future__ import annotations

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.genome.bct_transfer as bct
import cyclerfinder.search.cislunar_bct_search as search
from cyclerfinder.search.literature_check import SearchResult


def test_classifies_transfer_vs_chain() -> None:
    """A single BCT classifies as `transfer`; a return-leg-re-acquires-W as chain.

    A single backward-constructed BCT (no W re-acquisition on the return leg) is
    a `transfer` (precursor_mga). A synthetic case where the return leg
    re-acquires W within the budget classifies as `quasi_cycler_candidate`.
    """
    system = bcr4bp.andreu_default()
    target = bct.BCTTarget(theta2=0.70, branch="retrograde")
    arc = bct.construct_bct_backward(target, system, back_days=70.0)

    # A plain BCT (no W re-acquisition) -> transfer.
    label_transfer = search.classify_bct(arc, system, reacquires_w=False)
    assert label_transfer == "transfer"

    # A synthetic return-leg-re-acquires-W -> quasi_cycler_candidate.
    label_chain = search.classify_bct(arc, system, reacquires_w=True)
    assert label_chain == "quasi_cycler_candidate"


def test_search_emits_transfer_capability_records() -> None:
    """run_cislunar_bct_search over a tiny grid emits classified BCT records."""
    system = bcr4bp.andreu_default()
    grid = search.BCTSearchGrid(
        theta2_values=(0.70, 0.75),
        e2_values=(0.95,),
        back_days_values=(70.0,),
        branches=("retrograde",),
    )
    result = search.run_cislunar_bct_search(grid, system)
    assert result.n_evaluated == 2
    # Every record is classified; in the incoherent model the expected default
    # is `transfer` (capability), not a quasi-cycler.
    assert all(r.classification in ("transfer", "quasi_cycler_candidate") for r in result.records)
    # At least one reached the Hiten apoapsis band (the theta2=0.70 family).
    assert any(2.7 <= r.apoapsis_ld <= 5.1 for r in result.records)


def test_novel_candidate_runs_litcheck_nonnovel_suppressed() -> None:
    """A non-novel (Hiten-matching) candidate is suppressed; a novel one emitted.

    The Belbruno/Hiten corpus stub flags the BCT signature `published` ->
    suppressed from the novel-emission set. A stub returning no corpus hit leaves
    the candidate un-suppressed (eligible for human/gauntlet routing).
    """
    sig = bct.bct_candidate_signature(sequence=("E", "M"))

    def hiten_search(_q: str) -> list[SearchResult]:
        return [
            SearchResult(
                title="Belbruno Hiten ballistic capture weak stability boundary transfer",
                url="https://doi.org/10.2514/3.21079",
                snippet="Hiten MUSES-A low energy lunar transfer ballistic capture Moon.",
            )
        ]

    def empty_search(_q: str) -> list[SearchResult]:
        return [SearchResult(title="unrelated orbital mechanics", url="x", snippet="no match")]

    assert search.is_novel_emittable(sig, search_fn=hiten_search) is False
    assert search.is_novel_emittable(sig, search_fn=empty_search) is True
