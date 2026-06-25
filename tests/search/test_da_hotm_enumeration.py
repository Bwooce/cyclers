"""#450 Task 6: DA/HOTM enumeration driver + novelty routing.

The driver loops (C, n), enumerates fixed points, runs the filter cascade
(section-residual -> dedup -> reproduction screen -> corrector), and classifies
each survivor into {reproduction | known-family | novel-PO | novel-cycler-candidate}.
It reuses the existing modules unchanged (corrector, prefilter, lit-check). The
reproduction/known-family screen and the lit-check are INJECTABLE so the test runs
offline.

The base-family triangulation (n=1 DRO routed to reproduction) is the
anti-circularity guard: the SAME untuned machinery that surfaces a hybrid family
must also recover the known base families.
"""

from __future__ import annotations

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_enumerator import DomainBox
from cyclerfinder.search.da_hotm_enumeration import (
    EnumerationResult,
    classify_candidate,
    run_enumeration,
)


def _em() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Earth", "Moon")


def test_n1_dro_routes_to_reproduction() -> None:
    """Base-family triangulation: the n=1 DRO is recovered and routed to repro."""
    system = _em()
    box = DomainBox(x_lo=0.86, x_hi=0.91, xdot_lo=-0.03, xdot_hi=0.03)

    # Injected reproduction oracle: the DRO section point IS a known family.
    def known(x0: float, xdot0: float, n: int) -> bool:
        return abs(x0 - 0.88500968) < 0.01 and abs(xdot0) < 0.01

    results = run_enumeration(
        system,
        c_band=(3.00022,),
        n_range=(1,),
        domain_box=box,
        is_known_family=known,
        residual_tol=3e-2,
        grid=(17, 13),
    )
    dro = [r for r in results if abs(r.x0 - 0.88500968) < 0.02]
    assert dro, [(round(r.x0, 4), r.classification) for r in results]
    assert dro[0].classification == "reproduction", dro[0]


def test_classification_buckets() -> None:
    """classify_candidate routes by the known-family screen + cycler structure."""
    # Known family -> reproduction.
    assert classify_candidate(is_known=True, has_cycler_structure=False) == "reproduction"
    # Unknown PO, no cycler leg -> novel-PO (the Png' case).
    assert classify_candidate(is_known=False, has_cycler_structure=False) == "novel-PO"
    # Unknown with a cycler transfer leg -> novel-cycler-candidate.
    assert classify_candidate(is_known=False, has_cycler_structure=True) == "novel-cycler-candidate"


def test_result_ledger_carries_provenance() -> None:
    """Each result carries IC, residual, period, classification, n, C."""
    system = _em()
    box = DomainBox(x_lo=0.87, x_hi=0.90, xdot_lo=-0.02, xdot_hi=0.02)
    results = run_enumeration(
        system,
        c_band=(3.00022,),
        n_range=(1,),
        domain_box=box,
        is_known_family=lambda x0, xdot0, n: True,
        residual_tol=3e-2,
        grid=(13, 11),
    )
    assert results
    r = results[0]
    assert isinstance(r, EnumerationResult)
    assert r.n == 1
    assert r.c_target == 3.00022
    assert r.section_residual >= 0.0
    assert r.classification in {
        "reproduction",
        "known-family",
        "novel-PO",
        "novel-cycler-candidate",
        "uncertified",
    }


def test_empty_band_yields_no_results() -> None:
    """An off-family band yields an empty ledger (a legitimate negative)."""
    system = _em()
    box = DomainBox(x_lo=0.60, x_hi=0.62, xdot_lo=0.30, xdot_hi=0.32)
    results = run_enumeration(
        system,
        c_band=(3.00022,),
        n_range=(1,),
        domain_box=box,
        is_known_family=lambda x0, xdot0, n: False,
        residual_tol=1e-3,
        grid=(9, 9),
    )
    assert results == []
