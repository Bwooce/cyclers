"""Shared fixtures for the n-body SILVER rung tests (plan Phase B).

The two USER-HELD SILVER candidates' runtime ``data/review_queue.jsonl`` artefact
was NOT present on disk at execution time (concurrent agents; same situation the
diligence note records, ``docs/notes/2026-06-06-silver-candidates-russell-
diligence.md`` lines 22-25, 239-241). The candidates are therefore RECONSTRUCTED
from their recorded provenance:

* candidate 1 — ``OUTSTANDING.md`` + the diligence note: sequence E-M-E-E,
  2-synodic (period_k=2), per-encounter V∞ **[E 9.75, M 13.01, E 9.76, E 9.75]**
  km/s, ``match=novel``, bend-feasible, panel-survived.
* candidate 2 — same sources: sequence E-M-E-E, 2-synodic, V∞ **[E 9.62, M 12.06,
  …]** km/s (the trailing two E encounters were "not reported in the section"; we
  carry the two reported magnitudes and mirror the home-Earth value for the
  unreported return encounters, FLAGGED).

The ``tof_days`` per leg are NOT recorded anywhere on disk (neither OUTSTANDING,
the diligence note, nor the gauntlet ledger carry them). We use a representative
2-synodic E-M-E-E leg structure derived from the period: a ~150 d outbound E→M
transit, a ~330 d M→E return, and the balance of the 2-synodic period (~4.5 yr)
on the E→E loop leg. This is a REPRESENTATIVE reconstruction for the rung; the
rung's verdict is regime-level (correction-ΔV band), so it is robust to the exact
leg split — but the absence of the on-disk seed is a recorded limitation.
"""

from __future__ import annotations

import pytest

from cyclerfinder.data.review_queue import ReviewQueueEntry

# 2-synodic Earth-Mars repeat period ~ 2 * 2.135 yr synodic ~ 4.27 yr; we use a
# representative leg split summing to ~4.5 yr (the Russell 2-synodic class span).
_OUTBOUND_EM_DAYS = 150.0
_RETURN_ME_DAYS = 330.0
_LOOP_EE_DAYS = 4.27 * 365.25 - _OUTBOUND_EM_DAYS - _RETURN_ME_DAYS


def _silver_entry(candidate_id: str, vinf: tuple[float, ...], max_vinf: float) -> ReviewQueueEntry:
    return ReviewQueueEntry(
        candidate_id=candidate_id,
        signature_hash=f"reconstructed-{candidate_id}",
        verdict_tier="silver",
        match_outcome="novel",
        known_id=None,
        superseded_by=(),
        vinf_per_encounter_kms=vinf,
        tof_days=(_OUTBOUND_EM_DAYS, _RETURN_ME_DAYS, _LOOP_EE_DAYS),
        bend_feasible=True,
        max_vinf_kms=max_vinf,
        sequence=("E", "M", "E", "E"),
        period_k=2,
        model_assumption="real-de440",
        verdict_audit={"reconstructed_from": "OUTSTANDING.md + diligence note #126"},
        panel={},
        t_added="2026-06-06T00:00:00Z",
    )


@pytest.fixture
def silver_fixture() -> ReviewQueueEntry:
    """Held SILVER candidate 1 (E-M-E-E, V∞ [9.75, 13.01, 9.76, 9.75] km/s)."""
    return _silver_entry("forge-silver-1", (9.75, 13.01, 9.76, 9.75), 13.01)


@pytest.fixture
def silver_candidate_2() -> ReviewQueueEntry:
    """Held SILVER candidate 2 (E-M-E-E, V∞ [9.62, 12.06, ...] km/s).

    The two return-E encounter magnitudes were not reported; mirror the home-Earth
    9.62 value for them (FLAGGED — see module docstring).
    """
    return _silver_entry("forge-silver-2", (9.62, 12.06, 9.62, 9.62), 12.06)
