"""Phase 5 — the adversarial panel (Forge plan §34).

N independent re-verifications per candidate: re-solve from perturbed seeds,
independent code-path agreement, falsification probes. Majority-refute kills.

The falsification e2e (a fabricated candidate -> REJECTED) is the gate's teeth.
GOLDEN DISCIPLINE: the probes assert physical self-consistency (bend-feasibility,
V_inf cap, closure residual, seed-perturbation robustness), never a sourced value.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import BallisticClosureResult, ballistic_correct
from cyclerfinder.verify.adversarial import (
    PanelResult,
    adversarial_panel,
)

_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
_PERIOD_DAYS = (1.4612 + 2.8096) * 365.25
_PERIOD_SEC = _PERIOD_DAYS * 86400.0
_SEQ = ("E", "M", "E", "E")
_REVS = (0, 0, 1)
_BRANCH = ("single", "single", "low")


def _real_closure(ephem: Ephemeris) -> BallisticClosureResult:
    t0 = ((datetime(2030, 3, 22, tzinfo=UTC) + timedelta(days=-20)) - _J2000).total_seconds()
    return ballistic_correct(
        sequence=_SEQ,
        per_leg_revs=_REVS,
        per_leg_branch=_BRANCH,
        t0_seed_sec=t0,
        tof_seed_days=(154.0, 379.0),
        period_sec=_PERIOD_SEC,
        ephem=ephem,
        vinf_cap=14.0,
        slack_leg=2,
    )


@pytest.mark.slow
def test_panel_confirms_real_closure() -> None:
    """A genuine, demonstrated DE440 closure survives the adversarial panel:
    not majority-refuted, and seed-perturbation finds the family again."""
    ephem = Ephemeris("astropy")
    closure = _real_closure(ephem)
    assert closure.converged

    panel = adversarial_panel(
        closure,
        sequence=_SEQ,
        per_leg_revs=_REVS,
        per_leg_branch=_BRANCH,
        period_sec=_PERIOD_SEC,
        slack_leg=2,
        vinf_cap=14.0,
        ephem=ephem,
        n_verifiers=3,
        seed=0,
    )
    assert isinstance(panel, PanelResult)
    assert panel.n_verifiers == 3
    assert not panel.majority_refute
    assert panel.n_refuted < 2


@pytest.mark.slow
def test_panel_refutes_fabricated_impossible_bend() -> None:
    """FALSIFICATION E2E: a fabricated candidate with an impossible bend (a
    perturbed V_inf that violates the flyby turn limit) is majority-refuted ->
    the panel kills it. This is the gate's teeth."""
    ephem = Ephemeris("astropy")
    closure = _real_closure(ephem)

    # Fabricate: claim the closure is bend-feasible with an absurd V_inf cap that
    # the real geometry cannot satisfy, AND mark it as not actually converged.
    fabricated = BallisticClosureResult(
        t0_sec=closure.t0_sec,
        tof_days=closure.tof_days,
        max_residual_kms=50.0,  # nowhere near closed
        vinf_per_encounter_kms=tuple(v + 30.0 for v in closure.vinf_per_encounter_kms),
        converged=False,
        bend_feasible=False,
        vinf_cap_ok=False,
    )
    panel = adversarial_panel(
        fabricated,
        sequence=_SEQ,
        per_leg_revs=_REVS,
        per_leg_branch=_BRANCH,
        period_sec=_PERIOD_SEC,
        slack_leg=2,
        vinf_cap=8.0,  # the fabricated 30+ km/s V_inf busts this
        ephem=ephem,
        n_verifiers=3,
        seed=0,
    )
    assert panel.majority_refute
    assert panel.n_refuted >= 2


def test_panel_result_is_serialisable() -> None:
    """PanelResult exposes a dict for the review-queue audit trail."""
    pr = PanelResult(
        n_verifiers=3,
        n_refuted=0,
        majority_refute=False,
        verifier_verdicts=("confirmed", "confirmed", "confirmed"),
        notes="",
    )
    d = pr.as_dict()
    assert d["n_verifiers"] == 3
    assert d["majority_refute"] is False
