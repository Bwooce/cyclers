"""#826: bend feasibility + published-turn-ratio cross-check on #820's probe closures.

THE FINDING (no promotion; see docs/notes/2026-08-11-826-russell12-closure-adjudication.md).
------------------------------------------------------------------------------------------
#820 reported that ``probe_at_truth`` converges to a residual-zero closure at the
sourced geometry on 8 of 12 rows. Its ``stayed`` verdict is computed from
``BallisticClosureResult.converged`` ALONE — but the same dataclass also carries
``bend_feasible`` / ``vinf_cap_ok`` and a ``constraints_satisfied`` property
(= all three). Spec §14 **V0** requires the hard constraints (V_inf cap, r_p >=
r_p_min, **bend <= max**) on top of the closure residual, so "converged" is NOT
by itself even a V0-grade statement. Re-checking the 8 closures found **two are
bend-INFEASIBLE**: ``russell-ch4-5.30ggF3`` and ``russell-ch4-6.44Gg3``.

The two failures are NOT the same kind of thing, and the discriminator is a
SOURCED value — Russell's own published **turn ratio**:

    "Turn Ratio TR = max physically allowable turn angle / max required turn
     angle (delta_MAX)"; TR > 1 => all flybys physically attainable; max
     allowable based on a 200 km altitude Earth flyby.
    -- Russell-Ocampo 2003 p.13 (digest docs/notes/2026-06-17-digest-russell-ocampo-2003.md)

Our ``_max_bend_deg`` uses ``PLANETS[body].safe_alt_km = 200.0`` km for Earth and
Mars (itself sourced from Russell 2004 p.165 r_p,min), so ``min`` over flybys of
``max_bend / required_bend`` is a LIKE-FOR-LIKE reproduction of Russell's TR.

- ``russell-ch4-6.44Gg3`` measures TR = 0.956 against a published **0.95** — its
  bend infeasibility is a FAITHFUL reproduction of Russell's own Table 4.13
  NEAR-BALLISTIC classification (TR < 1: the cycler needs a small powered nudge),
  not a defect of the closure.
- ``russell-ch4-5.30ggF3`` measures TR = 0.613 against a published **1.27** — its
  third node reproduces the published value (1.255) but its SECOND node demands a
  141-degree turn where only ~87 degrees are available. That node is spurious: the
  closure is OFF-FAMILY there. This is the row #820 singled out as its strongest
  result (the sole STAYED-AT-TRUTH verdict, lowest truth residual 0.074) — the
  "it closed!" danger signal per [[feedback_orbit_closure_discipline]], caught only
  by a constraint the closure residual never sees.

GOLDEN DISCIPLINE: the EXPECTED side of every turn-ratio assertion below is the
row's own published ``invariants.turn_ratio`` (Russell 2004 Tables 4.9-4.13),
never a value this project computed. The bend angles EMERGE from the converged
Lambert chain.

Assertions record the qualitative finding with teeth (feasible/infeasible, and TR
agreement to a stated tolerance) without pinning brittle magnitudes. A fix under
``#829`` (gate the campaign on ``constraints_satisfied``) or ``#835``
(``5.30ggF3``'s off-family node) may legitimately change the ``5.30ggF3`` picture;
``6.44Gg3``'s TR < 1 is a published property and must NOT be "fixed" away.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import (
    _bend_deg,
    _max_bend_deg,
    _vinf_nodes,
    ballistic_correct,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
DAY_S = 86400.0
PHASE_EPOCHS = 64

# The 8 rows #820 reported as closing at truth. Values are #826's independent
# re-run (both 64 and 256 phase epochs agree).
BEND_FEASIBLE: dict[str, bool] = {
    "mcconaghy-2006-em-k2": True,
    "russell-ch4-4.991gG2": True,
    "russell-ch4-9.353Gg2": True,
    "russell-ch4-3.78Gg3": True,
    "russell-ch4-9.94Gg3": True,
    "russell-ch4-5.30ggF3": False,  # off-family node 2 (see module docstring)
    "russell-ch4-5.75ggF3": True,
    "russell-ch4-6.44Gg3": False,  # published TR = 0.95 < 1 (near-ballistic)
}

# Rows whose measured TR must reproduce the published TR. 5.30ggF3 is EXCLUDED
# (it is the off-family finding, asserted separately below).
TR_REPRODUCING = [r for r in BEND_FEASIBLE if r != "russell-ch4-5.30ggF3"]

# Tolerance on the TR reproduction. Russell prints TR to 2 decimals; the observed
# spread is 0.001-0.024, so 0.05 has teeth (it would catch a wrong node, a wrong
# flyby body, or a mis-posed genome) without pinning the exact value.
TR_TOL = 0.05


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _probe_closure(rid: str) -> tuple[Any, dict[str, np.ndarray], tuple[str, ...]]:
    """Replicate ``probe_at_truth``'s solve and return the full result + nodes.

    ``probe_at_truth`` returns a dict that drops ``bend_feasible`` entirely, so
    this reproduces its exact steps (select leg-1 topology at truth, then
    ``ballistic_correct`` seeded at the truth ToFs and best-phase t0) and keeps
    the ``BallisticClosureResult`` itself.
    """
    mod = _load_campaign()
    row = _row(rid)
    sel = mod.select_leg1_topology(
        mod.build_genome(row),
        model="circular",
        phase_epochs=PHASE_EPOCHS,
        t0_center=mod._t0_center(row),
    )
    genome = sel["genome"]
    ephem = Ephemeris("circular")
    solved = ballistic_correct(
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        t0_seed_sec=float(sel["best_t0_sec"]),
        tof_seed_days=mod._truth_seed(genome),
        period_sec=genome["period_sec"],
        ephem=ephem,
        vinf_cap=mod.VINF_CAP_KMS,
        slack_leg=genome["slack_leg"],
        tol_kms=mod.CORRECTOR_TOL_KMS,
        residual_mode="magnitude",
    )
    slack = genome["slack_leg"]
    nodes = _vinf_nodes(
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        t0_sec=float(solved.t0_sec),
        free_tof_days=tuple(float(v) for i, v in enumerate(solved.tof_days) if i != slack),
        slack_leg=slack,
        period_days=genome["period_sec"] / DAY_S,
        ephem=ephem,
    )
    return solved, nodes, genome["sequence"]


def _turn_ratios(nodes: dict[str, np.ndarray], sequence: tuple[str, ...]) -> list[float]:
    """Per-intermediate-flyby ``max_bend / required_bend`` (Russell's TR)."""
    out: list[float] = []
    for i in range(1, len(sequence) - 1):
        v_in, v_out = nodes[f"b{i}_in"], nodes[f"b{i}_out"]
        required = _bend_deg(v_in, v_out)
        max_turn = _max_bend_deg(float(np.linalg.norm(v_in)), sequence[i], None)
        # A node the closure turns through by ~0 degrees is unconstrained; it can
        # never be the binding flyby Russell's delta_MAX refers to.
        out.append(float("inf") if required < 1e-9 else max_turn / required)
    return out


@pytest.mark.parametrize("rid", sorted(BEND_FEASIBLE))
def test_probe_closure_bend_feasibility_is_recorded(rid: str) -> None:
    """#826: every #820 truth closure converges, but two are bend-INFEASIBLE.

    Guards the gap that ``probe_at_truth``'s ``stayed`` verdict (converged-only)
    hides: spec §14 V0 requires ``bend <= max`` as a hard constraint, so a
    converged-but-infeasible closure is not even V0-admissible evidence.
    """
    solved, _, _ = _probe_closure(rid)
    assert solved.converged, f"{rid}: #820's closure no longer converges"
    assert solved.vinf_cap_ok, f"{rid}: V_inf cap breached"
    assert solved.bend_feasible is BEND_FEASIBLE[rid], (
        f"{rid}: bend feasibility changed from #826's recorded finding "
        f"({BEND_FEASIBLE[rid]}); re-adjudicate before editing this expectation"
    )
    # The composite the campaign should have been gating on all along (#829).
    assert solved.constraints_satisfied is (
        solved.converged and solved.bend_feasible and solved.vinf_cap_ok
    )


@pytest.mark.parametrize("rid", sorted(TR_REPRODUCING))
def test_probe_closure_reproduces_published_turn_ratio(rid: str) -> None:
    """#826: the closure's binding flyby reproduces Russell's PUBLISHED turn ratio.

    Golden discipline: the expected side is the row's own sourced
    ``invariants.turn_ratio`` (Russell 2004 Tables 4.9-4.13, TR = max allowable /
    max required turn at a 200 km flyby); the bend angles emerge from the
    converged Lambert chain. This is an INDEPENDENT sourced check on the #820
    geometry — it constrains flyby turn angles, which the V_inf anchors do not.
    """
    published = _row(rid)["invariants"]["turn_ratio"]
    assert published is not None, f"{rid}: no published turn_ratio to check against"
    _, nodes, sequence = _probe_closure(rid)
    measured = min(_turn_ratios(nodes, sequence))
    assert abs(measured - float(published)) <= TR_TOL, (
        f"{rid}: measured TR {measured:.3f} vs published {published} "
        f"(tol {TR_TOL}) — the closure is not reproducing the sourced flyby geometry"
    )


def test_6_44gg3_bend_infeasibility_is_the_published_near_ballistic_property() -> None:
    """#826: 6.44Gg3's infeasible bend is Russell's own TR = 0.95, not a defect.

    Russell 2004 Table 4.13 is the NEAR-BALLISTIC table; this row is tagged
    "NEAR-BALLISTIC: TR = 0.95 < 1.0" (docs/notes/multi-arc-classification.md).
    Reproducing TR < 1 is therefore evidence the closure is ON the sourced family,
    even though it means the trajectory is not strictly ballistic. Do not "fix".
    """
    rid = "russell-ch4-6.44Gg3"
    published = float(_row(rid)["invariants"]["turn_ratio"])
    assert published < 1.0, "row's published TR is no longer sub-unity"
    solved, nodes, sequence = _probe_closure(rid)
    measured = min(_turn_ratios(nodes, sequence))
    assert not solved.bend_feasible
    assert measured < 1.0, "a TR < 1 row must measure an infeasible binding turn"
    assert abs(measured - published) <= TR_TOL


def test_5_30ggf3_closure_is_off_family_at_its_second_node() -> None:
    """#826: 5.30ggF3's closure demands a turn its published TR does not.

    #820's strongest-looking row (sole STAYED-AT-TRUTH, truth residual 0.074).
    Its LAST Earth node reproduces the published TR = 1.27, but an earlier Earth
    node demands roughly twice the available turn, dragging the binding ratio to
    ~0.61. A low closure residual with matching V_inf is therefore NOT sufficient
    evidence the sourced geometry was found. See #835.
    """
    rid = "russell-ch4-5.30ggF3"
    published = float(_row(rid)["invariants"]["turn_ratio"])
    solved, nodes, sequence = _probe_closure(rid)
    ratios = _turn_ratios(nodes, sequence)
    assert not solved.bend_feasible
    # The binding node is far BELOW the published TR — the off-family finding.
    assert min(ratios) < published - 0.5, (
        "5.30ggF3's binding turn ratio no longer contradicts its published TR; "
        "if a genome fix (#835) achieved this, re-adjudicate the row"
    )
    # ...while some node DOES reproduce it, which is why the residual looks good.
    assert any(abs(r - published) <= TR_TOL for r in ratios if r != float("inf")), (
        "no node reproduces the published TR any more"
    )
