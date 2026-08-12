"""#833 — the measured-turn-ratio instrument (``search/turn_ratio_check.py``).

`#826` used ``min(max_bend / required_bend)`` as throwaway adjudication code and
it caught an off-family closure (``russell-ch4-5.30ggF3``) that both the closure
residual and the V_inf anchors passed. `#833` promotes it to a named, tested,
reusable instrument; these tests pin its mechanics.

GOLDEN DISCIPLINE: the integration test's EXPECTED side is the row's own
published ``invariants.turn_ratio`` (Russell 2004 Tables 4.9-4.13, TR = max
allowable / max required turn at a 200 km flyby). The unit tests below assert
only the instrument's own algebra against angles CONSTRUCTED in the test, never
against a value our physics code produced.
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
from cyclerfinder.search.correct import _max_bend_deg, ballistic_correct
from cyclerfinder.search.turn_ratio_check import (
    TURN_RATIO_TOL,
    closure_turn_ratio,
    measure_turn_ratio,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
PHASE_EPOCHS = 64


def _nodes_from_turns(vinf_kms: float, turns_deg: list[float]) -> dict[str, np.ndarray]:
    """Synthetic node vectors: ``len(turns_deg)`` intermediate flybys, each with
    the requested in->out turn at a fixed V_inf magnitude."""
    nodes: dict[str, np.ndarray] = {"b0_out": np.array([vinf_kms, 0.0, 0.0])}
    for i, turn in enumerate(turns_deg, start=1):
        nodes[f"b{i}_in"] = np.array([vinf_kms, 0.0, 0.0])
        t = np.radians(turn)
        nodes[f"b{i}_out"] = vinf_kms * np.array([np.cos(t), np.sin(t), 0.0])
    nodes[f"b{len(turns_deg) + 1}_in"] = np.array([vinf_kms, 0.0, 0.0])
    return nodes


def test_required_bend_and_ratio_reproduce_the_constructed_turn() -> None:
    """The measured required bend IS the constructed in->out angle, and the ratio
    is ``max_bend / required``."""
    seq = ("E", "E", "E")
    report = measure_turn_ratio(_nodes_from_turns(5.3, [30.0]), seq)
    (f,) = report.flybys
    assert f.body == "E" and f.index == 1
    assert f.required_bend_deg == pytest.approx(30.0, abs=1e-9)
    assert f.max_bend_deg == pytest.approx(_max_bend_deg(5.3, "E"), abs=1e-12)
    assert f.ratio == pytest.approx(f.max_bend_deg / 30.0, rel=1e-12)
    assert report.turn_ratio == pytest.approx(f.ratio, rel=1e-12)
    assert report.binding_index == 1


def test_only_intermediate_encounters_are_flybys() -> None:
    """The chain's ends carry a single leg each — they are joined by the
    periodicity residual, not by a turn, so they are never measured."""
    seq = ("E", "M", "E", "E")
    report = measure_turn_ratio(_nodes_from_turns(5.3, [10.0, 20.0]), seq)
    assert [f.index for f in report.flybys] == [1, 2]
    assert [f.body for f in report.flybys] == ["M", "E"]


def test_binding_node_is_the_minimum_ratio() -> None:
    """``turn_ratio`` is the MIN over flybys — the node that binds first."""
    seq = ("E", "E", "E", "E")
    report = measure_turn_ratio(_nodes_from_turns(5.3, [10.0, 80.0]), seq)
    assert report.binding_index == 2  # the bigger turn has the smaller ratio
    assert report.turn_ratio == pytest.approx(min(f.ratio for f in report.flybys), rel=1e-12)


def test_zero_turn_node_is_unconstrained_and_never_binds() -> None:
    """A node the closure turns through by ~0 deg imposes no turn at all, so it
    cannot be the binding flyby Russell's delta_MAX refers to.

    Under the `#820` designated-arc posing the Mars node is exactly this: legs 0
    and 1 are one conic split at the Mars encounter, so its required turn is 0.
    """
    seq = ("E", "M", "E", "E")
    report = measure_turn_ratio(_nodes_from_turns(5.3, [0.0, 40.0]), seq)
    mars, earth = report.flybys
    assert mars.unconstrained and mars.ratio == float("inf")
    assert not earth.unconstrained
    assert report.binding_index == 2
    assert "M1:inf" in report.summary()


def test_feasibility_is_ratio_ge_one() -> None:
    """``feasible`` <=> required <= max <=> ratio >= 1, and ``all_feasible``
    is the conjunction — the same statement as ``bend_feasible``."""
    seq = ("E", "E", "E")
    max_e = _max_bend_deg(5.3, "E")
    ok = measure_turn_ratio(_nodes_from_turns(5.3, [max_e * 0.5]), seq)
    assert ok.flybys[0].feasible and ok.all_feasible and ok.turn_ratio > 1.0
    bad = measure_turn_ratio(_nodes_from_turns(5.3, [max_e * 1.5]), seq)
    assert not bad.flybys[0].feasible and not bad.all_feasible and bad.turn_ratio < 1.0


def test_all_unconstrained_reports_infinite_ratio_and_no_binding_node() -> None:
    report = measure_turn_ratio(_nodes_from_turns(5.3, [0.0]), ("E", "E", "E"))
    assert report.turn_ratio == float("inf")
    assert report.binding_index is None
    assert report.all_feasible


def test_agrees_with_published_uses_the_documented_tolerance() -> None:
    """The published-TR comparison is an absolute-difference test at
    :data:`TURN_RATIO_TOL` (0.05; `#826` measured a 0.001-0.024 spread)."""
    seq = ("E", "E", "E")
    max_e = _max_bend_deg(5.3, "E")
    report = measure_turn_ratio(_nodes_from_turns(5.3, [max_e / 1.27]), seq)
    assert report.turn_ratio == pytest.approx(1.27, abs=1e-9)
    assert report.agrees_with_published(1.27)
    assert report.agrees_with_published(1.27 + 0.9 * TURN_RATIO_TOL)
    assert not report.agrees_with_published(1.27 + 1.1 * TURN_RATIO_TOL)


# ---------------------------------------------------------------------------
# Integration: the instrument on a real converged closure.
# ---------------------------------------------------------------------------


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def test_closure_turn_ratio_matches_bend_feasible_and_published_tr() -> None:
    """On a real `#820` closure the instrument must (a) agree with the
    corrector's OWN ``bend_feasible`` verdict — a cross-check that it is reading
    the same nodes — and (b) reproduce the row's PUBLISHED turn ratio.

    ``russell-ch4-9.353Gg2`` is one of `#820`'s two CLOSE-AND-MATCH rows and is
    bend-feasible, so it exercises both halves at once.
    """
    rid = "russell-ch4-9.353Gg2"
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
    assert solved.converged
    report = closure_turn_ratio(
        solved,
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        slack_leg=genome["slack_leg"],
        period_sec=genome["period_sec"],
        ephem=ephem,
    )
    assert report.all_feasible is solved.bend_feasible
    published = float(row["invariants"]["turn_ratio"])
    assert report.agrees_with_published(published), (
        f"{rid}: measured TR {report.turn_ratio:.3f} vs published {published} ({report.summary()})"
    )
