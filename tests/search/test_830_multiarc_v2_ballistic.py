"""#830 — multi-arc V2-ballistic: the structural obstacle is gone, the GATE is degenerate.

FINDING (no promotion; see docs/notes/2026-08-12-829-833-835-830-russell12-bend-gate.md)
-----------------------------------------------------------------------------------------
``tests/search/test_free_return_v2_ballistic.py`` declined V2-ballistic on the
STRUCTURAL ground that a single-ellipse slice of a multi-arc cycler gives "no
continuous >=3-lap trajectory to propagate". #820's re-posed genome tiles the full
period, and :func:`cyclerfinder.search.multiarc_cycler.build_multiarc_cycler` now
feeds it to the §12 machinery — so that ground IS obsolete.

But the gate that machinery implements cannot decide the question for this
construction, and this module pins WHY with a negative control:
``verify_long_term_stability`` rebuilds each lap from the cycler's leg TEMPLATE at
lap-shifted planet positions instead of integrating across the wrap, so on an
exactly commensurate period with the circular ephemeris ANY template repeats and
the drift collapses to ~1e-5 km. Breaking the closure — leaving a multi-km/s
V_inf discontinuity at a flyby, a chain nothing could fly — still measures the
same ~1e-5 km and still reports ``stable=True``.

A ``stable`` verdict is therefore a statement about PERIOD COMMENSURABILITY, not
about ballistic periodicity, and must never be cited as V2 evidence for this
construction. Equally, the ~1e7 km drift measured on the CATALOGUE (rounded)
period is the period rounding, not a property of the trajectory.

``[[feedback_verify_gauntlet_with_positive_control]]``: "0/N all-fail" and "it
passed!" are the same danger; run the control and judge by the right criterion.

Not marked ``slow``: a V-gauntlet evidence test that CI skips is an unverified
claim (``[[feedback_delegation_fresh_agent_not_fork]]``).
"""

from __future__ import annotations

import dataclasses
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct
from cyclerfinder.search.multiarc_cycler import build_multiarc_cycler
from cyclerfinder.search.turn_ratio_check import closure_turn_ratio
from cyclerfinder.verify.propagate import DRIFT_TOLERANCE_KM, verify_long_term_stability

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
PHASE_EPOCHS = 64
DAY_S = 86400.0


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _closed(rid: str, period_mode: str = "exact-synodic") -> tuple[Any, dict[str, Any], Ephemeris]:
    mod = _load_campaign()
    sel = mod.select_topology(
        mod.build_genome(_row(rid), period_mode=period_mode),
        model="circular",
        phase_epochs=PHASE_EPOCHS,
        t0_center=mod._t0_center(_row(rid)),
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
    return solved, genome, ephem


def _cycler(solved: Any, genome: dict[str, Any], ephem: Ephemeris) -> Any:
    return build_multiarc_cycler(
        solved,
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        period_sec=genome["period_sec"],
        ephem=ephem,
    )


def test_exact_synodic_period_differs_from_the_printed_one_by_a_secular_amount() -> None:
    """The printed ``period.years`` is rounded to 2 dp; on the k=3 rows that is
    ~1.5 d, absorbed by the slack leg and secular over laps."""
    mod = _load_campaign()
    for rid, min_delta_days in (("russell-ch4-9.353Gg2", 0.2), ("russell-ch4-5.30ggF3", 1.4)):
        row = _row(rid)
        cat = mod.period_sec_for_row(row, "catalogue")
        exact = mod.period_sec_for_row(row, "exact-synodic")
        assert abs(cat - exact) / DAY_S >= min_delta_days
        # ...and the exact value really is k * T_syn of the row's own pair.
        assert exact / DAY_S == pytest.approx(
            float(row["period"]["k"]) * 779.9286472363009, rel=1e-9
        )


def test_multiarc_adapter_reconstructs_the_closure() -> None:
    """The adapter must reproduce the closure it is handed: one leg per ToF, one
    encounter per body, and V_inf MAGNITUDE continuity at every flyby (which is
    exactly what the corrector's residual drove to zero)."""
    solved, genome, ephem = _closed("russell-ch4-5.30ggF3")
    assert solved.converged
    cycler = _cycler(solved, genome, ephem)
    assert len(cycler.legs) == len(solved.tof_days)
    assert len(cycler.encounters) == len(genome["sequence"])
    assert cycler.period == pytest.approx(genome["period_sec"])
    for enc in cycler.encounters:
        assert float(np.linalg.norm(enc.vinf_in)) == pytest.approx(
            float(np.linalg.norm(enc.vinf_out)), abs=1e-6
        )
    # Boundary convention: the wrap endpoints carry vinf_in == vinf_out, so the
    # raw rotation sum does not charge the frame rotation as a manoeuvre.
    for enc in (cycler.encounters[0], cycler.encounters[-1]):
        assert np.allclose(enc.vinf_in, enc.vinf_out)


def test_lap_drift_gate_is_degenerate_on_a_commensurate_period() -> None:
    """NEGATIVE CONTROL — the load-bearing result of #830.

    A chain broken by 100 d of ToF (multi-km/s V_inf discontinuity at a flyby)
    measures the SAME lap-to-lap drift as the converged closure and is likewise
    reported ``stable``. So ``stable`` here cannot be cited as V2-ballistic
    evidence: it tests period commensurability, not the trajectory.
    """
    solved, genome, ephem = _closed("russell-ch4-9.353Gg2")
    assert solved.converged

    def drift_of(res: Any) -> float:
        report = verify_long_term_stability(
            _cycler(res, genome, ephem),
            n_laps=3,
            ephem=ephem,
            t_start=float(res.t0_sec),
            cycler_id="control",
        )
        assert report.n_laps_propagated == 3
        assert report.stable, "the metric is expected to report stable even when broken"
        return float(report.max_drift_km)

    good = drift_of(solved)
    assert good < DRIFT_TOLERANCE_KM

    n = len(solved.tof_days)
    tofs = list(solved.tof_days)
    tofs[n - 2] += 100.0
    tofs[n - 1] -= 100.0  # period preserved; the CHAIN is destroyed
    broken = dataclasses.replace(solved, tof_days=tuple(tofs))
    broken_cycler = _cycler(broken, genome, ephem)
    discontinuity = max(
        abs(float(np.linalg.norm(e.vinf_in)) - float(np.linalg.norm(e.vinf_out)))
        for e in broken_cycler.encounters
    )
    assert discontinuity > 1.0, "the control must actually break V_inf continuity"
    bad = drift_of(broken)
    assert bad < DRIFT_TOLERANCE_KM, "an unflyable chain also passes — the gate is degenerate"
    # Same order of magnitude: the metric does not even resolve the difference.
    assert 0.1 < bad / good < 10.0


def test_all_nodes_including_the_wrap_are_feasible_under_the_exact_period() -> None:
    """What IS load-bearing for a multi-lap claim in an idealized-flyby model.

    ``BallisticClosureResult.bend_feasible`` checks the INTERMEDIATE encounters
    only; the periodicity wrap is a real flyby the trajectory must fly every lap.
    Under the exact-synodic posing both #830 candidates measured here are
    feasible at every node INCLUDING the wrap, and the binding ratio still
    reproduces the row's PUBLISHED turn ratio (golden: the expected side is the
    published value).
    """
    for rid in ("russell-ch4-9.353Gg2", "russell-ch4-5.30ggF3"):
        solved, genome, ephem = _closed(rid)
        assert solved.converged and solved.constraints_satisfied
        report = closure_turn_ratio(
            solved,
            sequence=genome["sequence"],
            per_leg_revs=genome["per_leg_revs"],
            per_leg_branch=genome["per_leg_branch"],
            slack_leg=genome["slack_leg"],
            period_sec=genome["period_sec"],
            ephem=ephem,
            include_wrap=True,
        )
        assert report.all_feasible, f"{rid}: {report.summary()}"
        published = float(_row(rid)["invariants"]["turn_ratio"])
        assert report.agrees_with_published(published), (
            f"{rid}: measured {report.turn_ratio:.4f} vs published {published}"
        )
