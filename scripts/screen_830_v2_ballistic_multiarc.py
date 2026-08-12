"""#830 — V2-ballistic on #820's multi-arc genome, re-posed on the EXACT synodic period.

Question
--------
``tests/search/test_free_return_v2_ballistic.py`` declined V2-ballistic promotion
on a STRUCTURAL ground — a single-ellipse slice of a multi-arc cycler has no
continuous >=3-lap trajectory to propagate. #820's re-posed genome tiles the full
period (designated leg + every E-E phasing loop), so that ground is obsolete and
the spec §14 V2-ballistic gate (>=3 continuous laps, bounded rotating-frame
drift, in the row's defining model) is answerable for the first time.

Two prerequisites are honoured, both of which change the answer:

1. **Exact synodic period.** ``period.years`` is printed to 2 decimals; the
   rounding is +0.24 d (k=2) / -1.47 d (k=3), absorbed by the slack leg, and it
   is SECULAR — of order 3e6 km of rotating-frame drift per lap against a
   50,000 km tolerance. Every run below uses ``k * T_syn(pair)`` from the bodies'
   own mean motions. The catalogue-period result is measured alongside so the
   size of the rounding effect is reported, not asserted.
2. **Admissibility (#829).** Only rows whose seed-at-truth closure satisfies the
   spec §14 V0 hard constraints are eligible; a bend-infeasible closure cannot
   support a HIGHER tier.

THE GATE IS DEGENERATE FOR THIS CONSTRUCTION (measured, not assumed)
--------------------------------------------------------------------
``verify_long_term_stability`` reconstructs each lap from the cycler's LEG
TEMPLATE at lap-shifted planet positions (``propagate_lap`` step 3) rather than
integrating continuously across the wrap. With an exactly commensurate period on
the circular ephemeris the planets return to the same configuration, so ANY
template repeats itself and the lap-to-lap drift collapses — whether or not the
chain is a ballistic trajectory at all. A NEGATIVE CONTROL confirms it: breaking
the closure by +/-100 d of ToF (leaving a **35 km/s** V_inf discontinuity at a
flyby, i.e. a chain no spacecraft could fly) still measures ~1e-5 km of drift and
still reports ``stable=True``. A ``stable`` verdict here is therefore NOT evidence
of ballistic periodicity — it is a statement that the period is commensurate.
(``[[feedback_verify_gauntlet_with_positive_control]]``.)

By the same token the CATALOGUE-period drift (~5e6 km, ~100x tolerance) is not
evidence about the trajectory either: it is the period ROUNDING, i.e. #820's
predicted secular slack-leg error, now measured.

What IS load-bearing in an idealized-flyby model: once single-period closure
holds and every flyby's required turn is deliverable, the multi-lap trajectory is
the single lap repeated, so >=3-lap boundedness follows. That makes ALL-NODE bend
feasibility — INCLUDING the periodicity wrap node, which
``BallisticClosureResult.bend_feasible`` never checks — the real question, and it
is what this script reports per row.

This script COMPUTES; it decides nothing. Any validation-tier writeback is a
separate adjudication task (the #820->#826 / #822->#828 split).

Usage::

    uv run python scripts/screen_830_v2_ballistic_multiarc.py [--laps 3] [--out FILE]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct
from cyclerfinder.search.multiarc_cycler import build_multiarc_cycler
from cyclerfinder.search.turn_ratio_check import closure_turn_ratio
from cyclerfinder.verify.propagate import DRIFT_TOLERANCE_KM, verify_long_term_stability

REPO_ROOT = Path(__file__).resolve().parent.parent
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
CATALOGUE = REPO_ROOT / "data" / "catalogue.yaml"
DAY_S = 86400.0
PHASE_EPOCHS = 256

# The rows whose #829-gated seed-at-truth closure is ADMISSIBLE (spec §14 V0 hard
# constraints satisfied) under the corrected posing — the only V2 candidates.
# russell-ch4-6.44Gg3 is deliberately EXCLUDED: its closure is faithfully
# NEAR-BALLISTIC-AS-PUBLISHED (Russell Table 4.13, TR 0.95 < 1), so a
# *ballistic* gate is the wrong instrument for it.
ADMISSIBLE_ROWS = (
    "mcconaghy-2006-em-k2",
    "russell-ch4-4.991gG2",
    "russell-ch4-9.353Gg2",
    "russell-ch4-3.78Gg3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-5.30ggF3",
    "russell-ch4-5.75ggF3",
)


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def screen_row(
    mod: ModuleType, row: dict[str, Any], *, period_mode: str, n_laps: int, model: str = "circular"
) -> dict[str, Any]:
    """Close the row at truth under ``period_mode`` and run the V2 lap gate."""
    rid = row["id"]
    sel = mod.select_topology(
        mod.build_genome(row, period_mode=period_mode),
        model=model,
        phase_epochs=PHASE_EPOCHS,
        t0_center=mod._t0_center(row),
    )
    genome = sel["genome"]
    ephem = Ephemeris(model)
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
    out: dict[str, Any] = {
        "id": rid,
        "period_mode": period_mode,
        "period_days": round(genome["period_sec"] / DAY_S, 4),
        "loop_order": sel["loop_order"],
        "leg1": sel["leg1"],
        "truth_residual_kms": round(float(sel["best_truth_residual_kms"]), 4),
        "converged": bool(solved.converged),
        "solved_residual_kms": round(float(solved.max_residual_kms), 6),
        "constraints_satisfied": bool(solved.constraints_satisfied),
        "bend_feasible": bool(solved.bend_feasible),
        "vinf_cap_ok": bool(solved.vinf_cap_ok),
        "tof_days": [round(float(t), 3) for t in solved.tof_days],
        "vinf_per_encounter_kms": [round(float(v), 3) for v in solved.vinf_per_encounter_kms],
    }
    if not solved.converged:
        out["v2_ballistic"] = "NOT-CLOSED"
        return out

    tr = closure_turn_ratio(
        solved,
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        slack_leg=genome["slack_leg"],
        period_sec=genome["period_sec"],
        ephem=ephem,
    )
    out["measured_turn_ratio"] = None if tr.turn_ratio == float("inf") else round(tr.turn_ratio, 4)
    out["published_turn_ratio"] = (row.get("invariants") or {}).get("turn_ratio")
    out["turn_ratio_summary"] = tr.summary()

    # The PERIODICITY WRAP flyby, which BallisticClosureResult.bend_feasible
    # omits (it checks intermediate encounters only). For a MULTI-LAP question
    # this is a flyby the trajectory actually has to fly, so all-node
    # feasibility — not the drift number (see the module docstring / #830's
    # note) — is the load-bearing V2 evidence in an idealized-flyby model.
    with_wrap = closure_turn_ratio(
        solved,
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        slack_leg=genome["slack_leg"],
        period_sec=genome["period_sec"],
        ephem=ephem,
        include_wrap=True,
    )
    wrap = with_wrap.flybys[-1]
    out["wrap_node"] = {
        "required_bend_deg": round(wrap.required_bend_deg, 3),
        "max_bend_deg": round(wrap.max_bend_deg, 3),
        "ratio": None if wrap.unconstrained else round(wrap.ratio, 4),
        "feasible": bool(wrap.feasible),
    }
    out["all_nodes_feasible_incl_wrap"] = bool(with_wrap.all_feasible)
    out["turn_ratio_incl_wrap"] = (
        None if with_wrap.turn_ratio == float("inf") else round(with_wrap.turn_ratio, 4)
    )

    cycler = build_multiarc_cycler(
        solved,
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        period_sec=genome["period_sec"],
        ephem=ephem,
        sense=row.get("sense") or "n/a",
    )
    report = verify_long_term_stability(
        cycler, n_laps=n_laps, ephem=ephem, t_start=float(solved.t0_sec), cycler_id=rid
    )
    out.update(
        {
            "n_laps": int(report.n_laps_propagated),
            "max_drift_km": float(report.max_drift_km),
            "drift_tolerance_km": float(DRIFT_TOLERANCE_KM),
            "drift_over_tolerance": float(report.max_drift_km) / float(DRIFT_TOLERANCE_KM),
            "per_lap_drift_km": [float(d) for d in report.per_lap_drift_km],
            "stable": bool(report.stable),
            # NOT a delta-v: Cycler.maintenance_dv() is the raw sum of
            # |vinf_out - vinf_in|, which a ballistic flyby DELIVERS by turning.
            # Reported for completeness only; the ballistic question is the
            # per-node bend ratio above.
            "raw_vinf_rotation_sum_kms": float(cycler.maintenance_dv()),
            # DELIBERATELY NOT "PASS": the drift metric is degenerate for this
            # construction (see the module docstring's control result). The
            # honest per-row statement is the pair (drift, all-node feasibility).
            "drift_verdict": "BOUNDED" if report.stable else "UNBOUNDED",
            "v2_ballistic": "GATE-DEGENERATE",
        }
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laps", type=int, default=3)
    ap.add_argument("--rows", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    mod = _load_campaign()
    rows = yaml.safe_load(CATALOGUE.read_text())
    byid = {r["id"]: r for r in rows}
    ids = tuple(s.strip() for s in args.rows.split(",")) if args.rows else ADMISSIBLE_ROWS

    results: list[dict[str, Any]] = []
    for rid in ids:
        for period_mode in ("exact-synodic", "catalogue"):
            res = screen_row(mod, byid[rid], period_mode=period_mode, n_laps=args.laps)
            results.append(res)
            print(
                f"{rid:24s} {period_mode:14s} P={res['period_days']:.3f}d "
                f"conv={res['converged']} adm={res['constraints_satisfied']} "
                f"TR={res.get('measured_turn_ratio')}/{res.get('published_turn_ratio')} "
                f"drift={res.get('max_drift_km', float('nan')):.4g} km "
                f"({res.get('drift_over_tolerance', float('nan')):.4g}x tol) "
                f"{res.get('drift_verdict', '-')} "
                f"all_nodes_feasible={res.get('all_nodes_feasible_incl_wrap')} "
                f"wrap={res.get('wrap_node', {}).get('ratio')}",
                flush=True,
            )

    print("\n=== SUMMARY (>=3 laps, circular model; the drift gate is DEGENERATE) ===")
    counts: dict[str, int] = {}
    for r in results:
        key = f"{r['period_mode']}:{r.get('drift_verdict', r['v2_ballistic'])}"
        counts[key] = counts.get(key, 0) + 1
    feasible = sorted(
        r["id"]
        for r in results
        if r["period_mode"] == "exact-synodic" and r.get("all_nodes_feasible_incl_wrap")
    )
    print("drift counts:", counts)
    print("all-node feasible (incl. wrap), exact period:", feasible)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({"laps": args.laps, "results": results}, indent=2))
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
