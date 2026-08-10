"""#818 -- passive-node self-consistency gate, re-run over #816's stored roots.

PURE POST-PROCESSING: reads the already-stored per-root numbers from
``data/found/816_unequal_tof_asymmetric_roots/roots.json`` and applies the new
:func:`cyclerfinder.search.physical_sanity.passive_node_is_self_consistent`
gate to every ``anchor_passive`` / ``flyby_passive`` root. NO physics re-run,
no Lambert solves, no propagation — the gate is a closed-form function of the
stored (body, V_inf, required-turn) triple.

Purpose (from the #818 registration in data/OUTSTANDING.md): test #817's
explicitly-stated-as-unverified hypothesis that ALL 785 passive-node roots
fail the same way its hand-adjudicated Ariel-Oberon object did (parasitic
patched-conic deflection at the passive body's own SOI boundary is a
non-negligible fraction of the working turn budget).

Per-root inputs (stored by scan_816_unequal_tof_discrete_roots.py):

* ``class``: ``flyby_passive`` -> passive body is ``flyby``, working budget is
  ``required_turn_anchor_deg``; ``anchor_passive`` -> passive body is
  ``anchor``, working budget is ``required_turn_flyby_deg``.
* ``vinf_per_encounter_kms`` is ``[anchor, flyby, anchor]`` order (verified
  against the stored max_bend_deg_per_encounter values). For an anchor-passive
  root the two anchor entries differ only at FP noise; the LARGER V_inf is
  used (higher V_inf -> SMALLER parasitic deflection -> the most lenient
  bound, so a rejection is robust).

Output: ``data/found/818_passive_node_gate/gate_results.json`` — per-root gate
records + summary counts (including threshold-sensitivity counts at 5% / 10%
and at the more lenient Hill-radius bound), plus a console summary.

Usage::

    uv run python scripts/postprocess_818_passive_node_gate.py
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any

from cyclerfinder.search.physical_sanity import (
    DEFAULT_MAX_PARASITIC_TURN_FRACTION,
    PASSIVE_NODE_TURN_MAX_DEG,
    passive_node_is_self_consistent,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ROOTS = REPO_ROOT / "data/found/816_unequal_tof_asymmetric_roots/roots.json"
DEFAULT_OUT = REPO_ROOT / "data/found/818_passive_node_gate/gate_results.json"

PASSIVE_CLASSES = ("anchor_passive", "flyby_passive")


def passive_node_inputs(root: dict[str, Any]) -> tuple[str, float, float, float]:
    """Extract ``(passive_body, vinf_kms, working_budget_deg, passive_turn_deg)``.

    Pure function of the stored per-root record (see module docstring for the
    conventions). Raises ``ValueError`` on a non-passive class.
    """
    vinfs = root["vinf_per_encounter_kms"]
    if root["class"] == "flyby_passive":
        return (
            str(root["flyby"]),
            float(vinfs[1]),
            float(root["required_turn_anchor_deg"]),
            float(root["required_turn_flyby_deg"]),
        )
    if root["class"] == "anchor_passive":
        # Lenient bound: the larger of the two stored anchor V_infs.
        return (
            str(root["anchor"]),
            max(float(vinfs[0]), float(vinfs[2])),
            float(root["required_turn_flyby_deg"]),
            float(root["required_turn_anchor_deg"]),
        )
    raise ValueError(f"not a passive-node root: class={root['class']!r}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--roots", type=Path, default=DEFAULT_ROOTS)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    payload = json.loads(args.roots.read_text())
    passive = [r for r in payload["roots"] if r["class"] in PASSIVE_CLASSES]
    print(f"{len(payload['roots'])} stored roots, {len(passive)} passive-node roots")

    records: list[dict[str, Any]] = []
    for root in passive:
        body, vinf, budget, passive_turn = passive_node_inputs(root)
        # Consistency with the #816 classifier's own passivity criterion.
        assert passive_turn < PASSIVE_NODE_TURN_MAX_DEG, (passive_turn, root)
        verdict = passive_node_is_self_consistent(body, vinf, budget)
        frac_hill = verdict.parasitic_deflection_hill_deg / budget if budget > 0.0 else float("inf")
        records.append(
            {
                "anchor": root["anchor"],
                "flyby": root["flyby"],
                "q": root["q"],
                "n_rev": root["n_rev"],
                "beta_deg": root["beta_deg"],
                "class": root["class"],
                "passive_turn_deg": passive_turn,
                "turn_feasible": root["turn_feasible"],
                "passes_physical_gates": root["passes_physical_gates"],
                "parasitic_fraction_hill": frac_hill,
                "verdict": asdict(verdict),
            }
        )

    # --- summary -----------------------------------------------------------
    n = len(records)
    n_reject = sum(1 for r in records if not r["verdict"]["is_self_consistent"])
    fracs = [r["verdict"]["parasitic_fraction"] for r in records]
    tab: Counter[tuple[bool, bool, bool]] = Counter(
        (
            bool(r["turn_feasible"]),
            bool(r["passes_physical_gates"]),
            bool(r["verdict"]["is_self_consistent"]),
        )
        for r in records
    )
    feasible = [r for r in records if r["turn_feasible"]]
    feasible_rejected = [r for r in feasible if not r["verdict"]["is_self_consistent"]]
    # The joint battery: a passive root survives only if it passes BOTH the
    # #816 physical gates AND this gate.
    joint_survivors = [
        r for r in records if r["passes_physical_gates"] and r["verdict"]["is_self_consistent"]
    ]
    sensitivity = {
        f"n_reject_at_{int(100 * t)}pct": sum(1 for f in fracs if f > t) for t in (0.02, 0.05, 0.10)
    }
    sensitivity["n_reject_at_2pct_hill_bound"] = sum(
        1 for r in records if r["parasitic_fraction_hill"] > DEFAULT_MAX_PARASITIC_TURN_FRACTION
    )
    # Complementary ABSOLUTE reading (reported, not verdict-bearing): the
    # closure modelled the node as turning < PASSIVE_NODE_TURN_MAX_DEG, so if
    # the body's unavoidable SOI-boundary deflection exceeds even that, the
    # stored root is self-inconsistent regardless of budget.
    n_parasitic_exceeds_passivity = sum(
        1 for r in records if r["verdict"]["parasitic_deflection_deg"] > PASSIVE_NODE_TURN_MAX_DEG
    )
    parasitic_degs = [r["verdict"]["parasitic_deflection_deg"] for r in records]

    summary = {
        "task": 818,
        "source": str(args.roots.relative_to(REPO_ROOT)),
        "threshold_max_parasitic_fraction": DEFAULT_MAX_PARASITIC_TURN_FRACTION,
        "n_passive_roots": n,
        "n_rejected_by_gate": n_reject,
        "n_admitted_by_gate": n - n_reject,
        "parasitic_fraction_min": min(fracs),
        "parasitic_fraction_max": max(fracs),
        "n_turn_feasible": len(feasible),
        "n_turn_feasible_rejected_by_gate": len(feasible_rejected),
        "n_joint_survivors_physical_and_818": len(joint_survivors),
        "crosstab_turnfeasible_physgates_selfconsistent": {
            f"turn_feasible={k[0]},phys={k[1]},consistent={k[2]}": v for k, v in sorted(tab.items())
        },
        "threshold_sensitivity": sensitivity,
        "parasitic_deflection_deg_min": min(parasitic_degs),
        "parasitic_deflection_deg_max": max(parasitic_degs),
        "n_parasitic_exceeds_passivity_threshold_absolute": n_parasitic_exceeds_passivity,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"summary": summary, "records": records}, indent=1))
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {args.out}")

    print("\nturn-feasible passive roots (the only ones any other gate could not kill):")
    for r in feasible:
        v = r["verdict"]
        print(
            f"  {r['anchor']}-{r['flyby']} q={r['q']} n_rev={tuple(r['n_rev'])} "
            f"beta={r['beta_deg']:.4f} class={r['class']}: passive {v['body']} at "
            f"{v['vinf_kms']:.4f} km/s -> parasitic {v['parasitic_deflection_deg']:.4f} deg "
            f"= {100 * v['parasitic_fraction']:.2f}% of {v['working_turn_budget_deg']:.4f} deg "
            f"budget -> {'ADMIT' if v['is_self_consistent'] else 'REJECT'}"
        )


if __name__ == "__main__":
    main()
