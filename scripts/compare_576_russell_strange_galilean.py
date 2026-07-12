"""#576 step 5 -- Russell & Strange (2009) Table 3 literature golden comparison.

Mandatory per the #576 Fable plan review (C2-lit, load-bearing): Russell &
Strange, "Cycler Trajectories in Planetary Moon Systems," JGCD 32(1), 2009,
DOI 10.2514/1.36610, Table 3 ("Characteristics of promising Jovian ballistic
ideal model cyclers") is a SOURCED PUBLISHED census covering exactly this
territory -- transcribed here directly from the paper's OWN text layer
(``pdftotext -raw``, row order verified self-consistent against Table 2's
sourced synodic periods and the Fig. 6 caption's independently-legible
"BodyA-BodyB(v_inf_A#ID)" identifiers -- not hand-copied from the digest
note, which only paraphrased 3 of the 10 rows).

R-S ARCHITECTURE NOTE (why only "legs=1" rows are directly comparable)
------------------------------------------------------------------------
R-S's free-return cycler is built from "flyby-body free-return trajectories"
(Body A back to Body A), where "Number of legs" = how many such A-to-A
revolutions make up one full repeating cycle, and EXACTLY ONE of those legs
passes near the (massless-in-the-ideal-model) target Body B. For
``legs == 1``, the R-S cycler IS structurally identical to this project's
Anchor-Flyby-Anchor 2-leg symmetric construction (#563/#576): one A->B leg,
one B->A leg, closing at A. For ``legs > 1``, the R-S cycler makes multiple
A-to-A resonant loops with only a SINGLE B encounter total -- a genuinely
different topology (more like a resonant-hop + single-flyby tour) that this
project's 2-leg construction cannot represent and should not be compared to
on a "did we find the same object" basis. Both classes are reported below,
but only ``legs == 1`` rows are treated as apples-to-apples reproduction
targets.

CRITICAL ORDERING POINT (the #576 Fable correction, load-bearing)
--------------------------------------------------------------------
R-S's ideal-model search treats the TARGET body as massless -- it never
requires the target to bend at all. This project's #324 two-sided physical
gate requires EVERY node in the sequence (both the anchor's return bend and
the flyby/target's mid-tour bend) to clear 5 deg of useful bend. A genuine
R-S member could therefore fail this project's #324 gate on the target-body
side alone and never appear in ``enumerate_576_jupiter_galilean_symmetric_
closures.jsonl`` (which only records GATE-PASSING candidates), which would
misread a real R-S reproduction as a coverage failure -- exactly the #571
empty-region-stamp failure mode this project has hit before. So this script
compares against the UNGATED candidate set: every (pair, n, n_rev, rel_offset)
point that clears the RESIDUAL (closure) gate alone, re-derived here directly
from ``residual_at_point``/``GATE_RESIDUAL_KMS`` (step 3's own construction
machinery, imported verbatim) -- NOT read back from the enumerate_576 output
file, which only ever wrote the physically-gated survivors.

Comparison criterion: PERIOD MATCH. Every ungated symmetric candidate's cycle
period (``2 * tof_days = n * T_syn``) is compared against every R-S Table 3
Jovian row's ``Period, day`` column value, tolerance 0.5 days (loose enough
to absorb R-S's 1-decimal rounding and any patched-conic-vs-Lambert model
delta, tight enough that a spurious match at this system's short synodic
periods -- 2.35-12.52 d -- would be a coincidence, not noise).

Discipline: NO catalogue writeback, NO novelty claim (a period match here is
a POSITIVE reproduction result against a named published member -- it makes
the candidate MORE V0-known, not less; a non-match is reported honestly and
is NOT itself evidence of novelty -- literature_check.py + Opus/Fable
adjudication is the separate, later gate for that, per the #576 dispatch's
explicit scope limit).

Run as::

    uv run python scripts/compare_576_russell_strange_galilean.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from enumerate_563_symmetric_closures import (  # noqa: E402
    N_REV_VALUES,
    REL_OFFSETS_DEG,
    pair_n_max,
)
from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    gate_candidate,
    residual_at_point,
)

DATA_DIR = ROOT / "data"
OUT_PATH = DATA_DIR / "compare_576_russell_strange_galilean.jsonl"

PRIMARY = "Jupiter"
GALILEAN_MOONS: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")
TOF_SCALE_MAX = 2.0  # same bound as step 3 (verified against scan_433's own _meta record)

# Russell & Strange 2009, Table 3 ("Characteristics of promising Jovian
# ballistic ideal model cyclers"), transcribed verbatim from the paper's text
# layer (pdftotext -raw on the acquired PDF,
# /home/bruce/dev/cyclers_pdf/papers/russell-strange-2009-...pdf, page 7 of
# the extracted PDF / journal p.149), row order confirmed self-consistent
# against Table 2's sourced synodic periods (Io/Europa/Ganymede/Callisto mean
# orbital periods) and the independently-legible Fig. 6 caption identifiers.
# body_a = flyby body (bends every cycle); body_b = target body (passive,
# massless in R-S's ideal model); legs = number of body_a free-return
# revolutions per repeating cycle (ONE of which passes near body_b).
RUSSELL_STRANGE_TABLE3: tuple[dict[str, Any], ...] = (
    {
        "body_a": "Europa",
        "body_b": "Ganymede",
        "id": 131,
        "synodic_days": 7.05,
        "vinf_a_kms": 2.40,
        "vinf_b_kms": 4.10,
        "legs": 2,
        "period_days": 21.2,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Europa",
        "body_b": "Ganymede",
        "id": 159,
        "synodic_days": 7.05,
        "vinf_a_kms": 2.45,
        "vinf_b_kms": 4.11,
        "legs": 3,
        "period_days": 28.2,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Callisto",
        "id": 1,
        "synodic_days": 12.52,
        "vinf_a_kms": 3.18,
        "vinf_b_kms": 3.26,
        "legs": 3,
        "period_days": 37.6,
        "petal_yr": 0.41,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Callisto",
        "id": 5,
        "synodic_days": 12.52,
        "vinf_a_kms": 3.24,
        "vinf_b_kms": 3.34,
        "legs": 2,
        "period_days": 37.6,
        "petal_yr": 0.41,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Europa",
        "id": 5,
        "synodic_days": 7.05,
        "vinf_a_kms": 1.66,
        "vinf_b_kms": 2.57,
        "legs": 1,
        "period_days": 35.3,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Europa",
        "id": 43,
        "synodic_days": 7.05,
        "vinf_a_kms": 1.87,
        "vinf_b_kms": 3.89,
        "legs": 1,
        "period_days": 14.1,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Europa",
        "id": 316,
        "synodic_days": 7.05,
        "vinf_a_kms": 3.20,
        "vinf_b_kms": 3.81,
        "legs": 4,
        "period_days": 49.4,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Io",
        "id": 53,
        "synodic_days": 2.35,
        "vinf_a_kms": 3.90,
        "vinf_b_kms": 9.85,
        "legs": 2,
        "period_days": 21.2,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Io",
        "id": 185,
        "synodic_days": 2.35,
        "vinf_a_kms": 3.97,
        "vinf_b_kms": 9.90,
        "legs": 6,
        "period_days": 49.4,
        "petal_yr": 1.33,
    },
    {
        "body_a": "Ganymede",
        "body_b": "Io",
        "id": 403,
        "synodic_days": 2.35,
        "vinf_a_kms": 4.29,
        "vinf_b_kms": 4.34,
        "legs": 2,
        "period_days": 56.4,
        "petal_yr": 1.33,
    },
)

PERIOD_TOL_DAYS = 0.5


def enumerate_ungated_direction(anchor: str, flyby: str) -> list[dict[str, Any]]:
    """Every (n, n_rev, rel_offset) point that clears ONLY the residual
    (closure) gate -- the physical #324 bend gate is NOT applied here. This
    is the UNGATED set the #576 Fable correction mandates comparing R-S
    against, re-derived from step 3's own construction machinery verbatim
    (never read back from the gated-only enumerate_576 output file)."""
    t_syn, p_a, p_b, n_max = pair_n_max(anchor, flyby, primary=PRIMARY, tof_scale_max=TOF_SCALE_MAX)
    sqrt_papb = math.sqrt(p_a * p_b)
    ungated: list[dict[str, Any]] = []
    for n in range(1, n_max + 1):
        target_tof_days = n * t_syn / 2.0
        target_tof_scale = target_tof_days / sqrt_papb
        for n0, n1 in itertools.product(N_REV_VALUES, N_REV_VALUES):
            for rel in REL_OFFSETS_DEG:
                pt = residual_at_point(
                    anchor,
                    flyby,
                    rel_offset_deg=rel,
                    tof_scale=target_tof_scale,
                    n_rev=(n0, n1),
                    primary=PRIMARY,
                )
                if pt is None or pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                    continue
                gated = gate_candidate(anchor, flyby, pt, primary=PRIMARY)
                ungated.append(
                    {
                        "anchor": anchor,
                        "flyby": flyby,
                        "n_commensurate_int": n,
                        "t_syn_days": t_syn,
                        "rel_offset_deg": rel,
                        "n_rev": [n0, n1],
                        "tof_days": target_tof_days,
                        "cycle_period_days": 2.0 * target_tof_days,
                        "residual_kms": pt["residual_kms"],
                        "vinf_per_encounter_kms": gated["vinf_per_encounter_kms"],
                        "max_bend_deg_per_encounter": gated["max_bend_deg_per_encounter"],
                        "physical_gate_passed": gated["physical_gate_passed"],
                        "dop853_passed": gated["dop853_cross_check"]["passed"],
                        "all_gates_passed": gated["all_gates_passed"],
                    }
                )
    return ungated


def main() -> int:
    directions = []
    for a, b in itertools.combinations(GALILEAN_MOONS, 2):
        directions.append((a, b))
        directions.append((b, a))

    print(
        f"[576-RS] re-deriving UNGATED (residual-only) candidates for "
        f"{len(directions)} directions (tof_scale_max={TOF_SCALE_MAX})...",
        flush=True,
    )
    all_ungated: list[dict[str, Any]] = []
    for anchor, flyby in directions:
        pts = enumerate_ungated_direction(anchor, flyby)
        all_ungated.extend(pts)
        n_physical_pass = sum(1 for p in pts if p["all_gates_passed"])
        print(
            f"[576-RS]   {anchor}-{flyby}-{anchor}: {len(pts)} residual-only "
            f"candidates ({n_physical_pass} also clear the full physical gate)",
            flush=True,
        )
    print(f"[576-RS] {len(all_ungated)} total ungated (residual-only) candidates.", flush=True)

    # Compare against R-S Table 3: for each R-S row, is there an UNGATED
    # candidate on the SAME undirected pair whose cycle period matches
    # within tolerance? Report matches for ALL legs values (informative),
    # but flag legs==1 rows as the only architecturally comparable class.
    print("\n[576-RS] Russell-Strange Table 3 comparison (period match, tol=0.5d):", flush=True)
    rs_results = []
    n_legs1_matched = 0
    n_legs1_total = 0
    for rs in RUSSELL_STRANGE_TABLE3:
        pair = {rs["body_a"], rs["body_b"]}
        candidates_on_pair = [p for p in all_ungated if {p["anchor"], p["flyby"]} == pair]
        matches = [
            p
            for p in candidates_on_pair
            if abs(p["cycle_period_days"] - rs["period_days"]) <= PERIOD_TOL_DAYS
        ]
        is_legs1 = rs["legs"] == 1
        if is_legs1:
            n_legs1_total += 1
            if matches:
                n_legs1_matched += 1
        rec = {
            "rs_body_a": rs["body_a"],
            "rs_body_b": rs["body_b"],
            "rs_id": rs["id"],
            "rs_legs": rs["legs"],
            "rs_period_days": rs["period_days"],
            "rs_vinf_a_kms": rs["vinf_a_kms"],
            "rs_vinf_b_kms": rs["vinf_b_kms"],
            "architecturally_comparable": is_legs1,
            "n_ungated_candidates_on_pair": len(candidates_on_pair),
            "n_period_matches": len(matches),
            "matches": matches,
        }
        rs_results.append(rec)
        tag = "COMPARABLE" if is_legs1 else "different topology (legs>1)"
        match_str = f"{len(matches)} period match(es)" if matches else "no period match"
        print(
            f"[576-RS]   {rs['body_a']}-{rs['body_b']} #{rs['id']} (legs={rs['legs']}, "
            f"{tag}): R-S period={rs['period_days']:.1f}d, v_inf=({rs['vinf_a_kms']:.2f},"
            f"{rs['vinf_b_kms']:.2f}) km/s -> {match_str}",
            flush=True,
        )
        for m in matches:
            print(
                f"[576-RS]       MATCH: our n={m['n_commensurate_int']} n_rev={m['n_rev']} "
                f"rel_offset={m['rel_offset_deg']:.0f}deg period={m['cycle_period_days']:.3f}d "
                f"vinf={[f'{v:.2f}' for v in m['vinf_per_encounter_kms']]} km/s "
                f"physical_gate_passed={m['physical_gate_passed']}",
                flush=True,
            )

    print(
        f"\n[576-RS] DONE: {n_legs1_matched}/{n_legs1_total} architecturally-comparable "
        f"(legs=1) R-S Table 3 Jovian members have a period-matching ungated candidate in "
        "our own construction.",
        flush=True,
    )
    print(
        "[576-RS] Honest framing: a period match is a REPRODUCTION result (candidate is "
        "MORE V0-known), not evidence either way about novelty of any OTHER candidate. "
        "No novelty claim is made here for any candidate, matched or unmatched.",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#576 step 5 -- Russell-Strange 2009 Table 3 literature golden",
                    "source": "Russell & Strange, JGCD 32(1) 2009, DOI 10.2514/1.36610, "
                    "Table 3, transcribed via pdftotext -raw from the acquired PDF "
                    "(cyclers_pdf/papers), row order cross-validated against Table 2's "
                    "sourced synodic periods and the Fig.6 caption identifiers.",
                    "period_tol_days": PERIOD_TOL_DAYS,
                    "n_rs_rows": len(RUSSELL_STRANGE_TABLE3),
                    "n_rs_legs1_comparable": n_legs1_total,
                    "n_rs_legs1_matched": n_legs1_matched,
                    "n_ungated_candidates_total": len(all_ungated),
                }
            )
            + "\n"
        )
        for rec in rs_results:
            fh.write(json.dumps({"kind": "rs_row_comparison", **rec}) + "\n")
    print(f"[576-RS] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
