"""#1047: the literature step for the #1000/#1048 Earth-Moon families with the v3 gate
(`literature_check.check_literature`, offline corpus backend), using LABELLED and
TOPOLOGIED signatures only (#1045 ruling), controls first.

Signature for a family: primary "Earth", sequence ("Moon",), resonances ("p:q",), period_k = q,
period_band_tu = (T_min, T_max) of its computed members, topology_label {"resonant"} (a p:q
Earth-resonant periodic orbit with one lunar pass per period), working_bodies "one" (the Moon is the
only flyby body; the Earth is the primary). A sensitivity run repeats every signature with
topology {"repeated-moon"} (the Ross & Roberts-Tsoukkas cycler class).

Controls (run first; a failed control stops the run):
- positive: casoliva-7-3b (7:3, T 18.8496), labelled and topologied as above, must return
  "published" citing an Earth-Moon anchor (Casoliva 2010 or another Moon-body Earth anchor);
- negative (deliberately mislabelled): the same 7-3b signature with primary "Jupiter" (wrong
  system), and with topology {"pump-tour"} (a tour label), must NOT return "published".
Output: data/1047_literature/gate_results.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    _candidate_anchors,
    check_literature,
    offline_corpus_search,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1047_literal as lit

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1047_literature" / "gate_results.json"
EM_KEYS = ("Earth-Moon", "Casoliva", "Ross", "Braik", "Kumar", "Rawat", "Rosengren", "Antoniadou")


def sig(
    pq: tuple[int, int],
    t_band: tuple[float, float],
    topology: str,
    primary: str = "Earth",
) -> CandidateSignature:
    return CandidateSignature(
        primary=primary,
        sequence=("Moon",),
        period_k=pq[1],
        resonances=(f"{pq[0]}:{pq[1]}",),
        period_band_tu=t_band,
        topology_label=frozenset({topology}),
        working_bodies="one",
    )


def run(s: CandidateSignature) -> dict[str, Any]:
    r = check_literature(s, search=offline_corpus_search)
    return {
        "status": r.status,
        "citation": r.citation,
        "confidence": r.confidence,
        "notes": r.notes,
    }


def cites_em(res: dict[str, Any]) -> bool:
    return any(k in (res.get("citation") or "") for k in EM_KEYS)


def main() -> None:
    preflight_search(
        task_no=1047,
        region_id="literature-gate-v3-earth-moon-families",
        method=MethodCapability(
            genome="literature_check v3, offline corpus, labelled + topologied signatures",
            corrector="none (lookup)",
            capability_tags=frozenset({"literature"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    t73 = 18.8495559215
    controls = {
        "positive_casoliva_7-3b": run(sig((7, 3), (t73 * 0.999, t73 * 1.001), "resonant")),
        "positive_casoliva_7-3b_repeated-moon": run(
            sig((7, 3), (t73 * 0.999, t73 * 1.001), "repeated-moon")
        ),
        "negative_wrong_system_jupiter": run(
            sig((7, 3), (t73 * 0.999, t73 * 1.001), "resonant", primary="Jupiter")
        ),
        "negative_tour_label": run(sig((7, 3), (t73 * 0.999, t73 * 1.001), "pump-tour")),
    }
    for k, v in controls.items():
        print(k, v["status"], v["citation"], round(v["confidence"], 3), flush=True)
    pos_ok = controls["positive_casoliva_7-3b"]["status"] == "published" and cites_em(
        controls["positive_casoliva_7-3b"]
    )
    neg1 = controls["negative_wrong_system_jupiter"]
    neg2 = controls["negative_tour_label"]
    # The lead's criterion (2026-10-11): a deliberately mislabelled signature must NOT return
    # "published" (any citation). The first run used a looser test (tour label must not cite
    # Casoliva) and passed a control that fails this criterion; see the #1047 note sec. 3.
    neg_ok = neg1["status"] != "published" and neg2["status"] != "published"
    s73 = sig((7, 3), (t73 * 0.999, t73 * 1.001), "pump-tour")
    unlabelled = [a.name for a in _candidate_anchors(s73) if not a.topology_label]
    out: dict[str, Any] = {
        "controls": controls,
        "positive_pass": bool(pos_ok),
        "negative_pass": bool(neg_ok),
        "controls_pass": bool(pos_ok and neg_ok),
        "earth_moon_anchors_in_scope_with_no_topology_label": unlabelled,
    }
    if not out["controls_pass"]:
        OUT.write_text(json.dumps(out, indent=1) + "\n")
        print("CONTROLS FAILED: stop", flush=True)
        return
    fams = lit.families()
    res: dict[str, Any] = {}
    for key, f in sorted(fams.items()):
        ts = [m["T"] for m in f["members"]]
        band = (min(ts) * 0.999, max(ts) * 1.001)
        pq = tuple(f["pq"])
        res[key] = {
            "computed_pq": pq,
            "T_band": band,
            "resonant": run(sig(pq, band, "resonant")),
            "repeated-moon": run(sig(pq, band, "repeated-moon")),
        }
        r = res[key]["resonant"]
        print(f"{key}: {r['status']} {r['citation']} {r['confidence']:.3f}", flush=True)
    out["families"] = res
    OUT.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
