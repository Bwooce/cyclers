"""#942/#943: literature_check with the declared-scope tags — positive controls and candidates.

Pre-registration: docs/notes/2026-10-07-942-943-literature-gate-scope-preregistration.md. Every
signature is built the same way: the encounter sequence in time order, the period k, the V_inf per
encounter, topology {"repeated-moon"}, and the two MECHANICAL labels:
- working_bodies: "two" if every body has a demanded turn >= 0.05 deg at some encounter, else
  "one" (amendment A1: half the sec. 6.1 0.1-deg turn resolution);
- return_types: per same-body leg of the cycle key: R<b>/1:1 -> FR, other R -> FR-n:m, H -> HR,
  L<b>><b> with flight time / body period in (1, 2) -> SY (Menning's symmetric return), else GEN.
The catalogued H&M rows have no cycle key, so their labels are left undeclared.

Usage: uv run python scripts/litcheck_942_943_scope.py
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.search.literature_check import (
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)
from cyclerfinder.search.two_working_body import (
    HalfRevLeg,
    LambertLeg,
    ResonantLeg,
    correct_dates,
    eval_lambert_legs,
)
from cyclerfinder.search.two_working_body_enum import Zero, assess, flyby_table

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
#: A body "works" if it turns by at least this somewhere. Amendment A1 of the pre-registration
#: (after the EurGan#131 control failed at 1e-6 deg on a 1.2e-6-deg round-off turn): half the
#: 0.1-deg turn resolution of the results note's sec. 6.1 dedupe rule, which defines "turn 0".
TURN_EPS_DEG = 0.05


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")


def labels(cell: str, key: str, x_days: list[float]) -> tuple[Any, ...]:
    """(primary, sequence, k, vinf per encounter, working_bodies, return_types)."""
    system, a, b = ENUM.cell_system(cell)
    k, cyc = ENUM.parse_cycle_key(key, system, a, b)
    x = np.asarray(x_days) * DAY
    ass = assess(system, Zero(cyc, x, 0.0))
    table = sorted(flyby_table(system, ass), key=lambda f: float(f["t_days"]))
    seq = tuple(str(f["body"]) for f in table)
    vinf = tuple(round(float(f["vinf_kms"]), 3) for f in table)
    working = (
        "two"
        if all(
            any(float(f["turn_deg"]) >= TURN_EPS_DEG for f in table if f["body"] == c)
            for c in (a, b)
        )
        else "one"
    )
    legs = eval_lambert_legs(system, cyc, x)
    assert legs is not None
    rtypes: set[str] = set()
    lam_i = 0
    for leg in cyc.legs:
        if isinstance(leg, ResonantLeg):
            rtypes.add(
                "FR"
                if (leg.body_revs, leg.sc_revs) == (1, 1)
                else f"FR-{leg.body_revs}:{leg.sc_revs}"
            )
        elif isinstance(leg, HalfRevLeg):
            rtypes.add("HR")
        elif isinstance(leg, LambertLeg):
            ev = legs[lam_i]
            lam_i += 1
            if leg.frm == leg.to:
                ratio = (ev.t_arr - ev.t_dep) / system.period_s(leg.frm)
                rtypes.add("SY" if leg.nrev == 1 and 1.0 < ratio < 2.0 else "GEN")
    primary = "Jupiter" if cell in ("gc", "ge") else "Sun"
    return primary, seq, k, vinf, working, frozenset(rtypes)


def run(name: str, sig: CandidateSignature) -> dict[str, Any]:
    r = check_literature(sig, search=offline_corpus_search)
    rec = {
        "name": name,
        "sequence": list(sig.sequence),
        "working_bodies": sig.working_bodies,
        "return_types": sorted(sig.return_types) if sig.return_types is not None else None,
        "status": r.status,
        "confidence": r.confidence,
        "citation": r.citation,
        "notes": r.notes,
    }
    print(
        f"{name}: {r.status} ({r.confidence}) labels={sig.working_bodies} "
        f"{rec['return_types']} cit={str(r.citation)[:70]}"
    )
    return rec


def from_gauntlet(
    path: str, pick: dict[str, float], key: str | None = None
) -> tuple[str, str, list[float]]:
    g = json.loads((REPO / path).read_text())
    rows = [
        c
        for c in g["candidates"]
        if all(abs(c["vinf_kms"][q] - v) < 0.01 for q, v in pick.items())
        and (key is None or c["key"] == key)
    ]
    assert len(rows) == 1, (path, pick, len(rows))
    return g["cell"], rows[0]["key"], rows[0]["x_days"]


def sig_of(cell: str, key: str, x_days: list[float]) -> CandidateSignature:
    primary, seq, k, vinf, working, rtypes = labels(cell, key, x_days)
    return CandidateSignature(
        primary=primary,
        sequence=seq,
        period_k=k,
        vinf_per_encounter_kms=vinf,
        topology_label=frozenset({"repeated-moon"}),
        working_bodies=working,
        return_types=rtypes,
    )


def main() -> None:
    out: dict[str, list[dict[str, Any]]] = {"controls": [], "candidates": []}
    controls: list[tuple[str, str, str, list[float]]] = []
    gc = "data/943_cell_gc_gauntlet.json"
    ge = "data/943_cell_ge_gauntlet.json"
    ev = "data/942_cell_ev_gauntlet.json"
    controls.append(("GanCal#5", *from_gauntlet(gc, {"Ganymede": 3.238})))
    # GanCal#1: re-solved from its sec. 5 dates in cell gc (not in the gc cell's scope).
    s_gc, a, b = ENUM.cell_system("gc")
    key1 = "k3|LGanymede>Ganymede/1l|RGanymede/2:1|LGanymede>Callisto/0s|LCallisto>Ganymede/0s"
    _, cyc1 = ENUM.parse_cycle_key(key1, s_gc, a, b)
    sol = correct_dates(s_gc, cyc1, np.array([0.991426, 26.049946, 35.954310]) * DAY, tol_kms=1e-8)
    controls.append(("GanCal#1", "gc", key1, [float(v) / DAY for v in sol.x]))
    controls.append(("GanEur#43", *from_gauntlet(ge, {"Ganymede": 1.872, "Europa": 3.893})))
    controls.append(("EurGan#131", *from_gauntlet(ge, {"Ganymede": 4.103, "Europa": 2.402})))
    r316 = next(
        r
        for r in json.loads((REPO / "data/943_ganeur316_recall.json").read_text())
        if r.get("status") == "pass"
    )
    controls.append(("GanEur#316", "ge", r316["key"], r316["x_days"]))
    controls.append(("VenMar#45", *from_gauntlet("data/942_cell_vm2_gauntlet.json", {"V": 8.22})))
    controls.append(
        (
            "Hollister 1H",
            *from_gauntlet(ev, {"E": 2.994, "V": 3.19}, "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"),
        )
    )
    controls.append(("Hollister 2H", *from_gauntlet(ev, {"E": 5.599, "V": 6.022})))
    ro = next(
        r
        for r in json.loads((REPO / "data/942_em_ro_recall.json").read_text())
        if r["row"] == "2.5.1.+0" and r["status"] == "pass"
    )
    controls.append(("R-O 2.5.1.+0", "em", ro["key"], ro["x_days"]))
    for name, cell, key, xd in controls:
        out["controls"].append(run(name, sig_of(cell, key, xd)))
    # catalogued H&M rows: labels undeclared (no cycle key)
    cat = yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())
    rows = cat if isinstance(cat, list) else cat.get("cyclers", [])
    for r in rows:
        if str(r.get("id", "")).startswith("hollister-menning-1970-ev-orbit-"):
            enc = r.get("vinf_kms_at_encounters") or []
            sig = CandidateSignature(
                primary="Sun",
                sequence=tuple(str(e["body"]) for e in enc) or ("E", "V"),
                period_k=(r.get("period") or {}).get("k"),
                vinf_per_encounter_kms=tuple(
                    float(e["vinf_kms"]) for e in enc if e.get("vinf_kms") is not None
                ),
                topology_label=frozenset({"repeated-moon"}),
            )
            out["controls"].append(run(r["id"], sig))
    for name, path, pick in (
        ("gc-1", gc, {"Ganymede": 2.397}),
        ("gc-2", gc, {"Ganymede": 3.617}),
        ("ev-A", ev, {"E": 4.893, "V": 10.364}),
        ("ev-B", ev, {"E": 8.012, "V": 10.932}),
        ("ev-C", ev, {"E": 9.075, "V": 13.166}),
    ):
        out["candidates"].append(run(name, sig_of(*from_gauntlet(path, pick))))
    (REPO / "data" / "942_943_litcheck_scope.json").write_text(json.dumps(out, indent=1))
    bad = [c["name"] for c in out["controls"] if c["status"] != "published"]
    print("CONTROLS NOT PUBLISHED:", bad if bad else "none")


if __name__ == "__main__":
    main()
