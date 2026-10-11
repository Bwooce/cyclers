"""#1052 (Pluto part): controls C1-C3, diagnostics and the 224 #998 strong signatures through the
v3 gate (offline corpus backend). Pre-registration:
docs/notes/2026-10-11-1052-pluto-anchor-labels-prereg.md (secs. 5, 7, 8, 10).

Signatures are the #998 type: primary Pluto, time-ordered encounters (the Charon flybys of the
block, then the passive moon), period_k, V_inf per encounter, topology {"repeated-moon"},
working_bodies "one", closing_repeat False. The 224 are the strong members of
data/998_pluto_smallmoons/<cell>_gauntlet.json.

Usage: uv run python scripts/litcheck_1052_pluto.py --out FILE
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)

REPO = Path(__file__).resolve().parents[1]
MOON = {"ps": "Styx", "pn": "Nix", "pk": "Kerberos", "ph": "Hydra"}


def run(sig: CandidateSignature) -> dict[str, Any]:
    r = check_literature(sig, search=offline_corpus_search)
    return {"status": r.status, "citation": (r.citation or "")[:90], "confidence": r.confidence}


def sig(
    primary: str = "Pluto",
    topology: str | None = "repeated-moon",
    sequence: tuple[str, ...] = ("Charon", "Nix"),
    working: str | None = "one",
) -> CandidateSignature:
    return CandidateSignature(
        primary=primary,
        sequence=sequence,
        period_k=4,
        vinf_per_encounter_kms=(0.178, 0.130),
        topology_label=frozenset({topology}) if topology else frozenset(),
        working_bodies=working,
        closing_repeat=False,
    )


def probes() -> dict[str, Any]:
    return {
        "C1 positive: #998 type {repeated-moon}, one": run(sig()),
        "C2 negative: same, primary Jupiter": run(sig(primary="Jupiter")),
        "C3 negative: untopologied Pluto (Charon, Nix)": run(sig(topology=None)),
        "diag: untopologied Pluto (Nix, Hydra)": run(sig(topology=None, sequence=("Nix", "Hydra"))),
        "diag: {pump-tour} at Pluto": run(sig(topology="pump-tour")),
        "diag: {mga-tour} (Charon, Nix), one": run(sig(topology="mga-tour")),
        "diag: {mga-tour} (Charon, Nix), two": run(sig(topology="mga-tour", working="two")),
        "diag: {resonant} (Charon, Nix)": run(sig(topology="resonant")),
        "diag: {resonant} (Nix, Hydra)": run(sig(topology="resonant", sequence=("Nix", "Hydra"))),
        "diag: {repeated-moon} (Charon, Nix), two": run(sig(working="two")),
        "diag: {repeated-moon} (Nix, Hydra)": run(sig(sequence=("Nix", "Hydra"))),
        "diag: {binary-coorbital} (Charon,)": run(
            sig(topology="binary-coorbital", sequence=("Charon",))
        ),
    }


def members() -> list[dict[str, Any]]:
    out = []
    for cell, moon in MOON.items():
        g = json.loads(
            (REPO / "data" / "998_pluto_smallmoons" / f"{cell}_gauntlet.json").read_text()
        )
        for c in g["candidates"]:
            if not c["strong"]:
                continue
            n = int(c["n_charon_flybys"])
            s = CandidateSignature(
                primary="Pluto",
                sequence=("Charon",) * n + (moon,),
                period_k=int(c["k"]),
                vinf_per_encounter_kms=(round(c["vinf_charon_kms"], 3),) * n
                + (round(c["vinf_moon_kms"], 3),),
                topology_label=frozenset({"repeated-moon"}),
                working_bodies="one",
                closing_repeat=False,
            )
            r = run(s)
            out.append({"cell": cell, "key": c["key"], **r})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    preflight_search(
        task_no=1052,
        region_id="literature-gate-v3-pluto-anchor-labels",
        method=MethodCapability(
            genome="literature_check v3, offline corpus, labelled signatures",
            corrector="none (lookup)",
            capability_tags=frozenset({"literature"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    pr = probes()
    for k, v in pr.items():
        print(f"{k}: {v['status']} {v['confidence']:.2f} {v['citation'][:50]}", flush=True)
    mem = members()
    counts = collections.Counter((m["status"], m["citation"][:40]) for m in mem)
    print(dict(counts), flush=True)
    args.out.write_text(json.dumps({"probes": pr, "members": mem}, indent=1))


if __name__ == "__main__":
    main()
