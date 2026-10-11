"""#1052 rule U probe: the gate at a baseline commit against the working tree, over many signatures.

Signatures: every catalogue row with a sequence (primary inferred, unlabelled) and every 1-3-body
subset (with repeat forms) of every anchor's body set; each under the topologies in use plus none.
STOP RULE (docs/notes/2026-10-11-1052-rule-u-untopologied-signatures.md): no signature in any
topology may move toward "not-found" or "known-architecture-new-system"; topologied signatures must
not change at all; untopologied ones may move only published -> inconclusive.

Usage: uv run python scripts/probe_1052_rule_u.py --baseline <commit> --out <json>
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.search import literature_check as new

REPO = Path(__file__).resolve().parents[1]
TOPOLOGIES: list[frozenset[str]] = [
    frozenset(),
    *(
        frozenset({t})
        for t in (
            "repeated-moon",
            "resonant",
            "halo",
            "nrho",
            "tulip",
            "mga-tour",
            "pump-tour",
            "transfer",
            "quasi-satellite",
            "binary-coorbital",
        )
    ),
]
UNLABELLED_CITATIONS: set[str] = set()


def _baseline(commit: str) -> Any:
    src = subprocess.run(
        ["git", "show", f"{commit}:src/cyclerfinder/search/literature_check.py"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    path = Path(tempfile.mkdtemp()) / "lc_baseline.py"
    path.write_text(src)
    spec = importlib.util.spec_from_file_location("lc_baseline", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["lc_baseline"] = mod
    spec.loader.exec_module(mod)
    return mod


def _shapes() -> list[tuple[str, tuple[str, ...]]]:
    out: set[tuple[str, tuple[str, ...]]] = set()
    rows = yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())
    for r in rows:
        if r.get("sequence_canonical"):
            out.add((r.get("primary") or "Sun", tuple(str(r["sequence_canonical"]).split("-"))))
    for a in new.KNOWN_CORPUS:
        bodies = sorted(a.body_set)
        # 1-3 bodies (the pre-registration said 1-4; a 26-body anchor makes k = 4 about 36k
        # shapes, about 3 h; recorded as a deviation in the results)
        for k in range(1, min(3, len(bodies)) + 1):
            for comb in itertools.combinations(bodies, k):
                out.add((a.primary, comb))
                out.add((a.primary, (*comb, comb[0])))
                if k == 1:
                    out.add((a.primary, (comb[0], comb[0])))
    return sorted(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    old = _baseline(args.baseline)
    global UNLABELLED_CITATIONS
    UNLABELLED_CITATIONS = {x.citation for x in old.KNOWN_CORPUS if not x.topology_label}
    tally: Counter[tuple[str, str, str]] = Counter()
    bad: list[dict[str, Any]] = []
    n = 0
    shapes = _shapes()
    t0 = time.time()
    print(f"shapes: {len(shapes)} x {len(TOPOLOGIES)} topologies", flush=True)
    for i_shape, (primary, seq) in enumerate(shapes):
        if i_shape % 500 == 0:
            el = time.time() - t0
            eta = el / max(i_shape, 1) * (len(shapes) - i_shape)
            print(f"  {i_shape}/{len(shapes)} shapes, {n} checks, {len(bad)} violations, "
                  f"elapsed {el:.0f} s, ETA {eta:.0f} s", flush=True)  # fmt: skip
        for topo in TOPOLOGIES:
            kw = {
                "primary": primary,
                "sequence": seq,
                "period_k": None,
                "vinf_per_encounter_kms": (),
                "topology_label": topo,
            }
            a = old.check_literature(old.CandidateSignature(**kw), search=old.offline_corpus_search)
            b = new.check_literature(new.CandidateSignature(**kw), search=new.offline_corpus_search)
            n += 1
            if a.status == b.status:
                continue
            tname = ",".join(sorted(topo)) or "none"
            tally[(tname if topo else "none", a.status, b.status)] += 1
            moved_down = a.status == "published" and b.status == "inconclusive"
            single_sun_em = primary == "Sun" and len(set(seq)) == 1 and seq[0] in ("E", "M")
            via_unlabelled = a.citation in UNLABELLED_CITATIONS
            # note sec. A2: untopologied (rule U); single-body Sun E/M (the Aldrin fix); a
            # baseline citation from an unlabelled anchor (U2)
            allowed = moved_down and (not topo or single_sun_em or via_unlabelled)
            if not allowed or new.is_literature_fresh(b.status):
                bad.append({"primary": primary, "sequence": seq, "topology": sorted(topo),
                            "before": a.status, "after": b.status})  # fmt: skip
    print(f"signatures: {n}", flush=True)
    for k, v in sorted(tally.items()):
        print(" ", k, v, flush=True)
    print(f"STOP-RULE VIOLATIONS: {len(bad)}", flush=True)
    for x in bad[:20]:
        print("  ", x, flush=True)
    args.out.write_text(json.dumps({"n": n, "changes": {"|".join(k): v for k, v in tally.items()},
                                    "violations": bad}, indent=1))  # fmt: skip


if __name__ == "__main__":
    main()
