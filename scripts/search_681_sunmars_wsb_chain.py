#!/usr/bin/env python3
"""Task #681 -- Sun-Mars WSB repeating-capture (quasi-cycler) chain search.

The interplanetary analog of task #378's cislunar BCT-chain search
(`cislunar-bct-wsb-quasicycler-2026-06-26`, clean negative). The object class:
a *repeating* capture<->escape chain whose return leg re-acquires the Mars
weak-stability-boundary set W each cycle -- a `quasi_cycler`-class object, NOT a
one-shot ballistic-capture transfer (which is a precursor_mga, and is all the
positive-control paper Topputo & Belbruno 2015 itself constructs).

Method. For a grid of Mars WSB capture seeds -- periapsis states at radius r_p
on osculating ellipses (eccentricity e, orientation theta), for a set of Mars
initial true anomalies f0 and both Mars-relative senses -- integrate each seed
50 Mars revolutions (~94 yr, matching the paper's own Sect. 4.2 long-integration
test) in BOTH time directions in the planar Sun-Mars elliptic restricted problem
(`cyclerfinder.core.sunmars_wsb`). A seed is a repeating-chain CANDIDATE iff,
after escaping the Mars vicinity (a periapsis with E_2 > 0), the trajectory
undergoes >= 2 distinct SUSTAINED recapture episodes -- maximal runs of >= 2
consecutive bound Mars periapses (E_2 <= 0) within the Hill sphere, i.e. genuine
temporary captures that RECUR (`StabilityResult.n_recapture_episodes >= 2`). A
single such episode is an isolated temporary capture, not a cycler; a lone
bound-energy dip during a co-orbital conjunction (`recaptured_after_escape`,
kept only as a loose diagnostic) is not a capture at all. The expected default,
per Belbruno 2004 Thm 3.58 (capture on W is chaotic), task #378's cislunar clean
negative, and the paper's own Sect. 4.2 (no second capture in 50 revs), is ZERO
candidates.

PROVISIONAL admissibility choices (flagged, not user-decided -- see the #681
bullet's open-decision point). These are generous (bias toward FINDING a chain):
  * minimum captures per cycle: >= 1 re-acquisition after an escape (any bound
    Mars periapsis post-escape counts);
  * periodicity tolerance: none imposed at this stage -- ANY re-acquisition
    within the 50-rev horizon is flagged (a true quasi-cycler needs a repeating
    *period*, a strictly stronger condition; a zero count here forecloses it);
  * dv_band ceiling for per-cycle correction: not applied (a purely ballistic
    re-acquisition is sought first; deterministic small-correction chains are a
    strictly larger, later question).
A nonzero count would NOT be a catalogue admission -- it would be handed back for
scrutiny (dual Opus+Fable adjudication) per this project's standing practice.

Checkpoint/resume/assemble pattern (cf. scripts/certify_678_*). Results stream
to a JSONL checkpoint; --resume skips completed cells; --assemble summarizes.

Run:
  uv run python scripts/search_681_sunmars_wsb_chain.py            # full run
  uv run python scripts/search_681_sunmars_wsb_chain.py --assemble # summary
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.sunmars_wsb as sw

_CKPT = Path("docs/notes/scratch/681_sunmars_chain_ckpt.jsonl")

# Search grid (provisional; generous coverage of the capture family).
_RP_KM = (49896.0, 73896.0, 91897.0, 113897.0)  # Table 3 r_p span
_E_VALUES = (0.90, 0.99)  # paper's two eccentricities
_F0_VALUES = (0.0, math.pi / 4.0, math.pi / 2.0)  # paper's f0 set
_BRANCHES = ("prograde", "retrograde")
_N_THETA = 24  # periapsis orientation family selector
_DIRECTIONS = ("forward", "backward")
_HORIZON_REVS = 50.0  # Mars revolutions (~94 yr), the Sect. 4.2 span


def _cells() -> list[dict[str, object]]:
    cells: list[dict[str, object]] = []
    thetas = np.linspace(0.0, 2.0 * math.pi, _N_THETA, endpoint=False)
    for rp in _RP_KM:
        for ecc in _E_VALUES:
            for f0 in _F0_VALUES:
                for branch in _BRANCHES:
                    for theta in thetas:
                        for direction in _DIRECTIONS:
                            cells.append(
                                {
                                    "rp": rp,
                                    "e": ecc,
                                    "f0": f0,
                                    "branch": branch,
                                    "theta": float(theta),
                                    "direction": direction,
                                }
                            )
    return cells


def _key(cell: dict[str, object]) -> str:
    return (
        f"{cell['rp']}|{cell['e']}|{cell['f0']:.5f}|{cell['branch']}|"
        f"{cell['theta']:.6f}|{cell['direction']}"
    )


def _done_keys() -> set[str]:
    if not _CKPT.exists():
        return set()
    keys: set[str] = set()
    for line in _CKPT.read_text().splitlines():
        if line.strip():
            keys.add(json.loads(line)["key"])
    return keys


def _evaluate(cell: dict[str, object]) -> dict[str, object]:
    st = sw.capture_periapsis_state(
        r_p_km=float(cell["rp"]),
        ecc=float(cell["e"]),
        theta=float(cell["theta"]),
        f0=float(cell["f0"]),
        branch=cell["branch"],  # type: ignore[arg-type]
    )
    res = sw.integrate_stability(
        st,
        float(cell["f0"]),
        direction=cell["direction"],  # type: ignore[arg-type]
        horizon_revs=_HORIZON_REVS,
    )
    return {
        "key": _key(cell),
        **{k: cell[k] for k in ("rp", "e", "f0", "branch", "theta", "direction")},
        "escaped": res.escaped,
        "recaptured_after_escape": res.recaptured_after_escape,
        "n_recapture_episodes": res.n_recapture_episodes,
        "max_sustained_bound_revs": res.max_sustained_bound_revs,
        "min_recapture_dist_km": res.min_recapture_dist_km,
        "collided": res.collided,
        "n_captured_revs": res.n_captured_revs,
        "n_periapses": len(res.periapses),
    }


def run() -> None:
    _CKPT.parent.mkdir(parents=True, exist_ok=True)
    cells = _cells()
    done = _done_keys()
    todo = [c for c in cells if _key(c) not in done]
    print(f"total cells={len(cells)}  done={len(done)}  todo={len(todo)}")
    t0 = time.time()
    with _CKPT.open("a") as fh:
        for i, cell in enumerate(todo, 1):
            rec = _evaluate(cell)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if i % 25 == 0 or i == len(todo):
                el = time.time() - t0
                eta = el / i * (len(todo) - i)
                print(
                    f"  [{i}/{len(todo)}] {el:6.1f}s elapsed  ETA {eta / 60:5.1f}min  "
                    f"repeating_chains_so_far={_count_repeaters()}"
                )
    assemble()


def _count_repeaters() -> int:
    if not _CKPT.exists():
        return 0
    n = 0
    for line in _CKPT.read_text().splitlines():
        if line.strip() and json.loads(line).get("n_recapture_episodes", 0) >= 2:
            n += 1
    return n


def assemble() -> None:
    if not _CKPT.exists():
        print("no checkpoint")
        return
    recs = [json.loads(ln) for ln in _CKPT.read_text().splitlines() if ln.strip()]
    n = len(recs)
    n_esc = sum(1 for r in recs if r["escaped"])
    n_loose = sum(1 for r in recs if r["recaptured_after_escape"])
    n_sustained = sum(1 for r in recs if r.get("n_recapture_episodes", 0) >= 1)
    n_repeat = sum(1 for r in recs if r.get("n_recapture_episodes", 0) >= 2)
    n_coll = sum(1 for r in recs if r["collided"])
    print("\n=== #681 Sun-Mars WSB repeating-capture chain search -- summary ===")
    print(f"cells evaluated                        : {n}")
    print(f"escaped Mars vicinity                  : {n_esc}")
    print(f"collided with Mars                     : {n_coll}")
    print(f"loose post-escape energy dip (E2<0)    : {n_loose}  (co-orbital, NOT a capture)")
    print(f">=1 SUSTAINED recapture episode        : {n_sustained}  (isolated temporary captures)")
    print(f">=2 recapture episodes (REPEATING)     : {n_repeat}   <-- chain candidates")
    if n_repeat:
        print("  CANDIDATES (hand back for Opus+Fable adjudication, do NOT self-admit):")
        for r in recs:
            if r.get("n_recapture_episodes", 0) >= 2:
                print(
                    f"    rp={r['rp']} e={r['e']} f0={r['f0']:.3f} {r['branch']} "
                    f"theta={r['theta']:.3f} {r['direction']} "
                    f"episodes={r['n_recapture_episodes']}"
                )
    else:
        print("  VERDICT: EMPTY -- no REPEATING capture chain. The sustained")
        print("  recaptures are all single, isolated, chaotic re-encounters")
        print("  (episodes=1), reproducing the paper's Sect. 4.2 non-recurrence")
        print("  + Belbruno Thm 3.58 chaos + #378's cislunar clean negative.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemble", action="store_true", help="summarize checkpoint only")
    args = ap.parse_args()
    if args.assemble:
        assemble()
    else:
        run()


if __name__ == "__main__":
    main()
