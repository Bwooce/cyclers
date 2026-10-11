"""#1054: the #1000 seed-design fix.

Pre-registration: docs/notes/2026-10-11-1054-seed-design-fix.md.

Two checks:

* reproduction (regression against the #1000 records, not a golden): 20 fixed old seeds, 10 from
  the 6,578 km perigee row and 10 whose skeleton passes the Moon within 5,000 km (the two bins
  with 93-94 % impacts), re-solved on the default path give the status and iteration count
  recorded in data/1000_complement/shooting/grid.jsonl, seed by seed;
* the 93-94 % mechanism (skeleton arcs that hit a body before any Newton step) on 40 old
  floor-row seeds: at least 0.80 (pre-registered). The F1+F2 design's rate on 40 matched seeds
  (0.300) is recorded in the note (sec. 5.1) with no bound: that design was falsified (lead ruling
  2026-10-11 (a)); a bound for a new design is set only by its own pre-registration.

F3 (the osculating-perigee step limit, opt-in, default off): the helper's osculating perigee of a
skeleton node equals the skeleton's r_p, and three arm-B seeds re-solved with F3 on give their
recorded status, iteration count and halving count (regression against
data/1054_seed_design/arm_B.jsonl, not a golden).
"""

from __future__ import annotations

import importlib
import json
import random
import sys
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
# imported by name at run time (driver scripts, not a package; mypy does not follow them)
sh = importlib.import_module("run_1000_shooting")
sd = importlib.import_module("run_1054_seed_design")

# (p, q, sense, ir, iw): 10 from the floor row (ir 0), 10 with skeleton lunar pass < 5,000 km
# (seeded choice, random.Random(1054); docs/notes/2026-10-11-1054-seed-design-fix.md sec. 3)
REPRO_SEEDS = [
    (1, 2, -1, 0, 20),
    (5, 3, 1, 0, 32),
    (2, 1, 1, 0, 3),
    (7, 3, 1, 0, 32),
    (5, 3, -1, 0, 29),
    (8, 3, -1, 0, 5),
    (5, 2, -1, 0, 1),
    (3, 2, 1, 0, 32),
    (1, 3, -1, 0, 5),
    (4, 3, 1, 0, 28),
    (8, 3, -1, 2, 1),
    (5, 3, -1, 3, 8),
    (8, 3, -1, 1, 6),
    (5, 2, 1, 4, 12),
    (4, 3, 1, 2, 32),
    (5, 2, -1, 5, 19),
    (2, 3, 1, 3, 23),
    (7, 3, 1, 7, 6),
    (7, 3, -1, 1, 22),
    (7, 3, -1, 1, 19),
]


def _records() -> dict[tuple[int, int, int, int, int], dict[str, Any]]:
    path = REPO / "data" / "1000_complement" / "shooting" / "grid.jsonl"
    recs = [json.loads(x) for x in path.read_text().splitlines()]
    return {(r["p"], r["q"], r["sense"], r["ir"], r["iw"]): r for r in recs}


@pytest.mark.parametrize("seed", REPRO_SEEDS)
def test_old_seed_reproduces_record(seed: tuple[int, int, int, int, int]) -> None:
    rec = _records()[seed]
    out = sh.solve_seed((*seed, None))
    assert out["status"] == rec["status"]
    assert out["it"] == rec["it"]


def test_osculating_perigee_of_skeleton_node_is_rp() -> None:
    rp = sd.RP_NEW[2]
    nodes, _ = sh.skeleton_nodes(7, 3, -1, rp, 0.7)
    for node in nodes:
        assert sh.osc_perigee(node) == pytest.approx(rp, rel=1e-12)


def test_old_floor_row_skeleton_arcs_impact() -> None:
    rng = random.Random(1054)
    old_floor = [s for s in sd.old_seeds() if s["ir"] == 0]
    old = rng.sample(old_floor, 40)
    old_rate = sum(sd.skeleton_arc_impact(s) for s in old) / len(old)
    print(f"skeleton-arc impact rate, old floor row: {old_rate:.3f}")
    assert old_rate >= 0.80


# (p, q, sense, ir, iw) from arm B (F3 on) with at least one halved step and a short solve
F3_SEEDS = [(2, 3, 1, 2, 6), (4, 3, 1, 1, 7), (3, 2, 1, 1, 10)]


@pytest.mark.parametrize("seed", F3_SEEDS)
def test_f3_opt_in_reproduces_arm_b(seed: tuple[int, int, int, int, int]) -> None:
    path = REPO / "data" / "1054_seed_design" / "arm_B.jsonl"
    recs = {
        (r["p"], r["q"], r["sense"], r["ir"], r["iw"]): r
        for r in map(json.loads, path.read_text().splitlines())
    }
    rec = recs[seed]
    p, q, sense, ir, iw = seed
    out = sh.solve_values({}, p, q, sense, sh.RP[ir], sh.OMEGA[iw], None, sd.DROP)
    assert (out["status"], out["it"], out["n_halved"]) == (
        rec["status"],
        rec["it"],
        rec["n_halved"],
    )
    assert out["n_halved"] > 0
