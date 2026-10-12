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

Sec. 10 (amended 2026-10-12, sec. 10.1): the legacy path reproduces three #1048 gaps records and an
F1 seed its arm E record: status, iteration count, k_min, sheet and symmetry exactly; C and T of a
converged seed within a DERIVED rounding bound, not bit-for-bit. With C free the shooting system has
a family of solutions, so the stopping rule does not fix the landing point along the family; the
smallest singular value at the last Newton step (about 4e-8) amplifies the arc end points' rounding
noise into C. After the macOS 27.0.1 update (2026-10-12 09:38 AEDT) the same code lands 1.9e-7
(gaps) and 1.6e-9 (arm E) away from the records in C. Bound and positive control:
scripts/probe_1054_float_sensitivity.py.
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


# --- sec. 10: F1 is the default seed design; the legacy path reproduces the committed records ---

GAPS_SEEDS = [(3, 1, 1, 0, 0), (3, 1, 1, 0, 1), (3, 1, 1, 7, 0)]
ARM_E_SEEDS = [(2, 1, 1, 0, 24), (7, 3, -1, 3, 11)]

# Rounding-drift bounds (|C - C_rec|, |T - T_rec|) for the converged seeds, note sec. 10.1, from
# `scripts/probe_1054_float_sensitivity.py bound` (2026-10-12, rounded up). Derivation:
# |dx| <= sum over Newton steps of sqrt(4n) eps / sigma_min,k (minimum-norm step, error along the
# solution family never removed), with eps the arc end point's ROUNDING noise (1-ulp probe at the
# converged nodes: 2.55e-14 gaps, 6.21e-14 arm E); truncation error (rtol = atol = 1e-12) is common
# to both runs and cancels. |dC| <= |grad C(s0)| (|dx| + eps), |dT| <= sqrt(n) |dx| (k_min = p).
# Fixed constants, never recomputed from the solve under test. Positive control: the arm-E seed on
# the legacy row 4 converges with the same status and it at |dC| = 1.8e-4 and is rejected here; the
# truncation-scale bound (5.4e-4) would have accepted it.
DRIFT_TOL = {
    (3, 1, 1, 7, 0): (3.3e-6, 4.0e-6),
    (7, 3, -1, 3, 11): (1.4e-5, 1.3e-5),
}


def _assert_same_orbit(
    out: dict[str, Any], rec: dict[str, Any], seed: tuple[int, int, int, int, int]
) -> None:
    for key in ("k_min", "sheet", "symmetric"):
        assert out[key] == rec[key], key
    assert rec["k_min"] == seed[0]  # T = sum of the p durations, the bound's premise
    tol_c, tol_t = DRIFT_TOL[seed]
    assert abs(out["C"] - rec["C"]) <= tol_c
    assert abs(out["T"] - rec["T"]) <= tol_t


@pytest.mark.parametrize("seed", GAPS_SEEDS)
def test_legacy_gaps_seed_reproduces_1048_record(seed: tuple[int, int, int, int, int]) -> None:
    path = REPO / "data" / "1000_complement" / "shooting" / "gaps.jsonl"
    recs = {
        (r["p"], r["q"], r["sense"], r["ir"], r["iw"]): r
        for r in map(json.loads, path.read_text().splitlines())
    }
    rec = recs[seed]
    out = sh.solve_seed((*seed, None, "legacy"))
    assert (out["status"], out["it"]) == (rec["status"], rec["it"])
    if rec["status"] == "converged":
        _assert_same_orbit(out, rec, seed)
    assert "seed_design" not in out


def test_legacy_seven_tuple_equals_six_tuple() -> None:
    seed = REPRO_SEEDS[0]
    assert sh.solve_seed((*seed, None)) == sh.solve_seed((*seed, None, "legacy"))


def test_f1_rows_are_the_1054_rows() -> None:
    assert sh.RP_F1 is sd.RP_NEW
    assert sh.SEED_DESIGNS["f1"] is sh.RP_F1
    assert sh.SEED_DESIGNS["legacy"] is sh.RP


@pytest.mark.parametrize("seed", ARM_E_SEEDS)
def test_f1_seed_reproduces_arm_e_record(seed: tuple[int, int, int, int, int]) -> None:
    path = REPO / "data" / "1054_seed_design" / "arm_E.jsonl"
    recs = {
        (r["p"], r["q"], r["sense"], r["irn"], r["k"]): r
        for r in map(json.loads, path.read_text().splitlines())
    }
    rec = recs[seed]
    out = sh.solve_seed((*seed, None, "f1"))
    assert out["seed_design"] == "f1"
    assert (out["status"], out["it"]) == (rec["status"], rec["it"])
    if rec["status"] == "converged":
        _assert_same_orbit(out, rec, seed)
