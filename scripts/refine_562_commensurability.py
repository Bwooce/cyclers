"""#562 -- per-basin continuous-tof synodic-commensurability refinement.

Fixes the exact grid-aliasing flaw a Fable adversarial second-opinion found in
#561's Stage-2 filter (see ``data/OUTSTANDING.md`` #558/#561/#562): thresholding
the *grid's* ``tof_scale`` (0.05-step) against a 0.01-0.02 commensurability
tolerance is finer than the grid can resolve, and the literal 0.02-cut +
integer>=2 restriction structurally EXCLUDED the Titania-Oberon/Oberon-Titania
directions (the two largest basin counts and the two with real adjacent-
published-literature risk) from ever being tested.

This script does NOT re-sweep the grid. It reads #558's already-computed,
already-gated 3672-record raw candidate set
(``data/scan_558_uranus_all_pairs_index.jsonl``), clusters it into distinct
basins per (anchor, flyby, n_rev) direction, and for EACH basin:

1. Computes the pair's synodic period ``T_syn`` (generic
   ``1/|1/T_a - 1/T_b|`` formula, matching the #561/#562 write-ups' own
   reference values, e.g. T_syn(Umbriel,Oberon) = 5.987 d).
2. Finds the nearest integer ``n >= 1`` such that the basin's full cycle
   (``2*tof``) is closest to ``n * T_syn`` -- explicitly ALLOWING n=1 (the
   Fable review's own note: "the re-closure condition 2*tof/T_syn in Z is
   satisfied by 1 as much as by 5").
3. Sets ``tof`` to the EXACT commensurate value (``n * T_syn / 2`` per leg) --
   a direct formula, not a search -- then does a small bounded 1-D local
   search over ``rel_offset_deg`` (holding this exact tof fixed) to find the
   true nearby residual minimum, since the basin's natural (rel_offset, tof)
   zero is generically an ISOLATED point in the 2-parameter plane and will
   not sit exactly on both the basin's original rel_offset AND the exact
   commensurate tof simultaneously (per the #561 Opus adjudication's own
   "isolated solutions" finding). This keeps the refinement grid-independent:
   the commensurate tof is fixed by physics (T_syn), not by the sweep's 0.05
   step.
4. Re-evaluates residual + the #324 physical bend gate + an independent
   DOP853 cross-check AT the refined point, reusing #558's own gate functions
   verbatim (``gate_candidate`` / ``candidate_passes_physical_gate`` /
   ``dop853_cross_check_leg``) -- no new gate logic.
5. A basin PASSES Stage 2 iff the refined, exactly-commensurate point closes
   sub-gate (residual < 0.05 km/s, #558's own threshold).

Discipline: NO catalogue writeback, NO V1-V4-strict gauntlet run here (that is
Stage 3, a separate follow-up). This is Sonnet-tier mechanical work behind
strong deterministic gates, per #558's own successful precedent.

Run as::

    uv run python scripts/refine_562_commensurability.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

# Reuse #558's own gate machinery verbatim -- residual_at_point (exact-point
# evaluation, not grid-sampled), gate_candidate (physical bend + DOP853), and
# GATE_RESIDUAL_KMS (the same 0.05 km/s sub-gate threshold).
from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    gate_candidate,
    residual_at_point,
)

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day  # noqa: E402

DATA_DIR = ROOT / "data"
INDEX_PATH = DATA_DIR / "scan_558_uranus_all_pairs_index.jsonl"
OUT_PATH = DATA_DIR / "scan_562_commensurability_refinement.jsonl"

# Basins within a direction/n_rev group are clustered when adjacent raw
# candidates are within this many degrees of relative offset -- the sweep's
# own grid resolution is 1 deg (#558 n_offset=360), so a 2 deg gap tolerates
# one missing grid sample without merging genuinely distinct basins (basins
# are typically separated by tens of degrees; see the manual inspection in
# the task's own build notes).
CLUSTER_GAP_DEG = 2.0

# Local rel_offset search window (deg) around the basin's representative
# point, at the FIXED exact-commensurate tof -- a coarse-then-refined 1-D
# scan, not a new grid sweep (see module docstring point 3).
SEARCH_HALF_WINDOW_DEG = 3.0
SEARCH_COARSE_N = 61  # 0.1 deg steps across the +/-3 deg window
SEARCH_REFINE_N = 41  # second pass, 0.005 deg steps across +/-0.1 deg


def synodic_period_days(
    mu: float, sma_a_km: float, sma_b_km: float, *, opposite_sense: bool = False
) -> float:
    """Generic two-body synodic period (days).

    Same orbital sense (``opposite_sense=False``, the default -- every
    Uranian/Jovian/Saturnian pair this was originally written for):
    ``1 / |1/T_a - 1/T_b|``. Matches the reference value quoted in the
    #561/#562 OUTSTANDING.md write-ups (T_syn(Umbriel,Oberon) = 5.987 d) --
    verified below in ``_self_check_synodic_reference``.

    Opposite orbital sense (``opposite_sense=True``, #599): one body prograde
    and the other retrograde (e.g. Neptune's Triton vs Proteus). The two
    bodies then move toward each other rather than one catching up to the
    other, so their relative angular rate is the SUM of the two mean
    motions, not the difference: ``1 / (1/T_a + 1/T_b)``. This is the
    textbook synodic-period relation for counter-orbiting bodies (derivable
    the same way as the same-sense case: relative angular rate magnitude is
    ``|n_a - (-n_b)| = n_a + n_b`` instead of ``|n_a - n_b|``).
    """
    n_a = _mean_motion_rad_day(mu, sma_a_km)
    n_b = _mean_motion_rad_day(mu, sma_b_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    if opposite_sense:
        return 1.0 / (1.0 / p_a + 1.0 / p_b)
    return 1.0 / abs(1.0 / p_a - 1.0 / p_b)


def _self_check_synodic_reference() -> None:
    mu = PRIMARIES["Uranus"]
    t_syn = synodic_period_days(mu, SATELLITES["Umbriel"].sma_km, SATELLITES["Oberon"].sma_km)
    ref = 5.987
    assert abs(t_syn - ref) < 0.01, f"T_syn(Umbriel,Oberon)={t_syn:.4f}, expected ~{ref}"
    print(f"[562] self-check: T_syn(Umbriel,Oberon) = {t_syn:.4f} d (reference ~5.987 d)")

    # Synthetic counter-orbiting check (#599): two arbitrary sma values about
    # an arbitrary primary, checked directly against the textbook relation
    # (not just eyeballed) -- same-sense uses the difference of mean
    # motions, opposite-sense uses the sum. Also sanity-checks the direction
    # of the effect: counter-orbiting bodies conjunct MORE often (shorter
    # T_syn) than co-orbiting bodies at the same two periods.
    mu_synth = PRIMARIES["Neptune"]
    sma_x_km, sma_y_km = 100_000.0, 250_000.0
    n_x = _mean_motion_rad_day(mu_synth, sma_x_km)
    n_y = _mean_motion_rad_day(mu_synth, sma_y_km)
    p_x, p_y = 2.0 * math.pi / n_x, 2.0 * math.pi / n_y
    expect_same = 1.0 / abs(1.0 / p_x - 1.0 / p_y)
    expect_opposite = 1.0 / (1.0 / p_x + 1.0 / p_y)
    got_same = synodic_period_days(mu_synth, sma_x_km, sma_y_km, opposite_sense=False)
    got_opposite = synodic_period_days(mu_synth, sma_x_km, sma_y_km, opposite_sense=True)
    assert abs(got_same - expect_same) < 1e-9, (got_same, expect_same)
    assert abs(got_opposite - expect_opposite) < 1e-9, (got_opposite, expect_opposite)
    assert got_opposite < got_same, "counter-orbiting pair must conjunct more often"
    print(
        f"[562] self-check: synthetic counter-orbiting T_syn same-sense={got_same:.4f}d "
        f"vs opposite-sense={got_opposite:.4f}d (opposite < same: {got_opposite < got_same})"
    )


@dataclass
class Basin:
    anchor: str
    flyby: str
    n_rev: tuple[int, int]
    members: list[dict[str, Any]]

    @property
    def representative(self) -> dict[str, Any]:
        return min(self.members, key=lambda c: c["record"]["residual_kms"])


def load_raw_candidates() -> list[dict[str, Any]]:
    with INDEX_PATH.open(encoding="utf-8") as fh:
        lines = fh.readlines()
    last = json.loads(lines[-1])
    assert last.get("kind") == "candidates_needing_adjudication"
    cands = last["candidates"]
    assert len(cands) == 3672, f"expected 3672 raw candidates, got {len(cands)}"
    return cands


def cluster_basins(cands: list[dict[str, Any]]) -> list[Basin]:
    """Cluster the raw candidate list into distinct basins.

    Grouped by (anchor, flyby, n_rev) -- a DIFFERENT revolution-count branch
    is a physically distinct trajectory family, never merged -- then, within
    each group, adjacent-in-rel_offset candidates (gap <= CLUSTER_GAP_DEG)
    are merged into one basin. This is a lightweight 1-D connectivity
    clustering, not a claim of rigorous basin topology; it exists only to
    avoid refining near-duplicate grid points that already describe the same
    underlying zero as separate "basins."
    """
    import collections

    groups: dict[tuple[str, str, tuple[int, int]], list[dict[str, Any]]] = collections.defaultdict(
        list
    )
    for c in cands:
        r = c["record"]
        key = (c["anchor"], c["flyby"], tuple(r["n_rev"]))
        groups[key].append(c)

    basins: list[Basin] = []
    for (anchor, flyby, n_rev), items in groups.items():
        items.sort(key=lambda c: c["record"]["rel_offset_deg"])
        cur: list[dict[str, Any]] = [items[0]]
        for prev, nxt in itertools.pairwise(items):
            gap = nxt["record"]["rel_offset_deg"] - prev["record"]["rel_offset_deg"]
            if gap <= CLUSTER_GAP_DEG:
                cur.append(nxt)
            else:
                basins.append(Basin(anchor, flyby, n_rev, cur))
                cur = [nxt]
        basins.append(Basin(anchor, flyby, n_rev, cur))
    return basins


def refine_basin(basin: Basin) -> dict[str, Any]:
    """Refine one basin's tof to the exact commensurate value; re-gate there.

    Returns a full report dict: the target integer n, T_syn, the refined
    point's record + gate result, and the ORIGINAL (unrefined) grid
    representative for comparison.
    """
    mu = PRIMARIES["Uranus"]
    sat_a = SATELLITES[basin.anchor]
    sat_b = SATELLITES[basin.flyby]
    t_syn = synodic_period_days(mu, sat_a.sma_km, sat_b.sma_km)
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b

    rep = basin.representative
    rep_rec = rep["record"]
    tof_days_grid = rep_rec["tof_days"]
    ratio = 2.0 * tof_days_grid / t_syn
    n_int = max(1, round(ratio))
    target_tof_days = n_int * t_syn / 2.0
    target_tof_scale = target_tof_days / math.sqrt(p_a * p_b)

    rel0 = rep_rec["rel_offset_deg"]
    n_rev = tuple(rep_rec["n_rev"])

    def scan(center: float, half_window: float, n_pts: int) -> dict[str, Any] | None:
        best: dict[str, Any] | None = None
        for i in range(n_pts):
            d = -half_window + 2.0 * half_window * i / (n_pts - 1)
            rel = (center + d) % 360.0
            pt = residual_at_point(
                basin.anchor,
                basin.flyby,
                rel_offset_deg=rel,
                tof_scale=target_tof_scale,
                n_rev=n_rev,
            )
            if pt is not None and (best is None or pt["residual_kms"] < best["residual_kms"]):
                best = pt
        return best

    coarse = scan(rel0, SEARCH_HALF_WINDOW_DEG, SEARCH_COARSE_N)
    refined = coarse
    if coarse is not None:
        fine = scan(coarse["rel_offset_deg"], 0.1, SEARCH_REFINE_N)
        if fine is not None and fine["residual_kms"] < coarse["residual_kms"]:
            refined = fine

    result: dict[str, Any] = {
        "anchor": basin.anchor,
        "flyby": basin.flyby,
        "n_rev": list(n_rev),
        "n_basin_members": len(basin.members),
        "t_syn_days": t_syn,
        "n_commensurate_int": n_int,
        "target_tof_days": target_tof_days,
        "grid_rep_tof_days": tof_days_grid,
        "grid_rep_rel_offset_deg": rel0,
        "grid_rep_residual_kms": rep_rec["residual_kms"],
        "grid_frac": abs(ratio - n_int),
    }

    if refined is None:
        result["refined_feasible"] = False
        result["stage2_pass"] = False
        return result

    result["refined_feasible"] = True
    result["refined_rel_offset_deg"] = refined["rel_offset_deg"]
    result["refined_residual_kms"] = refined["residual_kms"]

    if refined["residual_kms"] >= GATE_RESIDUAL_KMS:
        result["stage2_pass"] = False
        result["reason"] = "refined residual above 0.05 km/s sub-gate"
        return result

    gated = gate_candidate(basin.anchor, basin.flyby, refined)
    result["vinf_per_encounter_kms"] = gated["vinf_per_encounter_kms"]
    result["max_bend_deg_per_encounter"] = gated["max_bend_deg_per_encounter"]
    result["physical_gate_passed"] = gated["physical_gate_passed"]
    result["dop853_cross_check"] = gated["dop853_cross_check"]
    result["stage2_pass"] = bool(gated["all_gates_passed"])
    if not result["stage2_pass"]:
        result["reason"] = "refined point failed physical bend or DOP853 gate"
    return result


def main() -> int:
    t0 = time.time()
    _self_check_synodic_reference()

    cands = load_raw_candidates()
    basins = cluster_basins(cands)
    print(f"[562] {len(cands)} raw candidates -> {len(basins)} distinct basins", flush=True)

    by_dir: dict[tuple[str, str], int] = {}
    for b in basins:
        by_dir[(b.anchor, b.flyby)] = by_dir.get((b.anchor, b.flyby), 0) + 1
    for k, v in sorted(by_dir.items(), key=lambda kv: -kv[1]):
        print(f"[562]   {k[0]}-{k[1]}: {v} basins", flush=True)

    results: list[dict[str, Any]] = []
    n_pass = 0
    n_int1 = 0
    n_int1_pass = 0
    for i, basin in enumerate(basins):
        res = refine_basin(basin)
        results.append(res)
        if res["n_commensurate_int"] == 1:
            n_int1 += 1
            if res["stage2_pass"]:
                n_int1_pass += 1
        if res["stage2_pass"]:
            n_pass += 1
        if (i + 1) % 60 == 0 or (i + 1) == len(basins):
            print(
                f"[562] {i + 1}/{len(basins)} basins refined, "
                f"{n_pass} PASS so far  (elapsed {time.time() - t0:.1f}s)",
                flush=True,
            )

    elapsed = time.time() - t0
    print(
        f"[562] DONE: {len(basins)} basins refined -> {n_pass} pass Stage 2 "
        f"(residual + bend + DOP853 at the exact commensurate tof); "
        f"{n_int1} basins targeted n=1 ({n_int1_pass} pass)  ({elapsed:.1f}s)",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#562 per-basin continuous-tof commensurability refinement",
                    "n_raw_candidates": len(cands),
                    "n_basins": len(basins),
                    "n_stage2_pass": n_pass,
                    "n_targeting_int1": n_int1,
                    "n_targeting_int1_pass": n_int1_pass,
                    "elapsed_s": elapsed,
                    "cluster_gap_deg": CLUSTER_GAP_DEG,
                    "gate_residual_kms": GATE_RESIDUAL_KMS,
                }
            )
            + "\n"
        )
        for res in results:
            fh.write(json.dumps(res) + "\n")
    print(f"[562] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
