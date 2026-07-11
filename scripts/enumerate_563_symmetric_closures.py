"""#563 -- direct symmetric-closure enumeration for the #312 Uranian moon-pair family.

Replaces grid-search-plus-refinement (#558 -> #562) with a direct, exhaustive
CONSTRUCTION of the finite symmetric-closure family, closing the coverage gap
#562 could not close: #562 only refines basins that #558's ORIGINAL 1-deg /
0.05-tof_scale discrete grid happened to trip its initial residual gate on. A
basin narrower than that grid's own resolution can fall entirely between grid
points and never surface at all -- #558 would never find it, so #562 (which
only refines what #558 already surfaced) could never recover it either.

#562's own result data shows the dominant closure signature is a SYMMETRIC
periodic orbit: ``rel_offset`` in {0 deg, 180 deg}, matched leg revolution
counts, exactly-commensurate ``tof`` -- closing to machine precision (1e-14 to
1e-16 km/s). This is the classical Miele-style "perpendicular crossing"
symmetric-orbit condition, not a numerical coincidence (#562 makes this
argument explicitly for #312 itself, and reproduces #312 to 8.9e-16 km/s at
EXACTLY this condition). Because of this, the full symmetric-closure family is
a FINITE, ENUMERABLE set -- for each of the 12 physically-viable (non-Miranda;
see below) anchor-flyby directions, each ``n_rev=(n0,n1)`` in {0,1,2,3}^2 (16
combinations), each commensurability integer ``n>=1`` up to a ceiling bound by
the SAME max-tof range #558's own sweep covered, and ``rel_offset`` in
{0 deg, 180 deg} -- the candidate is directly CONSTRUCTED at the exact
commensurate ``tof = n*T_syn/2`` and evaluated. Nothing is searched for, so
there is no grid-resolution risk of missing a basin: every candidate in the
finite enumeration is checked directly.

Scope of "12 physically-viable directions": the 8 Miranda-involving
directions are NOT re-tested here -- #558's own data already shows these fail
the #324 physical bend gate at every grid point (confirmed real physics, not
a coverage gap; see the #563 OUTSTANDING.md entry). Only the
C(4,2)*2 = 12 directions among {Ariel, Umbriel, Titania, Oberon} are
enumerated.

n_max bound (per pair, direction-independent since it depends only on the
synodic period and the geometric mean of the two periods, both symmetric in
the pair): ``n_max = floor(2 * tof_max / T_syn)`` where
``tof_max = 3.0 * sqrt(P_a * P_b)`` -- 3.0 is the literal maximum
``tof_scale`` #558's production sweep used (verified against every
``scan_558_uranus_*.jsonl`` file's own ``_meta`` record: all 20 direction
files report ``tof_scales`` maxing at 3.0). This is a STRICT bound -- it keeps
every constructed candidate within the tof range #558 actually tested, so
this enumeration cannot claim to search anywhere #558 didn't. (One #562
survivor, Oberon-Titania n=3 tof=36.95d, sits just outside this strict bound
and at a non-symmetric rel_offset=114.15 deg -- it is a genuinely asymmetric
family member, correctly out of this task's scope; see the module's honest
scope-limit note in OUTSTANDING.md.)

Discipline: NO catalogue writeback, NO V1-V4-strict gauntlet run here. Reuses
#558/#562's own gate functions verbatim (``residual_at_point``,
``gate_candidate``, ``GATE_RESIDUAL_KMS``) -- no new gate logic.

Run as::

    uv run python scripts/enumerate_563_symmetric_closures.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

# Reuse #558's own gate machinery verbatim.
# Reuse #562's own synodic-period helper verbatim.
from refine_562_commensurability import synodic_period_days  # noqa: E402
from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    N_REV_MAX,
    gate_candidate,
    residual_at_point,
)

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day  # noqa: E402

DATA_DIR = ROOT / "data"
OUT_PATH = DATA_DIR / "enumerate_563_symmetric_closures.jsonl"

# The 4 non-Miranda regular moons -- the 8 Miranda-involving directions are
# confirmed-real #324 bend-gate failures (per #558's own data), not re-tested.
NON_MIRANDA_MOONS: tuple[str, ...] = ("Ariel", "Umbriel", "Titania", "Oberon")

# #558's actual production max tof_scale, verified against every
# scan_558_uranus_*.jsonl _meta record (all 20 direction files agree).
TOF_SCALE_MAX = 3.0

REL_OFFSETS_DEG: tuple[float, ...] = (0.0, 180.0)
N_REV_VALUES: tuple[int, ...] = tuple(range(N_REV_MAX + 1))  # 0..3


def pair_n_max(anchor: str, flyby: str) -> tuple[float, float, float, int]:
    """(T_syn, P_a, P_b, n_max) for one pair -- n_max is direction-independent."""
    mu = PRIMARIES["Uranus"]
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    t_syn = synodic_period_days(mu, sat_a.sma_km, sat_b.sma_km)
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof_max = TOF_SCALE_MAX * math.sqrt(p_a * p_b)
    n_max = math.floor(2.0 * tof_max / t_syn)
    return t_syn, p_a, p_b, n_max


def enumerate_direction(anchor: str, flyby: str) -> dict[str, Any]:
    t_syn, p_a, p_b, n_max = pair_n_max(anchor, flyby)
    sqrt_papb = math.sqrt(p_a * p_b)

    n_evaluated = 0
    n_infeasible = 0
    n_subgate = 0
    passes: list[dict[str, Any]] = []

    for n in range(1, n_max + 1):
        target_tof_days = n * t_syn / 2.0
        target_tof_scale = target_tof_days / sqrt_papb
        for n0, n1 in itertools.product(N_REV_VALUES, N_REV_VALUES):
            for rel in REL_OFFSETS_DEG:
                n_evaluated += 1
                pt = residual_at_point(
                    anchor,
                    flyby,
                    rel_offset_deg=rel,
                    tof_scale=target_tof_scale,
                    n_rev=(n0, n1),
                )
                if pt is None:
                    n_infeasible += 1
                    continue
                if pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                    continue
                n_subgate += 1
                gated = gate_candidate(anchor, flyby, pt)
                if gated["all_gates_passed"]:
                    passes.append(
                        {
                            "anchor": anchor,
                            "flyby": flyby,
                            "n_commensurate_int": n,
                            "t_syn_days": t_syn,
                            "rel_offset_deg": rel,
                            "n_rev": [n0, n1],
                            "tof_days": target_tof_days,
                            "residual_kms": pt["residual_kms"],
                            "vinf_per_encounter_kms": gated["vinf_per_encounter_kms"],
                            "max_bend_deg_per_encounter": gated["max_bend_deg_per_encounter"],
                            "dop853_cross_check": gated["dop853_cross_check"],
                        }
                    )

    return {
        "anchor": anchor,
        "flyby": flyby,
        "t_syn_days": t_syn,
        "n_max": n_max,
        "n_evaluated": n_evaluated,
        "n_infeasible": n_infeasible,
        "n_subgate_residual_only": n_subgate,
        "n_all_gates_passed": len(passes),
        "passes": passes,
    }


def main() -> int:
    t0 = time.time()
    directions = []
    for a, b in itertools.combinations(NON_MIRANDA_MOONS, 2):
        directions.append((a, b))
        directions.append((b, a))
    assert len(directions) == 12, f"expected 12 non-Miranda directions, got {len(directions)}"

    print(f"[563] {len(directions)} directions x {len(N_REV_VALUES) ** 2} n_rev combos", flush=True)

    all_results: list[dict[str, Any]] = []
    total_evaluated = 0
    total_passes = 0
    for anchor, flyby in directions:
        res = enumerate_direction(anchor, flyby)
        all_results.append(res)
        total_evaluated += res["n_evaluated"]
        total_passes += res["n_all_gates_passed"]
        print(
            f"[563] {anchor}-{flyby}-{anchor}: T_syn={res['t_syn_days']:.4f}d n_max={res['n_max']} "
            f"evaluated={res['n_evaluated']} sub_gate={res['n_subgate_residual_only']} "
            f"all_gates_pass={res['n_all_gates_passed']}",
            flush=True,
        )

    elapsed = time.time() - t0
    print(
        f"[563] DONE: {total_evaluated} candidates directly evaluated across "
        f"{len(directions)} directions, {total_passes} pass ALL gates (residual+bend+DOP853) "
        f"({elapsed:.1f}s)",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#563 direct symmetric-closure enumeration",
                    "directions": len(directions),
                    "n_rev_combos": len(N_REV_VALUES) ** 2,
                    "rel_offsets_deg": list(REL_OFFSETS_DEG),
                    "tof_scale_max_bound": TOF_SCALE_MAX,
                    "total_evaluated": total_evaluated,
                    "total_all_gates_passed": total_passes,
                    "elapsed_s": elapsed,
                    "gate_residual_kms": GATE_RESIDUAL_KMS,
                }
            )
            + "\n"
        )
        for res in all_results:
            fh.write(
                json.dumps(
                    {
                        "kind": "direction_summary",
                        "anchor": res["anchor"],
                        "flyby": res["flyby"],
                        "t_syn_days": res["t_syn_days"],
                        "n_max": res["n_max"],
                        "n_evaluated": res["n_evaluated"],
                        "n_infeasible": res["n_infeasible"],
                        "n_subgate_residual_only": res["n_subgate_residual_only"],
                        "n_all_gates_passed": res["n_all_gates_passed"],
                    }
                )
                + "\n"
            )
            for p in res["passes"]:
                fh.write(json.dumps({"kind": "pass", **p}) + "\n")
    print(f"[563] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
