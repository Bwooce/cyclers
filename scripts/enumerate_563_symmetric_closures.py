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

Genericization (#575 C1)
-------------------------
``--primary``/``--moons``/``--out`` let this same construction be pointed at
ANY primary body + moon list (e.g. Saturn + Titan/Iapetus for #575) without
touching the construction logic above -- ``enumerate_direction``/
``pair_n_max`` simply forward ``primary`` to #558's already-genericized
``residual_at_point``/``gate_candidate`` (which have taken a ``primary=``
kwarg, defaulting to ``"Uranus"``, since #571). Running with NO arguments
reproduces the original Uranian, 4-non-Miranda-moon, 12-direction behavior
byte-for-byte -- this is the #575 C1 golden check: the genericization must
not alter a single Uranian result.
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


def pair_n_max(
    anchor: str, flyby: str, *, primary: str = "Uranus", tof_scale_max: float = TOF_SCALE_MAX
) -> tuple[float, float, float, int]:
    """(T_syn, P_a, P_b, n_max) for one pair -- n_max is direction-independent.

    ``primary``/``tof_scale_max`` (#575 genericization): default to the
    original Uranian values so every pre-#575 caller (which never passes
    these) is byte-for-byte unaffected.

    ``opposite_sense`` (#599) is derived from the pair's own registry
    ``retrograde`` flags (XOR: exactly one of the two orbits retrograde
    relative to the other), NOT a new parameter here -- every existing
    Uranian/Jovian/Saturnian moon has ``retrograde=False``, so this XOR is
    always False for every pre-#599 pair and ``synodic_period_days`` falls
    back to its original same-sense formula byte-for-byte. Only a pair like
    Neptune's Triton/Proteus (Triton retrograde, Proteus prograde) trips it.
    """
    mu = PRIMARIES[primary]
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    opposite_sense = sat_a.retrograde != sat_b.retrograde
    t_syn = synodic_period_days(mu, sat_a.sma_km, sat_b.sma_km, opposite_sense=opposite_sense)
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof_max = tof_scale_max * math.sqrt(p_a * p_b)
    n_max = math.floor(2.0 * tof_max / t_syn)
    return t_syn, p_a, p_b, n_max


def enumerate_direction(
    anchor: str, flyby: str, *, primary: str = "Uranus", tof_scale_max: float = TOF_SCALE_MAX
) -> dict[str, Any]:
    """Construct + gate every symmetric candidate for one ``anchor->flyby->anchor``
    direction. ``primary``/``tof_scale_max`` default to the original Uranian values
    (#575 genericization, construction logic below is UNCHANGED)."""
    t_syn, p_a, p_b, n_max = pair_n_max(anchor, flyby, primary=primary, tof_scale_max=tof_scale_max)
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
                    primary=primary,
                )
                if pt is None:
                    n_infeasible += 1
                    continue
                if pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                    continue
                n_subgate += 1
                gated = gate_candidate(anchor, flyby, pt, primary=primary)
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


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--primary",
        type=str,
        default="Uranus",
        help="Central body (PRIMARIES key). Default 'Uranus' reproduces the original "
        "#563 run byte-for-byte (the #575 C1 golden check).",
    )
    parser.add_argument(
        "--moons",
        type=str,
        default=",".join(NON_MIRANDA_MOONS),
        help="Comma-separated moon list; every ordered pair (both directions) among "
        "them is enumerated. Default: the 4 non-Miranda Uranian moons (12 directions).",
    )
    parser.add_argument(
        "--tof-scale-max",
        type=float,
        default=TOF_SCALE_MAX,
        help="Max tof_scale bound (must match the source discovery sweep's own tested "
        "max -- verify against that sweep's _meta record, do not assume). Default 3.0 "
        "(#558's own Uranian production sweep bound).",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="",
        help="Output JSONL path (default: data/enumerate_563_symmetric_closures.jsonl "
        "for the default Uranian args; otherwise required).",
    )
    args = parser.parse_args(argv)

    primary = args.primary
    moons = tuple(m.strip() for m in args.moons.split(",") if m.strip())
    tof_scale_max = args.tof_scale_max
    is_default_uranian = (
        primary == "Uranus" and moons == NON_MIRANDA_MOONS and tof_scale_max == TOF_SCALE_MAX
    )
    if args.out:
        out_path = Path(args.out)
    elif is_default_uranian:
        out_path = OUT_PATH
    else:
        raise SystemExit("--out is required when --primary/--moons/--tof-scale-max are non-default")

    t0 = time.time()
    directions = []
    for a, b in itertools.combinations(moons, 2):
        directions.append((a, b))
        directions.append((b, a))
    expected_directions = len(moons) * (len(moons) - 1)
    assert len(directions) == expected_directions, (
        f"expected {expected_directions} directions for {len(moons)} moons, got {len(directions)}"
    )

    print(
        f"[563] primary={primary} moons={moons} {len(directions)} directions x "
        f"{len(N_REV_VALUES) ** 2} n_rev combos, tof_scale_max={tof_scale_max}",
        flush=True,
    )

    all_results: list[dict[str, Any]] = []
    total_evaluated = 0
    total_passes = 0
    for anchor, flyby in directions:
        res = enumerate_direction(anchor, flyby, primary=primary, tof_scale_max=tof_scale_max)
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

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#563 direct symmetric-closure enumeration",
                    "primary": primary,
                    "moons": list(moons),
                    "directions": len(directions),
                    "n_rev_combos": len(N_REV_VALUES) ** 2,
                    "rel_offsets_deg": list(REL_OFFSETS_DEG),
                    "tof_scale_max_bound": tof_scale_max,
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
    print(f"[563] written: {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
