"""#583 (stage 3b of #581): widened-domain Sun-Earth ER3BP bounded-drift search.

Widens #581 stage 2's niching-GA reproduction beyond Gurfil-Kasdin (2002)'s
own 12 published optimization-set boxes, reusing BOTH
``search/niching_ga.py::run_deterministic_crowding`` and
``core/er3bp_geocentric.py::gurfil_kasdin_fitness`` UNMODIFIED (per #583's
own scope) -- the only new work is the widened bounds/genome + the long-span
bounded-vs-divergent drift classifier that gates any survivor before it is
trusted (``data/validation/er3bp_drift_classifier.py``).

Explicit, written-down decisions (per #583's own mandate -- do not leave to
the reader to infer):

1. **Free vs fixed state dims.** Genome is 7-dimensional: the 6 interleaved
   geocentric state slots ``[x, x', y, y', z, z']`` (same convention as
   stage 2's ``TABLE2``/``table_interleaved_to_state``) PLUS ``theta0`` as a
   free 7th gene. Each of the 7 :data:`PARTITIONS` below fixes a SUBSET of
   the 6 state slots at 0 (mirroring which published families it targets)
   and frees the rest -- same "fix most components, vary a few" structure
   Gurfil-Kasdin's own Table 2 uses, generalized across their per-set
   groupings rather than redone per set.
2. **theta0 IS free** (recommended by the task spec, done here): every
   partition searches ``theta0 in [0, 2*pi)`` as a genome variable, instead
   of stage 2's per-set FIXED theta0 in {0, pi}. Leaving it fixed would not
   genuinely explore beyond the 12 published sets (theta0 is the true-
   anomaly launch phase; stage 2 only ever tested Earth-perihelion and
   Earth-aphelion launches).
3. **n_rev stays 1**, matching stage 2 and ``gurfil_kasdin_fitness``'s own
   default -- not widened here (no motivating reason given in the spec to
   widen it, and doing so would conflate two independent axes of widening
   in one run).
4. **escape_radius stays 0.5 AU** (the ``gurfil_kasdin_fitness`` default,
   ~50x Earth's Hill radius, ~1.5e6 km) -- already implicitly admits the
   full heliocentric co-orbital regime (quasi-satellites, horseshoes,
   Sun-Earth DROs) that the #583 corpus-anchor prerequisite exists to
   police. This is a stated decision, not a surprise: a bounded-drift
   survivor near that boundary is exactly the territory the 3 new
   ``literature_check.py`` anchors (Gurfil-Kasdin 2002, Sun-Earth
   co-orbital, Henon family-f) were filed to catch.
5. **Widened bounds** (the actual numbers, justified below) cover the space
   BETWEEN and AROUND the union of the 12 published sets' boxes -- NOT a
   redo of stage 2's own per-set boxes, and NOT the paper's own
   sqrt(2)/sqrt(3) norm-budgeting convention (that convention existed only
   to keep a MULTI-COMPONENT vector's overall magnitude inside
   [``LO_R``, ``HI_R``] when 2-3 components shared one budget; the widened
   search instead gives every free position/velocity component its OWN
   independent bound, which is itself part of the widening -- a published
   family can no longer "hide" behind the norm-sharing constraint).

   * ``LO_WIDE = LO_R`` (1e6 km, the paper's own floor). CORRECTED
     2026-07-13 after a Fable review of a failed positive control: an
     earlier choice of ``0.002`` AU (~299,200 km, closer-in than ``LO_R``)
     re-admitted a physically trivial, strictly fitness-dominant deep-Hill-
     sphere quasi-circular basin that Eq. 15 cannot distinguish from a
     genuine family (it rewards only annulus thinness, with no periodicity
     content) -- deterministic crowding collapsed the entire population onto
     it. See the constant's own inline comment and docs/notes/2026-07-13-
     583-corpus-anchors-and-drift-classifier.md for the full diagnosis.
   * ``HI_WIDE = 0.15`` AU (~2.244e7 km): ~2.24x FARTHER than the paper's
     own ``HI_R`` (1e7 km) -- stays at 30% of ``escape_radius`` (0.5 AU),
     leaving room for a genuinely bounded family to breathe without
     immediately tripping the escape death-penalty.
   * ``V_WIDE = 0.5`` (AU/rad, ~15 km/s): numerically the SAME magnitude as
     the paper's own ``V``, but every partition below grants it to EACH free
     velocity component INDEPENDENTLY and over the FULL SIGNED range
     ``[-V_WIDE, V_WIDE]`` -- stage 2's sets each searched only one half
     (``[-V, 0]`` or ``[0, V]``) per set and shared the budget via
     sqrt(2)/sqrt(3) when >1 component was free. This is the actual
     "between and around" widening for the velocity axis.

Positive control (this dispatch's actual deliverable, per #583's
Fable-corrected gate -- NOT the full novel-territory sweep, which is a
longer job for a future dispatch): run partition ``P1`` (the richest and
cheapest -- 2 free state dims + theta0 -- covering 6 of the 14 published
families A-F) to completion and confirm it still recovers ITS neighboring
known families under stage 2's own pre-registered non-bit-exact match
criterion (``IC_PROXIMITY_TOL``/``FEATURE_FACTOR_TOL``, imported unmodified
from ``run_581_gurfil_reproduction.py``). See
``docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md`` for the
result.

Usage (chunked, same convention as #581's stage-2 script):
    uv run python scripts/run_583_widened_bounded_drift_search.py \\
        --partition P1 [--max-gens 100] --workers 8
    uv run python scripts/run_583_widened_bounded_drift_search.py --analyze --partition P1

Checkpoints/results: data/found/583_widened_search/
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from cyclerfinder.core.er3bp_geocentric import (
    SUN_EARTH_ER3BP,
    gurfil_kasdin_fitness,
    table_interleaved_to_state,
)
from cyclerfinder.data.preflight import MethodCapability, preflight_search
from cyclerfinder.data.validation.er3bp_drift_classifier import (
    N_REVS_DEFAULT,
    classify_bounded_drift,
    spot_check_theta0_robustness,
)
from cyclerfinder.search.niching_ga import (
    DeterministicCrowdingConfig,
    run_deterministic_crowding,
)

# Reuse stage 2's own match criterion + dynamics characterizer UNMODIFIED --
# the widened search must be judged by the SAME non-bit-exact bar, not a new
# one invented for this task. Same repo-root sys.path convention as other
# scripts/-to-scripts/ reuse in this project (see pyproject.toml's mypy
# overrides comment, e.g. run_392_v4_annual_sweep.py).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.run_581_gurfil_reproduction import (
    FEATURE_FACTOR_TOL,
    IC_PROXIMITY_TOL,
    LO_R,
    TABLE34,
    characterize,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "found" / "583_widened_search"

_REGION_ID = "sun-earth-er3bp-widened-bounded-drift-583-positive-control-2026-07-13"
_METHOD = MethodCapability(
    genome=(
        "Gurfil & Kasdin (2002) geocentric ER3BP objective (Eq. 15, "
        "core/er3bp_geocentric.py) evaluated by deterministic-crowding niching GA "
        "(search/niching_ga.py) over a genuinely WIDENED 7-dim genome (6 "
        "interleaved state slots + free theta0), bounds covering the space "
        "between/around the union of Gurfil-Kasdin's own 12 published sets"
    ),
    corrector=(
        "no corrector -- this dispatch runs ONLY the P1 positive-control "
        "partition (confirming the widened machinery still recovers stage "
        "2's own known families); the full novel-territory sweep across all "
        "7 partitions is explicitly deferred to a future dispatch"
    ),
    capability_tags=frozenset(
        {
            "er3bp",
            "geocentric",
            "niching-ga",
            "widened-domain",
            "positive-control",
            "sun-earth",
        }
    ),
    git_sha="working-tree",
)

# ---------------------------------------------------------------------------
# Widened bound constants (see module docstring §5 for the derivation).
# ---------------------------------------------------------------------------
LO_WIDE = LO_R  # AU (~1e6 km): CORRECTED 2026-07-13 (Fable review) -- do NOT
# go closer-in than the paper's own LO_R. A 0.002 AU floor (~299,200 km,
# ~0.2 Earth Hill radii) re-admits a physically trivial, strictly-dominant
# quasi-circular deep-Hill-sphere basin: Eq. 15 rewards only 1yr annulus
# thinness with no periodicity/family content, and deep inside the Hill
# sphere essentially ANY near-circular orbit scores ~1 (deficit ~5.6e-9 for
# the trivial basin vs ~8.7e-7 for the best genuine target family, a 150x
# fitness gap) -- deterministic crowding cannot protect a lower-fitness niche
# against a larger, smoothly-reachable, strictly-higher-fitness basin under
# child>=parent replacement, so the whole population collapsed onto the
# trivial solution (empirically confirmed: 75% of the gen-400 population
# within 0.7% of the 0.002 AU floor, zero members within 20% of any
# published family). See docs/notes/2026-07-13-583-corpus-anchors-and-drift-
# classifier.md for the full diagnosis. A future full sweep MAY explore the
# deep-Hill slice deliberately as ITS OWN partition, judged by clustering +
# drift classification + the corpus anchors filed below (NOT by Eq. 15 rank,
# which is uninformative there) -- not attempted in this dispatch.
HI_WIDE = 0.15  # AU (~2.244e7 km): ~2.24x farther than the paper's HI_R; 30% of escape_radius.
V_WIDE = 0.5  # AU/rad (~15 km/s): same magnitude as paper's V, full signed range, per-component.
THETA0_LO, THETA0_HI = 0.0, 2.0 * math.pi

_POS = (LO_WIDE, HI_WIDE)
_VEL = (-V_WIDE, V_WIDE)
_FIX = (0.0, 0.0)


def _bounds7(state6_bounds: list[tuple[float, float]]) -> list[tuple[float, float]]:
    return [*state6_bounds, (THETA0_LO, THETA0_HI)]


# name -> (bounds7 [x,x',y,y',z,z',theta0], families this partition targets,
# stage-2 free-dim-index signature it generalizes).
PARTITIONS: dict[str, tuple[list[tuple[float, float]], tuple[str, ...], tuple[int, ...]]] = {
    "P1": (
        _bounds7([_POS, _FIX, _FIX, _VEL, _FIX, _FIX]),
        ("A", "B", "C", "D", "E", "F"),
        (0, 3),
    ),
    "P2": (
        _bounds7([_POS, _FIX, _POS, _FIX, _FIX, _FIX]),
        ("G",),
        (0, 2),
    ),
    "P3": (
        _bounds7([_POS, _VEL, _FIX, _FIX, _FIX, _FIX]),
        ("H",),
        (0, 1),
    ),
    "P4": (
        _bounds7([_POS, _VEL, _POS, _VEL, _FIX, _FIX]),
        ("I",),
        (0, 1, 2, 3),
    ),
    "P5": (
        _bounds7([_POS, _VEL, _POS, _VEL, _POS, _VEL]),
        ("J", "K", "L"),
        (0, 1, 2, 3, 4, 5),
    ),
    "P6": (
        _bounds7([_POS, _FIX, _POS, _FIX, _POS, _FIX]),
        ("M",),
        (0, 2, 4),
    ),
    "P7": (
        _bounds7([_POS, _FIX, _FIX, _FIX, _POS, _FIX]),
        ("N",),
        (0, 4),
    ),
}


def fitness_widened(vec7: np.ndarray) -> float:
    """GA fitness: decode the 7-gene [state6_interleaved, theta0] genome."""
    vec7 = np.asarray(vec7, dtype=float)
    state = table_interleaved_to_state(vec7[:6])
    theta0 = float(vec7[6])
    return gurfil_kasdin_fitness(
        state, theta0, SUN_EARTH_ER3BP, n_rev=1.0, rtol=1e-9, atol=1e-9, n_samples=2000
    )


def run_partition(name: str, max_gens: int | None, workers: int) -> None:
    bounds, families, _sig = PARTITIONS[name]
    config = DeterministicCrowdingConfig(seed=583000 + sorted(PARTITIONS).index(name))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ckpt = OUT_DIR / f"{name}_checkpoint.npz"
    runlog = OUT_DIR / f"{name}_runlog.jsonl"
    t_start = time.monotonic()

    def progress(stats: dict[str, float]) -> None:
        rec = {
            "ts": datetime.now(UTC).isoformat(timespec="seconds"),
            "partition": name,
            **stats,
            "elapsed_s": round(time.monotonic() - t_start, 1),
        }
        with runlog.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        gen = int(stats["generation"])
        if gen % 20 == 0 or gen == config.generations:
            print(
                f"[{name}] gen {gen}/{config.generations} "
                f"mean={stats['fitness_mean']:.6f} max={stats['fitness_max']:.6f} "
                f"std={stats['fitness_std']:.6f} elapsed={rec['elapsed_s']}s",
                flush=True,
            )

    result = run_deterministic_crowding(
        fitness_widened,
        bounds,
        config,
        workers=workers,
        checkpoint_path=ckpt,
        max_generations_this_call=max_gens,
        progress_fn=progress,
    )
    print(
        f"[{name}] done through generation {result.generations_run} (target families: {families})",
        flush=True,
    )
    if result.generations_run >= config.generations:
        np.savez(
            OUT_DIR / f"{name}_final.npz",
            phenotypes=result.phenotypes,
            fitness=result.fitness,
        )
        print(f"[{name}] final population saved", flush=True)


# ---------------------------------------------------------------------------
# Analysis: stage 2's own match criterion, generalized to the 7-gene genome.
# ---------------------------------------------------------------------------


def match_family_in_widened_population(
    fam: str,
    partition: str,
    phen: np.ndarray,
    fitness: np.ndarray,
) -> dict[str, object]:
    """Family-directed match against the 7-gene widened population.

    Same three-part criterion as stage 2's own
    ``run_581_gurfil_reproduction.py::match_family_in_population``
    (IC proximity + type match + rmin/rmax feature ratio), generalized to
    include theta0 as a compared gene (stage 2 fixed theta0 per set so it
    never entered the distance; here it is free, so a genuine match must
    also land near the published theta0).
    """
    bounds, _families, _sig = PARTITIONS[partition]
    icv, ic_theta0, ftype, rmin_km, rmax_km = TABLE34[fam]
    target7 = np.array([*icv, ic_theta0])
    free = [k for k, (lo, hi) in enumerate(bounds) if hi > lo]
    span = np.array([bounds[k][1] - bounds[k][0] for k in free])
    lo_v = np.array([bounds[k][0] for k in free])
    target_norm = (target7[free] - lo_v) / span
    norm = (phen[:, free] - lo_v) / span
    diff = norm - target_norm
    if 6 in free:
        # theta0 (index 6) is periodic over its own bounds span -- a linear
        # normalized difference spuriously penalizes a genuine near-target
        # candidate that wrapped past 2*pi back toward 0 (found by Fable
        # review 2026-07-13). Wrap the normalized difference into [-0.5, 0.5]
        # (circular distance in unit-period coordinates) before combining.
        theta_col = free.index(6)
        diff[:, theta_col] -= np.round(diff[:, theta_col])
    dist = np.sqrt(np.mean(diff**2, axis=1))
    i = int(np.argmin(dist))
    cand_theta0 = float(phen[i, 6])
    cand = characterize(phen[i, :6], cand_theta0)
    got_rmin = float(cand["rmin_km_1yr"])
    got_rmax = float(cand["rmax_km_1yr"])
    ok_ic = float(dist[i]) < IC_PROXIMITY_TOL
    ok_type = cand["type"] == ftype
    ft = FEATURE_FACTOR_TOL
    ok_feat = 1.0 / ft <= got_rmin / rmin_km <= ft and 1.0 / ft <= got_rmax / rmax_km <= ft
    return {
        "family": fam,
        "partition": partition,
        "matched": bool(ok_ic and ok_type and ok_feat),
        "ic_rms_distance": float(dist[i]),
        "ic_ok": bool(ok_ic),
        "type_expected": ftype,
        "type_got": cand["type"],
        "type_ok": bool(ok_type),
        "rmin_ratio": got_rmin / rmin_km,
        "rmax_ratio": got_rmax / rmax_km,
        "features_ok": bool(ok_feat),
        "member_fitness": float(fitness[i]),
        "candidate": cand,
    }


def analyze_partition(name: str) -> None:
    _bounds, families, _sig = PARTITIONS[name]
    fpath = OUT_DIR / f"{name}_final.npz"
    if not fpath.exists():
        print(f"[{name}] final population missing, run to completion first")
        return
    data = np.load(fpath)
    phen, fitness = data["phenotypes"], data["fitness"]

    summary: dict[str, object] = {"partition": name, "families": {}}
    n_matched = 0
    for fam in families:
        rec = match_family_in_widened_population(fam, name, phen, fitness)
        summary["families"][fam] = rec  # type: ignore[index]
        n_matched += int(rec["matched"])
        print(
            f"[{name}] family {fam}: {'MATCH' if rec['matched'] else 'MISS '} "
            f"ic_dist={rec['ic_rms_distance']:.4f}({'ok' if rec['ic_ok'] else 'FAR'}) "
            f"type={rec['type_got']} vs {rec['type_expected']} "
            f"rmin_ratio={rec['rmin_ratio']:.3f} rmax_ratio={rec['rmax_ratio']:.3f}",
            flush=True,
        )
    summary["reproduction_rate"] = f"{n_matched}/{len(families)}"

    # Bounded-vs-divergent drift classification (#583's own new gate) on
    # every population member above a modest fitness floor -- the positive-
    # control question here is only "does the widened machinery still see
    # the known families", so this is diagnostic breadth, not a novelty claim.
    bounded_ids: list[int] = []
    checked = 0
    for idx, (row, fit) in enumerate(zip(phen, fitness, strict=True)):
        if fit < 0.9:  # skip low-fitness (unbounded/collided) genome noise
            continue
        checked += 1
        state = table_interleaved_to_state(row[:6])
        theta0 = float(row[6])
        v = classify_bounded_drift(state, theta0, n_revs=N_REVS_DEFAULT)
        if v.bounded:
            bounded_ids.append(idx)
    summary["drift_classified"] = {"checked": checked, "bounded": len(bounded_ids)}
    print(
        f"[{name}] drift classification: {len(bounded_ids)}/{checked} high-fitness "
        f"members bounded at {N_REVS_DEFAULT}yr",
        flush=True,
    )

    # Cheap theta0-robustness spot-check (per #583's own mandate) on up to 5
    # bounded survivors: re-test at 2 other phases and flag any flip. This is
    # diagnostic-only (does not itself gate anything) -- capped at 5 to keep
    # this analysis pass bounded regardless of how many members pass drift
    # classification.
    spot_checks: dict[str, object] = {}
    for idx in bounded_ids[:5]:
        row = phen[idx]
        state = table_interleaved_to_state(row[:6])
        theta0 = float(row[6])
        out = spot_check_theta0_robustness(state, theta0, n_revs=N_REVS_DEFAULT)
        flipped = not all(out.values())
        spot_checks[str(idx)] = {"per_phase_bounded": {str(k): v for k, v in out.items()}}
        print(
            f"[{name}] theta0 spot-check member {idx}: "
            f"{'FLIPPED' if flipped else 'stable'} across phases {list(out.values())}",
            flush=True,
        )
    summary["theta0_spot_checks"] = spot_checks

    out = OUT_DIR / f"{name}_analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[{name}] reproduction: {summary['reproduction_rate']}; summary -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--partition", choices=sorted(PARTITIONS), help="run/analyze one partition")
    ap.add_argument("--max-gens", type=int, default=None, help="generation cap this call")
    ap.add_argument("--workers", type=int, default=8, help="process-pool workers (M3: 8 cores)")
    ap.add_argument("--analyze", action="store_true", help="analyze a finished partition")
    args = ap.parse_args()
    preflight_search(
        task_no=583,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(PARTITIONS),
        override_reason=(
            "this dispatch runs ONLY the P1 positive-control partition to confirm "
            "the widened-bounds GA machinery still recovers #581 stage 2's own "
            "known families, per #583's own Fable-corrected gate; the full "
            "novel-territory sweep across all 7 partitions is explicitly deferred "
            "to a future dispatch, not an unbudgeted discovery sweep run here"
        ),
    )
    if args.partition is None:
        ap.error("specify --partition NAME (with --analyze or --max-gens)")
    if args.analyze:
        analyze_partition(args.partition)
    else:
        run_partition(args.partition, args.max_gens, args.workers)


if __name__ == "__main__":
    main()
