"""#581 stage 2: Gurfil-Kasdin (2002) niching-GA positive-control reproduction.

Runs the paper's 12 GA optimization sets (Table 2 bounds, Table 1 constants)
with the deterministic-crowding layer (``search/niching_ga.py``) on the
geocentric pulsating ER3BP (``core/er3bp_geocentric.py``), then compares the
final populations against the published 14 families A-N (Tables 3/4).

The paper's Table 2/3 vectors are printed in INTERLEAVED order
[x, x', y, y', z, z'] (established by r0/v0 cross-checks against Table 4; see
``core/er3bp_geocentric.py`` module docstring). The GA genome uses the same
interleaved slot order as the paper's tables.

Usage (chunked so every foreground call stays bounded):
    uv run python scripts/run_581_gurfil_reproduction.py --set 1 [--max-gens 100]
    uv run python scripts/run_581_gurfil_reproduction.py --analyze

Checkpoints/results: data/found/581_niching_ga/
"""

from __future__ import annotations

import argparse
import functools
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import cast

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core.er3bp import er3bp_eom
from cyclerfinder.core.er3bp_geocentric import (
    A_AU_KM_GURFIL_KASDIN,
    E_SUN_EARTH_GURFIL_KASDIN,
    MU_SUN_EARTH_GURFIL_KASDIN,
    R_EARTH_NORM,
    SUN_EARTH_ER3BP,
    geocentric_to_barycentric,
    gurfil_kasdin_fitness,
    table_interleaved_to_state,
)
from cyclerfinder.data.preflight import MethodCapability, preflight_search
from cyclerfinder.search.niching_ga import (
    DeterministicCrowdingConfig,
    run_deterministic_crowding,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "found" / "581_niching_ga"

_REGION_ID = "gurfil-kasdin-2002-geocentric-er3bp-581-stage2-positive-control-2026-07-12"
_METHOD = MethodCapability(
    genome=(
        "Gurfil & Kasdin (2002) geocentric pulsating ER3BP objective (Eq. 15, "
        "core/er3bp_geocentric.py) evaluated by deterministic-crowding niching GA "
        "(search/niching_ga.py) at the paper's own Table 1/2 constants and bounds"
    ),
    corrector=(
        "no corrector -- reproduces the paper's own published Table 3/4 families "
        "as a positive control for the niching mechanism, not a novel-orbit search"
    ),
    capability_tags=frozenset(
        {"er3bp", "geocentric", "niching-ga", "positive-control", "sun-earth"}
    ),
    git_sha="working-tree",
)

S2 = math.sqrt(2.0)
S3 = math.sqrt(3.0)
LO_R, HI_R = 0.0066845, 0.066845  # 1e6 km .. 1e7 km in AU (paper p. 5688)
V = 0.5  # velocity bound, AU/rad = 15 km/s (paper p. 5689)

# Table 2 (p. 5689), interleaved [x, x', y, y', z, z'] as printed.
# set -> (min6, max6, theta0, expected families)
TABLE2: dict[int, tuple[list[float], list[float], float, str]] = {
    1: ([LO_R, 0, 0, -V, 0, 0], [HI_R, 0, 0, 0, 0, 0], 0.0, "ABC"),
    2: ([LO_R, 0, 0, 0, 0, 0], [HI_R, 0, 0, V, 0, 0], 0.0, "D"),
    3: ([LO_R, 0, 0, 0, 0, 0], [HI_R, 0, 0, V, 0, 0], math.pi, "E"),
    4: ([LO_R, 0, 0, -V, 0, 0], [HI_R, 0, 0, 0, 0, 0], math.pi, "F"),
    5: (
        [v / S2 for v in [LO_R, 0, LO_R, 0, 0, 0]],
        [v / S2 for v in [HI_R, 0, HI_R, 0, 0, 0]],
        0.0,
        "G",
    ),
    6: (
        [v / S2 for v in [LO_R, -V, 0, 0, 0, 0]],
        [v / S2 for v in [HI_R, 0, 0, 0, 0, 0]],
        0.0,
        "H",
    ),
    7: (
        [v / S2 for v in [LO_R, -V, LO_R, -V, 0, 0]],
        [v / S2 for v in [HI_R, 0, HI_R, 0, 0, 0]],
        0.0,
        "I",
    ),
    8: (
        [v / S3 for v in [LO_R, -V, LO_R, -V, LO_R, -V]],
        [v / S3 for v in [HI_R, 0, HI_R, 0, HI_R, 0]],
        0.0,
        "J",
    ),
    9: (
        # Printed max slots 3/5 are "0.06684" (5 digits) in the paper; kept verbatim.
        [v / S3 for v in [LO_R, 0, LO_R, -V, -LO_R, 0]],
        [v / S3 for v in [HI_R, V, 0.06684, 0, 0.06684, V]],
        0.0,
        "K",
    ),
    10: (
        [v / S3 for v in [LO_R, 0, LO_R, 0, LO_R, 0]],
        [v / S3 for v in [HI_R, V, HI_R, V, HI_R, V]],
        0.0,
        "L",
    ),
    11: (
        [v / S3 for v in [LO_R, 0, LO_R, 0, LO_R, 0]],
        [v / S3 for v in [HI_R, 0, HI_R, 0, HI_R, 0]],
        0.0,
        "M",
    ),
    12: (
        [v / S2 for v in [LO_R, 0, 0, 0, LO_R, 0]],
        [v / S2 for v in [HI_R, 0, 0, 0, HI_R, 0]],
        0.0,
        "N",
    ),
}

# Table 3 published representative ICs (interleaved) + Table 4 features, for
# the --analyze comparison. family -> (icv, theta0, type, rmin_km, rmax_km)
TABLE34: dict[str, tuple[list[float], float, str, float, float]] = {
    "A": ([0.03894355345084, 0, 0, -0.07775997558556, 0, 0], 0.0, "DRO", 5769577, 11740892),
    "B": ([0.01482889291716, 0, 0, -0.03537804226749, 0, 0], 0.0, "DRO", 2191130, 3270091),
    "C": ([0.00680199460712, 0, 0, -0.02902265964752, 0, 0], 0.0, "DRO", 986996, 1126563),
    "D": ([0.00669917980717, 0, 0, 0.01073472190433, 0, 0], 0.0, "DPO", 317151, 1002197),
    "E": ([0.00669917980717, 0, 0, 0.01073472190433, 0, 0], math.pi, "DPO", 311893, 1002197),
    "F": ([0.00668449197861, 0, 0, -0.02343022812238, 0, 0], math.pi, "DRO", 224900, 1138177),
    "G": ([0.00472664960686, 0, 0.00472664960686, 0, 0, 0], 0.0, "ERO", 11532, 1084448),
    "H": ([0.00501336953102, -0.00097613025014, 0, 0, 0, 0], 0.0, "ERO", 15883, 750000),
    "I": (
        [0.01489874333797, -0.02624588842484, 0.01513377799728, -0.00995205804453, 0, 0],
        0.0,
        "DEO",
        2328213,
        3962617,
    ),
    "J": (
        [
            0.03348096835548,
            -0.00046191606162,
            0.00774766945226,
            -0.06652559750991,
            0.03675673393090,
            -0.00902011692574,
        ],
        0.0,
        "3D DRO",
        6892060,
        8510975,
    ),
    "K": (
        [
            0.00583817607709,
            0.00021213377191,
            0.00005468845617,
            -0.02914895441397,
            -0.00090918109209,
            0.00123524681998,
        ],
        0.0,
        "3D DRO",
        854531,
        939889,
    ),
    "L": (
        [
            0.00385929404386,
            0.00845656433011,
            0.00385933090466,
            0.00001057634437,
            0.00949726128195,
            0.00000748585020,
        ],
        0.0,
        "3D DEO",
        755223,
        1662761,
    ),
    "M": (
        [0.00386035674766, 0, 0.00385944561146, 0, 0.00386063960861, 0],
        0.0,
        "3D ERO",
        17282,
        1035854,
    ),
    "N": ([0.00766121767444, 0, 0, 0, 0.00668449197861, 0], 0.0, "3D DEO", 409758, 1521047),
}

# Family F: Table 4 prints rmax 1,002,197 km (suspected duplication of the D/E
# row value; inconsistent with F's own r0). The value above (1,138,177 km) is
# our high-accuracy re-integration of the published Table 3 IC, used only for
# candidate matching, and the discrepancy is reported in the results note.


def fitness_interleaved(vec6: np.ndarray, theta0: float) -> float:
    """GA fitness: decode the paper-order genome and score Eq. 15."""
    state = table_interleaved_to_state(np.asarray(vec6, dtype=float))
    return gurfil_kasdin_fitness(
        state, theta0, SUN_EARTH_ER3BP, n_rev=1.0, rtol=1e-9, atol=1e-9, n_samples=2000
    )


def run_set(set_no: int, max_gens: int | None, workers: int) -> None:
    lo, hi, theta0, expected = TABLE2[set_no]
    bounds = list(zip(lo, hi, strict=True))
    config = DeterministicCrowdingConfig(seed=581000 + set_no)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ckpt = OUT_DIR / f"set{set_no:02d}_checkpoint.npz"
    runlog = OUT_DIR / f"set{set_no:02d}_runlog.jsonl"
    fit = functools.partial(fitness_interleaved, theta0=theta0)
    t_start = time.monotonic()

    def progress(stats: dict[str, float]) -> None:
        rec = {
            "ts": datetime.now(UTC).isoformat(timespec="seconds"),
            "set": set_no,
            **stats,
            "elapsed_s": round(time.monotonic() - t_start, 1),
        }
        with runlog.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        gen = int(stats["generation"])
        if gen % 20 == 0 or gen == config.generations:
            print(
                f"[set {set_no}] gen {gen}/{config.generations} "
                f"mean={stats['fitness_mean']:.6f} max={stats['fitness_max']:.6f} "
                f"std={stats['fitness_std']:.6f} elapsed={rec['elapsed_s']}s",
                flush=True,
            )

    result = run_deterministic_crowding(
        fit,
        bounds,
        config,
        workers=workers,
        checkpoint_path=ckpt,
        max_generations_this_call=max_gens,
        progress_fn=progress,
    )
    print(
        f"[set {set_no}] done through generation {result.generations_run} "
        f"(expected families: {expected})",
        flush=True,
    )
    if result.generations_run >= config.generations:
        np.savez(
            OUT_DIR / f"set{set_no:02d}_final.npz",
            phenotypes=result.phenotypes,
            fitness=result.fitness,
            theta0=theta0,
        )
        print(f"[set {set_no}] final population saved", flush=True)


# ---------------------------------------------------------------------------
# Analysis: cluster final populations, characterize, match to Tables 3/4
# ---------------------------------------------------------------------------


def characterize(interleaved: np.ndarray, theta0: float) -> dict[str, object]:
    """High-accuracy 5-year characterization of one candidate IC."""
    state0 = table_interleaved_to_state(interleaved)
    bary0 = geocentric_to_barycentric(state0, MU_SUN_EARTH_GURFIL_KASDIN)
    offset = 1.0 - MU_SUN_EARTH_GURFIL_KASDIN

    def collision(f: float, s: np.ndarray, *_a: float) -> float:
        dx = s[0] - offset
        return float(dx * dx + s[1] ** 2 + s[2] ** 2 - R_EARTH_NORM**2)

    collision.terminal = True  # type: ignore[attr-defined]

    def escape(f: float, s: np.ndarray, *_a: float) -> float:
        dx = s[0] - offset
        return float(dx * dx + s[1] ** 2 + s[2] ** 2 - 0.25)  # 0.5 AU

    escape.terminal = True  # type: ignore[attr-defined]

    spans = {"1yr": 2.0 * math.pi, "5yr": 10.0 * math.pi}
    out: dict[str, object] = {}
    for label, dspan in spans.items():
        sol = solve_ivp(
            er3bp_eom,
            (theta0, theta0 + dspan),
            bary0,
            args=(MU_SUN_EARTH_GURFIL_KASDIN, E_SUN_EARTH_GURFIL_KASDIN),
            method="DOP853",
            rtol=1e-11,
            atol=1e-11,
            max_step=0.01,
            events=(collision, escape),
        )
        dx = sol.y[0] - offset
        r = np.sqrt(dx * dx + sol.y[1] ** 2 + sol.y[2] ** 2)
        out[f"rmin_km_{label}"] = float(r.min() * A_AU_KM_GURFIL_KASDIN)
        out[f"rmax_km_{label}"] = float(r.max() * A_AU_KM_GURFIL_KASDIN)
        out[f"terminated_{label}"] = bool(sol.status == 1)
        out[f"theta_end_{label}"] = float(sol.t[-1] - theta0)

    _x0, xp0, _y0, yp0, z0, zp0 = (float(v) for v in interleaved)
    spatial = abs(z0) > 1e-9 or abs(zp0) > 1e-9
    v0 = math.sqrt(xp0**2 + yp0**2 + zp0**2)
    # PS: never beyond 0.1 AU from Earth within 5 years (paper p. 5689).
    ps = not out["terminated_5yr"] and cast(float, out["rmax_km_5yr"]) <= (
        0.1 * A_AU_KM_GURFIL_KASDIN
    )
    if v0 < 2e-3:  # < ~0.06 km/s: the paper's "zero initial velocity" class
        kind = "ERO"
    elif yp0 < 0:
        kind = "DRO"
    else:
        kind = "DPO"
    if not ps:
        kind = "DEO"  # paper: PUS orbits are DEOs regardless of other traits
    out.update(
        {
            "ic_interleaved": [float(v) for v in interleaved],
            "theta0": theta0,
            "spatial": spatial,
            "v0_norm": v0,
            "practically_stable": ps,
            "type": ("3D " if spatial else "") + kind,
        }
    )
    return out


def cluster_population(
    phen: np.ndarray, fitness: np.ndarray, bounds: list[tuple[float, float]], tol: float = 0.03
) -> list[tuple[np.ndarray, float, int]]:
    """Greedy fitness-ranked clustering in bounds-normalized gene space.

    Returns [(representative_phenotype, fitness, member_count)] sorted by
    descending fitness.
    """
    free = [k for k, (lo, hi) in enumerate(bounds) if hi > lo]
    span = np.array([bounds[k][1] - bounds[k][0] for k in free])
    lo_v = np.array([bounds[k][0] for k in free])
    norm = (phen[:, free] - lo_v) / span
    order = np.argsort(-fitness)
    reps: list[tuple[np.ndarray, float, int]] = []
    rep_norms: list[np.ndarray] = []
    counts: list[int] = []
    for i in order:
        placed = False
        for ci, rn in enumerate(rep_norms):
            if np.linalg.norm(norm[i] - rn) < tol * math.sqrt(len(free)):
                counts[ci] += 1
                placed = True
                break
        if not placed:
            reps.append((phen[i], float(fitness[i]), 0))
            rep_norms.append(norm[i])
            counts.append(1)
    return [(r[0], r[1], c) for r, c in zip(reps, counts, strict=True)]


# Match criteria for "recognizably reproduced" (stochastic GA, not bit-exact):
# the population member nearest the published IC must (a) lie within
# IC_PROXIMITY_TOL of it in bounds-normalized RMS gene distance (same niche),
# (b) re-integrate to the same paper type, and (c) have 1-year rmin AND rmax
# within a factor of FEATURE_FACTOR_TOL of Table 4's values.
IC_PROXIMITY_TOL = 0.10
FEATURE_FACTOR_TOL = 2.0


def match_family_in_population(
    fam: str,
    set_no: int,
    phen: np.ndarray,
    fitness: np.ndarray,
) -> dict[str, object]:
    """Family-directed match: nearest niche member to the published IC."""
    lo, hi, theta0, _ = TABLE2[set_no]
    icv, _ic_theta0, ftype, rmin_km, rmax_km = TABLE34[fam]
    bounds = list(zip(lo, hi, strict=True))
    free = [k for k, (l_, h_) in enumerate(bounds) if h_ > l_]
    span = np.array([bounds[k][1] - bounds[k][0] for k in free])
    lo_v = np.array([bounds[k][0] for k in free])
    target = (np.array([icv[k] for k in free]) - lo_v) / span
    norm = (phen[:, free] - lo_v) / span
    dist = np.sqrt(np.mean((norm - target) ** 2, axis=1))  # normalized RMS
    i = int(np.argmin(dist))
    cand = characterize(phen[i], theta0)
    got_rmin = cast(float, cand["rmin_km_1yr"])
    got_rmax = cast(float, cand["rmax_km_1yr"])
    ok_ic = float(dist[i]) < IC_PROXIMITY_TOL
    ok_type = cand["type"] == ftype
    ft = FEATURE_FACTOR_TOL
    ok_feat = 1.0 / ft <= got_rmin / rmin_km <= ft and 1.0 / ft <= got_rmax / rmax_km <= ft
    return {
        "family": fam,
        "set": set_no,
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


def analyze() -> None:
    summary: dict[str, object] = {"sets": {}, "family_matches": {}}
    matches: dict[str, dict[str, object]] = {}
    for set_no, (lo, hi, _theta0, expected) in TABLE2.items():
        fpath = OUT_DIR / f"set{set_no:02d}_final.npz"
        if not fpath.exists():
            print(f"[set {set_no}] final population missing, skipping")
            continue
        data = np.load(fpath)
        phen, fitness = data["phenotypes"], data["fitness"]
        bounds = list(zip(lo, hi, strict=True))
        clusters = cluster_population(phen, fitness, bounds)
        summary["sets"][str(set_no)] = {  # type: ignore[index]
            "expected": list(expected),
            "n_clusters": len(clusters),
            "n_nonzero_fitness": int(np.sum(fitness > 0)),
        }
        for fam in expected:
            rec = match_family_in_population(fam, set_no, phen, fitness)
            matches[fam] = rec
            print(
                f"[set {set_no}] family {fam}: {'MATCH' if rec['matched'] else 'MISS '} "
                f"ic_dist={rec['ic_rms_distance']:.4f}({'ok' if rec['ic_ok'] else 'FAR'}) "
                f"type={rec['type_got']} vs {rec['type_expected']} "
                f"rmin_ratio={rec['rmin_ratio']:.3f} rmax_ratio={rec['rmax_ratio']:.3f}",
                flush=True,
            )
    summary["family_matches"] = {k: matches[k] for k in sorted(matches)}
    n = sum(1 for r in matches.values() if r["matched"])
    summary["reproduction_rate"] = f"{n}/14"
    out = OUT_DIR / "analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    reproduced = sorted(k for k, r in matches.items() if r["matched"])
    print(f"\nReproduction: {n}/14 families matched: {reproduced}")
    print(f"Missing: {sorted(set(TABLE34) - set(reproduced))}")
    print(f"Summary written to {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", type=int, choices=sorted(TABLE2), help="run one optimization set")
    ap.add_argument("--max-gens", type=int, default=None, help="generation cap this call")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--analyze", action="store_true", help="analyze finished sets")
    args = ap.parse_args()
    preflight_search(
        task_no=581,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(TABLE2),
        override_reason=(
            "positive-control reproduction of Gurfil & Kasdin (2002)'s own published "
            "Table 3/4 families -- validates the niching-GA mechanism against a known "
            "answer, not an unbudgeted discovery sweep; #581 stage 2 explicitly gated "
            "from stage 3 (novel-target search)"
        ),
    )
    if args.analyze:
        analyze()
    elif args.set is not None:
        run_set(args.set, args.max_gens, args.workers)
    else:
        ap.error("specify --set N or --analyze")


if __name__ == "__main__":
    main()
