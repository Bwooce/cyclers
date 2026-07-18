"""Task #633: fixed-Saturn-Titan-mu 2D (x0, C, hc) corridor grid search.

`#629`'s design read (docs/notes/2026-07-18-629-design-read-titan-kk-grid.md)
found `#627`'s 1D mu-continuation of the Ross-Roberts-Tsoukkas 2026
(k1,k2)=(1,1)/(3,3) ballistic-cycler families down to Saturn-Titan mu was
suspect due to an unscaled margin parameter. Phase A (`scripts/
run_629_titan_scaled_margin_phase_a.py`) reran the walk with a properly
SCALED margin across 6 alpha values and STILL found the target topology is
lost or the walk fails to converge -- decisive evidence that a 1D
continuation PATH (regardless of scaling) is unsuited to this problem: a
genuine hard fold exists at the (1,1) anchor itself, and both families lose
their target topology via different failure modes along the walk.

This script (Phase B) sidesteps path-dependence entirely: instead of
continuing FROM a distant Earth-Moon-scale anchor, it searches DIRECTLY at
the target mu via a 2D grid in (x0, C) (plus a discrete half-crossings
index), reusing the SAME machinery already validated for the Pluto-Charon
(2,1)/(2,2) grid searches
(:func:`cyclerfinder.search.pluto_charon_kk_sweep._grid_seed_search`,
:func:`c_sweep_find_nu_zero`) -- the new content here is coordinate scaling
(C expressed as a corridor-fraction rho, not a new method) plus the
per-point sweep loop, which (unlike ``_grid_seed_search``) does NOT early-
exit at the first topology match: this is a CENSUS over the whole grid, not
a single-seed search, so every point is evaluated and logged.

Grid (per the #629 design read's own spec, ~16,375 points):
  x0  in [-0.95, -0.30], step 0.005                       (131 points)
  rho in [0.6, 0.99], ~18 coarse + a FINE [0.955, 0.995]   (25 points)
      sub-band 7 points (without the fine sub-band, thin-corridor
      families like (3,3) are missed by construction -- the (3,3) anchor's
      own rho is 0.974, ABOVE the coarse band's ceiling)
      C = 3 + rho * (C_L1(Titan) - 3)
  hc  in {1, 3, 5, 7, 9}                                    (5 points)

For every grid point that clears the initial residual/spatial-closure gate
(the fixed-Jacobi corrector converges), the resulting orbit's winding
topology is classified. Running the further stability (Barden nu, escalating
to a c_sweep_find_nu_zero window search if the direct point is not already
stable) and encounter-relevance (perimoon_passage.py) gates on EVERY
converged point would be enormously more expensive (each c_sweep call is up
to 60 more corrector calls) for no benefit -- those extra gates are only
run for points whose topology ALREADY matches one of the two targets
((1,1) or (3,3)) AND reaches_secondary, exactly mirroring the escalation
order #627/#629's own scripts already use (try the direct point first, only
escalate to a C-sweep if that fails). This choice is a documented judgment
call (see the #633 OUTSTANDING.md bullet), not part of the design read's
literal text, made because a literal per-grid-point c_sweep would multiply
the total corrector-call budget by ~60x for no plausible benefit -- matching
topology is a hard, cheap-to-check necessary condition, and the design
read's own analysis is that ANY match at all would already be a rare event.

Timing pilot (see the #633 OUTSTANDING.md bullet for the measured numbers):
per-point cost (corrector + topology classification) is ~0.1-0.5s typically,
worst case ~4.4s observed in an edge-case probe near the L1 neck -- NOT the
1994s/355s Phase-A per-continuation-STEP costs (a single grid-point
correction is vastly cheaper than a 200-step continuation walk with
per-step C-sub-walks). A 20s per-call SIGALRM timeout is generous (>4x the
observed worst case) while still bounding any pathological point.

Runs in CHUNKS via --start-idx/--end-idx (a flattened index into the
deterministic (hc, rho, x0) grid ordering) so each invocation fits
comfortably inside a single foreground shell call, appending to the SAME
JSONL checkpoint file across chunks (per
[[feedback_incremental_progress_reports]] -- #520's silent-hours-zero-output
failure mode). Use --summarize to read back the accumulated JSONL and print
the aggregate census (converged fraction, topology histogram, candidates)
without re-running anything.

No catalogue writeback under any circumstances (see the #633 OUTSTANDING.md
bullet, step 10): any genuine on-topology/stable/encounter-relevant find
must stop here and be flagged for a Phase C Opus/Fable adjudication pass.

Usage
-----
  uv run python scripts/run_633_titan_corridor_grid.py --start-idx 0 --end-idx 1200
  uv run python scripts/run_633_titan_corridor_grid.py --start-idx 1200 --end-idx 2400
  ...
  uv run python scripts/run_633_titan_corridor_grid.py --summarize
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

_SRC = Path(__file__).parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.cr3bp_periodic as cp  # noqa: E402
from cyclerfinder.core.cr3bp import cr3bp_system  # noqa: E402
from cyclerfinder.core.satellites import SATELLITES  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.search.binary_star_search import winding_topology  # noqa: E402
from cyclerfinder.search.literature_check import (  # noqa: E402
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)
from cyclerfinder.search.perimoon_passage import find_perimoon_passage  # noqa: E402
from cyclerfinder.search.physical_sanity import flyby_is_useful  # noqa: E402
from cyclerfinder.search.pluto_charon_kk_sweep import (  # noqa: E402
    _c_l1,
    _run_with_timeout,
    c_sweep_find_nu_zero,
)

_REGION_ID = "rrt-kk-titan-corridor-grid-2026-07-18"
OUT_PATH = (
    Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "633_titan_corridor_grid.jsonl"
)

# --- Grid definition (per the #629 design read's own spec) ------------------
X0_LO, X0_HI, X0_STEP = -0.95, -0.30, 0.005
RHO_COARSE_LO, RHO_COARSE_HI, RHO_COARSE_N = 0.60, 0.95, 18
RHO_FINE_LO, RHO_FINE_HI, RHO_FINE_N = 0.955, 0.995, 7
HC_LIST: tuple[int, ...] = (1, 3, 5, 7, 9)

# Timing-pilot-measured constants (see #633 OUTSTANDING.md bullet for the raw
# pilot output): mean ~0.4s/point including topology classification, worst
# observed 4.4s near the L1 neck. Generous 20s SIGALRM margin (>4x worst
# case) bounds any pathological point without materially inflating the
# expected total runtime.
PER_CALL_TIMEOUT_S = 20
PERIOD_GUESS = 20.0  # TU; generous t_hi=1.25*this comfortably spans hc<=9
TIMING_PILOT_SECONDS_PER_POINT = 0.45  # measured, see bullet; conservative round-up

TARGETS: dict[str, tuple[int, int]] = {"11": (1, 1), "33": (3, 3)}


def build_grids() -> tuple[np.ndarray, np.ndarray, tuple[int, ...]]:
    x0_grid = np.round(np.arange(X0_LO, X0_HI + 1e-9, X0_STEP), 6)
    rho_coarse = np.linspace(RHO_COARSE_LO, RHO_COARSE_HI, RHO_COARSE_N)
    rho_fine = np.linspace(RHO_FINE_LO, RHO_FINE_HI, RHO_FINE_N)
    rho_grid = np.unique(np.round(np.concatenate([rho_coarse, rho_fine]), 6))
    return x0_grid, rho_grid, HC_LIST


def flatten_points(
    x0_grid: np.ndarray, rho_grid: np.ndarray, hc_list: tuple[int, ...]
) -> list[tuple[int, float, float, int]]:
    """Deterministic flattened (idx, x0, rho, hc) list -- hc outer, rho, x0 inner."""
    out = []
    idx = 0
    for hc in hc_list:
        for rho in rho_grid:
            for x0 in x0_grid:
                out.append((idx, float(x0), float(rho), int(hc)))
                idx += 1
    return out


def _append_jsonl(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(json.dumps(record) + "\n")


def _direct_gate(
    system: cr3bp.CR3BPSystem, orbit: cp.SymmetricOrbit, k1: int, k2: int
) -> dict[str, Any]:
    """Barden stability + topology at the AS-LANDED grid point, no C-walk."""
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo = winding_topology(system.mu, state0, orbit.period)
    nu, _ = cp.barden_stability(system, orbit, rtol=1e-13, atol=1e-13)
    stability_pass = bool(
        abs(nu) < 1.0 and topo.k1 == k1 and topo.k2 == k2 and topo.reaches_secondary
    )
    return {
        "topology": (topo.k1, topo.k2),
        "prograde": topo.prograde,
        "reaches_secondary": topo.reaches_secondary,
        "nu_direct": float(nu),
        "gate_stability_pass": stability_pass,
        "state0": state0,
    }


def _escalate_c_sweep(
    system: cr3bp.CR3BPSystem,
    orbit: cp.SymmetricOrbit,
    k1: int,
    k2: int,
    c_l1: float,
) -> dict[str, Any] | None:
    """C-sweep escalation for a topology-matching-but-unstable direct point.

    Band scaled to the Titan corridor's own (tiny, mu^(2/3)-shrunk) width
    rather than the Earth-Moon-scale +/-0.3 ``sweep_family_grid`` uses --
    that fixed offset would walk C far below C=3 (outside the ballistic-
    capture band entirely) for this corridor, per the #629 design read's own
    mu^(2/3)-scaling finding.
    """
    band = 0.5 * (c_l1 - 3.0)
    c_lo = max(3.0 + 1e-9, orbit.jacobi - band)
    c_hi = min(c_l1 - 1e-4, orbit.jacobi + band)
    if c_lo >= c_hi:
        return None
    stable = c_sweep_find_nu_zero(
        system,
        orbit.x0,
        orbit.jacobi,
        orbit.period,
        hc=None,
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    if stable is None:
        return None
    state0 = np.array([stable.x0, 0.0, 0.0, 0.0, stable.ydot0, 0.0])
    topo = winding_topology(system.mu, state0, stable.period)
    nu, _ = cp.barden_stability(system, stable, rtol=1e-13, atol=1e-13)
    if not (abs(nu) < 1.0 and topo.k1 == k1 and topo.k2 == k2 and topo.reaches_secondary):
        return None
    return {
        "orbit": stable,
        "topology": (topo.k1, topo.k2),
        "prograde": topo.prograde,
        "reaches_secondary": topo.reaches_secondary,
        "nu_direct": float(nu),
        "state0": state0,
    }


def offline_corpus_search(query: str) -> list[SearchResult]:
    """Deterministic offline literature backend (mirrors #627/#629's own use)."""
    q = query.lower()
    out: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        hit = any(a.lower() in q for a in anchor.authors) or any(
            kw.lower() in q for kw in anchor.keywords
        )
        bodies_named = sum(1 for b in anchor.body_set if b.lower() in q)
        if hit or (bodies_named >= 2 and "cycler" in q):
            bodies = " ".join(sorted(anchor.body_set))
            out.append(
                SearchResult(
                    title=anchor.name,
                    url=anchor.doi or "offline-corpus",
                    snippet=f"{anchor.citation} ({bodies})",
                )
            )
    return out


def evaluate_point(
    system: cr3bp.CR3BPSystem, idx: int, x0: float, rho: float, hc: int, c_l1: float
) -> dict[str, Any]:
    c_val = 3.0 + rho * (c_l1 - 3.0)
    record: dict[str, Any] = {"idx": idx, "x0": x0, "rho": rho, "c": c_val, "hc": hc}

    def _fn(_x0: float = x0, _c: float = c_val, _hc: int = hc) -> cp.SymmetricOrbit:
        return cp.correct_symmetric_fixed_jacobi(
            system, _x0, _c, PERIOD_GUESS, ydot0_sign=-1.0, half_crossings=_hc, tol=1e-10
        )

    try:
        o = _run_with_timeout(_fn, seconds=PER_CALL_TIMEOUT_S)
    except (ValueError, RuntimeError) as e:
        record["converged"] = False
        record["error"] = str(e)
        return record

    if o is None or not o.converged:
        record["converged"] = False
        if o is not None:
            record["crossing_residual"] = o.crossing_residual
        return record

    record["converged"] = True
    record["crossing_residual"] = o.crossing_residual
    record["period"] = o.period

    state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
    topo = winding_topology(system.mu, state0, o.period)
    record["k1"] = topo.k1
    record["k2"] = topo.k2
    record["prograde"] = topo.prograde
    record["reaches_secondary"] = topo.reaches_secondary

    target_match = None
    for key, (k1, k2) in TARGETS.items():
        if topo.k1 == k1 and topo.k2 == k2 and topo.reaches_secondary:
            target_match = key
            break
    record["target_match"] = target_match

    if target_match is None:
        return record

    # --- Candidate: chain the stability gate (direct point first, then a
    # C-sweep escalation if the direct point isn't already stable) ---
    k1, k2 = TARGETS[target_match]
    direct = _direct_gate(system, o, k1, k2)
    record["nu_direct"] = direct["nu_direct"]
    record["gate_stability_pass_direct"] = direct["gate_stability_pass"]

    final_orbit = o
    stable_info = direct
    if not direct["gate_stability_pass"]:
        esc = _escalate_c_sweep(system, o, k1, k2, c_l1)
        record["escalated_c_sweep"] = True
        record["escalation_found_stable"] = esc is not None
        if esc is None:
            record["gate_stability_pass"] = False
            return record
        final_orbit = esc["orbit"]
        stable_info = esc
        record["c_sweep_x0"] = final_orbit.x0
        record["c_sweep_jacobi"] = final_orbit.jacobi
        record["c_sweep_nu"] = stable_info["nu_direct"]
    else:
        record["escalated_c_sweep"] = False

    record["gate_stability_pass"] = True

    # --- Gate (b): encounter-relevance (perimoon-passage geometry) ---
    titan_radius_km = SATELLITES["Titan"].radius_eq_km
    passage = find_perimoon_passage(
        system, stable_info["state0"], final_orbit.period, titan_radius_km, rtol=1e-12, atol=1e-12
    )
    verdict = flyby_is_useful("Titan", passage.speed_rel_kms)
    record["perimoon"] = {
        "altitude_km": passage.altitude_km,
        "speed_rel_kms": passage.speed_rel_kms,
        "below_surface": passage.below_surface,
        "max_bend_deg": verdict.max_bend_deg,
        "is_useful_bend": verdict.is_useful,
    }
    record["gate_encounter_pass"] = bool((not passage.below_surface) and verdict.is_useful)

    if record["gate_encounter_pass"]:
        sig = CandidateSignature(
            primary="Saturn",
            sequence=("Titan",),
            topology_label=frozenset({"repeated-moon"}),
            vinf_per_encounter_kms=(passage.speed_rel_kms,),
        )
        lit = check_literature(sig, search=offline_corpus_search)
        record["lit_check_offline"] = {
            "status": lit.status,
            "citation": lit.citation,
            "confidence": lit.confidence,
        }

    return record


def run_chunk(start_idx: int, end_idx: int, out_path: Path) -> None:
    x0_grid, rho_grid, hc_list = build_grids()
    points = flatten_points(x0_grid, rho_grid, hc_list)
    n_total = len(points)
    end_idx = min(end_idx, n_total)
    system = cr3bp_system("Saturn", "Titan")
    c_l1 = _c_l1(system.mu)

    print(
        f"Task #633 Titan corridor grid: chunk [{start_idx}, {end_idx}) of {n_total} total points "
        f"(mu={system.mu:.8g}, C_L1={c_l1:.6f})",
        flush=True,
    )

    n_converged = 0
    n_candidates = 0
    n_stable_ontopology = 0
    t0 = time.time()
    for i in range(start_idx, end_idx):
        _idx, x0, rho, hc = points[i]
        rec = evaluate_point(system, i, x0, rho, hc, c_l1)
        _append_jsonl(out_path, rec)
        if rec.get("converged"):
            n_converged += 1
        if rec.get("target_match") is not None:
            n_candidates += 1
            print(
                f"  [{i}] CANDIDATE x0={x0:.4f} rho={rho:.4f} hc={hc} "
                f"target={rec['target_match']} gate_stability_pass="
                f"{rec.get('gate_stability_pass')} gate_encounter_pass="
                f"{rec.get('gate_encounter_pass')}",
                flush=True,
            )
        if rec.get("gate_stability_pass") and rec.get("target_match") is not None:
            n_stable_ontopology += 1
        if (i - start_idx + 1) % 200 == 0 or i == end_idx - 1:
            elapsed = time.time() - t0
            done = i - start_idx + 1
            rate = elapsed / done if done else 0.0
            print(
                f"  progress {done}/{end_idx - start_idx} "
                f"(global idx {i}/{n_total - 1}) elapsed={elapsed:.1f}s "
                f"rate={rate:.3f}s/pt converged={n_converged} candidates={n_candidates} "
                f"stable_ontopology={n_stable_ontopology}",
                flush=True,
            )

    elapsed = time.time() - t0
    print(
        f"Chunk done: {end_idx - start_idx} points in {elapsed:.1f}s "
        f"({elapsed / max(end_idx - start_idx, 1):.3f}s/pt). "
        f"converged={n_converged} candidates={n_candidates} "
        f"stable_ontopology={n_stable_ontopology}",
        flush=True,
    )


def summarize(out_path: Path) -> None:
    x0_grid, rho_grid, hc_list = build_grids()
    n_total = len(x0_grid) * len(rho_grid) * len(hc_list)
    if not out_path.exists():
        print(f"No checkpoint file at {out_path} yet.")
        return
    n_seen = 0
    n_converged = 0
    topo_hist: dict[str, int] = {}
    candidates: list[dict[str, Any]] = []
    with out_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            n_seen += 1
            if rec.get("converged"):
                n_converged += 1
                k = f"({rec.get('k1')},{rec.get('k2')})"
                topo_hist[k] = topo_hist.get(k, 0) + 1
            if rec.get("target_match") is not None:
                candidates.append(rec)

    pct_seen = 100 * n_seen / n_total
    print(f"Checkpoint {out_path}: {n_seen}/{n_total} points evaluated ({pct_seen:.1f}%)")
    pct_conv = 100 * n_converged / max(n_seen, 1)
    print(f"Converged: {n_converged}/{n_seen} ({pct_conv:.1f}%)")
    print("Topology histogram (top 15):")
    for k, v in sorted(topo_hist.items(), key=lambda kv: -kv[1])[:15]:
        print(f"  {k}: {v}")
    print(
        f"\nTarget-topology candidates ((1,1) or (3,3), reaches_secondary=True): {len(candidates)}"
    )
    stable_ontopology = [c for c in candidates if c.get("gate_stability_pass")]
    print(f"Of those, stable + on-topology (gate_stability_pass): {len(stable_ontopology)}")
    encounter_pass = [c for c in stable_ontopology if c.get("gate_encounter_pass")]
    print(f"Of those, encounter-relevant (gate_encounter_pass): {len(encounter_pass)}")
    for c in candidates:
        print(f"  CANDIDATE: {c}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-idx", type=int, default=None)
    parser.add_argument("--end-idx", type=int, default=None)
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--out", type=str, default=str(OUT_PATH))
    args = parser.parse_args()

    out_path = Path(args.out)

    if args.summarize:
        summarize(out_path)
        return

    if args.start_idx is None or args.end_idx is None:
        parser.error("--start-idx and --end-idx are required unless --summarize is given")

    x0_grid, rho_grid, hc_list = build_grids()
    n_total = len(x0_grid) * len(rho_grid) * len(hc_list)

    preflight_search(
        task_no=633,
        region_id=_REGION_ID,
        method=MethodCapability(
            genome=(
                "Ross-Roberts-Tsoukkas 2026 (k1,k2)-cycler CR3BP genome (#494's fixed-Jacobi "
                "symmetric corrector + Barden stability + winding-topology classifier), "
                "searched DIRECTLY at Saturn-Titan mu=2.36695e-4 via a 2D (x0, C) grid in "
                "corridor-scaled coordinates (C=3+rho*(C_L1-3)) plus a discrete "
                "half-crossings index -- #629 Phase B, sidesteps the continuation-path-"
                "dependence #627/#629-Phase-A found"
            ),
            corrector=(
                "correct_symmetric_fixed_jacobi grid scan (reusing "
                "pluto_charon_kk_sweep._grid_seed_search's primitive, generalized to a "
                "full-grid CENSUS rather than a first-hit search) + c_sweep_find_nu_zero "
                "stability escalation + perimoon_passage.py encounter-relevance gate"
            ),
            capability_tags=frozenset(
                {
                    "cr3bp",
                    "binary-cycler",
                    "k1k2-genome",
                    "real-moon",
                    "2d-grid",
                    "fixed-mu-direct-search",
                }
            ),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=n_total,
        timing_pilot_seconds_per_point=TIMING_PILOT_SECONDS_PER_POINT,
    )

    run_chunk(args.start_idx, args.end_idx, out_path)


if __name__ == "__main__":
    main()
