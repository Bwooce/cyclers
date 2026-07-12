"""#573 -- widened Titan-Iapetus 3D-closure population probe.

Generalizes ``scripts/probe_572_titan_iapetus_3d_closure.py`` (already
validated: correct inc=0 coplanar-reduction smoke test, correct
node-alignment search structure, correct DOP853 cross-check baked into the
upstream coplanar sweep) from the 2 hand-picked #572 candidates to the FULL
population of 69 Titan-anchored, coplanar-gate-passing candidates in
``data/scan_571_saturn_titan_iapetus.jsonl`` (``physical_gate_passed=True``
records). See ``data/OUTSTANDING.md`` #573 (including the Fable
plan-review corrections appended under it) for the authoritative scope --
this docstring only summarizes.

This script REUSES the core 3D Lambert-closure machinery from #572
unmodified (``iapetus_state_3d``, ``_leg_best``, the rotation algebra, the
node-alignment idea) and adds the at-scale guard logic the Fable review
found load-bearing:

1. Branch dedup -- the 69 records are grid samples along ~10 continuous
   solution branches, not 69 independent candidates. This script clusters
   CLOSING basins found across all 69 candidate runs by proximity in
   (n_rev, Omega*, tof_scale*, V_inf_mid) via union-find, and reports
   results in DISTINCT BRANCHES. It also groups the raw 69 seeds by
   (n_rev, tof_scale) for the singleton-flip check (item 6).
2. Adaptive Nelder-Mead Omega window: min(15 deg, half-distance to the
   nearest adjacent grid-minima Omega, periodic) instead of a fixed +-15,
   to avoid branch-jumping when basins are closely spaced. Post-refinement,
   verify the solution's V_inf stayed within ~0.5 km/s of its seed basin's
   grid V_inf (flag, don't silently drop, on drift).
3. Basin refinement capped at 12 deepest basins/candidate (matches #572);
   logs explicitly whenever a candidate has >12 grid basins (cap binds).
4. Nelder-Mead runs that pin at a box edge (Omega or tof_scale bound) are
   recorded as "formulation-conditional non-closures" -- excluded from the
   genuine-closure branch count but tracked and reported separately.
5. LambertGeometryError is tracked SEPARATELY from LambertConvergenceError
   (the #572 run never distinguished them because neither was ever raised
   at only 2 candidates x 3600 samples) with a small perturbed-offset retry
   before giving up on a grid point.
6. Isolated closure "flips" within an otherwise-contiguous (n_rev,
   tof_scale) rel_offset run are flagged via
   ``cyclerfinder.data.sweep_diagnostics.detect_isolated_singleton_anomalies``
   rather than counted at face value, per the project's standing
   sweep-singleton-artifact rule.
7. Eccentricity-robust sub-count: among distinct closing branches, count
   how many have a realized Iapetus (flyby-body, the MIDDLE value of the
   3-element ``max_bend_deg_per_encounter`` array) bend >= 6.0 deg.

Run as::

    uv run python scripts/run_573_titan_iapetus_population_closure.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from probe_572_titan_iapetus_3d_closure import (  # noqa: E402
    ANCHOR,
    FLYBY,
    INCLINATION_DEG,
    PRIMARY,
    _leg_best,
    _moon_state,
    _smoke_test_reduction,
    iapetus_state_3d,
)
from scan_558_uranus_all_pairs_offset_sweep import GATE_RESIDUAL_KMS  # noqa: E402

from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.data.sweep_diagnostics import detect_isolated_singleton_anomalies  # noqa: E402
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day  # noqa: E402
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

DATA_DIR = ROOT / "data"
SCAN_571_PATH = DATA_DIR / "scan_571_saturn_titan_iapetus.jsonl"
OUT_PATH = DATA_DIR / "probe_573_titan_iapetus_population_closure.jsonl"

N_OMEGA = 3600
BASIN_CAP = 12
ECC_ROBUST_BEND_DEG = 6.0

# Branch-dedup proximity thresholds (union-find merge criterion across
# closing basins found from different seed candidates). Chosen to be a few
# grid steps wide in each axis (grid steps: tof_scale ~0.05, rel_offset
# ~1-5 deg) while staying well inside "two distinct physical solutions
# would differ by more than this."
MERGE_OMEGA_DEG = 5.0
MERGE_TOF_SCALE = 0.05
MERGE_VINF_KMS = 0.10

# Basin-support drift flag: refined V_inf (mid, Iapetus) must stay within
# this of the *seed grid point's* V_inf to be trusted as "still in the same
# basin" rather than having wandered onto a different Lambert branch.
VINF_DRIFT_FLAG_KMS = 0.5

# Boundary-pin tolerance (NM landing within this of a box edge counts as
# "pinned", i.e. formulation-conditional, not a genuine free-optimum).
BOUNDARY_PIN_TOL_DEG = 0.05
BOUNDARY_PIN_TOL_TOF = 0.002


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def load_population() -> list[dict[str, Any]]:
    """Load the 69 Titan-anchored, physical_gate_passed=True candidates."""
    recs: list[dict[str, Any]] = []
    with SCAN_571_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            row = json.loads(line)
            if row.get("physical_gate_passed") is True and row.get("anchor") == ANCHOR:
                rec = row["record"]
                recs.append(
                    {
                        "label": (
                            f"rel{rec['rel_offset_deg']:.0f}_tof{rec['tof_scale']:.2f}_"
                            f"n{rec['n_rev'][0]}{rec['n_rev'][1]}"
                        ),
                        "rel_offset_deg": float(rec["rel_offset_deg"]),
                        "tof_scale": float(rec["tof_scale"]),
                        "n_rev": tuple(rec["n_rev"]),
                        "coplanar_residual_kms": float(rec["residual_kms"]),
                        "coplanar_iapetus_vinf_kms": float(row["vinf_per_encounter_kms"][1]),
                        "coplanar_bend_deg": float(row["max_bend_deg_per_encounter"][1]),
                    }
                )
    if len(recs) != 69:
        print(
            f"[573] WARNING: expected 69 Titan-anchored gate-passers, found {len(recs)}",
            flush=True,
        )
    return recs


def evaluate_point_tracked(
    rel_offset_deg: float,
    tof_scale: float,
    n_rev: tuple[int, int],
    omega_deg: float,
    inc_deg: float,
) -> tuple[dict[str, Any] | None, str | None]:
    """Same math as #572's ``evaluate_point``, but with
    LambertGeometryError tracked SEPARATELY from LambertConvergenceError
    (item 5 of the Fable correction -- never exercised at scale in #572's
    2-candidate run). Returns ``(result_or_None, error_kind)`` where
    ``error_kind`` is one of ``None`` (success), ``"geometry"``,
    ``"convergence"``, or ``"infeasible_n_rev"`` (Lambert solved but not at
    the requested revolution count -- not an exception, just no solution).
    """
    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    sat_b = SATELLITES[FLYBY]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    v_circ_b = math.sqrt(mu / sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof = tof_scale * math.sqrt(p_a * p_b)
    tof_s = tof * DAY_S

    omega = math.radians(omega_deg)
    inc = math.radians(inc_deg)
    u0 = math.radians(rel_offset_deg)

    r0, v0 = _moon_state(0.0, n_a, 0.0, sat_a.sma_km, mu)
    r1, v1 = iapetus_state_3d(u0 + n_b * tof, v_circ_b, sat_b.sma_km, omega, inc)
    r2, v2 = _moon_state(0.0, n_a, 2.0 * tof, sat_a.sma_km, mu)

    n0, n1 = n_rev
    try:
        leg0 = _leg_best(r0, v0, r1, v1, tof_s, mu, n0)
        leg1 = _leg_best(r1, v1, r2, v2, tof_s, mu, n1)
    except LambertGeometryError:
        return None, "geometry"
    except LambertConvergenceError:
        return None, "convergence"
    if leg0 is None or leg1 is None:
        return None, "infeasible_n_rev"

    r_mid = abs(leg0["vinf_in"] - leg1["vinf_out"])
    r_periodic = abs(leg0["vinf_out"] - leg1["vinf_in"])
    residual = max(r_mid, r_periodic)

    vinf0 = leg0["vinf_out"]
    vinf1 = max(leg0["vinf_in"], leg1["vinf_out"])
    vinf2 = leg1["vinf_in"]

    return (
        {
            "rel_offset_deg": rel_offset_deg,
            "tof_scale": tof_scale,
            "n_rev": list(n_rev),
            "omega_deg": omega_deg,
            "inc_deg": inc_deg,
            "residual_kms": residual,
            "vinf_kms": [vinf0, vinf1, vinf2],
        },
        None,
    )


def _circ_dist_deg(a: float, b: float) -> float:
    d = abs((a - b + 180.0) % 360.0 - 180.0)
    return d


def sweep_node_alignment(
    cand: dict[str, Any], err_counts: dict[str, int], *, n_omega: int = N_OMEGA
) -> dict[str, Any]:
    """Fine grid-search of Omega at fixed (rel_offset, tof_scale, n_rev),
    enumerate ALL local-minimum basins, then refine each of the (capped)
    deepest 12 with an ADAPTIVE-window bounded Nelder-Mead. Mirrors #572's
    ``sweep_node_alignment`` with the Fable at-scale guards layered on.
    """
    from scipy.optimize import minimize

    rel_offset_deg = cand["rel_offset_deg"]
    tof_scale = cand["tof_scale"]
    n_rev = tuple(cand["n_rev"])

    grid: list[tuple[float, dict[str, Any] | None]] = []
    n_geom = 0
    n_conv = 0
    n_geom_retried_ok = 0
    for i in range(n_omega):
        omega_deg = 360.0 * i / n_omega
        pt, err = evaluate_point_tracked(
            rel_offset_deg, tof_scale, n_rev, omega_deg, INCLINATION_DEG
        )
        if pt is None:
            if err == "geometry":
                n_geom += 1
                err_counts["geometry_total"] += 1
                for eps in (0.02, -0.02, 0.05, -0.05, 0.15, -0.15):
                    pt, _err2 = evaluate_point_tracked(
                        rel_offset_deg, tof_scale, n_rev, omega_deg + eps, INCLINATION_DEG
                    )
                    if pt is not None:
                        n_geom_retried_ok += 1
                        break
            elif err == "convergence":
                n_conv += 1
                err_counts["convergence_total"] += 1
        grid.append((omega_deg, pt))

    n_feasible = sum(1 for _o, p in grid if p is not None)

    basins_grid: list[dict[str, Any]] = []
    for i in range(n_omega):
        _o, p = grid[i]
        if p is None:
            continue
        _op, pp = grid[i - 1]
        _on, pn = grid[(i + 1) % n_omega]
        if pp is None or pn is None:
            continue
        if p["residual_kms"] <= pp["residual_kms"] and p["residual_kms"] <= pn["residual_kms"]:
            basins_grid.append(p)
    basins_grid.sort(key=lambda r: r["residual_kms"])

    cap_bound = len(basins_grid) > BASIN_CAP
    basins_to_refine = basins_grid[:BASIN_CAP]

    all_basin_omegas = [b["omega_deg"] for b in basins_grid]

    refined_basins: list[dict[str, Any]] = []
    for seed in basins_to_refine:
        # Adaptive window: min(15 deg, half-distance to nearest OTHER
        # basin's Omega, periodic wraparound) -- item 2 of the Fable
        # correction, prevents branch-jumping when basins sit close.
        other_omegas = [o for o in all_basin_omegas if o != seed["omega_deg"]]
        if other_omegas:
            nearest_dist = min(_circ_dist_deg(seed["omega_deg"], o) for o in other_omegas)
            half_window = nearest_dist / 2.0
        else:
            half_window = 15.0
        window_deg = min(15.0, half_window)

        omega_lo, omega_hi = seed["omega_deg"] - window_deg, seed["omega_deg"] + window_deg
        tof_lo, tof_hi = tof_scale - 0.1, tof_scale + 0.1

        def _obj(x: np.ndarray, _n_rev: tuple[int, int] = n_rev) -> float:
            omega_deg_x, tof_scale_x = float(x[0]), float(x[1])
            pt, _err = evaluate_point_tracked(
                rel_offset_deg, tof_scale_x, _n_rev, omega_deg_x, INCLINATION_DEG
            )
            if pt is None:
                return 1.0e3
            return pt["residual_kms"]

        x0 = np.array([seed["omega_deg"], tof_scale])
        res = minimize(
            _obj,
            x0,
            method="Nelder-Mead",
            bounds=[(omega_lo, omega_hi), (tof_lo, tof_hi)],
            options={"xatol": 1e-5, "fatol": 1e-9, "maxiter": 200, "maxfev": 200},
        )
        refined_omega = float(res.x[0]) % 360.0
        refined_tof = float(res.x[1])
        refined_pt, _err = evaluate_point_tracked(
            rel_offset_deg, refined_tof, n_rev, refined_omega, INCLINATION_DEG
        )
        if refined_pt is None or refined_pt["residual_kms"] > seed["residual_kms"]:
            refined_pt = dict(seed)
            refined_omega = seed["omega_deg"]
            refined_tof = tof_scale

        # Boundary-pin check (item 4): landed within tolerance of a box edge.
        boundary_pinned = (
            abs(float(res.x[0]) - omega_lo) < BOUNDARY_PIN_TOL_DEG
            or abs(float(res.x[0]) - omega_hi) < BOUNDARY_PIN_TOL_DEG
            or abs(refined_tof - tof_lo) < BOUNDARY_PIN_TOL_TOF
            or abs(refined_tof - tof_hi) < BOUNDARY_PIN_TOL_TOF
        )

        # Basin-support drift check (item 2): refined V_inf(mid) should stay
        # close to the SEED grid point's V_inf(mid); large drift means NM
        # wandered onto a different Lambert branch despite the bounded Omega
        # window (tof_scale is only loosely bounded at +-0.1).
        vinf_drift = abs(refined_pt["vinf_kms"][1] - seed["vinf_kms"][1])
        vinf_drift_flag = vinf_drift > VINF_DRIFT_FLAG_KMS

        refined_pt = dict(refined_pt)
        refined_pt["seed_omega_deg"] = seed["omega_deg"]
        refined_pt["seed_vinf_kms"] = seed["vinf_kms"]
        refined_pt["nm_window_deg"] = window_deg
        refined_pt["boundary_pinned"] = boundary_pinned
        refined_pt["vinf_drift_flag"] = vinf_drift_flag
        refined_pt["vinf_drift_kms"] = vinf_drift
        refined_basins.append(refined_pt)

    refined_basins.sort(key=lambda r: r["residual_kms"])

    return {
        "label": cand["label"],
        "rel_offset_deg": rel_offset_deg,
        "tof_scale_seed": tof_scale,
        "n_rev": list(n_rev),
        "coplanar_residual_kms": cand["coplanar_residual_kms"],
        "coplanar_iapetus_vinf_kms": cand["coplanar_iapetus_vinf_kms"],
        "coplanar_bend_deg": cand["coplanar_bend_deg"],
        "inclination_deg": INCLINATION_DEG,
        "n_omega_grid": n_omega,
        "n_geometry_errors_encountered": n_geom,
        "n_geometry_errors_resolved_by_retry": n_geom_retried_ok,
        "n_convergence_errors_encountered": n_conv,
        "n_feasible_omega_points": n_feasible,
        "n_basins_found": len(basins_grid),
        "basin_cap_bound": cap_bound,
        "basins": refined_basins,
    }


def evaluate_basins_against_gates(
    sweep: dict[str, Any], seq: tuple[str, str, str]
) -> list[dict[str, Any]]:
    basin_evals: list[dict[str, Any]] = []
    n_rev = sweep["n_rev"]
    for b in sweep["basins"]:
        residual_near = b["residual_kms"] < GATE_RESIDUAL_KMS
        gate_pass, gate_verdicts = candidate_passes_physical_gate(
            seq, tuple(b["vinf_kms"]), min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
        )
        bends = [v.max_bend_deg for v in gate_verdicts]
        formulation_conditional = bool(b["boundary_pinned"])
        genuine_closure = bool(residual_near and gate_pass and not formulation_conditional)
        basin_evals.append(
            {
                "n_rev": n_rev,
                "omega_deg": b["omega_deg"],
                "tof_scale": b["tof_scale"],
                "residual_kms": b["residual_kms"],
                "vinf_kms": b["vinf_kms"],
                "residual_near_coplanar": residual_near,
                "physical_gate_pass": gate_pass,
                "max_bend_deg_per_encounter": bends,
                "boundary_pinned": b["boundary_pinned"],
                "vinf_drift_flag": b["vinf_drift_flag"],
                "vinf_drift_kms": b["vinf_drift_kms"],
                "nm_window_deg": b["nm_window_deg"],
                "formulation_conditional_non_closure": formulation_conditional
                and residual_near
                and gate_pass,
                "closure": genuine_closure,
            }
        )
    return basin_evals


def cluster_branches(all_closing: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Union-find dedup of closing basins across ALL 69 candidate runs by
    proximity in (n_rev, Omega*, tof_scale*, V_inf_mid) -- item 1 of the
    Fable correction. Two closing-basin instances merge iff same n_rev AND
    within MERGE_OMEGA_DEG (periodic) AND MERGE_TOF_SCALE AND
    MERGE_VINF_KMS of each other.
    """
    n = len(all_closing)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            a, b = all_closing[i], all_closing[j]
            if tuple(a["n_rev"]) != tuple(b["n_rev"]):
                continue
            if _circ_dist_deg(a["omega_deg"], b["omega_deg"]) > MERGE_OMEGA_DEG:
                continue
            if abs(a["tof_scale"] - b["tof_scale"]) > MERGE_TOF_SCALE:
                continue
            if abs(a["vinf_kms"][1] - b["vinf_kms"][1]) > MERGE_VINF_KMS:
                continue
            union(i, j)

    groups: dict[int, list[dict[str, Any]]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(all_closing[i])
    return list(groups.values())


def main() -> int:
    t_start = time.time()
    sha = _git_sha()
    print(f"[573] Titan-Iapetus 3D-closure POPULATION probe -- sha={sha}", flush=True)

    print("[573] smoke test: iapetus_state_3d reduces to _moon_state at inc=0 ...", flush=True)
    smoke_ok = _smoke_test_reduction()
    print(f"[573] smoke test PASS: {smoke_ok}", flush=True)
    if not smoke_ok:
        print("[573] ABORTING -- 3D state generator does not reduce correctly.", flush=True)
        return 1

    population = load_population()
    print(f"[573] loaded {len(population)} Titan-anchored gate-passing candidates", flush=True)

    seq = (ANCHOR, FLYBY, ANCHOR)
    err_counts = {"geometry_total": 0, "convergence_total": 0}

    out_fh = OUT_PATH.open("w", encoding="utf-8")

    def _write(rec: dict[str, Any]) -> None:
        out_fh.write(json.dumps(rec, default=str) + "\n")
        out_fh.flush()

    _write(
        {
            "_meta": True,
            "task": "#573 widened Titan-Iapetus 3D-closure population probe",
            "git_sha": sha,
            "inclination_deg": INCLINATION_DEG,
            "gate_residual_kms": GATE_RESIDUAL_KMS,
            "min_useful_bend_deg": DEFAULT_MIN_USEFUL_BEND_DEG,
            "ecc_robust_bend_deg": ECC_ROBUST_BEND_DEG,
            "n_omega_grid": N_OMEGA,
            "basin_cap": BASIN_CAP,
            "merge_omega_deg": MERGE_OMEGA_DEG,
            "merge_tof_scale": MERGE_TOF_SCALE,
            "merge_vinf_kms": MERGE_VINF_KMS,
            "vinf_drift_flag_kms": VINF_DRIFT_FLAG_KMS,
            "smoke_test_reduction_pass": smoke_ok,
            "n_population": len(population),
        }
    )

    candidate_verdicts: list[dict[str, Any]] = []
    all_closing: list[dict[str, Any]] = []
    all_formulation_conditional: list[dict[str, Any]] = []
    cap_bound_candidates: list[str] = []
    n_basin_cap_binds = 0

    for idx, cand in enumerate(population):
        t0 = time.time()
        sweep = sweep_node_alignment(cand, err_counts, n_omega=N_OMEGA)
        elapsed = time.time() - t0
        sweep["elapsed_s"] = elapsed
        if sweep["basin_cap_bound"]:
            n_basin_cap_binds += 1
            cap_bound_candidates.append(cand["label"])

        basin_evals = evaluate_basins_against_gates(sweep, seq)
        closing = [be for be in basin_evals if be["closure"]]
        formulation_cond = [be for be in basin_evals if be["formulation_conditional_non_closure"]]

        for be in closing:
            be2 = dict(be)
            be2["seed_label"] = cand["label"]
            be2["seed_rel_offset_deg"] = cand["rel_offset_deg"]
            be2["seed_tof_scale"] = cand["tof_scale"]
            all_closing.append(be2)
        for be in formulation_cond:
            be2 = dict(be)
            be2["seed_label"] = cand["label"]
            all_formulation_conditional.append(be2)

        verdict = {
            "label": cand["label"],
            "n_rev": list(cand["n_rev"]),
            "rel_offset_deg": cand["rel_offset_deg"],
            "tof_scale": cand["tof_scale"],
            "closure_found": len(closing) > 0,
            "n_basins_found": sweep["n_basins_found"],
            "basin_cap_bound": sweep["basin_cap_bound"],
            "n_closing_basins": len(closing),
            "n_formulation_conditional_non_closures": len(formulation_cond),
            "n_geometry_errors": sweep["n_geometry_errors_encountered"],
            "n_geometry_errors_resolved": sweep["n_geometry_errors_resolved_by_retry"],
            "n_convergence_errors": sweep["n_convergence_errors_encountered"],
            "elapsed_s": elapsed,
        }
        candidate_verdicts.append(verdict)
        _write({"kind": "candidate_result", **sweep})
        _write({"kind": "candidate_verdict", **verdict})

        print(
            f"[573] [{idx + 1:2d}/{len(population)}] {cand['label']:24s} "
            f"basins={sweep['n_basins_found']:2d}{'(CAP)' if sweep['basin_cap_bound'] else '':5s} "
            f"closing={len(closing)} formcond={len(formulation_cond)} "
            f"geom_err={sweep['n_geometry_errors_encountered']} "
            f"conv_err={sweep['n_convergence_errors_encountered']} "
            f"({elapsed:.2f}s, total {time.time() - t_start:.1f}s)",
            flush=True,
        )

    # --- Item 6: singleton-flip guard, per contiguous (n_rev, tof_scale)
    # rel_offset run (only genuinely contiguous windows, per
    # detect_isolated_singleton_anomalies' own caution). ---
    from collections import defaultdict

    groups: dict[tuple[tuple[int, int], float], list[dict[str, Any]]] = defaultdict(list)
    for v in candidate_verdicts:
        key = (tuple(v["n_rev"]), v["tof_scale"])
        groups[key].append(v)

    singleton_flags: list[dict[str, Any]] = []
    for key, members in groups.items():
        if len(members) < 3:
            continue
        members_sorted = sorted(members, key=lambda m: m["rel_offset_deg"])
        vals = [m["closure_found"] for m in members_sorted]
        labels = [m["label"] for m in members_sorted]
        anomalies = detect_isolated_singleton_anomalies(vals, labels)
        for a in anomalies:
            singleton_flags.append(
                {
                    "n_rev": list(key[0]),
                    "tof_scale": key[1],
                    "label": a.label,
                    "closure_value": a.value,
                }
            )

    # --- Item 1: branch dedup via union-find on closing basins. ---
    branch_groups = cluster_branches(all_closing)
    branch_groups.sort(key=lambda g: min(b["residual_kms"] for b in g))

    branches_summary: list[dict[str, Any]] = []
    for gi, g in enumerate(branch_groups):
        rep = min(g, key=lambda b: b["residual_kms"])
        seed_labels = sorted({b["seed_label"] for b in g})
        n_revs_in_branch = sorted({tuple(b["n_rev"]) for b in g})
        iapetus_bend = rep["max_bend_deg_per_encounter"][1]
        any_drift_flag = any(b["vinf_drift_flag"] for b in g)
        branches_summary.append(
            {
                "branch_id": gi,
                "n_member_closures": len(g),
                "seed_labels": seed_labels,
                "n_rev": [list(nr) for nr in n_revs_in_branch],
                "representative_omega_deg": rep["omega_deg"],
                "representative_tof_scale": rep["tof_scale"],
                "representative_residual_kms": rep["residual_kms"],
                "representative_vinf_kms": rep["vinf_kms"],
                "representative_iapetus_bend_deg": iapetus_bend,
                "eccentricity_robust": bool(iapetus_bend >= ECC_ROBUST_BEND_DEG),
                "any_vinf_drift_flag_in_branch": any_drift_flag,
            }
        )

    n_distinct_branches = len(branch_groups)
    n_rev_classes_spanned = len(
        {tuple(nr) for g in branch_groups for nr in {tuple(b["n_rev"]) for b in g}}
    )
    n_ecc_robust = sum(1 for b in branches_summary if b["eccentricity_robust"])

    if n_distinct_branches >= 5 and n_rev_classes_spanned >= 2 and n_ecc_robust >= 3:
        bucket = "REAL_FAMILY_SIGNAL"
        phase2_mandatory = False
    elif n_distinct_branches <= 2:
        bucket = "ADDS_NOTHING_BEYOND_572"
        phase2_mandatory = False
    else:
        bucket = "MARGINAL"
        phase2_mandatory = True

    summary = {
        "kind": "population_summary",
        "n_population_candidates": len(population),
        "n_candidates_with_any_closure": sum(1 for v in candidate_verdicts if v["closure_found"]),
        "n_raw_closing_basin_instances": len(all_closing),
        "n_distinct_branches": n_distinct_branches,
        "n_rev_classes_spanned_by_branches": n_rev_classes_spanned,
        "n_eccentricity_robust_branches": n_ecc_robust,
        "ecc_robust_bend_threshold_deg": ECC_ROBUST_BEND_DEG,
        "decision_bucket": bucket,
        "phase2_iapetus_anchored_mandatory": phase2_mandatory,
        "n_basin_cap_binds": n_basin_cap_binds,
        "basin_cap_bound_candidates": cap_bound_candidates,
        "n_formulation_conditional_non_closures": len(all_formulation_conditional),
        "n_geometry_errors_total": err_counts["geometry_total"],
        "n_convergence_errors_total": err_counts["convergence_total"],
        "n_singleton_flip_flags": len(singleton_flags),
        "singleton_flip_flags": singleton_flags,
        "branches": branches_summary,
        "total_elapsed_s": time.time() - t_start,
    }
    _write(summary)
    out_fh.close()

    print("[573] === SUMMARY ===", flush=True)
    print(f"[573]   candidates run: {len(population)}", flush=True)
    print(
        f"[573]   candidates with >=1 closing basin: {summary['n_candidates_with_any_closure']}",
        flush=True,
    )
    print(f"[573]   raw closing-basin instances: {len(all_closing)}", flush=True)
    print(f"[573]   DISTINCT BRANCHES (deduped): {n_distinct_branches}", flush=True)
    print(f"[573]   n_rev classes spanned by branches: {n_rev_classes_spanned}", flush=True)
    print(
        f"[573]   eccentricity-robust branches (bend>={ECC_ROBUST_BEND_DEG}): {n_ecc_robust}",
        flush=True,
    )
    print(f"[573]   DECISION BUCKET: {bucket}", flush=True)
    print(f"[573]   phase2 (Iapetus-anchored) mandatory: {phase2_mandatory}", flush=True)
    print(f"[573]   basin cap binds: {n_basin_cap_binds} ({cap_bound_candidates})", flush=True)
    print(
        f"[573]   formulation-conditional non-closures: {len(all_formulation_conditional)}",
        flush=True,
    )
    print(f"[573]   geometry errors total: {err_counts['geometry_total']}", flush=True)
    print(f"[573]   convergence errors total: {err_counts['convergence_total']}", flush=True)
    print(f"[573]   singleton flip flags: {len(singleton_flags)}", flush=True)
    for b in branches_summary:
        print(
            f"[573]     branch {b['branch_id']}: n_rev={b['n_rev']} "
            f"members={b['n_member_closures']} "
            f"Omega*={b['representative_omega_deg']:.2f} tof*={b['representative_tof_scale']:.4f} "
            f"vinf={[f'{v:.3f}' for v in b['representative_vinf_kms']]} "
            f"iapetus_bend={b['representative_iapetus_bend_deg']:.2f} "
            f"ecc_robust={b['eccentricity_robust']}",
            flush=True,
        )
    print(f"[573] results written to {OUT_PATH}", flush=True)
    print(f"[573] total wall clock: {time.time() - t_start:.1f}s", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
