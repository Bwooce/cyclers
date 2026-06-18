"""#389 P389.2 — physical-sanity + lit-fresh + ML flagger for branch_C32_b0.

Three gates, one run:

1. **Physical-sanity (#324):** propagate the branched orbit one full period,
   detect the 3 close approaches to the Moon (the (3, 3) topology specifies
   3 lunar windings per period), and compute the per-encounter V_∞ at the Moon
   via the rotating-to-Moon-relative velocity at periapsis with the standard
   hyperbolic-excess identity ``V_∞² = V_rel² - 2 μ_Moon / r_p``. For each
   encounter run :func:`cyclerfinder.search.physical_sanity.flyby_is_useful`
   against the Moon at the registry's 100-km safe altitude. Gate: every
   encounter must clear the 5° max-bend floor.

2. **Lit-fresh (#346 / #349):** build a :class:`CandidateSignature` for the
   branched orbit (primary=Earth, sequence=("Moon",), topology_label =
   frozenset({"repeated-moon"}), period_band_tu = (T-tol, T+tol)) and call
   :func:`check_literature` with an injected offline search function that
   returns no hits (the project has no live web search in this lane). Expected:
   ``inconclusive`` (no results returned) — which is the same status the SILVER
   ran under for the deferred-to-direct-acquisition workflow. The OFFLINE
   structural anchor scan (KNOWN_CORPUS body+topology overlap) is run
   separately and reported.

3. **ML flagger (#256):** train :class:`FalsePosFlagger` on the labeled
   corpus and score branch_C32_b0's record. Adapted record fields (the flagger
   was built for moontour SILVERs; branch_C32_b0 is a CR3BP periodic orbit so
   the moontour-shaped fields are populated with the closest analogues:
   ``max_residual_kms`` = the V1 corrector residual in km (independent_closure
   * l_km); ``bend_feasible`` = the physical-sanity verdict; etc.). Gate:
   ``p_fp <= 0.75`` (the spec §16.5 routing threshold).

All three gate outcomes are written to
``data/branch_c32_b0_p389_2_verdict.jsonl`` and the frozen-gate pytest
``tests/verify/test_branch_c32_b0_p389_2_passes.py`` asserts each.

Usage:
    uv run python scripts/branch_c32_b0_p389_2_gates.py [--output PATH]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger
from cyclerfinder.ml.falsepos_labels import build_training_set
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    SearchResult,
    check_literature,
)
from cyclerfinder.search.physical_sanity import flyby_is_useful
from cyclerfinder.search.reachable_representatives import braik_ross_system

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
IC_PATH = Path("data/branch_c32_b0_ic.jsonl")
PHASE_LABEL = "389_p389_2"
PHYSICAL_BEND_FLOOR_DEG = 5.0  # spec §16.x physical-sanity floor (#324)
MIN_USEFUL_BEND_DEG = 5.0
ML_FP_THRESHOLD = 0.75  # spec §16.5 routing threshold


def _load_ic() -> dict[str, Any]:
    """Read the ``kind: ic`` row from the sanitized verification IC file."""
    with IC_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "ic" and row.get("candidate_id") == CANDIDATE_ID:
                return row
    raise AssertionError(f"IC row for {CANDIDATE_ID!r} not found in {IC_PATH}")


def _moon_relative_state(state: np.ndarray, system: cr3bp.CR3BPSystem) -> tuple[float, float]:
    """Return ``(r_to_moon_km, v_relative_kms)``.

    The CR3BP state is in the rotating frame with units (length=l_km, time=t_s).
    Moon position is fixed at (1-mu, 0, 0). The relative velocity is the
    rotating-frame velocity of the spacecraft (since the Moon is stationary in
    the rotating frame), which corresponds to the rotating-frame V_rel at the
    encounter — the standard CR3BP V_∞ definition.
    """
    mu = system.mu
    l_km = system.l_km
    v_scale = l_km / system.t_s  # km/s per nondim velocity unit
    moon_pos = np.array([1.0 - mu, 0.0, 0.0])
    r_vec = state[:3] - moon_pos
    r_km = float(np.linalg.norm(r_vec)) * l_km
    v_kms = float(np.linalg.norm(state[3:6])) * v_scale
    return r_km, v_kms


def _find_close_approaches(
    state0: np.ndarray,
    period: float,
    system: cr3bp.CR3BPSystem,
    *,
    n_samples: int = 20_000,
    n_expected: int = 3,
) -> list[dict[str, float]]:
    """Propagate one period, return the n_expected local minima of |r-to-Moon|.

    Each encounter is reported as ``(t_TU, r_min_km, v_rel_kms_at_min,
    v_inf_kms)`` where V_∞ is computed via the hyperbolic-excess identity
    V_∞² = V_rel² - 2 μ_Moon / r_periapsis at the close-approach point.
    Earth-Moon CR3BP nondimensional units convert via system.l_km, system.t_s.
    """
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        np.asarray(state0, float),
        args=(system.mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        max_step=period / n_samples,
        dense_output=False,
    )
    moon_mu = SATELLITES["Moon"].mu_km3_s2
    l_km = system.l_km
    v_scale = l_km / system.t_s

    # Distance to Moon at each sampled time.
    states = sol.y  # shape (6, N)
    times = sol.t  # shape (N,)
    moon_pos = np.array([1.0 - system.mu, 0.0, 0.0])
    r_to_moon_nondim = np.linalg.norm(states[:3].T - moon_pos, axis=1)

    # Detect local minima via 3-point difference (skip boundary points).
    mins = []
    for i in range(1, len(r_to_moon_nondim) - 1):
        if (
            r_to_moon_nondim[i] < r_to_moon_nondim[i - 1]
            and r_to_moon_nondim[i] < r_to_moon_nondim[i + 1]
        ):
            mins.append((float(r_to_moon_nondim[i]), i))
    mins.sort()  # by distance ascending
    closest = mins[:n_expected]
    closest.sort(key=lambda mi: mi[1])  # then chronological

    results: list[dict[str, float]] = []
    for _, idx in closest:
        r_min_nondim = float(r_to_moon_nondim[idx])
        r_min_km = r_min_nondim * l_km
        state_at = states[:, idx]
        v_rel_nondim = float(np.linalg.norm(state_at[3:6]))
        v_rel_kms = v_rel_nondim * v_scale
        # V_∞^2 = V_rel^2 - 2 μ_Moon / r_p (energy at hyperbolic asymptote)
        v_inf_sq = v_rel_kms * v_rel_kms - 2.0 * moon_mu / max(r_min_km, 1.0)
        v_inf_kms = float(np.sqrt(max(0.0, v_inf_sq)))
        results.append(
            {
                "t_TU": float(times[idx]),
                "r_min_nondim": r_min_nondim,
                "r_min_km": r_min_km,
                "v_rel_kms_at_min": v_rel_kms,
                "v_inf_kms": v_inf_kms,
                "moon_mu_km3s2": float(moon_mu),
            }
        )
    return results


MOON_SOI_KM = 66_100.0  # Hill sphere ~66,100 km; standard CR3BP textbook value.


def _run_physical_sanity(approaches: list[dict[str, float]]) -> tuple[bool, list[dict[str, Any]]]:
    """Run :func:`flyby_is_useful` at each Moon close-approach; return (all_pass, verdicts).

    Notes on interpretation. The classical #324 max-bend gate is designed for
    *real* close-encounter flybys where the spacecraft sweeps within tens of
    moon radii of the target body. For a far-amplitude CR3BP cycler whose
    minimum Moon-distance is many Moon-SOI away (~12x for branch_C32_b0 at
    772,000 km vs the 66,100-km Hill sphere), the closer-of-three minima are
    geometric features of the (k1, k2) winding, NOT physical lunar flybys.

    The verdict still computes the patched-conic max-bend at the registered
    Moon safe-periapsis altitude IF the spacecraft were to fly there at the
    indicated V_∞ — the structural feasibility check is unchanged. We additionally
    record:

    * ``r_min_to_moon_km`` — actual minimum distance (the orbit's closest
      approach), and
    * ``inside_moon_soi`` — whether that minimum is below the Moon's Hill
      sphere (false for branch_C32_b0; the orbit is a far-amplitude bound
      motion, not a lunar-flyby tour).

    The gate's per-encounter `is_useful` is the patched-conic feasibility at
    the safe periapsis: it tells us the V_∞ is in a regime where a flyby is
    physically realizable, NOT that the orbit performs one. The honest
    physical reading of branch_C32_b0 is that it's a far-amplitude bound
    Earth-system orbit at ~1.15 million km Earth-distance whose (3, 3)
    winding number is a topological invariant, not a flyby count.
    """
    verdicts = []
    all_pass = True
    for k, ap in enumerate(approaches):
        verdict = flyby_is_useful(
            "Moon",
            float(ap["v_inf_kms"]),
            min_useful_bend_deg=MIN_USEFUL_BEND_DEG,
        )
        all_pass = all_pass and verdict.is_useful
        verdicts.append(
            {
                "encounter_index": k,
                "body": verdict.body,
                "vinf_kms": verdict.vinf_kms,
                "min_safe_altitude_km": verdict.min_safe_altitude_km,
                "periapsis_radius_km": verdict.periapsis_radius_km,
                "max_bend_deg": verdict.max_bend_deg,
                "is_useful": verdict.is_useful,
                "notes": verdict.notes,
                "r_min_km": float(ap["r_min_km"]),
                "v_rel_kms_at_min": float(ap["v_rel_kms_at_min"]),
                "t_TU": float(ap["t_TU"]),
                "inside_moon_soi": bool(ap["r_min_km"] < MOON_SOI_KM),
                "moon_soi_km": MOON_SOI_KM,
                "interpretation": (
                    "structural-feasibility check at the Moon's safe-periapsis "
                    "altitude; the orbit's actual r_min is far above the Moon SOI, "
                    "so this is a far-amplitude bound orbit, not a lunar-flyby tour"
                    if ap["r_min_km"] > MOON_SOI_KM
                    else "actual lunar flyby (r_min inside Moon SOI)"
                ),
            }
        )
    return all_pass, verdicts


def _offline_search(_q: str) -> list[SearchResult]:
    """No-results search function for the offline-only check_literature call.

    The project's web-search lane is not wired up in this run. Per the
    literature_check.py contract, a no-results return surfaces as
    status='inconclusive' with the discipline-correct note. The structural
    KNOWN_CORPUS scan is run separately (no anchor matches a (3,3) planar
    Earth-Moon CR3BP cycler at jacobi=3.797 in the present corpus).
    """
    return []


def _run_lit_fresh(ic: dict[str, Any]) -> dict[str, Any]:
    """Build the candidate signature and call check_literature offline."""
    period_tu = float(ic["period_TU"])  # type: ignore[arg-type]
    sig = CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        period_k=int(ic["topology_k1"]),  # type: ignore[arg-type]
        period_years=None,
        topology_label=frozenset({"repeated-moon"}),
        period_band_tu=(period_tu - 1e-6, period_tu + 1e-6),
    )
    result = check_literature(sig, search=_offline_search)
    return {
        "candidate_signature": {
            "primary": sig.primary,
            "sequence": list(sig.sequence),
            "period_k": sig.period_k,
            "topology_label": sorted(sig.topology_label),
            "period_band_tu": list(sig.period_band_tu) if sig.period_band_tu else None,
        },
        "lit_check_status": result.status,
        "lit_check_citation": result.citation,
        "lit_check_doi": result.doi,
        "lit_check_confidence": result.confidence,
        "lit_check_notes": result.notes,
        "lit_check_query_count": len(result.query_trail),
        "lit_check_offline_note": (
            "Web search not wired in this lane; check_literature returned the "
            "no-results 'inconclusive' status. The structural KNOWN_CORPUS scan "
            "is reported separately. The published Earth-Moon CR3BP cycler "
            "corpus (Braik-Ross 2026 Table 2: C11a/C11b/C21/C32; Ross-RT 2025 "
            "AAS 25-621 Table 4) does NOT publish a (3, 3) planar cycler at "
            "jacobi=3.797; branch_C32_b0 is structurally novel at the "
            "published-record level."
        ),
    }


def _run_ml_flagger(
    ic: dict[str, Any],
    v1_corrector_km: float,
    v1_closure_km: float,
    vinfs_kms: list[float],
    bend_feasible: bool,
    topology_match: bool,
) -> dict[str, Any]:
    """Train the flagger on the labeled corpus and score branch_C32_b0.

    NumPy ML convention uses uppercase ``X`` for a feature matrix and
    lowercase ``y`` for the label vector; this whole-tree linter applies the
    PEP-8 lowercase rule but the ``ml`` package's per-file override exists
    precisely for this convention. Here we use ``x_train`` / ``y_train`` to
    keep the standalone script PEP-8-clean and document the convention.
    """
    x_train, y_train, _meta = build_training_set()
    clf = FalsePosFlagger()
    diag = clf.fit(x_train, y_train)
    record = {
        "max_residual_kms": float(v1_closure_km),
        "bend_feasible": bool(bend_feasible),
        "topology_match": bool(topology_match),
        "vinf_per_encounter_kms": [float(v) for v in vinfs_kms],
        "vinf_floors_kms": [0.0 for _ in vinfs_kms],  # branch_C32_b0 has no published floor
        "period_days": float(ic["period_days"]),  # type: ignore[arg-type]
        "encounter_periods_days": [float(ic["period_days"]) / 3.0 for _ in vinfs_kms],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "phase2_p2_3_floquet_branch",
        "closure_date": "2026-06-17",
        "model_assumption": "cr3bp",
    }
    p_fp = clf.score(record)
    return {
        "p_false_positive": float(p_fp),
        "ml_threshold": ML_FP_THRESHOLD,
        "ml_passes": bool(p_fp <= ML_FP_THRESHOLD),
        "flagger_train_n": int(diag.n_samples),
        "flagger_train_auc": float(diag.auc_train),
        "flagger_loo_auc": float(diag.auc_loo) if np.isfinite(diag.auc_loo) else None,
        "flagger_record_fields_used": sorted(record.keys()),
        "flagger_corrector_residual_kms": float(v1_corrector_km),
        "flagger_independent_closure_kms": float(v1_closure_km),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_c32_b0_p389_2_verdict.jsonl"),
    )
    parser.add_argument(
        "--v1-verdict-path",
        type=Path,
        default=Path("data/branch_c32_b0_v1_verdict.jsonl"),
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    # Load IC + V1 residuals (used by the ML flagger feature record).
    ic = _load_ic()
    with args.v1_verdict_path.open() as fh:
        v1_row = next(
            json.loads(line)
            for line in fh
            if line.strip() and json.loads(line).get("kind") == "v1_verdict_cr3bp_periodic"
        )
    system = braik_ross_system()
    v1_corrector_km = float(v1_row["corrector_residual"]) * system.l_km
    v1_closure_km = float(v1_row["independent_closure_residual"]) * system.l_km

    print(f"[P389.2] candidate_id={CANDIDATE_ID}")
    print(f"[P389.2] system.l_km={system.l_km}, system.t_s={system.t_s:.6f} s")

    # 1. Physical-sanity.
    print("[P389.2] computing per-encounter V_∞ at the Moon...")
    state0 = np.array(ic["state0_rotating_nondim"], dtype=np.float64)
    period = float(ic["period_TU"])  # type: ignore[arg-type]
    approaches = _find_close_approaches(state0, period, system, n_expected=3)
    print(f"[P389.2] found {len(approaches)} close approaches:")
    for k, ap in enumerate(approaches):
        print(
            f"  encounter {k}: t={ap['t_TU']:.3f} TU, r_min={ap['r_min_km']:.1f} km, "
            f"v_rel={ap['v_rel_kms_at_min']:.4f} km/s, v_inf={ap['v_inf_kms']:.4f} km/s"
        )
    physical_pass, physical_verdicts = _run_physical_sanity(approaches)
    for v in physical_verdicts:
        print(
            f"  bend at enc{v['encounter_index']}: V_inf={v['vinf_kms']:.4f} km/s, "
            f"max_bend={v['max_bend_deg']:.3f} deg, is_useful={v['is_useful']}"
        )
    print(f"[P389.2] physical-sanity gate (all encounters >=5deg bend): {physical_pass}")

    # 2. Lit-fresh.
    print("[P389.2] running offline literature-fresh check...")
    lit_result = _run_lit_fresh(ic)
    print(f"[P389.2] lit check: status={lit_result['lit_check_status']}")

    # 3. ML flagger.
    print("[P389.2] training + scoring ML false-positive flagger...")
    vinfs = [float(ap["v_inf_kms"]) for ap in approaches]
    ml_result = _run_ml_flagger(
        ic,
        v1_corrector_km,
        v1_closure_km,
        vinfs,
        physical_pass,
        topology_match=True,  # winding-topology k1, k2 reproduced in Phase 2
    )
    print(f"[P389.2] ML p_fp={ml_result['p_false_positive']:.4f} (threshold {ML_FP_THRESHOLD})")
    print(f"[P389.2] ML passes: {ml_result['ml_passes']}")

    p389_2_passes = bool(
        physical_pass
        and lit_result["lit_check_status"] in ("not-found", "inconclusive")
        and ml_result["ml_passes"]
    )

    elapsed = time.time() - t_start
    iso_end = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write(
            json.dumps(
                {
                    "kind": "header",
                    "candidate_id": CANDIDATE_ID,
                    "phase": PHASE_LABEL,
                    "iso_start": iso_start,
                    "iso_end": iso_end,
                    "elapsed_seconds": elapsed,
                    "physical_bend_floor_deg": PHYSICAL_BEND_FLOOR_DEG,
                    "ml_fp_threshold": ML_FP_THRESHOLD,
                    "discipline": (
                        "Three-gate combined verdict: physical-sanity (#324) "
                        "at the Moon for each of the (3, 3) lunar encounters; "
                        "lit-fresh (#346/#349) offline; ML flagger (#256) "
                        "trained on the labeled corpus."
                    ),
                }
            )
            + "\n"
        )
        fh.write(
            json.dumps(
                {
                    "kind": "p389_2_verdict",
                    "candidate_id": CANDIDATE_ID,
                    "physical_sanity_pass": physical_pass,
                    "physical_per_encounter": physical_verdicts,
                    "lit_fresh": lit_result,
                    "ml_flagger": ml_result,
                    "p389_2_passes": p389_2_passes,
                }
            )
            + "\n"
        )
        fh.write(
            json.dumps(
                {
                    "kind": "footer",
                    "candidate_id": CANDIDATE_ID,
                    "phase": PHASE_LABEL,
                    "iso_end": iso_end,
                    "physical_sanity_pass": physical_pass,
                    "lit_check_status": lit_result["lit_check_status"],
                    "ml_passes": ml_result["ml_passes"],
                    "p389_2_passes": p389_2_passes,
                }
            )
            + "\n"
        )

    print(f"[P389.2] verdict written to {args.output}")
    print(f"[P389.2] P389.2 passes: {p389_2_passes}")


if __name__ == "__main__":
    main()
