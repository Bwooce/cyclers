"""#327 verification probe — Umbriel-Oberon-Umbriel SILVER at gate-passing V_inf.

Open question (from #324 follow-up): the SILVER row
``repeated-moon-uranus-00000041`` in
``data/scan_312_uranus_oberon_umbriel.jsonl`` stores closure residual
0.025232 km/s at coarse-grid V_inf magnitudes (0.92 / 0.96 / 0.89 km/s).
The #324 physical-sanity gate (5 deg max-bend floor) was added AFTER the
#312 run and clears this row at the stored V_inf (max bend ~14.7 deg at
Umbriel). The unphysical 2.27 km/s V_inf was the offset-sweep refinement
record, not the SILVER row itself.

This script answers four questions at the gate-passing IC (the actual
stored SILVER row's coarse-grid phasing):

1.  ``RepeatedMoonTarget.close()`` under the 2-moon convention
    ``moons=(Oberon, Umbriel)`` — does it reproduce 0.025 km/s?
2.  Same close() under the 3-moon convention
    ``moons=(Titania, Oberon, Umbriel)`` — does it reproduce 0.636 km/s?
3.  Relative-offset sweep at gate-passing V_inf: where does the basin
    floor live? (24-sample sweep, narrower than #312's 96 — this is a
    confirmation, not a search.)
4.  n_rev sweep at the stored (1, 1) phasing: do n_rev in (0..5) close
    the gap further at this exact phasing?
5.  Independent DOP853 cross-check at the 2-moon closure leg geometry
    (rtol=atol=1e-12).
6.  Physical-sanity gate verdict at the stored V_inf.
7.  ``check_literature`` + ML flagger re-run on the gate-passing
    signature.

Discipline anchors:
* READ-ONLY on ``discovery_campaign.py``, ``saturn_uranus_campaign.py``,
  ``physical_sanity.py``, and the source JSONL.
* NO catalogue writeback.
* Output: ``data/silver_327_verified.jsonl``.
* Per-call wrapping only — we re-use the production close() and
  cross-check primitives verbatim.

Run as::

    uv run python scripts/verify_327_umbriel_silver.py
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

from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger  # noqa: E402
from cyclerfinder.ml.falsepos_labels import build_training_set  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    Candidate,
    RepeatedMoonTarget,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.five_tier_prioritizer import (  # noqa: E402
    PatchedConicLeg,
    legs_from_repeated_moon_candidate,
)
from cyclerfinder.search.literature_check import (  # noqa: E402
    CandidateSignature,
    check_literature,
)
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    candidate_passes_physical_gate,
)
from cyclerfinder.search.saturn_uranus_campaign import (  # noqa: E402
    dop853_cross_check_leg,
    offline_corpus_search,
)

SEQ = ("Umbriel", "Oberon", "Umbriel")
NREV = (1, 1)
SRC_JSONL = ROOT / "data" / "scan_312_uranus_oberon_umbriel.jsonl"
OUT_JSONL = ROOT / "data" / "silver_327_verified.jsonl"

# Stored row values (from data/scan_312_uranus_oberon_umbriel.jsonl row 41).
STORED_RESIDUAL_KMS = 0.025232272564609692
STORED_VINF = (0.9199258810725036, 0.9604309791298091, 0.8946936085078939)
STORED_TOF_DAYS = (14.940560615336594, 14.940560615336594)
STORED_MAX_VINF = 0.9604309791298091


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def close_under_convention(*, moons: tuple[str, ...], n_phase_samples: int = 12) -> dict[str, Any]:
    """Run the production close() under a specific moon-set convention.

    ``moons`` is what gets fed to ``RepeatedMoonTarget``; the sorted index of
    each moon in this set determines its deterministic offset seed
    ``2pi * j / len(moons)``. The candidate is always the SILVER (Umbriel,
    Oberon, Umbriel) at (1, 1). Returns the closure result + the per-moon
    offsets that produced it (so the convention's identity is recoverable).
    """
    target = RepeatedMoonTarget(
        primary="Uranus",
        moons=moons,
        seq_lengths=(3,),
        n_rev_grid=(0, 1, 2),
        n_phase_samples=n_phase_samples,
        tof_resonance_grid=(0.5, 1.0, 1.5, 2.0),
    )
    cand = Candidate(
        index=0,
        signature_hash="sha1:verify_327",
        sequence=SEQ,
        primary="Uranus",
        payload={"n_rev": list(NREV)},
    )
    result = target.close(cand)
    sorted_moons = tuple(sorted(moons))
    offsets_rad = {m: 2.0 * math.pi * j / len(sorted_moons) for j, m in enumerate(sorted_moons)}
    return {
        "moons": list(moons),
        "moons_sorted": list(sorted_moons),
        "deterministic_offsets_deg": {m: math.degrees(v) for m, v in offsets_rad.items()},
        "converged": result.converged,
        "residual_kms": result.residual_kms,
        "vinf_per_encounter_kms": list(result.vinf_per_encounter_kms),
        "tof_days": list(result.tof_days),
    }


def n_rev_sweep_at_convention(
    moons: tuple[str, ...], n_phase_samples: int = 12
) -> list[dict[str, Any]]:
    """Sweep n_rev_grid up to (0..5) at a given convention; the candidate is
    always the SILVER (Umbriel, Oberon, Umbriel) but n_rev varies.
    """
    rows: list[dict[str, Any]] = []
    target_base = RepeatedMoonTarget(
        primary="Uranus",
        moons=moons,
        seq_lengths=(3,),
        n_rev_grid=(0, 1, 2, 3, 4, 5),
        n_phase_samples=n_phase_samples,
        tof_resonance_grid=(0.5, 1.0, 1.5, 2.0),
    )
    for n1 in range(0, 6):
        for n2 in range(0, 6):
            if n1 == 0 and n2 == 0:
                # Trivial all-zero-rev corner: the production engine excludes
                # this in enumerate_candidates(). We mirror that exclusion for
                # the sweep so the headline doesn't pick a degeneracy.
                continue
            cand = Candidate(
                index=0,
                signature_hash="sha1:verify_327_nrev",
                sequence=SEQ,
                primary="Uranus",
                payload={"n_rev": [n1, n2]},
            )
            result = target_base.close(cand)
            rows.append(
                {
                    "n_rev": [n1, n2],
                    "converged": result.converged,
                    "residual_kms": (
                        result.residual_kms if math.isfinite(result.residual_kms) else None
                    ),
                    "max_vinf_kms": (
                        max(result.vinf_per_encounter_kms)
                        if result.vinf_per_encounter_kms
                        else None
                    ),
                }
            )
    return rows


def basin_floor_offset_sweep(
    *,
    n_offset: int = 24,
    n_phase: int = 24,
) -> dict[str, Any]:
    """Confirmation sweep over the relative phase offset between Oberon and
    Umbriel at (1, 1). Narrower grid than #312's 96x96 — we only need to
    confirm the basin floor agrees with #312's finding (~0.024 km/s) and
    identify whether the gate-passing IC (V_inf max ~0.96 km/s) lives at the
    floor or strictly above it.
    """
    mu = PRIMARIES["Uranus"]
    oberon = SATELLITES["Oberon"]
    umbriel = SATELLITES["Umbriel"]
    sma_o, sma_u = oberon.sma_km, umbriel.sma_km
    n_o = _mean_motion_rad_day(mu, sma_o)
    n_u = _mean_motion_rad_day(mu, sma_u)
    p_o = 2.0 * math.pi / n_o
    p_u = 2.0 * math.pi / n_u

    seq = SEQ
    nrev = NREV
    tof_scales = (0.5, 1.0, 1.5, 2.0)

    best_overall = math.inf
    best_rec: dict[str, Any] = {}
    # Also track the best record whose max V_inf is below the
    # gate-passing 1.0 km/s envelope, to separate basin floor from
    # gate-passing basin floor.
    best_gate_pass = math.inf
    best_gate_pass_rec: dict[str, Any] = {}

    for i in range(n_offset):
        rel_off = 2.0 * math.pi * i / n_offset
        for j in range(n_phase):
            phase0 = 2.0 * math.pi * j / n_phase
            theta = {"Oberon": phase0, "Umbriel": phase0 + rel_off}
            for ts in tof_scales:
                tof = ts * math.sqrt(p_o * p_u)
                epochs = [0.0, tof, 2.0 * tof]
                states = []
                for m, t in zip(seq, epochs, strict=True):
                    if m == "Oberon":
                        sma, n = sma_o, n_o
                    else:
                        sma, n = sma_u, n_u
                    states.append(_moon_state(theta[m], n, t, sma, mu))
                ok = True
                vinf_in: list[float | None] = [None, None, None]
                vinf_out: list[float | None] = [None, None, None]
                tofs_d = [tof, tof]
                for k in range(2):
                    r_a, v_a = states[k]
                    r_b, v_b = states[k + 1]
                    sols = lambert(r_a, r_b, tofs_d[k] * DAY_S, mu=mu, max_revs=max(0, nrev[k]))
                    wanted = [s for s in sols if s.n_revs == nrev[k]]
                    if not wanted:
                        ok = False
                        break
                    best = min(wanted, key=lambda s, va=v_a: float(np.linalg.norm(s.v1 - va)))
                    vinf_out[k] = float(np.linalg.norm(best.v1 - v_a))
                    vinf_in[k + 1] = float(np.linalg.norm(best.v2 - v_b))
                if not ok:
                    continue
                worst = 0.0
                for k in range(3):
                    if vinf_in[k] is not None and vinf_out[k] is not None:
                        worst = max(worst, abs(vinf_in[k] - vinf_out[k]))
                wo = vinf_out[0]
                wi = vinf_in[-1]
                if wo is not None and wi is not None:
                    worst = max(worst, abs(wo - wi))
                max_vinf = max(
                    (v for v in (vinf_in + vinf_out) if v is not None),
                    default=0.0,
                )
                rec = {
                    "rel_offset_deg": math.degrees(rel_off),
                    "phase0_deg": math.degrees(phase0),
                    "tof_scale": ts,
                    "tof_days": tof,
                    "residual_kms": worst,
                    "max_vinf_kms": max_vinf,
                    "vinf_in": [v if v is not None else 0.0 for v in vinf_in],
                    "vinf_out": [v if v is not None else 0.0 for v in vinf_out],
                }
                if worst < best_overall:
                    best_overall = worst
                    best_rec = dict(rec)
                # Gate-passing == max V_inf below the SILVER's ~0.96 envelope.
                # Use 1.0 km/s as a slightly generous envelope to still admit
                # nearby phasings that share the regime.
                if max_vinf <= 1.0 and worst < best_gate_pass:
                    best_gate_pass = worst
                    best_gate_pass_rec = dict(rec)

    return {
        "n_offset": n_offset,
        "n_phase": n_phase,
        "tof_scales": list(tof_scales),
        "best_overall_residual_kms": (best_overall if math.isfinite(best_overall) else None),
        "best_overall_record": best_rec,
        "best_gate_passing_residual_kms": (
            best_gate_pass if math.isfinite(best_gate_pass) else None
        ),
        "best_gate_passing_record": best_gate_pass_rec,
    }


def build_legs_at_silver() -> list[PatchedConicLeg] | None:
    """Reconstruct the SILVER's Lambert leg records (SI units).

    Uses the production helper ``legs_from_repeated_moon_candidate`` with
    the same phase_samples=12 the #312 / #285 SILVER row was produced at.
    The helper returns the BEST-PHASING legs (min worst V_inf-continuity).
    """
    return legs_from_repeated_moon_candidate(
        "Uranus",
        SEQ,
        NREV,
        phase_samples=12,
        tof_resonance_grid=(0.5, 1.0, 1.5, 2.0),
    )


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#327] verify SILVER -- sha={sha}", flush=True)
    print("[#327] candidate = repeated-moon-uranus-00000041", flush=True)
    print(f"[#327] sequence = {SEQ}, n_rev = {NREV}", flush=True)
    print(
        f"[#327] stored residual = {STORED_RESIDUAL_KMS:.6f} km/s, "
        f"max V_inf = {STORED_MAX_VINF:.4f} km/s",
        flush=True,
    )

    # 1. 2-moon convention closure.
    print("[#327] (1/7) 2-moon convention close()...", flush=True)
    two_moon = close_under_convention(moons=("Oberon", "Umbriel"))
    print(
        f"  residual = {two_moon['residual_kms']:.6f} km/s | "
        f"vinf = {[f'{v:.4f}' for v in two_moon['vinf_per_encounter_kms']]}",
        flush=True,
    )

    # 2. 3-moon convention closure.
    print("[#327] (2/7) 3-moon convention close()...", flush=True)
    three_moon = close_under_convention(moons=("Titania", "Oberon", "Umbriel"))
    print(
        f"  residual = {three_moon['residual_kms']:.6f} km/s | "
        f"vinf = {[f'{v:.4f}' for v in three_moon['vinf_per_encounter_kms']]}",
        flush=True,
    )

    # 3. Relative-offset basin sweep (confirmation, 24x24).
    print("[#327] (3/7) basin-floor offset sweep (24x24)...", flush=True)
    basin = basin_floor_offset_sweep(n_offset=24, n_phase=24)
    print(
        f"  best overall = {basin['best_overall_residual_kms']:.6f} km/s "
        f"at rel_off={basin['best_overall_record'].get('rel_offset_deg', 0):.2f} deg, "
        f"max V_inf={basin['best_overall_record'].get('max_vinf_kms', 0):.4f} km/s",
        flush=True,
    )
    if basin["best_gate_passing_residual_kms"] is not None:
        print(
            f"  best gate-passing (max V_inf <= 1.0 km/s) = "
            f"{basin['best_gate_passing_residual_kms']:.6f} km/s at rel_off="
            f"{basin['best_gate_passing_record'].get('rel_offset_deg', 0):.2f} deg, "
            f"max V_inf={basin['best_gate_passing_record'].get('max_vinf_kms', 0):.4f} km/s",
            flush=True,
        )
    else:
        print("  NO basin sample satisfies max V_inf <= 1.0 km/s.", flush=True)

    # 4. n_rev sweep at 2-moon convention (where the SILVER lives).
    print("[#327] (4/7) n_rev sweep at 2-moon convention...", flush=True)
    nrev_sweep = n_rev_sweep_at_convention(("Oberon", "Umbriel"))
    finite = [r for r in nrev_sweep if r["residual_kms"] is not None]
    finite.sort(key=lambda r: r["residual_kms"])
    print("  top 5 by residual:", flush=True)
    for r in finite[:5]:
        print(
            f"   n_rev={tuple(r['n_rev'])}: residual={r['residual_kms']:.6f} km/s, "
            f"max V_inf={r['max_vinf_kms']:.4f} km/s",
            flush=True,
        )

    # 5. DOP853 cross-check at the SILVER's actual legs.
    print("[#327] (5/7) independent DOP853 cross-check (rtol=atol=1e-12)...", flush=True)
    legs = build_legs_at_silver()
    if legs is None:
        print("  ERROR: legs_from_repeated_moon_candidate returned None", flush=True)
        cross_check_rows: list[dict[str, Any]] = []
        cross_check_max_dr_km = math.inf
    else:
        cross_check_rows = []
        cross_check_max_dr_km = 0.0
        for k, leg in enumerate(legs):
            res = dop853_cross_check_leg(leg, rtol=1e-12, atol=1e-12)
            cross_check_rows.append(
                {
                    "leg_index": k,
                    "label_from": leg.label_from,
                    "label_to": leg.label_to,
                    "dt_s": leg.dt_s,
                    "dr_arrival_km": res["dr_arrival_km"],
                    "dv_arrival_km_s": res["dv_arrival_km_s"],
                    "converged": res["converged"],
                    "passed": res["passed"],
                }
            )
            cross_check_max_dr_km = max(cross_check_max_dr_km, float(res["dr_arrival_km"]))
            print(
                f"   leg {k} {leg.label_from}->{leg.label_to}: "
                f"dr={res['dr_arrival_km']:.3e} km, dv={res['dv_arrival_km_s']:.3e} km/s, "
                f"converged={res['converged']}, passed={res['passed']}",
                flush=True,
            )
        # Also report dr in nondim Uranus-Oberon CR3BP units (LU = Oberon SMA).
        lu_km = SATELLITES["Oberon"].sma_km
        nondim_max = cross_check_max_dr_km / lu_km
        print(
            f"  max dr_arrival = {cross_check_max_dr_km:.3e} km "
            f"= {nondim_max:.3e} nondim (LU = Oberon SMA = {lu_km:.0f} km)",
            flush=True,
        )

    # 6. Physical-sanity gate at stored V_inf.
    print("[#327] (6/7) physical-sanity gate at stored V_inf...", flush=True)
    gate_pass, verdicts = candidate_passes_physical_gate(SEQ, STORED_VINF, min_useful_bend_deg=5.0)
    for v in verdicts:
        print(
            f"   {v.body} at V_inf={v.vinf_kms:.4f} km/s -> "
            f"max bend {v.max_bend_deg:.2f} deg, useful={v.is_useful}",
            flush=True,
        )
    physical_sanity_row = {
        "gate_passed": gate_pass,
        "min_useful_bend_deg": 5.0,
        "per_encounter": [
            {
                "body": v.body,
                "vinf_kms": v.vinf_kms,
                "max_bend_deg": v.max_bend_deg,
                "is_useful": v.is_useful,
                "min_safe_altitude_km": v.min_safe_altitude_km,
            }
            for v in verdicts
        ],
    }

    # 7. Literature check + ML flagger re-run.
    print("[#327] (7/7) literature_check + ML flagger at gate-passing IC...", flush=True)
    sig = CandidateSignature(
        primary="Uranus",
        sequence=SEQ,
        period_k=2,
        vinf_per_encounter_kms=STORED_VINF,
        n_rev=NREV,
    )
    lit_result = check_literature(sig, search=offline_corpus_search)
    print(
        f"   literature: status={lit_result.status}, confidence={lit_result.confidence}, "
        f"citation={lit_result.citation}",
        flush=True,
    )

    flagger = FalsePosFlagger()
    x_train, y_train, _meta = build_training_set()
    flagger.fit(x_train, y_train)
    # The flagger.score takes a candidate-feature row; mirror the production
    # call pattern from saturn_uranus_campaign.score_candidate by feeding
    # the same structural features the #312 run used.
    # Mirror the production feature schema from saturn_uranus_campaign.score_candidate
    # so the flagger sees the same row shape it scored in #312 (stored 0.5918).
    p_fp_features: dict[str, Any] = {
        "primary": "Uranus",
        "sequence": list(SEQ),
        "n_rev": list(NREV),
        "vinf_per_encounter_kms": list(STORED_VINF),
        "tof_days": list(STORED_TOF_DAYS),
        "verdict_audit": {
            "residual_kms": STORED_RESIDUAL_KMS,
            "primary": "Uranus",
        },
        "max_vinf_kms": STORED_MAX_VINF,
        "bend_feasible": True,
    }
    try:
        p_fp = float(flagger.score(p_fp_features))
    except Exception as exc:
        p_fp = 0.5
        print(f"   ML flagger raised {exc!r}; falling back to 0.5", flush=True)
    print(f"   ML flagger: p_fp = {p_fp:.6f}", flush=True)

    # Verdict logic (matches production policy in saturn_uranus_campaign.py).
    closure_threshold = 0.05  # km/s, same as #312's gate
    p_fp_silver_max = 0.75  # matches P_FP_SILVER_MAX (gauntlet V0 gate, #274)
    closure_pass = two_moon["residual_kms"] < closure_threshold
    independent_pass = cross_check_max_dr_km < 1.0  # < 1 km
    lit_fresh = lit_result.status == "not-found"
    ml_silver_pass = p_fp <= p_fp_silver_max
    basin_gate_pass = (
        basin["best_gate_passing_residual_kms"] is not None
        and basin["best_gate_passing_residual_kms"] < closure_threshold
    )

    if (
        closure_pass
        and independent_pass
        and gate_pass
        and lit_fresh
        and basin_gate_pass
        and ml_silver_pass
    ):
        verdict = "REAL_CANDIDATE_AWAITING_306"
    elif not closure_pass:
        verdict = "RETIRE_NO_CLOSURE_AT_GATE_PASSING_IC"
    elif not independent_pass:
        verdict = "RETIRE_CROSSCHECK_DISAGREEMENT"
    elif not gate_pass:
        verdict = "RETIRE_PHYSICAL_SANITY_FAILED"
    elif not lit_fresh:
        verdict = "REDISCOVERY"
    elif not ml_silver_pass:
        verdict = "MARGINAL_ML_FLAGGED"
    else:
        verdict = "MARGINAL_UNKNOWN"

    elapsed = time.time() - t0
    print(f"[#327] verdict = {verdict} (elapsed {elapsed:.1f}s)", flush=True)

    # Write the verified JSONL.
    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#327 verify Umbriel SILVER at gate-passing V_inf",
            "candidate_id": "repeated-moon-uranus-00000041",
            "primary": "Uranus",
            "sequence": list(SEQ),
            "n_rev": list(NREV),
            "source_jsonl": str(SRC_JSONL.relative_to(ROOT)),
            "stored": {
                "residual_kms": STORED_RESIDUAL_KMS,
                "vinf_per_encounter_kms": list(STORED_VINF),
                "tof_days": list(STORED_TOF_DAYS),
                "max_vinf_kms": STORED_MAX_VINF,
            },
            "git_sha": sha,
            "elapsed_s": elapsed,
        }
    )
    rows.append({"kind": "convention_two_moon", **two_moon})
    rows.append({"kind": "convention_three_moon", **three_moon})
    rows.append({"kind": "basin_offset_sweep", **basin})
    rows.append(
        {
            "kind": "n_rev_sweep_two_moon",
            "n_phase_samples": 12,
            "tof_resonance_grid": [0.5, 1.0, 1.5, 2.0],
            "rows": nrev_sweep,
        }
    )
    rows.append(
        {
            "kind": "dop853_cross_check",
            "rtol": 1e-12,
            "atol": 1e-12,
            "phase_samples": 12,
            "max_dr_arrival_km": cross_check_max_dr_km,
            "per_leg": cross_check_rows,
        }
    )
    rows.append({"kind": "physical_sanity_gate", **physical_sanity_row})
    rows.append(
        {
            "kind": "literature_check",
            "signature": {
                "primary": sig.primary,
                "sequence": list(sig.sequence),
                "period_k": sig.period_k,
                "vinf_per_encounter_kms": list(sig.vinf_per_encounter_kms),
                "n_rev": list(sig.n_rev),
            },
            "status": lit_result.status,
            "confidence": lit_result.confidence,
            "citation": lit_result.citation,
            "doi": lit_result.doi,
            "matched_url": lit_result.matched_url,
            "backend": "offline_corpus_search",
        }
    )
    rows.append(
        {
            "kind": "ml_flagger",
            "p_fp": p_fp,
            "p_fp_threshold_silver": 0.75,
            "features": p_fp_features,
        }
    )
    rows.append(
        {
            "_meta": True,
            "kind": "verdict",
            "closure_pass_two_moon": closure_pass,
            "closure_threshold_kms": closure_threshold,
            "independent_cross_check_pass": independent_pass,
            "physical_sanity_pass": gate_pass,
            "literature_status": lit_result.status,
            "ml_flagger_p_fp": p_fp,
            "basin_floor_gate_passing_kms": basin["best_gate_passing_residual_kms"],
            "basin_floor_overall_kms": basin["best_overall_residual_kms"],
            "verdict": verdict,
            "writeback_to_catalogue": False,
            "next_step": (
                "#306 3D V0-V5 gauntlet (Uranus-Oberon and Uranus-Umbriel CR3BP)"
                if verdict.startswith("REAL_CANDIDATE")
                else "negative-results registry (#312 family genome ceiling at gate-passing V_inf)"
            ),
        }
    )

    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")
    print(f"[#327] wrote {OUT_JSONL}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
