"""#372 P372.3 — STM-mode cross-check at the #347 Phase 1 (3,2) cluster point.

Computes the 6 Floquet multipliers of the monodromy at the (3,2) Earth-Moon
C32 family's post-saddle-center anchor (C = 3.14180; the parent state in
``data/floquet_phase1_reproduction.jsonl``) under BOTH STM-integration modes:

  * ``stm_mode='variable'`` — the legacy variable-step augmented variational
    path (the method #347 Phase 1 actually used).
  * ``stm_mode='fixed_path'`` — Pellegrini-Russell 2016 (JGCD,
    DOI 10.2514/1.G001920) mitigation: state-only DOP853 records the step
    grid; augmented state+STM replays one DOP853 step per recorded
    sub-interval with first_step=max_step=h_i, killing the eq. 17
    step-size IC-dependence.

Writes ``data/372_stm_mode_crosscheck.jsonl`` with one record containing both
multiplier sets, their |lambda_max|, and the relative disagreement at the
dominant eigenvalue. The Phase 0 design doc (line 116) flagged this as
"Unknown whether the eigenvalue separation will survive at the cluster
point"; this script delivers the quantified answer.

Acceptance per the plan (#372 spec):
  * <1% relative disagreement at |lambda_max| -- variable-step substrate is
    acceptable for the (3,2) family at this tolerance.
  * >=1% -- STM contamination confirmed; Phase 2 discovery sweep should
    default to ``stm_mode='fixed_path'``.

Usage:
    uv run python scripts/372_stm_mode_crosscheck.py

Output:
    data/372_stm_mode_crosscheck.jsonl
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy

# Sourced anchor: the (3,2) C32 post-saddle-center parent member emitted by
# scripts/floquet_phase1_p1_5_end_to_end.py (saved to floquet_phase1_reproduction.jsonl
# on 2026-06-17; first record). This C value is the post-bifurcation cluster
# point characterised in #347 Phase 1.
PARENT_STATE0 = np.array(
    [
        -0.28434935291806585,
        0.0,
        0.0,
        0.0,
        -2.0533998110844403,
        0.0,
    ],
    dtype=np.float64,
)
PARENT_PERIOD = 17.4642064064264  # TU
PARENT_JACOBI_TARGET = 3.1418000000000266
# Earth-Moon CR3BP mu per Ross & Roberts-Tsoukkas 2025 (AAS 25-621) p.3 -- the
# project-standard sourced value (1.2150584270572e-2).
EM_MU = 1.2150584270572e-2

OUTPUT = Path("data/372_stm_mode_crosscheck.jsonl")


def _eig_metrics(eigs: np.ndarray) -> dict[str, float | list[float]]:
    """Summarise a multiplier set: |lambda_max|, real/imag parts, ordered."""
    abs_sorted = np.sort(np.abs(eigs))[::-1]
    return {
        "lambda_max_abs": float(abs_sorted[0]),
        "lambda_min_abs": float(abs_sorted[-1]),
        "abs_sorted": [float(v) for v in abs_sorted],
        "real_parts": [float(v.real) for v in eigs],
        "imag_parts": [float(v.imag) for v in eigs],
    }


def main() -> None:
    sysm = cr3bp.CR3BPSystem(mu=EM_MU, primary="Earth", secondary="Moon", l_km=1.0, t_s=1.0)

    # Sanity: Jacobi at the anchor state should match the recorded parent C.
    c_actual = cr3bp.jacobi_constant(PARENT_STATE0, EM_MU)
    print(f"anchor C = {c_actual:.15g}  (target {PARENT_JACOBI_TARGET:.15g})")
    assert abs(c_actual - PARENT_JACOBI_TARGET) < 1e-9, "anchor C off"

    # Variable-step (legacy) monodromy.
    t0 = time.time()
    mono_var = monodromy(sysm, PARENT_STATE0, PARENT_PERIOD, stm_mode="variable")
    dt_var = time.time() - t0
    eigs_var = floquet_multipliers(mono_var)
    metrics_var = _eig_metrics(eigs_var)

    # Fixed-path monodromy (Pellegrini-Russell 2016 mitigation).
    t0 = time.time()
    mono_fix = monodromy(sysm, PARENT_STATE0, PARENT_PERIOD, stm_mode="fixed_path")
    dt_fix = time.time() - t0
    eigs_fix = floquet_multipliers(mono_fix)
    metrics_fix = _eig_metrics(eigs_fix)

    # Disagreement metrics.
    lam_max_var = metrics_var["lambda_max_abs"]
    lam_max_fix = metrics_fix["lambda_max_abs"]
    assert isinstance(lam_max_var, float) and isinstance(lam_max_fix, float)
    rel_disagree_lam_max = abs(lam_max_fix - lam_max_var) / max(lam_max_var, 1e-300)
    abs_sorted_var = np.asarray(metrics_var["abs_sorted"], dtype=float)
    abs_sorted_fix = np.asarray(metrics_fix["abs_sorted"], dtype=float)
    rel_disagree_per_eig = np.abs(abs_sorted_fix - abs_sorted_var) / np.maximum(
        abs_sorted_var, 1e-300
    )
    rel_disagree_max = float(rel_disagree_per_eig.max())
    frob_diff = float(np.linalg.norm(mono_fix - mono_var))
    frob_var = float(np.linalg.norm(mono_var))
    stm_frob_rel = frob_diff / max(frob_var, 1e-300)

    # Acceptance gate per #372 plan: <1% at |lambda_max|.
    acceptance_threshold = 1e-2
    passes_acceptance = bool(rel_disagree_lam_max < acceptance_threshold)

    record = {
        "kind": "stm_mode_crosscheck",
        "issue": "#372",
        "subissue": "P372.3",
        "phase": "347_phase1_cluster_point",
        "sourced_anchor": "braik_ross_2026_table2_C32_post_saddle_center",
        "reference": ("Pellegrini, E. & Russell, R.P. (2016), JGCD, DOI 10.2514/1.G001920"),
        "anchor": {
            "state0": PARENT_STATE0.tolist(),
            "period_TU": PARENT_PERIOD,
            "jacobi": float(c_actual),
            "mu_EM": EM_MU,
        },
        "variable_step": {
            "metrics": metrics_var,
            "wall_seconds": dt_var,
            "rtol": 1e-12,
            "atol": 1e-12,
        },
        "fixed_path": {
            "metrics": metrics_fix,
            "wall_seconds": dt_fix,
            "rtol": 1e-12,
            "atol": 1e-12,
        },
        "disagreement": {
            "rel_disagree_lam_max": float(rel_disagree_lam_max),
            "rel_disagree_max_over_eigs": rel_disagree_max,
            "stm_matrix_rel_frob": stm_frob_rel,
        },
        "acceptance": {
            "threshold": acceptance_threshold,
            "metric": "rel_disagree_lam_max",
            "passes": passes_acceptance,
            "verdict": (
                "variable-step substrate accepted for the (3,2) family at this tolerance"
                if passes_acceptance
                else "STM contamination confirmed; switch Phase 2 discovery "
                "to stm_mode='fixed_path'"
            ),
        },
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w") as f:
        f.write(json.dumps(record) + "\n")

    print(f"\nwritten: {OUTPUT}")
    print(f"variable |lambda_max| = {lam_max_var:.6e}")
    print(f"fixed_path |lambda_max| = {lam_max_fix:.6e}")
    print(f"rel_disagree_lam_max = {rel_disagree_lam_max:.6e}")
    print(f"rel_disagree_max_over_eigs = {rel_disagree_max:.6e}")
    print(f"STM Frobenius relative diff = {stm_frob_rel:.6e}")
    print(f"acceptance (threshold {acceptance_threshold}): {passes_acceptance}")


if __name__ == "__main__":
    main()
