"""#389 P389.1 — V1 re-confirm for branch_C32_b0.

Re-runs the Phase 2 branch corrector against the sanitized verification IC at
``data/branch_c32_b0_ic.jsonl``. The branched orbit's IC is passed through
:func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
in full-asymmetric mode (free vars = (x0, y0, z0, xdot0, ydot0, zdot0, T);
residual = full 6D state closure at T), and the resulting
:class:`Periodic3DOrbit` carries both:

* the corrector residual (DOP853 inside the Newton loop), and
* the independent Radau closure check (different integrator).

The V1 verdict is the compound: corrector residual < 1e-10 AND independent
closure < 1e-6 (the Periodic3DOrbit.converged compound gate).

The verdict JSONL row is written to ``data/branch_c32_b0_v1_verdict.jsonl``;
the frozen-gate pytest at ``tests/verify/test_branch_c32_b0_v1_passes.py``
asserts the recorded values against the V1 floor.

Discipline notes:

* The branch IC came from the Phase 2 sweep already. Re-running the corrector
  here is the structural re-confirmation: the IC must close to V1 spec the same
  way it did under the Phase 2 driver. If it doesn't, P389.1 FAILS at the
  closure gate and we HALT.
* "Closure" here is the CR3BP periodic-orbit closure (state at T - state at 0),
  NOT the moontour-flavoured per-leg V_inf-continuity floor the #327 SILVER V1
  used. branch_C32_b0 is a closed CR3BP orbit, not a patched-Lambert tour.

Usage:
    uv run python scripts/branch_c32_b0_v1_verify.py [--output PATH]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    RESIDUAL_FULL_STATE_AT_T,
    correct_general_periodic_3d,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
IC_PATH = Path("data/branch_c32_b0_ic.jsonl")
V1_CORRECTOR_RESIDUAL_FLOOR = 1.0e-10
V1_INDEPENDENT_CLOSURE_FLOOR = 1.0e-6
PHASE_LABEL = "389_p389_1"


def _load_ic() -> dict[str, object]:
    """Read the ``kind: ic`` row from the sanitized verification IC file."""
    with IC_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, object] = json.loads(line)
            if row.get("kind") == "ic" and row.get("candidate_id") == CANDIDATE_ID:
                return row
    raise AssertionError(f"IC row for {CANDIDATE_ID!r} not found in {IC_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_c32_b0_v1_verdict.jsonl"),
        help="Path to write the V1 verdict JSONL row.",
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    ic = _load_ic()
    state0 = np.array(ic["state0_rotating_nondim"], dtype=np.float64)
    period_guess = float(ic["period_TU"])  # type: ignore[arg-type]
    system = braik_ross_system()

    print(f"[V1-verify] candidate_id={CANDIDATE_ID}")
    print(f"[V1-verify] state0={state0}")
    print(f"[V1-verify] period_guess={period_guess:.15f} TU")
    print(f"[V1-verify] system.mu={system.mu}")
    print("[V1-verify] running correct_general_periodic_3d (full asymmetric, 6D at T)...")

    result = correct_general_periodic_3d(
        system,
        state0,
        period_guess,
        free_vars=FREE_VARS_FULL_ASYMMETRIC,
        residual_indices=RESIDUAL_FULL_STATE_AT_T,
        is_half_period_residual=False,
        tol=1e-12,
        max_iter=80,
        rtol=1e-12,
        atol=1e-12,
        independent_tol=V1_INDEPENDENT_CLOSURE_FLOOR,
        require_monotone_decrease=False,
    )
    elapsed = time.time() - t_start

    print(f"[V1-verify] corrector residual: {result.corrector_residual:.3e}")
    print(f"[V1-verify] independent closure: {result.independent_closure_residual:.3e}")
    print(f"[V1-verify] converged (compound gate): {result.converged}")
    print(f"[V1-verify] n_iter: {result.n_iter}")
    print(f"[V1-verify] T_TU corrected: {result.T_TU:.15f}")
    print(f"[V1-verify] jacobi corrected: {result.jacobi:.15f}")
    print(f"[V1-verify] degenerate_planar: {result.degenerate_planar}")
    print(f"[V1-verify] elapsed: {elapsed:.2f}s")

    passes_v1 = bool(
        result.converged
        and result.corrector_residual < V1_CORRECTOR_RESIDUAL_FLOOR
        and result.independent_closure_residual < V1_INDEPENDENT_CLOSURE_FLOOR
    )

    iso_end = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        header = {
            "kind": "header",
            "candidate_id": CANDIDATE_ID,
            "phase": PHASE_LABEL,
            "iso_start": iso_start,
            "iso_end": iso_end,
            "elapsed_seconds": elapsed,
            "v1_corrector_residual_floor": V1_CORRECTOR_RESIDUAL_FLOOR,
            "v1_independent_closure_floor": V1_INDEPENDENT_CLOSURE_FLOOR,
            "discipline": (
                "Spec §14 V1 for a CR3BP periodic orbit: corrector residual "
                "< 1e-10 AND independent (Radau) closure < 1e-6. The compound "
                "Periodic3DOrbit.converged gate enforces both."
            ),
        }
        fh.write(json.dumps(header) + "\n")
        verdict = {
            "kind": "v1_verdict_cr3bp_periodic",
            "candidate_id": CANDIDATE_ID,
            "system": "Earth-Moon CR3BP (Braik-Ross / Ross-RT mass ratio)",
            "mu": float(system.mu),
            "state0_input_nondim": [float(x) for x in state0],
            "period_input_TU": period_guess,
            "state0_corrected_nondim": [float(x) for x in result.state0],
            "period_corrected_TU": float(result.T_TU),
            "jacobi_corrected": float(result.jacobi),
            "corrector_residual": float(result.corrector_residual),
            "independent_closure_residual": float(result.independent_closure_residual),
            "n_iter": int(result.n_iter),
            "degenerate_planar": bool(result.degenerate_planar),
            "converged": bool(result.converged),
            "passes_v1": passes_v1,
            "v1_corrector_residual_floor": V1_CORRECTOR_RESIDUAL_FLOOR,
            "v1_independent_closure_floor": V1_INDEPENDENT_CLOSURE_FLOOR,
        }
        fh.write(json.dumps(verdict) + "\n")
        footer = {
            "kind": "footer",
            "candidate_id": CANDIDATE_ID,
            "phase": PHASE_LABEL,
            "iso_end": iso_end,
            "passes_v1": passes_v1,
        }
        fh.write(json.dumps(footer) + "\n")

    print(f"[V1-verify] verdict written to {args.output}")
    print(f"[V1-verify] passes_v1: {passes_v1}")


if __name__ == "__main__":
    main()
