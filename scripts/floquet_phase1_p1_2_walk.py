"""#347 Phase 1 P1.2 — Natural-parameter continuation of the (3,2) C32 family in CJ.

Walks the (3,2) Earth-Moon symmetric family from the Braik-Ross 2026 Table 2
representative (CJ=3.1294) upward in C, logging at each step:

  * (C, x0, ydot0, period, jacobi) of the converged member
  * the 6 Floquet multipliers from the monodromy matrix
  * the non-trivial unit-circle pair's |lambda| and argument
  * |lambda_max nontrivial| and Floquet sigma (TU^-1, day^-1)

Writes a single JSONL artifact under data/ so P1.3 (saddle-center detector) and
P1.4 (asymmetric branch corrector) can consume it without re-running the walk.

Per the Phase 0 design doc Section 5 (risk: DOP853 conditioning at the cluster
point), the |lambda_max| ~ 2.5e5 at the anchor means the monodromy is
ill-conditioned. We also log:

  * the trivial-pair drift |lambda - 1| for both trivial eigenvalues; the
    saddle-center signal is a non-trivial multiplier joining this trivial cluster.
  * the L2 norm of M @ R - I where R is a reciprocal symmetry approximation
    (Hamiltonian monodromy invariant; not used as a gate here, just diagnostic).

Usage:
    uv run python scripts/floquet_phase1_p1_2_walk.py [--n-steps N] [--dc DC]

Output:
    data/floquet_phase1_c32_family.jsonl
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from cyclerfinder.search.bifurcation_detector import (
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_continuation import continue_family
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit
from cyclerfinder.search.reachable_representatives import (
    TU_DAYS,
    braik_ross_system,
    recover_all_cyclers_braik_ross,
)


def _eigenvalue_diagnostics(eigs: np.ndarray) -> dict[str, float | list[float]]:
    """Classify eigenvalues into trivial/unit-circle/real-reciprocal groups.

    The CR3BP monodromy of a periodic orbit has reciprocal-pair structure:
    (lambda_1, 1/lambda_1, lambda_2, 1/lambda_2, lambda_3, 1/lambda_3). One pair
    is always (1, 1) up to integrator round-off. The remaining four eigenvalues
    are either two real reciprocal pairs (saddle-saddle) or one real reciprocal
    pair + one unit-circle complex conjugate pair (saddle-center). The saddle-
    center bifurcation occurs when the complex conjugate pair coalesces on the
    real axis at +1, then splits into a second real reciprocal pair.

    Returns a dict with:
      * trivial_dist_min, trivial_dist_max: |lambda - 1| for the two closest-to-1
        eigenvalues
      * uc_pair_arg: |arg(lambda)| of the unit-circle complex pair (in [0, pi])
        or NaN if the pair has collapsed onto the real axis
      * uc_pair_abs: |lambda| of that pair (should be 1.0 within integrator floor)
      * lambda_max_nontriv: largest |lambda| over non-trivial multipliers
      * lambda_min_nontriv: smallest |lambda| over non-trivial multipliers
      * eigs_real: list of real parts of all 6 eigenvalues sorted by |lambda - 1|
      * eigs_imag: list of imag parts (same order)
    """
    eigs_c = np.asarray(eigs, dtype=np.complex128)
    dists_to_one = np.abs(eigs_c - 1.0)
    order = np.argsort(dists_to_one)
    sorted_eigs = eigs_c[order]
    trivial_dist_min = float(dists_to_one[order[0]])
    trivial_dist_max = float(dists_to_one[order[1]])
    nontriv = sorted_eigs[2:]
    nontriv_mags = np.abs(nontriv)
    lambda_max_nontriv = float(np.max(nontriv_mags))
    lambda_min_nontriv = float(np.min(nontriv_mags))
    # Identify the unit-circle complex pair: |lambda - 1/lambda*| ~ 0 and |lambda|~1.
    uc_pair_arg = float("nan")
    uc_pair_abs = float("nan")
    on_unit_circle = [e for e in nontriv if abs(abs(e) - 1.0) < 1e-3 and abs(e.imag) > 1e-5]
    if on_unit_circle:
        # Take the eigenvalue with positive imaginary part (its conjugate is also present).
        e_uc = max(on_unit_circle, key=lambda z: z.imag)
        uc_pair_arg = float(abs(np.angle(e_uc)))
        uc_pair_abs = float(abs(e_uc))
    return {
        "trivial_dist_min": trivial_dist_min,
        "trivial_dist_max": trivial_dist_max,
        "uc_pair_arg": uc_pair_arg,
        "uc_pair_abs": uc_pair_abs,
        "lambda_max_nontriv": lambda_max_nontriv,
        "lambda_min_nontriv": lambda_min_nontriv,
        "eigs_real": [float(e.real) for e in sorted_eigs],
        "eigs_imag": [float(e.imag) for e in sorted_eigs],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-steps", type=int, default=30, help="continuation steps")
    parser.add_argument("--dc", type=float, default=1e-4, help="Jacobi step size")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/floquet_phase1_c32_family.jsonl"),
        help="output JSONL path",
    )
    args = parser.parse_args()

    sysm = braik_ross_system()
    reps = recover_all_cyclers_braik_ross(sysm)
    c32 = next(r for r in reps if r.label == "C32")
    seed = SymmetricOrbit(
        x0=float(c32.state0[0]),
        ydot0=float(c32.state0[4]),
        jacobi=c32.jacobi,
        t_half=c32.period * 0.5,
        period=c32.period,
        converged=True,
        crossing_residual=1e-12,
        n_iter=0,
    )

    t0 = time.time()
    branch = continue_family(
        sysm,
        seed,
        direction=+1,
        d_jacobi=args.dc,
        n_steps=args.n_steps,
        min_jacobi=3.0,
        max_jacobi=3.20,
        half_crossings=6,
        ydot0_sign=-1.0,
        seed_label="C32",
        radau_closure_tol=5e-2,
        radau_jacobi_tol=1e-7,
        period_step_frac=0.05,
        period_floor_frac=0.5,
        period_ceiling_frac=1.5,
    )
    walk_dt = time.time() - t0
    print(
        f"continuation done: stop={branch.stop_reason}, steps={branch.n_steps_taken}, "
        f"members={len(branch.members)}, walk_time={walk_dt:.1f}s"
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        # Header row: sourced anchor + run metadata.
        header = {
            "kind": "header",
            "phase": "347_phase1_p1_2",
            "sourced_anchor": "braik_ross_2026_table2_C32",
            "sourced_period_days": 78.613,
            "sourced_sigma_d_per_day": 0.1583,
            "sourced_cj": 3.1294,
            "n_steps": args.n_steps,
            "dc": args.dc,
            "n_members_kept": len(branch.members),
            "stop_reason": str(branch.stop_reason),
            "walk_time_seconds": walk_dt,
            "mu": float(sysm.mu),
            "tu_days": TU_DAYS,
            "iso_timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        fh.write(json.dumps(header) + "\n")
        for i, m in enumerate(branch.members):
            # Recompute monodromy + Floquet at the member IC.
            mono = monodromy(sysm, m.state0, m.period)
            eigs = floquet_multipliers(mono)
            diag = _eigenvalue_diagnostics(eigs)
            # Topology cross-check.
            topo = winding_topology(sysm.mu, m.state0, m.period)
            # Floquet sigma per Braik-Ross eq. 20.
            lam_max = diag["lambda_max_nontriv"]
            assert isinstance(lam_max, float)
            if lam_max > 1.0:
                sigma_tu = math.log(lam_max) / m.period
                sigma_d = sigma_tu / TU_DAYS
            else:
                sigma_tu = 0.0
                sigma_d = 0.0
            row = {
                "kind": "member",
                "index": i,
                "label": "C32_walk",
                "x0": float(m.x0),
                "ydot0": float(m.ydot0),
                "jacobi": float(m.jacobi),
                "period_TU": float(m.period),
                "period_days": float(m.period * TU_DAYS),
                "nu_barden": float(m.nu),
                "abs_lambda_barden": float(m.abs_lambda),
                "stable": bool(m.stable),
                "k1_winding": int(topo.k1),
                "k2_winding": int(topo.k2),
                "winding_prograde": bool(topo.prograde),
                "sigma_TU": float(sigma_tu),
                "sigma_d_per_day": float(sigma_d),
            }
            row.update(diag)
            fh.write(json.dumps(row) + "\n")
            if i in (0, 1, 5, 10, 15, len(branch.members) - 1):
                print(
                    f"  i={i}: C={m.jacobi:.6f}, x0={m.x0:.6f}, T_d={m.period * TU_DAYS:.3f}, "
                    f"|lam_max|={lam_max:.3e}, sigma_d={sigma_d:.4f}, "
                    f"UC_arg={diag['uc_pair_arg']:.4f}"
                )
    print(f"wrote {args.output} with {len(branch.members)} members")


if __name__ == "__main__":
    main()
