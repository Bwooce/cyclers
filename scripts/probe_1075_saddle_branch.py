"""#1075 probe: what does every rung of branch_at_saddle_center's epsilon ladder do today?

Parent = data/floquet_phase1_c32_family.jsonl member i=124 (the test's parent).
Prints the picked eigenvalue/eigenvector, then for each (eps, sign) the corrector's
converged flag, residual, independent closure, n_iter, final state, Jacobi and the
winding topology when converged. Also prints the committed reproduction record
(data/floquet_phase1_reproduction.jsonl) for the identity comparison.

Usage: uv run python scripts/probe_1075_saddle_branch.py [--eps 5e-4,7e-4,...] [--max-iter 80]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.asymmetric_branch import _select_saddle_center_eigenvector
from cyclerfinder.search.bifurcation_detector import monodromy
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    RESIDUAL_FULL_STATE_AT_T,
    correct_general_periodic_3d,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", default="5e-4,7e-4,3e-4,1e-3,1e-4,1.5e-3")
    ap.add_argument("--max-iter", type=int, default=80)
    ap.add_argument("--index", type=int, default=124)
    ap.add_argument(
        "--svd", action="store_true", help="print the first Newton Jacobian's SVD and step"
    )
    args = ap.parse_args()
    rows = [json.loads(x) for x in (ROOT / "data/floquet_phase1_c32_family.jsonl").open()]
    m = [r for r in rows if r.get("kind") == "member"][args.index]
    system = braik_ross_system()
    s0 = np.array([float(m["x0"]), 0, 0, 0, float(m["ydot0"]), 0], dtype=np.float64)
    T = float(m["period_TU"])  # noqa: N806
    print(f"parent i={args.index} x0={s0[0]!r} ydot0={s0[4]!r} T={T!r} C={m.get('jacobi')}")
    mono = monodromy(system, s0, T, rtol=1e-12, atol=1e-12)
    eig = np.linalg.eigvals(mono)
    print("multipliers:", [f"{complex(e):.12g}" for e in sorted(eig, key=lambda z: abs(z))])
    tangent = cr3bp.cr3bp_eom(0.0, s0, system.mu)
    pick = _select_saddle_center_eigenvector(mono, parent_tangent=tangent)
    assert pick is not None
    v, lam = pick
    print(f"picked lambda={lam!r}")
    print("picked v=", [repr(float(x)) for x in v])
    rec = json.loads((ROOT / "data/floquet_phase1_reproduction.jsonl").open().readline())
    print(
        f"committed record: lambda_re={rec['eigenvalue_used_real']!r} "
        f"topo=({rec['branch_k1']},{rec['branch_k2']}) C={rec['branch_jacobi']!r} "
        f"T={rec['branch_T_TU']!r} eps={rec['epsilon_used']!r} sign={rec['eigenvector_sign_used']}"
    )
    print("committed v=", rec["eigenvector_used"])
    if args.svd:
        from cyclerfinder.search.cr3bp_general_periodic_3d import _newton_step, _propagate_with_stm

        for sign in (+1, -1):
            x = s0 + sign * 5e-4 * v
            sf, stm = _propagate_with_stm(system, x, T, rtol=1e-12, atol=1e-12)
            res = sf - x
            jac = np.zeros((6, 7))
            jac[:, :6] = stm - np.eye(6)
            jac[:, 6] = cr3bp.cr3bp_eom(T, sf, system.mu)
            u, sv, _vt = np.linalg.svd(jac)
            raw, *_ = np.linalg.lstsq(jac, -res, rcond=None)
            capped = _newton_step(jac, res, tuple(range(7)), state_cap=0.1, period_cap=0.2 * T)
            print(f"sign={sign:+d} |res|={np.linalg.norm(res):.3e} sv={[f'{x_:.3e}' for x_ in sv]}")
            print(f"  res projected on left singular vectors: {[f'{x_:.3e}' for x_ in u.T @ res]}")
            print(f"  raw lstsq step={[f'{x_:.3e}' for x_ in raw]}")
            print(f"  capped step={[f'{x_:.3e}' for x_ in capped]}", flush=True)
        return
    parent_topo = winding_topology(system.mu, s0, T)
    print(f"parent topology ({parent_topo.k1},{parent_topo.k2})", flush=True)
    for eps in [float(x) for x in args.eps.split(",")]:
        for sign in (+1, -1):
            t0 = time.time()
            r = correct_general_periodic_3d(
                system,
                s0 + sign * eps * v,
                T,
                free_vars=FREE_VARS_FULL_ASYMMETRIC,
                residual_indices=RESIDUAL_FULL_STATE_AT_T,
                is_half_period_residual=False,
                tol=1e-10,
                max_iter=args.max_iter,
                rtol=1e-12,
                atol=1e-12,
                independent_tol=1e-6,
            )
            topo = "-"
            if r.converged or r.corrector_residual < 1e-6:
                try:
                    tp = winding_topology(system.mu, r.state0, r.T_TU)
                    topo = f"({tp.k1},{tp.k2})"
                except (RuntimeError, ValueError) as exc:
                    topo = f"err:{exc}"
            print(
                f"eps={eps:g} sign={sign:+d} conv={r.converged} res={r.corrector_residual:.3e} "
                f"indep={r.independent_closure_residual:.3e} it={r.n_iter} T={r.T_TU:.10f} "
                f"C={r.jacobi:.10f} topo={topo} planar={r.degenerate_planar} "
                f"x0={[round(float(x), 10) for x in r.state0]} [{time.time() - t0:.0f}s]",
                flush=True,
            )


if __name__ == "__main__":
    main()
