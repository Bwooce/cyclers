"""#311 Phase 1 Part C -- 3D CR3BP existence probe at Saturn-Rhea / Saturn-Dione.

Direct analog of the #312 Phase 1 Part C 3D probe on Uranus. The brief
calls out a "critical insight from #312": 0/24 seeds converged at Uranus
3D from a small-z0 L1-Lyapunov-tail seed. Saturn may behave similarly OR
differently because Saturn's larger primary mass (mu_Saturn-system =
3.79e7 vs mu_Uranus = 5.79e6) and the denser ring of regular moons
produces different out-of-plane resonant structure.

The best (k1, k2) near-miss from Parts A+B is Rhea-Dione-Rhea (1, 1) at
0.107 km/s -- the relevant CR3BPs are CR3BP(Saturn, Rhea) and
CR3BP(Saturn, Dione). The repeated-moon Lambert genome is planet-frame
coplanar by construction, so the probe asks the SAME narrower question
as #312-C: do 3D periodic orbits exist AT ALL in either CR3BP from a
small-z0 L1-Lyapunov-tail seed? A clean negative is OK; a clean positive
seeds Phase 2 (3D cycler maintenance / inclination trade).

The 3D corrector (#291 Phase 1+2) at
``src/cyclerfinder/search/cr3bp_general_periodic_3d.py`` is called
READ-ONLY. Seed grid + per-seed SIGALRM timeout match #312 Part C
exactly so the comparison is direct.

Output: ``data/scan_311_saturn_3d_probe.jsonl``.

Run as::

    uv run python scripts/scan_311_saturn_3d_probe.py
"""

from __future__ import annotations

import json
import math
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


class _SeedTimeoutError(Exception):
    """Raised by the per-seed SIGALRM handler -- bounded budget exceeded."""


def _alarm_handler(_signum: int, _frame: Any) -> None:  # pragma: no cover - signal path
    raise _SeedTimeoutError()


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core import cr3bp  # noqa: E402
from cyclerfinder.search.cr3bp_general_periodic_3d import (  # noqa: E402
    FREE_VARS_FULL_ASYMMETRIC,
    correct_general_periodic_3d,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _l1_distance(mu: float, *, max_iter: int = 80, tol: float = 1e-14) -> float:
    """Hill-cube + Newton refine of the L1 distance from the secondary.

    See #312 Part C script for derivation (Szebehely eq 4.4-16). Standard
    fifth-order polynomial; Newton converges in <10 iterations from the
    Hill seed.
    """
    gamma = (mu / 3.0) ** (1.0 / 3.0)  # Hill-cube starting guess
    for _ in range(max_iter):
        f = (
            gamma**5
            - (3.0 - mu) * gamma**4
            + (3.0 - 2.0 * mu) * gamma**3
            - mu * gamma**2
            + 2.0 * mu * gamma
            - mu
        )
        fp = (
            5.0 * gamma**4
            - 4.0 * (3.0 - mu) * gamma**3
            + 3.0 * (3.0 - 2.0 * mu) * gamma**2
            - 2.0 * mu * gamma
            + 2.0 * mu
        )
        if abs(fp) < 1e-30:
            break
        delta = f / fp
        gamma -= delta
        if abs(delta) < tol:
            break
    return float(gamma)


def _hill_period_guess(mu: float, gamma: float) -> float:
    """Hill-frame linearized vertical-Lyapunov period guess (nondim TU).

    See #312 Part C derivation: omega_z = sqrt(c_2), c_2 = 1 + 2 mu / gamma**3.
    """
    c_2 = 1.0 + 2.0 * mu / gamma**3
    omega_z = math.sqrt(c_2)
    return 2.0 * math.pi / omega_z


def probe_3d_at_system(
    primary: str,
    secondary: str,
    *,
    z0_grid: tuple[float, ...] = (1e-4, 1e-3, 5e-3, 1e-2),
    ydot0_grid: tuple[float, ...] = (-0.5, -0.2, -0.1, 0.1, 0.2, 0.5),
    tol: float = 1e-10,
    independent_tol: float = 1e-6,
    max_iter: int = 30,
    per_seed_timeout_s: int = 120,
) -> list[dict[str, Any]]:
    """Probe for 3D periodic orbits in CR3BP(primary, secondary).

    Seed structure: x0 just inside L1 (halo-family tail), z0 in z0_grid,
    ydot0 in ydot0_grid; perpendicular convention (y0 = xdot0 = zdot0 = 0).
    Period guess from Hill linearization. Full asymmetric 3D corrector.
    """
    system = cr3bp.cr3bp_system(primary, secondary)
    mu = system.mu
    gamma = _l1_distance(mu)
    x_l1 = 1.0 - mu - gamma
    x0 = x_l1 - 0.95 * gamma
    t_guess = _hill_period_guess(mu, gamma)
    rows: list[dict[str, Any]] = []

    prev_handler = signal.signal(signal.SIGALRM, _alarm_handler)
    try:
        for z0 in z0_grid:
            for ydot0 in ydot0_grid:
                seed = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
                t0 = time.time()
                signal.alarm(int(per_seed_timeout_s))
                try:
                    orb = correct_general_periodic_3d(
                        system,
                        seed,
                        t_guess,
                        free_vars=FREE_VARS_FULL_ASYMMETRIC,
                        tol=tol,
                        independent_tol=independent_tol,
                        max_iter=max_iter,
                    )
                    row = {
                        "primary": primary,
                        "secondary": secondary,
                        "mu": mu,
                        "x_l1": x_l1,
                        "gamma_l1": gamma,
                        "seed_x0": x0,
                        "seed_z0": z0,
                        "seed_ydot0": ydot0,
                        "seed_period_guess_TU": t_guess,
                        "converged": orb.converged,
                        "corrector_residual": orb.corrector_residual,
                        "independent_closure_residual": orb.independent_closure_residual,
                        "n_iter": orb.n_iter,
                        "degenerate_planar": orb.degenerate_planar,
                        "state0": list(orb.state0),
                        "T_TU": orb.T_TU,
                        "jacobi": orb.jacobi,
                        "elapsed_s": time.time() - t0,
                        "error": None,
                    }
                except _SeedTimeoutError:
                    row = {
                        "primary": primary,
                        "secondary": secondary,
                        "mu": mu,
                        "x_l1": x_l1,
                        "gamma_l1": gamma,
                        "seed_x0": x0,
                        "seed_z0": z0,
                        "seed_ydot0": ydot0,
                        "seed_period_guess_TU": t_guess,
                        "converged": False,
                        "corrector_residual": None,
                        "independent_closure_residual": None,
                        "n_iter": None,
                        "degenerate_planar": None,
                        "state0": None,
                        "T_TU": None,
                        "jacobi": None,
                        "elapsed_s": time.time() - t0,
                        "error": f"per-seed timeout {per_seed_timeout_s}s exceeded",
                    }
                except Exception as exc:
                    row = {
                        "primary": primary,
                        "secondary": secondary,
                        "mu": mu,
                        "x_l1": x_l1,
                        "gamma_l1": gamma,
                        "seed_x0": x0,
                        "seed_z0": z0,
                        "seed_ydot0": ydot0,
                        "seed_period_guess_TU": t_guess,
                        "converged": False,
                        "corrector_residual": None,
                        "independent_closure_residual": None,
                        "n_iter": None,
                        "degenerate_planar": None,
                        "state0": None,
                        "T_TU": None,
                        "jacobi": None,
                        "elapsed_s": time.time() - t0,
                        "error": repr(exc),
                    }
                finally:
                    signal.alarm(0)
                rows.append(row)
                print(
                    f"[311-C] {primary}-{secondary} z0={z0:.0e} "
                    f"ydot0={ydot0:+.2f} converged={row['converged']} "
                    f"corr_res={row['corrector_residual']!r} "
                    f"ind_res={row['independent_closure_residual']!r} "
                    f"planar={row['degenerate_planar']} "
                    f"jac={row['jacobi']!r} T_TU={row['T_TU']!r} "
                    f"elapsed={row['elapsed_s']:.1f}s err={row['error']}",
                    flush=True,
                )
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev_handler)
    return rows


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_311_saturn_3d_probe.jsonl"
    print(f"[311-C] Saturn 3D CR3BP existence probe -- sha={sha}", flush=True)
    print(f"[311-C] out={out_path}", flush=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, Any]] = []
    secondaries = ("Rhea", "Dione")  # best near-miss pair from Parts A+B
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#311 Phase 1 Part C -- Saturn 3D CR3BP existence probe",
                    "primary": "Saturn",
                    "secondaries": list(secondaries),
                    "z0_grid": [1e-4, 1e-3, 5e-3, 1e-2],
                    "ydot0_grid": [-0.5, -0.2, -0.1, 0.1, 0.2, 0.5],
                    "free_vars": list(FREE_VARS_FULL_ASYMMETRIC),
                    "tol": 1e-10,
                    "independent_tol": 1e-6,
                    "git_sha": sha,
                    "note": (
                        "Best (k1,k2) near-miss from Parts A+B is "
                        "Rhea-Dione-Rhea (1, 1) at 0.107 km/s; the "
                        "relevant CR3BPs are CR3BP(Saturn, Rhea) and "
                        "CR3BP(Saturn, Dione). Identical seed strategy + "
                        "grid as #312 Part C for direct Uranus-vs-Saturn "
                        "comparison."
                    ),
                }
            )
            + "\n"
        )
        fh.flush()
        for secondary in secondaries:
            print(f"[311-C] Probing Saturn-{secondary} CR3BP...", flush=True)
            rows = probe_3d_at_system("Saturn", secondary)
            for row in rows:
                fh.write(json.dumps(row) + "\n")
                fh.flush()
                all_rows.append(row)

        # Top-level verdict per system.
        verdict: dict[str, dict[str, Any]] = {}
        for secondary in secondaries:
            sys_rows = [r for r in all_rows if r["secondary"] == secondary]
            converged = [r for r in sys_rows if r["converged"] and not r["degenerate_planar"]]
            converged_planar = [r for r in sys_rows if r["converged"] and r["degenerate_planar"]]
            verdict[secondary] = {
                "n_seeds": len(sys_rows),
                "n_converged_3d": len(converged),
                "n_converged_planar_collapse": len(converged_planar),
                "n_failed": len(sys_rows) - len(converged) - len(converged_planar),
                "best_3d_closure": (
                    min(
                        (r["independent_closure_residual"] for r in converged),
                        default=None,
                    )
                ),
                "best_3d_jacobi_T_state0": [
                    {
                        "jacobi": r["jacobi"],
                        "T_TU": r["T_TU"],
                        "state0": r["state0"],
                        "z0": r["seed_z0"],
                        "ydot0": r["seed_ydot0"],
                        "ind_res": r["independent_closure_residual"],
                    }
                    for r in sorted(converged, key=lambda x: x["independent_closure_residual"])[:5]
                ],
            }
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "summary",
                    "verdict": verdict,
                    "git_sha": sha,
                }
            )
            + "\n"
        )
        fh.flush()
        print(f"[311-C] DONE -- verdict: {json.dumps(verdict, indent=2)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
