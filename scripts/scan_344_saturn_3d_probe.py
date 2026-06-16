"""#344 Phase 1 Part C -- 3D CR3BP existence probe at Saturn-moon systems.

Mirrors the #312 Phase 1 Part C / #341 Phase 1 Part C workflow. Asks the
narrower Phase-1 question: do 3D periodic orbits exist AT ALL in
CR3BP(Saturn, Titan) and CR3BP(Saturn, Rhea) from a small-z0 Lyapunov-
tail seed? A clean negative is OK; a clean positive seeds Phase 2.

Critical motivation from the brief: Saturn moons orbit close to Saturn's
equatorial plane -- out-of-plane orbits are structurally less natural
than Uranus's where the planet itself is tipped. We test anyway to keep
the verdict honest. Saturn's pole obliquity to the ecliptic is large
(~26.7 deg), but in the Saturn-EQUATORIAL frame the inner regular moons
are tightly co-planar (Titan inclination 0.35 deg). Iapetus is the
outlier at 14.7 deg, and is not probed here because the Lyapunov-tail
seed convention is for the planar collinear libration orbit at the L1
point of the secondary -- Iapetus's CR3BP exists like any other but the
3D probe seed convention does not couple to Iapetus's eccentric/inclined
orbit relative to Saturn's equator (this is a frame choice; see #341 doc
on the parallel point at Triton's retrograde inclination).

Seeds: x0 just inside L1 (x_L1 - 0.95 * gamma), sweep z0 and ydot0 on
a small log-spaced grid. Run the full asymmetric 3D corrector with
FREE_VARS_FULL_ASYMMETRIC.

NO catalogue writeback. The 3D corrector is called READ-ONLY (no
modifications to the module).

Output: ``data/scan_344_saturn_3d_probe.jsonl``.

Run as::

    uv run python scripts/scan_344_saturn_3d_probe.py
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
    """Per-seed SIGALRM-driven timeout."""


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
    """Newton refinement of L1's distance from the secondary (Szebehely 4.4-16)."""
    gamma = (mu / 3.0) ** (1.0 / 3.0)
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
    """Hill-linearized vertical-Lyapunov period guess (nondim TU)."""
    c_2 = 1.0 + 2.0 * mu / gamma**3
    omega_z = math.sqrt(c_2)
    return 2.0 * math.pi / omega_z


def probe_3d_at_system(
    primary: str,
    secondary: str,
    *,
    z0_grid: tuple[float, ...] = (1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1),
    ydot0_grid: tuple[float, ...] = (-0.5, -0.2, -0.1, 0.1, 0.2, 0.5),
    tol: float = 1e-10,
    independent_tol: float = 1e-6,
    max_iter: int = 30,
    per_seed_timeout_s: int = 120,
) -> list[dict[str, Any]]:
    """Probe for 3D periodic orbits in CR3BP(primary, secondary).

    Identical seed convention to #312 / #341 Part C.
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
                    f"[344-C] {primary}-{secondary} z0={z0:.0e} "
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
    out_path = ROOT / "data" / "scan_344_saturn_3d_probe.jsonl"
    print(f"[344-C] Saturn 3D CR3BP existence probe -- sha={sha}", flush=True)
    print(f"[344-C] out={out_path}", flush=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    z0_grid = (1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
    ydot0_grid = (-0.5, -0.2, -0.1, 0.1, 0.2, 0.5)
    all_rows: list[dict[str, Any]] = []
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#344 Phase 1 Part C -- Saturn 3D CR3BP existence probe",
                    "primary": "Saturn",
                    "secondaries": ["Titan", "Rhea"],
                    "z0_grid": list(z0_grid),
                    "ydot0_grid": list(ydot0_grid),
                    "free_vars": list(FREE_VARS_FULL_ASYMMETRIC),
                    "tol": 1e-10,
                    "independent_tol": 1e-6,
                    "git_sha": sha,
                    "note": (
                        "Saturn-Titan mu ~2.367e-4 (intermediate between Earth-Moon "
                        "1.215e-2 and Neptune-Triton 2.09e-4). Saturn-Rhea mu "
                        "~4.06e-6 (very small -- tiny gamma_L1, Hill radius "
                        "dominated by Saturn). The Lyapunov-tail seed at "
                        "x_L1 - 0.95*gamma probes the deep L1 cylinder; Phase 1 "
                        "asks generic 3D existence at this seed, not Saturn-Titan "
                        "halo families specifically (those are catalogued in "
                        "Brinckerhoff/Lo/Marsden AAS 09-129 and need their "
                        "documented continuation paths to surface)."
                    ),
                }
            )
            + "\n"
        )
        fh.flush()
        for secondary in ("Titan", "Rhea"):
            print(f"[344-C] Probing Saturn-{secondary} CR3BP...", flush=True)
            rows = probe_3d_at_system(
                "Saturn",
                secondary,
                z0_grid=z0_grid,
                ydot0_grid=ydot0_grid,
            )
            for row in rows:
                fh.write(json.dumps(row) + "\n")
                fh.flush()
                all_rows.append(row)

        verdict: dict[str, dict[str, Any]] = {}
        for secondary in ("Titan", "Rhea"):
            sys_rows = [r for r in all_rows if r["secondary"] == secondary]
            converged = [r for r in sys_rows if r["converged"] and not r["degenerate_planar"]]
            converged_planar = [r for r in sys_rows if r["converged"] and r["degenerate_planar"]]
            best_3d = None
            if converged:
                best = min(
                    converged,
                    key=lambda r: r["independent_closure_residual"] or math.inf,
                )
                best_3d = {
                    "seed_z0": best["seed_z0"],
                    "seed_ydot0": best["seed_ydot0"],
                    "independent_closure_residual": best["independent_closure_residual"],
                    "corrector_residual": best["corrector_residual"],
                    "T_TU": best["T_TU"],
                    "jacobi": best["jacobi"],
                }
            verdict[secondary] = {
                "n_seeds": len(sys_rows),
                "n_converged_3d": len(converged),
                "n_converged_planar_collapse": len(converged_planar),
                "n_failed": len(sys_rows) - len(converged) - len(converged_planar),
                "best_3d_brief": best_3d,
            }

        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "verdict",
                    "primary": "Saturn",
                    "verdict": verdict,
                    "interpretation": (
                        "n_converged_3d > 0 -> 3D periodic orbits exist at the "
                        "L1-Lyapunov-tail seed; track best vs the Part A 0.032 "
                        "km/s baseline. n_converged_3d == 0 -> negative at this "
                        "seed (does NOT rule out the family broadly; bifurcation-"
                        "continuation from a known planar member is Phase 2)."
                    ),
                    "git_sha": sha,
                }
            )
            + "\n"
        )
        fh.flush()
    print("[344-C] DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
