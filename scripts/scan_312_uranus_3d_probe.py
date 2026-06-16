"""#312 Phase 1 Part C -- 3D CR3BP existence probe at Uranus-moon systems.

The repeated-moon multi-rev genome (#254) used for Parts A and B is a
patched-conic Lambert search in the planet-inertial frame. The moons orbit
in Uranus's equatorial plane, so the planet-frame search IS coplanar by
construction. Adding ``z0 != 0`` to a planet-frame Lambert IC would push
the trajectory OUT of the moon orbital plane, and it would no longer
encounter the moons -- so the brief's "small z0 != 0" extension does NOT
map onto the repeated-moon genome itself.

What DOES map is a 3D CR3BP existence probe at the relevant moon-pair
system. The cycler that emerged from Part B (Umbriel-Oberon-Umbriel (1,1)
at residual 0.025 km/s) spends most of its cycle in the planet-frame
around Uranus, but the moon flybys themselves happen in the Uranus-Oberon
or Uranus-Umbriel CR3BP. If 3D (halo / NRHO / vertical Lyapunov) orbits
exist in either CR3BP, they offer a Phase-2 maintenance / inclination-
trade lever for the planar cycler (Phase 2 territory, deliberately
out-of-scope here).

This script answers a narrower Phase-1 question: at the fresh-mu
Uranus-Oberon and Uranus-Umbriel CR3BPs, can the 3D corrector
(#291's ``correct_general_periodic_3d``) close ANY 3D periodic orbit
from a Lyapunov-like + small-z0 seed? A clean negative is OK; a clean
positive is the seed for #312 Phase 2.

NO catalogue writeback. The 3D corrector is called READ-ONLY (no
modifications to the module). Seeds are drawn from the L1 Lyapunov family
geometry: x0 = 1 - mu - r_L1 (just inside L1), small ydot0 from a Hill
approximation, sweep z0 in {1e-4, 1e-3, 1e-2}.

Output: ``data/scan_312_uranus_3d_probe.jsonl``.

Run as::

    uv run python scripts/scan_312_uranus_3d_probe.py
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

    L1 is the equilibrium between primary and secondary; gamma = distance
    from secondary, satisfies
    ``gamma**5 - (3 - mu) gamma**4 + (3 - 2 mu) gamma**3 - mu gamma**2 +
    2 mu gamma - mu = 0`` (Szebehely eq 4.4-16). Standard fifth-order
    polynomial; Newton converges in <10 iterations from the Hill seed.
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

    At L1, the linearized vertical motion has angular frequency
    ``omega_z = sqrt(c_2)`` where ``c_2`` is a Legendre coefficient of the
    expansion of the CR3BP potential about L1 (Jorba-Masdemont). Use the
    leading approximation ``c_2 = 1 + 2 mu / gamma**3``; period = 2 pi /
    omega_z.
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

    Seed structure: place ``x0`` JUST INSIDE L1 (a halo-family-tail
    starting point), set ``y0 = xdot0 = zdot0 = 0`` (the perpendicular
    seed convention), sweep ``z0`` (small out-of-plane offset) and
    ``ydot0`` (the velocity that determines the orbit's energy /
    Jacobi). Period guess from the Hill linearisation. Run the full
    asymmetric 3D corrector and report converged closures.

    Each row records the converged state, period, Jacobi, and closure
    residuals (Newton + independent Radau). ``degenerate_planar`` is
    flagged if the corrector collapsed back to ``z0 == 0`` (i.e. no 3D
    family at this seed).
    """
    system = cr3bp.cr3bp_system(primary, secondary)
    mu = system.mu
    gamma = _l1_distance(mu)
    x_l1 = 1.0 - mu - gamma
    # x0 just INSIDE L1 (toward the primary), so the orbit envelops L1.
    # Use 0.95 * gamma offset from L1 -- close enough to be in the L1
    # halo basin, far enough to give the Newton iteration room.
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
                    f"[312-C] {primary}-{secondary} z0={z0:.0e} "
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
    out_path = ROOT / "data" / "scan_312_uranus_3d_probe.jsonl"
    print(f"[312-C] Uranus 3D CR3BP existence probe -- sha={sha}", flush=True)
    print(f"[312-C] out={out_path}", flush=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, Any]] = []
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#312 Phase 1 Part C -- Uranus 3D CR3BP existence probe",
                    "primary": "Uranus",
                    "secondaries": ["Oberon", "Umbriel"],
                    "z0_grid": [1e-4, 1e-3, 5e-3, 1e-2],
                    "ydot0_grid": [-0.5, -0.2, -0.1, 0.1, 0.2, 0.5],
                    "free_vars": list(FREE_VARS_FULL_ASYMMETRIC),
                    "tol": 1e-10,
                    "independent_tol": 1e-6,
                    "git_sha": sha,
                    "note": (
                        "The repeated-moon Lambert genome is planet-frame "
                        "coplanar by construction. This probe asks the "
                        "narrower Phase-1 question: do 3D periodic orbits "
                        "exist AT ALL in CR3BP(Uranus, Oberon) and "
                        "CR3BP(Uranus, Umbriel) from a small-z0 Lyapunov-"
                        "tail seed? A clean negative is OK; a clean positive "
                        "seeds Phase 2 (3D cycler maintenance / inclination)."
                    ),
                }
            )
            + "\n"
        )
        fh.flush()
        for secondary in ("Oberon", "Umbriel"):
            print(f"[312-C] Probing Uranus-{secondary} CR3BP...", flush=True)
            rows = probe_3d_at_system("Uranus", secondary)
            for row in rows:
                fh.write(json.dumps(row) + "\n")
                fh.flush()
                all_rows.append(row)

        # Top-level verdict per system.
        verdict: dict[str, dict[str, Any]] = {}
        for secondary in ("Oberon", "Umbriel"):
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
        print(f"[312-C] DONE -- verdict: {json.dumps(verdict, indent=2)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
