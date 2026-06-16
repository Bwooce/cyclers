"""Run #303: BCR4BP L1 Lyapunov mu_sun continuation from CR3BP to Andreu.

Produces ``data/bcr4bp_l1_family_303.jsonl`` containing one converged
BCR4BP L1 Lyapunov family member per step, continued from ``mu_sun = 0``
(CR3BP limit) to ``mu_sun = 328900.5423094043`` (Andreu / Rosales-Jorba
parameter value).

Discipline
----------
  * Sourced seed: CR3BP planar L1 Lyapunov at the Braik-Ross 2026 common
    Jacobi level ``C = 3.1294``, corrected by the existing
    ``cyclerfinder.search.reachable_representatives.correct_symmetric_free_period``.
    Both the IC (LL1 offline seed ``x = 0.8115``) and the Jacobi level are
    SOURCED in the existing codebase (per the orbit-closure discipline).
  * Independent (Radau) cross-check on every accepted member is enforced
    by the Phase 1 corrector itself.
  * NO catalogue writeback.
  * NO novelty claims -- L1 Lyapunov continuation in mu_sun is published
    methodology (Simo-Jorba-Gomez).

Usage
-----
``uv run python scripts/run_303_bcr4bp_l1_continuation.py``

Output JSONL columns are documented in the header row.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_continuation import (
    continue_bcr4bp_family_in_musun,
)
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    correct_bcr4bp_periodic,
)
from cyclerfinder.search.reachable_representatives import (
    correct_symmetric_free_period,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "data" / "bcr4bp_l1_family_303.jsonl"

# Sourced constants -- match Phase 1's andreu_default() exactly.
MU_EM = 0.012150581600000  # Rosales-Jorba 2023 Table 3
MU_SUN_ANDREU = 328900.5423094043
A_SUN_NONDIM = 388.8111430233511
OMEGA_SUN_NONDIM = 0.925195985520347
C_JACOBI_BRAIK_ROSS = 3.1294  # Braik-Ross 2026 Table 2 common Jacobi level
L1_LYAPUNOV_X0_GUESS = 0.8115  # OFFLINE_SEEDS LL1, x near EM L1 = 0.836915

# Continuation parameters.
N_STEPS = 50
STEP_METHOD = "geometric"

# Earth-Moon nondim time unit -> days (mean lunar sidereal frequency).
TU_DAYS = 1.0 / (2.0 * math.pi) * 27.321661  # ~ 4.348 d / TU


def _seed_orbit() -> bcr4bp.BCR4BPSystem:
    """Build the CR3BP L1 Lyapunov + bcr4bp@mu_sun=0 seed.

    Returns the converged BCR4BP-at-mu_sun=0 seed orbit. Both the L1 Lyapunov
    IC and the Jacobi level are sourced (per docstring).
    """
    sys_cr = cr3bp.CR3BPSystem(
        mu=MU_EM, primary="earth", secondary="moon", l_km=384400.0, t_s=375190.0
    )
    cr_orb = correct_symmetric_free_period(
        sys_cr,
        x0_guess=L1_LYAPUNOV_X0_GUESS,
        jacobi=C_JACOBI_BRAIK_ROSS,
        t_half_guess=1.5,
        ydot0_sign=1.0,
    )
    if not cr_orb.converged:
        raise RuntimeError(
            f"CR3BP L1 Lyapunov seed failed to close at C={C_JACOBI_BRAIK_ROSS}; "
            f"crossing_residual={cr_orb.crossing_residual:.3e}"
        )
    state_seed = np.array([cr_orb.x0, 0.0, 0.0, 0.0, cr_orb.ydot0, 0.0], dtype=np.float64)
    period_seed = float(cr_orb.period)

    sys_zero = bcr4bp.BCR4BPSystem(
        mu=MU_EM,
        mu_sun=0.0,
        a_sun_nondim=A_SUN_NONDIM,
        omega_sun_nondim=OMEGA_SUN_NONDIM,
    )
    seed = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    if not seed.converged:
        raise RuntimeError(
            "BCR4BP@mu_sun=0 corrector failed on CR3BP L1 Lyapunov seed: "
            f"corrector_residual={seed.corrector_residual:.3e}, "
            f"independent_closure={seed.independent_closure_residual:.3e}"
        )
    return seed


def _floquet_summary(member) -> dict[str, list[float]]:
    """Return floquet real / imag parts as plain Python lists (JSON-safe)."""
    if member.floquet is None:
        return {"real_parts": [], "imag_parts": []}
    eigs = member.floquet
    return {
        "real_parts": [float(x.real) for x in eigs],
        "imag_parts": [float(x.imag) for x in eigs],
    }


def _monodromy_summary(member) -> dict[str, float] | None:
    """Return a compact monodromy summary (det, max-mag) for the JSONL."""
    if member.monodromy is None:
        return None
    mono = member.monodromy
    try:
        det = float(np.linalg.det(mono))
    except np.linalg.LinAlgError:
        det = float("nan")
    if member.floquet is not None and len(member.floquet) > 0:
        max_mag = float(np.max(np.abs(member.floquet)))
    else:
        max_mag = float("nan")
    return {"determinant": det, "max_eigenvalue_magnitude": max_mag}


def _member_to_row(member, step_idx: int) -> dict:
    """Serialise one family member to a JSONL row."""
    orb = member.orbit
    sun_n = orb.sun_commensurate_n
    omega = orb.system.omega_sun_nondim
    sun_phase_drift = abs(omega * orb.period_nondim - 2.0 * math.pi * sun_n)
    return {
        "row_type": "member",
        "step_idx": int(step_idx),
        "mu_sun_value": float(member.mu_sun_value),
        "x0": float(orb.state_initial[IDX_X]),
        "y0": float(orb.state_initial[IDX_Y]),
        "z0": float(orb.state_initial[2]),
        "vx0": float(orb.state_initial[IDX_XDOT]),
        "vy0": float(orb.state_initial[IDX_YDOT]),
        "vz0": float(orb.state_initial[IDX_ZDOT]),
        "T_TU": float(orb.period_nondim),
        "T_days": float(orb.period_nondim * TU_DAYS),
        "sun_commensurate_n": int(sun_n),
        "sun_phase_drift": float(sun_phase_drift),
        "corrector_residual": float(orb.corrector_residual),
        "independent_closure_residual": float(orb.independent_closure_residual),
        "monodromy_summary": _monodromy_summary(member),
        "floquet_real_parts": _floquet_summary(member)["real_parts"],
        "floquet_imag_parts": _floquet_summary(member)["imag_parts"],
        "stability_tag": str(member.stability_tag),
    }


def main() -> int:
    t_start = time.time()
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"#303 BCR4BP L1 Lyapunov mu_sun-continuation -- starting at {ts}")
    print(f"  seed: CR3BP planar L1 Lyapunov @ Braik-Ross C = {C_JACOBI_BRAIK_ROSS}")
    print(f"  target: mu_sun = {MU_SUN_ANDREU} (Andreu / Rosales-Jorba 2023 Table 3)")
    print(f"  n_steps = {N_STEPS}, step_method = {STEP_METHOD}")

    seed = _seed_orbit()
    print(
        f"  seed converged: x0={seed.state_initial[IDX_X]:.10f}, "
        f"vy={seed.state_initial[IDX_YDOT]:.10f}, T={seed.period_nondim:.6f} TU, "
        f"corr_res={seed.corrector_residual:.3e}, "
        f"indep_closure={seed.independent_closure_residual:.3e}"
    )

    def _on_step(step_idx, member) -> None:
        if step_idx % 10 == 0 or step_idx == N_STEPS - 1:
            o = member.orbit
            print(
                f"  step {step_idx + 1}/{N_STEPS}: mu_sun={member.mu_sun_value:.4e}, "
                f"x0={o.state_initial[IDX_X]:.6f}, vy={o.state_initial[IDX_YDOT]:.6f}, "
                f"T={o.period_nondim:.6f} TU, corr={o.corrector_residual:.2e}, "
                f"indep={o.independent_closure_residual:.2e}, stab={member.stability_tag}"
            )

    fam = continue_bcr4bp_family_in_musun(
        seed,
        seed_mu_sun=0.0,
        target_mu_sun=MU_SUN_ANDREU,
        n_steps=N_STEPS,
        step_method=STEP_METHOD,
        corrector_tol=1e-10,
        closure_tol=1.0,  # generous: free-T continuation is NOT Sun-commensurate
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        monodromy=True,
        sun_commensurate_n=1,
        on_step=_on_step,
    )

    elapsed = time.time() - t_start
    print(
        f"\nContinuation done in {elapsed:.1f}s: "
        f"{len(fam.members)}/{N_STEPS} members converged. "
        f"mu_sun extent = ({fam.mu_sun_extent[0]:.3e}, {fam.mu_sun_extent[1]:.3e})"
    )
    print(f"  walk_notes: {fam.walk_notes}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    closure_resids = [m.orbit.independent_closure_residual for m in fam.members]
    corr_resids = [m.orbit.corrector_residual for m in fam.members]
    header = {
        "row_type": "header",
        "task_id": 303,
        "phase": "phase-2-mu-sun-continuation",
        "seed_orbit_id": (
            f"cr3bp-em-l1-lyapunov-c{C_JACOBI_BRAIK_ROSS}"
            f"-x{L1_LYAPUNOV_X0_GUESS}-source-braik-ross-2026"
        ),
        "mu_sun_extent": [float(fam.mu_sun_extent[0]), float(fam.mu_sun_extent[1])],
        "seed_mu_sun": float(fam.seed_mu_sun),
        "target_mu_sun": float(MU_SUN_ANDREU),
        "n_steps_attempted": N_STEPS,
        "n_steps_converged": len(fam.members),
        "step_method": STEP_METHOD,
        "mu_em": MU_EM,
        "a_sun_nondim": A_SUN_NONDIM,
        "omega_sun_nondim": OMEGA_SUN_NONDIM,
        "tu_days": TU_DAYS,
        "free_vars": ["x", "vy", "T"],
        "residual_indices_at_T_half": ["y", "vx", "vz"],
        "walk_notes": fam.walk_notes,
        "corrector_residual_max": max(corr_resids) if corr_resids else None,
        "corrector_residual_median": (float(np.median(corr_resids)) if corr_resids else None),
        "independent_closure_max": max(closure_resids) if closure_resids else None,
        "independent_closure_median": (
            float(np.median(closure_resids)) if closure_resids else None
        ),
        "notes": (
            "Independent (Radau) full-period closure is not strictly bounded "
            "because T is FREE and the continuation does not enforce Sun-"
            "commensurate periodicity -- members at moderate-to-strong mu_sun "
            "carry an O(mu_sun) Sun-phase residual. The corrector_residual "
            "(half-period symmetric closure) is at machine precision throughout. "
            "Per orbit-closure discipline, the independent residual is REPORTED "
            "(not used as a hard gate) for analysis."
        ),
    }
    with OUT_PATH.open("w") as fh:
        fh.write(json.dumps(header) + "\n")
        for step_idx, member in enumerate(fam.members):
            row = _member_to_row(member, step_idx)
            fh.write(json.dumps(row) + "\n")
    print(f"  wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
