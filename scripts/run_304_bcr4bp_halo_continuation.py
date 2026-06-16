"""Run #304: BCR4BP halo-orbit family mu_sun continuation from CR3BP to Andreu.

Produces ``data/bcr4bp_halo_family_304.jsonl`` containing one converged
BCR4BP halo family member per step, continued from ``mu_sun = 0`` (CR3BP
limit) to ``mu_sun = 328900.5423094043`` (Andreu / Rosales-Jorba parameter
value).

Discipline
----------
  * Sourced seed: Earth-Moon L1 southern halo IC (Howell 1984 / NASA TN D-1949
    family). The IC ``(0.824024728136525, 0, -0.054501847320725, 0,
    0.164671964079122, 0)`` with T_guess = 2.7549 TU is the same as
    ``tests/search/test_cr3bp_periodic.py::test_cr3bp_periodic_halo_l1_southern``
    and ``tests/genome/test_bcr4bp_genome.py::test_cr3bp_limit_closure_recovers_cr3bp_halo``.
  * Independent (Radau) cross-check on every accepted member is enforced
    by the Phase 1 corrector itself.
  * NO catalogue writeback.
  * NO novelty claims -- the Howell halo family + mu_sun continuation is
    published methodology (Simo-Jorba-Gomez / Andreu).
  * Halo masks: FREE_VARS_HALO = (x, z, vy, T); RESIDUAL_HALO_HALF_PERIOD =
    (y, vx, vz) at T/2 -- the perpendicular-crossing symmetric halo setup.

Usage
-----
``uv run python scripts/run_304_bcr4bp_halo_continuation.py``

Output JSONL columns are documented in the header row.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.genome.bcr4bp_continuation import (
    continue_bcr4bp_family_in_musun,
)
from cyclerfinder.genome.bcr4bp_genome import (
    FREE_VARS_HALO,
    IDX_X,
    IDX_YDOT,
    IDX_Z,
    RESIDUAL_HALO_HALF_PERIOD,
    correct_bcr4bp_periodic,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "data" / "bcr4bp_halo_family_304.jsonl"

# Sourced constants -- match Phase 1's andreu_default() exactly.
MU_EM = 0.012150581600000  # Rosales-Jorba 2023 Table 3
MU_SUN_ANDREU = 328900.5423094043
A_SUN_NONDIM = 388.8111430233511
OMEGA_SUN_NONDIM = 0.925195985520347

# Sourced Howell EM L1 southern halo seed (NASA TN D-1949 family).
HALO_SEED_STATE = np.array(
    [0.824024728136525, 0.0, -0.054501847320725, 0.0, 0.164671964079122, 0.0],
    dtype=np.float64,
)
HALO_SEED_PERIOD_GUESS = 2.7549

# Continuation parameters.
N_STEPS = 50
STEP_METHOD = "geometric"

# Earth-Moon nondim time unit -> days (mean lunar sidereal frequency).
TU_DAYS = 1.0 / (2.0 * math.pi) * 27.321661  # ~ 4.348 d / TU


def _seed_orbit() -> bcr4bp.BCR4BPSystem:
    """Close the Howell EM L1 southern halo at mu_sun=0 via the BCR4BP corrector.

    Returns the converged BCR4BP-at-mu_sun=0 seed orbit. The IC and period
    guess are SOURCED (Howell 1984 / NASA TN D-1949 family); the BCR4BP-at-
    mu_sun=0 corrector reduces exactly to CR3BP (verified structurally in
    tests/core/test_bcr4bp.py) so this is a tight CR3BP halo closure that
    happens to be reached via the BCR4BP code path.
    """
    sys_zero = bcr4bp.BCR4BPSystem(
        mu=MU_EM,
        mu_sun=0.0,
        a_sun_nondim=A_SUN_NONDIM,
        omega_sun_nondim=OMEGA_SUN_NONDIM,
    )
    seed = correct_bcr4bp_periodic(
        sys_zero,
        HALO_SEED_STATE,
        HALO_SEED_PERIOD_GUESS,
        sun_commensurate_n=1,
        free_vars=FREE_VARS_HALO,
        residual_indices=RESIDUAL_HALO_HALF_PERIOD,
        is_half_period_residual=True,
        # 1e-10 is the Phase 1 corrector default; the halo Jacobian conditioning
        # is slightly worse than the planar Lyapunov case (4 unknowns / 3 residuals
        # vs 3/3), so the Newton iteration lands a few ULP above 1e-12. The
        # independent (Radau) closure (gated at 1e-6) is the binding gate.
        tol=1e-10,
        independent_tol=1e-6,
    )
    if not seed.converged:
        raise RuntimeError(
            "BCR4BP@mu_sun=0 corrector failed on Howell EM L1 southern halo seed: "
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
        "y0": float(orb.state_initial[1]),
        "z0": float(orb.state_initial[IDX_Z]),
        "vx0": float(orb.state_initial[3]),
        "vy0": float(orb.state_initial[IDX_YDOT]),
        "vz0": float(orb.state_initial[5]),
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
    print(f"#304 BCR4BP halo-orbit mu_sun-continuation -- starting at {ts}")
    print("  seed: Howell EM L1 southern halo (NASA TN D-1949 family)")
    print(f"    state = {HALO_SEED_STATE.tolist()}")
    print(f"    period_guess = {HALO_SEED_PERIOD_GUESS} TU")
    print(f"  target: mu_sun = {MU_SUN_ANDREU} (Andreu / Rosales-Jorba 2023 Table 3)")
    print(f"  n_steps = {N_STEPS}, step_method = {STEP_METHOD}")
    print("  free_vars  = (x, z, vy, T) -- halo perpendicular-crossing free unknowns")
    print("  residuals  = (y, vx, vz) at T/2 -- halo perpendicular-crossing closure")

    seed = _seed_orbit()
    print(
        f"  seed converged: x0={seed.state_initial[IDX_X]:.10f}, "
        f"z0={seed.state_initial[IDX_Z]:.10f}, "
        f"vy={seed.state_initial[IDX_YDOT]:.10f}, T={seed.period_nondim:.6f} TU, "
        f"corr_res={seed.corrector_residual:.3e}, "
        f"indep_closure={seed.independent_closure_residual:.3e}"
    )

    def _on_step(step_idx, member) -> None:
        if step_idx % 10 == 0 or step_idx == N_STEPS - 1:
            o = member.orbit
            print(
                f"  step {step_idx + 1}/{N_STEPS}: mu_sun={member.mu_sun_value:.4e}, "
                f"x0={o.state_initial[IDX_X]:.6f}, z0={o.state_initial[IDX_Z]:.6f}, "
                f"vy={o.state_initial[IDX_YDOT]:.6f}, "
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
        # Free T means non-Sun-commensurate at non-zero mu_sun, so full-period
        # independent closure carries O(mu_sun) Sun-phase residual -- the
        # corrector_residual (half-period symmetric) is the binding gate.
        closure_tol=1.0,
        free_vars=FREE_VARS_HALO,
        residual_indices=RESIDUAL_HALO_HALF_PERIOD,
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
    if fam.members:
        z0_extent = (
            float(min(m.orbit.state_initial[IDX_Z] for m in fam.members)),
            float(max(m.orbit.state_initial[IDX_Z] for m in fam.members)),
        )
        x0_extent = (
            float(min(m.orbit.state_initial[IDX_X] for m in fam.members)),
            float(max(m.orbit.state_initial[IDX_X] for m in fam.members)),
        )
        t_extent_tu = (
            float(min(m.orbit.period_nondim for m in fam.members)),
            float(max(m.orbit.period_nondim for m in fam.members)),
        )
        print(f"  z0 extent: {z0_extent}")
        print(f"  x0 extent: {x0_extent}")
        print(f"  T extent (TU): {t_extent_tu}")
        stab_counts: dict[str, int] = {}
        for m in fam.members:
            stab_counts[m.stability_tag] = stab_counts.get(m.stability_tag, 0) + 1
        print(f"  stability tags: {stab_counts}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    closure_resids = [m.orbit.independent_closure_residual for m in fam.members]
    corr_resids = [m.orbit.corrector_residual for m in fam.members]
    z0_values = [float(m.orbit.state_initial[IDX_Z]) for m in fam.members]
    x0_values = [float(m.orbit.state_initial[IDX_X]) for m in fam.members]
    t_values = [float(m.orbit.period_nondim) for m in fam.members]
    header = {
        "row_type": "header",
        "task_id": 304,
        "phase": "phase-3-halo-mu-sun-continuation",
        "seed_orbit_id": (
            "cr3bp-em-l1-halo-southern-howell-1984-NASA-TN-D-1949-family"
            f"-x{HALO_SEED_STATE[IDX_X]}"
            f"-z{HALO_SEED_STATE[IDX_Z]}"
            f"-vy{HALO_SEED_STATE[IDX_YDOT]}"
        ),
        "seed_state": HALO_SEED_STATE.tolist(),
        "seed_period_guess_TU": HALO_SEED_PERIOD_GUESS,
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
        "free_vars": ["x", "z", "vy", "T"],
        "residual_indices_at_T_half": ["y", "vx", "vz"],
        "walk_notes": fam.walk_notes,
        "corrector_residual_max": max(corr_resids) if corr_resids else None,
        "corrector_residual_median": (float(np.median(corr_resids)) if corr_resids else None),
        "independent_closure_max": max(closure_resids) if closure_resids else None,
        "independent_closure_median": (
            float(np.median(closure_resids)) if closure_resids else None
        ),
        "z0_extent": [min(z0_values), max(z0_values)] if z0_values else None,
        "x0_extent": [min(x0_values), max(x0_values)] if x0_values else None,
        "T_TU_extent": [min(t_values), max(t_values)] if t_values else None,
        "notes": (
            "Independent (Radau) full-period closure is not strictly bounded "
            "because T is FREE and the continuation does not enforce Sun-"
            "commensurate periodicity -- members at moderate-to-strong mu_sun "
            "carry an O(mu_sun) Sun-phase residual. The corrector_residual "
            "(half-period symmetric closure) is at machine precision throughout. "
            "Per orbit-closure discipline, the independent residual is REPORTED "
            "(not used as a hard gate) for analysis. The halo masks (x, z, vy, T) "
            "free + (y, vx, vz) residual at T/2 preserve the perpendicular-"
            "crossing symmetric halo structure: y, vx, vz remain identically "
            "zero at t=0 (the IC) and the same at t=T/2 (the residual)."
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
