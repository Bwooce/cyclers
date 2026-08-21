"""#861: mu-continuation cross-check -- Neptune-Triton "4:5-saddle" descended to
Uranus-Oberon's own mass ratio, as an INDEPENDENT identity check on whatever the
fold-turning gate (`scripts/run_861_oberon_gate.py`) recovers (or fails to
recover) for Oberon's own 4:5 family.

`mu_continuation.continue_in_mu` (pre-existing, `#249`'s sibling) starts from a
TABLE-VERIFIED saddle (`neptune_triton_resonant_families.ESM_GATE_ROWS["4:5-saddle"]`,
Miceli & Bosanac 2026 ESM data, C=2.987089791658, `|lambda|~=105` per `#771`'s own
survey) and pseudo-arclength-continues it in mu down from Neptune-Triton's
2.089503183689124e-04 to Oberon's own 3.54326e-5 (~5.9x descent, `#860` Sec. 4(d)'s
own "modest span, safe" assessment) -- the source member's topology is verified by
construction (a real published row), so whatever arrives at Oberon's mu carries a
traceable identity, independent of BOTH `#861`'s own two Oberon-side fixes.

DIRECTIONALITY NOTE (found this task, not a pre-existing bug fixed here -- out of
`#861`'s own scope to change the shared module): every prior use of
`mu_continuation.continue_in_mu` in this project's history continues mu UPWARD
(Earth-Moon 0.0122 -> 0.1/0.3/0.5, Roberts-Tsoukkas reproduction). Its own initial
tangent (`_tangent(..., prev=None, ...)`) unconditionally orients toward
INCREASING mu (the module's own "orient toward increasing mu at the start"
convention) regardless of where the target actually is -- so calling the public
`continue_in_mu` toward a mu_target BELOW mu_start makes the very first step head
the WRONG way (confirmed empirically this task: a naive call walked mu from
2.0895e-4 UP to 2.45e-4, away from Oberon's 3.54e-5). This script works around it
locally (reusing the module's own private stepping helpers with a
DOWNWARD-oriented initial tangent) rather than editing the shared module -- flagged
here for a future task's own registration, not fixed as part of `#861`.

Foreground only, chunked (this project's own standing lesson,
`feedback_subagent_background_is_fatal`).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.mu_continuation as mc
import cyclerfinder.search.neptune_triton_resonant_families as ntrf

OBERON_MU = 3.54326e-5


def continue_in_mu_downward(
    seed: cp.SymmetricOrbit,
    mu_start: float,
    *,
    half_crossings: int,
    ydot0_sign: float,
    mu_target: float,
    ds0: float = 5e-6,
    ds_max: float = 3e-5,
    ds_min: float = 1e-9,
    max_steps: int = 400,
    corrector_tol: float = 1e-10,
    corrector_max_iter: int = 80,
    t_hi_frac: float = 1.8,
    radau_closure_tol: float = 1e-6,
    radau_jacobi_tol: float = 1e-7,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    record_every: int = 5,
    progress: bool = True,
) -> mc.MuBranch:
    """Same walk as ``mu_continuation.continue_in_mu`` but with the INITIAL
    tangent forced toward DECREASING mu (see module docstring's directionality
    note) -- required because ``mu_target < mu_start`` here (Neptune-Triton down
    to Oberon), a direction the shared module's own first-tangent convention was
    never exercised in.
    """
    mu_start = float(mu_start)
    z = np.array([seed.x0, seed.jacobi, mu_start])
    per = seed.period
    t_hi = t_hi_frac * per

    branch = mc.MuBranch(
        label="neptune_triton_4_5_saddle_to_oberon",
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
        mu_start=mu_start,
        mu_target=mu_target,
    )
    first = mc._make_member(
        z,
        ydot0_sign,
        half_crossings,
        t_hi,
        rtol=rtol,
        atol=atol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
    )
    if first is None:
        branch.stop_reason = mc.MuStopReason.NO_MEMBER
        return branch
    branch.members.append(first)

    # Force the INITIAL tangent toward decreasing mu (prev with a negative
    # mu-component -- see module docstring).
    tan = mc._tangent(
        z, ydot0_sign, half_crossings, t_hi, prev=np.array([0.0, 0.0, -1.0]), rtol=rtol, atol=atol
    )
    if tan is None:
        branch.stop_reason = mc.MuStopReason.NO_MEMBER
        return branch

    ds = ds0
    for _step in range(max_steps):
        if z[2] <= mu_target + 1e-14:
            branch.stop_reason = mc.MuStopReason.TARGET_REACHED
            break
        if tan[2] != 0.0:
            ds_to_target = (mu_target - z[2]) / tan[2]
            if 0.0 < ds_to_target <= ds:
                z_pred = z + ds_to_target * tan
                zland = mc._land_at_mu(
                    z_pred,
                    mu_target,
                    ydot0_sign,
                    half_crossings,
                    t_hi_frac * max(per, seed.period),
                    tol=corrector_tol,
                    max_iter=corrector_max_iter,
                    rtol=rtol,
                    atol=atol,
                )
                if zland is not None:
                    z = zland
                    branch.n_steps += 1
                    branch.stop_reason = mc.MuStopReason.TARGET_REACHED
                    break
        z_pred = z + ds * tan
        t_hi = t_hi_frac * max(per, seed.period)
        zc = mc._correct(
            z_pred,
            tan,
            ydot0_sign,
            half_crossings,
            t_hi,
            tol=corrector_tol,
            max_iter=corrector_max_iter,
            rtol=rtol,
            atol=atol,
        )
        if zc is None:
            ds *= 0.5
            if ds < ds_min:
                branch.stop_reason = mc.MuStopReason.STEP_UNDERFLOW
                break
            continue
        ntan = mc._tangent(zc, ydot0_sign, half_crossings, t_hi, prev=tan, rtol=rtol, atol=atol)
        if ntan is None:
            ds *= 0.5
            if ds < ds_min:
                branch.stop_reason = mc.MuStopReason.STEP_UNDERFLOW
                break
            continue
        z, tan, per = zc, ntan, per
        branch.n_steps += 1
        ds = min(ds * 1.3, ds_max)
        if branch.n_steps % record_every == 0:
            mem = mc._make_member(
                z,
                ydot0_sign,
                half_crossings,
                t_hi_frac * per,
                rtol=rtol,
                atol=atol,
                radau_closure_tol=radau_closure_tol,
                radau_jacobi_tol=radau_jacobi_tol,
            )
            if mem is not None:
                branch.members.append(mem)
                if progress:
                    print(
                        f"  step {branch.n_steps}: mu={mem.mu:.6e} C={mem.jacobi:.6f} "
                        f"|lambda|={mem.abs_lambda:.4g} stable={mem.stable}",
                        flush=True,
                    )
    else:
        branch.stop_reason = mc.MuStopReason.MAX_STEPS

    final = mc._make_member(
        z,
        ydot0_sign,
        half_crossings,
        t_hi_frac * per,
        rtol=rtol,
        atol=atol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
    )
    if final is not None and (not branch.members or abs(branch.members[-1].mu - final.mu) > 1e-15):
        branch.members.append(final)
        if progress:
            print(
                f"  FINAL: mu={final.mu:.6e} C={final.jacobi:.6f} "
                f"|lambda|={final.abs_lambda:.4g} stable={final.stable}",
                flush=True,
            )
    return branch


def main() -> None:
    out_path = Path("data/found/861_resonant_seeding_oberon_gate/mu_continuation_crosscheck.json")
    system = ntrf.neptune_triton_system()
    row = ntrf.ESM_GATE_ROWS["4:5-saddle"]
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        row.x0,
        row.jacobi,
        row.period,
        ydot0_sign=row.ydot0_sign,
        half_crossings=row.half_crossings,
        tol=1e-11,
        rtol=1e-13,
        atol=1e-13,
        x0_bounds=row.x0_bounds,
    )
    assert orbit.converged, "table-verified 4:5-saddle seed failed to converge -- regression"
    print(f"[861-mu] seed: C={orbit.jacobi} period={orbit.period} mu={ntrf.MICELI_MU}")

    t0 = time.time()
    branch = continue_in_mu_downward(
        orbit,
        ntrf.MICELI_MU,
        half_crossings=row.half_crossings,
        ydot0_sign=row.ydot0_sign,
        mu_target=OBERON_MU,
    )
    elapsed = time.time() - t0
    print(f"[861-mu] done in {elapsed:.1f}s, stop={branch.stop_reason}, ")
    print(f"[861-mu] n_members={len(branch.members)}")

    out = {
        "elapsed_s": elapsed,
        "stop_reason": str(branch.stop_reason),
        "n_steps": branch.n_steps,
        "mu_start": ntrf.MICELI_MU,
        "mu_target": OBERON_MU,
        "members": [
            {
                "mu": m.mu,
                "x0": m.x0,
                "ydot0": m.ydot0,
                "jacobi": m.jacobi,
                "period": m.period,
                "abs_lambda": m.abs_lambda,
                "stable": m.stable,
                "crossing_residual": m.crossing_residual,
                "radau_djacobi": m.radau_djacobi,
            }
            for m in branch.members
        ],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"[861-mu] wrote {out_path}")


if __name__ == "__main__":
    main()
