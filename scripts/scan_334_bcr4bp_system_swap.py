"""BCR4BP system-swap parameter sweep (#334 Phase 4).

Extends the #303/#304/#313 BCR4BP L1-Lyapunov mu_sun-continuation to the rest
of the Sun-primary-secondary triples in
``cyclerfinder.genome.bcr4bp_systems.REGISTRY`` and fits the geometric scaling
rule that #326 (commit ``c1896ef``) predicts from the SEM vs Sun-Jupiter
contrast:

    Δx0_target / Δx0_SEM  ≈  (mu_sun_target / mu_sun_SEM)
                              x (a_sun_SEM / a_sun_target) ** k

with ``k`` in ``[2, 3]``. The sweep records Δx0 (state-component drift from
the CR3BP-limit anchor to the full-mu_sun family endpoint), Δvy0, corrector
residual, independent (Radau) closure residual, and convergence rate per
system; then fits ``k`` from the SEM-anchored two-parameter geometric model.

For each system:

  1. Build the CR3BP from the registry GM / SMA values.
  2. Correct a planar L1 Lyapunov at C = C_L1 - 5e-4 via
     :func:`cyclerfinder.search.reachable_representatives.correct_symmetric_free_period`.
     (The seed grid is widened slightly vs scan_313 to cope with the broader
     mu range -- Pluto-Charon's mu ~ 0.108 is far outside the small-mu regime
     scan_313 was tuned for.)
  3. Build a BCR4BPSystem at mu_sun = 0 with the registry's a_sun / omega_sun
     and re-converge the CR3BP seed via
     :func:`cyclerfinder.genome.bcr4bp_genome.correct_bcr4bp_periodic`. This is
     the CR3BP-LIMIT ANCHOR (structural-correctness test for the registry
     constants -- per ``feedback_orbit_closure_discipline`` every system's
     anchor must pass at corrector precision).
  4. Continue the family in mu_sun from 0 to the registry mu_sun via
     :func:`cyclerfinder.genome.bcr4bp_continuation.continue_bcr4bp_family_in_musun`
     with 51 steps (1 attempted continuation step = #303 N_STEPS+1 includes
     the anchor; the per-system summary echoes ``n_steps_converged``).
  5. Compute the Δx0 / Δvy0 / ΔT signature: (final - anchor) in the converged
     state. This is the "Sun-perturbation displacement" measured by #326.

Then a SEM-anchored geometric-scaling fit:

    log(Δx0_target / Δx0_SEM) = log(mu_sun_target / mu_sun_SEM)
                                + k * log(a_sun_SEM / a_sun_target)

solved by ordinary least squares for ``k``. The fit reports per-system
predicted-vs-observed residuals so the caller can see which systems sit on
the rule and which deviate (suggesting a different dynamical regime).

Discipline
----------
  * Sourced constants only: every per-system mu / mu_sun / a_sun_nondim /
    omega_sun_nondim is a sourced/derived registry record.
  * Independent (Radau) cross-check on every CR3BP-limit anchor is enforced
    by the Phase 1 corrector itself.
  * NO catalogue writeback. JSONL + doc deliverables only.
  * NO novelty claims. This is a parameter-sweep / regime-map; nothing new
    is being asserted about the published BCR4BP literature.

Outputs
-------
  ``data/scan_334_bcr4bp_system_swap.jsonl``: one header record + one record
  per system + a final ``"row_type": "scaling_fit"`` record with the fitted k.

Usage
-----
``uv run python scripts/scan_334_bcr4bp_system_swap.py``
"""

from __future__ import annotations

import dataclasses
import json
import math
import subprocess
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_continuation import continue_bcr4bp_family_in_musun
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    correct_bcr4bp_periodic,
)
from cyclerfinder.genome.bcr4bp_systems import (
    REGISTRY,
    SEM_ANDREU,
    BCR4BPSystemConstants,
)
from cyclerfinder.search.reachable_representatives import (
    correct_symmetric_free_period,
    lagrange_collinear_x,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "data" / "scan_334_bcr4bp_system_swap.jsonl"
N_STEPS = 50  # 51 total members including the anchor; matches #303
# Seed-grid order: scan_313's (x_offset, t_half) grid at c_offset=5e-4 is the
# PRIMARY attempt, so SEM/SJE/SJI/Titan/Pluto-Charon land on the same family
# member as the scan_313 regression. The smaller c_offsets are tried ONLY as
# fallback for systems where 5e-4 overshoots the L1 Lyapunov family (the
# narrow-C regime at very small mu -- Phobos, Enceladus).
SEED_C_OFFSETS: Sequence[float] = (5e-4, 5e-5, 5e-6, 5e-7, 1e-7)
SEED_X_OFFSETS: Sequence[float] = (1e-4, 5e-4, 1e-3, 5e-3)
SEED_T_HALFS: Sequence[float] = (1.0, 1.5, 2.0)


# ---------------------------------------------------------------------------
# Result records.
# ---------------------------------------------------------------------------


@dataclass
class PerSystemResult:
    """Per-system Δx0 / convergence summary."""

    system_name: str
    primary: str
    secondary: str
    mu: float
    mu_sun: float
    a_sun_nondim: float
    omega_sun_nondim: float
    # CR3BP-limit anchor.
    anchor_converged: bool
    anchor_x0: float | None
    anchor_vy0: float | None
    anchor_period_tu: float | None
    anchor_corrector_residual: float | None
    anchor_independent_closure: float | None
    # mu_sun continuation outcome.
    n_steps_attempted: int
    n_steps_converged: int
    mu_sun_extent: tuple[float, float]
    final_mu_sun: float
    final_x0: float | None
    final_vy0: float | None
    final_period_tu: float | None
    final_corrector_residual: float | None
    final_independent_closure: float | None
    # Sun-perturbation signature.
    delta_x0: float | None
    delta_vy0: float | None
    delta_period_tu: float | None
    # Phase status.
    phase: str  # "seed_failed" | "anchor_failed" | "continuation_complete"
    walk_notes: str
    verdict: str

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def _seed_l1_lyapunov_cr3bp(
    system: cr3bp.CR3BPSystem,
) -> tuple[np.ndarray, float, float, float] | tuple[None, None, None, None]:
    """Correct a planar L1 Lyapunov via a small (c_offset, x_offset, t_half) grid.

    Returns ``(state6, period_nondim, jacobi_C, x_L1)`` on success, or
    ``(None, None, None, None)`` if the seed grid fails to find a converged
    member in the L1 Lyapunov basin.

    The grid scans several c_offsets so the same routine handles both the
    mu ~ 1e-2 regime (SEM, Pluto-Charon) and the mu ~ 1e-8 regime (Phobos,
    Enceladus) -- at very small mu the L1 Lyapunov family is very narrow in
    Jacobi space and a standard 5e-4 offset overshoots the family.
    """
    try:
        x_l1 = lagrange_collinear_x(system.mu, "L1")
    except Exception:
        return None, None, None, None
    c_l1 = cr3bp.jacobi_constant(np.array([x_l1, 0.0, 0.0, 0.0, 0.0, 0.0]), system.mu)

    # L1 Lyapunov basin acceptance window: x within 0.01 of x_L1 by default
    # (scan_313 setting). For large mu (e.g. Pluto-Charon 0.108) the basin
    # shifts substantially; widen proportionally to mu but never tighter than
    # 0.01 so small-mu seeds keep the strict scan_313 acceptance.
    x_window = max(0.01, system.mu)
    for c_offset in SEED_C_OFFSETS:
        jacobi = c_l1 - c_offset
        for x_offset in SEED_X_OFFSETS:
            for t_half in SEED_T_HALFS:
                try:
                    orb = correct_symmetric_free_period(
                        system,
                        x0_guess=x_l1 - x_offset,
                        jacobi=jacobi,
                        t_half_guess=t_half,
                        ydot0_sign=1.0,
                        tol=1e-11,
                        max_iter=100,
                    )
                except RuntimeError:
                    continue
                if not orb.converged:
                    continue
                if abs(orb.x0 - x_l1) <= x_window and 1.0 < orb.period < 6.0:
                    state = np.array([orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0], dtype=np.float64)
                    return state, float(orb.period), float(jacobi), float(x_l1)
    return None, None, None, None


def _sweep_one_system(consts: BCR4BPSystemConstants) -> PerSystemResult:
    """Run the L1 Lyapunov seed -> anchor -> mu_sun continuation for one entry."""
    name = consts.name
    print(f"\n[scan-334] === {name} ===", flush=True)
    print(
        f"[scan-334] {name}: mu={consts.mu:.4e} mu_sun={consts.mu_sun:.4e} "
        f"a_sun={consts.a_sun_nondim:.4f} omega_sun={consts.omega_sun_nondim:.6f}",
        flush=True,
    )

    cr3bp_sys = cr3bp.CR3BPSystem(
        mu=consts.mu,
        primary=consts.primary,
        secondary=consts.secondary,
        l_km=consts.l_km,
        t_s=consts.tu_seconds,
    )

    result = PerSystemResult(
        system_name=name,
        primary=consts.primary,
        secondary=consts.secondary,
        mu=consts.mu,
        mu_sun=consts.mu_sun,
        a_sun_nondim=consts.a_sun_nondim,
        omega_sun_nondim=consts.omega_sun_nondim,
        anchor_converged=False,
        anchor_x0=None,
        anchor_vy0=None,
        anchor_period_tu=None,
        anchor_corrector_residual=None,
        anchor_independent_closure=None,
        n_steps_attempted=0,
        n_steps_converged=0,
        mu_sun_extent=(0.0, 0.0),
        final_mu_sun=0.0,
        final_x0=None,
        final_vy0=None,
        final_period_tu=None,
        final_corrector_residual=None,
        final_independent_closure=None,
        delta_x0=None,
        delta_vy0=None,
        delta_period_tu=None,
        phase="seed_failed",
        walk_notes="",
        verdict="seed_failed",
    )

    # 1) CR3BP seed.
    print(f"[scan-334] {name}: seeding CR3BP L1 Lyapunov", flush=True)
    state_cr, period_cr, jacobi_cr, x_l1 = _seed_l1_lyapunov_cr3bp(cr3bp_sys)
    if state_cr is None:
        print(
            f"[scan-334] {name}: L1 LYAPUNOV CR3BP SEED FAILED -- skipping system.",
            flush=True,
        )
        result.walk_notes = "cr3bp seed grid exhausted without converging in L1 basin"
        result.verdict = "seed_failed"
        return result

    print(
        f"[scan-334] {name}: CR3BP seed converged "
        f"x0={state_cr[0]:.6f} vy={state_cr[4]:.6f} T={period_cr:.4f} TU "
        f"C={jacobi_cr:.6f} x_L1={x_l1:.6f}",
        flush=True,
    )

    # 2) BCR4BP @ mu_sun=0 anchor.
    print(f"[scan-334] {name}: BCR4BP@mu_sun=0 anchor", flush=True)
    sys_zero = bcr4bp.BCR4BPSystem(
        mu=consts.mu,
        mu_sun=0.0,
        a_sun_nondim=consts.a_sun_nondim,
        omega_sun_nondim=consts.omega_sun_nondim,
    )
    anchor = correct_bcr4bp_periodic(
        sys_zero,
        state_cr,
        period_cr,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    result.anchor_converged = bool(anchor.converged)
    result.anchor_x0 = float(anchor.state_initial[IDX_X])
    result.anchor_vy0 = float(anchor.state_initial[IDX_YDOT])
    result.anchor_period_tu = float(anchor.period_nondim)
    result.anchor_corrector_residual = float(anchor.corrector_residual)
    result.anchor_independent_closure = float(anchor.independent_closure_residual)
    print(
        f"[scan-334] {name}: anchor converged={anchor.converged} "
        f"corr_res={anchor.corrector_residual:.3e} "
        f"indep_closure={anchor.independent_closure_residual:.3e}",
        flush=True,
    )
    if not anchor.converged:
        # Anchor failed => registry constants inconsistent OR L1 Lyapunov basin
        # not reachable from this seed. Per the orbit-closure discipline this
        # is a STRUCTURAL FAILURE, surfaced explicitly rather than papered over.
        print(
            f"[scan-334] {name}: CR3BP-LIMIT ANCHOR FAILED -- aborting "
            "Sun-perturbation continuation for this system.",
            flush=True,
        )
        result.phase = "anchor_failed"
        result.walk_notes = (
            f"BCR4BP@mu_sun=0 corrector failed at CR3BP seed: "
            f"corrector_residual={anchor.corrector_residual:.3e}, "
            f"independent_closure={anchor.independent_closure_residual:.3e}"
        )
        result.verdict = "anchor_failed"
        return result

    # 3) mu_sun continuation: 0 -> target.
    target_mu_sun = consts.mu_sun
    print(
        f"[scan-334] {name}: continuing in mu_sun from 0 to {target_mu_sun:.4e} "
        f"with n_steps={N_STEPS}",
        flush=True,
    )

    def _on_step(step_idx: int, member) -> None:
        if step_idx % 10 == 0 or step_idx == N_STEPS - 1:
            o = member.orbit
            print(
                f"[scan-334] {name} step {step_idx + 1}/{N_STEPS}: "
                f"mu_sun={member.mu_sun_value:.4e} "
                f"x0={o.state_initial[0]:.6f} vy={o.state_initial[4]:.6f} "
                f"T={o.period_nondim:.4f} TU "
                f"corr={o.corrector_residual:.2e} "
                f"indep={o.independent_closure_residual:.2e} "
                f"stab={member.stability_tag}",
                flush=True,
            )

    t_start = time.time()
    family = continue_bcr4bp_family_in_musun(
        anchor,
        seed_mu_sun=0.0,
        target_mu_sun=target_mu_sun,
        n_steps=N_STEPS,
        step_method="geometric",
        corrector_tol=1e-10,
        # Generous closure_tol: free-T continuation is NOT Sun-commensurate;
        # the independent (Radau) closure is REPORTED, not enforced. Matches
        # #303's pattern.
        closure_tol=1.0,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        monodromy=False,  # not needed for the Δx0 measurement
        sun_commensurate_n=1,
        on_step=_on_step,
    )
    elapsed = time.time() - t_start

    result.n_steps_attempted = N_STEPS
    result.n_steps_converged = len(family.members)
    result.mu_sun_extent = (float(family.mu_sun_extent[0]), float(family.mu_sun_extent[1]))
    result.walk_notes = family.walk_notes
    result.phase = "continuation_complete"

    if family.members:
        last = family.members[-1].orbit
        result.final_mu_sun = float(family.members[-1].mu_sun_value)
        result.final_x0 = float(last.state_initial[IDX_X])
        result.final_vy0 = float(last.state_initial[IDX_YDOT])
        result.final_period_tu = float(last.period_nondim)
        result.final_corrector_residual = float(last.corrector_residual)
        result.final_independent_closure = float(last.independent_closure_residual)
        result.delta_x0 = float(last.state_initial[IDX_X] - anchor.state_initial[IDX_X])
        result.delta_vy0 = float(last.state_initial[IDX_YDOT] - anchor.state_initial[IDX_YDOT])
        result.delta_period_tu = float(last.period_nondim - anchor.period_nondim)
        survival_frac = result.final_mu_sun / target_mu_sun if target_mu_sun > 0 else 0.0
        if survival_frac >= 0.999:
            result.verdict = "family_survives_full_sun_perturbation"
        elif survival_frac >= 0.5:
            result.verdict = "family_partial_survival"
        else:
            result.verdict = "family_breaks_early"
    else:
        result.verdict = "no_members_converged"

    print(
        f"[scan-334] {name}: continuation done in {elapsed:.1f}s; "
        f"{len(family.members)}/{N_STEPS} members converged; "
        f"delta_x0={result.delta_x0}; verdict={result.verdict}",
        flush=True,
    )
    return result


# ---------------------------------------------------------------------------
# Geometric-scaling fit.
# ---------------------------------------------------------------------------


def _fit_geometric_scaling(
    results: Sequence[PerSystemResult],
    sem_consts: BCR4BPSystemConstants,
) -> dict:
    """SEM-anchored OLS fit of the (mu_sun, a_sun) geometric scaling rule.

    Model:

        log(|Δx0_target| / |Δx0_SEM|)  =  log(mu_sun_target / mu_sun_SEM)
                                          + k * log(a_sun_SEM / a_sun_target)

    Rearranged:

        y_i  =  log(|Δx0_i| / |Δx0_SEM|) - log(mu_sun_i / mu_sun_SEM)
             =  k * log(a_sun_SEM / a_sun_i)
             =  k * X_i

    OLS estimator:

        k_hat  =  Σ X_i y_i / Σ X_i^2

    SEM itself contributes ``X = 0`` and ``y = 0`` so it has no leverage on
    the fit (the SEM point is the anchor, by construction).

    Returns a dict with the fitted k, per-system observed-vs-predicted ratios,
    and a residual diagnostic. Only systems with a valid ``delta_x0`` AND
    ``phase == "continuation_complete"`` AND ``|delta_x0| > 0`` are included
    in the fit (silently-zero deltas at high a_sun would log-fail).
    """
    sem_match = [r for r in results if r.system_name == sem_consts.name and r.delta_x0 is not None]
    if not sem_match:
        return {
            "fitted_k": None,
            "k_se": None,
            "n_systems_in_fit": 0,
            "per_system": [],
            "notes": "SEM anchor missing or did not produce a delta_x0; cannot fit.",
        }
    sem = sem_match[0]
    if sem.delta_x0 is None or sem.delta_x0 == 0.0:
        return {
            "fitted_k": None,
            "k_se": None,
            "n_systems_in_fit": 0,
            "per_system": [],
            "notes": "SEM delta_x0 is zero or None; cannot anchor the geometric fit.",
        }

    sem_delta_x0 = abs(sem.delta_x0)
    sem_mu_sun = sem.mu_sun
    sem_a_sun = sem.a_sun_nondim

    xs: list[float] = []
    ys: list[float] = []
    per_system: list[dict] = []
    for r in results:
        # Exclude systems whose continuation broke before reaching the
        # registry mu_sun -- the recorded Δx0 reflects the break-point, not
        # the full-mu_sun perturbation the scaling rule predicts.
        survived_full = r.final_mu_sun >= 0.999 * r.mu_sun if r.mu_sun > 0 else False
        if (
            r.phase != "continuation_complete"
            or r.delta_x0 is None
            or r.delta_x0 == 0.0
            or r.mu_sun <= 0.0
            or r.a_sun_nondim <= 0.0
            or not survived_full
        ):
            per_system.append(
                {
                    "system_name": r.system_name,
                    "included_in_fit": False,
                    "reason": (
                        f"phase={r.phase}; delta_x0={r.delta_x0}; "
                        f"mu_sun={r.mu_sun}; a_sun={r.a_sun_nondim}; "
                        f"final_mu_sun={r.final_mu_sun}; "
                        f"survival_frac={(r.final_mu_sun / r.mu_sun) if r.mu_sun > 0 else 0.0:.3f}"
                    ),
                }
            )
            continue
        if r.system_name == sem.system_name:
            per_system.append(
                {
                    "system_name": r.system_name,
                    "included_in_fit": True,
                    "is_anchor": True,
                    "delta_x0": r.delta_x0,
                    "log_ratio_x0": 0.0,
                    "log_mu_sun_ratio": 0.0,
                    "log_a_sun_ratio": 0.0,
                }
            )
            continue
        ratio_x0 = abs(r.delta_x0) / sem_delta_x0
        mu_sun_ratio = r.mu_sun / sem_mu_sun
        a_sun_ratio = sem_a_sun / r.a_sun_nondim
        x_i = math.log(a_sun_ratio)
        y_i = math.log(ratio_x0) - math.log(mu_sun_ratio)
        xs.append(x_i)
        ys.append(y_i)
        per_system.append(
            {
                "system_name": r.system_name,
                "included_in_fit": True,
                "is_anchor": False,
                "delta_x0": r.delta_x0,
                "log_ratio_x0": math.log(ratio_x0),
                "log_mu_sun_ratio": math.log(mu_sun_ratio),
                "log_a_sun_ratio": x_i,
                "ols_residual_term": y_i,
            }
        )

    if not xs:
        return {
            "fitted_k": None,
            "k_se": None,
            "n_systems_in_fit": 1,  # SEM only
            "per_system": per_system,
            "notes": "only SEM available; need >=1 non-SEM system for OLS.",
        }

    xs_arr = np.asarray(xs, dtype=np.float64)
    ys_arr = np.asarray(ys, dtype=np.float64)
    denom = float(np.dot(xs_arr, xs_arr))
    k_hat = float(np.dot(xs_arr, ys_arr) / denom) if denom > 0 else float("nan")
    # OLS standard error on k (no-intercept simple regression).
    resids = ys_arr - k_hat * xs_arr
    n = len(xs_arr)
    if n > 1 and denom > 0:
        sigma2 = float(np.dot(resids, resids) / max(n - 1, 1))
        k_se = float(math.sqrt(sigma2 / denom))
    else:
        k_se = float("nan")

    # Per-system predicted-vs-observed delta_x0.
    for entry in per_system:
        if not entry.get("included_in_fit") or entry.get("is_anchor"):
            continue
        log_pred = entry["log_mu_sun_ratio"] + k_hat * entry["log_a_sun_ratio"]
        entry["predicted_log_ratio_x0"] = log_pred
        entry["observed_minus_predicted"] = entry["log_ratio_x0"] - log_pred
        entry["predicted_delta_x0"] = sem_delta_x0 * math.exp(log_pred)

    return {
        "fitted_k": k_hat,
        "k_se": k_se,
        "n_systems_in_fit": n + 1,  # +1 for SEM anchor
        "sem_anchor": {
            "system_name": sem.system_name,
            "delta_x0": sem.delta_x0,
            "mu_sun": sem.mu_sun,
            "a_sun_nondim": sem.a_sun_nondim,
        },
        "per_system": per_system,
        "notes": (
            "OLS k = sum(X_i * y_i) / sum(X_i^2) where "
            "X_i = log(a_sun_SEM / a_sun_i) and "
            "y_i = log(|dx_i|/|dx_SEM|) - log(mu_sun_i/mu_sun_SEM). "
            "SEM is the anchor; its (X, y) = (0, 0). "
            "Standard error uses n-1 df."
        ),
    }


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


def main() -> int:
    t_start = time.time()
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    sha = _git_sha()
    print(f"#334 BCR4BP system-swap sweep -- starting at {ts}, sha={sha}", flush=True)
    print(f"  registry: {len(REGISTRY)} systems", flush=True)

    results: list[PerSystemResult] = []
    for consts in REGISTRY:
        try:
            r = _sweep_one_system(consts)
        except Exception as exc:
            print(
                f"[scan-334] {consts.name}: UNHANDLED EXCEPTION {type(exc).__name__}: {exc}",
                flush=True,
            )
            r = PerSystemResult(
                system_name=consts.name,
                primary=consts.primary,
                secondary=consts.secondary,
                mu=consts.mu,
                mu_sun=consts.mu_sun,
                a_sun_nondim=consts.a_sun_nondim,
                omega_sun_nondim=consts.omega_sun_nondim,
                anchor_converged=False,
                anchor_x0=None,
                anchor_vy0=None,
                anchor_period_tu=None,
                anchor_corrector_residual=None,
                anchor_independent_closure=None,
                n_steps_attempted=0,
                n_steps_converged=0,
                mu_sun_extent=(0.0, 0.0),
                final_mu_sun=0.0,
                final_x0=None,
                final_vy0=None,
                final_period_tu=None,
                final_corrector_residual=None,
                final_independent_closure=None,
                delta_x0=None,
                delta_vy0=None,
                delta_period_tu=None,
                phase="exception",
                walk_notes=f"{type(exc).__name__}: {exc}",
                verdict="exception",
            )
        results.append(r)

    fit = _fit_geometric_scaling(results, SEM_ANDREU)

    elapsed = time.time() - t_start
    print(
        f"\n[scan-334] DONE in {elapsed:.1f}s; {len(results)} systems characterised; "
        f"fitted k = {fit['fitted_k']}, k_se = {fit.get('k_se')}",
        flush=True,
    )

    # Write output.
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "row_type": "header",
        "task_id": 334,
        "phase": "phase-4-bcr4bp-system-swap",
        "git_sha": sha,
        "ts_started": ts,
        "n_systems": len(REGISTRY),
        "n_steps_per_system": N_STEPS,
        "step_method": "geometric",
        "registry_source": "cyclerfinder.genome.bcr4bp_systems.REGISTRY",
        "sem_anchor": SEM_ANDREU.name,
        "notes": (
            "Per-system Δx0 / Δvy0 = (final BCR4BP-family member at full mu_sun) "
            "- (CR3BP-limit anchor at mu_sun=0). SEM anchor delta_x0 ~ 1.055e-4 "
            "per #326; non-SEM systems span ~5 orders of magnitude. Geometric "
            "scaling fit anchored on SEM."
        ),
    }
    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps(header, separators=(",", ":"), sort_keys=True) + "\n")
        for r in results:
            row = {"row_type": "system", **r.to_dict()}
            fh.write(json.dumps(row, separators=(",", ":"), sort_keys=True) + "\n")
        fh.write(
            json.dumps(
                {"row_type": "scaling_fit", **fit},
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        )
    print(f"[scan-334] WROTE {OUT_PATH}", flush=True)

    # Console summary.
    print()
    print(
        f"{'system':<55} {'mu_sun':>14} {'a_sun':>14} "
        f"{'|delta_x0|':>14} {'phase':<24} {'verdict':<40}"
    )
    print("-" * 168)
    for r in results:
        dx = "n/a" if r.delta_x0 is None else f"{abs(r.delta_x0):.3e}"
        print(
            f"{r.system_name:<55} {r.mu_sun:>14.4e} {r.a_sun_nondim:>14.4f} "
            f"{dx:>14} {r.phase:<24} {r.verdict:<40}"
        )
    print()
    print(f"fitted geometric scaling k = {fit['fitted_k']}, se = {fit.get('k_se')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
