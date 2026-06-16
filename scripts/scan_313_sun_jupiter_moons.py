"""Sun-Jupiter-Europa / Sun-Jupiter-Io BCR4BP L1-Lyapunov scout (#313 Phase 1 Part B).

Repurposes the #292 Phase 1 BCR4BPSystem (which is parameterised by 4 constants
and so can represent any Sun-primary-secondary bicircular system, not just
Sun-Earth-Moon) to test whether the CR3BP L1 Lyapunov substitute family
survives the Sun perturbation at Jupiter's distance for the Europa and Io
moons. Read-only over the BCR4BP modules.

Methodology (mirrors #303 ``run_303_bcr4bp_l1_continuation.py``):

  1. Build the Jupiter-moon CR3BP from the registry GM values.
  2. Correct a planar L1 Lyapunov at a Jacobi level just below C_L1 via
     ``correct_symmetric_free_period`` (same routine as #303 for Earth-Moon).
  3. Build a BCR4BPSystem with mu_sun=0 and the Sun-Jupiter-moon constants
     (mu_sun computed as GM_Sun/GM_Jupiter_sys, a_sun_nondim = Jupiter SMA /
     moon SMA, omega_sun = 1 - n_Jupiter_about_Sun * TU_moon).
  4. Re-converge the CR3BP seed via ``correct_bcr4bp_periodic`` at mu_sun=0 —
     this is the CR3BP-LIMIT ANCHOR (structural-correctness test). If the
     CR3BP seed does not round-trip through the BCR4BP corrector with the same
     IC and period, the BCR4BPSystem constants are wrong.
  5. Continue the family in mu_sun from 0 to the Sun-Jupiter-moon mu_sun via
     ``continue_bcr4bp_family_in_musun``.
  6. Report the converged family extent and the qualitative survival verdict.

Sourced constants:
  * GM_Sun = 1.32712440018e11 km^3/s^2 (IAU 2015 / JPL DE440 — already in
    ``core/constants.py`` MU_SUN_KM3_S2).
  * GM_Jupiter_sys = 1.26686534e8 km^3/s^2 (JPL SSD gm_de440 — already in
    ``core/satellites.py`` PRIMARIES).
  * Jupiter SMA = 5.20288700 AU = 7.7822e8 km (JPL DE440, ``core/constants.py``
    _JUPITER_SMA_AU).
  * AU = 1.49597870700e8 km (IAU 2012 Resolution B2 exact).
  * Europa SMA = 671100 km, GM = 3202.739 km^3/s^2 (JPL SSD, satellites.py).
  * Io SMA = 421800 km, GM = 5959.916 km^3/s^2 (JPL SSD, satellites.py).

Computed BCR4BP constants:

  Sun-Jupiter-Europa:
    mu        = 2.528e-5
    L         = 671100 km
    TU        = 48844.5 s = 0.5653 d
    a_sun_LU  = 1159.80
    mu_sun    = 1047.57
    omega_sun = 0.99918 (Sun synodic period ~6.29 TU ~3.55 d)

  Sun-Jupiter-Io:
    mu        = 4.704e-5
    L         = 421800 km
    TU        = 24338.5 s = 0.2817 d
    a_sun_LU  = 1845.28
    mu_sun    = 1047.57
    omega_sun = 0.99959 (Sun synodic period ~6.29 TU ~1.77 d)

Discipline:
  * NO catalogue writeback. JSONL deliverable + doc.
  * NO novelty claims. The Jupiter-moon L1 Lyapunov + BCR4BP framework are
    published; we are scouting whether the EM-style continuation result
    (family survives Sun perturbation) carries over qualitatively.
  * CR3BP-limit anchor (mu_sun=0 -> CR3BP) is the sourced-golden discipline
    test for our BCR4BPSystem constants. If the round-trip residual blows
    up, the constants are wrong — the script raises rather than continues.
"""

from __future__ import annotations

import json
import math
import subprocess
import time
from dataclasses import asdict, dataclass
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
from cyclerfinder.search.reachable_representatives import (
    correct_symmetric_free_period,
    lagrange_collinear_x,
)

# Sourced constants (see module docstring for provenance).
MU_SUN_KM3_S2 = 1.32712440018e11
GM_JUP_SYS = 1.26686534e8
AU_KM = 1.49597870700e8
A_JUP_KM = 5.20288700 * AU_KM
GM_EUROPA = 3202.739
GM_IO = 5959.916
A_EUROPA = 671100.0
A_IO = 421800.0


@dataclass
class SystemConstants:
    moon_name: str
    mu: float
    l_km: float
    t_s: float
    a_sun_nondim: float
    mu_sun: float
    omega_sun_nondim: float

    def to_dict(self) -> dict:
        return asdict(self)


def _build_system_constants(moon_name: str, gm_moon: float, a_moon: float) -> SystemConstants:
    """Compute the BCR4BP-Sun-Jupiter-moon constants from sourced physical inputs."""
    mu = gm_moon / GM_JUP_SYS
    l_km = a_moon
    t_s = math.sqrt(l_km**3 / GM_JUP_SYS)
    a_sun_nondim = A_JUP_KM / l_km
    mu_sun = MU_SUN_KM3_S2 / GM_JUP_SYS
    # omega_sun in the Jupiter-moon synodic frame: Sun appears to move at
    # 1 - n_Jupiter_about_Sun (in nondim TU units). Same formula as the
    # Sun-Earth-Moon case: 1 - n_Earth_about_Sun_in_TU_EM ~ 1 - 0.0748 ~ 0.9252.
    # At Jupiter the moon orbital frequency is huge vs Jupiter's orbital
    # frequency, so omega_sun is very close to 1 (synodic period ~ 2 pi TU).
    n_jupiter_about_sun = math.sqrt(MU_SUN_KM3_S2 / A_JUP_KM**3)
    omega_sun_nondim = 1.0 - n_jupiter_about_sun * t_s
    return SystemConstants(
        moon_name=moon_name,
        mu=mu,
        l_km=l_km,
        t_s=t_s,
        a_sun_nondim=a_sun_nondim,
        mu_sun=mu_sun,
        omega_sun_nondim=omega_sun_nondim,
    )


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


@dataclass
class FamilyRow:
    """One JSONL row representing one BCR4BP continuation member."""

    moon_name: str
    row_type: str  # "constants" | "cr3bp_seed" | "bcr4bp_anchor" | "member" | "summary"
    step_idx: int | None
    mu_sun_value: float | None
    x0: float | None
    y0: float | None
    z0: float | None
    vx0: float | None
    vy0: float | None
    vz0: float | None
    T_TU: float | None
    T_seconds: float | None
    corrector_residual: float | None
    independent_closure_residual: float | None
    sun_phase_drift: float | None
    stability_tag: str | None
    extras: dict

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self), separators=(",", ":"), sort_keys=True)


def _empty_row(moon_name: str, row_type: str) -> FamilyRow:
    return FamilyRow(
        moon_name=moon_name,
        row_type=row_type,
        step_idx=None,
        mu_sun_value=None,
        x0=None,
        y0=None,
        z0=None,
        vx0=None,
        vy0=None,
        vz0=None,
        T_TU=None,
        T_seconds=None,
        corrector_residual=None,
        independent_closure_residual=None,
        sun_phase_drift=None,
        stability_tag=None,
        extras={},
    )


def _seed_cr3bp_l1_lyapunov(
    system: cr3bp.CR3BPSystem, *, c_offset: float = 5e-4
) -> tuple[np.ndarray, float, float, float]:
    """Correct a planar L1 Lyapunov at C = C_L1 - c_offset.

    Returns (state6, period_nondim, jacobi_C, x_L1). Sweeps a small grid of
    (x_offset, t_half_guess) parameters and accepts the first converged orbit
    that lands in the L1 Lyapunov basin (x0 within 0.01 of x_L1 inside L1 and
    period in [1.5, 4.0] TU, the expected range at small mu where T ~ 2*pi).
    The corrector's basin shape varies with mu; the Earth-Moon defaults from
    #303 don't necessarily transfer to Jupiter-moon mu ~ 1e-5.
    """
    x_l1 = lagrange_collinear_x(system.mu, "L1")
    # Compute C at L1 via the standard pseudo-potential identity (v=0 at L1).
    c_l1 = cr3bp.jacobi_constant(np.array([x_l1, 0.0, 0.0, 0.0, 0.0, 0.0]), system.mu)
    jacobi = c_l1 - c_offset

    # Search a small grid of seed parameters for the L1 Lyapunov basin.
    last_resid = None
    last_orb = None
    for x_offset in (1e-4, 5e-4, 1e-3, 5e-3):
        for t_half_guess in (1.0, 1.5, 2.0):
            try:
                orb = correct_symmetric_free_period(
                    system,
                    x0_guess=x_l1 - x_offset,
                    jacobi=jacobi,
                    t_half_guess=t_half_guess,
                    ydot0_sign=1.0,
                    tol=1e-11,
                    max_iter=80,
                )
            except RuntimeError:
                continue
            if not orb.converged:
                last_resid = orb.crossing_residual
                last_orb = orb
                continue
            # Accept the L1 Lyapunov basin: x0 just inside L1, period 1.5-4 TU.
            if abs(orb.x0 - x_l1) < 0.01 and 1.5 < orb.period < 4.0:
                state = np.array([orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0], dtype=np.float64)
                return state, float(orb.period), float(jacobi), float(x_l1)
    raise RuntimeError(
        f"L1 Lyapunov seed failed at {system.secondary}: x_L1={x_l1:.6f}, "
        f"C_L1={c_l1:.6f}, jacobi={jacobi:.6f}, last_residual="
        f"{last_resid!r}, last_x0="
        f"{last_orb.x0 if last_orb is not None else None!r}, "
        f"last_T={last_orb.period if last_orb is not None else None!r}"
    )


def _scout_one_system(consts: SystemConstants) -> tuple[list[FamilyRow], dict]:
    """Run the full BCR4BP scout at one Sun-Jupiter-moon system."""
    moon = consts.moon_name
    print(f"\n[scan-313-B] === Sun-Jupiter-{moon} ===", flush=True)
    print(f"[scan-313-B] {moon}: constants = {consts.to_dict()}", flush=True)

    cr3bp_sys = cr3bp.CR3BPSystem(
        mu=consts.mu,
        primary="Jupiter",
        secondary=moon,
        l_km=consts.l_km,
        t_s=consts.t_s,
    )
    rows: list[FamilyRow] = []

    # Constants row.
    c_row = _empty_row(moon, "constants")
    c_row.extras = consts.to_dict()
    rows.append(c_row)

    # 1) CR3BP seed.
    print(f"[scan-313-B] {moon}: seeding CR3BP L1 Lyapunov", flush=True)
    state_cr, period_cr, jacobi_cr, x_l1 = _seed_cr3bp_l1_lyapunov(cr3bp_sys)
    print(
        f"[scan-313-B] {moon}: CR3BP seed converged "
        f"x0={state_cr[0]:.6f} vy={state_cr[4]:.6f} T={period_cr:.4f} TU "
        f"C={jacobi_cr:.6f} (C_L1-{jacobi_cr:.6f})",
        flush=True,
    )
    cr_row = _empty_row(moon, "cr3bp_seed")
    cr_row.x0 = float(state_cr[0])
    cr_row.y0 = float(state_cr[1])
    cr_row.z0 = float(state_cr[2])
    cr_row.vx0 = float(state_cr[3])
    cr_row.vy0 = float(state_cr[4])
    cr_row.vz0 = float(state_cr[5])
    cr_row.T_TU = period_cr
    cr_row.T_seconds = period_cr * consts.t_s
    cr_row.extras = {"jacobi_C": jacobi_cr, "x_L1": x_l1}
    rows.append(cr_row)

    # 2) BCR4BP @ mu_sun=0 anchor: the CR3BP seed must round-trip through the
    # BCR4BP corrector at mu_sun=0 with the same IC and period (the
    # CR3BP-limit anchor — structural-correctness test for our constants).
    print(f"[scan-313-B] {moon}: BCR4BP @ mu_sun=0 anchor", flush=True)
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
    a_row = _empty_row(moon, "bcr4bp_anchor")
    a_row.mu_sun_value = 0.0
    a_row.x0 = float(anchor.state_initial[0])
    a_row.y0 = float(anchor.state_initial[1])
    a_row.z0 = float(anchor.state_initial[2])
    a_row.vx0 = float(anchor.state_initial[3])
    a_row.vy0 = float(anchor.state_initial[4])
    a_row.vz0 = float(anchor.state_initial[5])
    a_row.T_TU = float(anchor.period_nondim)
    a_row.T_seconds = float(anchor.period_nondim * consts.t_s)
    a_row.corrector_residual = float(anchor.corrector_residual)
    a_row.independent_closure_residual = float(anchor.independent_closure_residual)
    a_row.extras = {"converged": bool(anchor.converged), "n_iter": int(anchor.n_iter)}
    rows.append(a_row)
    print(
        f"[scan-313-B] {moon}: anchor converged={anchor.converged} "
        f"corr_res={anchor.corrector_residual:.3e} "
        f"indep_closure={anchor.independent_closure_residual:.3e}",
        flush=True,
    )
    if not anchor.converged:
        # CR3BP-limit anchor failed — the BCR4BPSystem constants are wrong
        # for this system. Bail out loudly per the sourced-golden discipline.
        print(
            f"[scan-313-B] {moon}: CR3BP-LIMIT ANCHOR FAILED — aborting "
            "Sun-perturbation continuation.",
            flush=True,
        )
        summary = {
            "moon_name": moon,
            "phase": "anchor_failed",
            "n_members_converged": 0,
            "mu_sun_extent": (0.0, 0.0),
            "verdict": "anchor_failed",
        }
        s_row = _empty_row(moon, "summary")
        s_row.extras = summary
        rows.append(s_row)
        return rows, summary

    # 3) mu_sun continuation. Geometric stepping, generous closure_tol because
    # the free-T continuation lets period drift (Sun-commensurability is not
    # strictly enforced — same as the #303 Earth-Moon run).
    target_mu_sun = consts.mu_sun
    n_steps = 50
    print(
        f"[scan-313-B] {moon}: continuing in mu_sun from 0 to {target_mu_sun:.4f} "
        f"with n_steps={n_steps}",
        flush=True,
    )

    def _on_step(step_idx: int, member) -> None:
        if step_idx % 10 == 0 or step_idx == n_steps - 1:
            o = member.orbit
            print(
                f"[scan-313-B] {moon} step {step_idx + 1}/{n_steps}: "
                f"mu_sun={member.mu_sun_value:.4e} "
                f"x0={o.state_initial[0]:.6f} vy={o.state_initial[4]:.6f} "
                f"T={o.period_nondim:.4f} TU corr={o.corrector_residual:.2e} "
                f"indep={o.independent_closure_residual:.2e} "
                f"stab={member.stability_tag}",
                flush=True,
            )

    t_start = time.time()
    family = continue_bcr4bp_family_in_musun(
        anchor,
        seed_mu_sun=0.0,
        target_mu_sun=target_mu_sun,
        n_steps=n_steps,
        step_method="geometric",
        corrector_tol=1e-10,
        closure_tol=1.0,  # generous (free-T not Sun-commensurate, same as #303)
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        monodromy=True,
        sun_commensurate_n=1,
        on_step=_on_step,
    )
    elapsed = time.time() - t_start

    # Per-member rows.
    for i, mem in enumerate(family.members):
        orb = mem.orbit
        omega = orb.system.omega_sun_nondim
        sun_n = orb.sun_commensurate_n
        sun_drift = abs(omega * orb.period_nondim - 2.0 * math.pi * sun_n)
        m_row = _empty_row(moon, "member")
        m_row.step_idx = i
        m_row.mu_sun_value = float(mem.mu_sun_value)
        m_row.x0 = float(orb.state_initial[0])
        m_row.y0 = float(orb.state_initial[1])
        m_row.z0 = float(orb.state_initial[2])
        m_row.vx0 = float(orb.state_initial[3])
        m_row.vy0 = float(orb.state_initial[4])
        m_row.vz0 = float(orb.state_initial[5])
        m_row.T_TU = float(orb.period_nondim)
        m_row.T_seconds = float(orb.period_nondim * consts.t_s)
        m_row.corrector_residual = float(orb.corrector_residual)
        m_row.independent_closure_residual = float(orb.independent_closure_residual)
        m_row.sun_phase_drift = float(sun_drift)
        m_row.stability_tag = str(mem.stability_tag)
        rows.append(m_row)

    final_mu_sun = family.members[-1].mu_sun_value if family.members else 0.0
    survival_frac = final_mu_sun / target_mu_sun if target_mu_sun > 0 else 0.0
    if survival_frac >= 0.999:
        verdict = "family_survives_full_sun_perturbation"
    elif survival_frac >= 0.5:
        verdict = "family_partial_survival"
    elif survival_frac > 0.0:
        verdict = "family_breaks_early"
    else:
        verdict = "no_members_converged"
    summary = {
        "moon_name": moon,
        "phase": "continuation_complete",
        "n_members_converged": len(family.members),
        "n_steps_attempted": n_steps,
        "mu_sun_extent": list(family.mu_sun_extent),
        "target_mu_sun": target_mu_sun,
        "survival_fraction": survival_frac,
        "verdict": verdict,
        "walk_notes": family.walk_notes,
        "elapsed_seconds": elapsed,
    }
    s_row = _empty_row(moon, "summary")
    s_row.extras = summary
    rows.append(s_row)
    print(f"[scan-313-B] {moon}: SUMMARY {summary}", flush=True)
    return rows, summary


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    data_dir = repo_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    sha = _git_sha()
    t_start = time.time()

    targets = [
        ("Europa", GM_EUROPA, A_EUROPA, data_dir / "scan_313_sun_jupiter_europa.jsonl"),
        ("Io", GM_IO, A_IO, data_dir / "scan_313_sun_jupiter_io.jsonl"),
    ]
    summaries = []
    for moon, gm, a, out_path in targets:
        consts = _build_system_constants(moon, gm, a)
        rows, summary = _scout_one_system(consts)
        with out_path.open("w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(r.to_jsonl() + "\n")
        n_member = sum(1 for r in rows if r.row_type == "member")
        print(
            f"[scan-313-B] WROTE {out_path} : {len(rows)} rows ({n_member} members), sha={sha}",
            flush=True,
        )
        summaries.append(summary)

    total = time.time() - t_start
    print(f"\n[scan-313-B] PART B DONE : total wall = {total:.1f}s, sha={sha}", flush=True)
    for s in summaries:
        print(
            f"[scan-313-B] {s['moon_name']}: phase={s['phase']} verdict={s['verdict']} "
            f"members={s['n_members_converged']}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
