"""Phase 0 acceptance run: reproduce the Ross & Roberts-Tsoukkas 2025 (AAS 25-621)
stable Earth-Moon (k1,k2) cycler families with the fixed-Jacobi symmetric
corrector + Barden stability (#212 Part B).

For each family this emits the recovered (state0, T, C, nu) tuple, the per-family
residuals (period vs published T^stable; |nu|; crossing residual; independent
Radau closure), and the stable/unstable verdict. This is the artifact #216 will
ingest (proposed rows; NO catalogue writeback here).

SOURCED-GOLDEN DISCIPLINE: EXPECTED values are Ross's PRINTED mu (p. 3) and the
C^stable / T^stable / stability columns of Table 3 (p. 11). The recovered x0 /
ydot0 are DERIVED (the 1-D solve of the paper's symmetric-orbit structure, §5),
stored with provenance, never themselves goldens.

DATA GAPS (do NOT silently pick one): two C_(k1,k2) *bound*-column values in
Table 3 are internally inconsistent with Eq. 8 + Table 1 / Table 4 (up to
7.7e-3 for (3,1)); see the results note. Those bound columns are NOT used here
(only the internally-consistent C^stable / T^stable columns are).

Usage:
    uv run python scripts/cr3bp_ross_reproduce.py
    uv run python scripts/cr3bp_ross_reproduce.py --report /tmp/ross_reproduce.txt
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# Ross & Roberts-Tsoukkas 2025, AAS 25-621, p. 3.
ROSS_MU = 1.2150584270572e-2
ROSS_L_KM = 384400.0
ROSS_TU_DAYS = 27.321661 / (2.0 * np.pi)  # 1 TU in days (p. 3)
ROSS_T_S = 375699.8  # 1 TU in seconds (p. 3, ~ sidereal month / 2pi)


@dataclass(frozen=True)
class FamilySeed:
    label: str
    c_stable: float  # Table 3 C^stable (p. 11)
    t_stable: float  # Table 3 T^stable, TU (p. 11)
    x0_seed: float  # derived seed for the 1-D solve (NOT a golden)
    ydot0_sign: float
    half_crossings: int
    note: str


# The five reproduced families (seeds discovered by the §3-step-1/2 construction
# region + a half-period-crossing-index scan; held here as DERIVED seeds).
SEEDS = [
    FamilySeed(
        "(1,1)",
        3.151175879508174,
        10.29206921007976,
        -0.7682140805,
        -1.0,
        3,
        "Delta_p_m 0.13 km; ~2:3 synodic resonance (p. 11)",
    ),
    FamilySeed(
        "(2,1)",
        3.129389531088256,
        19.44043166795154,
        0.7237335857,
        1.0,
        4,
        "Delta_p_m 4.23 km; 3.09 sidereal months (p. 11); razor-thin window",
    ),
    FamilySeed(
        "(3,1)",
        3.161784147013429,
        14.78849241668140,
        -0.3209891696,
        -1.0,
        3,
        "window 750-1000 km perilune alt; Delta_p_m 253.70 km (Fig. 8)",
    ),
    FamilySeed(
        "(3,2)",
        3.182762663084288,
        17.90058010350006,
        -0.3210000000,
        -1.0,
        6,
        "Delta_p_m 42.08 km; 2 windows; near 2:5 synodic (p. 13); half-period "
        "= 6th x-axis crossing",
    ),
    FamilySeed(
        "(3,3)",
        3.177224018696528,
        18.14546057589189,
        -0.3217380626,
        -1.0,
        5,
        "window 4200-6200 km perilune alt; Delta_p_m 2041.34 km; 5 windows (p. 14)",
    ),
]


@dataclass(frozen=True)
class Reproduction:
    seed: FamilySeed
    orbit: cp.SymmetricOrbit
    nu: float
    abs_lambda: float
    radau_ok: bool
    radau_dj: float


def _system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=ROSS_MU, primary="Earth", secondary="Moon", l_km=ROSS_L_KM, t_s=ROSS_T_S
    )


def reproduce(seed: FamilySeed) -> Reproduction:
    sysm = _system()
    orbit = cp.correct_symmetric_fixed_jacobi(
        sysm,
        seed.x0_seed,
        seed.c_stable,
        seed.t_stable,
        ydot0_sign=seed.ydot0_sign,
        half_crossings=seed.half_crossings,
        tol=1e-10,
    )
    nu, lam = cp.barden_stability(sysm, orbit)
    po = cp.PeriodicOrbit(
        state0=np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0]),
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=orbit.converged,
        closure_residual=orbit.crossing_residual,
    )
    radau_ok, radau_dj = cp.crosscheck_periodic(sysm, po, closure_tol=1e-3, jacobi_tol=1e-8)
    return Reproduction(seed, orbit, float(nu), abs(lam), radau_ok, float(radau_dj))


def _fmt_state(orbit: cp.SymmetricOrbit) -> str:
    s = [orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0]
    return "[" + ", ".join(f"{v:.13g}" for v in s) + "]"


def build_report(reps: list[Reproduction], elapsed_s: float) -> str:
    ts = datetime.now(tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines: list[str] = []
    lines.append("Ross & Roberts-Tsoukkas 2025 (AAS 25-621) family reproduction")
    lines.append(f"run: {ts}   elapsed: {elapsed_s:.1f} s   mu (p.3): {ROSS_MU}")
    lines.append("=" * 78)
    n_stable = 0
    for r in reps:
        s = r.seed
        o = r.orbit
        stable = abs(r.nu) < 1.0
        n_stable += int(stable)
        lines.append(f"\nFamily {s.label}   ({s.note})")
        lines.append(f"  converged       : {o.converged}  (iters {o.n_iter})")
        lines.append(f"  state0 (nd)     : {_fmt_state(o)}   <- DERIVED, not a golden")
        lines.append(f"  C  (enforced)   : {o.jacobi:.15g}  (pub C^stable {s.c_stable:.15g})")
        lines.append(f"     dC           : {o.jacobi - s.c_stable:+.2e}")
        lines.append(f"  T = 2*t_half    : {o.period:.15g}  (pub T^stable {s.t_stable:.15g} TU)")
        lines.append(
            f"     dT           : {o.period - s.t_stable:+.2e} TU"
            f"  ({(o.period - s.t_stable) * ROSS_TU_DAYS:+.2e} d)"
        )
        lines.append(f"     T in days    : {o.period * ROSS_TU_DAYS:.6f}")
        lines.append(f"  crossing resid  : {o.crossing_residual:.2e}")
        lines.append(f"  nu = 1/2(l+1/l) : {r.nu:+.9f}   |lambda| {r.abs_lambda:.6f}")
        lines.append(
            f"  VERDICT         : {'STABLE (|nu|<1)' if stable else 'UNSTABLE'}"
            "  [Table-3 member is the nu=0 midpoint]"
        )
        lines.append(f"  Radau crosscheck: ok={r.radau_ok}  dJ={r.radau_dj:.2e}")
    lines.append("\n" + "=" * 78)
    lines.append(f"Reproduced families: {len(reps)} ; stable (|nu|<1): {n_stable}")
    lines.append("")
    lines.append("DATA GAPS (Table 3 vs Eq.8/Table1 + Table 4 -- do NOT pick one):")
    lines.append(
        "  * C_(2,1): Table 3 3.1297495000000 vs Eq.8 min 3.129751730201047 (Delta 2.24e-6)."
    )
    lines.append(
        "  * C_(3,1): Table 3 3.1833333078762 vs Table-4-implied 3.1756140 (Delta 7.7e-3)."
    )
    lines.append(
        "  These are the C_(k1,k2) BOUND columns; the C^stable/T^stable columns "
        "used above are internally consistent."
    )
    lines.append("\nNO catalogue writeback. Proposed (state0,T,C,nu) tuples for #216.")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproduce Ross & Roberts-Tsoukkas 2025 stable EM cycler families"
    )
    parser.add_argument("--report", type=Path, default=None, help="Write the runlog here")
    args = parser.parse_args()

    t0 = datetime.now(tz=UTC)
    print(f"cr3bp_ross_reproduce: mu={ROSS_MU} (Ross p.3)")
    reps: list[Reproduction] = []
    for seed in SEEDS:
        print(f"  reproducing {seed.label} ...", flush=True)
        reps.append(reproduce(seed))
    elapsed = (datetime.now(tz=UTC) - t0).total_seconds()

    report = build_report(reps, elapsed)
    print()
    print(report)
    if args.report is not None:
        args.report.write_text(report + "\n", encoding="utf-8")
        print(f"\nRunlog written: {args.report}")


if __name__ == "__main__":
    main()
