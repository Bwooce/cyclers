"""Discovery run: continue held Earth-Moon stable cyclers in the mass parameter
mu up into the binary-star range, recovering the families the
Roberts-Tsoukkas & Ross 2026 journal paper depicts in FIGURES ONLY (Fig. 3:
(1,3) at mu=0.1, (3,1) at mu=0.3, (1,1) at mu=0.5 -- all drawn stable, no
numbers printed).

Two families are continued (same-class as a journal target):
  * (1,1)  mu = 0.01215 -> 0.5  (equal-mass target, Fig. 3 bottom)
  * (3,1)  mu = 0.01215 -> 0.3  (Fig. 3 middle)

For each: pseudo-arclength continuation in (x0, C, mu) of the held nu=0 member,
then -- because the paper's depicted orbit is STABLE and the held branch may
turn unstable before the target mu -- a fixed-mu C-family scan at the target mu
to locate the stable subfamily (the paper's method) and report a figure-match
candidate.

DISCOVERY DISCIPLINE: the recovered (mu, C, T, IC, nu) are OUR OWN computed
values (the paper prints none). They are DISCOVERIES requiring the full V0-V5
gauntlet, NOT sourced rows. NO catalogue writeback here; this emits evidence for
a review-gated note.

Usage:
    uv run python scripts/cr3bp_mu_continuation.py
    uv run python scripts/cr3bp_mu_continuation.py --family 31 --report /tmp/mu31.txt
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.mu_continuation as mc

ROSS_MU = 1.2150584270572e-2
ROSS_L_KM = 384400.0
ROSS_T_S = 375699.8


def _system(mu: float) -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=mu, primary="Earth", secondary="Moon", l_km=ROSS_L_KM, t_s=ROSS_T_S)


@dataclass(frozen=True)
class FamilyTarget:
    label: str
    c_stable: float
    t_stable: float
    x0_seed: float
    ydot0_sign: float
    half_crossings: int
    mu_target: float
    journal_class: str


TARGETS = {
    "11": FamilyTarget(
        "(1,1)",
        3.151175879508174,
        10.29206921007976,
        -0.7682140805,
        -1.0,
        3,
        0.5,
        "Fig. 3 bottom: (1,1) equal-mass cycler at mu=0.5",
    ),
    "31": FamilyTarget(
        "(3,1)",
        3.161784147013429,
        14.78849241668140,
        -0.3209891696,
        -1.0,
        3,
        0.3,
        "Fig. 3 middle: (3,1) cycler at mu=0.3",
    ),
}


def _seed(t: FamilyTarget) -> cp.SymmetricOrbit:
    return cp.correct_symmetric_fixed_jacobi(
        _system(ROSS_MU),
        t.x0_seed,
        t.c_stable,
        t.t_stable,
        ydot0_sign=t.ydot0_sign,
        half_crossings=t.half_crossings,
        tol=1e-11,
    )


def run_family(t: FamilyTarget, lines: list[str]) -> None:
    lines.append(f"\n{'=' * 78}")
    lines.append(f"FAMILY {t.label}  ->  mu = {t.mu_target}   [{t.journal_class}]")
    lines.append("=" * 78)
    seed = _seed(t)
    lines.append(
        f"seed (mu={ROSS_MU:.10f}): x0={seed.x0:.10f} C={seed.jacobi:.12f} "
        f"T={seed.period:.10f} conv={seed.converged}"
    )

    branch = mc.continue_in_mu(
        seed,
        ROSS_MU,
        half_crossings=t.half_crossings,
        ydot0_sign=t.ydot0_sign,
        mu_target=t.mu_target,
        label=t.label,
        record_every=4,
    )
    lines.append(
        f"\nARCLENGTH CONTINUATION ({len(branch.members)} recorded members, "
        f"{branch.n_steps} steps, stop={branch.stop_reason}):"
    )
    lines.append(f"  {'mu':>9} {'C':>14} {'T':>13} {'nu':>12} {'res':>9} {'radauDJ':>9}  verdict")
    for m in branch.members:
        lines.append(
            f"  {m.mu:9.5f} {m.jacobi:14.9f} {m.period:13.8f} {m.nu:+12.5f} "
            f"{m.crossing_residual:9.1e} {m.radau_djacobi:9.1e}  "
            f"{'STABLE' if m.stable else 'unstable'}"
        )
    landed = branch.members[-1] if branch.members else None
    if landed is not None:
        lines.append(
            f"\nLANDED at mu={landed.mu:.6f}: x0={landed.x0:.10f} ydot0={landed.ydot0:.10f}"
        )
        lines.append(
            f"  C={landed.jacobi:.12f} T={landed.period:.10f} nu={landed.nu:+.6f} "
            f"({'STABLE' if landed.stable else 'UNSTABLE'})  "
            f"cross_res={landed.crossing_residual:.2e} radau_dJ={landed.radau_djacobi:.2e}"
        )

    # Stable-subfamily search at the target mu (the paper's depicted orbit is stable).
    if landed is not None and abs(landed.mu - t.mu_target) < 1e-4:
        lines.append(
            f"\nC-FAMILY SCAN at mu={t.mu_target} (locate the stable subfamily, paper's method):"
        )
        members = mc.scan_c_family_at_mu(
            t.mu_target,
            landed.x0,
            landed.jacobi,
            landed.period,
            half_crossings=t.half_crossings,
            ydot0_sign=t.ydot0_sign,
            dc=3e-3,
            n_each=25,
        )
        stable = [m for m in members if m.stable]
        lines.append(
            f"  scanned {len(members)} converged members; {len(stable)} linearly stable (|nu|<1)."
        )
        for m in members:
            flag = "  <== STABLE" if m.stable else ""
            lines.append(
                f"  mu={m.mu:.4f} C={m.jacobi:13.8f} T={m.period:12.7f} "
                f"nu={m.nu:+11.5f} x0={m.x0:+.7f}{flag}"
            )
        if stable:
            best = min(stable, key=lambda m: abs(m.nu))
            lines.append(f"\n  FIGURE-MATCH CANDIDATE (most-stable {t.label} at mu={t.mu_target}):")
            lines.append(f"    state0 = [{best.x0:.12g}, 0, 0, 0, {best.ydot0:.12g}, 0]  (DERIVED)")
            lines.append(
                f"    mu={best.mu:.10f} C={best.jacobi:.12f} T={best.period:.10f} "
                f"nu={best.nu:+.6f} |lambda|={best.abs_lambda:.6f}"
            )
            lines.append(
                f"    cross_res={best.crossing_residual:.2e} "
                f"radau_dJ={best.radau_djacobi:.2e}  -> genuine STABLE periodic orbit"
            )
        else:
            lines.append(
                f"\n  NO stable {t.label} member found in the scanned C-window at "
                f"mu={t.mu_target} (clean negative; widen dc/n_each or re-seed)."
            )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", choices=["11", "31", "all"], default="all")
    ap.add_argument("--report", type=Path, default=None)
    args = ap.parse_args()

    t0 = datetime.now(tz=UTC)
    lines: list[str] = []
    lines.append("CR3BP mass-parameter (mu) continuation -- binary-star cycler discovery")
    lines.append(f"run: {datetime.now(tz=UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}  mu0={ROSS_MU}")
    keys = ["11", "31"] if args.family == "all" else [args.family]
    for k in keys:
        print(f"continuing family {k} ...", flush=True)
        run_family(TARGETS[k], lines)
    lines.append(f"\nelapsed: {(datetime.now(tz=UTC) - t0).total_seconds():.1f} s")
    lines.append("\nNO catalogue writeback. DISCOVERED (mu,C,T,IC,nu) tuples are review-gated.")

    report = "\n".join(lines)
    print(report)
    if args.report is not None:
        args.report.write_text(report + "\n", encoding="utf-8")
        print(f"\nRunlog written: {args.report}")


if __name__ == "__main__":
    main()
