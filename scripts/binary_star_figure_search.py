"""Discovery run (#255): figure-seeded direct search for the binary-star cyclers
Roberts-Tsoukkas & Ross 2026 depict in Fig. 3 (FIGURES ONLY, no printed numbers):

  * (1,3) exterior cycler at mu = 0.1
  * (3,1) cycler          at mu = 0.3
  * (1,1) equal-mass      at mu = 0.5   (all drawn STABLE)

#252 continued the held Earth-Moon members in mu and branch-switched off the
cycler family before these targets. Here we instead seed the fixed-mu symmetric
corrector across the figure-read (x0, C) box and keep only orbits whose winding
topology matches the depicted (k1, k2), are prograde, and are linearly stable
(Barden |nu| < 1), each confirmed by an independent-Radau Jacobi cross-check.

DISCOVERY DISCIPLINE: a survivor is OUR OWN (mu, C, T, IC, nu) -- the paper
prints none. It is a discovery CANDIDATE requiring the literature-novelty check
(#261) + the V0-V5 gauntlet, never a sourced row. NO catalogue writeback.

Usage:
    uv run python scripts/binary_star_figure_search.py --report /tmp/bstar.txt
    uv run python scripts/binary_star_figure_search.py --mu 0.3
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from cyclerfinder.search.binary_star_search import (
    BinaryStarCandidate,
    figure_seeded_search,
)


@dataclass(frozen=True)
class FigureTarget:
    mu: float
    k1: int
    k2: int
    x0_lo: float
    x0_hi: float
    c_lo: float
    c_hi: float
    note: str


# Figure-read windows off Roberts-Tsoukkas & Ross 2026 Fig. 3 (rotating frame, DU).
# Jacobi bands bracket the #252 continuation's stable points (~3.0-3.18).
TARGETS: dict[str, FigureTarget] = {
    "0.1": FigureTarget(
        0.1, 1, 3, -3.0, 3.2, 3.00, 3.18, "Fig.3 top: (1,3) exterior, primaries at -0.1/+0.9"
    ),
    "0.3": FigureTarget(
        0.3, 3, 1, -0.85, 0.90, 3.00, 3.16, "Fig.3 middle: (3,1), primaries at -0.3/+0.7"
    ),
    "0.5": FigureTarget(
        0.5, 1, 1, -0.65, 0.65, 3.05, 3.16, "Fig.3 bottom: (1,1) equal-mass, primaries at +-0.5"
    ),
}


def run_target(key: str, t: FigureTarget, lines: list[str]) -> list[BinaryStarCandidate]:
    lines.append(f"\n{'=' * 78}")
    lines.append(f"mu = {t.mu}   target ({t.k1},{t.k2})   [{t.note}]")
    lines.append("=" * 78)
    cands = figure_seeded_search(
        t.mu,
        t.k1,
        t.k2,
        x0_lo=t.x0_lo,
        x0_hi=t.x0_hi,
        c_lo=t.c_lo,
        c_hi=t.c_hi,
        n_x0=26,
        n_c=14,
        half_crossings_set=(1, 2, 3, 4, 5, 6),
        ydot0_signs=(-1.0, 1.0),
        require_stable=True,
        require_prograde=True,
    )
    if not cands:
        lines.append(
            f"  NO stable prograde ({t.k1},{t.k2}) cycler found in the figure-read box "
            f"(x0 in [{t.x0_lo},{t.x0_hi}], C in [{t.c_lo},{t.c_hi}]) -- clean negative."
        )
        return cands
    cands.sort(key=lambda c: abs(c.nu))
    lines.append(f"  {len(cands)} stable prograde ({t.k1},{t.k2}) candidate(s):")
    for c in cands:
        lines.append(
            f"    x0={c.x0:+.10f} C={c.jacobi:.10f} T={c.period:.8f} nu={c.nu:+.5f} "
            f"|lam|={c.abs_lambda:.4f} (k1,k2)=({c.k1},{c.k2}) "
            f"xrange=[{c.x_min:+.3f},{c.x_max:+.3f}] "
            f"cross_res={c.crossing_residual:.1e} radau_dJ={c.radau_djacobi:.1e}"
        )
    best = cands[0]
    lines.append(f"\n  MOST-STABLE CANDIDATE (mu={t.mu}, DISCOVERY -- review-gated, NOT sourced):")
    lines.append(f"    state0 = [{best.x0:.12g}, 0, 0, 0, {best.ydot0:.12g}, 0]  mu={best.mu}")
    lines.append(
        f"    C={best.jacobi:.12f} T={best.period:.10f} nu={best.nu:+.6f} "
        f"(k1,k2)=({best.k1},{best.k2}) prograde={best.prograde}"
    )
    return cands


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mu", choices=[*sorted(TARGETS), "all"], default="all")
    ap.add_argument("--report", type=Path, default=None)
    args = ap.parse_args()

    t0 = datetime.now(tz=UTC)
    lines: list[str] = [
        "Binary-star figure-seeded cycler search (#255)",
        f"run: {t0.strftime('%Y-%m-%dT%H:%M:%SZ')}",
    ]
    keys = sorted(TARGETS) if args.mu == "all" else [args.mu]
    total = 0
    for k in keys:
        print(f"searching mu={k} ...", flush=True)
        cands = run_target(k, TARGETS[k], lines)
        total += len(cands)
    lines.append(f"\nelapsed: {(datetime.now(tz=UTC) - t0).total_seconds():.1f} s")
    lines.append(
        f"TOTAL stable prograde figure-topology candidates: {total} "
        "(each a DISCOVERY requiring literature-novelty #261 + V0-V5; NO writeback)."
    )
    report = "\n".join(lines)
    print(report)
    if args.report is not None:
        args.report.write_text(report + "\n", encoding="utf-8")
        print(f"\nRunlog written: {args.report}")


if __name__ == "__main__":
    main()
