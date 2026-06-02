"""Exploratory derivation of the S1L1 cycler's multi-rev encounter topology.

The S1L1 cycler is a 2-synodic-period Earth-Mars cycler. Its exact encounter
topology (the Earth/Mars flyby sequence and the heliocentric revolution count
per leg) is provisional in our catalogue. This script searches candidate
2-synodic E-M multi-rev topologies and reports which reproduce the published
V-infinity anchors.

PROVENANCE
----------
The only SOURCED (published) values used here are the V-infinity anchors
``5.65 km/s at Earth`` and ``3.05 km/s at Mars``, plus the ~154-day outbound
Earth->Mars time-of-flight. Every interval ToF and revolution count this
script prints is COMPUTED by the optimiser, not published.

# FINDING (2026-06-02):
#   no circular-coplanar hit; model limitation.
#   Neither ("E","M","E") nor ("E","M","E","E") at period_k=2, over all per-leg
#   rev combos in {0,1,2} and branches {single,low,high}, produced a result
#   satisfying constraints AND |vinf_E - 5.65| < 0.3 AND |vinf_M - 3.05| < 0.3
#   under Ephemeris(model="circular"), even after differential_evolution
#   confirmation of the closest shortlisted candidates. See the cheap-scan
#   table in the run log for the achieved vinf pairs. As with the Aldrin
#   cycler, S1L1's published anchors appear not to be hostable in the
#   circular-coplanar idealised model.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import OptimisationResult, optimise_cell_idealized
from cyclerfinder.search.sequence import Cell

# --- Sourced anchors (published) -------------------------------------------
VINF_E_KMS = 5.65
VINF_M_KMS = 3.05
VINF_CAP_KMS = VINF_E_KMS + 2.5  # 8.15

# --- Search configuration --------------------------------------------------
PERIOD_K = 2  # 2-synodic cycler
SEQUENCES: tuple[tuple[str, ...], ...] = (
    ("E", "M", "E"),
    ("E", "M", "E", "E"),
)
REV_CHOICES = (0, 1, 2)
SECONDS_PER_DAY = 86400.0

# Acceptance / shortlist tolerances (km/s).
STRICT_TOL = 0.3
SHORTLIST_TOL = 1.0

# Cheap-scan per-call soft time budget (seconds). The cheap scan uses
# use_de=False so calls are short; this is a guard, logged not enforced.
CHEAP_CALL_SOFT_BUDGET_S = 90.0


@dataclass(frozen=True)
class Candidate:
    """One topology to evaluate."""

    sequence: tuple[str, ...]
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]

    def cell(self) -> Cell:
        return Cell(
            bodies=("E", "M"),
            sequence=self.sequence,
            period_k=PERIOD_K,
            per_leg_revs=self.per_leg_revs,
            per_leg_branch=self.per_leg_branch,
        )


def _branches_for(revs: int) -> tuple[str, ...]:
    """Valid branch labels for a leg with ``revs`` heliocentric revolutions."""
    if revs == 0:
        return ("single",)
    return ("low", "high")


def enumerate_candidates() -> list[Candidate]:
    """All (sequence, per_leg_revs, per_leg_branch) topologies to try."""
    candidates: list[Candidate] = []
    for sequence in SEQUENCES:
        n_legs = len(sequence) - 1
        for rev_combo in itertools.product(REV_CHOICES, repeat=n_legs):
            branch_options = [_branches_for(r) for r in rev_combo]
            for branch_combo in itertools.product(*branch_options):
                candidates.append(
                    Candidate(
                        sequence=sequence,
                        per_leg_revs=tuple(rev_combo),
                        per_leg_branch=tuple(branch_combo),
                    )
                )
    return candidates


def _max_vinf_for_body(result: OptimisationResult, body: str) -> float:
    """Max over the body's encounters of max(|vinf_in|, |vinf_out|), km/s."""
    best = 0.0
    found = False
    for enc in result.best_cycler.encounters:
        if enc.body != body:
            continue
        found = True
        v_in = float(np.linalg.norm(enc.vinf_in))
        v_out = float(np.linalg.norm(enc.vinf_out))
        best = max(best, v_in, v_out)
    return best if found else float("nan")


def _leg_tofs_days(result: OptimisationResult) -> list[float]:
    """Per-leg time-of-flight in days for the optimised cycler."""
    return [(leg.t_arrive - leg.t_depart) / SECONDS_PER_DAY for leg in result.best_cycler.legs]


def cheap_scan(
    candidates: list[Candidate],
    ephem: Ephemeris,
) -> list[tuple[Candidate, float, float]]:
    """Pass 1: fast SLSQP-only scan over all candidates.

    Prints one line per candidate and returns the shortlist of candidates
    whose achieved vinf at both bodies is within ``SHORTLIST_TOL`` of the
    anchors.
    """
    print("=== CHEAP SCAN (n_starts=3, use_de=False) ===")
    print(f"{'cell.id':<46} {'vinf_E':>8} {'vinf_M':>8} {'residual':>10} {'feasible':>8}")
    shortlist: list[tuple[Candidate, float, float]] = []
    for cand in candidates:
        cell = cand.cell()
        t0 = time.monotonic()
        try:
            result = optimise_cell_idealized(
                cell,
                ephem,
                vinf_cap=VINF_CAP_KMS,
                n_starts=3,
                seed=0,
                use_de=False,
            )
        except Exception as exc:
            print(f"{cell.id:<46} ERROR: {exc!r}")
            continue
        elapsed = time.monotonic() - t0
        if elapsed > CHEAP_CALL_SOFT_BUDGET_S:
            print(
                f"  (note: {cell.id} took {elapsed:.0f}s, "
                f"over the {CHEAP_CALL_SOFT_BUDGET_S:.0f}s soft budget)"
            )
        vinf_e = _max_vinf_for_body(result, "E")
        vinf_m = _max_vinf_for_body(result, "M")
        print(
            f"{cell.id:<46} {vinf_e:>8.3f} {vinf_m:>8.3f} "
            f"{result.closure_residual_kms:>10.4f} "
            f"{result.constraints_satisfied!s:>8}"
        )
        if (
            np.isfinite(vinf_e)
            and np.isfinite(vinf_m)
            and abs(vinf_e - VINF_E_KMS) < SHORTLIST_TOL
            and abs(vinf_m - VINF_M_KMS) < SHORTLIST_TOL
        ):
            # Rank metric: combined distance to both anchors.
            shortlist.append((cand, vinf_e, vinf_m))
    return shortlist


def _anchor_distance(vinf_e: float, vinf_m: float) -> float:
    return abs(vinf_e - VINF_E_KMS) + abs(vinf_m - VINF_M_KMS)


@dataclass(frozen=True)
class Hit:
    """A candidate that passed the strict acceptance test."""

    candidate: Candidate
    vinf_e: float
    vinf_m: float
    residual: float
    leg_tofs_days: tuple[float, ...]


def confirm(
    shortlist: list[tuple[Candidate, float, float]],
    ephem: Ephemeris,
) -> list[Hit]:
    """Pass 2: DE-confirm at most the 3 best shortlisted candidates."""
    print()
    print("=== CONFIRM (n_starts=5, use_de=True) ===")
    if not shortlist:
        print("(shortlist empty; nothing to confirm)")
        return []

    ranked = sorted(shortlist, key=lambda t: _anchor_distance(t[1], t[2]))
    to_confirm = ranked[:3]
    hits: list[Hit] = []
    for cand, scan_e, scan_m in to_confirm:
        cell = cand.cell()
        print(f"confirming {cell.id} (cheap vinf_E={scan_e:.3f}, vinf_M={scan_m:.3f})")
        try:
            result = optimise_cell_idealized(
                cell,
                ephem,
                vinf_cap=VINF_CAP_KMS,
                n_starts=5,
                seed=0,
                use_de=True,
            )
        except Exception as exc:
            print(f"  ERROR: {exc!r}")
            continue
        vinf_e = _max_vinf_for_body(result, "E")
        vinf_m = _max_vinf_for_body(result, "M")
        feasible = result.constraints_satisfied
        accept = (
            feasible
            and np.isfinite(vinf_e)
            and np.isfinite(vinf_m)
            and abs(vinf_e - VINF_E_KMS) < STRICT_TOL
            and abs(vinf_m - VINF_M_KMS) < STRICT_TOL
        )
        print(
            f"  vinf_E={vinf_e:.3f} vinf_M={vinf_m:.3f} "
            f"residual={result.closure_residual_kms:.4f} "
            f"feasible={feasible} accept={accept}"
        )
        if accept:
            hits.append(
                Hit(
                    candidate=cand,
                    vinf_e=vinf_e,
                    vinf_m=vinf_m,
                    residual=result.closure_residual_kms,
                    leg_tofs_days=tuple(_leg_tofs_days(result)),
                )
            )
    return hits


def report_hits(hits: list[Hit]) -> None:
    print()
    print("=== HITS ===")
    if not hits:
        print(
            "NO HITS: zero candidates passed the strict 0.3 km/s acceptance "
            "test. S1L1's published anchors appear not to be hostable in the "
            "circular-coplanar idealised model (model limitation)."
        )
        return
    for hit in sorted(hits, key=lambda h: h.residual):
        cell = hit.candidate.cell()
        tofs = ", ".join(f"{t:.2f}" for t in hit.leg_tofs_days)
        print(f"cell.id          = {cell.id}")
        print(f"  per_leg_revs   = {hit.candidate.per_leg_revs}")
        print(f"  per_leg_branch = {hit.candidate.per_leg_branch}")
        print(f"  vinf_E         = {hit.vinf_e:.3f} km/s (anchor 5.65)")
        print(f"  vinf_M         = {hit.vinf_m:.3f} km/s (anchor 3.05)")
        print(f"  residual       = {hit.residual:.4f} km/s")
        print(f"  leg ToFs (days)= [{tofs}]")


def main() -> None:
    ephem = Ephemeris(model="circular")
    candidates = enumerate_candidates()
    print(f"enumerated {len(candidates)} candidate topologies")
    shortlist = cheap_scan(candidates, ephem)
    print()
    print(f"shortlist size: {len(shortlist)}")
    hits = confirm(shortlist, ephem)
    report_hits(hits)


if __name__ == "__main__":
    main()
