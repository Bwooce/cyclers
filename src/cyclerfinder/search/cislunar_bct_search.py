"""Cislunar BCT discovery driver (#378 Phase 3).

Sweeps the Belbruno ballistic-capture-transfer (BCT) family over its
parameters, constructs each via `genome/bct_transfer.py`, classifies it as a
`transfer` (precursor_mga capability output) or a `quasi_cycler_candidate` (the
return leg re-acquires the WSB surface W -- the only catalogue-row candidate,
design §6 / Table 1.1), runs the BCT novelty self-test, and produces records.

Honest verdict (design §6)
--------------------------
This is CAPABILITY-FIRST. A single BCT is a one-shot transfer, NOT a cycler
(Belbruno 2004: "WSB orbits are not cyclers"). The genuine catalogue-row
candidate is a *repeating* cislunar capture<->escape chain whose return leg
re-acquires W -- unproven, probably absent (Theorem 3.58: capture on W is a
CHAOTIC process, cutting against clean periodicity). The driver makes that
hypothesis *searchable*; it does NOT assert the object exists. A standalone BCT
that clears lit-check is routed as a capability record (precursor_mga), never
self-admitted to the catalogue (gauntlet/human gated).

Reuses
------
* `genome/bct_transfer.py` -- construct + the W-surface energy check;
* `core/wsb.py` -- the on-W re-acquisition predicate;
* `genome/bct_transfer.check_bct_novelty` -- the BCT lit-check seam;
* `data/empty_regions.py` -- the clean-negative registry (Phase 4 driver).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.wsb as wsb
import cyclerfinder.genome.bct_transfer as bct
from cyclerfinder.search.literature_check import CandidateSignature, SearchResult

BCTClassification = Literal["transfer", "quasi_cycler_candidate"]


@dataclass(frozen=True)
class BCTSearchGrid:
    """The (theta_2, e_2, back_days, branch) sweep grid.

    theta_2 is the family selector (which W manifold the backward arc escapes
    on); e_2 the capture eccentricity; back_days the arc-II horizon; branch the
    Moon-relative capture sense.
    """

    theta2_values: tuple[float, ...] = (0.65, 0.70, 0.75, 1.25, 1.60)
    e2_values: tuple[float, ...] = (0.95,)
    back_days_values: tuple[float, ...] = (70.0,)
    branches: tuple[bct.Branch, ...] = ("retrograde",)

    def points_total(self) -> int:
        """Total grid points (for the empty-region bound)."""
        return (
            len(self.theta2_values)
            * len(self.e2_values)
            * len(self.back_days_values)
            * len(self.branches)
        )


@dataclass(frozen=True)
class BCTRecord:
    """One classified BCT from the sweep."""

    theta2: float
    e2: float
    back_days: float
    branch: bct.Branch
    apoapsis_ld: float
    capture_e2: float
    classification: BCTClassification
    reacquires_w: bool


@dataclass(frozen=True)
class BCTSearchResult:
    """The aggregate sweep outcome."""

    records: tuple[BCTRecord, ...] = field(default_factory=tuple)
    n_evaluated: int = 0
    n_transfers: int = 0
    n_quasi_cycler_candidates: int = 0


def _return_leg_reacquires_w(
    arc: bct.BCTArc,
    system: bcr4bp.BCR4BPSystem,
    *,
    forward_days: float = 90.0,
    n_samples: int = 200,
) -> bool:
    """Does the FORWARD continuation past QF re-acquire W within the budget?

    A quasi-cycler must repeat: continue the trajectory FORWARD from QF and check
    whether it returns to a state ON W (bound to the Moon at a periapsis,
    E_2 <= 0, r-dot_23 ~ 0) after leaving the immediate capture -- a
    re-acquisition. The expected default is FALSE (a single capture, no clean
    repeat -- Theorem 3.58).
    """
    qf = arc.qf_state
    left_capture = False
    moon = np.array([1.0 - system.mu, 0.0, 0.0])
    r0 = float(np.linalg.norm(qf[:3] - moon))
    cur = qf.copy()
    t_prev = 0.0
    for d in np.linspace(forward_days / n_samples, forward_days, n_samples):
        t_nd = d / bct.TU_DAYS
        try:
            a = bcr4bp.propagate_bcr4bp(system, cur, t_nd - t_prev, t0=t_prev)
        except RuntimeError:
            return False
        cur = a.state_f
        t_prev = t_nd
        r = float(np.linalg.norm(cur[:3] - moon))
        # Must first LEAVE the capture vicinity (escape the immediate periapsis).
        if r > 3.0 * r0:
            left_capture = True
        # Then RE-ACQUIRE W: back on a Moon periapsis, bound, after leaving.
        if (
            left_capture
            and wsb.kepler_energy_moon(cur, system) <= 0.0
            and wsb.is_periapsis(cur, system, tol=1e-3)
            and r < 5.0 * r0
        ):
            return True
    return False


def classify_bct(
    arc: bct.BCTArc,
    system: bcr4bp.BCR4BPSystem,
    *,
    reacquires_w: bool | None = None,
) -> BCTClassification:
    """Classify a constructed BCT: `transfer` vs `quasi_cycler_candidate`.

    If ``reacquires_w`` is given it is used directly (test/override hook);
    otherwise it is computed via :func:`_return_leg_reacquires_w`. A re-acquiring
    return leg => `quasi_cycler_candidate` (the catalogue-row candidate);
    otherwise `transfer` (precursor_mga capability output).
    """
    if reacquires_w is None:
        reacquires_w = _return_leg_reacquires_w(arc, system)
    return "quasi_cycler_candidate" if reacquires_w else "transfer"


def run_cislunar_bct_search(
    grid: BCTSearchGrid,
    system: bcr4bp.BCR4BPSystem,
    *,
    check_reacquisition: bool = True,
) -> BCTSearchResult:
    """Sweep the grid, construct + classify each BCT, return the aggregate.

    For each grid point: build the backward arc (the Sun-shaped apoapsis + the
    on-W capture target), classify it, and record. ``check_reacquisition``
    controls whether the (expensive) return-leg W re-acquisition is computed; set
    False for a fast apoapsis-only census.
    """
    records: list[BCTRecord] = []
    n_transfer = 0
    n_chain = 0
    for branch in grid.branches:
        for theta2 in grid.theta2_values:
            for e2 in grid.e2_values:
                for back_days in grid.back_days_values:
                    target = bct.BCTTarget(e2=e2, theta2=theta2, branch=branch)
                    arc = bct.construct_bct_backward(target, system, back_days=back_days)
                    reacq = _return_leg_reacquires_w(arc, system) if check_reacquisition else False
                    classification = classify_bct(arc, system, reacquires_w=reacq)
                    if classification == "transfer":
                        n_transfer += 1
                    else:
                        n_chain += 1
                    records.append(
                        BCTRecord(
                            theta2=theta2,
                            e2=e2,
                            back_days=back_days,
                            branch=branch,
                            apoapsis_ld=arc.max_earth_apoapsis_ld,
                            capture_e2=wsb.kepler_energy_moon(arc.qf_state, system),
                            classification=classification,
                            reacquires_w=reacq,
                        )
                    )
    return BCTSearchResult(
        records=tuple(records),
        n_evaluated=len(records),
        n_transfers=n_transfer,
        n_quasi_cycler_candidates=n_chain,
    )


def is_novel_emittable(
    sig: CandidateSignature,
    *,
    search_fn: Callable[[str], list[SearchResult]],
) -> bool:
    """Is the BCT candidate novel-emittable (NOT a published rediscovery)?

    Runs the BCT novelty self-test (`bct.check_bct_novelty`). Returns False if the
    signature is flagged `published` (a Hiten rediscovery -> suppressed); True if
    `not-found` (eligible for human + gauntlet routing -- NECESSARY-not-sufficient
    for novelty, never a self-admission). An `inconclusive` search is treated as
    NOT-emittable (a not-found we cannot trust).
    """
    result = bct.check_bct_novelty(sig, search=search_fn)
    return result.status == "not-found"
