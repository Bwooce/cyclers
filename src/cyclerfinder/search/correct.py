"""N-arc ballistic differential corrector on the real ephemeris (spec §2.1).

Generalises scripts/correct_s1l1_twoarc.py: free vars x = [t0, leg ToFs] with
one leg pinned by the sourced period; residuals = flyby V_inf-magnitude
continuity + periodicity closure, driven to zero with least_squares; bend
feasibility checked post-hoc, never in the residual. Pure: depends only on
core/lambert, core/ephemeris, core/constants.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BallisticClosureResult:
    t0_sec: float
    tof_days: tuple[float, ...]
    max_residual_kms: float
    vinf_per_encounter_kms: tuple[float, ...]
    converged: bool
    bend_feasible: bool
    vinf_cap_ok: bool = True

    @property
    def constraints_satisfied(self) -> bool:
        return self.converged and self.bend_feasible and self.vinf_cap_ok
