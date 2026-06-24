"""#432 ER3BP discovery campaign: continue rotating-frame cycler families into
e>0, Floquet-monitor for survival/death and bifurcations. Report-only; no
catalogue writeback. See docs/superpowers/specs/2026-06-24-er3bp-discovery-campaign-design.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem


@dataclass(frozen=True)
class Er3bpSeed:
    label: str
    system: ER3BPSystem  # e=0 CR3BP μ lives in system.mu; system.e is the target
    state0: NDArray[np.float64]  # rotating-frame IC at e=0, shape (6,)
    period_f: float  # true-anomaly period (multiple of 2π for full period)
    is_half_period_residual: bool
    target_e: float
    source: str  # provenance string


# Broucke (1969) TR 32-1360 Table 12, Family 7P, Earth-Moon mu=0.0121550, Orbit 1.
_BROUCKE_EM_MU = 0.0121550
_BROUCKE_EM_ORBIT1 = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])


def standard_family_seeds(*, target_e: float = 0.0549) -> list[Er3bpSeed]:
    """Guaranteed seed floor from ICs already encoded in the repo (Earth-Moon).

    Currently the sourced Broucke-1969 Earth-Moon family. Additional Earth-Moon
    seeds (e.g. the Koblick NRHO table) are added by Task 4's catalogue provider.
    """
    sys = ER3BPSystem(mu=_BROUCKE_EM_MU, e=target_e, primary_name="E", secondary_name="M")
    return [
        Er3bpSeed(
            label="broucke-1969-em-7P-orbit1",
            system=sys,
            state0=_BROUCKE_EM_ORBIT1.copy(),
            period_f=2.0 * np.pi,
            is_half_period_residual=True,
            target_e=target_e,
            source="Broucke 1969 TR 32-1360 Table 12 Family 7P Orbit 1 (mu=0.0121550)",
        )
    ]
