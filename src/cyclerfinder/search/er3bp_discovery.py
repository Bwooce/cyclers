"""#432 ER3BP discovery campaign: continue rotating-frame cycler families into
e>0, Floquet-monitor for survival/death and bifurcations. Report-only; no
catalogue writeback. See docs/superpowers/specs/2026-06-24-er3bp-discovery-campaign-design.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_continuation import continue_er3bp_family_in_e
from cyclerfinder.search.er3bp_floquet import er3bp_monodromy, floquet_classify


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


@dataclass(frozen=True)
class Er3bpStep:
    e: float
    corrector_residual: float
    stability_tag: str
    on_unit_circle: bool


@dataclass(frozen=True)
class Er3bpContinuationTrace:
    seed_label: str
    outcome: str  # "survives" | "dies" | "bifurcates"
    steps: tuple[Er3bpStep, ...]
    e_max_reached: float
    e_star: float | None  # first bifurcation eccentricity, if any
    target_e: float


def continue_and_monitor(seed: Er3bpSeed, *, n_steps: int = 20) -> Er3bpContinuationTrace:
    """Continue ``seed`` e=0->target_e, Floquet-monitor each step, classify."""
    # The corrector's ``period_f`` is the integration span: the half-period when
    # the symmetry flag is set. The seed stores the FULL period, so halve it.
    integration_f = seed.period_f / 2.0 if seed.is_half_period_residual else seed.period_f
    try:
        family = continue_er3bp_family_in_e(
            ER3BPSystem(
                mu=seed.system.mu,
                e=0.0,
                primary_name=seed.system.primary_name,
                secondary_name=seed.system.secondary_name,
            ),
            seed.state0,
            integration_f,
            seed.target_e,
            n_steps,
            is_half_period_residual=seed.is_half_period_residual,
        )
    except Exception:
        family = []
    steps: list[Er3bpStep] = []
    e_star: float | None = None
    # Bifurcation in a Hamiltonian family is an elliptic<->hyperbolic stability
    # TRANSITION along the continuation (an eigenvalue leaving the unit circle),
    # NOT "an eigenvalue sits on the unit circle" (which is true for EVERY stable
    # orbit). We track the regime across steps and flag the first flip.
    prev_regime: str | None = None
    for orb in family:
        sys_e = ER3BPSystem(
            mu=orb.mu,
            e=orb.e,
            primary_name=seed.system.primary_name,
            secondary_name=seed.system.secondary_name,
        )
        try:
            mono = er3bp_monodromy(orb.state0, orb.period_f, sys_e)
            fl = floquet_classify(mono)
            tag, on_uc = fl.stability_tag, fl.on_unit_circle
        except Exception:
            tag, on_uc = "unknown", False
        steps.append(
            Er3bpStep(
                e=orb.e,
                corrector_residual=orb.corrector_residual,
                stability_tag=tag,
                on_unit_circle=on_uc,
            )
        )
        regime = (
            "hyperbolic"
            if tag == "unstable"
            else ("elliptic" if tag in ("stable", "marginal") else "unknown")
        )
        if (
            regime != "unknown"
            and prev_regime is not None
            and prev_regime != "unknown"
            and regime != prev_regime
            and e_star is None
        ):
            e_star = orb.e  # first genuine elliptic<->hyperbolic transition
        if regime != "unknown":
            prev_regime = regime
    e_max = steps[-1].e if steps else 0.0
    reached_target = bool(steps) and abs(e_max - seed.target_e) <= 1e-6
    if e_star is not None:
        outcome = "bifurcates"
    elif reached_target:
        outcome = "survives"
    else:
        outcome = "dies"
    return Er3bpContinuationTrace(
        seed_label=seed.label,
        outcome=outcome,
        steps=tuple(steps),
        e_max_reached=e_max,
        e_star=e_star,
        target_e=seed.target_e,
    )
