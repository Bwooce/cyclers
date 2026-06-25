"""V3 independent-integrator BCR4BP long-span gauntlet (#305 Part D).

Spec reinterpretation for BCR4BP periodic orbits
------------------------------------------------
Spec §14 V3 is "confirmed on an integrator independent of the finding solver;
bounded over 3-5 laps". V2-BCR4BP (:mod:`v2_bcr4bp`) measured the multi-lap
position drift under **DOP853** (the corrector's integrator). V3 re-runs the
SAME V2 long span (``k*T`` for ``k=1..n_cycles``) under an **independent
integrator family** and checks the per-lap terminal positions agree with V2's
within a floor — so the bounded-drift signature is integrator-independent, not
a single-integrator artefact.

Choice of independent integrator
--------------------------------
The independent family here is **scipy LSODA** — an Adams / BDF *multistep*
predictor-corrector (Hindmarsh ODEPACK), a genuinely different family from the
DOP853 single-step explicit Runge-Kutta the corrector and V2 use (and from the
Radau implicit-RK V1 used). It integrates the EXACT SAME ``bcr4bp_eom`` so no
re-derivation of the rotating-frame dynamics is introduced.

Why NOT REBOUND IAS15 (design draft §4 / risk register): driving REBOUND
through the BCR4BP would require re-expressing the full synodic rotating-frame
EOM (Coriolis + centrifugal + both primaries + Sun direct+indirect) as a
Python ``additional_forces`` callback — a substantial re-derivation and exactly
the bug surface the draft's `reference_rebound_variation_custom_force_gotcha`
flags. V3 needs only STATE re-propagation, so a different scipy family gives an
honest integrator-independent cross-check without that risk. (The stroboscopic
Option C upgrade can revisit a force-callback integrator later.)

The PERTURBED parity discipline
-------------------------------
The gotcha the draft warns about is that a cross-check can PASS even when the
Sun term is silently inert (a Sun-only-correct bug passes a mu_sun=0 test). V3
must therefore be exercised WITH the perturber active — the agreement number is
only meaningful on an orbit whose dynamics genuinely include the Sun term. The
frozen test includes a parity check: at the Andreu mu_sun the LSODA and DOP853
spans agree (Sun term wired into BOTH), and a deliberately Sun-term-stripped
re-propagation (mu_sun=0) DIVERGES from the mu_sun-on V2 span — confirming the
cross-check actually feels the Sun.

Floor rationale (labelled judgment call)
-----------------------------------------
``V3_BCR4BP_AGREEMENT_FLOOR_KMS = 100 km`` mirrors v3_3d's agreement floor —
two independent integrators of the same dynamics should agree to well under
this over 3 laps; a larger disagreement means one integrator is wrong.

Discipline
----------
* READ-ONLY on the BCR4BP genome (wrap, never re-solve).
* Chains on V2: requires ``len(v2_verdict.per_cycle_drift_km) >= n_cycles``.
* NO catalogue writeback.

References
----------
* spec §14 V3; design draft #305 §4 (V3), risk register (REBOUND gotcha).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.data.validation.v1_bcr4bp import SEM_L_KM
from cyclerfinder.data.validation.v2_bcr4bp import V2VerdictBCR4BP
from cyclerfinder.genome.bcr4bp_genome import BCR4BPPeriodicOrbit

V3_BCR4BP_N_CYCLES_MIN: Final[int] = 3
"""V3 minimum cycle count. Mirrors spec §14 ">=3 continuous laps"."""

V3_BCR4BP_AGREEMENT_FLOOR_KMS: Final[float] = 100.0
"""V3 per-lap integrator-agreement floor (km). LABELLED judgment call mirroring
v3_3d; two correct integrators of the same dynamics agree far below this."""


@dataclass(frozen=True)
class V3VerdictBCR4BP:
    """Frozen V3 verdict for a BCR4BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id :
        Identifier carried for the audit trail.
    n_cycles_requested :
        Lap count checked (must be ``>= V3_BCR4BP_N_CYCLES_MIN`` and
        ``<= len(v2_verdict.per_cycle_drift_km)``).
    per_cycle_agreement_km :
        Per-lap ``||X(kT)_pos[LSODA] - X(kT)_pos[DOP853]||`` (km). The DOP853
        terminal positions are recomputed here (V2 stored only the drift norm,
        not the endpoint), so both arms are produced fresh under their own
        integrators.
    max_agreement_km :
        Max of ``per_cycle_agreement_km`` — gated against ``agreement_floor_km``.
    agreement_floor_km :
        The labelled floor this verdict was held against.
    n_cycles_min :
        Spec minimum (3).
    independent_integrator :
        Name of the independent integrator family used ("LSODA").
    converged_each_cycle :
        Whether every lap propagated successfully under BOTH integrators.
    passes_v3_bcr4bp :
        ``converged_each_cycle AND n_cycles >= n_cycles_min AND
        max_agreement_km <= agreement_floor_km AND the V2 verdict it chains on
        passed``. Headline boolean.
    notes :
        Free-form audit string.
    """

    candidate_id: str
    n_cycles_requested: int
    per_cycle_agreement_km: tuple[float, ...]
    max_agreement_km: float
    agreement_floor_km: float
    n_cycles_min: int
    independent_integrator: str
    converged_each_cycle: bool
    passes_v3_bcr4bp: bool
    notes: str = ""


def _propagate_lsoda(
    system: bcr4bp.BCR4BPSystem,
    state0: np.ndarray,
    t: float,
    *,
    rtol: float,
    atol: float,
) -> np.ndarray | None:
    """Propagate ``state0`` over ``[0, t]`` under scipy LSODA. Returns the
    terminal 6-state, or ``None`` on integrator failure."""
    try:
        sol = solve_ivp(
            bcr4bp.bcr4bp_eom,
            (0.0, t),
            np.asarray(state0, dtype=np.float64),
            args=(system,),
            method="LSODA",
            rtol=max(rtol, 1e-12),
            atol=max(atol, 1e-12),
        )
    except (RuntimeError, ValueError):
        return None
    if not sol.success:
        return None
    yf = sol.y[:, -1]
    if not np.all(np.isfinite(yf)):
        return None
    return np.asarray(yf, dtype=np.float64)


def run_v3_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    v2_verdict: V2VerdictBCR4BP,
    n_cycles: int = V3_BCR4BP_N_CYCLES_MIN,
    agreement_floor_km: float = V3_BCR4BP_AGREEMENT_FLOOR_KMS,
    l_km: float = SEM_L_KM,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    notes: str = "",
) -> V3VerdictBCR4BP:
    """Run the V3 independent-integrator BCR4BP gauntlet on a periodic orbit.

    Pipeline:
      For ``k=1..n_cycles``:
        1. Propagate ``k*T`` under DOP853 (V2's integrator) — fresh endpoint.
        2. Propagate ``k*T`` under LSODA (independent multistep family).
        3. Record the terminal-position agreement (km).
      Gate the max agreement against ``agreement_floor_km``. Chains on the V2
      verdict: requires it passed AND has enough laps.

    Parameters
    ----------
    candidate_id :
        Identifier carried into the verdict.
    orbit :
        A :class:`BCR4BPPeriodicOrbit` from the genome.
    v2_verdict :
        The V2 verdict this V3 chains on. Must have passed and carry at least
        ``n_cycles`` per-cycle drift entries.
    n_cycles :
        Number of laps to cross-check. ``>= V3_BCR4BP_N_CYCLES_MIN`` and
        ``<= len(v2_verdict.per_cycle_drift_km)``.
    agreement_floor_km :
        Per-lap integrator-agreement floor (default
        :data:`V3_BCR4BP_AGREEMENT_FLOOR_KMS`).
    l_km :
        Length unit (km) for the agreement conversion (default :data:`SEM_L_KM`).
    rtol, atol :
        Integrator tolerances for BOTH arms.
    notes :
        Free-form audit note.

    Returns
    -------
    V3VerdictBCR4BP
        ``passes_v3_bcr4bp`` is the headline boolean.

    Notes
    -----
    A V3 PASS does NOT admit to the catalogue. The independent integrator is a
    different scipy family (LSODA multistep) integrating the exact same
    ``bcr4bp_eom`` — see the module docstring for why not REBOUND IAS15.
    """
    if not isinstance(orbit, BCR4BPPeriodicOrbit):
        raise TypeError(f"orbit must be a BCR4BPPeriodicOrbit; got {type(orbit).__name__}")
    if not isinstance(v2_verdict, V2VerdictBCR4BP):
        raise TypeError(f"v2_verdict must be a V2VerdictBCR4BP; got {type(v2_verdict).__name__}")
    if n_cycles < V3_BCR4BP_N_CYCLES_MIN:
        raise ValueError(f"V3_bcr4bp requires n_cycles >= {V3_BCR4BP_N_CYCLES_MIN}; got {n_cycles}")
    if len(v2_verdict.per_cycle_drift_km) < n_cycles:
        raise ValueError(
            f"V3 chains on V2: v2_verdict has {len(v2_verdict.per_cycle_drift_km)} laps "
            f"< requested n_cycles={n_cycles}"
        )
    if agreement_floor_km <= 0.0:
        raise ValueError(f"agreement_floor_km must be > 0; got {agreement_floor_km}")

    state0 = np.asarray(orbit.state_initial, dtype=np.float64)
    period = float(orbit.period_nondim)

    per_cycle: list[float] = []
    converged = True
    for k in range(1, n_cycles + 1):
        try:
            arc = bcr4bp.propagate_bcr4bp(
                orbit.system, state0, k * period, with_stm=False, rtol=rtol, atol=atol
            )
        except RuntimeError:
            converged = False
            break
        pos_dop = arc.state_f[:3]
        yf_lsoda = _propagate_lsoda(orbit.system, state0, k * period, rtol=rtol, atol=atol)
        if yf_lsoda is None or not np.all(np.isfinite(pos_dop)):
            converged = False
            break
        agree_km = float(np.linalg.norm(yf_lsoda[:3] - pos_dop)) * l_km
        per_cycle.append(agree_km)

    max_agree = max(per_cycle) if per_cycle else float("inf")
    converged_each_cycle = bool(converged and len(per_cycle) == n_cycles)
    passes = bool(
        converged_each_cycle
        and n_cycles >= V3_BCR4BP_N_CYCLES_MIN
        and max_agree <= agreement_floor_km
        and v2_verdict.passes_v2_bcr4bp
    )

    return V3VerdictBCR4BP(
        candidate_id=candidate_id,
        n_cycles_requested=int(n_cycles),
        per_cycle_agreement_km=tuple(per_cycle),
        max_agreement_km=float(max_agree),
        agreement_floor_km=float(agreement_floor_km),
        n_cycles_min=V3_BCR4BP_N_CYCLES_MIN,
        independent_integrator="LSODA",
        converged_each_cycle=converged_each_cycle,
        passes_v3_bcr4bp=passes,
        notes=notes,
    )


__all__ = [
    "V3_BCR4BP_AGREEMENT_FLOOR_KMS",
    "V3_BCR4BP_N_CYCLES_MIN",
    "V3VerdictBCR4BP",
    "run_v3_bcr4bp",
]
