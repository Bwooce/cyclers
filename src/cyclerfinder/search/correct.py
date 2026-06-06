"""N-arc ballistic differential corrector on the real ephemeris (spec §2.1).

Generalises scripts/correct_s1l1_twoarc.py: free vars x = [t0, leg ToFs] with
one leg pinned by the sourced period; residuals = flyby V_inf-magnitude
continuity + periodicity closure, driven to zero with least_squares; bend
feasibility checked post-hoc, never in the residual. Pure: depends only on
core/lambert, core/ephemeris, core/constants.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    LambertSolution,
    lambert,
)

DAY_S = 86400.0


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


# S1L1-prototype aliases for the first two-arc chain (E-M-E-E). The generic
# ``b{i}_in/b{i}_out`` keys are authoritative; these convenience aliases keep the
# prototype's residual/test vocabulary working unchanged.
_S1L1_ALIASES: dict[str, str] = {
    "e0": "b0_out",
    "m_in": "b1_in",
    "m_out": "b1_out",
    "e1_in": "b2_in",
    "e1_out": "b2_out",
    "e2_in": "b3_in",
}


def _pick(sols: list[LambertSolution], n_revs: int, branch: str) -> LambertSolution:
    """Select the requested (n_revs, branch) Lambert solution; fall back to the
    first (single-rev) solution if the exact branch is absent (prototype
    ``correct_s1l1_twoarc.py:48-52``)."""
    for s in sols:
        if s.n_revs == n_revs and s.branch == branch:
            return s
    return sols[0]


def _reconstruct_tofs(
    free_tof_days: Sequence[float], slack_leg: int, period_days: float
) -> list[float]:
    """Re-insert the eliminated slack leg ToF (``period - sum(free legs)``)."""
    slack_tof = period_days - float(sum(free_tof_days))
    tofs = list(free_tof_days)
    tofs.insert(slack_leg, slack_tof)
    return tofs


def _vinf_nodes(
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    t0_sec: float,
    free_tof_days: Sequence[float],
    slack_leg: int,
    period_days: float,
    ephem: Ephemeris,
) -> dict[str, np.ndarray]:
    """Per-encounter V_inf vectors for a closed N-arc chain (spec §2.1).

    Generalises ``_legs`` / ``_state_vinf`` (``correct_s1l1_twoarc.py:55-93``):
    reconstruct the slack leg ToF, walk cumulative encounter epochs, solve each
    leg's Lambert with its ``(n_revs, branch)``, and return ``V_inf =
    v_sc - v_planet`` for the outbound leg of every encounter (``b{i}_out``) and
    the inbound leg of every encounter (``b{i}_in``). End nodes ``b0`` and
    ``bn`` carry only the closure pair (``b0_out`` / ``bn_in``).
    """
    tofs = _reconstruct_tofs(free_tof_days, slack_leg, period_days)
    n_legs = len(sequence) - 1
    if len(tofs) != n_legs:
        raise ValueError(f"expected {n_legs} leg ToFs, got {len(tofs)}")

    # Cumulative epoch (seconds) at each encounter.
    epochs = [t0_sec]
    for tof in tofs:
        epochs.append(epochs[-1] + tof * DAY_S)

    # Heliocentric body states at each encounter.
    states = [ephem.state(body, t) for body, t in zip(sequence, epochs, strict=True)]

    nodes: dict[str, np.ndarray] = {}
    for i in range(n_legs):
        r1, v1_pl = states[i]
        r2, v2_pl = states[i + 1]
        sols = lambert(r1, r2, tofs[i] * DAY_S, max_revs=per_leg_revs[i])
        sol = _pick(sols, per_leg_revs[i], per_leg_branch[i])
        # V_inf leaving encounter i and arriving at encounter i+1.
        nodes[f"b{i}_out"] = np.asarray(sol.v1) - np.asarray(v1_pl)
        nodes[f"b{i + 1}_in"] = np.asarray(sol.v2) - np.asarray(v2_pl)

    for alias, key in _S1L1_ALIASES.items():
        if key in nodes:
            nodes[alias] = nodes[key]
    return nodes


def _residual_vector(
    nodes: dict[str, np.ndarray] | dict[str, tuple[float, float, float]],
    *,
    n_encounters: int,
) -> list[float]:
    """Ballistic-closure residuals (spec §2.1, ``correct_s1l1_twoarc.py:96-106``).

    For each intermediate encounter ``Bi`` (``1 <= i <= n-2``) a flyby conserves
    V_inf magnitude: ``|V_inf_in(Bi)| - |V_inf_out(Bi)|``. Plus the periodicity
    closure term ``|V_inf_in(Bn-1)| - |V_inf_out(B0)|``.
    """
    norm = np.linalg.norm
    res: list[float] = []
    for i in range(1, n_encounters - 1):
        res.append(
            float(norm(np.asarray(nodes[f"b{i}_in"]))) - float(norm(np.asarray(nodes[f"b{i}_out"])))
        )
    last = n_encounters - 1
    res.append(
        float(norm(np.asarray(nodes[f"b{last}_in"]))) - float(norm(np.asarray(nodes["b0_out"])))
    )
    return res


def _residuals(
    x: Sequence[float],
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    slack_leg: int,
    period_days: float,
    ephem: Ephemeris,
) -> list[float]:
    """least_squares residual callable: V_inf-continuity + closure, with Lambert
    pathologies mapped to a large finite penalty (``correct_s1l1_twoarc.py:100``).
    """
    n_encounters = len(sequence)
    n_res = n_encounters - 1
    try:
        nodes = _vinf_nodes(
            sequence=sequence,
            per_leg_revs=per_leg_revs,
            per_leg_branch=per_leg_branch,
            t0_sec=float(x[0]),
            free_tof_days=tuple(float(v) for v in x[1:]),
            slack_leg=slack_leg,
            period_days=period_days,
            ephem=ephem,
        )
    except (LambertConvergenceError, LambertGeometryError, ValueError):
        return [1e3] * n_res
    return _residual_vector(nodes, n_encounters=n_encounters)


def _max_bend_deg(vinf_kms: float, body: str, rp_factors: dict[str, float] | None = None) -> float:
    """Maximum single-flyby turn angle (deg) for ``vinf_kms`` at ``body``
    (``correct_s1l1_twoarc.py:109-114``). ``rp_factors`` optionally scales the
    body's ``safe_alt_km`` (spec §2.1 ``r_p_safe``)."""
    pl = PLANETS[body]
    mu = pl.mu_km3_s2
    safe_alt = pl.safe_alt_km
    if rp_factors is not None and body in rp_factors:
        safe_alt *= rp_factors[body]
    r_p = pl.radius_eq_km + safe_alt
    e = 1.0 + r_p * vinf_kms * vinf_kms / mu
    return float(np.degrees(2.0 * np.arcsin(1.0 / e)))


def _bend_deg(v_in: Sequence[float] | np.ndarray, v_out: Sequence[float] | np.ndarray) -> float:
    """Angle (deg) between the in/out V_inf vectors (``correct_s1l1_twoarc.py:117-120``)."""
    a, b = np.asarray(v_in), np.asarray(v_out)
    c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, c)))))
