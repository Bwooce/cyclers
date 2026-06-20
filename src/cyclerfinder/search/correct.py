"""N-arc ballistic differential corrector on the real ephemeris (spec §2.1).

Generalises scripts/correct_s1l1_twoarc.py: free vars x = [t0, leg ToFs] with
one leg pinned by the sourced period; residuals = flyby V_inf-magnitude
continuity + periodicity closure, driven to zero with least_squares. Two
residual modes (task #122):

* ``residual_mode="magnitude"`` (DEFAULT, unchanged) -- magnitude continuity
  only; bend feasibility checked post-hoc, never in the residual.
* ``residual_mode="vector"`` -- the Jones-method full v-inf vector residual:
  magnitude continuity PLUS a per-flyby bend-feasibility hinge INSIDE the
  residual, so the solve is steered toward bend-feasible (ballistic) families
  instead of the magnitude-continuous-but-powered basin.

Pure: depends only on core/lambert, core/ephemeris, core/constants.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.constants import MU_SUN_KM3_S2, PLANETS, VINF_CEILING_KMS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    LambertSolution,
    lambert,
)
from cyclerfinder.core.satellites import SATELLITES

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
    # LOUD physics flag (never a filter): any encounter V_inf exceeds the
    # elliptic-periodicity ceiling v_esc_sun(r_B) + v_B (~71.9 km/s at Earth).
    # A True here marks a degenerate / unit-error / off-family solve — a
    # periodic heliocentric cycler is physically incapable of it. Legitimate
    # high-energy rows (Russell 20.3 km/s) leave this False; it is surfaced, not
    # acted on, so a caller decides how to handle it (the publication layer
    # refuses; the search loop merely records it).
    hyperbolic_impossible: bool = False

    @property
    def constraints_satisfied(self) -> bool:
        return self.converged and self.bend_feasible and self.vinf_cap_ok


def _hyperbolic_impossible(
    sequence: tuple[str, ...], vinf_per_encounter_kms: tuple[float, ...]
) -> bool:
    """True if ANY encounter V_inf breaches its body's elliptic-periodicity
    ceiling (:data:`~cyclerfinder.core.constants.VINF_CEILING_KMS`).

    The per-encounter magnitudes align positionally with ``sequence``. A body
    absent from the ceiling table is skipped (cannot be assessed), never treated
    as a breach.
    """
    for body, vinf in zip(sequence, vinf_per_encounter_kms, strict=False):
        ceiling = VINF_CEILING_KMS.get(body)
        if ceiling is not None and vinf > ceiling:
            return True
    return False


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
    mu_central: float = MU_SUN_KM3_S2,
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
        sols = lambert(r1, r2, tofs[i] * DAY_S, mu=mu_central, max_revs=per_leg_revs[i])
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
    mode: str = "magnitude",
    sequence: tuple[str, ...] | None = None,
    rp_factors: dict[str, float] | None = None,
) -> list[float]:
    """Ballistic-closure residuals (spec §2.1, ``correct_s1l1_twoarc.py:96-106``).

    ``mode="magnitude"`` (DEFAULT, unchanged): for each intermediate encounter
    ``Bi`` (``1 <= i <= n-2``) a flyby conserves V_inf magnitude:
    ``|V_inf_in(Bi)| - |V_inf_out(Bi)|``. Plus the periodicity closure term
    ``|V_inf_in(Bn-1)| - |V_inf_out(B0)|``. This is necessary for ballistic
    closure but blind to whether the required in->out rotation is achievable.

    ``mode="vector"`` (Jones-method, task #122 Phase 1): bend feasibility is
    moved INSIDE the residual. Each intermediate flyby contributes TWO terms:
      1. the magnitude-continuity term (as above), and
      2. a feasibility hinge ``max(0, required_bend - max_bend)`` in degrees,
         scaled to km/s-comparable units -- the part of the required rotation
         that exceeds the body's V_inf-limited single-flyby turn. A
         magnitude-continuous-but-over-bent chain (the #110/#120 powered basin)
         now carries a non-zero residual, steering least_squares toward the
         bend-feasible (Jones) family. The closure term is magnitude-only (the
         wrap encounter is the home body, not an intermediate flyby).

    ``mode="vector"`` requires ``sequence`` (to look up each body's max bend).
    """
    norm = np.linalg.norm
    if mode == "magnitude":
        res: list[float] = []
        for i in range(1, n_encounters - 1):
            res.append(
                float(norm(np.asarray(nodes[f"b{i}_in"])))
                - float(norm(np.asarray(nodes[f"b{i}_out"])))
            )
        last = n_encounters - 1
        res.append(
            float(norm(np.asarray(nodes[f"b{last}_in"]))) - float(norm(np.asarray(nodes["b0_out"])))
        )
        return res

    if mode != "vector":
        raise ValueError(f"unknown residual mode {mode!r}")
    if sequence is None:
        raise ValueError("vector residual mode requires the encounter sequence")

    # km/s per degree of infeasible bend: the feasibility hinge is in degrees;
    # express it on the same scale as the magnitude (km/s) terms so least_squares
    # weighs them comparably. A scale of |V_inf|/deg-of-max-bend would be exact;
    # the V_inf magnitude itself is a stable, well-conditioned proxy (a fully
    # infeasible 180-deg turn at |V_inf| then costs ~|V_inf| km/s of residual,
    # the same order as a total magnitude loss).
    res = []
    for i in range(1, n_encounters - 1):
        v_in = np.asarray(nodes[f"b{i}_in"])
        v_out = np.asarray(nodes[f"b{i}_out"])
        mag_in = float(norm(v_in))
        mag_out = float(norm(v_out))
        res.append(mag_in - mag_out)
        required = _bend_deg(v_in, v_out)
        max_turn = _max_bend_deg(mag_in, sequence[i], rp_factors)
        excess_deg = max(0.0, required - max_turn)
        # Hinge -> km/s: fraction of a half-turn (180 deg) scaled by |V_inf|.
        res.append(mag_in * excess_deg / 180.0)
    last = n_encounters - 1
    res.append(
        float(norm(np.asarray(nodes[f"b{last}_in"]))) - float(norm(np.asarray(nodes["b0_out"])))
    )
    return res


def _residuals(
    x: Sequence[float] | np.ndarray,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    slack_leg: int,
    period_days: float,
    ephem: Ephemeris,
    residual_mode: str = "magnitude",
    rp_factors: dict[str, float] | None = None,
    mu_central: float = MU_SUN_KM3_S2,
) -> list[float]:
    """least_squares residual callable: V_inf-continuity + closure, with Lambert
    pathologies mapped to a large finite penalty (``correct_s1l1_twoarc.py:100``).

    ``residual_mode`` selects the magnitude-only (default) or the full v-inf
    vector (bend-feasibility-aware) residual (task #122 Phase 1). In vector mode
    each intermediate flyby contributes two residual terms, so the penalty
    length on a Lambert pathology grows accordingly.
    """
    n_encounters = len(sequence)
    # FRAGILITY (task #122/#137 review I2): ``n_res`` is the penalty-path residual
    # length returned on a Lambert pathology (the ``[1e3] * n_res`` branch below).
    # It MUST stay byte-for-byte equal to ``len(_residual_vector(..., mode=...))``
    # for the SAME mode, or least_squares sees a ragged residual vector (the
    # pathology evaluations and the successful ones disagree on length, which
    # corrupts the Jacobian / raises). The two formulas below mirror
    # ``_residual_vector`` exactly: magnitude -> ``(n-2)`` intermediate continuity
    # terms + 1 closure = ``n-1``; vector -> each of the ``(n-2)`` intermediate
    # flybys emits TWO terms (continuity + bend hinge) + 1 closure = ``2*(n-2)+1``.
    # Any change to the residual layout in either mode MUST be made in BOTH places.
    n_res = 2 * (n_encounters - 2) + 1 if residual_mode == "vector" else n_encounters - 1
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
    return _residual_vector(
        nodes,
        n_encounters=n_encounters,
        mode=residual_mode,
        sequence=sequence,
        rp_factors=rp_factors,
    )


def _max_bend_deg(vinf_kms: float, body: str, rp_factors: dict[str, float] | None = None) -> float:
    """Maximum single-flyby turn angle (deg) for ``vinf_kms`` at ``body``
    (``correct_s1l1_twoarc.py:109-114``). ``rp_factors`` optionally scales the
    body's ``safe_alt_km`` (spec §2.1 ``r_p_safe``).

    Resolves a heliocentric planet code via ``PLANETS`` and, failing that, a
    moon code via ``SATELLITES`` (moon-tour Tier-1): both expose
    ``mu_km3_s2``/``radius_eq_km``/``safe_alt_km``, so the V_inf-limited turn
    formula is centre-blind once the right body record is found."""
    pl = PLANETS.get(body) or SATELLITES[body]
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


def _per_encounter_vinf(nodes: dict[str, np.ndarray], n_encounters: int) -> tuple[float, ...]:
    """Per-encounter V_inf magnitude (mean of the in/out legs available at each
    encounter; ends carry only one leg). ``correct_s1l1_twoarc.py:139-140``."""
    norm = np.linalg.norm
    out: list[float] = []
    for i in range(n_encounters):
        mags = [float(norm(nodes[key])) for key in (f"b{i}_in", f"b{i}_out") if key in nodes]
        out.append(float(np.mean(mags)) if mags else 0.0)
    return tuple(out)


def _bend_feasible(
    nodes: dict[str, np.ndarray],
    sequence: tuple[str, ...],
    rp_factors: dict[str, float] | None,
) -> bool:
    """Every intermediate flyby's required turn must fit within its V_inf-limited
    maximum (``correct_s1l1_twoarc.py:141,151``)."""
    norm = np.linalg.norm
    for i in range(1, len(sequence) - 1):
        v_in = nodes[f"b{i}_in"]
        v_out = nodes[f"b{i}_out"]
        required = _bend_deg(v_in, v_out)
        max_turn = _max_bend_deg(float(norm(v_in)), sequence[i], rp_factors)
        if required > max_turn:
            return False
    return True


def ballistic_correct(
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    t0_seed_sec: float,
    tof_seed_days: Sequence[float],
    period_sec: float,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    rp_factors: dict[str, float] | None = None,
    slack_leg: int | None = None,
    tol_kms: float = 0.1,
    residual_mode: str = "magnitude",
    method: Literal["trf", "dogbox", "lm"] = "lm",
    mu_central: float = MU_SUN_KM3_S2,
) -> BallisticClosureResult:
    """N-arc ballistic differential corrector (spec §2.1; generalises
    ``correct_s1l1_twoarc.py:_solve``).

    Free vars ``x = [t0_sec, *tof_seed_days]`` (the slack leg is eliminated and
    reconstructed as ``period - sum(free legs)``). Drives the V_inf-continuity +
    closure residuals to zero with ``least_squares(method=method)``; converged iff
    the max residual is below ``tol_kms`` (default 0.1, the prototype threshold
    ``correct_s1l1_twoarc.py:169``). Bend feasibility and the V_inf cap are
    evaluated post-hoc.

    ``method`` selects the least_squares algorithm; the default ``"lm"``
    (Levenberg-Marquardt) is the historical, well-conditioned choice when the
    residual count exceeds the free-var count (m > n). Callers whose problem can
    be under-determined (m <= n; e.g. a short 2-encounter chain in vector mode)
    must pass ``method="trf"``, which handles m<n, m=n and m>n. ``lm`` raises a
    ``ValueError`` for m<n.
    """
    period_days = period_sec / DAY_S
    n_encounters = len(sequence)
    if slack_leg is None:
        # Default: pin the longest seed leg (most slack to absorb the period).
        slack_leg = int(np.argmax(tof_seed_days)) if len(tof_seed_days) else 0

    def _res(x: np.ndarray) -> list[float]:
        return _residuals(
            x,
            sequence=sequence,
            per_leg_revs=per_leg_revs,
            per_leg_branch=per_leg_branch,
            slack_leg=slack_leg,
            period_days=period_days,
            ephem=ephem,
            residual_mode=residual_mode,
            rp_factors=rp_factors,
            mu_central=mu_central,
        )

    x0 = np.array([t0_seed_sec, *tof_seed_days], dtype=np.float64)
    sol = least_squares(_res, x0, method=method, max_nfev=80, xtol=1e-9, ftol=1e-9)
    x = sol.x
    # least_squares stores the residual vector evaluated at the returned x in
    # sol.fun (LM's final fvec); reuse it instead of a redundant full re-eval of
    # the residual (each call solves every leg's Lambert). Byte-identical to
    # _res(sol.x): same deterministic residual at the same point.
    res = [float(r) for r in sol.fun]
    max_res = max(abs(r) for r in res)

    full_tofs = _reconstruct_tofs(tuple(float(v) for v in x[1:]), slack_leg, period_days)

    # Post-solve node extraction is NOT inside the residual guard, so a converged
    # x that still lands a Lambert pathology (e.g. a multi-rev VEM leg the solver
    # walked into) would raise here. Treat that as a non-converged outcome —
    # surface honestly rather than crash the caller (the headline-gate finding
    # depends on a result object, not an exception).
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
            mu_central=mu_central,
        )
    except (LambertConvergenceError, LambertGeometryError, ValueError):
        return BallisticClosureResult(
            t0_sec=float(x[0]),
            tof_days=tuple(full_tofs),
            max_residual_kms=float("inf"),
            vinf_per_encounter_kms=(),
            converged=False,
            bend_feasible=False,
            vinf_cap_ok=False,
        )
    vinf_per_encounter = _per_encounter_vinf(nodes, n_encounters)

    converged = max_res < tol_kms
    bend_feasible = _bend_feasible(nodes, sequence, rp_factors)
    vinf_cap_ok = max(vinf_per_encounter) <= vinf_cap

    return BallisticClosureResult(
        t0_sec=float(x[0]),
        tof_days=tuple(full_tofs),
        max_residual_kms=float(max_res),
        vinf_per_encounter_kms=vinf_per_encounter,
        converged=converged,
        bend_feasible=bend_feasible,
        vinf_cap_ok=vinf_cap_ok,
        hyperbolic_impossible=_hyperbolic_impossible(sequence, vinf_per_encounter),
    )
