"""One-DSM-per-leg genome (the Takao eta-coordinate) — multi-arc-return primitive.

This is the architectural sibling of :mod:`cyclerfinder.search.free_return` (the
#137 radial-crossing genome). Where the free-return module models an
Earth->Mars->Earth transfer as a *single* heliocentric ellipse (one shape ``(a, e)``
shared by both legs), this module supplies the primitive the single-ellipse genome
provably lacks: an *interior impulse* on each leg. That interior impulse is what
lets a leg follow a different ballistic arc on its front fraction than on its back
fraction — i.e. a leg that is NOT a piece of one repeating ellipse. The S1L1 /
Jones multi-arc closure blocker (``docs/notes/multi-arc-classification.md`` §7/§12;
the #137 / MBH Gate-3 negative result) is exactly the case where the sourced
geometry is *two generic-return arcs*, not one ellipse, so a single-ellipse genome
has no sourced-anchor basin. The one-DSM-per-leg genome can represent that.

The transcription
-----------------
Y. Takao, "Mission Analysis for the First-Ever Saturn Trojan 2019 UO14,"
arXiv:2501.06586 (astro-ph.EP), 2025. Per-leg genome
``Y_i = [V_inf, alpha, beta, tau, eta]_i`` (Eqs.1-2). The DSM leg evaluator
(Eqs.6-7, transcribed in ``docs/notes/2026-06-07-takao-2025-mpga-1dsm-mining.md``):

* propagate the outgoing state ballistically (heliocentric 2-body) for ``eta*tau``,
  giving the DSM position ``r_12`` and velocity ``v_12`` (the velocity arriving at
  the DSM, BEFORE the impulse);
* solve Lambert from the DSM position to the next body over the remaining
  ``(1-eta)*tau``, giving the post-impulse departure ``v_21`` and the arrival
  ``v_22``;
* the DSM impulse magnitude is the velocity mismatch at the DSM point,
  ``dV_DSM = ||v_21 - v_12||`` (Eq.6);
* the incoming hyperbolic velocity at the next body is ``v_22 - v_planet`` (Eq.7).

``eta in [0, 1]`` is the only genuinely new genome coordinate vs the free-return /
ballistic correctors; everything else (2-body propagate, Lambert) is reused from
:mod:`cyclerfinder.core`.

Conventions mirrored from ``free_return.py``
--------------------------------------------
* A frozen dataclass result carrying both the CONSTRAINED quantity (the residual /
  total dV the corrector drives) and the FREE / EMERGED evidence (per-leg dV
  breakdown, per-body emerged V_inf). The sourced V_inf is NEVER imposed; it emerges
  from the converged genome and is comparable as evidence (the golden-rule
  separation).
* ``converged`` decided by residual MAGNITUDE alone (``< tol_kms``), by design; the
  ``least_squares`` ``success`` flag is kept only for the audit trail.
* A full audit trail (per-leg dV, DSM states, solver nfev) on the result.
* RNG-free / deterministic: no randomness anywhere in this module (MBH supplies the
  hops; this module only evaluates and corrects).

Pure: depends only on core/constants, core/kepler (propagate), core/lambert.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import KeplerError, propagate
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.mbh import MBHStep

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64
DAY_S = SECONDS_PER_DAY


# ---------------------------------------------------------------------------
# The leg primitive (Takao Eq.6-7)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DsmLegResult:
    """Outcome of a single propagate-then-Lambert (one-DSM) leg.

    Attributes
    ----------
    v_arrive:
        Heliocentric velocity at the target endpoint (``v_22`` in Takao Eq.6),
        km/s. The post-flyby caller forms the incoming V_inf as
        ``v_arrive - v_planet_target`` (Eq.7).
    v_depart_post_dsm:
        Heliocentric departure velocity from the DSM point AFTER the impulse
        (``v_21``), km/s. The Lambert start velocity.
    v_arrive_pre_dsm:
        Heliocentric velocity arriving AT the DSM point BEFORE the impulse
        (``v_12``), km/s. The ballistic propagation endpoint velocity.
    dv_dsm_kms:
        The DSM impulse magnitude ``||v_21 - v_12||`` (Eq.6), km/s.
    r_dsm:
        Heliocentric position of the DSM point (``r_12 = r_21``), km. Audit/viz.
    t_dsm_sec:
        Time of the DSM along the leg, ``eta * tof`` seconds (relative to leg
        start). Audit/viz.
    n_revs_chosen:
        The revolution count of the back-arc Lambert branch actually selected,
        ``0`` for the single-revolution (direct) branch. When ``max_revs == 0``
        this is always ``0`` (the legacy single-rev path); when ``max_revs > 0``
        it records which multi-rev branch minimised :attr:`dv_dsm_kms`. Audit.
    branch_chosen:
        The branch label of the selected back-arc Lambert solution
        (``"single"`` / ``"low"`` / ``"high"``). Audit; pairs with
        :attr:`n_revs_chosen`.
    """

    v_arrive: Vec3
    v_depart_post_dsm: Vec3
    v_arrive_pre_dsm: Vec3
    dv_dsm_kms: float
    r_dsm: Vec3
    t_dsm_sec: float
    n_revs_chosen: int = 0
    branch_chosen: str = "single"


# eta is clamped this far from the singular endpoints {0, 1}: at exactly eta=0 the
# ballistic front arc has zero duration (DSM == leg start) and at exactly eta=1 the
# Lambert back arc has zero duration (singular t->0 Lambert). The degeneracy gate
# probes eta APPROACHING these endpoints (e.g. 1e-4), which is well-posed; this
# guard only rejects the exact-endpoint singular calls.
_ETA_EPS: float = 1.0e-9


def dsm_leg(
    r0: Vec3,
    v0: Vec3,
    tof: float,
    eta: float,
    target_r: Vec3,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
    max_revs: int = 0,
    rev_branch: tuple[int, str] | None = None,
) -> DsmLegResult:
    """One interior-impulse (DSM) leg: propagate ``eta*tof`` then Lambert the rest.

    Takao Eq.6-7. Starting from the outgoing heliocentric state ``(r0, v0)`` at the
    leg's first body, the trajectory is propagated ballistically (heliocentric
    2-body) for ``eta*tof`` to the DSM point; Lambert is then solved from the DSM
    position to ``target_r`` over the remaining ``(1-eta)*tof``. The DSM impulse is
    the velocity mismatch at the DSM point.

    Parameters
    ----------
    r0, v0:
        Heliocentric inertial state at the leg start, ``(3,)`` float64, km and km/s.
    tof:
        Total leg time of flight, seconds. Must be strictly positive.
    eta:
        DSM timing fraction in ``[0, 1]`` (Takao). ``eta`` of the ToF is flown
        ballistically before the impulse; ``(1-eta)`` is the Lambert arc. Must be
        strictly inside ``(0, 1)`` up to :data:`_ETA_EPS` (the exact endpoints are
        singular -- the degeneracy gate probes eta APPROACHING them).
    target_r:
        Heliocentric position of the target body at leg arrival, ``(3,)`` km.
    mu:
        Central-body gravitational parameter, km^3/s^2 (heliocentric default).
    prograde:
        Lambert transfer sense for the back arc (default prograde / short-way).
    max_revs:
        Maximum number of full revolutions to consider for the back-arc Lambert
        (passed straight through to :func:`cyclerfinder.core.lambert.lambert`).
        The default ``0`` keeps the historical single-revolution-only behaviour
        (bit-identical results). When ``max_revs > 0`` the solver enumerates the
        single-rev branch plus every feasible multi-rev ``low``/``high`` branch
        up to ``max_revs`` and selects the one that MINIMISES the DSM impulse
        ``dV_DSM`` (the leg objective). The resonant loop arcs of a multi-arc
        cycler are >1-period transfers, for which the single-rev branch is forced
        onto a degenerate near-radial high-energy solution; allowing multi-rev
        branches is what makes those arcs representable (the #153 diagnosis).
    rev_branch:
        Optional ``(n_revs, branch)`` selector. When given, the back arc uses
        exactly that branch (e.g. ``(2, "low")``) instead of minimising over the
        enumerated set; ``max_revs`` is widened to ``n_revs`` if needed so the
        requested branch is produced. A :class:`LambertError` is raised if the
        requested branch is infeasible for this geometry/ToF. ``None`` (default)
        selects the dV-minimising branch as described above.

    Returns
    -------
    DsmLegResult
        Arrival velocity, post/pre-DSM velocities, the DSM impulse magnitude, the
        DSM state, and the chosen back-arc revolution/branch for the audit trail.

    Raises
    ------
    ValueError
        On non-positive ``tof`` or ``eta`` at/outside the singular endpoints, or a
        negative ``max_revs``.
    LambertError
        If the back-arc Lambert solve fails (degenerate geometry / no
        convergence), or the explicitly requested ``rev_branch`` is infeasible.
    """
    if tof <= 0.0:
        raise ValueError(f"tof must be positive, got {tof}")
    if not (_ETA_EPS <= eta <= 1.0 - _ETA_EPS):
        raise ValueError(
            f"eta must lie in [{_ETA_EPS}, {1.0 - _ETA_EPS}] (the exact endpoints "
            f"0/1 are singular ballistic/Lambert degeneracies), got {eta}"
        )
    if max_revs < 0:
        raise ValueError(f"max_revs must be non-negative, got {max_revs}")

    r0_arr = np.asarray(r0, dtype=np.float64)
    v0_arr = np.asarray(v0, dtype=np.float64)
    target_arr = np.asarray(target_r, dtype=np.float64)

    t_front = eta * tof
    t_back = (1.0 - eta) * tof

    # Front arc: ballistic 2-body propagation to the DSM point (Eq.6, v_12).
    r_dsm, v12 = propagate(r0_arr, v0_arr, t_front, mu)

    # Back arc: Lambert from the DSM position to the target over (1-eta)*tof.
    # max_revs=0 returns ONLY the single-rev (direct) branch -> the historical
    # path is bit-identical (sols[0] is the single-rev solution). max_revs>0 also
    # returns the feasible multi-rev low/high branches, among which the chosen
    # branch is the one minimising the DSM impulse (the resonant loop arcs need
    # this; the single-rev branch on a >1-period leg is degenerate -- #153).
    enumerate_revs = max_revs
    if rev_branch is not None:
        # Widen the enumeration so the explicitly requested branch is produced.
        enumerate_revs = max(max_revs, int(rev_branch[0]))
    sols = lambert(r_dsm, target_arr, t_back, mu=mu, prograde=prograde, max_revs=enumerate_revs)

    if rev_branch is not None:
        want_n, want_b = int(rev_branch[0]), str(rev_branch[1])
        chosen = next(
            (s for s in sols if s.n_revs == want_n and s.branch == want_b),
            None,
        )
        if chosen is None:
            raise LambertError(
                f"requested rev_branch ({want_n}, {want_b!r}) is infeasible for this "
                f"geometry/ToF (available: "
                f"{[(s.n_revs, s.branch) for s in sols]})"
            )
    elif max_revs == 0:
        # Legacy path: the single-rev (direct) branch, identical to before.
        chosen = sols[0]
    else:
        # Pick the branch minimising the DSM impulse ||v21 - v12||.
        chosen = min(sols, key=lambda s: float(np.linalg.norm(s.v1 - v12)))

    v21 = chosen.v1  # post-impulse departure from the DSM (Eq.6, v_21)
    v22 = chosen.v2  # arrival at the target (Eq.6, v_22)

    dv_dsm = float(np.linalg.norm(v21 - v12))

    return DsmLegResult(
        v_arrive=np.asarray(v22, dtype=np.float64),
        v_depart_post_dsm=np.asarray(v21, dtype=np.float64),
        v_arrive_pre_dsm=np.asarray(v12, dtype=np.float64),
        dv_dsm_kms=dv_dsm,
        r_dsm=np.asarray(r_dsm, dtype=np.float64),
        t_dsm_sec=float(t_front),
        n_revs_chosen=int(chosen.n_revs),
        branch_chosen=str(chosen.branch),
    )


# ---------------------------------------------------------------------------
# Chained multi-leg evaluator + sequence-keyed bounds
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DsmChainResult:
    """Outcome + audit of a chained one-DSM-per-leg evaluation/correction.

    Mirrors :class:`cyclerfinder.search.free_return.FreeReturnClosureResult`'s
    constraint-vs-evidence separation.

    Attributes
    ----------
    total_dv_kms:
        The CONSTRAINED objective: ``sum(dV_DSM) + dV_arrive`` (Takao Eq.15 without
        the powered-flyby term, which the catalogue's flyby surrogate owns), km/s.
        This is the residual the corrector / MBH drives.
    max_residual_kms:
        Alias of ``total_dv_kms`` exposed under the name the MBH adapters expect
        (``free_return`` / ``ballistic`` correctors both report ``max_residual_kms``
        as the objective). Keeps :func:`make_dsm_chain_step` uniform with the
        existing adapters.
    dv_dsm_per_leg_kms:
        EMERGED per-leg DSM impulse magnitudes (the audit breakdown), km/s.
    dv_arrive_kms:
        EMERGED terminal-arrival V_inf magnitude (rendezvous) or 0 (flyby), km/s.
    vinf_in_kms, vinf_out_kms:
        EMERGED per-body incoming / outgoing V_inf magnitudes, keyed by leg index.
        Evidence -- never imposed (Takao Eq.7).
    eta_per_leg:
        The converged DSM fractions (the new genome coordinate), per leg.
    tof_days_per_leg:
        The converged per-leg ToFs (days).
    n_revs_per_leg:
        EMERGED back-arc revolution count selected on each leg (audit). All zeros
        when the chain was evaluated single-rev (``max_revs == 0``); otherwise the
        per-leg multi-rev branch the dV-minimising selection landed on.
    branch_per_leg:
        EMERGED back-arc Lambert branch label per leg (audit), pairs with
        :attr:`n_revs_per_leg`.
    t0_sec:
        Converged departure epoch (seconds).
    vinf_out0_kms, alpha0, beta0:
        The converged departure-V_inf genome (magnitude + azimuth + elevation,
        Takao Eq.5). These ENTER the dV objective through the departure state
        ``v_planet,0 + V_inf_out0`` and are moved by the corrector; the adapter
        echoes them so the landed decision vector is in the seed's coordinates.
    dsm_states:
        Per-leg ``(r_dsm_km, t_dsm_sec)`` for the audit trail / viz.
    converged:
        ``total_dv_kms < tol_kms`` -- residual-magnitude acceptance, BY DESIGN
        (mirrors ``free_return_correct``; the solver flag is secondary).
    solver_success, solver_nfev:
        DIAGNOSTIC ``least_squares`` outcome (audit trail only, never gates).
    """

    total_dv_kms: float
    max_residual_kms: float
    dv_dsm_per_leg_kms: tuple[float, ...]
    dv_arrive_kms: float
    vinf_in_kms: dict[int, float]
    vinf_out_kms: dict[int, float]
    eta_per_leg: tuple[float, ...]
    tof_days_per_leg: tuple[float, ...]
    t0_sec: float
    vinf_out0_kms: float = 0.0
    alpha0: float = 0.0
    beta0: float = 0.0
    n_revs_per_leg: tuple[int, ...] = field(default_factory=tuple)
    branch_per_leg: tuple[str, ...] = field(default_factory=tuple)
    dsm_states: tuple[tuple[Vec3, float], ...] = field(default_factory=tuple)
    converged: bool = False
    solver_success: bool = True
    solver_nfev: int = 0


def _vinf_out_dir(v_inf: float, alpha: float, beta: float) -> Vec3:
    """Outgoing V_inf vector from magnitude + azimuth/elevation (Takao Eq.5)."""
    return np.array(
        [
            v_inf * np.cos(alpha) * np.cos(beta),
            v_inf * np.sin(alpha) * np.cos(beta),
            v_inf * np.sin(beta),
        ],
        dtype=np.float64,
    )


def evaluate_dsm_chain(
    *,
    sequence: tuple[str, ...],
    t0_sec: float,
    vinf_out0_kms: float,
    alpha0: float,
    beta0: float,
    tof_days_per_leg: tuple[float, ...],
    eta_per_leg: tuple[float, ...],
    ephem: Ephemeris,
    mu: float = MU_SUN_KM3_S2,
    rendezvous: bool = False,
    tol_kms: float = 0.1,
    max_revs: int = 0,
) -> DsmChainResult:
    """Evaluate a chained one-DSM-per-leg trajectory (Takao Eqs.3-7, 14-15).

    The genome strings ``len(sequence) - 1`` legs. The departure state at body 0 is
    ``v_planet,0 + V_inf_out0`` with ``V_inf_out0`` from (``vinf_out0_kms``,
    ``alpha0``, ``beta0``) (Eq.4-5). Each leg runs :func:`dsm_leg`; the heliocentric
    arrival velocity is the spacecraft's incoming velocity at the next body. The
    departing velocity for the NEXT leg is that same incoming heliocentric velocity
    (a ballistic flyby preserves heliocentric speed continuity; the powered-flyby
    surrogate, owned by ``core/flyby.py``, is applied separately at scoring time and
    is NOT charged here -- this evaluator's objective is the sum of DSM impulses
    plus the terminal arrival, Eq.15 without the P-FB term).

    Body epochs chain as ``t_j = t0 + sum(tau_i)`` (Eq.3); body states pulled from
    ``ephem`` at those epochs.

    Returns
    -------
    DsmChainResult
        ``total_dv_kms`` is the objective; the per-leg DSM impulses, emerged V_inf,
        and DSM states are the evidence/audit trail.
    """
    n_legs = len(sequence) - 1
    if n_legs < 1:
        raise ValueError("sequence must name at least two bodies (one leg)")
    if len(tof_days_per_leg) != n_legs or len(eta_per_leg) != n_legs:
        raise ValueError(
            f"tof_days_per_leg and eta_per_leg must each have {n_legs} entries "
            f"(one per leg), got {len(tof_days_per_leg)} / {len(eta_per_leg)}"
        )

    r0, v_planet0 = ephem.state(sequence[0], t0_sec)
    v_out0 = _vinf_out_dir(vinf_out0_kms, alpha0, beta0)
    v_depart = np.asarray(v_planet0, dtype=np.float64) + v_out0

    dv_dsm: list[float] = []
    dsm_states: list[tuple[Vec3, float]] = []
    n_revs_legs: list[int] = []
    branch_legs: list[str] = []
    vinf_in: dict[int, float] = {}
    vinf_out: dict[int, float] = {0: float(vinf_out0_kms)}

    t_cursor = t0_sec
    r_curr = np.asarray(r0, dtype=np.float64)
    feasible = True
    for i in range(n_legs):
        tof_s = tof_days_per_leg[i] * DAY_S
        t_arrive = t_cursor + tof_s
        target_body = sequence[i + 1]
        r_target, v_planet_target = ephem.state(target_body, t_arrive)
        try:
            leg = dsm_leg(
                r_curr,
                v_depart,
                tof_s,
                eta_per_leg[i],
                np.asarray(r_target),
                mu=mu,
                max_revs=max_revs,
            )
        except (LambertError, KeplerError, ValueError):
            # A hop into a too-energetic departure can drive the ballistic
            # propagation hyperbolic past Newton convergence, or the back-arc
            # Lambert degenerate; report the chain INFEASIBLE rather than crashing
            # the whole MBH search (the optimiser simply rejects this hop).
            feasible = False
            break
        dv_dsm.append(leg.dv_dsm_kms)
        dsm_states.append((leg.r_dsm, leg.t_dsm_sec))
        n_revs_legs.append(leg.n_revs_chosen)
        branch_legs.append(leg.branch_chosen)
        v_inf_in_vec = leg.v_arrive - np.asarray(v_planet_target, dtype=np.float64)
        vinf_in[i + 1] = float(np.linalg.norm(v_inf_in_vec))

        # Ballistic-flyby heliocentric continuity for the next leg's departure.
        # (The powered-flyby bend cost is the catalogue scorer's job, not this
        # evaluator's; here the next leg simply departs on the arrival velocity.)
        v_depart = np.asarray(leg.v_arrive, dtype=np.float64)
        vinf_out[i + 1] = float(np.linalg.norm(v_inf_in_vec))
        r_curr = np.asarray(r_target, dtype=np.float64)
        t_cursor = t_arrive

    if not feasible:
        return DsmChainResult(
            total_dv_kms=float("inf"),
            max_residual_kms=float("inf"),
            dv_dsm_per_leg_kms=tuple(dv_dsm),
            dv_arrive_kms=float("inf"),
            vinf_in_kms=vinf_in,
            vinf_out_kms=vinf_out,
            eta_per_leg=tuple(eta_per_leg),
            tof_days_per_leg=tuple(tof_days_per_leg),
            t0_sec=float(t0_sec),
            vinf_out0_kms=float(vinf_out0_kms),
            alpha0=float(alpha0),
            beta0=float(beta0),
            n_revs_per_leg=tuple(n_revs_legs),
            branch_per_leg=tuple(branch_legs),
            dsm_states=tuple(dsm_states),
            converged=False,
        )

    dv_arrive = vinf_in[n_legs] if rendezvous else 0.0
    total_dv = float(sum(dv_dsm) + dv_arrive)
    return DsmChainResult(
        total_dv_kms=total_dv,
        max_residual_kms=total_dv,
        dv_dsm_per_leg_kms=tuple(dv_dsm),
        dv_arrive_kms=float(dv_arrive),
        vinf_in_kms=vinf_in,
        vinf_out_kms=vinf_out,
        eta_per_leg=tuple(eta_per_leg),
        tof_days_per_leg=tuple(tof_days_per_leg),
        t0_sec=float(t0_sec),
        vinf_out0_kms=float(vinf_out0_kms),
        alpha0=float(alpha0),
        beta0=float(beta0),
        n_revs_per_leg=tuple(n_revs_legs),
        branch_per_leg=tuple(branch_legs),
        dsm_states=tuple(dsm_states),
        converged=total_dv < tol_kms,
    )


# ---------------------------------------------------------------------------
# Sequence-keyed automatic bounds (Takao Appendix A.1-A.3)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DsmBounds:
    """Box bounds for the DSM-chain decision vector (Takao Appendix A.1-A.3).

    Layout matches :func:`dsm_chain_decision_vector`:
    ``[t0_sec, vinf_out0, alpha0, beta0, *tof_days_per_leg, *eta_per_leg]``.
    """

    lower: NDArray[np.float64]
    upper: NDArray[np.float64]


def _synodic_period_days(inner: str, outer: str) -> float:
    """Synodic period of two bodies (days), from their sidereal periods."""
    p_in = 2.0 * np.pi * np.sqrt((PLANETS[inner].sma_au * AU_KM) ** 3 / MU_SUN_KM3_S2) / DAY_S
    p_out = 2.0 * np.pi * np.sqrt((PLANETS[outer].sma_au * AU_KM) ** 3 / MU_SUN_KM3_S2) / DAY_S
    return float(abs(1.0 / (1.0 / p_in - 1.0 / p_out)))


def _hohmann_tof_days(inner: str, outer: str) -> float:
    """Hohmann (half-period of the transfer ellipse) ToF inner<->outer (days)."""
    a_t_km = 0.5 * (PLANETS[inner].sma_au + PLANETS[outer].sma_au) * AU_KM
    return float(np.pi * np.sqrt(a_t_km**3 / MU_SUN_KM3_S2) / DAY_S)


def sequence_keyed_bounds(
    *,
    sequence: tuple[str, ...],
    t0_window_sec: tuple[float, float],
    vinf_out0_bounds_kms: tuple[float, float] = (1.0, 5.1),
) -> DsmBounds:
    """Automatic box bounds from the body sequence (Takao Appendix A.1-A.3).

    Only the departure epoch window and the departure V_inf magnitude bound need to
    be supplied (Takao: "users need to specify only the epoch and the V_inf
    magnitude at launch"); the rest are sequence-keyed:

    * ``alpha in [-pi, pi]``, ``beta in [-pi/2, pi/2]``, ``eta in [0, 1]`` (A.1).
    * Per-leg ToF: an inner pair (both bodies inner planets, here the E<->M class)
      gets ``[30 d, P_s + P_H]`` (A.2); any leg touching an outer body gets
      ``[0.3 P_H, 1.3 P_H]`` (A.3). ``P_s`` synodic, ``P_H`` Hohmann for that leg's
      pair.

    The departure V_inf default ``[1, 5.1]`` km/s is Takao's Earth-departure window
    (5.1 km/s = the 1:2 Earth-resonant cap; 1 km/s prevents low-velocity flybys).
    """
    n_legs = len(sequence) - 1
    if n_legs < 1:
        raise ValueError("sequence must name at least two bodies (one leg)")

    t0_lo, t0_hi = t0_window_sec
    vinf_lo, vinf_hi = vinf_out0_bounds_kms

    tof_lo: list[float] = []
    tof_hi: list[float] = []
    # An "inner" pair: both endpoints are at/inside Mars (the E<->M resonant class
    # Takao's A.2 covers). Anything else (Jupiter/Saturn etc.) uses the A.3 window.
    inner_codes = {"Me", "V", "E", "M"}
    for i in range(n_legs):
        a, b = sequence[i], sequence[i + 1]
        # Order the pair by semi-major axis for the synodic/Hohmann helpers.
        inner, outer = (a, b) if PLANETS[a].sma_au <= PLANETS[b].sma_au else (b, a)
        p_h = _hohmann_tof_days(inner, outer)
        if a in inner_codes and b in inner_codes:
            p_s = _synodic_period_days(inner, outer)
            tof_lo.append(30.0)
            tof_hi.append(p_s + p_h)
        else:
            tof_lo.append(0.3 * p_h)
            tof_hi.append(1.3 * p_h)

    lower = np.array(
        [t0_lo, vinf_lo, -np.pi, -0.5 * np.pi, *tof_lo, *([0.0] * n_legs)],
        dtype=np.float64,
    )
    upper = np.array(
        [t0_hi, vinf_hi, np.pi, 0.5 * np.pi, *tof_hi, *([1.0] * n_legs)],
        dtype=np.float64,
    )
    return DsmBounds(lower=lower, upper=upper)


def dsm_chain_decision_vector(
    *,
    t0_sec: float,
    vinf_out0_kms: float,
    alpha0: float,
    beta0: float,
    tof_days_per_leg: tuple[float, ...],
    eta_per_leg: tuple[float, ...],
) -> NDArray[np.float64]:
    """Pack the chain genome into the flat decision vector the corrector consumes.

    Layout ``[t0_sec, vinf_out0, alpha0, beta0, *tof_days_per_leg, *eta_per_leg]``
    (matches :class:`DsmBounds`).
    """
    return np.array(
        [t0_sec, vinf_out0_kms, alpha0, beta0, *tof_days_per_leg, *eta_per_leg],
        dtype=np.float64,
    )


def _unpack(x: NDArray[np.float64], n_legs: int) -> dict[str, Any]:
    """Inverse of :func:`dsm_chain_decision_vector`."""
    return {
        "t0_sec": float(x[0]),
        "vinf_out0_kms": float(x[1]),
        "alpha0": float(x[2]),
        "beta0": float(x[3]),
        "tof_days_per_leg": tuple(float(v) for v in x[4 : 4 + n_legs]),
        "eta_per_leg": tuple(float(v) for v in x[4 + n_legs : 4 + 2 * n_legs]),
    }


def dsm_chain_correct(
    x0: NDArray[np.float64],
    *,
    sequence: tuple[str, ...],
    ephem: Ephemeris,
    bounds: DsmBounds | None = None,
    mu: float = MU_SUN_KM3_S2,
    rendezvous: bool = False,
    tol_kms: float = 0.1,
    max_nfev: int = 200,
    max_revs: int = 0,
) -> DsmChainResult:
    """Drive the chained-DSM total-dV to a minimum with bounded least-squares.

    Free variables = the full decision vector
    ``[t0_sec, vinf_out0, alpha0, beta0, *tof, *eta]`` (see
    :func:`dsm_chain_decision_vector`). The single scalar residual is the chain's
    ``total_dv_kms`` (Takao Eq.15 minus the P-FB term); ``least_squares`` minimises
    it inside ``bounds`` (Takao's box). Converged iff ``total_dv_kms < tol_kms``
    (residual-magnitude only, BY DESIGN -- mirrors ``free_return_correct``).

    This is the local-solve primitive the MBH wrapper drives via
    :func:`make_dsm_chain_step`; MBH supplies the basin hops, this supplies the
    refinement.
    """
    n_legs = len(sequence) - 1

    def _res(x: NDArray[np.float64]) -> NDArray[np.float64]:
        params = _unpack(np.asarray(x, dtype=np.float64), n_legs)
        r = evaluate_dsm_chain(
            sequence=sequence,
            ephem=ephem,
            mu=mu,
            rendezvous=rendezvous,
            tol_kms=tol_kms,
            max_revs=max_revs,
            **params,
        )
        val = r.total_dv_kms
        if not np.isfinite(val):
            val = 1.0e3
        return np.array([val], dtype=np.float64)

    x0_arr = np.asarray(x0, dtype=np.float64)
    if bounds is not None:
        # An MBH hop can perturb the seed just outside the box; ``trf`` rejects an
        # out-of-bounds x0 outright. Clip into the box (strictly interior, to avoid
        # the ``lb >= ub`` / on-edge degeneracies trf also rejects) so the hop's
        # intent is preserved rather than crashing the whole search.
        span = bounds.upper - bounds.lower
        eps = np.where(span > 0.0, 1.0e-9 * span, 0.0)
        x0_arr = np.clip(x0_arr, bounds.lower + eps, bounds.upper - eps)
        sol = least_squares(
            _res,
            x0_arr,
            bounds=(bounds.lower, bounds.upper),
            method="trf",
            max_nfev=max_nfev,
            xtol=1e-12,
            ftol=1e-12,
        )
    else:
        sol = least_squares(_res, x0_arr, method="trf", max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)

    params = _unpack(np.asarray(sol.x, dtype=np.float64), n_legs)
    result = evaluate_dsm_chain(
        sequence=sequence,
        ephem=ephem,
        mu=mu,
        rendezvous=rendezvous,
        tol_kms=tol_kms,
        max_revs=max_revs,
        **params,
    )
    # The evaluator already set ``converged`` from total_dv < tol_kms; carry the
    # solver diagnostics onto the otherwise-complete result (replace keeps every
    # field, including the converged departure-V_inf genome, in sync).
    return replace(result, solver_success=bool(sol.success), solver_nfev=int(sol.nfev))


# ---------------------------------------------------------------------------
# MBH adapter: wrap dsm_chain_correct as an objective_and_solve closure.
# This imports and CALLS the MBH machinery; it never edits mbh.py (the wrapper
# is generic and takes this closure -- see docs/notes/2026-06-07-mbh-wrapper.md).
# ---------------------------------------------------------------------------


def make_dsm_chain_step(
    *,
    sequence: tuple[str, ...],
    ephem: Ephemeris,
    bounds: DsmBounds | None = None,
    mu: float = MU_SUN_KM3_S2,
    rendezvous: bool = False,
    tol_kms: float = 0.1,
    max_nfev: int = 200,
    max_revs: int = 0,
) -> Callable[[np.ndarray, np.random.Generator], MBHStep]:
    """Adapter: :func:`dsm_chain_correct` as an MBH local-solve closure.

    Genome ``x = [t0_sec, vinf_out0, alpha0, beta0, *tof_days, *eta]`` (see
    :func:`dsm_chain_decision_vector`). Objective = ``total_dv_kms`` (Takao Eq.15
    minus the P-FB term); feasible = ``converged``. The emerged per-body V_inf is
    carried in ``info`` as the EVIDENCE the sourced-anchor gate compares against
    (it is never the objective, which would impose it). Deterministic: ``rng``
    accepted for the generic signature and ignored.
    """
    n_legs = len(sequence) - 1

    def step(x: np.ndarray, rng: np.random.Generator) -> MBHStep:
        r = dsm_chain_correct(
            np.asarray(x, dtype=np.float64),
            sequence=sequence,
            ephem=ephem,
            bounds=bounds,
            mu=mu,
            rendezvous=rendezvous,
            tol_kms=tol_kms,
            max_nfev=max_nfev,
            max_revs=max_revs,
        )
        # The corrector moves the WHOLE genome -- including the departure V_inf
        # (vinf_out0/alpha0/beta0), which enters the dV objective through the
        # departure state -- so the landed vector takes the corrector's converged
        # values, not the seed's. This is what lets MBH hop in the departure-energy
        # direction too (the 6.44Gg3 probe needs that lever).
        landed_x = dsm_chain_decision_vector(
            t0_sec=r.t0_sec,
            vinf_out0_kms=r.vinf_out0_kms,
            alpha0=r.alpha0,
            beta0=r.beta0,
            tof_days_per_leg=r.tof_days_per_leg,
            eta_per_leg=r.eta_per_leg,
        )
        return MBHStep(
            x=landed_x,
            objective=float(r.total_dv_kms),
            feasible=bool(r.converged),
            info={
                "total_dv_kms": float(r.total_dv_kms),
                "dv_dsm_per_leg_kms": tuple(r.dv_dsm_per_leg_kms),
                "dv_arrive_kms": float(r.dv_arrive_kms),
                "vinf_in_kms": dict(r.vinf_in_kms),
                "vinf_out_kms": dict(r.vinf_out_kms),
                "eta_per_leg": tuple(r.eta_per_leg),
                "tof_days_per_leg": tuple(r.tof_days_per_leg),
                "n_revs_per_leg": tuple(r.n_revs_per_leg),
                "branch_per_leg": tuple(r.branch_per_leg),
                "solver_nfev": int(r.solver_nfev),
            },
        )

    _ = n_legs  # documented genome arity; kept for adapter symmetry/readability
    return step
