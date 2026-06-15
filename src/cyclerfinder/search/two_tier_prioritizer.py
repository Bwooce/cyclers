"""Two-tier prioritizer: Braik-Ross heading screen + Zhou-Armellin impulse confirm (#263).

Composes the two existing reachable-set scorers into ONE accessibility API for the
discovery program's Track B prioritizer (docs/notes/2026-06-13-discovery-program-spec.md).

The two tiers are NOT redundant; each fills the other's blind spot:

* **Tier 1 -- Braik-Ross energy-PRESERVING reachable-set heading-fan**
  (:mod:`cyclerfinder.search.reachable_network`).
  A pure rotation of the rotating-frame velocity at fixed speed (Braik-Ross Eq. 26,
  ``dV_turn = 2*v*sin(|d|/2)``); the post-maneuver state stays on the SAME ``C_J``
  manifold. The reduced ``(x, y, theta)`` reachable atlas overlaps two families'
  forward/backward sets on a shared voxel grid; an overlap voxel is a heading-only
  proxy-dV bridge between them at one energy. Necessary-not-sufficient screening
  (paper Sec. 3): cheap, blind to energy-changing axes.

* **Tier 2 -- Zhou-Armellin energy-CHANGING single-impulse reachable footprint**
  (:mod:`cyclerfinder.search.reachable_impulsive`).
  An arbitrary-direction impulse at a fixed epoch on the max-magnitude sphere
  (Zhou Eqs. 4-6); the maneuver CHANGES the energy (Jacobi constant). Propagated
  to a horizon and tested for spatial proximity to the target representative's
  orbit. Lights up bridges that cross the energy axis -- exactly what tier 1
  cannot see by construction.

The composition: tier 1 filters the N^2 pair space cheaply (any nonzero overlap
is a candidate); tier 2 quantifies the minimum ΔV bridge for the survivors. Pairs
that tier 1 misses entirely (because the two families sit at different ``C_J``
manifolds and never overlap on the energy-preserving grid) get a tier-2-only
"impulsive only" verdict if their footprints reach with ΔV ≤ budget.

Pure composition: this module does NOT edit either tier's internals. It only
calls the existing public APIs (``build_reachable_set``, ``mirror_reachable_set``,
``pair_proxy``, ``reachable_cloud``, ``apply_impulse``, ``cr3bp.propagate``).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.reachable_impulsive as ri
import cyclerfinder.search.reachable_network as rn

# ---------------------------------------------------------------------------
# Unit conversions: ΔV budgets are user-facing in km/s; both scorers internally
# work in nondimensional VU. We reuse the same EM VU constant the Braik-Ross
# scorer uses (rn.VU_MS = 1023.16 m/s) so tier-1 and tier-2 budgets are
# directly comparable on one scale.
# ---------------------------------------------------------------------------


#: Convert km/s to nondimensional velocity at the Earth-Moon VU.
def kms_to_nd(dv_kms: float) -> float:
    """km/s -> nondimensional VU (Earth-Moon, VU ~ 1.0232 km/s)."""
    return dv_kms * 1000.0 / rn.VU_MS


def nd_to_kms(dv_nd: float) -> float:
    """Nondimensional VU -> km/s (Earth-Moon)."""
    return dv_nd * rn.VU_MS / 1000.0


# ---------------------------------------------------------------------------
# Representative protocol: anything with .state0 (6-vec) and .period (nd float).
# We accept the existing :class:`cyclerfinder.search.reachable_representatives.Representative`
# (which has both), or any duck-typed object exposing those attributes.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RepView:
    """Minimal view of a representative orbit used by the prioritizer.

    Wrappers around any object exposing ``state0`` and ``period`` (e.g. the
    existing :class:`cyclerfinder.search.reachable_representatives.Representative`)
    can be passed directly; this dataclass exists so callers without a
    Representative-typed object (e.g. raw scan output) can still drive the API.
    The ``label`` is a free-form identifier used in cache keys and reports.
    """

    label: str
    state0: NDArray[np.float64]
    period: float


def _as_repview(rep: object, fallback_label: str) -> RepView:
    """Normalize anything with ``state0`` + ``period`` to a :class:`RepView`.

    Accepts the existing
    :class:`cyclerfinder.search.reachable_representatives.Representative` (which
    has a ``label``), bare :class:`RepView`, or any duck-typed object; falls back
    to ``fallback_label`` if no ``.label`` attribute exists.
    """
    if isinstance(rep, RepView):
        return rep
    state0 = np.asarray(getattr(rep, "state0"), dtype=np.float64)  # noqa: B009
    period = float(getattr(rep, "period"))  # noqa: B009
    label = str(getattr(rep, "label", fallback_label))
    return RepView(label=label, state0=state0, period=period)


# ---------------------------------------------------------------------------
# Tier 2 minimum-impulse search.
# ---------------------------------------------------------------------------


def _orbit_position_samples(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    n_samples: int,
) -> NDArray[np.float64]:
    """Sample ``n_samples`` positions along one period of an orbit.

    Used as the target proximity set for tier-2's "did the impulsive cloud reach
    rep_to's orbit?" test. Returns an (n_samples, 3) array of rotating-frame
    positions.
    """
    if n_samples < 2:
        return np.asarray(state0, float)[:3].reshape(1, 3)
    ts = np.linspace(0.0, period, n_samples)
    positions = np.empty((n_samples, 3))
    positions[0] = np.asarray(state0, float)[:3]
    cur = np.asarray(state0, float)
    for i in range(1, n_samples):
        arc = cr3bp.propagate(system, cur, ts[i] - ts[i - 1])
        cur = arc.state_f
        positions[i] = cur[:3]
    return positions


def _impulsive_min_dv_for_pair(
    system: cr3bp.CR3BPSystem,
    rep_from: RepView,
    rep_to: RepView,
    *,
    dv_budget_nd: float,
    n_seeds: int,
    n_alpha: int,
    n_beta: int,
    n_dv_steps: int,
    horizon_frac: float,
    proximity_tol: float,
    n_target_samples: int,
) -> tuple[float, tuple[float, float, float] | None]:
    """Search the minimum single-impulse ΔV that reaches ``rep_to`` from ``rep_from``.

    Sweeps impulse magnitude on a log-spaced ladder ``[dv_budget/N_steps,
    dv_budget]`` (early-exit at the first magnitude whose cloud comes within
    ``proximity_tol`` of any sampled position on ``rep_to``'s orbit). For each
    magnitude, an ``(alpha, beta)`` grid of impulses is launched from ``n_seeds``
    arc-length-distributed seeds on ``rep_from`` and propagated for
    ``horizon_frac * rep_to.period``. The minimum reaching ΔV (nondimensional VU)
    is returned, or ``math.inf`` if no rung reaches.

    Returns
    -------
    (dv_min_nd, hit_info)
        ``dv_min_nd`` is the smallest reaching ΔV (or inf). ``hit_info`` is the
        ``(t_f, alpha, beta)`` of the first hit (or ``None`` if none).
    """
    target_positions = _orbit_position_samples(
        system, rep_to.state0, rep_to.period, n_target_samples
    )
    # Seeds along rep_from: reuse the Braik-Ross arc-length helper for consistency
    # with tier 1 (so both tiers query the same physical seed set).
    seeds = rn._seeds_on_orbit(system, rep_from.state0, rep_from.period, n_seeds)
    t_f = horizon_frac * rep_to.period
    # Log-spaced rung from a small fraction of the budget up to the budget.
    if dv_budget_nd <= 0.0 or n_dv_steps < 1:
        return math.inf, None
    rungs = np.geomspace(max(dv_budget_nd / max(n_dv_steps, 1), 1e-12), dv_budget_nd, n_dv_steps)
    best: float = math.inf
    best_hit: tuple[float, float, float] | None = None
    for dv_nd in rungs:
        for seed in seeds:
            try:
                cloud = ri.reachable_cloud(
                    system, seed, float(dv_nd), t_f, n_alpha=n_alpha, n_beta=n_beta
                )
            except (RuntimeError, ValueError):
                # Degenerate seed (rectilinear / outside Hill / integrator failure);
                # try the next seed rather than abort the search.
                continue
            if cloud.final_states.size == 0:
                continue
            # Min distance from any cloud final position to any target sample.
            dists = np.linalg.norm(
                cloud.final_states[:, None, :3] - target_positions[None, :, :], axis=-1
            )
            min_d = float(dists.min())
            if min_d <= proximity_tol:
                # Find which (alpha, beta) produced the hit for reporting.
                ij = np.unravel_index(int(np.argmin(dists)), dists.shape)
                # Reconstruct alpha/beta from grid index. ri.reachable_cloud uses
                # ``np.linspace(-pi/2, pi/2, n_alpha)`` x ``np.linspace(-pi, pi,
                # n_beta)`` with ``indexing='ij'`` raveled; row index = i*n_beta+j.
                idx = int(ij[0])
                ai = idx // n_beta
                bi = idx % n_beta
                alpha = -math.pi / 2.0 + ai * math.pi / max(n_alpha - 1, 1)
                beta = -math.pi + bi * (2.0 * math.pi) / max(n_beta - 1, 1)
                if float(dv_nd) < best:
                    best = float(dv_nd)
                    best_hit = (t_f, alpha, beta)
                # Found a hit at this rung; no smaller rung on this ladder can
                # do better -- break out of seeds + magnitude sweep.
                return best, best_hit
    return best, best_hit


# ---------------------------------------------------------------------------
# The two-tier prioritizer.
# ---------------------------------------------------------------------------


@dataclass
class TwoTierPrioritizer:
    """Two-tier accessibility prioritizer composing Braik-Ross + Zhou-Armellin.

    Tier 1 (energy-preserving heading-fan screen) is the cheap pre-filter; only
    pairs whose tier-1 overlap clears ``tier1_min_overlap`` enter the tier-2
    impulsive (energy-changing) ΔV-budget confirm. The defaults are tuned for the
    discovery program's Track B usage: every nonzero tier-1 overlap proceeds to
    tier 2 (``tier1_min_overlap = 0.0``, max recall), and the tier-2 budget is
    1 km/s (the Braik-Ross common-energy reference regime, slightly above the
    paper's 409.3 m/s cap to admit cross-energy bridges).

    Parameters
    ----------
    system :
        The CR3BP system both tiers propagate in (must match the ``C_J`` used in
        :meth:`score_pair`).
    grid :
        Tier-1 reduced-state voxel grid (Braik-Ross :class:`VoxelGrid`).
    tier1_min_overlap :
        Tier-1 admission threshold on the *overlap score* (a positive scalar:
        higher = more overlap, ``0.0`` = no overlap). The score is
        ``1.0 / (1.0 + proxy_dv_nd)`` (so a zero-cost overlap voxel gives 1.0
        and ``proxy_dv = inf`` gives 0.0). A pair is "tier-1 admitted" iff
        ``overlap_score > tier1_min_overlap``. Default ``0.0`` means every
        strictly-nonzero overlap proceeds to tier 2.
    tier1_n_seeds, tier1_n_fan, tier1_delta_max, tier1_horizon_frac :
        Tier-1 atlas knobs (passed through to
        :func:`reachable_network.build_reachable_set`).
    tier2_n_seeds, tier2_n_alpha, tier2_n_beta, tier2_n_dv_steps,
    tier2_horizon_frac, tier2_proximity_tol, tier2_n_target_samples :
        Tier-2 minimum-impulse-search knobs (passed through to
        :func:`_impulsive_min_dv_for_pair`).
    """

    system: cr3bp.CR3BPSystem
    grid: rn.VoxelGrid
    tier1_min_overlap: float = 0.0
    tier1_n_seeds: int = 10
    tier1_n_fan: int = 9
    tier1_delta_max: float = math.radians(30.0)
    tier1_horizon_frac: float = 1.0
    tier2_n_seeds: int = 4
    tier2_n_alpha: int = 7
    tier2_n_beta: int = 13
    tier2_n_dv_steps: int = 4
    tier2_horizon_frac: float = 0.5
    tier2_proximity_tol: float = 0.05  # ~19 200 km in LU
    tier2_n_target_samples: int = 24
    _atlas_cache: dict[tuple[str, float], tuple[rn.ReachableSet, rn.ReachableSet]] = field(
        default_factory=dict, repr=False
    )

    # ----- tier 1 -----

    def _atlas(self, rep: RepView, c_j: float) -> tuple[rn.ReachableSet, rn.ReachableSet]:
        """Build (or return cached) tier-1 forward + backward atlas for ``rep`` at ``c_j``."""
        key = (rep.label, float(c_j))
        if key in self._atlas_cache:
            return self._atlas_cache[key]
        horizon = self.tier1_horizon_frac * rep.period
        forward = rn.build_reachable_set(
            self.system,
            rep.state0,
            rep.period,
            self.grid,
            c_j,
            n_seeds=self.tier1_n_seeds,
            n_fan=self.tier1_n_fan,
            delta_max=self.tier1_delta_max,
            horizon=horizon,
        )
        backward = rn.mirror_reachable_set(forward, self.grid)
        self._atlas_cache[key] = (forward, backward)
        return forward, backward

    def _tier1_score(self, rep_from: RepView, rep_to: RepView, c_j: float) -> rn.PairProxy:
        """Tier-1 Braik-Ross proxy overlap (forward A + backward B) on the shared grid."""
        forward_a, _ = self._atlas(rep_from, c_j)
        _, backward_b = self._atlas(rep_to, c_j)
        return rn.pair_proxy(forward_a, backward_b, self.grid)

    # ----- API -----

    def score_pair(
        self,
        rep_from: object,
        rep_to: object,
        *,
        c_j: float,
        dv_budget_kms: float = 1.0,
    ) -> dict[str, object]:
        """Two-tier accessibility score for the ordered pair ``rep_from -> rep_to``.

        Returns a dict with at least:

        ``tier1_heading_overlap``
            The Braik-Ross energy-preserving overlap score (``0.0`` to ``1.0``;
            ``0.0`` = no overlap, larger = cheaper bridge). Defined as
            ``1.0 / (1.0 + proxy_dv_nd)`` so an in-budget overlap is bounded
            away from zero and a missing edge is exactly zero.

        ``tier1_proxy_dv_kms``
            The underlying Braik-Ross proxy ΔV in km/s (``inf`` if no overlap);
            exposed for users who want the raw cost rather than the score.

        ``tier1_admitted``
            ``True`` iff the overlap score exceeds ``tier1_min_overlap``. If
            ``False``, tier 2 is still run (the energy-CHANGING complement may
            light up bridges tier 1 misses by construction).

        ``tier2_impulsive_min_dv_kms``
            Smallest single-impulse ΔV (km/s) that drives the propagated cloud
            within ``tier2_proximity_tol`` of any sampled position on ``rep_to``'s
            orbit, searched on the rung ladder up to ``dv_budget_kms``. ``inf`` if
            no rung reaches.

        ``tier2_within_budget``
            ``True`` iff ``tier2_impulsive_min_dv_kms <= dv_budget_kms``.

        ``accessible``
            ``True`` iff EITHER tier admits the pair (tier 1 overlap below budget,
            OR tier 2 within budget). Tier 1 is the strict "shared C_J manifold"
            bridge; tier 2 is the cross-manifold complement.

        ``dominant_tier``
            ``"tier1"`` if tier 1 overlap (in km/s) is the cheaper binding bridge,
            ``"tier2"`` if the impulsive bridge is, ``"neither"`` if nothing is
            within budget, ``"both"`` if equal.
        """
        from_view = _as_repview(rep_from, "rep_from")
        to_view = _as_repview(rep_to, "rep_to")
        dv_budget_nd = kms_to_nd(dv_budget_kms)
        # ---- tier 1 ----
        # Degenerate: identical reps -> the pair is trivially accessible at zero
        # ΔV (no maneuver needed). Short-circuit before any propagation.
        if from_view.label == to_view.label and np.allclose(from_view.state0, to_view.state0):
            return {
                "rep_from": from_view.label,
                "rep_to": to_view.label,
                "tier1_heading_overlap": 1.0,
                "tier1_proxy_dv_kms": 0.0,
                "tier1_admitted": True,
                "tier2_impulsive_min_dv_kms": 0.0,
                "tier2_within_budget": True,
                "accessible": True,
                "dominant_tier": "both",
                "tier2_hit": None,
            }
        try:
            tier1 = self._tier1_score(from_view, to_view, c_j)
        except (RuntimeError, ValueError):
            tier1 = rn.PairProxy(
                accessible=False, proxy_dv=math.inf, proxy_time=math.inf, voxel=None
            )
        tier1_dv_nd = tier1.proxy_dv if tier1.accessible else math.inf
        tier1_dv_kms = nd_to_kms(tier1_dv_nd) if math.isfinite(tier1_dv_nd) else math.inf
        tier1_overlap_score = 1.0 / (1.0 + tier1_dv_nd) if math.isfinite(tier1_dv_nd) else 0.0
        tier1_admitted = tier1_overlap_score > self.tier1_min_overlap
        # ---- tier 2 (always run; cheap when the budget is small and grid coarse) ----
        if dv_budget_nd > 0.0:
            tier2_nd, hit = _impulsive_min_dv_for_pair(
                self.system,
                from_view,
                to_view,
                dv_budget_nd=dv_budget_nd,
                n_seeds=self.tier2_n_seeds,
                n_alpha=self.tier2_n_alpha,
                n_beta=self.tier2_n_beta,
                n_dv_steps=self.tier2_n_dv_steps,
                horizon_frac=self.tier2_horizon_frac,
                proximity_tol=self.tier2_proximity_tol,
                n_target_samples=self.tier2_n_target_samples,
            )
        else:
            tier2_nd, hit = math.inf, None
        tier2_kms = nd_to_kms(tier2_nd) if math.isfinite(tier2_nd) else math.inf
        tier2_within = tier2_kms <= dv_budget_kms
        # ---- compose ----
        tier1_within = tier1_admitted and tier1_dv_kms <= dv_budget_kms
        accessible = tier1_within or tier2_within
        dominant = self._dominant_tier(tier1_dv_kms, tier2_kms, dv_budget_kms)
        return {
            "rep_from": from_view.label,
            "rep_to": to_view.label,
            "tier1_heading_overlap": tier1_overlap_score,
            "tier1_proxy_dv_kms": tier1_dv_kms,
            "tier1_admitted": tier1_admitted,
            "tier2_impulsive_min_dv_kms": tier2_kms,
            "tier2_within_budget": tier2_within,
            "accessible": accessible,
            "dominant_tier": dominant,
            "tier2_hit": hit,
        }

    @staticmethod
    def _dominant_tier(t1_kms: float, t2_kms: float, budget_kms: float) -> str:
        """Which tier is the cheaper binding bridge (within ``budget_kms``)?"""
        t1_ok = math.isfinite(t1_kms) and t1_kms <= budget_kms
        t2_ok = math.isfinite(t2_kms) and t2_kms <= budget_kms
        if not t1_ok and not t2_ok:
            return "neither"
        if t1_ok and not t2_ok:
            return "tier1"
        if t2_ok and not t1_ok:
            return "tier2"
        # Both within budget -- the cheaper bridge wins, ties broken to "both".
        if math.isclose(t1_kms, t2_kms, rel_tol=1e-9, abs_tol=1e-9):
            return "both"
        return "tier1" if t1_kms < t2_kms else "tier2"

    def rank_destinations(
        self,
        rep_from: object,
        candidates: Sequence[object] | Iterable[object],
        *,
        c_j: float,
        dv_budget_kms: float = 1.0,
    ) -> list[dict[str, object]]:
        """Score every ``rep_from -> candidate`` pair and return them sorted.

        The sort key is the minimum of ``(tier1_proxy_dv_kms,
        tier2_impulsive_min_dv_kms)`` (the cheapest binding bridge in km/s), so
        the best-accessible destinations come first; inaccessible pairs (both
        tiers ``inf``) sort to the end.
        """
        scored: list[dict[str, object]] = []
        for cand in candidates:
            scored.append(self.score_pair(rep_from, cand, c_j=c_j, dv_budget_kms=dv_budget_kms))

        def _key(d: dict[str, object]) -> float:
            t1 = float(d["tier1_proxy_dv_kms"])  # type: ignore[arg-type]
            t2 = float(d["tier2_impulsive_min_dv_kms"])  # type: ignore[arg-type]
            return min(t1, t2)

        return sorted(scored, key=_key)


__all__ = [
    "RepView",
    "TwoTierPrioritizer",
    "kms_to_nd",
    "nd_to_kms",
]
