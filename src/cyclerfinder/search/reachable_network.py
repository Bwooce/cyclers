"""Reduced (x, y, theta) heading-fan reachable-set accessibility-network scorer.

Implements the Braik & Ross 2026 ("Orbital Networks in the Three-Body Problem",
arXiv:2605.31543) reachable-set family-accessibility method for the *planar*
Earth-Moon CR3BP, framed as a family-selection / continuation prioritizer (NOT a
transfer designer, NOT a new-family generator). See the mining note
``docs/notes/2026-06-13-braik-ross-2026-orbital-networks-mining.md`` for the full
spec.

Method (Secs. 2-5 of the paper):

1. At one fixed Jacobi constant ``C_J``, reduce the planar state to ``(x, y, t)``
   where ``t`` (theta) is the rotating-frame velocity heading and the speed
   ``v = sqrt(-2*Ubar - C_J)`` is fixed by position (Eq. 8). Time-reversal
   symmetry ``R(x, y, t) = (x, -y, pi - t)`` (Eq. 14).
2. From arc-length-spaced seeds on each representative orbit, apply an
   energy-preserving instantaneous heading-change maneuver (a pure rotation of
   the rotating-frame velocity at fixed speed -> stays on the same ``C_J``
   manifold; ``dV_turn = 2*v*sin(|d|/2)``, Eq. 26), ``|d| <= delta_max``;
   propagate one-sided to horizon ``T_a`` and voxel-log every crossed
   ``(x, y, t)`` cell. The backward reachable set is the time-reversal mirror of
   the forward set (free).
3. A -> B is directly accessible iff ``forward(A) & backward(B) != {}`` on the
   shared voxel grid. Each overlap voxel gets a proxy ``dV`` (min source-side +
   min target-side turning cost + a voxel-scale heading-mismatch patch term) and
   a proxy time; the pair value is the min-proxy-dV voxel. This is explicitly
   NECESSARY-NOT-SUFFICIENT screening.
4. Assemble the N x N symmetric proxy-dV matrix into a weighted undirected graph;
   compute cost-aware centralities: strength (hub), harmonic closeness (gateway),
   betweenness (relay).

Limitations carried throughout: planar only / single ``C_J`` / heading-only
maneuver / screening-not-transfer. Proxy costs are a conservative upper bound
(the paper's Sec. 7 validation: corrected dV < proxy dV in every tested case),
not optimized transfers.

Pure: math / numpy / scipy + ``cyclerfinder.core.cr3bp``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp

# ---------------------------------------------------------------------------
# Reduced (x, y, theta) model: speed relation (Eq. 8), turn cost (Eq. 26),
# time-reversal symmetry (Eq. 14).
# ---------------------------------------------------------------------------


def ubar(x: float, y: float, mu: float) -> float:
    """Planar pseudo-potential Ubar = -1/2 (x^2 + y^2) - (1-mu)/r1 - mu/r2.

    Braik-Ross Eq. 3 (planar; z = 0). With this sign convention the Jacobi
    constant is ``C_J = -2*Ubar - v^2`` (Eq. 5), matching our
    :func:`cyclerfinder.core.cr3bp.jacobi_constant`.
    """
    r1 = math.sqrt((x + mu) ** 2 + y * y)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y)
    return -0.5 * (x * x + y * y) - (1.0 - mu) / r1 - mu / r2


def reduced_speed(x: float, y: float, mu: float, c_j: float) -> float:
    """Rotating-frame speed at position ``(x, y)`` on the ``C_J`` manifold (Eq. 8).

    ``v = sqrt(-2*Ubar(x, y) - C_J)``. The radicand is non-negative only inside
    the energetically admissible (Hill) region; outside it the position is not
    reachable at this energy.

    Raises
    ------
    ValueError
        If the radicand ``-2*Ubar - C_J`` is negative (position outside the Hill
        region at this ``C_J``: the zero-velocity curve has been crossed).
    """
    rad = -2.0 * ubar(x, y, mu) - c_j
    if rad < 0.0:
        raise ValueError(
            f"reduced_speed: position ({x:.6f}, {y:.6f}) is forbidden at C_J={c_j:.6f} "
            f"(radicand {rad:.3e} < 0; outside the Hill region)"
        )
    return math.sqrt(rad)


def is_admissible(x: float, y: float, mu: float, c_j: float) -> bool:
    """True iff ``(x, y)`` is energetically admissible (inside the Hill region)."""
    return (-2.0 * ubar(x, y, mu) - c_j) >= 0.0


def dv_turn(v: float, delta: float) -> float:
    """Energy-preserving heading-change maneuver cost (Braik-Ross Eq. 26).

    A pure rotation of the rotating-frame velocity by angle ``delta`` at fixed
    speed ``v`` costs ``dV_turn = 2*v*sin(|delta|/2)`` (the chord of the turn on
    the constant-speed circle). Because speed is unchanged the post-maneuver
    state stays on the same ``C_J`` manifold.
    """
    return 2.0 * v * abs(math.sin(0.5 * delta))


def heading(vx: float, vy: float) -> float:
    """Velocity heading theta = atan2(vy, vx), in (-pi, pi]."""
    return math.atan2(vy, vx)


def velocity_from_heading(v: float, theta: float) -> tuple[float, float]:
    """Inverse of :func:`heading`: (vx, vy) = (v cos theta, v sin theta)."""
    return v * math.cos(theta), v * math.sin(theta)


def wrap_angle(theta: float) -> float:
    """Wrap an angle to (-pi, pi]."""
    a = math.fmod(theta + math.pi, 2.0 * math.pi)
    if a <= 0.0:
        a += 2.0 * math.pi
    return a - math.pi


def angular_diff(a: float, b: float) -> float:
    """Signed smallest difference a - b wrapped to (-pi, pi]."""
    return wrap_angle(a - b)


def time_reversal(x: float, y: float, theta: float) -> tuple[float, float, float]:
    """Time-reversal symmetry R(x, y, theta) = (x, -y, pi - theta) (Eq. 14).

    The CR3BP is invariant under ``(t, x, y, vx, vy) -> (-t, x, -y, vx, -vy)``.
    In reduced coordinates a velocity ``(vx, vy)`` at heading ``theta`` maps to
    ``(vx, -vy)`` at heading ``-theta``; the reversed *direction of travel*
    (forward arc run backward in time) is the opposite heading, i.e.
    ``pi - theta``. The returned theta is wrapped to (-pi, pi].
    """
    return x, -y, wrap_angle(math.pi - theta)


# ---------------------------------------------------------------------------
# Voxel grid (x, y, theta).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VoxelGrid:
    """Uniform (x, y, theta) voxel grid over the admissible planar domain.

    ``dx``, ``dy`` in length units (Braik-Ross use 0.001 LU ~ 384 km);
    ``dtheta`` in radians (paper uses 1 deg). The grid origin is implicit at
    (0, 0, -pi); a cell index is ``(floor((x-x0)/dx), floor((y-y0)/dy),
    floor((theta+pi)/dtheta))`` with theta wrapped to (-pi, pi]. The theta axis
    wraps modulo ``n_theta``.
    """

    dx: float
    dy: float
    dtheta: float
    x0: float = -2.0
    y0: float = -2.0

    @property
    def n_theta(self) -> int:
        return round(2.0 * math.pi / self.dtheta)

    def index(self, x: float, y: float, theta: float) -> tuple[int, int, int]:
        """Voxel index containing ``(x, y, theta)`` (theta wrapped, axis modular)."""
        ix = math.floor((x - self.x0) / self.dx)
        iy = math.floor((y - self.y0) / self.dy)
        it = math.floor((wrap_angle(theta) + math.pi) / self.dtheta) % self.n_theta
        return ix, iy, it

    def center(self, idx: tuple[int, int, int]) -> tuple[float, float, float]:
        """Geometric center of voxel ``idx`` (theta in (-pi, pi])."""
        ix, iy, it = idx
        cx = self.x0 + (ix + 0.5) * self.dx
        cy = self.y0 + (iy + 0.5) * self.dy
        ct = wrap_angle(-math.pi + (it + 0.5) * self.dtheta)
        return cx, cy, ct


# ---------------------------------------------------------------------------
# Reachable-set atlas for a single representative orbit.
# ---------------------------------------------------------------------------


@dataclass
class ReachableSet:
    """Forward reachable atlas for one representative orbit on a shared grid.

    ``voxels`` maps voxel index -> minimum source-side turning cost ``dV_turn``
    used to first reach that voxel (the cheapest maneuver across all seeds/fans
    whose propagated arc crossed it). ``times`` maps voxel index -> the minimum
    proxy time (seconds-free; nondimensional propagation time) at which the voxel
    was first crossed by the cheapest-cost arc reaching it.
    """

    voxels: dict[tuple[int, int, int], float] = field(default_factory=dict)
    times: dict[tuple[int, int, int], float] = field(default_factory=dict)

    def _log(self, idx: tuple[int, int, int], cost: float, t: float) -> None:
        prev = self.voxels.get(idx)
        if prev is None or cost < prev:
            self.voxels[idx] = cost
            self.times[idx] = t
        elif cost == prev and t < self.times[idx]:
            self.times[idx] = t


def _seeds_on_orbit(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    n_seeds: int,
) -> list[NDArray[np.float64]]:
    """Sample ``n_seeds`` states evenly in *arc length* along one period.

    Arc length is approximated by integrating speed over a fine time sampling of
    the orbit and inverting the cumulative-arc function; this spaces seeds by
    geometry rather than by time (so fast and slow portions of the orbit get
    proportionate coverage), matching the paper's "arc-length-spaced seeds".
    """
    n_fine = max(200, 20 * n_seeds)
    ts = np.linspace(0.0, period, n_fine)
    states = np.empty((n_fine, 4))
    states[0] = np.asarray(state0, float)[[0, 1, 3, 4]]
    cur = np.asarray(state0, float)
    for i in range(1, n_fine):
        arc = cr3bp.propagate(system, cur, ts[i] - ts[i - 1])
        cur = arc.state_f
        states[i] = cur[[0, 1, 3, 4]]
    speeds = np.hypot(states[:, 2], states[:, 3])
    seg = 0.5 * (speeds[1:] + speeds[:-1]) * np.diff(ts)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(cum[-1])
    if total <= 0.0:
        # Degenerate (should not happen for a real orbit); fall back to time.
        targets = np.linspace(0.0, period, n_seeds, endpoint=False)
        axis = ts
    else:
        targets = np.linspace(0.0, total, n_seeds, endpoint=False)
        axis = cum
    seed_ts = np.interp(targets, axis, ts)
    seeds: list[NDArray[np.float64]] = []
    cur = np.asarray(state0, float)
    t_prev = 0.0
    for tt in seed_ts:
        arc = cr3bp.propagate(system, cur, float(tt) - t_prev)
        cur = arc.state_f
        t_prev = float(tt)
        seeds.append(cur.copy())
    return seeds


def build_reachable_set(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    grid: VoxelGrid,
    c_j: float,
    *,
    n_seeds: int = 12,
    n_fan: int = 9,
    delta_max: float = math.radians(30.0),
    horizon: float | None = None,
    n_log: int = 60,
) -> ReachableSet:
    """Forward heading-fan reachable atlas for one representative orbit.

    From ``n_seeds`` arc-length-spaced seeds, fan the rotating-frame velocity by
    ``n_fan`` headings in ``[-delta_max, +delta_max]`` (Eq. 26 turn cost),
    propagate each one-sided to ``horizon`` (default one period), and voxel-log
    every crossed cell with the source-side turn cost and crossing time. The
    backward set of any orbit is the time-reversal mirror of its forward set
    (:func:`mirror_reachable_set`), so only the forward set is propagated.
    """
    rs = ReachableSet()
    horizon = float(period) if horizon is None else float(horizon)
    deltas = np.linspace(-delta_max, delta_max, n_fan)
    seeds = _seeds_on_orbit(system, state0, period, n_seeds)
    log_ts = np.linspace(horizon / n_log, horizon, n_log)
    for seed in seeds:
        x0, y0 = float(seed[0]), float(seed[1])
        try:
            v = reduced_speed(x0, y0, system.mu, c_j)
        except ValueError:
            continue
        theta0 = heading(float(seed[3]), float(seed[4]))
        for delta in deltas:
            cost = dv_turn(v, float(delta))
            vx, vy = velocity_from_heading(v, theta0 + float(delta))
            man_state = np.array([x0, y0, 0.0, vx, vy, 0.0])
            cur = man_state
            t_prev = 0.0
            ok = True
            for tt in log_ts:
                try:
                    arc = cr3bp.propagate(system, cur, float(tt) - t_prev)
                except RuntimeError:
                    ok = False
                    break
                cur = arc.state_f
                t_prev = float(tt)
                th = heading(float(cur[3]), float(cur[4]))
                idx = grid.index(float(cur[0]), float(cur[1]), th)
                rs._log(idx, cost, float(tt))
            if not ok:
                continue
    return rs


def mirror_reachable_set(forward: ReachableSet, grid: VoxelGrid) -> ReachableSet:
    """Backward reachable atlas = time-reversal mirror of the forward atlas.

    Each forward voxel ``(ix, iy, it)`` with center ``(x, y, theta)`` maps to its
    time-reversal image ``R(x, y, theta) = (x, -y, pi - theta)`` (Eq. 14); the
    cost and time are preserved (a forward arc reaching a voxel becomes, run
    backward in time, an arc *arriving from* the mirror voxel at the same turn
    cost and elapsed time). This makes the backward set free of extra
    propagation.
    """
    back = ReachableSet()
    for idx, cost in forward.voxels.items():
        cx, cy, ct = grid.center(idx)
        mx, my, mt = time_reversal(cx, cy, ct)
        midx = grid.index(mx, my, mt)
        back._log(midx, cost, forward.times[idx])
    return back


# ---------------------------------------------------------------------------
# Pairwise overlap -> proxy (dV, time).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PairProxy:
    """Min-proxy-dV overlap value for an ordered pair (source A -> target B)."""

    accessible: bool
    proxy_dv: float
    proxy_time: float
    voxel: tuple[int, int, int] | None


def pair_proxy(
    forward_a: ReachableSet,
    backward_b: ReachableSet,
    grid: VoxelGrid,
) -> PairProxy:
    """Screening proxy for A -> B: min over shared voxels of the proxy cost.

    A -> B is directly accessible iff ``forward(A) & backward(B) != {}``. For each
    shared voxel the proxy ``dV`` is the source-side turn cost (to reach the
    voxel from A) plus the target-side turn cost (to arrive at B, from the mirror
    set) plus a *voxel-scale heading-mismatch patch term*: the cost of the small
    residual turn needed to reconcile the two arcs' headings within the voxel,
    bounded by half a voxel of heading, ``dV_turn(v, dtheta/2)``. The pair value
    is the minimum-proxy-dV voxel; its proxy time is the sum of the two arcs'
    crossing times. NECESSARY-NOT-SUFFICIENT screening.
    """
    shared = forward_a.voxels.keys() & backward_b.voxels.keys()
    if not shared:
        return PairProxy(accessible=False, proxy_dv=math.inf, proxy_time=math.inf, voxel=None)
    best_dv = math.inf
    best_t = math.inf
    best_idx: tuple[int, int, int] | None = None
    for idx in shared:
        src = forward_a.voxels[idx]
        tgt = backward_b.voxels[idx]
        # Voxel-scale heading-mismatch patch: a half-voxel residual turn. The
        # local speed is reconstructed from the voxel center on the C_J manifold;
        # it is independent of the pair, so it only breaks ties geometrically.
        patch = _patch_cost(idx, grid)
        dv = src + tgt + patch
        t = forward_a.times[idx] + backward_b.times[idx]
        if dv < best_dv or (dv == best_dv and t < best_t):
            best_dv = dv
            best_t = t
            best_idx = idx
    return PairProxy(accessible=True, proxy_dv=best_dv, proxy_time=best_t, voxel=best_idx)


def _patch_cost(idx: tuple[int, int, int], grid: VoxelGrid) -> float:
    """Voxel-scale heading-mismatch patch term (a half-voxel residual turn).

    Uses a unit reference speed so the term is purely a geometric tie-breaker of
    order ``dtheta`` and never dominates the physical source/target turn costs;
    it encodes that two arcs landing in the same voxel may differ by up to one
    voxel of heading and need a small patch maneuver to reconcile.
    """
    return dv_turn(1.0, 0.5 * grid.dtheta)


def proxy_matrix(
    forward: list[ReachableSet],
    backward: list[ReachableSet],
    grid: VoxelGrid,
) -> NDArray[np.float64]:
    """Symmetric N x N proxy-dV matrix from per-family forward/backward atlases.

    ``M[i, j] = min(proxy(i -> j), proxy(j -> i))`` so the result is symmetric
    (the paper's undirected graph). Missing edges (no overlap either direction)
    are ``inf``; the diagonal is 0.
    """
    n = len(forward)
    mat = np.full((n, n), math.inf)
    np.fill_diagonal(mat, 0.0)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            p = pair_proxy(forward[i], backward[j], grid)
            if p.accessible and p.proxy_dv < mat[i, j]:
                mat[i, j] = p.proxy_dv
    # Symmetrize: undirected edge weight = the cheaper direction.
    sym: NDArray[np.float64] = np.minimum(mat, mat.T)
    return sym


# ---------------------------------------------------------------------------
# Cost-aware graph centralities (strength / harmonic closeness / betweenness).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Centralities:
    """Per-node cost-aware centralities on the proxy-dV graph."""

    strength: NDArray[np.float64]  # sum of reciprocal edge costs (hub)
    harmonic_closeness: NDArray[np.float64]  # sum of 1/shortest-path-cost (gateway)
    betweenness: NDArray[np.float64]  # fraction of shortest paths through node (relay)


def _floyd_warshall(weights: NDArray[np.float64]) -> NDArray[np.float64]:
    """All-pairs shortest-path cost on a symmetric weighted graph (inf = no edge)."""
    n = weights.shape[0]
    dist: NDArray[np.float64] = weights.copy()
    np.fill_diagonal(dist, 0.0)
    for k in range(n):
        dik = dist[:, k][:, None]
        dkj = dist[k, :][None, :]
        dist = np.minimum(dist, dik + dkj)
    return dist


def _betweenness(weights: NDArray[np.float64]) -> NDArray[np.float64]:
    """Brandes-style betweenness on a symmetric weighted graph.

    Counts, for each node, the fraction of all shortest source-target paths
    (over unordered pairs, both endpoints excluded) that pass through it,
    handling ties by even split. Edge weight = proxy dV (lower = preferred).
    """
    n = weights.shape[0]
    bet = np.zeros(n)
    dist = _floyd_warshall(weights)
    # Count shortest paths via DP over nodes sorted by distance from each source.
    for s in range(n):
        order = np.argsort(dist[s])
        sigma = np.zeros(n)
        sigma[s] = 1.0
        preds: list[list[int]] = [[] for _ in range(n)]
        finite = [int(t) for t in order if math.isfinite(dist[s, t])]
        for t in finite:
            for u in range(n):
                if u == t or not math.isfinite(weights[u, t]):
                    continue
                via = dist[s, u] + weights[u, t]
                if math.isclose(via, dist[s, t], rel_tol=1e-9, abs_tol=1e-12):
                    sigma[t] += sigma[u]
                    preds[t].append(u)
        delta = np.zeros(n)
        for t in reversed(finite):
            if t == s:
                continue
            for u in preds[t]:
                if sigma[t] > 0:
                    delta[u] += (sigma[u] / sigma[t]) * (1.0 + delta[t])
            if t != s:
                bet[t] += delta[t]
    # Undirected graph: each unordered pair counted twice.
    return bet / 2.0


def centralities(weights: NDArray[np.float64]) -> Centralities:
    """Cost-aware strength / harmonic closeness / betweenness (Braik-Ross Sec. 5).

    ``weights`` is the symmetric proxy-dV matrix (inf = no edge, 0 diagonal).
    - strength: sum over neighbours of ``1/weight`` (reciprocal-cost hub score).
    - harmonic closeness: sum over other nodes of ``1/shortest-path-cost``
      (gateway/staging score; robust to disconnected components).
    - betweenness: fraction of all-pairs shortest paths through the node (relay).
    """
    n = weights.shape[0]
    strength = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j and math.isfinite(weights[i, j]) and weights[i, j] > 0.0:
                strength[i] += 1.0 / weights[i, j]
    dist = _floyd_warshall(weights)
    harmonic = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j and math.isfinite(dist[i, j]) and dist[i, j] > 0.0:
                harmonic[i] += 1.0 / dist[i, j]
    bet = _betweenness(weights)
    return Centralities(strength=strength, harmonic_closeness=harmonic, betweenness=bet)


# ---------------------------------------------------------------------------
# Budget-capped network (Braik-Ross Sec. 5.1, Eqs. 53-62).
# ---------------------------------------------------------------------------

#: Earth-Moon velocity unit VU = LU / TU = 384400 km / (4.34837740 d) -> m/s.
#: Converts a nondimensional proxy turn-cost (Eq. 26) to m/s so the budget cap
#: (Braik-Ross max-budget reference dV_cap = 409.3 m/s) can be applied.
VU_MS = 384400.0 / (4.34837740 * 86400.0) * 1000.0

#: Braik-Ross maximum-budget reference maneuver cap (Sec. 5.1): dV_cap = 409.3 m/s.
#: At this cap the 13-node network retains 75 of 78 edges; the three dropped edges
#: all involve the 2:1 stable resonant R21-S (the persistent hard-access family).
DV_CAP_MS = 409.3


def apply_budget_cap(
    weights_nd: NDArray[np.float64],
    *,
    dv_cap_ms: float = DV_CAP_MS,
) -> NDArray[np.float64]:
    """Convert a nondimensional proxy-dV matrix to m/s and drop over-budget edges.

    Braik-Ross edge-retention rule (Eq. 54): an undirected edge ``(A, B)`` is kept
    iff its direct proxy cost ``dV_A,B <= dV_cap``. Edges above the cap (and the
    already-missing ``inf`` edges) are set to ``inf`` (no edge). The returned
    matrix is in m/s with a zero diagonal -- ready for :func:`centralities`.
    """
    mat = weights_nd * VU_MS
    np.fill_diagonal(mat, 0.0)
    over = mat > dv_cap_ms
    mat[over] = math.inf
    np.fill_diagonal(mat, 0.0)
    return mat


def normalized_centralities(weights_ms: NDArray[np.float64], n_families: int) -> Centralities:
    """Braik-Ross normalized centralities (Eqs. 59-62) on a capped m/s matrix.

    Applies the paper's per-metric normalizations to :func:`centralities`:
      - strength  S(A) = (1/(Nf-1)) * sum 1/dV  (Eq. 60);
      - harmonic  H(A) = (1/(Nf-1)) * sum 1/d   (Eq. 61);
      - betweenness B(A) = (2/((Nf-1)(Nf-2))) * sum sigma_PQ(A)/sigma_PQ (Eq. 62).
    ``n_families`` is ``Nf`` (the full representative-set size the paper normalizes
    against). The normalizations are constant per metric, so node *ranks* are
    identical to :func:`centralities`; they are applied so the reported values are
    directly comparable to Table 4.
    """
    raw = centralities(weights_ms)
    nf = float(n_families)
    s_norm = 1.0 / (nf - 1.0)
    b_norm = 2.0 / ((nf - 1.0) * (nf - 2.0))
    return Centralities(
        strength=raw.strength * s_norm,
        harmonic_closeness=raw.harmonic_closeness * s_norm,
        betweenness=raw.betweenness * b_norm,
    )
