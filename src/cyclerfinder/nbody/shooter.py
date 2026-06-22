"""Multiple-shooting differential corrector in restricted n-body (design §3, #133).

Consumer 2 of the harness: the **Jones flyby-propagation shooter**. Multiple
shooting over encounter nodes — full-state continuity in real (rails) dynamics
between nodes propagated by :class:`cyclerfinder.nbody.propagator.RestrictedNBody`
— that corrects a patched-conic seed to an n-body-ballistic cycler. This is the
SNOPT analogue of Jones, Hernandez & Jesick (AAS 17-577, 2017): their stage 1
accepted conic v∞-mismatches <= 200 m/s, then ran SNOPT n-body correction to
ballistic; here the conic seed is the START and the multiple-shooting solve drives
the full-state defects to zero.

**SEEDING — the #135 verdict (binding).** The #135 like-for-like diagnostic
(``docs/notes/2026-06-06-russell12-likeforlike.md``) found 0 CLOSE-AND-MATCH
coplanar-vs-coplanar on known-solvable instances: the corrector closes
geometrically but lands OFF-ANCHOR (our V∞ 9-28 vs sourced 3-10 km/s). The
verdict is **seeding/basin, not solver deficiency**. The implication, carried into
this shooter, is absolute: a single-shoot from a naive equispaced / coplanar /
blind-scan seed falls into the same high-V∞ basin. **This shooter MUST be seeded
from the #133 near-miss conic survey** — the lowest-residual low-V∞ conic chains
near the Jones anchors collected by :func:`near_miss_survey` — NEVER from a blind
scan. The near-miss survey (Phase 3a) is the seeding source; the multiple-shooting
solve refines from there.

**GOLDEN DISCIPLINE.** The shooter OUTPUT (n-body-ballistic V∞ at the converged
nodes) is the side under test; the EXPECTED side is the SOURCED Jones multiset
(AAS 17-577 Tables 2/3) only. No self-computed value is ever EXPECTED. Divergence
is a first-class non-converged record (mirror ``correct.py``'s honest
non-converged result), never an exception.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS, SECONDS_PER_DAY

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris

Vec3 = NDArray[np.float64]

# Full Cartesian state dimension per node (r[3] + v[3]).
_STATE_DIM = 6


def defect_count(*, n_encounters: int) -> int:
    """Number of scalar full-state continuity defects for ``n_encounters`` nodes.

    A cycler with ``n`` encounter nodes has ``n - 1`` interior legs, each
    contributing a full 6-component state defect (the propagated arc end vs the
    next node state). E.g. E-M-E-V-V-E (6 encounters) -> 5 legs -> 30 defects.
    """
    if n_encounters < 2:
        raise ValueError("a cycler needs at least 2 encounter nodes")
    return (n_encounters - 1) * _STATE_DIM


def build_shooting_vector(
    nodes: Mapping[str, Vec3],
    epochs: Sequence[float],
    tofs: Sequence[float],
    *,
    slack_leg: int,
    period_days: float,
) -> NDArray[np.float64]:
    """Pack the multiple-shooting free-variable vector ``x``.

    Layout: ``[ {node states, 6 each, in key order}, {node epochs}, {free ToFs} ]``.
    The **slack leg** ToF is eliminated from the free vector (it is reconstructed
    from the period pin ``period_days - sum(other ToFs)``, the
    ``correct.py:_reconstruct_tofs`` convention, ``correct.py:77-86``) so the
    period is held exactly. Node states are stacked in ``sorted`` key order for a
    stable, reproducible packing (the same ``b{i}`` vocabulary as
    ``correct._vinf_nodes``).
    """
    if not 0 <= slack_leg < len(tofs):
        raise ValueError(f"slack_leg {slack_leg} out of range for {len(tofs)} ToFs")
    parts: list[NDArray[np.float64]] = []
    for key in sorted(nodes):
        parts.append(np.asarray(nodes[key], dtype=np.float64).ravel())
    parts.append(np.asarray(epochs, dtype=np.float64).ravel())
    free_tofs = [t for i, t in enumerate(tofs) if i != slack_leg]
    parts.append(np.asarray(free_tofs, dtype=np.float64))
    return np.concatenate(parts)


@dataclass(frozen=True)
class ShootingSeed:
    """A multiple-shooting seed: per-node full states + epochs + leg ToFs.

    ``node_states`` is one ``(6,)`` heliocentric J2000-ecliptic state ``[r|v]``
    (km, km/s) per encounter, in encounter order ``b0..b{n-1}``. ``epochs`` are the
    encounter epochs on the **TDB-J2000-seconds axis** (§0). ``tofs`` are the
    inter-node flight times (days). ``sequence`` is the body code at each node.
    ``slack_leg`` / ``period_days`` carry the period pin (``correct.py``
    convention). The per-node V∞ vectors (``vinf_out`` for the outgoing leg of each
    node) are retained so the flyby-bend hinge and the correction-ΔV accounting can
    be computed without re-solving Lambert.
    """

    node_states: list[Vec3]
    epochs: list[float]
    tofs: list[float]
    sequence: tuple[str, ...]
    slack_leg: int
    period_days: float
    vinf_in: list[Vec3]
    vinf_out: list[Vec3]


def seed_from_conic(
    *,
    sequence: tuple[str, ...],
    vinf_nodes: Mapping[str, Vec3],
    t0_sec: float,
    tofs_days: Sequence[float],
    slack_leg: int,
    period_days: float,
    ephem: Ephemeris,
) -> ShootingSeed:
    """Map a patched-conic chain to a multiple-shooting seed (design §3 + drift).

    Reads node V∞ **vectors** from a ``correct._vinf_nodes``-shaped mapping
    (``b{i}_in`` / ``b{i}_out`` keys) — NOT from a ``best_cycler.encounters``
    attribute, which the live ``BallisticClosureResult`` does not carry (the design
    §3 DRIFT flagged in the plan: the result holds only scalar
    ``vinf_per_encounter_kms``; the vectors live in ``_vinf_nodes``). Each node's
    spacecraft full state is the encounter-body state plus the **outgoing** V∞ at
    that node (``v_sc = v_planet + vinf_out``, the ``verify/propagate.py``
    reconstruction): the state that *begins* the leg the shooter propagates. The
    wrap node ``b{n-1}`` (home arrival) uses its inbound V∞.
    """
    n = len(sequence)
    epochs = [t0_sec]
    for tof in tofs_days:
        epochs.append(epochs[-1] + float(tof) * SECONDS_PER_DAY)
    node_states: list[Vec3] = []
    vinf_in: list[Vec3] = []
    vinf_out: list[Vec3] = []
    for i, body in enumerate(sequence):
        r_pl, v_pl = ephem.state(body, epochs[i])
        r_pl = np.asarray(r_pl, dtype=np.float64)
        v_pl = np.asarray(v_pl, dtype=np.float64)
        v_out = np.asarray(vinf_nodes.get(f"b{i}_out", np.zeros(3)), dtype=np.float64)
        v_inn = np.asarray(vinf_nodes.get(f"b{i}_in", np.zeros(3)), dtype=np.float64)
        # The node state begins the leg the shooter propagates -> outgoing V∞.
        # The wrap node has no outgoing leg; use its inbound V∞ for the state.
        vinf_state = v_inn if (i == n - 1) else v_out
        v_sc = v_pl + vinf_state
        node_states.append(np.concatenate([r_pl, v_sc]))
        vinf_in.append(v_inn)
        vinf_out.append(v_out)
    return ShootingSeed(
        node_states=node_states,
        epochs=epochs,
        tofs=[float(t) for t in tofs_days],
        sequence=tuple(sequence),
        slack_leg=slack_leg,
        period_days=period_days,
        vinf_in=vinf_in,
        vinf_out=vinf_out,
    )


def _flyby_periapsis_hinge_km(v_in: Vec3, v_out: Vec3, body: str) -> float:
    """Bend-feasibility hinge: shortfall of the achievable periapsis below safe.

    An unpowered flyby turns the V∞ vector by ``delta = angle(v_in, v_out)`` at
    constant magnitude; the periapsis radius needed for that turn is
    ``r_p = (mu / v∞^2) (1/sin(delta/2) - 1)``. If ``r_p`` is below the safe flyby
    radius (``radius_eq_km + safe_alt_km``) the turn is infeasible (the spacecraft
    would have to pass inside the body). The hinge returns ``max(0, r_safe - r_p)``
    in km — zero when feasible, positive (a residual the solver drives down) when
    the required bend is too sharp for the body to deliver ballistically.
    """
    p = PLANETS[body]
    r_safe = p.radius_eq_km + p.safe_alt_km
    vin = float(np.linalg.norm(v_in))
    vout = float(np.linalg.norm(v_out))
    vmag = 0.5 * (vin + vout)
    if vmag <= 0.0:
        return 0.0
    cos_d = float(np.dot(v_in, v_out) / (vin * vout)) if vin > 0 and vout > 0 else 1.0
    cos_d = max(-1.0, min(1.0, cos_d))
    delta = float(np.arccos(cos_d))
    half = 0.5 * delta
    sin_half = float(np.sin(half))
    if sin_half <= 1e-12:
        return 0.0  # negligible bend -> any periapsis works
    r_p = (p.mu_km3_s2 / (vmag * vmag)) * (1.0 / sin_half - 1.0)
    return max(0.0, r_safe - r_p)


def defect_residual(
    seed: ShootingSeed,
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
    vinf_anchors: Mapping[str, float] | None = None,
    vinf_weight: float = 0.0,
) -> NDArray[np.float64]:
    """Full-state multiple-shooting residual in restricted n-body (design §3).

    Components, concatenated:
      1. **Leg defects** — for each interior leg ``i``, propagate node ``i``'s
         full state from ``epoch_i`` for ``tof_i`` in the rails n-body and subtract
         node ``i+1``'s state: 6 components per leg (km on position, km/s on
         velocity). Self-consistent nodes (sampled from one arc) give ~0.
      2. **Flyby hinges** — one per interior encounter: the
         :func:`_flyby_periapsis_hinge_km` shortfall (km), keeping the corrected
         flybys bend-feasible (the ``correct.py`` bend lesson, moved into the
         residual).
      3. **Periodicity wrap** — the position+velocity gap between the wrap node and
         the period-shifted home node (the cycler must repeat).

    Divergence is a first-class outcome: a leg whose propagation does not converge
    contributes a large finite sentinel defect (never a NaN/raise), so
    ``least_squares`` sees it as a bad region, not a crash (mirror ``correct.py``).

    ``max_wall_sec`` is the per-leg propagation wall budget (the propagator's own
    divergence budget). The default (90 s) matches the propagator; a multiple-
    shooting solve that re-propagates many legs over a high-V∞ / step-collapsing
    seed (the #135 high-V∞ basin) needs a tighter budget to stay bounded — every
    leg that blows the budget returns the divergence sentinel instead of grinding
    on. Pass a small value (e.g. 8 s) for a wall-capped bounded solve.
    """
    from cyclerfinder.nbody.forces import RailsEphemerisCache
    from cyclerfinder.nbody.propagator import RestrictedNBody

    prop = RestrictedNBody("rebound")
    n = len(seed.sequence)
    res: list[float] = []

    # Build ONE rails cache spanning the whole itinerary and reuse it across every
    # leg propagation — the spline build (~0.5 s/call) is the dominant per-call
    # cost, and a multiple-shooting Jacobian re-propagates the same window many
    # times. Without reuse a single 3-node shoot does not finish in 8 min.
    bodies = tuple(bodies)
    cache = (
        RailsEphemerisCache(bodies, ephem, min(seed.epochs), max(seed.epochs)) if bodies else None
    )

    # 1. Leg continuity defects.
    for i in range(n - 1):
        s_i = seed.node_states[i]
        r0 = np.asarray(s_i[:3], dtype=np.float64)
        v0 = np.asarray(s_i[3:], dtype=np.float64)
        t0 = seed.epochs[i]
        t1 = seed.epochs[i + 1]
        arc = prop.propagate(
            r0,
            v0,
            t0_sec=t0,
            t1_sec=t1,
            bodies=bodies,
            accuracy=accuracy,
            ephem=ephem,
            cache=cache,
            max_wall_sec=max_wall_sec,
        )
        s_next = seed.node_states[i + 1]
        if arc.converged and np.all(np.isfinite(arc.r_km)) and np.all(np.isfinite(arc.v_km_s)):
            dr = arc.r_km - np.asarray(s_next[:3], dtype=np.float64)
            dv = arc.v_km_s - np.asarray(s_next[3:], dtype=np.float64)
            res.extend(float(x) for x in dr)
            res.extend(float(x) for x in dv)
        else:
            # Honest divergence sentinel (large, finite) — never a NaN/exception.
            res.extend([1e9] * _STATE_DIM)

    # 2. Flyby bend-feasibility hinges (interior encounters only).
    for i in range(1, n - 1):
        res.append(_flyby_periapsis_hinge_km(seed.vinf_in[i], seed.vinf_out[i], seed.sequence[i]))

    # 3. Periodicity wrap: home node repeats one period later. Compared in the
    #    home-body-relative frame so a pure ephemeris shift of the home body is not
    #    spuriously charged as a defect.
    r_home, v_home = ephem.state(seed.sequence[0], seed.epochs[0])
    r_wrap_pl, v_wrap_pl = ephem.state(seed.sequence[-1], seed.epochs[-1])
    s0 = seed.node_states[0]
    sn = seed.node_states[-1]
    rel0_r = np.asarray(s0[:3], dtype=np.float64) - np.asarray(r_home, dtype=np.float64)
    rel0_v = np.asarray(s0[3:], dtype=np.float64) - np.asarray(v_home, dtype=np.float64)
    reln_r = np.asarray(sn[:3], dtype=np.float64) - np.asarray(r_wrap_pl, dtype=np.float64)
    reln_v = np.asarray(sn[3:], dtype=np.float64) - np.asarray(v_wrap_pl, dtype=np.float64)
    res.extend(float(x) for x in (reln_r - rel0_r))
    res.extend(float(x) for x in (reln_v - rel0_v))

    # 4. V∞-anchor penalty (family-pinned homotopy, #388). Opt-in: only when
    #    anchors are supplied AND the weight is positive. Each penalty row pulls a
    #    node's V∞ magnitude toward the SOURCED anchor for its body; the homotopy
    #    driver ramps vinf_weight to zero so the recorded V∞ emerges from the
    #    unpenalized (weight=0) solve. Appended after leg/hinge/wrap, in node order.
    if vinf_anchors and vinf_weight > 0.0:
        sw = float(np.sqrt(vinf_weight))
        for i, body in enumerate(seed.sequence):
            anchor = vinf_anchors.get(body)
            if anchor is None:
                continue
            s_i = seed.node_states[i]
            _, v_pl = ephem.state(body, seed.epochs[i])
            vinf_mag = float(
                np.linalg.norm(
                    np.asarray(s_i[3:], dtype=np.float64) - np.asarray(v_pl, dtype=np.float64)
                )
            )
            res.append(sw * (vinf_mag - float(anchor)))

    return np.asarray(res, dtype=np.float64)


@dataclass(frozen=True)
class NearMissSeed:
    """A low-V∞ near-miss conic chain — a shooting seed (Phase 3a, #133).

    The #135 verdict mandates the shooter be seeded from the near-miss survey, NOT
    a blind scan: the lowest-V∞ conic chains found near the Jones anchors, even
    when they carry a small V∞-continuity mismatch (Jones's own stage 1 accepted
    <= 200 m/s conic mismatches, then SNOPT-corrected to ballistic — the deepdive
    note §2.2). ``max_residual_kms`` is the conic-continuity residual; the seed is
    retained even when it exceeds the strict 0.1 km/s closure floor (down to a
    relaxed near-miss tolerance), because the shooter's job is to absorb that
    residual into n-body-ballistic continuity.
    """

    t0_sec: float
    tof_days: tuple[float, ...]
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    slack_leg: int
    max_residual_kms: float
    vinf_per_encounter_kms: tuple[float, ...]
    max_vinf_kms: float
    bend_feasible: bool


def near_miss_survey(
    *,
    sequence: tuple[str, ...],
    period_sec: float,
    t0_base_sec: float,
    tof_seed_days: Sequence[float],
    ephem: Ephemeris,
    n_epochs: int = 64,
    vinf_cap: float = 8.0,
    near_miss_tol_kms: float = 0.5,
    branch_topologies: Sequence[tuple[tuple[int, ...], tuple[str, ...]]] | None = None,
) -> list[NearMissSeed]:
    """Phase 3a near-miss survey: lowest-V∞ conic chains near the Jones anchors.

    Scans an epoch grid over one repeat period (centred on ``t0_base_sec``) x the
    given rev/branch topologies, running the conic ``ballistic_correct`` at each
    point, and collects every chain whose continuity residual is within
    ``near_miss_tol_kms`` (the relaxed near-miss tolerance — Jones's <=200 m/s
    stage-1 analogue, here 0.5 km/s), ranked by max per-encounter V∞ (lowest first
    = closest to the Jones 2.5-3.9 basin). These are the shooter seeds (#135
    verdict: near-miss seeding, never a blind scan).
    """
    from cyclerfinder.core.constants import SECONDS_PER_DAY as _SPD
    from cyclerfinder.search.correct import ballistic_correct

    n_legs = len(sequence) - 1
    slack_leg = int(np.argmax(tof_seed_days)) if len(tof_seed_days) else 0
    free_seed = [t for i, t in enumerate(tof_seed_days) if i != slack_leg]
    if branch_topologies is None:
        branch_topologies = [((0,) * n_legs, ("single",) * n_legs)]

    period_days = period_sec / _SPD
    step = period_days / max(1, n_epochs)
    t0_seeds = [t0_base_sec + i * step * _SPD for i in range(n_epochs)]

    seeds: list[NearMissSeed] = []
    for revs, branch in branch_topologies:
        for t0 in t0_seeds:
            result = ballistic_correct(
                sequence,
                revs,
                branch,
                t0,
                free_seed,
                period_sec,
                ephem,
                vinf_cap=vinf_cap,
                slack_leg=slack_leg,
                tol_kms=near_miss_tol_kms,
            )
            if not result.vinf_per_encounter_kms:
                continue
            if result.max_residual_kms > near_miss_tol_kms:
                continue
            seeds.append(
                NearMissSeed(
                    t0_sec=result.t0_sec,
                    tof_days=result.tof_days,
                    per_leg_revs=revs,
                    per_leg_branch=branch,
                    slack_leg=slack_leg,
                    max_residual_kms=result.max_residual_kms,
                    vinf_per_encounter_kms=result.vinf_per_encounter_kms,
                    max_vinf_kms=max(result.vinf_per_encounter_kms),
                    bend_feasible=result.bend_feasible,
                )
            )
    seeds.sort(key=lambda s: s.max_vinf_kms)
    return seeds


def shooting_seed_from_near_miss(
    near_miss: NearMissSeed,
    sequence: tuple[str, ...],
    period_sec: float,
    ephem: Ephemeris,
) -> ShootingSeed:
    """Build the multiple-shooting seed for a near-miss chain (reuses _vinf_nodes)."""
    from cyclerfinder.core.constants import SECONDS_PER_DAY as _SPD
    from cyclerfinder.search.correct import _vinf_nodes

    period_days = period_sec / _SPD
    free_tofs = [t for i, t in enumerate(near_miss.tof_days) if i != near_miss.slack_leg]
    nodes = _vinf_nodes(
        sequence=sequence,
        per_leg_revs=near_miss.per_leg_revs,
        per_leg_branch=near_miss.per_leg_branch,
        t0_sec=near_miss.t0_sec,
        free_tof_days=free_tofs,
        slack_leg=near_miss.slack_leg,
        period_days=period_days,
        ephem=ephem,
    )
    return seed_from_conic(
        sequence=sequence,
        vinf_nodes=nodes,
        t0_sec=near_miss.t0_sec,
        tofs_days=near_miss.tof_days,
        slack_leg=near_miss.slack_leg,
        period_days=period_days,
        ephem=ephem,
    )


@dataclass(frozen=True)
class ShootResult:
    """Frozen multiple-shooting result (design §3 / §5 honest-record discipline).

    ``converged`` follows ``scipy``'s success flag AND a defect-norm acceptance
    floor; a non-converged solve is recorded honestly (mirror ``correct.py``), not
    raised. ``vinf_per_encounter_kms`` is the per-node V∞ magnitude of the
    *corrected* solution (the side the Jones gate compares to the SOURCED multiset).
    ``correction_dv`` is the node-impulse ΔV from seed to corrected (Task B.2).
    """

    converged: bool
    defect_norm: float
    seed_defect_norm: float
    corrected_states: list[Vec3]
    vinf_per_encounter_kms: list[float]
    correction_dv_kms: float
    bend_feasible: bool
    sequence: tuple[str, ...]
    integrator_accuracy: float
    n_iterations: int


# SNOPT-analogue continuity acceptance (Jones AAS 17-577 method deep-dive §2.5:
# "continuity to 1.0E-3 km position and 1.0E-6 km/s velocity"). We accept a solve
# as ballistic-converged when its full-state defect L2 norm falls below a scaled
# multiple of these published per-component tolerances. Documented provenance.
_POS_CONTINUITY_KM = 1.0e-3
_VEL_CONTINUITY_KMS = 1.0e-6


def _states_to_x(states: Sequence[Vec3]) -> NDArray[np.float64]:
    return np.concatenate([np.asarray(s, dtype=np.float64).ravel() for s in states])


def _x_to_states(x: NDArray[np.float64], n: int) -> list[Vec3]:
    return [
        np.asarray(x[i * _STATE_DIM : (i + 1) * _STATE_DIM], dtype=np.float64) for i in range(n)
    ]


def _seed_with_states(seed: ShootingSeed, states: Sequence[Vec3]) -> ShootingSeed:
    """Clone ``seed`` with new node states (epochs/ToFs/V∞ vectors carried)."""
    return ShootingSeed(
        node_states=[np.asarray(s, dtype=np.float64) for s in states],
        epochs=list(seed.epochs),
        tofs=list(seed.tofs),
        sequence=seed.sequence,
        slack_leg=seed.slack_leg,
        period_days=seed.period_days,
        vinf_in=list(seed.vinf_in),
        vinf_out=list(seed.vinf_out),
    )


def _node_vinf(state: Vec3, body: str, epoch: float, ephem: Ephemeris) -> Vec3:
    _, v_pl = ephem.state(body, epoch)
    return np.asarray(state[3:], dtype=np.float64) - np.asarray(v_pl, dtype=np.float64)


# --- Parallel finite-difference Jacobian (perf lever 1, #159) ------------------
#
# A multiple-shooting Jacobian is (3*n_nodes+1) residual evaluations, each of
# which re-propagates the full multi-segment trajectory; serially this dominates
# the LM solve cost (~0.5-5 CPU-h/member, the #159 cost model). The columns are
# independent, so they parallelise across a process pool. We use the SAME
# primitives-only worker pattern the scan rung established
# (``search/scan.py``): only the small ``ShootingSeed`` (numpy arrays + floats,
# fully picklable) and the perturbed free-vector cross the process boundary; each
# worker constructs its own ``Ephemeris(model=...)`` once via the pool
# initialiser, never pickling a live ``Ephemeris``. Threads are NOT used: the
# rails additional-force callback is a Python CFUNCTYPE invoked on every IAS15
# force evaluation, so REBOUND re-acquires the GIL inside the inner loop and
# threads would serialise there (measured; see the perf note).

# Per-worker pinned context (one binding per worker process — process-safe, NOT
# thread-safe), set by :func:`_shoot_init_worker`. The parent passes the seed and
# the ephemeris MODEL STRING (never a live Ephemeris) through the initialiser.
_WORKER_SEED: ShootingSeed | None = None
_WORKER_EPHEM: Ephemeris | None = None
_WORKER_BODIES: tuple[str, ...] = ()
_WORKER_ACCURACY: float = 1e-10
_WORKER_MAX_WALL_SEC: float = 90.0


def _shoot_init_worker(
    seed: ShootingSeed,
    ephem_model: str,
    bodies: tuple[str, ...],
    accuracy: float,
    max_wall_sec: float,
) -> None:
    """Pool initialiser: pin the seed + a fresh per-process Ephemeris."""
    global _WORKER_SEED, _WORKER_EPHEM, _WORKER_BODIES, _WORKER_ACCURACY, _WORKER_MAX_WALL_SEC
    from cyclerfinder.core.ephemeris import Ephemeris

    _WORKER_SEED = seed
    _WORKER_EPHEM = Ephemeris(model=ephem_model)
    _WORKER_BODIES = tuple(bodies)
    _WORKER_ACCURACY = float(accuracy)
    _WORKER_MAX_WALL_SEC = float(max_wall_sec)


def _residual_worker(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Worker entry point: evaluate the full-state residual at free-vector ``x``.

    Runs in a worker process against the pinned seed/ephemeris. Only the
    perturbed free-vector ``x`` crossed the boundary; the result is the small
    residual vector. Mirrors :func:`_evaluate_point` in ``search/scan.py``.
    """
    assert _WORKER_SEED is not None and _WORKER_EPHEM is not None
    n = len(_WORKER_SEED.sequence)
    states = _x_to_states(x, n)
    trial = _seed_with_states(_WORKER_SEED, states)
    return defect_residual(
        trial,
        ephem=_WORKER_EPHEM,
        bodies=_WORKER_BODIES,
        accuracy=_WORKER_ACCURACY,
        max_wall_sec=_WORKER_MAX_WALL_SEC,
    )


def _fd_step(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Forward-difference step per element (MINPACK ``lmdif`` convention).

    ``h_j = sqrt(eps) * |x_j|``, with ``sqrt(eps)`` substituted when ``x_j`` is
    zero — the same relative step ``scipy``'s ``method="lm"`` uses internally, so
    the parallel Jacobian matches the serial reference to working precision.
    """
    rel = float(np.sqrt(np.finfo(np.float64).eps))
    h = rel * np.abs(x)
    h[h == 0.0] = rel
    return h


def _fd_jacobian(
    residual_of_x: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    x: NDArray[np.float64],
    f0: NDArray[np.float64],
    *,
    column_eval: Callable[
        [Callable[[NDArray[np.float64]], NDArray[np.float64]], list[NDArray[np.float64]]],
        list[NDArray[np.float64]],
    ],
) -> NDArray[np.float64]:
    """Forward-difference Jacobian ``df/dx`` (columns optionally parallelised).

    ``column_eval`` maps the residual over the per-column perturbed vectors;
    passing a serial map or a pool map selects the execution mode. The arithmetic
    (step size, forward difference) is identical either way, so the serial and
    parallel Jacobians agree to working precision — asserted in the tests.
    """
    h = _fd_step(x)
    perturbed = []
    for j in range(x.size):
        xj = x.copy()
        xj[j] += h[j]
        perturbed.append(xj)
    fcols = column_eval(residual_of_x, perturbed)
    jac = np.empty((f0.size, x.size), dtype=np.float64)
    for j in range(x.size):
        jac[:, j] = (fcols[j] - f0) / h[j]
    return jac


def _serial_columns(
    residual_of_x: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    perturbed: list[NDArray[np.float64]],
) -> list[NDArray[np.float64]]:
    """Serial column evaluator — the determinism oracle for the parallel path."""
    return [residual_of_x(xj) for xj in perturbed]


# --- Analytic block-bidiagonal STM Jacobian (perf lever, #388) -----------------
#
# The FD Jacobian is (6*n_nodes+1) full residual re-propagations per LM step; on
# the multi-year, multi-rev SnLm cyclers that is the compute wall (a single FD
# shoot >400 s no-return). The defects are full-state continuity, so the Jacobian
# is exactly known from per-leg state-transition matrices — ONE augmented
# (variational) propagation per leg replaces the (6*n_nodes+1) FD re-propagations.
# FD stays the default + the parity oracle; this is the opt-in `jacobian="stm"`.


def _stm_jacobian(
    seed: ShootingSeed,
    x: NDArray[np.float64],
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
    vinf_anchors: Mapping[str, float] | None = None,
    vinf_weight: float = 0.0,
) -> NDArray[np.float64]:
    """Analytic Jacobian of :func:`defect_residual` w.r.t. the node-state vector.

    The free variables are the per-node full Cartesian states (6 each, the same
    packing as :func:`_states_to_x`); the residual is, in order, the ``n-1`` leg
    continuity defects (6 each), the ``n-2`` interior flyby hinges (1 each), and
    the 6-component periodicity wrap. The Jacobian is **block-bidiagonal**:

    - Leg defect ``c_i = propagate(node_i, leg_i) - node_{i+1}`` gives
      ``dc_i/dnode_i = Phi_i`` (the 6x6 per-leg STM, co-integrated in ONE
      augmented propagation via ``propagate(with_stm=True)``) and
      ``dc_i/dnode_{i+1} = -I_6``; every other block is zero.
    - The flyby hinges read ``seed.vinf_in/out`` (carried constants, never
      recomputed from the node states), so their rows are identically zero —
      matching what FD sees (a small node-state perturbation leaves them fixed).
    - The wrap residual ``(node_{n-1} - r_wrap_pl) - (node_0 - r_home)`` gives
      ``-I_6`` on node 0 and ``+I_6`` on node ``n-1`` (the planet states are
      epoch-fixed constants).

    A leg that does not converge (divergence sentinel in the residual) leaves its
    ``Phi_i`` block zero — the sentinel is a constant in ``x`` there, so a zero
    derivative is the honest local linearisation (and the FD oracle agrees: the
    constant sentinel differences to ~0). ``-I_6`` on the next node is kept so the
    defect still couples the nodes.
    """
    from cyclerfinder.nbody.forces import RailsEphemerisCache
    from cyclerfinder.nbody.propagator import RestrictedNBody

    n = len(seed.sequence)
    states = _x_to_states(x, n)
    bodies = tuple(bodies)
    prop = RestrictedNBody("rebound")
    cache = (
        RailsEphemerisCache(bodies, ephem, min(seed.epochs), max(seed.epochs)) if bodies else None
    )

    n_leg = (n - 1) * _STATE_DIM
    n_hinge = max(0, n - 2)
    # Penalty node indices (those whose body has an anchor), in node order — must
    # match the residual's penalty-row order in defect_residual (block 4).
    pen_nodes: list[int] = []
    if vinf_anchors and vinf_weight > 0.0:
        pen_nodes = [i for i, b in enumerate(seed.sequence) if b in vinf_anchors]
    n_pen = len(pen_nodes)
    n_rows = n_leg + n_hinge + _STATE_DIM + n_pen
    n_cols = n * _STATE_DIM
    jac = np.zeros((n_rows, n_cols), dtype=np.float64)
    eye6 = np.eye(_STATE_DIM)

    for i in range(n - 1):
        s_i = states[i]
        r0 = np.asarray(s_i[:3], dtype=np.float64)
        v0 = np.asarray(s_i[3:], dtype=np.float64)
        arc = prop.propagate(
            r0,
            v0,
            t0_sec=seed.epochs[i],
            t1_sec=seed.epochs[i + 1],
            bodies=bodies,
            accuracy=accuracy,
            ephem=ephem,
            cache=cache,
            max_wall_sec=max_wall_sec,
            with_stm=True,
        )
        rows = slice(i * _STATE_DIM, (i + 1) * _STATE_DIM)
        if arc.converged and arc.stm is not None and np.all(np.isfinite(arc.stm)):
            jac[rows, i * _STATE_DIM : (i + 1) * _STATE_DIM] = arc.stm
        jac[rows, (i + 1) * _STATE_DIM : (i + 2) * _STATE_DIM] = -eye6

    # Hinge rows (n_leg : n_leg + n_hinge) are constant in x -> left zero.

    # Periodicity wrap rows.
    wrap = slice(n_leg + n_hinge, n_leg + n_hinge + _STATE_DIM)
    jac[wrap, 0:_STATE_DIM] = -eye6
    jac[wrap, (n - 1) * _STATE_DIM : n * _STATE_DIM] = eye6

    # Penalty rows: d(sqrt(w)*|v∞_i|)/d(node_i velocity) = sqrt(w) * v̂∞_i, zero on
    # the position block and on every other node. Matches the residual's penalty
    # rows (defect_residual block 4). A near-zero-V∞ node -> zero row.
    if n_pen:
        sw = float(np.sqrt(vinf_weight))
        base = n_leg + n_hinge + _STATE_DIM
        for k, i in enumerate(pen_nodes):
            s_i = states[i]
            _, v_pl = ephem.state(seed.sequence[i], seed.epochs[i])
            dv = np.asarray(s_i[3:], dtype=np.float64) - np.asarray(v_pl, dtype=np.float64)
            mag = float(np.linalg.norm(dv))
            if mag > 1e-12:
                vhat = dv / mag
                jac[base + k, i * _STATE_DIM + 3 : i * _STATE_DIM + 6] = sw * vhat

    return jac


def _parallel_columns_for_test(
    seed: ShootingSeed,
    ephem: Ephemeris,
    x: NDArray[np.float64],
    f0: NDArray[np.float64],
    *,
    n_jobs: int,
    bodies: Sequence[str] = (),
    accuracy: float = 1e-11,
    max_wall_sec: float = 90.0,
) -> NDArray[np.float64]:
    """Build the FD Jacobian via the process pool (correctness-test entry point).

    Exercises the exact pool path :func:`shoot` uses for ``n_jobs > 1`` in
    isolation (one Jacobian build, no LM loop), so the parallel-vs-serial
    equality gate runs in milliseconds. Not part of the solve API.
    """
    from concurrent.futures import ProcessPoolExecutor

    n = len(seed.sequence)

    def residual_of_x(xv: NDArray[np.float64]) -> NDArray[np.float64]:
        trial = _seed_with_states(seed, _x_to_states(xv, n))
        return defect_residual(
            trial, ephem=ephem, bodies=bodies, accuracy=accuracy, max_wall_sec=max_wall_sec
        )

    workers = max(1, min(int(n_jobs), x.size))
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_shoot_init_worker,
        initargs=(seed, ephem.model, tuple(bodies), accuracy, max_wall_sec),
    ) as pool:

        def _cols(
            _res: Callable[[NDArray[np.float64]], NDArray[np.float64]],
            perturbed: list[NDArray[np.float64]],
        ) -> list[NDArray[np.float64]]:
            return list(pool.map(_residual_worker, perturbed))

        return _fd_jacobian(residual_of_x, x, f0, column_eval=_cols)


def shoot(
    seed: ShootingSeed,
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    accuracy: float = 1e-10,
    max_nfev: int = 200,
    max_wall_sec: float = 90.0,
    n_jobs: int = 1,
    jacobian: Literal["fd", "stm"] = "fd",
    progress: Callable[[str, int, float, float], None] | None = None,
    vinf_anchors: Mapping[str, float] | None = None,
    vinf_weight: float = 0.0,
) -> ShootResult:
    """Multiple-shooting differential correction (the SNOPT analogue, design §3).

    Free variables: the per-node full Cartesian states (6 each). Node epochs are
    held at the (near-miss-seeded) seed epochs — a documented choice: the conic
    seed already fixes good encounter dates near the Jones anchors, and freeing
    epochs multiplies the finite-difference Jacobian cost without changing the
    family selection the #135 verdict identified as the real lever (seeding, not
    DOF count). Solver: ``scipy.optimize.least_squares(method="lm")`` (mirror
    ``correct.py``), Jacobian by finite difference over REBOUND (Q1 baseline).

    A solve is ``converged`` when ``least_squares`` reports success AND the
    corrected full-state defect norm clears the published SNOPT continuity floor
    (§2.5: 1e-3 km / 1e-6 km/s). Divergence is recorded honestly, not raised.

    ``max_wall_sec`` is the per-leg propagation wall budget, forwarded to
    :func:`defect_residual`. A bounded (tractable-gate) solve combines a small
    ``max_nfev`` with a small ``max_wall_sec`` so the worst case is
    ``O(max_nfev * (n-1) * max_wall_sec)`` wall — divergent legs short-circuit on
    the budget instead of collapsing IAS15 steps for the full 90 s default.

    ``n_jobs`` (perf lever 1, #159): worker-process count for the
    finite-difference Jacobian columns. ``n_jobs == 1`` (the default) is today's
    serial path EXACTLY — ``least_squares`` builds its own internal FD Jacobian,
    byte-identical to the pre-#159 behaviour. ``n_jobs > 1`` supplies an explicit
    forward-difference Jacobian whose (3*n_nodes+1) independent columns are
    evaluated across a :class:`~concurrent.futures.ProcessPoolExecutor` (same
    primitives-only worker contract as ``search/scan.py``); the FD arithmetic
    matches ``scipy``'s ``lm`` step to working precision. The serial and parallel
    Jacobians agree on a fixture to ~1e-9 relative (asserted in the tests). The
    speedup is bounded by ``min(n_jobs, 3*n_nodes+1)``.

    ``jacobian`` (perf lever, #388): ``"fd"`` (default) is the finite-difference
    path above (serial scipy-internal FD for ``n_jobs == 1``, the parallel FD pool
    for ``n_jobs > 1``) and stays the parity oracle. ``"stm"`` supplies the analytic
    block-bidiagonal Jacobian assembled from per-leg state-transition matrices
    (:func:`_stm_jacobian`) — ONE co-integrated variational propagation per leg
    instead of the ``6*n_nodes+1`` FD re-propagations, the lever that makes the
    multi-year multi-rev shoot tractable. The STM and FD Jacobians agree to ~5e-3
    relative on the fixture (asserted in the tests); ``"stm"`` ignores ``n_jobs``
    (it is already one propagation per leg, nothing to fan out).

    ``progress`` (optional heartbeat): if given, it is called as
    ``progress(kind, count, defect_norm, elapsed_sec)`` after each residual
    evaluation (``kind="f"``) and each STM-Jacobian build (``kind="J"``, with
    ``defect_norm`` carried from the last residual), where ``elapsed_sec`` is the
    wall time that single evaluation took. This gives a long detached solve a
    per-evaluation intermediate-output stream (every ~residual/Jacobian, i.e.
    seconds-to-minutes) instead of only a per-row record at the very end.
    """
    import time as _time

    from scipy.optimize import least_squares

    n = len(seed.sequence)
    if vinf_anchors and vinf_weight > 0.0 and n_jobs > 1:
        raise ValueError(
            "vinf penalty (family-pinned homotopy) is supported only with n_jobs=1 "
            "(the parallel FD worker path does not carry anchors); use jacobian='stm'."
        )
    _f_count = [0]
    _j_count = [0]
    _last_norm = [float("inf")]

    def residual_of_x(x: NDArray[np.float64]) -> NDArray[np.float64]:
        t_eval = _time.monotonic()
        states = _x_to_states(x, n)
        trial = _seed_with_states(seed, states)
        res = defect_residual(
            trial,
            ephem=ephem,
            bodies=bodies,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
            vinf_anchors=vinf_anchors,
            vinf_weight=vinf_weight,
        )
        if progress is not None:
            _f_count[0] += 1
            norm = float(np.linalg.norm(res))
            _last_norm[0] = norm
            progress("f", _f_count[0], norm, _time.monotonic() - t_eval)
        return res

    x0 = _states_to_x(seed.node_states)
    seed_res = residual_of_x(x0)
    seed_defect_norm = float(np.linalg.norm(seed_res))

    if jacobian == "stm":

        def jac_of_x_stm(x: NDArray[np.float64]) -> NDArray[np.float64]:
            t_eval = _time.monotonic()
            jac = _stm_jacobian(
                seed,
                x,
                ephem=ephem,
                bodies=bodies,
                accuracy=accuracy,
                max_wall_sec=max_wall_sec,
                vinf_anchors=vinf_anchors,
                vinf_weight=vinf_weight,
            )
            if progress is not None:
                _j_count[0] += 1
                progress("J", _j_count[0], _last_norm[0], _time.monotonic() - t_eval)
            return jac

        sol = least_squares(residual_of_x, x0, jac=jac_of_x_stm, method="lm", max_nfev=max_nfev)
    elif n_jobs <= 1:
        sol = least_squares(residual_of_x, x0, method="lm", max_nfev=max_nfev)
    else:
        from concurrent.futures import ProcessPoolExecutor

        n_cols = x0.size
        workers = max(1, min(int(n_jobs), n_cols))
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_shoot_init_worker,
            initargs=(seed, ephem.model, tuple(bodies), accuracy, max_wall_sec),
        ) as pool:

            def _parallel_columns(
                _res: Callable[[NDArray[np.float64]], NDArray[np.float64]],
                perturbed: list[NDArray[np.float64]],
            ) -> list[NDArray[np.float64]]:
                # The residual is evaluated in the workers against the pinned seed;
                # only the perturbed free-vectors cross the boundary.
                return list(pool.map(_residual_worker, perturbed))

            def jac_of_x(x: NDArray[np.float64]) -> NDArray[np.float64]:
                f0 = residual_of_x(x)
                return _fd_jacobian(residual_of_x, x, f0, column_eval=_parallel_columns)

            sol = least_squares(residual_of_x, x0, jac=jac_of_x, method="lm", max_nfev=max_nfev)
    corrected_states = _x_to_states(np.asarray(sol.x, dtype=np.float64), n)
    final_res = residual_of_x(np.asarray(sol.x, dtype=np.float64))
    defect_norm = float(np.linalg.norm(final_res))

    # Only the leg-continuity block governs the ballistic-continuity acceptance
    # (the trailing hinge + wrap terms are constraints, not continuity defects).
    n_leg_defects = (n - 1) * _STATE_DIM
    leg_block = final_res[:n_leg_defects]
    continuity_floor = float(
        np.linalg.norm([_POS_CONTINUITY_KM] * 3 + [_VEL_CONTINUITY_KMS] * 3) * np.sqrt(n - 1)
    )
    converged = bool(sol.success) and float(np.linalg.norm(leg_block)) < continuity_floor

    corrected_vinf = [
        _node_vinf(s, body, ep, ephem)
        for s, body, ep in zip(corrected_states, seed.sequence, seed.epochs, strict=True)
    ]
    vinf_mag = [float(np.linalg.norm(v)) for v in corrected_vinf]

    # Correction ΔV: per-node velocity change seed -> corrected (Task B.2).
    seed_nodes = {f"b{i}": np.asarray(seed.node_states[i][3:]) for i in range(n)}
    corr_nodes = {f"b{i}": np.asarray(corrected_states[i][3:]) for i in range(n)}
    from cyclerfinder.nbody.correction_dv import node_impulse_correction_dv

    cdv = node_impulse_correction_dv(seed_nodes, corr_nodes)

    # Bend feasibility of the interior flybys, evaluated with the Jones B-plane
    # kernel (AAS 17-577 Eq. 7, #142): the targeted flyby turns the SEED incoming
    # asymptote (vinf_in, end of the prior conic leg) into the SEED outgoing
    # asymptote (vinf_out, start of the next), and that powered flyby must clear
    # the published v∞-mismatch tolerance + 100 km-100,000 km altitude window.
    # (The corrected n-body node carries a single v∞; the conic seed's distinct
    # in/out asymptotes are the physically meaningful bend the flyby must deliver.)
    from cyclerfinder.nbody.bplane import interior_flyby_feasible

    bend_feasible = all(
        interior_flyby_feasible(seed.vinf_in[i], seed.vinf_out[i], seed.sequence[i])
        for i in range(1, n - 1)
    )

    return ShootResult(
        converged=converged,
        defect_norm=defect_norm,
        seed_defect_norm=seed_defect_norm,
        corrected_states=corrected_states,
        vinf_per_encounter_kms=vinf_mag,
        correction_dv_kms=cdv.total_kms,
        bend_feasible=bend_feasible,
        sequence=seed.sequence,
        integrator_accuracy=accuracy,
        n_iterations=int(sol.nfev),
    )


__all__ = [
    "NearMissSeed",
    "ShootResult",
    "ShootingSeed",
    "build_shooting_vector",
    "defect_count",
    "defect_residual",
    "near_miss_survey",
    "seed_from_conic",
    "shoot",
    "shooting_seed_from_near_miss",
]
