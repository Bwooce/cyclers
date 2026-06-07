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

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

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
    """
    from cyclerfinder.nbody.propagator import RestrictedNBody

    prop = RestrictedNBody("rebound")
    n = len(seed.sequence)
    res: list[float] = []

    # 1. Leg continuity defects.
    for i in range(n - 1):
        s_i = seed.node_states[i]
        r0 = np.asarray(s_i[:3], dtype=np.float64)
        v0 = np.asarray(s_i[3:], dtype=np.float64)
        t0 = seed.epochs[i]
        t1 = seed.epochs[i + 1]
        arc = prop.propagate(
            r0, v0, t0_sec=t0, t1_sec=t1, bodies=bodies, accuracy=accuracy, ephem=ephem
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

    return np.asarray(res, dtype=np.float64)


__all__ = [
    "ShootingSeed",
    "build_shooting_vector",
    "defect_count",
    "defect_residual",
    "seed_from_conic",
]
