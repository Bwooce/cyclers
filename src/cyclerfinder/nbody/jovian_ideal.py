"""Ideal circular-coplanar Galilean model — the λ=0 endpoint of the EGGIE homotopy.

This is the make-or-break positive control for the #480 EGGIE reproduction (Stage 2
of ``docs/superpowers/plans/2026-06-29-480-eggie-resonant-conic-generator-plan.md``).
The Stage-1 resonant-conic generator (:mod:`cyclerfinder.search.resonant_conic`)
proved the single conic at e≈0.62 puts all three Galilean V∞ on the Table-4 targets
(Europa 9.12, Ganymede 7.07, Io 8.38 km/s) — the diagnosis spike-4 crux. The
question Stage 2 answers: does the Jupiter-central multiple-shooting corrector
(:func:`cyclerfinder.nbody.jovian.jovian_defect_residual`) close the ideal-model
EGGIE tour near-ballistic when seeded from that family-correct conic guess?

To answer it without the real-ephemeris off-basin confound, this module supplies the
paper's IDEAL model (p.3): the Galilean moons are circular + coplanar at the
:func:`cyclerfinder.search.resonant_conic.ideal_moon_smas` radii, with inertial
phases set so each moon sits exactly at the Stage-1 guess's node position at the
guess epoch. :class:`IdealJovianEphemeris` is duck-typed like
:class:`cyclerfinder.nbody.jovian.JovianEphemeris` (``state(moon, t_sec)``) so it
drops straight into :class:`~cyclerfinder.nbody.jovian.JovianRailsCache`,
:class:`~cyclerfinder.nbody.jovian.JovianRestrictedNBody`, and
:func:`~cyclerfinder.nbody.jovian.jovian_defect_residual` unchanged.

Provenance (never from our own code; ``feedback_golden_tests_sourced_only``):

* Ideal model: Hernandez-Jones-Jesick 2017 (AAS 17-608) p.3 (moons circular +
  coplanar; ideal smas via the 5.2°-per-synodic resonance factors).
* Table-4 V∞ golden targets: Europa 9.12, Ganymede 7.07, Io 8.38 km/s
  (``docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md``).
"""

from __future__ import annotations

import math
import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2

if TYPE_CHECKING:
    from cyclerfinder.nbody.jovian import JovianEphemeris, JovianRailsCache
    from cyclerfinder.nbody.shooter import ShootingSeed, ShootResult
    from cyclerfinder.search.resonant_conic import EggieGuess

Vec3 = NDArray[np.float64]

#: The three Galilean moons of the EGGIE tour (the ideal model has no Callisto).
EGGIE_MOONS: tuple[str, ...] = ("Io", "Europa", "Ganymede")


class IdealJovianEphemeris:
    """Circular-coplanar Galilean moon ephemeris, duck-typed like ``JovianEphemeris``.

    Each moon orbits Jupiter on a circle of radius ``smas[moon]`` (km) in the z=0
    plane, prograde (counter-clockwise). Its inertial angle is
    ``theta = phase[moon] + n*(t_sec - anchor[moon])`` with ``n = sqrt(mu/a^3)`` —
    the exact convention of :func:`cyclerfinder.search.resonant_conic._tour_states`,
    so a moon placed from a Stage-1 guess sits on the guess's node positions at the
    guess epochs. This is the λ=0 endpoint of the ideal→real homotopy.
    """

    def __init__(
        self,
        smas: dict[str, float],
        phases: dict[str, float],
        anchors: dict[str, float],
        *,
        mu: float = MU_JUPITER_KM3_S2,
    ) -> None:
        self.smas = dict(smas)
        self.phases = dict(phases)
        self.anchors = dict(anchors)
        self.mu = mu
        self._n = {m: math.sqrt(mu / smas[m] ** 3) for m in smas}

    def state(self, moon: str, t_sec: float) -> tuple[Vec3, Vec3]:
        """Moon ``(r_km, v_km_s)`` relative to Jupiter, J2000 in-plane (z=0)."""
        a = self.smas[moon]
        theta = self.phases[moon] + self._n[moon] * (float(t_sec) - self.anchors[moon])
        c, s = math.cos(theta), math.sin(theta)
        r_vec = a * np.array([c, s, 0.0], dtype=np.float64)
        v_mag = math.sqrt(self.mu / a)
        v_vec = v_mag * np.array([-s, c, 0.0], dtype=np.float64)
        return r_vec, v_vec


def ideal_ephemeris_from_guess(
    guess: EggieGuess, *, t0_sec: float = 0.0, mu: float = MU_JUPITER_KM3_S2
) -> IdealJovianEphemeris:
    """Build the ideal ephemeris whose moons sit on ``guess``'s node positions.

    Phases come straight from the guess (each moon's inertial angle at its FIRST
    encounter); anchors are that first-encounter epoch on the seed's absolute axis
    (``t0_sec`` + the guess's cumulative encounter time). The smas are the ideal
    Galilean radii (paper p.3).
    """
    from cyclerfinder.search.resonant_conic import ideal_moon_smas

    smas = ideal_moon_smas()
    anchors: dict[str, float] = {}
    for k, moon in enumerate(guess.sequence):
        if moon not in anchors:
            anchors[moon] = t0_sec + guess.epochs_days[k] * SECONDS_PER_DAY
    return IdealJovianEphemeris(smas, dict(guess.moon_phases_rad), anchors, mu=mu)


def build_eggie_shooting_seed(
    guess: EggieGuess, *, t0_sec: float = 0.0, mu: float = MU_JUPITER_KM3_S2
) -> ShootingSeed:
    """Map a Stage-1 :class:`EggieGuess` to a :class:`ShootingSeed` (ideal model).

    Each node's spacecraft state is the conic node position (= the moon position by
    construction) with the guess's per-node spacecraft velocity. Epochs are the
    guess's cumulative encounter times on the ``t0_sec`` axis; ToFs and V∞ vectors
    are carried from the guess. The moon positions the corrector integrates against
    (via :func:`ideal_ephemeris_from_guess`) match these node positions by design.
    """
    from cyclerfinder.nbody.shooter import ShootingSeed

    n = len(guess.sequence)
    epochs = [t0_sec + d * SECONDS_PER_DAY for d in guess.epochs_days]
    node_states = [
        np.concatenate(
            [
                np.asarray(guess.node_positions[k], dtype=np.float64),
                np.asarray(guess.node_sc_velocities[k], dtype=np.float64),
            ]
        )
        for k in range(n)
    ]
    tofs = [float(d) for d in guess.tofs_days]
    slack_leg = int(np.argmax(tofs))
    return ShootingSeed(
        node_states=node_states,
        epochs=epochs,
        tofs=tofs,
        sequence=tuple(guess.sequence),
        slack_leg=slack_leg,
        period_days=float(sum(tofs)),
        vinf_in=[np.asarray(v, dtype=np.float64) for v in guess.vinf_in],
        vinf_out=[np.asarray(v, dtype=np.float64) for v in guess.vinf_out],
    )


def build_eggie_periapsis_seed(
    guess: EggieGuess, *, t0_sec: float = 0.0, mu: float = MU_JUPITER_KM3_S2
) -> ShootingSeed:
    """Map a Stage-1 guess to a ShootingSeed with nodes at flyby PERIAPSIS (not centers).

    The zero-SOI patched conic places each encounter AT the moon centre; propagating an
    n-body leg from a moon centre is a zero-altitude plunge (a ~1e6-km artefact that
    dominates the residual and is not a wrong-family signal — the 2-body legs close to
    ~1e-8 km). This builder instead converts each encounter to the moon-centred-hyperbola
    periapsis state (:func:`cyclerfinder.nbody.jovian.periapsis_node`), so the corrector
    starts on a physical flyby. Interior nodes use ``(vinf_in, vinf_out)``; the boundary
    Europa nodes (one V∞ is zero) use the single available V∞ as a degenerate zero-turn
    encounter — :func:`periapsis_node` then offsets them ~0.6 SOI off the moon where its
    gravity is negligible (a clean periodic boundary node).
    """
    from cyclerfinder.nbody.jovian import periapsis_node
    from cyclerfinder.nbody.shooter import ShootingSeed

    ephem = ideal_ephemeris_from_guess(guess, t0_sec=t0_sec, mu=mu)
    epochs = [t0_sec + d * SECONDS_PER_DAY for d in guess.epochs_days]
    vinf_in = [np.asarray(v, dtype=np.float64) for v in guess.vinf_in]
    vinf_out = [np.asarray(v, dtype=np.float64) for v in guess.vinf_out]

    node_states: list[Vec3] = []
    for k, moon in enumerate(guess.sequence):
        vin = vinf_in[k]
        vout = vinf_out[k]
        if float(np.linalg.norm(vin)) <= 0.0:  # departure boundary node
            vin = vout
        if float(np.linalg.norm(vout)) <= 0.0:  # arrival boundary node
            vout = vin
        # IdealJovianEphemeris is duck-typed (.state) to the JovianEphemeris API.
        r_p, v_p, _d = periapsis_node(moon, epochs[k], vin, vout, ephem)  # type: ignore[arg-type]
        node_states.append(np.concatenate([r_p, v_p]))

    tofs = [float(d) for d in guess.tofs_days]
    return ShootingSeed(
        node_states=node_states,
        epochs=epochs,
        tofs=tofs,
        sequence=tuple(guess.sequence),
        slack_leg=int(np.argmax(tofs)),
        period_days=float(sum(tofs)),
        vinf_in=vinf_in,
        vinf_out=vinf_out,
    )


@dataclass(frozen=True)
class SubarcSeed:
    """A sub-arc multiple-shooting seed: encounter nodes PLUS interior continuity nodes.

    Stage 4 of the #480 EGGIE corrector. Each original leg ``i`` (encounter ``i`` ->
    encounter ``i+1``) is split into ``n_subarcs`` sub-arcs by inserting
    ``n_subarcs - 1`` interior shooting nodes, sampled along the seed leg's actual
    propagated trajectory. The interior nodes are PURE continuity points (no flyby
    hinge, no V∞): only the original encounter nodes carry flyby hinges / participate
    in the periodicity wrap. Splitting each defect into ``n_subarcs`` shorter sub-arcs
    distributes the perijove sensitivity so each sub-arc's STM is well-conditioned —
    the lever for the ~0.1 km/s velocity-continuity plateau characterised in
    ``docs/notes/2026-06-29-480-eggie-stage3-stm-verdict.md``.

    ``node_states`` / ``epochs`` run over ALL ``M = (n-1)*n_subarcs + 1`` nodes in
    encounter order with interior nodes interleaved. ``encounter_idx`` are the ``n``
    indices into ``node_states`` that are real encounters; ``sequence`` / ``vinf_in``
    / ``vinf_out`` are the per-encounter constants carried from the base seed. With
    ``n_subarcs == 1`` it reduces to the encounter nodes themselves (``encounter_idx
    == range(n)``), so the residual/Jacobian match the one-node-per-leg path.
    """

    node_states: list[Vec3]
    epochs: list[float]
    encounter_idx: tuple[int, ...]
    sequence: tuple[str, ...]
    vinf_in: list[Vec3]
    vinf_out: list[Vec3]
    n_subarcs: int


def build_subarc_seed(
    seed: ShootingSeed,
    n_subarcs: int,
    *,
    ephem: IdealJovianEphemeris,
    moons: Sequence[str] = EGGIE_MOONS,
    rtol: float = 1e-11,
    atol: float = 1e-9,
) -> SubarcSeed:
    """Insert ``n_subarcs - 1`` interior continuity nodes per leg, on the seed arc.

    Each leg is propagated from its encounter node with
    :func:`cyclerfinder.nbody.jovian_stm.propagate_with_stm` (the same Jupiter-central
    dynamics the corrector integrates) and sampled at ``n_subarcs - 1`` evenly spaced
    interior epochs. The interior node states therefore lie on the seed leg's actual
    trajectory, so the interior sub-arc continuity defects start ~0 and the leg's
    full (conic-vs-n-body) defect remains concentrated at the original encounter-node
    boundary — which the corrector then redistributes across the now-shorter sub-arcs.
    """
    from cyclerfinder.nbody.jovian_stm import propagate_with_stm

    if n_subarcs < 1:
        raise ValueError(f"n_subarcs must be >= 1, got {n_subarcs}")
    n = len(seed.sequence)
    moons = tuple(moons)

    node_states: list[Vec3] = []
    epochs: list[float] = []
    encounter_idx: list[int] = []

    for i in range(n - 1):
        encounter_idx.append(len(node_states))
        node_states.append(np.asarray(seed.node_states[i], dtype=np.float64))
        epochs.append(float(seed.epochs[i]))
        if n_subarcs > 1:
            t_a = float(seed.epochs[i])
            t_b = float(seed.epochs[i + 1])
            dt = (t_b - t_a) / n_subarcs
            r = np.asarray(seed.node_states[i][:3], dtype=np.float64)
            v = np.asarray(seed.node_states[i][3:], dtype=np.float64)
            t_cur = t_a
            for j in range(1, n_subarcs):
                t_next = t_a + j * dt
                rf, vf, _ = propagate_with_stm(
                    r, v, t_cur, t_next, ephem=ephem, moons=moons, rtol=rtol, atol=atol
                )
                node_states.append(np.concatenate([rf, vf]))
                epochs.append(t_next)
                r, v, t_cur = rf, vf, t_next

    # Final encounter node.
    encounter_idx.append(len(node_states))
    node_states.append(np.asarray(seed.node_states[n - 1], dtype=np.float64))
    epochs.append(float(seed.epochs[n - 1]))

    return SubarcSeed(
        node_states=node_states,
        epochs=epochs,
        encounter_idx=tuple(encounter_idx),
        sequence=tuple(seed.sequence),
        vinf_in=[np.asarray(v, dtype=np.float64) for v in seed.vinf_in],
        vinf_out=[np.asarray(v, dtype=np.float64) for v in seed.vinf_out],
        n_subarcs=n_subarcs,
    )


def subarc_defect_residual(
    sub: SubarcSeed,
    states: Sequence[Vec3],
    *,
    ephem: JovianEphemeris,
    cache: JovianRailsCache,
    moons: Sequence[str] = EGGIE_MOONS,
    accuracy: float = 1e-11,
    max_wall_sec: float = 30.0,
) -> NDArray[np.float64]:
    """Sub-arc multiple-shooting residual (generalises :func:`jovian_defect_residual`).

    Layout, concatenated (matches :func:`cyclerfinder.nbody.jovian_stm.subarc_stm_jacobian`):

    1. ``M-1`` sub-arc continuity defects (6 each; velocity x ``_W_VEL``) — every
       consecutive node pair, encounter and interior alike, is a continuity point.
    2. ``n-2`` interior flyby hinges (one per interior ENCOUNTER only; interior
       continuity nodes have none) — read from the carried ``vinf_in/out`` constants.
    3. The 6-component periodicity wrap between the first and last encounter nodes
       (which are node 0 and node ``M-1``), in the home-moon-relative frame.

    With ``sub.n_subarcs == 1`` (``encounter_idx == range(n)``) this is term-for-term
    the same residual as :func:`jovian_defect_residual` on the encounter nodes.
    """
    from cyclerfinder.nbody.jovian import (
        _W_VEL,
        JovianRestrictedNBody,
        _jovian_flyby_hinge_km,
    )
    from cyclerfinder.nbody.shooter import _STATE_DIM

    prop = JovianRestrictedNBody()
    states = [np.asarray(s, dtype=np.float64) for s in states]
    m = len(states)
    res: list[float] = []

    # 1. Sub-arc continuity defects (every consecutive node pair).
    for j in range(m - 1):
        r0 = states[j][:3]
        v0 = states[j][3:]
        arc = prop.propagate(
            r0,
            v0,
            sub.epochs[j],
            sub.epochs[j + 1],
            moons=moons,
            cache=cache,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
        )
        s_next = states[j + 1]
        if arc.converged and np.all(np.isfinite(arc.r_km)) and np.all(np.isfinite(arc.v_km_s)):
            dr = arc.r_km - s_next[:3]
            dv = arc.v_km_s - s_next[3:]
            res.extend(float(x) for x in dr)
            res.extend(float(x) * _W_VEL for x in dv)
        else:
            res.extend([1e9] * _STATE_DIM)

    # 2. Flyby hinges at interior ENCOUNTERS only.
    n_enc = len(sub.encounter_idx)
    for k in range(1, n_enc - 1):
        res.append(_jovian_flyby_hinge_km(sub.vinf_in[k], sub.vinf_out[k], sub.sequence[k]))

    # 3. Periodicity wrap between the first and last encounter nodes.
    i0 = sub.encounter_idx[0]
    i_last = sub.encounter_idx[-1]
    r_home, v_home = ephem.state(sub.sequence[0], sub.epochs[i0])
    r_wrap_pl, v_wrap_pl = ephem.state(sub.sequence[-1], sub.epochs[i_last])
    s0 = states[i0]
    sn = states[i_last]
    rel0_r = s0[:3] - np.asarray(r_home, dtype=np.float64)
    rel0_v = s0[3:] - np.asarray(v_home, dtype=np.float64)
    reln_r = sn[:3] - np.asarray(r_wrap_pl, dtype=np.float64)
    reln_v = sn[3:] - np.asarray(v_wrap_pl, dtype=np.float64)
    res.extend(float(x) for x in (reln_r - rel0_r))
    res.extend(float(x) * _W_VEL for x in (reln_v - rel0_v))

    return np.asarray(res, dtype=np.float64)


def ideal_eggie_shoot(
    guess: EggieGuess,
    *,
    t0_sec: float = 0.0,
    mu: float = MU_JUPITER_KM3_S2,
    accuracy: float = 1e-11,
    max_nfev: int = 20,
    max_wall_sec: float = 20.0,
    periapsis_nodes: bool = True,
    method: str = "lm",
    jacobian: str = "fd",
    n_subarcs: int = 1,
    runlog_path: str | None = None,
) -> ShootResult:
    """Ideal-model multiple-shooting corrector for the EGGIE tour (Stage 2 gate).

    Mirrors :func:`cyclerfinder.nbody.jovian.jovian_shoot`'s least-squares loop but
    injects the ideal circular-coplanar ephemeris + a matching
    :class:`~cyclerfinder.nbody.jovian.JovianRailsCache` (instead of jup365), so the
    corrector runs in the paper's own ideal model. Each residual evaluation is the
    full Jupiter-central multiple-shooting residual
    (:func:`~cyclerfinder.nbody.jovian.jovian_defect_residual`). The result records
    the seed/final defect norms, the continuity-floor convergence flag, per-encounter
    V∞ (compare to Table 4), and the node-impulse correction ΔV. Divergence is honest
    (never raised). Writes an incremental runlog (append+flush per call) if
    ``runlog_path`` is given — the residual is FD-Jacobian-bound, so never a black box.

    ``jacobian`` (#480 Stage 3): ``"fd"`` (default) is the original finite-difference
    path (``least_squares`` builds its own internal FD Jacobian — byte-unchanged) and
    stays the parity oracle. ``"stm"`` supplies the analytic block-bidiagonal Jacobian
    (:func:`cyclerfinder.nbody.jovian_stm.jovian_stm_jacobian`) assembled from per-leg
    analytic state-transition matrices — ONE co-integrated variational propagation per
    leg instead of the ``6*n_nodes+1`` FD re-propagations, the lever built to break the
    Stage-2 FD-noise plateau (``docs/notes/2026-06-29-480-eggie-stage2-nbody-verdict``).

    ``n_subarcs`` (#480 Stage 4): ``1`` (default) keeps the byte-unchanged
    one-node-per-leg path (the existing :func:`jovian_defect_residual` /
    :func:`cyclerfinder.nbody.jovian_stm.jovian_stm_jacobian`). ``>1`` inserts
    ``n_subarcs - 1`` interior continuity nodes per leg (:func:`build_subarc_seed`) and
    closes the expanded :func:`subarc_defect_residual` with
    :func:`cyclerfinder.nbody.jovian_stm.subarc_stm_jacobian` — splitting each defect
    into shorter sub-arcs so the perijove sensitivity is distributed (the Stage-3
    verdict's highest-leverage next lever). Each iter co-integrates
    ``(n-1) * n_subarcs`` STMs.
    """
    from scipy.optimize import least_squares

    from cyclerfinder.nbody.correction_dv import node_impulse_correction_dv
    from cyclerfinder.nbody.jovian import (
        JovianRailsCache,
        _jovian_flyby_hinge_km,
        jovian_defect_residual,
    )
    from cyclerfinder.nbody.shooter import (
        _POS_CONTINUITY_KM,
        _STATE_DIM,
        _VEL_CONTINUITY_KMS,
        ShootResult,
        _seed_with_states,
        _states_to_x,
        _x_to_states,
    )

    if n_subarcs < 1:
        raise ValueError(f"n_subarcs must be >= 1, got {n_subarcs}")

    w_vel = 1.0e3  # velocity-residual weight (mirror jovian_defect_residual _W_VEL)

    ephem = ideal_ephemeris_from_guess(guess, t0_sec=t0_sec, mu=mu)
    if periapsis_nodes:
        seed = build_eggie_periapsis_seed(guess, t0_sec=t0_sec, mu=mu)
    else:
        seed = build_eggie_shooting_seed(guess, t0_sec=t0_sec, mu=mu)
    moons = EGGIE_MOONS
    n = len(seed.sequence)

    # Sub-arc node set (n_subarcs == 1 -> the encounter nodes themselves).
    if n_subarcs == 1:
        sub = None
        states0 = [np.asarray(s, dtype=np.float64) for s in seed.node_states]
        epochs = list(seed.epochs)
        encounter_idx = tuple(range(n))
    else:
        sub = build_subarc_seed(seed, n_subarcs, ephem=ephem, moons=moons)
        states0 = list(sub.node_states)
        epochs = list(sub.epochs)
        encounter_idx = sub.encounter_idx
    m = len(states0)
    # IdealJovianEphemeris is duck-typed (.state) to the JovianEphemeris API.
    cache = JovianRailsCache(moons, ephem, min(epochs), max(epochs))  # type: ignore[arg-type]

    log = open(runlog_path, "a", buffering=1) if runlog_path else None  # noqa: SIM115
    nfev = {"count": 0}

    def residual_of_x(x: NDArray[np.float64]) -> NDArray[np.float64]:
        states = _x_to_states(x, m)
        if sub is None:
            trial = _seed_with_states(seed, states)
            r = jovian_defect_residual(
                trial,
                ephem=ephem,  # type: ignore[arg-type]  # duck-typed ideal ephemeris
                cache=cache,
                moons=moons,
                accuracy=accuracy,
                max_wall_sec=max_wall_sec,
            )
        else:
            r = subarc_defect_residual(
                sub,
                states,
                ephem=ephem,  # type: ignore[arg-type]  # duck-typed ideal ephemeris
                cache=cache,
                moons=moons,
                accuracy=accuracy,
                max_wall_sec=max_wall_sec,
            )
        if log is not None:
            nfev["count"] += 1
            log.write(
                f"{time.strftime('%Y-%m-%dT%H:%M:%S')} nfev={nfev['count']} "
                f"defect={float(np.linalg.norm(r)):.6e}\n"
            )
        return r

    x0 = _states_to_x(states0)
    seed_res = residual_of_x(x0)
    seed_defect_norm = float(np.linalg.norm(seed_res))

    ls_kwargs = {"x_scale": "jac"} if method != "lm" else {}
    if jacobian == "stm":
        from cyclerfinder.nbody.jovian_stm import jovian_stm_jacobian, subarc_stm_jacobian

        def jac_of_x(x: NDArray[np.float64]) -> NDArray[np.float64]:
            t_eval = time.monotonic()
            if sub is None:
                jac = jovian_stm_jacobian(seed, x, ephem=ephem, moons=moons)
            else:
                jac = subarc_stm_jacobian(sub, x, ephem=ephem, moons=moons)
            if log is not None:
                log.write(
                    f"{time.strftime('%Y-%m-%dT%H:%M:%S')} jac_stm "
                    f"dt={time.monotonic() - t_eval:.2f}s\n"
                )
            return jac

        ls_kwargs["jac"] = jac_of_x  # type: ignore[assignment]
    sol = least_squares(residual_of_x, x0, method=method, max_nfev=max_nfev, **ls_kwargs)  # type: ignore[arg-type]
    corrected_states = _x_to_states(np.asarray(sol.x, dtype=np.float64), m)
    final_res = residual_of_x(np.asarray(sol.x, dtype=np.float64))
    defect_norm = float(np.linalg.norm(final_res))
    if log is not None:
        log.close()

    n_leg_defects = (m - 1) * _STATE_DIM
    leg_block = final_res[:n_leg_defects]
    floor_vec = ([_POS_CONTINUITY_KM] * 3 + [_VEL_CONTINUITY_KMS * w_vel] * 3) * (m - 1)
    continuity_floor = float(np.linalg.norm(floor_vec))
    converged = bool(sol.success) and float(np.linalg.norm(leg_block)) < continuity_floor

    vinf_mag: list[float] = []
    for enc, body in zip(encounter_idx, seed.sequence, strict=True):
        _, v_m = ephem.state(body, epochs[enc])
        state = corrected_states[enc]
        vinf_mag.append(float(np.linalg.norm(np.asarray(state[3:], dtype=np.float64) - v_m)))

    seed_nodes = {f"b{i}": np.asarray(states0[enc][3:]) for i, enc in enumerate(encounter_idx)}
    corr_nodes = {
        f"b{i}": np.asarray(corrected_states[enc][3:]) for i, enc in enumerate(encounter_idx)
    }
    cdv = node_impulse_correction_dv(seed_nodes, corr_nodes)

    bend_feasible = all(
        _jovian_flyby_hinge_km(seed.vinf_in[i], seed.vinf_out[i], seed.sequence[i]) <= 0.0
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
    "EGGIE_MOONS",
    "IdealJovianEphemeris",
    "SubarcSeed",
    "build_eggie_periapsis_seed",
    "build_eggie_shooting_seed",
    "build_subarc_seed",
    "ideal_eggie_shoot",
    "ideal_ephemeris_from_guess",
    "subarc_defect_residual",
]
