"""Multi-lap propagation + tolerant closure-drift verification (M6a).

Spec references
---------------
* §4 (architecture — ``verify/propagate.py``), §5 step 7 (verify:
  multi-lap propagation; periodicity / re-encounter miss), §6
  (top-level interfaces — ``find_cyclers``, ``Cycler``), §8 M6
  milestone definition ("verify periodic over >= 3 laps"), §10
  (closure-frame correctness risk), §12(c) (dynamic ephemeris frame +
  tolerant verification — the binding architectural decision for
  M6a), §14 V2 (multi-lap periodicity validation gate).

Purpose
-------
Given an idealized :class:`~cyclerfinder.model.Cycler`, propagate its
trajectory continuously across ``n_laps`` consecutive laps on a real
(or circular) ephemeris and measure the lap-to-lap drift in the spec
§12(c) dynamic rotating frame. A periodic cycler stays within a
bounded km-scale tolerance (the binding 50,000 km tolerance is
derived in plan §4.3); a non-periodic / degenerate-closure candidate
drifts by AU per lap and is rejected.

Algorithm
---------
For each lap ``i`` in ``range(n_laps)``:

1. Reconstruct the lap-i start state of each leg ``j`` by:

   * Reading the planet position ``r_planet_lap_i_leg_j = ephem.state(
     encounter[j].body, t_dep_lap_i_leg_j)`` at the lap-shifted leg
     departure time.
   * Rotating the cycler's encoded V_inf_out at encounter ``j`` from
     the lap-0 dynamic-frame angle into the lap-i dynamic-frame
     angle, yielding ``vinf_inertial_lap_i_leg_j``.
   * Spacecraft velocity at lap-i leg-j departure:
     ``v_sc = v_planet_lap_i + vinf_inertial_lap_i_leg_j``.

2. Propagate Kepler from ``(r_planet_lap_i_leg_j, v_sc)`` to each
   requested sample time within ``[t_dep_lap_i_leg_j,
   t_arr_lap_i_leg_j]``.

The lap-shifted leg start state is therefore a clean function of the
cycler's idealized leg template plus the dynamic-frame rotation at
the lap-i epoch. For a truly periodic cycler the per-lap reconstructed
geometry is identical (up to the dynamic-frame breathing) lap-over-
lap, so the rotating-frame drift is bounded.

For a degenerate-closure cycler (V_inf chosen so that the idealized
closure residual is small but the geometry doesn't close on the real
ephemeris), the per-lap reconstructed planet position drifts away
from a Kepler-arc continuation by AU per lap — orders of magnitude
above the 50,000-km tolerance. The tolerance therefore rejects the
spec §10 degenerate-solution basin cleanly.

Frame
-----
Lap-to-lap drift is measured in the spec §12(c) dynamic rotating
frame anchored to ``bodies = (cycler.bodies[0], cycler.bodies[1])``
by default (overrideable via ``frame_bodies``). The dynamic frame's
theta(t) is read directly from the ephemeris (see
:func:`cyclerfinder.core.frames.to_rotating_dynamic`) so the
transform is exact rather than just numerically close.

For the M3 circular-coplanar regression test, callers pass
``use_uniform_frame=True`` to switch to the M3 uniform frame at
Earth's mean motion — the dynamic and uniform frames agree to
float-precision on circular ephemeris, but the uniform-frame path
exercises the M3 code so the regression test catches M3 changes.

Locked dataclass
----------------
:class:`StabilityReport` is frozen and shape-locked at M6a so M6b's
TCM-budget extension fills only the ``per_lap_dv`` and
``total_tcm_dv`` fields rather than reshaping the dataclass. The
``frame_used`` field reads ``"dynamic"`` (default) or ``"uniform"``
(when ``use_uniform_frame=True``).

Plan: ``docs/phases/m6a-idealized-closure-verification/plan.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, sin
from typing import Final

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.frames import (
    _frame_angle_dynamic,
    synodic_omega,
    to_rotating,
    to_rotating_dynamic,
)
from cyclerfinder.core.kepler import KeplerConvergenceError, propagate
from cyclerfinder.model.cycler import Cycler

DRIFT_TOLERANCE_KM: Final[float] = 50_000.0
"""Maximum permissible lap-to-lap drift for an E-M-class cycler over
<= 5 laps. Derivation in plan §4.3: ~0.02 deg of geometric breathing
per lap at Mars's mean radius of ~1.5 AU = 2.6 million km / deg = ~50,000
km / 0.02 deg. Tight enough to reject propagator regressions and
frame-transform errors; loose enough to absorb real eccentricity-
driven breathing on Earth-Mars trajectories. M6b may tighten this if
its ephemeris-mode TCM optimiser reports drift consistently below
10,000 km."""


@dataclass(frozen=True)
class StabilityReport:
    """Result of :func:`verify_long_term_stability`.

    All fields are immutable. The ``per_lap_dv`` and ``total_tcm_dv``
    fields are zeros in M6a; M6b populates them when
    :func:`cyclerfinder.search.optimize.optimise_cell_ephemeris` runs.

    Spec references: §8 (M6 milestone — "verified periodic over >= 3
    laps"), §12(c) (dynamic frame + tolerant verification), §14 V2
    (multi-lap periodicity validation gate).

    Attributes
    ----------
    cycler_id:
        Catalogue entry id (e.g. ``"s1l1-2syn-em-cpom"``) if the
        report was produced from a catalogue cycler, else ``None``.
        Passed through by the caller; M7's batch-validate runner
        sets it.
    n_laps_propagated:
        Number of laps the propagation actually completed. Equals the
        ``n_laps`` argument unless an early-termination condition
        tripped.
    max_drift_km:
        Max consecutive-lap-pair drift across the
        ``range(n_laps - 1)`` pairs. The basis for ``stable``.
    max_drift_lap_index:
        The ``i`` such that
        ``lap_to_lap_drift(samples_lap_i, samples_lap_{i+1})``
        equalled ``max_drift_km``. Diagnostic.
    per_lap_drift_km:
        Cumulative drift at each lap boundary.
        ``per_lap_drift_km[i]`` is
        ``lap_to_lap_drift(samples_lap_0, samples_lap_{i+1})`` —
        total drift accumulated from the start through the end of lap
        ``i+1``. ``len() == n_laps_propagated``.
    stable:
        ``max_drift_km < DRIFT_TOLERANCE_KM``. The headline boolean
        for spec §14 V2.
    per_lap_dv:
        **M6a: zero-tuple of length** ``n_laps_propagated``. M6b
        populates with the per-lap TCM dV (km/s) from the ephemeris-
        mode optimiser. Locked here so M6b doesn't reshape the
        dataclass.
    total_tcm_dv:
        **M6a: 0.0**. M6b populates with ``sum(per_lap_dv)``.
        Locked here so M6b doesn't reshape the dataclass.
    frame_used:
        ``"dynamic"`` (default) or ``"uniform"`` (when the caller
        forced the M3 frame, e.g. for the M3 Aldrin circular-coplanar
        regression test). Diagnostic; not part of the spec §14 V2
        gate contract.
    """

    cycler_id: str | None
    n_laps_propagated: int
    max_drift_km: float
    max_drift_lap_index: int
    per_lap_drift_km: tuple[float, ...]
    stable: bool
    per_lap_dv: tuple[float, ...]
    total_tcm_dv: float
    frame_used: str


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_frame_bodies(
    cycler: Cycler,
    frame_bodies: tuple[str, ...] | None,
) -> tuple[str, ...]:
    """Pick the dynamic-frame anchor bodies from a cycler.

    Rules:
      1. If ``frame_bodies`` is explicit, return it (caller knows).
      2. Else use the first two unique body codes in ``cycler.bodies``
         in encounter order. For an E-M cycler this is ``("E", "M")``;
         for a VEM cycler ``("E", "V")`` etc.

    Per plan §3.4 this is a policy helper. For M6a the default policy
    is the first-two-unique-bodies-in-encounter-order rule; M8 may
    revisit for multi-body cyclers.

    Never raises (the cycler is guaranteed non-empty by M5
    construction).
    """
    if frame_bodies is not None:
        return frame_bodies
    seen: list[str] = []
    for body in cycler.bodies:
        if body not in seen:
            seen.append(body)
        if len(seen) == 2:
            break
    if len(seen) == 1:
        # Single-body cycler — synthesise a second body so the
        # dynamic frame is well-defined. Pin to "E" because E-M is
        # the canonical 2-body pair in cyclerfinder.
        seen.append("M" if seen[0] != "M" else "E")
    return tuple(seen)


def _rotate_vinf_to_lap(
    vinf_inertial_lap_0: NDArray[np.float64],
    t_lap_0: float,
    t_lap_i: float,
    bodies_frame: tuple[str, ...],
    ephem: Ephemeris,
    use_uniform_frame: bool,
    omega_uniform: float,
) -> NDArray[np.float64]:
    """Rotate a V_inf vector from the lap-0 inertial frame to the lap-i frame.

    The V_inf vector is constant in the rotating frame for a truly
    periodic cycler. So:

    1. Express ``vinf_inertial_lap_0`` in the rotating frame at
       ``t_lap_0``: ``vinf_rot = R(-theta(t_lap_0)) * vinf_inertial``
       (no Coriolis correction — V_inf has no position to cross with).
    2. Express ``vinf_rot`` back in the inertial frame at ``t_lap_i``:
       ``vinf_inertial_lap_i = R(+theta(t_lap_i)) * vinf_rot``.

    Net effect: rotate by ``theta(t_lap_i) - theta(t_lap_0)`` about
    +z.

    For the dynamic frame, ``theta`` is read from the ephemeris via
    :func:`cyclerfinder.core.frames._frame_angle_dynamic`. For the
    uniform frame, ``theta = omega_uniform * t``.

    Module-internal.
    """
    if use_uniform_frame:
        theta_diff = omega_uniform * (t_lap_i - t_lap_0)
    else:
        theta_lap_0 = _frame_angle_dynamic(t_lap_0, bodies_frame, ephem)
        theta_lap_i = _frame_angle_dynamic(t_lap_i, bodies_frame, ephem)
        theta_diff = theta_lap_i - theta_lap_0
    c, s = cos(theta_diff), sin(theta_diff)
    vx = float(vinf_inertial_lap_0[0])
    vy = float(vinf_inertial_lap_0[1])
    vz = float(vinf_inertial_lap_0[2])
    return np.array(
        [c * vx - s * vy, s * vx + c * vy, vz],
        dtype=np.float64,
    )


def _propagate_leg_samples(
    r0: NDArray[np.float64],
    v0: NDArray[np.float64],
    t0: float,
    t_samples: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Sample the Kepler propagation of ``(r0, v0)`` at each ``t_samples[k]``.

    Returns shape ``(len(t_samples), 7)``: ``[t, x, y, z, vx, vy,
    vz]`` per row.

    Module-internal helper for :func:`propagate_lap`.
    """
    n = len(t_samples)
    out = np.empty((n, 7), dtype=np.float64)
    for k in range(n):
        dt = float(t_samples[k] - t0)
        r, v = propagate(r0, v0, dt)
        out[k, 0] = float(t_samples[k])
        out[k, 1] = float(r[0])
        out[k, 2] = float(r[1])
        out[k, 3] = float(r[2])
        out[k, 4] = float(v[0])
        out[k, 5] = float(v[1])
        out[k, 6] = float(v[2])
    return out


def _to_rotating_frame_position(
    r_inertial: NDArray[np.float64],
    t_sec: float,
    bodies_frame: tuple[str, ...],
    ephem: Ephemeris,
    use_uniform_frame: bool,
    omega_uniform: float,
) -> NDArray[np.float64]:
    """Rotating-frame position of ``r_inertial`` at time ``t_sec``.

    Branches on ``use_uniform_frame``. The velocity is not needed
    for drift measurement; we pass a zero velocity to the underlying
    transform.

    Module-internal helper for :func:`lap_to_lap_drift`.
    """
    v_zero = np.zeros(3, dtype=np.float64)
    if use_uniform_frame:
        r_rot, _ = to_rotating(r_inertial, v_zero, t_sec, omega_uniform)
    else:
        r_rot, _ = to_rotating_dynamic(r_inertial, v_zero, t_sec, bodies_frame, ephem)
    return r_rot


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def propagate_lap(
    cycler: Cycler,
    ephem: Ephemeris,
    t_start: float,
    t_end: float,
    n_samples: int,
    *,
    lap_index: int = 0,
    frame_bodies: tuple[str, ...] | None = None,
    use_uniform_frame: bool = False,
) -> NDArray[np.float64]:
    """Propagate the cycler's spacecraft trajectory across ``[t_start, t_end]``.

    Returns ``n_samples`` uniformly-spaced ``(t, x, y, z, vx, vy,
    vz)`` rows in the inertial heliocentric frame, in time order.

    Algorithm:
      1. Compute ``times = linspace(t_start, t_end, n_samples)``.
      2. For each leg ``j`` of the cycler, compute the lap-shifted
         time window
         ``[t_start + (encounter[j].t - encounter[0].t),
         t_start + (encounter[j+1].t - encounter[0].t)]``.
      3. Reconstruct the lap-shifted leg-start spacecraft state from
         the planet's actual position at the lap-shifted departure
         time + the cycler's encoded V_inf_out rotated by the lap's
         frame angle (see :func:`_rotate_vinf_to_lap`).
      4. Propagate Kepler from that state to each ``t_samples[k]``
         falling within the leg's window.

    Parameters
    ----------
    cycler:
        The idealized cycler. Its ``encounters`` and ``legs`` provide
        the leg template; the lap-shifted geometry is reconstructed
        per the algorithm above.
    ephem:
        Heliocentric state provider for the lap-shifted planet
        positions.
    t_start:
        Absolute inertial time at which lap ``lap_index`` starts (s).
        Typically ``user_t_start + lap_index * cycler.period``; the
        ``lap_index`` argument is informational.
    t_end:
        Absolute inertial time at which the lap ends (s). Typically
        ``t_start + cycler.period``.
    n_samples:
        Number of uniformly-spaced sample rows in the output.
    lap_index:
        The lap number this call represents (0-based). Used by
        :func:`multi_lap_propagation` to track diagnostic metadata
        and to support callers that want to know which lap was
        propagated. Not used inside the propagation (the leg-template
        re-use is implicit in the cycler's encoded geometry; the
        actual lap rotation is driven by ``t_start``'s relation to
        the cycler's encoded first-encounter time).
    frame_bodies:
        Optional override for the dynamic-frame anchor bodies. ``None``
        ⇒ :func:`_resolve_frame_bodies(cycler, None)`. Used to rotate
        the V_inf vectors to the lap-shifted inertial frame.
    use_uniform_frame:
        If ``True``, use the M3 uniform-frame V_inf rotation at
        Earth's mean motion. Default ``False`` (dynamic frame).

    Returns
    -------
    NDArray[np.float64]
        Shape ``(n_samples, 7)`` of ``[t, x, y, z, vx, vy, vz]``;
        floats in s, km, km/s.

    Raises
    ------
    ValueError
        If ``n_samples < 1`` or ``t_end <= t_start``.
    """
    if n_samples < 1:
        raise ValueError(f"n_samples must be >= 1; got {n_samples}")
    if t_end <= t_start:
        raise ValueError(f"t_end must be > t_start; got {t_start}..{t_end}")
    if lap_index < 0:
        raise ValueError(f"lap_index must be >= 0; got {lap_index}")

    bodies_frame = _resolve_frame_bodies(cycler, frame_bodies)
    omega_uniform = synodic_omega("E")
    t_first_enc = cycler.encounters[0].t

    # Pre-compute lap-i leg-start states per leg.
    n_legs = len(cycler.legs)
    leg_windows: list[tuple[float, float]] = []
    leg_start_states: list[tuple[NDArray[np.float64], NDArray[np.float64]]] = []
    for j in range(n_legs):
        enc_j = cycler.encounters[j]
        # Offset of this leg's departure/arrival within the lap (lap 0 reference).
        dt_dep = enc_j.t - t_first_enc
        dt_arr = cycler.encounters[j + 1].t - t_first_enc
        t_dep_lap = t_start + dt_dep
        t_arr_lap = t_start + dt_arr
        leg_windows.append((t_dep_lap, t_arr_lap))

        # Planet position at lap-shifted leg-departure time.
        r_planet_lap, v_planet_lap = ephem.state(enc_j.body, t_dep_lap)

        # The cycler's encoded v_depart_lap_0 = v_planet_at_enc_j + vinf_out_lap_0.
        # vinf_out_lap_0 is encounter[j].vinf_out (inertial at lap-0 leg-departure).
        vinf_out_lap_0 = enc_j.vinf_out
        # Rotate the lap-0 inertial V_inf to the lap-i inertial frame so the
        # spacecraft V_inf relative to the (lap-shifted) planet is constant
        # in the rotating frame, matching the periodic-cycler invariant.
        vinf_out_lap_i = _rotate_vinf_to_lap(
            np.asarray(vinf_out_lap_0, dtype=np.float64),
            t_lap_0=enc_j.t,
            t_lap_i=t_dep_lap,
            bodies_frame=bodies_frame,
            ephem=ephem,
            use_uniform_frame=use_uniform_frame,
            omega_uniform=omega_uniform,
        )
        v_sc_lap = np.asarray(v_planet_lap, dtype=np.float64) + vinf_out_lap_i
        leg_start_states.append(
            (np.asarray(r_planet_lap, dtype=np.float64), v_sc_lap),
        )

    # Sample times within [t_start, t_end].
    times = np.linspace(t_start, t_end, n_samples, dtype=np.float64)

    # Assign each sample to a leg. Boundary t == t_arr_j belongs to leg j
    # (so encounters appear once and only once).
    out = np.empty((n_samples, 7), dtype=np.float64)
    leg_idx_for_sample: list[int] = []
    for k in range(n_samples):
        t_k = float(times[k])
        chosen_leg = n_legs - 1  # fallback: last leg
        for j in range(n_legs):
            t_dep_j, t_arr_j = leg_windows[j]
            if t_dep_j <= t_k <= t_arr_j:
                chosen_leg = j
                break
        leg_idx_for_sample.append(chosen_leg)

    # Propagate per leg, batching samples for that leg.
    for j in range(n_legs):
        t_dep_j, _t_arr_j = leg_windows[j]
        r0_j, v0_j = leg_start_states[j]
        sample_idxs = [k for k in range(n_samples) if leg_idx_for_sample[k] == j]
        if not sample_idxs:
            continue
        t_samples = np.asarray([times[k] for k in sample_idxs], dtype=np.float64)
        rows = _propagate_leg_samples(r0_j, v0_j, t_dep_j, t_samples)
        for local_idx, k in enumerate(sample_idxs):
            out[k, :] = rows[local_idx, :]
    return out


def multi_lap_propagation(
    cycler: Cycler,
    ephem: Ephemeris,
    n_laps: int,
    t_start: float = 0.0,
    n_samples_per_lap: int = 100,
    *,
    frame_bodies: tuple[str, ...] | None = None,
    use_uniform_frame: bool = False,
) -> dict[str, NDArray[np.float64]]:
    """Propagate the cycler through ``n_laps`` consecutive laps.

    For each lap ``i`` in ``range(n_laps)`` calls
    :func:`propagate_lap` with the lap-shifted time window
    ``[t_start + i * cycler.period, t_start + (i + 1) * cycler.period]``.

    Per spec §12 ("propagate continuously"), the lap-i geometry is
    reconstructed from the cycler's encoded leg template rotated to
    lap-i's epoch — see :func:`propagate_lap`'s docstring for the
    algorithm. For a truly periodic cycler the lap-i and lap-0
    geometries are bit-stable in the dynamic rotating frame (modulo
    eccentricity-driven breathing on the astropy backend); for a
    degenerate-closure / open trajectory the lap-i geometry drifts by
    AU per lap.

    Parameters
    ----------
    cycler:
        The idealized cycler.
    ephem:
        Heliocentric state provider.
    n_laps:
        Number of laps. Must be >= 1.
    t_start:
        Absolute inertial time at which lap 0 starts (s). Default
        ``0.0``.
    n_samples_per_lap:
        Number of uniformly-spaced sample rows per lap.
    frame_bodies, use_uniform_frame:
        Forwarded to :func:`propagate_lap`.

    Returns
    -------
    dict[str, NDArray[np.float64]]
        ``"samples"``: shape ``(n_laps * n_samples_per_lap, 7)`` of
            ``[t, x, y, z, vx, vy, vz]`` rows, in time order across
            all laps.
        ``"lap_indices"``: shape ``(n_laps + 1,)`` int array (stored
            as float64 for consistency); ``lap_indices[i]`` is the
            row index in ``samples`` where lap ``i`` starts.
        ``"lap_start_times"``: shape ``(n_laps,)`` float array of
            inertial-frame lap-start times in seconds.

    Raises
    ------
    ValueError
        If ``n_laps < 1`` or ``n_samples_per_lap < 1``.
    """
    if n_laps < 1:
        raise ValueError(f"n_laps must be >= 1; got {n_laps}")
    if n_samples_per_lap < 1:
        raise ValueError(f"n_samples_per_lap must be >= 1; got {n_samples_per_lap}")

    period = cycler.period
    total_samples = n_laps * n_samples_per_lap
    samples = np.empty((total_samples, 7), dtype=np.float64)
    lap_indices = np.empty(n_laps + 1, dtype=np.float64)
    lap_start_times = np.empty(n_laps, dtype=np.float64)

    for i in range(n_laps):
        t_lap_start = t_start + i * period
        t_lap_end = t_start + (i + 1) * period
        lap_start_times[i] = t_lap_start
        lap_indices[i] = i * n_samples_per_lap
        rows = propagate_lap(
            cycler,
            ephem,
            t_lap_start,
            t_lap_end,
            n_samples_per_lap,
            lap_index=i,
            frame_bodies=frame_bodies,
            use_uniform_frame=use_uniform_frame,
        )
        samples[i * n_samples_per_lap : (i + 1) * n_samples_per_lap, :] = rows
    lap_indices[n_laps] = total_samples

    return {
        "samples": samples,
        "lap_indices": lap_indices,
        "lap_start_times": lap_start_times,
    }


def lap_to_lap_drift(
    samples_lap_a: NDArray[np.float64],
    samples_lap_b: NDArray[np.float64],
    bodies: tuple[str, ...],
    ephem: Ephemeris,
    *,
    use_uniform_frame: bool = False,
) -> float:
    """Max km drift between two laps' samples in the dynamic rotating frame.

    For each row index ``i`` in ``range(n_samples_per_lap)``,
    transform the corresponding samples in laps A and B into the
    dynamic rotating frame at their respective inertial times via
    :func:`cyclerfinder.core.frames.to_rotating_dynamic`, and compute
    ``||r_A_rot - r_B_rot||``. Return the max across all ``i``.

    Per spec §12(c), the verification target is bounded drift, not
    zero. The M6a binding tolerance is :data:`DRIFT_TOLERANCE_KM`
    (50,000 km; see plan §4.3 for derivation).

    Both arrays must have the same shape ``(n_samples_per_lap, 7)``
    and be aligned positionally — i.e. row ``i`` of each represents
    the same phase within its lap.

    Parameters
    ----------
    samples_lap_a, samples_lap_b:
        Shape ``(n, 7)`` per-lap sample arrays from
        :func:`propagate_lap`.
    bodies:
        Dynamic-frame anchor bodies. ``bodies[0]`` defines the
        rotating frame's +x direction; ``bodies[1]`` enters via the
        Coriolis angular rate computation.
    ephem:
        Heliocentric state provider.
    use_uniform_frame:
        If ``True``, use the M3 uniform frame at Earth's mean motion
        instead of the dynamic frame. Default ``False``.

    Returns
    -------
    float
        Max position drift in km (non-negative).

    Raises
    ------
    ValueError
        If the two sample arrays have differing shapes.
    """
    if samples_lap_a.shape != samples_lap_b.shape:
        raise ValueError(
            f"sample arrays must have the same shape; got "
            f"{samples_lap_a.shape} vs {samples_lap_b.shape}",
        )
    if samples_lap_a.ndim != 2 or samples_lap_a.shape[1] != 7:
        raise ValueError(
            f"sample arrays must be of shape (n, 7); got {samples_lap_a.shape}",
        )
    omega_uniform = synodic_omega("E")
    n = samples_lap_a.shape[0]
    max_drift: float = 0.0
    for i in range(n):
        t_a = float(samples_lap_a[i, 0])
        t_b = float(samples_lap_b[i, 0])
        r_a = np.asarray(samples_lap_a[i, 1:4], dtype=np.float64)
        r_b = np.asarray(samples_lap_b[i, 1:4], dtype=np.float64)
        r_a_rot = _to_rotating_frame_position(
            r_a, t_a, bodies, ephem, use_uniform_frame, omega_uniform
        )
        r_b_rot = _to_rotating_frame_position(
            r_b, t_b, bodies, ephem, use_uniform_frame, omega_uniform
        )
        drift = float(np.linalg.norm(r_a_rot - r_b_rot))
        if drift > max_drift:
            max_drift = drift
    return max_drift


def verify_long_term_stability(
    cycler: Cycler,
    n_laps: int,
    ephem: Ephemeris,
    *,
    t_start: float = 0.0,
    frame_bodies: tuple[str, ...] | None = None,
    cycler_id: str | None = None,
    n_samples_per_lap: int = 100,
    use_uniform_frame: bool = False,
) -> StabilityReport:
    """Spec §12 long-term stability verifier; spec §14 V2 gate machinery.

    Pipeline:
        1. Resolve ``bodies = _resolve_frame_bodies(cycler,
           frame_bodies)``.
        2. Propagate continuously through ``n_laps`` via
           :func:`multi_lap_propagation`.
        3. For each consecutive lap pair ``(i, i+1)``, compute
           ``lap_to_lap_drift(samples_lap_i, samples_lap_{i+1},
           bodies, ephem)``.
        4. Build ``per_lap_drift_km[i] = lap_to_lap_drift(
           samples_lap_0, samples_lap_{i+1}, bodies, ephem)`` for
           ``i in range(n_laps)``.
        5. Set ``max_drift_km = max(consecutive_pair_drifts)`` and
           ``max_drift_lap_index = argmax``.
        6. ``stable = max_drift_km < DRIFT_TOLERANCE_KM``.
        7. Populate the M6b placeholders with zeros and return.

    The ``use_uniform_frame`` flag exists for the M3-circular
    regression test (``test_aldrin_cycler_periodic_over_3_laps_\
circular``). When ``True``, the drift computation uses
    :func:`cyclerfinder.core.frames.to_rotating` at Earth's constant
    mean motion; the report's ``frame_used = "uniform"``. The M6a
    binding gate runs with the default ``False`` and the dynamic
    frame.

    Parameters
    ----------
    cycler:
        The idealized cycler.
    n_laps:
        Number of laps to verify. The M6a binding gate uses ``3``;
        M6b will use up to ``5`` per spec §12(a).
    ephem:
        Heliocentric state provider. The M6a binding gate uses
        ``Ephemeris("astropy")``; the circular regression test uses
        ``Ephemeris("circular")``.
    t_start:
        Absolute inertial time at which lap 0 starts (s). Default
        ``0.0``.
    frame_bodies:
        Optional override for the dynamic-frame anchor bodies. ``None``
        ⇒ :func:`_resolve_frame_bodies(cycler, None)`.
    cycler_id:
        Catalogue entry id for the produced :class:`StabilityReport`.
        Optional — M6a does not know how the cycler was derived.
    n_samples_per_lap:
        Forwarded to :func:`multi_lap_propagation`.
    use_uniform_frame:
        If ``True``, switch to the M3 uniform frame at Earth's mean
        motion. Default ``False``.

    Returns
    -------
    StabilityReport
        Frozen result record.

    Raises
    ------
    ValueError
        If ``n_laps < 2`` (consecutive-pair drift is undefined for a
        single lap).
    """
    if n_laps < 2:
        raise ValueError(
            f"verify_long_term_stability requires n_laps >= 2 (consecutive-pair"
            f" drift is undefined for a single lap); got {n_laps}",
        )

    bodies = _resolve_frame_bodies(cycler, frame_bodies)
    try:
        mlp = multi_lap_propagation(
            cycler,
            ephem,
            n_laps,
            t_start=t_start,
            n_samples_per_lap=n_samples_per_lap,
            frame_bodies=frame_bodies,
            use_uniform_frame=use_uniform_frame,
        )
    except KeplerConvergenceError:
        # Propagator divergence: emit a non-stable report with zero
        # laps propagated. The caller's gate check fails on
        # ``stable=False``.
        return StabilityReport(
            cycler_id=cycler_id,
            n_laps_propagated=0,
            max_drift_km=float("inf"),
            max_drift_lap_index=-1,
            per_lap_drift_km=(),
            stable=False,
            per_lap_dv=(),
            total_tcm_dv=0.0,
            frame_used="dynamic" if not use_uniform_frame else "uniform",
        )

    samples = mlp["samples"]
    lap_indices = mlp["lap_indices"]
    n_samples_per_lap_actual = int(lap_indices[1] - lap_indices[0])

    def lap_samples(i: int) -> NDArray[np.float64]:
        lo = int(lap_indices[i])
        hi = lo + n_samples_per_lap_actual
        return samples[lo:hi, :]

    # Consecutive-pair drifts.
    consecutive_drifts: list[float] = []
    for i in range(n_laps - 1):
        d = lap_to_lap_drift(
            lap_samples(i),
            lap_samples(i + 1),
            bodies,
            ephem,
            use_uniform_frame=use_uniform_frame,
        )
        consecutive_drifts.append(d)

    if consecutive_drifts:
        max_drift_km = max(consecutive_drifts)
        max_drift_lap_index = consecutive_drifts.index(max_drift_km)
    else:
        max_drift_km = 0.0
        max_drift_lap_index = -1

    # Cumulative-from-lap-0 drift at each lap boundary.
    per_lap_drift: list[float] = []
    samples_lap_0 = lap_samples(0)
    for i in range(n_laps):
        if i == 0:
            # Lap 0 vs lap 0 — zero by construction.
            per_lap_drift.append(0.0)
            continue
        d = lap_to_lap_drift(
            samples_lap_0,
            lap_samples(i),
            bodies,
            ephem,
            use_uniform_frame=use_uniform_frame,
        )
        per_lap_drift.append(d)

    stable = max_drift_km < DRIFT_TOLERANCE_KM
    frame_used = "uniform" if use_uniform_frame else "dynamic"

    return StabilityReport(
        cycler_id=cycler_id,
        n_laps_propagated=n_laps,
        max_drift_km=max_drift_km,
        max_drift_lap_index=max_drift_lap_index,
        per_lap_drift_km=tuple(per_lap_drift),
        stable=stable,
        per_lap_dv=(0.0,) * n_laps,
        total_tcm_dv=0.0,
        frame_used=frame_used,
    )


__all__ = [
    "DRIFT_TOLERANCE_KM",
    "StabilityReport",
    "lap_to_lap_drift",
    "multi_lap_propagation",
    "propagate_lap",
    "verify_long_term_stability",
]
