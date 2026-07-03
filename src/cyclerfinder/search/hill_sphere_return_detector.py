"""Repeated-Hill-sphere-encounter detector for transient-drift quasi_cyclers (#535).

Implements the admission criterion settled in writing before any sweep code
was built (see ``docs/notes/2026-07-03-535-quasi-cycler-transient-drift-
admission-criterion.md`` -- the single source of truth this module must
match exactly):

1. **Encounter**: distance to the target body < ``r_hill``.
2. **Return**: one MAXIMAL CONTINUOUS Hill-sphere-residency interval (not
   sub-counting periapsis wiggles inside one episode).
3. **Distinctness**: two returns count as separate only if the gap between
   them (time OUTSIDE the Hill sphere) is >= ``min_separation``.
4. **Admission window**: a sliding window of length in
   ``[window_lo, window_hi]`` containing between ``n_returns_lo`` and
   ``n_returns_hi`` (inclusive) distinct returns is ADMISSIBLE; report the
   window bounds and the actual return epochs used, never a silently
   cherry-picked window.
5. **Bounded geometry**: within an admissible window, the loosest return's
   closest-approach distance must be <= ``geometry_factor`` times the
   window's own tightest closest-approach distance.

This module does NOT decide velocity/flyby-quality (``dv_band``) or seed
selection -- those are separate, later pipeline stages (see the criterion
note's "what this does NOT decide" section).

Pure: numpy only, no CR3BP-specific code -- operates on a plain time series
of ``(t, position)`` samples relative to the target body. A caller wanting
sub-sample crossing precision should sample densely; this module linearly
interpolates crossing TIMES from the given samples but does not re-propagate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class Return:
    """One distinct Hill-sphere-residency episode."""

    t_enter: float
    t_exit: float
    t_closest: float
    closest_distance: float
    closest_position: NDArray[np.float64]


@dataclass(frozen=True)
class AdmissionWindow:
    """One admissible sliding window (criterion §4)."""

    t_start: float
    t_end: float
    returns: tuple[Return, ...]
    geometry_ratio: float  # loosest / tightest closest-approach distance in the window
    geometry_ok: bool


def _interp_crossing(t0: float, t1: float, d0: float, d1: float, r_hill: float) -> float:
    """Linear-interpolated crossing time of ``d(t) = r_hill`` between two samples."""
    if d1 == d0:
        return 0.5 * (t0 + t1)
    frac = (r_hill - d0) / (d1 - d0)
    frac = min(max(frac, 0.0), 1.0)
    return t0 + frac * (t1 - t0)


def find_returns(
    t: NDArray[np.float64],
    positions: NDArray[np.float64],
    r_hill: float,
    *,
    min_separation: float,
) -> list[Return]:
    """Extract distinct returns (criterion §§2-3) from a sampled trajectory.

    ``t`` is a strictly increasing 1D array of sample times; ``positions``
    is ``(len(t), d)`` for any dimension ``d`` (2D or 3D), the position
    RELATIVE TO the target body at each sample.
    """
    dist = np.linalg.norm(positions, axis=1)
    inside = dist < r_hill

    raw_episodes: list[tuple[float, float]] = []
    episode_start: float | None = None
    for i in range(len(t)):
        if inside[i] and episode_start is None:
            episode_start = (
                float(t[i])
                if i == 0
                else _interp_crossing(t[i - 1], t[i], dist[i - 1], dist[i], r_hill)
            )
        elif not inside[i] and episode_start is not None:
            t_exit = _interp_crossing(t[i - 1], t[i], dist[i - 1], dist[i], r_hill)
            raw_episodes.append((episode_start, t_exit))
            episode_start = None
    if episode_start is not None:
        raw_episodes.append((episode_start, float(t[-1])))

    # Merge episodes separated by less than min_separation (criterion §3).
    merged: list[tuple[float, float]] = []
    for start, end in raw_episodes:
        if merged and start - merged[-1][1] < min_separation:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))

    returns: list[Return] = []
    for start, end in merged:
        mask = (t >= start) & (t <= end)
        if not np.any(mask):
            # Degenerate: episode entirely between two samples (very short
            # dip). Use the nearer sample as the closest-approach proxy.
            idx = int(np.argmin(np.abs(t - 0.5 * (start + end))))
            mask = np.zeros_like(t, dtype=bool)
            mask[idx] = True
        sub_dist = dist[mask]
        sub_t = t[mask]
        sub_pos = positions[mask]
        i_min = int(np.argmin(sub_dist))
        returns.append(
            Return(
                t_enter=start,
                t_exit=end,
                t_closest=float(sub_t[i_min]),
                closest_distance=float(sub_dist[i_min]),
                closest_position=np.asarray(sub_pos[i_min], dtype=np.float64),
            )
        )
    return returns


def find_admission_windows(
    returns: list[Return],
    t_span_start: float,
    t_span_end: float,
    *,
    window_lo: float,
    window_hi: float,
    n_returns_lo: int,
    n_returns_hi: int,
    geometry_factor: float,
    window_step: float | None = None,
) -> list[AdmissionWindow]:
    """Slide windows of length in ``[window_lo, window_hi]`` and report every
    admissible one (criterion §§4-5), rather than the first found.

    ``t_span_start``/``t_span_end`` are the trajectory's ACTUAL propagated
    time bounds (not derived from ``returns`` -- a window must stay within
    data that was genuinely checked for returns, an honesty requirement:
    claiming a window is "quiet" past where the propagation stopped would
    be an unverified claim). ``window_step`` defaults to ``window_lo / 20``
    -- fine enough to not miss a qualifying window at the return-epoch scale
    this criterion targets (multi-year returns), coarse enough to keep the
    scan cheap.
    """
    if not returns:
        return []
    if window_step is None:
        window_step = window_lo / 20.0

    windows: list[AdmissionWindow] = []

    for length in (window_lo, window_hi):
        start = t_span_start
        while start + length <= t_span_end + 1e-9:
            end = start + length
            in_window = [r for r in returns if start <= r.t_closest <= end]
            n = len(in_window)
            if n_returns_lo <= n <= n_returns_hi:
                dists = [r.closest_distance for r in in_window]
                tightest, loosest = min(dists), max(dists)
                ratio = loosest / tightest if tightest > 0 else float("inf")
                windows.append(
                    AdmissionWindow(
                        t_start=start,
                        t_end=end,
                        returns=tuple(in_window),
                        geometry_ratio=ratio,
                        geometry_ok=ratio <= geometry_factor,
                    )
                )
            start += window_step
    return windows


__all__ = ["AdmissionWindow", "Return", "find_admission_windows", "find_returns"]
