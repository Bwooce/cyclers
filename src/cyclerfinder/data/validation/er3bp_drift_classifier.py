"""Bounded-vs-divergent long-span drift classifier for the geocentric ER3BP (#583).

#581 stage 2 validated the Gurfil-Kasdin niching-GA layer against the paper's
own 12 published optimization sets, using the paper's own 1-revolution (1
Earth year) fitness window (``core/er3bp_geocentric.py::gurfil_kasdin_fitness``,
Eq. 15). #583 widens the search DOMAIN past those 12 sets. A GA fitness peak
under a 1-rev window is a basin *indicator*, not evidence the orbit stays
bounded over any longer span -- exactly the concern
[[feedback_orbit_closure_discipline]] and
[[feedback_verify_gauntlet_with_positive_control]] exist to catch. This module
is the missing long-span gate: propagate a candidate well past the 1-rev
fitness window and classify it BOUNDED (stationary geocentric r-band /
recurrence) or DIVERGENT (secular rmax growth or escape/collision).

Closest existing template is
:mod:`cyclerfinder.data.validation.v2_3d` (propagate-and-gate-drift), but that
module's design assumes a corrector-refined PERIOD (``n_cycles`` consecutive
returns to a known ``T``). The geocentric ER3BP quasi-orbits this module
judges have no such period -- Gurfil-Kasdin's own families are QUASI-periodic
bounded-oscillation orbits (DRO/DPO/ERO/DEO), not closed periodic orbits. So
instead of "drift at each period return", this classifier windows the long
propagation into ``n_revs`` consecutive 1-year (1-revolution) blocks --
matching :func:`~cyclerfinder.core.er3bp_geocentric.gurfil_kasdin_fitness`'s
own window convention exactly -- and tests whether the per-window r-band
(``rmin_i``, ``rmax_i``) is STATIONARY (bounded quasi-periodic oscillation) or
exhibits a SECULAR TREND / escape (divergent).

Horizon choice (``N_REVS_DEFAULT = 50``): the task spec calls for "N ~= 50-100,
justify the exact number". 50 (Earth years) is chosen as the LOWER end of that
range because:
  * It is 50x the 1-rev GA fitness window -- two orders of magnitude beyond
    what the GA ever saw, comfortably long enough for any secular instability
    with an e-folding time of a few years (the timescale relevant at these
    geocentric distances, 1e6-1e7+ km, where solar perturbation dominates
    Earth's gravity) to become visible as either escape or an unmistakable
    r-band trend.
  * It stays computationally tractable for a BATCH of checks (11 known-good
    families + 1 escaping control + N partition survivors), each ~50 x 500 =
    25,000-sample DOP853 integrations at rtol=1e-9 -- the same tolerance
    stage 2's own 5-year extended check already used successfully
    (``run_581_gurfil_reproduction.py::characterize``).
  * :func:`spot_check_theta0_robustness` and callers running a slower
    hardening pass MAY re-run at ``n_revs=100`` for extra margin (the module
    constant is a default, not a hard ceiling); a doubled horizon is the
    natural next check for any BOUNDED verdict a caller wants extra
    confidence in, per [[project_388_wall_energy_selective]]'s
    epoch/horizon-fragility lesson.

Numeric thresholds (written into the code per the task's own mandate, not
left to a builder default):
  * ``RMAX_GROWTH_RATIO_THRESHOLD = 3.0`` -- if the mean rmax of the LAST
    quarter of windows exceeds 3x the mean rmax of the FIRST quarter, the
    orbit is DIVERGENT (unambiguous secular growth). 3x mirrors the same
    order-of-magnitude tolerance #581's own family-matching criterion used
    for feature ratios (``FEATURE_FACTOR_TOL = 2.0`` in
    ``run_581_gurfil_reproduction.py``, loosened slightly here because a
    50-year windowed trend is noisier than a single 1-year measurement).
  * ``RBAND_TREND_FRACTION_THRESHOLD = 0.30`` -- a linear least-squares fit
    of rmax_i vs window index; if the fitted total drift over the run
    (``|slope| * n_windows``) exceeds 30% of the mean rmax, the orbit is
    DIVERGENT. This catches slow, still-secular drift that a first/last
    quartile ratio alone might miss (e.g. a shallow but consistent trend that
    hasn't yet tripled by window 50).
  * Collision (``r <= R_EARTH_NORM``) or escape (``r >= escape_radius``) at
    ANY point during the run is an immediate DIVERGENT verdict (mirrors
    ``gurfil_kasdin_fitness``'s own death-penalty convention, Eq. 17).
  * Fewer than ``n_revs`` COMPLETE windows (i.e. the propagation terminated
    early on collision/escape, or produced too few sample points near the
    tail) is treated conservatively as DIVERGENT -- a truncated run is never
    upgraded to "bounded" by default.

This classifier needs its OWN positive control (the paper's own
"not-found is necessary-not-sufficient" discipline applies symmetrically to
a new gate): :func:`classify_bounded_drift` must classify all 11 of #581
stage 2's known-good reproduced families as BOUNDED, and a deliberately
escaping test IC as DIVERGENT, BEFORE it is trusted on anything new. See
``tests/data/test_er3bp_drift_classifier.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.er3bp import ER3BPSystem, er3bp_eom
from cyclerfinder.core.er3bp_geocentric import (
    A_AU_KM_GURFIL_KASDIN,
    R_EARTH_NORM,
    SUN_EARTH_ER3BP,
    geocentric_to_barycentric,
)

N_REVS_DEFAULT: Final[int] = 50
"""Default long-span horizon in revolutions (1 rev = 1 Earth year in this
frame's true-anomaly independent variable). See module docstring for the
justification (50x past the 1-rev GA fitness window; tractable in batch)."""

SAMPLES_PER_REV: Final[int] = 500
"""Samples per 1-year window used to resolve rmin/rmax. The natural orbital
timescale at Gurfil-Kasdin's own 1e6-1e7 km geocentric distances is of order
weeks-to-months (well under a year), so 500 samples/year comfortably
resolves multiple turning points per window without the O(1e4-1e5)-point
memory cost of matching the 1-rev fitness function's own finer 2000-4000
sample density over a 50x longer span."""

ESCAPE_RADIUS_NORM_DEFAULT: Final[float] = 0.5
"""Normalized (AU) escape radius. Matches
:func:`cyclerfinder.core.er3bp_geocentric.gurfil_kasdin_fitness`'s own
default exactly, so the same "beyond this we don't even try to classify a
trend, it's just gone" boundary applies at both the 1-rev fitness stage and
the long-span drift-classification stage."""

RMAX_GROWTH_RATIO_THRESHOLD: Final[float] = 3.0
"""DIVERGENT if mean(rmax, last quarter of windows) / mean(rmax, first
quarter of windows) exceeds this. See module docstring for the derivation
(loosened from #581's own 2.0 feature-ratio tolerance to absorb 50-window
noise)."""

RBAND_TREND_FRACTION_THRESHOLD: Final[float] = 0.30
"""DIVERGENT if a linear fit of rmax_i vs window index predicts a total
drift, over the whole run, exceeding this fraction of the mean rmax. Catches
slow secular trends the quartile-ratio test alone might miss."""


@dataclass(frozen=True)
class DriftVerdict:
    """Frozen bounded-vs-divergent verdict for one long-span propagation.

    ``bounded`` is the headline boolean. All other fields are audit trail.
    """

    candidate_id: str
    theta0: float
    n_revs_requested: int
    n_windows_complete: int
    terminated_early: bool
    termination_reason: str
    rmin_per_window_km: tuple[float, ...]
    rmax_per_window_km: tuple[float, ...]
    growth_ratio: float
    trend_fraction: float
    growth_ratio_threshold: float
    trend_fraction_threshold: float
    bounded: bool
    notes: str = ""


def classify_bounded_drift(
    state6_geo: NDArray[np.float64],
    theta0: float,
    sys: ER3BPSystem = SUN_EARTH_ER3BP,
    *,
    n_revs: int = N_REVS_DEFAULT,
    samples_per_rev: int = SAMPLES_PER_REV,
    escape_radius: float = ESCAPE_RADIUS_NORM_DEFAULT,
    rtol: float = 1e-9,
    atol: float = 1e-9,
    growth_ratio_threshold: float = RMAX_GROWTH_RATIO_THRESHOLD,
    trend_fraction_threshold: float = RBAND_TREND_FRACTION_THRESHOLD,
    candidate_id: str = "",
) -> DriftVerdict:
    """Propagate ``n_revs`` revolutions past the 1-rev GA fitness window and
    classify the geocentric r-band as bounded (stationary / recurrent) or
    divergent (secular growth or escape/collision).

    Parameters mirror :func:`~cyclerfinder.core.er3bp_geocentric.gurfil_kasdin_fitness`
    where they overlap (``escape_radius``, ``rtol``/``atol``) so a candidate's
    fitness-stage and drift-classification-stage runs use the same physical
    boundaries.
    """
    if n_revs < 4:
        raise ValueError(
            f"n_revs must be >= 4 (need at least first/last quartile windows); got {n_revs}"
        )
    mu, e = sys.mu, sys.e
    bary0 = geocentric_to_barycentric(np.asarray(state6_geo, dtype=float), mu)
    offset = 1.0 - mu

    def geo_r2(_f: float, s: NDArray[np.float64]) -> float:
        dx = s[0] - offset
        return float(dx * dx + s[1] * s[1] + s[2] * s[2])

    def collision(f: float, s: NDArray[np.float64], *_a: float) -> float:
        return geo_r2(f, s) - R_EARTH_NORM**2

    def escape(f: float, s: NDArray[np.float64], *_a: float) -> float:
        return geo_r2(f, s) - escape_radius**2

    collision.terminal = True  # type: ignore[attr-defined]
    escape.terminal = True  # type: ignore[attr-defined]

    f_end = theta0 + n_revs * 2.0 * math.pi
    n_points = n_revs * samples_per_rev + 1
    t_eval = np.linspace(theta0, f_end, n_points)

    sol = solve_ivp(
        er3bp_eom,
        (theta0, f_end),
        bary0,
        args=(mu, e),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        events=(collision, escape),
    )

    terminated_early = sol.status == 1
    termination_reason = ""
    t_events = sol.t_events
    if terminated_early and t_events is not None:
        if len(t_events[0]) > 0:
            termination_reason = "collision"
        elif len(t_events[1]) > 0:
            termination_reason = "escape"
        else:
            termination_reason = "unknown-terminal-event"
    elif terminated_early:
        termination_reason = "unknown-terminal-event"
    elif not sol.success:
        termination_reason = "integration-failure"

    if sol.y.shape[1] < 2:
        return DriftVerdict(
            candidate_id=candidate_id,
            theta0=theta0,
            n_revs_requested=n_revs,
            n_windows_complete=0,
            terminated_early=True,
            termination_reason=termination_reason or "no-samples",
            rmin_per_window_km=(),
            rmax_per_window_km=(),
            growth_ratio=float("inf"),
            trend_fraction=float("inf"),
            growth_ratio_threshold=growth_ratio_threshold,
            trend_fraction_threshold=trend_fraction_threshold,
            bounded=False,
            notes="propagation produced <2 samples -- immediate divergent",
        )

    t = sol.t
    dx = sol.y[0] - offset
    r = np.sqrt(dx * dx + sol.y[1] ** 2 + sol.y[2] ** 2)
    window_idx = np.floor((t - theta0) / (2.0 * math.pi)).astype(int)

    rmin_per_window: list[float] = []
    rmax_per_window: list[float] = []
    for i in range(n_revs):
        mask = window_idx == i
        if not np.any(mask):
            break  # ran out of data: early termination mid-run
        rmin_per_window.append(float(r[mask].min()))
        rmax_per_window.append(float(r[mask].max()))

    n_complete = len(rmax_per_window)
    rmax_arr = np.array(rmax_per_window)

    if n_complete < n_revs or n_complete < 4:
        # Truncated run (collision/escape event fired, or pathologically few
        # samples near the tail): conservative default is DIVERGENT, never
        # upgraded to bounded on partial data.
        ratio = float("inf")
        trend_fraction = float("inf")
        bounded = False
        notes = (
            f"only {n_complete}/{n_revs} windows complete "
            f"({termination_reason or 'truncated'}) -- treated as divergent"
        )
    else:
        q = max(1, n_complete // 4)
        first_q_mean = float(np.mean(rmax_arr[:q]))
        last_q_mean = float(np.mean(rmax_arr[-q:]))
        ratio = last_q_mean / first_q_mean if first_q_mean > 0 else float("inf")

        idx = np.arange(n_complete, dtype=float)
        slope, _intercept = np.polyfit(idx, rmax_arr, 1)
        mean_rmax = float(np.mean(rmax_arr))
        trend_fraction = (
            abs(float(slope) * n_complete) / mean_rmax if mean_rmax > 0 else float("inf")
        )

        bounded = (
            not terminated_early
            and ratio <= growth_ratio_threshold
            and trend_fraction <= trend_fraction_threshold
        )
        notes = (
            f"growth_ratio={ratio:.3f} (thr {growth_ratio_threshold}), "
            f"trend_fraction={trend_fraction:.3f} (thr {trend_fraction_threshold})"
        )

    return DriftVerdict(
        candidate_id=candidate_id,
        theta0=theta0,
        n_revs_requested=n_revs,
        n_windows_complete=n_complete,
        terminated_early=terminated_early,
        termination_reason=termination_reason,
        rmin_per_window_km=tuple(v * A_AU_KM_GURFIL_KASDIN for v in rmin_per_window),
        rmax_per_window_km=tuple(v * A_AU_KM_GURFIL_KASDIN for v in rmax_per_window),
        growth_ratio=ratio,
        trend_fraction=trend_fraction,
        growth_ratio_threshold=growth_ratio_threshold,
        trend_fraction_threshold=trend_fraction_threshold,
        bounded=bounded,
        notes=notes,
    )


def spot_check_theta0_robustness(
    state6_geo: NDArray[np.float64],
    theta0_base: float,
    sys: ER3BPSystem = SUN_EARTH_ER3BP,
    *,
    n_revs: int = N_REVS_DEFAULT,
    phase_offsets: tuple[float, ...] = (
        2.0 * math.pi / 3.0,
        4.0 * math.pi / 3.0,
    ),
    **kwargs: object,
) -> dict[float, bool]:
    """Cheap theta0-robustness spot-check for a BOUNDED verdict.

    Re-runs :func:`classify_bounded_drift` at ``theta0_base`` plus each of
    ``phase_offsets`` (default: two other epochs at +/-120 deg true-anomaly
    phase, i.e. testing 3 roughly evenly-spaced points around the orbit) and
    returns ``{theta0: bounded}`` for all three. Per
    [[project_388_wall_energy_selective]]'s epoch-fragility lesson, a
    "bounded" verdict that flips at a different phase is NOT trustworthy
    without noting the fragility -- this is diagnostic-only, it does not
    itself gate anything; the caller decides how to weigh a flip.
    """
    thetas = (theta0_base, *(theta0_base + off for off in phase_offsets))
    out: dict[float, bool] = {}
    for th in thetas:
        verdict = classify_bounded_drift(
            state6_geo,
            th,
            sys,
            n_revs=n_revs,
            **kwargs,  # type: ignore[arg-type]
        )
        out[th] = verdict.bounded
    return out


__all__ = [
    "ESCAPE_RADIUS_NORM_DEFAULT",
    "N_REVS_DEFAULT",
    "RBAND_TREND_FRACTION_THRESHOLD",
    "RMAX_GROWTH_RATIO_THRESHOLD",
    "SAMPLES_PER_REV",
    "DriftVerdict",
    "classify_bounded_drift",
    "spot_check_theta0_robustness",
]
