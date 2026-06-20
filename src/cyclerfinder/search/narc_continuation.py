"""N-arc real-ephemeris continuation lane (#388): Russell parent -> N-arc seed.

This module is the bridge from a Russell idealized parent cycler (the planar,
synodic-normalized :class:`~cyclerfinder.search.generic_return.Cycler` produced by
:func:`~cyclerfinder.search.cycler_assembly.assemble_cycler`) to the per-leg seed
description an N-arc real-ephemeris corrector consumes.

It pairs three sourced inputs:

- the catalogue row's ``sequence_canonical`` (the body chain),
- the SOURCED published per-leg ToFs from
  :func:`~cyclerfinder.search.dsm_descriptor_seed.seed_dsm_chain_from_descriptor`
  (Russell-table arc ToFs + transit times, in days),
- the parent cycler's resonant-return shape (``generic_return.n_revs`` /
  ``generic_return.branch``) to classify each leg.

The output :class:`NarcSeed` carries the per-leg revolution count and Lambert
branch label, a ballistic period (in seconds, scaled by the REAL Earth-Mars
synodic year, not the idealized 2.0 used in the synodic-normalized parent), and
the row's sourced V_inf anchors at Earth and Mars. Pure; raises only when the row
carries no descriptor.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, pi, tau
from typing import TYPE_CHECKING, Any

from cyclerfinder.search.dsm_descriptor_seed import seed_dsm_chain_from_descriptor
from cyclerfinder.search.generic_return import RussellModel

if TYPE_CHECKING:
    from collections.abc import Iterable

    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.cycler_search import Cycler

# Real Earth-Mars synodic-year scaling. The synodic-normalized Russell parent
# uses an idealized 2.0-year synodic period (model units); the real-ephemeris seed
# must scale the parent period ``p`` (in synodic units) back to physical years by
# the REAL synodic year. Mars sidereal period = 1.8808 yr (Earth = 1.0 yr):
#   synodic = 1 / (1/T_E - 1/T_M).
MARS_SIDEREAL_YEARS = 1.8808
DAY_S = 86400.0
YEAR_DAYS = 365.25


@dataclass(frozen=True)
class NarcSeed:
    """Per-leg seed for the N-arc real-ephemeris corrector.

    Attributes
    ----------
    sequence:
        Body chain, e.g. ``("E", "E", "M", "M")``.
    per_leg_revs:
        Lambert revolution count per leg (``nlegs = len(sequence) - 1`` entries).
        0 for a cross-body transit leg; the parent's resonant ``n_revs`` for a
        same-body resonant leg.
    per_leg_branch:
        Lambert branch label per leg: ``"single"`` (single-rev / transit),
        ``"low"`` (slow resonant branch) or ``"high"`` (fast resonant branch).
    tof_seed_days:
        SOURCED published per-leg seed ToFs (days), in sequence order.
    period_sec:
        Ballistic cycler period (seconds), scaled by the REAL Earth-Mars
        synodic year.
    vinf_anchor_e_kms:
        Sourced V_inf at Earth (km/s) from the row's encounter table.
    vinf_anchor_m_kms:
        Sourced V_inf at Mars (km/s) from the row's encounter table.
    """

    sequence: tuple[str, ...]
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    tof_seed_days: tuple[float, ...]
    period_sec: float
    vinf_anchor_e_kms: float
    vinf_anchor_m_kms: float


def _vinf_anchor(row: dict[str, Any], body: str) -> float:
    """Sourced V_inf (km/s) for *body* from the row's encounter table."""
    for e in row.get("vinf_kms_at_encounters") or []:
        if e.get("body") == body and e.get("vinf_kms") is not None:
            return float(e["vinf_kms"])
    raise ValueError(f"row has no V_inf anchor for body {body!r}")


def russell_parent_to_ballistic_seed(
    model: RussellModel,
    cycler: Cycler,
    row: dict[str, Any],
) -> NarcSeed:
    """Bridge a Russell idealized parent cycler to an N-arc real-eph seed.

    Parameters
    ----------
    model:
        The Russell model (carries the ``tu_days`` time-unit; accepted for API
        symmetry with the assembly path).
    cycler:
        The assembled parent cycler (provides ``p`` and the resonant-return
        shape ``generic_return.n_revs`` / ``generic_return.branch``).
    row:
        Raw catalogue row dict (``sequence_canonical``,
        ``vinf_kms_at_encounters``, and a g/G descriptor).

    Returns
    -------
    NarcSeed

    Raises
    ------
    ValueError
        If the row carries no descriptor (no sourced per-leg ToFs).
    """
    sequence = tuple(row["sequence_canonical"].split("-"))

    dsm_seed = seed_dsm_chain_from_descriptor(row)
    if dsm_seed is None:
        raise ValueError("row has no descriptor")
    tof_seed_days = tuple(float(t) for t in dsm_seed.per_leg_tof_days)

    nlegs = len(sequence) - 1
    revs_list: list[int] = []
    branch_list: list[str] = []
    for i in range(nlegs):
        a, b = sequence[i], sequence[i + 1]
        if a != b:
            # Cross-body transit leg.
            revs_list.append(0)
            branch_list.append("single")
        else:
            # Same-body resonant leg: take the parent's resonant return shape.
            revs = int(cycler.generic_return.n_revs)
            if revs == 0:
                branch_list.append("single")
            else:
                branch_list.append("low" if cycler.generic_return.branch == "slow" else "high")
            revs_list.append(revs)

    real_syn_yr = 1.0 / (1.0 / 1.0 - 1.0 / MARS_SIDEREAL_YEARS)
    period_sec = cycler.p * real_syn_yr * YEAR_DAYS * DAY_S

    return NarcSeed(
        sequence=sequence,
        per_leg_revs=tuple(revs_list),
        per_leg_branch=tuple(branch_list),
        tof_seed_days=tof_seed_days,
        period_sec=period_sec,
        vinf_anchor_e_kms=_vinf_anchor(row, "E"),
        vinf_anchor_m_kms=_vinf_anchor(row, "M"),
    )


# Real Earth-Mars synodic period (years). Russell §5.3 derives candidate launch
# epochs by scanning whole-synodic LaunchWindow intervals from the J2000 epoch.
EARTH_MARS_SYNODIC_YEARS = 2.1354


def _wrap_pi(x: float) -> float:
    """Wrap an angle (radians) to the half-open interval ``(-pi, pi]``."""
    w = (x + pi) % tau - pi
    # ``(x + pi) % tau`` lands in ``[0, tau)`` so the result is in ``[-pi, pi)``;
    # map the -pi endpoint to +pi to match the (-pi, pi] convention.
    if w <= -pi:
        w += tau
    return w


def parent_phase_angle(ephem: Ephemeris, t0_sec: float) -> float:
    """Signed in-ecliptic Earth->Mars relative phase angle (radians) at ``t0_sec``.

    The angle from Earth's heliocentric direction to Mars's, measured in the
    ecliptic plane (xy components only), wrapped to ``(-pi, pi]``. Requires a
    real-ephemeris backend (``Ephemeris("astropy")``) for an absolute-calendar
    meaning; the circular backend places both planets at theta=0 at t=0.
    """
    r_e, _ = ephem.state("E", t0_sec)
    r_m, _ = ephem.state("M", t0_sec)
    ang = atan2(float(r_m[1]), float(r_m[0])) - atan2(float(r_e[1]), float(r_e[0]))
    return _wrap_pi(ang)


def candidate_epochs(
    ephem: Ephemeris,
    target_phase: float,
    *,
    launch_window_synodics: Iterable[int] = range(1, 22),
    grid: int = 100,
) -> list[float]:
    """Real-ephemeris epochs (seconds since J2000) at a target Earth-Mars phase.

    Russell §5.3: scan each whole-synodic LaunchWindow interval
    ``[w*T_syn, (w+1)*T_syn)`` from the J2000 epoch on a ``grid``-point grid,
    and record the epoch minimising the wrapped phase error against
    ``target_phase``. One best epoch is returned per window; the list is sorted
    ascending by phase error (best-first).

    Parameters
    ----------
    ephem:
        Real-ephemeris provider (``Ephemeris("astropy")``).
    target_phase:
        Desired signed Earth->Mars relative phase angle (radians).
    launch_window_synodics:
        Whole-synodic LaunchWindow indices to scan (default ``1..21``).
    grid:
        Number of equally-spaced epochs per window interval (default 100).

    Returns
    -------
    list[float]
        Per-window best epochs (seconds since J2000), best-first.
    """
    t_syn = EARTH_MARS_SYNODIC_YEARS * YEAR_DAYS * DAY_S

    best: list[tuple[float, float]] = []  # (phase_error, epoch)
    for w in launch_window_synodics:
        t_start = w * t_syn
        win_best_err = float("inf")
        win_best_t = t_start
        for k in range(grid):
            t = t_start + (k / grid) * t_syn
            err = abs(_wrap_pi(parent_phase_angle(ephem, t) - target_phase))
            if err < win_best_err:
                win_best_err = err
                win_best_t = t
        best.append((win_best_err, win_best_t))

    best.sort(key=lambda pe: pe[0])
    return [t for _, t in best]
