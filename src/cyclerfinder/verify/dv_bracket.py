"""Per-leg ΔV bracket diagnostic - phase-free base floor ≤ DSM ≤ Lambert ceiling.

A cheap sanity band for a single transfer leg's ΔV, after Şaloğlu & Taheri
(2025) §2 (the "ΔV-optimality certificate") and the 2025 mining note
(``docs/notes/2026-06-10-saloglu-2025-iso-impulse-3d-mining.md`` §5-6).

Theory (citation)
-----------------
K. Şaloğlu, E. Taheri, "Classification and Feasibility Assessment of Infinitely
Many Iso-Impulse Three-Dimensional Trajectories," *J. Astronautical Sciences*
(2025), DOI 10.1007/s40295-025-00528-0 (arXiv:2501.01583). For a fixed-time
rendezvous between two orbits the paper establishes a two-sided bound on the
achievable ΔV (their §2 / p. 32):

    base ΔV (phase-free, time-free min of the 2- and 3-impulse base solutions)
        ≤  any fixed-time solution's ΔV
        ≤  the two-impulse Lambert ΔV at the fixed TOF.

The **lower** bound is the phase-free *base solution* - the minimum-ΔV
two-/three-impulse transfer when departure / arrival phases are FREE (it cannot
be beaten by any time-fixed schedule, since freeing the phases only enlarges the
feasible set). The **upper** bound is the ordinary ballistic two-impulse Lambert
arc at the leg's actual TOF (a feasible - usually non-optimal - solution).

What this module is (and is not)
--------------------------------
This is the **cheap diagnostic** of the certificate, NOT the base-solution
optimiser. The base-solution ΔV floor is taken as an INPUT (a sourced number for
goldens, or a caller-supplied floor); reproducing it from orbital elements needs
the Eq.-1 / Eq.-7 phase-free global optimiser the paper used (particleswarm) -
deliberately out of scope here (see the 2025 mining note §5.1 caveat). What this
module does is:

* order-check a leg: ``base ≤ dsm ≤ lambert`` (within tolerance), flagging an
  inversion (a bug - a DSM solution below the phase-free floor is impossible in
  the frozen-element model, and one above the ballistic Lambert means the DSM
  made the leg *worse*);
* flag a leg whose DSM solution sits suspiciously close to the Lambert ceiling
  - i.e. the interior impulse bought little, so there is unexploited DSM
  headroom worth a second look (the headroom fraction toward the base floor).

Honesty (binding, per the 2025 mining note §6.2)
------------------------------------------------
With real ephemerides the departure / arrival "orbits" differ between epochs, so
the frozen-element base bound is exact only in the frozen-element approximation
- a DIAGNOSTIC band, never a hard constraint. A bracket violation is a signal to
investigate (unit error, off-family solve, or a base floor computed under a
different model), not an automatic reject.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.lambert import lambert

Vec3 = NDArray[np.float64]


class BracketVerdict(enum.Enum):
    """Per-leg ΔV bracket verdict."""

    WITHIN_BRACKET = "within_bracket"
    """``base ≤ dsm ≤ lambert`` holds (within tolerance) and the DSM solution is
    not anomalously close to the Lambert ceiling."""

    NEAR_LAMBERT_CEILING = "near_lambert_ceiling"
    """The bracket holds, but the DSM ΔV sits within ``near_ceiling_frac`` of the
    Lambert upper bound (measured across the base→Lambert span): the interior
    impulse bought little, so there may be unexploited DSM headroom on this leg."""

    BELOW_BASE_FLOOR = "below_base_floor"
    """``dsm < base`` beyond tolerance: physically impossible in the
    frozen-element model - a bug signal (unit error / off-family / mismatched
    base model)."""

    ABOVE_LAMBERT_CEILING = "above_lambert_ceiling"
    """``dsm > lambert`` beyond tolerance: the DSM made the leg *worse* than a
    plain ballistic Lambert arc - a bug or a badly-seeded solve."""


@dataclass(frozen=True)
class DvBracket:
    """Outcome of the per-leg ΔV bracket check.

    Attributes
    ----------
    base_floor_kms:
        Phase-free base-solution ΔV lower bound (km/s) - the certificate floor
        (sourced / caller-supplied; NOT computed here).
    dsm_dv_kms:
        The leg's actual (DSM / MBH) ΔV being bracketed, km/s.
    lambert_ceiling_kms:
        Ballistic two-impulse Lambert ΔV at the leg's TOF, km/s (the upper
        bound).
    verdict:
        :class:`BracketVerdict`.
    headroom_fraction:
        Where the DSM solution sits in the band, ``(lambert - dsm) /
        (lambert - base)`` clipped to ``[0, 1]``: ``1`` ≈ at the base floor
        (excellent), ``0`` ≈ at the Lambert ceiling (no DSM benefit). ``0`` when
        the band is degenerate (``lambert ≈ base``).
    """

    base_floor_kms: float
    dsm_dv_kms: float
    lambert_ceiling_kms: float
    verdict: BracketVerdict
    headroom_fraction: float


def lambert_ceiling_dv_kms(
    r1: Vec3,
    v1_before: Vec3,
    r2: Vec3,
    v2_after: Vec3,
    tof_s: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
) -> float:
    """Two-impulse ballistic Lambert ΔV for a leg (the bracket upper bound).

    Solves the zero-revolution Lambert arc ``r1 → r2`` over ``tof_s`` and sums
    the departure burn ``‖v_dep - v1_before‖`` and arrival burn
    ``‖v2_after - v_arr‖`` against the externally-imposed boundary velocities
    (the states the leg must match on each side). This is a feasible - generally
    non-optimal - two-impulse solution, so it upper-bounds the achievable ΔV.

    Parameters
    ----------
    r1, r2:
        Departure and arrival positions, km.
    v1_before:
        Velocity the spacecraft has just BEFORE the departure impulse (km/s) -
        e.g. the planet velocity for a ballistic match, or the incoming leg's
        arrival velocity.
    v2_after:
        Velocity required just AFTER the arrival impulse, km/s.
    tof_s:
        Leg time of flight, seconds (> 0).
    mu:
        Central-body gravitational parameter, km³/s².
    prograde:
        Lambert transfer sense (passed through).
    """
    if tof_s <= 0.0:
        raise ValueError(f"leg TOF must be positive, got {tof_s}")
    sol = lambert(
        np.asarray(r1, dtype=np.float64),
        np.asarray(r2, dtype=np.float64),
        float(tof_s),
        mu=mu,
        prograde=prograde,
    )[0]
    v_dep = np.asarray(sol.v1, dtype=np.float64)
    v_arr = np.asarray(sol.v2, dtype=np.float64)
    dv_dep = float(np.linalg.norm(v_dep - np.asarray(v1_before, dtype=np.float64)))
    dv_arr = float(np.linalg.norm(np.asarray(v2_after, dtype=np.float64) - v_arr))
    return dv_dep + dv_arr


def bracket_leg_dv(
    base_floor_kms: float,
    dsm_dv_kms: float,
    lambert_ceiling_kms: float,
    *,
    tol_kms: float = 1.0e-6,
    near_ceiling_frac: float = 0.05,
) -> DvBracket:
    """Order-check a leg's ΔV against its phase-free floor and Lambert ceiling.

    Parameters
    ----------
    base_floor_kms:
        Phase-free base-solution ΔV lower bound (sourced / supplied), km/s.
    dsm_dv_kms:
        The leg's actual ΔV (DSM / MBH solution), km/s.
    lambert_ceiling_kms:
        Ballistic two-impulse Lambert ΔV at the leg's TOF, km/s - see
        :func:`lambert_ceiling_dv_kms`.
    tol_kms:
        Slack on the order checks (km/s), absorbing numerical noise and the
        frozen-element approximation.
    near_ceiling_frac:
        A leg is flagged ``NEAR_LAMBERT_CEILING`` when its DSM ΔV sits within
        this fraction of the base→Lambert span below the ceiling (i.e.
        ``headroom_fraction < near_ceiling_frac``).

    Returns
    -------
    DvBracket
    """
    span = lambert_ceiling_kms - base_floor_kms
    if span > tol_kms:
        headroom = (lambert_ceiling_kms - dsm_dv_kms) / span
        headroom = float(min(1.0, max(0.0, headroom)))
    else:
        # Degenerate band (Lambert ≈ base): no headroom to speak of.
        headroom = 0.0

    if dsm_dv_kms < base_floor_kms - tol_kms:
        verdict = BracketVerdict.BELOW_BASE_FLOOR
    elif dsm_dv_kms > lambert_ceiling_kms + tol_kms:
        verdict = BracketVerdict.ABOVE_LAMBERT_CEILING
    elif span > tol_kms and headroom < near_ceiling_frac:
        verdict = BracketVerdict.NEAR_LAMBERT_CEILING
    else:
        verdict = BracketVerdict.WITHIN_BRACKET

    return DvBracket(
        base_floor_kms=float(base_floor_kms),
        dsm_dv_kms=float(dsm_dv_kms),
        lambert_ceiling_kms=float(lambert_ceiling_kms),
        verdict=verdict,
        headroom_fraction=headroom,
    )


@dataclass(frozen=True)
class IsoDvSplitFeasibility:
    """Result of the Eq.-16 iso-ΔV impulse-split feasibility probe.

    Attributes
    ----------
    feasible:
        ``True`` when at least one same-direction impulse split (one phasing
        orbit, ``Σ N_k = 1``) fits in the available time without ΔV penalty.
    required_phasing_period_s:
        The phasing-orbit period the split would need, ``surplus_tof / Σ N_k``
        seconds, for the probed revolution count; ``inf`` when ``Σ N_k = 0``.
    period_window_s:
        The admissible phasing-period window ``(T(alpha=0), T(alpha=1))`` seconds - the
        initial-orbit period and the phase-free-arc period. ``feasible`` iff the
        required period lies (inclusively) inside this window.
    max_phasing_revs:
        Upper bound on the number of phasing revolutions, ``floor(surplus_tof /
        T(alpha=0))`` (the §III.E ``n_p < TOF/T(alpha=0)`` count bound); ``0`` when no
        single phasing orbit of the smallest admissible period fits.
    """

    feasible: bool
    required_phasing_period_s: float
    period_window_s: tuple[float, float]
    max_phasing_revs: int


def iso_dv_split_feasible(
    surplus_tof_s: float,
    initial_orbit_period_s: float,
    phase_free_arc_period_s: float,
    *,
    n_phasing_revs: int = 1,
) -> IsoDvSplitFeasibility:
    """Eq.-16 one-line test: does a long multi-rev leg admit iso-ΔV impulse splits?

    Şaloğlu, Taheri & Landau (2023), Sec. III.E, Eq. 16 (mining note
    ``docs/notes/2026-06-10-saloglu-2023-iso-impulse-mining.md`` §3): a same-
    direction anchor-impulse split onto ``Σ N_k`` phasing revolutions is feasible
    **iff** the required phasing-orbit period lands inside the period window
    bounded below by the initial-orbit period ``T(alpha=0)`` and above by the
    phase-free-arc period ``T(alpha=1)``::

        T(alpha=0)  ≤  surplus_tof / Σ N_k  ≤  T(alpha=1)

    The split is **iso-ΔV by construction** - it redistributes the anchor impulse
    across phasing orbits without changing total ΔV, trading only mission time.
    A short (< 1 rev) leg has too little surplus for any phasing orbit and fails
    for every ``Σ N_k ≥ 1`` (note §3): a single interior impulse there is not
    leaving anchor-type savings on the table.

    Our S/L Earth-to-Earth resonant intervals (the S1L1 nomenclature) ARE this
    paper's phasing orbits, so on a long multi-rev cycler leg this probe answers
    whether the encounter-anchor impulse could be split across phasing
    revolutions without a ΔV penalty - useful for impulse-magnitude-capped
    variants and for explaining iso-ΔV plateaus the MBH lane finds. Caveat: our
    "impulses" at an Earth encounter are FLYBYS with bounded turn, not free
    same-direction burns, so the realisable allocation is narrower than the
    free-burn bound this test reports (note §3).

    Parameters
    ----------
    surplus_tof_s:
        Surplus time available for phasing (total leg/mission time minus the
        common terminal coasts and the phase-free arc), seconds.
    initial_orbit_period_s:
        ``T(alpha=0)``, the initial-orbit period (lower period bound), seconds.
    phase_free_arc_period_s:
        ``T(alpha=1)``, the phase-free-arc orbit period (upper period bound),
        seconds.
    n_phasing_revs:
        ``Σ N_k``, the total phasing-revolution count to probe (≥ 1).

    Returns
    -------
    IsoDvSplitFeasibility
    """
    if n_phasing_revs < 1:
        raise ValueError(f"n_phasing_revs must be >= 1, got {n_phasing_revs}")
    lo = float(initial_orbit_period_s)
    hi = float(phase_free_arc_period_s)
    required = float(surplus_tof_s) / float(n_phasing_revs)
    feasible = lo <= required <= hi
    max_revs = int(np.floor(surplus_tof_s / lo)) if lo > 0.0 else 0
    return IsoDvSplitFeasibility(
        feasible=feasible,
        required_phasing_period_s=required,
        period_window_s=(lo, hi),
        max_phasing_revs=max(0, max_revs),
    )


__all__ = [
    "BracketVerdict",
    "DvBracket",
    "IsoDvSplitFeasibility",
    "bracket_leg_dv",
    "iso_dv_split_feasible",
    "lambert_ceiling_dv_kms",
]
