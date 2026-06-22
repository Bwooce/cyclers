"""Band-aware maintenance-ΔV acceptance thresholds (task #420).

Couples the v4.8 Axis-B ``dv_band`` taxonomy (task #415/#417,
``docs/notes/2026-06-22-dv-band-definitions.md``) to the §14 V2/V3 closure
gauntlet so reproduction acceptance is *honest per class*: a strictly-ballistic
row must close to < 1 m/s of deterministic maintenance over 7 cycles, whereas a
powered-DSM row must close to within its STATED powered budget (a non-zero
target — it is not a failure for being powered).

What lives here
---------------
This module is the single mapping from a catalogue row's ``dv_band`` to the
**per-7-cycle deterministic maintenance-ΔV acceptance window** (in m/s), plus a
helper that judges a measured total-over-N-cycles maintenance ΔV against the
band's window. The generic (band-blind) criterion stays exactly where it was —
this module only *parameterises* it when a row carries a band, and falls back to
the caller's generic threshold when ``dv_band`` is ``None``.

Sourcing (all thresholds trace to the band-definitions note)
------------------------------------------------------------
The three ballistic ceilings are the Russell-Ocampo 2006 JGCD 29(2) census
(< 1 / < 10 / < 300 m/s total deterministic maneuver over a 7-cycle
real-ephemeris propagation at the best launch window) and scale pro-rata with
the propagated cycle count. The powered-DSM window is ``[300 m/s, 3.5 km/s x
7]`` per 7 cycles — its **lower** bound is the top of Russell's net (a powered
cycler is, by definition, above the low-maintenance ceiling and must close to a
*strictly-positive* budget, not zero), and its **upper** bound is the existing
V2-powered sanity ceiling (3.5 km/s/cycle =
:data:`cyclerfinder.verify.v2_powered._MAINTENANCE_DV_SANITY_MAX_KMS`, expressed
in the 7-cycle basis) that rejects the off-family ~55 km/s degenerate basin
(finding #114). The powered upper bound is a project-convention degenerate
rejector, NOT a sourced tight tier: the V2-powered gate's MEASURED ΔV is an
over-estimating turn-deficit surrogate (~2.9 km/s/cycle) that must not be
asserted against the sourced Aldrin budget (~1.73-2.04 km/s / 15 yr,
Byrnes-Longuski-Aldrin 1993) — see the ``_POWERED_DSM_*`` basis note (conflict
C4). The powered band's binding criterion is therefore its strictly-positive
floor (rejects the ΔV≈0 ballistic neighbour), with the sanity ceiling on top.

``low_thrust_sep`` carries no impulsive m/s ceiling — it is a regime, not a
magnitude (band-definitions §"Proposed bands"). It is reported as a propellant /
Isp-converted budget elsewhere (:mod:`cyclerfinder.search.lowthrust_maintenance`)
and is intentionally *not* band-accepted here; it returns ``None`` (no impulsive
window) so the caller keeps the generic path rather than silently passing.

Golden discipline
-----------------
The *measured* maintenance ΔV the caller passes in is OUR computed value
(McConaghy 2002 defers the magnitude); these thresholds are the published
band *boundaries*, never asserted as a rediscovered per-row number. The powered
upper bound and the ``None`` low-thrust handling are project-convention,
flagged as such above.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

# Russell's native basis: the m/s tiers are TOTAL deterministic maneuver over a
# 7-cycle real-ephemeris propagation. A gate that runs a different cycle count
# scales the window pro-rata (band-definitions §"Mandatory basis").
RUSSELL_BASIS_CYCLES: Final[int] = 7

# --- Russell-Ocampo 2006 JGCD 29(2) census ceilings (total over 7 cycles, m/s).
_STRICTLY_BALLISTIC_MAX_MPS: Final[float] = 1.0
_ESSENTIALLY_BALLISTIC_MAX_MPS: Final[float] = 10.0
_LOW_MAINTENANCE_MAX_MPS: Final[float] = 300.0

# --- powered-DSM window (total over 7 cycles, m/s).
#
# IMPORTANT basis note (task #420, conflict C4 of the band-definitions doc).
# The powered band differs from the three ballistic tiers in TWO ways:
#   1. Its acceptance is a strictly-POSITIVE floor, not an upper ceiling — a
#      powered cycler must close to a real budget above Russell's net (≥ 300
#      m/s / 7 cycles), so a ΔV≈0 ballistic neighbour does NOT satisfy it.
#   2. Its SOURCED budget (Aldrin ~1.73-2.04 km/s / 15 yr ≈ per 7 synodic
#      cycles, Byrnes-Longuski-Aldrin 1993) is NOT directly comparable to the
#      V2-powered gate's MEASURED maintenance ΔV, which is the over-estimating
#      turn-deficit surrogate (~2.9 km/s PER cycle, McConaghy 2002 defers the
#      true magnitude). Asserting the gate's surrogate against the sourced
#      budget would wrongly fail the recorded-V2 Aldrin (golden-discipline
#      violation — a computed surrogate must never be matched to a sourced
#      number). So the powered UPPER bound is the existing project-convention
#      V2-powered sanity ceiling (3.5 km/s/cycle; rejects the off-family
#      ~55 km/s basin, finding #114), expressed in the 7-cycle basis so the
#      pro-rata scaling returns it to 3.5 km/s/cycle unchanged. It is a degenerate
#      rejector, NOT a sourced tight tier.
_POWERED_DSM_LOWER_MPS: Final[float] = 300.0
# 3.5 km/s/cycle * 7 cycles (project-convention sanity ceiling, not a sourced tier).
_POWERED_DSM_UPPER_MPS: Final[float] = 3_500.0 * RUSSELL_BASIS_CYCLES


@dataclass(frozen=True)
class BandThreshold:
    """The per-7-cycle deterministic maintenance-ΔV acceptance window for a band.

    Attributes
    ----------
    band:
        The ``dv_band`` enum value this window is for.
    lower_mps:
        Inclusive lower bound on the total deterministic maintenance ΔV over the
        Russell 7-cycle basis (m/s). ``0.0`` for the ballistic bands (a cheaper
        close is always acceptable); the powered floor (300 m/s) for
        ``powered_dsm`` — a powered cycler must close to a *positive* budget, so
        a ΔV≈0 ballistic neighbour does NOT satisfy the powered band.
    upper_mps:
        Inclusive upper bound on the total deterministic maintenance ΔV over the
        7-cycle basis (m/s). For the three ballistic bands this is the **sourced**
        Russell tier ceiling (< 1 / < 10 / < 300 m/s). For ``powered_dsm`` it is
        the **project-convention** V2-powered sanity ceiling (3.5 km/s/cycle in
        the 7-cycle basis) — a degenerate-basin rejector, NOT a sourced tier;
        the powered band's binding criterion is its strictly-positive ``lower_mps``
        floor (see the module-level basis note, conflict C4).
    """

    band: str
    lower_mps: float
    upper_mps: float


# The single band → per-7-cycle window map. ``low_thrust_sep`` is intentionally
# absent (regime, not impulsive magnitude — handled as a propellant/Isp budget
# elsewhere); ``dv_band_threshold`` returns ``None`` for it and for any unknown
# / null band, which is the caller's signal to keep the generic criterion.
_BAND_WINDOWS: Final[dict[str, BandThreshold]] = {
    "strictly_ballistic": BandThreshold("strictly_ballistic", 0.0, _STRICTLY_BALLISTIC_MAX_MPS),
    "essentially_ballistic": BandThreshold(
        "essentially_ballistic", 0.0, _ESSENTIALLY_BALLISTIC_MAX_MPS
    ),
    "low_maintenance": BandThreshold("low_maintenance", 0.0, _LOW_MAINTENANCE_MAX_MPS),
    "powered_dsm": BandThreshold("powered_dsm", _POWERED_DSM_LOWER_MPS, _POWERED_DSM_UPPER_MPS),
}


def dv_band_threshold(dv_band: str | None) -> BandThreshold | None:
    """Map a ``dv_band`` to its per-7-cycle maintenance-ΔV acceptance window.

    Parameters
    ----------
    dv_band:
        The catalogue row's ``dv_band`` (one of the five v4.8 enum values, or
        ``None``). ``None`` is the common, valid default (most rows carry no
        band) and returns ``None`` — the caller's signal to fall back to the
        existing generic criterion.

    Returns
    -------
    BandThreshold | None
        The per-7-cycle acceptance window for the band, or ``None`` when the
        band is ``None``, unknown, or ``"low_thrust_sep"`` (no impulsive m/s
        ceiling — a propellant/Isp regime, not a magnitude). A row therefore
        **cannot be band-promoted on a band it does not carry** — the generic
        path always applies in the ``None`` case.
    """
    if dv_band is None:
        return None
    return _BAND_WINDOWS.get(dv_band)


def accept_maintenance_dv(
    total_maintenance_dv_mps: float,
    *,
    dv_band: str | None,
    n_cycles: int,
    generic_max_mps: float,
) -> bool:
    """Judge a measured total maintenance ΔV against the band-aware threshold.

    When ``dv_band`` carries an impulsive window (the four magnitude bands), the
    measured total deterministic maintenance ΔV over ``n_cycles`` must fall
    within the band's window **scaled to ``n_cycles``** (Russell's native basis
    is 7 cycles; a gate running ``n_cycles`` compares against
    ``window * n_cycles / 7``). When ``dv_band`` is ``None`` /
    ``"low_thrust_sep"`` / unknown, the band path does not apply and the
    measured ΔV is judged against the caller's existing ``generic_max_mps`` (the
    band-blind criterion preserved verbatim).

    Parameters
    ----------
    total_maintenance_dv_mps:
        Measured TOTAL deterministic maintenance ΔV over ``n_cycles`` (m/s),
        non-negative.
    dv_band:
        The catalogue row's ``dv_band`` (or ``None``).
    n_cycles:
        Number of consecutive cycles the gate actually propagated (>= 1). Used
        to pro-rata the 7-cycle band window.
    generic_max_mps:
        The existing band-blind upper acceptance bound (m/s) over the same
        ``n_cycles`` — the fallback applied when no band window exists. A lower
        bound of 0 is implied for the generic path.

    Returns
    -------
    bool
        ``True`` when the measured ΔV is accepted for the row's class.
    """
    if total_maintenance_dv_mps < 0.0:
        raise ValueError(
            f"total_maintenance_dv_mps must be non-negative, got {total_maintenance_dv_mps}"
        )
    if n_cycles < 1:
        raise ValueError(f"n_cycles must be >= 1, got {n_cycles}")

    window = dv_band_threshold(dv_band)
    if window is None:
        # Generic (band-blind) fallback: preserved verbatim. A row cannot be
        # promoted on a band it does not carry.
        return total_maintenance_dv_mps <= generic_max_mps

    scale = n_cycles / RUSSELL_BASIS_CYCLES
    lower = window.lower_mps * scale
    upper = window.upper_mps * scale
    return lower <= total_maintenance_dv_mps <= upper


__all__ = [
    "RUSSELL_BASIS_CYCLES",
    "BandThreshold",
    "accept_maintenance_dv",
    "dv_band_threshold",
]
