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


def classify_dv_band(
    total_maintenance_dv_mps: float, *, n_cycles: int = RUSSELL_BASIS_CYCLES
) -> str:
    """Assign a ``dv_band`` from a MEASURED real-ephemeris maintenance ΔV (task #422).

    This is the *output* direction of the band↔threshold map: where
    :func:`dv_band_threshold` answers "what window does this *literature* band
    imply?", this answers "what band does this *measured* real-ephemeris
    maintenance ΔV fall in?". It is what lets the ~215 null-band rows acquire a
    band by being reproduced (rather than only carrying a literature label).

    The classification is on Russell's native basis — TOTAL deterministic
    maintenance ΔV over a 7-cycle real-ephemeris propagation, best launch window
    (band-definitions §"Mandatory basis"). A gate that propagated a different
    cycle count passes its real ``n_cycles``; the measured total is scaled
    pro-rata to the 7-cycle basis *before* binning, so the returned band is
    always expressed in Russell's basis regardless of how many cycles were run.

    Bins (the sourced Russell-Ocampo 2006 JGCD 29(2) ceilings, on the 7-cycle
    basis):

    * ``< 1 m/s``   → ``"strictly_ballistic"``
    * ``< 10 m/s``  → ``"essentially_ballistic"``
    * ``< 300 m/s`` → ``"low_maintenance"``
    * ``>= 300 m/s`` → ``"powered_dsm"``

    Note the bins are half-open at the *top* (``< ceiling``), matching the
    literature wording ("less than 1 m/s", etc.) and the inclusive-upper windows
    of :func:`dv_band_threshold` (a value exactly ON a ballistic ceiling, e.g.
    1.0 m/s, falls into the *next* band — it is not ``< 1``). ``low_thrust_sep``
    is never returned: it is a propulsion regime, not a measured impulsive
    magnitude, so it cannot be inferred from a ΔV number.

    Parameters
    ----------
    total_maintenance_dv_mps:
        Measured TOTAL deterministic maintenance ΔV over ``n_cycles`` (m/s),
        non-negative.
    n_cycles:
        Number of consecutive cycles the gate actually propagated (>= 1).
        Defaults to the 7-cycle Russell basis (the common case where the
        measurement is already a 7-cycle total — pass it through unscaled).

    Returns
    -------
    str
        One of ``"strictly_ballistic"`` / ``"essentially_ballistic"`` /
        ``"low_maintenance"`` / ``"powered_dsm"``.
    """
    if total_maintenance_dv_mps < 0.0:
        raise ValueError(
            f"total_maintenance_dv_mps must be non-negative, got {total_maintenance_dv_mps}"
        )
    if n_cycles < 1:
        raise ValueError(f"n_cycles must be >= 1, got {n_cycles}")

    # Scale the measured total to Russell's 7-cycle basis before binning.
    dv_7cycle = total_maintenance_dv_mps * RUSSELL_BASIS_CYCLES / n_cycles
    if dv_7cycle < _STRICTLY_BALLISTIC_MAX_MPS:
        return "strictly_ballistic"
    if dv_7cycle < _ESSENTIALLY_BALLISTIC_MAX_MPS:
        return "essentially_ballistic"
    if dv_7cycle < _LOW_MAINTENANCE_MAX_MPS:
        return "low_maintenance"
    return "powered_dsm"


# Band ordering by maintenance cost (cheapest → most powered). Used to decide
# whether a MEASURED band is *worse* (more powered) than a SOURCED band, which
# is the human-review-flag condition (per orbit-closure discipline: a measured
# band that costs MORE than the sourced literature label is a real conflict;
# a measured band that is cheaper is consistent with — and bounded by — the
# sourced ceiling). ``low_thrust_sep`` is intentionally absent (a regime, not a
# point on this cost scale); a sourced ``low_thrust_sep`` never participates in
# the impulsive measured-vs-sourced comparison.
_BAND_COST_ORDER: Final[tuple[str, ...]] = (
    "strictly_ballistic",
    "essentially_ballistic",
    "low_maintenance",
    "powered_dsm",
)


@dataclass(frozen=True)
class DvBandClassification:
    """The dv_band assignment for a reproduced row + its provenance (task #422).

    Attributes
    ----------
    dv_band:
        The band assigned to the row. When the row carried a *sourced* band this
        is the sourced band UNCHANGED (a sourced literature value is never
        silently overwritten — orbit-closure discipline); otherwise it is the
        :func:`classify_dv_band` of the measured ΔV.
    dv_band_source:
        Provenance marker. ``"computed-v3"`` when the band was assigned from the
        measured real-ephemeris maintenance ΔV; ``"sourced"`` when the row
        already carried a literature band that was kept.
    measured_band:
        The band the MEASURED ΔV classifies into (always populated). Equals
        ``dv_band`` in the computed case; may differ from it in the sourced case
        (that difference is what ``mismatch`` reports).
    measured_total_dv_mps:
        The measured total deterministic maintenance ΔV (m/s) that was
        classified, echoed for provenance.
    n_cycles:
        The cycle count the measurement spanned.
    sourced_band:
        The row's pre-existing sourced band, or ``None`` if it carried none.
    mismatch:
        ``True`` when the row had a sourced band AND the measured band is
        STRICTLY WORSE (more powered) than it — a real-ephemeris reproduction
        that costs more than the literature claim. This is surfaced for human
        review; it never auto-overwrites the sourced value.
    detail:
        Human-readable summary of the assignment / mismatch.
    """

    dv_band: str
    dv_band_source: str
    measured_band: str
    measured_total_dv_mps: float
    n_cycles: int
    sourced_band: str | None
    mismatch: bool
    detail: str


def assign_dv_band_from_measurement(
    total_maintenance_dv_mps: float,
    *,
    n_cycles: int,
    sourced_dv_band: str | None,
) -> DvBandClassification:
    """Assign a row's dv_band from a measured ΔV, respecting any sourced band.

    The #422 "dv_band as a reproduction OUTPUT" rule, with the orbit-closure
    discipline conflict-surfacing guarantee:

    * **No sourced band** (the ~215 null rows): assign the measured band
      (:func:`classify_dv_band`) with ``dv_band_source="computed-v3"``. The row
      acquires a band *by being reproduced*.
    * **Has a sourced band**: keep the sourced band (``dv_band_source="sourced"``,
      ``dv_band`` UNCHANGED — a literature value is never silently overwritten).
      If the measured band is strictly WORSE (more powered) than the sourced one,
      set ``mismatch=True`` so a human reviews the conflict. A measured band that
      is *cheaper* than (or equal to) the sourced band is consistent and not
      flagged. A sourced ``low_thrust_sep`` (a regime, off the impulsive cost
      scale) never produces a mismatch from an impulsive measurement.

    Parameters
    ----------
    total_maintenance_dv_mps:
        Measured TOTAL deterministic maintenance ΔV over ``n_cycles`` (m/s).
    n_cycles:
        Cycle count the measurement spanned (>= 1).
    sourced_dv_band:
        The row's pre-existing literature band, or ``None``.

    Returns
    -------
    DvBandClassification
    """
    measured = classify_dv_band(total_maintenance_dv_mps, n_cycles=n_cycles)

    if sourced_dv_band is None:
        return DvBandClassification(
            dv_band=measured,
            dv_band_source="computed-v3",
            measured_band=measured,
            measured_total_dv_mps=float(total_maintenance_dv_mps),
            n_cycles=n_cycles,
            sourced_band=None,
            mismatch=False,
            detail=(
                f"computed-v3: assigned '{measured}' from measured "
                f"{total_maintenance_dv_mps:.4g} m/s over {n_cycles} cycle(s)"
            ),
        )

    # Row carries a sourced band: keep it, compare, flag a worse measurement.
    mismatch = _measured_band_is_worse(measured=measured, sourced=sourced_dv_band)
    if mismatch:
        detail = (
            f"MISMATCH: measured band '{measured}' "
            f"({total_maintenance_dv_mps:.4g} m/s over {n_cycles} cycle(s)) is more "
            f"powered than sourced band '{sourced_dv_band}' — flagged for human "
            f"review; sourced band NOT overwritten"
        )
    else:
        detail = (
            f"sourced: kept '{sourced_dv_band}'; measured '{measured}' "
            f"({total_maintenance_dv_mps:.4g} m/s over {n_cycles} cycle(s)) is "
            f"consistent (not more powered)"
        )
    return DvBandClassification(
        dv_band=sourced_dv_band,
        dv_band_source="sourced",
        measured_band=measured,
        measured_total_dv_mps=float(total_maintenance_dv_mps),
        n_cycles=n_cycles,
        sourced_band=sourced_dv_band,
        mismatch=mismatch,
        detail=detail,
    )


def _measured_band_is_worse(*, measured: str, sourced: str) -> bool:
    """True when ``measured`` is strictly more powered (costlier) than ``sourced``.

    Both must be impulsive magnitude bands for the comparison to apply; a band
    off the impulsive cost scale (``low_thrust_sep`` or any unknown band) is not
    comparable and yields ``False`` (no mismatch flagged — the conflict, if any,
    is not an impulsive-magnitude one).
    """
    if measured not in _BAND_COST_ORDER or sourced not in _BAND_COST_ORDER:
        return False
    return _BAND_COST_ORDER.index(measured) > _BAND_COST_ORDER.index(sourced)


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


# --- #424: formalise the §14 / #175 V3 ballistic-vs-powered class-split ---------
#
# spec §14 V3 has a class-split (#175): a BALLISTIC row's reproduced maintenance
# ΔV is held to a generic 120 m/s budget; a POWERED row is held to its OWN
# documented ΔV budget (x a 1.10 margin) — NOT failed for being non-ballistic.
# This was applied by hand in catalogue evidence text; the function below makes it
# a single tested decision.
#
# CRITICAL basis distinction (do NOT "improve" this to use the dv_band tier):
# the V3 bar is REGIME-driven, not Axis-B-tier-driven. The dv_band is the SOURCED
# best-window maintenance ΔV (e.g. McConaghy's S1L1 ~10 m/s -> essentially_ballistic);
# the V3 input is the REPRODUCED continuous-from-one-seed horizon TCM (#169), a
# different quantity in a different basis (S1L1 reproduces at 62 m/s). 62 m/s is a
# valid V3-ballistic close (< 120) even though the row is sourced essentially_ballistic
# (< 10). Judging the 62 m/s reproduction against the <10 tier would WRONGLY demote a
# legitimately-reproduced row. So dv_band is used ONLY to identify the powered class.
V3_BALLISTIC_BUDGET_MPS: Final[float] = 120.0  # spec §14 generic V3-ballistic bar (7-cycle basis)
V3_POWERED_MARGIN: Final[float] = 1.10  # #175: continuous TCM <= documented budget x this


@dataclass(frozen=True)
class V3ClassVerdict:
    """Outcome of the §14/#175 V3 ballistic/powered class-split.

    Attributes
    ----------
    passed:
        ``True`` when the reproduced continuous TCM is within the class bar.
    cls:
        ``"ballistic"`` or ``"powered"`` — which bar was applied.
    bar_mps:
        The acceptance bar actually compared against, scaled to ``n_cycles``.
    basis:
        Human-readable note on which bar + why.
    """

    passed: bool
    cls: str
    bar_mps: float
    basis: str


def v3_class_split_verdict(
    measured_tcm_mps: float,
    *,
    n_cycles: int,
    dv_band: str | None = None,
    trajectory_regime: str | None = None,
    documented_budget_mps: float | None = None,
    powered_margin: float = V3_POWERED_MARGIN,
    ballistic_budget_mps: float = V3_BALLISTIC_BUDGET_MPS,
) -> V3ClassVerdict:
    """Judge a reproduced V3 continuous-TCM ΔV under the §14/#175 class-split.

    A row is POWERED iff ``dv_band == "powered_dsm"`` or
    ``trajectory_regime == "powered"``; then the bar is its own
    ``documented_budget_mps * powered_margin``. Otherwise it is BALLISTIC and the
    bar is ``ballistic_budget_mps`` (the generic 120 m/s, NOT the dv_band tier —
    see the module note above). Both bars are stated in the Russell 7-cycle basis
    and pro-rata scaled to ``n_cycles``.

    Reproduces the recorded #175 decisions: S1L1 (62 m/s, ballistic) PASS at < 120;
    russell-ch4-8.049gGf2 / App-C #188 (163.6 m/s, powered, budget 420) PASS; its
    sibling #192 (2040.6 m/s, powered, budget 1678) FAIL.
    """
    if measured_tcm_mps < 0.0:
        raise ValueError(f"measured_tcm_mps must be non-negative, got {measured_tcm_mps}")
    if n_cycles < 1:
        raise ValueError(f"n_cycles must be >= 1, got {n_cycles}")

    scale = n_cycles / RUSSELL_BASIS_CYCLES
    is_powered = dv_band == "powered_dsm" or trajectory_regime == "powered"
    if is_powered:
        if documented_budget_mps is None or documented_budget_mps <= 0.0:
            raise ValueError(
                "powered V3 class-split (#175) requires a positive documented_budget_mps "
                f"(the row's own published maintenance budget), got {documented_budget_mps!r}"
            )
        bar = documented_budget_mps * powered_margin * scale
        return V3ClassVerdict(
            passed=measured_tcm_mps <= bar,
            cls="powered",
            bar_mps=bar,
            basis=(
                f"V3-powered (#175): continuous TCM <= documented "
                f"{documented_budget_mps:g} m/s x {powered_margin:g}"
            ),
        )
    bar = ballistic_budget_mps * scale
    return V3ClassVerdict(
        passed=measured_tcm_mps <= bar,
        cls="ballistic",
        bar_mps=bar,
        basis=f"V3-ballistic (#175): continuous TCM <= {ballistic_budget_mps:g} m/s generic budget",
    )


__all__ = [
    "RUSSELL_BASIS_CYCLES",
    "V3_BALLISTIC_BUDGET_MPS",
    "V3_POWERED_MARGIN",
    "BandThreshold",
    "DvBandClassification",
    "V3ClassVerdict",
    "accept_maintenance_dv",
    "assign_dv_band_from_measurement",
    "classify_dv_band",
    "dv_band_threshold",
    "v3_class_split_verdict",
]
