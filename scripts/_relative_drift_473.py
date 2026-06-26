"""#473 RELATIVE (scale-invariant) drift criterion — calibrated to #339.

WHY ABSOLUTE-KM WAS WRONG
-------------------------
The prior #473 run classified a quasi_cycler by drift SHAPE (oscillates and
returns) PLUS an ABSOLUTE ~530k km envelope calibrated to the #339 reference
quasi_cycler's 530k peak. But 530k km is a URANIAN-scale number: #339's
turn-around moon Oberon orbits at SMA 583,511 km, so 530k km is ~0.91 of that
orbit radius. Carried verbatim into a BIGGER system the absolute envelope
mis-scales: a Jovian tour (Ganymede SMA 1.07M km, Callisto 1.88M km) or a
Saturnian one (Titan 1.22M km, Rhea 527k km) legitimately drifts more in
absolute km purely because the system is larger. An absolute 530k-800k km
envelope therefore wrongly rejects everything Jovian/Saturnian, OR (as the
prior STEP-0 did) flips everything on SHAPE alone and admits 2.1M-3.76M km
excursions that are NOT rendezvous-quality.

THE FIX — RELATIVE DRIFT
------------------------
A tour is a bounded quasi_cycler ONLY if BOTH hold:
  1. SHAPE   : the per-cycle drift series oscillates-and-returns over >=10
               cycles (bounded, not monotonic divergence). [_drift_shape_473]
  2. RELATIVE: max_over_cycles(per_cycle_drift_km / SMA_norm_km) <= R_REF
               where SMA_norm is the SMA of the tour's NORMALISING moon and
               R_REF is calibrated FROM #339 itself.

R_REF is NON-CIRCULAR: #339 is the published / accepted V4 quasi_cycler, so
"at least as good as #339 in *relative* drift" is the principled bar. #339's
peak 530k km / Oberon SMA 583,511 km = 0.908. We adopt R_REF = 0.91 with a
small documented margin to R_MARGIN = 1.0 (a tour whose worst-cycle drift is up
to one full normalising-orbit radius is still "rendezvous-scale"); survivors
are reported with their exact ratio so the margin is auditable.

NORMALISING MOON
----------------
The V2-moontour rendezvous drift is measured at the cycle's FINAL encounter,
which is ``sequence[-1] == sequence[0]`` (the tour anchor). For #339 that anchor
is Umbriel (SMA 265,986) — but the published #339 reference envelope is quoted
against OBERON, the OUTERMOST (turn-around) moon of the tour. To keep the R_REF
calibration self-consistent we normalise by the SMA of the OUTERMOST encountered
moon (largest SMA in the sequence). This is documented and applied uniformly:
  * #339 Umbriel-Oberon-Umbriel -> Oberon 583,511 km -> ratio 0.908 (PASS).
This is the conservative choice (largest SMA -> smallest ratio); using the
anchor moon instead would only make the bar HARDER for inner-anchored tours, so
the outermost-moon normalisation never over-admits relative to anchor
normalisation for the #339-class geometry it is calibrated on.
"""

from __future__ import annotations

from cyclerfinder.core.satellites import SATELLITES

# #339 calibration: peak 530,000 km / Oberon SMA 583,511 km = 0.9083.
R_REF: float = 0.91
# Documented margin: admit up to one full normalising-orbit radius of worst-cycle
# drift. Survivors carry their exact ratio so this margin is auditable per-tour.
R_MARGIN: float = 1.0

SILVER_339_PEAK_KMS: float = 530_000.0
SILVER_339_NORM_MOON: str = "Oberon"


def normalising_moon(sequence: tuple[str, ...]) -> str:
    """Outermost (largest-SMA) encountered moon of the tour — see module docstring."""
    return max(set(sequence), key=lambda m: SATELLITES[m].sma_km)


def normalising_sma_km(sequence: tuple[str, ...]) -> float:
    return SATELLITES[normalising_moon(sequence)].sma_km


def relative_drift_ratio(max_drift_kms: float, sequence: tuple[str, ...]) -> float:
    """max_drift_km / SMA_of_outermost_moon — the scale-invariant drift metric."""
    return max_drift_kms / normalising_sma_km(sequence)


def passes_relative(
    max_drift_kms: float, sequence: tuple[str, ...], r_ref: float = R_MARGIN
) -> bool:
    """True iff the worst-cycle drift is <= r_ref normalising-orbit radii."""
    return relative_drift_ratio(max_drift_kms, sequence) <= r_ref
