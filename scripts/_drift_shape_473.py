"""#473 drift-series SHAPE classifier (shared by posctl / #470 re-adjudication / search).

The V2-moontour gauntlet's strict ``passes_v2`` only compares max drift to the
50,000 km floor. A QUASI-cycler (the #339 class) deliberately fails that floor
but is still a valid catalogue row because its per-cycle drift OSCILLATES within
a bounded envelope and RETURNS toward its low band rather than diverging
monotonically (see tests/verify/test_silver_327_v2_quasi_cycler.py).

This module classifies a per-cycle drift series (cycle >= 1, the seed cycle 0 is
0 by construction) by SHAPE so a search can score bounded-vs-divergent directly
rather than only against the strict floor. It is calibrated against the #339
positive control: drift 0 -> 326k -> 515k -> 487k -> 254k -> 86k -> 390k -> 530k
-> 447k -> 176k must classify as ``bounded-oscillating-and-returns``.
"""

from __future__ import annotations

import itertools


def classify_drift_shape(drifts_after_seed: list[float]) -> dict:
    """Classify a per-cycle drift series (cycle >= 1) by SHAPE.

    Returns a dict with ``shape`` in:
      * ``monotonic-divergence`` — nearly monotone rise, ends at the peak,
        never returns to its low band (genuine reject).
      * ``bounded-oscillating-and-returns`` — rises to a peak then revisits its
        low band, ending off-peak (the #339 quasi-cycler signature).
      * ``bounded-oscillating`` — returns to its low band but ends near the peak.
      * ``bounded-flat`` — max within 1.3x of min (essentially stable).
      * ``bounded-non-returning`` — spread but never revisits the low band.
    """
    if not drifts_after_seed:
        return {"shape": "no-data", "max": None, "min": None}
    mx = max(drifts_after_seed)
    mn = min(drifts_after_seed)
    last = drifts_after_seed[-1]
    n = len(drifts_after_seed)
    peak_idx = drifts_after_seed.index(mx)
    rises = sum(1 for a, b in itertools.pairwise(drifts_after_seed) if b > a)
    rise_frac = rises / max(1, n - 1)
    # "returns" = after rising above 2x the minimum, the series comes back down
    # to within 1.3x of its minimum at least once (a genuine peak-then-return,
    # anywhere in the horizon — not necessarily after the single global max).
    returns = False
    risen = False
    for d in drifts_after_seed:
        if d >= 2.0 * mn:
            risen = True
        elif risen and d <= 1.3 * mn:
            returns = True
            break
    ends_at_peak = last >= 0.95 * mx
    divergent = rise_frac >= 0.85 and ends_at_peak and not returns
    spread = mx > 1.3 * mn
    if divergent:
        shape = "monotonic-divergence"
    elif not spread:
        shape = "bounded-flat"
    elif returns and not ends_at_peak:
        shape = "bounded-oscillating-and-returns"
    elif returns:
        shape = "bounded-oscillating"
    else:
        shape = "bounded-non-returning"
    return {
        "shape": shape,
        "max": mx,
        "min": mn,
        "last": last,
        "peak_cycle": peak_idx + 1,
        "rise_frac": round(rise_frac, 3),
        "returns": returns,
        "ends_at_peak": ends_at_peak,
    }


# Shapes a quasi-cycler may legitimately take (bounded, non-divergent).
BOUNDED_SHAPES = frozenset(
    {"bounded-oscillating-and-returns", "bounded-oscillating", "bounded-flat"}
)
