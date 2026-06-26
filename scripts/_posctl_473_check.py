"""#473 POSITIVE-CONTROL gate — reproduce the #339 quasi-cycler bounded signature.

This is the HARD GATE before any search. It reconstructs the #339 known-good
Umbriel-Oberon-Umbriel quasi-cycler with its REAL geometry (NOT the killed
run's wrong 2.993-d skeleton) and runs the OFFICIAL run_v2_moontour at
n_cycles=10. It then classifies the per-cycle drift series by SHAPE and asserts
the bounded-oscillating-and-returns signature documented in
tests/verify/test_silver_327_v2_quasi_cycler.py:

  * 10/10 cycles complete (every Lambert leg converges)
  * drift oscillates within ~86k-530k km and RETURNS toward the ~86k floor
  * NOT monotonic divergence
  * strict passes_v2 is FALSE (it is a QUASI-cycler, meant to fail the 50k floor)

If this does not reproduce, the apparatus is wrong and the search must NOT run.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from cyclerfinder.data.validation.v2_moontour import run_v2_moontour

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _drift_shape_473 import classify_drift_shape

# REAL #339 geometry — verified in tests/data/test_v2_moontour.py,
# tests/data/test_v4_uranus.py, catalogue row
# umbriel-oberon-1-1-uranian-quasi-cycler-2026.
SILVER_ID = "repeated-moon-uranus-00000041"
SILVER_SEQ = ("Umbriel", "Oberon", "Umbriel")
SILVER_VINF = (0.9199258810725036, 0.9604309791298091, 0.8946936085078939)
SILVER_TOF = (14.940560615336594, 14.940560615336594)  # days per leg — NOT 2.993
SILVER_REL_OFF = 180.0
SILVER_NREV = (1, 1)
SILVER_PHASE0 = 29.999999999999996


def main() -> None:
    verdict = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=10,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    per_cycle_drifts = [float(c.rendezvous_drift_kms) for c in verdict.per_cycle]
    drifts_after_seed = per_cycle_drifts[1:]
    cls = classify_drift_shape(drifts_after_seed)

    out = {
        "candidate_id": verdict.candidate_id,
        "sequence": list(verdict.sequence),
        "leg_tofs_days": list(SILVER_TOF),
        "n_cycles_requested": verdict.n_cycles_requested,
        "n_cycles_completed": verdict.n_cycles_completed,
        "passes_v2_strict": verdict.passes_v2,
        "max_drift_kms": round(verdict.max_drift_kms, 1),
        "per_cycle_drift_kms": [round(d, 1) for d in per_cycle_drifts],
        "shape_classification": cls,
    }
    print(json.dumps(out, indent=2, ensure_ascii=True))

    # Gate assertions: the documented #339 bounded signature.
    assert verdict.n_cycles_completed == 10, "POSCTL must complete 10/10 cycles"
    assert verdict.passes_v2 is False, "POSCTL is a QUASI-cycler; strict V2 must FAIL"
    assert cls["shape"] == "bounded-oscillating-and-returns", (
        f"POSCTL drift shape is {cls['shape']!r}, expected bounded-oscillating-and-returns"
    )
    assert 60_000.0 < cls["min"] < 120_000.0, f"POSCTL min drift {cls['min']} off the ~86k floor"
    assert cls["max"] < 600_000.0, f"POSCTL max drift {cls['max']} exceeds the ~530k band"
    print("\nPOSITIVE CONTROL: PASS (bounded-oscillating-and-returns reproduced).")


if __name__ == "__main__":
    main()
