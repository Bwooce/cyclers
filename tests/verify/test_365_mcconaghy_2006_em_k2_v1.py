"""§14 V1 frozen-gate evidence for #365 mcconaghy-2006-em-k2 (S1L1) honest-negative.

This module is the FROZEN-GATE wrapper for the McConaghy 2006 / 2004 S1L1
honest-negative — the single sink for the McConaghy 2004 JSR Table 6 (22-
encounter S1L1 DE405 itinerary) AND the eight McConaghy 2006 JSR Tables 2-9
itineraries (4 cycler-Δv-only + 4 cycler+taxi vehicles, all S1L1 family at
different launch epochs).

The §14 V1 like-for-like ask: does the closer_sweep_v1 substrate (single
circular-coplanar free-return ellipse) reproduce the row's abstract V_inf
anchor (E=4.7, M=5.0 km/s, transit 153 d) within 0.5 km/s AND clear the
V_inf-continuity gate?

VERDICT: the corrector closes with EMERGED V_inf E=4.771, M=5.036 km/s
(|ΔE|=0.07, |ΔM|=0.04, BOTH below the 0.5 km/s V1 floor) — but the §14 V1
mechanics gate REJECTS the row: V_inf-continuity break of 24 km/s at the
Mars flyby, the reconstructed return leg does not close to Earth. **The
S1L1 cycler is intrinsically MULTI-ARC** (McConaghy-Russell-Longuski 2005
Table 2 formal label: ``2g(2.8277, 657.97°, U) g(1.4508, 522.29°, L)``, two
generic Earth-Earth arcs joined at the Mars flyby) and a single radial-
crossing ellipse fundamentally cannot represent it.

HONEST NEGATIVE per ``feedback_orbit_closure_discipline``: V0 stays V0; the
multi-arc topology is the verdict, not a corrector deficiency. The row's
REAL elevated-tier evidence lives at the sibling ``russell-ch4-4.991gG2``
row (V3, real-eph DE440 closure with REBOUND/IAS15 cross-check per #167/
#94) — em-k2's V0 is by **publication gap at this row's metric** (a
circular-coplanar single-ellipse V_inf anchor), not by validation
infrastructure deficit.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
ROW_ID = "mcconaghy-2006-em-k2"
V1_FLOOR_KMS = 0.5


def _load_verdict() -> dict[str, Any]:
    path = REPO_ROOT / "data" / f"{ROW_ID}_v1_verdict.jsonl"
    assert path.exists(), f"frozen verdict file missing: {path}"
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            return json.loads(line)  # type: ignore[no-any-return]
    raise AssertionError(f"verdict JSONL empty: {path}")


def test_mcconaghy_2006_em_k2_v1_honest_negative() -> None:
    """em-k2 stays V0: V_inf matches sourced anchor within 0.5 km/s, BUT the
    §14 V1 mechanics gate rejects the row on V_inf-continuity (multi-arc).

    The single-ellipse vinf-continuity break IS the honest-negative verdict;
    this row's elevated-tier evidence lives at sibling row russell-ch4-4.991gG2
    (V3, real-eph closure), not via single-ellipse V1.
    """
    row = _load_verdict()
    assert row["candidate_id"] == ROW_ID
    assert row["passes_v1"] is False, (
        f"V1 verdict reports passes_v1={row['passes_v1']!r} — em-k2 is a "
        "documented honest-negative (multi-arc); registry would be invalidated."
    )
    assert row["verdict"] == "FAIL"
    assert "CLOSE-NOT-V1" in row["fail_mode"], (
        f"expected CLOSE-NOT-V1 fail_mode, got {row['fail_mode']!r}"
    )
    # Both V_inf legs ARE within the 0.5 km/s V1 floor — the FAIL is the
    # V_inf-continuity break, not the magnitude.
    assert row["vinf_gate_E_passed"] is True
    assert row["vinf_gate_M_passed"] is True
    # The discriminating fail is the V_inf-continuity gate (multi-arc).
    assert row["corrector"]["converged"] is True
    assert row["v1_mechanics_passed"] is False
    assert row["vinf_continuous"] is False
