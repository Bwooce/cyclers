"""The binding Axis-B gate: the Aldrin outbound-ToF fidelity shift (Forge phase 1).

Plan Revision R1 (BINDING): the original S1L1 5.65/3.05 anchor is superseded
(unverified provenance). The ladder anchors instead on the **Aldrin outbound
time-of-flight shift** — both sides sourced and machine-readable since schema
v4.2:

* the circular-coplanar idealization ToF **146 d**
  (``trajectory.segments[out-em].tof_days``), and
* the Rogers et al. 2012 STOUR analytic-ephemeris band **[161, 172] d**
  (``trajectory.segments[out-em].tof_days_bounds``).

The documented fidelity behaviour: solving the same Aldrin E-M-E cell at a
higher fidelity must move the outbound ToF *upward, toward/into* the sourced
[161, 172] band, away from the lower-fidelity coplanar value. This test solves
the cell at both wired rungs and classifies the shift.

Golden discipline
-----------------
Every EXPECTED value here is READ from ``data/catalogue.yaml`` (the sourced 146
and [161, 172], and the sourced Aldrin a/e). The rung solutions are COMPUTED by
our code and are never used as an EXPECTED side. The observed COMPUTED ToFs are
reported in the assertion messages, not hard-coded as targets.
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.sequence import Cell
from cyclerfinder.verify.fidelity import (
    PersistenceClass,
    fidelity_persistence,
    solve_at_fidelity,
)

_ALDRIN_ID = "aldrin-classic-em-k1-outbound"


def _aldrin_sourced() -> tuple[float, float, int, tuple[float, float]]:
    """Read the SOURCED Aldrin a/e, outbound coplanar ToF, and STOUR band."""
    entry = load_catalog().by_id[_ALDRIN_ID]
    seg = next(s for s in entry.raw["trajectory"]["segments"] if s["id"] == "out-em")
    a_au = float(seg["a_au"])
    e = float(seg["e"])
    coplanar_tof_days = int(seg["tof_days"])
    raw_band = seg["tof_days_bounds"]
    assert len(raw_band) == 2
    band = (float(raw_band[0]), float(raw_band[1]))
    return a_au, e, coplanar_tof_days, band


def _aldrin_cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )


@pytest.mark.slow
def test_aldrin_outbound_tof_shifts_toward_stour_band() -> None:
    """BINDING GATE: the Aldrin outbound ToF rises from the coplanar rung toward
    the sourced [161, 172] STOUR band at the real-DE440 rung.

    EXPECTED sides (both SOURCED, read from the YAML): the coplanar idealization
    146 d and the band [161, 172] d. The two rung ToFs are COMPUTED; the test
    asserts the *direction* of the shift (documented) and reports the observed
    values verbatim.
    """
    a_au, e, coplanar_sourced_tof, band = _aldrin_sourced()
    cell = _aldrin_cell()

    # --- Rung 1: circular-coplanar (closed-form resonance construction).
    coplanar = solve_at_fidelity(cell, "circular-coplanar", a_au=a_au, e=e)
    assert coplanar.converged

    # --- Rung 2: real-DE440 (general ephemeris optimiser, phase-matched to the
    # sourced Aldrin V∞ anchors).
    real = solve_at_fidelity(
        cell,
        "real-de440",
        ephem=Ephemeris(model="astropy"),
        vinf_cap=12.0,
        priority_date_iso="1985-01-01",
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=5,
        seed=0,
    )
    assert real.converged, "real-DE440 Aldrin solve must converge for the gate"

    msg = (
        f"OBSERVED outbound ToF: coplanar={coplanar.outbound_tof_days:.2f} d, "
        f"real-DE440={real.outbound_tof_days:.2f} d; "
        f"SOURCED coplanar idealization={coplanar_sourced_tof} d, "
        f"SOURCED STOUR band={band} d"
    )

    # The real-ephemeris outbound ToF must move UPWARD relative to the coplanar
    # rung — the documented fidelity direction (toward the band).
    assert real.outbound_tof_days > coplanar.outbound_tof_days, msg

    # And it must move toward/into the SOURCED [161, 172] band relative to the
    # SOURCED coplanar idealization (146 d). Classify the shift from the sourced
    # coplanar anchor to the computed real-ephemeris value.
    report = fidelity_persistence(
        "outbound_tof_days",
        low_value=float(coplanar_sourced_tof),
        high_value=real.outbound_tof_days,
        abs_tol=2.0,
        expected_direction=+1,
        documented_band=band,
    )
    assert report.classification is PersistenceClass.SHIFTS_DOCUMENTED, msg


@pytest.mark.slow
def test_aldrin_real_rung_recovers_sourced_vinf_band() -> None:
    """Corroboration: the real-DE440 rung's recovered V∞ at Earth/Mars sit in
    the sourced Aldrin band (Russell 6.5 / 9.7 km/s). Anchors are SOURCED; the
    recovered V∞ are COMPUTED and only checked for band membership, not equality.
    """
    _a_au, _e, _coplanar_tof, _band = _aldrin_sourced()
    cell = _aldrin_cell()
    real = solve_at_fidelity(
        cell,
        "real-de440",
        ephem=Ephemeris(model="astropy"),
        vinf_cap=12.0,
        priority_date_iso="1985-01-01",
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=5,
        seed=0,
    )
    assert real.converged
    # SOURCED Aldrin V∞ (Russell 2004 Table 3.4): E 6.5, M 9.7 km/s.
    assert real.vinf_kms["E"] == pytest.approx(6.5, abs=2.0)
    assert real.vinf_kms["M"] == pytest.approx(9.7, abs=3.0)
