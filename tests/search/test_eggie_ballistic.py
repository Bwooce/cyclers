"""#480 regression — true ballistic-cycler construction of the Hernandez 2017 EGGIE.

Two assertions, both deterministic (fixed converged seeds refined by the corrector):

* GATE A — a genuine *feasible* ballistic EGGIE exists in the ideal circular-coplanar
  model: equal-in/out |V∞| at all four flybys, the two Ganymede flybys at equal V∞,
  the cycle closed, every flyby altitude inside the paper's 25-70,000 km window, total
  flyby ΔV ~ 0. These are model-intrinsic PHYSICAL properties (not arbitrary goldens).
* GATE B — the construction reaches the *sourced* Table-4 V∞ levels (Europa 9.12, both
  Ganymede 7.07, Io 8.38 km/s) on the ballistic manifold, ballistically, but with
  SUB-SURFACE (infeasible) flyby altitudes. This documents the binding constraint of
  the strict 2D ideal model (the paper's gentle 653-6263 km altitudes require the
  real-ephemeris/3D conversion). The EXPECTED V∞ trace to the paper, never to our code
  (``feedback_golden_tests_sourced_only``).

Verdict note: ``docs/notes/2026-06-29-480-eggie-ballistic-construction-verdict.md``.
"""

from __future__ import annotations

from cyclerfinder.search.eggie_ballistic import (
    ALT_MAX_KM,
    ALT_MIN_KM,
    feasible_ballistic_eggie,
    interior_table4_eggie,
    table4_vinf_eggie,
)

# SOURCED from Hernandez-Jones-Jesick 2017, AAS 17-608, Table 4 (digest
# docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md).
# These EXPECTED values trace to the paper, never to our own computation.
TABLE4_VINF_EUROPA_KMS = 9.12
TABLE4_VINF_GANYMEDE_KMS = 7.07  # both Ganymede flybys, equal
TABLE4_VINF_IO_KMS = 8.38


def test_gate_a_feasible_ballistic_eggie_exists() -> None:
    """A closed, equal-V∞, bend-feasible ballistic EGGIE exists in the ideal model."""
    e = feasible_ballistic_eggie()

    # Ballistic: equal in/out |V∞| at every flyby (the defining property).
    assert e.ballistic_resnorm_kms < 1.0e-3, e.ballistic_resnorm_kms
    # The two Ganymede flybys occur at equal V∞ (resonant return).
    assert e.ganymede_equal_resid_kms < 1.0e-3, e.ganymede_equal_resid_kms
    # The cycle closes at the Europa seam (periodicity magnitude).
    assert e.seam_defect_kms < 1.0e-3, e.seam_defect_kms

    # Every flyby altitude lies inside the paper's 25-70,000 km window.
    assert e.all_feasible, e.flyby_alt_km
    for key, alt in e.flyby_alt_km.items():
        assert ALT_MIN_KM <= alt <= ALT_MAX_KM, (key, alt)

    # Near-ballistic total ΔV (the flybys carry essentially no maneuver).
    assert e.total_dv_ms < 1.0, e.total_dv_ms

    # It is a genuine EGGIE-topology member at lower excess speed than Table 4
    # (reported, not a Table-4 reproduction) — Ganymede well below 7.07.
    assert 6.0 < e.vinf_kms["Ganymede1"] < 7.0, e.vinf_kms


def test_gate_b_table4_vinf_reached_but_subsurface() -> None:
    """Table-4 V∞ is reproduced ballistically but with sub-surface (infeasible) flybys.

    The 2D ideal model puts the Table-4 V∞ levels exactly on the ballistic manifold
    (ΔV ~ 0, equal Ganymede, seam closed) yet forces ~180deg reversal bends ->
    sub-surface periapsides. The feasible gentle-bend altitudes the paper prints need
    the real-ephemeris conversion. This asserts that characterised result.
    """
    e = table4_vinf_eggie()

    # Ballistic and equal-Ganymede on the manifold.
    assert e.ballistic_resnorm_kms < 1.0e-3, e.ballistic_resnorm_kms
    assert e.ganymede_equal_resid_kms < 1.0e-3, e.ganymede_equal_resid_kms
    assert e.total_dv_ms < 1.0, e.total_dv_ms

    # V∞ match the SOURCED Table-4 levels (within 0.05 km/s).
    assert abs(e.vinf_kms["Europa_dep"] - TABLE4_VINF_EUROPA_KMS) < 0.05, e.vinf_kms
    assert abs(e.vinf_kms["Ganymede1_out"] - TABLE4_VINF_GANYMEDE_KMS) < 0.05, e.vinf_kms
    assert abs(e.vinf_kms["Ganymede2"] - TABLE4_VINF_GANYMEDE_KMS) < 0.05, e.vinf_kms
    assert abs(e.vinf_kms["Io_out"] - TABLE4_VINF_IO_KMS) < 0.05, e.vinf_kms

    # ...but the flybys are NOT bend-feasible: at least the Io flyby is sub-surface.
    assert not e.all_feasible
    assert e.flyby_alt_km["Io"] < ALT_MIN_KM, e.flyby_alt_km
    assert e.flyby_alt_km["G1"] < ALT_MIN_KM, e.flyby_alt_km


def test_gate_b_interior_subtour_reproduces_table4_feasibly() -> None:
    """The interior G->G->I sub-tour reproduces Table-4 V∞ with FEASIBLE altitudes.

    Pinpoints the binding constraint: with the Europa periodicity seam dropped, the 3
    interior flybys are exactly ballistic at the sourced Table-4 V∞ and ALL their
    altitudes lie in the 25-70000 km window; only the seam stays open (Europa arrival
    != departure, Europa flyby sub-surface). So the interior is feasible at Table-4 V∞
    in the 2D ideal model — it is full periodic closure (the seam) that fails.
    """
    e = interior_table4_eggie()

    # Interior flybys exactly ballistic (equal in/out |V∞|) and Ganymede equal.
    assert abs(e.vinf_kms["Ganymede1"] - e.vinf_kms["Ganymede1_out"]) < 1.0e-3, e.vinf_kms
    assert abs(e.vinf_kms["Ganymede2"] - e.vinf_kms["Ganymede2_out"]) < 1.0e-3, e.vinf_kms
    assert abs(e.vinf_kms["Io"] - e.vinf_kms["Io_out"]) < 1.0e-3, e.vinf_kms
    assert e.ganymede_equal_resid_kms < 1.0e-3, e.ganymede_equal_resid_kms

    # Interior V∞ match the SOURCED Table-4 levels.
    assert abs(e.vinf_kms["Europa_dep"] - TABLE4_VINF_EUROPA_KMS) < 0.05, e.vinf_kms
    assert abs(e.vinf_kms["Ganymede1_out"] - TABLE4_VINF_GANYMEDE_KMS) < 0.05, e.vinf_kms
    assert abs(e.vinf_kms["Io_out"] - TABLE4_VINF_IO_KMS) < 0.05, e.vinf_kms

    # The 3 interior flyby altitudes are all in-window.
    for key in ("G1", "G2", "Io"):
        assert ALT_MIN_KM <= e.flyby_alt_km[key] <= ALT_MAX_KM, (key, e.flyby_alt_km)

    # ...but the periodicity seam is open (the binding constraint).
    assert e.seam_defect_kms > 0.2, e.seam_defect_kms
    assert e.flyby_alt_km["E"] < ALT_MIN_KM, e.flyby_alt_km
