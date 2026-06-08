"""Tests for the GMAT V4 report parser + two-part V4 predicate (#171).

Fixture-string only — NO GMAT is installed or invoked. The reference ΔVs
(Aldrin 2.9138 km/s; S1L1 62 m/s) are OUR values under external check; the ±5%
band is declared up front. The convergence bar (Jones 1e-3 km / 1e-6 km/s) is the
SOURCED tolerance.
"""

from __future__ import annotations

from scripts.parse_gmat_report import (
    parse_convergence,
    parse_maintenance_dv,
    parse_report,
    v4_pass,
)

# --- Phase 3.1: maintenance ΔV sums per-flyby TCMs ----------------------------

_REPORT_CONVERGED_MULTI = """\
GMAT Mission Run
Target FlybyTCM_Mars3: targeting converged
TCM_Mars3.Magnitude = 0.0210
Target FlybyTCM_Mars6: targeting converged
TCM_Mars6.Magnitude = 0.0085
Target FlybyTCM_Mars9: targeting converged
TCM_Mars9.Magnitude = 0.0325
*** GMAT run complete ***
"""

_REPORT_NOT_CONVERGED = """\
GMAT Mission Run
Target FlybyTCM_Mars3: targeting converged
TCM_Mars3.Magnitude = 0.0210
Target FlybyTCM_Mars6: Targeter did not converge
*** GMAT run complete ***
"""

_REPORT_ALDRIN_OK = """\
GMAT Mission Run
Target FlybyTCM_MarsReturn: optimization converged
TCM_MarsReturn.TotalDV = 2.9500
*** GMAT run complete ***
"""


def test_parse_maintenance_dv_sums_tcm() -> None:
    dv = parse_maintenance_dv(_REPORT_CONVERGED_MULTI)
    assert dv == 0.0210 + 0.0085 + 0.0325


def test_parse_report_extracts_per_flyby() -> None:
    parsed = parse_report(_REPORT_CONVERGED_MULTI)
    assert parsed.tcm_magnitudes_kms == (0.0210, 0.0085, 0.0325)
    assert parsed.converged is True
    assert parsed.n_failed_blocks == 0
    assert parsed.maintenance_dv_kms == 0.0210 + 0.0085 + 0.0325


# --- Phase 3.2: convergence predicate ----------------------------------------


def test_parse_convergence_true_when_all_converged() -> None:
    assert parse_convergence(_REPORT_CONVERGED_MULTI) is True


def test_parse_convergence_false_on_failure_line() -> None:
    assert parse_convergence(_REPORT_NOT_CONVERGED) is False


def test_parse_convergence_false_when_no_evidence() -> None:
    assert parse_convergence("GMAT Mission Run\nno targeting blocks here\n") is False


# --- Phase 3.3: two-part predicate (band + convergence) ----------------------


def test_v4_pass_within_band_and_converged() -> None:
    # Aldrin: 2.95 within 5% of 2.9138 (band 2.768-3.060), converged.
    assert v4_pass(2.95, ref_dv=2.9138, converged=True) is True


def test_v4_pass_fails_when_not_converged() -> None:
    assert v4_pass(2.95, ref_dv=2.9138, converged=False) is False


def test_v4_pass_fails_outside_band() -> None:
    # 3.20 is outside the 2.768-3.060 band even though converged.
    assert v4_pass(3.20, ref_dv=2.9138, converged=True) is False


def test_v4_band_edges() -> None:
    ref = 2.9138
    band = 0.05 * ref
    # Just inside the band passes; just outside fails (FP-robust margins).
    assert v4_pass(ref - band + 1e-9, ref_dv=ref, converged=True) is True
    assert v4_pass(ref + band - 1e-9, ref_dv=ref, converged=True) is True
    assert v4_pass(ref - band - 1e-6, ref_dv=ref, converged=True) is False
    assert v4_pass(ref + band + 1e-6, ref_dv=ref, converged=True) is False


# --- Phase 3.4: convergence-only fallback (no reference) ----------------------


def test_v4_pass_convergence_only_when_no_reference() -> None:
    # The Mars-perturbed continuous arm: no prior reference dV.
    assert v4_pass(0.137, ref_dv=None, converged=True) is True
    assert v4_pass(0.137, ref_dv=None, converged=False) is False


_REPORT_COLUMNAR = """\
Sat.A1ModJulian           dv_Mars3
31393.638                 0.3340549804286069
"""

_REPORT_COLUMNAR_MULTI = """\
*** GMAT Mission Run
*** The Targeter converged!
Sat.A1ModJulian           dv_Mars3
31393.638                 0.3340549804286069
*** The Targeter converged!
Sat.A1ModJulian           dv_Mars6
31393.648                 1.273223947801619
"""


def test_parse_columnar_reportfile() -> None:
    """The real GMAT ReportFile form (header + numeric data row)."""
    parsed = parse_report(_REPORT_COLUMNAR)
    assert parsed.tcm_magnitudes_kms == (0.3340549804286069,)
    assert parsed.maintenance_dv_kms == 0.3340549804286069


def test_parse_columnar_multiblock_sums_all() -> None:
    """Concatenated per-flyby reports (the combined-log case): sum all dv columns."""
    parsed = parse_report(_REPORT_COLUMNAR_MULTI)
    assert parsed.maintenance_dv_kms == 0.3340549804286069 + 1.273223947801619
    assert parsed.converged is True


def test_aldrin_report_end_to_end() -> None:
    parsed = parse_report(_REPORT_ALDRIN_OK)
    assert parsed.converged is True
    assert parsed.maintenance_dv_kms == 2.9500
    assert v4_pass(parsed.maintenance_dv_kms, ref_dv=2.9138, converged=parsed.converged) is True
