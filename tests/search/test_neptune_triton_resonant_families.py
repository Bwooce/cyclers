"""Tests for the Neptune-Triton resonant/Lyapunov/DPO/LPO family module
(`#776`).

Reproduce-before-trust: sourced constants (mu, l*, t*, the ten ESM_GATE_ROWS)
are checked verbatim against the raw ESM text files' own printed values and
the JAS-2026 paper's/Miceli 2025 dissertation's own text layers (see the
module docstring). The gate is HONEST: all ten rows pass every one of the
three criteria (periodicity self-consistency, reproduction, internal
Barden-vs-planar_floquet cross-check) -- a genuinely clean sweep, unlike the
Saturn-Titan module's own two honest partial fails, explained by this
data's own already-nondimensional 12-decimal precision (no unit-conversion
loss). The two-body-seed-lineage attempt and one of the two 4:3-family
continuation branches are documented, well-characterized NEGATIVES, tested
here as standing regression evidence, not swept under the rug. See
``docs/notes/2026-08-01-776-neptune-triton-resonant-families-gate.md`` for
the full evidentiary writeup.
"""

from __future__ import annotations

import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.jovian_resonant_families as jrf
import cyclerfinder.search.neptune_triton_resonant_families as ntf


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return ntf.neptune_triton_system()


# ---------------------------------------------------------------------------
# (1) Sourced-constant reproduce-before-trust gate.
# ---------------------------------------------------------------------------


def test_mu_matches_every_esm_header_verbatim() -> None:
    """Every one of the three ESM text files' own header states, verbatim:
    'Mass ratio: 2.089503183689124e-04'."""
    assert ntf.MICELI_MU == 2.089503183689124e-04


def test_l_star_t_star_match_jas2026_p5_body_text() -> None:
    """JAS-2026 p.5, verbatim: 'l* = 354, 760 km' / 't* ~= 8.081353 x 10^4 s'."""
    assert ntf.MICELI_L_KM == 354760.0
    assert pytest.approx(8.081353e4) == ntf.MICELI_T_S


def test_esm_gate_rows_match_raw_esm_files_verbatim() -> None:
    """Ten rows, each (x0, ydot0, period, jacobi) copied verbatim from the
    raw ESM2/ESM3/ESM4 text files this task -- see module docstring for the
    per-row source_file/source_row_label provenance."""
    expected = {
        "1:7": (
            "ESM3",
            "Res17+x+h",
            7.011583577132,
            -6.902215957345,
            43.981049667607,
            1.806962818639,
        ),
        "3:2-start": (
            "ESM4",
            "Res32-x+h",
            -0.901451717937,
            -0.051900595659,
            12.554472960470,
            3.028835529717,
        ),
        "L2-lyapunov-target": (
            "ESM4",
            "L2LyapunovTarget",
            1.056460480396,
            -0.110092617115,
            4.110480863772,
            3.003706759619,
        ),
        "L1-lyapunov": (
            "ESM2",
            "L1Lyapunov",
            0.963347425337,
            -0.026829573549,
            2.968384422151,
            3.013995763057,
        ),
        "L2-lyapunov": (
            "ESM2",
            "L2Lyapunov",
            1.045065938781,
            -0.024387884231,
            3.139197034268,
            3.013770826075,
        ),
        "DPO": ("ESM4", "DPO", 1.018236356947, 0.094629509354, 1.522062471575, 3.013873926461),
        "LPO": ("ESM4", "LPO", 1.010469990495, 0.130296567487, 0.527889248185, 3.021659863278),
        "4:5-saddle": (
            "ESM4",
            "Res45+x+h",
            -0.969056422016,
            -0.126767119783,
            30.398418802755,
            2.987089791658,
        ),
        "4:7-stress": (
            "ESM4",
            "Res47-x+h",
            1.787757094694,
            -1.147924599743,
            44.514684550531,
            2.997230642137,
        ),
        "4:3-saddle": (
            "ESM4",
            "Res43+x+h",
            -0.597950460396,
            -0.828492541088,
            12.795970900869,
            3.016635194282,
        ),
    }
    assert set(ntf.ESM_GATE_ROWS) == set(expected)
    for label, (source_file, source_row_label, x0, ydot0, period, jacobi) in expected.items():
        row = ntf.ESM_GATE_ROWS[label]
        assert row.source_file == source_file, label
        assert row.source_row_label == source_row_label, label
        assert row.x0 == x0, label
        assert row.ydot0 == ydot0, label
        assert row.period == period, label
        assert row.jacobi == jacobi, label


def test_esm_gate_rows_span_all_three_source_files() -> None:
    files = {row.source_file for row in ntf.ESM_GATE_ROWS.values()}
    assert files == {"ESM2", "ESM3", "ESM4"}


def test_mu_differs_from_registry_by_small_known_amount() -> None:
    """Registry Neptune-Triton mu (GM_Triton / (GM_Neptune-system + GM_Triton))
    is NOT the paper's own value -- a documented ~2e-4 relative GM-vintage
    difference, the same class of delta the Jovian/Saturn-Titan modules
    document for their own systems."""
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    registry_mu = SATELLITES["Triton"].mu_km3_s2 / (
        PRIMARIES["Neptune"] + SATELLITES["Triton"].mu_km3_s2
    )
    rel = abs(registry_mu - ntf.MICELI_MU) / ntf.MICELI_MU
    assert 1e-5 < rel < 1e-2, f"unexpected mu delta {rel:.2e} (expected ~2e-4)"


def test_system_uses_source_mu_by_default(system: cr3bp.CR3BPSystem) -> None:
    assert system.mu == ntf.MICELI_MU
    assert system.primary == "Neptune"
    assert system.secondary == "Triton"


# ---------------------------------------------------------------------------
# (2)(a) Periodicity self-consistency: propagate each row's OWN raw IC for
# its OWN stated Integration Time, independent of this module's own
# corrector.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(ntf.ESM_GATE_ROWS))
def test_periodicity_self_consistency(label: str, system: cr3bp.CR3BPSystem) -> None:
    row = ntf.ESM_GATE_ROWS[label]
    resid = ntf.periodicity_residual(row, system)
    assert resid < ntf.PERIODICITY_GATE_TOL, f"{label}: periodicity residual {resid:.3e}"


# ---------------------------------------------------------------------------
# (2)(b) Reproduction: this module's own corrector, seeded at the row's own
# printed IC, must reproduce x0/ydot0/period.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(ntf.ESM_GATE_ROWS))
def test_sourced_seed_converges_and_reproduces(label: str, system: cr3bp.CR3BPSystem) -> None:
    row = ntf.ESM_GATE_ROWS[label]
    cand = ntf.recover_esm_candidate(label, system)
    assert cand.crossing_residual < 1e-9
    assert cand.jacobi == pytest.approx(row.jacobi, abs=1e-8)
    x0_err = abs(cand.x0 - row.x0) / abs(row.x0) if row.x0 != 0 else abs(cand.x0)
    ydot0_err = abs(cand.ydot0 - row.ydot0) / abs(row.ydot0)
    period_err = abs(cand.period - row.period) / abs(row.period)
    assert x0_err < ntf.REPRODUCTION_GATE_REL_TOL, f"{label}: x0_rel_err={x0_err:.2e}"
    assert ydot0_err < ntf.REPRODUCTION_GATE_REL_TOL, f"{label}: ydot0_rel_err={ydot0_err:.2e}"
    assert period_err < ntf.REPRODUCTION_GATE_REL_TOL, f"{label}: period_rel_err={period_err:.2e}"


# ---------------------------------------------------------------------------
# (2)(c) Barden vs planar_floquet -- INTERNAL cross-check only (no published
# eigenvalue target exists for this system).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(ntf.ESM_GATE_ROWS))
def test_barden_matches_planar_floquet_internal_crosscheck(
    label: str, system: cr3bp.CR3BPSystem
) -> None:
    cand = ntf.recover_esm_candidate(label, system)
    rel = abs(cand.max_eigenvalue - cand.planar_floquet_eigenvalue) / abs(cand.max_eigenvalue)
    assert rel < ntf.CROSSCHECK_GATE_REL_TOL, (
        f"{label}: barden={cand.max_eigenvalue} pf={cand.planar_floquet_eigenvalue}"
    )


def test_eigenvalue_magnitudes_match_771_scoping_survey(system: cr3bp.CR3BPSystem) -> None:
    """Cross-check against `#771`'s own scoping-pass eigenvalue survey
    (module docstring): 4:5-saddle ~105, 4:7-stress ~1.5e4, both genuine
    real saddles; the near-unit rows are genuinely NOT saddles."""
    saddle45 = ntf.recover_esm_candidate("4:5-saddle", system)
    assert saddle45.is_real_unstable
    assert saddle45.max_eigenvalue == pytest.approx(-105.05, abs=0.5)

    stress47 = ntf.recover_esm_candidate("4:7-stress", system)
    assert stress47.is_real_unstable
    assert stress47.max_eigenvalue == pytest.approx(14624.0, rel=1e-2)

    for label in ("1:7", "3:2-start", "LPO"):
        cand = ntf.recover_esm_candidate(label, system)
        assert not cand.is_real_unstable, label
        assert cand.max_eigenvalue == pytest.approx(1.0, abs=1e-6), label


# ---------------------------------------------------------------------------
# (3) The honest full gate itself -- reported row-by-row.
# ---------------------------------------------------------------------------


def test_gate_report_all_ten_rows_pass(system: cr3bp.CR3BPSystem) -> None:
    """A genuinely clean sweep: 10/10 rows pass all three criteria. Unlike
    the Saturn-Titan module's own two honest partial fails, explained by
    this data's own already-nondimensional 12-decimal precision (no
    unit-conversion loss like Vaquero's dimensional km/s table)."""
    rows = {r.label: r for r in ntf.gate_report(system)}
    assert set(rows) == set(ntf.ESM_GATE_ROWS)
    for label, row in rows.items():
        assert row.periodicity_confirmed, (label, row)
        assert row.reproduction_confirmed, (label, row)
        assert row.crosscheck_confirmed, (label, row)
        assert row.passed, (label, row)
    assert sum(r.passed for r in rows.values()) == 10


# ---------------------------------------------------------------------------
# (4) Equilibrium-point positive control (Miceli 2025 dissertation Table 2.1).
# ---------------------------------------------------------------------------


def test_dissertation_table21_equilibria_verbatim() -> None:
    assert ntf.DISSERTATION_TABLE21_EQUILIBRIA == {
        "L1": 0.959217,
        "L2": 1.041493,
        "L3": -1.000087,
        "L4": 0.499791,
        "L5": 0.499791,
    }


def test_equilibrium_gate_report_all_confirmed() -> None:
    rows = ntf.equilibrium_gate_report()
    assert {r.point for r in rows} == {"L1", "L2", "L3", "L4", "L5"}
    for row in rows:
        assert row.confirmed, row
        assert row.rel_err < 1e-6, row  # far inside EQUILIBRIUM_GATE_REL_TOL


# ---------------------------------------------------------------------------
# (5) Two-body-resonant-seed lineage: documented, well-characterized
# negative -- neither naive seed lands on its own labeled resonance.
# ---------------------------------------------------------------------------


def test_two_body_seed_4_3_lands_on_wrong_topology(system: cr3bp.CR3BPSystem) -> None:
    seed = jrf.two_body_resonant_seed(4, 3, x0_sign=-1)
    state0 = [seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0]
    import numpy as np

    c_nat = cr3bp.jacobi_constant(np.array(state0), system.mu)
    cand = jrf.converge_candidate(
        system, "4:3_2body", seed.x0, c_nat, seed.period_full, ydot0_sign=1.0, half_crossings=None
    )
    assert cand is not None
    # Lands near period/2pi ~= 4.0, NOT the seed's own 3.0 (4:3's own q).
    assert cand.period_over_2pi == pytest.approx(4.0, abs=0.05)


def test_two_body_seed_2_3_lands_on_wrong_topology(system: cr3bp.CR3BPSystem) -> None:
    seed = jrf.two_body_resonant_seed(2, 3, x0_sign=-1)
    state0 = [seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0]
    import numpy as np

    c_nat = cr3bp.jacobi_constant(np.array(state0), system.mu)
    cand = jrf.converge_candidate(
        system, "2:3_2body", seed.x0, c_nat, seed.period_full, ydot0_sign=-1.0, half_crossings=None
    )
    assert cand is not None
    # Lands near period/2pi ~= 7.0, NOT the seed's own 3.0 (2:3's own q).
    assert cand.period_over_2pi == pytest.approx(7.0, abs=0.05)


def test_two_body_seed_lineage_note_is_documented() -> None:
    assert "period/2pi ~= 4.0" in ntf.TWO_BODY_SEED_LINEAGE_NOTE
    assert "period/2pi ~= 7.0" in ntf.TWO_BODY_SEED_LINEAGE_NOTE


# ---------------------------------------------------------------------------
# (6) Continuation-in-C_J onto two multi-member families (item (d) of the
# `#776` gate). Two clean confirmations; one honest fold-reversal negative.
# ---------------------------------------------------------------------------


def test_continue_23_family_reaches_jacobi_bound_and_matches_printed_members(
    system: cr3bp.CR3BPSystem,
) -> None:
    branch = ntf.continue_23_family(system)
    assert branch.stop_reason is cc.StopReason.JACOBI_BOUND
    assert branch.n_rejected == 0
    assert len(branch.members) > 100
    for target in ntf.FAMILY_23[1:]:
        closest = min(branch.members, key=lambda m: abs(m.jacobi - target.jacobi))
        x0_err = abs(closest.x0 - target.x0) / abs(target.x0)
        t_err = abs(closest.period - target.period) / abs(target.period)
        assert x0_err < 1e-3, (target.jacobi, x0_err)
        assert t_err < 1e-3, (target.jacobi, t_err)


def test_continue_43_saddle_family_reaches_jacobi_bound_and_matches_printed_member(
    system: cr3bp.CR3BPSystem,
) -> None:
    branch = ntf.continue_43_saddle_family(system)
    assert branch.stop_reason is cc.StopReason.JACOBI_BOUND
    assert branch.n_rejected == 0
    assert len(branch.members) > 10
    target = ntf.FAMILY_43_HC2[-1]
    last = branch.members[-1]
    x0_err = abs(last.x0 - target.x0) / abs(target.x0)
    t_err = abs(last.period - target.period) / abs(target.period)
    assert x0_err < 1e-3, x0_err
    assert t_err < 1e-3, t_err


def test_continue_43_near_unit_family_hits_fold_reversal_honest_negative(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Documented negative (module docstring "family-mixing" finding): this
    branch does NOT cleanly reach its own target C -- an honest, well-
    characterized negative, not a bug to work around."""
    branch = ntf.continue_43_near_unit_family_fold_reversal(system)
    assert branch.stop_reason is cc.StopReason.FOLD_REVERSAL
    target = ntf.FAMILY_43_NEAR_UNIT[-1]
    assert branch.members[-1].jacobi < target.jacobi - 1e-4  # stopped well short


# ---------------------------------------------------------------------------
# (7) Connection stage explicitly out of scope for `#776` (its own Task-B
# analog, a later task) -- this module deliberately exposes no manifold/
# homoclinic/heteroclinic machinery.
# ---------------------------------------------------------------------------


def test_module_exposes_no_connection_machinery() -> None:
    connection_names = {"find_homoclinic", "find_heteroclinic", "correct_connection"}
    assert connection_names.isdisjoint(set(ntf.__all__))
