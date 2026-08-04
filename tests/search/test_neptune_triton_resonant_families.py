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

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_families as jrf
import cyclerfinder.search.neptune_triton_resonant_families as ntf


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return ntf.neptune_triton_system()


@pytest.fixture(scope="module")
def gate_777_rows(system: cr3bp.CR3BPSystem) -> list[ntf.GateRow]:
    """`#777`'s own full gate sweep over :data:`ntf.ESM_ROWS_777` (the 48
    rows `#776` did not vendor), computed ONCE per test session -- reused by
    every `#777` gate-row assertion below rather than re-converging all 48
    rows per test (the existing per-row parametrized tests above each
    independently re-converge all ten :data:`ntf.ESM_GATE_ROWS`; at 48 rows
    that pattern would mean hundreds of redundant corrector/STM runs)."""
    return ntf.gate_report_777(system)


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


# ---------------------------------------------------------------------------
# `#777`: vendor the remaining ~48 canonical ESM rows `#776` did not cover.
# Same honesty discipline: 47/48 pass all three criteria; one honest
# crosscheck-only negative, plus three honest continuation negatives
# alongside two clean continuation confirmations. See
# ``docs/notes/2026-08-04-777-neptune-triton-remaining-rows.md`` for the
# full evidentiary writeup.
# ---------------------------------------------------------------------------


def test_esm_rows_777_count_and_disjoint_from_esm_gate_rows() -> None:
    """48 new rows -- NOT `#776`'s own dispatch-note "~34" estimate, which
    this task's own full line-by-line audit of all three ESM files
    supersedes (see the `#777` results note)."""
    assert len(ntf.ESM_ROWS_777) == 48
    assert set(ntf.ESM_ROWS_777).isdisjoint(set(ntf.ESM_GATE_ROWS))


def test_esm_rows_777_span_both_esm3_and_esm4() -> None:
    files = {row.source_file for row in ntf.ESM_ROWS_777.values()}
    assert files == {"ESM3", "ESM4"}


def test_esm_rows_777_sample_values_match_raw_esm_files_verbatim() -> None:
    """Spot-check (not all 48, for readability) against the raw ESM3/ESM4
    text files this task -- one row per interesting case: a large-|x0|
    high-energy ESM3 row, the lowest-C DPO member, the mixed-half_crossings
    "3:2-x-esm4-4-hc2" outlier, and the topology-outlier "2:1+x-esm4-1"."""
    expected = {
        "1:2+x-esm3": (
            "ESM3",
            "Res12+x+h",
            -2.881985569172,
            2.629061067518,
            12.565322956698,
            2.087857684887,
        ),
        "dpo-esm4-8": (
            "ESM4",
            "DPO",
            1.004676477572,
            0.289504961848,
            7.065417702085,
            3.000962711997,
        ),
        "3:2-x-esm4-4-hc2": (
            "ESM4",
            "Res32-x+h",
            -0.925033394133,
            0.001786689533,
            12.557846095824,
            3.018021426630,
        ),
        "2:1+x-esm4-1": (
            "ESM4",
            "Res21+x+h",
            -0.422381198297,
            -1.116798990950,
            2.375767706138,
            3.667872679873,
        ),
    }
    for key, (source_file, source_row_label, x0, ydot0, period, jacobi) in expected.items():
        row = ntf.ESM_ROWS_777[key]
        assert row.source_file == source_file, key
        assert row.source_row_label == source_row_label, key
        assert row.x0 == x0, key
        assert row.ydot0 == ydot0, key
        assert row.period == period, key
        assert row.jacobi == jacobi, key


def test_half_crossings_none_auto_detection_matches_esm_gate_rows_hand_picked_values(
    system: cr3bp.CR3BPSystem,
) -> None:
    """The assumption `#777`'s own automatic ``half_crossings`` determination
    for all 48 :data:`ntf.ESM_ROWS_777` rows rests on: the SAME logic
    :func:`cp.correct_symmetric_fixed_jacobi` uses internally when
    ``half_crossings=None`` (the x-axis crossing nearest ``T/2`` on the raw
    seed) recovers the EXACT integer `#776` hand-picked for every one of
    the ten :data:`ntf.ESM_GATE_ROWS`."""
    for label, row in ntf.ESM_GATE_ROWS.items():
        orbit = cp.correct_symmetric_fixed_jacobi(
            system,
            row.x0,
            row.jacobi,
            row.period,
            ydot0_sign=row.ydot0_sign,
            half_crossings=None,
            tol=1e-11,
            rtol=1e-13,
            atol=1e-13,
            x0_bounds=row.x0_bounds,
        )
        assert orbit.converged, label
        x0_err = abs(orbit.x0 - row.x0) / abs(row.x0) if row.x0 != 0 else abs(orbit.x0)
        assert x0_err < 1e-4, (label, x0_err)


def test_gate_report_777_has_48_rows(gate_777_rows: list[ntf.GateRow]) -> None:
    assert len(gate_777_rows) == 48
    assert {r.label for r in gate_777_rows} == set(ntf.ESM_ROWS_777)


def test_gate_report_777_periodicity_and_reproduction_all_pass(
    gate_777_rows: list[ntf.GateRow],
) -> None:
    """(a)/(b) are a clean 48/48 -- every row's own raw seed self-closes and
    every row's own corrected candidate reproduces its printed x0/ydot0/
    period. Only (c), the internal crosscheck, has one honest near-miss
    (see the two tests below)."""
    for row in gate_777_rows:
        assert row.periodicity_confirmed, row
        assert row.reproduction_confirmed, row


def test_gate_report_777_47_of_48_pass_all_three_criteria(
    gate_777_rows: list[ntf.GateRow],
) -> None:
    n_pass = sum(r.passed for r in gate_777_rows)
    assert n_pass == 47


def test_gate_report_777_single_honest_crosscheck_negative(
    gate_777_rows: list[ntf.GateRow],
) -> None:
    """The one honest negative: ``"dpo-esm4-2"`` (ESM4 DPO, 2nd printed
    member) reproduces cleanly on (a)/(b) but its Barden-vs-planar_floquet
    internal cross-check misses :data:`ntf.CROSSCHECK_GATE_REL_TOL` by
    ~2x -- a genuine near-unit-circle numerical-sensitivity near-miss on a
    very short-period (T~1.63 nondim) orbit, `#777`'s own analog of
    `#776`'s own 4:3 fold-reversal finding: a real, well-characterized
    negative, not a bug."""
    failing = [r for r in gate_777_rows if not r.passed]
    assert len(failing) == 1
    bad = failing[0]
    assert bad.label == "dpo-esm4-2"
    assert bad.periodicity_confirmed
    assert bad.reproduction_confirmed
    assert not bad.crosscheck_confirmed
    assert bad.crosscheck_rel_err == pytest.approx(2.145766139904934e-05, rel=1e-3)
    assert bad.crosscheck_rel_err > ntf.CROSSCHECK_GATE_REL_TOL


# ---------------------------------------------------------------------------
# `#777`'s own two-body-seed-lineage extension: (1,2) and (2,1), two
# resonance ratios `#776` never tried.
# ---------------------------------------------------------------------------


def test_two_body_seed_1_2_period_index_hits_but_wrong_orbit(system: cr3bp.CR3BPSystem) -> None:
    """Converges to the RIGHT period index (q=2) -- unlike every other
    naive attempt in this family of checks -- but at a hugely different C
    and x0 from the paper's own printed "1:2+x-esm3" row: still not a
    genuine identification of the paper's own labeled family member."""
    seed = jrf.two_body_resonant_seed(1, 2, x0_sign=-1)
    state0 = [seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0]
    c_nat = cr3bp.jacobi_constant(np.array(state0), system.mu)
    cand = jrf.converge_candidate(
        system, "1:2_2body", seed.x0, c_nat, seed.period_full, ydot0_sign=-1.0, half_crossings=None
    )
    assert cand is not None
    assert cand.period_over_2pi == pytest.approx(2.0, abs=0.05)
    row = ntf.ESM_ROWS_777["1:2+x-esm3"]
    assert abs(cand.jacobi - row.jacobi) > 0.5  # same period index, very different orbit


def test_two_body_seed_2_1_lands_on_wrong_topology(system: cr3bp.CR3BPSystem) -> None:
    seed = jrf.two_body_resonant_seed(2, 1, x0_sign=-1)
    state0 = [seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0]
    c_nat = cr3bp.jacobi_constant(np.array(state0), system.mu)
    cand = jrf.converge_candidate(
        system, "2:1_2body", seed.x0, c_nat, seed.period_full, ydot0_sign=1.0, half_crossings=None
    )
    assert cand is not None
    # Lands near period/2pi ~= 2.0, NOT the seed's own 1.0 (2:1's own q).
    assert cand.period_over_2pi == pytest.approx(2.0, abs=0.05)


def test_two_body_seed_lineage_note_777_is_documented() -> None:
    assert "period/2pi ~= 2.0" in ntf.TWO_BODY_SEED_LINEAGE_NOTE_777
    assert "1:2" in ntf.TWO_BODY_SEED_LINEAGE_NOTE_777
    assert "2:1" in ntf.TWO_BODY_SEED_LINEAGE_NOTE_777


# ---------------------------------------------------------------------------
# `#777`'s own continuation-in-C_J checks (item (e) of its gate): two clean
# confirmations, three honest negatives.
# ---------------------------------------------------------------------------


def test_continue_32_esm4_hc1_family_777_reaches_jacobi_bound(system: cr3bp.CR3BPSystem) -> None:
    branch = ntf.continue_32_esm4_hc1_family_777(system)
    assert branch.stop_reason is cc.StopReason.JACOBI_BOUND
    assert branch.n_rejected == 0
    assert len(branch.members) > 100
    for target in ntf.FAMILY_32_ESM4_HC1_777[1:]:
        closest = min(branch.members, key=lambda m: abs(m.jacobi - target.jacobi))
        x0_err = abs(closest.x0 - target.x0) / abs(target.x0)
        t_err = abs(closest.period - target.period) / abs(target.period)
        assert x0_err < 1e-3, (target.jacobi, x0_err)
        assert t_err < 1e-3, (target.jacobi, t_err)


def test_continue_47_esm4_pair_777_reaches_jacobi_bound(system: cr3bp.CR3BPSystem) -> None:
    branch = ntf.continue_47_esm4_pair_777(system)
    assert branch.stop_reason is cc.StopReason.JACOBI_BOUND
    assert branch.n_rejected == 0
    assert len(branch.members) > 50
    target = ntf.FAMILY_47_ESM4_HC3_777[-1]
    last = branch.members[-1]
    x0_err = abs(last.x0 - target.x0) / abs(target.x0)
    t_err = abs(last.period - target.period) / abs(target.period)
    assert x0_err < 1e-3, x0_err
    assert t_err < 1e-3, t_err


def test_continue_dpo_family_777_hits_gauntlet_reject_honest_negative(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Documented negative (module docstring): unlike a fold or a topology
    jump, the DPO family walk is stopped by the gauntlet's own physical-
    plausibility rejection partway through -- a genuine, well-characterized
    negative, not a bug to work around."""
    branch = ntf.continue_dpo_family_gauntlet_reject_777(system)
    assert branch.stop_reason is cc.StopReason.GAUNTLET_REJECT
    assert branch.n_rejected >= 1
    assert len(branch.members) < len(ntf.FAMILY_DPO_777) * 10  # stopped well short of a full sweep


def test_continue_21_esm4_family_777_hits_fold_reversal_honest_negative(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Documented negative (module docstring): the walk folds back before
    reaching the 6th printed member, a genuine outlier ~0.5 higher in C
    than its own 5 siblings."""
    branch = ntf.continue_21_esm4_family_fold_reversal_777(system)
    assert branch.stop_reason is cc.StopReason.FOLD_REVERSAL
    target = ntf.FAMILY_21_ESM4_777[-1]
    assert branch.members[-1].jacobi < target.jacobi - 1e-2  # stopped well short


def test_continue_25m_esm4_family_777_hits_topology_jump_honest_negative(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Documented negative (module docstring): the two printed "Res25-x+h"
    rows are NOT two points on one continuous branch -- the walk fails at
    (or before) the very first continuation step."""
    branch = ntf.continue_25m_esm4_family_topology_jump_777(system)
    assert branch.stop_reason is cc.StopReason.TOPOLOGY_JUMP
    assert len(branch.members) <= 1
