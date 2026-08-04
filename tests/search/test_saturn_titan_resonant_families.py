"""Tests for the Saturn-Titan resonant-orbit family module (#765).

Reproduce-before-trust: sourced constants (mu, C, Table 4.1 targets) are
checked verbatim against Vaquero 2013's own stated values (verified directly
against the PDF text layer AND a rendered page-image read of Table 4.1,
p.109 -- see the module docstring). The gate is HONEST: 3:4 and L1 pass
outright; 6:5 fails on eigenvalue only (a genuine, small, well-characterized
near-miss); L2 passes on eigenvalue but fails on period (a documented,
likely source-errata period value) -- none of this is fudged. See
``docs/notes/2026-07-29-765-saturn-titan-resonant-families-vaquero-gate.md``
for the full evidentiary writeup.
"""

from __future__ import annotations

import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_families as jrf
import cyclerfinder.search.saturn_titan_resonant_families as stf


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return stf.saturn_titan_system()


# ---------------------------------------------------------------------------
# (1) Sourced-constant reproduce-before-trust gate.
# ---------------------------------------------------------------------------


def test_mu_matches_thesis_p132() -> None:
    """Vaquero 2013 p.132: mu ~ 2.3658e-4, verbatim (her own text uses '~',
    an explicit approximate-display signal -- see module docstring)."""
    assert stf.VAQUERO_MU == 2.3658e-4


def test_c_matches_table41_caption() -> None:
    """Table 4.1 caption (p.109): C = 3.010000, verbatim."""
    assert stf.VAQUERO_C == 3.010000


def test_table41_targets_match_p109_verbatim() -> None:
    assert stf.TABLE41_TARGETS == {
        "3:4": 2129.81,
        "6:5": 191.641,
        "L1": 1004.72,
        "L2": 892.850,
    }


def test_table41_dimensional_matches_p109_verbatim() -> None:
    assert stf.TABLE41_DIMENSIONAL == {
        "3:4": (1.25869e6, 0.477301, 66.3312),
        "6:5": (1.14214e6, 0.545759, 71.2638),
        "L1": (1.15897e6, 0.447315, 8.2829),
        "L2": (1.25231e6, 0.549329, 79.7260),
    }


def test_mu_differs_from_registry_by_small_known_amount() -> None:
    """Registry Saturn-Titan mu (GM_Titan / (GM_Saturn-system + GM_Titan))
    is NOT the thesis's own rounded display value -- a documented ~0.03%
    GM-vintage difference, the same class of delta the Jovian module
    documents for its own Jupiter-Europa gap."""
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    registry_mu = SATELLITES["Titan"].mu_km3_s2 / (
        PRIMARIES["Saturn"] + SATELLITES["Titan"].mu_km3_s2
    )
    rel = abs(registry_mu - stf.VAQUERO_MU) / stf.VAQUERO_MU
    assert 1e-5 < rel < 1e-2, f"unexpected mu delta {rel:.2e} (expected ~3e-4)"


def test_system_uses_thesis_mu_by_default(system: cr3bp.CR3BPSystem) -> None:
    assert system.mu == stf.VAQUERO_MU
    assert system.primary == "Saturn"
    assert system.secondary == "Titan"


# ---------------------------------------------------------------------------
# (2) l*/t* self-consistency: nondimensionalizing Table 4.1's own printed
# (x, ydot) with this module's derived l*/t* and the thesis's own mu must
# reproduce the STATED Jacobi constant, C=3.010000 -- independent evidence
# the l*/t* choice matches whatever the thesis used internally.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(stf.TABLE41_DIMENSIONAL))
def test_table41_ic_reproduces_stated_jacobi_constant(
    label: str, system: cr3bp.CR3BPSystem
) -> None:
    import numpy as np

    x_km, ydot_kms, _t_days = stf.TABLE41_DIMENSIONAL[label]
    x0 = x_km / system.l_km
    vel_star = system.l_km / system.t_s
    ydot0 = ydot_kms / vel_star
    state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    c_computed = cr3bp.jacobi_constant(state0, system.mu)
    rel = abs(c_computed - stf.VAQUERO_C) / stf.VAQUERO_C
    assert rel < 2e-5, f"{label}: computed C={c_computed}, rel_err={rel:.2e}"


def test_naive_two_body_seed_does_not_converge_at_vaquero_c(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Matches Anderson & Lo's own documented finding for the analogous
    Jupiter-Europa attempt (#753 module docstring item 1): the literal
    periapsis-at-secondary two-body seed does not even produce a valid
    ydot0 at this Jacobi constant -- an expected, DOCUMENTED negative, not
    a bug. Confirms this module leads with the sourced Table 4.1 seeds
    rather than the blind two-body construction, per the #765 dispatch
    note's own instruction."""
    seed = jrf.two_body_resonant_seed(3, 4, x0_sign=-1)
    with pytest.raises(ValueError, match="negative radicand"):
        import cyclerfinder.search.cr3bp_periodic as cp

        cp.ydot0_from_jacobi(seed.x0, stf.VAQUERO_C, system.mu, sign=1.0)


# ---------------------------------------------------------------------------
# (3) Sourced seeds: standing regression (still converge, tight residual).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(stf.TABLE41_TARGETS))
def test_sourced_seed_converges(label: str, system: cr3bp.CR3BPSystem) -> None:
    cand = stf.recover_table41_candidate(label, system)
    assert cand.crossing_residual < 1e-9
    assert cand.jacobi == pytest.approx(stf.VAQUERO_C, abs=1e-8)
    assert cand.is_real_unstable  # every Table 4.1 row is a real saddle, not
    # a complex unit-modulus pair (see module docstring finding 4)


@pytest.mark.parametrize("label", sorted(stf.TABLE41_TARGETS))
def test_barden_matches_planar_floquet(label: str, system: cr3bp.CR3BPSystem) -> None:
    """Independent cross-check (feedback_orbit_closure_discipline): Barden's
    half-period identity and a direct full-period monodromy eigendecomposition
    must agree closely for a trustworthy classification."""
    cand = stf.recover_table41_candidate(label, system)
    rel = abs(cand.max_eigenvalue - cand.planar_floquet_eigenvalue) / cand.max_eigenvalue
    assert rel < 1e-6, f"{label}: barden={cand.max_eigenvalue} pf={cand.planar_floquet_eigenvalue}"


@pytest.mark.parametrize("label", sorted(stf.TABLE41_TARGETS))
def test_sourced_seed_is_basin_robust(label: str, system: cr3bp.CR3BPSystem) -> None:
    """Every one of 11 evenly-spaced seeds across a +-2e-4 window around the
    sourced Table 4.1 seed converges to the SAME eigenvalue -- not an
    isolated numerical fluke (basin_robustness_scan, reused directly from
    the Jovian module)."""
    x0_guess, period_guess = stf._table41_seed_nondim(label, system)
    hc = stf._HALF_CROSSINGS[label]
    results = jrf.basin_robustness_scan(
        system,
        x0_lo=x0_guess - 2e-4,
        x0_hi=x0_guess + 2e-4,
        n_seeds=11,
        jacobi=stf.VAQUERO_C,
        ydot0_sign=1.0,
        half_crossings=hc,
        period_guess=period_guess,
    )
    converged = [c for _, c in results if c is not None]
    assert len(converged) == len(results), f"{label}: not all seeds converged"
    eigs = {round(c.max_eigenvalue, 2) for c in converged}
    assert len(eigs) == 1, f"{label}: basin scan found multiple distinct eigenvalues: {eigs}"


# ---------------------------------------------------------------------------
# (4) The honest Table 4.1 gate itself -- reported family-by-family, no
# fudged tolerance. 2/4 rows fully pass (3:4, L1); 6:5 fails on eigenvalue
# only (real, small, well-characterized near-miss); L2 fails on period only
# (documented, likely source-errata value -- see TABLE41_L2_PERIOD_ERRATA_NOTE).
# ---------------------------------------------------------------------------


def test_table41_gate_honest_report(system: cr3bp.CR3BPSystem) -> None:
    rows = {r.label: r for r in stf.gate_report(system)}
    assert set(rows) == set(stf.TABLE41_TARGETS)

    # 3:4: near-machine-precision reproduction on every axis.
    assert rows["3:4"].passed, rows["3:4"]
    assert rows["3:4"].eigenvalue_rel_err < 1e-5
    assert rows["3:4"].period_rel_err < 1e-3
    assert rows["3:4"].ic_confirmed

    # L1: same, near-machine-precision.
    assert rows["L1"].passed, rows["L1"]
    assert rows["L1"].eigenvalue_rel_err < 1e-5
    assert rows["L1"].period_rel_err < 1e-3
    assert rows["L1"].ic_confirmed

    # 6:5: HONEST FAIL on eigenvalue (2.34e-3, just outside the 1e-3 gate)
    # despite an excellent IC/period match -- not fudged, not silently
    # loosened. A future fix that finds a better match should update this
    # assertion, not work around it.
    row_65 = rows["6:5"]
    assert not row_65.eigenvalue_confirmed, row_65
    assert 1e-3 < row_65.eigenvalue_rel_err < 1e-2, row_65
    assert row_65.period_confirmed, row_65
    assert row_65.ic_confirmed, row_65
    assert not row_65.passed

    # L2: eigenvalue confirms to near-machine precision, IC confirms too,
    # but the row's OWN printed period is off by a factor of ~9.27 -- the
    # documented, likely-errata anomaly (TABLE41_L2_PERIOD_ERRATA_NOTE).
    row_l2 = rows["L2"]
    assert row_l2.eigenvalue_confirmed, row_l2
    assert row_l2.eigenvalue_rel_err < 1e-5, row_l2
    assert row_l2.ic_confirmed, row_l2
    assert not row_l2.period_confirmed, row_l2
    assert row_l2.period_rel_err > 0.5, row_l2  # a huge, not-borderline miss
    assert not row_l2.passed

    n_passed = sum(r.passed for r in rows.values())
    assert n_passed == 2


def test_l2_period_errata_note_is_documented() -> None:
    """The L2 period anomaly must be documented in an importable, testable
    place -- not just prose in the module docstring."""
    assert "79.7260" in stf.TABLE41_L2_PERIOD_ERRATA_NOTE
    assert "8.603" in stf.TABLE41_L2_PERIOD_ERRATA_NOTE


def test_eigenvalue_sensitivity_to_mu_is_measured_not_assumed(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Measures (not assumes) how sensitive each row's recovered eigenvalue
    is to a small mu perturbation. A +0.1% mu shift moves every row's
    eigenvalue by a comparable RELATIVE amount (not absolute -- see below).

    `#769` (``docs/notes/2026-08-05-769-saturn-titan-65-eigenvalue.md``) used
    exactly this comparable-relative-sensitivity result, plus 3:4's shift
    running in the OPPOSITE direction from 6:5/L1/L2's, to EXCLUDE mu
    imprecision as the cause of 6:5's own near-miss (module docstring
    finding 2) -- a fixed mu precision floor would produce a comparable
    RELATIVE eigenvalue error in ALL four rows, not the observed 1e-6-level
    match in three rows and 2.3e-3 miss in the fourth. See the module
    docstring finding 2: a follow-up C-sensitivity check was ALSO
    adversarially reviewed and found not to reconcile all four rows either
    (L1/L2's own baselines pin C far tighter than 6:5 would need) -- the
    miss is best characterized as row-specific/source-side, not a single
    shared-parameter correction of any kind."""
    mu_perturbed = stf.VAQUERO_MU * 1.001
    perturbed_system = cr3bp.CR3BPSystem(
        mu=mu_perturbed,
        primary=system.primary,
        secondary=system.secondary,
        l_km=system.l_km,
        t_s=system.t_s,
    )
    rel_shifts: dict[str, float] = {}
    for label, target in stf.TABLE41_TARGETS.items():
        baseline = stf.recover_table41_candidate(label, system)
        perturbed = stf.recover_table41_candidate(label, perturbed_system)
        shift = abs(perturbed.max_eigenvalue - baseline.max_eigenvalue) / baseline.max_eigenvalue
        # Every row must show SOME measurable sensitivity (not numerically
        # inert to a real 0.1% mu change) but stay bounded/sane (not a
        # chaotic order-of-magnitude jump) -- both directions matter: a
        # `shift` of exactly 0 would mean this test isn't measuring anything
        # real, and a huge shift would mean the seed sits in a genuinely
        # unstable/chaotic basin, invalidating the whole comparison.
        assert 1e-5 < shift < 0.5, f"{label}: mu-sensitivity shift={shift:.2e} (target {target})"
        rel_shifts[label] = shift
    # Regression guard for #769's own finding: relative mu-sensitivity is
    # comparable ACROSS rows (this is the fact that excludes mu-imprecision
    # as 6:5's root cause -- see docstring above). Generous bounds (not
    # tuned tight) so this doesn't become brittle to minor corrector
    # changes; it only needs to keep catching a return to the "6:5 is
    # uniquely mu-sensitive" hypothesis this module's docstring used to make.
    max_shift, min_shift = max(rel_shifts.values()), min(rel_shifts.values())
    assert max_shift / min_shift < 5.0, (
        f"relative mu-sensitivity spread across rows grew unexpectedly large: {rel_shifts}"
    )


# ---------------------------------------------------------------------------
# (5) Connection stage explicitly out of scope for #765 (its own Task-B
# analog, a later task) -- this module deliberately exposes no manifold/
# homoclinic/heteroclinic machinery.
# ---------------------------------------------------------------------------


def test_module_exposes_no_connection_machinery() -> None:
    connection_names = {"find_homoclinic", "find_heteroclinic", "correct_connection"}
    assert connection_names.isdisjoint(set(stf.__all__))
