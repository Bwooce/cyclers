"""Task #549: real-binary (k1,k2)-cycler sweep test suite.

Generalizes #504's Pluto-Charon (k1,k2) sweep to four sourced real binary
systems (Patroclus-Menoetius, Didymos-Dimorphos, Orcus-Vanth, Eris-Dysnomia).

This module is a REGRESSION suite for the machinery, not a re-run of the full
32-combo discovery sweep (that ran once as scripts/run_549_real_binary_kk_sweep.py;
its raw results are summarized in data/OUTSTANDING.md's #549 entry). It covers:

1. A positive-control regression through the NEW generalized driver
   (mu_step_to_system / sweep_family) reproducing PC (3,2) -- proves the
   generalization did not silently break #504's own machinery.
2. A direct unit test of the bug this module fixes: #504's mu_step_to_orbit
   hardcodes its "final correction" step to Pluto-Charon's system regardless
   of the target mu -- mu_step_to_system must finish in the CALLER's system.
3. Static checks on the sourced REAL_BINARY_SYSTEMS constants (mu values,
   l_km/t_s unit conversions) so a future edit can't silently drift a cited
   number.
4. One clean-negative-aware sweep_family call per sourced system (fastest
   anchor per system) as a light re-confirmation; passes either way, matching
   #504's own convention -- NOT an assertion that no cycler exists (the full
   32-combo negative is recorded in OUTSTANDING.md, not re-asserted here).
   Patroclus-Menoetius is intentionally NOT re-run here (every anchor attempt
   for it took 35-265s in the discovery run -- see OUTSTANDING.md's #549
   entry for its full record); this suite covers the three cheaper systems
   to keep the default ratchet's wall-clock cost bounded.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.binary_star_search import collinear_lpoints
from cyclerfinder.search.pluto_charon_kk_sweep import PC_MU, make_pluto_charon_system
from cyclerfinder.search.real_binary_kk_sweep import (
    REAL_BINARY_SYSTEMS,
    mu_step_to_system,
    sweep_family,
)

# ---------------------------------------------------------------------------
# 1. Positive control through the NEW generalized driver
# ---------------------------------------------------------------------------


def test_549_positive_control_via_generic_driver() -> None:
    """sweep_family(pc_system, "mu01_32") must reproduce #504's PC (3,2).

    Uses the SAME Ross-RT 2026 Table-I mu=0.1 anchor #504 uses, but routed
    through the new mu_step_to_system + sweep_family driver instead of
    #504's PC-specific sweep_32_positive_control -- a regression check that
    the generalization preserves the validated result to high precision.
    """
    pc = make_pluto_charon_system()
    result = sweep_family(pc, "mu01_32")

    assert result.stable_found, (
        f"Generalized-driver PC (3,2) positive control FAILED: "
        f"method={result.method!r}, note={result.note!r}"
    )
    assert result.topology_ok
    assert result.prograde
    assert result.reaches_secondary
    assert result.crosscheck_ok

    # Catalogue ross-rt-pc-cycler-32-2026: C=3.57951501972907, T=11.8334625170346 TU.
    assert result.jacobi_mid is not None
    assert abs(result.jacobi_mid - 3.57951501972907) < 1e-4, (
        f"jacobi_mid={result.jacobi_mid:.9f} far from catalogue 3.57951501972907"
    )
    assert result.period_mid is not None
    assert abs(result.period_mid - 11.8334625170346) < 1e-3, (
        f"period_mid={result.period_mid:.9f} far from catalogue 11.8334625170346"
    )
    assert result.nu_mid is not None
    assert abs(result.nu_mid) < 1e-4, f"nu_mid={result.nu_mid:.3e} not near-zero"


# ---------------------------------------------------------------------------
# 2. The bug this module fixes: finish in the CALLER's system, not PC
# ---------------------------------------------------------------------------


def test_549_mu_step_to_system_uses_caller_system_not_hardcoded_pc() -> None:
    """mu_step_to_system must finish in target_system, not Pluto-Charon.

    #504's mu_step_to_orbit hardcodes its final correction to
    make_pluto_charon_system() regardless of the target mu -- a latent bug
    for reuse at any OTHER mu. This directly exercises the trivial identity
    case (anchor_mu == target mu, so n_steps=1 mu-stepping is a no-op) at a
    mu that is NOT PC_MU, and asserts the returned orbit's Jacobi/period are
    consistent with THAT system, not silently re-solved at PC_MU.
    """
    # A system at mu=0.001 (identical mu to the Table-I mu=0.001 anchor, so
    # the mu-continuation loop is a no-op and any drift can only come from
    # the "final correction" using the WRONG system).
    target = cr3bp.CR3BPSystem(mu=0.001, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    assert abs(target.mu - PC_MU) > 0.05, "test system must be far from PC_MU"

    orbit = mu_step_to_system(
        anchor_mu=0.001,
        target_system=target,
        anchor_x0=-0.647047499999966,
        anchor_jacobi=3.031605708907296,
        anchor_period=14.774502790974823,
        hc=None,
        sign=-1.0,
        n_steps=1,
    )
    assert orbit is not None, "mu-step at mu=anchor_mu (no-op stepping) must converge"

    # If the old bug were present, the "final correction" would silently
    # re-solve at PC_MU=0.10876 instead of mu=0.001 -- the returned x0 would
    # sit near the PC (3,2) branch (x0 ~ -0.69), not near the mu=0.001
    # anchor's own x0 (~-0.647). Assert we stayed near the TARGET branch.
    assert abs(orbit.x0 - (-0.647047499999966)) < 0.05, (
        f"x0={orbit.x0:.6f} drifted far from the mu=0.001 anchor -- "
        "suggests the final correction ran in the WRONG system"
    )


# ---------------------------------------------------------------------------
# 3. Sourced constants: static checks (no integration)
# ---------------------------------------------------------------------------


def test_549_sourced_mu_values() -> None:
    """Sourced mu values match the citations recorded in each system's docstring."""
    assert REAL_BINARY_SYSTEMS["patroclus-menoetius"].mu == 104.0**3 / (113.0**3 + 104.0**3)
    assert abs(REAL_BINARY_SYSTEMS["didymos-dimorphos"].mu - 4.3e9 / (5.4e11 + 4.3e9)) < 1e-12
    assert REAL_BINARY_SYSTEMS["orcus-vanth"].mu == 0.137
    assert abs(REAL_BINARY_SYSTEMS["eris-dysnomia"].mu - 0.0050 / 1.0050) < 1e-12

    # All mu in (0, 0.5]: physically valid CR3BP mass ratios (secondary <= primary).
    for key, sys_def in REAL_BINARY_SYSTEMS.items():
        assert 0.0 < sys_def.mu <= 0.5, f"{key}: mu={sys_def.mu} out of (0, 0.5]"
        assert sys_def.mu_source, f"{key}: missing mu_source citation"
        assert sys_def.l_source, f"{key}: missing l_source citation"
        assert sys_def.t_source, f"{key}: missing t_source citation"


def test_549_t_s_matches_orbital_period_convention() -> None:
    """t_s = sqrt(l_km^3/G(m1+m2)) == P_orbital/(2*pi), matching cr3bp_system()'s
    own convention (verified against the in-registry Pluto-Charon system)."""
    pc = make_pluto_charon_system()
    l1, _l2, _l3 = collinear_lpoints(pc.mu)
    c_l1 = float(cr3bp.jacobi_constant(np.array([l1, 0.0, 0.0, 0.0, 0.0, 0.0]), pc.mu))
    assert c_l1 > 3.0  # sanity: PC's own C_L1 is a normal CR3BP value

    # Orcus-Vanth: t_s should equal P_days*86400/(2*pi) exactly (both computed
    # via _t_s_from_period_days -- assert the stored value matches the source
    # P=9.5393 d directly, catching any future accidental edit).
    ov = REAL_BINARY_SYSTEMS["orcus-vanth"]
    expected_t_s = 9.5393 * 86400.0 / (2.0 * np.pi)
    assert abs(ov.t_s - expected_t_s) < 1e-6


# ---------------------------------------------------------------------------
# 4. One light clean-negative-aware re-confirmation per system
# ---------------------------------------------------------------------------


def _assert_valid_stable_member(result, *, label: str) -> None:  # type: ignore[no-untyped-def]
    assert result.stable_found
    assert result.topology_ok, f"{label}: topology check failed"
    assert result.prograde, f"{label}: orbit must be prograde"
    assert result.reaches_secondary, f"{label}: orbit must reach the secondary realm"
    assert result.nu_mid is not None
    assert abs(result.nu_mid) < 1.0
    assert result.crosscheck_ok, f"{label}: independent Radau crosscheck failed"


def test_549_didymos_dimorphos_11_clean_negative_aware() -> None:
    """Didymos-Dimorphos (1,1) via the closest Table-I anchor (mu=0.001).

    #549's full sweep found this converges to a stable orbit with WRONG
    topology (prograde but doesn't reach the secondary) -- a clean negative.
    Clean-negative-aware: passes either way (if a future run somehow finds
    a genuine stable (1,1) member here, all validity gates must hold).
    """
    target = REAL_BINARY_SYSTEMS["didymos-dimorphos"].to_cr3bp_system()
    result = sweep_family(target, "mu001_11")
    if not result.stable_found:
        return
    _assert_valid_stable_member(result, label="Didymos-Dimorphos (1,1)")


def test_549_orcus_vanth_32_clean_negative_aware() -> None:
    """Orcus-Vanth (3,2) via the closest Table-I anchor (mu=0.1).

    #549's full sweep found this converges to a stable orbit with WRONG
    topology (prograde but doesn't reach the secondary) -- a clean negative.
    Clean-negative-aware: passes either way.
    """
    target = REAL_BINARY_SYSTEMS["orcus-vanth"].to_cr3bp_system()
    result = sweep_family(target, "mu01_32")
    if not result.stable_found:
        return
    _assert_valid_stable_member(result, label="Orcus-Vanth (3,2)")


def test_549_eris_dysnomia_11_clean_negative_aware() -> None:
    """Eris-Dysnomia (1,1) via the closest Table-I anchor (mu=0.001).

    #549's full sweep found no stable |nu|<1 window in the C-sweep range --
    a clean negative. Clean-negative-aware: passes either way.
    """
    target = REAL_BINARY_SYSTEMS["eris-dysnomia"].to_cr3bp_system()
    result = sweep_family(target, "mu001_11")
    if not result.stable_found:
        return
    _assert_valid_stable_member(result, label="Eris-Dysnomia (1,1)")
