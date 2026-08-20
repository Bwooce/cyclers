"""Tests for `#839`: targeted Wu(3:1) <-> Ws(2:1) Earth-Moon heteroclinic
search AT C=3.13, touching `vaquero-31-c313-em-resonant-po-2013`
(docs/notes/2026-08-2*-839-c313-targeted-search.md).

Reproduce-before-trust discipline, mirroring `#827`'s own test module
(``tests/search/test_kumar_em_resonant_heteroclinics.py``): the connection-
evidence tests below rebuild a result from scratch through the library
primitives -- fresh node builds (:func:`build_node31_c313` re-derives the
CATALOGUE row's own state and rejects disagreement;
:func:`build_node21_c313` step-continues a Kumar Table-6 bracket and rejects
a branch jump), a fresh Newton connection refinement seeded ONLY by the
recorded phase indices (``tau_u``, ``tau_s``, ``k_u``, ``k_s``, ``branch_u``,
``branch_s`` -- ``refine_connection`` ignores a seed's ``x``/``xdot``/
``ydot``), and a fresh full `#822` verification battery. Nothing here replays
the committed JSON directly -- the achieved values are pinned against
``data/found/839_c313_targeted_search/results.json``'s own recorded values as
regression targets, following the `#810`/`#827` precedent. None are
``@pytest.mark.slow`` -- an evidence test CI never runs is an unverified
claim (``feedback_delegation_fresh_agent_not_fork``).
"""

from __future__ import annotations

import math

import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.kumar_em_resonant_heteroclinics as keh
import cyclerfinder.search.vaquero_c313_targeted_search as m
import cyclerfinder.search.vaquero_em_cycler_connections as vcc
import cyclerfinder.search.vaquero_em_cyclers as vem


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return m.em_system()


@pytest.fixture(scope="module")
def node31(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    return m.build_node31_c313(system)


@pytest.fixture(scope="module")
def node21_and_prov(system: cr3bp.CR3BPSystem) -> tuple[jrc.ResonantNode, m.Node21Provenance]:
    return m.build_node21_c313(system)


# ---------------------------------------------------------------------------
# (1) Sourced constants and mu identity.
# ---------------------------------------------------------------------------


def test_em_system_uses_this_projects_registry_mu_not_kumars() -> None:
    """The catalogue row's own ``mass_ratio`` is this project's registry mu,
    NOT Kumar's printed mu -- the two differ at the ~1.24e-10 absolute
    level (module docstring, "MU")."""
    system = m.em_system()
    assert system.mu == pytest.approx(0.01215058439469525, abs=0.0)
    assert abs(system.mu - keh.KUMAR_MU) > 1e-11


def test_c313_matches_the_catalogue_rows_own_recorded_jacobi_constant() -> None:
    assert m.C313 == 3.129999999999993


def test_epsilon_is_kumar_epsilon_reused_verbatim() -> None:
    """Same manifold-offset magnitude `#827` derived for this SAME family
    branch's extreme 2:1 saddle (module docstring, "EPSILON")."""
    assert m.EPSILON == keh.KUMAR_EPSILON
    assert m.EPSILON == 1e-4


# ---------------------------------------------------------------------------
# (2) Node re-derivation (never a raw trust of the catalogue/Kumar prints).
# ---------------------------------------------------------------------------


def test_node31_rederives_the_catalogue_rows_own_state(node31: jrc.ResonantNode) -> None:
    """node31 IS the catalogue row `vaquero-31-c313-em-resonant-po-2013`
    itself, re-converged independently and checked against its own recorded
    ``state_nd``/``jacobi_constant``."""
    assert node31.converged
    assert node31.jacobi == pytest.approx(m.C313, abs=1e-9)
    assert abs(node31.state0[0] - m.CATALOGUE_STATE_ND_31[0]) < m.NODE31_IC_ABS_TOL
    assert abs(node31.state0[4] - m.CATALOGUE_STATE_ND_31[4]) < m.NODE31_IC_ABS_TOL


def test_node21_continuation_does_not_jump_branches(
    node21_and_prov: tuple[jrc.ResonantNode, m.Node21Provenance],
) -> None:
    """The C=3.13 2:1 member, reached by 0.001-step continuation from the
    Kumar Table-6 C=3.10 bracket, must land within the interpolation-check
    tolerance of a linear interpolation of the two own-mu bracket endpoints
    -- a 0.01 step was measured (module docstring) to jump onto a different
    branch entirely; this asserts the actually-used 0.001 step does not."""
    node21, prov = node21_and_prov
    assert node21.converged
    assert prov.interp_distance < m.NODE21_INTERP_CHECK_TOL
    assert prov.interp_distance == pytest.approx(0.0037025615644105656, abs=1e-9)
    assert len(prov.trace) == 31  # 30 x 0.001 steps + the C=3.10 anchor


def test_node21_at_c310_matches_827s_own_kumar_module_node(
    system: cr3bp.CR3BPSystem, node21_and_prov: tuple[jrc.ResonantNode, m.Node21Provenance]
) -> None:
    """node21's own construction anchor (the C=3.10 bracket, built at THIS
    project's registry mu) must be the SAME physical orbit `#827`'s own
    ``build_kumar_node(kumar_system(), 2, 3.10)`` reproduces at Kumar's mu --
    period and |lambda| agreeing to ~1e-7 relative despite the ~1.24e-10
    absolute mu difference confirms node21 is not a wrong branch (advisor
    review finding, see the results note)."""
    _node21, prov = node21_and_prov
    ksys = keh.kumar_system()
    kumar_node21_310 = keh.build_kumar_node(ksys, 2, 3.10)
    assert prov.bracket_310_lambda == pytest.approx(333.04100533262795, rel=1e-6)
    from cyclerfinder.genome.heteroclinic_cycle import _planar_floquet_pair

    lam_kumar, _, _, _ = _planar_floquet_pair(
        ksys, kumar_node21_310.state0, kumar_node21_310.period
    )
    assert prov.bracket_310_lambda == pytest.approx(lam_kumar, rel=1e-6)


def test_node21_rejects_a_branch_jumping_continuation_step(system: cr3bp.CR3BPSystem) -> None:
    """The interpolation-check gate (``build_node21_c313``'s
    ``interp_check_tol``) is a real gate: an absurdly coarse continuation
    step (0.5, jumping straight from C=3.10 to C=3.13 in one shot -- the same
    failure mode measured for the 0.01 step, module docstring) lands the
    fixed-``half_crossings``-index corrector on a different branch, and the
    resulting disagreement with the own-mu bracket interpolation MUST raise
    a ``ValueError``, not silently degrade."""
    with pytest.raises(ValueError, match="disagrees with the linear own-mu bracket interpolation"):
        m.build_node21_c313(
            system,
            step=0.5,  # absurdly coarse: forces a branch jump before C313 is even reached
        )


def test_kumar_table6_31_interpolation_lands_near_the_catalogue_state() -> None:
    """Pure arithmetic: linearly interpolating Kumar's OWN printed Table-6
    3:1 rows (C=3.10, 3.15) to C=3.13 lands close to the catalogue row's own
    independently-DERIVED ``state_nd`` -- the family-identity evidence this
    module's docstring cites ("~2e-3" scale), made into an assertion."""
    x310, ydot310 = keh.KUMAR_TABLE6_ICS[(3, 3.10)]
    x315, ydot315 = keh.KUMAR_TABLE6_ICS[(3, 3.15)]
    frac = (3.13 - 3.10) / (3.15 - 3.10)
    interp_x0 = x310 + frac * (x315 - x310)
    interp_ydot0 = ydot310 + frac * (ydot315 - ydot310)
    dist = math.hypot(
        interp_x0 - m.CATALOGUE_STATE_ND_31[0], interp_ydot0 - m.CATALOGUE_STATE_ND_31[4]
    )
    assert dist < 0.005
    assert dist == pytest.approx(0.0017902178701488518, abs=1e-9)


# ---------------------------------------------------------------------------
# (3) Independent from-scratch connection reconstruction, pinned against the
#     committed run record (data/found/839_c313_targeted_search/results.json)
#     as a regression target, not as a value to replay.
# ---------------------------------------------------------------------------


def test_connection_reconverges_at_n_tau_48(
    system: cr3bp.CR3BPSystem,
    node31: jrc.ResonantNode,
    node21_and_prov: tuple[jrc.ResonantNode, m.Node21Provenance],
) -> None:
    """The n_tau=48 phase-grid's converged, verified Wu(3:1) -> Ws(2:1) hit
    at C=3.13, seeded ONLY by the recorded phase indices."""
    node21, _prov = node21_and_prov
    seed = vcc.ConnectionSeed(
        distance=0.0,
        branch_u=-1,
        branch_s=1,
        k_u=17,
        k_s=14,
        tau_u=5.8726639534049845,
        tau_s=6.473865290398373,
        x=0.0,
        xdot=0.0,
        ydot=0.0,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=m.EPSILON)
    assert conn.converged
    assert conn.residual < 1e-9
    assert abs(conn.crossing_xv[0] - (-0.21975612699700106)) < 1e-6
    assert abs(conn.crossing_xv[1] - 1.0625086172960396) < 1e-6

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=m.EPSILON)
    assert ev.passed
    assert ev.ydot_signs_match
    assert ev.full_state_gap < 1e-6
    assert ev.full_state_gap == pytest.approx(1.652810462522925e-08, abs=1e-11)
    assert ev.ghost_distance_from > vcc.GHOST_GUARD_DELTA
    assert ev.ghost_distance_to > vcc.GHOST_GUARD_DELTA


def test_connection_reconverges_at_n_tau_64_a_different_crossing(
    system: cr3bp.CR3BPSystem,
    node31: jrc.ResonantNode,
    node21_and_prov: tuple[jrc.ResonantNode, m.Node21Provenance],
) -> None:
    """The n_tau=64 phase-grid's converged, verified hit is a DIFFERENT
    specific manifold crossing (different k_u/k_s) from the n_tau=48 hit --
    both pass the same unmodified battery, corroborating EXISTENCE of a
    genuine connection at this C via two independent phase grids, without
    claiming the two hits are the same object (results note, "ROBUSTNESS")."""
    node21, _prov = node21_and_prov
    seed = vcc.ConnectionSeed(
        distance=0.0,
        branch_u=-1,
        branch_s=1,
        k_u=20,
        k_s=13,
        tau_u=4.586306798124203,
        tau_s=6.421078095155724,
        x=0.0,
        xdot=0.0,
        ydot=0.0,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=m.EPSILON)
    assert conn.converged
    assert conn.residual < 1e-9
    assert abs(conn.crossing_xv[0] - 0.5638883664836785) < 1e-6
    assert abs(conn.crossing_xv[1] - (-0.5802586614639252)) < 1e-6

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=m.EPSILON)
    assert ev.passed
    assert ev.ydot_signs_match
    assert ev.full_state_gap < 1e-6
    assert ev.full_state_gap == pytest.approx(7.673612548784641e-08, abs=1e-11)

    # Confirm it is genuinely a distinct crossing from the n_tau=48 hit.
    assert (conn.k_u, conn.k_s) != (17, 14)


def test_half_crossings_uses_vaquero_convention_for_node31(system: cr3bp.CR3BPSystem) -> None:
    """node31 is built with ``half_crossings=3``
    (``VAQUERO_HALF_CROSSINGS``), the SAME topology convention `#799`/`#811`
    used for this exact catalogue row -- not an independent guess."""
    assert vem.VAQUERO_HALF_CROSSINGS == 3
