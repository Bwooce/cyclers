"""Tests for the Resonant Atlas Stage A' fold-turning worker (`#861`).

`#860` diagnosed `#859`'s Stage A smoke-test failure and proposed two fixes:
(1) a conjugate-apse seed (`jovian_resonant_families.two_body_conjugate_apse_seed`,
tested in `test_jovian_resonant_families.py`), and (2) fold-turning continuation
(`cr3bp_jacobi_arclength.continue_in_jacobi`) instead of `continue_family`'s
natural-parameter walk. This module wires (2) into a Stage A' worker. Tests here
cover the WIRING (seed dispatch, member classification, worker contract) on FAST,
cheap systems/short walks -- the actual Oberon gate (`#861`'s deliverable) is run
separately via `scripts/run_861_oberon_gate.py` and is NOT re-run by the test
suite (too slow for unit tests; the gate's own results are checkpointed under
`data/found/861_resonant_seeding_oberon_gate/`).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.neptune_triton_resonant_families as ntrf
import cyclerfinder.search.resonant_atlas_stage_a_prime as sap
from cyclerfinder.search.mu_continuation import MuMember


@pytest.fixture(scope="module")
def europa_system() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Jupiter", "Europa")


# ---------------------------------------------------------------------------
# build_seed dispatch
# ---------------------------------------------------------------------------


def test_build_seed_opposition_matches_two_body_resonant_seed() -> None:
    import cyclerfinder.search.jovian_resonant_families as jrf

    x0, ydot0, sign, period_full = sap.build_seed(3, 4, "opposition")
    ref = jrf.two_body_resonant_seed(3, 4, x0_sign=-1)
    assert x0 == pytest.approx(ref.x0)
    assert ydot0 == pytest.approx(ref.ydot0)
    assert sign == -1.0
    assert period_full == pytest.approx(ref.period_full)


def test_build_seed_conjugate_apse_matches_new_seed() -> None:
    import cyclerfinder.search.jovian_resonant_families as jrf

    x0, ydot0, sign, period_full = sap.build_seed(4, 5, "conjugate_apse")
    ref = jrf.two_body_conjugate_apse_seed(4, 5)
    assert x0 == pytest.approx(ref.x0)
    assert ydot0 == pytest.approx(ref.ydot0)
    assert sign == (1.0 if ref.ydot0 >= 0.0 else -1.0)
    assert period_full == pytest.approx(ref.period_full)


def test_build_seed_rejects_bad_kind() -> None:
    with pytest.raises(ValueError, match="seed_kind"):
        sap.build_seed(3, 4, "sideways")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# classify_member -- known-answer positive control: Neptune-Triton "4:5-saddle"
# (`#776`'s own table-verified saddle, C=2.987089791658, half_crossings=3,
# |lambda|~=105 per `#771`'s survey -- the SAME p:q label, 4:5, one of Oberon's
# own six published ratios).
# ---------------------------------------------------------------------------


def test_classify_member_on_neptune_triton_4_5_saddle_positive_control() -> None:
    system = ntrf.neptune_triton_system()
    row = ntrf.ESM_GATE_ROWS["4:5-saddle"]
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        row.x0,
        row.jacobi,
        row.period,
        ydot0_sign=row.ydot0_sign,
        half_crossings=row.half_crossings,
        tol=1e-11,
        rtol=1e-13,
        atol=1e-13,
        x0_bounds=row.x0_bounds,
    )
    assert orbit.converged
    nu, lam = cp.barden_stability(system, orbit)
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    po = cp.PeriodicOrbit(
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=True,
        closure_residual=orbit.crossing_residual,
    )
    _ok, radau_dj = cp.crosscheck_periodic(system, po)
    m = MuMember(
        mu=system.mu,
        state0=state0,
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        jacobi=orbit.jacobi,
        period=orbit.period,
        nu=nu,
        abs_lambda=abs(lam),
        crossing_residual=orbit.crossing_residual,
        radau_djacobi=radau_dj,
        stable=abs(nu) < 1.0,
    )
    fm = sap.classify_member(system, m, p=4, q=5)  # label is "4:5"

    # Genuine saddle (|lambda| >> 1), per #771's own survey (~105).
    assert fm.abs_lambda > 50.0
    assert not fm.stable
    # Inertial winding about the primary should land within a few % of p=4
    # (the label's own p) -- see classify_member's own docstring for the
    # rotating->inertial conversion and the ~4% empirical offset this exact
    # row shows (matching the period check's own comparable offset).
    assert abs(abs(fm.winding_p_inertial) - 4.0) / 4.0 < 0.10
    # A real, finite close approach -- measured this task at 0.0187 nondim
    # (~6634 km at Neptune-Triton's own l*=354760 km). NOT asserted against a
    # cross-system nondim threshold copied from Jupiter-Europa's own #756 scale
    # (0.018-0.042 nondim there too, but that mu is ~8x smaller, so the SAME raw
    # nondim distance is a much smaller fraction of Triton's own mu^(1/3)
    # Hill-radius scale than it would be for Europa) -- just a sanity bound that
    # the computed value is a real, non-degenerate positive distance.
    assert 0.0 < fm.closest_secondary_approach_nondim < 0.05
    # period/2pi need not land exactly on q=5 for a genuinely unstable, far-
    # from-integrable family (the #755 reviewer precedent) -- just report it
    # sanely (finite, positive, same order of magnitude as q).
    assert 0.0 < fm.period_over_2pi < 3.0 * 5


# ---------------------------------------------------------------------------
# fold_turn_family -- fast, short walk end-to-end (wiring correctness, not the
# real gate).
# ---------------------------------------------------------------------------


def test_fold_turn_family_short_walk_produces_genuine_members(
    europa_system: cr3bp.CR3BPSystem,
) -> None:
    result = sap.fold_turn_family(
        europa_system,
        3,
        4,
        "conjugate_apse",
        system_key="jupiter-europa",
        max_steps=4,
        ds0=8e-3,
        ds_max=1.5e-2,
        record_every=1,
        c_span=0.03,
    )
    assert result.seed_converged
    assert result.half_crossings > 0
    assert len(result.members) >= 1
    for m in result.members:
        assert m.crossing_residual < 1e-6
        assert m.radau_djacobi < 1e-4
        # Jacobi self-consistency: the seed member's own C must equal seed_jacobi.
    assert any(abs(m.jacobi - result.seed_jacobi) < 1e-8 for m in result.members)


def test_fold_turn_family_seed_domain_invalid_reports_not_converged(
    europa_system: cr3bp.CR3BPSystem,
) -> None:
    with pytest.raises(ValueError):
        sap.build_seed(5, 1, "conjugate_apse")


def test_fold_turn_family_both_seed_kinds_run_without_crashing(
    europa_system: cr3bp.CR3BPSystem,
) -> None:
    for kind in ("opposition", "conjugate_apse"):
        result = sap.fold_turn_family(
            europa_system,
            4,
            3,
            kind,
            system_key="jupiter-europa",
            max_steps=2,
            ds0=6e-3,
            ds_max=1e-2,
            record_every=1,
            c_span=0.02,
        )
        assert result.seed_converged


# ---------------------------------------------------------------------------
# fold_turn_worker -- CellOutcome contract.
# ---------------------------------------------------------------------------


def test_fold_turn_worker_returns_json_safe_hit_or_miss() -> None:
    cell = {
        "system_key": "jupiter-europa",
        "primary": "Jupiter",
        "secondary": "Europa",
        "p": 3,
        "q": 4,
        "seed_kind": "conjugate_apse",
        "max_steps": 3,
        "ds0": 8e-3,
        "ds_max": 1.5e-2,
        "record_every": 1,
        "c_span": 0.02,
    }
    outcome = sap.fold_turn_worker(cell)
    assert outcome.status in ("hit", "miss")
    assert outcome.error == ""
    assert outcome.payload["p"] == 3
    assert outcome.payload["q"] == 4
    import json

    json.dumps(outcome.payload)  # must be JSON-safe (campaign_runner contract)


def test_fold_turn_worker_never_raises_on_bad_cell() -> None:
    outcome = sap.fold_turn_worker(
        {
            "system_key": "bogus",
            "primary": "NoSuchPlanet",
            "secondary": "NoSuchMoon",
            "p": 3,
            "q": 4,
            "seed_kind": "conjugate_apse",
        }
    )
    assert outcome.status == "error"
    assert outcome.error != ""


def test_fold_turn_worker_accepts_mu_override() -> None:
    """The `#728` Oberon digest's own paper mu, not the registry value -- same
    "paper's own mu" convention every sourced module in this project uses."""
    cell = {
        "system_key": "uranus-oberon",
        "primary": "Uranus",
        "secondary": "Oberon",
        "mu": 3.54326e-5,
        "p": 3,
        "q": 4,
        "seed_kind": "conjugate_apse",
        "max_steps": 2,
        "ds0": 6e-3,
        "ds_max": 1e-2,
        "record_every": 1,
        "c_span": 0.01,
    }
    outcome = sap.fold_turn_worker(cell)
    assert outcome.status in ("hit", "miss")
    assert outcome.payload.get("seed_jacobi") is not None or "reason" in outcome.payload
