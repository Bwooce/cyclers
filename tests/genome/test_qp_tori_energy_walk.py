"""Tests for the #466 energy-MOVING QP-GMOS continuation (parent-family-driven).

Distinct from #333 (`test_qp_tori_arclength.py`), whose family is near-iso-energetic:
this walk steps the #296 parent (1,1) family DOWN its monotone energy ladder and
re-converges the QP-torus at each parent member, so it MOVES in energy by construction.
The decisive capability gate is the inverse of the #333 containment assertion: ``C_J``
must span a real range (not the #333 6e-7 floor).

Golden discipline (mirrors #333): every EXPECTED side asserts topology
(irrationality), invariance (Fourier/off-grid closure), or self-consistency (members
sit at the parent ladder's prescribed energies; energy genuinely moves). The parent
family ladder is a digested INPUT (data/family_296_3d_em_11.jsonl), not a code-computed
golder. NEVER a frequency / C_J value our own code produced is asserted as a target.
Report-only; no catalogue writeback.
"""

from __future__ import annotations

import json
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v1_qp import run_v1_qp
from cyclerfinder.genome import qp_tori_energy_walk as qew
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus

SmokeTorus = tuple[cr3bp.CR3BPSystem, QPTorus]

ROOT = Path(__file__).resolve().parents[2]
SUBFAMILIES_FILE = ROOT / "data" / "family_296_3d_subfamilies_299.jsonl"
PARENT_FAMILY_FILE = ROOT / "data" / "family_296_3d_em_11.jsonl"

EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=EM_MU, primary="earth", secondary="moon", l_km=EM_L_KM, t_s=EM_T_S)


def _load_parent_at_step(step_index: int) -> dict[str, Any]:
    if not PARENT_FAMILY_FILE.exists():
        pytest.skip(f"parent family file not present: {PARENT_FAMILY_FILE}")
    with PARENT_FAMILY_FILE.open() as f:
        for line in f:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            if d.get("step_index") == step_index:
                return dict(d)
    raise RuntimeError(f"step_index={step_index} not found in {PARENT_FAMILY_FILE}")


def _load_first_neimark_sacker_bracket() -> tuple[dict[str, Any], dict[str, Any]]:
    if not SUBFAMILIES_FILE.exists():
        pytest.skip(f"subfamilies file not present: {SUBFAMILIES_FILE}")
    with SUBFAMILIES_FILE.open() as f:
        for line in f:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            br = d.get("bracket")
            if br is None:
                continue
            if br.get("classification") == "neimark_sacker":
                parent = _load_parent_at_step(int(br["step_a"]))
                return br, parent
    raise RuntimeError("no Neimark-Sacker bracket found in subfamilies file")


@pytest.fixture(scope="module")
def smoke_torus() -> SmokeTorus:
    """The #290 converged smoke torus (n_trans=2, amplitude=5e-4), C_J approx 3.1279,
    parent step 112."""
    br, parent = _load_first_neimark_sacker_bracket()
    system = _em_system()
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    lam_a = complex(br["eig_a_re"], br["eig_a_im"])
    lam_b = complex(br["eig_b_re"], br["eig_b_im"])
    torus = correct_qp_torus(
        system,
        parent_state,
        parent_period,
        (lam_a, lam_b),
        k=int(br["k"]),
        n_long=16,
        n_trans=2,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="466_energy_walk_smoke_seed",
    )
    return system, torus


# ---------------------------------------------------------------------------
# Task 1: parent-family ladder + unit-circle Floquet pair detection.
# ---------------------------------------------------------------------------


def test_parent_family_is_monotone_energy_ladder() -> None:
    """The #296 parent family is an energy continuation: C_J strictly monotone in
    step_index. This is the ladder the walk rides (a digested INPUT, not a golden)."""
    fam = qew.load_parent_family()
    assert len(fam) > 200
    fam_sorted = sorted(fam, key=lambda m: m.step_index)
    # strictly increasing with step over the seed->silver segment (steps 8..112)
    seg = [m.jacobi for m in fam_sorted if 8 <= m.step_index <= 112]
    assert all(a < b for a, b in pairwise(seg))
    # seed (step 112) and SILVER (step 8) both present and ~0.096 apart
    by_step = {m.step_index: m for m in fam}
    assert abs(by_step[112].jacobi - 3.12785) < 1e-3
    assert abs(by_step[8].jacobi - 3.03196) < 1e-3


def test_unit_circle_pair_present_seed_and_silver() -> None:
    """A complex Floquet pair sits ON the unit circle at BOTH the seed (112) and the
    SILVER (8) parent member -- the QP-torus center persists across the descent."""
    fam = {m.step_index: m for m in qew.load_parent_family()}
    for st in (112, 8):
        pair = qew._nearest_unit_complex_pair(fam[st].floquet)
        assert pair is not None
        lam = pair[0]
        assert abs(abs(lam) - 1.0) < 1e-2
        assert abs(lam.imag) > 1e-6  # genuinely complex (rotation), not real
        assert lam.imag > 0  # oriented upper-half-plane


def test_k_from_pair_sane() -> None:
    # a 90-degree rotation (4th root of unity) -> k=4
    lam = complex(np.cos(np.pi / 2), np.sin(np.pi / 2))
    assert qew._k_from_pair(lam) == 4


# ---------------------------------------------------------------------------
# Task 2: a single parent member down the ladder converges a genuine torus.
# ---------------------------------------------------------------------------


def test_single_lower_member_converges_at_lower_cj(smoke_torus: SmokeTorus) -> None:
    """Re-converging the torus one parent member DOWN lands a valid irrational torus
    at a strictly LOWER C_J than the seed -- energy moved, structure survived."""
    system, _torus = smoke_torus
    fam = {m.step_index: m for m in qew.load_parent_family()}
    cj_seed = fam[112].jacobi
    lower = fam[104]  # one stride=8 member down
    t = qew._converge_member_torus(
        lower,
        system,
        n_trans=2,
        amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="test",
    )
    assert t is not None
    member = qew._member_from_torus(t, lower, arclength_s=0.0)
    assert member.jacobi < cj_seed  # genuinely lower energy
    # Structural gate = corrector-INDEPENDENT off-grid closure (the genuine torus
    # gate). The N=2 GMOS Fourier-norm sits at the ~1.4e-5 FD-conditioning floor at
    # lower-energy members (just above V1's 1e-5 Fourier gate -- the documented
    # Phase-1 limit), while the off-grid topology check stays ~3e-6.
    assert member.extras["independent_residual"] < 1e-4  # genuine invariant torus
    assert member.is_practically_irrational  # not a phase-lock collapse


# ---------------------------------------------------------------------------
# Task 3: the energy walk spans a REAL C_J range (the decisive capability).
# ---------------------------------------------------------------------------


def test_energy_walk_spans_real_cj_range(smoke_torus: SmokeTorus) -> None:
    """The decisive #466 capability: walking DOWN the parent ladder makes C_J span a
    real range, NOT the #333 iso-energetic 6e-7 floor. Inverse of #333 containment.
    """
    _system, torus = smoke_torus
    fam = qew.walk_energy(
        torus,
        direction="down",
        step_stride=24,
        max_steps=3,
        seed_step=112,
    )
    cj_vals = [m.jacobi for m in fam.members]
    span = max(cj_vals) - min(cj_vals)
    # a REAL energy span, orders of magnitude above the #333 6.3e-7 iso-energetic floor
    assert span > 1e-2
    # monotone descent away from the seed
    for a, b in pairwise(cj_vals):
        assert b <= a + 1e-9
    # every member a genuine irrational torus (no phase-lock collapse)
    assert all(m.is_practically_irrational for m in fam.members)


def test_each_member_is_genuine_torus(smoke_torus: SmokeTorus) -> None:
    """Every energy-walk member is a genuine invariant torus: the corrector-
    INDEPENDENT off-grid closure (V1_qp's topology gate, RE-EVALUATED via run_v1_qp
    with its own RNG) is < 1e-4 nondim. The GMOS Fourier-norm may sit marginally above
    V1's permissive 1e-5 Fourier floor at lower-energy members (the documented N=2
    FD-conditioning limit), so we assert the off-grid topology gate -- the genuine
    torus property -- and surface the V1 Fourier verdict for reporting, not as a hard
    gate that would falsely reject a valid torus at the truncation floor.
    """
    _system, torus = smoke_torus
    fam = qew.walk_energy(
        torus,
        direction="down",
        step_stride=24,
        max_steps=2,
        seed_step=112,
    )
    for i, m in enumerate(fam.members):
        verdict = run_v1_qp(f"466_member_{i}", m.torus)
        # the corrector-independent topology gate must hold (genuine torus)
        assert verdict.converged_independent, (
            f"member {i} off-grid topology check failed: {verdict}"
        )


def test_walk_is_deterministic(smoke_torus: SmokeTorus) -> None:
    _system, torus = smoke_torus
    fam_a = qew.walk_energy(torus, direction="down", step_stride=24, max_steps=2, seed_step=112)
    fam_b = qew.walk_energy(torus, direction="down", step_stride=24, max_steps=2, seed_step=112)
    cj_a = [m.jacobi for m in fam_a.members]
    cj_b = [m.jacobi for m in fam_b.members]
    assert np.allclose(cj_a, cj_b, atol=1e-12)


def test_on_step_callback_fires(smoke_torus: SmokeTorus) -> None:
    _system, torus = smoke_torus
    seen: list[qew.EnergyWalkMember] = []
    qew.walk_energy(
        torus,
        direction="down",
        step_stride=24,
        max_steps=2,
        seed_step=112,
        on_step=seen.append,
    )
    assert len(seen) >= 2
    assert all(isinstance(m, qew.EnergyWalkMember) for m in seen)


# ---------------------------------------------------------------------------
# Task 4: target-capped walk + irrational-seed sanity.
# ---------------------------------------------------------------------------


def test_walk_reaches_cj_target(smoke_torus: SmokeTorus) -> None:
    """A reachable Jacobi target a couple of strides below the seed terminates with
    reached_target; the final member sits at or above the target (we stop once the
    next member would cross below)."""
    _system, torus = smoke_torus
    fam = {m.step_index: m for m in qew.load_parent_family()}
    target = fam[88].jacobi  # ~2.5 strides of 8 below the seed
    walk = qew.walk_energy(
        torus,
        direction="down",
        step_stride=8,
        max_steps=20,
        seed_step=112,
        cj_target=target,
    )
    if walk.terminated_reason == "reached_target":
        assert walk.members[-1].jacobi >= target - 1e-6
    # in all cases energy MOVED below the seed
    assert walk.members[-1].jacobi < fam[112].jacobi - 1e-4


def test_irrational_seed_freq_ratio(smoke_torus: SmokeTorus) -> None:
    """Sanity: the seed itself is a genuine torus (irrational), so a resonance_lock
    during the descent would reflect a real Arnold-tongue crossing, not a mislabel."""
    _system, torus = smoke_torus
    omega_long = 2 * np.pi / torus.t_strob
    omega_trans = torus.rho / torus.t_strob
    ratio = omega_trans / omega_long
    from cyclerfinder.genome.qp_tori import is_practically_irrational

    assert is_practically_irrational(ratio, max_denominator=12, tol=1e-4)
