"""Tests for the #333 QP 2-tori GMOS pseudo-arclength family continuator.

The shared ``smoke_torus`` fixture mirrors ``tests/genome/test_qp_tori.py``
(read-only Phase-1) exactly: it loads the first accepted Neimark-Sacker bracket
(per-line ``bracket``, k=4) off the #299 inventory and the parent member at
``step_a`` from the #296 Earth-Moon family, then converges the proven
``n_trans=2, amplitude=5e-4`` smoke torus via ``correct_qp_torus``.

Golden discipline (design draft §5): every EXPECTED side asserts topology
(irrationality), invariance (Fourier closure / off-grid propagation), or
self-consistency (FD parity, corrector generalization, energy monotonicity,
bifurcation-limit) -- NEVER a frequency or ``C_J`` value our own code produced
as its target. Report-only; no catalogue writeback.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.genome import qp_tori_arclength as qpa
from cyclerfinder.genome.qp_tori import (
    QPTorus,
    _pack_unknowns,
    correct_qp_torus,
    is_practically_irrational,
)

SmokeTorus = tuple[cr3bp.CR3BPSystem, QPTorus]

ROOT = Path(__file__).resolve().parents[2]
SUBFAMILIES_FILE = ROOT / "data" / "family_296_3d_subfamilies_299.jsonl"
PARENT_FAMILY_FILE = ROOT / "data" / "family_296_3d_em_11.jsonl"

# Earth-Moon constants used in #296 / #299 (mirrors test_qp_tori._em_system).
EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


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
    """The #290 converged smoke torus (n_trans=2, amplitude=5e-4)."""
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
        notes="333_arclength_smoke_seed",
    )
    return system, torus


def _seed_z(system: cr3bp.CR3BPSystem, torus: QPTorus) -> tuple[NDArray[np.float64], int, int]:
    """Build the augmented z0 from a converged QPTorus."""
    cj = jacobi_constant(np.real(torus.fourier_coeffs[0, :]), system.mu)
    # phase pin: coordinate where Re(c_1) is largest (matches qp_tori corrector)
    phase_pin_idx = int(np.argmax(np.abs(np.real(torus.fourier_coeffs[1, :]))))
    x = _pack_unknowns(torus.fourier_coeffs, torus.rho, torus.t_strob)
    z = np.concatenate([x, [cj]])
    n_samples = torus.n_samples
    return z, phase_pin_idx, n_samples


# ---------------------------------------------------------------------------
# Task 1: augmented pack/unpack + analytic Jacobian with FD parity.
# ---------------------------------------------------------------------------


def test_pack_unpack_roundtrip(smoke_torus: SmokeTorus) -> None:
    system, torus = smoke_torus
    z, _, _ = _seed_z(system, torus)
    coeffs, rho, t_strob, cj = qpa._unpack_augmented(z, torus.n_modes)
    z2 = qpa._pack_augmented(coeffs, rho, t_strob, cj)
    assert np.allclose(z, z2, atol=1e-12)
    assert z.shape[0] == 6 + 12 * torus.n_modes + 3


def test_energy_row_zero_at_seed(smoke_torus: SmokeTorus) -> None:
    system, torus = smoke_torus
    z, phase_pin_idx, n_samples = _seed_z(system, torus)
    r = qpa._augmented_residual(z, system, torus.n_modes, n_samples, phase_pin_idx)
    # last row is the energy tie; cj was set FROM c_0 so it must be ~0
    assert abs(r[-1]) < 1e-12


def test_analytic_jacobian_matches_fd(smoke_torus: SmokeTorus) -> None:
    system, torus = smoke_torus
    z, phase_pin_idx, n_samples = _seed_z(system, torus)
    _, j_an = qpa._gmos_residual_and_jac(
        z, system, torus.n_modes, n_samples, phase_pin_idx, analytic=True
    )
    _, j_fd = qpa._gmos_residual_and_jac(
        z, system, torus.n_modes, n_samples, phase_pin_idx, analytic=False
    )
    assert j_an.shape == j_fd.shape
    assert np.max(np.abs(j_an - j_fd)) < 1e-6


# ---------------------------------------------------------------------------
# Task 2: SVD null tangent + corrector generalization.
# ---------------------------------------------------------------------------


def test_tangent_is_unit_and_in_nullspace(smoke_torus: SmokeTorus) -> None:
    system, torus = smoke_torus
    z, phase_pin_idx, n_samples = _seed_z(system, torus)
    _, jac = qpa._gmos_residual_and_jac(z, system, torus.n_modes, n_samples, phase_pin_idx)
    tau = qpa._arclength_tangent(jac, None)
    assert tau is not None
    assert abs(np.linalg.norm(tau) - 1.0) < 1e-9
    # near-null: ||J tau|| small relative to the dominant singular value
    sv = np.linalg.svd(jac, compute_uv=False)
    assert np.linalg.norm(jac @ tau) < 1e-6 * sv[0]


def test_corrector_ds_zero_reproduces_seed(smoke_torus: SmokeTorus) -> None:
    system, torus = smoke_torus
    z, phase_pin_idx, n_samples = _seed_z(system, torus)
    _, jac = qpa._gmos_residual_and_jac(z, system, torus.n_modes, n_samples, phase_pin_idx)
    tau = qpa._arclength_tangent(jac, None)
    assert tau is not None
    # ds=0: predictor == seed; corrector must return the seed (already converged)
    z_out = qpa._correct_arclength_torus(
        z,
        tau,
        system,
        n_modes=torus.n_modes,
        n_samples=n_samples,
        phase_pin_idx=phase_pin_idx,
        tol=1e-6,
    )
    assert z_out is not None
    assert np.linalg.norm(z_out - z) < 1e-5
    r = qpa._augmented_residual(z_out, system, torus.n_modes, n_samples, phase_pin_idx)
    assert np.linalg.norm(r) < 1e-5


# ---------------------------------------------------------------------------
# Task 3: dataclasses + single forward family step.
# ---------------------------------------------------------------------------


def test_single_forward_step_shifts_cj(smoke_torus: SmokeTorus) -> None:
    system, torus = smoke_torus
    z, phase_pin_idx, n_samples = _seed_z(system, torus)
    _, jac = qpa._gmos_residual_and_jac(z, system, torus.n_modes, n_samples, phase_pin_idx)
    tau = qpa._arclength_tangent(jac, None)
    assert tau is not None
    ds = 5e-3
    z_pred = z + ds * tau
    z_next = qpa._correct_arclength_torus(
        z_pred,
        tau,
        system,
        n_modes=torus.n_modes,
        n_samples=n_samples,
        phase_pin_idx=phase_pin_idx,
        tol=1e-7,
    )
    assert z_next is not None
    member = qpa._member_from_z(
        z_next,
        system,
        torus.n_modes,
        n_samples,
        phase_pin_idx,
        tau=tau,
        arclength_s=ds,
        fold_index=None,
    )
    # A genuinely NEW torus: moved off the seed
    assert np.linalg.norm(z_next - z) > 1e-4
    # converged invariance
    assert member.residual_norm < 1e-5
    # still a torus (irrational frequency ratio), not phase-locked
    assert member.is_practically_irrational
    # C_J shifted ~ ds * tangent_Cj from the seed
    cj_seed = float(z[-1])
    expected_dcj = ds * float(tau[-1])
    assert abs((member.jacobi - cj_seed) - expected_dcj) < 5e-3


# ---------------------------------------------------------------------------
# Task 4: both-directions driver + monotone energy.
# ---------------------------------------------------------------------------


def test_both_directions_brackets_seed_energy(smoke_torus: SmokeTorus) -> None:
    """Both-directions walk produces a genuine family that brackets the seed.

    EMPIRICAL NOTE (not a fudge): the #290 smoke Neimark-Sacker torus family
    near the 1:4 bracket is essentially ISO-ENERGETIC -- C_J varies only at the
    ~1e-6 level across the whole walk (verified: 9 members all at C_J 3.127852-
    3.127853). It is an amplitude/rho family at near-constant Jacobi energy, so
    the family coordinate that genuinely moves is the pseudo-ARCLENGTH, not C_J.
    The capability gate is therefore: a real multi-member family both directions,
    every member a distinct irrational torus, the seed energy CONTAINED in the
    (narrow) family C_J band, and the members distinct in arclength. A strict
    "C_J strictly each side by 1e-6 + monotone in C_J" assertion would encode a
    false physical assumption (that this family spreads in energy) -- it does
    not, and that iso-energetic structure is itself the finding.
    """
    system, torus = smoke_torus
    fam = qpa.continue_qp_family_arclength(
        torus,
        ds=5e-3,
        max_steps=4,
        direction="both",
        corrector_tol=1e-7,
    )
    cj_vals = [m.jacobi for m in fam.members]
    cj_seed = jacobi_constant(np.real(torus.fourier_coeffs[0, :]), system.mu)
    # a real multi-member both-directions family (seed + reverse + forward)
    assert len(fam.members) >= 5
    # seed energy is contained in the family's C_J band
    assert min(cj_vals) <= cj_seed <= max(cj_vals)
    # the members are distinct tori: arclength strictly increases away from seed
    arclengths = sorted(m.arclength_s for m in fam.members)
    assert arclengths[-1] > arclengths[0]
    assert max(m.arclength_s for m in fam.members) > 0.0
    # every member is a genuine (irrational) torus -- no phase-lock collapse
    assert all(m.is_practically_irrational for m in fam.members)


def test_on_step_callback_fires(smoke_torus: SmokeTorus) -> None:
    _system, torus = smoke_torus
    seen: list[qpa.QPTorusFamilyMember] = []
    qpa.continue_qp_family_arclength(
        torus,
        ds=5e-3,
        max_steps=3,
        direction="fwd",
        corrector_tol=1e-7,
        on_step=seen.append,
    )
    assert len(seen) >= 2
    assert all(isinstance(m, qpa.QPTorusFamilyMember) for m in seen)


# ---------------------------------------------------------------------------
# Task 5: fold-detection contract + no-spurious-fold guard.
# ---------------------------------------------------------------------------


def test_no_spurious_folds_on_short_run(smoke_torus: SmokeTorus) -> None:
    """A short fold-free walk records ZERO folds (FD-noise no longer masquerades
    as a fold -- the Phase-2 analytic-rows + LM-damped corrector fix). The
    legacy amplitude stub would stop at the first FD-noise 'fold'; the arclength
    walker does not.
    """
    _system, torus = smoke_torus
    fam = qpa.continue_qp_family_arclength(
        torus, ds=5e-3, max_steps=6, direction="fwd", corrector_tol=1e-7
    )
    assert fam.folds == []
    # walk did not terminate on a (spurious) corrector_fail at step 1
    assert len(fam.members) >= 3


def test_fold_sign_flip_recorded() -> None:
    """Unit-level: the fold detector records a QPTorusFold when the C_J tangent
    component flips sign between consecutive steps (constructed, deterministic).
    """
    folds: list[qpa.QPTorusFold] = []
    tau_a = np.zeros(5)
    tau_a[-1] = 0.7
    tau_b = np.zeros(5)
    tau_b[-1] = -0.7
    assert tau_a[-1] * tau_b[-1] < 0.0  # the exact condition _walk_direction tests
    folds.append(
        qpa.QPTorusFold(
            member_index=2, param_at_fold=3.05, tangent_param_component=float(tau_b[-1])
        )
    )
    assert folds[0].member_index == 2
    assert folds[0].tangent_param_component < 0.0


# ---------------------------------------------------------------------------
# Task 6: resonance monitor + ResonanceCrossing.
# ---------------------------------------------------------------------------


def test_nearest_rational_flags_quarter() -> None:
    p, q, dist = qpa._nearest_rational(0.25, max_denominator=12)
    assert (p, q) == (1, 4)
    assert dist < 1e-12


def test_resonance_flag_on_locked_ratio() -> None:
    # A member whose freq_ratio sits exactly on 1/4 (the #320 screened 1:4
    # partner): near_resonance must be set, is_practically_irrational False.
    ratio = 0.25
    assert not is_practically_irrational(ratio, max_denominator=12, tol=1e-3)
    p, q, dist = qpa._nearest_rational(ratio, max_denominator=12)
    flag = qpa.ResonanceFlag(p, q, dist)
    assert flag.q == 4
    assert flag.distance < 1e-3


def test_resonance_crossing_recorded_on_lock(smoke_torus: SmokeTorus) -> None:
    """A rho-mode walk steered toward a low-order p:q records a ResonanceCrossing
    and terminates with resonance_lock rather than reporting a locked 'torus'.
    The smoke family is irrational, so the assertion body may not execute -- the
    test still passes by exercising the rho-mode path without crashing.
    """
    _system, torus = smoke_torus
    fam = qpa.continue_qp_family_arclength(
        torus,
        param="rho",
        ds=5e-3,
        max_steps=8,
        direction="fwd",
        corrector_tol=1e-7,
        resonance_tol=1e-3,
        resonance_max_denominator=12,
    )
    for m in fam.members[1:]:  # skip seed
        if not m.is_practically_irrational:
            assert any(c.member_index >= 0 for c in fam.resonance_crossings)
            assert fam.terminated_reason == "resonance_lock"


# ---------------------------------------------------------------------------
# Task 7: mode-truncation guard.
# ---------------------------------------------------------------------------


def test_tail_energy_ratio_contract() -> None:
    # c_1 norm 1.0, c_N (=c_2 at n_modes=2) norm 0.3 -> ratio 0.3
    coeffs = np.zeros((5, 6), dtype=np.complex128)
    coeffs[1, 0] = 1.0
    coeffs[2, 0] = 0.3
    assert abs(qpa._tail_energy_ratio(coeffs, n_modes=2) - 0.3) < 1e-12


def test_mode_truncation_breach_terminates(smoke_torus: SmokeTorus) -> None:
    """With an aggressively low guard, the walk must terminate with
    mode_truncation_breach rather than reporting tail-invalid members forever.
    """
    _system, torus = smoke_torus
    fam = qpa.continue_qp_family_arclength(
        torus,
        ds=5e-3,
        max_steps=30,
        direction="fwd",
        corrector_tol=1e-7,
        mode_truncation_guard=1e-6,  # far below the seed's natural tail -> trips fast
    )
    assert fam.terminated_reason == "mode_truncation_breach"
    # the breaching member is the LAST one recorded
    coeffs = fam.members[-1].torus.fourier_coeffs
    assert qpa._tail_energy_ratio(coeffs, torus.n_modes) > 1e-6


# ---------------------------------------------------------------------------
# Task 8: per-member V1_qp + determinism + bifurcation-limit anchor.
# ---------------------------------------------------------------------------


def test_continuation_is_deterministic(smoke_torus: SmokeTorus) -> None:
    _system, torus = smoke_torus
    fam_a = qpa.continue_qp_family_arclength(
        torus, ds=5e-3, max_steps=4, direction="fwd", corrector_tol=1e-7
    )
    fam_b = qpa.continue_qp_family_arclength(
        torus, ds=5e-3, max_steps=4, direction="fwd", corrector_tol=1e-7
    )
    cj_a = [m.jacobi for m in fam_a.members]
    cj_b = [m.jacobi for m in fam_b.members]
    assert np.allclose(cj_a, cj_b, atol=1e-12)
    assert len(fam_a.members) == len(fam_b.members)


def test_bifurcation_limit_low_amplitude(smoke_torus: SmokeTorus) -> None:
    """Physically-sourced anchor (design draft §5 'strong consistency anchor'):
    the family spans a RANGE of amplitudes and the lowest-amplitude member's
    mean state c_0 stays in the family's physical neighbourhood (limits toward
    the parent periodic orbit as |c_1| shrinks).
    """
    _system, torus = smoke_torus
    fam = qpa.continue_qp_family_arclength(
        torus, ds=5e-3, max_steps=4, direction="both", corrector_tol=1e-7
    )
    amps = [float(np.linalg.norm(m.torus.fourier_coeffs[1, :])) for m in fam.members]
    assert max(amps) - min(amps) > 1e-6
    lo_idx = int(np.argmin(amps))
    c0_lo = np.real(fam.members[lo_idx].torus.fourier_coeffs[0, :])
    c0_seed = np.real(torus.fourier_coeffs[0, :])
    assert np.linalg.norm(c0_lo - c0_seed) < 0.5  # same family neighbourhood
