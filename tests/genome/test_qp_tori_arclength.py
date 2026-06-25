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

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.genome import qp_tori_arclength as qpa
from cyclerfinder.genome.qp_tori import _pack_unknowns, correct_qp_torus

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
def smoke_torus():
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


def _seed_z(system, torus):
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


def test_pack_unpack_roundtrip(smoke_torus):
    system, torus = smoke_torus
    z, _, _ = _seed_z(system, torus)
    coeffs, rho, t_strob, cj = qpa._unpack_augmented(z, torus.n_modes)
    z2 = qpa._pack_augmented(coeffs, rho, t_strob, cj)
    assert np.allclose(z, z2, atol=1e-12)
    assert z.shape[0] == 6 + 12 * torus.n_modes + 3


def test_energy_row_zero_at_seed(smoke_torus):
    system, torus = smoke_torus
    z, phase_pin_idx, n_samples = _seed_z(system, torus)
    r = qpa._augmented_residual(z, system, torus.n_modes, n_samples, phase_pin_idx)
    # last row is the energy tie; cj was set FROM c_0 so it must be ~0
    assert abs(r[-1]) < 1e-12


def test_analytic_jacobian_matches_fd(smoke_torus):
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
