"""Tests for the V1 same-model QP-torus gauntlet (#319 Phase 1 Part A).

Sourced golden discipline (per ``feedback_golden_tests_sourced_only``)
---------------------------------------------------------------------
The tests assert on CLOSURE QUALITY + TOPOLOGY, NOT on specific torus
state numbers our own code produced. The acceptable torus inputs are
seeded from a sourced #299 Neimark-Sacker bracket (the same seed
:mod:`tests.genome.test_qp_tori` uses); treating the resulting torus as
the V1 input is honest because V1's whole job is to RE-VERIFY closure
on a candidate, not to re-solve it from scratch.

The floors (V1_QP_FOURIER_FLOOR = 1e-5, V1_QP_INDEPENDENT_FLOOR = 1e-4)
are sourced from Olikara-Howell 2014 / Olikara 2016 -- module constants,
NOT test-tunable.

Test cases
----------
  1. Sourced #299 smoke-test torus: corrector converges, V1_qp passes both
     gates with a clean margin.
  2. Periodic-orbit limit (zero-amplitude torus): the seed at
     ``initial_torus_amplitude=0`` is just the parent fixed point. V1_qp
     should pass trivially -- the torus is "invariant" (it's a fixed
     point of the stroboscopic map).
  3. Corrupted Fourier coefficients: take a converged torus and bit-flip
     its Fourier modes. V1_qp must FAIL (Fourier residual blows up).
  4. Fabrication-guard: the V1_qp floors equal the module constants.
  5. Audit-trail fields: candidate_id and notes survive into the verdict.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v1_qp import (
    V1_QP_FOURIER_FLOOR,
    V1_QP_INDEPENDENT_FLOOR,
    V1_QP_INDEPENDENT_N_SAMPLES,
    V1VerdictQP,
    run_v1_qp,
)
from cyclerfinder.genome.qp_tori import (
    QPTorus,
    _seed_invariant_circle,
    correct_qp_torus,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
SUBFAMILIES_FILE = DATA_DIR / "family_296_3d_subfamilies_299.jsonl"
PARENT_FAMILY_FILE = DATA_DIR / "family_296_3d_em_11.jsonl"

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


def _build_smoke_torus() -> QPTorus:
    """Corrector-converged torus from the sourced #299 Neimark-Sacker bracket.

    Mirrors ``tests.genome.test_qp_tori.test_sourced_neimark_sacker_smoke``
    so the V1 input is exactly the same torus the genome-tier smoke test
    accepts. NO numerical novelty is asserted -- the torus is the
    project's own computation.
    """
    br, parent = _load_first_neimark_sacker_bracket()
    system = _em_system()
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    k = int(br["k"])
    lam_a = complex(br["eig_a_re"], br["eig_a_im"])
    lam_b = complex(br["eig_b_re"], br["eig_b_im"])
    return correct_qp_torus(
        system,
        parent_state,
        parent_period,
        (lam_a, lam_b),
        k=k,
        n_long=16,
        n_trans=2,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="phase1_smoke_test_from_299_bracket (V1_qp input)",
    )


def _build_zero_amplitude_torus() -> QPTorus:
    """Construct a degenerate "torus" at zero amplitude -- the parent
    periodic orbit dressed as a constant invariant circle.

    The corrector sees the constant-mode seed and the amplitude pin at 0
    is degenerate, so we build the QPTorus directly via the seed helper
    rather than running Newton. This is what a phase-locked /
    fully-collapsed torus looks like -- it IS invariant under the
    stroboscopic flow (it's a fixed point).
    """
    br, parent = _load_first_neimark_sacker_bracket()
    system = _em_system()
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    k = int(br["k"])
    n_modes = 2
    arc = cr3bp.propagate(system, parent_state, parent_period, with_stm=True)
    monod = arc.stm
    assert monod is not None
    coeffs_seed, lam_seed, _ = _seed_invariant_circle(
        parent_state, monod, k=k, n_modes=n_modes, amplitude=0.0
    )
    rho_seed = math.atan2(np.imag(lam_seed), np.real(lam_seed))
    t_strob = float(parent_period)
    n_samples = 2 * n_modes + 3
    return QPTorus(
        system=system,
        omega_long=2 * math.pi / t_strob,
        omega_trans=rho_seed / t_strob,
        rho=rho_seed,
        t_strob=t_strob,
        fourier_coeffs=coeffs_seed,
        n_modes=n_modes,
        n_samples=n_samples,
        invariance_residual=0.0,
        independent_closure_residual=0.0,
        converged=True,
        n_iter=0,
        notes="zero_amplitude_periodic_orbit_limit",
    )


# ---------------------------------------------------------------------------
# Test 1: floor fabrication guard.
# ---------------------------------------------------------------------------


def test_v1_qp_floors_are_sourced_constants() -> None:
    """The Fourier-norm floor and independent floor match the module's
    sourced constants -- not silently tunable.

    Per Olikara-Howell 2014, the Fourier-norm floor of 1e-5 sits in the
    middle of the published 1e-4 to 1e-6 truncation-error band for
    moderate ``N``. The independent floor of 1e-4 matches the QPTorus
    default ``independent_tol``.
    """
    assert V1_QP_FOURIER_FLOOR == 1.0e-5
    assert V1_QP_INDEPENDENT_FLOOR == 1.0e-4
    assert V1_QP_INDEPENDENT_N_SAMPLES == 100


# ---------------------------------------------------------------------------
# Test 2: sourced #299 smoke-test torus passes V1_qp.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v1_qp_sourced_smoke_torus_passes() -> None:
    """V1_qp on the #299 Neimark-Sacker-seeded torus.

    The torus is the same one ``test_sourced_neimark_sacker_smoke``
    accepts at the genome tier. V1_qp re-evaluates the Fourier-norm
    invariance and runs an INDEPENDENT off-grid sampler. Both gates
    must hold.
    """
    torus = _build_smoke_torus()
    verdict = run_v1_qp(
        "qp-torus-smoke-299",
        torus,
        notes="sourced #299 Neimark-Sacker bracket; n_modes=2 amp=5e-4",
    )
    assert isinstance(verdict, V1VerdictQP)
    print(
        f"\n[V1_qp smoke] fourier_norm={verdict.invariance_residual_fourier_norm:.3e} "
        f"| off_grid_nondim={verdict.independent_invariance_residual_nondim:.3e} "
        f"| off_grid_km={verdict.independent_residual_km:.3e} "
        f"| n_modes={verdict.n_modes} | passes={verdict.passes_v1_qp}"
    )
    assert verdict.converged_corrector, (
        f"Fourier-norm residual {verdict.invariance_residual_fourier_norm:.3e} "
        f"exceeds V1_qp floor {verdict.v1_floor_invariance:.0e}"
    )
    assert verdict.converged_independent, (
        f"independent off-grid residual {verdict.independent_invariance_residual_nondim:.3e} "
        f"exceeds V1_qp independent floor {verdict.independent_floor_nondim:.0e}"
    )
    assert verdict.passes_v1_qp


# ---------------------------------------------------------------------------
# Test 3: periodic-orbit limit (zero amplitude) passes V1_qp.
# ---------------------------------------------------------------------------


def test_v1_qp_zero_amplitude_limit_passes() -> None:
    """A zero-amplitude "torus" (collapsed to the parent fixed point) is
    invariant under the stroboscopic flow trivially -- the GMOS residual
    at the constant mode is the parent's periodicity residual, which is
    O(integrator tolerance) = ~1e-12 for the project's defaults.

    Degenerate but ADMISSIBLE: this is the limit where a QP-torus
    family approaches its parent periodic orbit. V1_qp does NOT reject
    it; the orbit_class question (cycler vs quasi_cycler) is a Phase 2
    discriminator, not a V1 question.
    """
    torus = _build_zero_amplitude_torus()
    verdict = run_v1_qp(
        "qp-torus-zero-amplitude-limit",
        torus,
        notes="degenerate parent-orbit limit",
    )
    print(
        f"\n[V1_qp zero-amp] fourier_norm={verdict.invariance_residual_fourier_norm:.3e} "
        f"| off_grid_nondim={verdict.independent_invariance_residual_nondim:.3e} "
        f"| passes={verdict.passes_v1_qp}"
    )
    # At zero amplitude the GMOS residual is the integrator's
    # round-trip error on the parent orbit -- well below 1e-5.
    assert verdict.passes_v1_qp


# ---------------------------------------------------------------------------
# Test 4: corrupted Fourier coefficients fail V1_qp.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v1_qp_rejects_corrupted_fourier_coefficients() -> None:
    """Take a converged torus and apply a gross perturbation to its
    Fourier modes -- V1_qp must catch the corruption.

    This is the negative control: a corrector that converged at one
    point in time but whose state has been mutated since (bit-flip in
    storage, deserialization bug, deliberate tamper) emits an inflated
    Fourier-norm residual and V1_qp REFUSES it.
    """
    torus = _build_smoke_torus()
    # Add a 0.05-nondim displacement to ALL n=1 mode components -- huge
    # vs the seed amplitude of 5e-4; guaranteed to blow the Fourier
    # residual past 1e-5.
    bad_coeffs = torus.fourier_coeffs.copy()
    bad_coeffs[1, :] += 0.05 + 0.05j
    # Maintain the reality constraint c_{-n} = conj(c_n)
    n_total = bad_coeffs.shape[0]
    bad_coeffs[n_total - 1, :] = np.conj(bad_coeffs[1, :])
    corrupted = replace(torus, fourier_coeffs=bad_coeffs)
    verdict = run_v1_qp(
        "qp-torus-corrupted",
        corrupted,
        notes="negative control: deliberately corrupted Fourier modes",
    )
    print(
        f"\n[V1_qp corrupted] fourier_norm={verdict.invariance_residual_fourier_norm:.3e} "
        f"| off_grid_nondim={verdict.independent_invariance_residual_nondim:.3e} "
        f"| passes={verdict.passes_v1_qp}"
    )
    assert not verdict.passes_v1_qp, (
        "V1_qp false-positive on corrupted torus: "
        f"fourier={verdict.invariance_residual_fourier_norm:.3e}, "
        f"off_grid={verdict.independent_invariance_residual_nondim:.3e}"
    )


# ---------------------------------------------------------------------------
# Test 5: audit-trail fields survive.
# ---------------------------------------------------------------------------


def test_v1_qp_verdict_carries_audit_fields() -> None:
    """The verdict carries the floors it was held against AND echoes the
    candidate_id and notes -- the audit trail must let a later reader
    reconstruct exactly what passed/failed.
    """
    torus = _build_zero_amplitude_torus()
    verdict = run_v1_qp(
        "audit-trail-check",
        torus,
        notes="audit-trail check",
    )
    assert verdict.candidate_id == "audit-trail-check"
    assert verdict.notes == "audit-trail check"
    assert verdict.v1_floor_invariance == V1_QP_FOURIER_FLOOR
    assert verdict.independent_floor_nondim == V1_QP_INDEPENDENT_FLOOR
    assert verdict.n_off_grid_samples == V1_QP_INDEPENDENT_N_SAMPLES
    assert verdict.n_modes == torus.n_modes


# ---------------------------------------------------------------------------
# Test 6: malformed system rejected.
# ---------------------------------------------------------------------------


def test_v1_qp_rejects_malformed_system() -> None:
    """A torus carrying a CR3BPSystem with l_km=0 must fail loudly --
    the km conversion of the off-grid residual would otherwise emit a
    meaningless zero.
    """
    torus = _build_zero_amplitude_torus()
    bad_system = cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=0.0,
        t_s=EM_T_S,
    )
    bad_torus = replace(torus, system=bad_system)
    with pytest.raises(ValueError, match="invalid CR3BP system"):
        run_v1_qp("bad-system", bad_torus)


# ---------------------------------------------------------------------------
# Test 7: bad caller args.
# ---------------------------------------------------------------------------


def test_v1_qp_rejects_bad_caller_args() -> None:
    """Invalid floors or sample counts fail loudly."""
    torus = _build_zero_amplitude_torus()
    with pytest.raises(ValueError, match="invariance_floor"):
        run_v1_qp("bad", torus, invariance_floor=0.0)
    with pytest.raises(ValueError, match="independent_tol_nondim"):
        run_v1_qp("bad", torus, independent_tol_nondim=-1.0)
    with pytest.raises(ValueError, match="n_off_grid_samples"):
        run_v1_qp("bad", torus, n_off_grid_samples=0)
    with pytest.raises(TypeError, match="QPTorus"):
        run_v1_qp("bad", "not-a-torus")  # type: ignore[arg-type]
