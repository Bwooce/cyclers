"""Permanent numerical math PROPERTY tests (promoted from the 2026-06-11/12 review probes).

PROPERTY / CONSISTENCY tests, NOT sourced golden validation
-----------------------------------------------------------
Every assertion here checks a *mathematical identity* between two pieces of
this repo's own code (or between code and a finite-difference / round-trip
reconstruction of itself). Nothing on any EXPECTED side traces to a published
source, so these tests can never substitute for the sourced golden suites —
they catch internal inconsistencies (a wrong variational A-matrix, a Lambert
velocity-reconstruction bug, a non-orthogonal frame rotation), not shared
systematic errors. See ``docs/notes/2026-06-11-project-review-results.md``
(the math-verification probe suites these are promoted from).

The four probe families promoted to permanence (task #204):

1. **FD-vs-STM** — columns of the CR3BP state-transition matrix from
   :func:`cyclerfinder.core.cr3bp.propagate(..., with_stm=True)` must match
   central finite differences of the plain flow map (cross-checks the
   variational A-matrix in ``cr3bp_stm_eom`` against ``cr3bp_eom``).
2. **∇C·f = 0** — the Jacobi-constant gradient is orthogonal to the CR3BP
   flow field at random states (C is conserved along the flow). Checked two
   ways: an FD gradient of ``jacobi_constant`` (pure cross-module check, no
   test-local algebra) and the analytic gradient (exact identity, tight tol).
3. **Lambert-vs-Kepler truth** — solve Lambert between two points generated
   by propagating a known Kepler orbit; the recovered velocities must equal
   the Kepler-propagated truth velocities (short way, long way, and a
   multi-revolution branch).
4. **Frame orthonormality** — the rotation operators implied by the frames /
   ephemeris layer satisfy R Rᵀ = I and det R = +1 across random epochs
   (uniform synodic frame, dynamic frame on the circular and
   inclined-circular backends, the tilted Rodrigues omega-vector path, and
   the inclined backend's constant per-planet plane rotations). The DE440
   astropy backend's constant obliquity rotation is exercised by
   ``tests/verify/test_ephemeris_crosscheck.py`` and is not duplicated here.

All randomness is seeded (``np.random.default_rng`` with fixed seeds) so the
suite is deterministic; the whole file runs in a few seconds.
"""

from __future__ import annotations

from math import pi, sqrt
from typing import Protocol

import numpy as np
import pytest
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris, _InclinedCircularBackend, inclined_planets
from cyclerfinder.core.frames import (
    synodic_omega,
    to_rotating,
    to_rotating_dynamic,
    to_rotating_omega_vec,
)
from cyclerfinder.core.kepler import coe_to_rv
from cyclerfinder.core.kepler import propagate as kepler_propagate
from cyclerfinder.core.lambert import lambert

Vec6 = NDArray[np.float64]

# Test-local mass ratio of Earth-Moon magnitude. This is a PROPERTY suite:
# the identities below hold for any mu in (0, 1/2), so no registry / published
# value is needed (and none is claimed).
_MU_EM_LIKE: float = 0.0121505856

_AU_KM: float = 1.495978707e8


def _sample_cr3bp_states(rng: np.random.Generator, n: int, mu: float) -> list[Vec6]:
    """Random nondimensional CR3BP states staying >= 0.3 from both primaries.

    Positions in a [-1.5, 1.5]^2 x [-0.45, 0.45] box (z compressed), bounded
    away from the r1/r2 singularities so DOP853 and the finite differences
    stay well-conditioned; velocities O(0.5) nondimensional.
    """
    states: list[Vec6] = []
    p1 = np.array([-mu, 0.0, 0.0])
    p2 = np.array([1.0 - mu, 0.0, 0.0])
    while len(states) < n:
        pos = rng.uniform(-1.5, 1.5, size=3)
        pos[2] *= 0.3
        if np.linalg.norm(pos - p1) < 0.3 or np.linalg.norm(pos - p2) < 0.3:
            continue
        vel = rng.uniform(-0.8, 0.8, size=3)
        states.append(np.concatenate([pos, vel]).astype(np.float64))
    return states


# ---------------------------------------------------------------------------
# 1. FD-vs-STM: variational equations vs central finite differences
# ---------------------------------------------------------------------------


def test_cr3bp_stm_columns_match_central_finite_differences() -> None:
    """STM columns == central FD of the flow map, rel tol 1e-6 per column.

    Step h = 3e-6 nondim sits at the truncation/roundoff sweet spot for the
    rtol=1e-12 DOP853 flow map (probe-measured worst column error ~4e-10;
    the 1e-6 tolerance leaves >1000x headroom while still catching any term
    error in the variational A-matrix, which would show up at O(1)).
    """
    rng = np.random.default_rng(20260612)
    system = cr3bp.CR3BPSystem(
        mu=_MU_EM_LIKE,
        primary="test-primary",
        secondary="test-secondary",
        l_km=384400.0,
        t_s=375190.0,
    )
    t_nd = 1.0
    h = 3.0e-6
    worst = 0.0
    for s0 in _sample_cr3bp_states(rng, 5, system.mu):
        arc = cr3bp.propagate(system, s0, t_nd, with_stm=True)
        assert arc.stm is not None
        for j in range(6):
            sp = s0.copy()
            sm = s0.copy()
            sp[j] += h
            sm[j] -= h
            fp = cr3bp.propagate(system, sp, t_nd).state_f
            fm = cr3bp.propagate(system, sm, t_nd).state_f
            fd_col = (fp - fm) / (2.0 * h)
            rel = float(
                np.linalg.norm(fd_col - arc.stm[:, j])
                / max(float(np.linalg.norm(arc.stm[:, j])), 1.0)
            )
            worst = max(worst, rel)
    assert worst < 1.0e-6, f"STM column vs central FD: worst rel error {worst:.3e} >= 1e-6"


# ---------------------------------------------------------------------------
# 2. grad(C) . f = 0: Jacobi gradient orthogonal to the CR3BP flow field
# ---------------------------------------------------------------------------


def _fd_jacobi_gradient(state: Vec6, mu: float, h: float = 1.0e-6) -> Vec6:
    grad = np.zeros(6, dtype=np.float64)
    for j in range(6):
        sp = state.copy()
        sm = state.copy()
        sp[j] += h
        sm[j] -= h
        grad[j] = (cr3bp.jacobi_constant(sp, mu) - cr3bp.jacobi_constant(sm, mu)) / (2.0 * h)
    return grad


def _analytic_jacobi_gradient(state: Vec6, mu: float) -> Vec6:
    """grad of C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2 (test-local algebra)."""
    x, y, z, vx, vy, vz = (float(v) for v in state)
    r1 = sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    r1c, r2c = r1**3, r2**3
    om1 = 1.0 - mu
    gx = 2.0 * x - 2.0 * om1 * (x + mu) / r1c - 2.0 * mu * (x - 1.0 + mu) / r2c
    gy = 2.0 * y - 2.0 * om1 * y / r1c - 2.0 * mu * y / r2c
    gz = -2.0 * om1 * z / r1c - 2.0 * mu * z / r2c
    return np.array([gx, gy, gz, -2.0 * vx, -2.0 * vy, -2.0 * vz], dtype=np.float64)


@pytest.mark.parametrize("mu", [_MU_EM_LIKE, 0.3, 1.0e-4])
def test_jacobi_gradient_orthogonal_to_flow_fd(mu: float) -> None:
    """FD grad of jacobi_constant dotted with cr3bp_eom vanishes (tol 1e-8).

    Pure cross-module consistency: the gradient is built from
    ``jacobi_constant`` itself (no test-local algebra), so a sign/term
    mismatch between the conserved quantity and the EOM cannot hide.
    Normalisation: |grad . f| / max(|grad| |f|, 1). FD floor ~1e-10.
    """
    rng = np.random.default_rng(42)
    worst = 0.0
    for state in _sample_cr3bp_states(rng, 8, mu):
        grad = _fd_jacobi_gradient(state, mu)
        flow = cr3bp.cr3bp_eom(0.0, state, mu)
        dot = float(np.dot(grad, flow))
        scale = max(float(np.linalg.norm(grad)) * float(np.linalg.norm(flow)), 1.0)
        worst = max(worst, abs(dot) / scale)
    assert worst < 1.0e-8, f"FD grad(C).f: worst normalised |dot| {worst:.3e} >= 1e-8"


@pytest.mark.parametrize("mu", [_MU_EM_LIKE, 0.3, 1.0e-4])
def test_jacobi_gradient_orthogonal_to_flow_analytic(mu: float) -> None:
    """Analytic grad(C) . f == 0 to near machine precision (exact identity).

    The Coriolis accelerations are workless and the pseudo-potential terms
    cancel exactly, so the only residual is float roundoff (probe-measured
    worst ~1.4e-16 normalised; tolerance 1e-13 leaves cross-platform slack).
    """
    rng = np.random.default_rng(42)
    worst = 0.0
    for state in _sample_cr3bp_states(rng, 8, mu):
        grad = _analytic_jacobi_gradient(state, mu)
        flow = cr3bp.cr3bp_eom(0.0, state, mu)
        dot = float(np.dot(grad, flow))
        scale = max(float(np.linalg.norm(grad)) * float(np.linalg.norm(flow)), 1.0)
        worst = max(worst, abs(dot) / scale)
    assert worst < 1.0e-13, f"analytic grad(C).f: worst normalised |dot| {worst:.3e} >= 1e-13"


# ---------------------------------------------------------------------------
# 3. Lambert-vs-Kepler round-trip truth
# ---------------------------------------------------------------------------


def _transfer_angle_prograde(r1: Vec6, r2: Vec6) -> float:
    """Transfer angle measured CCW about +z (prograde), in [0, 2*pi)."""
    r1n = float(np.linalg.norm(r1))
    r2n = float(np.linalg.norm(r2))
    cos_dnu = float(np.clip(np.dot(r1, r2) / (r1n * r2n), -1.0, 1.0))
    dnu = float(np.arccos(cos_dnu))
    cross_z = float(r1[0] * r2[1] - r1[1] * r2[0])
    if cross_z < 0.0:
        dnu = 2.0 * pi - dnu
    return dnu


def _lambert_kepler_roundtrip_worst(
    rng: np.random.Generator, frac_lo: float, frac_hi: float, want_long_way: bool
) -> float:
    """Worst per-component |v - v_truth| (km/s) over 5 seeded Kepler legs.

    Truth: an elliptic heliocentric orbit (e <= 0.4, prograde, planar) is
    propagated by the repo's Kepler propagator; Lambert between the two
    positions must recover both endpoint velocities. Time fractions
    [0.10, 0.22] of the period guarantee transfer angle < pi for e <= 0.4
    (the fastest pi-sweep, periapsis-centred at e = 0.4, takes 0.252 of a
    period); fractions [0.78, 0.90] guarantee > pi by symmetry. Each case's
    actual geometry is asserted as a precondition.
    """
    worst = 0.0
    n_cases = 0
    while n_cases < 5:
        a_km = _AU_KM * rng.uniform(0.8, 1.8)
        ecc = rng.uniform(0.05, 0.4)
        nu0 = rng.uniform(0.0, 2.0 * pi)
        argp = rng.uniform(0.0, 2.0 * pi)
        r1, v1_truth = coe_to_rv(a_km, ecc, nu0, MU_SUN_KM3_S2, arg_peri_rad=argp)
        period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
        tof_s = rng.uniform(frac_lo, frac_hi) * period_s
        r2, v2_truth = kepler_propagate(r1, v1_truth, tof_s, MU_SUN_KM3_S2)
        dnu = _transfer_angle_prograde(r1, r2)
        # Geometry precondition (also keeps clear of the 180-deg singularity).
        if want_long_way:
            assert dnu > pi + 0.05, f"case not long-way: dnu={dnu:.4f}"
        else:
            assert dnu < pi - 0.05, f"case not short-way: dnu={dnu:.4f}"
        sol = lambert(r1, r2, tof_s, mu=MU_SUN_KM3_S2, prograde=True)[0]
        err = max(
            float(np.max(np.abs(sol.v1 - v1_truth))),
            float(np.max(np.abs(sol.v2 - v2_truth))),
        )
        worst = max(worst, err)
        n_cases += 1
    return worst


def test_lambert_recovers_kepler_velocities_short_way() -> None:
    """Lambert on Kepler-generated endpoints, transfer angle < pi (tol 1e-6 km/s)."""
    worst = _lambert_kepler_roundtrip_worst(
        np.random.default_rng(7), 0.10, 0.22, want_long_way=False
    )
    assert worst < 1.0e-6, f"short-way Lambert vs Kepler truth: worst |dv| {worst:.3e} km/s"


def test_lambert_recovers_kepler_velocities_long_way() -> None:
    """Lambert on Kepler-generated endpoints, transfer angle > pi (tol 1e-6 km/s)."""
    worst = _lambert_kepler_roundtrip_worst(
        np.random.default_rng(11), 0.78, 0.90, want_long_way=True
    )
    assert worst < 1.0e-6, f"long-way Lambert vs Kepler truth: worst |dv| {worst:.3e} km/s"


def test_lambert_multirev_branch_recovers_kepler_velocities() -> None:
    """A 1-rev Lambert branch recovers the Kepler truth for tof = period + arc.

    Deterministic fixture: a = 1.3 AU, e = 0.2 ellipse propagated one full
    revolution plus 0.3 of a period. The direct (n_revs=0) conic connecting
    the same endpoints in that tof is a *different* orbit, so only a 1-rev
    branch can match the truth velocities — this pins the multi-rev branch
    construction, not just the root solve.
    """
    a_km = _AU_KM * 1.3
    ecc = 0.2
    r1, v1_truth = coe_to_rv(a_km, ecc, 0.7, MU_SUN_KM3_S2)
    period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
    tof_s = 1.3 * period_s
    r2, v2_truth = kepler_propagate(r1, v1_truth, tof_s, MU_SUN_KM3_S2)
    sols = lambert(r1, r2, tof_s, mu=MU_SUN_KM3_S2, prograde=True, max_revs=2)
    one_rev = [s for s in sols if s.n_revs == 1]
    assert one_rev, f"no 1-rev branches returned (got {[(s.n_revs, s.branch) for s in sols]})"
    best = min(
        max(
            float(np.max(np.abs(s.v1 - v1_truth))),
            float(np.max(np.abs(s.v2 - v2_truth))),
        )
        for s in one_rev
    )
    assert best < 1.0e-6, f"multi-rev Lambert vs Kepler truth: best 1-rev |dv| {best:.3e} km/s"


# ---------------------------------------------------------------------------
# 4. Frame orthonormality: R R^T = I, det R = +1 across random epochs
# ---------------------------------------------------------------------------

_ORTHO_TOL: float = 1.0e-13


def _assert_rotation(rot: NDArray[np.float64], label: str) -> None:
    ortho_err = float(np.max(np.abs(rot @ rot.T - np.eye(3))))
    det_err = abs(float(np.linalg.det(rot)) - 1.0)
    assert ortho_err < _ORTHO_TOL, f"{label}: ||R R^T - I||_max = {ortho_err:.3e}"
    assert det_err < _ORTHO_TOL, f"{label}: |det R - 1| = {det_err:.3e}"


class _PositionTransform(Protocol):
    def __call__(
        self, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...


def _implied_position_rotation(
    transform: _PositionTransform,
) -> NDArray[np.float64]:
    """Extract the 3x3 position-rotation matrix of a frame transform.

    The position part of every rotating-frame transform here is linear
    (r_rot = R r), so applying the transform to the three basis vectors with
    zero velocity reconstructs R column by column. (The velocity part is
    affine in r via the Coriolis correction and is covered by the round-trip
    tests in tests/core/test_frames*.py.)
    """
    cols = []
    zero_v = np.zeros(3, dtype=np.float64)
    for j in range(3):
        e_j = np.zeros(3, dtype=np.float64)
        e_j[j] = 1.0
        r_rot, _ = transform(e_j, zero_v)
        cols.append(r_rot)
    return np.column_stack(cols)


def test_uniform_rotating_frame_rotation_orthonormal_across_epochs() -> None:
    """Implied R(-omega t) of to_rotating is orthonormal at random epochs."""
    rng = np.random.default_rng(99)
    omega = synodic_omega("E")
    for t_sec in rng.uniform(-50.0 * 365.25 * 86400.0, 50.0 * 365.25 * 86400.0, size=10):
        t_s = float(t_sec)
        rot = _implied_position_rotation(lambda r, v, t=t_s: to_rotating(r, v, t, omega))
        _assert_rotation(rot, f"uniform frame t={t_s:.6e}s")


def test_omega_vec_rodrigues_rotation_orthonormal_across_epochs() -> None:
    """Implied rotation of the tilted omega-vector (Rodrigues) path is orthonormal.

    A genuinely tilted omega vector forces the general Rodrigues branch of
    to_rotating_omega_vec (the pure-z case delegates to the scalar form,
    already covered above).
    """
    rng = np.random.default_rng(101)
    for _ in range(8):
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        if abs(axis[0]) < 1.0e-3 and abs(axis[1]) < 1.0e-3:
            axis[0] += 0.1  # ensure the tilted (non-pure-z) branch is taken
            axis /= np.linalg.norm(axis)
        omega_vec = (axis * synodic_omega("E")).astype(np.float64)
        t_sec = float(rng.uniform(-20.0, 20.0) * 365.25 * 86400.0)
        rot = _implied_position_rotation(
            lambda r, v, t=t_sec, w=omega_vec: to_rotating_omega_vec(r, v, t, w)
        )
        _assert_rotation(rot, f"omega-vec frame axis={axis}, t={t_sec:.6e}s")


@pytest.mark.parametrize("model", ["circular", "inclined-circular"])
def test_dynamic_frame_rotation_orthonormal_across_epochs(model: str) -> None:
    """Implied dynamic-frame rotation is orthonormal at random epochs."""
    rng = np.random.default_rng(303)
    ephem = Ephemeris(model)
    bodies = ("E", "M")
    for t_sec in rng.uniform(0.0, 30.0 * 365.25 * 86400.0, size=8):
        t_s = float(t_sec)
        rot = _implied_position_rotation(
            lambda r, v, t=t_s: to_rotating_dynamic(r, v, t, bodies, ephem)
        )
        _assert_rotation(rot, f"dynamic frame ({model}) t={t_s:.6e}s")


def test_inclined_backend_plane_rotations_orthonormal() -> None:
    """Per-planet orbital-plane -> ecliptic rotations are orthonormal, det +1.

    These are constant in epoch (built once from inc/Omega), so 'across
    epochs' reduces to checking every cached per-body matrix.
    """
    backend = _InclinedCircularBackend(inclined_planets())
    rotations = backend._rotations
    assert rotations, "inclined backend produced no rotation matrices"
    for code, rot in rotations.items():
        _assert_rotation(rot, f"inclined-circular plane rotation for {code!r}")
