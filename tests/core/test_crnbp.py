"""CRNBP (N=5 Jupiter-Europa-Io-Ganymede) core EOM / STM / equilibrium tests (#717).

`#717` (shortlist item 1 of 3 from `#714`'s CCR5BP/CRNBP scoping pass). Mirrors
``tests/core/test_ccr4bp.py``'s gate structure exactly, extended with the new
N=5-specific gates: coupling-term sign lock and the equilibrium-point
("dynamical substitute") positive control against Gilliam's thesis Ch. V
Tables 5-6.

Sign reconciliation (`#717` Step 0) locked in by Gate 2 below: two prior
digests of Negri & Prado (2022)'s Eq. 11 disagreed on the sign of the mutual
coupling term between extra perturbers. Read directly from both source PDFs
(Negri & Prado's own arXiv:2307.10881 Eqs. 7-9, unambiguous; Gilliam's thesis
Eqs. 25-27 AND an independently-appended Gilliam/Bettinger/et al. *Icarus* 429
(2025) 116455 reprint's Eqs. 9-11 -- three independent sources, all agreeing):
the coupling term is ADDED to the direct term with the SAME sign (an inner
PLUS), not subtracted. The
``docs/notes/2026-07-26-digest-negri-prado-2022-crnbp.md`` digest had the
transcription error (an OCR/multi-line-bracket-rendering artifact); the
Gilliam digest was correct. See ``src/cyclerfinder/core/crnbp.py``'s module
docstring for the full derivation.

Equilibrium-point positive control (`#717` Step 2): Gilliam's thesis Ch. V
computes CRNBP "dynamical substitute" equilibrium points (Newton-Raphson
force-balance loci -- NOT periodic orbits or connections) for the
Jupiter-Europa system with Io+Ganymede as the two extra perturbers -- exactly
this module's N=5 configuration. Her Tables 5/6 report ``dx``, ``dy``,
Lagrange-Box-area, and "average [substitute location] vs CR3BP [point]"
(read from her own wording: the DISTANCE between the time-averaged substitute
location and the base CR3BP point, not the mean of pointwise displacements --
the latter reading was checked and does NOT reproduce her numbers).

Precision caveat (documented, not hidden -- the project's own "Gilliam-table
mismatch policy"): our registry-sourced ``mu_io``/``mu_gan`` values differ
from Gilliam's own stated mass-ratio values (from an appended Icarus 2025
reprint's Table 5: mu_io=4.7055e-5, mu_gan=7.8070e-5) at the ~3-5e-4 relative
level -- both trace to the SAME underlying JPL moon masses/SMAs, but through
different mass-vs-GM rounding conventions (her Table 2 uses rounded body
masses in kg; this project's registry uses GM_de440 gravitational
parameters). Additionally, her Ch. V equilibrium tables do not state the
initial synodic phase angles used for Io/Ganymede at the propagation epoch --
the tests below explicitly request ``theta_io0=0.0`` (an "all perturbers
aligned at t=0" convention, the simplest choice and the one her Ch. IV
prose describes as "test case 1"). NOTE (`#723`): the module DEFAULT is no
longer this all-aligned phase -- :func:`crnbp.jupiter_europa_io_ganymede_default`
now defaults to the physically-real Laplace-libration phase ``theta_io0 =
pi - 2*theta_gan0`` (see ``crnbp.py``'s module docstring); the all-aligned
"test case 1" configuration used below is obtained by passing
``theta_io0=0.0`` explicitly. Because the N=5 configuration here is
EXACTLY periodic (the Laplace-projected omega_io = -2*omega_gan lock),
propagating over exactly one Ganymede synodic period vs. her stated 100 TU
window gives the SAME box/average to < 1% (checked directly -- 100 TU is not
an exact multiple of the period, so this is a confirmed close match, not a
mathematical identity). Shifting Io's phase RELATIVE to Ganymede's (a
genuinely different perturber configuration, not a time-shift along one fixed
solution) is found EMPIRICALLY below to move the "avg vs CR3BP" metric by
only ~3% across a full sweep while moving the Lagrange-Box area by up to
~40% -- an observed asymmetry (not a proven invariance), used here to justify
asserting box-area and avg-vs-CR3BP tightly at the specific default
(all-aligned) phase, while leaving dx/dy individually (phase-split-sensitive)
only sanity-bounded (same order of magnitude), not tightly asserted.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

_SAMPLE_STATE_PLANAR = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)
_SAMPLE_STATE_3D = np.array([0.6, 0.2, 0.04, 0.05, 0.35, -0.02], dtype=np.float64)
_L_KM_EUROPA = 671100.0  # this project's established Jupiter-Europa DU (tests/core/test_ccr4bp.py)


# ---------------------------------------------------------------------------
# Gate 1: mu_io -> 0 reduces the N=5 CRNBP EXACTLY to ccr4bp's N=4 (Ganymede only).
# ---------------------------------------------------------------------------


def _ccr4bp_default() -> ccr4bp.CCR4BPSystem:
    return ccr4bp.jupiter_europa_ganymede_default()


def _crnbp_matching_ccr4bp(io_mu: float = 0.0) -> crnbp.CRNBPSystem:
    """A CRNBPSystem whose Ganymede perturber matches ccr4bp's default exactly,
    with Io's mass set to ``io_mu`` (0.0 for the reduction test)."""
    d = _ccr4bp_default()
    io = crnbp.CRNBPPerturber(mu=io_mu, a=0.6285203397407242, omega=1.0, theta0=0.0)
    gan = crnbp.CRNBPPerturber(mu=d.mu_gan, a=d.a_gan, omega=d.omega_gan, theta0=d.theta_gan0)
    return crnbp.CRNBPSystem(mu=d.mu, perturbers=(io, gan))


def test_mu_io_zero_reduces_to_ccr4bp_eom_planar() -> None:
    sys5 = _crnbp_matching_ccr4bp(io_mu=0.0)
    d = _ccr4bp_default()
    rhs5 = crnbp.crnbp_eom(0.0, _SAMPLE_STATE_PLANAR, sys5)
    rhs4 = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE_PLANAR, d)
    assert np.allclose(rhs5, rhs4, rtol=0.0, atol=1e-14), (
        f"CRNBP@mu_io=0 deviates from CCR4BP: max|delta|={np.max(np.abs(rhs5 - rhs4)):.3e}"
    )


def test_mu_io_zero_reduces_to_ccr4bp_eom_3d() -> None:
    sys5 = _crnbp_matching_ccr4bp(io_mu=0.0)
    d = _ccr4bp_default()
    for t in (0.0, 1.7, 4.3):
        rhs5 = crnbp.crnbp_eom(t, _SAMPLE_STATE_3D, sys5)
        rhs4 = ccr4bp.ccr4bp_eom(t, _SAMPLE_STATE_3D, d)
        assert np.allclose(rhs5, rhs4, rtol=0.0, atol=1e-14), (
            f"CRNBP@mu_io=0 deviates from CCR4BP at t={t}: "
            f"max|delta|={np.max(np.abs(rhs5 - rhs4)):.3e}"
        )


def test_mu_io_zero_reduces_to_ccr4bp_stm() -> None:
    sys5 = _crnbp_matching_ccr4bp(io_mu=0.0)
    d = _ccr4bp_default()
    y42 = np.concatenate([_SAMPLE_STATE_3D, np.eye(6).reshape(36)])
    rhs5 = crnbp.crnbp_stm_eom(0.5, y42, sys5)
    rhs4 = ccr4bp.ccr4bp_stm_eom(0.5, y42, d)
    assert np.allclose(rhs5, rhs4, rtol=0.0, atol=1e-13), (
        f"CRNBP STM@mu_io=0 deviates from CCR4BP STM: max|delta|={np.max(np.abs(rhs5 - rhs4)):.3e}"
    )


def test_mu_io_zero_reduces_to_ccr4bp_propagator() -> None:
    sys5 = _crnbp_matching_ccr4bp(io_mu=0.0)
    d = _ccr4bp_default()
    arc5 = crnbp.propagate_crnbp(sys5, _SAMPLE_STATE_3D, 1.0)
    arc4 = ccr4bp.propagate_ccr4bp(d, _SAMPLE_STATE_3D, 1.0)
    assert np.allclose(arc5.state_f, arc4.state_f, rtol=1e-11, atol=1e-11), (
        f"CRNBP@mu_io=0 propagation deviates from CCR4BP: "
        f"max|delta|={np.max(np.abs(arc5.state_f - arc4.state_f)):.3e}"
    )


def test_both_perturbers_zero_reduces_to_cr3bp() -> None:
    """A fully degenerate check: both extra perturbers massless -> plain CR3BP."""
    d = _ccr4bp_default()
    io = crnbp.CRNBPPerturber(mu=0.0, a=0.6285, omega=1.0, theta0=0.0)
    gan = crnbp.CRNBPPerturber(mu=0.0, a=d.a_gan, omega=d.omega_gan, theta0=0.0)
    sys5 = crnbp.CRNBPSystem(mu=d.mu, perturbers=(io, gan))
    rhs5 = crnbp.crnbp_eom(0.7, _SAMPLE_STATE_3D, sys5)
    rhs3 = cr3bp.cr3bp_eom(0.7, _SAMPLE_STATE_3D, d.mu)
    assert np.allclose(rhs5, rhs3, rtol=0.0, atol=1e-14)


def test_empty_perturbers_tuple_reduces_to_cr3bp() -> None:
    d = _ccr4bp_default()
    sys5 = crnbp.CRNBPSystem(mu=d.mu, perturbers=())
    rhs5 = crnbp.crnbp_eom(1.3, _SAMPLE_STATE_3D, sys5)
    rhs3 = cr3bp.cr3bp_eom(1.3, _SAMPLE_STATE_3D, d.mu)
    assert np.allclose(rhs5, rhs3, rtol=0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# Gate 2: coupling-term sign lock (`#717` Step 0).
# ---------------------------------------------------------------------------


def test_coupling_term_sign_matches_hand_derivation() -> None:
    """A hand-computed regression lock for the corrected Eq. 11 coupling sign,
    isolated at the PER-PERTURBER (not summed-over-all-j) level.

    Important subtlety found while building this test (a genuine structural
    finding of `#717`, documented fully in ``crnbp.py``'s module docstring):
    summing the coupling contribution over ALL perturbers ``j`` cancels to
    EXACTLY ZERO regardless of which sign convention is used (both the
    corrected inner-PLUS and the erroneous digest's inner-MINUS give a
    zero total, by the same antisymmetry argument under ``j<->k`` swap --
    the overall sign in front does not affect whether the sum cancels). So
    testing the TOTAL acceleration cannot detect a sign error here; this test
    isolates perturber "a"'s OWN bracket contribution via
    :func:`crnbp._coupling_term_for_perturber`, which DOES carry real sign
    information (it is not itself symmetrized away).
    """
    mu_a, a_a, theta_a = 0.001, 1.0, 0.0  # perturber "j" on +x axis
    mu_b, a_b, theta_b = 0.0007, 1.6, math.pi / 3.0  # perturber "k"
    pert_a = crnbp.CRNBPPerturber(mu=mu_a, a=a_a, omega=0.0, theta0=theta_a)
    pert_b = crnbp.CRNBPPerturber(mu=mu_b, a=a_b, omega=0.0, theta0=theta_b)
    t = 0.0
    sys_both = crnbp.CRNBPSystem(mu=0.0, perturbers=(pert_a, pert_b))
    positions = crnbp._perturber_positions(t, sys_both)
    coupling_a_numeric = crnbp._coupling_term_for_perturber(0, positions)

    gxa, gya = a_a * math.cos(theta_a), a_a * math.sin(theta_a)
    gxb, gyb = a_b * math.cos(theta_b), a_b * math.sin(theta_b)
    rab_x, rab_y = gxa - gxb, gya - gyb
    rab3 = (rab_x * rab_x + rab_y * rab_y) ** 1.5
    # coupling_a = -mu_a * mu_b * (r_a - r_b) / |r_a - r_b|^3 (the corrected sign)
    coupling_hand_x = -mu_a * mu_b * rab_x / rab3
    coupling_hand_y = -mu_a * mu_b * rab_y / rab3
    # what the ERRONEOUS digest's inner-MINUS would have given instead
    # (coupling_a = +mu_a*mu_b*(r_a-r_b)/|r_a-r_b|^3, the opposite sign):
    coupling_wrong_x = +mu_a * mu_b * rab_x / rab3

    assert coupling_a_numeric[0] == pytest.approx(coupling_hand_x, abs=1e-15)
    assert coupling_a_numeric[1] == pytest.approx(coupling_hand_y, abs=1e-15)
    # Non-triviality + sign discrimination: not accidentally zero, and not
    # accidentally equal to the WRONG-sign value either.
    assert abs(coupling_a_numeric[0]) > 1e-8
    assert coupling_a_numeric[0] != pytest.approx(coupling_wrong_x, abs=1e-8)


def test_coupling_term_net_contribution_is_exactly_zero() -> None:
    """A genuinely new structural finding of `#717` (beyond what `#714`
    assessed): summed over ALL perturbers, the coupling term contributes
    EXACTLY ZERO to the total spacecraft acceleration -- for N=5 (two extra
    perturbers) and, per the algebraic proof in ``crnbp.py``'s module
    docstring, for any N. Verified here against the actual physical Io+
    Ganymede default system (not just the toy hand-example above): the
    full ``_perturbers_acceleration`` (WITH the coupling term, as coded)
    must equal the sum of two INDEPENDENT single-perturber accelerations
    (each with mu_io or mu_gan alone, so the coupling call is a no-op)."""
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    io, gan = sys5.perturbers
    sys_io_only = crnbp.CRNBPSystem(mu=sys5.mu, perturbers=(io,))
    sys_gan_only = crnbp.CRNBPSystem(mu=sys5.mu, perturbers=(gan,))

    for x, y, z, t in [(0.6, 0.2, 0.03, 1.7), (0.98, -0.1, 0.0, 4.2), (1.05, 0.05, 0.01, 9.9)]:
        a_full = np.array(crnbp._perturbers_acceleration(x, y, z, t, sys5))
        a_io = np.array(crnbp._perturbers_acceleration(x, y, z, t, sys_io_only))
        a_gan = np.array(crnbp._perturbers_acceleration(x, y, z, t, sys_gan_only))
        assert np.allclose(a_full, a_io + a_gan, rtol=0.0, atol=1e-16), (
            f"coupling term did not cancel to zero at (x,y,z,t)=({x},{y},{z},{t}): "
            f"full={a_full}, superposition={a_io + a_gan}"
        )


def test_coupling_term_net_zero_holds_at_n6_independent_transcription() -> None:
    """From-scratch (NOT reusing any ``crnbp.py`` code) literal transcription
    of Gilliam's Eq. 25 x-component at N=6 (THREE extra perturbers), to check
    the cancellation is a general-N property and not an N=5-specific
    coincidence or an artifact of this module's own implementation."""

    def xddot_literal(
        x: float, y: float, mu: float, extras: dict[int, tuple[float, float, float]]
    ) -> float:
        r1 = math.sqrt((x + mu) ** 2 + y * y)
        r2 = math.sqrt((x - 1 + mu) ** 2 + y * y)
        mu1, mu2 = 1.0 - mu, mu
        bodies = {1: (mu1, 0.0, 0.0), 2: (mu2, 1.0, 0.0)}
        bodies.update(extras)
        term = x - mu1 * (x + mu) / r1**3 - mu2 * (x - 1 + mu) / r2**3
        for j, (mu_j, r_j, psi_j) in extras.items():
            rj = math.sqrt((x - r_j * math.cos(psi_j)) ** 2 + (y - r_j * math.sin(psi_j)) ** 2)
            direct = (x + mu - r_j * math.cos(psi_j)) / rj**3
            inner = 0.0
            for k, (mu_k, r_k, psi_k) in bodies.items():
                if k == j:
                    continue
                rkj = math.sqrt(r_k**2 + r_j**2 - 2 * r_k * r_j * math.cos(psi_k - psi_j))
                inner += mu_k * (r_j * math.cos(psi_j) - r_k * math.cos(psi_k)) / rkj**3
            term -= mu_j * (direct + inner)
        return term

    def xddot_superposition(
        x: float, y: float, mu: float, extras: dict[int, tuple[float, float, float]]
    ) -> float:
        """Same but each extra body's inner sum uses ONLY the primaries
        (k=1,2), i.e. no cross-extra-body coupling at all."""
        r1 = math.sqrt((x + mu) ** 2 + y * y)
        r2 = math.sqrt((x - 1 + mu) ** 2 + y * y)
        mu1, mu2 = 1.0 - mu, mu
        term = x - mu1 * (x + mu) / r1**3 - mu2 * (x - 1 + mu) / r2**3
        for _j, (mu_j, r_j, psi_j) in extras.items():
            rj = math.sqrt((x - r_j * math.cos(psi_j)) ** 2 + (y - r_j * math.sin(psi_j)) ** 2)
            direct = (x + mu - r_j * math.cos(psi_j)) / rj**3
            r1j = r_j
            r2j = math.sqrt(1.0 + r_j**2 - 2.0 * r_j * math.cos(psi_j))
            inner = (
                mu1 * (r_j * math.cos(psi_j)) / r1j**3
                + mu2 * (r_j * math.cos(psi_j) - 1.0) / r2j**3
            )
            term -= mu_j * (direct + inner)
        return term

    mu = 2.528e-5
    extras = {3: (4.70e-5, 0.6285, 0.9), 4: (7.80e-5, 1.595, 2.1), 5: (5.67e-5, 2.804, 4.2)}
    x, y = 0.6, 0.2
    full = xddot_literal(x, y, mu, extras)
    superp = xddot_superposition(x, y, mu, extras)
    assert full == pytest.approx(superp, abs=1e-13), (
        f"N=6 (3 extra bodies) coupling term did not cancel: full={full}, superposition={superp}"
    )


def test_coupling_term_vanishes_with_one_perturber() -> None:
    """With only ONE nonzero-mass perturber the inner sum is vacuous (k has no
    valid range) -- the coupling contributes nothing, exactly ccr4bp's shape."""
    d = _ccr4bp_default()
    sys5 = _crnbp_matching_ccr4bp(io_mu=0.0)
    ax, ay, az = crnbp._perturbers_acceleration(0.6, 0.2, 0.0, 0.3, sys5)
    ax_g, ay_g, az_g = ccr4bp._ganymede_acceleration(0.6, 0.2, 0.0, 0.3, d)
    assert (ax, ay, az) == pytest.approx((ax_g, ay_g, az_g), abs=1e-15)


def test_coupling_term_state_independent_no_jacobian_contribution() -> None:
    """The coupling term (and indirect term) must not appear in the STM's
    Hessian block: perturbing (x, y, z) must change the perturber acceleration
    by EXACTLY the analytic direct-term-only Hessian block
    (:func:`crnbp._perturber_second_deriv_block`), confirming `#714`'s
    structural claim directly at the acceleration-function level (independent
    of, and more targeted than, the whole-propagation FD-vs-STM check)."""
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    t = 1.1
    eps = 1e-6
    base_pos = [0.6, 0.2, 0.03]

    def accel_from(x: float, y: float, z: float) -> np.ndarray:
        return np.array(crnbp._perturbers_acceleration(x, y, z, t, sys5))

    fd_jac = np.zeros((3, 3), dtype=np.float64)
    for coord in range(3):
        pos_p = list(base_pos)
        pos_p[coord] += eps
        pos_m = list(base_pos)
        pos_m[coord] -= eps
        fd_jac[:, coord] = (accel_from(*pos_p) - accel_from(*pos_m)) / (2 * eps)

    sxx, syy, szz, sxy, sxz, syz = crnbp._perturber_second_deriv_block(
        base_pos[0], base_pos[1], base_pos[2], t, sys5
    )
    analytic_jac = np.array([[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]], dtype=np.float64)
    assert np.allclose(fd_jac, analytic_jac, rtol=1e-5, atol=1e-8), (
        f"FD Jacobian of the perturber acceleration disagrees with the analytic "
        f"direct-term-only Hessian block: max|delta|="
        f"{np.max(np.abs(fd_jac - analytic_jac)):.3e} -- the coupling/indirect "
        f"terms leaked into the finite-difference derivative, contradicting "
        f"`#714`'s state-independence claim."
    )


# ---------------------------------------------------------------------------
# Gate 3: STM finite-difference consistency (full N=5 model, both perturbers).
# ---------------------------------------------------------------------------


def test_stm_finite_difference_consistency() -> None:
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    state0 = _SAMPLE_STATE_3D.copy()
    t_horizon = 0.3

    arc = crnbp.propagate_crnbp(sys5, state0, t_horizon, with_stm=True)
    assert arc.stm is not None
    stm = arc.stm

    eps = 1e-6
    fd_jac = np.zeros((6, 6), dtype=np.float64)
    for i in range(6):
        sp = state0.copy()
        sp[i] += eps
        sm = state0.copy()
        sm[i] -= eps
        arc_p = crnbp.propagate_crnbp(sys5, sp, t_horizon)
        arc_m = crnbp.propagate_crnbp(sys5, sm, t_horizon)
        fd_jac[:, i] = (arc_p.state_f - arc_m.state_f) / (2.0 * eps)

    scale = max(1.0, float(np.max(np.abs(stm))))
    rel_delta = float(np.max(np.abs(fd_jac - stm))) / scale
    assert rel_delta < 1e-6, (
        f"N=5 STM vs central-FD relative disagreement: {rel_delta:.3e} "
        f"(abs max|delta|={np.max(np.abs(fd_jac - stm)):.3e}, stm_scale={scale:.3e})"
    )


# ---------------------------------------------------------------------------
# Gate 4: sourced-constant cross-check + Laplace-lock projection.
# ---------------------------------------------------------------------------


def test_default_ganymede_matches_ccr4bp_exactly() -> None:
    """The N=5 default's Ganymede perturber must be BIT-IDENTICAL to ccr4bp's
    own default (same registry, same formula, no independent re-derivation)."""
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    d4 = ccr4bp.jupiter_europa_ganymede_default()
    _io, gan = sys5.perturbers
    assert sys5.mu == d4.mu
    assert gan.mu == d4.mu_gan
    assert gan.a == d4.a_gan
    assert gan.omega == d4.omega_gan


def test_io_constants_traced_to_registry() -> None:
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    io, _gan = sys5.perturbers
    gm_j = PRIMARIES["Jupiter"]
    gm_e = SATELLITES["Europa"].mu_km3_s2
    gm_io = SATELLITES["Io"].mu_km3_s2
    denom = gm_j + gm_e
    assert io.mu == pytest.approx(gm_io / denom, rel=0.0, abs=1e-18)
    assert io.a == pytest.approx(
        SATELLITES["Io"].sma_km / SATELLITES["Europa"].sma_km, rel=0.0, abs=1e-15
    )


def test_laplace_lock_projection_exact_and_close_to_observed() -> None:
    """omega_io = -2*omega_gan EXACTLY (the projection); the registry two-body
    rate ratio is close (~2e-4) but not exact -- confirming `#714`'s finding
    that the projection is a mild, well-justified idealisation."""
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    io, gan = sys5.perturbers
    assert io.omega == pytest.approx(-2.0 * gan.omega, rel=0.0, abs=1e-15)

    raw_omega_io = ccr4bp.two_body_synodic_rate(sys5.mu, io.mu, io.a)
    ratio = raw_omega_io / gan.omega
    assert ratio == pytest.approx(-2.0, abs=5e-4)
    # The projection changes the rate by the same small order as the raw-vs-
    # projected deviation (not by an unrelated, much larger amount).
    rel_change = abs(io.omega - raw_omega_io) / abs(raw_omega_io)
    assert rel_change < 1e-3


def test_default_phase_sits_at_physical_laplace_libration_center() -> None:
    """`#723` regression: :func:`jupiter_europa_io_ganymede_default`'s DEFAULT
    phase must be the physically-real Laplace-libration center, not the
    all-aligned antipode `#721` found shipped as the default.

    In this module's Europa-synodic phase convention (``theta_j(t) ==
    lambda_j(t) - lambda_Europa(t)`` exactly, since Europa itself is fixed at
    synodic angle 0 by construction), the classical Laplace resonance argument
    ``phi_L = lambda_Io - 3*lambda_Europa + 2*lambda_Ganymede`` becomes
    ``phi_L = theta_io + 2*theta_gan``. The real Galilean trio librates about
    ``phi_L = 180 deg`` (Sinclair 1975; Murray & Dermott; corroborated by
    Gilliam's thesis prose and by Baresi/Owen/Scheeres AAS 23-257 Table 1's
    Jupiter-Europa-frame ``(phi_Io, phi_Gan) = (pi, 0)``) -- NOT 0 deg, the
    old buggy default's value.
    """
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    io, gan = sys5.perturbers
    assert gan.theta0 == pytest.approx(0.0, abs=1e-15)
    phi_l = (io.theta0 + 2.0 * gan.theta0) % (2.0 * math.pi)
    assert phi_l == pytest.approx(math.pi, abs=1e-12), (
        f"default Laplace argument phi_L={math.degrees(phi_l):.4f} deg, expected 180 deg"
    )

    # The relation generalizes correctly for nonzero theta_gan0 -- not a
    # hardcoded special case that only happens to work at theta_gan0=0.
    for theta_gan0 in (0.4, -1.1, 2.7):
        sys_shifted = crnbp.jupiter_europa_io_ganymede_default(theta_gan0=theta_gan0)
        io_s, gan_s = sys_shifted.perturbers
        assert gan_s.theta0 == pytest.approx(theta_gan0, abs=1e-15)
        phi_l_s = (io_s.theta0 + 2.0 * gan_s.theta0) % (2.0 * math.pi)
        assert phi_l_s == pytest.approx(math.pi, abs=1e-12), (
            f"theta_gan0={theta_gan0}: phi_L={math.degrees(phi_l_s):.4f} deg, expected 180 deg"
        )

    # theta_io0 remains an explicit override for the Gilliam all-aligned
    # ("test case 1") equilibrium-point positive control.
    sys_aligned = crnbp.jupiter_europa_io_ganymede_default(theta_io0=0.0)
    io_a, gan_a = sys_aligned.perturbers
    assert io_a.theta0 == 0.0
    assert gan_a.theta0 == 0.0


# ---------------------------------------------------------------------------
# Gate 5: equilibrium-point ("dynamical substitute") solver correctness.
# ---------------------------------------------------------------------------


def test_cr3bp_collinear_points_match_gilliam_table4() -> None:
    """Base CR3BP L1/L2/L3 (independent of any perturber) must match Gilliam's
    thesis Table 4 Jupiter-Europa values (JPL Three-Body Periodic Orbit
    Catalog), confirming this project's registry-derived mu2 agrees with hers
    to high precision (2.528017724591319e-5 vs. her stated
    2.528017528540000e-5 -- ~7.8e-8 relative)."""
    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    x_l1 = crnbp.cr3bp_collinear_point(sys5.mu, "L1")
    x_l2 = crnbp.cr3bp_collinear_point(sys5.mu, "L2")
    x_l3 = crnbp.cr3bp_collinear_point(sys5.mu, "L3")
    assert x_l1 == pytest.approx(0.97976410, abs=2e-7)
    assert x_l2 == pytest.approx(1.02046139, abs=2e-7)
    assert x_l3 == pytest.approx(-1.00001053, abs=2e-7)


def test_equilibrium_solver_matches_base_cr3bp_at_zero_perturber_mass() -> None:
    """Structural sanity: with both perturbers massless the "equilibrium
    dynamical substitute" solver must sit exactly at the base CR3BP point for
    all times (no drift at all)."""
    d = _ccr4bp_default()
    io = crnbp.CRNBPPerturber(mu=0.0, a=0.6285, omega=1.0, theta0=0.0)
    gan = crnbp.CRNBPPerturber(mu=0.0, a=d.a_gan, omega=d.omega_gan, theta0=0.0)
    sys5 = crnbp.CRNBPSystem(mu=d.mu, perturbers=(io, gan))
    result = crnbp.equilibrium_dynamical_substitute(
        sys5, "L1", l_km=_L_KM_EUROPA, t_window=5.0, n_samples=20
    )
    assert result.dx_km == pytest.approx(0.0, abs=1e-6)
    assert result.dy_km == pytest.approx(0.0, abs=1e-6)
    assert result.avg_vs_cr3bp_km == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize(
    "which,gilliam_dx_km,gilliam_dy_km,gilliam_box_km2,gilliam_avg_km",
    [
        ("L1", 57.8603, 44.0709, 2549.59, 4.27893),
        ("L2", 57.8969, 44.3684, 2568.65, 3.72655),
    ],
)
def test_equilibrium_dynamical_substitute_matches_gilliam_thesis_table(
    which: str,
    gilliam_dx_km: float,
    gilliam_dy_km: float,
    gilliam_box_km2: float,
    gilliam_avg_km: float,
) -> None:
    """`#717` Step 2 positive control: reproduce Gilliam's thesis Table 5 (E1,
    ``which="L1"``) / Table 6 (E2, ``which="L2"``) Jupiter-Europa-Io-Ganymede
    dynamical-substitute Lagrange-Box numbers.

    Propagated over exactly ONE Ganymede synodic period (the system is exactly
    periodic under the Laplace-lock projection, so this equals any longer
    window's steady-state box/average -- confirmed separately by re-running at
    a 100 TU window, matching Gilliam's own stated propagation length, to
    < 1% difference from the one-period result).

    See the module docstring for the precision caveats: our mu_io/mu_gan
    differ from Gilliam's own stated values by ~3-5e-4 relative (different
    mass-rounding convention, same underlying JPL data), and her epoch phase
    angles are not stated (we use the "all perturbers aligned at t=0"
    convention, requested here via ``theta_io0=0.0`` -- `#723` corrected
    ``jupiter_europa_io_ganymede_default``'s own DEFAULT to the physically-real
    Laplace-libration phase, so this control test must now request the
    all-aligned idealisation explicitly rather than relying on the default).
    The average-vs-CR3BP distance is empirically found (see
    ``test_equilibrium_average_is_phase_choice_invariant_over_one_period``
    below) to change little (~3%) across different Io/Ganymede relative-phase
    choices, and the Lagrange-Box area agrees closely at this specific
    (all-aligned) phase -- both are asserted tightly here; dx/dy
    individually (phase-split-sensitive, varying up to ~40% across phase
    choices) are only sanity-bounded.
    """
    sys5 = crnbp.jupiter_europa_io_ganymede_default(theta_io0=0.0)
    _io, gan = sys5.perturbers
    result = crnbp.equilibrium_dynamical_substitute(
        sys5, which, l_km=_L_KM_EUROPA, t_window=gan.synodic_period, n_samples=800
    )

    assert result.box_area_km2 == pytest.approx(gilliam_box_km2, rel=0.10), (
        f"{which} Lagrange Box area {result.box_area_km2:.2f} km^2 vs. "
        f"Gilliam's {gilliam_box_km2:.2f} km^2"
    )
    assert result.avg_vs_cr3bp_km == pytest.approx(gilliam_avg_km, rel=0.15), (
        f"{which} avg-vs-CR3BP {result.avg_vs_cr3bp_km:.4f} km vs. "
        f"Gilliam's {gilliam_avg_km:.4f} km"
    )
    # dx/dy: same order of magnitude, sum-of-both within a documented loose
    # bound (the individual x/y split depends on the undocumented epoch
    # phase; the TOTAL excursion is a milder, still-informative check).
    assert 0.3 * gilliam_dx_km < result.dx_km < 3.0 * gilliam_dx_km
    assert 0.3 * gilliam_dy_km < result.dy_km < 3.0 * gilliam_dy_km


def test_equilibrium_average_is_phase_choice_invariant_over_one_period() -> None:
    """Empirical basis for the tight ``avg_vs_cr3bp_km`` tolerance above:
    shifting Io's initial phase RELATIVE to Ganymede's is a genuinely
    different perturber configuration (not a time-shift along one fixed
    solution, since the -2 Laplace lock ties a pure time-shift's phase changes
    together), yet the time-averaged substitute location over one full period
    is found to move only slightly. This is an OBSERVED property (checked
    here), not a proven theorem -- the Lagrange-Box area, checked the same
    way, is NOT similarly stable (it varies by up to ~40% across the same
    phase sweep; only asserted at the specific default phase in the
    Gilliam-comparison test above)."""
    from dataclasses import replace

    base = crnbp.jupiter_europa_io_ganymede_default()
    io, gan = base.perturbers
    avgs = []
    boxes = []
    for frac in (0.0, 0.25, 0.5, 0.75):
        io_shifted = replace(io, theta0=2.0 * math.pi * frac)
        sys5 = crnbp.CRNBPSystem(mu=base.mu, perturbers=(io_shifted, gan))
        result = crnbp.equilibrium_dynamical_substitute(
            sys5, "L1", l_km=_L_KM_EUROPA, t_window=gan.synodic_period, n_samples=400
        )
        avgs.append(result.avg_vs_cr3bp_km)
        boxes.append(result.box_area_km2)
    avgs_arr = np.array(avgs)
    boxes_arr = np.array(boxes)
    assert (np.max(avgs_arr) - np.min(avgs_arr)) / np.mean(avgs_arr) < 0.10, (
        f"avg_vs_cr3bp varied more than expected across phase shifts: {avgs}"
    )
    # Contrast: the Lagrange-Box area is NOT similarly stable (the asymmetry
    # this test's docstring and the module docstring both describe).
    assert (np.max(boxes_arr) - np.min(boxes_arr)) / np.mean(boxes_arr) > 0.20, (
        f"expected the Lagrange-Box area to vary noticeably across phase shifts "
        f"(unlike avg_vs_cr3bp); got: {boxes}"
    )
