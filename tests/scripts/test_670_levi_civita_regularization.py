"""#670 -- tests for Levi-Civita close-approach regularization of the CR3BP.

Stage 4 of the Wilczak-Zgliczynski proof-machinery build (``#668`` Stages 1-2,
``#669`` Stage 3).  ``#669`` proved the rigorous Oterma L1->L2 arc enclosure's
real horizon is a PHYSICAL close flyby of the secondary (Jupiter perijove
``r2 ~ 0.004`` at ``t ~ 0.466``): a genuine ``1/r2`` singularity in the equations
of motion that no coordinate REFRAMING (QR/Lohner-C1) can remove.  This adds the
classical fix -- **Levi-Civita regularization** (``z = w^2`` complex-square map +
fictitious-time ``dt = r2 dtau``), which turns the near-collision into a REGULAR,
analytic arc (:func:`_validated_taylor_integrator.make_cr3bp_lc_secondary_jet`).

Load-bearing properties under test, in the same containment discipline as the
``#668``/``#669`` suites (contains-true-value AND excludes-nearby-wrong-value,
non-vacuous):

  * the coordinate map ``physical <-> regularized`` is a rigorous identity on
    round-trip;
  * **closed-form anchor**: a rectilinear (radial) fall to a point mass -- a
    genuine two-body *collision* orbit -- is a harmonic oscillator in
    Levi-Civita coordinates, so the rigorous enclosure traverses ``r2 = 0``
    SMOOTHLY and matches the closed-form ``sqrt(a) cos(omega tau)`` and Kepler
    free-fall time;
  * the regularized CR3BP enclosure carries the real Oterma arc THROUGH the
    Jupiter perijove that stopped ``#668``/``#669`` (``t > 0.466``), still
    rigorously conserving the Jacobi constant and non-vacuous;
  * the Stage 2/3 positive controls (closed-form exp / harmonic goldens) are
    unbroken by the addition.

Note on dyadic steps: the closed-form-comparison tests use ``t_final / n_steps``
a power of two so the integrator's accumulated time is *exact* -- otherwise the
(rigorous) enclosure encloses the solution at ``n_steps * float(h)``, which
differs from the round ``t_final`` by ~1e-16 and (correctly) would not contain a
closed form evaluated at the round value.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath", reason="mpmath is an optional 'interval' extra (task #610/#625/#668/#669/#670)"
)

import scripts._validated_taylor_integrator as vti  # noqa: E402

MU_OTERMA = 0.0009537
C_OTERMA = 3.03
H_OTERMA = -C_OTERMA / 2  # fixed Jacobi energy h = -C/2
# L1->L2 heteroclinic section point 0 (data/golden/wz_oterma_heteroclinic.yaml).
X0, VX0 = 0.9522928423486199945, 1.23e-05
# Jupiter perijove of this arc, from #669 (float DOP853 diagnostic): the wall.
PERIJOVE_T = 0.466


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _contains(enc: Any, value: Any) -> bool:
    return bool(enc.a <= value) and bool(value <= enc.b)


def _oterma_ic(iv: Any) -> tuple[list[Any], Any]:
    """Regularized IC + Jacobi enclosure for the Oterma L1->L2 section point 0."""
    r1 = iv.sqrt((iv.mpf(X0) + iv.mpf(MU_OTERMA)) ** 2)
    r2 = iv.sqrt((iv.mpf(X0) + iv.mpf(MU_OTERMA - 1.0)) ** 2)
    vy2 = (
        iv.mpf(X0) ** 2
        + 2 * (1 - iv.mpf(MU_OTERMA)) / r1
        + 2 * iv.mpf(MU_OTERMA) / r2
        - iv.mpf(VX0) ** 2
        - iv.mpf(C_OTERMA)
    )
    vy0 = iv.sqrt(vy2)
    y0 = vti.lc_secondary_from_physical(iv, _pt(X0), _pt(0.0), _pt(VX0), vy0, MU_OTERMA)
    phys0 = vti.lc_secondary_to_physical(iv, y0, MU_OTERMA)
    return y0, vti.cr3bp_planar_jacobi(iv, phys0, MU_OTERMA)


# --- closed-form controls (unchanged from the #668/#669 machinery) ------------
def _jet_exp(iv: Any, state: list[Any], mu: float) -> list[Any]:
    return [state[0]]  # y' = y


def _jet_harmonic(iv: Any, state: list[Any], mu: float) -> list[Any]:
    y1, y2 = state
    return [y2, vti.ts_scale(iv, y1, iv.mpf(-1.0))]


# --- pure-Kepler-in-Levi-Civita jet: the closed-form collision anchor ---------
def _kepler_lc_jet_factory(h: float) -> Any:
    """Rectilinear/planar two-body Kepler around a point mass, regularized.

    No rotation, no primary: ``Gamma = (Pu^2+Pv^2)/8 - mu - r2 h`` gives
    ``u' = Pu/4``, ``Pu' = 2 h u`` (and same for v), i.e. ``u'' = (h/2) u`` -- a
    harmonic oscillator with ``omega = sqrt(-h/2)``.  ``T' = r2`` carries physical
    time.  A collision ``r2 = u^2 + v^2 -> 0`` is just ``u, v`` passing through 0.
    """

    def jet(iv: Any, state: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = state
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        return [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]

    return jet


def test_lc_coordinate_roundtrip_is_rigorous_identity() -> None:
    """physical -> regularized -> physical rigorously encloses the original, tight."""
    iv = mp.iv
    y0, _ = _oterma_ic(iv)
    back = vti.lc_secondary_to_physical(iv, y0, MU_OTERMA)
    r1 = iv.sqrt((iv.mpf(X0) + iv.mpf(MU_OTERMA)) ** 2)
    r2 = iv.sqrt((iv.mpf(X0) + iv.mpf(MU_OTERMA - 1.0)) ** 2)
    vy0 = iv.sqrt(
        iv.mpf(X0) ** 2
        + 2 * (1 - iv.mpf(MU_OTERMA)) / r1
        + 2 * iv.mpf(MU_OTERMA) / r2
        - iv.mpf(VX0) ** 2
        - iv.mpf(C_OTERMA)
    )
    for enc, true in zip(back, [iv.mpf(X0), iv.mpf(0.0), iv.mpf(VX0), vy0], strict=True):
        assert _contains(enc, true.a) and _contains(enc, true.b)
        assert float(enc.delta.b) < 1e-25  # tight, not vacuous


def test_radial_fall_traverses_collision_matches_closed_form() -> None:
    """Closed-form anchor: a two-body collision orbit is smooth in Levi-Civita.

    Radial fall from rest at distance ``a`` toward a point mass ``mu`` is a genuine
    ``r2 -> 0`` collision.  In L-C it becomes ``u(tau) = sqrt(a) cos(omega tau)``
    (a harmonic oscillator), so the *rigorous* enclosure sails straight through
    the collision at ``omega tau = pi/2`` and matches the closed form -- proving
    the regularization removes the singularity, not merely hides it.  Dyadic step
    (4.0 / 32 = 0.125) so the accumulated time is exact.
    """
    iv = mp.iv
    a = 1.0
    mu = 0.5
    h = -mu / a  # energy of rest-at-distance-a
    omega = mp.sqrt(mp.mpf(-h) / 2)
    jet = _kepler_lc_jet_factory(h)
    # IC: rest at r2 = a on the real axis -> u = sqrt(a), v = Pu = Pv = 0.
    y0 = [_pt(mp.sqrt(mp.mpf(a))), _pt(0.0), _pt(0.0), _pt(0.0), _pt(0.0)]
    tau_f = 4.0  # collision is at tau = pi ~ 3.14, INSIDE [0, 4]
    res = vti.integrate(iv, jet, y0, tau_f, mu, n_steps=32, order=16)
    assert res["validated"] is True  # never trips a singularity guard
    u_enc, t_enc = res["final"][0], res["final"][4]
    # (1) matches the closed-form harmonic solution, tightly, excludes wrong.
    u_cf = mp.sqrt(mp.mpf(a)) * mp.cos(omega * mp.mpf(tau_f))
    assert _contains(u_enc, u_cf)
    assert float(u_enc.delta.b) < 1e-25
    assert not _contains(u_enc, u_cf + mp.mpf("1e-18"))
    assert not _contains(u_enc, u_cf - mp.mpf("1e-18"))
    # (2) elapsed physical time = integral r2 dtau = a (tau/2 + sin(2 omega tau)/(4 omega));
    #     with a=1, omega=1/2 this is tau/2 + sin(tau)/2 -- the Kepler free-fall relation
    #     (collision free-fall time = pi/2 = (pi/2) sqrt(a^3/(2 mu))).
    t_cf = mp.mpf(a) * (mp.mpf(tau_f) / 2 + mp.sin(2 * omega * mp.mpf(tau_f)) / (4 * omega))
    assert _contains(t_enc, t_cf)
    assert float(t_enc.delta.b) < 1e-20


def test_regularized_arc_carries_through_jupiter_perijove() -> None:
    """THE Stage-4 result: the rigorous enclosure crosses the #668/#669 wall.

    ``#668``/``#669``'s naive-C0 and QR enclosures both stopped at the physical
    Jupiter perijove (``t ~ 0.466``).  In regularized coordinates the naive C0
    integrator carries the real Oterma arc rigorously THROUGH it: here the enclosed
    physical time's LOWER bound already exceeds the perijove time, the run stays
    validated (no collision-guard trip -- the field is analytic through the close
    approach), the Jacobi constant stays rigorously enclosed around ``C0``, and the
    enclosure is non-vacuous.  Dyadic step (19.0 / 152 = 0.125).
    """
    iv = mp.iv
    y0, c0 = _oterma_ic(iv)
    jet = vti.make_cr3bp_lc_secondary_jet(H_OTERMA)
    res = vti.integrate(iv, jet, y0, 19.0, MU_OTERMA, n_steps=152, order=10)
    assert res["validated"] is True
    box = res["final"]
    t_phys = box[4]
    assert bool(t_phys.a > PERIJOVE_T)  # rigorously PAST the perijove wall
    phys = vti.lc_secondary_to_physical(iv, box, MU_OTERMA)
    cf = vti.cr3bp_planar_jacobi(iv, phys, MU_OTERMA)
    assert bool(cf.a <= c0.a) and bool(c0.b <= cf.b)  # first integral conserved
    hw = res["final_max_halfwidth"]
    assert hw is not None and 0.0 < hw < 1.0  # non-vacuous, still tight


def test_stage2_3_positive_controls_intact() -> None:
    """Adding regularization did not break the Stage-2/3 machinery.

    Re-run the closed-form exp (``e^t``) and harmonic (rotation) goldens through
    the unchanged ``integrate`` -- contains-true-value AND excludes-nearby-wrong,
    non-vacuous, exactly as ``#668``/``#669`` established.
    """
    iv = mp.iv
    # exp: y' = y encloses e at t = 1
    res = vti.integrate(iv, _jet_exp, [_pt(1.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    enc = res["final"][0]
    true = iv.exp(iv.mpf(1.0))
    assert bool(enc.a <= true.a) and bool(true.b <= enc.b)
    assert float(enc.delta.b) < 1e-20
    assert not _contains(enc, mp.e + mp.mpf("1e-12"))
    # harmonic: encloses (cos 1, -sin 1) at t = 1
    res2 = vti.integrate(iv, _jet_harmonic, [_pt(1.0), _pt(0.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res2["validated"] is True
    y1, y2 = res2["final"]
    c1, s1 = iv.cos(iv.mpf(1.0)), -iv.sin(iv.mpf(1.0))
    assert bool(y1.a <= c1.a) and bool(c1.b <= y1.b)
    assert bool(y2.a <= s1.a) and bool(s1.b <= y2.b)
    assert float(y1.delta.b) < 1e-18
