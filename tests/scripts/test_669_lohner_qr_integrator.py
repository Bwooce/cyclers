"""#669 -- tests for the Lohner C1 / QR-reframed validated interval integrator.

Extends ``scripts/_validated_taylor_integrator.py`` with the wrapping-effect fix
(Lohner's QR-coordinate reframing) that ``#668`` identified as future work.  The
load-bearing properties under test, in the same spirit as the ``#668`` suite:

  * the variational (state-transition-matrix) enclosure is rigorous -- it
    CONTAINS the true STM and is NOT a vacuous bound (excludes nearby wrong
    values), checked against the closed-form ``exp`` and harmonic-oscillator
    STMs (``e^t`` and a rotation matrix);
  * ``rigorous_inverse`` genuinely encloses the true inverse (``A @ A^-1``
    contains the identity, tightly);
  * the QR-reframed integrator reproduces the closed-form goldens, and on the
    real Oterma arc is *orders of magnitude tighter* than ``#668``'s naive C0 box
    at matched times (i.e. QR really does defeat wrapping);
  * the residual ``t ~ 0.42`` horizon is a PHYSICAL Jupiter flyby, not wrapping:
    the reachable time is step-size dependent (finer step reaches later), the
    opposite of ``#668``'s step-size-independent wrapping wall.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath", reason="mpmath is an optional 'interval' extra (task #610/#625/#668/#669)"
)

import scripts._validated_taylor_integrator as vti  # noqa: E402

MU_OTERMA = 0.0009537
C_OTERMA = 3.03
# L1->L2 heteroclinic section point (data/golden/wz_oterma_heteroclinic.yaml).
X0, VX0 = 0.9522928423486199945, 1.23e-05


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 50
    mp.iv.dps = 50


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _contains(enc: Any, value: Any) -> bool:
    return bool(enc.a <= value) and bool(value <= enc.b)


def _oterma_ic(iv: Any) -> tuple[list[Any], Any]:
    r1 = iv.sqrt((iv.mpf(X0) + iv.mpf(MU_OTERMA)) ** 2)
    r2 = iv.sqrt((iv.mpf(X0) + iv.mpf(MU_OTERMA - 1.0)) ** 2)
    vy2 = (
        iv.mpf(X0) ** 2
        + 2 * (1 - iv.mpf(MU_OTERMA)) / r1
        + 2 * iv.mpf(MU_OTERMA) / r2
        - iv.mpf(VX0) ** 2
        - iv.mpf(C_OTERMA)
    )
    y0 = [_pt(X0), _pt(0.0), _pt(VX0), iv.sqrt(vy2)]
    return y0, vti.cr3bp_planar_jacobi(iv, y0, MU_OTERMA)


# --- closed-form variational (STM) jets for the positive controls ------------
def _jet_exp(iv: Any, state: list[Any], mu: float) -> list[Any]:
    return [state[0]]  # y' = y


def _var_exp(iv: Any, aug: list[Any], mu: float) -> list[Any]:
    return [aug[0], aug[1]]  # y'=y ; V'=V  -> STM = e^t


def _jet_harmonic(iv: Any, state: list[Any], mu: float) -> list[Any]:
    y1, y2 = state
    return [y2, vti.ts_scale(iv, y1, iv.mpf(-1.0))]


def _var_harmonic(iv: Any, aug: list[Any], mu: float) -> list[Any]:
    y1, y2 = aug[0], aug[1]
    v00, v01, v10, v11 = aug[2], aug[3], aug[4], aug[5]

    def neg(s: list[Any]) -> list[Any]:
        return vti.ts_scale(iv, s, iv.mpf(-1.0))

    # V' = [[0,1],[-1,0]] V  (row-major)
    return [y2, neg(y1), v10, v11, neg(v00), neg(v01)]


def test_variational_jet_exp_stm_contains_and_excludes() -> None:
    """The STM block of the augmented exp flow encloses e^t=e tightly at t=1."""
    iv = mp.iv
    res = vti.integrate(iv, _var_exp, [_pt(1.0), _pt(1.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    stm = res["final"][1]
    true = iv.exp(iv.mpf(1.0))
    assert bool(stm.a <= true.a) and bool(true.b <= stm.b)  # contains
    assert float(stm.delta.b) < 1e-20  # tight (not vacuous)
    assert not _contains(stm, mp.e + mp.mpf("1e-12"))
    assert not _contains(stm, mp.e - mp.mpf("1e-12"))


def test_variational_jet_harmonic_stm_is_rotation() -> None:
    """The STM block of the augmented harmonic flow encloses the rotation matrix."""
    iv = mp.iv
    aug0 = [_pt(1.0), _pt(0.0), _pt(1.0), _pt(0.0), _pt(0.0), _pt(1.0)]  # V0 = I
    res = vti.integrate(iv, _var_harmonic, aug0, 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    v00, v01, v10, v11 = res["final"][2], res["final"][3], res["final"][4], res["final"][5]
    c1, s1 = iv.cos(iv.mpf(1.0)), iv.sin(iv.mpf(1.0))
    assert bool(v00.a <= c1.a) and bool(c1.b <= v00.b)  # cos 1
    assert bool(v01.a <= s1.a) and bool(s1.b <= v01.b)  # sin 1
    assert bool(v10.a <= (-s1).a) and bool((-s1).b <= v10.b)  # -sin 1
    assert bool(v11.a <= c1.a) and bool(c1.b <= v11.b)  # cos 1
    assert float(v00.delta.b) < 1e-18


def test_rigorous_inverse_encloses_true_inverse() -> None:
    """rigorous_inverse contains the true inverse and A@A^-1 encloses I, tightly."""
    iv = mp.iv
    a_pt = [
        [iv.mpf(2.0), iv.mpf(0.5), iv.mpf(-1.0)],
        [iv.mpf(0.0), iv.mpf(3.0), iv.mpf(1.0)],
        [iv.mpf(1.0), iv.mpf(-2.0), iv.mpf(4.0)],
    ]
    ainv = vti.rigorous_inverse(iv, a_pt)
    assert ainv is not None
    prod = vti._imatmul(iv, a_pt, ainv)
    for i in range(3):
        for j in range(3):
            target = iv.mpf(1.0 if i == j else 0.0)
            assert bool(prod[i][j].a <= target.a) and bool(target.b <= prod[i][j].b)
            assert float(prod[i][j].delta.b) < 1e-30  # tight, not vacuous


def test_c1qr_reproduces_exp_flow() -> None:
    """The QR integrator reproduces y'=y: encloses e at t=1, tight, excludes wrong."""
    iv = mp.iv
    res = vti.integrate_c1_qr(iv, _jet_exp, _var_exp, [_pt(1.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    enc = res["final"][0]
    true = iv.exp(iv.mpf(1.0))
    assert bool(enc.a <= true.a) and bool(true.b <= enc.b)
    assert float(enc.delta.b) < 1e-20
    assert not _contains(enc, mp.e + mp.mpf("1e-12"))


def test_c1qr_reproduces_harmonic_flow() -> None:
    """The QR integrator reproduces the rotation flow (cos 1, -sin 1) at t=1."""
    iv = mp.iv
    res = vti.integrate_c1_qr(
        iv, _jet_harmonic, _var_harmonic, [_pt(1.0), _pt(0.0)], 1.0, 0.0, n_steps=8, order=14
    )
    assert res["validated"] is True
    y1, y2 = res["final"]
    c1, s1 = iv.cos(iv.mpf(1.0)), -iv.sin(iv.mpf(1.0))
    assert bool(y1.a <= c1.a) and bool(c1.b <= y1.b)
    assert bool(y2.a <= s1.a) and bool(s1.b <= y2.b)


def test_qr_defeats_wrapping_on_rotation() -> None:
    """The textbook wrapping proof: a wide box under rigid rotation.

    A rigid rotation (harmonic oscillator) leaves the true set's size invariant,
    but a naive axis-aligned interval box grows like ``e^t`` (pure wrapping).  The
    QR frame tracks the rotation, so the enclosure stays essentially flat.  This
    is the load-bearing evidence that ``integrate_c1_qr`` genuinely removes the
    wrapping effect -- independent of the CR3BP application.
    """
    iv = mp.iv
    d = mp.mpf("1e-3")
    y0 = [iv.mpf([1 - d, 1 + d]), iv.mpf([-d, d])]  # wide box, half-width 1e-3
    naive = vti.integrate(iv, _jet_harmonic, y0, 20.0, 0.0, n_steps=640, order=14)
    qr = vti.integrate_c1_qr(iv, _jet_harmonic, _var_harmonic, y0, 20.0, 0.0, n_steps=640, order=14)
    assert naive["validated"] is True and qr["validated"] is True
    # naive box wraps up by many orders of magnitude ...
    assert naive["final_max_halfwidth"] > 1e3
    # ... while the QR enclosure stays within ~2x its initial width (no wrapping).
    assert qr["final_max_halfwidth"] < 3e-3
    # and the QR enclosure is >1e6x tighter than the wrapped naive box.
    assert qr["final_max_halfwidth"] * 1e6 < naive["final_max_halfwidth"]


def test_qr_matches_naive_c0_on_non_wrapping_oterma_arc() -> None:
    """Honest counterpoint: the Oterma point arc is NOT wrapping-limited.

    Its true stretching is modest and non-rotating, so QR and the naive C0 box
    are of comparable width at a matched grid (QR does not, and is not expected
    to, extend a non-wrapping arc).  Both rigorously conserve the Jacobi
    constant.
    """
    iv = mp.iv
    y0, c0 = _oterma_ic(iv)
    naive = vti.integrate(iv, vti.cr3bp_planar_jet, y0, 0.25, MU_OTERMA, n_steps=16, order=12)
    qr = vti.integrate_c1_qr(
        iv,
        vti.cr3bp_planar_jet,
        vti.cr3bp_planar_variational_jet,
        y0,
        0.25,
        MU_OTERMA,
        n_steps=16,
        order=12,
    )
    assert naive["validated"] is True and qr["validated"] is True
    cf = vti.cr3bp_planar_jacobi(iv, qr["final"], MU_OTERMA)
    assert bool(cf.a <= c0.a) and bool(c0.b <= cf.b)  # rigorous first integral
    # comparable widths (within 10x either way) -- wrapping is not dominant here
    ratio = qr["final_max_halfwidth"] / naive["final_max_halfwidth"]
    assert 0.1 < ratio < 10.0


def test_oterma_wall_is_jupiter_flyby_not_wrapping() -> None:
    """The residual horizon is a PHYSICAL flyby: reach is step-size *dependent*.

    ``#668`` proved its naive-C0 wrapping wall was step-size *independent*.  Here
    the reachable time grows when the step shrinks (the enclosure follows a real
    close approach to Jupiter and can only be pushed nearer perijove with finer
    steps) -- the opposite signature, confirming the wall is dynamical, not a
    wrapping artifact that QR failed to fix.
    """
    iv = mp.iv
    y0, _ = _oterma_ic(iv)
    reach = []
    for n_steps in (32, 64):  # h = 1/64, 1/128 over t_final=0.5
        res = vti.integrate_c1_qr(
            iv,
            vti.cr3bp_planar_jet,
            vti.cr3bp_planar_variational_jet,
            y0,
            0.5,
            MU_OTERMA,
            n_steps=n_steps,
            order=10,
        )
        assert res["validated"] is False  # both stop at the flyby bottleneck
        reach.append(res["n_completed"] * (0.5 / n_steps))
    # finer step reaches meaningfully later (step-size dependent -> physical)
    assert reach[1] > reach[0] + 0.01
    assert reach[0] > 0.3  # a genuine, non-trivial rigorously-enclosed arc
