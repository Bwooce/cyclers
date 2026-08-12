"""#671 -- tests for the REGULARIZED variational jet + QR/Lohner-C1 composition.

Stage 5 of the Wilczak-Zgliczynski proof-machinery build (``#668`` Stages 1-2,
``#669`` Stage 3 QR reframing, ``#670`` Stage 4 Levi-Civita regularization).
``#670`` regularized the Jupiter close approach that stopped Stages 2/3, carrying
the rigorous Oterma arc THROUGH the perijove -- but then hit a NEW wall (``t~0.85``,
``tau~40``) characterized as classic exponential (``e^t``) WRAPPING, exactly the
failure mode ``#669``'s QR reframing defeats.  ``#670`` integrated the regularized
field with the plain (non-variational) jet, so QR could not be applied in the
regularized frame.  This stage supplies the missing piece --
:func:`_validated_taylor_integrator.make_cr3bp_lc_secondary_variational_jet`, the
state-transition-matrix flow ``V' = Dg(w) V`` of the Levi-Civita field -- so
:func:`integrate_c1_qr` can reframe the enclosure in the *regularized*
coordinates' own local stretching directions.  Both fixes then compose.

Load-bearing properties under test, in the ``#668``/``#669``/``#670`` containment
discipline (contains-true-value AND excludes-nearby-wrong-value, non-vacuous):

  * **the analytic regularized Jacobian ``Dg`` matches finite differences** of the
    exact regularized field :func:`make_cr3bp_lc_secondary_jet` -- the derivation
    gate (incl. at near-collision ``r2 -> 0`` points);
  * **the closed-form STM anchor**: the pure-Kepler-in-Levi-Civita field is a
    harmonic oscillator, so its STM is the matrix exponential ``cos/sin`` rotation
    -- the augmented flow encloses it, tightly, excluding wrong values;
  * **radial fall through the collision WITH QR active**: the QR-reframed enclosure
    still traverses ``r2 = 0`` smoothly and matches the closed form (both fixes
    compose on a case with a known answer);
  * the Stage 2/3 exp/harmonic goldens are intact through the QR path;
  * the Stage 4 plain-``integrate`` Levi-Civita anchor is unbroken;
  * the QR-regularized real Oterma arc crosses the Jupiter perijove, rigorously
    conserves the Jacobi constant, and is *dramatically tighter* than the naive C0
    box at a matched grid (QR really defeats the regularized-frame wrapping).
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath", reason="mpmath is an optional 'interval' extra (task #610/#625/#668/#669/#670/#671)"
)

import scripts._validated_taylor_integrator as vti  # noqa: E402

MU_OTERMA = 0.0009537
C_OTERMA = 3.03
H_OTERMA = -C_OTERMA / 2  # fixed Jacobi energy h = -C/2
X0, VX0 = 0.9522928423486199945, 1.23e-05
PERIJOVE_T = 0.466  # #669 float diagnostic: the physical wall Stages 2/3 hit


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


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
    vy0 = iv.sqrt(vy2)
    y0 = vti.lc_secondary_from_physical(iv, _pt(X0), _pt(0.0), _pt(VX0), vy0, MU_OTERMA)
    phys0 = vti.lc_secondary_to_physical(iv, y0, MU_OTERMA)
    return y0, vti.cr3bp_planar_jacobi(iv, phys0, MU_OTERMA)


# --- closed-form controls (unchanged from the #668/#669 machinery) ------------
def _jet_exp(iv: Any, state: list[Any], mu: float) -> list[Any]:
    return [state[0]]


def _var_exp(iv: Any, aug: list[Any], mu: float) -> list[Any]:
    return [aug[0], aug[1]]  # y'=y ; V'=V -> STM = e^t


def _jet_harmonic(iv: Any, state: list[Any], mu: float) -> list[Any]:
    y1, y2 = state
    return [y2, vti.ts_scale(iv, y1, iv.mpf(-1.0))]


def _var_harmonic(iv: Any, aug: list[Any], mu: float) -> list[Any]:
    y1, y2 = aug[0], aug[1]
    v00, v01, v10, v11 = aug[2], aug[3], aug[4], aug[5]

    def neg(s: list[Any]) -> list[Any]:
        return vti.ts_scale(iv, s, iv.mpf(-1.0))

    return [y2, neg(y1), v10, v11, neg(v00), neg(v01)]


# --- pure-Kepler-in-Levi-Civita: closed-form STM = harmonic rotation ----------
def _kepler_lc_jet(h: float) -> Any:
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


def _kepler_lc_var_jet(h: float) -> Any:
    """Augmented (5-state + 5x5 STM) pure-Kepler-LC jet.

    Df is constant except the T row ``[2u, 2v, 0, 0, 0]``:  u'=Pu/4, v'=Pv/4,
    Pu'=2h u, Pv'=2h v, T'=u^2+v^2.  The [u,Pu] (and [v,Pv]) 2x2 blocks are
    ``B = [[0, 1/4], [2h, 0]]``, whose exponential is the harmonic rotation.
    """

    def jet(iv: Any, aug: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = aug[0], aug[1], aug[2], aug[3], aug[4]
        vmat = aug[5:30]
        order = len(u) - 1
        zero = vti.ts_const(iv, iv.mpf(0), order)
        q = vti.ts_const(iv, iv.mpf(0.25), order)
        th = vti.ts_const(iv, iv.mpf(2.0 * h), order)
        two_u = vti.ts_scale(iv, u, iv.mpf(2.0))
        two_v = vti.ts_scale(iv, v, iv.mpf(2.0))
        df = [
            [zero, zero, q, zero, zero],
            [zero, zero, zero, q, zero],
            [th, zero, zero, zero, zero],
            [zero, th, zero, zero, zero],
            [two_u, two_v, zero, zero, zero],
        ]
        dV: list[Any] = [zero for _ in range(25)]
        for c in range(5):
            col = [vmat[r * 5 + c] for r in range(5)]
            for rr in range(5):
                acc = None
                for kk in range(5):
                    if df[rr][kk] is zero:
                        continue
                    term = vti.ts_mul(iv, df[rr][kk], col[kk])
                    acc = term if acc is None else vti.ts_add(iv, acc, term)
                dV[rr * 5 + c] = acc if acc is not None else zero
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        state = [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]
        return [*state, *dV]

    return jet


# ---------------------------------------------------------------------------- #
# 1. The derivation gate: analytic Dg vs finite differences of the field.       #
# ---------------------------------------------------------------------------- #
def _dg_analytic(iv: Any, w: list[float]) -> list[list[float]]:
    """Dg(w): order-0 STM-derivative block of the variational jet seeded at V=I."""
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_OTERMA)
    box = [iv.mpf([x, x]) for x in w]
    aug = list(box) + [iv.mpf(1.0 if r == c else 0.0) for r in range(5) for c in range(5)]
    aug0 = [vti.ts_const(iv, aug[i], 0) for i in range(30)]
    d = vjet(iv, aug0, MU_OTERMA)
    stm = d[5:30]
    return [[float(stm[r * 5 + c][0].mid) for c in range(5)] for r in range(5)]


def _g_at(iv: Any, w: list[float]) -> list[float]:
    sjet = vti.make_cr3bp_lc_secondary_jet(H_OTERMA)
    order0 = [vti.ts_const(iv, iv.mpf([x, x]), 0) for x in w]
    d = sjet(iv, order0, MU_OTERMA)
    return [float(di[0].mid) for di in d]


def test_regularized_variational_jacobian_matches_finite_difference() -> None:
    """The analytic Dg equals central finite differences of the regularized field.

    This is the correctness gate for the Stage-5 derivation: the hand-derived
    5x5 Jacobian of :func:`make_cr3bp_lc_secondary_jet` must agree with a central
    finite-difference Jacobian of that exact field at a spread of points -- crucially
    including near-collision ``r2 -> 0`` states, where the physical CR3BP would be
    singular but the regularized field is analytic.
    """
    iv = mp.iv
    eps = 1e-7
    points = [
        [0.7, 0.3, 0.2, -0.4, 0.0],
        [1.2, -0.5, 0.1, 0.6, 3.0],
        [0.05, 0.02, -0.3, 0.15, 1.0],  # near collision
        [0.001, -0.002, 0.05, 0.03, 0.5],  # very near collision
        [-0.9, 0.4, 0.5, -0.2, 0.0],
    ]
    for w in points:
        a = _dg_analytic(iv, w)
        for j in range(5):
            wp = list(w)
            wp[j] += eps
            wm = list(w)
            wm[j] -= eps
            gp, gm = _g_at(iv, wp), _g_at(iv, wm)
            for i in range(5):
                fd = (gp[i] - gm[i]) / (2 * eps)
                scale = max(1.0, abs(a[i][j]), abs(fd))
                assert abs(a[i][j] - fd) / scale < 1e-5, (w, i, j, a[i][j], fd)


# ---------------------------------------------------------------------------- #
# 2. Closed-form STM anchor: Kepler-LC harmonic oscillator matrix exponential.   #
# ---------------------------------------------------------------------------- #
def test_kepler_lc_stm_matches_closed_form_matrix_exponential() -> None:
    """The augmented Kepler-LC flow encloses the harmonic-rotation STM, tightly.

    ``B = [[0, 1/4], [2h, 0]]`` has ``B^2 = (h/2) I = -omega^2 I`` (``omega =
    sqrt(-h/2)``), so ``exp(B tau) = cos(omega tau) I + (sin(omega tau)/omega) B``:
    STM_uu = cos, STM_u,Pu = sin/(4 omega), STM_Pu,u = 2h sin/omega, STM_Pu,Pu =
    cos.  The rigorous STM block must contain these and exclude nearby wrong values.
    """
    iv = mp.iv
    h = -0.5
    omega = mp.sqrt(mp.mpf(-h) / 2)
    tau = 1.0
    aug0 = [_pt(0.8), _pt(0.3), _pt(0.1), _pt(-0.2), _pt(0.0)]
    aug0 += [_pt(1.0 if r == c else 0.0) for r in range(5) for c in range(5)]
    res = vti.integrate(iv, _kepler_lc_var_jet(h), aug0, tau, 0.0, n_steps=8, order=16)
    assert res["validated"] is True
    stm = res["final"][5:30]

    def s(r: int, c: int) -> Any:
        return stm[r * 5 + c]

    cos_t = mp.cos(omega * mp.mpf(tau))
    sin_t = mp.sin(omega * mp.mpf(tau))
    exp = {
        (0, 0): cos_t,
        (0, 2): sin_t / (4 * omega),
        (2, 0): 2 * mp.mpf(h) * sin_t / omega,
        (2, 2): cos_t,
        (1, 1): cos_t,
        (1, 3): sin_t / (4 * omega),
        (3, 1): 2 * mp.mpf(h) * sin_t / omega,
        (3, 3): cos_t,
    }
    for (r, c), val in exp.items():
        assert _contains(s(r, c), val), (r, c, val)
        assert float(s(r, c).delta.b) < 1e-20  # tight, not vacuous
        assert not _contains(s(r, c), val + mp.mpf("1e-15"))
        assert not _contains(s(r, c), val - mp.mpf("1e-15"))
    # cross-block coupling is exactly zero (u,Pu do not couple to v,Pv)
    for r, c in [(0, 1), (0, 3), (2, 1), (2, 3)]:
        assert _contains(s(r, c), mp.mpf(0))
        assert float(s(r, c).delta.b) < 1e-20


# ---------------------------------------------------------------------------- #
# 3. Radial fall through the collision, WITH QR active (both fixes compose).     #
# ---------------------------------------------------------------------------- #
def test_radial_fall_qr_traverses_collision_matches_closed_form() -> None:
    """The QR-reframed regularized enclosure sails through ``r2 = 0``, matches cf.

    Same closed-form collision anchor as ``#670`` (radial fall to a point mass ->
    harmonic oscillator in Levi-Civita), now run through :func:`integrate_c1_qr`
    with the Stage-5 regularized variational jet supplying the STM.  The QR machinery
    and the regularization compose: the enclosure still traverses the collision at
    ``omega tau = pi/2`` and encloses ``sqrt(a) cos(omega tau)``, tightly, excluding
    wrong values.  Dyadic step (4.0 / 32).
    """
    iv = mp.iv
    a, mu = 1.0, 0.5
    h = -mu / a
    omega = mp.sqrt(mp.mpf(-h) / 2)
    y0 = [_pt(mp.sqrt(mp.mpf(a))), _pt(0.0), _pt(0.0), _pt(0.0), _pt(0.0)]
    tau_f = 4.0
    res = vti.integrate_c1_qr(
        iv, _kepler_lc_jet(h), _kepler_lc_var_jet(h), y0, tau_f, mu, n_steps=32, order=16
    )
    assert res["validated"] is True  # never trips a singularity guard, through r2=0
    u_enc = res["final"][0]
    u_cf = mp.sqrt(mp.mpf(a)) * mp.cos(omega * mp.mpf(tau_f))
    assert _contains(u_enc, u_cf)
    assert float(u_enc.delta.b) < 1e-20  # QR keeps it tight through the collision
    assert not _contains(u_enc, u_cf + mp.mpf("1e-15"))
    assert not _contains(u_enc, u_cf - mp.mpf("1e-15"))


# ---------------------------------------------------------------------------- #
# 4. Prior positive controls intact under the addition.                         #
# ---------------------------------------------------------------------------- #
def test_qr_positive_controls_intact() -> None:
    """Stage 2/3 exp/harmonic goldens still hold through the QR integrator."""
    iv = mp.iv
    res = vti.integrate_c1_qr(iv, _jet_exp, _var_exp, [_pt(1.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    enc = res["final"][0]
    true = iv.exp(iv.mpf(1.0))
    assert bool(enc.a <= true.a) and bool(true.b <= enc.b)
    assert float(enc.delta.b) < 1e-20
    assert not _contains(enc, mp.e + mp.mpf("1e-12"))
    res2 = vti.integrate_c1_qr(
        iv, _jet_harmonic, _var_harmonic, [_pt(1.0), _pt(0.0)], 1.0, 0.0, n_steps=8, order=14
    )
    assert res2["validated"] is True
    y1, y2 = res2["final"]
    c1, s1 = iv.cos(iv.mpf(1.0)), -iv.sin(iv.mpf(1.0))
    assert bool(y1.a <= c1.a) and bool(c1.b <= y1.b)
    assert bool(y2.a <= s1.a) and bool(s1.b <= y2.b)


def test_stage4_plain_lc_anchor_intact() -> None:
    """The Stage-4 plain-``integrate`` Levi-Civita collision anchor is unbroken."""
    iv = mp.iv
    a, mu = 1.0, 0.5
    h = -mu / a
    omega = mp.sqrt(mp.mpf(-h) / 2)
    y0 = [_pt(mp.sqrt(mp.mpf(a))), _pt(0.0), _pt(0.0), _pt(0.0), _pt(0.0)]
    tau_f = 4.0
    res = vti.integrate(iv, _kepler_lc_jet(h), y0, tau_f, mu, n_steps=32, order=16)
    assert res["validated"] is True
    u_enc = res["final"][0]
    u_cf = mp.sqrt(mp.mpf(a)) * mp.cos(omega * mp.mpf(tau_f))
    assert _contains(u_enc, u_cf)
    assert float(u_enc.delta.b) < 1e-25


# ---------------------------------------------------------------------------- #
# 5. The Stage-5 result: QR-regularized crosses the perijove and beats naive.    #
# ---------------------------------------------------------------------------- #
# #837: 3x the suite-default 600s cap. This certificate replay measures ~79s
# serial on a quiet machine, but is single-threaded CPU-bound mpmath, so heavy
# external load (sibling sessions; #810 measured load avg 25 on 8 cores) can
# inflate WALL clock >6.5x past 600s on an otherwise-fine computation. The
# replay depth (n_steps=96, order=10, tau=12) is what the tightness assertions
# certify, so it must not be downgraded; a serial/xdist marker would not help
# (#810 saw the timeout in a -n0 isolation re-run under the same load). The
# cap stays finite as a belt-and-braces guard against genuine hangs.
@pytest.mark.timeout(1800)
def test_qr_regularized_oterma_defeats_wrapping_and_conserves_jacobi() -> None:
    """QR + regularization compose on the real Oterma arc: wrapping is defeated.

    At a matched grid (``h = 0.125``, order 10) both integrators run to ``tau_f = 12``
    -- into the regime where ``#670``'s naive C0 box has already ballooned by the
    ``e^t`` wrapping growth (``hw ~ 1e-8``).  The QR-regularized enclosure instead
    stays tight (``hw ~ 1e-14``), rigorously conserves the Jacobi constant, and is
    >=1e4x tighter than the naive box -- the wrapping ``#670`` diagnosed in the
    regularized frame is genuinely removed once QR operates on the regularized
    coordinates' own local stretching directions.  (The full perijove crossing and
    the ``tau = 80 / t ~ 2.05`` reach are the driver certificate's headline; a unit
    test keeps ``tau`` modest for runtime.)
    """
    iv = mp.iv
    y0, c0 = _oterma_ic(iv)
    jet = vti.make_cr3bp_lc_secondary_jet(H_OTERMA)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_OTERMA)
    qr = vti.integrate_c1_qr(iv, jet, vjet, y0, 12.0, MU_OTERMA, n_steps=96, order=10)
    naive = vti.integrate(iv, jet, y0, 12.0, MU_OTERMA, n_steps=96, order=10)
    assert qr["validated"] is True and naive["validated"] is True
    box = qr["final"]
    t_phys = box[4]
    assert bool(t_phys.a > 0.30)  # a genuine, non-trivial rigorously-enclosed arc
    phys = vti.lc_secondary_to_physical(iv, box, MU_OTERMA)
    cf = vti.cr3bp_planar_jacobi(iv, phys, MU_OTERMA)
    assert bool(cf.a <= c0.a) and bool(c0.b <= cf.b)  # first integral conserved
    hw = qr["final_max_halfwidth"]
    assert hw is not None and 0.0 < hw < 1e-6  # non-vacuous, tight
    # QR is >=1e4x tighter than the (exponentially-wrapping) naive box at tau=12
    assert qr["final_max_halfwidth"] * 1e4 < naive["final_max_halfwidth"]
