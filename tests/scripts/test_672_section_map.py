"""#672 -- tests for the rigorous Poincare SECTION-CROSSING MAP + its Jacobian.

Stage 6 of the Wilczak-Zgliczynski proof-machinery build (``#668`` Stages 1-2,
``#669`` Stage 3 QR, ``#670`` Stage 4 Levi-Civita, ``#671`` Stage 5 composition).
The W-Z existence proof is verified through the Poincare SECTION MAP ``P`` on
``Theta = {y = 0}`` (parameterized by ``(x, xdot)``), not the raw flow.  This stage
builds the precursor to h-sets: rigorously ISOLATE the section-crossing time (an
interval-Newton root of ``sigma(state(s)) = 0`` on the crossing step's rigorous
Taylor model -- NOT a floating-point event), evaluate the enclosed crossing state,
and form the rigorous section-map Jacobian

    DP = Dphi_tau*  -  f(P) (Dsigma(P) . Dphi_tau*) / (Dsigma(P) . f(P))

(the flow STM plus the implicit-function correction for the crossing time's IC
dependence).  Building h-sets / proving covering relations remains a future stage.

Load-bearing properties, in the ``#668``-``#671`` containment discipline (contains
the true value AND excludes nearby wrong values, non-vacuous):

  * **the section-map Jacobian formula is exact** -- the pure-linear-algebra
    ``section_map_jacobian`` matches the hand-computed ``DP`` for a rotation flow;
  * **the whole pipeline is correct** -- on the 2D harmonic oscillator (section
    ``{y1=0}``, start ``(1,1)``) the rigorous crossing time (``3 pi/4``), post-map
    state (``(0, -sqrt2)``) AND full ``DP = [[0,0],[-1/sqrt2,-1/sqrt2]]`` all match
    the independent closed form, tightly;
  * **a second, 5D control** -- pure-Kepler-in-Levi-Civita, section ``{u=0}``:
    ``tau* = pi``, ``pu* = -2``, ``Dsigma.f = -1/2``, on-section ``DP`` row vanishes;
  * ``integrate_c1_qr`` now also returns the accumulated flow STM (``Phi = A@B``),
    which encloses the harmonic rotation STM tightly;
  * **the real W-Z Oterma section point** is rigorously isolated: a tight
    crossing-time interval, a transversality certificate, and a tight non-vacuous
    ``DP`` -- the QR-accumulated STM staying tight through the first Jupiter perijove.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#672)",
)

import scripts._validated_taylor_integrator as vti  # noqa: E402

MU_OTERMA = 0.0009537
C_OTERMA = 3.03
H_OTERMA = -C_OTERMA / 2
X0, VX0 = 0.9522928423486199945, 1.23e-05


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _contains(enc: Any, value: Any) -> bool:
    return bool(enc.a <= value) and bool(value <= enc.b)


# --- flows ------------------------------------------------------------------- #
def _harm_jet(iv: Any, s: list[Any], mu: float) -> list[Any]:
    return [s[1], vti.ts_scale(iv, s[0], iv.mpf(-1.0))]


def _harm_var(iv: Any, aug: list[Any], mu: float) -> list[Any]:
    y1, y2 = aug[0], aug[1]
    v00, v01, v10, v11 = aug[2], aug[3], aug[4], aug[5]

    def neg(s: list[Any]) -> list[Any]:
        return vti.ts_scale(iv, s, iv.mpf(-1.0))

    return [y2, neg(y1), v10, v11, neg(v00), neg(v01)]


def _sig_first(iv: Any, st: list[Any]) -> Any:
    return st[0]


def _grad_first_2(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(1), iv.mpf(0)]


def _kep_jet(h: float) -> Any:
    def jet(iv: Any, s: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = s
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        return [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]

    return jet


def _kep_var(h: float) -> Any:
    def jet(iv: Any, aug: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = aug[0], aug[1], aug[2], aug[3], aug[4]
        vmat = aug[5:30]
        order = len(u) - 1
        zero = vti.ts_const(iv, iv.mpf(0), order)
        q = vti.ts_const(iv, iv.mpf(0.25), order)
        th = vti.ts_const(iv, iv.mpf(2.0 * h), order)
        tu = vti.ts_scale(iv, u, iv.mpf(2.0))
        tv = vti.ts_scale(iv, v, iv.mpf(2.0))
        df = [
            [zero, zero, q, zero, zero],
            [zero, zero, zero, q, zero],
            [th, zero, zero, zero, zero],
            [zero, th, zero, zero, zero],
            [tu, tv, zero, zero, zero],
        ]
        dV: list[Any] = [zero for _ in range(25)]
        for c in range(5):
            col = [vmat[r * 5 + c] for r in range(5)]
            for rr in range(5):
                acc = None
                for kk in range(5):
                    if df[rr][kk] is zero:
                        continue
                    t = vti.ts_mul(iv, df[rr][kk], col[kk])
                    acc = t if acc is None else vti.ts_add(iv, acc, t)
                dV[rr * 5 + c] = acc if acc is not None else zero
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        st = [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]
        return [*st, *dV]

    return jet


def _grad_first_5(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(1), iv.mpf(0), iv.mpf(0), iv.mpf(0), iv.mpf(0)]


# ---------------------------------------------------------------------------- #
# 1. The section-map Jacobian FORMULA (pure linear algebra, closed form).       #
# ---------------------------------------------------------------------------- #
def test_section_map_jacobian_matches_hand_computation() -> None:
    """``DP = Phi - f (Dsigma Phi)/(Dsigma f)`` for a quarter-turn rotation flow.

    Harmonic oscillator, section ``{y1=0}``, start ``(1,0)``: crossing at ``pi/2``,
    ``Phi = [[0,1],[-1,0]]``, ``f = [y2,-y1] = [-1,0]``, ``Dsigma = [1,0]``.
    Hand result ``DP = [[0,0],[-1,0]]`` (first row zero on-section; ``d(post y2)/
    d(init y1) = -1``).
    """
    iv = mp.iv
    phi = [[_pt(0.0), _pt(1.0)], [_pt(-1.0), _pt(0.0)]]
    fvec = [_pt(-1.0), _pt(0.0)]
    dsig = [_pt(1.0), _pt(0.0)]
    dp, dsf, _ = vti.section_map_jacobian(iv, phi, fvec, dsig)
    assert _contains(dsf, mp.mpf(-1))  # Dsigma.f = -1 (transversal, nonzero)
    expected = [[0.0, 0.0], [-1.0, 0.0]]
    for i in range(2):
        for j in range(2):
            assert _contains(dp[i][j], mp.mpf(expected[i][j])), (i, j)
            assert float(dp[i][j].delta.b) < 1e-30


# ---------------------------------------------------------------------------- #
# 2. flow_taylor_model rigor: encloses the exact harmonic solution over a step.  #
# ---------------------------------------------------------------------------- #
def test_flow_taylor_model_encloses_solution_over_step() -> None:
    """The step Taylor model encloses ``(cos s, -sin s)`` for all ``s in [0,h]``."""
    iv = mp.iv
    h = 0.3
    ybox = [_pt(1.0), _pt(0.0)]
    tm = vti.flow_taylor_model(iv, _harm_jet, ybox, h, 0.0, 18)
    assert tm is not None
    for s in [0.0, 0.1, 0.2, 0.3]:
        y1 = vti.eval_series_horner(iv, tm[0], _pt(s))
        y2 = vti.eval_series_horner(iv, tm[1], _pt(s))
        assert _contains(y1, mp.cos(mp.mpf(s)))
        assert _contains(y2, -mp.sin(mp.mpf(s)))
        assert float(y1.delta.b) < 1e-16  # tight, non-vacuous (whole-step Taylor model)


# ---------------------------------------------------------------------------- #
# 3. integrate_c1_qr now returns the accumulated flow STM (Phi = A @ B).         #
# ---------------------------------------------------------------------------- #
def test_integrate_c1_qr_returns_accumulated_stm() -> None:
    """The QR-accumulated STM encloses the harmonic rotation ``[[cos,sin],[-sin,cos]]``."""
    iv = mp.iv
    res = vti.integrate_c1_qr(
        iv, _harm_jet, _harm_var, [_pt(1.0), _pt(0.0)], 1.0, 0.0, n_steps=8, order=14
    )
    assert res["validated"] is True
    stm = res["final_stm"]
    assert stm is not None
    c1, s1 = mp.cos(mp.mpf(1.0)), mp.sin(mp.mpf(1.0))
    exp = [[c1, s1], [-s1, c1]]
    for i in range(2):
        for j in range(2):
            assert _contains(stm[i][j], exp[i][j]), (i, j)
            assert float(stm[i][j].delta.b) < 1e-20


# ---------------------------------------------------------------------------- #
# 4. End-to-end closed-form control: 2D harmonic, full DP validation.           #
# ---------------------------------------------------------------------------- #
def test_harmonic_section_map_matches_closed_form() -> None:
    """Crossing time, post-map state, AND the full ``DP`` all match the closed form.

    Start ``(1,1)``: ``tau* = 3 pi/4``, post state ``(0, -sqrt2)``,
    ``DP = [[0,0],[-1/sqrt2,-1/sqrt2]]``, ``Dsigma.f = -sqrt2``.
    """
    iv = mp.iv
    res = vti.rigorous_section_map(
        iv,
        _harm_jet,
        _harm_var,
        [_pt(1.0), _pt(1.0)],
        0.0,
        sigma_val=_sig_first,
        sigma_grad=_grad_first_2,
        tau_max=4.0,
        n_steps=64,
        order=16,
        tau_min=0.1,
    )
    assert res["found"] is True
    assert res["transversal"] is True
    ts = res["tau_star"]
    assert _contains(ts, 3 * mp.pi / 4)
    assert float(ts.delta.b) < 1e-20
    cs = res["crossing_state"]
    assert float(cs[0].delta.b) < 1e-20  # y1* = 0
    assert _contains(cs[1], -mp.sqrt(2))
    assert _contains(res["dsigma_dot_f"], -mp.sqrt(2))
    dp = res["section_jacobian"]
    off = -1 / mp.sqrt(2)
    assert _contains(dp[0][0], mp.mpf(0)) and _contains(dp[0][1], mp.mpf(0))
    assert _contains(dp[1][0], off) and _contains(dp[1][1], off)
    assert float(dp[1][0].delta.b) < 1e-20  # tight, non-vacuous
    assert not _contains(dp[1][0], off + mp.mpf("1e-12"))  # excludes nearby wrong
    assert not _contains(dp[1][1], off - mp.mpf("1e-12"))


# ---------------------------------------------------------------------------- #
# 5. Second control: pure-Kepler-in-Levi-Civita (5D regularized), section {u=0}. #
# ---------------------------------------------------------------------------- #
def test_kepler_lc_section_map_matches_closed_form() -> None:
    """``u(tau)=cos(omega tau)`` crosses 0 at ``tau=pi``; ``pu*=-2``; ``Dsigma.f=-1/2``."""
    iv = mp.iv
    h = -0.5
    y0 = [_pt(1.0), _pt(0.0), _pt(0.0), _pt(0.0), _pt(0.0)]
    res = vti.rigorous_section_map(
        iv,
        _kep_jet(h),
        _kep_var(h),
        y0,
        0.0,
        sigma_val=_sig_first,
        sigma_grad=_grad_first_5,
        tau_max=5.0,
        n_steps=64,
        order=16,
        tau_min=0.5,
    )
    assert res["found"] is True and res["transversal"] is True
    assert _contains(res["tau_star"], mp.pi)
    cs = res["crossing_state"]
    assert float(cs[0].delta.b) < 1e-20  # u* = 0
    assert _contains(cs[2], mp.mpf(-2))  # pu* = -2
    assert _contains(res["dsigma_dot_f"], mp.mpf("-0.5"))
    dp = res["section_jacobian"]
    for j in range(5):  # on-section: the u-row of DP vanishes
        assert _contains(dp[0][j], mp.mpf(0))
        assert float(dp[0][j].delta.b) < 1e-20


# ---------------------------------------------------------------------------- #
# 6. The real W-Z Oterma section point: rigorous, transversal, tight DP.         #
# ---------------------------------------------------------------------------- #
def _oterma_ic(iv: Any) -> list[Any]:
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
    return vti.lc_secondary_from_physical(iv, _pt(X0), _pt(0.0), _pt(VX0), vy0, MU_OTERMA)


def _sig_2uv(iv: Any, st: list[Any]) -> Any:
    return iv.mpf(2) * st[0] * st[1]


def _grad_2uv(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(2) * st[1], iv.mpf(2) * st[0], iv.mpf(0), iv.mpf(0), iv.mpf(0)]


# #837: 3x the suite-default 600s cap. This certificate replay measures ~92s
# serial on a quiet machine, but is single-threaded CPU-bound mpmath, so heavy
# external load (sibling sessions; #810 measured load avg 25 on 8 cores) can
# inflate WALL clock >6.5x past 600s on an otherwise-fine computation. The
# replay depth (n_steps=136, order=10, tau_max=17) is what the tightness
# assertions certify, so it must not be downgraded; a serial/xdist marker would
# not help (#810 saw the timeout in a -n0 isolation re-run under the same
# load). The cap stays finite as a belt-and-braces guard against genuine hangs.
@pytest.mark.timeout(1800)
def test_oterma_first_section_crossing_rigorously_isolated() -> None:
    """The first forward re-crossing of ``{y=0}`` is rigorously isolated + transversal.

    From the first published golden crossing (on ``{y=0}``), the machinery isolates a
    tight crossing-time interval, encloses the physical ``(x, xdot)``, certifies
    transversality, and returns a tight non-vacuous ``DP`` -- the QR-accumulated STM
    staying tight straight through the first Jupiter perijove (``T ~ 0.462``) the
    crossing sits just past.  (A rounded IC does not track the unstable connection,
    so this crossing is not expected to equal golden point 2; matching the published
    sequence is the future h-set / covering-relations stage.)
    """
    iv = mp.iv
    w0 = _oterma_ic(iv)
    jet = vti.make_cr3bp_lc_secondary_jet(H_OTERMA)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_OTERMA)
    res = vti.rigorous_section_map(
        iv,
        jet,
        vjet,
        w0,
        MU_OTERMA,
        sigma_val=_sig_2uv,
        sigma_grad=_grad_2uv,
        tau_max=17.0,
        n_steps=136,
        order=10,
        tau_min=1.0,
    )
    assert res["found"] is True
    assert res["transversal"] is True
    # 0 not in the transversality scalar (Dsigma.f) -- a genuine transversal crossing
    dsf = res["dsigma_dot_f"]
    assert bool(dsf.a > 0) or bool(dsf.b < 0)
    ts = res["tau_star"]
    assert bool(ts.a > 1.0)  # a genuine forward crossing, not the trivial start
    assert float(ts.delta.b) < 1e-6  # crossing time rigorously pinned down
    cs = res["crossing_state"]
    phys = vti.lc_secondary_to_physical(iv, cs, MU_OTERMA)
    # physical y is (rigorously) zero at the crossing
    assert _contains(phys[1], mp.mpf(0))
    assert float(phys[1].delta.b) < 1e-10
    # physical x is a tight, non-vacuous interval near Jupiter (x ~ 1.0035)
    assert float(phys[0].delta.b) < 1e-6
    assert bool(phys[0].a > 0.9) and bool(phys[0].b < 1.1)
    # the rigorous section-map Jacobian is non-vacuous and tight
    dp = res["section_jacobian"]
    dp_hw = max(float(dp[i][j].delta.b) for i in range(5) for j in range(5)) / 2
    assert 0.0 < dp_hw < 1e-6
    # the crossing lands just past the first Jupiter perijove (physical time T ~ 0.46)
    assert bool(cs[4].a > 0.4) and bool(cs[4].b < 0.5)
