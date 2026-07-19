"""#668 -- tests for the validated (rigorous) interval Taylor ODE integrator.

Covers ``scripts/_validated_taylor_integrator.py``, the rigorous CR3BP flow
enclosure primitive built for the Wilczak-Zgliczynski proof-machinery task.

The load-bearing property under test is *rigor*: every enclosure must provably
CONTAIN the true solution, and -- crucially, per
[[feedback_verify_gauntlet_with_positive_control]] -- must not be a vacuous
"contains everything" bound.  So the closed-form goldens (exp, harmonic
oscillator) assert BOTH containment of the true value AND exclusion of nearby
wrong values, and the integrator must honestly return ``None`` (not a silently
wrong enclosure) when a step is too large to validate.

Note on exactness: a validated step encloses the solution at the *dyadic* time
actually reached (the Python-float step size), so the goldens use exactly
representable steps (1/8) and reach an exactly representable time (t=1) so the
comparison target ``exp(1)`` / ``cos(1)`` is the true value at the reached time.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath", reason="mpmath is an optional 'interval' extra (task #610/#625/#668)"
)

import scripts._validated_taylor_integrator as vti  # noqa: E402

MU_OTERMA = 0.0009537
C_OTERMA = 3.03


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 50
    mp.iv.dps = 50


def _pt(x: float) -> Any:
    return mp.iv.mpf([x, x])


def _contains(enc: Any, value: Any) -> bool:
    """True iff the interval ``enc`` rigorously contains the point ``value``."""
    return bool(enc.a <= value) and bool(value <= enc.b)


def _jet_exp(iv: Any, state: list[Any], mu: float) -> list[Any]:
    return [state[0]]  # y' = y


def _jet_harmonic(iv: Any, state: list[Any], mu: float) -> list[Any]:
    y1, y2 = state
    return [y2, vti.ts_scale(iv, y1, iv.mpf(-1.0))]  # y1'=y2, y2'=-y1


def test_exp_enclosure_contains_and_excludes() -> None:
    """y'=y from 1 over t=1 (8 steps of 1/8) rigorously encloses e, tightly."""
    iv = mp.iv
    res = vti.integrate(iv, _jet_exp, [_pt(1.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    enc = res["final"][0]
    true = iv.exp(iv.mpf(1.0))  # e at the reached (exact) time t=1
    # (1) rigorous containment
    assert bool(enc.a <= true.a) and bool(true.b <= enc.b)
    # (2) NOT vacuous: the enclosure is tight and excludes clearly-wrong values
    assert float(enc.delta.b) < 1e-20
    assert not _contains(enc, mp.e + mp.mpf("1e-12"))
    assert not _contains(enc, mp.e - mp.mpf("1e-12"))


def test_harmonic_oscillator_enclosure() -> None:
    """Rotation: encloses (cos 1, -sin 1) at t=1 and conserves energy tightly."""
    iv = mp.iv
    res = vti.integrate(iv, _jet_harmonic, [_pt(1.0), _pt(0.0)], 1.0, 0.0, n_steps=8, order=14)
    assert res["validated"] is True
    y1, y2 = res["final"]
    c1 = iv.cos(iv.mpf(1.0))
    s1 = -iv.sin(iv.mpf(1.0))
    assert bool(y1.a <= c1.a) and bool(c1.b <= y1.b)
    assert bool(y2.a <= s1.a) and bool(s1.b <= y2.b)
    assert float(y1.delta.b) < 1e-18


def test_apriori_enclosure_satisfies_picard_inclusion() -> None:
    """The returned rough box must actually satisfy the Picard inclusion it claims."""
    iv = mp.iv
    y0 = [_pt(1.0), _pt(0.0)]
    h = 0.1
    W = vti.apriori_enclosure(iv, _jet_harmonic, y0, h, 0.0)
    assert W is not None
    # Re-verify independently: y0 + [0,h]*f(W)  subset of  W.
    hint = iv.mpf([0.0, h])
    order0 = [vti.ts_const(iv, W[i], 0) for i in range(len(W))]
    fW = [d[0] for d in _jet_harmonic(iv, order0, 0.0)]
    for i in range(len(W)):
        img = y0[i] + hint * fW[i]
        assert bool(W[i].a <= img.a) and bool(img.b <= W[i].b)


def test_step_too_large_returns_none() -> None:
    """An unvalidatable step returns None -- an honest failure, not a wrong bound."""
    iv = mp.iv
    # y'=y cannot be validated over a full unit step from 1 (solution ~e).
    assert vti.validated_step(iv, _jet_exp, [_pt(1.0)], 1.0, 0.0, order=12) is None


def test_cr3bp_jacobi_conserved_along_real_oterma_arc() -> None:
    """Rigorous: Jacobi C is enclosed around C0 along a real Oterma heteroclinic arc."""
    iv = mp.iv
    # L1->L2 heteroclinic section point, full C=3.03 state.
    x0, vx0 = 0.9522928423486199945, 1.23e-05
    r1 = iv.sqrt((iv.mpf(x0) + iv.mpf(MU_OTERMA)) ** 2)
    r2 = iv.sqrt((iv.mpf(x0) + iv.mpf(MU_OTERMA - 1.0)) ** 2)
    vy2 = (
        iv.mpf(x0) ** 2
        + 2 * (1 - iv.mpf(MU_OTERMA)) / r1
        + 2 * iv.mpf(MU_OTERMA) / r2
        - iv.mpf(vx0) ** 2
        - iv.mpf(C_OTERMA)
    )
    assert bool(vy2.a > 0)  # energetically admissible on the C=3.03 manifold
    y0 = [_pt(x0), _pt(0.0), _pt(vx0), iv.sqrt(vy2)]
    c0 = vti.cr3bp_planar_jacobi(iv, y0, MU_OTERMA)

    res = vti.integrate(iv, vti.cr3bp_planar_jet, y0, 0.25, MU_OTERMA, n_steps=16, order=12)
    assert res["validated"] is True  # short arc validates cleanly
    cf = vti.cr3bp_planar_jacobi(iv, res["final"], MU_OTERMA)
    # first integral: end-of-arc Jacobi enclosure rigorously contains C0
    assert bool(cf.a <= c0.a) and bool(c0.b <= cf.b)
    # and the arc stayed tight (a real, non-vacuous enclosure)
    assert res["final_max_halfwidth"] is not None
    assert res["final_max_halfwidth"] < 1e-6


def test_wrapping_effect_is_stepsize_independent() -> None:
    """Blow-up horizon barely moves when the step halves -> wrapping, not truncation.

    Documents the honest limitation: the naive C0 enclosure blows up near the L1
    neck at a step-size-independent time; only Lohner QR reframing (future work)
    would extend it.
    """
    iv = mp.iv
    x0, vx0 = 0.9522928423486199945, 1.23e-05
    r1 = iv.sqrt((iv.mpf(x0) + iv.mpf(MU_OTERMA)) ** 2)
    r2 = iv.sqrt((iv.mpf(x0) + iv.mpf(MU_OTERMA - 1.0)) ** 2)
    vy2 = (
        iv.mpf(x0) ** 2
        + 2 * (1 - iv.mpf(MU_OTERMA)) / r1
        + 2 * iv.mpf(MU_OTERMA) / r2
        - iv.mpf(vx0) ** 2
        - iv.mpf(C_OTERMA)
    )
    y0 = [_pt(x0), _pt(0.0), _pt(vx0), iv.sqrt(vy2)]
    t_reached = []
    for n_steps in (64, 128):
        res = vti.integrate(iv, vti.cr3bp_planar_jet, y0, 1.0, MU_OTERMA, n_steps=n_steps, order=10)
        assert res["validated"] is False  # naive C0 cannot reach t=1
        t_reached.append(res["n_completed"] * (1.0 / n_steps))
    # halving the step extends the reachable time by well under 2x (wrapping)
    assert max(t_reached) < 2.0 * min(t_reached)
    assert min(t_reached) > 0.3  # but a genuine, non-trivial arc IS certified
