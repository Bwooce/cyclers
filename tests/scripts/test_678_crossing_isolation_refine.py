"""#678 -- controls for the TWO-SIDED BRACKET REFINEMENT in ``isolate_section_crossing``.

Stage 12 of the Wilczak-Zgliczynski proof-machinery build (``#668``-``#677``).  #677's
nu5 N=3 box return reached the true crossing (tau~24.1) without the enclosure blowing
up, yet ``isolate_section_crossing`` could not certify a unique transversal crossing.
The #678 diagnosis (instrumented replay) showed WHY: it is NOT a within-step precision
issue -- ``dsigma/ds`` over the step cleanly excludes 0 (transversal) -- it is that the
box has inflated ~157x over the long tau~24 arc, so the crossing is SMEARED across ~5
integration steps; by the time the enclosure is strictly signed on the far side the
step-start box already straddles sigma=0, so no single step shows a strict endpoint
sign change.

This stage adds a sound, backward-compatible TWO-SIDED bracket refinement: when the
coarse whole-step ``dsigma/ds`` interval WRAPS through 0 (interval dependency over the
step) even though every IC crosses transversally within the step, the bracket is
bisected -- always keeping the half whose endpoints still show a *certified strict*
sign change -- until ``dsigma/ds`` is single-signed over the narrower bracket.  It is a
strict generalisation (a step already transversal over ``[0,h]`` is unchanged, so every
``#672``-``#677`` control is untouched), and it cannot admit a false crossing: a midpoint
whose ``sigma`` straddles 0 (genuine box-width smear) or a never-transversal bracket
(genuine tangency) both return ``None``.

These controls drive ``isolate_section_crossing`` directly on hand-built Taylor models
(state = a scalar polynomial in the step variable ``s``, ``sigma(state)=state[0]``) with
a KNOWN true crossing, covering: (A) a genuine transversal crossing whose whole-step
derivative interval wraps through 0 -- the refinement recovers it and brackets the true
root; (B) a near-degenerate tangency (``p'(root)=0``) -- correctly refused; (C) a
box-width-smeared crossing (a synthetic mirror of the real nu5 N=3 failure) -- correctly
refused.  Soundness of every accepted bracket is checked against the closed-form root.
"""

from __future__ import annotations

from typing import Any

import pytest

mp = pytest.importorskip(
    "mpmath",
    reason="mpmath is an optional 'interval' extra (task #610/#625/#668-#678)",
)

import scripts._validated_taylor_integrator as vti  # noqa: E402


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 40
    mp.iv.dps = 40


def _sig_first(iv: Any, st: list[Any]) -> Any:
    return st[0]


def _grad_first(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(1)]


def _series(iv: Any, coeffs: list[Any]) -> list[list[Any]]:
    """A one-component state Taylor model: state[0](s) = sum_k coeffs[k] s^k."""
    return [[c for c in coeffs]]


# --------------------------------------------------------------------------- #
# A. POSITIVE: whole-step dsigma/ds WRAPS through 0, but the refinement recovers #
#    the unique transversal crossing and brackets the true root.                #
#    p(s) = s^3 - 1.5 s^2 + s - 0.2 on [0,1];  p'(s) = 3s^2 - 3s + 1 > 0 always  #
#    (discriminant 9-12<0), yet its interval Horner eval over [0,1] is [-2,1].   #
# --------------------------------------------------------------------------- #
def test_refinement_recovers_transversal_crossing_when_fullstep_wraps() -> None:
    iv = mp.iv
    h = 1.0
    coeffs = [iv.mpf("-0.2"), iv.mpf(1), iv.mpf("-1.5"), iv.mpf(1)]
    state_tm = _series(iv, coeffs)

    # Document the wrapping: the coarse whole-step derivative interval contains 0,
    # so the pre-refinement (full-step-only) transversality test would have failed.
    dstate = [vti.deriv_series(iv, c) for c in state_tm]
    dfull = vti.eval_series_horner(iv, dstate[0], iv.mpf([0.0, h]))
    assert bool(dfull.a < 0) and bool(dfull.b > 0), "control must actually wrap through 0"

    s_star = vti.isolate_section_crossing(iv, state_tm, _sig_first, _grad_first, h)
    assert s_star is not None, "refinement should isolate the genuine transversal crossing"

    # closed-form root of p on [0,1]
    root = mp.findroot(lambda s: s**3 - mp.mpf("1.5") * s**2 + s - mp.mpf("0.2"), mp.mpf("0.3"))
    assert bool(s_star.a <= root) and bool(root <= s_star.b), "bracket must contain true root"
    assert bool(s_star.a >= 0) and bool(s_star.b <= h), "bracket stays within [0,h]"
    # non-vacuous: excludes a clearly-wrong location
    assert not (bool(s_star.a <= mp.mpf("0.8")) and bool(mp.mpf("0.8") <= s_star.b))
    # sigma at the bracketed root really is (rigorously) ~0
    sig = vti.eval_series_horner(iv, state_tm[0], s_star)
    assert bool(sig.a <= 0) and bool(sig.b >= 0)


# --------------------------------------------------------------------------- #
# B. NEGATIVE (near-degenerate tangency): p(s) = (s-0.4)^3 on [0,0.8].          #
#    Strict endpoint sign change (p(0)=-0.064<0, p(0.8)=+0.064>0) but p'(0.4)=0  #
#    -- a NON-transversal inflection crossing.  Must return None (no false pos). #
# --------------------------------------------------------------------------- #
def test_refinement_refuses_nontransversal_tangency() -> None:
    iv = mp.iv
    h = 0.8
    # (s-0.4)^3 = s^3 - 1.2 s^2 + 0.48 s - 0.064
    coeffs = [iv.mpf("-0.064"), iv.mpf("0.48"), iv.mpf("-1.2"), iv.mpf(1)]
    state_tm = _series(iv, coeffs)
    s0 = vti.eval_series_horner(iv, state_tm[0], iv.mpf([0.0, 0.0]))
    sh = vti.eval_series_horner(iv, state_tm[0], iv.mpf([h, h]))
    assert bool(s0.b < 0) and bool(sh.a > 0), "control must have a strict endpoint sign change"
    s_star = vti.isolate_section_crossing(iv, state_tm, _sig_first, _grad_first, h)
    assert s_star is None, "a non-transversal (tangency) crossing must not be certified"


# --------------------------------------------------------------------------- #
# C. NEGATIVE (box-width smear -- the synthetic mirror of the real nu5 N=3      #
#    failure): the same wrapping cubic but with a thick constant-term interval  #
#    (+-0.05 box uncertainty).  Strict endpoints, dsigma/ds wraps so refinement  #
#    engages, but sigma at the interior crossing STRADDLES 0 over the box -- the  #
#    crossing is not localizable within the step.  Must return None.            #
# --------------------------------------------------------------------------- #
def test_refinement_refuses_box_width_smeared_crossing() -> None:
    iv = mp.iv
    h = 1.0
    w = mp.mpf("0.05")
    c0 = iv.mpf([-mp.mpf("0.2") - w, -mp.mpf("0.2") + w])
    coeffs = [c0, iv.mpf(1), iv.mpf("-1.5"), iv.mpf(1)]
    state_tm = _series(iv, coeffs)
    s0 = vti.eval_series_horner(iv, state_tm[0], iv.mpf([0.0, 0.0]))
    sh = vti.eval_series_horner(iv, state_tm[0], iv.mpf([h, h]))
    assert bool(s0.b < 0) and bool(sh.a > 0), "endpoints still strictly straddle"
    dfull = vti.eval_series_horner(iv, vti.deriv_series(iv, state_tm[0]), iv.mpf([0.0, h]))
    assert bool(dfull.a < 0) and bool(dfull.b > 0), "dsigma/ds still wraps (refinement engages)"
    s_star = vti.isolate_section_crossing(iv, state_tm, _sig_first, _grad_first, h)
    assert s_star is None, "a box-width-smeared crossing must not be certified (mirror of nu5 N=3)"


# --------------------------------------------------------------------------- #
# D. BACKWARD COMPATIBILITY: a cleanly transversal whole step is isolated       #
#    exactly as before (refinement breaks immediately with [0,h]).              #
# --------------------------------------------------------------------------- #
def test_clean_transversal_step_isolates_as_before() -> None:
    iv = mp.iv
    h = 1.0
    # p(s) = s - 0.5, p' = 1 > 0 over the whole step: no refinement needed.
    coeffs = [iv.mpf("-0.5"), iv.mpf(1)]
    state_tm = _series(iv, coeffs)
    dfull = vti.eval_series_horner(iv, vti.deriv_series(iv, state_tm[0]), iv.mpf([0.0, h]))
    assert bool(dfull.a > 0), "control is transversal over the whole step already"
    s_star = vti.isolate_section_crossing(iv, state_tm, _sig_first, _grad_first, h)
    assert s_star is not None
    assert bool(s_star.a <= mp.mpf("0.5")) and bool(mp.mpf("0.5") <= s_star.b)
    assert float(s_star.delta.b) < 1e-25, "clean transversal crossing pinned tightly"
